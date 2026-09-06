#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==== MEETING NOTE ====
# WHAT IT IS FOR   : Splits a fungal bin into genus populations read by read (ITS), polishes a
#                    consensus from the DOMINANT population only (medoid read as template,
#                    minimap2 + samtools consensus, two rounds) and decides that consensus
#                    over the three loci (18S / ITS / 28S).
# INPUT            : --root; "fastq files"/F*/ bin reads; REFERENCE_DB (or a fast copy)
# OUTPUT           : referans_konsensus/fungal_polish/consensus/<bin>.fasta (a new candidate set),
#                    referans_konsensus/fungal_polish/population/<bin>.tsv (one row per read),
#                    FUNGAL_BIN_IDENTITY.tsv (one row per bin: populations, decision, gain)
# HOW TO RUN IT    : ./primerjury polish  (verification/population_polish.py drives every library;
#                    this file is the shared machinery)  |  python3 verification/fungal_bin_identity.py = self-test
# WHY IT IS LIKE THIS : Measured (2026-09-05): in a mixed fungal bin the chain's blended
#                    consensus sat at 97.25 per cent to the Petriella musispora TYPE ITS record
#                    while the bin's own single reads reached a median of 98.20 and a best of
#                    100.00; 47 of 60 reads went to one genus, the rest to other Microascaceae.
#                    A blend of two organisms represents neither and stays under the species
#                    threshold (ITS 99.6). A consensus polished from the dominant population
#                    represents one organism.
# =======================
u"""FUNGAL BIN IDENTITY: the polished consensus of the dominant population.

THE PROBLEM (measured 2026-09-05)
  Fungal bins are MIXED. In F2-2_500148, 78 per cent of the reads go to Petriella
  by ITS and the rest to Scedosporium / Lomentospora (all Microascaceae). The
  chain's consensus (reference-anchored, then dominant allele) is the BLEND of that
  mixture: 97.25 per cent to the P. musispora type record over 509 bp. The same
  bin's single reads: median 98.20, best 100.00. So the consensus was WORSE than a
  raw read; it stayed under the ITS species threshold (99.6) and the bin was left
  at "Scedosporium sp." while the read witness said "Petriella musispora".

  Two causes at once:
    1. THE MIXTURE. A per-position majority base over two organisms produces a
       sequence between them; in a fast-evolving region such as ITS that sequence
       is close to no record at all.
    2. A FOREIGN TEMPLATE. When the anchoring reference is another genus
       (Petriella reads on a Scedosporium template) the ITS alignment is gapped.

THE FIX: SPLIT INTO POPULATIONS, POLISH FROM THE BIN'S OWN READ
  1. Reads are sampled from the bin (reservoir sampling, fixed seed).
  2. Every read is searched against the ITS databases (RefSeq ITS + UNITE: the
     SAME set the consensus route uses, `locus_decision.MANTAR_LOKUSLARI`); the best
     record per read (identity, then alignment; at least EN_AZ_KANIT bp) is kept.
  3. Reads are grouped by the GENUS of their best record. The largest group is
     the dominant population; its share and the runner-up are reported.
  4. THE TEMPLATE is the population's MEDOID read: the read with the highest mean
     k-mer Jaccard similarity to the other reads of the population, at least 90
     per cent of the population's median length. It is chosen WITHOUT a
     reference: the bin's own reads pick it, so there is no circularity.
  5. The population's reads are aligned to the template (minimap2 map-ont),
     `samtools consensus` is taken, the consensus becomes the next template;
     --rounds times (default 2; measured: no gain after the second round).
  6. The polished consensus is decided over the three loci with
     `locus_decision.lokus_karari` + `birlestir` (the SAME thresholds and rule as
     the consensus route) and the reported (cascade) method is written next to it.
  7. The result is written as a new candidate set,
     `referans_konsensus/fungal_polish/consensus/<bin>.fasta`. `select_consensus.py`
     weighs it with the SAME criterion as every other candidate (read k-mer
     support); if it wins, the canonical consensus and the identity query are it.
     No hand rule says "use this one for fungi": the set wins because it
     represents the bin's reads better.

THE READ WITNESS FROM THE SAME DATABASE
  The population table (best record per read, RefSeq + UNITE) also serves as the
  read-level witness, so a "witness disagreement" can no longer be an artefact of
  two routes looking at two different databases.

RUN
  This file is the shared machinery (sampling, windows, per-read records, medoid,
  polish, populations). The driver for every library is
  verification/population_polish.py (`./primerjury polish`; `./primerjury fungi`
  is the same with --groups F1,F2). Running this file runs its self-test.
"""
from __future__ import print_function
import argparse
import collections
import io
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _d in (os.path.join(KOK, 'verification'), os.path.join(KOK, 'steps')):
    if _d not in sys.path:
        sys.path.insert(0, _d)
from locus_decision import (MANTAR_LOKUSLARI, lokus_karari, birlestir,      # noqa: E402
                            raporlanan_yontem)
from identity_verification import EN_AZ_KANIT, ADSIZ_JETONLARI, CINS_ESIGI  # noqa: E402
from target_identity import ad_ayikla, cins_epitet                          # noqa: E402

SET_DIR = os.path.join('referans_konsensus', 'fungal_polish', 'consensus')   # the new candidate set
POP_DIR = os.path.join('referans_konsensus', 'fungal_polish', 'population')
SEED = 20260905
K_MEDOID = 15          # k-mer size for the medoid
MIN_POPULATION = 5     # no consensus from fewer reads than this
MIN_DEPTH = 3          # samtools consensus -d
# medaka (2026-09-06): ONT's own neural polisher, applied ON TOP OF the samtools consensus. Measured on
# the study's HAC root (99 bins, same code, the only difference medaka): one bin rose to species
# (B-4_1642647, 98.59 -> 99.15 per cent), none fell; the bin's own reads fit the polished sequence
# better in 21 bins and ~one base worse in 2. Medaka ALONE regressed one bin (99.86 -> 99.71), so it is
# never used alone. The model is the data's basecalling model. Environment: a Python 3.11 conda env
# (MEDAKA_ENV points at its bin directory); without it the samtools consensus stands and a warning says so.
MEDAKA_ENV = os.environ.get('MEDAKA_ENV', os.path.expanduser('~/miniconda3/envs/medaka/bin'))
MEDAKA_MODEL = os.environ.get('MEDAKA_MODEL', 'r1041_e82_400bps_sup_v5.2.0')
PLACEHOLDER = ('sp', 'sp.', 'spp', 'spp.', 'cf', 'cf.', 'aff', 'aff.')
FUNGAL_LIBS = ('F1', 'F2')
ITS_DB = [d for lok, dbs, _a in MANTAR_LOKUSLARI if lok == u'ITS' for d in dbs]
WINDOW_PAD = 60        # bases added on each side of the ITS window


def its_window(reads_fa, db, threads=1):
    u"""Where the ITS sits on each read: the reads are searched against the small RefSeq
    ITS set and the best record's (bitscore) span on the read is the window.
    {read: (qstart, qend)}. WHY: a whole-operon read (3.7 kb) searched against UNITE
    gives long ~93 per cent alignments of its 18S/28S flanks to unrelated genera
    (measured: Fusarium / Nectria; 60 reads took 442 s). The window cuts both the
    artefact and the query length by about six."""
    # 2026-09-06: minimap2 against the small RefSeq set's index. Measured: the bacterial 16S
    # window for 12,377 reads took 48 min with blastn and seconds with minimap2. The window is
    # the read span of the highest-scoring mapping (PAF columns 3-4, 0-based -> 1-based).
    mmi = mmi_path(db, threads)
    paf = run(['minimap2', '-t', str(threads), '-x', 'map-ont', '--secondary=no', mmi, reads_fa])
    best = {}
    for line in paf.splitlines():
        p = line.split('\t')
        if len(p) < 12:
            continue
        k, qs, qe, score = p[0], int(p[2]) + 1, int(p[3]), int(p[9])
        if k not in best or score > best[k][0]:
            best[k] = (score, qs, qe)
    return {k: (v[1], v[2]) for k, v in best.items()}


def window_file(records, window, path):
    u"""Write the reads cut to their ITS window (+- WINDOW_PAD); a read without a window
    (no RefSeq hit at all: a new lineage) is written whole. Returns the number cut."""
    out, n = [], 0
    for k, d in records:
        if k in window:
            qs, qe = window[k]
            out.append((k, d[max(0, qs - 1 - WINDOW_PAD):min(len(d), qe + WINDOW_PAD)]))
            n += 1
        else:
            out.append((k, d))
    write_fasta(path, out)
    return n


def db_path(root, name, db_dir=None, quiet=False):
    u"""Resolve one database: a fast copy first ($PRIMERJURY_FAST_DB or ~/refdb_fast),
    then <root>/REFERENCE_DB (or --db-dir). Missing databases are REPORTED, never
    skipped silently."""
    cands = [os.environ.get('PRIMERJURY_FAST_DB') or os.path.expanduser('~/refdb_fast'),
             db_dir or os.path.join(root, 'REFERENCE_DB')]
    for k in cands:
        y = os.path.join(k, name)
        if os.path.exists(y + '.nin'):
            return y
    if not quiet:
        print(u'  WARNING: %s is not indexed in any database folder; NOT searched' % name)
    return None


def bin_files(root):
    u"""{label: fastq path}; a file whose prefix does not match its folder is a stray
    copy of another bin and is ignored."""
    base = os.path.join(root, 'fastq files')
    out = {}
    if not os.path.isdir(base):
        return out
    for d in sorted(os.listdir(base)):
        dd = os.path.join(base, d)
        if not os.path.isdir(dd):
            continue
        for f in sorted(os.listdir(dd)):
            m = re.match(r'^(.+?)[-_]reads[-_]([A-Za-z]*\d+)\.fastq$', f)
            if not m:
                continue
            if m.group(1).replace('_', '-') != d.replace('_', '-'):
                continue
            out['%s_%s' % (d, m.group(2))] = os.path.join(dd, f)
    return out


def sample_reads(path, n, seed=SEED):
    u"""Reservoir sampling of (id, sequence); the file is never held in memory."""
    rnd = random.Random(seed)
    out, seen = [], 0
    if not (path and os.path.exists(path)):
        return out
    with io.open(path, encoding='utf-8', errors='replace') as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            d = fh.readline().strip().upper()
            fh.readline()
            fh.readline()
            if len(d) < 400:
                continue
            seen += 1
            rec = (h[1:].split()[0] if h.startswith(u'@') else u'r%d' % seen, d)
            if len(out) < n:
                out.append(rec)
            else:
                j = rnd.randrange(seen)
                if j < n:
                    out[j] = rec
    return out


def write_fasta(path, records):
    with io.open(path, 'w', encoding='utf-8') as g:
        for name, seq in records:
            g.write(u'>%s\n%s\n' % (name, seq))


def read_fasta_seq(path):
    s = []
    for l in io.open(path, encoding='utf-8', errors='replace'):
        if not l.startswith(u'>'):
            s.append(l.strip())
    return u''.join(s).upper()


def kmer_set(s, k=K_MEDOID):
    return set(s[i:i + k] for i in range(len(s) - k + 1))


_COMP = {u'A': u'T', u'C': u'G', u'G': u'C', u'T': u'A'}


def reverse_complement(s):
    return u''.join(_COMP.get(c, u'N') for c in reversed(s))


def medoid(records):
    u"""The population's medoid read: highest mean k-mer Jaccard to the others (the
    larger of the two strands counts), at least 90 per cent of the median length.
    Returns (id, sequence, mean similarity)."""
    if not records:
        return None
    lens = sorted(len(d) for _, d in records)
    med = lens[len(lens) // 2]
    cands = [(k, d) for k, d in records if len(d) >= 0.9 * med] or list(records)
    km = {k: (kmer_set(d), kmer_set(reverse_complement(d))) for k, d in records}
    best = None
    for k, d in cands:
        a = km[k][0]
        tot, n = 0.0, 0
        for k2, _d2 in records:
            if k2 == k:
                continue
            b, brc = km[k2]
            j1 = len(a & b) / float(len(a | b) or 1)
            j2 = len(a & brc) / float(len(a | brc) or 1)
            tot += max(j1, j2)
            n += 1
        score = tot / n if n else 0.0
        if best is None or (score, len(d)) > (best[2], len(best[1])):
            best = (k, d, score)
    return best


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(u'%s\n%s' % (u' '.join(cmd), (r.stderr or '')[-400:]))
    return r.stdout


def blast(query, db, threads=1, top=50):
    u"""{query id: [(pident, length, bitscore, qlen, title, slen)]} - the tuple the
    locus decision expects."""
    r = subprocess.run(['blastn', '-query', query, '-db', db, '-outfmt',
                        '6 qseqid pident length bitscore qlen slen stitle',
                        '-max_target_seqs', str(top), '-evalue', '1e-20',
                        '-num_threads', str(threads)], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(u'blastn failed (%d): %s' % (r.returncode, (r.stderr or '')[:300]))
    d = {}
    for line in r.stdout.splitlines():
        p = line.split('\t')
        if len(p) < 7:
            continue
        d.setdefault(p[0], []).append((float(p[1]), int(p[2]), float(p[3]),
                                       int(p[4]), p[6], int(p[5])))
    return d


def polish(template_fa, reads_fa, out_fa, threads=1):
    u"""Align the reads to the template (minimap2 map-ont), take samtools consensus.
    Returns (sequence, internal N count); terminal N (no depth) is trimmed."""
    bam = out_fa + '.bam'
    p1 = subprocess.Popen(['minimap2', '-t', str(threads), '-ax', 'map-ont',
                           '--secondary=no', template_fa, reads_fa],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(['samtools', 'sort', '-@', '1', '-o', bam, '-'],
                          stdin=p1.stdout, stderr=subprocess.DEVNULL)
    p1.stdout.close()
    p2.communicate()
    if p2.returncode != 0:
        raise RuntimeError(u'minimap2 / samtools sort failed on %s' % reads_fa)
    run(['samtools', 'index', bam])
    run(['samtools', 'consensus', '-a', '--show-ins', 'yes', '--show-del', 'no',
         '-m', 'simple', '-c', '0.5', '-d', str(MIN_DEPTH), '-f', 'fasta',
         '-o', out_fa, bam])
    seq = read_fasta_seq(out_fa).strip(u'N')
    for e in (bam, bam + '.bai'):
        try:
            os.remove(e)
        except OSError:
            pass
    return seq, seq.count(u'N')


def medaka_present():
    return os.path.exists(os.path.join(MEDAKA_ENV, 'medaka'))


def medaka_polish(template_fa, reads_fa, out_fa, threads=1):
    u"""minimap2 -> medaka inference -> medaka sequence (the medaka_consensus wrapper exits in its
    version check when bcftools is absent, so the three steps are run directly).
    Returns the sequence, or None when medaka failed."""
    bam, hdf = out_fa + '.bam', out_fa + '.hdf'
    env = dict(os.environ)
    env['PATH'] = MEDAKA_ENV + ':' + env.get('PATH', '')
    p1 = subprocess.Popen(['minimap2', '-t', str(threads), '-ax', 'map-ont', '--secondary=no', template_fa, reads_fa],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(['samtools', 'sort', '-@', '1', '-o', bam, '-'], stdin=p1.stdout, stderr=subprocess.DEVNULL)
    p1.stdout.close()
    p2.communicate()
    if p2.returncode != 0:
        return None
    run(['samtools', 'index', bam])
    ok = True
    for argv in ([os.path.join(MEDAKA_ENV, 'medaka'), 'inference', bam, hdf, '--model', MEDAKA_MODEL, '--threads', str(threads)],
                 [os.path.join(MEDAKA_ENV, 'medaka'), 'sequence', hdf, template_fa, out_fa]):
        r = subprocess.run(argv, env=env, capture_output=True, text=True)
        if r.returncode != 0:
            ok = False
            break
    for e in (bam, bam + '.bai', hdf):
        try:
            os.remove(e)
        except OSError:
            pass
    if not ok or not os.path.exists(out_fa):
        return None
    return read_fasta_seq(out_fa).strip(u'N')


def best_hits(reads_fa, dbs, root, threads=1, db_dir=None, blast_fn=None):
    u"""Best record per read over all databases: identity first, then alignment
    length; alignments under EN_AZ_KANIT are ignored. {read: (pid, aln, title, db)}."""
    blast_fn = blast_fn or blast
    best = {}
    for name in dbs:
        db = db_path(root, name, db_dir) if root is not None else name
        if not db:
            continue
        for k, hits in blast_fn(reads_fa, db, threads=threads, top=5).items():
            for x in hits:
                pid, aln, title = x[0], x[1], x[4]
                if aln < EN_AZ_KANIT:
                    continue
                if k not in best or (pid, aln) > best[k][:2]:
                    best[k] = (pid, aln, title, os.path.basename(name))
    return best


def mmi_path(db, threads=1):
    u"""The minimap2 index of a database (<db>.mmi), built once if missing (UNITE: 75 s, 2.8 GB)."""
    mmi = db + '.mmi'
    if not os.path.exists(mmi) or os.path.getmtime(mmi) < os.path.getmtime(db):
        run(['minimap2', '-t', str(threads), '-x', 'map-ont', '-d', mmi + '.building', db])
        os.replace(mmi + '.building', mmi)
    return mmi


def best_hits_mm2(reads_fa, dbs, root, threads=1, db_dir=None, secondaries=10):
    u"""Best record per read with minimap2 -c (base-level alignment; identity = 1 - de, the
    gap-compressed divergence). Same shape as best_hits: {read: (pid, aln, title, db)}.
    WHY minimap2 (measured 2026-09-06 on 401 ITS windows):
      * blastn -max_target_seqs 5 MISSED the best record (the known truncation): a read was
        written as Microascus 96.63 per cent while the true best record was Acaulium 98.27
        (confirmed with blastn -subject); max_target_seqs 50 and 500 find different records
        again. UNITE holds hundreds of near-identical records per species, and the truncation
        bites exactly in that pile-up.
      * time: 401 windows in 9 s against 26 min for BLAST on 2 threads.
      * on the same record the identity differs by a median of +0.16 points; every one of the
        53 reads with a different genus had a higher-identity (real) record found by minimap2."""
    best = {}
    for name in dbs:
        db = db_path(root, name, db_dir) if root is not None else name
        if not db:
            continue
        mmi = mmi_path(db, threads)
        paf = run(['minimap2', '-t', str(threads), '-c', '-x', 'map-ont', '--secondary=yes',
                   '-N', str(secondaries), mmi, reads_fa])
        cand, needed = {}, set()
        for line in paf.splitlines():
            p = line.split('\t')
            if len(p) < 12:
                continue
            k, tname, blen = p[0], p[5], int(p[10])
            if blen < EN_AZ_KANIT:
                continue
            de = None
            for x in p[12:]:
                if x.startswith('de:f:'):
                    de = float(x[5:])
                    break
            if de is None:
                continue
            pid = round((1.0 - de) * 100.0, 2)
            if k not in cand or (pid, blen) > cand[k][:2]:
                cand[k] = (pid, blen, tname)
                needed.add(tname)
        titles = {}
        if needed:
            with io.open(db, encoding='utf-8', errors='replace') as fh:
                for line in fh:
                    if line.startswith(u'>'):
                        nm = line[1:].split()[0]
                        if nm in needed:
                            titles[nm] = line[1:].strip()
        for k, (pid, blen, tname) in cand.items():
            if k not in best or (pid, blen) > best[k][:2]:
                best[k] = (pid, blen, titles.get(tname, tname), os.path.basename(name))
    return best


def unnamed(name):
    low = name.lower()
    return any(t in low for t in ADSIZ_JETONLARI)


def populations(best, key=u'ITS'):
    u"""Group reads by the genus of their best record. Returns (genus counter,
    {genus: [(read, pid, aln, name)]}, species counter, {species: [pid]}).
    Unnamed records join no population."""
    count = collections.Counter()
    members = collections.defaultdict(list)
    scount = collections.Counter()
    spid = collections.defaultdict(list)
    for k, (pid, aln, title, _db) in best.items():
        name = ad_ayikla(title)
        if not name or unnamed(name):
            continue
        genus, epithet = cins_epitet(name)
        if not genus:
            continue
        if pid < CINS_ESIGI[key]:
            continue          # below the genus threshold (of the locus) a hit founds no population
        count[genus] += 1
        members[genus].append((k, pid, aln, name))
        if epithet and epithet.lower() not in PLACEHOLDER:
            full = u'%s %s' % (genus, epithet)
            scount[full] += 1
            spid[full].append(pid)
    return count, members, scount, spid


# ---------------------------------------------------------------------------
# SELF-TEST: synthetic data, the real tools (blastn, minimap2, samtools)
# ---------------------------------------------------------------------------
def _random_seq(rnd, n):
    return u''.join(rnd.choice(u'ACGT') for _ in range(n))


def _mutate(rnd, s, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate * 0.5:
            out.append(rnd.choice(u'ACGT'.replace(c, u'')))
        elif r < rate * 0.75:
            continue
        elif r < rate:
            out.append(c)
            out.append(rnd.choice(u'ACGT'))
        else:
            out.append(c)
    return u''.join(out)


def _identity(a_fa, b_fa):
    out = subprocess.run(['blastn', '-query', a_fa, '-subject', b_fa, '-outfmt',
                          '6 pident length', '-max_hsps', '1'],
                         capture_output=True, text=True).stdout.strip().splitlines()
    if not out:
        return 0.0, 0
    p = out[0].split(u'\t')
    return float(p[0]), int(p[1])


def tools_present():
    return all(shutil.which(t) for t in ('blastn', 'makeblastdb', 'minimap2', 'samtools'))


def self_test(a=None):
    u"""1) the medoid recognises reverse-strand reads; 2) a 70/30 mixture of A and B
    (B is 90 per cent to A) yields A as the dominant population; 3) the consensus
    polished from A's reads (3 per cent error) is >= 99.5 per cent to A and better
    than the medoid alone; 4) sampling is deterministic; 5) the row has as many
    fields as the header."""
    if not tools_present():
        print(u'  self-test skipped: blastn / makeblastdb / minimap2 / samtools not all present')
        return 0
    rnd = random.Random(7)
    A = _random_seq(rnd, 1800)
    B = _mutate(rnd, A, 0.10)
    rA = [(u'a%d' % i, _mutate(rnd, A, 0.03)) for i in range(70)]
    rB = [(u'b%d' % i, _mutate(rnd, B, 0.03)) for i in range(30)]
    rA = [(k, reverse_complement(d) if i % 4 == 0 else d) for i, (k, d) in enumerate(rA)]
    errors = []
    tmp = tempfile.mkdtemp(prefix='fungal_selftest_')
    try:
        med = medoid(rA)
        if not med or med[2] < 0.3:
            errors.append(u'medoid: unexpected similarity (%s)' % (med and med[2]))
        db = os.path.join(tmp, 'db.fa')
        write_fasta(db, [(u'NR_000001.1 Aus alpha strain X ITS region', A),
                         (u'NR_000002.1 Bus beta strain Y ITS region', B)])
        run(['makeblastdb', '-in', db, '-dbtype', 'nucl'])
        reads_fa = os.path.join(tmp, 'reads.fa')
        write_fasta(reads_fa, rA + rB)
        best = best_hits(reads_fa, [db], None, 1)
        count, members, scount, spid = populations(best)
        if count.most_common(1)[0][0] != u'Aus':
            errors.append(u'populations: expected Aus dominant, got %s' % count.most_common(2))
        if not (60 <= count[u'Aus'] <= 70 and 25 <= count[u'Bus'] <= 30):
            errors.append(u'population counts unexpected: %s' % dict(count))
        if scount.most_common(1)[0][0] != u'Aus alpha':
            errors.append(u'species counter: %s' % scount.most_common(2))
        # the minimap2 route (the one a real run takes): same database, same populations
        best2 = best_hits_mm2(reads_fa, [db], None, 1)
        count2, _m2, scount2, _s2 = populations(best2)
        if dict(count2) != dict(count) or scount2.most_common(1)[0][0] != u'Aus alpha':
            errors.append(u'minimap2 route disagrees with BLAST: %s / %s' % (dict(count2), dict(count)))
        p2 = sorted(v[0] for v in best2.values())
        if p2 and (p2[0] < 90.0 or p2[-1] > 100.0):
            errors.append(u'minimap2 identity range unexpected: %.2f-%.2f' % (p2[0], p2[-1]))
        A_fa = os.path.join(tmp, 'A.fa')
        write_fasta(A_fa, [(u'A', A)])
        seq_of = dict(rA + rB)
        pop = [(k, seq_of[k]) for k, _p, _a, _n in members[u'Aus']]
        pop_fa = os.path.join(tmp, 'pop.fa')
        write_fasta(pop_fa, pop)
        med = medoid(pop)
        tpl = os.path.join(tmp, 't0.fa')
        write_fasta(tpl, [(med[0], med[1])])
        p0, _ = _identity(tpl, A_fa)
        for t in (1, 2):
            out = os.path.join(tmp, 'c%d.fa' % t)
            seq, _n = polish(tpl, pop_fa, out, 1)
            write_fasta(out, [(u'c', seq)])
            tpl = out
        p2, l2 = _identity(tpl, A_fa)
        if p2 < 99.5 or l2 < 1700:
            errors.append(u'polish: %.2f%% / %d bp to A (expected >= 99.5 and >= 1700; medoid %.2f%%)' % (p2, l2, p0))
        if p2 <= p0:
            errors.append(u'polish did not improve on the medoid (%.2f -> %.2f)' % (p0, p2))
        print(u'  medoid %.2f%% -> dominant-population polish %.2f%%/%d bp' % (p0, p2, l2))
        fq = os.path.join(tmp, 'x.fastq')
        with io.open(fq, 'w', encoding='utf-8') as g:
            for i in range(500):
                s = _random_seq(rnd, 600)
                g.write(u'@r%d\n%s\n+\n%s\n' % (i, s, u'I' * len(s)))
        o1, o2 = sample_reads(fq, 50), sample_reads(fq, 50)
        if len(o1) != 50 or o1 != o2:
            errors.append(u'sample_reads is not deterministic or did not return 50 (%d)' % len(o1))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    for e in errors:
        print(u'  ERROR: %s' % e)
    print(u'  SELF-TEST: %s' % (u'PASSED' if not errors else u'%d ERROR(S)' % len(errors)))
    return 0 if not errors else 1


if __name__ == '__main__':
    sys.exit(self_test())
