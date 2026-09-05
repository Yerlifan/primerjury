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
# HOW TO RUN IT    : python3 verification/fungal_bin_identity.py --root . [--reads 150] [--rounds 2]
#                    python3 verification/fungal_bin_identity.py --self-test
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
  python3 verification/fungal_bin_identity.py --root . [--reads 150] [--rounds 2]
      [--threads 1] [--bins F2-2_500148,...] [--redo] [--db-dir REFERENCE_DB]
  python3 verification/fungal_bin_identity.py --self-test
An existing bin output is skipped unless --redo is given (a restarted chain
continues where it stopped).
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
from identity_verification import EN_AZ_KANIT, ADSIZ_JETONLARI              # noqa: E402
from target_identity import ad_ayikla, cins_epitet                          # noqa: E402

SET_DIR = os.path.join('referans_konsensus', 'fungal_polish', 'consensus')   # the new candidate set
POP_DIR = os.path.join('referans_konsensus', 'fungal_polish', 'population')
SEED = 20260905
K_MEDOID = 15          # k-mer size for the medoid
MIN_POPULATION = 5     # no consensus from fewer reads than this
MIN_DEPTH = 3          # samtools consensus -d
PLACEHOLDER = ('sp', 'sp.', 'spp', 'spp.', 'cf', 'cf.', 'aff', 'aff.')
FUNGAL_LIBS = ('F1', 'F2')
ITS_DB = [d for lok, dbs, _a in MANTAR_LOKUSLARI if lok == u'ITS' for d in dbs]


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


def unnamed(name):
    low = name.lower()
    return any(t in low for t in ADSIZ_JETONLARI)


def populations(best):
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
        count[genus] += 1
        members[genus].append((k, pid, aln, name))
        if epithet and epithet.lower() not in PLACEHOLDER:
            full = u'%s %s' % (genus, epithet)
            scount[full] += 1
            spid[full].append(pid)
    return count, members, scount, spid


def process_bin(root, label, path, a, tmp, blast_fn=None):
    reads = sample_reads(path, a.reads)
    if len(reads) < MIN_POPULATION:
        return dict(bin=label, reads=len(reads), note=u'too few reads (%d)' % len(reads))
    kd = os.path.join(tmp, label)
    os.makedirs(kd)
    reads_fa = os.path.join(kd, 'reads.fa')
    write_fasta(reads_fa, reads)
    best = best_hits(reads_fa, ITS_DB, root, a.threads, a.db_dir, blast_fn)
    count, members, scount, spid = populations(best)
    seq_of = dict(reads)
    pop_dir = os.path.join(root, POP_DIR)
    if not os.path.isdir(pop_dir):
        os.makedirs(pop_dir)
    rows = [u'read\tbp\tgenus\tname\tidentity %\talignment bp\tdatabase']
    for k, d in reads:
        if k in best:
            pid, aln, title, db = best[k]
            name = ad_ayikla(title) or u''
            rows.append(u'%s\t%d\t%s\t%s\t%.2f\t%d\t%s'
                        % (k, len(d), cins_epitet(name)[0] if name else u'', name, pid, aln, db))
        else:
            rows.append(u'%s\t%d\t\t\t\t\t' % (k, len(d)))
    io.open(os.path.join(pop_dir, label + '.tsv'), 'w', encoding='utf-8',
            newline='\n').write(u'\n'.join(rows) + u'\n')

    total = sum(count.values())
    r = dict(bin=label, reads=len(reads), assigned=total)
    if scount:
        s1, n1 = scount.most_common(1)[0]
        pids = sorted(spid[s1])
        r.update(read_species=s1, species_share=100.0 * n1 / sum(scount.values()),
                 species_n=n1, species_median=pids[len(pids) // 2])
    if total:
        g1, n1 = count.most_common(1)[0]
        g2, n2 = count.most_common(2)[1] if len(count) > 1 else (u'-', 0)
        r.update(genus=g1, share=100.0 * n1 / total, genus2=g2, share2=100.0 * n2 / total)
        pop = [(k, seq_of[k]) for k, _p, _a, _n in members[g1]]
        pop_note = u'dominant population %s %d/%d' % (g1, n1, total)
    else:
        pop = list(reads)
        r.update(genus=u'-', share=0.0, genus2=u'-', share2=0.0)
        pop_note = u'no ITS hit; all reads form one population'
    if len(pop) < MIN_POPULATION:
        r['note'] = u'%s: %d reads < %d, no consensus made' % (pop_note, len(pop), MIN_POPULATION)
        return r
    med = medoid(pop)
    r.update(template=med[0], template_bp=len(med[1]), template_similarity=med[2])
    pop_fa = os.path.join(kd, 'population.fa')
    write_fasta(pop_fa, pop)
    tpl = os.path.join(kd, 'template0.fa')
    write_fasta(tpl, [(med[0], med[1])])
    seq, n_in = med[1], 0
    for rnd in range(1, a.rounds + 1):
        out = os.path.join(kd, 'round%d.fa' % rnd)
        seq, n_in = polish(tpl, pop_fa, out, a.threads)
        if len(seq) < 200:
            r['note'] = u'%s; round %d left %d bp, previous round kept' % (pop_note, rnd, len(seq))
            seq = read_fasta_seq(tpl)
            break
        write_fasta(out, [(label, seq)])
        tpl = out
    r.update(bp=len(seq), n_internal=n_in, sequence=seq, pop_note=pop_note, pop_n=len(pop))
    query = os.path.join(kd, 'final.fa')
    write_fasta(query, [(label, seq)])
    hits = {}
    for lok, dbs, _an in MANTAR_LOKUSLARI:
        acc = []
        for name in dbs:
            db = db_path(root, name, a.db_dir) if root is not None else name
            if db:
                acc += (blast_fn or blast)(query, db, threads=a.threads).get(label, [])
        hits[lok] = acc
    decisions = {lok: lokus_karari(hits.get(lok), an, ad_ayikla, cins_epitet, ADSIZ_JETONLARI)
                 for lok, _d, an in MANTAR_LOKUSLARI}
    name, ident, note = birlestir(decisions)
    ra, rp, rl, rlok, rsecond, rgap = raporlanan_yontem(hits, ad_ayikla, cins_epitet, ADSIZ_JETONLARI)
    r.update(decisions=decisions, name=name, identity=ident, decision_note=note,
             reported=(u'%s (%.2f%%, %d bp, %s)' % (ra, rp, rl, rlok)) if ra else u'-',
             reported_second=(u'%s, gap %.2f' % (rsecond, rgap)) if rsecond else u'-')
    return r


HEADER = [u'bin', u'sampled reads', u'ITS-assigned reads', u'dominant genus', u'share %',
          u'second genus', u'second share %', u'read species (ITS)', u'species share %',
          u'reads at species', u'species median identity %', u'template read', u'template bp',
          u'template k-mer similarity', u'polishing reads', u'consensus bp', u'internal N',
          u'18S name', u'18S %', u'18S bp', u'ITS name', u'ITS %', u'ITS bp',
          u'28S name', u'28S %', u'28S bp', u'COMBINED NAME', u'identity %', u'decision note',
          u'reported method', u'reported second', u'note']


def make_row(r):
    def f(x, fmt=u'%s'):
        return (fmt % x) if x not in (None, u'') else u'-'
    dec = r.get('decisions', {})
    p = [r['bin'], f(r.get('reads')), f(r.get('assigned')), f(r.get('genus')),
         f(r.get('share'), u'%.0f'), f(r.get('genus2')), f(r.get('share2'), u'%.0f'),
         f(r.get('read_species')), f(r.get('species_share'), u'%.0f'), f(r.get('species_n')),
         f(r.get('species_median'), u'%.2f'), f(r.get('template')), f(r.get('template_bp')),
         f(r.get('template_similarity'), u'%.3f'), f(r.get('pop_n')), f(r.get('bp')),
         f(r.get('n_internal'))]
    for lok in (u'18S', u'ITS', u'28S'):
        k = dec.get(lok)
        p += [k['ad'], u'%.2f' % k['kimlik'], u'%d' % k['hizalama']] if k else [u'-', u'-', u'-']
    p += [f(r.get('name')), f(r.get('identity')), f(r.get('decision_note')), f(r.get('reported')),
          f(r.get('reported_second')), f(r.get('note'))]
    return u'\t'.join(p)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', default='.')
    ap.add_argument('--reads', type=int, default=150)
    ap.add_argument('--rounds', type=int, default=2)
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--bins', default='')
    ap.add_argument('--db-dir', default=None)
    ap.add_argument('--out', default='FUNGAL_BIN_IDENTITY.tsv')
    ap.add_argument('--redo', action='store_true')
    ap.add_argument('--self-test', action='store_true')
    a = ap.parse_args(argv)
    if a.self_test:
        return self_test(a)
    root = os.path.abspath(a.root)
    for name in ITS_DB:
        db_path(root, name, a.db_dir)
    files = bin_files(root)
    bins = [b for b in sorted(files) if b.split(u'-')[0] in FUNGAL_LIBS]
    if a.bins:
        want = set(a.bins.split(u','))
        bins = [b for b in bins if b in want]
    set_dir = os.path.join(root, SET_DIR)
    if not os.path.isdir(set_dir):
        os.makedirs(set_dir)
    out_path = os.path.join(root, a.out)
    previous = {}
    if os.path.exists(out_path) and not a.redo:
        for s in io.open(out_path, encoding='utf-8', errors='replace').read().splitlines()[1:]:
            if s.strip():
                previous[s.split(u'\t')[0]] = s
    print(u'  fungal bins: %d  (reads %d, rounds %d, threads %d)' % (len(bins), a.reads, a.rounds, a.threads))
    rows = {}
    tmp = tempfile.mkdtemp(prefix='fungal_polish_')
    levels = collections.Counter()
    try:
        for n, label in enumerate(bins, 1):
            fa = os.path.join(set_dir, label + '.fasta')
            if os.path.exists(fa) and not a.redo and label in previous:
                rows[label] = previous[label]
                print(u'  [%3d/%d] %-16s done, skipped' % (n, len(bins), label))
                continue
            t0 = time.time()
            r = process_bin(root, label, files[label], a, tmp)
            if r.get('sequence'):
                with io.open(fa, 'w', encoding='utf-8') as g:
                    g.write(u'>%s fungal_polish population=%s share=%.0f reads=%d template=%s rounds=%d\n%s\n'
                            % (label, r.get('genus'), r.get('share', 0.0), r.get('pop_n', 0),
                               r.get('template'), a.rounds, r['sequence']))
                nm = r['name'] or u''
                lvl = (u'species' if u' ' in nm and u'sp.' not in nm and u'cf.' not in nm
                       else (u'cf.' if u'cf.' in nm else (u'genus' if nm.endswith(u'sp.') else u'unnamed')))
                levels[lvl] += 1
                d = r['decisions']
                print(u'  [%3d/%d] %-16s %s %.0f%% | ITS %s %.2f%%/%d | 28S %s %.2f%% | -> %s (%s) %.0f s'
                      % (n, len(bins), label, r.get('genus'), r.get('share', 0.0),
                         d[u'ITS']['ad'][:24], d[u'ITS']['kimlik'], d[u'ITS']['hizalama'],
                         d[u'28S']['ad'][:22], d[u'28S']['kimlik'], r['name'], lvl, time.time() - t0))
            else:
                print(u'  [%3d/%d] %-16s %s' % (n, len(bins), label, r.get('note')))
            rows[label] = make_row(r)
            sys.stdout.flush()
            io.open(out_path, 'w', encoding='utf-8', newline='\n').write(
                u'\t'.join(HEADER) + u'\n' + u'\n'.join(rows[k] for k in sorted(rows)) + u'\n')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    io.open(out_path, 'w', encoding='utf-8', newline='\n').write(
        u'\t'.join(HEADER) + u'\n' + u'\n'.join(rows[k] for k in sorted(rows)) + u'\n')
    print(u'')
    print(u'  levels (this run): %s' % dict(levels))
    print(u'  set: %s' % set_dir)
    print(u'  report: %s' % out_path)
    return 0


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
        r = dict(bin=u'F9-9_1', reads=3, note=u'too few reads (3)')
        if len(make_row(r).split(u'\t')) != len(HEADER):
            errors.append(u'make_row field count does not match the header')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    for e in errors:
        print(u'  ERROR: %s' % e)
    print(u'  SELF-TEST: %s' % (u'PASSED' if not errors else u'%d ERROR(S)' % len(errors)))
    return 0 if not errors else 1


if __name__ == '__main__':
    sys.exit(main())
