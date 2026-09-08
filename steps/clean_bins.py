# -*- coding: utf-8 -*-
"""CLEAN BINS: one organism per bin, full-length amplicons only (bin_reads v4, 2026-09-07/08).

WHAT IT IS FOR
--------------
Groups the reads of one barcode by organism and writes the groups in the
layout the chain reads:

    <out>/<name>/<name>-reads_BIN<n>.fastq      full-length reads + concatemer segments
    <out>/<name>/<name>_centres.fasta            one centre read per bin
    <out>/<name>/BIN_TABLE.tsv                   counts, share, centre length, median identity
    <out>/<name>/_unassigned/                    everything that fitted nowhere (counted, never deleted)

WHY A NEW VERSION
-----------------
bin_reads.py (v3) opened a window per LENGTH PEAK and assigned at 90 per cent
identity. Measured on the study's from-scratch SUP root (2026-09-07, A1-1): of
28 bins, 11 were full-length amplicons, 10 were FRAGMENTS of the same
organisms (500-760 bp) and 7 were CONCATEMERS (2.9 kb: two 16S copies ligated
end to end, 30 of 30 reads). Fragment bins showed up in the primer check as
"members that do not amplify" (the primer site is not in the read);
concatemer bins were named after whichever half sorted first. A bin was a
length class, not an organism. The user's decision: "bins must be clean,
organisms must not mix".

THE METHOD
----------
0. LENGTH CLASS against the library's expected amplicon (AMPLICON_BP):
   full = 0.80-1.25x, fragment = 300 bp-0.80x, concatemer = >1.25x (the
   1.25-1.6x band is handled the concatemer way: a full segment is cut out
   if there is one).
1. SEEDS from the full reads, core length only (0.85-1.15x): a 1,753 bp seed
   (1.17x, 16S plus 300 bp of the next copy) once collected 11,261 segments
   and 24 full reads. minimap2 ava-ont, identity >= 0.97, >= 80 per cent of
   the shorter one covered, connected components; the centre is the member
   nearest the component's median length.
2. TIERED ASSIGNMENT of every full read (map-ont, secondaries on): >= 0.97
   identity over >= 85 per cent of the read goes straight in; 0.95-0.97 only
   when the nearest centre leads the second nearest by >= 0.02 (sister
   genera sit at 92-95 per cent and must not be forced). Measured: at a flat
   0.97 44 per cent of the full reads were left out, their identity to the
   centre had a median of 0.962 and a consensus centre did not help (0.961):
   the loss is read error, not a bad centre. A SECOND seed round is drawn
   from the reads left over, then assignment again.
3. CONCATEMERS: each read is aligned to all centres (secondaries on); every
   interval that covers >= 80 per cent of a centre at >= 0.97 identity AND is
   no longer than 1.25x the centre is CUT OUT as a read of its own
   (<id>_p1, _p2 ...) and joins that centre's bin. The length cap matters:
   minimap2 chains the shared end of two tandem copies into one alignment
   with a large insertion and returns a 1.6-2.0x "segment" (measured
   2026-09-08: 0.6 to 3 per cent of the segments, some small bins made of
   nothing else). A two-copy segment is not one amplicon and is skipped.
4. FRAGMENTS are counted, not mapped. They never enter a bin, the identity
   step and the primer check use full-length reads only, and mapping them
   is ruinously slow on rRNA operons: 2.3 million fragments of one fungal
   barcode against 524 bin centres produced 0 bytes of PAF in 3.7 hours,
   because every fragment chains to every centre. --fragments-max > 0 maps a
   seeded random sample into _fragments/ for diagnosis only.
5. Clusters with >= MIN_READS full-length reads (segments included) become
   bins, numbered largest first. Nothing is deleted: the unassigned reads of
   every class are written under _unassigned/ and counted in the table.

Two more measured traps are handled downstream, in select_bins.py, because
they need the finished table: a SATELLITE bin (the second seed round picks a
read with a ~310 bp library tail, and the tailed reads of an organism that
already has a bin form a second bin of the same organism) and the two-copy
segments of runs made before the length cap. See select_bins.py.

The PAF of every minimap2 call is streamed from a file, never held in memory:
a 1 GB PAF held as a Python string pushed the process into swap (measured).

RUN
---
    python3 steps/clean_bins.py --fastq barcode01.fastq --name A1-1 --out "fastq files_bins" \\
        [--amplicon 1500] [--identity 0.97] [--threads 3] [--seeds 1200] [--seed 7] [--fragments-max 0]

--amplicon overrides the AMPLICON_BP lookup by the name's library prefix.
"""
from __future__ import print_function
import argparse
import collections
import importlib.util
import io
import os
import random
import shutil
import subprocess
import sys
import tempfile

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, KOK)
try:
    from verification.identity_verification import EN_AZ_OKUMA as MIN_READS_DEFAULT, AMPLICON_BP
except Exception:                       # the module needs the aligner and the databases
    MIN_READS_DEFAULT = 5
    AMPLICON_BP = {'A1': 1500, 'B': 1500, 'A2': 4300, 'F1': 3700, 'F2': 3700}


def _bin_reads():
    """The v3 module: tool lookup, fastq stream and the seed clustering are shared."""
    path = os.path.join(KOK, 'steps', 'bin_reads.py')
    spec = importlib.util.spec_from_file_location('bin_reads', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


FULL_LOW, FULL_HIGH = 0.80, 1.25        # full-length class (x amplicon)
SEED_LOW, SEED_HIGH = 0.85, 1.15        # seed (centre) candidates: core length only
FRAGMENT_MIN = 300                      # bp; shorter reads are counted, not processed
ASSIGN_IDENTITY = 0.97                  # v3 had 0.90: Proteiniphilum and Petrimonas shared a bin
ASSIGN_IDENTITY_2 = 0.95                # second tier, only with a clear margin
ASSIGN_MARGIN = 0.02
ASSIGN_COVERAGE = 0.85                  # share of the READ that must be covered (v3: 0.80)
SEGMENT_COVERAGE = 0.80                 # a concatemer segment must cover this share of the centre


def length_class(length, expected):
    """'short' | 'fragment' | 'full' | 'concatemer' for one read length."""
    if length < FRAGMENT_MIN:
        return 'short'
    if length < FULL_LOW * expected:
        return 'fragment'
    if length <= FULL_HIGH * expected:
        return 'full'
    return 'concatemer'


def paf_to_file(cmd, path):
    """Run minimap2 with stdout in a file and return the file, opened for reading."""
    with io.open(path, 'wb') as g, io.open(path + '.err', 'wb') as e:
        rc = subprocess.call(cmd, stdout=g, stderr=e)
    if rc != 0:
        msg = io.open(path + '.err', encoding='utf-8', errors='replace').read()[:400]
        raise SystemExit('STOPPED: %s -> %d: %s' % (cmd[0], rc, msg))
    return io.open(path, encoding='utf-8', errors='replace')


def paf_rows(paf):
    """(read, read_len, qs, qe, target, target_len, ts, te, identity) per PAF line with a de tag.
    paf may be a string or a line stream (an open file)."""
    for line in (paf.splitlines() if isinstance(paf, type(u'')) else paf):
        p = line.split('\t')
        if len(p) < 12:
            continue
        de = None
        for tag in p[12:]:
            if tag.startswith('de:f:'):
                de = float(tag[5:])
                break
        if de is None:
            continue
        yield (p[0], int(p[1]), int(p[2]), int(p[3]), p[5], int(p[6]), int(p[7]), int(p[8]), 1.0 - de)


def assign_from_rows(rows, identity, coverage, tiered=True):
    """Tiered assignment from PAF rows: read -> (centre, identity). Pure; the aligner is outside."""
    best = collections.defaultdict(dict)            # read -> {centre: best identity}
    for rid, ql, qs, qe, centre, _tl, _ts, _te, ident in rows:
        cov = (qe - qs) / float(ql) if ql else 0.0
        if cov >= coverage and ident > best[rid].get(centre, 0.0):
            best[rid][centre] = ident
    out = {}
    for rid, d in best.items():
        order = sorted(d.items(), key=lambda x: -x[1])
        centre, ident = order[0]
        second = order[1][1] if len(order) > 1 else 0.0
        if ident >= identity or (tiered and ident >= ASSIGN_IDENTITY_2 and ident - second >= ASSIGN_MARGIN):
            out[rid] = (centre, ident)
    return out


def assign(mm2, threads, centres_fa, reads_fq, identity, coverage, tiered=True):
    paf = paf_to_file([mm2, '-x', 'map-ont', '-c', '--secondary=yes', '-N', '5', '-p', '0.8',
                       '-t', str(threads), centres_fa, reads_fq], reads_fq + '.paf')
    out = assign_from_rows(paf_rows(paf), identity, coverage, tiered)
    paf.close()
    return out


def segments_from_rows(rows, identity, seg_coverage=SEGMENT_COVERAGE, cap=FULL_HIGH):
    """Concatemer segments from PAF rows: read -> [(qs, qe, centre, identity)], non-overlapping,
    longest-and-most-identical first; an interval longer than cap x the centre is a two-copy
    chain and is refused. Pure."""
    cand = collections.defaultdict(list)
    for rid, _ql, qs, qe, centre, tl, ts, te, ident in rows:
        if ident >= identity and (te - ts) >= seg_coverage * tl and (qe - qs) <= cap * tl:
            cand[rid].append((qs, qe, centre, ident))
    out = {}
    for rid, L in cand.items():
        L.sort(key=lambda x: (-(x[1] - x[0]) * x[3], x[0]))
        chosen = []
        for qs, qe, centre, ident in L:
            if all(qe <= s or qs >= e for s, e, _c, _i in chosen):
                chosen.append((qs, qe, centre, ident))
        out[rid] = sorted(chosen)
    return out


def segment(mm2, threads, centres_fa, reads_fq, identity):
    paf = paf_to_file([mm2, '-x', 'map-ont', '-c', '--secondary=yes', '-N', '20', '-p', '0.5',
                       '-t', str(threads), centres_fa, reads_fq], reads_fq + '.seg.paf')
    out = segments_from_rows(paf_rows(paf), identity)
    paf.close()
    return out


def write_centres(centres, path):
    with io.open(path, 'w', encoding='ascii') as g:
        for name, seq in centres.items():
            g.write(u'>%s\n%s\n' % (name, seq))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--fastq', required=True)
    ap.add_argument('--name', required=True, help='library-year, e.g. A1-1; the prefix selects the amplicon')
    ap.add_argument('--out', required=True)
    ap.add_argument('--amplicon', type=int, default=0, help='expected amplicon bp (default: AMPLICON_BP by prefix)')
    ap.add_argument('--identity', type=float, default=ASSIGN_IDENTITY)
    ap.add_argument('--threads', type=int, default=3)
    ap.add_argument('--seeds', type=int, default=1200)
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--min-reads', type=int, default=MIN_READS_DEFAULT)
    ap.add_argument('--fragments-max', type=int, default=0,
                    help='0 = fragments are counted, not mapped (default); >0 = map a seeded sample into _fragments/')
    a = ap.parse_args(argv)
    br = _bin_reads()
    mm2 = br.tool('minimap2')
    expected = a.amplicon or AMPLICON_BP.get(a.name.split('-')[0], 0)
    if not expected:
        raise SystemExit('STOPPED: no expected amplicon for %s (AMPLICON_BP or --amplicon)' % a.name)
    out = os.path.join(os.path.abspath(a.out), a.name)
    for d in (out, os.path.join(out, '_unassigned')):
        if not os.path.isdir(d):
            os.makedirs(d)
    tmp = tempfile.mkdtemp(prefix='_tmp_clean_bins_', dir=os.path.dirname(out))
    try:
        # 0. length classes, streamed to class files
        path = {k: os.path.join(tmp, k + '.fq') for k in ('full', 'fragment', 'concatemer', 'short')}
        writers = {k: io.open(v, 'w', encoding='ascii') for k, v in path.items()}
        count = collections.Counter()
        for rid, seq, qual in br.fastq_stream(a.fastq):
            k = length_class(len(seq), expected)
            writers[k].write(u'@%s\n%s\n+\n%s\n' % (rid, seq, qual))
            count[k] += 1
        for w in writers.values():
            w.close()
        total = sum(count.values())
        print('  %s: %d reads; amplicon ~%d bp; full %d (%.1f%%), fragment %d, concatemer %d, <%d bp %d'
              % (a.name, total, expected, count['full'], 100.0 * count['full'] / max(1, total),
                 count['fragment'], count['concatemer'], FRAGMENT_MIN, count['short']))
        if count['full'] < a.min_reads:
            raise SystemExit('STOPPED: %s has no full-length reads' % a.name)
        # 1-2. seeds and tiered assignment, two rounds
        random.seed(a.seed)
        centres = collections.OrderedDict()
        assignment = {}
        pool = path['full']
        for rnd in (1, 2):
            n_pool = sum(1 for _ in br.fastq_stream(pool))
            if n_pool < a.min_reads:
                break
            candidates = [j for j, (_r, s, _q) in enumerate(br.fastq_stream(pool))
                          if SEED_LOW * expected <= len(s) <= SEED_HIGH * expected]
            if len(candidates) < a.min_reads:
                break
            pick = set(random.sample(candidates, min(a.seeds, len(candidates))))
            seeds = collections.OrderedDict((r, s) for j, (r, s, _q) in enumerate(br.fastq_stream(pool)) if j in pick)
            comp = br.cluster_seeds(seeds, a.identity, mm2, a.threads, tmp, 'seeds%d' % rnd)
            for c in comp:
                if c not in centres:
                    centres[c] = seeds[c]
            centres_fa = os.path.join(tmp, 'centres%d.fa' % rnd)
            write_centres(centres, centres_fa)
            assignment.update(assign(mm2, a.threads, centres_fa, pool, a.identity, ASSIGN_COVERAGE))
            left = os.path.join(tmp, 'left%d.fq' % rnd)
            n_left = 0
            with io.open(left, 'w', encoding='ascii') as g:
                for rid, seq, qual in br.fastq_stream(pool):
                    if rid not in assignment:
                        g.write(u'@%s\n%s\n+\n%s\n' % (rid, seq, qual))
                        n_left += 1
            print('  round %d: seeds %d -> centres %d (total %d); assigned %d, left %d'
                  % (rnd, len(seeds), len(comp), len(centres), len(assignment), n_left))
            pool = left
            if n_left < a.min_reads:
                break
        centres_fa = os.path.join(tmp, 'centres_final.fa')
        write_centres(centres, centres_fa)
        # 3. concatemer segments
        segs = segment(mm2, a.threads, centres_fa, path['concatemer'], a.identity) if count['concatemer'] else {}
        seg_count = collections.Counter(c for L in segs.values() for _qs, _qe, c, _i in L)
        # bins (before fragments: fragments never count)
        full_count = collections.Counter(c for c, _i in assignment.values())
        idents = collections.defaultdict(list)
        for c, ident in assignment.values():
            idents[c].append(ident)
        bin_count = collections.Counter(full_count)
        bin_count.update(seg_count)
        passing = sorted((c for c in centres if bin_count[c] >= a.min_reads), key=lambda c: (-bin_count[c], c))
        number = {c: i + 1 for i, c in enumerate(passing)}
        # 4. fragments: counted; mapped only on request, to the bin centres, as a seeded sample
        frag_sample = set()                                 # empty = not mapped
        frag_assignment = {}
        if a.fragments_max > 0 and count['fragment'] and passing:
            bin_centres_fa = os.path.join(tmp, 'centres_bins.fa')
            write_centres(collections.OrderedDict((c, centres[c]) for c in passing), bin_centres_fa)
            frag_fq = path['fragment']
            if count['fragment'] > a.fragments_max:
                frag_sample = set(random.Random(a.seed).sample(range(count['fragment']), a.fragments_max))
                frag_fq = os.path.join(tmp, 'fragment_sample.fq')
                with io.open(frag_fq, 'w', encoding='ascii') as g:
                    for j, (rid, seq, qual) in enumerate(br.fastq_stream(path['fragment'])):
                        if j in frag_sample:
                            g.write(u'@%s\n%s\n+\n%s\n' % (rid, seq, qual))
            else:
                frag_sample = set(range(count['fragment']))
            frag_assignment = assign(mm2, a.threads, bin_centres_fa, frag_fq, a.identity, ASSIGN_COVERAGE, tiered=False)
            if not os.path.isdir(os.path.join(out, '_fragments')):
                os.makedirs(os.path.join(out, '_fragments'))
        # write
        with io.open(os.path.join(out, '%s_centres.fasta' % a.name), 'w', encoding='ascii') as g:
            for c in passing:
                g.write(u'>%s BIN%d reads=%d concatemer_segments=%d\n%s\n'
                        % (c, number[c], bin_count[c], seg_count[c], centres[c]))
        bins = {c: io.open(os.path.join(out, '%s-reads_BIN%d.fastq' % (a.name, number[c])), 'w', encoding='ascii')
                for c in passing}
        un = {k: io.open(os.path.join(out, '_unassigned', '%s_%s_unassigned.fastq' % (a.name, k)), 'w', encoding='ascii')
              for k in ('full', 'concatemer')}
        small = 0
        for rid, seq, qual in br.fastq_stream(path['full']):
            c = assignment.get(rid, (None, 0))[0]
            if c in bins:
                bins[c].write(u'@%s\n%s\n+\n%s\n' % (rid, seq, qual))
            elif c is not None:
                small += 1
                un['full'].write(u'@%s small_cluster\n%s\n+\n%s\n' % (rid, seq, qual))
            else:
                un['full'].write(u'@%s\n%s\n+\n%s\n' % (rid, seq, qual))
        cut = not_cut = 0
        for rid, seq, qual in br.fastq_stream(path['concatemer']):
            written = False
            for j, (qs, qe, c, _i) in enumerate(segs.get(rid) or [], 1):
                if c in bins:
                    bins[c].write(u'@%s_p%d concatemer_segment=%d-%d\n%s\n+\n%s\n'
                                  % (rid, j, qs, qe, seq[qs:qe], qual[qs:qe]))
                    written = True
            if written:
                cut += 1
            else:
                not_cut += 1
                un['concatemer'].write(u'@%s\n%s\n+\n%s\n' % (rid, seq, qual))
        frag_written = collections.Counter()
        if frag_sample:
            fw = {c: io.open(os.path.join(out, '_fragments', '%s-reads_BIN%d_fragments.fastq' % (a.name, number[c])),
                             'w', encoding='ascii') for c in passing}
            for j, (rid, seq, qual) in enumerate(br.fastq_stream(path['fragment'])):
                c = frag_assignment.get(rid, (None, 0))[0] if j in frag_sample else None
                if c in fw:
                    fw[c].write(u'@%s\n%s\n+\n%s\n' % (rid, seq, qual))
                    frag_written[c] += 1
            for f in fw.values():
                f.close()
        for f in list(bins.values()) + list(un.values()):
            f.close()
        table = [u'\t'.join([u'bin', u'reads', u'concatemer segments', u'fragments (separate file)',
                             u'share % (barcode)', u'centre read', u'centre bp', u'median identity', u'file'])]
        for c in passing:
            ids = sorted(idents[c]) or [0.0]
            table.append(u'\t'.join([u'BIN%d' % number[c], str(full_count[c]), str(seg_count[c]),
                                     str(frag_written[c]), u'%.2f' % (100.0 * bin_count[c] / total), c,
                                     str(len(centres[c])), u'%.4f' % ids[len(ids) // 2],
                                     u'%s-reads_BIN%d.fastq' % (a.name, number[c])]))
        full_unassigned = count['full'] - len(assignment)
        for label, n in ((u'full unassigned', full_unassigned), (u'small cluster (<%d)' % a.min_reads, small),
                         (u'concatemer cut', cut), (u'concatemer unassigned', not_cut),
                         (u'fragments (length class, total)', count['fragment']),
                         (u'fragments mapped (0 = not mapped, diagnostic only)', len(frag_sample)),
                         (u'<%d bp' % FRAGMENT_MIN, count['short'])):
            table.append(u'\t'.join([label, str(n), u'-', u'-', u'%.2f' % (100.0 * n / total), u'-', u'-', u'-', u'-']))
        io.open(os.path.join(out, 'BIN_TABLE.tsv'), 'w', encoding='utf-8', newline='\n').write(u'\n'.join(table) + u'\n')
        print('  bins %d (full %d + concatemer segments %d); full unassigned %d (%.1f%%), small clusters %d; '
              'concatemers cut %d / unassigned %d; written: %s'
              % (len(passing), sum(full_count[c] for c in passing), sum(seg_count[c] for c in passing),
                 full_unassigned, 100.0 * full_unassigned / max(1, count['full']), small, cut, not_cut, out))
        for c in passing[:8]:
            ids = sorted(idents[c]) or [0.0]
            print('    BIN%-3d full %6d  segments %5d  centre %d bp  identity %.3f'
                  % (number[c], full_count[c], seg_count[c], len(centres[c]), ids[len(ids) // 2]))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
