# -*- coding: utf-8 -*-
"""The pure parts of the clean binning step (2026-09-08) and, when minimap2 is present, one end-to-end
run on synthetic reads.

What must stay true:

1. Length classes: <300 short; <0.80x fragment; 0.80-1.25x full; >1.25x concatemer.
2. Tiered assignment: >= 0.97 goes in; 0.95-0.97 only with a >= 0.02 margin over the second centre;
   coverage below 0.85 of the read refuses; without tiering the second tier is off.
3. Concatemer segments: intervals must cover >= 80 per cent of the centre, be >= 0.97 identical,
   be no longer than 1.25x the centre (a two-copy chain is refused) and not overlap; longer and more
   identical first.
4. Satellite rule (pure part): half the sample on one chosen centre at >= 0.97 over >= 80 per cent of
   the centre makes a satellite; spread hits do not.
5. End to end (needs minimap2): two organisms at 90 per cent identity, reads with 0.7 per cent error (SUP),
   tailed reads, tandem concatemers and fragments -> two bins, one organism each, every read of a bin
   within 0.80-1.25x, concatemers cut into segments, no segment longer than 1.25x, fragments counted and
   not in any bin, unassigned reads written and counted; then select_bins drops nothing wrongly and
   filters long reads.

RUN
    python3 tests/test_clean_bins.py
"""
from __future__ import print_function
import importlib.util
import io
import os
import random
import shutil
import subprocess
import sys
import tempfile

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load(name):
    path = os.path.join(KOK, 'steps', name + '.py')
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def row(q, ql, qs, qe, t, tl, ts, te, de):
    return '\t'.join(map(str, [q, ql, qs, qe, '+', t, tl, ts, te, 1000, 1200, 60, 'de:f:%s' % de]))


def test_length_class():
    m = _load('clean_bins')
    got = [m.length_class(L, 1500) for L in (250, 600, 1199, 1200, 1500, 1875, 1876, 2900)]
    want = ['short', 'fragment', 'fragment', 'full', 'full', 'full', 'concatemer', 'concatemer']
    ok = got == want
    print('  length classes: %s %s' % ('ok' if ok else 'FAIL', got))
    return ok


def test_tiered_assignment():
    m = _load('clean_bins')
    paf = '\n'.join([
        row('direct', 1500, 0, 1490, 'c1', 1500, 0, 1490, 0.02),          # 0.98 -> in
        row('margin', 1500, 0, 1480, 'c1', 1500, 0, 1480, 0.04),          # 0.96 vs 0.93: margin 0.03 -> in
        row('margin', 1500, 0, 1480, 'c2', 1500, 0, 1480, 0.07),
        row('close', 1500, 0, 1480, 'c1', 1500, 0, 1480, 0.04),           # 0.96 vs 0.95: margin 0.01 -> out
        row('close', 1500, 0, 1480, 'c2', 1500, 0, 1480, 0.05),
        row('lowcov', 1500, 0, 1200, 'c1', 1500, 0, 1200, 0.01),          # 80 per cent covered -> out
        row('weak', 1500, 0, 1490, 'c1', 1500, 0, 1490, 0.06),            # 0.94 -> out
    ])
    got = m.assign_from_rows(m.paf_rows(paf), 0.97, 0.85)
    ok1 = sorted(got) == ['direct', 'margin'] and got['margin'][0] == 'c1'
    flat = m.assign_from_rows(m.paf_rows(paf), 0.97, 0.85, tiered=False)
    ok2 = sorted(flat) == ['direct']
    ok = ok1 and ok2
    print('  tiered assignment: %s %s / flat %s' % ('ok' if ok else 'FAIL', sorted(got), sorted(flat)))
    return ok


def test_segments():
    m = _load('clean_bins')
    paf = '\n'.join([
        row('two', 3000, 0, 1450, 'c1', 1500, 0, 1450, 0.02),             # copy 1
        row('two', 3000, 1500, 2980, 'c1', 1500, 10, 1490, 0.02),         # copy 2
        row('two', 3000, 100, 1300, 'c2', 1500, 0, 1200, 0.01),           # overlaps copy 1, 80% -> refused by overlap
        row('chain', 3084, 0, 3084, 'c1', 1500, 0, 1430, 0.02),           # 2.06x the centre: two-copy chain -> refused
        row('chain', 3084, 0, 1437, 'c1', 1500, 0, 1429, 0.03),           # the honest half -> taken
        row('short', 3000, 0, 1000, 'c1', 1500, 0, 1000, 0.01),           # covers 67% of the centre -> refused
    ])
    got = m.segments_from_rows(m.paf_rows(paf), 0.97)
    ok = (got.get('two') == [(0, 1450, 'c1', 0.98), (1500, 2980, 'c1', 0.98)]
          and got.get('chain') == [(0, 1437, 'c1', 0.97)] and 'short' not in got)
    print('  segments: %s %s' % ('ok' if ok else 'FAIL', got))
    return ok


def test_satellite_rule():
    m = _load('select_bins')
    rows_sat = [(u'r%d' % i, 1750, 0, 1440, 'BIN1', 1430, 0, 1429, 0.99) for i in range(31)] + \
               [(u'r%d' % i, 1750, 0, 1440, 'BIN2', 1430, 0, 1429, 0.99) for i in range(31, 60)]
    sat = m.satellite_from_rows(rows_sat, 0.97, 0.8, 60)
    rows_spread = [(u'r%d' % i, 1750, 0, 1440, 'BIN%d' % (i % 4), 1430, 0, 1429, 0.99) for i in range(60)]
    spread = m.satellite_from_rows(rows_spread, 0.97, 0.8, 60)
    rows_low = [(u'r%d' % i, 1750, 0, 1440, 'BIN1', 1430, 0, 1000, 0.99) for i in range(60)]   # 70% of centre
    low = m.satellite_from_rows(rows_low, 0.97, 0.8, 60)
    ok = sat == ('BIN1', 31) and spread is None and low is None
    print('  satellite rule: %s %s %s %s' % ('ok' if ok else 'FAIL', sat, spread, low))
    return ok


def _mutate(seq, rate, rnd):
    out = []
    for b in seq:
        r = rnd.random()
        if r < rate / 3:
            continue                                        # deletion
        if r < 2 * rate / 3:
            out.append(rnd.choice('ACGT'))                  # substitution
            out.append(b) if rnd.random() < 0.5 else None
        elif r < rate:
            out.append(b)
            out.append(rnd.choice('ACGT'))                  # insertion
        else:
            out.append(b)
    return ''.join(out)


def test_end_to_end():
    if not shutil.which('minimap2'):
        print('  end to end: skipped (minimap2 not on PATH)')
        return True
    m = _load('clean_bins')
    rnd = random.Random(3)
    org1 = ''.join(rnd.choice('ACGT') for _ in range(1500))
    org2 = list(org1)                                       # 90 per cent identical sister
    for i in rnd.sample(range(1500), 150):
        org2[i] = rnd.choice([b for b in 'ACGT' if b != org2[i]])
    org2 = ''.join(org2)
    tail = ''.join(rnd.choice('ACGT') for _ in range(310))
    reads = []
    for k in range(400):
        reads.append(('o1_%d' % k, _mutate(org1, 0.007, rnd)))
    for k in range(200):
        reads.append(('o2_%d' % k, _mutate(org2, 0.007, rnd)))
    for k in range(40):
        reads.append(('tail_%d' % k, _mutate(org1, 0.007, rnd) + tail))               # 1.2x, still full class
    for k in range(60):
        reads.append(('cat_%d' % k, _mutate(org1, 0.007, rnd) + _mutate(org1, 0.007, rnd)))   # tandem concatemer
    for k in range(80):
        reads.append(('frag_%d' % k, _mutate(org1[200:800], 0.007, rnd)))              # 600 bp fragment
    for k in range(20):
        reads.append(('junk_%d' % k, ''.join(rnd.choice('ACGT') for _ in range(1500))))   # nothing
    rnd.shuffle(reads)
    tmp = tempfile.mkdtemp(prefix='clean_bins_test_')
    try:
        fq = os.path.join(tmp, 'barcode99.fastq')
        with io.open(fq, 'w', encoding='ascii') as g:
            for rid, seq in reads:
                g.write(u'@%s\n%s\n+\n%s\n' % (rid, seq, 'I' * len(seq)))
        out = os.path.join(tmp, 'bins')
        rc = m.main(['--fastq', fq, '--name', 'A1-9', '--out', out, '--threads', '2', '--seeds', '200'])
        binroot = os.path.join(out, 'A1-9')
        table = [l.rstrip('\n').split('\t') for l in io.open(os.path.join(binroot, 'BIN_TABLE.tsv'), encoding='utf-8')]
        bins = [r for r in table[1:] if r[0].startswith('BIN')]
        members = {}
        for r in bins:
            for rid, seq, _q in m._bin_reads().fastq_stream(os.path.join(binroot, r[8])):
                members.setdefault(r[0], []).append((rid, len(seq)))
        # one organism per bin: the o1/tail/cat prefixes belong together, o2 alone, no frag/junk in any bin
        pure = all(len(set(('o2' if rid.startswith('o2') else 'o1') for rid, _ in L)) == 1 for L in members.values())
        no_frag = all(not rid.startswith(('frag', 'junk')) for L in members.values() for rid, _ in L)
        in_class = all(1200 <= L <= 1875 for M in members.values() for _r, L in M)
        cut = any(rid.startswith('cat_') and '_p' in rid for L in members.values() for rid, _ in L)
        two_big = len([r for r in bins if int(r[1]) + int(r[2]) >= 100]) == 2
        frag_row = [r for r in table if r[0].startswith('fragments (length class')][0]
        unassigned = os.path.isfile(os.path.join(binroot, '_unassigned', 'A1-9_full_unassigned.fastq'))
        ok = rc == 0 and pure and no_frag and in_class and cut and two_big and int(frag_row[1]) == 80 and unassigned
        print('  end to end: %s bins %d (>=100 reads: %d), pure %s, no fragments in bins %s, in class %s, '
              'concatemers cut %s, fragments counted %s' % ('ok' if ok else 'FAIL', len(bins),
              sum(1 for r in bins if int(r[1]) + int(r[2]) >= 100), pure, no_frag, in_class, cut, frag_row[1]))
        # selection on top: the long-read filter and the table columns
        s = _load('select_bins')
        target = os.path.join(tmp, 'chosen')
        rc2 = s.main(['--source', out, '--target', target, '--name', 'A1-9'])
        sel = [l.rstrip('\n').split('\t') for l in io.open(os.path.join(target, 'SELECTION_TABLE.tsv'), encoding='utf-8')]
        chosen = [r for r in sel[1:] if r[5] == 'YES']
        copied = sum(1 for r in chosen for _ in m._bin_reads().fastq_stream(os.path.join(target, 'A1-9', 'A1-9-reads_%s.fastq' % r[1])))
        # a bin made of the tailed reads of organism 1 (round two seeds from the leftovers) must be
        # flagged as BIN1's satellite, not chosen
        tail_bins = [b for b, L in members.items() if all(rid.startswith('tail_') for rid, _ in L)]
        sat_rows = [r for r in sel[1:] if r[1] in tail_bins]
        sat_ok = all(r[5] == 'no' and r[6].startswith('satellite of') for r in sat_rows)
        ok2 = rc2 == 0 and len(chosen) >= 2 and copied > 0 and sel[0][-1].startswith('long dropped') and sat_ok
        print('  selection: %s chosen %d, copied reads %d, tail bins %s flagged as satellites %s'
              % ('ok' if ok2 else 'FAIL', len(chosen), copied, tail_bins, sat_ok))
        return ok and ok2
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    print('test_clean_bins')
    results = [test_length_class(), test_tiered_assignment(), test_segments(), test_satellite_rule(),
               test_end_to_end()]
    ok = all(results)
    print('  %s: %d/%d' % ('PASS' if ok else 'FAIL', sum(results), len(results)))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
