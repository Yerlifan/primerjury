# -*- coding: utf-8 -*-
"""The pure parts of the Kraken-free binning step (2026-09-05).

Two things that must stay true whatever minimap2 does:

1. Length peaks are non-overlapping windows of peak +- 15 per cent, ordered by
   share, and a peak below the share floor is not a window.
2. A PAF pair is accepted only when the alignment covers >= 80 per cent of the
   SHORTER sequence and 1 - de reaches the identity; the coverage is measured
   on the shorter one on purpose (a 4 kb read covering a 1.4 kb centre is a
   match for the centre, not a 35 per cent match).

RUN
    python3 tests/test_bin_reads.py
"""
from __future__ import print_function
import importlib.util
import os
import sys

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load():
    path = os.path.join(KOK, 'steps', 'bin_reads.py')
    spec = importlib.util.spec_from_file_location('bin_reads', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_peaks():
    m = _load()
    # 1300 bp would fall inside the 1450 window (+-15 per cent) and is folded into it;
    # 800 bp stands apart; 500 bp is 2 per cent and below the floor.
    lengths = [1450] * 540 + [3700] * 300 + [800] * 120 + [500] * 20 + [200] * 20
    w = m.find_length_peaks(lengths, min_share=0.03, max_peaks=4)
    ok = ([x[2] for x in w] == [1450, 3750, 850]            # by share
          and w[0][0] == int(1450 * 0.85) and w[0][1] == int(1450 * 1.15)
          and all(not (a[1] >= b[0] and b[1] >= a[0]) for i, a in enumerate(w) for b in w[i + 1:]))
    print('  peaks: %s %s' % ('ok' if ok else 'FAIL', [(x[2], round(x[3], 1)) for x in w]))
    return ok


def test_paf_pairs():
    m = _load()
    def row(q, ql, qs, qe, t, tl, ts, te, de):
        return '\t'.join(map(str, [q, ql, qs, qe, '+', t, tl, ts, te, 1000, 1200, 60, 'de:f:%s' % de]))
    paf = '\n'.join([
        row('r1', 4000, 100, 3900, 'c1', 1400, 10, 1390, 0.02),   # long read over short centre: covered
        row('r2', 1400, 0, 700, 'c1', 1400, 0, 700, 0.01),        # half covered: refused
        row('r3', 1400, 0, 1400, 'c1', 1400, 0, 1400, 0.05),      # identity 0.95 < 0.97: refused
        row('c1', 1400, 0, 1400, 'c1', 1400, 0, 1400, 0.0),       # self hit: ignored
        row('r4', 1400, 0, 1400, 'c2', 1400, 0, 1400, 0.03) + '\tcm:i:5',
    ])
    got = sorted((a, b, round(i, 2)) for a, b, i in m.paf_pairs(paf, 0.97, 0.80))
    ok = got == [('r1', 'c1', 0.98), ('r4', 'c2', 0.97)]
    print('  paf pairs: %s %s' % ('ok' if ok else 'FAIL', got))
    return ok


if __name__ == '__main__':
    r = [test_peaks(), test_paf_pairs()]
    print('bin_reads: %s' % ('PASS' if all(r) else 'FAIL'))
    sys.exit(0 if all(r) else 1)
