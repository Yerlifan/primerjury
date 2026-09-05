# -*- coding: utf-8 -*-
"""The read-level primer check and the ARMS trap (2026-09-05).

1. A product is found in both orientations, with at most one mismatch per
   primer and the two 3' bases exact; a mismatch in the 3' end is refused.
2. The verdict flags a SYSTEMATIC mismatch when the members amplify at mm1 but
   the pooled mm0 rate is zero, which is exactly what a deliberate -3 variant
   looks like on the reads.
3. The consensus chooser prefers read support over "fewest N" and marks a
   choice below the support floor as untrusted.

RUN
    python3 tests/test_primer_read_check.py
"""
from __future__ import print_function
import importlib.util
import os
import sys

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, KOK)


def _load(rel, name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(KOK, rel))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_products():
    m = _load('verification/read_level_primer_check.py', 'rlpc')
    F = 'AGCGCAACCCTCATCATTAG'
    R = 'GGCTGCTGGCACGGAGTT'
    insert = 'ACGT' * 30
    read = 'TTTT' + F + insert + m.rc(R) + 'GGGG'
    ok = m.products(read, F, R, 1) == (len(F) + len(insert) + len(R), 0)
    ok = ok and m.products(m.rc(read), F, R, 1) == (len(F) + len(insert) + len(R), 0)
    one = 'TTTT' + F[:5] + ('A' if F[5] != 'A' else 'C') + F[6:] + insert + m.rc(R) + 'GGGG'
    ok = ok and m.products(one, F, R, 1) == (len(F) + len(insert) + len(R), 1)
    three_prime = 'TTTT' + F[:-1] + ('A' if F[-1] != 'A' else 'C') + insert + m.rc(R) + 'GGGG'
    ok = ok and m.products(three_prime, F, R, 1) is None
    print('  products: %s' % ('ok' if ok else 'FAIL'))
    return ok


def test_arms_verdict():
    m = _load('verification/read_level_primer_check.py', 'rlpc')
    v, _, _, systematic = m.verdict([95.0, 90.0], 92.0, 0.0, 2.0, 1.0, False)
    ok = systematic and 'SYSTEMATIC MISMATCH' in v and v.startswith('PASSED')
    v2, _, _, s2 = m.verdict([95.0, 90.0], 92.0, 88.0, 2.0, 1.0, False)
    ok = ok and not s2 and v2 == 'PASSED (worst bin >=8x)'
    v3, _, _, _ = m.verdict([5.0], 5.0, 5.0, 0.0, 0.0, False)
    ok = ok and v3.startswith('NO / WEAK')
    print('  ARMS verdict: %s (%s)' % ('ok' if ok else 'FAIL', v))
    return ok


def test_consensus_choice():
    m = _load('steps/select_consensus.py', 'select_consensus')
    reads = ['ACGTTGCA' * 40]
    rk = set()
    for r in reads:
        rk |= m.kmers(r) | m.kmers(m.rc(r))
    right = 'ACGTTGCA' * 30                     # supported by the reads
    wrong = 'TTGACCAG' * 30                     # zero N, no support
    masked = ('ACGTTGCA' * 30)[:100] + 'N' * 20 + ('ACGTTGCA' * 30)[120:]
    best = m.choose([('zero_n', wrong), ('supported', right), ('masked', masked)], rk)
    ok = best is not None and best[0] == 'supported' and best[2] > m.MIN_SUPPORT
    best2 = m.choose([('zero_n', wrong)], rk)
    ok = ok and best2 is not None and best2[2] < m.MIN_SUPPORT
    print('  consensus choice: %s (%s %.0f%% / lone bad %.0f%%)'
          % ('ok' if ok else 'FAIL', best and best[0], best and best[2], best2 and best2[2]))
    return ok


if __name__ == '__main__':
    r = [test_products(), test_arms_verdict(), test_consensus_choice()]
    print('read-level primer check: %s' % ('PASS' if all(r) else 'FAIL'))
    sys.exit(0 if all(r) else 1)
