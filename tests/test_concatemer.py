# -*- coding: utf-8 -*-
"""A bin of ligated double amplicons is not an organism and must not be named.

Measured 2026-09-07 in the study's from-scratch SUP root: A1 bins with 2.9 kb reads
where every read carried two distinct 16S copies (Methanosarcina + Nitrosocosmicus).
The table had named such a bin from whichever half sorted first; 95 per cent of its
reads were the other organism. The rule lives in identity_verification and the
population driver applies it before polishing.

RUN
    python3 tests/test_concatemer.py
"""
from __future__ import print_function
import os
import sys

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(KOK, 'verification'))
from identity_verification import is_concatemer, AMPLICON_BP, CONCATEMER_FACTOR   # noqa: E402
import population_polish as pp                                                       # noqa: E402


def test_double_length_reads_are_a_concatemer():
    assert is_concatemer('A1', 2880) == (True, 1500)
    assert is_concatemer('B', 2980)[0]
    assert is_concatemer('F2', 7400)[0]
    return True


def test_normal_amplicons_are_not():
    for g, bp in (('A1', 1440), ('B', 1500), ('A2', 4300), ('F1', 3700), ('F2', 3695)):
        assert not is_concatemer(g, bp)[0], (g, bp)
    assert is_concatemer('ZZ', 99999) == (False, 0)      # unknown group: no rule, no claim
    return True


def test_factor_is_the_single_source():
    assert CONCATEMER_FACTOR == 1.6 and AMPLICON_BP['A1'] == 1500
    assert 'is_concatemer' in dir(pp)
    return True


def main():
    tests = [test_double_length_reads_are_a_concatemer, test_normal_amplicons_are_not, test_factor_is_the_single_source]
    failed = 0
    for t in tests:
        try:
            ok = t()
        except Exception as e:      # noqa: BLE001 - a test harness reports, it does not hide
            ok = False
            print('  %s raised %s: %s' % (t.__name__, type(e).__name__, e))
        print('  %-52s %s' % (t.__name__, 'PASS' if ok else 'FAIL'))
        failed += 0 if ok else 1
    print('  %d of %d tests passed' % (len(tests) - failed, len(tests)))
    return 0 if not failed else 1


if __name__ == '__main__':
    sys.exit(main())
