# -*- coding: utf-8 -*-
"""Function-gene primer design on an alignment (tools/function_gene_design.py, 2026-09-18).

1. The consensus keeps only columns filled in at least half of the sequences,
   and remembers which alignment column each consensus base came from.
2. A binding site is scored against the IUPAC primer: a 3' end mismatch or a
   site that is mostly gap is not testable (None), not "one more mismatch".
3. Degeneracy goes to the positions that win the most members first, stops at
   the variant cap and never touches the 3' end.
4. Coverage is counted per genus: products for sites within one mismatch,
   "testable" for any site that can be read at all.

Pure Python: neither primer3 nor MAFFT is needed.

RUN
    python3 tests/test_function_gene_design.py
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


m = _load('tools/function_gene_design.py', 'fgd')


def test_iupac():
    ok = m.rc('ACGK') == 'MCGT' and sorted(m.varyantlar('AKR')) == ['AGA', 'AGG', 'ATA', 'ATG']
    ok = ok and m.IKILI[frozenset('AG')] == 'R' and m.IKILI[frozenset('CT')] == 'Y'
    print('  IUPAC and reverse complement: %s' % ('ok' if ok else 'FAIL'))
    return ok


def test_consensus():
    kons, sut = m.uzlasi(['AC-T', 'ACGT', 'A--T'])
    ok = kons == 'ACT' and sut == [0, 1, 3]
    print('  consensus columns: %s (%s %s)' % ('ok' if ok else 'FAIL', kons, sut))
    return ok


def test_site_mismatch():
    p = 'ACGTACGTAC'
    ok = m.uyumsuzluk(p, p) == 0
    ok = ok and m.uyumsuzluk('TCGTACGTAC', p) == 1
    ok = ok and m.uyumsuzluk('ACGTACGTAG', p) is None          # 3' end
    ok = ok and m.uyumsuzluk('AC------AC', p) is None          # mostly gap
    ok = ok and m.uyumsuzluk('ACGTACGTAC', 'ACKTACGTAC') == 0  # K = G/T
    print('  site mismatch: %s' % ('ok' if ok else 'FAIL'))
    return ok


def test_degenerate():
    p = 'AAAAAAAAAA'
    uye = ['AAAAAAAAAA'] * 3 + ['GAAAAAAAAA'] * 2 + ['ACAAAAAAAA'] + ['AAAAAAAAAT'] * 4
    ok = m.dejenere_kur(p, uye, 2) == 'RAAAAAAAAA'
    ok = ok and m.dejenere_kur(p, uye, 4) == 'RMAAAAAAAA'
    ok = ok and m.dejenere_kur(p, [], 8) == p
    print('  degeneracy: %s' % ('ok' if ok else 'FAIL'))
    return ok


def test_coverage():
    F, yer_R = 'ACGTTGCAAC', 'GGATCCTTAG'
    R = m.rc(yer_R)
    iyi = F + 'TTTTT' + yer_R
    kotu = 'TTGTTGCAAC' + 'TTTTT' + yer_R                  # two mismatches in F
    kayit = [('a1', 'X', '', ''), ('a2', 'X', '', ''), ('b1', 'Y', '', '')]
    hiz = {'a1': iyi, 'a2': iyi, 'b1': kotu}
    say, boy = m.degerlendir(F, R, list(range(0, 10)), list(range(15, 25)), kayit, hiz)
    ok = say['X'] == [2, 2] and say['Y'] == [0, 1] and boy == [25, 25]
    print('  per-genus coverage: %s (%s %s)' % ('ok' if ok else 'FAIL', dict(say), boy))
    return ok


if __name__ == '__main__':
    r = [test_iupac(), test_consensus(), test_site_mismatch(), test_degenerate(), test_coverage()]
    print('function-gene design: %s' % ('PASS' if all(r) else 'FAIL'))
    sys.exit(0 if all(r) else 1)
