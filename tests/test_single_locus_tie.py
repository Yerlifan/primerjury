# -*- coding: utf-8 -*-
"""The single-locus ladder (steps/target_identity.py): the same three faults the
multi-locus decision had, measured in the study on 2026-09-06.

1. 'Candidatus Methanofastidiosum' (a two word genus) was written as a SPECIES,
   because "has a space" stood in for "has a species epithet".
2. 'Petrimonas sp. IBARAKI' on top left a bin at genus while 'Petrimonas
   sulfuriphila' sat 0.10 points behind; the placeholder 'sp' even became the
   epithet of a 'cf.' name.
3. A record such as 'Methanosarcina sp. X' counted as a rival species and turned a
   clean 'Methanosarcina mazei' into 'cf.'.

RUN
    python3 tests/test_single_locus_tie.py
"""
from __future__ import print_function
import os
import sys

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(KOK, 'steps'))
sys.path.insert(0, os.path.join(KOK, 'verification'))
import target_identity as ti   # noqa: E402


def hit(pid, aln, title, db='archaea.16S.fna', slen=1500):
    return (pid, aln, title, db, slen)


def ladder(hits):
    return ti.basamaktan_sec({'16S': hits}, ['16S'], True)


def test_candidatus_genus_is_not_a_species():
    assert ti.cins_epitet(u'Candidatus Methanofastidiosum') == (u'Candidatus Methanofastidiosum', u'')
    assert ti.cins_epitet(u'Candidatus Methanofastidiosum methylthiophilus')[1] == u'methylthiophilus'
    r = ladder([hit(99.72, 1400, u'NR_3 Candidatus Methanofastidiosum strain X 16S ribosomal RNA')])
    assert r and r[0] == u'Candidatus Methanofastidiosum sp.', r
    return True


def test_placeholder_top_record_falls_to_species_within_margin():
    r = ladder([hit(100.0, 1450, u'NR_1 Petrimonas sp. IBARAKI 16S ribosomal RNA'),
                hit(99.9, 1450, u'NR_2 Petrimonas sulfuriphila strain BN3 16S ribosomal RNA')])
    assert r and r[0] == u'Petrimonas sulfuriphila', r
    return True


def test_placeholder_record_is_not_a_rival():
    r = ladder([hit(100.0, 1450, u'NR_4 Methanosarcina mazei strain S-6 16S ribosomal RNA'),
                hit(99.9, 1450, u'NR_5 Methanosarcina sp. X 16S ribosomal RNA')])
    assert r and r[0] == u'Methanosarcina mazei', r
    r = ladder([hit(100.0, 1450, u'NR_4 Methanosarcina mazei strain S-6 16S ribosomal RNA'),
                hit(99.9, 1450, u'NR_6 Methanosarcina soligelidi strain Y 16S ribosomal RNA')])
    assert r and r[0] == u'Methanosarcina cf. mazei', r
    return True


def test_genus_only_record_leading_by_more_than_margin_stays_genus():
    r = ladder([hit(100.0, 1450, u'NR_1 Petrimonas sp. IBARAKI 16S ribosomal RNA'),
                hit(99.3, 1450, u'NR_2 Petrimonas sulfuriphila strain BN3 16S ribosomal RNA')])
    assert r and r[0] == u'Petrimonas sp.', r
    return True


def main():
    tests = [test_candidatus_genus_is_not_a_species,
             test_placeholder_top_record_falls_to_species_within_margin,
             test_placeholder_record_is_not_a_rival,
             test_genus_only_record_leading_by_more_than_margin_stays_genus]
    failed = 0
    for t in tests:
        try:
            ok = t()
        except Exception as e:      # noqa: BLE001
            ok = False
            print('  %s raised %s: %s' % (t.__name__, type(e).__name__, e))
        print('  %-64s %s' % (t.__name__, 'PASS' if ok else 'FAIL'))
        failed += 0 if ok else 1
    print('  %d of %d tests passed' % (len(tests) - failed, len(tests)))
    return 0 if not failed else 1


if __name__ == '__main__':
    sys.exit(main())
