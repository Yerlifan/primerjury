# -*- coding: utf-8 -*-
"""The top record cannot always name a species; the decision must not stop there.

Measured 2026-09-06 on a fungal bin: the ITS search put a UNITE "Petriella sp."
record at 100.00 per cent over 600 bp on top and the RefSeq TYPE record
"Petriella musispora" at 100.00 per cent over 500 bp right under it. The old
rule looked only at the top record and left the bin at genus. Now, when the top
record cannot name a species (unnamed, or genus only) and a species-named
record sits within the separation margin, the decision rests on that record
and the skipped one is written in the note. A genus-only record that leads by
more than the margin still means "the nearest reference is an undescribed
lineage" and the bin stays at genus.

RUN
    python3 tests/test_locus_tie.py
"""
from __future__ import print_function
import os
import sys

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(KOK, 'verification'))
sys.path.insert(0, os.path.join(KOK, 'steps'))
import locus_decision as ld                                   # noqa: E402
from identity_verification import ADSIZ_JETONLARI, AYRIM_PAYI  # noqa: E402
from target_identity import ad_ayikla, cins_epitet             # noqa: E402


def hit(pid, aln, title, slen=None):
    return (pid, aln, 1000.0, 3000, title, slen or aln)


def decide(hits, key='ITS'):
    return ld.lokus_karari(hits, key, ad_ayikla, cins_epitet, ADSIZ_JETONLARI)


def test_genus_only_top_record_within_margin_yields_species():
    r = decide([hit(100.0, 600, u'UDB1|k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Microascales;f__Microascaceae;g__Petriella;s__Petriella_sp|SH1'),
                hit(100.0, 500, u'NR_172285.1 Petriella musispora CBS 745.69 ITS region; from TYPE material'),
                hit(92.0, 500, u'NR_156504.1 Petriella guttulata CBS 362.61 ITS region')])
    assert r['duzey'] == u'tur' and r['ad'] == u'Petriella musispora', r
    assert r['notu'], r
    return True


def test_unnamed_top_record_within_margin_yields_species():
    r = decide([hit(99.9, 550, u'MN1 uncultured fungus clone X ITS'),
                hit(99.8, 500, u'NR_1 Petriella musispora CBS 1 ITS region')])
    assert r['duzey'] == u'tur' and r['ad'] == u'Petriella musispora', r
    return True


def test_genus_only_record_leading_by_more_than_margin_stays_genus():
    r = decide([hit(100.0, 600, u'UDB2|k__Fungi;g__Petriella;s__Petriella_sp|SH2'),
                hit(100.0 - AYRIM_PAYI - 0.1, 500, u'NR_1 Petriella musispora CBS 1 ITS region')])
    assert r['duzey'] == u'cins' and r['ad'] == u'Petriella sp.', r
    return True


def test_close_rival_species_still_gives_cf():
    r = decide([hit(100.0, 500, u'NR_1 Petriella musispora CBS 1 ITS region'),
                hit(99.8, 500, u'NR_2 Petriella setifera CBS 2 ITS region')])
    assert r['duzey'] == u'cins' and r['ad'].startswith(u'Petriella cf.'), r
    return True


def test_placeholder_epithet_is_not_a_species_and_not_a_rival():
    r = decide([hit(100.0, 600, u'NG_1 Petriella sp. CBS 3 28S rRNA')], 'LSU_MANTAR')
    assert r['duzey'] == u'cins' and r['ad'] == u'Petriella sp.', r
    r = decide([hit(100.0, 500, u'NR_1 Petriella musispora CBS 1 ITS region'),
                hit(99.9, 600, u'NG_2 Petriella sp. CBS 9 ITS region')])
    assert r['duzey'] == u'tur' and r['ad'] == u'Petriella musispora', r
    return True


def main():
    tests = [test_genus_only_top_record_within_margin_yields_species,
             test_unnamed_top_record_within_margin_yields_species,
             test_genus_only_record_leading_by_more_than_margin_stays_genus,
             test_close_rival_species_still_gives_cf,
             test_placeholder_epithet_is_not_a_species_and_not_a_rival]
    failed = 0
    for t in tests:
        try:
            ok = t()
        except Exception as e:      # noqa: BLE001
            ok = False
            print('  %s raised %s: %s' % (t.__name__, type(e).__name__, e))
        print('  %-60s %s' % (t.__name__, 'PASS' if ok else 'FAIL'))
        failed += 0 if ok else 1
    print('  %d of %d tests passed' % (len(tests) - failed, len(tests)))
    return 0 if not failed else 1


if __name__ == '__main__':
    sys.exit(main())
