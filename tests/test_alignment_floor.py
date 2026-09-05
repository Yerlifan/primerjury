# -*- coding: utf-8 -*-
"""The coverage-aware alignment floor and the short-hit veto (2026-09-05).

1. hizalama_yeterli(): the floor passes on length, or on >= 250 bases covering
   >= 90 per cent of the record; 200 bases of a 210-base record do not pass.
2. locus_decision: a 72 bp hit at 100 per cent no longer vetoes the genus that
   a 1,400 bp hit at 99 per cent gives; six-tuples (with the record length) and
   five-tuples (without) are both accepted.
3. target_identity.esik_uygula(): a species over 520 of a 540-base ITS record
   keeps its species name; the same 520 bases of a 1,900-base record come down
   to genus.

RUN
    python3 tests/test_alignment_floor.py
"""
from __future__ import print_function
import os
import sys

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, KOK)
sys.path.insert(0, os.path.join(KOK, 'verification'))
sys.path.insert(0, os.path.join(KOK, 'steps'))


def test_floor():
    from verification import identity_verification as iv
    ok = iv.hizalama_yeterli('SSU', 1250, 1500)[0]
    ok = ok and iv.hizalama_yeterli('ITS', 520, 540)[0]          # coverage exception
    ok = ok and iv.hizalama_yeterli('ITS', 520, 540)[1] != u''   # and it says so
    ok = ok and not iv.hizalama_yeterli('ITS', 520, 1900)[0]     # 27 per cent of the record
    ok = ok and not iv.hizalama_yeterli('ITS', 200, 210)[0]      # under 250 bases
    ok = ok and not iv.hizalama_yeterli('SSU', 1000, None)[0]    # no record length, floor rules
    print('  floor: %s' % ('ok' if ok else 'FAIL'))
    return ok


def test_short_hit_does_not_veto():
    import locus_decision as ld
    import target_identity as ti
    hits6 = [
        (100.0, 72, 130.0, 1500, 'KX000001.1 Methanosarcina mazei strain X 16S', 72),
        (99.0, 1400, 2500.0, 1500, 'NR_041956.1 Methanosarcina mazei DSM 2053 16S', 1450),
    ]
    r = ld.lokus_karari(hits6, 'SSU', ti.ad_ayikla, ti.cins_epitet, ti.ADSIZ_JETONLARI)
    ok = r['ad'] != u'cannot be named' and r['hizalama'] == 1400
    hits5 = [h[:5] for h in hits6]
    r5 = ld.lokus_karari(hits5, 'SSU', ti.ad_ayikla, ti.cins_epitet, ti.ADSIZ_JETONLARI)
    ok = ok and r5['hizalama'] == 1400
    only_short = [hits6[0]]
    r0 = ld.lokus_karari(only_short, 'SSU', ti.ad_ayikla, ti.cins_epitet, ti.ADSIZ_JETONLARI)
    ok = ok and r0['ad'] == u'cannot be named'
    print('  short-hit veto: %s (%s / %d bp)' % ('ok' if ok else 'FAIL', r['ad'], r['hizalama']))
    return ok


def test_species_over_short_its_record():
    import target_identity as ti
    ad1, note1 = ti.esik_uygula('Petriella musispora', 99.8, 'ITS', 520, None, 540)
    ad2, note2 = ti.esik_uygula('Petriella musispora', 99.8, 'ITS', 520, None, 1900)
    ok = ad1 == 'Petriella musispora' and 'covers' in note1 and ad2 == 'Petriella sp.'
    print('  ITS coverage: %s (%s | %s)' % ('ok' if ok else 'FAIL', ad1, ad2))
    return ok


if __name__ == '__main__':
    r = [test_floor(), test_short_hit_does_not_veto(), test_species_over_short_its_record()]
    print('alignment floor: %s' % ('PASS' if all(r) else 'FAIL'))
    sys.exit(0 if all(r) else 1)
