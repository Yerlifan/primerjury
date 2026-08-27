# -*- coding: utf-8 -*-
"""THE EVIDENCE for the unnamed environmental record trap, in the identity
verification's name resolver and its defensible level rule.

WHAT IS TESTED
--------------
Sets such as NCBI nt are full of unnamed records:
    KJ734864.1 Uncultured prokaryote clone D5 16S ribosomal RNA gene
    GQ503828.1 Bacterium enrichment culture clone R4-53B 16S ribosomal RNA

Before the fix the resolver took these for binomials ('Uncultured prokaryote',
'Bacterium enrichment'), and because the identity was 99 per cent the defensible
level rule produced a name at SPECIES level. The result: in the identity claims
table, having no answer stood there as a CONFIRMED identity.

Two things have to be true at once:
  1) an unnamed record MUST NOT PRODUCE A NAME (the false positive has to go)
  2) a real name MUST NOT BE BROKEN (no false negative may appear)

The second matters as much as the first: a filter that is too wide silences real
identities too, and that is the more insidious fault.

TO RUN IT
---------
    python3 tests/test_unnamed_records.py
"""
from __future__ import print_function
import argparse
import importlib.util
import os
import sys

VARSAYILAN_KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# (baslik, ad_uretmeli_mi, beklenen_tur)
ORNEKLER = [
    # --- ADSIZ: ad URETILMEMELI -------------------------------------------
    ('KJ734864.1 Uncultured prokaryote clone D5 16S ribosomal RNA gene, partial',
     False, None),
    ('KJ957653.1 Uncultured bacterium clone 4B-11 16S ribosomal RNA gene',
     False, None),
    ('GQ503828.1 Bacterium enrichment culture clone R4-53B 16S ribosomal RNA',
     False, None),
    ('MK123456.1 Unclassified organism isolate ABC 18S ribosomal RNA',
     False, None),
    ('XX000000.1 Environmental sample clone Z 16S ribosomal RNA',
     False, None),

    # --- GERCEK AD: BOZULMAMALI -------------------------------------------
    ('CP073276.1 Petrimonas sulfuriphila strain BN3 16S ribosomal RNA',
     True, 'Petrimonas sulfuriphila'),
    ('AF411468.1.1199 Archaea;Halobacteriota;Methanosarcinia;Methanosarcinales;'
     'Methanosarcina;Methanosarcina mazei',
     True, 'Methanosarcina mazei'),
    ('AY882347.1.484 Eukaryota;Amorphea;Obazoa;Opisthokonta;Fungi;Petriella;'
     'Petriella setifera',
     True, 'Petriella setifera'),
    ('NR_201921.1 Methanothermococcus jasoni strain Ax23 16S ribosomal RNA',
     True, 'Methanothermococcus jasoni'),
    # Cins adi 'Bacterium' ile BASLAMAYAN gercek bir ad susturulmamali
    ('NR_112345.1 Bacteroides fragilis strain NCTC 9343 16S ribosomal RNA',
     True, 'Bacteroides fragilis'),
]


def kd_yukle(kok):
    y = os.path.join(kok, 'verification', 'identity_verification.py')
    if not os.path.exists(y):
        sys.stderr.write(u'not found: %s\n' % y)
        sys.exit(2)
    spec = importlib.util.spec_from_file_location('kd', y)
    m = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(m)
    except SystemExit:
        pass
    return m


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', dest='kok', default=VARSAYILAN_KOK)
    a = p.parse_args()
    kd = kd_yukle(a.kok)

    gecti = True
    print(u'=== 1) ad_coz: an unnamed record MUST NOT PRODUCE a name, and a real name MUST NOT BE BROKEN ===')
    print('%-56s | %-24s | %s' % ('header', 'species resolved', 'result'))
    print('-' * 104)
    for baslik, ad_olmali, bek in ORNEKLER:
        _c, t, _tam = kd.ad_coz(baslik)
        ok = (bool(t) == ad_olmali) and (not ad_olmali or t == bek)
        gecti = gecti and ok
        print('%-56s | %-24s | %s'
              % (baslik[:56], str(t)[:24],
                 'RIGHT' if ok else 'WRONG (expected: %s)' % bek))

    print()
    print(u'=== 2) savunulabilir_duzey: 99 per cent identity MUST NOT PRODUCE a name (when the record is unnamed) ===')
    print('%-38s | %-38s | %s' % ('case', 'level', 'result'))
    print('-' * 104)
    senaryo = [
        ('an unnamed clone, identity 99 per cent',
         'KJ734864.1 Uncultured prokaryote clone D5 16S ribosomal RNA gene',
         99.0, False),
        ('a real name, identity 99 per cent',
         'CP073276.1 Petrimonas sulfuriphila strain BN3 16S ribosomal RNA',
         99.0, True),
    ]
    for ad, baslik, k, ad_olmali in senaryo:
        r = kd.savunulabilir_duzey(
            [dict(baslik=baslik, kimlik=k), dict(baslik='X 16S', kimlik=90.0)], 'SSU')
        duzey = r['duzey']
        adlandirildi = duzey in ('TUR', 'CINS') or duzey.startswith('CINS (')
        ok = (adlandirildi == ad_olmali)
        gecti = gecti and ok
        print('%-38s | %-38s | %s' % (ad, duzey[:38], 'RIGHT' if ok else 'WRONG'))

    print()
    print(u'RESULT: ' + (u'EVERY TEST PASSED' if gecti else 'BASARISIZ'))
    return 0 if gecti else 1


if __name__ == '__main__':
    sys.exit(main())
