# -*- coding: utf-8 -*-
"""ADSIZ CEVRE KAYDI tuzaginin KANITI (kimlik_dogrulama.ad_coz / savunulabilir_duzey).

NE SINANIYOR
------------
NCBI nt gibi kumeler adlandirilmamis kayitlarla doludur:
    KJ734864.1 Uncultured prokaryote clone D5 16S ribosomal RNA gene
    GQ503828.1 Bacterium enrichment culture clone R4-53B 16S ribosomal RNA

Duzeltme oncesi ad_coz bunlari ikili ad saniyordu ('Uncultured prokaryote',
'Bacterium enrichment') ve kimlik %99 oldugu icin savunulabilir_duzey TUR
duzeyinde bir ad uretiyordu. Sonuc: KIMLIK_SONUC/kimlik_iddialari.tsv icinde
CEVAPSIZLIK, DOGRULANDI damgali kimlik gibi duruyordu.

Iki sey ayni anda dogru olmali:
  1) adsiz kayit AD URETMEMELI  (yanlis pozitif gitmeli)
  2) gercek ad BOZULMAMALI      (yanlis negatif olusmamali)

Ikincisi birincisi kadar onemlidir: fazla genis bir suzgec gercek kimlikleri de
susturur ve bu daha sinsi bir hatadir.

KOSMA
-----
    python3 DENETIM_SINAMALARI/sina_adsiz_kayit.py
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
    y = os.path.join(kok, 'KURTARMA', 'kimlik_dogrulama.py')
    if not os.path.exists(y):
        sys.stderr.write('bulunamadi: %s\n' % y)
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
    p.add_argument('--kok', default=VARSAYILAN_KOK)
    a = p.parse_args()
    kd = kd_yukle(a.kok)

    gecti = True
    print('=== 1) ad_coz: adsiz kayit ad URETMEMELI, gercek ad BOZULMAMALI ===')
    print('%-56s | %-24s | %s' % ('baslik', 'cozulen tur', 'sonuc'))
    print('-' * 104)
    for baslik, ad_olmali, bek in ORNEKLER:
        _c, t, _tam = kd.ad_coz(baslik)
        ok = (bool(t) == ad_olmali) and (not ad_olmali or t == bek)
        gecti = gecti and ok
        print('%-56s | %-24s | %s'
              % (baslik[:56], str(t)[:24],
                 'DOGRU' if ok else 'YANLIS (beklenen: %s)' % bek))

    print()
    print('=== 2) savunulabilir_duzey: %99 kimlik AD URETMEMELI (kayit adsizsa) ===')
    print('%-30s | %-38s | %s' % ('durum', 'duzey', 'sonuc'))
    print('-' * 104)
    senaryo = [
        ('adsiz klon, kimlik %99',
         'KJ734864.1 Uncultured prokaryote clone D5 16S ribosomal RNA gene',
         99.0, False),
        ('gercek ad, kimlik %99',
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
        print('%-30s | %-38s | %s' % (ad, duzey[:38], 'DOGRU' if ok else 'YANLIS'))

    print()
    print('SONUC: ' + ('BUTUN SINAMALAR GECTI' if gecti else 'BASARISIZ'))
    return 0 if gecti else 1


if __name__ == '__main__':
    sys.exit(main())
