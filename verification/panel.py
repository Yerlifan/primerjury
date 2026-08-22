# -*- coding: utf-8 -*-
"""panel.py  -  TEK GIRIS NOKTASI.

NEDEN BIRLESTIRILDI
-------------------
Bir gecede sekiz ayri betik birikti (denetim, geometri, plaka, esik, guncel
durum, toplanti durumu, NCBI katmani, kutu uretimi). Sekiz ayri komut demek,
birini kosup otekini unutmak demek; bu projenin butun hatalari tam da boyle
cikti - bir sey degisti, ona bagli olan yeniden kosulmadi.

Burada hepsi TEK komutun alt komutlari. Kod KOPYALANMADI: her alt komut kendi
betigini cagirir, boylece iki kopya zamanla ayrisamaz.

Kullanim:
    python verification/panel.py hepsi          # sirasiyla hepsi (olcum, yazmaz)
    python verification/panel.py denetle        # 9 maddelik denetim kapisi
    python verification/panel.py geometri       # primer3 ile Tm ve kural denetimi
    python verification/panel.py geometri --write # Tm/dTm sutunlarini duzelt
    python verification/panel.py plaka          # jel cakismasi icin plaka onerisi
    python verification/panel.py esik           # iki esik kurali yan yana
    python verification/panel.py guncel         # GUNCEL_DURUM.md uret
    python verification/panel.py toplanti       # toplantida ne istendi, ne oldu
    python verification/panel.py ncbi4          # NCBI 4. katman (~30 dk)
    python verification/panel.py kutu --plan    # olculmeyen taksonlar icin plan
    python verification/panel.py referans       # hizli test referanslarini yenile

"hepsi" olcum yapar, HICBIR dosyayi degistirmez: geometri --write ve ncbi4 gibi
yan etkili adimlar hepsi'ye DAHIL DEGILDIR. Yan etkiyi insan ister.
"""
from __future__ import print_function

import os
import subprocess
import sys

BURA = os.path.dirname(os.path.abspath(__file__))
KOK = os.path.dirname(BURA)

# alt komut -> (betik yolu, aciklama, "hepsi"ye dahil mi)
KOMUT = {
    'denetle':  ('verification/audit_all.py',
                 u'18 maddelik denetim kapisi (tablolar, referans, muhur, '
                 u'geometri, plaka, belge sayilari)', True),
    'guncel':   ('verification/current_status.py',
                 u'GUNCEL_DURUM.md - panelin bugunku sayilari', True),
    'toplanti': ('verification/decision_status.py',
                 u'toplantida ne istendi, hangisi oldu', True),
    'esik':     ('verification/recompute_thresholds.py',
                 u'duz esik ile bolluga agirlikli esik yan yana', True),
    'geometri': ('verification/refresh_geometry.py',
                 u'primer3 ile Tm/GC/uzunluk/toka/dimer + urun boyu olcumu', True),
    'plaka':    ('verification/assign_plate.py',
                 u'plaka ici jel cakismasini azaltan dagilim onerisi', True),
    'referans': ('verification/refresh_reference.py',
                 u'hizli test referanslarini tam kosudan yenile', False),
    'ncbi4':    ('verification/ncbi_layer.py',
                 u'NCBI Primer-BLAST 4. katman (~30 dk, aga cikar)', False),
    'kutu':     ('verification/build_bins.py',
                 u'olculmeyen taksonlari kutuya cevir (kalibrasyon kapili)', False),
    'siparis':  ('verification/order_form.py',
                 u'tedarikciye gidecek TEK dogru oligo listesi (uretilir)', True),
    'excel':    ('verification/build_excel.py',
                 u'butun guncel veriyi TEK Excel e yaz (uretilir)', True),
    'lokus':    ('engine/target_full.py',
                 u'bir hedef icin BUTUN lokuslarda cift ara (primer3 gerekir)', False),
    'arsiv':    ('verification/archive.py',
                 u'eski dosyalari _SILINECEKLER e tasi (once PLAN basar)', False),
    'sinif':    ('verification/ncbi_reclassify.py',
                 u'NCBI sonuclarini siki ad kuraliyla yeniden say (aga cikmaz)', True),
    'kapsama':  ('screening/exclusion_coverage_check.py',
                 u'dislama taksonu uyeleri kapsiyor mu (NCBI Taxonomy)', False),
}

# "hepsi" sirasi: once durum uretilir, sonra denetim onu da gorur.
# "sabah" = hepsi + NCBI yeniden sayimi. Ikisi de OLCER, DEGISTIRMEZ.
SIRA = ['guncel', 'siparis', 'toplanti', 'esik', 'geometri', 'plaka',
        'sinif', 'excel', 'denetle']


def kos(ad, ek):
    yol = os.path.join(KOK, *KOMUT[ad][0].split('/'))
    if not os.path.exists(yol):
        print(u'  SKIPPED: %s is missing (%s)' % (ad, yol))
        return 127
    # target_full.py --root almaz; kokU _FL_KOK ortam degiskeninden okur.
    if ad == 'lokus':
        os.environ['_FL_KOK'] = KOK
        komut = [sys.executable, yol] + list(ek)
    else:
        komut = [sys.executable, yol, '--root', KOK] + list(ek)
    print()
    print('=' * 78)
    print('  >> %s   (%s)' % (ad, KOMUT[ad][1]))
    print('  $ %s' % ' '.join(komut[1:]))
    print('=' * 78)
    sys.stdout.flush()
    return subprocess.call(komut)


def yardim():
    print(__doc__)
    print('  alt komutlar:')
    for k in sorted(KOMUT):
        print('    %-10s %s%s' % (k, KOMUT[k][1],
                                  '' if KOMUT[k][2] else u'   [not included in \'all\']'))


def main(argv):
    if not argv or argv[0] in ('-h', '--help', 'yardim'):
        yardim()
        return 0
    ad, ek = argv[0], argv[1:]
    if ad in ('hepsi', 'sabah'):
        kodlar = {}
        for k in SIRA:
            kodlar[k] = kos(k, [])
        print()
        print('=' * 78)
        print(u'  ALL DONE')
        for k in SIRA:
            print('    %-10s cikis kodu %s' % (k, kodlar[k]))
        print()
        print('  Yan etkili adimlar BILEREK kosulmadi. Gerekiyorsa ayri ayri:')
        for k in sorted(KOMUT):
            if not KOMUT[k][2]:
                print('    python verification/panel.py %s' % k)
        print('=' * 78)
        return 1 if any(v not in (0, 2) for v in kodlar.values()) else 0
    if ad not in KOMUT:
        print('Bilinmeyen alt komut: %s' % ad)
        yardim()
        return 2
    return kos(ad, ek)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
