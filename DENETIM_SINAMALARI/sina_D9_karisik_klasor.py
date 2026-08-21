# -*- coding: utf-8 -*-
"""M5-D9 (karisik yonlu konsensus klasoru) kontrolunun KANITI.

NE SINANIYOR
------------
`capraz_kontrol.d9_karisik_klasor_yollari()` sunu dogru cevaplamali:
"bu dosya karisik yonlu 'consensus sequences' klasorunu KOD ICINDE okuyor mu?"

Bu soru onemli cunku o klasor karisik yonludur (olculen: 71 antisense /
27 sense) ve ters yonlu bir konsensuste in-silico PCR SESSIZCE 0 urun verir;
olculen kayip %100 (KAPSAMLI_ARAMA/yon_etki_testi.py).

NEDEN AYRI BIR SINAMA GEREKIYOR
-------------------------------
Kontrol eskiden dosyanin TAMAMINDA duz metin aramasi yapiyordu. 2026-08-09
kosusunda bes yanlis pozitif uretti: aciklamada klasor adini ANAN dosyalar da
RISKLI isaretlendi. Bir denetleyicinin yanlis pozitifi zararsiz degildir -
gercek bulgulari gurultuye gomer ve rapora guveni dusurur.

2026-08-21 duzeltmesi uc suzgec ekledi:
  1) docstring ve yorumlar atilir      -> aciklamada gecen ad kod sayilmaz
  2) cikti cagrilarinin argumanlari atilir -> ekrana basilan mesaj kod degil
  3) gorevi geregi okuyanlar muaf      -> D9_MESRU_OKUYUCULAR

KOSMA
-----
    python3 DENETIM_SINAMALARI/sina_D9_karisik_klasor.py
    python3 DENETIM_SINAMALARI/sina_D9_karisik_klasor.py --kok /baska/yol

Cikis kodu 0 = yedi sinamanin yedisi de beklendigi gibi.
Yol GOMULU DEGILDIR: varsayilan olarak bu dosyanin bir ust dizini kok sayilir.
"""
from __future__ import print_function
import argparse
import importlib.util
import io
import os
import sys

VARSAYILAN_KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

KARISIK = u'consensus sequences'
KANONIK = u'konsensus_kanonik'

# --- sentetik ornekler: her biri TEK bir ayrimi sinar --------------------
GERCEK_RISK = u'''
import os
def yukle(kok):
    return open(os.path.join(kok, "consensus sequences", "A1-1.fasta")).read()
'''

SADECE_ACIKLAMA = u'''
"""Bu betik eskiden "consensus sequences" klasorunu okurdu."""
# consensus sequences artik kullanilmiyor
import os
def yukle(kok):
    return open(os.path.join(kok, "kanonik", "A1-1.fasta")).read()
'''

SADECE_MESAJ = u'''
import os
def kos(yaz, kok):
    yaz("  Kaynak: \\"consensus sequences\\" (eski set).")
    return open(os.path.join(kok, "kanonik", "A1-1.fasta")).read()
'''

# --- gercek dosyalar: 2026-08-09'da yanlis pozitif verenler --------------
GERCEK_DOSYALAR = [
    (os.path.join('KAPSAMLI_ARAMA', 'hepsi.py'), False,
     u'yalniz ekrana basilan mesaj dizesi'),
    (os.path.join('KAPSAMLI_ARAMA', 'yon_denetimi.py'), False,
     u'gorevi geregi okuyor (muaf liste)'),
    (os.path.join('WSL_betikleri', '03_primer_aday_uret.py'), False,
     u'yol CLI argumani, docstring ornegi'),
    (os.path.join('WSL_betikleri', '04_grup_primer.py'), False,
     u'yol CLI argumani, docstring ornegi'),
]


def ck_yukle(kok):
    yol = os.path.join(kok, 'capraz_kontrol.py')
    if not os.path.exists(yol):
        sys.stderr.write('capraz_kontrol.py bulunamadi: %s\n' % yol)
        sys.exit(2)
    spec = importlib.util.spec_from_file_location('ck', yol)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--kok', default=VARSAYILAN_KOK)
    a = p.parse_args()
    ck = ck_yukle(a.kok)

    def riskli(metin, ad=u''):
        v = ck.d9_karisik_klasor_yollari(metin, ad)
        return bool(v) and KANONIK not in ck._kod_govdesi(metin)

    def duz_metin(metin):          # duzeltme ONCESI mantik, karsilastirma icin
        return KARISIK in metin and KANONIK not in metin

    gecti = True
    print(u'%-42s | %-8s | %-8s | %-8s | %s'
          % (u'sinama', u'ESKI', u'YENI', u'beklenen', u'sonuc'))
    print(u'-' * 100)

    for rel, beklenen, gerekce in GERCEK_DOSYALAR:
        yol = os.path.join(a.kok, rel)
        if not os.path.exists(yol):
            print(u'%-42s | DOSYA YOK - sinama kosulamadi' % rel)
            gecti = False
            continue
        m = io.open(yol, encoding='utf-8', errors='replace').read()
        y = riskli(m, yol)
        ok = (y == beklenen)
        gecti = gecti and ok
        print(u'%-42s | %-8s | %-8s | %-8s | %-6s (%s)'
              % (rel.replace(os.sep, '/'),
                 u'RISKLI' if duz_metin(m) else u'temiz',
                 u'RISKLI' if y else u'temiz',
                 u'RISKLI' if beklenen else u'temiz',
                 u'DOGRU' if ok else u'YANLIS', gerekce))

    print()
    for ad, metin, beklenen in [
            (u'sentetik: KODDA gercekten okuyor', GERCEK_RISK, True),
            (u'sentetik: yalniz aciklamada aniyor', SADECE_ACIKLAMA, False),
            (u'sentetik: yalniz ekrana basiyor', SADECE_MESAJ, False)]:
        y = riskli(metin, u'sentetik.py')
        ok = (y == beklenen)
        gecti = gecti and ok
        print(u'%-42s | %-8s | %-8s | %-8s | %s'
              % (ad, u'-', u'RISKLI' if y else u'temiz',
                 u'RISKLI' if beklenen else u'temiz',
                 u'DOGRU' if ok else u'YANLIS'))

    print()
    print(u'SONUC: ' + (u'YEDI SINAMANIN YEDISI DE BEKLENDIGI GIBI'
                        if gecti else u'BASARISIZ'))
    return 0 if gecti else 1


if __name__ == '__main__':
    sys.exit(main())
