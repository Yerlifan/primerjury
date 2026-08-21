# -*- coding: utf-8 -*-
"""GUNCEL_DURUM.md  -  panelin BUGUNKU sayilari, TEK yerden ve URETILEREK.

NEDEN VAR
---------
Ayni soruya uc belge uc farkli cevap veriyordu: OKU_ONCE "KESIN 16 cift",
NASIL_DEVAM_EDILIR "KESIN 16 cift", CALISTIRMA_KILAVUZU "siparis edilebilir
11 cift". Panelin bugunku hali 20 cift (15 hedef ozgul + 5 evrensel). Ucu de
yazildiklari gunde dogruydu; sorun sayinin belgeye ELLE yazilmis olmasi.

Bu dosya elle yazilmaz. Her kosuda panelin kendi ciktisindan uretilir. Belgeler
sayiyi tekrar etmek yerine buraya isaret eder.

Kosum:
    python verification/current_status.py --kok .
"""
from __future__ import print_function

import argparse
import csv
import io
import os
import sys
import time


def _tsv(yol):
    if not os.path.exists(yol):
        return []
    with io.open(yol, encoding='utf-8') as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith('#')),
                                   delimiter='\t'))


def sayilar(kok):
    sl = _tsv(os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'SIPARIS_LISTESI.tsv'))
    kesin = [r for r in sl if (r.get('SINIF') or '').strip().upper() == 'KESIN']
    evr = [r for r in sl if (r.get('SINIF') or '').strip().upper() == 'EVRENSEL']
    disi = [r for r in sl if (r.get('SINIF') or '').strip().upper()
            not in ('KESIN', 'EVRENSEL')]
    return sl, kesin, evr, disi


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', '--kok', dest='kok', default='.')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)
    sl, kesin, evr, disi = sayilar(kok)
    if not sl:
        sys.exit('HATA: TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv okunamadi.')

    y = os.path.join(kok, 'GUNCEL_DURUM.md')
    with io.open(y, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# Panelin bugünkü hâli\n\n')
        fh.write(u'Üretim: %s · kaynak `TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv`\n\n'
                 % time.strftime('%Y-%m-%d %H:%M'))
        fh.write(u'> Bu dosya **elle yazılmaz**, her koşuda panelin kendi çıktısından '
                 u'üretilir. Başka bir belgede bununla çelişen bir sayı görürseniz '
                 u'doğru olan budur; o belge yazıldığı günün sayısını taşıyordur.\n\n')
        fh.write(u'| | çift | oligo |\n|---|---|---|\n')
        fh.write(u'| **Sipariş edilecek** | **%d** | **%d** |\n'
                 % (len(kesin) + len(evr), 2 * (len(kesin) + len(evr))))
        fh.write(u'| bunun hedef özgül olanı | %d | %d |\n' % (len(kesin), 2 * len(kesin)))
        fh.write(u'| bunun evrensel/kontrol olanı | %d | %d |\n' % (len(evr), 2 * len(evr)))
        fh.write(u'| Sipariş dışı (silinmedi) | %d | — |\n\n' % len(disi))
        fh.write(u'## Sipariş edilecek çiftler\n\n')
        fh.write(u'| hedef | sınıf | şart | ürün (bp) | dCq |\n|---|---|---|---|---|\n')
        for r in kesin + evr:
            fh.write(u'| %s | %s | %s | %s | %s |\n'
                     % (r.get('hedef', ''), (r.get('SINIF') or '').strip(),
                        (r.get('siparis_sarti') or '—').strip(),
                        (r.get('urun_bp') or '—').strip(),
                        (r.get('dCq_karsiligi') or '—').strip()))
        if disi:
            fh.write(u'\n## Sipariş dışı\n\n| hedef | durum | dCq |\n|---|---|---|\n')
            for r in disi:
                fh.write(u'| %s | %s | %s |\n'
                         % (r.get('hedef', ''), (r.get('durum') or '—').strip(),
                            (r.get('dCq_karsiligi') or '—').strip()))
        fh.write(u'\n## Bu sayının anlamadığı şey\n\n')
        fh.write(u'Bu tablo **eşik ve kapsam** hükmüdür. Bir çiftin burada olması '
                 u'geometri kapısından geçtiği, plaka içi jel ayrımının temiz olduğu '
                 u'ya da NCBI katmanının onayladığı anlamına gelmez. Onlar ayrı '
                 u'dosyalarda ölçülür: `primer_final/geometri_denetimi_*.tsv`, '
                 u'`TEK_TUS_SONUC/DENETIM_RAPORU.md`, '
                 u'`DOGRULAMA_SONUC/NCBI_KATMAN4_RAPORU.md`.\n')
    print('yazildi: %s' % y)
    print('  siparis edilecek: %d cift (%d hedef ozgul + %d evrensel), %d oligo'
          % (len(kesin) + len(evr), len(kesin), len(evr), 2 * (len(kesin) + len(evr))))
    print('  siparis disi: %d' % len(disi))
    return 0


if __name__ == '__main__':
    sys.exit(main())
