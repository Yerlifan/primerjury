# -*- coding: utf-8 -*-
"""BIR HEDEFIN UYELIGINI OLCULEN KIMLIKTEN TANIMLA.

NEDEN (2026-08-11)
------------------
Petriella_cinsi cifti siparis listesinde SINIF=KESIN, sart=KOSULSUZ yazili ve
dCq'su 3,46. Ama uc sey eksikti:
  * uyelik tablosunda satiri YOK -> "G asamasinda uyelik tanimi YOK"
  * panel kaynaginda satiri YOK  -> plakasi ve Ta'si belirsiz
  * dolayisiyla dislama kapsamasi da denetlenemiyor
"Kosulsuz" sozcugu raporda "tartisilacak bir sey yok" diye okunur; oysa vardi.

NE YAPAR
--------
Uyeligi ELLE YAZMAZ: hedefin sinifindaki (F1/F2/A1/A2/B) butun kutulari gezer,
OLCULEN kimliginde verilen kalibi tasiyanlari UYE, kalanlari RAKIP yapar ve
gerekcesini satir satir yazar. Sonra:
  1) uyelik_yeniden_turetme_uyelik_*.tsv dosyasina satir ekler (P ve K bunu okur)
  2) screening/hedef_uyelik.tsv dosyasina acik tanimi ekler
  3) istenirse panel kaynagina plaka/Ta satiri ekler
Her dosyanin once yedegini alir ve ne yaptigini basar.

Kosum (once PLAN):
  python verification/uyelik_tanimla.py --kok . --hedef Petriella_cinsi --kalip Petriella
  ... --yaz            uyelik dosyalarina yaz
  ... --panel-satiri P1:55   panel kaynagina plaka/Ta satirini da ekle
"""
from __future__ import print_function

import argparse
import glob
import io
import os
import re
import shutil
import sys
import time


def kimlikler(kok):
    y = os.path.join(kok, 'TUM_KIMLIK_SONUC', 'tum_kutu_kimlikleri.tsv')
    sat = [x.rstrip('\n').split('\t') for x in io.open(y, encoding='utf-8')]
    bi = None
    for i, r in enumerate(sat):
        if r and r[0].strip() == 'kutu':
            bi = i
            break
    if bi is None:
        return []
    H = sat[bi]
    return [dict(zip(H, r)) for r in sat[bi + 1:]
            if len(r) >= len(H) - 1 and r[0].strip()]


def uyelik_dosyasi(kok):
    a = glob.glob(os.path.join(kok, 'uyelik_yeniden_turetme_uyelik_*.tsv'))
    a += glob.glob(os.path.join(kok, 'engine_SONUC', '*uyelik*.tsv'))
    if not a:
        return None
    a.sort(key=lambda p: (os.path.getmtime(p), os.path.basename(p)))
    return a[-1]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--kok', default='.')
    p.add_argument('--hedef', required=True)
    p.add_argument('--kalip', required=True,
                   help='olculen kimlikte aranacak metin, orn: Petriella')
    p.add_argument('--sinif', default='',
                   help='kutu sinifi (F2 gibi); bos ise siparis listesinden okunur')
    p.add_argument('--yaz', action='store_true')
    p.add_argument('--panel-satiri', default='',
                   help='panel kaynagina satir ekle: PLAKA:Ta  (orn P1:55)')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)

    # sinif
    sinif = a.sinif.strip().upper()
    if not sinif:
        import csv
        y = os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'panel_tek_protokol.tsv')
        if os.path.exists(y):
            with io.open(y, encoding='utf-8') as fh:
                for r in csv.DictReader((l for l in fh if not l.startswith('#')),
                                        delimiter='\t'):
                    if (r.get('hedef') or '').strip() == a.hedef:
                        sinif = (r.get('sinif') or '').strip().upper()
                        break
    if not sinif:
        sys.exit('HATA: sinif bulunamadi, --sinif ile verin.')

    K = kimlikler(kok)
    ayni_sinif = [r for r in K if r['kutu'].split('-')[0].upper() == sinif]
    if not ayni_sinif:
        sys.exit('HATA: %s sinifinda kutu bulunamadi.' % sinif)

    uye, rakip = [], []
    for r in ayni_sinif:
        metin = ' '.join((r.get('DOGRULANAN_KIMLIK') or '',
                          r.get('ONERILEN_AD') or '',
                          r.get('en_iyi_isabet') or ''))
        (uye if a.kalip.lower() in metin.lower() else rakip).append(r)

    print('=' * 78)
    print('  UYELIK TANIMI  hedef=%s  sinif=%s  kalip="%s"'
          % (a.hedef, sinif, a.kalip))
    print('=' * 78)
    print('  %s sinifinda kutu: %d' % (sinif, len(ayni_sinif)))
    print()
    print('  UYE (%d) - olculen kimliginde "%s" geciyor:' % (len(uye), a.kalip))
    for r in uye:
        print('    %-18s %-30s %s' % (r['kutu'],
                                      (r.get('ONERILEN_AD') or '')[:30],
                                      (r.get('SAVUNULABILIR_DUZEY') or '')[:18]))
    print()
    print('  RAKIP (%d):' % len(rakip))
    for r in rakip:
        print('    %-18s %s' % (r['kutu'], (r.get('ONERILEN_AD') or '')[:52]))
    if not uye:
        sys.exit('\nHATA: hicbir kutu kalibi tasimiyor - uyelik tanimlanamaz.')

    if not a.yaz:
        print()
        print('  Bu bir PLANDIR. Yazmak icin --yaz ekleyin.')
        return 0

    # ---- 1) uyelik_yeniden_turetme dosyasi ----
    uy = uyelik_dosyasi(kok)
    if not uy:
        sys.exit('HATA: uyelik_yeniden_turetme_uyelik_*.tsv bulunamadi.')
    sat = [l.rstrip('\n').split('\t') for l in io.open(uy, encoding='utf-8')]
    bas = sat[0]
    varsa = [i for i, r in enumerate(sat[1:], 1) if r and r[0].strip() == a.hedef]
    shutil.copy2(uy, uy + '.yedek_%s_uyelik' % time.strftime('%H%M'))
    yeni = [''] * len(bas)
    def koy(ad, deger):
        if ad in bas:
            yeni[bas.index(ad)] = deger
    koy(bas[0], a.hedef)
    koy('sinif', sinif)
    koy('eski_uye_kutular', '')
    koy('yeni_uye_kutular', ';'.join(r['kutu'] for r in uye))
    koy('rakip_kutular', ';'.join(r['kutu'] for r in rakip))
    koy('kanit', u'olculen kimlikte "%s" gecen %d kutu UYE, kalan %d kutu RAKIP '
                 u'(%s). Kraken etiketine BAKILMADI.'
                 % (a.kalip, len(uye), len(rakip), time.strftime('%Y-%m-%d')))
    if varsa:
        sat[varsa[0]] = yeni
        print('\n  uyelik satiri GUNCELLENDI: %s' % os.path.basename(uy))
    else:
        sat.append(yeni)
        print('\n  uyelik satiri EKLENDI: %s' % os.path.basename(uy))
    with io.open(uy, 'w', encoding='utf-8', newline='') as fh:
        for r in sat:
            fh.write(u'\t'.join(r) + u'\n')

    # ---- 2) hedef_uyelik.tsv (acik tanim) ----
    hy = os.path.join(kok, 'screening', 'hedef_uyelik.tsv')
    if os.path.exists(hy):
        shutil.copy2(hy, hy + '.yedek_%s_uyelik' % time.strftime('%H%M'))
        metin = io.open(hy, encoding='utf-8').read()
        if not re.search(r'(?m)^%s\t' % re.escape(a.hedef), metin):
            tx = sorted({r['kutu'].split('_')[-1] for r in uye})
            if not metin.endswith('\n'):
                metin += '\n'
            metin += (u'%s\t%s\t\tOLCULEN KIMLIK\tolculen kimlikte "%s" gecen '
                      u'%s sinifi kutular (%s). Taxidler Kraken etiketidir, '
                      u'uyelik OLCULEN kimlige gore secildi.\n'
                      % (a.hedef, ','.join(tx), a.kalip, sinif,
                         time.strftime('%Y-%m-%d')))
            io.open(hy, 'w', encoding='utf-8', newline='').write(metin)
            print('  acik tanim EKLENDI: screening/hedef_uyelik.tsv')
        else:
            print('  acik tanim zaten var: screening/hedef_uyelik.tsv')

    # ---- 3) panel kaynagina satir ----
    if a.panel_satiri:
        try:
            plaka, ta = a.panel_satiri.split(':')
        except ValueError:
            sys.exit('HATA: --panel-satiri PLAKA:Ta biciminde olmali (orn P1:55)')
        py = os.path.join(kok, 'primer_final',
                          'devir_ciftleri_20260802_sonrotus_TESLIM.tsv')
        psat = [l.rstrip('\n').split('\t') for l in io.open(py, encoding='utf-8')]
        pb = psat[0]
        if any(len(r) > pb.index('Hedef') and r[pb.index('Hedef')].strip() == a.hedef
               for r in psat[1:]):
            print('  panel satiri zaten var, dokunulmadi.')
        else:
            import csv as _csv
            sl = None
            with io.open(os.path.join(kok, 'TEK_PROTOKOL_SONUC',
                                      'SIPARIS_LISTESI.tsv'), encoding='utf-8') as fh:
                for r in _csv.DictReader((l for l in fh if not l.startswith('#')),
                                         delimiter='\t'):
                    if (r.get('hedef') or '').strip() == a.hedef:
                        sl = r
                        break
            if not sl:
                sys.exit('HATA: %s siparis listesinde yok, panel satiri yazilamaz.'
                         % a.hedef)
            shutil.copy2(py, py + '.yedek_%s_panelsatiri' % time.strftime('%H%M'))
            yeni_p = [''] * len(pb)
            def p_koy(ad, deger):
                if ad in pb:
                    yeni_p[pb.index(ad)] = deger
            p_koy('Plaka', plaka.strip())
            p_koy('Ta (C)', ta.strip())
            p_koy('Hedef', a.hedef)
            p_koy('Duzey', 'cins')
            p_koy('Amplikon sinifi', sinif)
            for k in pb:
                if k.startswith('Ileri primer'):
                    yeni_p[pb.index(k)] = (sl.get('F') or '').strip().upper()
                if k.startswith('Geri primer'):
                    yeni_p[pb.index(k)] = (sl.get('R') or '').strip().upper()
            p_koy('Urun (bp)', (sl.get('urun_bp') or '').strip())
            # Tm/GC/dTm BOS birakilir: bir sonraki geometri kosusu OLCUP yazar.
            # Buraya tahmin yazmak, tam da duzelttigimiz hatayi geri getirirdi.
            psat.append(yeni_p)
            with io.open(py, 'w', encoding='utf-8', newline='') as fh:
                for r in psat:
                    fh.write(u'\t'.join(r) + u'\n')
            print('  panel satiri EKLENDI: plaka %s, Ta %s (Tm/GC bos - geometri '
                  'kosusu olcup yazacak)' % (plaka, ta))
    print('=' * 78)
    return 0


if __name__ == '__main__':
    sys.exit(main())
