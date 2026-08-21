# -*- coding: utf-8 -*-
"""SIPARIS FORMU  -  oligo tedarikcisine gidecek TEK dogru liste.

NEDEN VAR (2026-08-10 gece, denetimde yakalandi)
------------------------------------------------
Sabah ozeti "diziler buradan kopyalanacak" diye
PrimerJury_PCR_Paneli_2026-08-02_TESLIM.xlsx dosyasinin "2 Panel" sayfasini
gosteriyordu. O sayfadaki ALTI ciftin dizisi ESKI:

    Bacteroidales_kumesi, Bakteri_universal, Mantar_universal (F1),
    Microascaceae_askomikot, Petriella_musispora, Petrimonas_cinsi

Yani o dosyadan siparis verilseydi yirmi ciftin altisi YANLIS oligo olarak
gelirdi. Ustelik geometri denetimi dosyasi da baska bir altili tasiyordu.
Uc ayri dosya, uc ayri "dogru" - hangisinin dogru oldugunu soyleyen bir sey
yoktu.

Bu betik o boslugu kapatir: siparis listesi PANELIN KENDI KAYNAGINDAN
uretilir, elle yazilmaz, ve her uretimde kaynak dosyanin ozeti (md5) yazilir.
Kaynak degisirse form da degisir; degismezse ayni cikar.

Kosum:
    python verification/siparis_formu.py --kok .
Cikti:
    SIPARIS_FORMU.tsv   (tedarikciye yapistirilacak)
    SIPARIS_FORMU.md    (insanin okuyacagi hali)
"""
from __future__ import print_function

import argparse
import csv
import hashlib
import io
import os
import re
import sys
import time

PANEL = os.path.join('primer_final', 'devir_ciftleri_20260802_sonrotus_TESLIM.tsv')


def _tsv(yol):
    if not os.path.exists(yol):
        return []
    with io.open(yol, encoding='utf-8') as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith('#')),
                                   delimiter='\t'))


def panel_oku(kok):
    y = os.path.join(kok, PANEL)
    sat = [l.rstrip('\n').split('\t') for l in io.open(y, encoding='utf-8')]
    b = sat[0]
    iP, iT, iH = b.index('Plaka'), b.index('Ta (C)'), b.index('Hedef')
    iU = b.index('Urun (bp)')
    iF = next(i for i, x in enumerate(b) if x.startswith('Ileri primer'))
    iR = next(i for i, x in enumerate(b) if x.startswith('Geri primer'))
    out = []
    for r in sat[1:]:
        if len(r) <= max(iU, iF, iR) or not r[iH].strip():
            continue
        F, R = r[iF].strip().upper(), r[iR].strip().upper()
        if not re.fullmatch(r'[ACGT]+', F or '') or not re.fullmatch(r'[ACGT]+', R or ''):
            continue
        out.append(dict(hedef=r[iH].strip(), plaka=r[iP].strip(), ta=r[iT].strip(),
                        F=F, R=R, urun=r[iU].strip()))
    return out, y


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--kok', default='.')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)

    pnl, panel_yolu = panel_oku(kok)
    ozet = hashlib.md5(io.open(panel_yolu, 'rb').read()).hexdigest()[:12]
    PN = {c['hedef']: c for c in pnl}

    # HANGI CIFT SIPARIS EDILIR sorusunun tek yetkilisi SIPARIS_LISTESI'dir;
    # dizinin kendisi de oradan alinir. Panel kaynagi ile KARSILASTIRILIR ama
    # yerine gecmez: 2026-08-10'da Petriella_cinsi siparis listesinde KESIN
    # yaziliyken panel kaynaginda SATIRI YOKTU - panel kaynagindan uretilen
    # bir form o cifti sessizce dusururdu.
    sl_yol = os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'SIPARIS_LISTESI.tsv')
    sl = _tsv(sl_yol)
    if not sl:
        sys.exit('HATA: %s okunamadi.' % sl_yol)

    girer, girmez, uyari = [], [], []
    for r in sl:
        ad = (r.get('hedef') or '').strip()
        F = (r.get('F') or '').strip().upper()
        R = (r.get('R') or '').strip().upper()
        if not ad:
            continue
        if not re.fullmatch(r'[ACGT]+', F or '') or not re.fullmatch(r'[ACGT]+', R or ''):
            uyari.append(u'%s: siparis listesinde dizi yok ya da A/C/G/T disi '
                         u'karakter var - forma ALINMADI' % ad)
            continue
        p = PN.get(ad)
        if p is None:
            uyari.append(u'%s: siparis listesinde VAR ama panel kaynaginda SATIRI '
                         u'YOK. Plaka ve Ta bilgisi bu cift icin bilinmiyor.' % ad)
        elif (p['F'], p['R']) != (F, R):
            uyari.append(u'%s: siparis listesi ile panel kaynagi FARKLI dizi '
                         u'soyluyor. Form siparis listesini kullandi.' % ad)
        sinif = (r.get('SINIF') or '').strip().upper()
        c = dict(hedef=ad, F=F, R=R, sinif=sinif or '?',
                 plaka=(p or {}).get('plaka', '?'), ta=(p or {}).get('ta', '?'),
                 urun=(r.get('urun_bp') or (p or {}).get('urun', '')).strip(),
                 sart=(r.get('siparis_sarti') or '').strip(),
                 adF=(r.get('oligo_adi_F') or (ad + '_F')).strip(),
                 adR=(r.get('oligo_adi_R') or (ad + '_R')).strip())
        (girer if sinif in ('KESIN', 'EVRENSEL') else girmez).append(c)

    ty = os.path.join(kok, 'SIPARIS_FORMU.tsv')
    with io.open(ty, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# SIPARIS FORMU - %s\n' % time.strftime('%Y-%m-%d %H:%M'))
        fh.write(u'# Kaynak: %s (md5 %s)\n' % (PANEL.replace('\\', '/'), ozet))
        fh.write(u'# Bu dosya ELLE YAZILMAZ, panelin kaynagindan uretilir.\n')
        fh.write(u'# xlsx dosyasindan KOPYALAMAYIN: 2026-08-10 itibariyla oradaki\n'
                 u'# alti ciftin dizisi ESKIDIR (Bacteroidales, Bakteri_universal,\n'
                 u'# Mantar F1, Microascaceae, Petriella_musispora, Petrimonas).\n')
        for u_ in uyari:
            fh.write(u'# UYARI: %s\n' % u_)
        fh.write(u'oligo_adi\tdizi_5_3\tuzunluk\thedef\tyon\tplaka\tTa_C\turun_bp\tsinif\n')
        for c in girer:
            for ad, d, yon in ((c['adF'], c['F'], 'ileri'), (c['adR'], c['R'], 'geri')):
                fh.write(u'%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n'
                         % (ad, d, len(d), c['hedef'], yon, c['plaka'], c['ta'],
                            c['urun'], c['sinif']))

    my = os.path.join(kok, 'SIPARIS_FORMU.md')
    with io.open(my, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# Sipariş formu\n\n')
        fh.write(u'Üretim: %s · kaynak `%s` (md5 `%s`)\n\n'
                 % (time.strftime('%Y-%m-%d %H:%M'), PANEL.replace('\\', '/'), ozet))
        fh.write(u'> **Excel dosyasından kopyalamayın.** '
                 u'`PrimerJury_PCR_Paneli_2026-08-02_TESLIM.xlsx` içindeki '
                 u'`2 Panel` sayfasında altı çiftin dizisi eskidir: Bacteroidales, '
                 u'Bakteri evrensel, Mantar F1, Microascaceae, Petriella musispora, '
                 u'Petrimonas. O sayfadan sipariş verilirse yirmi çiftin altısı '
                 u'yanlış oligo olarak gelir. Doğru diziler aşağıdadır ve panelin '
                 u'kendi kaynağından üretilmiştir.\n\n')
        fh.write(u'**%d çift = %d oligo**\n\n' % (len(girer), 2 * len(girer)))
        if uyari:
            fh.write(u'### Uyarılar\n\n')
            for u_ in uyari:
                fh.write(u'- %s\n' % u_)
            fh.write(u'\n')
        fh.write(u'| oligo adı | dizi (5→3) | uz | hedef | yön | plaka | Ta | ürün |\n'
                 u'|---|---|---|---|---|---|---|---|\n')
        for c in girer:
            for ad, d, yon in ((c['adF'], c['F'], 'ileri'), (c['adR'], c['R'], 'geri')):
                fh.write(u'| %s | `%s` | %d | %s | %s | %s | %s | %s |\n'
                         % (ad, d, len(d), c['hedef'], yon, c['plaka'], c['ta'], c['urun']))
        if girmez:
            fh.write(u'\n## Sipariş dışı (silinmedi, bilgi için)\n\n')
            fh.write(u'| hedef | sınıf | ürün |\n|---|---|---|\n')
            for c in girmez:
                fh.write(u'| %s | %s | %s |\n' % (c['hedef'], c['sinif'], c['urun']))
        fh.write(u'\n## Bu formun söylemediği şey\n\n')
        fh.write(u'Bu form **hangi dizinin sipariş edileceğini** söyler. Bu '
                 u'çiftlerin geometri kapısından geçtiğini, plaka içi jel ayrımının '
                 u'temiz olduğunu ya da eşik kuralının hangisi olduğunu söylemez. '
                 u'Onlar `TEK_TUS_SONUC/DENETIM_RAPORU.md` ve `GECE_BULGULARI.md` '
                 u'içindedir ve sipariş vermeden önce okunmalıdır.\n')

    print('yazildi: %s' % ty)
    print('yazildi: %s' % my)
    print('  siparise giren: %d cift = %d oligo' % (len(girer), 2 * len(girer)))
    print('  siparis disi  : %d cift' % len(girmez))
    print('  kaynak md5    : %s' % ozet)
    return 0


if __name__ == '__main__':
    sys.exit(main())
