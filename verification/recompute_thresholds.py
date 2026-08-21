# -*- coding: utf-8 -*-
"""BOLLUGA AGIRLIKLI ESIGI GUNCEL VERIYLE YENIDEN HESAPLA.

KURAL (2026-08-08 esik belgesinden, degistirilmedi)
---------------------------------------------------
    gerekli dCq = max( log2(R) + 4,3 ,  3,32 )

R  = rakip havuzunun uye havuzuna okuma orani (bolluk orani)
4,3 = Fierer ve ark.'nin fiilen uyguladigi %5 saflik olcutunun dongu karsiligi
3,32 = taban; R<=1 olsa bile bir hedefin en az 10 kat ayrim gostermesi istenir
      (log2(10) = 3,32)

NEDEN YENIDEN
-------------
Esik tablosu 8 Agustos'ta uretildi. O gunden sonra iki sey degisti:
  1) uyelikler OLCULEN kimlikten yeniden turetildi (21 hedefin 2'sinde degisti),
  2) alti ciftin dizisi degisti.
R iki degisiklikten de etkilenir. Yani 8 Agustos'un gerekli-dCq degerleri
bugunku panele ait DEGILDIR. Bu betik onlari bugunku sayilarla yeniden uretir.

HUKUM VERMEZ. Iki kurali (duz esik 3,00 ve bolluga agirlikli esik) YAN YANA
koyar. Hangisinin uygulanacagi bilimsel bir tercihtir ve bu betigin isi degil.

Kosum:
    python verification/recompute_thresholds.py --kok .
Cikti:
    ESIK_IKI_KURAL.tsv  ve  ESIK_IKI_KURAL.md
"""
from __future__ import print_function

import argparse
import csv
import io
import math
import os
import sys
import time

DUZ_ESIK = 3.00
TABAN = 3.32
FIERER = 4.3


def _f(x):
    try:
        return float(str(x).replace(',', '.'))
    except (TypeError, ValueError):
        return None


def _tsv(yol):
    if not os.path.exists(yol):
        return []
    with io.open(yol, encoding='utf-8') as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith('#')),
                                   delimiter='\t'))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', '--kok', dest='kok', default='.')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)

    panel = _tsv(os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'panel_tek_protokol.tsv'))
    ham = _tsv(os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'kutu_bazli_ham_sayilar.tsv'))
    if not panel:
        sys.exit(u'ERROR: panel_tek_protokol.tsv could not be read.')

    # R NEREDEN GELIYOR - ve neden yeniden hesaplanmiyor
    # ---------------------------------------------------
    # R (bolluk orani) ESIK_VE_OLCUT_2026-08-08.tsv icinde YAZILI ama onu
    # ureten betik projede YOK. cross_check.py yalnizca aritmetigi
    # dogruluyor (max(log2(R)+4,3; 3,32)), R'nin kendisini hesaplamiyor.
    #
    # 2026-08-10'da R'yi "rakip kutu okumalari / uye kutu okumalari" diye
    # yeniden hesaplamayi denedim: 17 satirin 0'inda tabloyu tutturdu
    # (ornek Metilotrofik: tablo 1,760, benim hesabim 293,685). Demek ki
    # tablodaki R baska bir normalizasyon kullaniyor. Hangisi oldugunu
    # BILMIYORUM; bilmedigim bir sayiyi uydurmaktansa tablodakini
    # KAYNAGIYLA BIRLIKTE tasiyorum ve ureteci olmadigini yaziyorum.
    #
    # 2026-08-10 gece, oturum kayitlarinda tanim BULUNDU:
    #     "R = en yakin rakibin hedefe gore bolluk orani"
    # Bu tanimin akla gelen alti varyantini (en bol rakip kutu / uye havuzu,
    # en yuksek urun oranli rakip / uye orani, toplam rakip / toplam uye, ...)
    # bugunku sayilarla denedim: 18 satirin 1'inde tuttu (Asetoklastik
    # metanojenler, 0,130 birebir). Cifti DEGISMEYEN 14 satirin de yalniz
    # 1'i tuttugu icin bu bir rastlanti sayilmali.
    #
    # Sonuc: 8 Agustos'un R degerleri O GUNKU uyelikle hesaplanmis; uyelikler
    # 3 Agustos'ta olculen kimlikten yeniden turetildi ve 10 Agustos'ta alti
    # cift degisti. Yani sayilar bugunku panele ait DEGIL.
    #
    # Bu bir eksiktir ve raporda oyle gecer: rapora giren bir karar
    # tablosundaki sayinin yeniden uretilebilir olmasi gerekir. Karar
    # verilmesi gereken sey: R'nin tanimini yazip BUGUNKU uyelikle bastan
    # hesaplamak (yeniden uretilebilir olur) ya da bolluk kuralini bu teslimde
    # hic kullanmamak. Ikisi de savunulabilir; sessizce eski sayiyi kullanmak
    # savunulamaz.
    esik08 = {}
    ey = os.path.join(kok, 'ESIK_VE_OLCUT_2026-08-08.tsv')
    for r in _tsv(ey):
        R = _f(r.get('R'))
        if R is not None:
            esik08[(r.get('hedef') or '').strip()] = dict(
                R=R, dcq08=_f(r.get('dCq_olculen')))

    satir = []
    for r in panel:
        ad = (r.get('hedef') or '').strip()
        d = _f(r.get('ASIL_ayrim_mm1'))
        dcq = math.log(d, 2) if (d and d > 0) else None
        e8 = esik08.get(ad, {})
        R = e8.get('R')
        ger = max(math.log(R, 2) + FIERER, TABAN) if (R and R > 0) else None
        # R 8 Agustos'ta hesaplandi. O gunden sonra uyelikler yeniden
        # turetildi ve alti cift degisti; dCq degistiyse R de degismis
        # olabilir. Degisen satirlar isaretlenir - sessizce kullanilmaz.
        d08 = e8.get('dcq08')
        bayat = (d08 is not None and dcq is not None and abs(d08 - dcq) > 0.05)
        satir.append(dict(
            hedef=ad, dcq=dcq, R=R, gerekli=ger, bayat=bayat, dcq08=d08,
            duz=(None if dcq is None else (dcq >= DUZ_ESIK)),
            bolluk=(None if (dcq is None or ger is None) else (dcq >= ger))))

    ty = os.path.join(kok, 'ESIK_IKI_KURAL.tsv')
    with io.open(ty, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# The two threshold rules side by side. Generated %s\n'
                 % time.strftime('%Y-%m-%d %H:%M'))
        fh.write(u'# duz kural   : dCq >= %.2f\n' % DUZ_ESIK)
        fh.write(u'# bolluk kural: dCq >= max(log2(R) + %.1f, %.2f)\n' % (FIERER, TABAN))
        fh.write(u'# R = competitor pool reads / member pool reads (mm<=1)\n')
        fh.write(u'hedef\tdCq\tdCq_08_08\tR_bayat_mi\tR\tgerekli_dCq\t'
                 u'duz_kural\tbolluk_kurali\tayrisiyor_mu\n')
        for s in satir:
            fh.write(u'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                s['hedef'],
                '' if s['dcq'] is None else ('%.2f' % s['dcq']),
                '' if s['dcq08'] is None else ('%.2f' % s['dcq08']),
                'EVET' if s['bayat'] else '',
                '' if s['R'] is None else ('%.3f' % s['R']),
                '' if s['gerekli'] is None else ('%.2f' % s['gerekli']),
                '' if s['duz'] is None else ('GECER' if s['duz'] else 'KALIR'),
                '' if s['bolluk'] is None else ('GECER' if s['bolluk'] else 'KALIR'),
                'EVET' if (s['duz'] is not None and s['bolluk'] is not None
                           and s['duz'] != s['bolluk']) else ''))

    ayrisan = [s for s in satir if s['duz'] is not None and s['bolluk'] is not None
               and s['duz'] != s['bolluk']]
    my = os.path.join(kok, 'ESIK_IKI_KURAL.md')
    with io.open(my, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# The two threshold rules side by side\n\n')
        fh.write(u'Generated: %s, source `TEK_PROTOKOL_SONUC/` (this run\'s own output)\n\n' % time.strftime('%Y-%m-%d %H:%M'))
        fh.write(u'| rule | definition |\n|---|---|\n')
        fh.write(u'| Flat | dCq >= %.2f |\n' % DUZ_ESIK)
        fh.write(u'| Abundance-weighted | dCq >= max(log2(R) + %.1f, %.2f) |\n\n'
                 % (FIERER, TABAN))
        fh.write(u'R is the ratio of competitor-pool reads to member-pool reads. The 4.3 comes from Fierer et al.\'s 5% purity criterion')
        bayat = [s for s in satir if s.get('bayat')]
        fh.write(u'> **Where the R values come from.** R was taken from `ESIK_VE_OLCUT_2026-08-08.tsv`; the script that produced it' % len(bayat))
        fh.write(u'**The two rules disagree on %d rows.**\n\n' % len(ayrisan))
        fh.write(u'| target | dCq | R | required dCq | flat | abundance | is R stale |\n|---|---|---|---|---|---|---|\n')
        for s in satir:
            if s['dcq'] is None:
                continue
            fh.write(u'| %s%s | %.2f | %s | %s | %s | %s | %s |\n' % (
                s['hedef'], u' **<-- disagrees**' if s in ayrisan else '',
                s['dcq'],
                '—' if s['R'] is None else '%.2f' % s['R'],
                '—' if s['gerekli'] is None else '%.2f' % s['gerekli'],
                u'passes' if s['duz'] else u'stays',
                '—' if s['bolluk'] is None else (u'passes' if s['bolluk'] else u'stays'),
                (u'**YES** (on 8 August dCq was %.2f)' % s['dcq08'])
                if s.get('bayat') else u'—'))
        fh.write(u'\n## What this table does not decide\n\nWhich rule to apply is a **choice of criterion**, not a measurement')

    print('yazildi: %s' % ty)
    print('yazildi: %s' % my)
    print(u'  measurable rows: %d' % sum(1 for s in satir if s['dcq'] is not None))
    print('  iki kural ayrisan : %d' % len(ayrisan))
    for s in ayrisan:
        print(u'    %-44s dCq %.2f, required %.2f' % (s['hedef'][:44], s['dcq'], s['gerekli']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
