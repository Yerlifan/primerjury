# -*- coding: utf-8 -*-
"""GEOMETRI KAPISINI SU ANKI PANEL DIZILERIYLE YENIDEN KOS.

NEDEN (2026-08-10, denetimde yakalandi)
---------------------------------------
Paneldeki ALTI ciftin dizisi 2 Agustos'taki geometri denetiminden bu yana
degisti; geometri dosyasi hala ESKI dizileri tasiyor. Yani bu alti cift
panelin kendi geometri kurallarindan (uzunluk 18-25, GC %40-60, Tm 58-62,
sac tokasi, dimer, 3' uc) HIC GECMEDI.

Dahasi, teslim tablosunun Tm sutunu bu bes yeni cift icin 48,9-51,1 yaziyor.
Bu deger panelin motoruyla (primer3, mv=50 dv=1,5 dntp=0,6 dna=50) uyusmuyor:
bagimsiz en yakin komsu hesabiyla ayni bes primer 55-58 cikiyor ve saglam 16
ciftte tablo ile bagimsiz hesap arasindaki fark sabit +3,72 ± 0,22 C iken bu
beste -6,38 C. Iki ayri kume; tek bir yontemle olusamaz.

Neden onemli: teslim tablosunda "Ta = min(Tm) - 3" yaziyor. Yazili Tm'lerle
P1 icin Ta 47 cikar; oysa plakada 55 yaziyor. Tabloyu okuyan biri Ta'yi
yeniden hesaplarsa panelin tamami yanlis sicaklikta kosar. Tutarli Tm ile
ayni kural bes plaka grubunun BESINDE de yazili Ta'yi 0,6 C icinde veriyor -
yani TASARIM DOGRU, YAZILAN SAYI YANLIS.

Bu betik hukum vermez, OLCER: panelin kendi motoruyla butun primerleri
yeniden hesaplar, kural ihlallerini listeler ve tabloda yazan degerle
karsilastirir.

Kosum:
    python verification/refresh_geometry.py --kok .
    python verification/refresh_geometry.py --kok . --yaz     (panel tablosunu da guncelle)
"""
from __future__ import print_function

import argparse
import io
import os
import re
import sys
import time

PANEL = os.path.join('primer_final', 'devir_ciftleri_20260802_sonrotus_TESLIM.tsv')


def geo_yukle(kok):
    """Panelin KENDI geometri modulunu getirir. Ayri bir kopya YAZILMAZ -
    iki kopya zamanla ayrisir ve hangisinin dogru oldugu bilinmez."""
    for aday in ('engine', 'engine'):
        d = os.path.join(kok, aday)
        if os.path.exists(os.path.join(d, 'geometry_core.py')):
            sys.path.insert(0, d)
            try:
                import geometry_core
            except ImportError as e:
                # primer3 yoksa YIGIN IZI basip kullaniciyi korkutma; ne
                # eksik oldugunu ve nasil kurulacagini soyle.
                return None, (u'geometry_core.py bulundu (%s) ama yuklenemedi: %s\n'
                              u'      Kurulum: pip3 install primer3-py '
                              u'--break-system-packages'
                              % (os.path.join(aday, 'geometry_core.py'), e))
            return geo, os.path.join(aday, 'geometry_core.py')
    return None, u'geometry_core.py bulunamadi (engine ya da engine)'


# --- URUN BOYU: yazilan degil OLCULEN --------------------------------------
# Tm gibi urun boyu da elle yazilabilen bir sayidir ve elle yazilan her sayi
# bir gun kayar. 2026-08-10: Proteolitik_Synergistaceae satirinda 173 bp
# yaziyordu, olcum 172 diyor. Bir baz, ama jel ayrimi ve bant sinifi
# hesaplarina giriyor. Olcut panelin kendi olcutu: mm<=1, 3' son iki baz tam.
def _rc(x):
    return x.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def _yer(p, d):
    n = len(p)
    out = []
    for i in range(len(d) - n + 1):
        f = 0
        for a, c in zip(p, d[i:i + n]):
            if a != c:
                f += 1
                if f > 1:
                    break
        else:
            if p[-2:] == d[i + n - 2:i + n]:
                out.append(i)
    return out


def konsensus_yukle(kok):
    """Kanonik konsensusler - GLOB DEGIL, INDEKS.tsv.

    konsensus_kanonik klasorunde 250 dosya var, yalnizca 100'u gecerli.
    Kalan 150'si bagli klasorde silinemeyen kalinti; 33 kutuda ayni kutunun
    FARKLI icerikli iki-uc surumu duruyor. Glob ile okumak hangi surumun
    kullanilacagini dosya adi sirasina birakir - panelin gormedigi bir diziyle
    olcum yapmis oluruz. Panelin kendi yukleyicisi de bu yuzden indeks okur.
    """
    import csv as _csv
    d = os.path.join(kok, 'konsensus_kanonik')
    ix = os.path.join(d, 'INDEKS.tsv')
    out = {}
    if not os.path.exists(ix):
        return out
    with io.open(ix, encoding='utf-8') as fh:
        for r in _csv.DictReader(fh, delimiter='\t'):
            f = os.path.join(d, (r.get('dosya') or '').strip())
            if os.path.exists(f):
                out[os.path.basename(f)] = ''.join(
                    l.strip() for l in io.open(f, encoding='utf-8', errors='replace')
                    if not l.startswith('>')).upper()
    return out


def urun_boylari(F, R, kons):
    Rrc = _rc(R)
    boylar = set()
    for d in kons.values():
        fs = _yer(F, d)
        if not fs:
            continue
        rs = _yer(Rrc, d)
        for i in fs:
            for j in rs:
                if j >= i:
                    L = j + len(R) - i
                    if 40 <= L <= 600:
                        boylar.add(L)
    return sorted(boylar)


def panel_oku(kok):
    y = os.path.join(kok, PANEL)
    sat = [l.rstrip('\n').split('\t') for l in io.open(y, encoding='utf-8')]
    bas = sat[0]
    iP, iT, iH = bas.index('Plaka'), bas.index('Ta (C)'), bas.index('Hedef')
    iF = next(i for i, b in enumerate(bas) if b.startswith('Ileri primer'))
    iR = next(i for i, b in enumerate(bas) if b.startswith('Geri primer'))
    iTf = bas.index('Ileri Tm')
    iTr = bas.index('Geri Tm')
    idT = bas.index('dTm') if 'dTm' in bas else None
    iU = bas.index('Urun (bp)') if 'Urun (bp)' in bas else None
    out = []
    for r in sat[1:]:
        if len(r) <= iTr or not r[iH].strip():
            continue
        if not re.match(r'^[A-Za-z]', r[iH].strip()):
            continue
        F, R = r[iF].strip().upper(), r[iR].strip().upper()
        if not re.fullmatch(r'[ACGT]+', F or '') or not re.fullmatch(r'[ACGT]+', R or ''):
            continue
        out.append(dict(satir=r, plaka=r[iP].strip(), ta=r[iT].strip(),
                        hedef=r[iH].strip(), F=F, R=R,
                        tmF=r[iTf].strip(), tmR=r[iTr].strip(),
                        urun=(r[iU].strip() if iU is not None and len(r) > iU else ''),
                        iTf=iTf, iTr=iTr, idT=idT, iU=iU))
    return out, bas, sat, y


def yalniz_urun_boyu(kok, yaz_mi):
    """primer3 yokken bile kosabilen bolum: urun boyu olcumu."""
    kons = konsensus_yukle(kok)
    if not kons:
        print(u'  konsensus_kanonik is empty or missing, so product length could not be measured either.')
        return 2
    ciftler, bas, sat, panel_yolu = panel_oku(kok)
    sapan = []
    print(u'  --- product length (mm<=1, exact match at the last two 3\' bases) ---')
    for c in ciftler:
        bl = urun_boylari(c['F'], c['R'], kons)
        try:
            yazili = int(re.sub(r'\D', '', c.get('urun') or '0'))
        except ValueError:
            yazili = 0
        if bl and yazili and yazili not in bl:
            sapan.append((c, yazili, bl))
            print(u'      %-44s tabloda %d bp, olculen %s'
                  % (c['hedef'][:44], yazili, bl[:5]))
    print(u'    %d pairs measured, %d deviations' % (len(ciftler), len(sapan)))
    if yaz_mi and sapan:
        yed = panel_yolu + '.yedek_%s_urun' % time.strftime('%H%M')
        io.open(yed, 'w', encoding='utf-8', newline='').write(
            io.open(panel_yolu, encoding='utf-8').read())
        n = 0
        for c, yazili, bl in sapan:
            if len(bl) == 1 and c.get('iU') is not None:
                r = c['satir']
                if len(r) > c['iU']:
                    print(u'  product length corrected: %s %s -> %d'
                          % (c['hedef'], r[c['iU']], bl[0]))
                    r[c['iU']] = str(bl[0])
                    n += 1
        if n:
            with io.open(panel_yolu, 'w', encoding='utf-8', newline='') as fh:
                for r in sat:
                    fh.write(u'\t'.join(r) + u'\n')
            print(u'  %d rows updated (backup: %s)' % (n, os.path.basename(yed)))
    return 1 if sapan else 0


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', '--kok', dest='kok', default='.')
    p.add_argument('--write', '--yaz', dest='yaz', action='store_true',
                   help='also update the Tm and dTm columns of the panel table '
                        '(once yedek alinir)')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)

    geo, geo_yol = geo_yukle(kok)
    # primer3 yoksa Tm olculemez ama URUN BOYU olculebilir: o olcum
    # konsensus dizileri uzerinde yapilir ve primer3'e bagli degildir.
    # Bir eksik yuzunden yapilabilecek olcumu de atlamak, denetimi
    # gereksizce kor birakir.
    if geo is None:
        print('=' * 78)
        print(u'  Tm OLCULEMEDI - %s' % geo_yol)
        print(u'  Without the panel\'s OWN engine no Tm is computed, and no other formula')
        print(u'  is INVENTED. Product length does not need primer3, so only that')
        print(u'  is measured below.')
        print('=' * 78)
        return yalniz_urun_boyu(kok, a.yaz)
    print('=' * 78)
    print(u'  GEOMETRY RE-MEASUREMENT   %s' % time.strftime('%Y-%m-%d %H:%M'))
    print('  motor: %s   ayarlar: %s' % (geo_yol, geo.KW))
    print('=' * 78)

    ciftler, bas, sat, panel_yolu = panel_oku(kok)
    print(u'  pairs to measure in the panel: %d' % len(ciftler))

    sonuc = []
    ihlal = []
    uyusmaz = []
    for c in ciftler:
        tF, tR = geo.tm(c['F']), geo.tm(c['R'])
        vF, vR = geo.viol(c['F']), geo.viol(c['R'])
        het = geo.het(c['F'], c['R'])
        c['yeni_tmF'], c['yeni_tmR'] = tF, tR
        c['yeni_dTm'] = round(abs(tF - tR), 2)
        sonuc.append((c, tF, tR, vF, vR, het))
        if vF:
            ihlal.append(u'%s Ileri: %s' % (c['hedef'], ', '.join(vF)))
        if vR:
            ihlal.append(u'%s Geri: %s' % (c['hedef'], ', '.join(vR)))
        for et, yazili, yeni in (('F', c['tmF'], tF), ('R', c['tmR'], tR)):
            try:
                f = abs(float(str(yazili).replace(',', '.')) - yeni)
            except ValueError:
                f = 99.0
            if f > 1.0:
                uyusmaz.append(u'%s %s: tabloda %s, motor %.2f (fark %.2f)'
                               % (c['hedef'], et, yazili, yeni, f))

    print()
    print('  %-44s %-15s %-15s %s' % ('hedef', 'tablo F/R', 'motor F/R', 'ihlal'))
    for c, tF, tR, vF, vR, _het in sonuc:
        print('  %-44s %6s/%-6s %6.2f/%-6.2f %s'
              % (c['hedef'][:44], c['tmF'], c['tmR'], tF, tR,
                 ('F:' + ','.join(vF) if vF else '') +
                 (' R:' + ','.join(vR) if vR else '')))

    # --- plaka basina Ta = min(Tm) - 3 sinavi -----------------------------
    print()
    print(u'  --- the Ta = min(Tm) - 3 rule (with the engine\'s values) ---')
    grup = {}
    for c, tF, tR, _a, _b, _h in sonuc:
        grup.setdefault((c['plaka'], c['ta']), []).extend([tF, tR])
    ta_sorun = []
    for k, v in sorted(grup.items()):
        bek = min(v) - 3
        try:
            yazili = float(k[1])
        except ValueError:
            continue
        fark = abs(yazili - bek)
        print('    plaka %-3s Ta %-4s | min Tm %.2f -> beklenen Ta %.2f | fark %.2f %s'
              % (k[0], k[1], min(v), bek, fark, '' if fark <= 1.2 else '<-- TUTMUYOR'))
        if fark > 1.2:
            ta_sorun.append(u'plaka %s Ta %s: motor min Tm %.2f, kural %.2f diyor'
                            % (k[0], k[1], min(v), bek))

    # --- urun boyu olcumu --------------------------------------------------
    kons = konsensus_yukle(kok)
    urun_sapan = []
    if kons:
        for c, _tF, _tR, _a, _b, _h in sonuc:
            bl = urun_boylari(c['F'], c['R'], kons)
            c['olculen_urun'] = bl
            try:
                yazili = int(re.sub(r'\D', '', c.get('urun') or '0'))
            except ValueError:
                yazili = 0
            if bl and yazili and yazili not in bl:
                urun_sapan.append(u'%s: tabloda %d bp, olculen %s'
                                  % (c['hedef'], yazili, bl[:5]))
                c['yeni_urun'] = bl[0] if len(bl) == 1 else None
            else:
                c['yeni_urun'] = None
        print()
        print(u'  --- product length (mm<=1, exact match at the last two 3\' bases) ---')
        print(u'    %d pairs measured, %d deviations' % (len(sonuc), len(urun_sapan)))
        for x in urun_sapan:
            print('      %s' % x)
    else:
        print()
        print(u'  product length NOT MEASURED: the konsensus_kanonik directory is empty or missing.')

    cy = os.path.join(kok, 'primer_final',
                      'geometri_denetimi_%s.tsv' % time.strftime('%Y%m%d'))
    with io.open(cy, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# The geometry audit, with the panel sequences AS THEY ARE NOW. Produced %s\n'
                 % time.strftime('%Y-%m-%d %H:%M'))
        fh.write(u'# Motor: %s  ayarlar: %s\n' % (geo_yol, geo.KW))
        fh.write(u'Hedef\tPrimer\tDizi\tUz\tGC%\tTm\tHairpin Tm\tHomodimer Tm\tIhlal\n')
        for c, tF, tR, vF, vR, het in sonuc:
            for et, d, t, v in (('Ileri', c['F'], tF, vF), ('Geri', c['R'], tR, vR)):
                fh.write(u'%s\t%s\t%s\t%d\t%.1f\t%.2f\t%.1f\t%.1f\t%s\n'
                         % (c['hedef'], et, d, len(d), geo.gc(d), t,
                            geo.hp(d), geo.hd(d), '; '.join(v)))
            fh.write(u'%s\tCIFT\tdTm=%.2f  het=%.1f\t\t\t\t\t\t\n'
                     % (c['hedef'], abs(tF - tR), het))
    print()
    print(u'  written: %s' % cy)

    print()
    if uyusmaz:
        print(u'  THE Tm IN THE TABLE AND THE ENGINE DISAGREE (%d):' % len(uyusmaz))
        for x in uyusmaz:
            print('    * %s' % x)
    else:
        print(u'  The Tm in the table agrees with the engine.')
    if ihlal:
        print('  PANEL KURALI IHLALI (%d):' % len(ihlal))
        for x in ihlal:
            print('    * %s' % x)
    else:
        print(u'  No violation of the panel geometry rules.')
    if ta_sorun:
        print('  Ta KURALI TUTMUYOR:')
        for x in ta_sorun:
            print('    * %s' % x)

    if a.yaz and (uyusmaz or urun_sapan):
        yed = panel_yolu + '.yedek_%s_geometri' % time.strftime('%H%M')
        io.open(yed, 'w', encoding='utf-8', newline='').write(
            io.open(panel_yolu, encoding='utf-8').read())
        for c, tF, tR, _a, _b, _h in sonuc:
            r = c['satir']
            r[c['iTf']] = '%.2f' % tF
            r[c['iTr']] = '%.2f' % tR
            if c['idT'] is not None and len(r) > c['idT']:
                r[c['idT']] = '%.2f' % abs(tF - tR)
            # Urun boyu yalniz TEK bir olculen deger varsa duzeltilir.
            # Birden cok boy cikiyorsa (evrensel primerler) hangisinin
            # yazilacagi bir KARARDIR ve betik karar vermez.
            if c.get('yeni_urun') and c.get('iU') is not None and len(r) > c['iU']:
                print(u'  product length corrected: %s %s -> %d'
                      % (c['hedef'], r[c['iU']], c['yeni_urun']))
                r[c['iU']] = str(c['yeni_urun'])
        with io.open(panel_yolu, 'w', encoding='utf-8', newline='') as fh:
            for r in sat:
                fh.write(u'\t'.join(r) + u'\n')
        print(u'  the panel table was updated (backup: %s)' % os.path.basename(yed))

    print('=' * 78)
    return 1 if (uyusmaz or ihlal or ta_sorun or urun_sapan) else 0


if __name__ == '__main__':
    sys.exit(main())
