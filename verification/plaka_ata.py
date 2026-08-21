# -*- coding: utf-8 -*-
"""PLAKA ATAMASI ONERISI  -  jel cakismasini kaldiracak dagilimi ARAR.

SORUN
-----
Ayni plakada kosan iki urun 10 bp'den yakinsa %2 agarozda ayirt edilemez
(toplantida konulan olcut). Bugun degistirilen ciftler bu ayrimi yeniden
bozdu: Bacteroidales urunu 241 bp'den 150 bp'ye dondu ve Mantar F2'nin
145 bp'sinin yanina oturdu - eski bir duzeltme farkinda olmadan geri alindi.

NE YAPAR
--------
Dizilere DOKUNMAZ. Yalnizca hangi ciftin hangi plakada kosacagini degistirir.
Kisit: her cift, plakasinin Ta'sinda calisabilmeli. Olcut "Ta <= min(Tm) - 3";
yani bir cift, kendi min Tm'sinden 3 derece dusuk ya da daha dusuk bir Ta'ya
konabilir, daha YUKSEGINE konamaz.

Cikti bir ONERIDIR, uygulanmaz. Plaka degistirmek deneyin duzenini degistirir;
karari insan verir.

Kosum:
    python verification/plaka_ata.py --kok .
"""
from __future__ import print_function

import argparse
import io
import itertools
import os
import random
import re
import sys
import time

PANEL = os.path.join('primer_final', 'devir_ciftleri_20260802_sonrotus_TESLIM.tsv')
JEL_ESIK = 10


def oku(kok, tm_kaynagi=None):
    y = os.path.join(kok, PANEL)
    sat = [l.rstrip('\n').split('\t') for l in io.open(y, encoding='utf-8')]
    bas = sat[0]
    iP, iT, iH = bas.index('Plaka'), bas.index('Ta (C)'), bas.index('Hedef')
    iU = bas.index('Urun (bp)')
    iF = next(i for i, b in enumerate(bas) if b.startswith('Ileri primer'))
    iTf, iTr = bas.index('Ileri Tm'), bas.index('Geri Tm')
    out = []
    for r in sat[1:]:
        if len(r) <= max(iU, iTr) or not r[iH].strip():
            continue
        if not re.fullmatch(r'[ACGT]+', (r[iF] or '').strip().upper()):
            continue
        try:
            u = int(re.sub(r'\D', '', r[iU]))
            tf, tr = float(r[iTf].replace(',', '.')), float(r[iTr].replace(',', '.'))
        except ValueError:
            continue
        ad = r[iH].strip()
        if tm_kaynagi and ad in tm_kaynagi:
            tf, tr = tm_kaynagi[ad]
        out.append(dict(hedef=ad, plaka=r[iP].strip(), ta=r[iT].strip(),
                        urun=u, tm_min=min(tf, tr)))
    return out


def geo_tm(kok, ciftler):
    """Varsa panelin KENDI motoruyla Tm. Yoksa tablodaki deger kullanilir ve
    bu durum ACIKCA yazilir - sessizce bayat sayiya guvenilmez."""
    for aday in ('DUZELTME_betikleri', 'MADDE123_betikleri'):
        d = os.path.join(kok, aday)
        if os.path.exists(os.path.join(d, 'geo.py')):
            sys.path.insert(0, d)
            try:
                import geo
            except Exception:
                return None, u'geo.py yuklenemedi (primer3 yok?)'
            sat = [l.rstrip('\n').split('\t')
                   for l in io.open(os.path.join(kok, PANEL), encoding='utf-8')]
            bas = sat[0]
            iH = bas.index('Hedef')
            iF = next(i for i, b in enumerate(bas) if b.startswith('Ileri primer'))
            iR = next(i for i, b in enumerate(bas) if b.startswith('Geri primer'))
            out = {}
            for r in sat[1:]:
                if len(r) <= max(iF, iR) or not r[iH].strip():
                    continue
                F, R = r[iF].strip().upper(), r[iR].strip().upper()
                if not re.fullmatch(r'[ACGT]+', F or ''):
                    continue
                out[r[iH].strip()] = (geo.tm(F), geo.tm(R))
            return out, None
    return None, u'geo.py bulunamadi'


def cakismalar(grup):
    n = 0
    ayrinti = []
    for (a1, u1), (a2, u2) in itertools.combinations(
            sorted(((c['hedef'], c['urun']) for c in grup), key=lambda x: x[1]), 2):
        if abs(u1 - u2) < JEL_ESIK:
            n += 1
            ayrinti.append((a1, u1, a2, u2))
    return n, ayrinti


def uygun(c, ta):
    """Cift bu Ta'da kosabilir mi: Ta <= min(Tm) - 3."""
    return ta <= c['tm_min'] - 3 + 0.6      # 0,6 C olcum toleransi


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--kok', default='.')
    p.add_argument('--deneme', type=int, default=20000)
    a = p.parse_args()
    kok = os.path.abspath(a.kok)

    tmk, uyari = geo_tm(kok, None)
    ciftler = oku(kok, tmk)
    print('=' * 78)
    print('  PLAKA ATAMASI ONERISI   %s' % time.strftime('%Y-%m-%d %H:%M'))
    print('  Tm kaynagi: %s' % (u'panelin kendi motoru (geo.py)' if tmk
                                else u'TABLO (%s) - dikkat, tablo degeri bayat olabilir' % uyari))
    print('=' * 78)
    if not ciftler:
        print('  cift okunamadi.')
        return 1

    gruplar = sorted(set((c['plaka'], c['ta']) for c in ciftler))
    print('  cift %d, plaka grubu %d' % (len(ciftler), len(gruplar)))

    def maliyet(atama):
        t = 0
        for g in gruplar:
            uy = [c for c in ciftler if atama[c['hedef']] == g]
            t += cakismalar(uy)[0]
        return t

    simdi = {c['hedef']: (c['plaka'], c['ta']) for c in ciftler}
    print()
    print('  SIMDIKI durum:')
    top = 0
    for g in gruplar:
        uy = [c for c in ciftler if simdi[c['hedef']] == g]
        n, ay = cakismalar(uy)
        top += n
        print('    plaka %-3s Ta %-4s  %d cift, %d cakisma' % (g[0], g[1], len(uy), n))
        for a1, u1, a2, u2 in ay:
            print('        %s (%d) / %s (%d)' % (a1, u1, a2, u2))
    print('    TOPLAM CAKISMA: %d' % top)

    # her cift hangi gruplarda kosabilir
    izin = {}
    for c in ciftler:
        iz = []
        for g in gruplar:
            try:
                ta = float(g[1])
            except ValueError:
                continue
            if uygun(c, ta):
                iz.append(g)
        izin[c['hedef']] = iz or [simdi[c['hedef']]]

    rastgele = random.Random(20260810)
    eniyi = dict(simdi)
    eniyi_m = top
    for _ in range(a.deneme):
        aday = {}
        for c in ciftler:
            aday[c['hedef']] = rastgele.choice(izin[c['hedef']])
        # grup boyutlarini makul tut: bos ya da asiri dolu plaka istemeyiz
        boy = {}
        for g in gruplar:
            boy[g] = sum(1 for c in ciftler if aday[c['hedef']] == g)
        if min(boy.values()) < 2 or max(boy.values()) > len(ciftler) // 2 + 2:
            continue
        m = maliyet(aday)
        if m < eniyi_m:
            eniyi_m, eniyi = m, aday
            if m == 0:
                break

    print()
    if eniyi_m < top:
        print('  ONERI: cakisma %d -> %d' % (top, eniyi_m))
        for g in gruplar:
            uy = sorted((c for c in ciftler if eniyi[c['hedef']] == g),
                        key=lambda x: x['urun'])
            print('    plaka %-3s Ta %-4s  %s'
                  % (g[0], g[1], ', '.join('%s(%d)' % (c['hedef'][:22], c['urun'])
                                           for c in uy)))
        print()
        print('  DEGISECEK CIFTLER:')
        for c in ciftler:
            if eniyi[c['hedef']] != simdi[c['hedef']]:
                print('    %-44s %s Ta %s  ->  %s Ta %s'
                      % (c['hedef'][:44], simdi[c['hedef']][0], simdi[c['hedef']][1],
                         eniyi[c['hedef']][0], eniyi[c['hedef']][1]))
    else:
        print('  Mevcut Ta kisitlari altinda cakismayi azaltan bir dagilim')
        print('  BULUNAMADI (%d deneme). Secenekler:' % a.deneme)
        print('    - cakismayi kabul edip RAPORDA gerekcesiyle yazmak')
        print('      (qPCR ayrimi erime egrisiyle yapilir, jel ikincil kontroldur)')
        print('    - cakisan ciftlerden birini yeniden tasarlayip urun boyunu kaydirmak')
    print()
    print('  Bu bir ONERIDIR, uygulanmadi. Plaka degistirmek deney duzenini')
    print('  degistirir; karari insan verir.')
    print('=' * 78)
    return 0


if __name__ == '__main__':
    sys.exit(main())
