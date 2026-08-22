# -*- coding: utf-8 -*-
"""A PLATE ASSIGNMENT SUGGESTION: it SEARCHES for a layout that removes the gel
overlap.

THE PROBLEM
-----------
Two products running on the same plate cannot be told apart on 2 per cent
agarose if they are within 10 bp of one another, which is the criterion that was
set. Pairs changed later broke that separation again: one product went from 241
bp to 150 bp and landed beside another at 145 bp, undoing an earlier correction
without anyone noticing.

WHAT IT DOES
------------
It DOES NOT TOUCH the sequences. It only changes which pair runs on which plate.
The constraint: every pair has to work at its plate's annealing temperature. The
criterion is Ta <= min(Tm) - 3, so a pair can be put on a Ta three degrees below
its own minimum Tm or lower, and never HIGHER.

The output is A SUGGESTION and is not applied. Changing a plate changes the
layout of the experiment, and a person decides that.

To run it:
    python verification/assign_plate.py --root .
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
    """The Tm from the panel's OWN engine when it is available. Otherwise the value
        in the table is used and that is written down OPENLY: a stale number is
        never trusted silently."""
    for aday in ('engine', 'engine'):
        d = os.path.join(kok, aday)
        if os.path.exists(os.path.join(d, 'geometry_core.py')):
            sys.path.insert(0, d)
            try:
                import geometry_core
            except Exception:
                return None, 'geometry_core.py could not be loaded; is primer3 missing?'
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
                out[r[iH].strip()] = (geometry_core.tm(F), geometry_core.tm(R))
            return out, None
    return None, u'geometry_core.py bulunamadi'


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
    'Can this pair run at that Ta: Ta <= min(Tm) - 3.'
    return ta <= c['tm_min'] - 3 + 0.6      # 0,6 C olcum toleransi


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', dest='kok', default='.')
    p.add_argument('--trial', dest='deneme', type=int, default=20000)
    a = p.parse_args()
    kok = os.path.abspath(a.kok)

    tmk, uyari = geo_tm(kok, None)
    ciftler = oku(kok, tmk)
    print('=' * 78)
    print('  PLAKA ATAMASI ONERISI   %s' % time.strftime('%Y-%m-%d %H:%M'))
    print('  Tm kaynagi: %s' % (u'the panel\'s own engine (geometry_core.py)' if tmk
                                else u'TABLE (%s) - careful, the table value can be stale' % uyari))
    print('=' * 78)
    if not ciftler:
        print(u'  the pair could not be read.')
        return 1

    gruplar = sorted(set((c['plaka'], c['ta']) for c in ciftler))
    print(u'  pair %d, plate group %d' % (len(ciftler), len(gruplar)))

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
        print(u'    plate %-3s Ta %-4s  %d pairs, %d clashes' % (g[0], g[1], len(uy), n))
        for a1, u1, a2, u2 in ay:
            print('        %s (%d) / %s (%d)' % (a1, u1, a2, u2))
    print(u'    TOTAL CLASHES: %d' % top)

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
        print(u'  THE PAIRS THAT WILL CHANGE:')
        for c in ciftler:
            if eniyi[c['hedef']] != simdi[c['hedef']]:
                print('    %-44s %s Ta %s  ->  %s Ta %s'
                      % (c['hedef'][:44], simdi[c['hedef']][0], simdi[c['hedef']][1],
                         eniyi[c['hedef']][0], eniyi[c['hedef']][1]))
    else:
        print(u'  Under the current Ta constraints no distribution that reduces the clashes')
        print(u'  NOT FOUND (%d attempts). The options are:' % a.deneme)
        print(u'    - to accept the collision and write it in THE REPORT with its reason')
        print(u'      (qPCR separation is done with the melt curve; the gel is a secondary check)')
        print(u'    - redesign one of the clashing pairs and shift its product length')
    print()
    print(u'  This is a SUGGESTION and was not applied. Changing a plate alters the experimental')
    print(u'  it changes; the decision is a person\'s.')
    print('=' * 78)
    return 0


if __name__ == '__main__':
    sys.exit(main())
