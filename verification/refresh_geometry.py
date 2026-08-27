# -*- coding: utf-8 -*-
"""RUN THE GEOMETRY GATE AGAIN WITH THE PANEL SEQUENCES AS THEY ARE NOW.

WHY (2026-08-10, caught in an audit)
------------------------------------
The sequence of SIX pairs in the panel has changed since the geometry audit of 2
August, and the geometry file still carries the OLD sequences. So those six pairs
have NEVER PASSED the panel's own geometry rules (length 18-25, GC 40-60 percent, Tm
58-62, hairpin, dimer, the 3' end).

On top of that, the Tm column of the delivery table says 48.9 to 51.1 for these five
new pairs. That does not agree with the panel's engine (primer3, mv=50 dv=1.5
dntp=0.6 dna=50): an independent nearest neighbour calculation gives 55 to 58 for the
same five primers, and while the difference between the table and the independent
calculation is a steady +3.72 +- 0.22 C across the 16 sound pairs, on these five it
is -6.38 C. Two separate populations; they cannot have come out of a single method.

Why it matters: the delivery table says "Ta = min(Tm) - 3". With the Tm values as
written, Ta comes out 47 for P1, while the plate says 55. Anyone reading the table
and recomputing Ta would run the whole panel at the wrong temperature. With a
consistent Tm the same rule gives the written Ta within 0.6 C on ALL FIVE plate
groups: THE DESIGN IS RIGHT, THE NUMBER WRITTEN IS WRONG.

This script gives no verdict, it MEASURES: it recomputes every primer with the
panel's own engine, lists the rule violations and compares them against the value
written in the table.

To run:
    python verification/refresh_geometry.py --root .
    python verification/refresh_geometry.py --root . --write  (update the panel table too)

"""
from __future__ import print_function

import argparse
import io
import os
import re
import sys
import time

PANEL = os.path.join('final_primers', 'devir_ciftleri_20260802_sonrotus_TESLIM.tsv')


def geo_yukle(kok):
    """Brings in the panel's OWN geometry module. A separate copy IS NOT WRITTEN: two
    copies drift apart in time and there is no telling which one is right.

    """
    for aday in ('engine', 'engine'):
        d = os.path.join(kok, aday)
        if os.path.exists(os.path.join(d, 'geometry_core.py')):
            sys.path.insert(0, d)
            try:
                import geometry_core
            except ImportError as e:
                # primer3 yoksa YIGIN IZI basip kullaniciyi korkutma; ne
                # eksik oldugunu ve nasil kurulacagini soyle.
                return None, ('geometry_core.py was found (%s) but could not '
                              'be loaded: %s       To install: pip3 install '
                              'primer3-py --break-system-packages'
                              % (os.path.join(aday, 'geometry_core.py'), e))
            return geo, os.path.join(aday, 'geometry_core.py')
    return None, 'geometry_core.py was not found under engine'


# --- THE PRODUCT LENGTH: measured, not written ------------------------------
# Like the Tm, the product length is a number that can be written by hand, and every
# number written by hand drifts one day. 2026-08-10: the
# Proteolitik_Synergistaceae row said 173 bp and the measurement says 172. One base,
# but it enters the gel separation and the band class calculations. The criterion is
# the panel's own: mm<=1, the last two bases at the 3' end exact.
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
    """The canonical consensuses, from INDEX.tsv and NOT from a glob.

    There are 250 files in the canonical_consensus directory and only 100 of them are
    valid. The other 150 are leftovers that cannot be deleted on the mounted directory;
    in 33 bins two or three versions of the same bin with DIFFERENT content sit side by
    side. Reading with a glob leaves the choice of version to file name order, so we
    would be measuring with a sequence the panel never saw. That is why the panel's own
    loader reads the index too.

    """
    import csv as _csv
    d = os.path.join(kok, 'canonical_consensus')
    ix = os.path.join(d, 'INDEX.tsv')
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
    'The part that runs even without primer3: measuring the product length.'
    kons = konsensus_yukle(kok)
    if not kons:
        print(u'  canonical_consensus is empty or missing, so product length could not be measured either.')
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
            print('      %-44s the table says %d bp and the measured value is '
                  '%s'
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
    p.add_argument('--root', dest='kok', default='.')
    p.add_argument('--write', dest='yaz', action='store_true',
                   help='also update the Tm and dTm columns of the panel '
                        'table; a backup is taken first')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)

    geo, geo_yol = geo_yukle(kok)
    # Without primer3 the Tm cannot be measured but THE PRODUCT LENGTH can:
    # that measurement is made on the consensus sequences and does not depend
    # on primer3. Skipping a measurement that can be made because of one
    # missing piece leaves the audit needlessly blind.
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
    print('  motor: %s   ayarlar: %s' % (geo_yol, geometry_core.KW))
    print('=' * 78)

    ciftler, bas, sat, panel_yolu = panel_oku(kok)
    print(u'  pairs to measure in the panel: %d' % len(ciftler))

    sonuc = []
    ihlal = []
    uyusmaz = []
    for c in ciftler:
        tF, tR = geometry_core.tm(c['F']), geometry_core.tm(c['R'])
        vF, vR = geometry_core.viol(c['F']), geometry_core.viol(c['R'])
        het = geometry_core.het(c['F'], c['R'])
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
                uyusmaz.append('%s %s: the table says %s and the engine %.2f, '
                               'a difference of %.2f'
                               % (c['hedef'], et, yazili, yeni, f))

    print()
    print('  %-44s %-15s %-15s %s'
          % ('target', 'table F/R', 'engine F/R', 'violation'))
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
                urun_sapan.append('%s: the table says %d bp and the measured '
                                  'value is %s'
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
        print(u'  product length NOT MEASURED: the canonical_consensus directory is empty or missing.')

    cy = os.path.join(kok, 'final_primers',
                      'geometri_denetimi_%s.tsv' % time.strftime('%Y%m%d'))
    with io.open(cy, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# The geometry audit, with the panel sequences AS THEY ARE NOW. Produced %s\n'
                 % time.strftime('%Y-%m-%d %H:%M'))
        fh.write(u'# Motor: %s  ayarlar: %s\n' % (geo_yol, geometry_core.KW))
        fh.write(u'Hedef\tPrimer\tDizi\tUz\tGC%\tTm\tHairpin Tm\tHomodimer Tm\tIhlal\n')
        for c, tF, tR, vF, vR, het in sonuc:
            for et, d, t, v in (('Ileri', c['F'], tF, vF), ('Geri', c['R'], tR, vR)):
                fh.write(u'%s\t%s\t%s\t%d\t%.1f\t%.2f\t%.1f\t%.1f\t%s\n'
                         % (c['hedef'], et, d, len(d), geometry_core.gc(d), t,
                            geometry_core.hp(d), geometry_core.hd(d), '; '.join(v)))
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
            # The product length is corrected only if there is ONE measured value.
            # If several lengths come out (the universal primers), which one to
            # write is A DECISION and the script does not decide.
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
