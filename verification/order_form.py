# -*- coding: utf-8 -*-
"""THE ORDER FORM: the ONE correct list to send the oligo supplier.

WHY IT EXISTS, caught in an audit
---------------------------------
The morning summary said "copy the sequences from here" and pointed at the panel
sheet of an older workbook. The sequences of SIX pairs on that sheet were OUT OF
DATE, so ordering from that file would have brought six of the twenty pairs as
THE WRONG oligo. On top of that the geometry audit file carried a different six.
Three separate files, three separate versions of "the right one", and nothing
that said which was right.

This script closes that gap: the order list is PRODUCED FROM THE PANEL'S OWN
SOURCE rather than written by hand, and the source file's md5 is written on every
run. If the source changes the form changes; if it does not, the form comes out
the same.

To run it:
    python verification/order_form.py --root .
The output:
    the order form as a TSV, to paste for the supplier
    the order form as markdown, for a person to read
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

PANEL = os.path.join('final_primers', 'devir_ciftleri_20260802_sonrotus_TESLIM.tsv')


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
    p.add_argument('--root', dest='kok', default='.')
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
    sl_yol = os.path.join(kok, 'ONE_PROTOCOL_RESULT', 'SIPARIS_LISTESI.tsv')
    sl = _tsv(sl_yol)
    if not sl:
        sys.exit(u'ERROR: %s could not be read.' % sl_yol)

    girer, girmez, uyari = [], [], []
    for r in sl:
        ad = (r.get('hedef') or '').strip()
        F = (r.get('F') or '').strip().upper()
        R = (r.get('R') or '').strip().upper()
        if not ad:
            continue
        if not re.fullmatch(r'[ACGT]+', F or '') or not re.fullmatch(r'[ACGT]+', R or ''):
            uyari.append(u'%s: there is no sequence in the order list, or it holds a character outside A/C/G/T, so it WAS NOT TAKEN into the form' % ad)
            continue
        p = PN.get(ad)
        if p is None:
            uyari.append(u'%s: it IS in the order list but it HAS NO ROW in the panel source. The plate and Ta of this pair are unknown.' % ad)
        elif (p['F'], p['R']) != (F, R):
            uyari.append(u'%s: the order list and the panel source give DIFFERENT sequences. The form used the order list.' % ad)
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
        fh.write(u'# THE ORDER FORM - %s\n' % time.strftime('%Y-%m-%d %H:%M'))
        fh.write("""# The source: %s (md5 %s)
""" % (PANEL.replace('\\', '/'), ozet))
        fh.write(u'# This file is NEVER WRITTEN BY HAND; it is generated from the panel source.\n')
        fh.write(u'# DO NOT COPY FROM THE xlsx FILE: as of 2026-08-10 the sequence of six\n# pairs there IS OUT OF DATE (Bacteroidales, Bakteri_universal,\n# Mantar F1, Microascaceae, Petriella_musispora, Petrimonas).\n')
        for u_ in uyari:
            fh.write(u'# WARNING: %s\n' % u_)
        fh.write(u'oligo_adi\tdizi_5_3\tuzunluk\thedef\tyon\tplaka\tTa_C\turun_bp\tsinif\n')
        for c in girer:
            for ad, d, yon in ((c['adF'], c['F'], 'ileri'), (c['adR'], c['R'], 'geri')):
                fh.write(u'%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n'
                         % (ad, d, len(d), c['hedef'], yon, c['plaka'], c['ta'],
                            c['urun'], c['sinif']))

    my = os.path.join(kok, 'SIPARIS_FORMU.md')
    with io.open(my, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# Order form\n\n')
        fh.write(u'Generated: %s, source `%s` (md5 `%s`)\n\n'
                 % (time.strftime('%Y-%m-%d %H:%M'), PANEL.replace('\\', '/'), ozet))
        fh.write(u'> **Do not copy from the Excel file.** The `2 Panel` sheet inside `PrimerJury_PCR_Paneli_2026-08-02_TESLIM.xlsx` ')
        fh.write(u'**%d pairs = %d oligos**\n\n' % (len(girer), 2 * len(girer)))
        if uyari:
            fh.write(u'### Warnings\n\n')
            for u_ in uyari:
                fh.write(u'- %s\n' % u_)
            fh.write(u'\n')
        fh.write(u'| oligo name | sequence (5->3) | len | target | orientation | plate | Ta | product |\n|---|---|---|---|---|---|---|---|\n')
        for c in girer:
            for ad, d, yon in ((c['adF'], c['F'], 'ileri'), (c['adR'], c['R'], 'geri')):
                fh.write(u'| %s | `%s` | %d | %s | %s | %s | %s | %s |\n'
                         % (ad, d, len(d), c['hedef'], yon, c['plaka'], c['ta'], c['urun']))
        if girmez:
            fh.write(u'\n## Not ordered (kept, not deleted, for reference)\n\n')
            fh.write(u'| target | class | product |\n|---|---|---|\n')
            for c in girmez:
                fh.write(u'| %s | %s | %s |\n' % (c['hedef'], c['sinif'], c['urun']))
        fh.write(u'\n## What this form does not tell you\n\n')
        fh.write(u'This form tells you **which sequences to order**. It does not tell you that these pairs passed the geometry gate, that the in-plate gel separation is clean, or which threshold rule was applied. Those live in `ONE_KEY_RESULT/DENETIM_RAPORU.md` and `GECE_BULGULARI.md`, and both should be read before an order is placed.\n')

    print(u'written: %s' % ty)
    print(u'written: %s' % my)
    print(u'  going into the order: %d pairs = %d oligos' % (len(girer), 2 * len(girer)))
    print(u'  not ordered        : %d pairs' % len(girmez))
    print('  the source md5 : %s' % ozet)
    return 0


if __name__ == '__main__':
    sys.exit(main())
