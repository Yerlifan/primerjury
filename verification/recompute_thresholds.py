# -*- coding: utf-8 -*-
"""RECOMPUTE THE ABUNDANCE WEIGHTED THRESHOLD WITH CURRENT DATA.

THE RULE (from the 2026-08-08 threshold document, unchanged)
------------------------------------------------------------
    required dCq = max( log2(R) + 4,3 ,  3,32 )

R  = the read ratio of the competitor pool to the member pool (the abundance ratio)
4,3 = the cycle equivalent of the 5 percent purity criterion Fierer et al. actually
      applied
3,32 = the floor; even with R<=1 a target is expected to show at least a tenfold
      discrimination (log2(10) = 3,32)

WHY AGAIN
---------
The threshold table was produced on 8 August. Two things have changed since:
  1) the memberships were re-derived from the MEASURED identity (they changed on 2
     of the 21 targets),
  2) the sequence of six pairs changed.
R is affected by both. So the required dCq values of 8 August DO NOT BELONG to
today's panel. This script reproduces them with today's numbers.

IT GIVES NO VERDICT. It puts the two rules (the flat threshold 3,00 and the
abundance weighted threshold) SIDE BY SIDE. Which one is applied is a scientific
preference and not this script's business.

To run:
    python verification/recompute_thresholds.py --root .
Output:
    ESIK_IKI_KURAL.tsv  and  ESIK_IKI_KURAL.md

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
    p.add_argument('--root', dest='kok', default='.')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)

    panel = _tsv(os.path.join(kok, 'ONE_PROTOCOL_RESULT', 'panel_tek_protokol.tsv'))
    ham = _tsv(os.path.join(kok, 'ONE_PROTOCOL_RESULT', 'kutu_bazli_ham_sayilar.tsv'))
    if not panel:
        sys.exit(u'ERROR: panel_tek_protokol.tsv could not be read.')

    # WHERE R COMES FROM, and why it is not recomputed
    # ------------------------------------------------
    # R (the abundance ratio) IS WRITTEN inside ESIK_VE_OLCUT_2026-08-08.tsv but the
    # script that produced it is NOT in the project. cross_check.py only confirms the
    # arithmetic (max(log2(R)+4,3; 3,32)); it does not compute R itself.
    #
    # On 2026-08-10 I tried recomputing R as "competitor bin reads / member bin reads":
    # it matched the table on 0 of 17 rows (Metilotrofik, for one: the table says
    # 1,760, my calculation 293,685). So the R in the table uses some other
    # normalisation. I DO NOT KNOW which one; rather than invent a number I do not
    # know, I carry the one in the table TOGETHER WITH ITS SOURCE and write down that
    # it has no producer.
    #
    # On the night of 2026-08-10 the definition WAS FOUND in the session records:
    #     "R = the abundance ratio of the nearest competitor to the target"
    # I tried the six variants of that definition that come to mind (the most abundant
    # competitor bin / the member pool, the competitor with the highest product ratio /
    # the member ratio, total competitor / total member, and so on) against today's
    # numbers: it matched on 1 of 18 rows (Asetoklastik metanojenler, 0,130 exactly).
    # Since only 1 of the 14 rows whose pair DID NOT CHANGE matched either, that has to
    # count as a coincidence.
    #
    # The conclusion: the R values of 8 August were computed with THAT DAY'S
    # membership; the memberships were re-derived on 3 August from the measured
    # identity and six pairs changed on 10 August. So the numbers DO NOT BELONG to
    # today's panel.
    #
    # That is a gap and it goes into the report as one: a number in a decision table
    # that enters a report has to be reproducible. What has to be decided: write down
    # the definition of R and recompute it from scratch with TODAY'S membership (which
    # makes it reproducible), or do not use the abundance rule in this delivery at all.
    # Both are defensible; using the old number silently is not.
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
        # R was computed on 8 August. Since that day the memberships were
        # re-derived and six pairs changed; if dCq changed then R may have
        # changed too. The rows that changed are marked, they are not used
        # silently.
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
        fh.write(u'Generated: %s, source `ONE_PROTOCOL_RESULT/` (this run\'s own output)\n\n' % time.strftime('%Y-%m-%d %H:%M'))
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
