# -*- coding: utf-8 -*-
"""REFRESH THE QUICK TEST REFERENCES FROM A FULL RUN.

WHY
---
The "reference" values the quick consistency test compares against were constants
embedded in the code, and nowhere did it say which primer pair they had been measured
from. On 2026-08-10 the Bacteroidales pair was changed; the test compared the old
pair's 0.74x against the new pair's 14.23x and said "ZINCIR TUTARSIZ". The chain was
not inconsistent, the reference was stale. Updating the constants by eye would let
the same thing happen again, so the reference is now produced FROM THE OUTPUT OF A
FULL RUN and every row carries the F and R sequences that measurement was made with.

The source: ONE_PROTOCOL_RESULT/panel_tek_protokol.tsv (a full depth run)
The output: QUICK_TEST/referans_degerler.tsv

To run:
    python verification/refresh_reference.py --root .

"""
from __future__ import print_function

import argparse
import csv
import io
import os
import sys
import time


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', dest='kok', default='.')
    p.add_argument('--source', dest='kaynak', default=None,
                   help='default: ONE_PROTOCOL_RESULT/panel_tek_protokol.tsv')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)
    kaynak = a.kaynak or os.path.join(kok, 'ONE_PROTOCOL_RESULT',
                                      'panel_tek_protokol.tsv')
    if not os.path.exists(kaynak):
        sys.exit(u'ERROR: there is no source: %s\n      the full run (stage P) has to be finished first.' % kaynak)

    with io.open(kaynak, encoding='utf-8') as fh:
        satirlar = list(csv.DictReader(
            (l for l in fh if not l.startswith('#')), delimiter='\t'))
    if not satirlar:
        sys.exit(u'ERROR: the source is empty: %s' % kaynak)

    hedef_dizin = os.path.join(kok, 'QUICK_TEST')
    if not os.path.isdir(hedef_dizin):
        os.makedirs(hedef_dizin)
    cikti = os.path.join(hedef_dizin, 'referans_degerler.tsv')

    n = 0
    atlanan = []
    with io.open(cikti, 'w', encoding='utf-8', newline='') as fh:
        fh.write("""# THE QUICK TEST REFERENCE VALUES
""")
        fh.write("""# Produced: %s
""" % time.strftime('%Y-%m-%d %H:%M'))
        fh.write("""# Source  : %s
""" % os.path.relpath(kaynak, kok))
        fh.write(u'# Every row carries the F/R SEQUENCE the measurement was made with. The quick\n# test compares the sequence first; if the pair changed it does NOT\n# compare against the old number, it says "reference stale" and does\n# not stop the chain. That is why hard coded numbers were abandoned.\n')
        fh.write(u'hedef\tF\tR\treferans_x\tkarar\tkapsam\n')
        for r in satirlar:
            ad = (r.get('hedef') or '').strip()
            F = (r.get('F') or '').strip().upper()
            R = (r.get('R') or '').strip().upper()
            x = (r.get('ASIL_ayrim_mm1') or '').strip()
            if not ad or not F or not R or not x:
                if ad:
                    atlanan.append(ad)
                continue
            fh.write(u'%s\t%s\t%s\t%s\t%s\t%s\n'
                     % (ad, F, R, x.replace(',', '.'),
                        (r.get('esik_gecti_mi') or '').strip(),
                        (r.get('ASIL_kapsam_mm1') or '').strip()))
            n += 1

    print('yazildi: %s' % cikti)
    print(u'  %d rows recorded as the reference' % n)
    if atlanan:
        print(u'  skipped (the sequence or the measurement is empty): %s' % ', '.join(atlanan))
    return 0 if n else 1


if __name__ == '__main__':
    sys.exit(main())
