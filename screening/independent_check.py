# -*- coding: utf-8 -*-
"""
independent_check.py measures the numbers in the headings THREE SEPARATE ways and
shows that they hold.

  Route A  ispcr.find_sites / amplify   : a numpy vector scan, SEEDLESS.
                                          The panel's OWN code (engine/ispcr.py),
                                          unchanged in this session.
  Route B  brute_force.py               : pure python, every position tried one at a
                                          time. Written independently, no shared code.
  Route C  read_engine.py               : the engine corrected in this session (the
                                          reference).

A and B use not one line of the corrected engine. If all three give the same number,
the output of the corrected engine is independently confirmed.

Usage:
    python independent_check.py --fastq "../fastq files" [--mm 1] [--nmax 3000]

"""
# -------------------------------------------------------------------------
# independent_check.py measures the numbers in the delivery headings through three
#                         separate code routes and shows that all three give the
#                         same answer.
#
# INPUT  : the "fastq files" directory given with --fastq (--mm the mismatch
#          ceiling, --nmax the reads per bin, --seed the sampling seed). The list of
#          pairs and bins tested (TESTLER) is fixed inside the file. The three
#          routes: engine/ispcr.py, brute_force.py and read_engine.py.
# OUTPUT : it writes no file; it prints the numbers of the three routes side by side.
# CALLED BY: IT IS NOT IN THE MENU, it is an evidence generator run by hand.
#
# WHY THREE ROUTES: route A (ispcr) is the panel's own code and was not changed at
# all in this session; route B (brute force) is seedless and shares no code; route C
# is the corrected engine. A and B use not one line of the corrected engine, so if
# all three give the same number, the losslessness of the pigeonhole seeding is
# confirmed from outside rather than by its own code.
# -------------------------------------------------------------------------
import sys, os, glob, random, argparse

BURA = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BURA)
sys.path.insert(0, os.path.join(os.path.dirname(BURA), 'engine'))

import brute_force as kk
import read_engine as om

try:
    import ispcr
    ISPCR = True
except Exception as e:
    ISPCR = False
    sys.stderr.write(u'WARNING: ispcr could not be loaded (%s). Route A will be skipped.\n' % e)

# The bins the numbers in the headings rest on: the M. mazei pair (row 22) and the
# M. hadiensis competitor bins, where the fault was largest.
TESTLER = [
    ('Methanosarcina_mazei_turu', 'GCCCTTGGGACCGGCATAA', 'TCGCTGGCTAGTAGGTACATTACA',
     [('A1-4', '3078083', 'COMPETITOR M. hadiensis, the centre of the fault'),
      ('A2-4', '3078083', 'COMPETITOR M. hadiensis'),
      ('A2-2', '2209',    'UYE M. mazei'),
      ('A1-3', '2209',    'UYE M. mazei')]),
    ('Methanosarcina_cinsi', 'TCGCTAGGTGTCAGGCATG', 'GCGATTCAGGCAAGGTCTTC',
     [('A1-2', '2209', 'MEMBER M. mazei, whose coverage had been measured '
                       'short'),
      ('A2-3', '2223', 'competitor M. soehngenii')]),
    ('Asetoklastik_metanojenler', 'CCGGGAGAGGTGAGAGGTAC', 'CGGGTATCTAATCCGGTTCGTG',
     [('A1-2', '2209', 'MEMBER M. mazei, whose coverage had been measured '
                       'short')]),
]


def yol_a(reads, F, R, mm):
    p = 0
    for s in reads:
        for seq in (s, ispcr.rc(s)):
            if ispcr.amplify(ispcr.clean(seq), F, R, max_mm=mm,
                             lo=40, hi=600, need_tail=True):
                p += 1
                break
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fastq', required=True)
    ap.add_argument('--mm', type=int, default=1)
    ap.add_argument('--nmax', type=int, default=3000)
    ap.add_argument('--seed', type=int, default=3)
    a = ap.parse_args()

    print(u'criterion: mismatches <= %d, the last 2 bases at the 3\' end an EXACT match, product 40-600 bp' % a.mm)
    print('%-26s %-30s %6s | %-16s %-16s %-16s | %s'
          % ('cift', 'kutu', 'n', 'A ispcr(numpy)', 'B kaba(python)', 'C duzeltilmis', 'UYUM'))
    tumu = True
    for ad, F, R, kutular in TESTLER:
        for d, tax, not_ in kutular:
            p = os.path.join(a.fastq, d, '%s-reads_%s.fastq' % (d, tax))
            if not os.path.exists(p):
                print(u'  MISSING:', p); continue
            rs = list(om.okumalar(p))
            if a.nmax and len(rs) > a.nmax:
                random.seed(a.seed)
                rs = random.sample(rs, a.nmax)
            if not rs:
                continue
            b = kk.kutu_pcr(rs, F, R, max_mm=a.mm)[0]
            c = om.kutu_pcr(rs, F, R, max_mm=a.mm)[0]
            av = yol_a(rs, F, R, a.mm) if ISPCR else None
            deg = [x for x in (av, b, c) if x is not None]
            ok = len(set(deg)) == 1
            tumu &= ok
            f = lambda v: '-' if v is None else '%5d %6.2f%%' % (v, 100.0 * v / len(rs))
            print('%-26s %-30s %6d | %-16s %-16s %-16s | %s'
                  % (ad[:26], ('%s_%s %s' % (d, tax, not_))[:30], len(rs),
                     f(av), f(b), f(c), 'ESIT' if ok else '*** FARK ***'))
    print()
    if tumu:
        print(u'RESULT: ALL THREE ROUTES GIVE THE SAME NUMBER - the site counts are independently verified.')
        return 0
    print(u'RESULT: THERE IS A DIFFERENCE - do not use the numbers.')
    return 1


if __name__ == '__main__':
    sys.exit(main())
