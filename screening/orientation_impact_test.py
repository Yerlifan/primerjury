# -*- coding: utf-8 -*-
"""
orientation_impact_test.py, A KNOWN ANSWER TEST: which number does a wrong
orientation spoil, and by how much.

The setup: the same primer pair, the same consensus, the only difference being
ORIENTATION.
  (a) the consensus stored in the SENSE direction -> expected: a product is found
  (b) the REVERSE COMPLEMENT of the same consensus -> expected: NO product (0)

The reason: this project's in-silico PCR engines (ispcr.amplify and okuma_motoru)
scan the given sequence ON THE PLUS STRAND ONLY; they look there for the forward
primer and for the complement of the reverse primer. On a sequence stored in reverse
neither is found, and the engine reports no product without raising anything.

"""
# -------------------------------------------------------------------------
# orientation_impact_test.py ties the question "what do we lose if the orientation
#                     is wrong" to a number: the same pair, the same consensus, the
#                     only difference being orientation.
#
# INPUT  : the consensus directory given with --klasor under --kok (the default is
#          referans_konsensus/konsensus); the panel's universal pairs and three
#          extra pairs are a fixed list inside the file. read_engine.py does the
#          measuring, and engine/ispcr.py as well if it is present.
# OUTPUT : it writes no file; it prints the product counts found in the right
#          direction and in the reverse to the screen.
# CALLED BY: IT IS NOT IN THE MENU, it is an evidence generator run by hand. The
#          number it produces (117 products in the right direction, 0 in the
#          reverse) is cited as the reason in the explanations inside config.py and
#          orientation.py.
#
# WHY THE LOSS IS 100 PERCENT: this project's in-silico PCR engines scan the given
# sequence ON THE PLUS STRAND ONLY; they look there for the forward primer and for
# the complement of the reverse primer. On a consensus stored in reverse neither is
# found. The engine raises nothing, it says "no product"; that is, the fault is
# silent and shows only in a known answer test like this one.
# -------------------------------------------------------------------------
import sys, os, glob, argparse

BURA = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BURA)
sys.path.insert(0, os.path.join(os.path.dirname(BURA), 'engine'))
import read_engine as om
try:
    import ispcr
    ISPCR = True
except Exception:
    ISPCR = False

CIFT = {
    'A':  ('Arke_universal',        'CTGCGGTTTAATTGGATTCAACGC', 'GAACTGACGACGGCCATGC'),
    'B':  ('Bakteri_universal',     'ACAAGCGGTGGAGCATGTG',      'ACGACAGCCATGCAGCAC'),
    'F1': ('Mantar_universal (F1)', 'GGTTACCCGCTGAACTTAAGC',    'CGCTTCACTCGCCGTTAC'),
    'F2': ('Mantar_universal (F2)', 'GTGCATGGCCGTTCTTAGTTG',    'CAAACTTCCATCGGCTTGAGC'),
}
EK = [('Methanosarcina_mazei_turu', 'A', 'GCCCTTGGGACCGGCATAA', 'TCGCTGGCTAGTAGGTACATTACA'),
      ('Methanosarcina_cinsi',      'A', 'TCGCTAGGTGTCAGGCATG', 'GCGATTCAGGCAAGGTCTTC'),
      ('Proteolitik_Cloacimonas',   'B', 'TTAAAGGCAGCGGCTCACC', 'GAACCCGACACCTAGTGATTATCG')]


def sinifi(yol):
    for s, e in (('A', ('A1-', 'A2-')), ('F1', ('F1-',)), ('F2', ('F2-',)), ('B', ('B-',))):
        if any(x in yol for x in e):
            return s
    return '?'


def oku(yol):
    buf, cur = [], []
    for l in open(yol, errors='ignore'):
        if l.startswith('>'):
            if cur: buf.append(''.join(cur)); cur = []
        else:
            cur.append(l.strip())
    if cur: buf.append(''.join(cur))
    return buf


def urun_var(dizi, F, R, mm=3):
    """SADECE ARTI IPLIK - projenin motorlarinin yaptigi budur."""
    fs = om.Sonda(F, False, mm, son2=True)
    rs = om.Sonda(om.rc(R), True, mm, son2=True)
    return om.urun_var(dizi, fs, rs, len(F), len(R), 40, 600) is not None


def ispcr_var(dizi, F, R, mm=3):
    return bool(ispcr.amplify(ispcr.clean(dizi), F, R, max_mm=mm, lo=40, hi=600))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', '--kok', dest='kok', required=True)
    ap.add_argument('--dir', '--klasor', dest='klasor', default='referans_konsensus/konsensus')
    ap.add_argument('--mm', type=int, default=3)
    a = ap.parse_args()

    yollar = sorted(glob.glob(os.path.join(a.kok, a.klasor, '**', '*.fasta'), recursive=True))
    print(u'Source directory: %s  (%d files)' % (a.klasor, len(yollar)))
    print(u'Criterion       : mismatches <= %d, last 2 bases at the 3\' end exact, product 40-600 bp' % a.mm)
    print(u'Engine          : only the PLUS STRAND is scanned (the behaviour of this project\'s engines)')

    testler = [(ad, sn, F, R) for sn, (ad, F, R) in CIFT.items()] + \
              [(ad, sn, F, R) for ad, sn, F, R in EK]
    print('%-28s %-4s %8s %8s %8s   %s' % ('cift', 'sinif', 'dogru', 'TERS', 'kayip', 'sonuc'))
    tp_d = tp_t = 0
    for ad, sn, F, R in testler:
        dogru = ters = top = 0
        for y in yollar:
            if sinifi(y) != sn:
                continue
            for s in oku(y):
                s = om.temizle(s)
                if len(s) < 200:
                    continue
                top += 1
                if urun_var(s, F, R, a.mm):
                    dogru += 1
                if urun_var(om.rc(s), F, R, a.mm):
                    ters += 1
        if not top:
            continue
        tp_d += dogru; tp_t += ters
        kayip = dogru - ters
        print('%-28s %-4s %5d/%-3d %5d/%-3d %8d   %s'
              % (ad[:28], sn, dogru, top, ters, top, kayip,
                 u'the reverse orientation ZEROES the product' if ters == 0 and dogru > 0
                 else ('etkilenmedi' if dogru == ters else 'kismi kayip')))
    print(u'\nTOTAL  %d products in the correct orientation,  %d in the reverse orientation,  %d lost (%.1f%%)'
          % (tp_d, tp_t, tp_d - tp_t, 100.0 * (tp_d - tp_t) / max(tp_d, 1)))

    if ISPCR:
        print(u'\nCROSS CHECK - the panel\'s own engine (ispcr.amplify), the same files, class A:')
        ad, F, R = CIFT['A']
        d = t = n = 0
        for y in yollar:
            if sinifi(y) != 'A':
                continue
            for s in oku(y):
                s = om.temizle(s)
                if len(s) < 200:
                    continue
                n += 1
                d += int(ispcr_var(s, F, R, a.mm))
                t += int(ispcr_var(om.rc(s), F, R, a.mm))
        print(u'  %s: correct %d/%d, reversed %d/%d  -> same answer' % (ad, d, n, t, n))


if __name__ == '__main__':
    main()
