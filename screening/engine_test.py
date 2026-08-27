# -*- coding: utf-8 -*-
"""
engine_test.py proves that read_engine.py is LOSSLESS.

The method: `okuma_motoru.Sonda` (pigeonhole seeding plus full verification) and
`kaba_kuvvet.yerler` (seedless, every position one at a time) are run on the SAME
sequences under the SAME criterion and the lists of binding sites are compared
EXACTLY. A single site of difference fails the test.

Three layers:
  T1  synthetic sequences: a mismatch is placed at each primer position in turn
      (a mismatch falling inside the seed is exactly the case the old engine missed)
  T2  real reads: the panel's 21 pairs, on a subset of the given fastq files
  T3  product level: are the kutu_pcr counts equal in the two implementations
It also reports HOW MANY sites the OLD engine misses in the same test (the evidence
for the fault).

Usage:
    python engine_test.py                          # T1 only (synthetic, no data needed)
    python engine_test.py "fastq files/A1-4/*.fastq" [--n 300] [--mm 1]
Exit code 0 = passed, 1 = failed.

"""
# -------------------------------------------------------------------------
# engine_test.py proves that read_engine.py's pigeonhole seeding is LOSSLESS by
#                 comparing it exactly against seedless brute force.
#
# INPUT  : T1 needs no data (it produces synthetic sequences); T2 and T3 take the
#          fastq files or a wildcard pattern given on the command line. It imports
#          read_engine.py and brute_force.py directly as modules; the panel's 21
#          pairs sit in the file as a fixed list.
# OUTPUT : it writes no file; it prints the result to the screen. Exit code 0 =
#          passed, 1 = failed (a single binding site difference fails the test).
# CALLED BY: IT IS NOT IN THE MENU, it is a test run by hand. A reduced form of the
#          same comparison is made automatically inside self_test.py on every run
#          (key 8 and the head of every measuring key).
#
# THE PIGEONHOLE PRINCIPLE IS TESTED HERE: if a primer of length k is split into
# (m + 1) non-overlapping pieces in a scan allowing at most m mismatches, at least
# one piece has to match exactly; that is a guarantee. T1 aims at exactly that
# guarantee: it puts the mismatch at EVERY position of the primer one at a time,
# because the case the old engine missed was the mismatch falling inside the seed.
# -------------------------------------------------------------------------
import sys, os, glob, random, argparse, itertools

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import read_engine as om
import brute_force as kk

CIFTLER = [
    ('Metanojen_universal',              'GTGGAGCTTGCGGTTTAATTG',    'CAGGATGCTTCACAGTACGAAC'),
    ('Methanothrix_cinsi',               'GAGAGGTACTTCAGGGGTAGG',    'CTAGCTTTCGTCCCTTGCC'),
    ('Petrimonas_cinsi',                 'AAGTCGCGTGAAGGATGAAG',     'AAAATTTCACCGCCGACTTAAC'),
    ('Metanomikrobiyales_hidrojenotrof', 'TGGGACCGCCTCTGCTAAAG',     'CATTGTAGCCCGCGTGTAGC'),
    ('Mantar_universal_F1',              'GGTTACCCGCTGAACTTAAGC',    'CGCTTCACTCGCCGTTAC'),
    ('Methanosarcina_cinsi',             'TCGCTAGGTGTCAGGCATG',      'GCGATTCAGGCAAGGTCTTC'),
    ('Metilotrofik_metanojen',           'CAATCCTGAAACCCGTCCATAG',   'ATATTCACCGCCTGATGTTGAC'),
    ('Nitrosocosmicus_AOA',              'ACTCTGAGTGATTTCCGTTAAGG',  'TGCTTTAGGCCCAATAAACGTC'),
    ('Proteiniphilum_cinsi',             'GGTTCCTTGAGTGTGGATGAGG',   'CTTGAGCGTCAGTTATGGCTTAG'),
    ('Proteolitik_Synergistaceae',       'AGCTAGTAGGTTGGGTAACGG',    'GATTTCTTCACCCACGCGG'),
    ('Bacteroidales_kumesi',             'GAAGCTAGGATTTGGTTGCTGTG',  'CTCCCCAGGTGGATAACTTATCG'),
    ('Mantar_universal_F2',              'GTGCATGGCCGTTCTTAGTTG',    'CAAACTTCCATCGGCTTGAGC'),
    ('Microascaceae_askomikot',          'ATCAATAAGCGGAGGAAAAGAAACC','CCTCTTCAAATTACAACTCGGACTG'),
    ('Sakarolitik_Sphaerochaeta',        'ATCTGGCCATGTACTGACGC',     'CTGGTGCACATCGTTTACTGTG'),
    ('Asetoklastik_metanojenler',        'CCGGGAGAGGTGAGAGGTAC',     'CGGGTATCTAATCCGGTTCGTG'),
    ('Bakteri_universal',                'ACAAGCGGTGGAGCATGTG',      'ACGACAGCCATGCAGCAC'),
    ('Petriella_musispora',              'GGAGTCGTCCTAATATGCGAGTG',  'CAAATCCATCCGAGAACATCAGG'),
    ('Proteolitik_Cloacimonas',          'TTAAAGGCAGCGGCTCACC',      'GAACCCGACACCTAGTGATTATCG'),
    ('Arke_universal',                   'CTGCGGTTTAATTGGATTCAACGC', 'GAACTGACGACGGCCATGC'),
    ('Methanothrix_soehngenii_turu',     'AATGTAGCAATACATGGCGAACTG', 'TTCCAGCAATCGAGACCTATCG'),
    ('Methanosarcina_mazei_turu',        'GCCCTTGGGACCGGCATAA',      'TCGCTGGCTAGTAGGTACATTACA'),
]


# --- an exact copy of the OLD (faulty) engine, for the comparison only -------
def eski_yerler(seq, primer, max_mm=1, uc5=False, SEED=13):
    """engine/reads.py -> Sonda.bul davranisinin birebir kopyasi."""
    sd = primer[:SEED] if uc5 else primer[-SEED:]
    off = 0 if uc5 else len(primer) - SEED
    tohumlar = [''.join(x) for x in itertools.product(
        *[om.IUPAC.get(c, 'ACGT') for c in sd])]
    L = len(primer)
    out = []
    for t in tohumlar:
        i = seq.find(t)
        while i != -1:
            st = i - off
            if 0 <= st and st + L <= len(seq):
                mm = 0
                iyi = True
                for a, b in zip(primer, seq[st:st + L]):
                    if b not in om.IUPAC.get(a, 'ACGT'):
                        mm += 1
                        if mm > max_mm:
                            iyi = False
                            break
                if iyi:
                    out.append((st, mm))
            i = seq.find(t, i + 1)
    return sorted(set(out))


def t1_sentetik(max_mm=1):
    """For each primer, synthetic reads are produced with a single mismatch placed at
    each position of the primer. The corrected engine has to match brute force exactly.

    """
    random.seed(11)
    fark = 0; toplam = 0; eski_kacan = 0
    for ad, F, R in CIFTLER:
        for primer, uc5 in ((F, False), (kk.rc(R), True)):
            for poz in range(len(primer)):
                for yeni in 'ACGT':
                    if yeni == primer[poz]:
                        continue
                    varyant = primer[:poz] + yeni + primer[poz + 1:]
                    dolgu1 = ''.join(random.choice('ACGT') for _ in range(120))
                    dolgu2 = ''.join(random.choice('ACGT') for _ in range(120))
                    seq = dolgu1 + varyant + dolgu2
                    a = om.Sonda(primer, uc5, max_mm, son2=True).bul(seq)
                    b = kk.yerler(seq, primer, max_mm, son2=True, uc5=uc5)
                    toplam += len(b)
                    if sorted(a) != sorted(b):
                        fark += 1
                        if fark <= 5:
                            print('  T1 FARK %s poz %d %s: yeni=%s kaba=%s'
                                  % (ad, poz, yeni, a, b))
                    e = eski_yerler(seq, primer, max_mm, uc5)
                    eski_kacan += len(b) - len([x for x in e if x in b])
    print('T1 sentetik : %d baglanma yeri, duzeltilmis motor farki = %d'
          % (toplam, fark))
    print(u'   in the same test the sites the OLD engine missed = %d (%.1f%%)'
          % (eski_kacan, 100.0 * eski_kacan / max(toplam, 1)))
    return fark == 0


def t2_t3_gercek(desen, n=300, max_mm=1):
    yollar = []
    for d in desen:
        g = sorted(glob.glob(d))
        yollar.extend(g if g else [d])
    yollar = [p for p in yollar if os.path.exists(p)]
    if not yollar:
        print(u'T2/T3 skipped (no fastq found)')
        return True
    fark_yer = 0; toplam_yer = 0; fark_urun = 0; eski_kacan = 0
    for p in yollar:
        rs = list(itertools.islice(om.okumalar(p), n))
        if not rs:
            continue
        for ad, F, R in CIFTLER:
            for primer, uc5 in ((F, False), (kk.rc(R), True)):
                s = om.Sonda(primer, uc5, max_mm, son2=True)
                for x in rs:
                    for seq in (x, om.rc(x)):
                        a = s.bul(seq)
                        b = kk.yerler(seq, primer, max_mm, son2=True, uc5=uc5)
                        toplam_yer += len(b)
                        if sorted(a) != sorted(b):
                            fark_yer += 1
                        e = eski_yerler(seq, primer, max_mm, uc5)
                        eski_kacan += len(b) - len([y for y in e if y in b])
            p1 = om.kutu_pcr(rs, F, R, max_mm=max_mm)[0]
            p2 = kk.kutu_pcr(rs, F, R, max_mm=max_mm)[0]
            if p1 != p2:
                fark_urun += 1
                print('  T3 FARK %s %s: yeni=%d kaba=%d' % (os.path.basename(p), ad, p1, p2))
    print(u'T2 real reads   : %d binding sites, difference from the corrected engine = %d'
          % (toplam_yer, fark_yer))
    print(u'   in the same test the sites the OLD engine missed = %d (%.1f%%)'
          % (eski_kacan, 100.0 * eski_kacan / max(toplam_yer, 1)))
    print(u'T3 product count : mismatching (file x pair) = %d' % fark_urun)
    return fark_yer == 0 and fark_urun == 0


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('fastq', nargs='*')
    ap.add_argument('--n', type=int, default=300)
    ap.add_argument('--mm', type=int, default=1)
    a = ap.parse_args()
    print(u'read engine version:', om.__version__, ' criterion mm<=%d' % a.mm)
    ok1 = t1_sentetik(a.mm)
    ok2 = t2_t3_gercek(a.fastq, a.n, a.mm) if a.fastq else True
    print()
    if ok1 and ok2:
        print(u'RESULT: PASSED - the corrected engine is IDENTICAL to brute force.')
        sys.exit(0)
    print(u'RESULT: FAILED - there is a difference, so do NOT write to the panel.')
    sys.exit(1)
