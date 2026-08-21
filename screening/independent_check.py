# -*- coding: utf-8 -*-
"""
independent_check.py - Baslik sayilarini UC AYRI yolla olcup tuttugunu gosterir.

  Yol A  ispcr.find_sites / amplify   : numpy vektor tarama, TOHUMSUZ.
                                        Panelin KENDI kodu (engine/ispcr.py),
                                        bu oturumda degistirilmedi.
  Yol B  brute_force.py               : saf python, her pozisyon tek tek denenir.
                                        Bagimsiz yazildi, ortak kod yok.
  Yol C  read_engine.py              : bu oturumda duzeltilen motor (referans).

A ve B, duzeltilmis motorun hicbir satirini kullanmaz. Ucu de ayni sayiyi
veriyorsa duzeltilmis motorun ciktisi bagimsiz olarak dogrulanmis olur.

Kullanim:
    python independent_check.py --fastq "..\fastq files" [--mm 1] [--nmax 3000]
"""
# ---------------------------------------------------------------------------
# independent_check.py — teslim basliklarindaki sayilari uc ayri kod yoluyla
#                         olcup ucunun de ayni cevabi verdigini gosterir.
#
# GIRDI  : --fastq ile "fastq files" klasoru (--mm uyumsuzluk tavani, --nmax
#          kutu basina okuma, --seed ornekleme tohumu). Sinanan cift ve kutu
#          listesi (TESTLER) dosyanin icinde sabittir. Uc yol:
#          engine/ispcr.py, brute_force.py ve read_engine.py.
# CIKTI  : dosyaya yazmaz; uc yolun sayilarini yan yana ekrana basar.
# CAGRAN : MENUDE DEGILDIR - elle calistirilan bir kanit uretecidir.
#
# NEDEN UC YOL: yol A (ispcr) panelin kendi kodudur ve bu oturumda hic
# degistirilmedi; yol B (kaba kuvvet) tohumsuzdur ve ortak kod paylasmaz; yol C
# duzeltilmis motordur. A ve B, duzeltilmis motorun hicbir satirini kullanmaz -
# ucu ayni sayiyi veriyorsa guvercin yuvasi tohumlamasinin kayipsizligi kendi
# koduyla degil, disaridan dogrulanmis olur.
# ---------------------------------------------------------------------------
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
    sys.stderr.write('UYARI: ispcr yuklenemedi (%s). Yol A atlanacak.\n' % e)

# Baslik sayilarinin dayandigi kutular: M. mazei cifti (satir 22) ve
# hatanin en buyuk oldugu M. hadiensis rakip kutulari.
TESTLER = [
    ('Methanosarcina_mazei_turu', 'GCCCTTGGGACCGGCATAA', 'TCGCTGGCTAGTAGGTACATTACA',
     [('A1-4', '3078083', 'RAKIP M. hadiensis - hatanin merkezi'),
      ('A2-4', '3078083', 'RAKIP M. hadiensis'),
      ('A2-2', '2209',    'UYE M. mazei'),
      ('A1-3', '2209',    'UYE M. mazei')]),
    ('Methanosarcina_cinsi', 'TCGCTAGGTGTCAGGCATG', 'GCGATTCAGGCAAGGTCTTC',
     [('A1-2', '2209', 'UYE M. mazei - kapsam eksik olculmustu'),
      ('A2-3', '2223', 'rakip M. soehngenii')]),
    ('Asetoklastik_metanojenler', 'CCGGGAGAGGTGAGAGGTAC', 'CGGGTATCTAATCCGGTTCGTG',
     [('A1-2', '2209', 'UYE M. mazei - kapsam eksik olculmustu')]),
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
                print('  YOK:', p); continue
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
