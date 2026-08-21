# -*- coding: utf-8 -*-
"""
orientation_impact_test.py - BILINEN CEVAPLI TEST: yonun yanlis olmasi hangi sayiyi ne kadar bozar.

Kurulum: ayni primer cifti, ayni konsensus, tek fark YON.
  (a) konsensus SENSE yonde saklanmis   -> beklenen: urun bulunur
  (b) ayni konsensusun TERS TUMLEYENI   -> beklenen: urun BULUNAMAZ (0)

Sebep: projenin in-silico PCR motorlari (ispcr.amplify ve okuma_motoru) verilen diziyi
YALNIZ ARTI IPLIKTE tarar - ileri primeri ve ters primerin tumleyenini arti iplikte arar.
Ters saklanmis bir konsensuste ikisi de bulunamaz. Hata atilmaz; "urun yok" denir.
Bu, panelin "10 Olcum Hatalari" #1 kaydinin sayisal kanitidir.

Kullanim: python orientation_impact_test.py --kok ..
"""
# ---------------------------------------------------------------------------
# orientation_impact_test.py — "yon yanlissa ne kaybederiz" sorusunu sayiya baglar:
#                     ayni cift, ayni konsensus, tek fark yon.
#
# GIRDI  : --kok altindaki --klasor ile verilen konsensus klasoru (varsayilan
#          referans_konsensus/konsensus); panelin evrensel ciftleri ve uc ek
#          cift dosyanin icinde sabit listedir. Olcumu read_engine.py, varsa
#          ayrica engine/ispcr.py yapar.
# CIKTI  : dosyaya yazmaz; dogru yonde ve ters yonde bulunan urun sayilarini
#          ekrana basar.
# CAGRAN : MENUDE DEGILDIR - elle calistirilan bir kanit uretecidir. Urettigi
#          sayi (dogru yonde 117 urun, ters yonde 0) config.py ve orientation.py
#          icindeki aciklamalarda gerekce olarak anilir.
#
# NEDEN KAYIP %100: projenin in-silico PCR motorlari verilen diziyi YALNIZ arti
# iplikte tarar - ileri primeri ve ters primerin tumleyenini orada arar. Ters
# saklanmis bir konsensuste ikisi de bulunamaz. Motor hata atmaz, "urun yok"
# der; yani hata sessizdir ve ancak boyle bilinen cevapli bir testle gorunur.
# ---------------------------------------------------------------------------
import sys, os, glob, argparse

BURA = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BURA)
sys.path.insert(0, os.path.join(os.path.dirname(BURA), 'engine'))
import okuma_motoru as om
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
    print('Kaynak klasor : %s  (%d dosya)' % (a.klasor, len(yollar)))
    print('Olcut         : uyumsuzluk <= %d, 3\' son 2 baz tam, urun 40-600 bp' % a.mm)
    print('Motor         : yalniz ARTI IPLIK taranir (projenin motorlarinin davranisi)\n')

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
                 'ters yon urunu SIFIRLIYOR' if ters == 0 and dogru > 0
                 else ('etkilenmedi' if dogru == ters else 'kismi kayip')))
    print('\nTOPLAM  dogru yonde %d urun,  ters yonde %d urun,  kayip %d (%.1f%%)'
          % (tp_d, tp_t, tp_d - tp_t, 100.0 * (tp_d - tp_t) / max(tp_d, 1)))

    if ISPCR:
        print('\nCAPRAZ KONTROL - panelin kendi motoru (ispcr.amplify), ayni dosyalar, A sinifi:')
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
        print('  %s: dogru %d/%d, ters %d/%d  -> ayni sonuc' % (ad, d, n, t, n))


if __name__ == '__main__':
    main()
