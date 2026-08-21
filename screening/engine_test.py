# -*- coding: utf-8 -*-
"""
engine_test.py - read_engine.py'nin KAYIPSIZ oldugunu kanitlar.

Yontem: `okuma_motoru.Sonda` (guvercin yuvasi tohumlamasi + tam dogrulama) ile
`kaba_kuvvet.yerler` (tohumsuz, her pozisyon tek tek) AYNI dizilerde, AYNI olcutle
kosulur ve baglanma yerlerinin listesi BIREBIR karsilastirilir. Tek bir yer
farki bile testi dusurur.

Uc katman:
  T1  sentetik diziler - her primer pozisyonuna tek tek uyumsuzluk yerlestirilir
      (tohumun icine dusen uyumsuzluk tam da eski motorun kacirdigi durumdur)
  T2  gercek okumalar - panelin 21 cifti, verilen fastq dosyalarindan alt kume
  T3  urun duzeyi - kutu_pcr sayilari iki uygulamada esit mi
Ayrica ESKI motorun ayni testte KAC yer kacirdigi da raporlanir (hatanin kaniti).

Kullanim:
    python engine_test.py                          # yalniz T1 (sentetik, veri gerekmez)
    python engine_test.py "fastq files/A1-4/*.fastq" [--n 300] [--mm 1]
Cikis kodu 0 = gecti, 1 = kaldi.
"""
# ---------------------------------------------------------------------------
# engine_test.py — read_engine.py'nin guvercin yuvasi tohumlamasinin KAYIPSIZ
#                 oldugunu, tohumsuz kaba kuvvetle birebir karsilastirarak
#                 kanitlar.
#
# GIRDI  : T1 icin veri gerekmez (sentetik diziler uretilir); T2 ve T3 icin
#          komut satirinda verilen fastq dosyalari ya da joker desen. Modul
#          olarak read_engine.py ve brute_force.py'yi dogrudan ice aktarir;
#          panelin 21 cifti dosyanin icinde sabit liste olarak durur.
# CIKTI  : dosyaya yazmaz; sonucu ekrana basar. Cikis kodu 0 = gecti,
#          1 = kaldi (tek bir baglanma yeri farki bile testi dusurur).
# CAGRAN : MENUDE DEGILDIR - elle calistirilan bir sinamadir. Ayni
#          karsilastirmanin kucultulmus hali her kosuda self_test.py
#          icinden otomatik yapilir (tus 8 ve butun olcum tuslarinin basi).
#
# GUVERCIN YUVASI ILKESI BURADA SINANIR: k uzunlugundaki bir primer, en fazla
# m uyumsuzluk aranan bir taramada (m + 1) ortusmeyen parcaya bolunurse en az
# bir parca tam eslesmek zorundadir - bu bir garantidir. T1 tam da bu garantiyi
# hedefler: uyumsuzlugu primerin HER pozisyonuna tek tek yerlestirir, cunku
# eski motorun kacirdigi durum uyumsuzlugun tohumun icine dusmesiydi.
# ---------------------------------------------------------------------------
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


# --- ESKI (hatali) motorun birebir kopyasi - yalniz karsilastirma icin -------
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
    """Her primer icin, primerin her pozisyonuna tek uyumsuzluk konulmus sentetik
    okumalar uretilir. Duzeltilmis motor kaba kuvvetle birebir tutmali."""
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
    print('   ayni testte ESKI motorun kacirdigi yer = %d (%.1f%%)'
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
    print('T2 gercek okuma : %d baglanma yeri, duzeltilmis motor farki = %d'
          % (toplam_yer, fark_yer))
    print('   ayni testte ESKI motorun kacirdigi yer = %d (%.1f%%)'
          % (eski_kacan, 100.0 * eski_kacan / max(toplam_yer, 1)))
    print(u'T3 product count : mismatching (file x pair) = %d' % fark_urun)
    return fark_yer == 0 and fark_urun == 0


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('fastq', nargs='*')
    ap.add_argument('--n', type=int, default=300)
    ap.add_argument('--mm', type=int, default=1)
    a = ap.parse_args()
    print(u'read engine version:', om.__version__, ' olcut mm<=%d' % a.mm)
    ok1 = t1_sentetik(a.mm)
    ok2 = t2_t3_gercek(a.fastq, a.n, a.mm) if a.fastq else True
    print()
    if ok1 and ok2:
        print('SONUC: GECTI - duzeltilmis motor kaba kuvvetle BIREBIR ayni.')
        sys.exit(0)
    print(u'RESULT: FAILED - there is a difference, so do NOT write to the panel.')
    sys.exit(1)
