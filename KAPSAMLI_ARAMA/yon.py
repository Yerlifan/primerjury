# -*- coding: utf-8 -*-
"""
yon.py - KANONIK YON TANIMI VE NORMALIZASYONU. Tek kaynak.

PROJENIN KANONIK YONU: SENSE (referans / artı iplik).
  SSU rRNA ve ITS konsensusleri, referans veritabanlarindaki (SILVA/RefSeq) yonde
  saklanir. Sebep: projenin butun in-silico PCR motorlari (ispcr.amplify,
  okuma_motoru.Sonda) verilen diziyi YALNIZ ARTI IPLIKTE tarar. Ters saklanmis bir
  konsensuste ileri primer de, ters primerin tumleyeni de bulunamaz; motor hata
  atmaz, "urun yok" der.
  OLCULEN ETKI: yon_etki_testi.py -> dogru yonde 117 urun, ters yonde 0. Kayip %100.

NEDEN BU MODUL VAR: yon hatasi gece boyunca EN AZ UC AYRI YERDE ayri ayri bulunup
ayri ayri yamandi. Uc ayri yama = tek kanonik cozum yok = bir sonraki degisiklikte
yine kacar. Bundan sonra yon TEK BIR YERDE, burada tanimlanir; her betik
`konsensus_kanonik/` klasorunu okur ve yon sorgusu icin bu modulu cagirir.

IKI BAGIMSIZ OLCUT (proje kurali: hicbir karar tek kod yoluna birakilmaz).
Ikisi ayrilirsa dosya BELIRSIZ sayilir ve normalize EDILMEZ, isaretlenir.

  Olcut 1  panelin kendi evrensel ciftleri: sense yonde F ve rc(R) dogrudan bulunur
  Olcut 2  literatur evrensel motifleri (panelden bagimsiz):
             SSU  515F, 806R-sense, 1100-sense
             ITS  ITS1, rc(ITS4)

API:
    yon.tespit(dizi, sinif)        -> ('SENSE'|'ANTISENSE'|'BELIRSIZ', ayrinti_dict)
    yon.kanonik(dizi, sinif)       -> (kanonik_dizi, karar, cevrildi_mi)
    yon.sinifi(yol_veya_ad)        -> 'A'|'B'|'F1'|'F2'|'?'
"""
# ---------------------------------------------------------------------------
# yon.py — kanonik yonun (SENSE) TEK tanim yeri: bir dizinin hangi yonde
#          saklandigini olcer ve gerekirse ters tumleyenini alarak kanonige
#          cevirir.
#
# GIRDI  : dogrudan dizi ve amplikon sinifi alir; dosya_kanonik() ise bir fasta
#          dosyasi okur. Baglanma yeri aramasi icin okuma_motoru.Sonda
#          kullanilir (kayipsiz, <=2 uyumsuzluk toleransli).
# CIKTI  : dosyaya yazmaz. tespit() (karar, ayrinti) ikilisi; kanonik()
#          (kanonik_dizi, karar, cevrildi_mi) uclusu; dosya_kanonik() kayit
#          listesi; kendini_sina() bos liste (gecti) ya da hata metinleri
#          dondurur. Dogrudan calistirilirsa sinav sonucunu ekrana basar.
# CAGRAN : kanonik_uret.py (konsensus_kanonik uretimi), hepsi.yon_kapisi (her
#          asamanin basindaki kapi), kendini_sina.yon_sinamasi ve
#          konsensus_uret.py. Yani KAPSAMLI_ARAMA.bat tuslari 1, 2, 3, 4, 5, 6,
#          7, 8 ve 9'un tamaminda dolayli olarak calisir.
#
# NEDEN AYRI VE TEK BIR MODUL: yon hatasi ayni gece uc ayri yerde ayri ayri
# bulunup ayri ayri yamandi. Uc yama demek, tek bir kanonik cozum olmamasi ve
# bir sonraki degisiklikte hatanin yeniden kacmasi demektir. Karar iki BAGIMSIZ
# olcute baglanmistir (panelin kendi evrensel ciftleri ve literatur motifleri);
# ikisi ayrilirsa dizi BELIRSIZ sayilir ve normalize EDILMEZ, isaretlenir -
# cunku yanlis yone cevrilmis bir dizi, hic cevrilmemis kadar sessiz zarar verir.
# ---------------------------------------------------------------------------
import os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import okuma_motoru as om

KANONIK_YON = 'SENSE'
TOLERANS_MM = 2          # nanopore konsensusunde tek tuk hata olur

PANEL_CIFT = {
    'A':  ('Arke_universal',        'CTGCGGTTTAATTGGATTCAACGC', 'GAACTGACGACGGCCATGC'),
    'B':  ('Bakteri_universal',     'ACAAGCGGTGGAGCATGTG',      'ACGACAGCCATGCAGCAC'),
    'F1': ('Mantar_universal (F1)', 'GGTTACCCGCTGAACTTAAGC',    'CGCTTCACTCGCCGTTAC'),
    'F2': ('Mantar_universal (F2)', 'GTGCATGGCCGTTCTTAGTTG',    'CAAACTTCCATCGGCTTGAGC'),
}
MOTIF = {
    'SSU': ['GTGYCAGCMGCCGCGGTAA', 'GGATTAGATACCC', 'AGTCCCGCAACGAGCGCAACCC'],
    'ITS': ['TCCGTAGGTGAACCTGCGG', 'GCATATCAATAAGCGGAGGA'],
}
rc = om.rc
temizle = om.temizle


def sinifi(ad):
    for s, onek in (('F1', ('F1-', 'F1_')), ('F2', ('F2-', 'F2_')),
                    ('A', ('A1-', 'A2-', 'A1_', 'A2_')), ('B', ('B-', 'B_'))):
        if any(x in ad for x in onek):
            return s
    return '?'


def _var(dizi, desen, mm=TOLERANS_MM):
    return len(om.Sonda(desen, uc5=False, max_mm=mm, son2=False).bul(dizi)) > 0


def _olcut1(dizi, sinif):
    if sinif not in PANEL_CIFT:
        return 'BELIRSIZ', 0, 0
    _, F, R = PANEL_CIFT[sinif]
    d = int(_var(dizi, F)) + int(_var(dizi, rc(R)))
    t = int(_var(dizi, rc(F))) + int(_var(dizi, R))
    return ('SENSE' if d > t else 'ANTISENSE' if t > d else 'BELIRSIZ'), d, t


def _olcut2(dizi, sinif):
    tip = 'ITS' if sinif in ('F1', 'F2') else 'SSU'
    d = sum(1 for m in MOTIF[tip] if _var(dizi, m))
    t = sum(1 for m in MOTIF[tip] if _var(rc(dizi), m))
    return ('SENSE' if d > t else 'ANTISENSE' if t > d else 'BELIRSIZ'), d, t


def tespit(dizi, sinif):
    """Donus: (karar, ayrinti). karar in {'SENSE','ANTISENSE','BELIRSIZ'}"""
    # Iki olcut BIRBIRINDEN BAGIMSIZ secilmistir: olcut 1 panelin kendi
    # evrensel primerlerini kullanir, olcut 2 panelden hic haberi olmayan
    # literatur motiflerini (515F, 806R-sense, ITS1, rc(ITS4)) kullanir. Ayni
    # kaynaktan turemis iki olcut ayni yonde yanilabilirdi.
    # Karar kurali muhafazakardir: ikisi uyusursa karar kesin, biri sessiz
    # kalirsa konusan gecerli, IKISI CELISIRSE karar BELIRSIZ olur ve dizi
    # cevrilmez. Yanlis yone cevirmek, cevirmemekten daha zararlidir cunku
    # sonraki asamalar dosyanin kanonik oldugunu varsayar.
    # 200 bp'den kisa diziler dogrudan BELIRSIZ sayilir: motiflerin hepsinin
    # sigmadigi bir parcada sayimlar rastlantisal olur.
    dizi = temizle(dizi)
    if len(dizi) < 200:
        return 'BELIRSIZ', dict(sebep='dizi 200 bp\'den kisa', uzunluk=len(dizi))
    o1, d1, t1 = _olcut1(dizi, sinif)
    o2, d2, t2 = _olcut2(dizi, sinif)
    ay = dict(olcut1=o1, olcut1_duz=d1, olcut1_ters=t1,
              olcut2=o2, olcut2_duz=d2, olcut2_ters=t2)
    if o1 == o2 and o1 != 'BELIRSIZ':
        ay['sebep'] = 'iki olcut de ayni yonu veriyor'
        return o1, ay
    if o1 != 'BELIRSIZ' and o2 == 'BELIRSIZ':
        ay['sebep'] = 'yalniz olcut 1 karar verdi (olcut 2 sessiz)'
        return o1, ay
    if o2 != 'BELIRSIZ' and o1 == 'BELIRSIZ':
        ay['sebep'] = 'yalniz olcut 2 karar verdi (olcut 1 sessiz)'
        return o2, ay
    if o1 != o2:
        ay['sebep'] = 'IKI OLCUT AYRILDI - normalize edilmez, maskelenir'
        return 'BELIRSIZ', ay
    ay['sebep'] = 'iki olcut de sessiz (motif bulunamadi)'
    return 'BELIRSIZ', ay


def kanonik(dizi, sinif):
    """Donus: (kanonik_dizi, karar, cevrildi_mi).
    BELIRSIZ ise dizi DEGISTIRILMEZ - cagiran taraf isaretlemelidir."""
    dizi = temizle(dizi)
    karar, ay = tespit(dizi, sinif)
    if karar == 'ANTISENSE':
        return rc(dizi), karar, True
    return dizi, karar, False


def dosya_kanonik(yol):
    """Bir fasta dosyasindaki her kaydi kanonige cevirir.
    Donus: [(baslik, kanonik_dizi, karar, cevrildi)]"""
    sn = sinifi(os.path.basename(yol)) or sinifi(yol)
    if sn == '?':
        sn = sinifi(yol)
    out, ad, buf = [], None, []
    with open(yol, errors='ignore') as fh:
        for l in fh:
            if l.startswith('>'):
                if ad is not None:
                    out.append((ad, ''.join(buf)))
                ad, buf = l[1:].rstrip('\n'), []
            else:
                buf.append(l.strip())
    if ad is not None:
        out.append((ad, ''.join(buf)))
    son = []
    for a, s in out:
        k, karar, cev = kanonik(s, sn)
        son.append((a, k, karar, cev))
    return son, sn


def kendini_sina():
    """Modulun kendi sinavi. Ana is baslamadan once kosar (proje kurali 2)."""
    hata = []
    # 1) bilinen sense bir SSU parcasi: 515F + 806R-sense icerir
    s = ('GG' * 60 + 'GTGCCAGCAGCCGCGGTAA' + 'AC' * 120 + 'GGATTAGATACCC' + 'TT' * 60)
    k, karar, cev = kanonik(s, 'A')
    if karar != 'SENSE' or cev:
        hata.append('sentetik sense SSU yanlis: %s' % karar)
    # 2) ayni dizinin tersi ANTISENSE bulunmali ve geri cevrilmeli
    k2, karar2, cev2 = kanonik(rc(s), 'A')
    if karar2 != 'ANTISENSE' or not cev2:
        hata.append('sentetik antisense SSU yanlis: %s' % karar2)
    elif k2 != temizle(s):
        hata.append('cevirme diziyi geri getirmedi')
    # 3) motifsiz dizi BELIRSIZ olmali
    _, karar3, _ = kanonik('A' * 400, 'A')
    if karar3 != 'BELIRSIZ':
        hata.append('motifsiz dizi BELIRSIZ degil: %s' % karar3)
    # 4) idempotans: kanonigi tekrar kanoniklestirmek degistirmemeli
    k4, _, cev4 = kanonik(k2, 'A')
    if cev4 or k4 != k2:
        hata.append('idempotans bozuk')
    # 5) sinif tespiti
    if sinifi('A1-4_2209_konsensus.fasta') != 'A' or sinifi('F2-1_101201.fasta') != 'F2':
        hata.append('sinif tespiti bozuk')
    return hata


if __name__ == '__main__':
    h = kendini_sina()
    if h:
        print('KENDINI SINAMA KALDI:')
        for x in h:
            print('  -', x)
        sys.exit(1)
    print('yon.py kendini sinama: GECTI (5/5)  | kanonik yon =', KANONIK_YON)
