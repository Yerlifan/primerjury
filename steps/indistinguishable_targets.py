#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
indistinguishable_targets.py
Aynı amplikon sınıfı içinde, FARKLI taksonlara atanmış ama dizi düzeyinde
birbirinden ayırt edilemeyen konsensüsleri bulur.

Neden gerekli: Kraken2 tek bir popülasyonun okumalarını kardeş tür
düğümlerine dağıtabiliyor. Bu durumda iki ayrı "takson" kutusu aynı diziyi
üretir. Böyle bir kutu rakip listesine girerse "rakipte ürün oluşmasın"
kuralı mantıken sağlanamaz hale gelir ve hedef sessizce sıfır aday verir.
Bu betik o çiftleri açıkça listeler; 08 aynı ölçümü kullanarak rakip
kümesini temizler ve her çıkarmayı log'lar.

Ölçüm iki bağımsız yoldan yapılır:
  1. minimap2 (mappy) hizalaması, iki yönde ayrı ayrı
  2. hizalamadan bağımsız kanonik k-mer kapsaması
İki ölçüm ayrışırsa çift "ayrisan_olcum" olarak işaretlenir ve karar
temkinli tarafta verilir (ayırt edilemez sayılır).

Kullanım:
  python3 indistinguishable_targets.py --kons <klasor> [--out ayirt_edilemez.tsv]
  python3 indistinguishable_targets.py --kons <klasor> --adlar taxid_adlari.tsv
"""
import argparse, csv, glob, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import alignment

MAPPY = hizalama.ARKA_UC is not None

TAMLAYICI = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")

# IUPAC kodunun temsil ettigi baz kumesi. Hizalayici bu kodlari sessizce
# eslesme sayabilir; o zaman iki konsensus "ayni" gorunur ama aslinda biri
# W (A ya da T), oteki T'dir. Bu yuzden ozdeslik ASAGIDA kendimiz sayilir.
IUPAC_KUME = {"A": "A", "C": "C", "G": "G", "T": "T",
              "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT",
              "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
              "N": "ACGT"}


def baz_kesisir(x, y):
    """Iki konumun alel kumeleri kesisiyor mu.

    N EŞLEŞME SAYILMAZ. N "her baz olabilir" degil, "veri yok" demektir;
    her bazla kesistirilirse N orani yuksek dosyalarda kesisimli ozdeslik
    yapay olarak yukselir ve iki takson gercekte ayrilabilirken ayirt
    edilemez sayilir. Ayni kural baglanma tarafinda da gecerlidir
    (04'teki base_match kalipta N'i uyumsuz sayar)."""
    if x == "N" or y == "N":
        return False
    return bool(set(IUPAC_KUME.get(x, "")) & set(IUPAC_KUME.get(y, "")))

# Eşikler. Veriden değil kuraldan gelir ve burada tek yerde tanımlıdır:
#   OZDESLIK_ESIK : bu özdeşliğin üstünde iki kutu ayırt edilemez sayılır
#   UZUNLUK_ESIK  : bu kadar bazdan kısa hizalama karar için yetersizdir
#   KAPSAMA_ESIK  : mappy yokken kullanılan kaba ölçüt
OZDESLIK_ESIK = 99.5
UZUNLUK_ESIK = 500
KAPSAM_ESIK = 0.60      # hizalama, kisa dizinin en az bu kadarini ortmeli
KAPSAMA_ESIK = 0.97


def oku(f):
    return "".join(l.strip() for l in open(f, encoding="utf-8", errors="replace")
                   if not l.startswith(">")).upper()


def rc(s):
    return s.translate(TAMLAYICI)[::-1]


def kanonik_kmer(s, k=21):
    out = set()
    for i in range(len(s) - k + 1):
        m = s[i:i + k]
        if all(c in "ACGT" for c in m):
            out.add(min(m, rc(m)))
    return out


def hizala(s1, s2):
    """(hizalanan_uzunluk, kesisimli_ozdeslik, kati_ozdeslik) doner.

    kesisimli_ozdeslik : alel kumeleri kesisen konumlar eslesme sayilir.
        Primer tasarimini ilgilendiren olcut budur; kumeler kesisiyorsa o
        konumda hicbir primer bazi iki kalibi ayiramaz.
    kati_ozdeslik      : karakterler birebir ayni mi. W ile T farkli sayilir.

    Hizalayicinin kendi mlen degeri kullanilmaz, cunku IUPAC kodlarini
    sessizce eslesme sayabilir ve iki konsensus oldugundan daha benzer
    gorunur."""
    if not MAPPY:
        return (0, 0.0, 0.0)
    try:
        A = hizalama.Hizalayici(seq=s1, preset="asm20")
    except RuntimeError:
        return (0, 0.0, 0.0)
    if not A:
        return (0, 0.0, 0.0)
    best = None
    for h in A.map(s2):
        if best is None or h.blen > best.blen:
            best = h
    if best is None:
        return (0, 0.0, 0.0)
    # Ters zincir hizalamasinda mappy'nin q_st/q_en degerleri OZGUN sorgu
    # koordinatlarindadir, CIGAR ise ters tumleyen uzerinde yurur. Baslangici
    # duzeltmezsek karsilastirma kayar ve ozdeslik anlamsiz cikar.
    if best.strand == 1:
        q, qp = s2, best.q_st
    else:
        q, qp = rc(s2), len(s2) - best.q_en
    rp = best.r_st
    kesisen = kati = toplam = 0
    for ln, op in best.cigar:
        if op == 0:
            for k in range(ln):
                if rp + k < len(s1) and qp + k < len(q):
                    x, y = s1[rp + k], q[qp + k]
                    toplam += 1
                    if x == y:
                        kati += 1
                    if baz_kesisir(x, y):
                        kesisen += 1
            rp += ln; qp += ln
        elif op == 1:
            qp += ln
        elif op == 2:
            rp += ln
        elif op == 4:
            qp += ln
    if not toplam:
        return (best.blen, 0.0, 0.0)
    return (best.blen, 100.0 * kesisen / toplam, 100.0 * kati / toplam)


def sinif_of(tag):
    return re.split(r"[-_]", tag)[0]


def envanter(kons, min_kapsanan=200):
    """Sinif ve takson basina TEK temsilci: en cok baz kapsayan konsensus."""
    temsil = {}
    for f in sorted(glob.glob(os.path.join(kons, "*_konsensus.fasta"))):
        b = os.path.basename(f)
        m = re.match(r"((?:A1|A2|B|F1|F2))-(\d+)_(\d+)_", b)
        if not m:
            continue
        s = oku(f)
        kapsanan = len(s.replace("N", ""))
        if kapsanan < min_kapsanan:
            continue
        key = (m.group(1), m.group(3))
        if key not in temsil or kapsanan > temsil[key][1]:
            temsil[key] = (b, kapsanan, s)
    return temsil


def ayirt_edilemezler(temsil, ozdeslik_esik=OZDESLIK_ESIK,
                      uzunluk_esik=UZUNLUK_ESIK, kapsama_esik=KAPSAMA_ESIK):
    """[(sinif, taxid1, taxid2, uzunluk, ozdeslik, kapsama, gerekce)] doner."""
    sonuc = []
    siniflar = sorted(set(s for s, _ in temsil))
    for sn in siniflar:
        txs = sorted([t for (s, t) in temsil if s == sn], key=int)
        kmer = {t: kanonik_kmer(temsil[(sn, t)][2]) for t in txs}
        for i, t1 in enumerate(txs):
            s1 = temsil[(sn, t1)][2]
            for t2 in txs[i + 1:]:
                s2 = temsil[(sn, t2)][2]
                # 1. olcum: iki yonde hizalama
                # Uzunluk, kapsam ve ozdeslik AYNI hizalamadan alinir;
                # farkli yonlerden karistirilirsa olcutler birbirini
                # tutmayan hizalamalara ait olur.
                u12, o12, k12 = hizala(s1, s2)
                u21, o21, k21 = hizala(s2, s1)
                if u12 >= u21:
                    u, o, kati = u12, o12, k12
                else:
                    u, o, kati = u21, o21, k21
                # Yerel hizalama yalnizca korunmus omurgayi ortebilir ve
                # ozdeslik boylece yukari sapar. Bu yuzden hizalamanin KISA
                # dizinin ne kadarini ortugu de kosula alinir.
                # Pay ve payda ayni birimde olmali: hizalamanin ortugu
                # KOLON sayisi (toplam) ile kisa dizinin ham uzunlugu.
                # Eski surumde pay hizalama blok uzunlugu, payda N'siz
                # uzunluktu; N orani yuksek dosyalarda kapsam 1'i asabiliyordu.
                kisa = min(len(s1), len(s2))
                kapsam = min(1.0, u / max(1, kisa))
                hiz_der = (u >= uzunluk_esik and kapsam >= KAPSAM_ESIK and
                           o >= ozdeslik_esik)
                # 2. olcum: hizalamadan bagimsiz k-mer kapsamasi
                K1, K2 = kmer[t1], kmer[t2]
                kap = len(K1 & K2) / max(1, min(len(K1), len(K2)))
                kap_der = kap >= kapsama_esik
                if MAPPY:
                    if hiz_der and kap_der:
                        ger = "iki_olcum_uyustu"
                    elif hiz_der or kap_der:
                        ger = "ayrisan_olcum"      # temkinli taraf
                    else:
                        continue
                else:
                    if not kap_der:
                        continue
                    ger = "yalniz_kmer_olcumu"
                sonuc.append((sn, t1, t2, u, round(o, 2),
                              round(kap, 4), ger, round(kapsam, 3),
                              round(kati, 2)))
    return sonuc


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--kons", required=True)
    p.add_argument("--adlar", default=None,
                   help="taxid<TAB>ad file, only raporu okunakli kilar")
    p.add_argument("--out", default=None)
    p.add_argument("--ozdeslik-esik", type=float, default=OZDESLIK_ESIK)
    p.add_argument("--uzunluk-esik", type=int, default=UZUNLUK_ESIK)
    p.add_argument("--kapsama-esik", type=float, default=KAPSAMA_ESIK)
    return p.parse_args()


def main():
    a = get_args()
    print(hizalama.durum())
    if not MAPPY:
        print("UYARI: hizalama arka ucu yok. Yalniz k-mer kapsamasi "
              "kullanilacak, olcum kabalasir.")
    ad = {}
    if a.adlar and os.path.exists(a.adlar):
        for line in open(a.adlar, encoding="utf-8"):
            p = line.rstrip("\n").split("\t")
            if len(p) > 1:
                ad[p[0]] = p[1]
    temsil = envanter(a.kons)
    print("temsilci konsensus: %d (sinif, takson) kutusu" % len(temsil))
    rows = ayirt_edilemezler(temsil, a.ozdeslik_esik, a.uzunluk_esik,
                             a.kapsama_esik)
    print("\nAYIRT EDILEMEZ TAKSON CIFTLERI "
          "(ozdeslik >= %%%.1f ve hizalanan >= %d bp, ya da k-mer kapsamasi >= %.2f)"
          % (a.ozdeslik_esik, a.uzunluk_esik, a.kapsama_esik))
    print("%-3s %-28s %-28s %8s %8s %8s %6s %8s  %s"
          % ("sn", "takson 1", "takson 2", "hizalanan", "kesisimli", "kati",
             "kapsam", "kmer", "gerekce"))
    for sn, t1, t2, u, o, k, g, kp, kt in sorted(rows, key=lambda r: (r[0], -r[4])):
        print("%-3s %-28s %-28s %8d  %%%.2f  %%%.2f %6.2f %8.4f  %s"
              % (sn, ad.get(t1, t1)[:27], ad.get(t2, t2)[:27], u, o, kt, kp, k, g))
    print("\ntoplam: %d cift" % len(rows))
    if a.out:
        os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
        with open(a.out, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["sinif", "taxid1", "ad1", "taxid2", "ad2",
                        "hizalanan_bp", "kesisimli_ozdeslik", "kati_ozdeslik",
                        "hizalama_kapsami", "kmer_kapsamasi", "gerekce"])
            for sn, t1, t2, u, o, k, g, kp, kt in rows:
                w.writerow([sn, t1, ad.get(t1, ""), t2, ad.get(t2, ""),
                            u, o, kt, kp, k, g])
        print("yazildi: %s" % a.out)


if __name__ == "__main__":
    main()
