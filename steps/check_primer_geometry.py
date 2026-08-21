#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_primer_geometry.py
İleri ve geri primerin doğru tasarlanıp tasarlanmadığını, tasarım kodunun
HİÇBİR fonksiyonunu kullanmadan, sıfırdan denetler.

Neden ayrı bir betik: tasarım kodu kendi ürettiği koordinatlarla kendini
doğrularsa, koordinat hatası ikisinde de aynı yönde olur ve görünmez. Bu
betik yalnızca kalıp dizisini, primer dizilerini ve bildirilen konumları
alır; ürünü kalıptan keser ve şu dört şartı bağımsız olarak sınar:

  1. Ürünün BAŞI ileri primerin kendisidir.
  2. Ürünün SONU geri primerin TERS TÜMLEYENİDİR.
  3. Geri primer, kalıbın eksi zincirinde okunduğunda 3' ucu ürünün içine
     bakar; yani artı zincirdeki karşılığının 5' ucuna denk gelir.
  4. Bildirilen ürün uzunluğu, kesilen ürünün uzunluğuna eşittir.

Ayrıca her primerin 3' ucunun kalıpta belirsiz (IUPAC ya da N) bir
pozisyona denk gelip gelmediği denetlenir; geri primerin 3' ucu pencerenin
BAŞINA düştüğü için burası kolayca gözden kaçar.

Kullanım:
  python3 check_primer_geometry.py --tsv primer_adaylari/X__A1.tsv \
      --kons referans_konsensus/baskin/konsensus --capa A1-1_2209
"""
import argparse, csv, glob, os, re, sys

TAM = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn",
                    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T",
         "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def uyar(oligo, kalip_parca):
    """Oligo, kalip parcasiyla bagdasiyor mu.

    Kalipta IUPAC kodu bulunabilir; primer ise salt ACGT'dir ve kodun
    temsil ettigi alellerden birine cozulmustur. Bu yuzden karakter
    esitligi degil, KUME UYUMU aranir. Kalipta N varsa o konum bilinmiyor
    demektir ve uyum sayilmaz."""
    if len(oligo) != len(kalip_parca):
        return False
    for o, k in zip(oligo, kalip_parca):
        if k == "N":
            return False
        if o not in IUPAC.get(k, ""):
            return False
    return True


def rc(s):
    """Ters tümleyen. Önce tümlenir, sonra ters çevrilir."""
    return s.translate(TAM)[::-1]


def oku(f):
    return "".join(l.strip() for l in open(f, encoding="utf-8",
                                           errors="replace")
                   if not l.startswith(">")).upper()


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--tsv", required=True)
    p.add_argument("--kons", required=True)
    p.add_argument("--capa", default=None,
                   help="capa etiketi; verilmezse TSV'deki capa sutunu kullanilir")
    p.add_argument("--en-fazla", type=int, default=200)
    return p.parse_args()


def main():
    a = get_args()
    rows = list(csv.DictReader(open(a.tsv, encoding="utf-8"), delimiter="\t"))
    if not rows:
        sys.exit("bos TSV")
    capa_ad = a.capa or rows[0].get("capa", "")
    aday = sorted(glob.glob(os.path.join(a.kons, "*%s*" % capa_ad))) or \
        sorted(glob.glob(os.path.join(a.kons, "*")))
    if not aday:
        sys.exit("capa konsensusu bulunamadi: %s" % capa_ad)
    kalip = oku(aday[0])
    print("capa dosyasi : %s" % os.path.basename(aday[0]))
    print("kalip uzunluk: %d" % len(kalip))
    print("sinanan satir: %d\n" % min(len(rows), a.en_fazla))

    say = dict(tamam=0, urun_basi=0, urun_sonu=0, yon=0, uzunluk=0,
               f_3p_belirsiz=0, r_3p_belirsiz=0, kalip_disi=0)
    ornek = []
    for r in rows[:a.en_fazla]:
        F, R = r["ileri_dizi"], r["geri_dizi"]
        try:
            fb = int(r["ileri_baslangic"]) - 1
            fl = int(r["ileri_uzunluk"])
            gb = int(r["geri_baslangic"]) - 1
            gl = int(r["geri_uzunluk"])
            bildirilen_min = int(r["urun_min"])
            bildirilen_maks = int(r.get("urun_maks") or r["urun_min"])
        except (KeyError, ValueError):
            sys.exit("TSV'de konum sutunlari yok; 04'un guncel surumu gerekli")
        if fb < 0 or gb + gl > len(kalip):
            say["kalip_disi"] += 1
            continue

        # Urun: ileri primerin 5' ucundan, geri primerin kalip uzerindeki
        # penceresinin sonuna kadar. Geri primer eksi zincirde okundugu icin
        # kalip uzerindeki penceresi [gb, gb+gl) araligidir.
        urun = kalip[fb:gb + gl]
        hata = []

        # 1. urunun basi = ileri primer
        if not uyar(F, urun[:fl]):
            say["urun_basi"] += 1
            hata.append("urun_basi")
        # 2. urunun sonu = geri primerin ters tumleyeni
        if not uyar(rc(R), urun[-gl:]):
            say["urun_sonu"] += 1
            hata.append("urun_sonu")
        # 3. yon: geri primer kalip penceresinin ters tumleyeni olmali.
        #    Kalip kodlarinin ters tumleyeni alinip oligo ona karsi sinanir.
        if not uyar(R, rc(kalip[gb:gb + gl])):
            say["yon"] += 1
            hata.append("yon")
        # 4. uzunluk: bildirilen aralik uyelere gore hesaplanir, capadaki
        #    urun bu araligin icinde olmali
        if not (bildirilen_min <= len(urun) <= bildirilen_maks):
            say["uzunluk"] += 1
            hata.append("uzunluk(%d, bildirilen %d-%d)"
                        % (len(urun), bildirilen_min, bildirilen_maks))
        # 5. 3' uclarin kalipta belirsiz pozisyona denk gelmesi
        #    ileri primerin 3' ucu  -> kalipta fb+fl-1
        #    geri primerin 3' ucu   -> kalipta gb (pencerenin BASI)
        if kalip[fb + fl - 1] not in "ACGT":
            say["f_3p_belirsiz"] += 1
            hata.append("F_3p_belirsiz(%s)" % kalip[fb + fl - 1])
        if kalip[gb] not in "ACGT":
            say["r_3p_belirsiz"] += 1
            hata.append("R_3p_belirsiz(%s)" % kalip[gb])

        if not hata:
            say["tamam"] += 1
        elif len(ornek) < 5:
            ornek.append((r, urun, hata))

    n = min(len(rows), a.en_fazla)
    print("SONUC")
    print("   dort geometri sartini da gecen satir : %d / %d" % (say["tamam"], n))
    print("   urunun basi ileri primer degil       : %d" % say["urun_basi"])
    print("   urunun sonu geri primerin rc'si degil: %d" % say["urun_sonu"])
    print("   geri primer kalibin rc'si degil      : %d" % say["yon"])
    print("   bildirilen uzunluk tutmuyor          : %d" % say["uzunluk"])
    print("   ileri primerin 3' ucu belirsiz bazda : %d" % say["f_3p_belirsiz"])
    print("   geri primerin 3' ucu belirsiz bazda  : %d" % say["r_3p_belirsiz"])
    print("   kalip disina tasan satir             : %d" % say["kalip_disi"])
    for r, urun, hata in ornek:
        print("\n   HATA %s" % ", ".join(hata))
        print("      F=%s  R=%s" % (r["ileri_dizi"], r["geri_dizi"]))
        print("      urun bas=%s ... son=%s" % (urun[:28], urun[-28:]))
        print("      rc(R)   =%s" % rc(r["geri_dizi"]))
    if say["tamam"] == n:
        print("\nButun satirlar dort geometri sartini da gecti.")


if __name__ == "__main__":
    main()
