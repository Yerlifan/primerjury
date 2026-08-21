#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
abundance_rank.py
BOLLUK TABLOSUNU, VERININ DESTEKLEDIGI RUTBEDE VERIR.

Neden Bracken calistirilmadi
----------------------------
Ilk plan, guven duzeltmesi uygulanmis Kraken2 raporlarini Bracken'a verip
cins duzeyinde bolluk uretmekti. Once olculdu ve plan birakildi.

Bracken, ust rutbelerde duran okumalari veritabanindaki k-mer dagilimini
kullanarak asagi dagitir. Bu, gercek organizmanin veritabaninda BULUNDUGU
varsayimina dayanir. Bu numunede o varsayim kurulamiyor:

  Guven esigi 0,02'den sonra 2 314 887 okumanin
     cins duzeyine inebilen  :  400 190  (%17,3)
     cins duzeyine inemeyen  : 1 723 684  (%77,7)
     siniflandirilmamis      :  115 343  (%5,0)

  Grup grup bakildiginda tablo daha da keskin:
     Arke  (barcode01-08) : okumalarin %72-96'si cins duzeyine iniyor
     Bakteri (barcode17-20): yalnizca %8-9
     Mantar F1 (13-16)    : %0,6-2,4
     Mantar F2 (09-12)    : %0,02-0,3

Yani arke tarafinda cins duzeyi bolluk saglam, bakteri ve mantar tarafinda
degil. Mantarda okumalarin %99,7'sini veritabani onceliklerine gore cinse
dagitmak, sayi uretmek olur. Ustelik bu numunede baskin mantarin
(Microascaceae, Petriella'ya yakin) ve bakteri soylarinin cogunun
veritabaninda bulunmadigi ayrica olculdu.

Bu betigin yaptigi
------------------
Her orneğin bollugunu, o ornekte VERININ DESTEKLEDIGI en dar rutbede verir.
Rutbe elle secilmez: okumalarin belirli bir orani (varsayilan %50) hangi
rutbede ya da altinda yerlesebiliyorsa o rutbe secilir. Yerlesemeyen
okumalar gizlenmez, ayri satirda "daha ust rutbede kaldi" olarak durur.

Kullanim:
  python3 abundance_rank.py --kraken kraken_c0.02 --out bolluk_c0.02
"""
import argparse, csv, glob, os, re, sys
from collections import defaultdict

# Kraken2 rutbe kodlari, genisten dara. Alt rutbeler (P1, C2 gibi) ana
# rutbeye katlanir; ara rutbeler taksonomi surumune gore degisir ve
# karsilastirmayi kirilgan yapar.
RUTBE_SIRA = ["D", "K", "P", "C", "O", "F", "G", "S"]
RUTBE_ADI = {"D": "alem üstü", "K": "alem", "P": "şube", "C": "sınıf",
             "O": "takım", "F": "aile", "G": "cins", "S": "tür"}

GRUP_ARALIK = [("A1", 1, "Arke, kısa amplikon"), ("A2", 5, "Arke, uzun amplikon"),
               ("F2", 9, "Mantar, uzun amplikon"), ("F1", 13, "Mantar, kısa amplikon"),
               ("B", 17, "Bakteri")]
YILLAR = (2021, 2023, 2024, 2025)
BARKOD_GRUP, BARKOD_YIL = {}, {}
for _g, _b0, _ac in GRUP_ARALIK:
    for _i, _y in enumerate(YILLAR):
        BARKOD_GRUP[_b0 + _i] = _g
        BARKOD_YIL[_b0 + _i] = _y


def ana_rutbe(kod):
    """P1 -> P, S2 -> S, U -> U. Bilinmeyen kod None doner."""
    if not kod:
        return None
    k = kod[0]
    return k if k in RUTBE_SIRA else None


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--kraken", required=True,
                   help="guven duzeltmesi uygulanmis rapor klasoru")
    p.add_argument("--out", required=True)
    p.add_argument("--desen", default="*_kraken2.report")
    p.add_argument("--kapsama", type=float, default=0.50,
                   help="bir rutbenin secilebilmesi icin siniflandirilmis "
                        "okumalarin en az bu orani o rutbede ya da altinda "
                        "yerlesmeli")
    p.add_argument("--ust", type=int, default=15, help="tabloda gosterilen takson")
    return p.parse_args()


def rapor_oku(yol):
    """Kraken2 raporunu okur. Girinti derinliginden ebeveyn zinciri de kurulur;
    ayni rutbedeki ic ice dugumleri ayirt etmek icin gerekli.
    Doner: (dogrudan{taxid:(rutbe, ad, sayi)}, klan{taxid:sayi},
            ebeveyn{taxid:taxid}, toplam)"""
    dogrudan, klan, ebeveyn = {}, {}, {}
    toplam = 0
    yol_yigin = []
    with open(yol, errors="replace") as fh:
        for line in fh:
            q = line.rstrip("\n").split("\t")
            if len(q) < 6:
                continue
            try:
                kl, dg, tx = int(q[1]), int(q[2]), int(q[4])
            except ValueError:
                continue
            rut, etiket = q[3], q[5]
            derinlik = (len(etiket) - len(etiket.lstrip(" "))) // 2
            yol_yigin = yol_yigin[:derinlik]
            ebeveyn[tx] = yol_yigin[-1] if yol_yigin else None
            yol_yigin.append(tx)
            dogrudan[tx] = (rut, etiket.strip(), dg)
            klan[tx] = kl
            toplam += dg
    return dogrudan, klan, ebeveyn, toplam


def main():
    a = get_args()
    dosyalar = sorted(glob.glob(os.path.join(a.kraken, a.desen))
                      + glob.glob(os.path.join(a.kraken, "*", a.desen)))
    if not dosyalar:
        sys.exit("rapor bulunamadi: %s" % a.kraken)
    os.makedirs(a.out, exist_ok=True)
    print("rapor: %d" % len(dosyalar))

    kapsama_satir, bolluk_satir, ozet_satir = [], [], []
    for yol in dosyalar:
        m = re.search(r"barcode(\d+)", os.path.basename(yol))
        if not m:
            print("   ATLANDI, barkod cozulemedi: %s" % os.path.basename(yol))
            continue
        bc = int(m.group(1))
        if bc not in BARKOD_GRUP:
            print("   ATLANDI, ornekleme haritasinda yok: barcode%02d" % bc)
            continue
        dogrudan, klan, ebeveyn, toplam = rapor_oku(yol)
        if toplam <= 0:
            print("   ATLANDI, bos rapor: %s" % os.path.basename(yol))
            continue
        sinifsiz = sum(d for r, _, d in dogrudan.values() if r == "U")
        siniflanan = toplam - sinifsiz

        # her ana rutbe icin: o rutbede YA DA ALTINDA dogrudan atanan okuma
        rut_dogrudan = defaultdict(int)
        for tx, (rut, ad, dg) in dogrudan.items():
            r = ana_rutbe(rut)
            if r:
                rut_dogrudan[r] += dg
        kumulatif = {}
        toplam_alt = 0
        for r in reversed(RUTBE_SIRA):          # S'den D'ye
            toplam_alt += rut_dogrudan.get(r, 0)
            kumulatif[r] = toplam_alt

        # secilen rutbe: kapsama esigini gecen EN DAR rutbe
        secilen = None
        for r in reversed(RUTBE_SIRA):          # once en dar
            if siniflanan and kumulatif[r] / siniflanan >= a.kapsama:
                secilen = r
                break
        if secilen is None:
            secilen = "D"

        for r in RUTBE_SIRA:
            kapsama_satir.append(dict(
                barkod="barcode%02d" % bc, grup=BARKOD_GRUP[bc], yil=BARKOD_YIL[bc],
                rutbe=r, rutbe_adi=RUTBE_ADI[r],
                bu_rutbede_ve_altinda=kumulatif[r],
                oran=round(100.0 * kumulatif[r] / siniflanan, 3) if siniflanan else 0.0,
                secilen=("EVET" if r == secilen else "")))

        # Secilen rutbede bolluk: o rutbedeki taksonlarin KLAN sayilari.
        #
        # DIKKAT: ayni ana rutbede IC ICE dugumler bulunabilir. Kraken2
        # raporunda gercek sube "P", alt sube "P1", onun alti "P2" diye
        # kodlanir ve hepsi ayni ana rutbeye katlanir. Ust dugumun klani
        # alt dugumun klanini zaten icerdigi icin ikisini birden saymak
        # cift sayimdir. Olculdu: duzeltmeden once barcode10'un yuzdeleri
        # toplami %368,54 cikiyordu.
        #
        # Cozum: yalnizca ayni ana rutbeden ATASI OLMAYAN dugumler alinir.
        # Bu kume tanim geregi ayriktir, klanlari ortusmez.
        def ust_ayni_rutbede(tx, r):
            p = ebeveyn.get(tx)
            gorulen = set()
            while p is not None and p not in gorulen:
                gorulen.add(p)
                if p in dogrudan and ana_rutbe(dogrudan[p][0]) == r:
                    return True
                p = ebeveyn.get(p)
            return False

        sec = []
        for tx, (rut, ad, dg) in dogrudan.items():
            if ana_rutbe(rut) == secilen and not ust_ayni_rutbede(tx, secilen):
                sec.append((ad, klan.get(tx, 0), tx))
        sec.sort(key=lambda z: -z[1])
        # ayrik kume denetimi: klan toplami, o rutbede ya da altinda yerlesen
        # okuma sayisini asamaz
        _t = sum(z[1] for z in sec)
        if _t > kumulatif[secilen] + 1:
            print("   UYARI barcode%02d: %s rutbesinde klan toplami (%d) "
                  "yerlesen okumayi (%d) asiyor, cift sayim olabilir"
                  % (bc, secilen, _t, kumulatif[secilen]))
        # o rutbenin ustunde kalanlar
        ust_kalan = siniflanan - kumulatif[secilen]
        for ad, n, tx in sec[:a.ust]:
            bolluk_satir.append(dict(
                barkod="barcode%02d" % bc, grup=BARKOD_GRUP[bc],
                yil=BARKOD_YIL[bc], rutbe=secilen, takson=ad, taxid=tx,
                okuma=n, yuzde=round(100.0 * n / toplam, 4)))
        if ust_kalan > 0:
            bolluk_satir.append(dict(
                barkod="barcode%02d" % bc, grup=BARKOD_GRUP[bc],
                yil=BARKOD_YIL[bc], rutbe=secilen,
                takson="[%s düzeyine inemeyen okuma]" % RUTBE_ADI[secilen],
                taxid="", okuma=ust_kalan,
                yuzde=round(100.0 * ust_kalan / toplam, 4)))
        if sinifsiz > 0:
            bolluk_satir.append(dict(
                barkod="barcode%02d" % bc, grup=BARKOD_GRUP[bc],
                yil=BARKOD_YIL[bc], rutbe=secilen, takson="[sınıflandırılmamış]",
                taxid="", okuma=sinifsiz,
                yuzde=round(100.0 * sinifsiz / toplam, 4)))

        ozet_satir.append(dict(
            barkod="barcode%02d" % bc, grup=BARKOD_GRUP[bc], yil=BARKOD_YIL[bc],
            toplam_okuma=toplam, siniflanan=siniflanan, sinifsiz=sinifsiz,
            secilen_rutbe=secilen, secilen_rutbe_adi=RUTBE_ADI[secilen],
            bu_rutbede_yerlesen=kumulatif[secilen],
            yerlesme_orani=round(100.0 * kumulatif[secilen] / siniflanan, 2)
            if siniflanan else 0.0,
            cins_orani=round(100.0 * kumulatif["G"] / siniflanan, 2)
            if siniflanan else 0.0,
            tur_orani=round(100.0 * kumulatif["S"] / siniflanan, 2)
            if siniflanan else 0.0))
        print("   barcode%02d %-3s toplam=%8d  secilen rutbe=%s (%s)  "
              "yerlesme=%.1f%%  cins=%.1f%%  tur=%.1f%%"
              % (bc, BARKOD_GRUP[bc], toplam, secilen, RUTBE_ADI[secilen],
                 100.0 * kumulatif[secilen] / siniflanan if siniflanan else 0,
                 100.0 * kumulatif["G"] / siniflanan if siniflanan else 0,
                 100.0 * kumulatif["S"] / siniflanan if siniflanan else 0))

    for adf, satirlar in (("rutbe_kapsamasi.tsv", kapsama_satir),
                          ("bolluk.tsv", bolluk_satir),
                          ("ozet.tsv", ozet_satir)):
        if not satirlar:
            continue
        with open(os.path.join(a.out, adf), "w", newline="",
                  encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(satirlar[0].keys()),
                               delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(satirlar)
    print("\nyazildi: %s" % a.out)

    # grup grup ozet
    print("\nGRUP BAZINDA VERININ DESTEKLEDIGI RUTBE")
    print("%-4s %-24s %-14s %10s %10s" % ("grup", "aciklama", "secilen rutbe",
                                          "cins orani", "tur orani"))
    for g, b0, ac in GRUP_ARALIK:
        ilgili = [x for x in ozet_satir if x["grup"] == g]
        if not ilgili:
            continue
        rutler = sorted(set(x["secilen_rutbe_adi"] for x in ilgili))
        co = sum(x["cins_orani"] for x in ilgili) / len(ilgili)
        to = sum(x["tur_orani"] for x in ilgili) / len(ilgili)
        print("%-4s %-24s %-14s %9.1f%% %9.1f%%"
              % (g, ac, ", ".join(rutler), co, to))
    print("\nCins orani, o gruptaki okumalarin yuzde kacinin cins ya da daha "
          "dar bir rutbeye yerlesebildigini gosterir. Dusuk oran, o gruptaki "
          "organizmalarin Kraken2 veritabaninda temsil edilmedigi anlamina "
          "gelir; sayilari cinse dagitmak veri uretmek olur.")


if __name__ == "__main__":
    main()
