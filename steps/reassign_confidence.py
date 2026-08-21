#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reassign_confidence.py
KRAKEN2 GUVEN ESIGINI VERITABANI OLMADAN, CIKTI DOSYALARINDAN UYGULAR.

Neden: Kraken2 varsayilan ayarinda (--confidence 0) cekimser kalmaz. Gercek
organizma veritabaninda yoksa okuma, hangi ayirt edici k-mer'leri saglam
kaldiysa en yuksek puanli kardes yaprakla etiketlenir. Olculdu: A2-4
ornegindeki dort Methanosarcina kutusunun okumalari ayni referanslara
gidiyor, hicbiri kendi atandigi turu tercih etmiyor. Cozum guven esigi
vermektir; ama bunun icin normalde 106 GB veritabaniyla saatler suren bir
yeniden siniflandirma gerekir.

Gerekmiyor. Kraken2'nin --output dosyasi her okuma icin k-mer LCA dizisini
(taxid:sayi ciftleri) zaten tasiyor. Guven puani tam olarak bu diziden
hesaplanir:
    guven = C / Q
    C = etiketin altindaki klanda yer alan k-mer sayisi
    Q = belirsiz baz icermeyen, yani veritabanina sorulan k-mer sayisi
        (eslesmeyenler, yani 0: girdileri Q'ya DAHILDIR; yalnizca A:
        girdileri haric tutulur)
Esigi gecemeyen okuma, gecene kadar takson agacinda yukari tasinir; hicbir
ata gecemezse siniflandirilmamis sayilir. Bu betik ayni islemi yapar.

Takson agaci rapor dosyalarindan kurulur (girinti derinligi ebeveyni verir);
butun raporlar birlestirilir. Olculdu: barcode01 ciktisindaki vurus
taksonlarinin %99,84'u bu agacta yer aliyor. Agacta bulunamayan vurus
KLANIN DISINDA sayilir; bu guveni dusurur, yani atamayi yukari iter, yani
guvenli yondedir. Bulunamayan orani her dosya icin ayrica raporlanir.

Cikti:
  <out>/<taban>_c<esik>_kraken2.report   Kraken2 rapor bicimi (Bracken okur)
  <out>/karsilastirma.tsv                okumanin hangi rutbeye tasindigi
  <out>/ozet.tsv

Kullanim:
  python3 reassign_confidence.py --kraken "kraken results" \
      --confidence 0.1 --out kraken_c0.1
"""
import argparse, csv, glob, os, re, sys
from collections import defaultdict, Counter

RUTBE_SIRA = ["U", "R", "D", "K", "P", "C", "O", "F", "G", "S"]


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--kraken", required=True,
                   help="'kraken results' directory (report and output files)")
    p.add_argument("--out", required=True)
    p.add_argument("--confidence", type=float, default=0.1)
    p.add_argument("--desen-rapor", default="*_kraken2.report")
    p.add_argument("--desen-cikti", default="*_output")
    p.add_argument("--max-okuma", type=int, default=0,
                   help="0 = all; deneme for kucultun")
    p.add_argument("--tarama", default=None,
                   help="comma-separated ayrilmis esikler; if given only tarama "
                        "yapilir, rapor yazilmaz. Ornek: 0,0.005,0.01,0.02,0.05")
    p.add_argument("--tarama-okuma", type=int, default=20000,
                   help="reads per file in scan mode")
    return p.parse_args()


def agac_kur(kraken_kok, desen):
    """Butun raporlardan takson agacini kurar.
    Doner: (ebeveyn, rutbe, ad, rapor_sayisi)"""
    ebeveyn, rutbe, ad = {}, {}, {}
    n = 0
    for rep in sorted(glob.glob(os.path.join(kraken_kok, "*", desen))
                      + glob.glob(os.path.join(kraken_kok, desen))):
        n += 1
        yol = []
        with open(rep, errors="replace") as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 6:
                    continue
                try:
                    tx = int(p[4])
                except ValueError:
                    continue
                etiket = p[5]
                derinlik = (len(etiket) - len(etiket.lstrip(" "))) // 2
                yol = yol[:derinlik]
                par = yol[-1] if yol else None
                # ilk gorulen ebeveyn tutulur; raporlar tutarli olmali
                if tx not in ebeveyn:
                    ebeveyn[tx] = par
                rutbe[tx] = p[3]
                ad[tx] = etiket.strip()
                yol.append(tx)
    return ebeveyn, rutbe, ad, n


def main():
    a = get_args()
    ebeveyn, rutbe, ad, nrep = agac_kur(a.kraken, a.desen_rapor)
    if not ebeveyn:
        sys.exit("rapor bulunamadi: %s" % a.kraken)
    print("rapor: %d, agactaki takson: %d" % (nrep, len(ebeveyn)))

    # atalar zinciri, bir kez hesaplanir
    ata_onbellek = {}

    def atalar(t):
        """t'den koke kadar olan zincir, t dahil. Bilinmeyen takson icin
        bos liste doner; boyle bir vurus hicbir klanda sayilmaz."""
        if t in ata_onbellek:
            return ata_onbellek[t]
        if t not in ebeveyn:
            ata_onbellek[t] = ()
            return ()
        zincir, k, gorulen = [], t, set()
        while k is not None and k not in gorulen:
            gorulen.add(k)
            zincir.append(k)
            k = ebeveyn.get(k)
        ata_onbellek[t] = tuple(zincir)
        return ata_onbellek[t]

    ciktilar = sorted(glob.glob(os.path.join(a.kraken, "*", a.desen_cikti))
                      + glob.glob(os.path.join(a.kraken, a.desen_cikti)))
    ciktilar = [c for c in ciktilar if not c.endswith(".report")]
    if not ciktilar:
        sys.exit("output dosyasi bulunamadi (desen: %s)" % a.desen_cikti)
    print("output dosyasi: %d" % len(ciktilar))
    os.makedirs(a.out, exist_ok=True)

    if a.tarama:
        esikler = [float(x) for x in a.tarama.split(",") if x.strip()]
        return tarama_yap(a, ciktilar, esikler, atalar, rutbe)

    ozet, karsilastirma = [], []
    for yol in ciktilar:
        taban = re.sub(r"_output$", "", os.path.basename(yol))
        dogrudan = Counter()          # yeni atama -> okuma sayisi
        eski_yeni = Counter()         # (eski_rutbe, yeni_rutbe) -> sayi
        n = tasinan = sinifsiz_olan = 0
        vurus_top = vurus_eksik = 0
        with open(yol, errors="replace") as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 5:
                    continue
                n += 1
                if a.max_okuma and n > a.max_okuma:
                    n -= 1
                    break
                eski = p[2]
                # --use-names kullanildiysa "Ad (taxid 1234)" bicimi gelir
                m = re.search(r"taxid\s+(\d+)", eski)
                if m:
                    eski_tx = int(m.group(1))
                else:
                    try:
                        eski_tx = int(eski)
                    except ValueError:
                        eski_tx = 0
                if eski_tx == 0:
                    dogrudan[0] += 1
                    eski_yeni[("U", "U")] += 1
                    continue
                # k-mer vuruslarini oku
                Q = 0
                klan = defaultdict(int)
                for tok in p[4].split():
                    t, _, c = tok.partition(":")
                    if not c:
                        continue
                    try:
                        say = int(c)
                    except ValueError:
                        continue
                    if t == "A":
                        continue          # belirsiz baz, Q'ya girmez
                    Q += say
                    if t == "0":
                        continue          # eslesme yok, Q'ya girer, klana girmez
                    try:
                        ti = int(t)
                    except ValueError:
                        continue
                    zincir = atalar(ti)
                    vurus_top += say
                    if not zincir:
                        vurus_eksik += say
                        continue
                    for k in zincir:
                        klan[k] += say
                if Q <= 0:
                    dogrudan[0] += 1
                    eski_yeni[(rutbe.get(eski_tx, "?"), "U")] += 1
                    sinifsiz_olan += 1
                    continue
                # eski atamanin kokten yaprağa yolunda, esigi gecen EN OZGUL
                # dugum secilir. atalar() yaprak-koke sirali doner.
                yeni_tx = 0
                for k in atalar(eski_tx):
                    if klan.get(k, 0) / Q >= a.confidence:
                        yeni_tx = k
                        break
                dogrudan[yeni_tx] += 1
                er = rutbe.get(eski_tx, "?")
                yr = "U" if yeni_tx == 0 else rutbe.get(yeni_tx, "?")
                eski_yeni[(er, yr)] += 1
                if yeni_tx != eski_tx:
                    tasinan += 1
                if yeni_tx == 0:
                    sinifsiz_olan += 1

        # klan sayimlari
        klan_say = Counter()
        for tx, c in dogrudan.items():
            if tx == 0:
                continue
            for k in atalar(tx):
                klan_say[k] += c
        toplam = n

        # Kraken2 rapor bicimi
        rap = os.path.join(a.out, "%s_c%g_kraken2.report" % (taban, a.confidence))
        cocuk = defaultdict(list)
        for tx in klan_say:
            par = ebeveyn.get(tx)
            if par is not None and par in klan_say:
                cocuk[par].append(tx)
        kokler = [tx for tx in klan_say if ebeveyn.get(tx) not in klan_say]
        with open(rap, "w", encoding="utf-8") as fh:
            u = dogrudan.get(0, 0)
            if toplam:
                fh.write("%6.2f\t%d\t%d\tU\t0\tunclassified\n"
                         % (100.0 * u / toplam, u, u))
            yigin = [(tx, 0) for tx in sorted(kokler,
                                              key=lambda x: -klan_say[x])]
            while yigin:
                tx, d = yigin.pop()
                fh.write("%6.2f\t%d\t%d\t%s\t%d\t%s%s\n"
                         % (100.0 * klan_say[tx] / toplam if toplam else 0.0,
                            klan_say[tx], dogrudan.get(tx, 0),
                            rutbe.get(tx, "-"), tx, "  " * d, ad.get(tx, str(tx))))
                for c in sorted(cocuk.get(tx, []), key=lambda x: klan_say[x]):
                    yigin.append((c, d + 1))

        ozet.append(dict(
            dosya=taban, okuma=toplam, tasinan=tasinan,
            tasinan_yuzde=round(100.0 * tasinan / toplam, 2) if toplam else 0,
            sinifsiz_olan=sinifsiz_olan,
            agacta_bulunamayan_vurus_yuzde=(
                round(100.0 * vurus_eksik / vurus_top, 3) if vurus_top else 0.0),
            rapor=os.path.basename(rap)))
        for (er, yr), c in sorted(eski_yeni.items(),
                                  key=lambda x: -x[1]):
            karsilastirma.append(dict(dosya=taban, eski_rutbe=er,
                                      yeni_rutbe=yr, okuma=c))
        print("   %-34s okuma=%8d tasinan=%7d (%5.2f%%) sinifsiz=%7d "
              "agac_disi_vurus=%.3f%%"
              % (taban, toplam, tasinan,
                 100.0 * tasinan / toplam if toplam else 0, sinifsiz_olan,
                 100.0 * vurus_eksik / vurus_top if vurus_top else 0))

    for adf, satirlar in (("ozet.tsv", ozet), ("karsilastirma.tsv", karsilastirma)):
        if not satirlar:
            continue
        with open(os.path.join(a.out, adf), "w", newline="",
                  encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(satirlar[0].keys()),
                               delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(satirlar)
    print("\nyazildi: %s" % a.out)
    print("guven esigi: %g" % a.confidence)
    tt = sum(x["okuma"] for x in ozet)
    ts = sum(x["tasinan"] for x in ozet)
    print("toplam okuma: %d, rutbesi degisen: %d (%.2f%%)"
          % (tt, ts, 100.0 * ts / tt if tt else 0))
    print("\nBu raporlar Bracken'a dogrudan verilebilir. Veritabani ve ham "
          "fastq gerekmez; hesap Kraken2'nin kendi k-mer ciktisindan yapilir.")


def tarama_yap(a, ciktilar, esikler, atalar, rutbe):
    """Esik secimini kural ezberinden degil VERIDEN yapmak icin.

    Olculdu (bu veri, 20 dosya): ONT okumalarinin k-mer'lerinin buyuk
    cogunlugu veritabaninda karsilik bulmuyor ve bu k-mer'ler guven
    puaninin PAYDASINA giriyor. Bu yuzden kisa okuma verisi icin sik
    onerilen 0,1 esigi burada okumalarin yarisindan cogunu
    siniflandirilmamis birakiyor. Dogru esik veriden secilmeli.

    Her esik icin uc sayi verilir:
      tur_okuma   tur (S) rutbesinde kalan okuma
      cins_okuma  cins (G) rutbesine tasinan ya da orada kalan okuma
      sinifsiz    hicbir atayi gecemeyen okuma
    Amac tur duzeyindeki sahte ayrimi cins duzeyinde toplamak, okuma
    kaybetmek degil; dolayisiyla cins_okuma artarken sinifsiz'in dusuk
    kaldigi en yuksek esik secilmelidir."""
    print("\nESIK TARAMASI (dosya basina en fazla %d okuma)" % a.tarama_okuma)
    print("%8s %10s %10s %10s %10s %10s"
          % ("esik", "okuma", "tur(S)", "cins(G)", "ust_rutbe", "sinifsiz"))
    satirlar = []
    for esik in esikler:
        top = Counter()
        n_top = 0
        for yol in ciktilar:
            n = 0
            with open(yol, errors="replace") as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 5:
                        continue
                    n += 1
                    if n > a.tarama_okuma:
                        break
                    n_top += 1
                    m = re.search(r"taxid\s+(\d+)", p[2])
                    try:
                        eski_tx = int(m.group(1)) if m else int(p[2])
                    except ValueError:
                        eski_tx = 0
                    if eski_tx == 0:
                        top["U"] += 1
                        continue
                    Q = 0
                    klan = defaultdict(int)
                    for tok in p[4].split():
                        t, _, c = tok.partition(":")
                        if not c:
                            continue
                        try:
                            say = int(c)
                        except ValueError:
                            continue
                        if t == "A":
                            continue
                        Q += say
                        if t == "0":
                            continue
                        try:
                            ti = int(t)
                        except ValueError:
                            continue
                        for k in atalar(ti):
                            klan[k] += say
                    if Q <= 0:
                        top["U"] += 1
                        continue
                    yeni = 0
                    for k in atalar(eski_tx):
                        if klan.get(k, 0) / Q >= esik:
                            yeni = k
                            break
                    if yeni == 0:
                        top["U"] += 1
                    else:
                        r = rutbe.get(yeni, "?")
                        if r == "S" or r.startswith("S"):
                            top["S"] += 1
                        elif r == "G" or r.startswith("G"):
                            top["G"] += 1
                        else:
                            top["ust"] += 1
        print("%8g %10d %10d %10d %10d %10d"
              % (esik, n_top, top["S"], top["G"], top["ust"], top["U"]))
        satirlar.append(dict(esik=esik, okuma=n_top, tur=top["S"],
                             cins=top["G"], ust_rutbe=top["ust"],
                             sinifsiz=top["U"]))
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "esik_taramasi.tsv"), "w", newline="",
              encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(satirlar[0].keys()),
                           delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(satirlar)
    print("\nyazildi: %s" % os.path.join(a.out, "esik_taramasi.tsv"))
    print("Esigi buradan secin: cins(G) sutunu yukselirken sinifsiz sutunu "
          "hala dusukse o esik uygundur. Sectiginiz degeri --confidence ile "
          "verip --tarama'siz calistirin.")
    return 0


if __name__ == "__main__":
    main()
