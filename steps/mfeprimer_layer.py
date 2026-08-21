#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mfeprimer_layer.py
DIS VERITABANI OZGULLUGUNUN IKINCI, BAGIMSIZ OLCUMU.

external_databases.py blastn kullanir: her primeri ayri ayri arar, sonra
vuruslari kendisi eslestirip urun olusup olusmadigina karar verir.
Bu betik mfeprimer'in 'spec' alt komutunu kullanir: mfeprimer amplikonu
kendi termodinamik modeliyle, k-mer indeksinden dogrudan tahmin eder.
Iki yontem birbirinden bagimsizdir; toplanti kararindaki "hicbir karar tek
koda birakilmasin" ilkesi dis ozgullukte ancak ikisi birden calisinca
saglanir.

Sonuc uc sutunda verilir:
  blast_urun     14'un buldugu hedef disi urun sayisi
  mfe_urun       mfeprimer'in buldugu hedef disi amplikon sayisi
  uyum           iki_olcum_uyustu | ayrisan_olcum | tek_olcum

ONEMLI: mfeprimer toplanti kararindaki "son iki baz birebir" kuralini
uygulamaz. Olculdu: 3' son bazi degistirilmis bir primer, bozulmamis
primerle ayni sayida amplikon veriyor. Bu yuzden mfeprimer sonuclari eleme
olcutu degil, ikinci bir bakis acisidir.

"Uyum" mutlak sayilarin esitligi demek DEGILDIR; iki yontemin esik ve
model farklari sayilari kacinilmaz olarak ayirir. Uyum, ikisinin de AYNI
KARARI vermesidir: sifir mi, sifir degil mi. Ayrisma sessizce gecilmez,
ayri bir sutunda ve ozet satirinda raporlanir.

Kullanim:
  python3 mfeprimer_layer.py --final primer_final --db REFERANS_DB \
      --mfe tools/mfeprimer --out primer_final/mfeprimer.tsv
"""
import argparse, csv, collections, os, re, shutil, subprocess, sys, tempfile

# 14 ile AYNI sinif-veritabani eslemesi; ikisi ayrisirsa karsilastirma
# anlamsiz olur.
SINIF_DB = {
    "A1": ["archaea.16S.fna"],
    "A2": ["archaea.16S.fna"],
    "B":  ["bacteria.16S.fna"],
    "F1": ["fungi.ITS.fna", "fungi.28SrRNA.fna"],
    "F2": ["fungi.ITS.fna", "fungi.28SrRNA.fna"],
}

# 14 ile ayni genis kume; iki olcumun ayni veritabanlarini gormesi sart,
# yoksa "ayrisan olcum" satirlari yontem farkini degil kapsam farkini
# gosterir. Liste external_databases.py'deki ile BIREBIR AYNI olmak
# zorundadir; 14'te ROD'un okaryot-yalniz oldugu olculdukten sonra
# (60320/60320 Eukaryota, 0 Bacteria, 0 Archaea) alan bazli elle secim
# kaldirildi ve tum siniflar tum rDNA veritabanlarini gorur oldu.
GENIS_ORTAK = [
    "SILVA_138.2_SSURef_NR99.fasta",
    "SILVA_138.2_LSURef_NR99.fasta",
    "ROD_v1.2_operon_variants.fasta",
    "UNITE_ITS.fasta",
    "PR2_SSU_taxo_long.fasta",
    "archaea.16S.fna",
    "bacteria.16S.fna",
    "fungi.ITS.fna",
    "fungi.18SrRNA.fna",
    "fungi.28SrRNA.fna",
]
SINIF_DB_GENIS = {
    s: [d for d in GENIS_ORTAK if d not in SINIF_DB[s]] for s in SINIF_DB
}


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--final", required=True, help="09'un output directory")
    p.add_argument("--db", required=True, help="REFERANS_DB directory")
    p.add_argument("--mfe", required=True, help="mfeprimer calistirilabiliri")
    p.add_argument("--blast", default=None,
                   help="output of the external-databases step; if given, the two measurements are compared")
    p.add_argument("--out", required=True)
    p.add_argument("--prod-min", type=int, default=50)
    p.add_argument("--prod-max", type=int, default=400)
    p.add_argument("--tm-min", type=float, default=30.0,
                   help="mfeprimer amplikon Tm alt siniri")
    p.add_argument("--mismatch", type=int, default=3,
                   help="toplanti kararindaki toplam mismatch siniri")
    p.add_argument("--mis-end", type=int, default=3,
                   help="mfeprimer'in mismatch penceresi. 3 hizli and "
                        "secici, 9 (mfeprimer varsayilani) asiri gevsek ve "
                        "cok yavas. Bu deger toplanti kararindaki 3' uc "
                        "kuralini UYGULAMAZ, betik basindaki nota bakin.")
    p.add_argument("--genis", action="store_true",
                   help="also scan the same wide database set as the external-databases step")
    p.add_argument("--cpu", type=int, default=4)
    p.add_argument("--zaman-asimi", type=int, default=3600)
    return p.parse_args()


def indeks_eksik(fna):
    """mfeprimer 'spec' icin gereken yardimci dosyalar. Yalnizca .primerqc.bin
    bakmak yetmez: .fai ve .json yoksa mfeprimer 'no valid db found' yazar
    ve CIKIS KODU 0 DONER, yani sessizce hicbir sey uretmez. Eksik dosya
    adlari doner, bos liste her seyin yerinde oldugunu gosterir."""
    eksik = []
    if not (os.path.exists(fna + ".primerqc.bin")
            or os.path.exists(fna + ".primerqc")):
        eksik.append(".primerqc.bin")
    for ek in (".fai", ".json"):
        if not os.path.exists(fna + ek):
            eksik.append(ek)
    return eksik


def main():
    a = get_args()
    if not os.path.exists(a.mfe):
        sys.exit("mfeprimer bulunamadi: %s" % a.mfe)
    if not os.access(a.mfe, os.X_OK):
        try:
            os.chmod(a.mfe, 0o755)
        except OSError:
            sys.exit("mfeprimer calistirilabilir degil: %s\n"
                     "   chmod +x '%s' deneyin" % (a.mfe, a.mfe))
    tsv = os.path.join(a.final, "primer_final.tsv")
    if not os.path.exists(tsv):
        sys.exit("bulunamadi: %s" % tsv)
    rows = [r for r in csv.DictReader(open(tsv, encoding="utf-8"), delimiter="\t")
            if r.get("ozgulluk_durum") == "GECTI"]
    if not rows:
        sys.exit("gecen aday yok")
    print("sinanacak cift: %d" % len(rows))

    # 14'un sonucu: (hedef, sinif, ileri, geri, veritabani) -> urun sayisi
    blast = {}
    byol = a.blast or os.path.join(a.final, "dis_veritabani.tsv")
    if os.path.exists(byol):
        for x in csv.DictReader(open(byol, encoding="utf-8"), delimiter="\t"):
            k = (x["hedef"], x["sinif"], x["ileri_dizi"], x["geri_dizi"],
                 x["veritabani"])
            blast[k] = blast.get(k, 0) + int(x.get("hedef_disi_urun", 0) or 0)
        print("14'un ciktisi okundu: %s (%d kayit)" % (byol, len(blast)))
    else:
        print("UYARI: 14'un ciktisi yok, karsilastirma yapilamayacak: %s" % byol)

    calisma = tempfile.mkdtemp(prefix="mfe_")
    sonuc = []
    hata = 0
    # sinif basina tek cagri: mfeprimer tsv girdisinde her satir bir cifttir
    sinif_cift = collections.defaultdict(list)
    for i, r in enumerate(rows):
        sinif_cift[r["sinif"]].append((i, r))

    for sinif, ciftler in sorted(sinif_cift.items()):
        dblist = list(SINIF_DB.get(sinif, []))
        if a.genis:
            dblist += SINIF_DB_GENIS.get(sinif, [])
        for dbad in dblist:
            fna = os.path.join(a.db, dbad)
            if not os.path.exists(fna):
                print("   veritabani yok, atlandi: %s" % fna)
                continue
            eks = indeks_eksik(fna)
            if eks:
                print("   mfeprimer indeksi eksik (%s): %s"
                      % (", ".join(eks), dbad))
                print("   Kurmak icin: %s index -i %s -c %d"
                      % (a.mfe, fna, a.cpu))
                hata += 1
                continue
            girdi = os.path.join(calisma, "%s_%s.tsv" % (sinif, dbad))
            with open(girdi, "w", encoding="utf-8") as fh:
                for i, r in ciftler:
                    fh.write("p%d\t%s\t%s\n" % (i, r["ileri_dizi"], r["geri_dizi"]))
            cikti = girdi + ".mfe.tsv"
            # mfeprimer'in KENDI modeli calisir; toplanti kararindaki
            # baglanma kurali burada TAKLIT EDILMEZ. Taklit edilseydi iki
            # olcum bagimsiz olmaz ve ikinci olcumun anlami kalmazdi.
            #
            # OLCULDU (archaea.16S, 2026-08-01): mfeprimer --misEnd 3 ile
            # calistirildiginda, 3' UCTAKI SON BAZI degistirilmis bir primer
            # bozulmamis primerle AYNI sayida amplikon veriyor (323'e 323).
            # Yani mfeprimer "son iki baz birebir uymali" kuralini
            # uygulamiyor; kendi k-mer cekirdegini baska yerde kuruyor.
            # --misEnd 9 (varsayilan) ise ayni primer icin 30 000'den fazla
            # amplikon veriyor ve pratikte kullanilamayacak kadar yavas.
            # Bu yuzden mfeprimer sayilari BIR ELEME OLCUTU DEGILDIR;
            # blastn'in bulmadigi bir urun bulursa bu, o cift icin
            # laboratuvarda ayrica bakilmasi gereken bir uyaridir.
            cmd = [a.mfe, "spec", "-i", girdi, "-o", cikti, "-d", fna,
                   "-s", str(a.prod_min), "-S", str(a.prod_max),
                   "-t", str(a.tm_min), "-c", str(a.cpu),
                   "--misMatch", str(a.mismatch),
                   "--misStart", "1", "--misEnd", str(a.mis_end)]
            print("   mfeprimer %-4s x %-22s (%d cift)"
                  % (sinif, dbad, len(ciftler)))
            try:
                p = subprocess.run(cmd, capture_output=True, text=True,
                                   timeout=a.zaman_asimi)
            except subprocess.TimeoutExpired:
                print("      ZAMAN ASIMI (%d sn), atlandi" % a.zaman_asimi)
                continue
            ciktisi = (p.stdout or "") + (p.stderr or "")
            if p.returncode != 0:
                print("      HATA: %s" % ciktisi.strip()[:220])
                hata += 1
                continue
            # mfeprimer bazi hatalarda CIKIS KODU 0 doner ve yalnizca
            # ekrana yazar. Cikis kodunu tek olcut saymak, bu durumda
            # "hedef disi urun yok" gibi okunan sahte bir temizlik uretir.
            if "no valid db" in ciktisi.lower() or "error" in ciktisi.lower():
                print("      HATA (cikis kodu 0 ama ileti var): %s"
                      % ciktisi.strip()[:220])
                hata += 1
                continue
            # mfeprimer ciktisini oku: amplikon satirlarini cift basina say
            say = collections.Counter()
            ornek = collections.defaultdict(list)
            # mfeprimer iki dosya yazar: <out> insan icin okunabilir rapor,
            # <out>.spec.tsv makine icin. Ikincisi okunur; ilkinin bicimi
            # surumler arasinda degisiyor ve sayim icin guvenilir degil.
            okunan = None
            for c in (cikti + ".spec.tsv", cikti + ".tsv", cikti):
                if os.path.exists(c) and c.endswith(".spec.tsv"):
                    okunan = c
                    break
            if not okunan:
                print("      .spec.tsv olusmadi, bu veritabani OLCULEMEDI. "
                      "mfeprimer iletisi: %s" % (ciktisi.strip()[:160] or "yok"))
                hata += 1
                continue
            with open(okunan, encoding="utf-8", errors="replace") as fh:
                basliklar = None
                for line in fh:
                    if not line.strip():
                        continue
                    p2 = line.rstrip("\n").split("\t")
                    if basliklar is None and p2[0].startswith("#1-based"):
                        continue          # dosyanin ilk aciklama satiri
                    if basliklar is None:
                        basliklar = [x.strip().lstrip("#") for x in p2]
                        continue
                    d = dict(zip(basliklar, p2))
                    # fpName ve rpName ayni cifte ait olmali; degilse bu
                    # amplikon iki farkli ciftin primerlerinden olusmustur
                    # ve o cifte yazilmaz.
                    mf = re.search(r"p(\d+)", d.get("fpName", ""))
                    mr = re.search(r"p(\d+)", d.get("rpName", ""))
                    if not mf or not mr or mf.group(1) != mr.group(1):
                        continue
                    idx = mf.group(1)
                    say[idx] += 1
                    if len(ornek[idx]) < 5:
                        ornek[idx].append("%s:%s bp" % (d.get("chrom", "")[:28],
                                                        d.get("ampSize", "")))
            for i, r in ciftler:
                k = (r["hedef"], sinif, r["ileri_dizi"], r["geri_dizi"], dbad)
                b = blast.get(k)
                m = say.get(str(i), 0)
                if b is None:
                    uyum = "tek_olcum"
                elif (b == 0) == (m == 0):
                    uyum = "iki_olcum_uyustu"
                else:
                    uyum = "ayrisan_olcum"
                sonuc.append(dict(
                    hedef=r["hedef"], sinif=sinif, veritabani=dbad,
                    ileri_dizi=r["ileri_dizi"], geri_dizi=r["geri_dizi"],
                    blast_urun=("" if b is None else b), mfe_urun=m,
                    uyum=uyum, ornekler=";".join(ornek.get(str(i), []))))
    shutil.rmtree(calisma, ignore_errors=True)

    if not sonuc:
        print("hicbir veritabani taranamadi; cikti yazilmadi")
        # Bos cikmasi sessizce 'temiz' okunmamali; bos ama basliklari olan
        # bir dosya yazilir ki bayat cikti hayatta kalmasin.
        if a.out:
            d = os.path.dirname(os.path.abspath(a.out))
            if d:
                os.makedirs(d, exist_ok=True)
            with open(a.out, "w", encoding="utf-8") as fh:
                fh.write("hedef\tsinif\tveritabani\tileri_dizi\tgeri_dizi\t"
                         "blast_urun\tmfe_urun\tuyum\tornekler\n")
        return 2

    d = os.path.dirname(os.path.abspath(a.out))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sonuc[0].keys()),
                           delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(sonuc)
    print("\nyazildi: %s" % a.out)

    say_uyum = collections.Counter(x["uyum"] for x in sonuc)
    print("kayit: %d" % len(sonuc))
    for k in ("iki_olcum_uyustu", "ayrisan_olcum", "tek_olcum"):
        if say_uyum.get(k):
            print("   %-20s %d" % (k, say_uyum[k]))
    ayrisan = [x for x in sonuc if x["uyum"] == "ayrisan_olcum"]
    if ayrisan:
        print("\nAYRISAN OLCUMLER (biri urun buldu, oteki bulmadi):")
        for x in ayrisan[:20]:
            print("   %-30s %-3s %-18s blast=%-5s mfe=%-5s"
                  % (x["hedef"][:29], x["sinif"], x["veritabani"][:17],
                     x["blast_urun"], x["mfe_urun"]))
        if len(ayrisan) > 20:
            print("   ... %d kayit daha" % (len(ayrisan) - 20))
        print("Bu satirlar elenmis sayilmaz; iki yontemin ayrildigi yerlerdir "
              "ve laboratuvarda ayri dikkat gerektirir.")
    temiz = sorted(set((x["hedef"], x["sinif"]) for x in sonuc
                       if x["mfe_urun"] == 0))
    print("\nmfeprimer'e gore hedef disi amplikon vermeyen hedef-sinif: %d"
          % len(temiz))
    if hata:
        print("\nDIKKAT: %d veritabani olculemedi. Bu veritabanlari icin "
              "ikinci olcum YAPILMAMISTIR; eksikligi sonuc dosyasinda "
              "'tek_olcum' olarak degil, hic satir olmayarak gorunur." % hata)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
