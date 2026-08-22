#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OZEL KRAKEN2 VERITABANI ICIN TAKSONOMI VE KUTUPHANE URETIMI  (EN SON CARE)

ONCE BUNU OKUYUN
PlusPFP kuruluysa BU BETIGE GEREK YOKTUR. PlusPFP zaten Standard'a protozoa,
mantar ve bitki ekler, yani eksik oldugunu olctugumuz gruplarin tamami.
Bu betik, PlusPFP hicbir sekilde kurulamiyorsa ya da marker gen (16S/18S/ITS)
duzeyinde ikinci bir gorus isteniyorsa vardir.

NE YAPAR
SILVA, UNITE ve PR2 dosyalarinin basliklarindaki soy dizgilerinden sentetik bir
NCBI benzeri taksonomi kurar (nodes.dmp, names.dmp) ve dizileri kraken2'nin
bekledigi `kraken:taxid|N` bicimiyle library/ altina yazar. Sonrasini
kraken2-build --build yapar.

TAXID'LER SENTETIKTIR. NCBI taxid'leriyle karsilastirilamaz. Karsilastirma bu
yuzden isim duzeyinde yapilir (comparison_table.py boyle davranir).

Calistirma:
  python3 custom_taxonomy.py --cikti <db klasoru> --kume silva_ssu=/yol/a.fasta ...
  python3 custom_taxonomy.py --selftest
"""
import argparse, os, re, sys

RANKLAR = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
# UNITE basliklarindaki rank onekleri: k__Fungi, p__Ascomycota, ...
ONEK = re.compile(r"^[kdpcofgs]__")

def soy_ayikla(baslik):
    """
    Fasta basligindan soy dizgisini cikarir ve duzeylere boler.

    Uc bicim de desteklenir:
      SILVA : AB000389.1.1487 Bacteria;Pseudomonadota;...;Genus species
      PR2   : AB000000.1.1000 Eukaryota;TSAR;Alveolata;...
      UNITE : Ad|KY123|SH123.08FU|reps|k__Fungi;p__Ascomycota;...;s__Petriella_x
    Doner: (dizi_kimligi, [duzey, duzey, ...])
    """
    b = baslik.lstrip(">").rstrip()
    if not b:
        return "", []
    # Soy dizgisi noktali virgul iceren parcadir. Once bosluga, sonra boruya bak.
    kimlik = b.split()[0]
    kalan = b[len(kimlik):].strip()
    if ";" not in kalan:
        parcalar = b.split("|")
        kalan = next((p for p in parcalar if ";" in p), "")
        if parcalar and parcalar[0]:
            kimlik = parcalar[0].split()[0] if parcalar[0].split() else kimlik
        # UNITE'ta daha kararli kimlik SH kodudur, varsa o kullanilir.
        sh = next((p for p in parcalar if p.startswith("SH")), "")
        if sh:
            kimlik = sh
    if ";" not in kalan:
        return kimlik, []
    duzeyler = []
    for p in kalan.split(";"):
        p = ONEK.sub("", p.strip()).replace("_", " ").strip()
        if not p or p.lower() in {"unidentified", "incertae sedis", "na", "unclassified"}:
            continue
        duzeyler.append(p)
    return kimlik, duzeyler

class Taksonomi:
    """
    Sentetik agac. Kok taxid 1. Her benzersiz soy YOLU bir dugumdur; ayni ad
    farkli soylarda gecerse ayri dugum olur, cunku 'Incertae sedis' gibi adlar
    birden cok yerde geciyor ve birlestirilirse agac catallanip bozulur.
    """
    def __init__(self):
        self.ebeveyn = {1: 1}
        self.ad = {1: "root"}
        self.rank = {1: "no rank"}
        self.yol = {(): 1}
        self.sonraki = 2

    def ekle(self, duzeyler):
        anahtar = ()
        tx = 1
        for i, d in enumerate(duzeyler):
            anahtar = anahtar + (d,)
            if anahtar in self.yol:
                tx = self.yol[anahtar]
                continue
            yeni = self.sonraki
            self.sonraki += 1
            self.yol[anahtar] = yeni
            self.ebeveyn[yeni] = tx
            self.ad[yeni] = d
            self.rank[yeni] = RANKLAR[i] if i < len(RANKLAR) else "no rank"
            tx = yeni
        return tx

    def yaz(self, klasor):
        os.makedirs(klasor, exist_ok=True)
        with open(os.path.join(klasor, "nodes.dmp"), "w", encoding="utf-8") as fh:
            for tx in sorted(self.ebeveyn):
                fh.write(f"{tx}\t|\t{self.ebeveyn[tx]}\t|\t{self.rank[tx]}\t|\t-\t|\n")
        with open(os.path.join(klasor, "names.dmp"), "w", encoding="utf-8") as fh:
            for tx in sorted(self.ad):
                fh.write(f"{tx}\t|\t{self.ad[tx]}\t|\t\t|\tscientific name\t|\n")
        return len(self.ebeveyn)

def kume_isle(yol, tak, cikti_fh):
    """Bir fasta dosyasini okur, taksonomiyi buyutur, kutuphane fastasini yazar."""
    okunan = soysuz = yazilan = 0
    tx = None
    with open(yol, errors="replace") as fh:
        for satir in fh:
            if satir.startswith(">"):
                okunan += 1
                kimlik, duzeyler = soy_ayikla(satir)
                if not duzeyler:
                    soysuz += 1
                    tx = None
                    continue
                tx = tak.ekle(duzeyler)
                yazilan += 1
                cikti_fh.write(f">{kimlik}|kraken:taxid|{tx} {duzeyler[-1]}\n")
            elif tx is not None:
                # SILVA RNA olarak gelir, kraken2 DNA bekler. U -> T cevrilir.
                cikti_fh.write(satir.upper().replace("U", "T"))
    return okunan, yazilan, soysuz

def selftest():
    print("=" * 72)
    print(u'THE CUSTOM TAXONOMY, A TEST WITH KNOWN ANSWERS')
    print("=" * 72)
    hata = 0
    def K(ad, bul, bek):
        nonlocal hata
        ok = bul == bek
        if not ok:
            hata += 1
        print(f"  {'PASS' if ok else 'FAIL'}  {ad:<58} {bul} / {bek}")

    k, d = soy_ayikla(">AB000389.1.1487 Bacteria;Pseudomonadota;Gammaproteobacteria")
    K(u'a SILVA header, the identity', k, "AB000389.1.1487")
    K("SILVA basligi, soy", d, ["Bacteria", "Pseudomonadota", "Gammaproteobacteria"])
    k2, d2 = soy_ayikla(">Petriella|KY123456|SH1234567.08FU|reps|k__Fungi;p__Ascomycota;g__Petriella")
    K(u'a UNITE header, the SH identity is chosen', k2, "SH1234567.08FU")
    K("UNITE rank onekleri atilir", d2, ["Fungi", "Ascomycota", "Petriella"])
    k3, d3 = soy_ayikla(">X1 Eukaryota;TSAR;Alveolata;unidentified;Genus x")
    K(u'empty and meaningless levels are skipped', d3, ["Eukaryota", "TSAR", "Alveolata", "Genus x"])
    K(u'a header with no lineage returns an empty list', soy_ayikla(">SADECE_KIMLIK")[1], [])
    K("bos baslik cokmez", soy_ayikla(">")[1], [])

    t = Taksonomi()
    a = t.ekle(["Bacteria", "Firmicutes", "Bacilli"])
    b = t.ekle(["Bacteria", "Firmicutes", "Clostridia"])
    K("ortak soy paylasilir, iki kez yaratilmaz", t.ebeveyn[a], t.ebeveyn[b])
    K("kok 1'dir", t.ebeveyn[t.ebeveyn[t.ebeveyn[a]]], 1)
    K(u'the rank is assigned by depth', t.rank[a], "class")
    c = t.ekle(["Archaea", "Firmicutes"])
    K(u'the same name in a different lineage is a SEPARATE node', t.ebeveyn[c] != t.ebeveyn[t.ebeveyn[a]], True)

    import tempfile
    with tempfile.TemporaryDirectory() as dizin:
        fa = os.path.join(dizin, "g.fasta")
        open(fa, "w").write(">S1 Bacteria;Firmicutes;Bacilli\nAUGCAUGC\n"
                            ">S2 SOYSUZ\nAAAA\n"
                            ">S3 Bacteria;Firmicutes;Clostridia\nGGGG\n")
        t2 = Taksonomi()
        cy = os.path.join(dizin, "lib.fna")
        with open(cy, "w") as fh:
            ok, yz, sz = kume_isle(fa, t2, fh)
        K(u'every header is read', ok, 3)
        K(u'a sequence with no lineage is skipped', sz, 1)
        K(u'the sequences with a lineage are written', yz, 2)
        icerik = open(cy).read()
        K("baslik kraken:taxid tasir", "|kraken:taxid|" in icerik, True)
        K("RNA'daki U, DNA'da T olur", "ATGCATGC" in icerik, True)
        K(u'the sequence of a record with no lineage is not written either', "AAAA" in icerik, False)
        n = t2.yaz(os.path.join(dizin, "taxonomy"))
        K(u'nodes.dmp is written', os.path.exists(os.path.join(dizin, "taxonomy", "nodes.dmp")), True)
        satir = open(os.path.join(dizin, "taxonomy", "nodes.dmp")).readline()
        K("nodes.dmp bicimi kraken2'nin bekledigi gibi", satir.startswith("1\t|\t1\t|\t"), True)
        # kok + Bacteria + Firmicutes + Bacilli + Clostridia = 5
        K(u'the node count including the root', n, 5)

    print("=" * 72)
    print("THE TEST PASSED" if hata == 0 else f"THE TEST FAILED, {hata} items")
    print("=" * 72)
    return 0 if hata == 0 else 1

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cikti", default="")
    ap.add_argument("--kume", action="append", default=[],
                    help="ad=/yol/file.fasta biciminde, birden cok kez verilebilir")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(selftest())
    if not a.cikti:
        ap.error("--cikti gerekli")
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        if selftest() != 0:
            print("THE TEST FAILED, stopped (project rule 2)")
            sys.exit(2)
    if not a.kume:
        ap.error("en az bir --kume gerekli")

    lib = os.path.join(a.cikti, "library", "ozel")
    os.makedirs(lib, exist_ok=True)
    tak = Taksonomi()
    hedef = os.path.join(lib, "library.fna")
    toplam_ok = toplam_yz = toplam_sz = 0
    with open(hedef, "w") as fh:
        for k in a.kume:
            if "=" not in k:
                print(f"ERROR: --kume must have the form name=/path, received: {k}")
                sys.exit(1)
            ad, yol = k.split("=", 1)
            if not os.path.exists(yol):
                print(f"ERROR: no such file: {yol}")
                sys.exit(1)
            print(f"  {ad}: {os.path.basename(yol)}")
            ok, yz, sz = kume_isle(yol, tak, fh)
            print(f"     {ok} headers, {yz} written, {sz} skipped for having no lineage")
            if sz > ok * 0.5:
                print(f"     WARNING: more than half the headers have no lineage. The header format")
                print(f"     may differ from what is expected, so the result is not trustworthy.")
            toplam_ok += ok; toplam_yz += yz; toplam_sz += sz
    n = tak.yaz(os.path.join(a.cikti, "taxonomy"))
    print(f"\n  total {toplam_ok} headers, {toplam_yz} sequences written, {toplam_sz} without a lineage")
    print(f"  taxonomy: {n} nodes")
    if toplam_yz == 0:
        print(u'ERROR: not one sequence could be written, so the build is pointless. Stopping.')
        sys.exit(1)
    print(f"  library: {hedef}")
    print(f"  The next step is kraken2-build --build (kraken_tool.sh does it)")

if __name__ == "__main__":
    main()
