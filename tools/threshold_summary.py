#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THE SUMMARY OF THE CONFIDENCE THRESHOLD SCAN AND THE COMPARISON OF TWO DATABASES

It reads the esik_<C>.report and esik_<C>.out files kraken_tool.sh produces.

WHAT IS BEING ASKED
Kraken2 assigns a read only if a certain fraction of its k-mers go to the same
clade. At threshold 0 a single k-mer is enough for an assignment. As the threshold
rises, weak assignments drop out. How many drop is Kraken's own measure of how weak
those assignments were to begin with.

THE REAL QUESTION (when two databases are given)
On the old database, when the threshold rises, do the assignments that collapse stay
standing on the new one? If they do, the problem was the database's COVERAGE: the
organism was not in the old database, Kraken assigned the read to the nearest
relative, and that assignment was weak.

"""
import argparse, csv, glob, os, re, sys
from collections import Counter, defaultdict

# The domain definitions. taxid is for NCBI; the name is a fallback for the custom
# database with a synthetic taxonomy. Fungi sit UNDER eukaryota, so the two are
# separate columns but nested.
ALANLAR = [
    ("arke",    "2157",  {"archaea"}),
    ("bakteri", "2",     {"bacteria"}),
    ("okaryot", "2759",  {"eukaryota", "eukarya"}),
    ("mantar",  "4751",  {"fungi"}),
    ("bitki",   "33090", {"viridiplantae", "plantae"}),
    ("virus",   "10239", {"viruses"}),
]

# ------------------------------------------------------------------ reading
def rapor_oku(yol):
    """
    Reads a Kraken2 report.
    The columns: percent, clade_reads, node_reads, rank, taxid, the indented name.
    Returns: (nodes, parent, total_reads, unclassified)

    The indentation gives a tree, two spaces per level. The parent map comes out of it,
    and measurement 2's summing is done with that map.

    """
    dugumler = []
    ebeveyn = {}
    yigin = []
    siniflandirilmayan = 0
    kok_klad = 0
    if not os.path.exists(yol):
        return [], {}, 0, 0
    with open(yol, errors="replace") as fh:
        for satir in fh:
            a = satir.rstrip("\n").split("\t")
            if len(a) < 6:
                continue
            try:
                klad = int(a[1]); dugum = int(a[2])
            except ValueError:
                continue
            rank, tx, ham = a[3].strip(), a[4].strip(), a[5]
            derinlik = (len(ham) - len(ham.lstrip(" "))) // 2
            ad = ham.strip()
            if rank == "U":
                siniflandirilmayan = klad
                continue
            if tx == "1":
                kok_klad = klad
            while yigin and yigin[-1][0] >= derinlik:
                yigin.pop()
            ebeveyn[tx] = yigin[-1][1] if yigin else ""
            yigin.append((derinlik, tx))
            dugumler.append(dict(klad=klad, dugum=dugum, rank=rank,
                                 taxid=tx, ad=ad, derinlik=derinlik))
    toplam = siniflandirilmayan + kok_klad
    return dugumler, ebeveyn, toplam, siniflandirilmayan

def alan_dugumu(dugumler, taxid, adlar):
    """Looks for the domain's node by taxid first and by name if that fails."""
    for d in dugumler:
        if d["taxid"] == taxid:
            return d
    for d in dugumler:
        if d["ad"].strip().lower() in adlar:
            return d
    return None

def alan_raporundan(dugumler, toplam):
    """MEASUREMENT 1: the domain percentages from the clade read counts in the report."""
    sonuc = {}
    for ad, tx, adlar in ALANLAR:
        d = alan_dugumu(dugumler, tx, adlar)
        sonuc[ad] = d["klad"] if d else 0
    return sonuc

def taxid_cek(kimlik):
    """'Methanosarcina mazei (taxid 2209)' -> '2209'. Bulamazsa bos."""
    if "(taxid" in kimlik:
        p = kimlik.rsplit("(taxid", 1)[1].strip().rstrip(")").strip()
        return p if p.isdigit() else ""
    return ""

def out_oku(yol):
    """
    The Kraken2 output: C/U, the read name, the identity, the length, the k-mer map.
    Returns: (counter_taxid, total, unclassified, bin_counter)
    bin_counter: the source taxid (the @tx<taxid>_ prefix) -> Counter(assigned taxid)

    """
    sayac = Counter()
    kutu = defaultdict(Counter)
    toplam = 0
    sinifsiz = 0
    if not os.path.exists(yol):
        return sayac, 0, 0, kutu
    with open(yol, errors="replace") as fh:
        for satir in fh:
            a = satir.rstrip("\n").split("\t")
            if len(a) < 3:
                continue
            toplam += 1
            ad = a[1]
            kaynak = ""
            if ad.startswith("tx") and "_" in ad:
                k = ad[2:ad.index("_")]
                if k.isdigit():
                    kaynak = k
            if a[0] != "C":
                sinifsiz += 1
                if kaynak:
                    kutu[kaynak]["U"] += 1
                continue
            tx = taxid_cek(a[2]) or a[2].strip()
            sayac[tx] += 1
            if kaynak:
                kutu[kaynak][tx] += 1
    return sayac, toplam, sinifsiz, kutu

def atalar(tx, ebeveyn, sinir=200):
    """tx and all its ancestors. If a cycle forms it stops at the bound rather than running forever."""
    out = []
    g = tx
    n = 0
    while g and n < sinir:
        out.append(g)
        g = ebeveyn.get(g, "")
        n += 1
    return out

def alan_outtan(sayac, ebeveyn, dugumler):
    """
    MEASUREMENT 2: it sums the per read assignments into domains over the tree.
    It is not a route wholly independent of the report (the tree comes from the report),
    but the COUNTING is independent; that catches clade count parsing faults.

    """
    hedef = {}
    for ad, tx, adlar in ALANLAR:
        d = alan_dugumu(dugumler, tx, adlar)
        if d:
            hedef[d["taxid"]] = ad
    sonuc = {ad: 0 for ad, _, _ in ALANLAR}
    onbellek = {}
    for tx, n in sayac.items():
        if tx not in onbellek:
            bulunan = [hedef[a] for a in atalar(tx, ebeveyn) if a in hedef]
            onbellek[tx] = bulunan
        for ad in onbellek[tx]:
            sonuc[ad] += n
    return sonuc

# ------------------------------------------------------------------ tarama
def esik_listesi(klasor):
    'Sorts the threshold report files into numeric order.'
    out = []
    for y in glob.glob(os.path.join(klasor, "esik_*.report")):
        m = re.match(r"esik_(.+)\.report$", os.path.basename(y))
        if not m:
            continue
        try:
            c = float(m.group(1))
        except ValueError:
            continue
        out.append((c, m.group(1), y))
    return sorted(out)

def tarama_oku(klasor):
    """
    Reads one scan directory.
    Returns: (rows, bins)
      rows: [dict(esik, toplam, sinifsiz, sinifsiz_oran, alanlar, alanlar2, ayrilik)]
      bins: {threshold_text: {source_taxid: Counter(assigned)}}

    """
    satirlar = []
    kutular = {}
    for c, metin, rap in esik_listesi(klasor):
        dugumler, ebeveyn, toplam, sinifsiz = rapor_oku(rap)
        a1 = alan_raporundan(dugumler, toplam)
        out_yol = os.path.join(klasor, f"esik_{metin}.out")
        sayac, toplam2, sinifsiz2, kutu = out_oku(out_yol)
        a2 = alan_outtan(sayac, ebeveyn, dugumler) if sayac else {}
        if kutu:
            kutular[metin] = kutu
        # The two measurements are compared. The report's total is used as the
        # denominator; if there is an out file the totals are expected to be equal too.
        ayrilik = []
        if toplam2 and toplam and toplam2 != toplam:
            ayrilik.append(f"read count report {toplam} / out {toplam2}")
        if a2:
            for ad, _, _ in ALANLAR:
                if a1.get(ad, 0) != a2.get(ad, 0):
                    ayrilik.append(f"{ad} rapor {a1.get(ad,0)} / out {a2.get(ad,0)}")
        satirlar.append(dict(
            esik=c, esik_metni=metin, toplam=toplam, sinifsiz=sinifsiz,
            alanlar=a1, alanlar2=a2, ayrilik="; ".join(ayrilik)))
    return satirlar, kutular

def yuzde(n, toplam):
    return (100.0 * n / toplam) if toplam else 0.0

# ------------------------------------------------------------------ yazim
def egri_metni(satirlar, baslik):
    g = []
    g.append("=" * 78)
    g.append(baslik)
    g.append("=" * 78)
    g.append(f"{'esik':>6} {'okuma':>8} {'sinifsiz':>9} " +
             " ".join(f"{a:>9}" for a, _, _ in ALANLAR))
    g.append("-" * 78)
    for s in satirlar:
        t = s["toplam"]
        g.append(f"{s['esik']:>6} {t:>8} {yuzde(s['sinifsiz'], t):>8.2f}% " +
                 " ".join(f"{yuzde(s['alanlar'].get(a,0), t):>8.2f}%" for a, _, _ in ALANLAR))
    g.append("")
    g.append(u'fungi sit UNDER eukaryota, so the two columns are nested.')
    ayr = [s for s in satirlar if s["ayrilik"]]
    if ayr:
        g.append("")
        g.append(u'THERE IS A DISAGREEMENT. The two measurements did not give the same number, so do not trust them:')
        for s in ayr:
            g.append(f"  esik {s['esik']}: {s['ayrilik']}")
    else:
        g.append(u'the two independent measurements gave the same result at every threshold.')
    return "\n".join(g)

def kutu_hakim(kutu_sayac):
    """A bin's most frequent assignment and its fraction, with the unclassified reads in the denominator."""
    toplam = sum(kutu_sayac.values())
    if not toplam:
        return "", 0.0
    en = [(tx, n) for tx, n in kutu_sayac.most_common() if tx != "U"]
    if not en:
        return "U", 0.0
    return en[0][0], en[0][1] / toplam

def cokme_esigi(kutular, kaynak, alt=0.20):
    """
    The threshold at which a bin's dominant assignment first falls below 20 percent.
    If it never falls, None is returned. That is exactly what "a collapsing assignment"
    means.

    """
    for metin in sorted(kutular, key=lambda m: float(m)):
        k = kutular[metin].get(kaynak)
        if not k:
            continue
        _, oran = kutu_hakim(k)
        if oran < alt:
            return float(metin)
    return None

def ayakta_kalma(kutular_a, kutular_b, ad_a, ad_b):
    """
    THE REAL QUESTION. Do the assignments that collapse on the old database stay
    standing on the new one?
    Returns: (rows, summary)

    """
    kaynaklar = set()
    for k in (kutular_a, kutular_b):
        for m in k.values():
            kaynaklar |= set(m.keys())
    satirlar = []
    for kay in sorted(kaynaklar, key=lambda x: (len(x), x)):
        ca = cokme_esigi(kutular_a, kay)
        cb = cokme_esigi(kutular_b, kay)
        if ca is None and cb is None:
            durum = "ikisinde de ayakta"
        elif ca is not None and cb is None:
            durum = "ESKIDE COKTU, YENIDE AYAKTA"
        elif ca is None and cb is not None:
            durum = "eskide ayakta, yenide coktu"
        elif cb > ca:
            durum = "yenide daha dayanikli"
        elif cb < ca:
            durum = 'more fragile in the new one'
        else:
            durum = 'both collapsed at the same threshold'
        satirlar.append(dict(kaynak=kay, cokme_a=ca, cokme_b=cb, durum=durum))
    ozet = Counter(s["durum"] for s in satirlar)
    return satirlar, ozet

# ------------------------------------------------------------------ selftest
def selftest():
    print("=" * 72)
    print(u'THRESHOLD SUMMARY, A TEST WITH KNOWN ANSWERS')
    print("=" * 72)
    hata = 0
    def K(ad, bul, bek):
        nonlocal hata
        ok = bul == bek
        if not ok:
            hata += 1
        print(f"  {'PASS' if ok else 'FAIL'}  {ad:<58} {bul} / {bek}")

    import tempfile
    # A report built by hand whose answer was worked out on paper.
    # 1000 reads: 100 unclassified, 900 root. Archaea 400, bacteria 300,
    # eukaryota 200, of which 150 fungi.
    rap = ("  10.00\t100\t100\tU\t0\tunclassified\n"
           "  90.00\t900\t10\tR\t1\troot\n"
           "  89.00\t890\t0\tR1\t131567\t  cellular organisms\n"
           "  40.00\t400\t50\tD\t2157\t    Archaea\n"
           "  35.00\t350\t350\tS\t2209\t      Methanosarcina mazei\n"
           "  30.00\t300\t100\tD\t2\t    Bacteria\n"
           "  20.00\t200\t200\tS\t1642647\t      Proteiniphilum saccharofermentans\n"
           "  20.00\t200\t50\tD\t2759\t    Eukaryota\n"
           "  15.00\t150\t0\tK\t4751\t      Fungi\n"
           "  15.00\t150\t150\tS\t101201\t        Trichoderma asperellum\n")
    with tempfile.TemporaryDirectory() as d:
        ry = os.path.join(d, "esik_0.report")
        open(ry, "w").write(rap)
        dug, eb, top, sifsiz = rapor_oku(ry)
        K(u'the total reads are the unclassified plus the root clade', top, 1000)
        K(u'the unclassified are read', sifsiz, 100)
        K("agac girintiden kurulur, mantar okaryotun altinda", eb.get("4751"), "2759")
        K("tur, alanin altinda dogru baglanir", eb.get("2209"), "2157")
        a1 = alan_raporundan(dug, top)
        K(u'the archaea clade count', a1["arke"], 400)
        K(u'the bacteria clade count', a1["bakteri"], 300)
        K(u'the eukaryota clade count', a1["okaryot"], 200)
        K(u'the fungi clade count', a1["mantar"], 150)
        K("bulunmayan alan sifirdir, cokmez", a1["bitki"], 0)
        K("yuzde hesabi", round(yuzde(a1["arke"], top), 2), 40.0)

        # An out file giving the same numbers for MEASUREMENT 2.
        oy = os.path.join(d, "esik_0.out")
        satir = []
        satir += ["C\ttx2209_r%d\tMethanosarcina mazei (taxid 2209)\t1500\t\n" % i for i in range(350)]
        satir += ["C\ttx2157_r%d\tArchaea (taxid 2157)\t1500\t\n" % i for i in range(50)]
        satir += ["C\ttx1642647_r%d\tP. saccharofermentans (taxid 1642647)\t1500\t\n" % i for i in range(200)]
        satir += ["C\ttx2_r%d\tBacteria (taxid 2)\t1500\t\n" % i for i in range(100)]
        satir += ["C\ttx101201_r%d\tTrichoderma asperellum (taxid 101201)\t1500\t\n" % i for i in range(150)]
        satir += ["C\ttx2759_r%d\tEukaryota (taxid 2759)\t1500\t\n" % i for i in range(50)]
        satir += ["U\ttx101201_u%d\tunclassified (taxid 0)\t1500\t\n" % i for i in range(100)]
        open(oy, "w").writelines(satir)
        sayac, top2, sifsiz2, kutu = out_oku(oy)
        K(u'out the total reads', top2, 1000)
        K("out siniflandirilamayan", sifsiz2, 100)
        a2 = alan_outtan(sayac, eb, dug)
        K(u'measurement 2 archaea, the same as measurement 1', a2["arke"], 400)
        K(u'measurement 2 bacteria, the same as measurement 1', a2["bakteri"], 300)
        K(u'measurement 2 fungi, the same as measurement 1', a2["mantar"], 150)
        K(u'measurement 2 eukaryota holds the fungi too', a2["okaryot"], 200)
        K(u'the bin is taken from the read name', sorted(kutu.keys())[0:2], ["101201", "1642647"])

        satirlar, kutular = tarama_oku(d)
        K(u'the scan reads a single threshold', len(satirlar), 1)
        K(u'no disagreement, the two measurements match', satirlar[0]["ayrilik"], "")

        # AYRILIK GERCEKTEN YAKALANIYOR MU. Bu madde, olcutun kor olmadigini
        # sinar: bilerek bozulmus bir out dosyasi AYRILIK vermeli.
        open(oy, "a").write("C\ttx2_x1\tBacteria (taxid 2)\t1500\t\n")
        s2, _ = tarama_oku(d)
        K(u'a corrupted file is caught as a DISAGREEMENT', s2[0]["ayrilik"] != "", True)

    # The collapse and survival of a bin.
    A = {"0":   {"K1": Counter({"9": 90, "U": 10})},
         "0.1": {"K1": Counter({"9": 10, "U": 90})}}
    B = {"0":   {"K1": Counter({"7": 95, "U": 5})},
         "0.1": {"K1": Counter({"7": 92, "U": 8})}}
    K("hakim atama bulunur", kutu_hakim(A["0"]["K1"])[0], "9")
    K("hakim oran paydaya sinifsizi katar", round(kutu_hakim(A["0"]["K1"])[1], 2), 0.9)
    K("cokme esigi bulunur", cokme_esigi(A, "K1"), 0.1)
    K(u'a bin that does not collapse returns None', cokme_esigi(B, "K1"), None)
    sat, ozet = ayakta_kalma(A, B, "eski", "yeni")
    K("eskide coken, yenide ayakta kalan yakalanir",
      sat[0]["durum"], "ESKIDE COKTU, YENIDE AYAKTA")

    K("dongulu agac sonsuza gitmez", len(atalar("a", {"a": "b", "b": "a"})), 200)
    K(u'an empty report does not crash it', rapor_oku(u'/no/such/file')[2], 0)

    print("=" * 72)
    print("THE TEST PASSED" if hata == 0 else f"THE TEST FAILED, {hata} items")
    print("=" * 72)
    return 0 if hata == 0 else 1

def selftest_sessiz():
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        return selftest()

# ------------------------------------------------------------------ ana
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--job", dest="is1", default="")
    ap.add_argument("--name", default="database 1")
    ap.add_argument("--job2", default="")
    ap.add_argument("--name2", default="database 2")
    ap.add_argument("--root", dest="kok", default="")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(selftest())
    if not a.is1:
        ap.error("--job is required")
    if selftest_sessiz() != 0:
        print("THE TEST FAILED, stopped (project rule 2)")
        sys.exit(2)

    s1, k1 = tarama_oku(a.is1)
    if not s1:
        print(f"ERROR: {a.is1} holds no esik_<C>.report file.")
        print(u'  The scan has to be run first:  bash kraken_tool.sh threshold-old')
        sys.exit(1)

    metin = egri_metni(s1, f"THE CONFIDENCE THRESHOLD CURVE, {a.ad}")
    print(metin)
    yaz_csv(os.path.join(a.is1, "esik_egrisi.csv"), s1, a.ad)
    open(os.path.join(a.is1, "esik_egrisi.txt"), "w", encoding="utf-8").write(metin + "\n")
    print(f"\nwritten: {os.path.join(a.is1, 'esik_egrisi.csv')}")

    if not a.is2:
        return
    s2, k2 = tarama_oku(a.is2)
    if not s2:
        print(f"\nWARNING: {a.is2} holds no scan, so the two databases were not compared.")
        print(u'  To run it:  bash kraken_tool.sh threshold-new')
        return
    m2 = egri_metni(s2, f"THE CONFIDENCE THRESHOLD CURVE, {a.ad2}")
    print("\n" + m2)
    yaz_csv(os.path.join(a.is2, "esik_egrisi.csv"), s2, a.ad2)
    open(os.path.join(a.is2, "esik_egrisi.txt"), "w", encoding="utf-8").write(m2 + "\n")

    birlesik = yan_yana(s1, s2, a.ad, a.ad2)
    print("\n" + birlesik)
    sat, ozet = ayakta_kalma(k1, k2, a.ad, a.ad2)
    kalma = kalma_metni(sat, ozet, a.ad, a.ad2, a.kok)
    print("\n" + kalma)
    hedef = os.path.join(a.is1, "esik_iki_veritabani.txt")
    open(hedef, "w", encoding="utf-8").write(birlesik + "\n\n" + kalma + "\n")
    with open(os.path.join(a.is1, "esik_ayakta_kalma.csv"), "w",
              newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["kaynak_taxid", f"cokme_esigi_{a.ad}", f"cokme_esigi_{a.ad2}", "durum"])
        for s in sat:
            w.writerow([s["kaynak"], s["cokme_a"] if s["cokme_a"] is not None else "cokmedi",
                        s["cokme_b"] if s["cokme_b"] is not None else "cokmedi", s["durum"]])
    print(f"\nwritten: {hedef}")

def yaz_csv(yol, satirlar, ad):
    with open(yol, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        basliklar = ["veritabani", "esik", "toplam_okuma", "siniflandirilamayan",
                     "siniflandirilamayan_yuzde"]
        for x, _, _ in ALANLAR:
            basliklar += [f"{x}_okuma", f"{x}_yuzde", f"{x}_yuzde_olcum2"]
        basliklar.append("ayrilik")
        w.writerow(basliklar)
        for s in satirlar:
            t = s["toplam"]
            sat = [ad, s["esik"], t, s["sinifsiz"], round(yuzde(s["sinifsiz"], t), 3)]
            for x, _, _ in ALANLAR:
                sat += [s["alanlar"].get(x, 0),
                        round(yuzde(s["alanlar"].get(x, 0), t), 3),
                        round(yuzde(s["alanlar2"].get(x, 0), t), 3) if s["alanlar2"] else ""]
            sat.append(s["ayrilik"])
            w.writerow(sat)

def yan_yana(s1, s2, ad1, ad2):
    g = ["=" * 96,
         f"THE TWO DATABASES SIDE BY SIDE   ({ad1}  and  {ad2})",
         "=" * 96,
         f"{'esik':>6} | " + " | ".join(
             f"{x:>10}" for x in ["unclassified", "archaea", "bacteria",
                                  "fungi"]),
         f"{'':>6} | " + " | ".join(
             f"{ad1[:4]:>4}/{ad2[:5]:>5}" for _ in range(4)),
         "-" * 96]
    h2 = {s["esik"]: s for s in s2}
    for s in s1:
        o = h2.get(s["esik"])
        hucre = []
        for anahtar in ["sinifsiz", "arke", "bakteri", "mantar"]:
            v1 = (yuzde(s["sinifsiz"], s["toplam"]) if anahtar == "sinifsiz"
                  else yuzde(s["alanlar"].get(anahtar, 0), s["toplam"]))
            if o:
                v2 = (yuzde(o["sinifsiz"], o["toplam"]) if anahtar == "sinifsiz"
                      else yuzde(o["alanlar"].get(anahtar, 0), o["toplam"]))
                hucre.append(f"{v1:>4.1f}/{v2:>5.1f}")
            else:
                hucre.append(f"{v1:>4.1f}/{'none':>5}")
        g.append(f"{s['esik']:>6} | " + " | ".join(f"{h:>10}" for h in hucre))
    eksik = [s["esik"] for s in s1 if s["esik"] not in h2]
    if eksik:
        g.append("")
        g.append(f"WARNING: these thresholds are MISSING on the {ad2} side, so they were "
                 f"not compared: {eksik}")
        g.append("  A missing threshold must not be read as a compared threshold.")
    return "\n".join(g)

def kalma_metni(sat, ozet, ad1, ad2, kok):
    isimler = isimleri_oku(kok)
    g = ["=" * 96,
         "THE REAL QUESTION: do the assignments that collapse on the old one stay standing on the new one",
         "=" * 96,
         "",
         f"The threshold at which a bin 'collapses' is the one where its dominant",
         f"assignment first falls below %20 of the reads. The left column is "
         f"{ad1}, the right {ad2}.",
         ""]
    for d, n in ozet.most_common():
        g.append(f"  {n:>3} bins   {d}")
    g.append("")
    g.append(f"{'bin':<40}{ad1[:12]:>12}{ad2[:12]:>12}   state")
    g.append("-" * 96)
    for s in sat:
        a = "cokmedi" if s["cokme_a"] is None else str(s["cokme_a"])
        b = "cokmedi" if s["cokme_b"] is None else str(s["cokme_b"])
        ad = isimler.get(s["kaynak"], f"taxid {s['kaynak']}")
        g.append(f"{ad[:39]:<40}{a:>12}{b:>12}   {s['durum']}")
    g.append("")
    n_dogrulayan = ozet.get("ESKIDE COKTU, YENIDE AYAKTA", 0) + ozet.get("yenide daha dayanikli", 0)
    n_toplam = sum(ozet.values())
    g.append("READING")
    if n_toplam == 0:
        g.append("  There is no bin to compare.")
    elif n_dogrulayan > n_toplam / 2:
        g.append(f"  In {n_dogrulayan} of {n_toplam} bins the assignments are more "
                 f"resilient on the new database.")
        g.append("  That supports the diagnosis that the problem was COVERAGE: the old")
        g.append("  database did not hold these organisms, Kraken assigned them to the")
        g.append("  nearest relative, and those assignments were weak. The new database")
        g.append("  holds the organism, so the assignment is strong.")
    else:
        g.append(f"  In only {n_dogrulayan} of {n_toplam} bins are the assignments more "
                 f"resilient on the new one.")
        g.append("  That DOES NOT SUPPORT the COVERAGE diagnosis. The weakness may come")
        g.append("  from the reads themselves rather than from the database coverage. "
                 "The diagnosis has to be reviewed.")
    return "\n".join(g)

def isimleri_oku(kok):
    """taxid -> ad. The same route as tools/kraken_summary.py, with no numpy dependency."""
    import ast
    adaylar = [os.path.join(kok, "tools", "blast_ispcr.py") if kok else "",
               os.path.join(os.path.dirname(os.path.abspath(__file__)), "blast_ispcr.py")]
    for yol in adaylar:
        if not yol or not os.path.exists(yol):
            continue
        try:
            agac = ast.parse(open(yol, encoding="utf-8", errors="replace").read())
        except SyntaxError:
            continue
        for d in agac.body:
            if not isinstance(d, ast.Assign):
                continue
            for h in d.targets:
                if isinstance(h, ast.Name) and h.id == "ISIMLER":
                    try:
                        return ast.literal_eval(d.value)
                    except (ValueError, SyntaxError):
                        return {}
    return {}

if __name__ == "__main__":
    main()
