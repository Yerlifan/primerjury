#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
KRAKEN2 RECLASSIFICATION SUMMARY

Reads the <taxid>.out files produced by rerun_kraken.sh and works out, for each
bin, what the reads actually look like.

WHAT IS BEING ASKED
Every bin was labelled as belonging to one taxon. Do most of its reads really go
to that taxon, or to something else? Seven eukaryotes were already known to be
mislabelled; this summary puts every bin through the same criterion rather than
only those seven, because a mistake in a place nobody looks still counts.

WHAT IT DOES NOT SAY
Kraken2 is a measurement too and it has its own blind spots: an organism that is
not represented in the database gets labelled as its nearest relative. So the
output is a second opinion, not a decision. For a bin that disagrees with its
label, the decision is made by hand.

The status values written to the CSV (uyusuyor / UYUSMUYOR / KARAR YOK) are DATA:
kraken_tool.sh reads column 9 and compares it against "uyusuyor". They are kept
as they are and the screen text is printed in English beside them.

Run:
  python3 tools/kraken_summary.py --job <kraken_yeniden directory>
  python3 tools/kraken_summary.py --selftest
"""
import argparse, csv, glob, os, sys
from collections import Counter

# The three states, as stored, with the wording used on screen.
UYUSUYOR = "uyusuyor"
UYUSMUYOR = "UYUSMUYOR"
KARAR_YOK = "KARAR YOK"
KARISIK = "KARISIK"
EKRAN = {UYUSUYOR: "agrees", UYUSMUYOR: "DISAGREES",
         KARAR_YOK: "NO DECISION", KARISIK: "MIXED"}


def read_names(toolkit):
    """
    taxid -> organism name.

    WHY THE MODULE IS NOT IMPORTED
    The names live in the ISIMLER dictionary inside blast_ispcr.py. An earlier
    version loaded that file as a module, but it wants numpy and calls sys.exit
    when numpy is missing. sys.exit raises SystemExit, which "except Exception"
    does NOT catch, so the summary would crash because of numpy. That happens in
    exactly the situation rerun_kraken.sh creates: a micromamba environment that
    has kraken2 but no numpy.

    The summary needs no numpy at all, only the name list. So the file is
    parsed with ast straight out of its source text and never executed.
    Zero dependencies.
    """
    import ast
    adaylar = [os.path.join(toolkit, "blast_ispcr.py") if toolkit else "",
               os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "blast_ispcr.py")]
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


def species_map(rapor):
    """
    Builds taxid -> (species_taxid, species_name) from tum.report.

    WHY IT IS NEEDED
    kraken2 can leave the LCA at a node BELOW the species. The report carries rank
    codes: S species, S1 and S2 strain. The reads of the Methanothrix soehngenii
    bin, for example, are labelled "Methanothrix soehngenii GP6 (taxid 990316)",
    and that taxid DIFFERS from the species taxid 2223. A criterion that compares
    taxids literally calls this "DISAGREES", when it is the same species.

    For the same reason two STRAINS are enough to make a bin look MIXED: in the
    Methanosarcina barkeri bin, "barkeri 3" comes out at 72 percent and "barkeri
    str. Wiesmoor" at 24 percent, and both are the same species. So the counting
    is aggregated at SPECIES level.

    Report columns: percent, clade_reads, node_reads, rank, taxid, indented name.
    The tree comes in as indentation, two spaces per level.
    """
    if not rapor or not os.path.exists(rapor):
        return {}
    harita = {}
    yigin = []
    with open(rapor, errors="replace") as fh:
        for satir in fh:
            a = satir.rstrip("\n").split("\t")
            if len(a) < 6:
                continue
            rank, tx, ham_ad = a[3].strip(), a[4].strip(), a[5]
            derinlik = (len(ham_ad) - len(ham_ad.lstrip(" "))) // 2
            ad = ham_ad.strip()
            while yigin and yigin[-1][0] >= derinlik:
                yigin.pop()
            yigin.append((derinlik, rank, tx, ad))
            # If the node is itself a species, that is the answer; otherwise the
            # nearest species ancestor. With no species ancestor (genus, family,
            # kingdom and other higher nodes) the node stays as itself, and it is
            # visible that it did.
            tur = next(((t, n) for _, r, t, n in reversed(yigin) if r == "S"), None)
            harita[tx] = tur if tur else (tx, ad)
    return harita


def read_rows(yol):
    """
    Kraken2 output is tab separated: C/U, read name, identity, length, LCA map.
    With --use-names the third field has the form "Name (taxid NNN)".
    Returns: [(was_classified, identity_text), ...]
    """
    out = []
    with open(yol, errors="replace") as fh:
        for satir in fh:
            a = satir.rstrip("\n").split("\t")
            if len(a) < 3:
                continue
            out.append((a[0] == "C", a[2].strip()))
    return out


def summarise_bin(kayitlar, en_fazla=4, tur_map=None):
    """
    Returns: (total, classified, [(identity, count, fraction, species_taxid), ...])

    The counting is aggregated at SPECIES level, so strains of the same species
    merge into one entry. Otherwise two strains of one species make a bin look
    both "DISAGREES" and "MIXED".

    Unclassified reads DO go into the denominator of the fractions; shrinking the
    denominator makes a bin look more consistent than it is.
    """
    toplam = len(kayitlar)
    if not toplam:
        return 0, 0, []
    tur_map = tur_map or {}
    sayac = Counter()
    adlar = {}
    sinif = 0
    for c, k in kayitlar:
        if not c:
            continue
        sinif += 1
        tx = pull_taxid(k)
        ttx, tad = tur_map.get(tx, (tx, k))
        anahtar = ttx or k
        sayac[anahtar] += 1
        adlar.setdefault(anahtar, tad if ttx else k)
    ilk = [(adlar[a], n, n / toplam, a) for a, n in sayac.most_common(en_fazla)]
    return toplam, sinif, ilk


def split(tum_out, klasor):
    """
    Splits tum.out into <taxid>.out files by the tx<taxid>_ prefix in the read
    name. Files are written in one go, so many files are never open at once.
    """
    from collections import defaultdict
    kova = defaultdict(list)
    adsiz = []
    with open(tum_out, errors="replace") as fh:
        for satir in fh:
            a = satir.split("\t")
            if len(a) < 3:
                continue
            ad = a[1]
            if ad.startswith("tx") and "_" in ad:
                tx = ad[2:ad.index("_")]
                if tx.isdigit():
                    kova[tx].append(satir)
                    continue
            adsiz.append(satir)
    for tx, satirlar in kova.items():
        with open(os.path.join(klasor, "%s.out" % tx), "w") as fh:
            fh.writelines(satirlar)
    if adsiz:
        with open(os.path.join(klasor, "ADSIZ.out"), "w") as fh:
            fh.writelines(adsiz)
        print("WARNING: the source taxon of %d reads could not be resolved, ADSIZ.out"
              % len(adsiz))
    print("split into %d taxon files" % len(kova))
    return kova


def pull_taxid(kimlik):
    """'Methanosarcina mazei (taxid 2209)' -> '2209'. Empty if there is none."""
    if "(taxid" in kimlik:
        p = kimlik.rsplit("(taxid", 1)[1].strip().rstrip(")").strip()
        return p if p.isdigit() else ""
    return ""


# ------------------------------------------------------------------ selftest
def selftest():
    print("=" * 72)
    print("KRAKEN SUMMARY, A TEST WITH KNOWN ANSWERS")
    print("=" * 72)
    hata = 0

    def K(ad, bul, bek):
        nonlocal hata
        ok = bul == bek
        if not ok:
            hata += 1
        print("  %s  %-56s %s / %s" % ('PASS' if ok else 'FAIL', ad, bul, bek))

    K("the taxid is parsed", pull_taxid("Methanosarcina mazei (taxid 2209)"), "2209")
    K("no taxid gives an empty string", pull_taxid("unclassified"), "")
    K("a name with parentheses does not break it",
      pull_taxid("Candidatus X (sp.) (taxid 44)"), "44")

    # 10 reads: 6 correct, 2 something else, 2 unclassified
    kay = ([(True, "A (taxid 1)")] * 6 + [(True, "B (taxid 2)")] * 2
           + [(False, "unclassified")] * 2)
    top, sin, ilk = summarise_bin(kay)
    K("total reads are counted", top, 10)
    K("classified reads are counted", sin, 8)
    K("the most frequent identity is right", ilk[0][0], "A (taxid 1)")
    # 6/10 = 0.6, NOT 6/8. A smaller denominator makes a bin look more consistent.
    K("the fraction divides by all reads, not by the classified ones",
      round(ilk[0][2], 3), 0.6)
    K("the second identity is reported too", ilk[1][0], "B (taxid 2)")

    # If blast_ispcr.py sits beside it, the names must be readable WITHOUT running
    # it. This item guarantees the summary still works where numpy is not installed.
    _ad = read_names("")
    if _ad:
        K("names are read without executing the module",
          _ad.get("2209"), "Methanosarcina mazei")
        K("the name count is plausible", len(_ad) > 30, True)
    else:
        print("  SKIPPED  blast_ispcr.py is not beside it, the name test was not run")

    # SPECIES LEVEL AGGREGATION. These items exist because of a real bug: in the
    # Methanosarcina barkeri bin the strains "barkeri 3" and "barkeri str.
    # Wiesmoor" were counted separately, and the bin looked both DISAGREES and
    # MIXED. Both are the same species.
    import tempfile as _tf
    with _tf.TemporaryDirectory() as _d:
        _r = os.path.join(_d, "tum.report")
        open(_r, "w").write(
            "  4.35\t3779\t49\tS\t2208\t                    Methanosarcina barkeri\n"
            "  3.22\t2796\t2796\tS1\t1434107\t                      Methanosarcina barkeri 3\n"
            "  1.07\t932\t932\tS1\t1434109\t                      Methanosarcina barkeri str. Wiesmoor\n"
            "  3.53\t3064\t0\tS\t2223\t                    Methanothrix soehngenii\n"
            "  3.53\t3064\t3064\tS1\t990316\t                      Methanothrix soehngenii GP6\n")
        _t = species_map(_r)
        K("a strain is attached to its species ancestor",
          _t.get("1434107"), ("2208", "Methanosarcina barkeri"))
        K("a species node is attached to itself",
          _t.get("2208"), ("2208", "Methanosarcina barkeri"))
        K("the second species is right too",
          _t.get("990316"), ("2223", "Methanothrix soehngenii"))
        _kay = ([(True, "Methanosarcina barkeri 3 (taxid 1434107)")] * 72
                + [(True, "Methanosarcina barkeri str. Wiesmoor (taxid 1434109)")] * 24
                + [(True, "Methanocorpusculum labreanum Z (taxid 410358)")] * 4)
        _t2, _s2, _i2 = summarise_bin(_kay, tur_map=_t)
        K("two strains of one species merge into ONE entry", _i2[0][1], 96)
        K("the merged entry is shown under the species name",
          _i2[0][0], "Methanosarcina barkeri")
        K("the species taxid comes back, comparable with the label", _i2[0][3], "2208")
        K("the second entry is below 20 percent, so the bin is not MIXED",
          _i2[1][2] < 0.20, True)

    top0, sin0, ilk0 = summarise_bin([])
    K("an empty bin does not crash", (top0, sin0, ilk0), (0, 0, []))

    print("=" * 72)
    print("TEST PASSED" if hata == 0 else "TEST FAILED, %d items" % hata)
    print("=" * 72)
    return 0 if hata == 0 else 1


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--job", "--is", dest="is_klasor", default="",
                    help="the kraken_yeniden directory holding tum.out and tum.report")
    ap.add_argument("--toolkit", default="",
                    help="directory holding blast_ispcr.py (default: beside this script)")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(selftest())
    if not a.is_klasor:
        ap.error("--job is required")

    if quiet_selftest() != 0:
        print("THE SELF TEST FAILED, stopping")
        sys.exit(2)

    ISIMLER = read_names(a.toolkit)
    if not ISIMLER:
        print("WARNING: the organism names could not be read, taxids will be shown instead")

    def isim(tx):
        return ISIMLER.get(tx, "taxid %s" % tx)

    # If tum.out exists it is split by taxon first. The split happens here because
    # the number of files awk can hold open at once is limited in some versions,
    # and past that limit lines are lost without an error.
    tum = os.path.join(a.is_klasor, "tum.out")
    if os.path.exists(tum):
        split(tum, a.is_klasor)
    TUR = species_map(os.path.join(a.is_klasor, "tum.report"))
    if TUR:
        print("species map built: %d nodes. Counts are aggregated at SPECIES level, "
              "so strains of one species merge into a single entry." % len(TUR))
    else:
        print("WARNING: there is no tum.report, so strain level labels could not be "
              "lifted to species level. Strains of one species can wrongly look "
              "DISAGREES and MIXED.")

    dosyalar = sorted(glob.glob(os.path.join(a.is_klasor, "*.out")))
    dosyalar = [d for d in dosyalar if os.path.basename(d)[:-4].isdigit()]
    if not dosyalar:
        print("there is no <taxid>.out file at all")
        sys.exit(1)

    # The expected taxon set is read from kaynak_sayim.tsv, and anything missing
    # is not skipped silently. "A taxon was never scanned" and "it was scanned and
    # came out clean" must never arrive at the same conclusion.
    bekleniyor = set()
    ks = os.path.join(a.is_klasor, "kaynak_sayim.tsv")
    if os.path.exists(ks):
        with open(ks, errors="replace") as fh:
            for satir in fh:
                p = satir.split("\t")
                if p and p[0].strip().isdigit():
                    bekleniyor.add(p[0].strip())
    bulunan = {os.path.basename(d)[:-4] for d in dosyalar}
    eksik = sorted(bekleniyor - bulunan)

    satirlar = []
    print("\n%-38s%7s%7s  most frequent identity" % ('bin', 'reads', 'class'))
    print("-" * 100)
    for d in dosyalar:
        tx = os.path.basename(d)[:-4]
        kay = read_rows(d)
        toplam, sinif, ilk = summarise_bin(kay, tur_map=TUR)
        birinci_ad, birinci_n, birinci_oran, birinci_tx = (
            ilk[0] if ilk else ("", 0, 0.0, ""))
        # The species taxid of the bin's own label. If the label is a strain taxid
        # it is lifted to its species, otherwise the two sides would be compared
        # at different levels.
        kutu_tur = TUR.get(tx, (tx, ""))[0]

        # Three states: the same taxid as the label, a different taxid, or no
        # decision possible. "No decision" is kept apart and never put in the same
        # sentence as a disagreement.
        if not ilk or birinci_oran < 0.20:
            durum = KARAR_YOK
        elif birinci_tx == kutu_tur:
            durum = UYUSUYOR
        else:
            durum = UYUSMUYOR

        ikinci_oran = ilk[1][2] if len(ilk) > 1 else 0.0
        # If a second identity also holds a large share, the bin is mixed.
        karisim = KARISIK if ikinci_oran >= 0.20 else ""

        print("%-38s%7d%7d  %-45s %6.1f%%  %s %s"
              % (isim(tx)[:37], toplam, sinif, birinci_ad[:44],
                 birinci_oran * 100, EKRAN[durum],
                 EKRAN[KARISIK] if karisim else ""))
        for k, n, o, _ in ilk[1:]:
            print("%-52s%-45s %6.1f%%" % ('', k[:44], o * 100))

        satirlar.append(dict(
            taxid=tx, etiket=isim(tx), okuma=toplam, siniflandirilan=sinif,
            kraken_1=birinci_ad, kraken_1_oran=round(birinci_oran, 4),
            kraken_2=(ilk[1][0] if len(ilk) > 1 else ""),
            kraken_2_oran=round(ikinci_oran, 4),
            durum=durum, karisim=karisim))

    yol = os.path.join(a.is_klasor, "kraken_ozet.csv")
    with open(yol, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(satirlar[0]))
        w.writeheader()
        w.writerows(satirlar)

    uyusmaz = [s for s in satirlar if s["durum"] == UYUSMUYOR]
    kararsiz = [s for s in satirlar if s["durum"] == KARAR_YOK]
    karisik = [s for s in satirlar if s["karisim"]]
    print("\n" + "=" * 100)
    print("%d bins scanned" % len(satirlar))
    print("  disagreeing with the label : %d" % len(uyusmaz))
    for s in uyusmaz:
        print("      %-36s -> %s (%.1f%%)"
              % (s['etiket'], s['kraken_1'], s['kraken_1_oran'] * 100))
    print("  no decision possible       : %d   (no identity reached 20%% of the reads)"
          % len(kararsiz))
    for s in kararsiz:
        print("      %-36s most frequent: %s (%.1f%%)"
              % (s['etiket'], s['kraken_1'] or 'none', s['kraken_1_oran'] * 100))
    print("  looking mixed              : %d   (the second identity is also above 20%%)"
          % len(karisik))
    for s in karisik:
        print("      %-36s %s + %s" % (s['etiket'], s['kraken_1'], s['kraken_2']))
    if eksik:
        print("\n  %d TAXA WERE NEVER SCANNED, and for those bins there is NO result:"
              % len(eksik))
        for tx in eksik:
            print("      %s (taxid %s)" % (isim(tx), tx))
    else:
        print("\n  every taxon seen in the source was scanned (%d)"
              % (len(bekleniyor) or len(satirlar)))
    print("\nwritten: %s" % yol)
    print("\nNOTE: kraken2 is a measurement too. An organism that is not represented")
    print("in the database gets labelled as its nearest relative. This table is a")
    print("second opinion, not a decision.")


def quiet_selftest():
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        return selftest()


if __name__ == "__main__":
    main()
