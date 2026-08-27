#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
abundance_rank.py
GIVES THE ABUNDANCE TABLE AT THE RANK THE DATA SUPPORTS.

Why Bracken was not run
-----------------------
The first plan was to feed the confidence corrected Kraken2 reports to Bracken and
produce abundance at genus level. It was measured first and the plan was dropped.

Bracken distributes the reads sitting at upper ranks downward, using the k-mer
distribution of the database. That rests on the assumption that the real organism IS
IN the database. In this sample that assumption cannot be made:

  After a confidence threshold of 0,02, of 2 314 887 reads
     able to reach genus level    :  400 190  (17.3 percent)
     unable to reach genus level  : 1 723 684  (77.7 percent)
     unclassified                 :  115 343  (5.0 percent)

Group by group the table is sharper still:
     Archaea  (barcode01-08) : 72 to 96 percent of the reads reach genus level
     Bacteria (barcode17-20) : only 8 to 9 percent
     Fungi F1 (13-16)        : 0.6 to 2.4 percent
     Fungi F2 (09-12)        : 0.02 to 0.3 percent

So genus level abundance is sound on the archaeal side and not on the bacterial and
fungal sides. Distributing 99.7 percent of the fungal reads to genera by the
database's priors would be manufacturing numbers. On top of that it was separately
measured that the dominant fungus in this sample (Microascaceae, close to
Petriella) and most of the bacterial lineages are not in the database.

What this script does
---------------------
It gives each sample's abundance at the narrowest rank THE DATA SUPPORTS in that
sample. The rank is not chosen by hand: whichever rank a given proportion of the
reads (50 percent by default) can settle at or below is the rank chosen. The reads
that cannot settle are not hidden; they stay on a separate row as "left at a higher
rank".

Usage:
  python3 abundance_rank.py --kraken kraken_c0.02 --out bolluk_c0.02

"""
import argparse, csv, glob, os, re, sys
from collections import defaultdict

# The Kraken2 rank codes, from broad to narrow. Sub-ranks (P1, C2 and the like) are
# folded into the main rank; the intermediate ranks change with the taxonomy version
# and make the comparison fragile.
RUTBE_SIRA = ["D", "K", "P", "C", "O", "F", "G", "S"]
RUTBE_ADI = {"D": "domain", "K": "kingdom", "P": "phylum", "C": "class",
             "O": "order", "F": "family", "G": "genus", "S": "species"}

GRUP_ARALIK = [("A1", 1, "Archaea, short amplicon"),
               ("A2", 5, "Archaea, long amplicon"),
               ("F2", 9, "Fungi, long amplicon"),
               ("F1", 13, "Fungi, short amplicon"),
               ("B", 17, "Bacteria")]
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
                   help="report directory with the confidence correction applied")
    p.add_argument("--out", required=True)
    p.add_argument("--pattern", default="*_kraken2.report")
    p.add_argument("--coverage", type=float, default=0.50,
                   help='the classified fraction required before a rank can '
                        'be chosen: at least this fraction of the reads has '
                        'to sit at that rank or below it')
    p.add_argument("--top", type=int, default=15, help="taxon shown in the table")
    return p.parse_args()


def rapor_oku(yol):
    """Reads a Kraken2 report. The parent chain is built from the indentation depth as
    well; that is needed to tell nested nodes at the same rank apart.
    Returns: (direct{taxid:(rank, name, count)}, clade{taxid:count},
              parent{taxid:taxid}, total)

    """
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
    dosyalar = sorted(glob.glob(os.path.join(a.kraken, a.pattern))
                      + glob.glob(os.path.join(a.kraken, "*", a.pattern)))
    if not dosyalar:
        sys.exit(u'the report was not found: %s' % a.kraken)
    os.makedirs(a.out, exist_ok=True)
    print(u'reports: %d' % len(dosyalar))

    kapsama_satir, bolluk_satir, ozet_satir = [], [], []
    for yol in dosyalar:
        m = re.search(r"barcode(\d+)", os.path.basename(yol))
        if not m:
            print(u'   SKIPPED, the barcode could not be resolved: %s' % os.path.basename(yol))
            continue
        bc = int(m.group(1))
        if bc not in BARKOD_GRUP:
            print(u'   SKIPPED, not in the sampling map: barcode%02d' % bc)
            continue
        dogrudan, klan, ebeveyn, toplam = rapor_oku(yol)
        if toplam <= 0:
            print(u'   SKIPPED, empty report: %s' % os.path.basename(yol))
            continue
        sinifsiz = sum(d for r, _, d in dogrudan.values() if r == "U")
        siniflanan = toplam - sinifsiz

        # for each main rank: the reads assigned directly AT or BELOW that rank
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
            if siniflanan and kumulatif[r] / siniflanan >= a.coverage:
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

        # The abundance at the chosen rank: the CLADE counts of the taxa at that rank.
        #
        # CAREFUL: NESTED nodes can appear at the same main rank. In a Kraken2
        # report the real phylum is coded "P", the subphylum "P1", the one under
        # it "P2", and all of them fold into the same main rank. Because the
        # upper node's clade already holds the lower node's clade, counting both
        # is double counting. Measured: before the fix the percentages of
        # barcode10 summed to 368.54 percent.
        #
        # The fix: only the nodes with NO ANCESTOR at the same main rank are taken.
        # That set is disjoint by definition and its clades do not overlap.
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
        # the disjoint set check: the clade total cannot exceed the number of reads that
        # settle at or below that rank
        _t = sum(z[1] for z in sec)
        if _t > kumulatif[secilen] + 1:
            print(u'   WARNING barcode%02d: at rank %s the clade total (%d) exceeds the placed reads (%d), which can be double counting'
                  % (bc, secilen, _t, kumulatif[secilen]))
        # o rutbenin ustunde kalanlar
        ust_kalan = siniflanan - kumulatif[secilen]
        for ad, n, tx in sec[:a.top]:
            bolluk_satir.append(dict(
                barkod="barcode%02d" % bc, grup=BARKOD_GRUP[bc],
                yil=BARKOD_YIL[bc], rutbe=secilen, takson=ad, taxid=tx,
                okuma=n, yuzde=round(100.0 * n / toplam, 4)))
        if ust_kalan > 0:
            bolluk_satir.append(dict(
                barkod="barcode%02d" % bc, grup=BARKOD_GRUP[bc],
                yil=BARKOD_YIL[bc], rutbe=secilen,
                takson="[reads that cannot reach %s level]" % RUTBE_ADI[secilen],
                taxid="", okuma=ust_kalan,
                yuzde=round(100.0 * ust_kalan / toplam, 4)))
        if sinifsiz > 0:
            bolluk_satir.append(dict(
                barkod="barcode%02d" % bc, grup=BARKOD_GRUP[bc],
                yil=BARKOD_YIL[bc], rutbe=secilen, takson="[unclassified]",
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
        print(u'   barcode%02d %-3s total=%8d  chosen rank=%s (%s)  placed=%.1f%%  genus=%.1f%%  species=%.1f%%'
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
    print("\nwritten: %s" % a.out)

    # a summary group by group
    print("\nTHE RANK THE DATA SUPPORTS, PER GROUP")
    print("%-4s %-24s %-14s %10s %10s"
          % ("group", "description", "rank chosen",
             "genus share", "species share"))
    for g, b0, ac in GRUP_ARALIK:
        ilgili = [x for x in ozet_satir if x["grup"] == g]
        if not ilgili:
            continue
        rutler = sorted(set(x["secilen_rutbe_adi"] for x in ilgili))
        co = sum(x["cins_orani"] for x in ilgili) / len(ilgili)
        to = sum(x["tur_orani"] for x in ilgili) / len(ilgili)
        print("%-4s %-24s %-14s %9.1f%% %9.1f%%"
              % (g, ac, ", ".join(rutler), co, to))
    print(u'\nThe genus fraction shows what percentage of the reads in that group could be placed at genus rank or narrower. A low fraction means the organisms in that group are not represented in the Kraken2 database, and spreading their counts over a genus would be inventing data.')


if __name__ == "__main__":
    main()
