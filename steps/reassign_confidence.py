#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reassign_confidence.py
APPLIES THE KRAKEN2 CONFIDENCE THRESHOLD FROM THE OUTPUT FILES, WITH NO DATABASE.

Why: under its default setting (--confidence 0) Kraken2 does not abstain. If the real
organism is not in the database, a read is labelled with the highest scoring sibling
leaf, whichever discriminating k-mers survived. Measured: in sample A2-4 the reads of
four Methanosarcina bins go to the same references and not one of them prefers the
species it was assigned to. The fix is to raise the confidence threshold, and this
script applies it from the output files themselves so that the whole classification
does not have to be run again.

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
    p.add_argument("--report-pattern", default="*_kraken2.report")
    p.add_argument("--output-pattern", default="*_output")
    p.add_argument("--max-reads", type=int, default=0,
                   help='0 means all of them; lower it for a trial run')
    p.add_argument("--scan", default=None,
                   help='comma separated thresholds; when it is given only '
                        'the scan runs and no report is written. For example: '
                        '0,0.005,0.01,0.02,0.05')
    p.add_argument("--scan-reads", type=int, default=20000,
                   help="reads per file in scan mode")
    return p.parse_args()


def agac_kur(kraken_kok, desen):
    """Builds the taxon tree from all the reports.
    Returns: (parent, rank, name, report_count)

    """
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
    ebeveyn, rutbe, ad, nrep = agac_kur(a.kraken, a.report_pattern)
    if not ebeveyn:
        sys.exit(u'the report was not found: %s' % a.kraken)
    print(u'reports: %d, taxa in the tree: %d' % (nrep, len(ebeveyn)))

    # atalar zinciri, bir kez hesaplanir
    ata_onbellek = {}

    def atalar(t):
        """The chain from t to the root, t included. For an unknown taxon it returns an empty
        list; such a hit is counted in no clade.

        """
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

    ciktilar = sorted(glob.glob(os.path.join(a.kraken, "*", a.output_pattern))
                      + glob.glob(os.path.join(a.kraken, a.output_pattern)))
    ciktilar = [c for c in ciktilar if not c.endswith(".report")]
    if not ciktilar:
        sys.exit(u'no output file was found (pattern: %s)' % a.output_pattern)
    print("output dosyasi: %d" % len(ciktilar))
    os.makedirs(a.out, exist_ok=True)

    if a.scan:
        esikler = [float(x) for x in a.scan.split(",") if x.strip()]
        return tarama_yap(a, ciktilar, esikler, atalar, rutbe)

    ozet, karsilastirma = [], []
    for yol in ciktilar:
        taban = re.sub(r"_output$", "", os.path.basename(yol))
        dogrudan = Counter()          # the new assignment -> the read count
        eski_yeni = Counter()         # (eski_rutbe, yeni_rutbe) -> sayi
        n = tasinan = sinifsiz_olan = 0
        vurus_top = vurus_eksik = 0
        with open(yol, errors="replace") as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 5:
                    continue
                n += 1
                if a.max_reads and n > a.max_reads:
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
                        continue          # no match; it enters Q but not the clade
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
                # on the path from the root to the leaf of the old assignment, the MOST
                # SPECIFIC node passing the threshold is chosen. atalar() returns them
                # ordered from leaf to root.
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

        # the Kraken2 report format
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
        print(u'   %-34s reads=%8d moved=%7d (%5.2f%%) unclassified=%7d off_tree_hits=%.3f%%'
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
    print(u'total reads: %d, rank changed: %d (%.2f%%)'
          % (tt, ts, 100.0 * ts / tt if tt else 0))
    print(u'\nThese reports can be handed straight to Bracken. No database and no raw fastq are needed; the calculation comes from Kraken2\'s own k-mer output.')


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
    print(u'\nTHRESHOLD SCAN (at most %d reads per file)' % a.scan_reads)
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
                    if n > a.scan_reads:
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
    print(u'Choose the threshold here: if the genus(G) column is rising while the unclassified column is still low, that threshold is suitable. Pass the value you chose with --confidence and run without --scan.')
    return 0


if __name__ == "__main__":
    main()
