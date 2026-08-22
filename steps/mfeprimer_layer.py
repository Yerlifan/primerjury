#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mfeprimer_layer.py
THE SECOND, INDEPENDENT MEASUREMENT OF EXTERNAL DATABASE SPECIFICITY.

external_databases.py uses blastn: it searches for each primer separately, then
matches the hits itself and decides whether a product forms.
This script uses mfeprimer's spec subcommand: mfeprimer predicts the amplicon
directly from a k-mer index with its own thermodynamic model.
The two methods are independent of one another; the principle from the meeting
decision that no decision be left to a single piece of code is met in external
specificity only when both of them run.

The result is given in three columns:
  blast_urun     the off target product count external_databases.py found
  mfe_urun       the off target amplicon count mfeprimer found
  uyum           iki_olcum_uyustu | ayrisan_olcum | tek_olcum

IMPORTANT: mfeprimer does not apply the rule from the meeting decision that the
last two bases match exactly. It was measured: a primer whose last base at the 3'
end has been changed gives the same number of amplicons as the intact primer. That
is why the mfeprimer results are not an elimination criterion but a second point of
view.

Agreement DOES NOT mean that the absolute numbers are equal; the threshold and
model differences of the two methods separate the numbers unavoidably. Agreement is
the two of them giving THE SAME DECISION: is it zero or is it not zero. A
divergence is not passed over silently, it is reported in a separate column and in
the summary row.

Usage:
  python3 mfeprimer_layer.py --final final_primers --db REFERENCE_DB       --mfe tools/mfeprimer --out final_primers/mfeprimer.tsv

"""
import argparse, csv, collections, os, re, shutil, subprocess, sys, tempfile

# THE SAME class to database mapping as external_databases.py; if the two drift
# apart the comparison becomes meaningless.
SINIF_DB = {
    "A1": ["archaea.16S.fna"],
    "A2": ["archaea.16S.fna"],
    "B":  ["bacteria.16S.fna"],
    "F1": ["fungi.ITS.fna", "fungi.28SrRNA.fna"],
    "F2": ["fungi.ITS.fna", "fungi.28SrRNA.fna"],
}

# The same broad set as external_databases.py; the two measurements must see the
# same databases, otherwise the diverging measurement rows show a difference of
# coverage rather than a difference of method. The list has to be EXACTLY THE SAME
# as the one in external_databases.py; after ROD was measured there to be eukaryote
# only (60320/60320 Eukaryota, 0 Bacteria, 0 Archaea) the hand made selection by
# domain was dropped and every class came to see every rDNA database.
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
    p.add_argument("--db", required=True, help='the reference database '
                                               'directory')
    p.add_argument("--mfe", required=True, help='the mfeprimer executable')
    p.add_argument("--blast", default=None,
                   help="output of the external-databases step; if given, the two measurements are compared")
    p.add_argument("--out", required=True)
    p.add_argument("--prod-min", type=int, default=50)
    p.add_argument("--prod-max", type=int, default=400)
    p.add_argument("--tm-min", type=float, default=30.0,
                   help="mfeprimer amplikon Tm alt siniri")
    p.add_argument("--mismatch", type=int, default=3,
                   help='the total mismatch bound from the design rules')
    p.add_argument("--mis-end", type=int, default=3,
                   help="mfeprimer's mismatch window. 3 is fast and "
                        'selective, while 9, the mfeprimer default, is far '
                        'too loose and very slow. This value DOES NOT APPLY '
                        "the 3' end rule of the design decisions; see the "
                        'note at the head of the script.')
    p.add_argument("--wide", action="store_true",
                   help="also scan the same wide database set as the external-databases step")
    p.add_argument("--cpu", type=int, default=4)
    p.add_argument("--timeout", type=int, default=3600)
    return p.parse_args()


def indeks_eksik(fna):
    """The helper files mfeprimer spec needs. Looking only at .primerqc.bin
    is not enough: without .fai and .json mfeprimer writes no valid db found and
    RETURNS EXIT CODE 0, that is, it silently produces nothing. It returns the
    missing file names; an empty list shows that everything is in place.

    """
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
        sys.exit(u'mfeprimer was not found: %s' % a.mfe)
    if not os.access(a.mfe, os.X_OK):
        try:
            os.chmod(a.mfe, 0o755)
        except OSError:
            sys.exit(u'mfeprimer is not executable: %s\n   try chmod +x \'%s\'' % (a.mfe, a.mfe))
    tsv = os.path.join(a.final, "final_primers.tsv")
    if not os.path.exists(tsv):
        sys.exit(u'not found: %s' % tsv)
    rows = [r for r in csv.DictReader(open(tsv, encoding="utf-8"), delimiter="\t")
            if r.get("ozgulluk_durum") == "GECTI"]
    if not rows:
        sys.exit(u'no candidate passed')
    print(u'pairs to test: %d' % len(rows))

    # what external_databases.py produced: (target, class, forward, reverse, database) -> product count
    blast = {}
    byol = a.blast or os.path.join(a.final, "dis_veritabani.tsv")
    if os.path.exists(byol):
        for x in csv.DictReader(open(byol, encoding="utf-8"), delimiter="\t"):
            k = (x["hedef"], x["sinif"], x["ileri_dizi"], x["geri_dizi"],
                 x["veritabani"])
            blast[k] = blast.get(k, 0) + int(x.get("hedef_disi_urun", 0) or 0)
        print(u'read the external-databases output: %s (%d records)' % (byol, len(blast)))
    else:
        print(u'WARNING: the external-databases output is missing, no comparison will be possible: %s' % byol)

    calisma = tempfile.mkdtemp(prefix="mfe_")
    sonuc = []
    hata = 0
    # one call per class: in mfeprimer's tsv input every row is a pair
    sinif_cift = collections.defaultdict(list)
    for i, r in enumerate(rows):
        sinif_cift[r["sinif"]].append((i, r))

    for sinif, ciftler in sorted(sinif_cift.items()):
        dblist = list(SINIF_DB.get(sinif, []))
        if a.wide:
            dblist += SINIF_DB_GENIS.get(sinif, [])
        for dbad in dblist:
            fna = os.path.join(a.db, dbad)
            if not os.path.exists(fna):
                print(u'   no database, skipped: %s' % fna)
                continue
            eks = indeks_eksik(fna)
            if eks:
                print(u'   the mfeprimer index is missing (%s): %s'
                      % (", ".join(eks), dbad))
                print(u'   To build it: %s index -i %s -c %d'
                      % (a.mfe, fna, a.cpu))
                hata += 1
                continue
            girdi = os.path.join(calisma, "%s_%s.tsv" % (sinif, dbad))
            with open(girdi, "w", encoding="utf-8") as fh:
                for i, r in ciftler:
                    fh.write("p%d\t%s\t%s\n" % (i, r["ileri_dizi"], r["geri_dizi"]))
            cikti = girdi + ".mfe.tsv"
            # mfeprimer's OWN model runs; the binding rule from the meeting
            # decision IS NOT IMITATED here. Had it been imitated the two
            # measurements would not be independent and the second one would
            # mean nothing.
            #
            # MEASURED (archaea.16S, 2026-08-01): when it is run with
            # mfeprimer --misEnd 3, a primer whose LAST BASE AT THE 3' END
            # has been changed gives THE SAME number of amplicons as the
            # intact primer (323 against 323). So mfeprimer does not apply
            # the rule that the last two bases must match exactly; it builds
            # its own k-mer seed somewhere else. --misEnd 9 (the default)
            # gives more than 30 000 amplicons for the same primer and is
            # too slow to use in practice. That is why the mfeprimer counts
            # ARE NOT AN ELIMINATION CRITERION; if it finds a product blastn
            # did not, that is a warning that the pair needs a separate look
            # in the laboratory.
            cmd = [a.mfe, "spec", "-i", girdi, "-o", cikti, "-d", fna,
                   "-s", str(a.prod_min), "-S", str(a.prod_max),
                   "-t", str(a.tm_min), "-c", str(a.cpu),
                   "--misMatch", str(a.mismatch),
                   "--misStart", "1", "--misEnd", str(a.mis_end)]
            print(u'   mfeprimer %-4s x %-22s (%d pairs)'
                  % (sinif, dbad, len(ciftler)))
            try:
                p = subprocess.run(cmd, capture_output=True, text=True,
                                   timeout=a.timeout)
            except subprocess.TimeoutExpired:
                print(u'      TIMEOUT (%d s), skipped' % a.timeout)
                continue
            ciktisi = (p.stdout or "") + (p.stderr or "")
            if p.returncode != 0:
                print(u'      ERROR: %s' % ciktisi.strip()[:220])
                hata += 1
                continue
            # On some errors mfeprimer RETURNS EXIT CODE 0 and only writes
            # to the screen. Counting the exit code as the only criterion
            # produces a false cleanliness that reads as there being no off
            # target product.
            if "no valid db" in ciktisi.lower() or "error" in ciktisi.lower():
                print(u'      ERROR (exit code 0 but there is a message): %s'
                      % ciktisi.strip()[:220])
                hata += 1
                continue
            # mfeprimer ciktisini oku: amplikon satirlarini cift basina say
            say = collections.Counter()
            ornek = collections.defaultdict(list)
            # mfeprimer writes two files: <out> is a human readable report
            # and <out>.spec.tsv is for the machine. The second one is read;
            # the first one's format changes between versions and is not
            # reliable for counting.
            okunan = None
            for c in (cikti + ".spec.tsv", cikti + ".tsv", cikti):
                if os.path.exists(c) and c.endswith(".spec.tsv"):
                    okunan = c
                    break
            if not okunan:
                print(u'      no .spec.tsv was produced, so this database was NOT MEASURED. mfeprimer said: %s' % (ciktisi.strip()[:160] or "yok"))
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
                    # fpName and rpName must belong to the same pair; if they do not,
                    # this amplicon is made of the primers of two different pairs and
                    # it is not written to that pair.
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
        print(u'no database could be scanned; no output was written')
        # Coming out empty must not be read silently as clean; an empty file WITH
        # its headers is written so that stale output does not survive.
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
    print(u'records: %d' % len(sonuc))
    for k in ("iki_olcum_uyustu", "ayrisan_olcum", "tek_olcum"):
        if say_uyum.get(k):
            print("   %-20s %d" % (k, say_uyum[k]))
    ayrisan = [x for x in sonuc if x["uyum"] == "ayrisan_olcum"]
    if ayrisan:
        print(u'\nDIVERGING MEASUREMENTS (one found a product, the other did not):')
        for x in ayrisan[:20]:
            print("   %-30s %-3s %-18s blast=%-5s mfe=%-5s"
                  % (x["hedef"][:29], x["sinif"], x["veritabani"][:17],
                     x["blast_urun"], x["mfe_urun"]))
        if len(ayrisan) > 20:
            print(u'   ... %d more records' % (len(ayrisan) - 20))
        print(u'These rows do not count as eliminated. They are where the two methods disagree, and they need separate attention in the laboratory.')
    temiz = sorted(set((x["hedef"], x["sinif"]) for x in sonuc
                       if x["mfe_urun"] == 0))
    print(u'\ntarget classes with no off-target amplicon according to mfeprimer: %d'
          % len(temiz))
    if hata:
        print(u'\nCAUTION: %d databases could not be measured. No second measurement was made for them, and the gap is recorded in the result file' % hata)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
