#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
external_databases.py
Tests the primer pairs that passed verification against EXTERNAL reference
databases. It shows whether a product forms in relatives that are absent from the
sample but could be present in the environment.

The method: every primer is searched against the database with blastn (task
blastn-short). If two hits on the same reference sequence, on opposite strands and
with their 3' ends facing one another, meet within the product length range, that
counts as an off-target product. BLAST alone does not answer "will the primer
bind"; so every hit is additionally checked against the panel's own criteria.

"""
import argparse, csv, os, re, shutil, subprocess, sys, tempfile, collections

TAMLAYICI = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T",
         "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}

# Which class is tested against which database.
#
# THE NARROW SET (the default): NCBI RefSeq's curated 16S and ITS sets.
# These are made of type strains and representative sequences; uncultured
# environmental lineages are largely ABSENT. A significant part of an anaerobic
# digester community is exactly that group, so the narrow set on its own is not
# enough to say "there is no off-target product".
SINIF_DB = {
    "A1": ["archaea.16S.fna"],
    "A2": ["archaea.16S.fna"],
    "B":  ["bacteria.16S.fna"],
    "F1": ["fungi.ITS.fna", "fungi.28SrRNA.fna"],
    "F2": ["fungi.ITS.fna", "fungi.28SrRNA.fna"],
}

# THE WIDE SET (with --genis): the large databases that also hold environmental
# sequences. SILVA SSU/LSU NR99 and UNITE carry uncultured lineages that are not in
# RefSeq; ROD covers full rRNA operon variants and PR2 the eukaryotic SSU. The
# running time grows considerably, which is why it is not the default.
#
# NO PER-CLASS MANUAL SELECTION IS MADE. In the earlier version the wide set was
# written out by class by hand, and that led to A SILENT LOSS OF MEASUREMENT:
#   MEASURED (2026-08-01): all 60320 of the 60320 records in
#   ROD_v1.2_operon_variants.fasta are Eukaryota; it holds 0 Bacteria and 0
#   Archaea (counted from the lineage in the headers). ROD had been assigned to
#   classes A1/A2/B. The result: "no off-target product in ROD" was written for 71
#   archaeal and bacterial pairs. That is not evidence of specificity; the database
#   holds no sequence from that domain at all. The same mistake had a mirror image:
#   ROD's 9753 full fungal operons had never been scanned in the fungal classes.
#
# The new rule: EVERY CLASS SEES EVERY rDNA DATABASE. In a real PCR a primer meets
# not only the rDNA of its own domain but all the DNA in the environment, and an
# archaeal primer mis-binding in a bacterial 23S is exactly the kind of error we are
# looking for. There is no scientific reason to limit by domain, only a speed
# reason, and the speed reason is already met by the narrow set.
#
# SILVA_138.2_LSUParc.fasta is DELIBERATELY OUT: the Parc set also holds partial and
# low quality records, while LSURef_NR99 is the 99% dereplicated curated
# representative of the same coverage. NR99 rather than Parc is used on the SSU side
# too; to treat the two the same way, NR99 is taken for the LSU as well.
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

# ref_all.fna and ref_all2.fna ARE DELIBERATELY ABSENT.
#   MEASURED (2026-08-01, by comparing the identity sets in the .fai files):
#   ref_all2.fna = archaea.16S + bacteria.16S + fungi.ITS + fungi.18SrRNA
#                + fungi.28SrRNA, exactly 65358 records, zero difference in either
#                direction.
#   ref_all.fna  = archaea.16S + bacteria.16S + fungi.ITS, 48431 records.
# So both are subsets of the list above; adding them to the scan brings not one new
# sequence, it only scans the same records a second time.
# (Which is also why there is no need to build a BLAST index for those two files.)

SINIF_DB_GENIS = {
    s: [d for d in GENIS_ORTAK if d not in SINIF_DB[s]] for s in SINIF_DB
}


def rc(s):
    return s.translate(TAMLAYICI)[::-1]


def uyar(p, t):
    return bool(set(IUPAC.get(p, "")) & set(IUPAC.get(t, "")))


def baglanma_uygun(oligo, hedef, son_tam=2, son_pencere=5, son_mm=1, toplam_mm=3):
    """A panel decision: because extension starts from the 3' end, the last bases are critical."""
    if len(oligo) != len(hedef):
        return False, None
    mm = [i for i, (p, t) in enumerate(zip(oligo, hedef)) if not uyar(p, t)]
    n = len(oligo)
    if any(i >= n - son_tam for i in mm):
        return False, None
    if sum(1 for i in mm if i >= n - son_pencere) > son_mm:
        return False, None
    if len(mm) > toplam_mm:
        return False, None
    return True, len(mm)


def db_hazirla(fna, calisma):
    """Uses an existing BLAST index, or builds one in the working directory if there is none."""
    if os.path.exists(fna + ".nin") or os.path.exists(fna + ".00.nin"):
        return fna
    hedef = os.path.join(calisma, os.path.basename(fna))
    if not os.path.exists(hedef + ".nin"):
        os.symlink(os.path.abspath(fna), hedef) if not os.path.exists(hedef) else None
        r = subprocess.run(["makeblastdb", "-in", hedef, "-dbtype", "nucl"],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print("   makeblastdb basarisiz: %s" % r.stderr.strip()[:200])
            return None
    return hedef


KONS_ETIKET = re.compile(r"((?:A1|A2|B|F1|F2))-\d+_(\d+)")


def sinif_konsensuslari(kons_klasoru, sinif):
    """Verilen sinifin konsensus dosyalarindan dizileri doner."""
    diziler = []
    if not kons_klasoru or not os.path.isdir(kons_klasoru):
        return diziler
    for ad in sorted(os.listdir(kons_klasoru)):
        if not ad.endswith(".fasta"):
            continue
        m = KONS_ETIKET.match(ad)
        if not m or m.group(1) != sinif:
            continue
        ad_, dizi = None, []
        with open(os.path.join(kons_klasoru, ad), encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(">"):
                    if ad_ and dizi:
                        diziler.append((ad_, "".join(dizi)))
                    ad_, dizi = line[1:].strip().split()[0], []
                else:
                    dizi.append(line.strip())
        if ad_ and dizi:
            diziler.append((ad_, "".join(dizi)))
    return diziler


def kapsam_olc(db, diziler, calisma, etiket, threads, zaman_asimi):
    """Does this database hold this class's organisms AT ALL?

        Added after the ROD mistake. If a database holds not one sequence from the
        relevant domain, a result of 'no off-target product found' IS NOT evidence of
        specificity; it came out empty because there was nothing to measure. Those two
        situations have to be distinguishable in the output.

        The measure: the class's own consensus sequences are searched against the
        database with megablast and the LONGEST alignment is taken. The threshold is not
        invented, it comes from the data: the product looked for is at most prod_max
        bases long.

    """
    if not diziler:
        return ("KAPSAM_OLCULMEDI", 0, 0.0)
    sorgu = os.path.join(calisma, "kapsam_%s.fa" % etiket)
    with open(sorgu, "w") as fh:
        for ad, dizi in diziler:
            fh.write(">%s\n%s\n" % (ad, dizi))
    cmd = ["blastn", "-query", sorgu, "-db", db,
           "-outfmt", "6 length pident", "-max_target_seqs", "5",
           "-num_threads", str(threads)]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=zaman_asimi)
    except subprocess.TimeoutExpired:
        return ("KAPSAM_OLCULMEDI", 0, 0.0)
    if r.returncode != 0:
        return ("KAPSAM_OLCULMEDI", 0, 0.0)
    en_uzun, en_kimlik = 0, 0.0
    for line in r.stdout.splitlines():
        p = line.split("\t")
        if len(p) < 2:
            continue
        try:
            uz, kim = int(p[0]), float(p[1])
        except ValueError:
            continue
        if uz > en_uzun:
            en_uzun, en_kimlik = uz, kim
    return ("", en_uzun, en_kimlik)


# --- telling the target's own taxon from a foreign taxon ----------------
#
# Measured on 2026-08-01: some of the records giving the HIGHEST numbers in the wide
# scan's "off-target product" column are in fact THE TARGET ITSELF.
#   Asetoklastik_metanojenler x archaea.16S = 306 products; example hits
#     NR_104707.1 Methanothrix soehngenii GP6
#     NR_028242.1 Methanothrix soehngenii Opfikon
#   Methanothrix_soehngenii_turu x SILVA = 308 products; every hit is
#     Methanosaetaceae;Methanothrix
# Those are not a specificity error, they are the primer doing its job: the database
# holds the sequences of the target taxon too. Against that, some records with LOWER
# counts are real errors:
#   Nitrosocosmicus_AOA x SILVA = 1119; the hits are Halobacteria,
#     Methanoperedenaceae, Thermoplasmata, Cenarchaeum
#   Petrimonas_cinsi x SILVA = 707; the hits are Clostridium, Bacteroides
# Ranking by the raw count was sending the user to fix the wrong primers. So every
# product is split, by the taxon of the reference it formed in, into ITS OWN TAXON
# and A FOREIGN TAXON.
#
# The taxon names are not written by hand; they come from two data sources:
#   1) the taxid list in hedefler.tsv plus taxid_adlari.tsv  (THE DECLARED name)
#   2) the olculen_kimlik column of hedef_kimlik.tsv         (THE MEASURED name)
# The two can diverge, and they do: the measured identity of the Trichoderma_cinsi
# target is Petriella musispora. Both count as "its own taxon"; which of them
# matched is visible in the output.

KUCUK_TOKEN = 4
ATLA_TOKEN = {"uncultured", "candidatus", "bacterium", "archaeon", "sp",
              "strain", "clone", "isolate", "unidentified", "environmental",
              "samples", "incertae", "sedis", "type", "material", "partial",
              "complete", "sequence", "ribosomal", "gene", "rrna", "genes",
              "fungal", "endophyte", "voucher", "culture", "enrichment"}


def _tokenlar(metin):
    """Reduces a header of any format to a common set of words.
        The SILVA (;), ROD (|;_), UNITE (k__/p__ and _), PR2 (|_) and RefSeq (space)
        formats all go through the same process.

    """
    return {t for t in re.split(r"[^A-Za-z]+", metin.lower())
            if len(t) >= KUCUK_TOKEN and t not in ATLA_TOKEN}


def _ad_cinsi(ad):
    """'Ca. Nitrosocosmicus hydrocola' -> nitrosocosmicus
       'uncultured Acetobacteroides sp.' -> acetobacteroides"""
    for kelime in re.split(r"[^A-Za-z]+", ad):
        k = kelime.lower()
        if len(k) >= KUCUK_TOKEN and k not in ATLA_TOKEN:
            return k
    return ""


def hedef_taksonlari(hedefler_tsv, adlar_tsv, kimlik_tsv):
    """{target: {'beyan': set, 'olculen': set, 'evrensel': bool}}"""
    adlar = {}
    if adlar_tsv and os.path.exists(adlar_tsv):
        for line in open(adlar_tsv, encoding="utf-8"):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0].strip().isdigit():
                adlar[p[0].strip()] = p[1].strip()
    out = {}
    if hedefler_tsv and os.path.exists(hedefler_tsv):
        for line in open(hedefler_tsv, encoding="utf-8"):
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or p[0] == "karar":
                continue
            ad, tidler = p[1], [t.strip() for t in p[3].split(",") if t.strip()]
            # '*A', '*B', '*F' evrensel hedeflerin isaretidir: bu hedeflerde
            # "yabanci takson" kavrami yoktur, cok sayida taksonu birden
            # cogaltmak zaten amactir.
            evrensel = any(t.startswith("*") for t in tidler)
            beyan = set()
            for t in tidler:
                c = _ad_cinsi(adlar.get(t, ""))
                if c:
                    beyan.add(c)
            out[ad] = {"beyan": beyan, "olculen": set(), "evrensel": evrensel}
    if kimlik_tsv and os.path.exists(kimlik_tsv):
        for r in csv.DictReader(open(kimlik_tsv, encoding="utf-8"),
                                delimiter="\t"):
            ad = r.get("hedef")
            if not ad:
                continue
            out.setdefault(ad, {"beyan": set(), "olculen": set(),
                                "evrensel": False})
            for parca in re.split(r"[;,]", r.get("olculen_kimlik") or ""):
                parca = re.sub(r"\(.*?\)", " ", parca)
                c = _ad_cinsi(parca)
                if c:
                    out[ad]["olculen"].add(c)
    return out


# --- the DISTANCE of a foreign hit --------------------------------------
#
# Measured after the second round, 2026-08-01: "a foreign taxon" is not a sufficient
# measure on its own either, because some of the targets are not A SINGLE TAXON but
# A FUNCTIONAL GROUP. Examples, from the real output:
#   Hidrojenotrofik_metanojenler -> the foreign hits are Methanobacterium
#     alcaliphilum and Methanosphaera stadtmanae. Neither is in the declared taxid
#     list, but both are hydrogenotrophic methanogens, and catching that group is
#     exactly what the target is for.
#   Nitrosocosmicus_AOA -> the foreign hits are Nitrosotalea and Nitrosopumilus.
#     Both are ammonia oxidising archaea, that is, AOA.
#   Zoopagomycota_mantari -> the foreign hit is Piptocephalis moniliformis, which is
#     itself a Zoopagomycota.
# Against that:
#   Petrimonas_cinsi -> Flavobacterium, Phocaeicola vulgatus
#   Trichoderma_cinsi -> Calonectria, Acremonium, Trichothecium, Aniptodera
# Those really are distant taxa.
#
# Making the distinction with a hand written "functional group" table would be
# exactly the thing to avoid. Instead the distance is MEASURED FROM THE DATA: the
# database headers carry the lineage (SILVA, ROD, PR2, UNITE). The target's own
# lineage is derived from the common prefix of the references that gave a product IN
# ITS OWN TAXON, and the depth each foreign hit shares with that lineage is measured.
# A hit sharing everything except the last two ranks counts as NEAR, the rest as FAR.
# The priority ordering is by the FAR count.

def _soyagaci(baslik):
    """Extracts the ordered lineage fields from a header. On formats that carry no
        lineage (RefSeq, for instance) it returns an empty list.

    """
    if not baslik:
        return []
    govde = baslik.split(None, 1)
    # SILVA: 'KIMLIK Archaea;Halobacteriota;...'
    if len(govde) > 1 and ";" in govde[1]:
        alanlar = govde[1].split(";")
    elif "|" in baslik:
        # ROD: 'GCA|kaynak|Eukaryota;...;Tur|size=1'
        # PR2: 'KIMLIK|18S_rRNA|nucleus||Eukaryota|TSAR|...'
        # UNITE: 'UDB|k__Fungi;p__...;s__Tur|SH...'
        parcalar = baslik.split("|")
        soy = [p for p in parcalar if ";" in p]
        alanlar = soy[0].split(";") if soy else parcalar[3:]
    else:
        return []
    out = []
    for x in alanlar:
        x = re.sub(r"^[a-z]__", "", x.strip())
        x = re.sub(r"[^A-Za-z ]+", " ", x).strip().lower()
        if x and x not in ("incertae sedis", "unclassified"):
            out.append(x)
    return out


def _ortak_derinlik(a, b):
    n = 0
    for x, y in zip(a, b):
        if x != y:
            break
        n += 1
    return n


def _baskin_soy(soylar):
    """The DOMINANT (most frequent) lineage of the hits in the target's own taxon.

        THE COMMON PREFIX IS NOT TAKEN. Had it been, the threshold would shift with the
        number and the diversity of the own hits: with one hit the prefix is the full
        lineage (depth 7), with two diverse hits only the order level (depth 4), and the
        same foreign hit could count as FAR on one record and NEAR on another. The
        dominant lineage keeps the depth fixed. Since the own hits are by definition all
        of the same genus, those lineages are already nearly identical.

    """
    soylar = [tuple(s) for s in soylar if s]
    if not soylar:
        return []
    return list(collections.Counter(soylar).most_common(1)[0][0])


def basliklari_coz(fna, kimlikler, onbellek):
    """Urun olusturan referanslarin tam basligini fasta'dan tek gecisle alir.
    blastn ciktisi yalniz kisa kimligi verir; takson bilgisi baslikta durur."""
    d = onbellek.setdefault(fna, {})
    eksik = {k for k in kimlikler if k not in d}
    if not eksik or not os.path.exists(fna):
        return d
    with open(fna, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            kimlik = line[1:].split(None, 1)[0]
            if kimlik in eksik:
                d[kimlik] = line[1:].rstrip("\n")
                eksik.discard(kimlik)
                if not eksik:
                    break
    for k in eksik:
        d[k] = ""          # cozulemedi; yabanci sayilmaz, bilinmiyor sayilir
    return d


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--final", required=True, help="09'un output directory")
    p.add_argument("--db", required=True, help="REFERANS_DB directory")
    p.add_argument("--hedefler", default=None, help="targets.tsv")
    p.add_argument("--adlar", default=None, help="taxid_adlari.tsv")
    p.add_argument("--kimlik", default=None,
                   help="hedef_kimlik.tsv from the target-identity step (measured identity)")
    p.add_argument("--kons", default=None,
                   help="consensus directory; if given, every (class, database) "
                        "ikilisi icin KAPSAM DENETIMI yapilir")
    p.add_argument("--out", default=None)
    p.add_argument("--prod-min", type=int, default=50)
    p.add_argument("--prod-max", type=int, default=400)
    p.add_argument("--evalue", type=float, default=1000.0)
    p.add_argument("--max-hedef", type=int, default=5000)
    p.add_argument("--is-parcacigi", type=int, default=4)
    p.add_argument("--genis", action="store_true",
                   help="also the large databases that include environmental sequences "
                        "tarar (SILVA, UNITE, ROD, PR2). Uzun surer.")
    p.add_argument("--yalniz-genis", action="store_true",
                   help="only genis kumeyi tarar")
    p.add_argument("--zaman-asimi", type=int, default=14400,
                   help="veritabani basina saniye siniri")
    return p.parse_args()


def main():
    a = get_args()
    tsv = os.path.join(a.final, "primer_final.tsv")
    if not os.path.exists(tsv):
        sys.exit(u'not found: %s' % tsv)
    rows = [r for r in csv.DictReader(open(tsv, encoding="utf-8"), delimiter="\t")
            if r.get("ozgulluk_durum") == "GECTI"]
    if not rows:
        sys.exit(u'no candidate passed')
    print(u'pairs to test: %d' % len(rows))

    # sinif -> primer kumesi
    sinif_primer = collections.defaultdict(dict)   # sinif -> dizi -> ad
    for i, r in enumerate(rows):
        sinif_primer[r["sinif"]]["%s_F%d" % (r["hedef"][:24], i)] = r["ileri_dizi"]
        sinif_primer[r["sinif"]]["%s_R%d" % (r["hedef"][:24], i)] = r["geri_dizi"]

    calisma = tempfile.mkdtemp(prefix="blastdb_")
    sonuc = []
    atlanan = []
    gecersiz = []          # (sinif, db, en_uzun_hizalama, kimlik)
    kons_onbellek = {}
    baslik_onbellek = {}
    taksonlar = hedef_taksonlari(a.hedefler, a.adlar, a.kimlik)
    if taksonlar:
        ev = sum(1 for v in taksonlar.values() if v["evrensel"])
        ad = sum(1 for v in taksonlar.values()
                 if not v["evrensel"] and not (v["beyan"] | v["olculen"]))
        print(u'targets whose taxon name resolved: %d (universal: %d, unnamed: %d)'
              % (len(taksonlar), ev, ad))
        if ad:
            print(u'   WARNING: own versus foreign cannot be separated for unnamed targets; those products count as \'taxon unknown\'')
    else:
        print(u'WARNING: --targets/--names/--identity were not given, so own taxon and foreign taxon will NOT be separated and raw product counts are reported')
    for sinif, primerler in sorted(sinif_primer.items()):
        dblist = [] if a.yalniz_genis else list(SINIF_DB.get(sinif, []))
        if a.genis or a.yalniz_genis:
            dblist += SINIF_DB_GENIS.get(sinif, [])
        for dbad in dblist:
            fna = os.path.join(a.db, dbad)
            if not os.path.exists(fna):
                print(u'   no database, skipped: %s' % fna)
                atlanan.append((sinif, dbad, "dosya yok"))
                continue
            db = db_hazirla(fna, calisma)
            if not db:
                atlanan.append((sinif, dbad, "indeks kurulamadi"))
                continue
            # THE COVERAGE AUDIT, before the scan. See kapsam_olc().
            kap_durum, kap_uz, kap_kim = kapsam_olc(
                db, kons_onbellek.setdefault(
                    sinif, sinif_konsensuslari(a.kons, sinif)),
                calisma, "%s_%s" % (sinif, dbad), a.is_parcacigi, a.zaman_asimi)
            if not kap_durum:
                if kap_uz < a.prod_max:
                    kap_durum = "KAPSAM_YOK"
                    print(u'      NO COVERAGE: the longest region in %s resembling class %s is %d bp (%.1f%%), while the product sought is at most %d bp.' % (dbad, sinif, kap_uz, kap_kim, a.prod_max))
                    gecersiz.append((sinif, dbad, kap_uz, kap_kim))
                else:
                    kap_durum = "KAPSANIYOR"
            sorgu = os.path.join(calisma, "%s_%s.fa" % (sinif, dbad))
            with open(sorgu, "w") as fh:
                for ad, dizi in primerler.items():
                    fh.write(">%s\n%s\n" % (ad, dizi))
            cikti = sorgu + ".tsv"
            cmd = ["blastn", "-task", "blastn-short", "-query", sorgu, "-db", db,
                   "-outfmt", "6 qseqid sseqid sstart send sstrand length "
                              "qstart qend qlen sseq qseq mismatch",
                   "-evalue", str(a.evalue), "-max_target_seqs", str(a.max_hedef),
                   "-num_threads", str(a.is_parcacigi), "-dust", "no",
                   "-out", cikti]
            print("   blastn %-14s x %-22s (%d primer)"
                  % (sinif, dbad, len(primerler)))
            try:
                r = subprocess.run(cmd, capture_output=True, text=True,
                                   timeout=a.zaman_asimi)
            except subprocess.TimeoutExpired:
                # Skipping silently would show that database as "clean".
                print("      ZAMAN ASIMI (%d sn): %s OLCULEMEDI"
                      % (a.zaman_asimi, dbad))
                atlanan.append((sinif, dbad, "zaman asimi"))
                continue
            if r.returncode != 0:
                print(u'      ERROR: %s' % r.stderr.strip()[:200])
                atlanan.append((sinif, dbad, "blastn hatasi"))
                continue
            # vuruslari referans basina topla
            vurus = collections.defaultdict(list)
            with open(cikti) as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 12:
                        continue
                    q, s, sst, sen, strand = p[0], p[1], int(p[2]), int(p[3]), p[4]
                    qs, qe, ql = int(p[6]), int(p[7]), int(p[8])
                    sseq, qseq = p[9].upper(), p[10].upper()
                    # yalniz 3' ucu kapsayan vurus uzama yapabilir
                    if qe != ql:
                        continue
                    # When BLAST returns a short alignment, the 5' side of the primer
                    # falls outside the alignment and the mismatches there are not
                    # counted. Since the unseen part is unknown, the worst case is
                    # assumed: every unaligned 5' base counts as a possible mismatch and
                    # the total limit is checked against that.
                    gorulmeyen = qs - 1
                    ok, mm = baglanma_uygun(qseq, sseq)
                    if not ok:
                        continue
                    if mm + gorulmeyen > 3:
                        continue
                    mm += gorulmeyen
                    uc = sen        # 3' ucun referanstaki konumu
                    # The primer length travels too; the product length is NOT the
                    # distance between the two 3' ends but the distance from the 5' end
                    # of the forward primer to the 5' end of the reverse primer. In the
                    # old version the difference (lenF + lenR - 2, typically 44 bases)
                    # meant that real 70-94 bp off-target products fell below the lower
                    # bound and were never counted at all.
                    vurus[s].append((q, uc, strand, mm, ql))
            # on the same reference, on opposite strands, and with their 3' ends facing each other
            ekstra = collections.Counter()
            detay = {}
            urunler = collections.defaultdict(list)   # cift -> [(ref, boy)]
            for s, v in vurus.items():
                arti = [x for x in v if x[2] == "plus"]
                eksi = [x for x in v if x[2] == "minus"]
                for qf, uf, _, mf, lf in arti:
                    for qr, ur, _, mr, lr in eksi:
                        if ur <= uf:
                            continue
                        # The same measure as stage 04's product_len: the 5' end of the
                        # forward primer is uf - lf + 1, of the reverse primer ur + lr - 1
                        boy = (ur + lr - 1) - (uf - lf + 1) + 1
                        if a.prod_min <= boy <= a.prod_max:
                            i1 = re.sub(r"_[FR]\d+$", "", qf)
                            i2 = re.sub(r"_[FR]\d+$", "", qr)
                            n1 = qf.rsplit("_", 1)[1][1:]
                            n2 = qr.rsplit("_", 1)[1][1:]
                            if n1 != n2:
                                continue          # the primers of different pairs
                            ekstra[n1] += 1
                            urunler[n1].append((s, boy))
                            detay.setdefault(n1, []).append("%s:%d bp" % (s, boy))

            # resolve the taxon of the references that formed a product and split own from foreign
            tum_ref = {s for lst in urunler.values() for s, _ in lst}
            bas = basliklari_coz(fna, tum_ref, baslik_onbellek) if tum_ref else {}
            ref_token = {s: _tokenlar(bas.get(s, "")) for s in tum_ref}

            for i, r2 in enumerate(rows):
                if r2["sinif"] != sinif:
                    continue
                k = str(i)
                tk = taksonlar.get(r2["hedef"], {})
                evrensel = bool(tk.get("evrensel"))
                kendi_adlar = set(tk.get("beyan", ())) | set(tk.get("olculen", ()))
                kendi = yabanci = bilinmiyor = 0
                yab_ornek = []
                kendi_soylar, yab_kayit = [], []
                for s, boy in urunler.get(k, []):
                    tok = ref_token.get(s, set())
                    if evrensel:
                        kendi += 1
                    elif not tok:
                        bilinmiyor += 1
                    elif kendi_adlar & tok:
                        kendi += 1
                        kendi_soylar.append(_soyagaci(bas.get(s, "")))
                    else:
                        yabanci += 1
                        yab_kayit.append((s, boy))
                        if len(yab_ornek) < 5:
                            yab_ornek.append("%s:%d bp"
                                             % (bas.get(s, s)[:70], boy))
                # yabanci vuruslarin hedefin soyagacina uzakligi
                ref_soy = _baskin_soy(kendi_soylar)
                yakin = uzak = soysuz = 0
                for s, _boy in yab_kayit:
                    soy = _soyagaci(bas.get(s, ""))
                    if not soy or not ref_soy:
                        soysuz += 1
                    elif _ortak_derinlik(ref_soy, soy) >= max(1, len(ref_soy) - 2):
                        yakin += 1
                    else:
                        uzak += 1
                sonuc.append(dict(
                    hedef=r2["hedef"], sinif=sinif, veritabani=dbad,
                    ileri_dizi=r2["ileri_dizi"], geri_dizi=r2["geri_dizi"],
                    hedef_disi_urun=ekstra.get(k, 0),
                    urun_kendi_taksonda=kendi,
                    urun_yabanci_taksonda=yabanci,
                    yabanci_yakin=yakin,
                    yabanci_uzak=uzak,
                    yabanci_soyagacsiz=soysuz,
                    urun_takson_bilinmiyor=bilinmiyor,
                    hedef_soyagaci=";".join(ref_soy),
                    hedef_turu="evrensel" if evrensel else "ozgul",
                    ornekler=";".join(detay.get(k, [])[:5]),
                    yabanci_ornekler=";".join(yab_ornek),
                    kapsam_durumu=kap_durum,
                    kapsam_en_uzun_bp=kap_uz,
                    kapsam_kimlik=("%.1f" % kap_kim) if kap_uz else ""))
    try:
        shutil.rmtree(calisma, ignore_errors=True)
    except Exception:
        pass
    if not sonuc:
        print(u'no database could be scanned; no output was written')
        return
    if a.out:
        d = os.path.dirname(os.path.abspath(a.out))
        if d:
            os.makedirs(d, exist_ok=True)
        with open(a.out, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(sonuc[0].keys()),
                               delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(sonuc)
        print("\nyazildi: %s" % a.out)
    if atlanan:
        print("\nOLCULEMEYEN VERITABANLARI (temiz sayilmazlar):")
        for sn, db, sb in atlanan:
            print("   %-4s %-34s %s" % (sn, db, sb))
    if gecersiz:
        print(u'\nOUT-OF-COVERAGE PAIRS (a \'no product\' result on these rows is not evidence of specificity):')
        for sn, db, uz, kim in gecersiz:
            print("   %-4s %-34s en uzun benzer bolge %5d bp  %%%.1f"
                  % (sn, db, uz, kim))
    temiz = sum(1 for x in sonuc if x["hedef_disi_urun"] == 0)
    temiz_gecerli = sum(1 for x in sonuc if x["hedef_disi_urun"] == 0
                        and x.get("kapsam_durumu") == "KAPSANIYOR")
    print(u'records: %d, giving no product at all: %d (of those, with verified coverage: %d)'
          % (len(sonuc), temiz, temiz_gecerli))

    # THE RANKING IS NOT BY THE RAW COUNT. The reasoning is in the note above
    # hedef_taksonlari(): some of the highest raw counts are the target's own taxon and
    # are not an error.
    ozgul = [x for x in sonuc
             if x.get("hedef_turu") == "ozgul"
             and x.get("kapsam_durumu") == "KAPSANIYOR"
             and x.get("yabanci_uzak", 0) > 0]
    print(u'\nSPECIFICITY FINDINGS, ordered by DISTANT TAXON')
    print(u'(near = shares the target\'s lineage except for the last two ranks, usually the same functional group)')
    if not ozgul:
        print(u'   no record gives a product in a distant taxon')
    for x in sorted(ozgul, key=lambda x: -x["yabanci_uzak"])[:15]:
        print(u'   %-28s %-3s %-28s far=%5d near=%5d own=%5d'
              % (x["hedef"][:27], x["sinif"], x["veritabani"][:27],
                 x["yabanci_uzak"], x["yabanci_yakin"], x["urun_kendi_taksonda"]))
        if x["yabanci_ornekler"]:
            print("        %s" % x["yabanci_ornekler"].split(";")[0][:100])

    # a summary at pair level
    cift = {}
    for x in sonuc:
        if x.get("kapsam_durumu") != "KAPSANIYOR":
            continue
        k = (x["hedef"], x["ileri_dizi"], x["geri_dizi"])
        d = cift.setdefault(k, {"uzak": 0, "yakin": 0, "kendi": 0,
                                "turu": x.get("hedef_turu")})
        d["uzak"] += x.get("yabanci_uzak", 0)
        d["yakin"] += x.get("yabanci_yakin", 0)
        d["kendi"] += x.get("urun_kendi_taksonda", 0)
    oz = [v for v in cift.values() if v["turu"] == "ozgul"]
    # A pair giving NO PRODUCT AT ALL IN ITS OWN TAXON cannot count as "clean". That is
    # the pair level counterpart of the ROD mistake: it came out empty because there was
    # nothing to measure. Measured (2026-08-01): three pairs of Sakarolitik F2 give no
    # product in any of the eight covered databases, not even in their own taxon; those
    # are the pairs stage 18 already marks as alan_karisimi.
    olculebilir = [v for v in oz if v["kendi"] > 0]
    print(u'\nSPECIFIC pairs measured against a covered database: %d' % len(oz))
    print(u'   of those, MEASURABLE because they give a product in their own taxon: %d'
          % len(olculebilir))
    print(u'   giving no product in any distant taxon (genuinely clean)          : %d'
          % sum(1 for v in olculebilir if v["uzak"] == 0))
    inert = [v for v in oz if v["kendi"] == 0]
    if inert:
        print(u'   giving no product even in their own taxon (MEASUREMENT INVALID)   : %d'
              % len(inert))


if __name__ == "__main__":
    main()
