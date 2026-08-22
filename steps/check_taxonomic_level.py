#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_taxonomic_level.py
TESTING THE PANEL DECISION DIRECTLY: does every target separate at the level the
decision asked for (species or genus)?

The decision:
  species specific : Methanosarcina mazei, Methanothrix soehngenii,
                     Methanosarcina barkeri, Podospora pseudopauciseta,
                     Dictyostelium discoideum, Trichoderma asperellum
  genus specific   : Bacteroides, Alistipes, Proteiniphilum, Petrimonas
This list is not written by hand; it is read from the "duzey" column of
targets.tsv.

WHY A SEPARATE MEASUREMENT IS NEEDED
steps/specificity.py tests against the competitors in the sample, and
steps/external_databases.py looks for off-target products in external databases.
Neither asks "does this pair separate its target from the OTHER SPECIES OF ITS OWN
GENUS". Species specificity is exactly that question, and it can only be answered
against a panel built from sibling species.

THE METHOD
Each target's genus is taken from the declared taxid's name. The reference
databases are scanned, every record of that genus WITH A KNOWN SPECIES NAME is
collected, and a panel is built. blastn is run against the panel, and the products
are counted with external_databases.py's product rule (the same reference, opposite
strands, 3' ends facing one another, the product length in range) and the same
binding rule. The records that give a product are grouped by species.

The panel uses ONLY records with a known species name. Environmental databases
dominated by "uncultured archaeon", such as SILVA, cannot enter a species
separation panel; a record carrying no species identity neither confirms nor
refutes a species separation.

THE VERDICT
  duzey=tur : if there is a product in the target species and none in any other
              species of the same genus, TUR_OZGUL. Since the panel decision
              tolerates 1-2 cross reacting species, when the number of cross
              reacting SPECIES is below the threshold it is TUR_OZGUL_ESIKLI (the
              threshold changes with --cross-species-tolerance, default 2). Above the
              threshold it is TUR_AYRIMI_YOK. The measure is the number of cross
              reacting SPECIES, not the number of products formed in them.
              If the target species is absent from the panel altogether,
              HEDEF_TUR_PANELDE_YOK; species specificity can then be neither shown
              nor refuted.
  duzey=cins: the within-genus coverage is reported (how many species are
              amplified). Whether a product forms OUTSIDE the genus is
              external_databases.py's job and is not repeated here.

Usage:
  python3 check_taxonomic_level.py --targets targets.tsv       --names taxid_names.tsv --final final_primers --db REFERENCE_DB       --identity final_primers/hedef_kimlik.tsv       --out final_primers/duzey_denetimi.tsv

"""
import argparse, collections, csv, os, re, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
# The dependency is DELIBERATE: the binding and product rule must be EXACTLY the
# same as in steps/external_databases.py. Writing a separate copy opens the door to
# the two measurements drifting apart silently.
import importlib.util as _il
_s = _il.spec_from_file_location("_dv", os.path.join(HERE, "external_databases.py"))
DV = _il.module_from_spec(_s)
_yedek, sys.argv = sys.argv, ["external_databases.py"]
try:
    _s.loader.exec_module(DV)
except SystemExit:
    pass
finally:
    sys.argv = _yedek

# The databases that carry species names. SILVA is left out: the great majority of
# NR99 records end in 'uncultured archaeon/bacterium' and carry no species identity.
# Had they been taken into the panel, the species separation measurement would have
# been diluted with records that have no species name.
PANEL_DB = [
    "archaea.16S.fna",
    "bacteria.16S.fna",
    "fungi.ITS.fna",
    "fungi.18SrRNA.fna",
    "fungi.28SrRNA.fna",
    "UNITE_ITS.fasta",
    "ROD_v1.2_operon_variants.fasta",
    "PR2_SSU_taxo_long.fasta",
]

# 'sp', 'spp', 'cf' and 'aff' ARE NOT SPECIES NAMES.
#   MEASURED (2026-08-01): records of the form 's__Trichoderma_sp' in UNITE were
#   being counted as species names and stood in the panel as a separate "species";
#   16910 records for Trichoderma alone, 1712 for Marasmius, 1326 for Podospora.
#   A product forming in a record whose species is not known DOES NOT SHOW that the
#   species separation failed; which species it is, is unknown. If these enter the
#   panel, the panel's species count inflates and a wrong TUR_AYRIMI_YOK verdict can
#   be produced.
#   (The RefSeq form 'Trichoderma sp.' already carries a full stop and so did not
#   catch on the regular expression; the ones escaping were only the underscored forms.)
GURULTU = {"uncultured", "unidentified", "environmental", "sample", "clone",
           "isolate", "strain", "voucher", "candidatus", "bacterium",
           "archaeon", "fungal", "endophyte", "culture", "enrichment",
           "sp", "spp", "cf", "aff", "indet", "incertae", "sedis"}


def tur_adi(baslik):
    """Extracts the 'Genus species' binomial from a header; empty if there is none.

        RefSeq : 'NR_104707.1 Methanothrix soehngenii GP6 16S ...'
        UNITE  : '...;s__Thelephora_albomarginata|SH...'
        ROD    : '...;Drosophila;Drosophila_melanogaster|size=1'
        PR2    : '...|Unruhdinium|Unruhdinium_kevei'
        SILVA  : '... ;Methanothrix;uncultured archaeon'  -> empty

    """
    if not baslik:
        return ""
    # once alt cizgili ikili bicimler (UNITE, ROD, PR2)
    for parca in re.split(r"[|;]", baslik):
        p = parca.strip()
        p = re.sub(r"^[a-z]__", "", p)
        m = re.match(r"^([A-Z][a-z]+)_([a-z][a-z\-]+)$", p)
        if (m and m.group(1).lower() not in GURULTU
                and m.group(2).lower() not in GURULTU):
            return "%s %s" % (m.group(1), m.group(2))
    # then the space separated RefSeq form: the two words after the accession
    #
    # A SEMICOLON IN THE HEADER IS NOT ON ITS OWN A REASON TO DISCARD.
    #   MEASURED (2026-08-01): RefSeq ITS records have the form
    #   'NR_172285.1 Petriella musispora CBS 745.69 ITS region; from TYPE
    #   material'. In the earlier version, if the header held a semicolon this
    #   branch never ran, so the species name of the TYPE MATERIAL records in
    #   fungi.ITS.fna and fungi.28SrRNA.fna could not be extracted and they were
    #   not taken into the panel. Those are exactly the records the panel needs
    #   most: type strain sequences are the gold standard of species separation.
    #   (Measured: 35 records fell out this way for Petriella alone.)
    # SILVA's lineage headers cannot pass this branch anyway, because their
    # second word holds a semicolon and catches on the regular expression below;
    # no separate guard is needed.
    kelime = baslik.split()
    if len(kelime) >= 3:
        g, t = kelime[1], kelime[2]
        if (re.match(r"^[A-Z][a-z]+$", g) and re.match(r"^[a-z][a-z\-]+$", t)
                and g.lower() not in GURULTU and t.lower() not in GURULTU):
            return "%s %s" % (g, t)
    return ""


def ad_parcala(ad):
    """'Ca. Nitrosocosmicus hydrocola' -> ('Nitrosocosmicus','hydrocola')"""
    kelime = [k for k in re.split(r"[^A-Za-z]+", ad or "") if k]
    kelime = [k for k in kelime if k.lower() not in GURULTU and len(k) > 2]
    if not kelime:
        return ("", "")
    cins = kelime[0]
    tur = kelime[1].lower() if len(kelime) > 1 else ""
    return (cins, tur)


def referans_esle(ref_hedef, hedef_adlari):
    """Links the target name in reference_primers.tsv to the name in targets.tsv.

        STRIPPING THE '_referans' SUFFIX ALONE IS NOT ENOUGH.
          MEASURED (2026-08-01): 'Methanosarcina_barkeri_referans' strips to
          'Methanosarcina_barkeri', while the name in targets.tsv is
          'Methanosarcina_barkeri_turu'. When the match fails, the target's ONLY
          primer set (it has no de novo pair at all) drops silently and the target
          appears as CIFT_YOK. A silent drop makes a target that was never measured
          look as though no pair could be found.
        So an exact match is tried first, then a prefix match; if neither holds it
        returns None and the caller reports that LOUDLY.

    """
    kok = re.sub(r"_referans$", "", ref_hedef)
    if kok in hedef_adlari:
        return kok
    adaylar = [h for h in hedef_adlari if h.startswith(kok + "_")]
    if len(adaylar) == 1:
        return adaylar[0]
    return None


def hedefleri_oku(hedefler_tsv, adlar_tsv, kimlik_tsv):
    adlar = {}
    if adlar_tsv and os.path.exists(adlar_tsv):
        for line in open(adlar_tsv, encoding="utf-8"):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0].strip().isdigit():
                adlar[p[0].strip()] = p[1].strip()
    olculen = {}
    if kimlik_tsv and os.path.exists(kimlik_tsv):
        for r in csv.DictReader(open(kimlik_tsv, encoding="utf-8"),
                                delimiter="\t"):
            olculen[r["hedef"]] = re.sub(r"\(.*?\)", "",
                                         r.get("olculen_kimlik") or "").strip()
    out = []
    for line in open(hedefler_tsv, encoding="utf-8"):
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 4 or p[0] == "karar":
            continue
        duzey = p[2].strip()
        if duzey not in ("tur", "cins"):
            continue
        tidler = [t.strip() for t in p[3].split(",") if t.strip()]
        beyan = [adlar.get(t, "") for t in tidler if adlar.get(t)]
        cinsler, turler = set(), set()
        for ad in beyan:
            c, t = ad_parcala(ad)
            if c:
                cinsler.add(c)
            if c and t:
                turler.add("%s %s" % (c, t))
        # AN OPTIONAL 7TH COLUMN: hedef_tur
        # On some targets the declared taxid does not correspond to the organism found in
        # the sample, and the right species may have no entry in our taxid_names.tsv. An
        # example: the measured identity of the target whose bins are labelled 101201
        # (Trichoderma asperellum) is Petriella musispora. A TAXID IS NEVER INVENTED; the
        # species name is written directly in this column and added to the target's own name
        # set. When the column is empty nothing changes and the old behaviour continues.
        # The column REPLACES rather than ADDS. It exists because the taxid's name is wrong,
        # and adding the name would mean keeping the wrong name too. Measured: on the
        # Petriella_musispora row, in_taxid holds the bins' Kraken2 labels,
        # 101201/2034170/63577 (the in-sample support is measured with those). Had it added,
        # Trichoderma asperellum, atroviride and breve would have entered the target species
        # set too, and a pair amplifying Trichoderma would have counted as "there is a
        # product in the target species".
        if len(p) > 6 and p[6].strip():
            cinsler, turler = set(), set()
            for parca in p[6].split(","):
                c2, t2 = ad_parcala(parca)
                if c2:
                    cinsler.add(c2)
                if c2 and t2:
                    turler.add("%s %s" % (c2, t2))
        olc = olculen.get(p[1], "")
        oc, ot = ad_parcala(olc.split(";")[0]) if olc else ("", "")
        out.append({"hedef": p[1], "duzey": duzey, "beyan_ad": "; ".join(beyan),
                    "cinsler": cinsler, "hedef_turler": turler,
                    "olculen_ad": olc, "olculen_cins": oc,
                    "olculen_tur": ("%s %s" % (oc, ot)) if oc and ot else ""})
    return out


def panelleri_topla(dbklasor, cinsler, en_fazla_tur_basina, gunluk):
    """Scans the databases in one pass, collecting species labelled records per genus.
        Returns: {genus: {species: [(label, sequence), ...]}}"""
    panel = {c: collections.defaultdict(list) for c in cinsler}
    kirpilan = collections.Counter()
    dusen_tursuz = collections.Counter()
    dusen_baska_cins = collections.Counter()
    kucuk = {c.lower() for c in cinsler}
    for dbad in PANEL_DB:
        yol = os.path.join(dbklasor, dbad)
        if not os.path.exists(yol):
            gunluk.append(u'panel: there is no %s, skipped' % dbad)
            continue
        alinan = 0
        with open(yol, encoding="utf-8", errors="replace") as fh:
            baslik, dizi, sec = None, [], None
            for line in fh:
                if line.startswith(">"):
                    if sec and dizi:
                        c, t = sec
                        if len(panel[c][t]) < en_fazla_tur_basina:
                            panel[c][t].append(
                                ("%s|%s" % (dbad.split(".")[0], baslik[:60]),
                                 "".join(dizi)))
                            alinan += 1
                        else:
                            kirpilan[(c, t)] += 1
                    baslik = line[1:].rstrip("\n")
                    dizi, sec = [], None
                    dusuk = baslik.lower()
                    for c in cinsler:
                        if c.lower() in dusuk:
                            ta = tur_adi(baslik)
                            if ta and ta.split()[0].lower() == c.lower():
                                sec = (c, ta)
                            elif ta:
                                # The genus name appears in the header but the binomial
                                # belongs to another genus (it may be a strain name, a note
                                # or host information). It is not taken into the panel; it
                                # DOES NOT DROP SILENTLY, it is counted and reported.
                                dusen_baska_cins[c] += 1
                            else:
                                dusen_tursuz[c] += 1
                            break
                elif sec is not None:
                    dizi.append(line.strip())
            if sec and dizi:
                c, t = sec
                if len(panel[c][t]) < en_fazla_tur_basina:
                    panel[c][t].append(
                        ("%s|%s" % (dbad.split(".")[0], baslik[:60]),
                         "".join(dizi)))
                    alinan += 1
                else:
                    kirpilan[(c, t)] += 1
        gunluk.append(u'panel: %-34s %d records taken' % (dbad, alinan))
    # NO SILENT TRIMMING: the records trimmed, and the ones dropped for having no
    # species name, are reported; otherwise a gap in the panel reads as full coverage.
    for (c, t), n in sorted(kirpilan.items(), key=lambda x: -x[1])[:10]:
        gunluk.append(u'the panel was TRIMMED: %d records for %s / %s were not taken'
                      % (n, c, t))
    for c, n in dusen_tursuz.items():
        if n:
            gunluk.append(u'panel: %d records carrying no species name for %s WERE NOT TAKEN into the panel' % (n, c))
    for c, n in dusen_baska_cins.items():
        if n:
            gunluk.append(u'panel: %d records where the name %s appears but the binomial belongs to another genus WERE NOT TAKEN into the panel' % (n, c))
    return panel


def urun_say(primerler, panel_fa, calisma, etiket, a):
    """blastn against the panel, applying the outside database step's binding and
        product rule. Returns: {pair_number: {reference: product_count}}"""
    db = os.path.join(calisma, "panel_%s.fa" % etiket)
    shutil.copyfile(panel_fa, db)
    r = subprocess.run(["makeblastdb", "-in", db, "-dbtype", "nucl"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return None
    sorgu = os.path.join(calisma, "sorgu_%s.fa" % etiket)
    with open(sorgu, "w") as fh:
        for ad, dizi in primerler.items():
            fh.write(">%s\n%s\n" % (ad, dizi))
    cikti = sorgu + ".tsv"
    cmd = ["blastn", "-task", "blastn-short", "-query", sorgu, "-db", db,
           "-outfmt", "6 qseqid sseqid sstart send sstrand length "
                      "qstart qend qlen sseq qseq mismatch",
           "-evalue", str(a.evalue), "-max_target_seqs", "100000",
           "-num_threads", str(a.threads), "-dust", "no", "-out", cikti]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        return None
    vurus = collections.defaultdict(list)
    with open(cikti) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            q, s, strand = p[0], p[1], p[4]
            qs, qe, ql = int(p[6]), int(p[7]), int(p[8])
            sen = int(p[3])
            sseq, qseq = p[9].upper(), p[10].upper()
            if qe != ql:
                continue
            gorulmeyen = qs - 1
            ok, mm = DV.baglanma_uygun(qseq, sseq)
            if not ok or mm + gorulmeyen > 3:
                continue
            vurus[s].append((q, sen, strand, ql))
    say = collections.defaultdict(collections.Counter)
    for s, v in vurus.items():
        arti = [x for x in v if x[2] == "plus"]
        eksi = [x for x in v if x[2] == "minus"]
        for qf, uf, _, lf in arti:
            for qr, ur, _, lr in eksi:
                if ur <= uf:
                    continue
                boy = (ur + lr - 1) - (uf - lf + 1) + 1
                if not (a.prod_min <= boy <= a.prod_max):
                    continue
                n1 = qf.rsplit("_", 1)[1][1:]
                n2 = qr.rsplit("_", 1)[1][1:]
                if n1 == n2:
                    say[n1][s] += 1
    return say


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--targets", required=True)
    p.add_argument("--names", required=True)
    p.add_argument("--final", required=True, help="final_primers directory")
    p.add_argument("--reference", default=None, help="reference_primers.tsv")
    p.add_argument("--db", required=True, help='the reference database '
                                               'directory')
    p.add_argument("--identity", default=None, help='the identity table produced by '
                                                    'the target identity step')
    p.add_argument("--out", required=True)
    p.add_argument("--prod-min", type=int, default=50)
    p.add_argument("--prod-max", type=int, default=400)
    p.add_argument("--evalue", type=float, default=1000.0)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--max-per-species", type=int, default=200)
    # The panel decision: "where there are 1-2 cross reacting species, it still counts
    # as species specific". The value is not hard coded here but stands as an option; if
    # it is changed, which threshold the run used is written at the top of the output.
    p.add_argument("--cross-species-tolerance", type=int, default=2,
                   help='the number of cross reacting SPECIES tolerated at '
                        'species level specificity, not the number of '
                        'products; the default is 2')
    return p.parse_args()


def main():
    a = get_args()
    gunluk = []
    hedefler = hedefleri_oku(a.targets, a.names, a.identity)
    if not hedefler:
        sys.exit(u'targets.tsv holds no row with duzey=tur or duzey=cins')
    print(u'targets with a declared decision level: %d (species: %d, genus: %d)'
          % (len(hedefler),
             sum(1 for h in hedefler if h["duzey"] == "tur"),
             sum(1 for h in hedefler if h["duzey"] == "cins")))
    print(u'cross reacting SPECIES tolerated under species specificity: %d'
          % a.cross_species_tolerance)

    # ciftleri topla
    ciftler = collections.defaultdict(list)
    tsv = os.path.join(a.final, "final_primers.tsv")
    for r in csv.DictReader(open(tsv, encoding="utf-8"), delimiter="\t"):
        if r.get("ozgulluk_durum") == "GECTI":
            ciftler[r["hedef"]].append((r["ileri_dizi"], r["geri_dizi"],
                                        "de novo"))
    if a.reference and os.path.exists(a.reference):
        # ALL the names in targets.tsv (without distinguishing the level), because the
        # reference set also holds targets at level=group
        tum_ad = []
        for line in open(a.targets, encoding="utf-8"):
            if line.startswith("#") or not line.strip():
                continue
            pp = line.rstrip("\n").split("\t")
            if len(pp) >= 2 and pp[0] != "karar":
                tum_ad.append(pp[1])
        eslesmeyen = collections.Counter()
        for r in csv.DictReader(open(a.reference, encoding="utf-8"),
                                delimiter="\t"):
            ham = r.get("hedef", "")
            ad = referans_esle(ham, tum_ad)
            if ad is None:
                eslesmeyen[ham] += 1
                continue
            ciftler[ad].append((r["ileri_dizi"], r["geri_dizi"], "referans"))
        for ham, n in eslesmeyen.items():
            print(u'   WARNING: reference target \'%s\' (%d pairs) matched no name in targets.tsv and was LEFT OUT OF THE MEASUREMENT' % (ham, n))

    tum_cins = set()
    for h in hedefler:
        tum_cins |= h["cinsler"]
        if h["olculen_cins"]:
            tum_cins.add(h["olculen_cins"])
    print("panel kurulacak cins: %s" % ", ".join(sorted(tum_cins)))
    panel = panelleri_topla(a.db, tum_cins, a.max_per_species, gunluk)
    for g in gunluk:
        print("   " + g)

    calisma = tempfile.mkdtemp(prefix="duzey_")
    sonuc = []
    for h in hedefler:
        hd = h["hedef"]
        cf = ciftler.get(hd, [])
        # panel: beyan edilen cins(ler) + olculen cins
        kendi_cins = set(h["cinsler"])
        if h["olculen_cins"]:
            kendi_cins.add(h["olculen_cins"])
        kayitlar = []
        for c in sorted(kendi_cins):
            for t, lst in sorted(panel.get(c, {}).items()):
                for i, (etiket, dizi) in enumerate(lst):
                    kayitlar.append(("%s__%s__%d" % (c, t.replace(" ", "_"), i),
                                     t, dizi))
        panel_turler = sorted({t for _, t, _ in kayitlar})
        if not cf:
            sonuc.append(dict(hedef=hd, duzey=h["duzey"], cift_no="",
                              kaynak="", beyan_ad=h["beyan_ad"],
                              olculen_ad=h["olculen_ad"],
                              panel_tur_sayisi=len(panel_turler),
                              hedef_turde_urun="", diger_turde_urun="",
                              cogaltilan_hedef_turler="",
                              capraz_tur_sayisi="", cogaltilan_turler="",
                              karar="CIFT_YOK"))
            continue
        if not kayitlar:
            for i, (F, R, kaynak) in enumerate(cf):
                sonuc.append(dict(hedef=hd, duzey=h["duzey"], cift_no=str(i),
                                  kaynak=kaynak, beyan_ad=h["beyan_ad"],
                                  olculen_ad=h["olculen_ad"],
                                  panel_tur_sayisi=0, hedef_turde_urun="",
                                  diger_turde_urun="",
                                  cogaltilan_hedef_turler="",
                                  capraz_tur_sayisi="",
                                  cogaltilan_turler="", karar="PANEL_YOK"))
            continue
        panel_fa = os.path.join(calisma, "%s_panel.fa" % hd)
        with open(panel_fa, "w") as fh:
            for kimlik, _t, dizi in kayitlar:
                fh.write(">%s\n%s\n" % (kimlik, dizi))
        primerler = {}
        for i, (F, R, _k) in enumerate(cf):
            primerler["%s_F%d" % (hd[:20], i)] = F
            primerler["%s_R%d" % (hd[:20], i)] = R
        print(u'   %-32s %2d pairs x %3d panel records (%d species)'
              % (hd, len(cf), len(kayitlar), len(panel_turler)))
        say = urun_say(primerler, panel_fa, calisma, hd[:16], a)
        if say is None:
            print(u'      blastn/makeblastdb failed, skipped')
            continue
        kimlik_tur = {k: t for k, t, _ in kayitlar}
        # the target species set: the declared species plus the measured species
        hedef_turler = set(h["hedef_turler"])
        if h["olculen_tur"]:
            hedef_turler.add(h["olculen_tur"])
        for i, (F, R, kaynak) in enumerate(cf):
            per = say.get(str(i), collections.Counter())
            tur_urun = collections.Counter()
            for ref, n in per.items():
                tur_urun[kimlik_tur.get(ref, "?")] += n
            # WHICH target species the product formed in is written out. The target
            # species set holds the declared name AND the measured identity together,
            # and the two can diverge. An example: the declared species of
            # Zoopagomycota_mantari is Dictyostelium discoideum while its measured
            # identity is a below-threshold Marasmius. Without writing which one the
            # product formed in, a TUR_OZGUL_ESIKLI verdict could be read as "specific
            # to Dictyostelium", when it may be specific to Marasmius.
            cogaltilan_hedef = sorted(t for t in tur_urun if t in hedef_turler)
            hedefte = sum(n for t, n in tur_urun.items() if t in hedef_turler)
            digerde = sum(n for t, n in tur_urun.items()
                          if t not in hedef_turler)
            digerler = sorted({t for t in tur_urun if t not in hedef_turler})
            hedef_panelde = bool(hedef_turler & set(panel_turler))
            if h["duzey"] == "tur":
                # THE NUMBER OF CROSS REACTING SPECIES, not the number of products. The
                # panel decision says "where there are 1-2 cross reacting species it still
                # counts as species specific"; the measure is HOW MANY DIFFERENT SPECIES
                # were amplified, not how many products formed in them. Zero cross reactions
                # and a within-threshold situation are written as separate verdicts; folding
                # the two into one label would equate the weaker pair with the stronger one.
                if not hedef_panelde:
                    karar = "HEDEF_TUR_PANELDE_YOK"
                elif hedefte == 0:
                    karar = "HEDEF_TURDE_URUN_YOK"
                elif len(digerler) == 0:
                    karar = "TUR_OZGUL"
                elif len(digerler) <= a.cross_species_tolerance:
                    karar = "TUR_OZGUL_ESIKLI"
                else:
                    karar = "TUR_AYRIMI_YOK"
            else:
                # GENUS LEVEL. The earlier version counted only "in how many species of
                # the panel does a product form", which hid products forming OUTSIDE the
                # genus.
                #   MEASURED (2026-08-01): two of Proteiniphilum_cinsi's five pairs also
                #   amplify Fermentimonas caenicola, which is another genus. The old count
                #   wrote both as CINS_ICI_3_3, so a pair violating genus specificity
                #   looked like the pair with the widest coverage.
                # Since the acceptance criterion is "genus specific", no product may form
                # outside the DECLARED genus. The measured identity does not count as the
                # target here: if Proteiniphilum was asked for, Fermentimonas is a cross
                # reaction, even when the organism in the bins is that one.
                tum_cogaltilan = sorted(t for t in tur_urun if t)
                ici = [t for t in tum_cogaltilan
                       if t.split()[0] in h["cinsler"]]
                disi = [t for t in tum_cogaltilan
                        if t.split()[0] not in h["cinsler"]]
                if not ici:
                    karar = "CINS_ICINDE_URUN_YOK"
                elif not disi:
                    karar = "CINS_OZGUL"
                else:
                    karar = "CINS_AYRIMI_YOK"
                cogaltilan_hedef = ici
                digerler = disi
                hedefte = sum(n for t, n in tur_urun.items() if t in ici)
                digerde = sum(n for t, n in tur_urun.items() if t in disi)
                karar += "_%d_%d" % (len(ici), len(panel_turler))
            sonuc.append(dict(
                hedef=hd, duzey=h["duzey"], cift_no=str(i), kaynak=kaynak,
                beyan_ad=h["beyan_ad"], olculen_ad=h["olculen_ad"],
                panel_tur_sayisi=len(panel_turler),
                hedef_turde_urun=hedefte, diger_turde_urun=digerde,
                cogaltilan_hedef_turler="; ".join(cogaltilan_hedef),
                capraz_tur_sayisi=len(digerler),
                cogaltilan_turler="; ".join(digerler[:8]), karar=karar))
    shutil.rmtree(calisma, ignore_errors=True)

    d = os.path.dirname(os.path.abspath(a.out))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sonuc[0].keys()),
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(sonuc)
    print(u'\nwritten: %s' % a.out)

    print("\n%-32s %-5s %-6s %-26s %s"
          % ("hedef", "duzey", "cift", "karar", "cogaltilan diger turler"))
    print("-" * 118)
    for h in hedefler:
        satir = [x for x in sonuc if x["hedef"] == h["hedef"]]
        # The species list is gathered PER VERDICT. Gathered across the target as a whole,
        # the species amplified by another pair were being written beside the TUR_OZGUL row
        # that did achieve the separation, and that pair read as if it separated nothing.
        gruplu = collections.defaultdict(list)
        for x in satir:
            gruplu[x["karar"]].append(x)
        for karar, lst in sorted(gruplu.items(), key=lambda x: -len(x[1])):
            digerler = sorted({d for x in lst for d in
                               x["cogaltilan_turler"].split("; ") if d})
            hedefteki = sorted({d for x in lst for d in
                                x["cogaltilan_hedef_turler"].split("; ") if d})
            print("%-32s %-5s %-6d %-22s %s"
                  % (h["hedef"][:31], h["duzey"], len(lst), karar,
                     ("; ".join(digerler))[:40]))
            if hedefteki:
                print("%-32s %s"
                      % ("", u'  amplified in the target species: '
                         + "; ".join(hedefteki)[:60]))
    tur_hedef = {h["hedef"] for h in hedefler if h["duzey"] == "tur"}
    kati = {x["hedef"] for x in sonuc if x["karar"] == "TUR_OZGUL"}
    esikli = {x["hedef"] for x in sonuc
              if x["karar"] in ("TUR_OZGUL", "TUR_OZGUL_ESIKLI")}
    print(u'\ntargets where species specificity was requested: %d' % len(tur_hedef))
    print(u'   at least one pair with NO CROSS-REACTION (SPECIES-SPECIFIC) : %d'
          % len(kati & tur_hedef))
    print(u'   at least one pair within the threshold (<= %d cross-reacting species) : %d'
          % (a.cross_species_tolerance, len(esikli & tur_hedef)))


if __name__ == "__main__":
    main()
