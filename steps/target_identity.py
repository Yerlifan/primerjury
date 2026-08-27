#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
target_identity.py
COMPARES THE TARGET NAME WITH THE ORGANISM THE DATA SHOWS.

The target names in the meeting decisions rest on Kraken2's taxon assignments.
Kraken2 does not abstain under its default setting: if the real organism is not in
the database, a read is labelled with the highest scoring leaf. So a bin's label and
what the SEQUENCE in that bin really is can come apart.

This script queries the bin consensuses against the reference databases with blastn
and puts three things side by side for each target:
  toplanti_adi     the target name written in the decision
  kraken_etiketi   the names of that target's taxids
  olculen_kimlik   the consensus's best match in the reference databases

EVERY CLASS SEES EVERY rDNA DATABASE. The discriminating region of the class comes
first, but the rest of the databases are asked as well, and the script writes on
every row which region and which database the name came from, so the level of
confidence is not hidden.

Usage:
  python3 target_identity.py --consensus reference_consensus/dominant/consensus \
      --db REFERENCE_DB --targets targets.tsv --names taxid_names.tsv \
      --out final_primers/hedef_kimlik.tsv
"""
import argparse, csv, collections, glob, os, re, subprocess, sys, tempfile, shutil

_KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if os.path.join(_KOK, 'verification') not in sys.path:
    sys.path.insert(0, os.path.join(_KOK, 'verification'))

# The checkpoint ledger stores the expensive blastn outputs and rescues an
# interrupted run. If it cannot be loaded the script still works, only without the
# rescue.
try:
    from checkpoint import Defter
except Exception as _e2:                       # pragma: no cover
    print(u'WARNING: the checkpoint module could not be loaded (%s); an interrupted '
          u'run will start from the beginning' % type(_e2).__name__)

    class Defter(object):
        def __init__(self, *a, **k):
            pass

        def var(self, _a):
            return False

        def al(self, _a, varsayilan=None):
            return varsayilan

        def yaz(self, _a, _d):
            pass

        def kapat(self):
            pass


# --------------------------------------------------------------- THRESHOLDS
# THE 2026-08-25 FIX. This script used to apply only A FLOOR (90 per cent) through
# --min-identity; on any match above the floor, if the name held, the verdict became
# "tur_uyusuyor" outright. Yet in the same project verification/identity_verification.py
# takes the species threshold for the SSU to be 98.70 per cent. One project was
# using two different yardsticks, and the concrete result was this (the 2026-08-25
# run):
#     Methanosarcina_mazei_turu -> "tur_uyusuyor", the evidence 98.59 and 98.66
# that is, a match BELOW its own species threshold was producing a species name.
#
# The answer: the thresholds are read FROM ONE SOURCE. TUR_ESIGI, CINS_ESIGI and
# AYRIM_PAYI are imported from identity_verification.py, so when that file changes
# this script follows on its own. If the import fails (the directory was moved) the
# same values are used as a backup and A WARNING is printed while running; there is
# no silent return to the old behaviour.
try:
    from identity_verification import (TUR_ESIGI, CINS_ESIGI, AYRIM_PAYI,
                                       EN_AZ_HIZALAMA, ADSIZ_JETONLARI)
    _ESIK_KAYNAK = 'verification/identity_verification.py'
except Exception as _e:
    TUR_ESIGI = {'SSU': 98.7, 'LSU': 98.7, 'LSU_MANTAR': 99.8,
                 'ITS': 99.6, 'OPERON': 98.7, 'KARISIK': 98.7}
    CINS_ESIGI = {'SSU': 94.5, 'LSU': 94.5, 'LSU_MANTAR': 98.2,
                  'ITS': 94.3, 'OPERON': 94.5, 'KARISIK': 94.5}
    EN_AZ_HIZALAMA = {'SSU': 1200, 'LSU': 600, 'LSU_MANTAR': 600,
                      'ITS': 600, 'OPERON': 1200, 'KARISIK': 1200}
    ADSIZ_JETONLARI = ('uncultured', 'unclassified', 'unidentified',
                       'environmental', 'metagenome', 'enrichment', 'clone')
    AYRIM_PAYI = 0.5
    _ESIK_KAYNAK = 'A BACKUP COPY (%s)' % type(_e).__name__


# The fungal combination rule lives in its own module: a bin of the F classes
# carries a whole operon and its three loci can be combined by a vote rather than
# a ladder. Losing the module only costs that one option.
try:
    import locus_decision as _LD
except Exception:                              # pragma: no cover
    _LD = None


def cins_epitet(ad):
    """(genus, epithet) out of a name; the epithet is empty when there is none."""
    cins = cins_ayikla(ad)
    p = (ad or '').split()
    return cins, (p[-1] if len(p) > 1 and p[-1] != cins else '')


# --------------------------------------------------- WHICH DATABASE IS ASKED
# EVERY CLASS SEES EVERY rDNA DATABASE.
#
# The earlier version asked each class only of its own hand picked shortlist:
# archaea and bacteria were asked ONLY of RefSeq 16S, fungi only of ITS and 28S.
# Two separate faults came out of that.
#
#   1. A SILENT LOSS OF COVERAGE. archaea.16S.fna holds 1,160 records and
#      bacteria.16S.fna 26,877, while SILVA SSU NR99 sits in the same directory
#      with 510,495. MEASURED (2026-08-26): once SILVA SSU was added, three bins
#      changed name and all three came into line with Kraken or with the claim
#      itself: an unnameable bin became "Synergistaceae", "Nitrososphaera sp."
#      became "Candidatus Nitrosocosmicus", and "Fermentimonas" became
#      "Proteiniphilum".
#
#   2. LIMITING BY DOMAIN ASSUMES THE LABEL IS RIGHT. This script exists precisely
#      because a bin's label and its sequence can come apart. A bin labelled
#      bacterial that is really fungal can only be caught by asking the fungal
#      databases. Asking a 16S consensus of an ITS set costs a few seconds and
#      returns nothing when the label is right; when the label is wrong it is the
#      only thing that shows it.
#
# The same rule was already adopted for the specificity scan, for the same reason,
# and is documented in steps/external_databases.py: there is no scientific reason
# to limit by domain, only a reason of speed. --databases narrow is there for when
# speed matters.
#
# THE ORDER IS THE LADDER OF THE CLASS: the discriminating region of the class
# comes first. Which locus a name is taken from is decided by SINIF_BASAMAK below,
# not by this order, but asking the discriminating region first means the
# checkpoint fills up in a useful order.
#
# NOTE: SILVA is full of environmental clones; unnamed records produce no name at
# all (see ad_ayikla and ADSIZ_JETONLARI), they only ever say "the nearest record
# is this one".
TUM_VTB = [
    ("archaea.16S.fna",                "16S"),
    ("bacteria.16S.fna",               "16S"),
    ("SILVA_138.2_SSURef_NR99.fasta",  "16S"),
    ("SILVA_138.2_LSURef_NR99.fasta",  "23S"),
    ("ROD_v1.2_operon_variants.fasta", "OPERON"),
    ("fungi.ITS.fna",                  "ITS"),
    ("UNITE_ITS.fasta",                "ITS"),
    ("fungi.28SrRNA.fna",              "28S"),
    ("fungi.18SrRNA.fna",              "18S"),
    ("PR2_SSU_taxo_long.fasta",        "18S"),
]

# The discriminating database of each class, asked first.
AYIRT_EDICI = {
    "A1": [("archaea.16S.fna", "16S"), ("SILVA_138.2_SSURef_NR99.fasta", "16S")],
    "A2": [("archaea.16S.fna", "16S"), ("SILVA_138.2_SSURef_NR99.fasta", "16S")],
    "B":  [("bacteria.16S.fna", "16S"), ("SILVA_138.2_SSURef_NR99.fasta", "16S")],
    "F1": [("fungi.ITS.fna", "ITS"), ("UNITE_ITS.fasta", "ITS"),
           ("fungi.28SrRNA.fna", "28S"), ("fungi.18SrRNA.fna", "18S")],
    "F2": [("fungi.ITS.fna", "ITS"), ("UNITE_ITS.fasta", "ITS"),
           ("fungi.28SrRNA.fna", "28S"), ("fungi.18SrRNA.fna", "18S")],
}

# --databases narrow: the discriminating set alone. Faster, and it is what the
# earlier version did; it is kept so an old run can be reproduced.
SINIF_DB_DAR = dict(AYIRT_EDICI)

# --databases all (the default): the discriminating set first, then everything else.
SINIF_DB_GENIS = {
    s: AYIRT_EDICI[s] + [d for d in TUM_VTB if d not in AYIRT_EDICI[s]]
    for s in AYIRT_EDICI
}

SINIF_DB = dict(SINIF_DB_GENIS)

# The region asked of blastn -> the locus key in the threshold table.
BOLGE_LOKUS = {'16S': 'SSU', '18S': 'SSU', '23S': 'LSU',
               '28S': 'LSU_MANTAR', 'ITS': 'ITS', 'OPERON': 'OPERON'}

# THE LADDER OF EACH CLASS, which is where the name is taken from.
#
# This is the method as reported: "the discriminating region of the class is looked
# at first: 16S in archaea and bacteria, ITS in fungi. If nothing comes back from
# ITS then 28S is used, and if that gives nothing either then 18S."
#
# WHY A LADDER AND NOT THE HIGHEST PERCENTAGE: 99 per cent in 28S and 99 per cent
# in ITS are not the same thing. Racing the raw identities of different loci
# against one another produces a wrong winner. MEASURED (2026-08-26): with the
# fungal loci asked together and compared on raw identity, the Petriella targets
# fell to unnameable. The loci are therefore walked in order and the first one
# that gives a usable answer supplies the name.
#
# The loci of the other domain sit at the end of every ladder, so that a bin whose
# label has the domain wrong is still named instead of being left empty.
SINIF_BASAMAK = {
    "A1": ["16S", "23S", "OPERON", "ITS", "28S", "18S"],
    "A2": ["16S", "23S", "OPERON", "ITS", "28S", "18S"],
    "B":  ["16S", "23S", "OPERON", "ITS", "28S", "18S"],
    "F1": ["ITS", "28S", "18S", "OPERON", "16S", "23S"],
    "F2": ["ITS", "28S", "18S", "OPERON", "16S", "23S"],
}
VARSAYILAN_BASAMAK = ["16S", "ITS", "28S", "18S", "23S", "OPERON"]


def ad_ayikla(baslik):
    """Pull 'Genus species' (or the genus alone) out of a reference header.

    Three shapes are supported:
      RefSeq : "NR_041956.1 Methanosarcina mazei strain DSM 2053 ..."
      SILVA  : "CP009514.x.y Archaea;Halobacteriota;...;Methanosarcina mazei"
      UNITE  : "AY882347|k__Fungi;p__Ascomycota;...;s__Petriella_setifera"

    On the SILVA and UNITE path the LAST element IS NOT always a binomial; it can
    be a two word GENUS such as "Candidatus Nitrosocosmicus", a string carrying no
    information such as "uncultured bacterium", or a family name. So the path is
    walked FROM THE END BACKWARDS and the first element THAT CARRIES INFORMATION is
    taken. The earlier version looked only at the last element; once SILVA was
    added, names were being written out as raw taxonomy strings such as
    "Bacteria;Bacteroidota;..." (seen on 2026-08-26).
    """
    b = (baslik or "").strip()
    if not b:
        return None
    # The RefSeq shape is tried FIRST: "ACCESSION Genus species ...".
    # WHY FIRST: RefSeq headers can carry a SEMICOLON
    #   "NR_172285.1 Petriella musispora CBS 745.69 ITS region; from TYPE material"
    # That was taken for a taxonomy path and walked from the end backwards; with
    # "from TYPE material" as the last element no name came out and the bin was
    # left unnamed (seen on the Petriella targets on 2026-08-26). Patterns anchored
    # at the start are not affected by a semicolon, so they are tried first.
    m = re.match(r"^\S+\s+(Candidatus [A-Z][a-z]+\s+[a-z][a-z-]+)", b)
    if m:
        return m.group(1)
    m = re.match(r"^\S+\s+([A-Z][a-z]+\s+[a-z][a-z-]+)", b)
    if m:
        return m.group(1)
    if ";" in b:
        yol = [x.strip() for x in b.split(";") if x.strip()]
        for oge in reversed(yol):
            ad = re.sub(r"^[kpcofgs]__", "", oge).replace("_", " ").strip()
            # the first element can carry an accession: "CP009514.x.y Archaea"
            ad = re.sub(r"^\S*\d\S*\s+", "", ad).strip()
            if not ad or any(j in ad.lower() for j in
                             ("uncultured", "unidentified", "unclassified",
                              "environmental", "metagenome", "incertae",
                              "enrichment", "unknown")):
                continue
            m = re.match(r"^(Candidatus [A-Z][a-z]+ [a-z][a-z-]+)$", ad)
            if m:
                return m.group(1)
            m = re.match(r"^([A-Z][a-z]+ [a-z][a-z-]+)$", ad)
            if m:
                return m.group(1)
            m = re.match(r"^(Candidatus [A-Z][a-z]{2,})$", ad)
            if m:
                return m.group(1)
            m = re.match(r"^([A-Z][a-z]{2,})$", ad)
            if m:
                return m.group(1)
        return None
    m = re.match(r"\S+\s+([A-Z][a-z]+\s+[a-z][a-z-]+)", b)
    if m:
        return m.group(1)
    m = re.match(r"\S+\s+(Candidatus [A-Z][a-z]+\s+[a-z][a-z-]+)", b)
    if m:
        return m.group(1)
    m = re.match(r"\S+\s+(\S+\s+\S+)", b)
    return m.group(1) if m else None


# The suffixes of ranks above genus. A SILVA or UNITE lineage whose deeper
# elements are all unnamed ("uncultured bacterium") leaves the deepest element
# that carries information sitting at family, order or phylum level, and the
# name that comes out is then "Bacteroidota" or "Synergistaceae" rather than a
# genus. That IS worth reporting: it says the taxon is real and this is as deep
# as the reference goes. But it is not a genus, so no species threshold and no
# "cf." may be applied to it, and it must not be silently read as one.
UST_RANK_SONU = ('ota', 'aceae', 'ales', 'ineae', 'mycota', 'mycetes',
                 'bacteria', 'archaeota', 'phyta', 'idae')


def ust_rank_mi(ad):
    """Is this name a rank ABOVE genus? A binomial never is."""
    if not ad or ' ' in ad:
        return False
    return ad.endswith(UST_RANK_SONU)


def tur_adi(baslik):
    """The 'Genus species' pair out of a reference header, or None."""
    ad = ad_ayikla(baslik)
    return ad if ad and " " in ad else None


def cins_ayikla(ad):
    """The GENUS out of a name. 'Candidatus X' is a two word genus; taking the
    first token was producing meaningless names such as 'Candidatus cf. exaquare'."""
    p = (ad or "").split()
    if not p:
        return ""
    if p[0] == "Candidatus" and len(p) > 1:
        return "Candidatus " + p[1]
    return p[0]


def lokus_duzelt(bolge, baslik):
    """The locus key, corrected by THE DOMAIN OF THE HIT.

    The locus is a property of the database, not of the query: SILVA LSU NR99
    carries both prokaryotic and eukaryotic records and arrives under '23S'. A
    fungal hit in that set would have the prokaryotic threshold (98.7 per cent)
    applied to it, while the fungal LSU needs 99.8. The lineage in the header says
    which domain it is.
    """
    lokus = BOLGE_LOKUS.get(bolge, 'SSU')
    if lokus == 'LSU' and re.search(r'(Eukaryota|Fungi|Dikarya)', baslik or ''):
        return 'LSU_MANTAR'
    return lokus


def esik_uygula(adi, pid, lokus, aln=0, rakip=None):
    """A species name is given only if the SPECIES THRESHOLD of the locus is passed.

    Returns (name, note). If the threshold is not passed the name comes down to
    genus level ("Methanosarcina sp."), and if the genus threshold is not passed
    either no name is given. THE NAME IS NOT THROWN AWAY; what happened stands in
    the note, because "the nearest record is this one" is worth having even when it
    is not an identity.
    """
    if adi and adi.startswith(u"cannot be named"):
        return adi, u""
    te = TUR_ESIGI.get(lokus, 98.7)
    ce = CINS_ESIGI.get(lokus, 94.5)
    enaz = EN_AZ_HIZALAMA.get(lokus, 600)
    if ust_rank_mi(adi):
        # A rank above genus. Neither the species threshold nor "cf." applies to
        # it: there is no species claim here to qualify. The genus threshold is
        # still asked, because below it even the lineage is not evidence.
        if pid >= ce:
            return adi, (u'a rank above genus: the record matched carries no genus '
                         u'name, so this is as deep as the reference goes')
        return u'cannot be named', (u'below the genus threshold of %.2f per cent '
                                    u'(%s), no name can be given' % (ce, lokus))
    if pid >= te:
        # A short alignment is not evidence of a species. 100 per cent over 484
        # bases and 100 per cent over 2900 bases are not the same evidence; the
        # hit is not thrown away, it is only kept from rising to species level.
        if aln and aln < enaz:
            cins = cins_ayikla(adi)
            return (cins + ' sp.' if cins else adi,
                    u'the alignment is only %d bases and the %s locus wants at least '
                    u'%d for a species name, so it came down to genus level'
                    % (aln, lokus, enaz))
        # The species threshold was passed. But if ANOTHER SPECIES sits at nearly
        # the same closeness, a species name is indefensible. MEASURED (A1-2_2209):
        # M. mazei at 100.000 while M. soligelidi is at 99.925, which is 0.075
        # apart against a separation margin of 0.5.
        if rakip:
            rpid, rad = rakip
            if rad and rpid is not None and (pid - rpid) < AYRIM_PAYI:
                cins = cins_ayikla(adi)
                _p = (adi or '').split()
                epitet = _p[-1] if len(_p) > 1 else ''
                return ('%s cf. %s' % (cins, epitet) if cins and epitet else adi,
                        u'the second species %s is only %.2f per cent behind at '
                        u'%.2f per cent, against a separation margin of %.2f, so the '
                        u'species name was given with "cf."'
                        % (rad, pid - rpid, rpid, AYRIM_PAYI))
        return adi, ''
    cins = cins_ayikla(adi)
    if pid >= ce and cins:
        return cins + ' sp.', (u'below the species threshold of %.2f per cent (%s), '
                               u'so it came down to genus level' % (te, lokus))
    return u'cannot be named', (u'below the genus threshold of %.2f per cent (%s), '
                                u'no name can be given' % (ce, lokus))


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--consensus", required=True)
    p.add_argument("--db", required=True)
    p.add_argument("--targets", default="targets.tsv")
    p.add_argument("--names", default="taxid_names.tsv")
    p.add_argument("--out", required=True)
    p.add_argument("--min-alignment", type=int, default=250,
                   help='an alignment shorter than this does not count as an '
                        'identity')
    p.add_argument("--min-identity", type=float, default=90.0,
                   help='a hit below this counts as WEAK; the naming threshold is '
                        'a separate matter')
    p.add_argument("--species-threshold", choices=["apply", "ignore"],
                   default="apply", dest="species_threshold",
                   help='apply the locus specific threshold to the species name '
                        '(default: apply)')
    p.add_argument("--threads", type=int, default=4)
    # --databases: all (the default) asks every rDNA database of every class;
    # narrow asks only the discriminating set of the class, which is faster and
    # reproduces what the earlier version did.
    p.add_argument("--databases", choices=["all", "narrow"], default="all",
                   help='ask every database, or only the discriminating set of the '
                        'class (default: all)')
    # --class-db sets by hand which database a class is asked of.
    # The form: "A2=SILVA_138.2_LSURef_NR99.fasta:23S,B=..."
    # WHY IT IS THERE (2026-08-25): 16S cannot separate close species in principle
    # (M. mazei against M. soligelidi, 0.08 per cent apart). The reads of the A2
    # library are a whole operon of 4,300 bp and the LSU fragment is far more
    # discriminating. This option is for asking that fragment.
    p.add_argument("--class-db", default=None, dest="class_db",
                   help='a list of class=file:region, which overrides the default')
    # --fungal-rule: how the three loci of a fungal bin are combined.
    #   ladder (the default) is the method as reported: ITS first, then 28S, then
    #     18S, and the first locus that answers supplies the name.
    #   vote uses locus_decision.py: every locus is decided on its own and the
    #     genus is settled by a vote across them, with the species taken only from
    #     ITS or 28S. MEASURED (2026-08-27, 44 fungal bins): asking the three loci
    #     separately produced information ITS alone could not give in 16 of them.
    # The ladder stays the default because it is the method that was reported;
    # changing the default would make the produced result and the reported result
    # part company in silence.
    p.add_argument("--fungal-rule", choices=["ladder", "vote"], default="ladder",
                   dest="fungal_rule",
                   help='how the three loci of a fungal bin are combined '
                        '(default: ladder, the method as reported)')
    p.add_argument("--checkpoint", default=None,
                   help='the checkpoint file (default: <out directory>/'
                        'checkpoint/target_identity.json)')
    return p.parse_args()


def oku_fasta(f):
    return "".join(l.strip() for l in open(f, encoding="utf-8", errors="replace")
                   if not l.startswith(">")).upper()


def db_hazirla(fna, calisma):
    if os.path.exists(fna + ".nin") or os.path.exists(fna + ".00.nin"):
        return fna
    hedef = os.path.join(calisma, os.path.basename(fna))
    if not os.path.exists(hedef):
        os.symlink(os.path.abspath(fna), hedef)
    if not os.path.exists(hedef + ".nin"):
        r = subprocess.run(["makeblastdb", "-in", hedef, "-dbtype", "nucl"],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(u'   makeblastdb failed: %s' % r.stderr.strip()[:180])
            return None
    return hedef


def sinif_db_coz(metin):
    """A2=file.fasta:23S,B=... -> {class: [(file, region)]}"""
    out = {}
    for parca in (metin or "").split(","):
        parca = parca.strip()
        if not parca or "=" not in parca:
            continue
        sinif, kalan = parca.split("=", 1)
        dosya, _, bolge = kalan.partition(":")
        out.setdefault(sinif.strip(), []).append(
            (dosya.strip(), (bolge or "16S").strip()))
    return out


def adsiz_mi(ad):
    return bool(ad) and any(j in ad.lower() for j in ADSIZ_JETONLARI)


def basamaktan_sec(lokus_isabet, basamak, esik_uygula_mi):
    """Walk the ladder of the class and take the name from the first locus that
    gives a usable answer.

    lokus_isabet: {region: [(pid, aln, tit), ...]} already filtered and sorted by
    (-pid, -aln).

    Returns (name, note, identity, alignment, region, database) or None.
    """
    for bolge in basamak:
        isb = lokus_isabet.get(bolge) or []
        if not isb:
            continue
        pid, aln, tit, dbad = isb[0]
        adi = ad_ayikla(tit)
        if adi and adsiz_mi(adi):
            adi = None
        if not adi:
            # 2026-08-26: IF THE BEST HIT GIVES NO NAME, fall back to the first
            # NAMED hit within the separation margin. After SILVA was added, the
            # best hit in many bins became an unnamed environmental record and the
            # bin was left with no name at all, while a named record sat at almost
            # the same closeness. The same rule is in identity_verification.py.
            yakin = None
            for rpid, raln, rtit, rdb in isb[1:]:
                rad = ad_ayikla(rtit)
                if rad and not adsiz_mi(rad) and (pid - rpid) <= AYRIM_PAYI:
                    yakin = (rpid, raln, rad, rtit, rdb)
                    break
            if not yakin:
                continue
            rpid, raln, rad, tit, dbad = yakin
            adi, pid, aln = rad, rpid, raln
            on_not = (u'the best hit gives no name, so a named record at %.2f per '
                      u'cent was used' % rpid)
        else:
            on_not = u''
        lokus = lokus_duzelt(bolge, tit)
        # A rival counts only when it is A SPECIES DIFFERENT from ours. Comparing
        # a binomial against a family name says nothing about whether the species
        # can be defended, and putting "cf." on the strength of one would be
        # inventing an uncertainty that was not measured.
        rakip = None
        if ' ' in (adi or ''):
            for rpid, _raln, rtit, _rdb in isb:
                rad = ad_ayikla(rtit)
                if (rad and not adsiz_mi(rad) and rad != adi and ' ' in rad
                        and not ust_rank_mi(rad)):
                    rakip = (rpid, rad)
                    break
        if esik_uygula_mi:
            adi, aciklama = esik_uygula(adi, pid, lokus, aln, rakip)
        else:
            aciklama = u''
        if on_not:
            aciklama = (on_not + u'; ' + aciklama) if aciklama else on_not
        if adi == u'cannot be named':
            # This locus could not name it; the ladder carries on to the next one.
            # The reason is kept in case no locus manages it.
            continue
        # The ladder stops at the first locus that answers, which is the method as
        # reported. But when a LATER locus reaches species level and this one did
        # not, the reader has to know: the answer taken is the shallower of the
        # two, on purpose, and hiding that would make the row look more settled
        # than the data is.
        if ' ' not in adi or adi.endswith(' sp.'):
            derin = []
            for b2 in basamak:
                if b2 == bolge:
                    continue
                i2 = lokus_isabet.get(b2) or []
                if not i2:
                    continue
                p2, a2, t2, _d2 = i2[0]
                ad2 = ad_ayikla(t2)
                if not ad2 or adsiz_mi(ad2) or ust_rank_mi(ad2) or ' ' not in ad2:
                    continue
                l2 = lokus_duzelt(b2, t2)
                if (p2 >= TUR_ESIGI.get(l2, 98.7)
                        and a2 >= EN_AZ_HIZALAMA.get(l2, 600)):
                    derin.append(u'%s says %s at %.2f per cent over %d bp'
                                 % (b2, ad2, p2, a2))
            if derin:
                ek = (u'a deeper answer exists at another locus but the ladder of the '
                      u'class takes %s first: %s' % (bolge, '; '.join(derin)))
                aciklama = (aciklama + u'; ' + ek) if aciklama else ek
        return adi, aciklama, pid, aln, bolge, dbad
    return None


def oyla_mantar(lokus_isabet):
    """The three loci of a fungal bin, combined by locus_decision.py.

    Returns the same shape basamaktan_sec does, or None when no locus answered.
    lokus_karari sorts by (identity, then the third field); the third field is the
    alignment length here, because that is the tie break this project uses, and
    the bitscore is not carried this far.
    """
    kararlar = {}
    for bolge in (u'ITS', u'28S', u'18S'):
        isb = lokus_isabet.get(bolge) or []
        if not isb:
            continue
        kararlar[bolge] = _LD.lokus_karari(
            [(pid, aln, aln, 0, tit) for pid, aln, tit, _db in isb],
            lokus_duzelt(bolge, isb[0][2]), ad_ayikla, cins_epitet,
            ADSIZ_JETONLARI)
    if not kararlar:
        return None
    adi, _kimlik, notu = _LD.birlestir(kararlar)
    if not adi or adi == u'cannot be named':
        return None
    # Report the locus the name came from, not the one that happened to be first.
    en = max((v for v in kararlar.values() if v['cins']),
             key=lambda v: v['hizalama'], default=None)
    bolge = next((k for k, v in kararlar.items() if v is en), u'ITS')
    isb = lokus_isabet.get(bolge) or [(0.0, 0, '', '')]
    return (adi, notu, en['kimlik'] if en else isb[0][0],
            en['hizalama'] if en else isb[0][1], bolge, isb[0][3])


def main():
    a = get_args()
    if not shutil.which("blastn"):
        sys.exit(u'blastn was not found. To install it: sudo apt-get install -y ncbi-blast+')
    global SINIF_DB
    SINIF_DB = dict(SINIF_DB_DAR if a.databases == "narrow" else SINIF_DB_GENIS)
    print(u'thresholds: %s' % _ESIK_KAYNAK)
    print(u'databases: %s' % ('the discriminating set of each class'
                              if a.databases == "narrow"
                              else 'every rDNA database, for every class'))
    if a.class_db:
        ezme = sinif_db_coz(a.class_db)
        for k, v in ezme.items():
            SINIF_DB[k] = v
        print(u'   the class databases were overridden: %s'
              % ", ".join("%s -> %s" % (k, ",".join(d for d, _ in v))
                          for k, v in sorted(ezme.items())))
    ad = {}
    if os.path.exists(a.names):
        for l in open(a.names, encoding="utf-8"):
            q = l.rstrip("\n").split("\t")
            if len(q) > 1:
                ad[q[0]] = q[1]
    hedef_taxid = {}
    for l in open(a.targets, encoding="utf-8"):
        if l.startswith("#") or not l.strip():
            continue
        q = l.rstrip("\n").split("\t")
        if len(q) < 4 or q[0] == "karar":
            continue
        hedef_taxid[q[1]] = [t.strip() for t in q[3].split(",") if t.strip()]

    # the bin inventory
    kutular = {}
    for p in sorted(glob.glob(os.path.join(a.consensus, "*.fasta"))):
        et = re.sub(r"_(baskin|ref|self)?_?konsensus\.fasta$", "",
                    os.path.basename(p))
        m = re.match(r"((?:A1|A2|B|F1|F2))-\d+_(\d+)$", et)
        if m:
            kutular[et] = (m.group(1), m.group(2), oku_fasta(p))
    if not kutular:
        sys.exit(u'no consensus found: %s' % a.consensus)
    print(u'bins: %d' % len(kutular))

    calisma = tempfile.mkdtemp(prefix="kimlik_")
    # The checkpoint holds the raw blastn output per class and database. The
    # verdict and the naming are DERIVED AGAIN on every run, so a change of
    # threshold does not repeat the scan. The signature carries the database
    # mapping and the blast parameters.
    kp = a.checkpoint or os.path.join(
        os.path.dirname(os.path.abspath(a.out)), "checkpoint",
        "target_identity.json")
    defter = Defter(kp, imza={"sinif_db": SINIF_DB, "db": a.db,
                              "blast": "blastn -evalue 1e-20 -max_target_seqs 50",
                              "min_hizalama": a.min_alignment, "surum": "v2"})

    # one blastn call per class and database
    # kutu_lokus: {bin: {region: [(pid, aln, title, database), ...]}}
    kutu_lokus = collections.defaultdict(lambda: collections.defaultdict(list))
    zayif_en = {}
    for sinif in sorted(set(v[0] for v in kutular.values())):
        etler = [e for e, v in kutular.items() if v[0] == sinif]
        sorgu = os.path.join(calisma, "%s.fa" % sinif)
        with open(sorgu, "w", encoding="utf-8") as fh:
            for e in etler:
                fh.write(">%s\n%s\n" % (e, kutular[e][2].replace("N", "")))
        for dbad, bolge in SINIF_DB.get(sinif, []):
            fna = os.path.join(a.db, dbad)
            if not os.path.exists(fna):
                print(u'   no such database: %s' % dbad)
                continue
            cikti = os.path.join(calisma, "%s.%s.tsv" % (sinif, dbad))
            # THE CHECKPOINT IS CONSULTED FIRST: when the hits are cached both
            # blastn and makeblastdb are skipped. makeblastdb is more expensive
            # than blastn on some databases, and doing it the other way round
            # would eat most of the benefit of the rescue.
            ck = u"%s|%s|%s" % (sinif, dbad, a.min_alignment)
            onbellekten = defter.al(ck) if defter.var(ck) else None
            if onbellekten is not None:
                open(cikti, "w", encoding="utf-8").write(onbellekten)
                print(u'   cache  %-3s x %-32s (%d bins)' % (sinif, dbad, len(etler)))
                r = None
            else:
                db = db_hazirla(fna, calisma)
                if not db:
                    continue
                r = subprocess.run(
                    ["blastn", "-query", sorgu, "-db", db, "-outfmt",
                     "6 qseqid pident length bitscore qlen stitle",
                     # 5 -> 50 (2026-08-25). max_target_seqs cuts by BLAST's OWN
                     # bitscore ordering. Because bitscore grows with length, a
                     # short but exactly matching type strain record could not get
                     # into the first 5 and never appeared in the output file at
                     # all. MEASURED (A1-2_2209): the NR_041956.1 M. mazei DSM 2053
                     # record at 100.000 per cent over 1333 bp was not in the first
                     # 5, and a record at 98.591 per cent over 1419 bp was there
                     # instead. Choosing well is only useful if the right candidate
                     # is in hand, so the number of candidates fetched was raised.
                     "-max_target_seqs", "50",
                     "-evalue", "1e-20", "-num_threads", str(a.threads),
                     "-out", cikti], capture_output=True, text=True)
                print(u'   blastn %-3s x %-32s (%d bins)' % (sinif, dbad, len(etler)))
            if r is not None:
                if r.returncode != 0:
                    print(u'      ERROR: %s' % r.stderr.strip()[:160])
                    continue
                # Only a SUCCESSFUL run is stored. Storing the empty output of a
                # failed run would silently answer "no hits" on later runs.
                try:
                    defter.yaz(ck, open(cikti, encoding="utf-8").read())
                except (IOError, OSError) as e:
                    print(u'      the checkpoint could not be written to: %s'
                          % type(e).__name__)
            for line in open(cikti, encoding="utf-8"):
                q = line.rstrip("\n").split("\t")
                if len(q) < 6:
                    continue
                et, pid, aln, bit, qlen, tit = (q[0], float(q[1]), int(q[2]),
                                                float(q[3]), int(q[4]), q[5])
                if aln < a.min_alignment or pid < a.min_identity:
                    # The best hit that fails the filter is kept too. "No match"
                    # and "the nearest relative is at 88 per cent" are not the same
                    # thing; the first is an absence of data, the second an
                    # organism with no close relative in the database.
                    if et not in zayif_en or bit > zayif_en[et][0]:
                        zayif_en[et] = (bit, pid, aln, tit, bolge, dbad)
                    continue
                kutu_lokus[et][bolge].append((pid, aln, tit, dbad))

    shutil.rmtree(calisma, ignore_errors=True)

    # WITHIN a locus, order by identity and break ties by the longer alignment.
    # ACROSS loci nothing is ordered: the ladder decides which locus is used.
    for et in kutu_lokus:
        for bolge in kutu_lokus[et]:
            kutu_lokus[et][bolge].sort(key=lambda x: (-x[0], -x[1]))

    # collect per target
    sonuc = []
    for hedef, tl in sorted(hedef_taxid.items()):
        ilgili = [e for e, v in kutular.items() if v[1] in tl]
        if not ilgili:
            continue
        kraken = sorted(set(ad.get(v, v) for v in tl
                            if any(kutular[e][1] == v for e in ilgili)))
        say = collections.Counter()
        detay = {}
        lokus_of = {}
        for e in ilgili:
            sinif = kutular[e][0]
            basamak = SINIF_BASAMAK.get(sinif, VARSAYILAN_BASAMAK)
            if a.fungal_rule == "vote" and sinif in ("F1", "F2") and _LD:
                secim = oyla_mantar(kutu_lokus.get(e) or {})
            else:
                secim = None
            if secim is None:
                secim = basamaktan_sec(kutu_lokus.get(e) or {}, basamak,
                                       a.species_threshold == "apply")
            if secim:
                adi, aciklama, pid, aln, bolge, dbad = secim
                guclu = True
            else:
                z = zayif_en.get(e)
                if not z:
                    say[u"no hit at all in the databases"] += 1
                    continue
                _bit, pid, aln, tit, bolge, dbad = z
                # The name is NOT replaced by the RAW HEADER when it cannot be
                # extracted. It used to read "adi = ad_ayikla(tit) or tit[:40]";
                # with no name the header stood in for one, and then the genus was
                # taken as adi.split()[0], so AN ACCESSION NUMBER was mistaken for
                # a genus ("NG_064074.1 sp."). Such a record now counts as
                # unnameable; the header is in the evidence text anyway, so nothing
                # is lost.
                adi = ad_ayikla(tit)
                if not adi or adsiz_mi(adi):
                    adi = u'cannot be named'
                adi = adi + " (below the threshold)"
                aciklama = u''
                guclu = False
            say[adi] += 1
            lokus_of.setdefault(adi, (bolge, dbad))
            detay.setdefault(adi, []).append(
                "%s %.2f per cent/%dbp/%s/%s%s%s"
                % (e, pid, aln, bolge, dbad,
                   "" if guclu else " WEAK",
                   (" [" + aciklama + "]") if aciklama else ""))
        if not say:
            continue
        baskin, n = say.most_common(1)[0]
        bolge, dbad = lokus_of.get(baskin, ("", ""))
        # The agreement is reported IN LEVELS rather than as a yes or no: the genus
        # agreeing while the species differs is not the same as no agreement at all.
        if "no hit" in baskin:
            uyum = "vurus_yok"
        elif baskin.endswith("(below the threshold)"):
            uyum = "YAKIN_AKRABA_YOK"
        elif baskin.startswith(u"cannot be named"):
            uyum = "ADLANDIRILAMIYOR"
        elif baskin.endswith(" sp."):
            # We stayed at genus level because the species threshold was not
            # passed. Kraken may well be giving a species name, but WE CANNOT CLAIM
            # A SPECIES; writing that down as "uyusuyor" would say more than was
            # measured.
            uyum = ("cins_uyusuyor_tur_ATANAMAZ"
                    if any(baskin.split()[0] == k.split()[0] for k in kraken if k)
                    else "CINS_FARKLI")
        elif any(baskin == k for k in kraken):
            uyum = "tur_uyusuyor"
        elif any(baskin.split()[0] == k.split()[0] for k in kraken if k):
            uyum = "cins_uyusuyor_tur_farkli"
        else:
            uyum = "CINS_FARKLI"
        sonuc.append(dict(
            hedef=hedef,
            kraken_etiketi="; ".join(kraken),
            olculen_kimlik=baskin,
            kutu_sayisi=len(ilgili),
            destekleyen_kutu=n,
            uyum=uyum,
            lokus=bolge,
            veritabani=dbad,
            kanit="; ".join(detay.get(baskin, [])[:4]),
            diger=("; ".join("%s(%d)" % (k, v)
                             for k, v in say.most_common()[1:4]) or "")))

    if not sonuc:
        sys.exit(u'no identity could be measured for any target')
    d = os.path.dirname(os.path.abspath(a.out))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sonuc[0].keys()),
                           delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(sonuc)
    defter.kapat()
    print(u'\nwritten: %s' % a.out)
    say_uyum = collections.Counter(x["uyum"] for x in sonuc)
    print(u'targets: %d' % len(sonuc))
    for k in ("tur_uyusuyor", "cins_uyusuyor_tur_farkli",
              "cins_uyusuyor_tur_ATANAMAZ", "ADLANDIRILAMIYOR", "CINS_FARKLI",
              "YAKIN_AKRABA_YOK", "vurus_yok"):
        if say_uyum.get(k):
            print("   %-26s %d" % (k, say_uyum[k]))
    print("\n%-34s %-30s %-34s %-8s %s"
          % ("TARGET", "KRAKEN2 LABEL", "MEASURED IDENTITY", "LOCUS", "AGREEMENT"))
    for x in sonuc:
        print("%-34s %-30s %-34s %-8s %s"
              % (x["hedef"][:33], x["kraken_etiketi"][:29],
                 x["olculen_kimlik"][:33], x["lokus"], x["uyum"]))


if __name__ == "__main__":
    main()
