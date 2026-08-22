# -*- coding: utf-8 -*-
"""Configuration: every path and every numeric constant, in one place.

THIS IS THE FILE YOU EDIT to run the pipeline on your own data. It is the only
place paths are defined; nothing else hard-codes a location.

Contents, in order:
  * paths     panel tables, canonical consensus, reads, reference databases,
              outputs, checkpoints, cache
  * qPCR      product size windows, annealing temperature, structure limits.
              FIXED for QIAGEN Rotor-Gene Q + QuantiNova SYBR Green. Change
              them only if your chemistry differs, and change them BEFORE a
              run, never in the middle of one.
  * oligo     invariant primer rules (length, GC, 3' end, base runs)
  * grid      axes of the 144-cell parameter grid
  * sampling  in-sample measurement settings and funnel capacities

WARNING: one threshold here affects every stage at once. ENKOTU_ASGARI_OKUMA
and KAPSAM_ESIGI in particular carry comparability BETWEEN stages; changing
either mid-run makes stages incomparable without any error being raised.

Turkish identifiers are kept deliberately: several are also TSV column names
and checkpoint keys, so renaming them would change output schemas.

--- ozgun aciklama ---
Butun sabitler ve yollar tek yerde. Kullanici bu dosyayi duzenleyebilir.
qPCR kisitlari QIAGEN Rotor-Gene Q + QuantiNova SYBR Green icin SABIT tutuldu.
"""
# -------------------------------------------------------------------------
# config.py is the one place every file path and numeric constant is defined.
#
# INPUT  : it reads only its own location; KOK is the directory above this package.
#          It reads no other file and makes no measurement.
# OUTPUT : it writes no file. Constants are exposed at module level: the paths
#          (PANEL_TSV, HEDEFLER_TSV, KONSENSUS_KANONIK, FASTQ, REFDB, CIKTI,
#          KONTROL, ONBELLEK), the qPCR constraints, the fixed primer rules, the
#          axes of the 144 cell parameter grid, the sample measurement settings and
#          the funnel capacities.
# CALLED BY: nearly every module in the package imports it with
#          "from . import config as C"; from outside,
#          protocol/single_protocol_measure.py (key P),
#          verification/recovery_round.py (key K) and
#          verification/specificity_round.py (key D) use it too. So it is loaded
#          indirectly on every menu key.
#
# Changing a threshold here affects every stage at once. The ENKOTU_ASGARI_OKUMA
# and KAPSAM_ESIGI values in particular carry comparability between the stages and
# must not be changed in the middle of a run.
# -------------------------------------------------------------------------
import os

# ---------------------------------------------------------------- yollar
# KOK = proje klasoru (bu paketin bir ust dizini)
KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def y(*p):
    return os.path.join(KOK, *p)

PANEL_TSV   = y('primer_final', 'devir_ciftleri_20260802_sonrotus_TESLIM.tsv')
PANEL_TSV_YEDEK = y('primer_final', 'devir_ciftleri_20260802_sonrotus.tsv')
HEDEFLER_TSV= y('steps', 'hedefler.tsv')
TAXID_ADLARI= y('steps', 'taxid_adlari.tsv')
# ORIENTATION NORMALISATION (2026-08-02): consensuses are now read from ONE
# CANONICAL directory. 'consensus sequences' is MIXED orientation (71 antisense /
# 27 sense, measured) and reading it directly is forbidden; on a reversed consensus
# in-silico PCR SILENTLY gives 0 products (the measured loss is 100%,
# screening/orientation_impact_test.py).
# The canonical directory is produced by screening/build_canonical.py. The
# definition: orientation.py
KONSENSUS_KANONIK = y('konsensus_kanonik')
KONSENSUS_INDEKS  = os.path.join(KONSENSUS_KANONIK, 'INDEKS.tsv')
KONSENSUS_HAM     = y('consensus sequences')      # for canonical generation ONLY
KONSENSUS         = KONSENSUS_KANONIK
FASTQ       = y('fastq files')
REFDB       = y('REFERANS_DB')
SILVA_SSU   = os.path.join(REFDB, 'SILVA_138.2_SSURef_NR99.fasta')
SILVA_LSU   = os.path.join(REFDB, 'SILVA_138.2_LSURef_NR99.fasta')
UNITE_ITS   = os.path.join(REFDB, 'UNITE_ITS.fasta')
PR2         = os.path.join(REFDB, 'PR2_SSU_taxo_long.fasta')

# the directories holding the existing measurement code (IMPORTED, not rewritten)
# The core engine modules (ispcr, okuma, tarayici, cift). screening/engine_gateway.py
# loads them from a file path; they are kept outside the package because in the
# source study they lived in separate script directories. In this repository they
# were brought into one place.
BETIK_YOLLARI = [y('engine')]

CIKTI       = y('SCREENING_RESULT')
KONTROL     = os.path.join(CIKTI, 'kontrol')      # checkpoint
ONBELLEK    = os.path.join(CIKTI, 'onbellek')     # cache

# ---------------------------------------------------------------- qPCR kisitlari (SABIT)
# QuantiNova SYBR Green + Rotor-Gene Q
URUN_IDEAL      = (60, 150)   # tercih edilen
URUN_KABUL      = (150, 250)  # kabul edilebilir, protokolde 30 sn annealing/extension
URUN_MUTLAK_UST = 400         # aramada uretilen en buyuk urun (400 ustu hic denenmez)
URUN_ONERILMEZ  = 250         # above this it is marked "not recommended"

TA_HEDEF        = 60.0        # the aim is for the whole panel to run at one Ta, 60 C preferred
TA_KURALI       = 3.0         # Ta = min(Tm) - 3  (panelin kurali)

# SYBR Green -> the structural criteria are ELIMINATING (not a warning)
HAIRPIN_TM_UST     = 45.0
HOMODIMER_TM_UST   = 45.0
HETERODIMER_TM_UST = 45.0
# the delta-G thresholds (kcal/mol, computed at 60 C); more negative = worse
DG_YAPI_ALT        = -9.0     # herhangi bir yapi
DG_UC_ALT          = -5.0     # 3' uc kaynakli dimer

# the primer3 salt and concentration conditions - THE SAME VALUES AS THE PANEL'S
# (identical to geometry_core.py)
P3 = dict(mv_conc=50, dv_conc=1.5, dntp_conc=0.6, dna_conc=50)

# ---------------------------------------------------------------- degismez primer kurallari
UZUNLUK    = (18, 25)   # kullanicinin istegi: her pozisyondan 18-25 arasi her oligo
TEKRAR_UST = 4          # a run of 4 identical bases is forbidden (the same as ara.py)
DTM_UST    = 1.5        # cift ici Tm farki

# ---------------------------------------------------------------- THE PARAMETER GRID
# It starts strict and loosens step by step. THE ORDER MATTERS: the report answers
# "at which setting did an answer appear" against this order.
IZGARA_GC    = [(40, 60), (37, 63), (35, 65)]
IZGARA_TM    = [(58, 62), (57, 63), (56, 64)]
IZGARA_URUN  = [(60, 150), (60, 200), (60, 250), (60, 400)]
IZGARA_UC_GC = [True, False]    # 3' son baz G/C sart mi
IZGARA_SON5  = [True, False]    # 3' son 5 bazda en cok 3 G/C sart mi
# 3*3*4*2*2 = 144 cells in total

# ---------------------------------------------------------------- numune olcumu
NUMUNE_OKUMA_MIN, NUMUNE_OKUMA_MAX = 200, 6000   # the corrected filter (see the measurement-bug note #2)
NUMUNE_OKUMA_SAYISI = 300      # sample reads per bin (a fixed seed)
NUMUNE_TOHUM        = 20260802
NUMUNE_MAX_MM       = 1        # <=1 uyumsuzluk + 3' son 2 baz TAM (panel numune olcutu)
KURESEL_MAX_MM      = 5        # the global criterion (panel: total <=5, F and R written separately)
REFERANS_MAX_MM     = 1
# This threshold has to be ABSOLUTE. When a relative threshold was tried ("half the
# largest bin"), at full depth the largest bin is ~46 000 reads, so the threshold
# rose to ~23 000 and only 1-5 of the 10-33 competitor bins entered the measurement.
# In other words "the worst competitor bin" came to mean "the deepest bin" and the
# real worst competitor stayed outside the measurement. Measuring the same pair with
# 300 reads let every bin in, and a spurious difference of up to 40 fold appeared
# between the two stages.
ENKOTU_ASGARI_OKUMA = 150      # to enter the "worst single competitor bin" measure
                               # the MINIMUM reads required. It has to be ABSOLUTE: if it
                               # shifts with depth the stages stop being comparable.
# THE COVERAGE axis stays the only meaningful measure when the discrimination ratio
# becomes undefined (on universal targets the competitor set approaches empty and
# the denominator goes to zero). Besides, a single member bin coming out empty zeroes
# the discrimination ratio; coverage is not affected by that.
KAPSAM_ESIGI        = 0.20     # a member bin counts as 'covered' at >=20% product

# ---------------------------------------------------------------- huni (funnel) kapasiteleri
# Her asama bir sonrakine kac aday gecirir. Buyutmek suresi uzatir.
HUNI = dict(
    cift_ust        = 400000,   # asama B: sayilan en fazla cift (izgara tablosu bunun uzerinden)
    numuneye_giden  = 1200,     # stage C: the candidate pairs going to the raw read scan
    referansa_giden = 120,      # stage D: the ones going to reference coverage and competitor separation
    arms_taban      = 25,       # the number of 'best' candidates ARMS is tried on
    arms_ust        = 400,      # numunede olculen en fazla ARMS varyanti
    kusele_giden    = 12,       # asama E: kuresel taramaya giden (EN PAHALI)
)

# the number of bases processed at a time in the global scan (the memory ceiling is
# ~6 times this in bytes)
KURESEL_PARCA = 40_000_000


# -------------------------------------------------------------------------
# THE DISCRIMINATION THRESHOLD - ONE SOURCE
#
# 2026-08-06: the threshold is now defined not as a FOLD but as a dCq (delta Cq).
# The user's decision: dCq = 3.
#
# WHY dCq: that is the language the laboratory speaks. What qPCR measures is a cycle
# difference; the fold difference is DERIVED from it. With the threshold embedded as
# a fold ("10x"), the number looked arbitrary and could not be compared against the
# literature.
#
# THE CONVERSION: at 100% efficiency every cycle DOUBLES the product, so
#     fold = 2 ** dCq        and      dCq = log2(fold)
# dCq 3 -> 2**3 = 8.00 fold. (The old tool threshold, 10x, was dCq 3.32.)
#
# THE SOURCE: dCq >= 3 is the accepted floor in the literature for a specificity or
# NTC passing criterion (NEB's high efficiency qPCR data analysis guide and the MIQE
# reporting language). That number IS NO LONGER A TOOL THRESHOLD WE INVENTED; it is
# marked in the outputs as a "literature criterion".
#
# THE EFFICIENCY WARNING - TO KEEP IN MIND WHEN READING THE NUMBER:
# The conversion above assumes 100% efficiency. If the real efficiency is lower, THE
# SAME dCq corresponds to a SMALLER fold difference:
#     fold = (1 + E) ** dCq        E = the efficiency (0-1)
#     100% efficiency -> 3 cycles = 8.00 fold
#      90% efficiency -> 3 cycles = 6.86 fold
#      80% efficiency -> 3 cycles = 5.83 fold
# Once a calibration curve is produced in the laboratory, the real efficiency should
# be substituted and the rows past the threshold RE-EVALUATED. This note is printed
# into every output (ESIK_VERIM_NOTU).
# -------------------------------------------------------------------------
ESIK_DCQ = 3.0                      # <-- DEGISTIRILECEK TEK YER
AYRIM_ESIK = 2.0 ** ESIK_DCQ        # kat karsiligi: 8,00

ESIK_KOKENI = (u'LITERATUR OLCUTU dCq >= %.1f (ozgulluk/NTC gecme tabani; '
               u'NEB qPCR veri analizi, MIQE raporlama dili)' % ESIK_DCQ)

ESIK_VERIM_NOTU = (
    u'VERIM %%100 VARSAYILDI. dCq %.1f = %.2f kat donusumu her dongunun urunu '
    u'ikiye katladigini kabul eder. Gercek verim daha dusukse ayni dCq daha '
    u'kucuk kat farkina denk gelir: %%90 verimde %.2f kat, %%80 verimde %.2f '
    u'kat. Kalibrasyon egrisi cikinca gercek verimle yeniden degerlendirin.'
    % (ESIK_DCQ, 2.0 ** ESIK_DCQ, 1.9 ** ESIK_DCQ, 1.8 ** ESIK_DCQ))


def esik_metni(kat=None):
    """THE SAME format in every output: 'dCq 3,0 (8,00x)'."""
    k = AYRIM_ESIK if kat is None else kat
    return u'dCq %s (%sx)' % (('%.1f' % ESIK_DCQ).replace('.', ','),
                              ('%.2f' % k).replace('.', ','))


def kat_dcq(kat, verim=1.0):
    """Fold -> dCq. An efficiency of 1.0 is 100%. This is what gets used once the real efficiency is measured."""
    import math as _m
    try:
        k = float(kat)
    except (TypeError, ValueError):
        return None
    if k <= 0:
        return None
    return round(_m.log(k) / _m.log(1.0 + verim), 2)


def kat_ve_dcq(kat):
    """'8,45x (dCq 3,08)' - the fold and the dCq SIDE BY SIDE. Used on measurement rows."""
    d = kat_dcq(kat)
    if d is None:
        return u'-'
    return u'%sx (dCq %s)' % (('%.2f' % float(kat)).replace('.', ','),
                              ('%.2f' % d).replace('.', ','))
