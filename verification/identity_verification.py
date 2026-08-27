# -*- coding: utf-8 -*-
"""Identity verification, tests reported identity claims INDEPENDENTLY.

WHY THE METHOD IS DELIBERATELY DIFFERENT
    This round repeats NONE of the methods that produced the claims:

      * Kraken2       k-mers + lowest common ancestor on a taxonomy tree.
                      This module never touches a taxonomy tree.
      * our round 1   within-class consensus alignment + discriminating 21-mers.
                      This module compares consensus sequences against NAMED
                      records in external databases, not against each other.
      * our round 2   in-silico PCR. There are no primers in this module.

    Sharing a mechanism would mean sharing its blind spots: a wrong call would
    be confirmed rather than caught.

METHOD: SEED + ALIGNMENT (BLAST-like, taxonomy-free)
    1  k-mer seeds are taken from the query consensus; the database is streamed
       and short-listed by seed count.
    2  EVERY short-listed record is fully aligned (Levenshtein DP, infix).
    3  Identity is measured TWICE, over the whole overlap, and over the
       DISCRIMINATING WINDOW: the columns where the best reference records
       differ from one another. Conserved regions (18S, 5.8S, LSU core) fall
       outside that window, so a claim resting on falsely high conserved-region
       identity becomes visible instead of passing.
    4  AGREEMENT BETWEEN DATABASES IS REQUIRED. A single database's best hit is
       never an identification: deduplicated sets delete rare genera (measured:
       0 Petriella records in SILVA LSURef NR99, 82 in the Parc set of the same
       release).

VERDICTS
    DOGRULANDI      >=2 independent databases support the claim
    DUZELTILMELI    >=2 databases agree on a DIFFERENT answer; the corrected
                    wording is written out
    DOGRULANAMADI   evidence insufficient, contradictory, or single-source

    No verdict is invented. When evidence is thin the answer is
    DOGRULANAMADI, not a guess.

UNNAMED RECORDS CANNOT BECOME A NAME
    A 99% match to "Uncultured bacterium clone 4B-11" says your sequence
    overlaps environmental clones. It is not a species. Such records return
    ADLANDIRILAMIYOR (referans adsiz) and never a taxon, see ad_coz().

    Output also lists the five nearest organisms, deduplicated by ORGANISM
    rather than by record, so the list shows what else is close instead of the
    same species repeated from five databases.

    Writes to IDENTITY_RESULT/ only. Never touches panel files.

"""

# -------------------------------------------------------------------------
# identity_verification.py tests the reported identity claims independently
# against external reference databases, with SEEDING plus ALIGNMENT.
#
# INPUT  : the local FASTA sets under REFERENCE_DB/ (the VTB list; twins and
#          subsets have been taken out of the vote),
#          canonical_consensus/ (through screening.targets.konsensusler),
#          NCBI nt (a separate layer, over the network) and the hand filled
#          IDENTITY_RESULT/nt_elle/NT_SONUC_SABLONU.tsv.
# OUTPUT : IDENTITY_RESULT/kimlik_iddialari.tsv (the main table),
#          IDENTITY_RESULT/KIMLIK_DOGRULAMA_RAPORU.md,
#          IDENTITY_RESULT/VERITABANI_ENVANTERI.md,
#          IDENTITY_RESULT/LITERATUR_ELLE_KONTROL.tsv,
#          IDENTITY_RESULT/nt_ham/, nt_elle/, kontrol/ .
#          It WRITES NOTHING into the panel files.
# CALLED BY: verification/full_chain.py -> key I
#          (python3 verification/identity_verification.py --root .)
#          This file's functions are also imported as a module by stage G
#          (all_bin_identities.py) and stage E (access_check.py).
#
# THREE RULES, each marked in the code below:
#   1) THE SHORT LIST IS 500 and ALL of it is aligned, so the cut off is not binding.
#   2) A claim counts as VERIFIED only when AT LEAST TWO independent databases agree.
#   3) Identity is also measured in a DISCRIMINATING WINDOW, with the conserved
#      regions left out.
# -------------------------------------------------------------------------
import os, sys, csv, json, time, re, argparse

VERSIYON = '1.1 (2026-08-04)'

# --------------------------------------------------- THE SETTINGS SIGNATURE
# This is MIXED INTO the checkpoint key. Raise this line whenever the short list
# size or the ranking and selection logic changes: vtb_*.json files left over from
# older runs then become invalid AUTOMATICALLY, instead of silently handing back
# the old result produced with the smaller list.
# 2026-08-26: v3 -> v4. The hit ordering changed (on equal identity the LONGER
# alignment now comes first) and so did the naming rules (a minimum alignment
# length, the fungal thresholds, the unnamed second hit). Because the ordering
# affects the SCAN stage, the signature was raised: the old vtb_*.json files
# become invalid on their own and the scan runs from scratch. That way THE
# VERDICT is re-derived as well, not the naming alone.
# v4 -> v5 (2026-08-26): the independence measurement moved into the verdict
# path, and the nt search mode was changed from MEGABLAST to ordinary blastn.
# Both are visible only in A FULL SCAN; a re-derivation takes the verdict from
# the cache.
AYAR_IMZASI = 'idfbm25-kl1200-v5'

K_TOHUM = 16            # seed length

# ------------------------------------------- THE RANKING CRITERION (2026-08-05)
# THE OLD CRITERION: score = "how many distinct seeds of the query occur in this
# record" (a flat count). It broke in two ways at once, and both got worse as the
# database grew:
#
#   (a) CONSERVED REGION NOISE. Seeds coming from the 18S/5.8S/LSU core occur in
#       almost the whole database (measured: one seed occurred in 1 612 663 of
#       2 069 188 records). Those seeds discriminate nothing, yet under a flat
#       count they carry the SAME weight as a rare seed that is real evidence of
#       relatedness.
#   (b) LENGTH BIAS - THIS WAS THE REAL CAUSE. Because the score is a SUM, the
#       longer record always wins. The target record AY882347 is 484 bp in UNITE
#       and can hold only 19 of the query's seeds. The records that take the top
#       places instead are 892-2000 bp and hold 34-59 seeds. A short record that
#       matches EXACTLY gets buried purely for being SHORT.
#
# THE NEW CRITERION = inverse frequency weighting plus BM25 length normalisation:
#
#       score(record) = SUM_{matching seeds} ln(N / (1 + df(seed)))
#                       ----------------------------------------------
#                          1 - b + b * length(record) / MEAN_LENGTH
#
#   df(seed) : IN HOW MANY RECORDS that seed occurs (counted free in the same pass)
#   N        : the number of records scanned,  b = BM25_B = 0.75,  MEAN = mean length
#
# MEASURED (the rank of record AY882347, all against a FULL database scan):
#   UNITE ITS (2 069 188 records)  F2-1_101201    1869 -> 19      (idf alone: 1320)
#                                  F2-1_2034170   2037 -> 18
#                                  F1-4_101201     136 -> 17
#                                  F1-4_2093779    162 -> 18
#                                  F1-4_2093780    153 -> 17
#                                  F1-2_101201     189 -> 30
#                                  F1-1_2093779    225 -> 35
#   SILVA LSU Parc (1 312 521)     F2-1_101201    6794 ->  2
#                                  F2-1_2034170   5749 ->  3
#   RefSeq ref_all2 (65 358)       F2-1_101201    1320 ->  1
#   RefSeq fungal ITS (20 394)     1 -> 1 in all seven queries  (NO DEGRADATION)
#
# Inverse frequency ALONE WAS NOT ENOUGH (1869 -> 1320). Length normalisation is
# required; the two work together. Detail: SIRALAMA_COZUMU.md.
BM25_B = 0.75           # the BM25 length normalisation coefficient (0 = no normalisation)
ADAY_HAVUZU = 3000      # candidates kept during the pass; re-ranked with idf
                        # ranked, and the first KISA_LISTE of them go to alignment.
                        # The worst measured pre-ranking position is 45, so there is 66x of headroom.

# THE SHORT LIST SIZE - 2026-08-04, the second correction.
#
# THE PROBLEM (measured): the short list was ranked and cut by SEED COUNT, but the
# decision is made by ALIGNMENT IDENTITY. Because the two criteria differ, the real
# best match could fail to enter the list at all:
#   * F2-1_101201: the list showed Parascedosporium putredinis as best (98.16%),
#     while the same database held Petriella guttulata at 99.65%.
#   * SILVA LSU NR99, a query degraded by 8%: 47 of the target record's seeds were
#     intact and it still did not reach the first 60, because it was competing with
#     thousands of records sharing the same conserved regions.
# The "guarantee the expected taxon" patch only closes that hole WHEN WE ALREADY
# KNOW WHAT WE ARE LOOKING FOR, and in an identity question we do not.
#
# THE FIX: make the cut off NON-BINDING. The list went 60 -> 500 and ALL 500
# candidates are aligned; the decision is left entirely to alignment.
# The cost: the vectorised aligner takes ~0.02 s for 1.5 kb, so 500 candidates take
# ~10 s (the database pass of the short list stage already takes minutes, so the
# alignment cost is a small part of the total).
#
# SELF CALIBRATION: for every query the seed rank of the WINNING hit inside the
# short list is recorded (see kazanan_sira). If the winners always come from the
# first 100, the cut off is not binding, and we prove that WITHOUT A SEPARATE
# MEASUREMENT.
# 2026-08-05 (morning): 500 -> 1000. Self calibration showed 500 WAS BINDING (in 13
# of 118 queries the winner came from beyond position 400; the highest winning rank
# was 4171). 1000 recovered part of that but DID NOT RECOVER 4171.
#
# 2026-08-05 (evening): 1000 -> 500 AGAIN. Growing the list did not solve the
# problem, it only dug deeper. The problem was IN THE RANKING CRITERION (see "THE
# RANKING CRITERION" above: conserved region noise plus length bias). Under the new
# criterion the worst measured winning rank is 35; 500 is 14 times that and it
# halves the running time of stages I and G. Growing it is no longer necessary.
KISA_LISTE = 500        # how many records from each database are fully aligned
SIRA_UYARI_ESIGI = 200  # if the winner comes from beyond this, the criterion is still broken
SIRA_GUVENLI_BOLGE = 50   # if the winners stay inside this, the cut off is not binding
GARANTI_UST = 40        # the most records taken through the "expected taxon guarantee"
AYIRT_EDICI_UST = 8     # how many best records the discriminating window uses
UYUM_TOLERANS = 1.0     # the acceptable gap between the claimed percentage and the measured one

# (label, file, type, use_in_identity_stage, note)
# The rule for the IDENTITY stage: a DEREPLICATED set IS NOT ENOUGH. Sets such as
# NR99 delete rare genera outright. Measured: SILVA LSURef NR99 holds ZERO
# Petriella records, while the Parc set of the same release holds 82. For an
# identity question, Parc is REQUIRED.
VTB = [
    ('SILVA SSU NR99',     'SILVA_138.2_SSURef_NR99.fasta', 'SSU',    True,
     '510,495 records; SSU, deduplicated'),
    ('SILVA LSU NR99',     'SILVA_138.2_LSURef_NR99.fasta', 'LSU',    True,
     '95,279 records; LSU, deduplicated, which DELETES RARE GENERA'),
    ('SILVA LSU Parc',     'SILVA_138.2_LSUParc.fasta',     'LSU',    True,
     '1,312,521 records; NOT deduplicated. Petriella: 0 in NR99 and 82 in '
     'Parc, measured. REQUIRED at the identity stage.'),
    ('UNITE ITS',          'UNITE_ITS.fasta',               'ITS',    True,
     '2,069,189 records; fungal ITS. Petriella: 113 records, measured.'),
    ('PR2 SSU',            'PR2_SSU_taxo_long.fasta',       'SSU',    True,
     '240,201 records; eukaryotic 18S'),
    ('ROD operon',         'ROD_v1.2_operon_variants.fasta','OPERON', True,
     '60,320 records; rRNA operon variants'),
    ('RefSeq bakteri 16S', 'bacteria.16S.fna',              'SSU',    True,
     '26,877 records; weighted towards named type material'),
    ('RefSeq arke 16S',    'archaea.16S.fna',               'SSU',    True,
     '1,160 records'),
    ('RefSeq mantar ITS',  'fungi.ITS.fna',                 'ITS',    True,
     '20,394 records'),
    # 2026-08-26: the locus went from 'LSU' to 'LSU_MANTAR'. This set is fungal
    # throughout and the fungal LSU species threshold is 99.8 per cent (Vu et al.
    # 2018); the 98.7 per cent of the prokaryotic 23S would be too loose here.
    ('RefSeq mantar 28S',  'fungi.28SrRNA.fna',             'LSU_MANTAR', True,
     '12,890 records; Petriella: 2 records, measured'),
    ('RefSeq mantar 18S',  'fungi.18SrRNA.fna',             'SSU',    True,
     '4,037 records'),
    ('RefSeq ref_all2',    'ref_all2.fna',                  'KARISIK', True,
     '65,358 records; the merged RefSeq set, a superset of ref_all'),
    ('RefSeq ref_all',     'ref_all.fna',                   'KARISIK', False,
     '48,431 records; ref_all2 is a SUPERSET of this one, so it was taken out '
     'of the identity vote to stop the same records being counted twice. It '
     'is not an independent source.'),
]

# NCBI nt is a SEPARATE layer (not a local file, it goes over the network).
NT_ETIKET = 'NCBI nt'

# --------------------------------------------------------------- THE CLAIMS
# type: 'kimlik'  -> an identity claim between a bin and a named taxon
#       'kutu2'   -> the identity of two bins to one another
#       'dagilim' -> a claim that one organism is spread over several bins
#       'adsiz'   -> a "cannot be named" claim
#       'ayrilmaz'-> a claim that two species cannot be told apart
#       'gecici'  -> a "identity cannot be separated in this class" claim
IDDIALAR = [
 dict(no=1, oncelik=2, tip='kimlik', kutu=['F2-1_101201'], sinif='F2',
      beklenen_cins='Petriella', beklenen_yuzde=None,
      metin='the taxid 101201 bin is from the genus Petriella'),
 dict(no=2, oncelik=1, tip='kutu2', kutu=['F1-4_2093780'], karsi=['F2-1_101201'],
      beklenen_yuzde=99.58, beklenen_cins='Petriella',
      metin='the Podospora pseudopauciseta bin (F1-4_2093780) is 99.58 per cent identical to Petriella',
      not_='THE NEWEST AND LEAST TESTED CLAIM, so a high priority'),
 dict(no=3, oncelik=1, tip='dagilim', beklenen_cins='Petriella',
      kutu=['F2-1_101201','F2-2_101201','F2-3_101201','F2-4_101201',
            'F2-1_2034170','F2-4_2034170','F2-1_500148','F2-2_500148','F2-4_500148'],
      supheli=['F2-1_2034170','F2-4_2034170','F2-1_500148','F2-2_500148','F2-4_500148'],
      beklenen_alt=76.0, beklenen_ust=86.0,
      # THE 2026-08-06 CORRECTION: the text said "nine bins" and that number came from
      # no measurement; it was the length of the claim's own bin list.
      # Stage G scanned ALL 96 bins independently and measured this:
      #   * EIGHT of this claim's nine bins came out Petriella,
      #   * the ninth (F2-4_500148) DID NOT. It is an unnameable lineage whose nearest
      #     record is Lomentospora prolifica at 83.68%; the claim's own evidence row
      #     already had it diverging at 52.21%,
      #   * against that, FOUR F1 bins the claim did not count (F1-2_101201,
      #     F1-4_101201, F1-4_2093779, F1-4_2093780) also came out Petriella.
      # Twelve bins in total. The text now carries the measured number, so that no
      # contradiction between "nine" and the evidence reaches the report.
      metin='Petriella is spread over 12 bins in the sample: eight of this '
            "claim's nine bins plus four F1 bins, while F2-4_500148 is NOT "
            'Petriella. Of the reads in the bins Kraken called T. breve and '
            'M. brunneum, 76 to 86 per cent are Petriella.',
      not_='it rests on one round, so it takes priority'),
 dict(no=4, oncelik=3, tip='kimlik', kutu=['A1-4_2208'], sinif='A1',
      beklenen_tur='Methanosarcina vacuolata', beklenen_yuzde=97.6,
      metin='the M. barkeri bin is 97.6 per cent to M. vacuolata'),
 dict(no=5, oncelik=3, tip='kimlik', kutu=['B-2_818','B-3_818','B-2_214856'], sinif='B',
      beklenen_tur='Alistipes putredinis', beklenen_yuzde=85.0, tolerans=2.0,
      adsiz_bekleniyor=True,
      metin='the Bacteroides bins are an unnameable Bacteroidales whose '
            'nearest record, Alistipes putredinis, is around 85 per cent'),
 dict(no=6, oncelik=3, tip='kimlik', kutu=['B-2_1197717'], sinif='B',
      beklenen_yuzde=99.4, beklenen_aile='Synergistaceae',
      rakip_cins='Cloacibacillus', rakip_yuzde=90.0, adsiz_bekleniyor=True,
      metin='the Cloacibacillus target is an unnameable Synergistaceae at '
            '99.4 per cent, with the nearest Cloacibacillus at 90 per cent'),
 dict(no=7, oncelik=2, tip='adsiz', kutu=['F1-1_44689','F1-2_44689','F1-3_44689','F1-4_44689'],
      metin=u'taxid 44689 etiketi curutuldu ve yerine isim konulamiyor'),
 dict(no=8, oncelik=2, tip='ayrilmaz', kutu=['A1-2_2209'], sinif='A1',
      turler=[('Methanosarcina soligelidi', 99.93), ('Methanosarcina mazei', 99.85)],
      metin='the taxid 2209 bin is 99.93 per cent to M. soligelidi and 99.85 per cent to M. mazei; the two do not separate'),
 dict(no=9, oncelik=2, tip='cins_duzeyi', kutu=['F2-1_101201'], beklenen_cins='Petriella',
      metin=u'Petriella cins duzeyinde kalmali, tur adi verilemez, cf. setifera'),
 # --- The targets stage K calls "HETEROJEN": to decide anything, the bin
 # identities have to be tested AGAINST AN EXTERNAL REFERENCE. K cannot do that; I can.
 dict(no=11, oncelik=1, tip='kimlik', kutu=['B-4_285070'], sinif='B',
      beklenen_cins='Petrimonas',
      metin='Are the member bins of the Petrimonas target the same organism? '
            'The recovery stage called them HETEROGENEOUS: not one of the 3 '
            'bins is 99 per cent or more identical to another.',
      not_='The recovery stage COULD NOT narrow this target and said the bin '
           'identities have to be confirmed against a reference first. This '
           'claim is exactly that check, and without its result nothing can '
           'be decided about Petrimonas.'),
 dict(no=12, oncelik=1, tip='kimlik', kutu=['B-2_818'], sinif='B',
      beklenen_tur='Alistipes putredinis', beklenen_yuzde=85.0, tolerans=3.0,
      adsiz_bekleniyor=True,
      metin='Are the member bins of the Bacteroidales set target the same '
            'organism? The recovery stage called them HETEROGENEOUS: not one '
            'of the 12 bins is 99 per cent or more identical to another.',
      not_='The same reasoning: the recovery stage could not narrow it, so confirmation against an outside reference is required.'),
 dict(no=10, oncelik=2, tip='gecici', sinif='B',
      metin='the class B bin identities are provisional, because consensus '
            'identity cannot separate anything in that class'),
]


# --------------------------------------------------------------- basics
_C = str.maketrans('ACGTUNacgtun', 'TGCAANtgcaan')


# Reverse complement. The query is searched in both directions: the orientation of
# reference records varies from set to set, and searching one direction loses half.
def rc(s):
    return s.translate(_C)[::-1]


# Everything outside ACGT becomes N (U -> T). Ns always count as a mismatch in the
# alignment, so an ambiguous base is never read in the candidate's favour.
def temizle(s):
    return re.sub(r'[^ACGT]', 'N', s.upper().replace('U', 'T'))


def sure_metni(sn):
    sn = int(sn)
    return ('%d saniye' % sn) if sn < 90 else ('%d dakika' % round(sn / 60.0)) \
        if sn < 5400 else ('%.1f saat' % (sn / 3600.0))


# The decimal separator. None and values that cannot be converted become "-", so
# that "0" and "not measured" are never confused inside a report.
def vir(x, b=2):
    if x is None:
        return '-'
    try:
        return ('%.*f' % (b, float(x))).replace('.', ',')
    except (TypeError, ValueError):
        return str(x)


_M = None


def enc(s):
    global _M
    import numpy as np
    if _M is None:
        m = np.full(256, 4, dtype=np.uint8)
        for i, c in enumerate('ACGT'):
            m[ord(c)] = i
        _M = m
    return _M[np.frombuffer(s.encode(), dtype=np.uint8)]


# -------------------------------------------------------------------------
# INFIX (HW) LEVENSHTEIN - THE MEASURE THAT MAKES THE DECISION.
#
# The short query is aligned INSIDE the long target, so the overhang left at the
# start and the end of the target is not penalised. Reference records come in very
# different lengths (some a full operon, some a single locus); global alignment
# counts that length difference as mismatch and rules the right record out.
#
# THIS FUNCTION MAKES THE DECISION, not the seed count. Seed ranking is only the
# criterion for COLLECTING CANDIDATES; ALL 500 candidates in the short list are put
# through here and the best hit is chosen by alignment identity. The two criteria
# being separate is what caused the real best matches to be missed back when the
# list was 60.
# -------------------------------------------------------------------------
def hizala(q, t):
    """Infix (HW) Levenshtein: aligns the short query INSIDE the long target.
        Returns: (percent_identity, distance). Row by row DP with numpy.

    """
    import numpy as np
    if not q or not t:
        return (0.0, len(q or t or ' '))
    Q, T = enc(q), enc(t)
    onceki = np.zeros(len(T) + 1, dtype=np.int32)
    for i in range(len(Q)):
        simdi = np.empty_like(onceki)
        simdi[0] = i + 1
        farkli = (T != Q[i]) | (T == 4) | (Q[i] == 4)
        # The left neighbour dependency (insertion) is VECTORISED:
        #   now[j] = min(cand[j], now[j-1]+1)
        # Setting a[j] = now[j]-j gives a[j] = min(cand[j]-j, a[j-1]), which is a running
        # minimum, so np.minimum.accumulate does it in one pass. The Python inner loop is
        # gone: a 1.5 kb x 1.5 kb alignment drops from minutes to seconds.
        aday = np.minimum(onceki[:-1] + farkli, onceki[1:] + 1)
        aday = np.concatenate(([i + 1], aday))
        idx = np.arange(len(aday))
        simdi = np.minimum.accumulate(aday - idx) + idx
        onceki = simdi
    d = int(onceki.min())
    return (round(100.0 * (1 - d / float(len(q))), 2), d)


def fasta_akisi(yol):
    """Yields (header, sequence). A stream, for files that do not fit in memory.

        THERE IS NO CAP; it reads to the END of the file. (There used to be an
        unused 'parca=200000' parameter. It did nothing anywhere, but it gave the
        impression that a cap existed, so it was removed. A real cap problem did
        happen in the access test: the first run cut off at 120 001 records, and
        SILVA SSU NR99 (510 495), LSU Parc (1 312 521) and UNITE ITS (2 069 189)
        were effectively being scanned truncated. The identity stage does NOT have
        that bug; the number of records scanned is printed on every line and
        compared against the expected record count.)

    """
    bas, par = None, []
    with open(yol, encoding='utf-8', errors='ignore') as fh:
        for satir in fh:
            if satir.startswith('>'):
                if bas is not None:
                    yield bas, temizle(''.join(par))
                bas = satir[1:].strip(); par = []
            else:
                par.append(satir.strip())
    if bas is not None:
        yield bas, temizle(''.join(par))


# THE SEED STEP - chosen by measurement (2026-08-04).
# The short list stage scans the whole database pass, and its cost is proportional
# to THE SEED COUNT. Measured (SILVA LSU Parc, 2.1 GB / 1.31 M records):
#   step= 7 -> 410 seeds -> ~21 minutes
#   step=25 -> 116 seeds -> ~ 6 minutes   <- CHOSEN
#   step=60 ->  48 seeds -> ~ 3 minutes
# Whether the short list degrades was tested: at step=25, 23/25 and 17/25 of the
# first 25 stayed the same, and THE BEST HIT CAME OUT THE SAME IN BOTH TESTS. Since
# every verdict rests on hit[0], step=25 is safe; step=60 loses more.
def tohumlar(q, k=K_TOHUM, adim=25):
    return {q[i:i + k] for i in range(0, max(1, len(q) - k + 1), adim) if 'N' not in q[i:i + k]}


# -------------------------------------------------------------------------
# THE MEAN RECORD LENGTH - the denominator of the BM25 normalisation
#
# Because the normalisation rests on L/MEAN, the MEAN has to be known BEFORE THE
# PASS STARTS: otherwise the first records in the file are scored wrongly and are
# dropped from the candidate pool for the wrong reason (tried with a running mean;
# the ranking moves in the first 5%). It is computed once and written to a .ortuz
# file beside the database; later runs read it.
# -------------------------------------------------------------------------
def ortalama_uzunluk(yol):
    yan = yol + '.ortuz'
    try:
        if os.path.exists(yan) and os.path.getmtime(yan) >= os.path.getmtime(yol):
            n, t = open(yan).read().split()
            if int(n) > 0:
                return float(t) / int(n)
    # If the side file is corrupt or unreadable the length is counted again, so
    # we fall onto the right path anyway. But WHICH error is expected should be
    # written down: this used to be "except Exception", and a fault in the code
    # would land here and stay hidden.
    except (IOError, OSError, ValueError) as e:
        print(u'  [note] %s could not be read (%s), the length is being counted again'
              % (os.path.basename(yan), type(e).__name__))
    n = 0
    t = 0
    fai = yol + '.fai'
    if os.path.exists(fai):                  # samtools indeksi varsa bedava
        # THE A6 CORRECTION (2026-08-21): a malformed line used to be skipped silently
        # while the t/n accumulation carried on, so the mean was computed from PARTIAL
        # data. That value feeds the short list scoring. If there is a malformed line the
        # .fai is not trusted at all: an inconsistent index is thrown away entirely and the
        # full scan branch below takes over. A visible full measurement instead of a silent
        # partial one.
        fai_bozuk = 0
        for satir in open(fai, encoding='utf-8', errors='ignore'):
            p = satir.split('\t')
            if len(p) > 1:
                try:
                    t += int(p[1]); n += 1
                except ValueError:
                    fai_bozuk += 1
        if fai_bozuk:
            sys.stderr.write(
                u'WARNING: %s contains %d malformed lines. The index was ignored and the average length is measured by scanning the FASTA.\n' % (fai, fai_bozuk))
            n = 0
            t = 0
    if not n:                                # if absent, a fast pass that counts length alone
        u = 0
        with open(yol, 'rb') as f:
            for ham in f:
                if ham[:1] == b'>':
                    if u:
                        t += u; n += 1
                    u = 0
                else:
                    u += len(ham.strip())
        if u:
            t += u; n += 1
    ort = float(t) / max(1, n)
    try:
        open(yan, 'w').write('%d %d' % (n, t))
    except (IOError, OSError):
        # A read only directory or a full disk. The calculation has already been
        # made; it simply cannot be kept for the next run. Not worth stopping for.
        pass
    return ort if ort > 0 else 1.0


class _TersBaslik(object):
    """Baslik ARTAN sirali olsun diye min-heap'te ters karsilastirilan sarmalayici."""
    __slots__ = ('s',)

    def __init__(self, s):
        self.s = s

    def __lt__(self, o):
        return self.s > o.s

    def __eq__(self, o):
        return self.s == o.s

    def __str__(self):
        return self.s


# -------------------------------------------------------------------------
# THE SHORT LIST - SELECTION IN TWO STAGES
#
# The list used to be 60, then 500, then 1000. Growing it did not solve the
# problem: in one measured case the real relative (AY882347) sat at position 1869
# in UNITE ITS, and no reasonable list size would have caught it. The problem was
# not the list SIZE but the RANKING CRITERION (for the reasoning see the "THE
# RANKING CRITERION" section at the top of this file: conserved region noise plus
# length bias).
#
# There are now two stages inside one pass:
#   1) PRE-FILTER - the best ADAY_HAVUZU (3000) records are kept by RAW seed count
#      normalised for length. At this stage df is not yet known; but in the same
#      pass, how many records each seed occurs in is counted FOR FREE.
#   2) RE-RANKING - once the pass ends, idf is computed from the real df values,
#      the 3000 candidates are re-ranked by the idf+BM25 score, and the first
#      'ust' of them go to alignment. The worst measured pre-filter position is 45,
#      and since the pool is 3000, the pre-filter cut off is NOT binding.
#
# The "expected taxon guarantee" patch (the garanti parameter) stays where it is,
# but IT IS NO LONGER NEEDED: in all seven measured queries the winner enters the
# list on its seeds. Where a record entered through the patch, its position under
# normal ranking is still computed, and if the winner came in through the patch the
# report prints a WARNING.
# -------------------------------------------------------------------------
def kisa_liste(yol, q, ust=KISA_LISTE, ilerle=None, garanti=(), havuz=None,
               suzgec=None):
    """Pick the most promising records (BLAST's seeding step), with idf plus BM25.

        Returns: [dict(tohum, skor, baslik, dizi, sira, kaynak), ...]
          tohum  : how many distinct query seeds occurred in this record
                   (the raw count, for information)
          skor   : the idf weighted, length normalised ranking score
          sira   : the place in the FINAL ranking (1 = best). The winner's value
                   is what gets recorded.
          kaynak : 'tohum' -> it came in on the criterion | 'garanti' -> through
                   the taxon patch

    """
    import heapq, math
    th_l = sorted(tohumlar(q) | tohumlar(rc(q)))
    if not th_l:
        return []
    nt = len(th_l)
    ORT = ortalama_uzunluk(yol)
    # ust=0 KEEPS ITS OLD MEANING: no cut, EVERY record holding a seed is aligned.
    # In that case the pre-filter pool must also be unlimited, otherwise "no cut"
    # would quietly turn into "cut at 3000".
    havuz = None if not ust else max(int(havuz or ADAY_HAVUZU), ust * 3, 500)
    df = [0] * nt
    N = 0
    yigin = []                     # min-heap: (on_skor, _TersBaslik, sirano, ...)
    zorunlu = []
    gar = [g.lower() for g in (garanti or ()) if g]
    n = 0
    kisa_liste.son_taranan = 0
    for bas, diz in fasta_akisi(yol):
        n += 1
        if ilerle and n % 20000 == 0:
            ilerle(n)
        L = len(diz)
        if L < 100:
            continue
        # THE FILTER (optional, off by default, so the old behaviour is preserved exactly).
        # 2026-08-11: nearly all 500 slots of the short list were filling with unnamed
        # environmental clones, and the "nearest named species" column was coming out
        # empty. With the filter on, the list fills only with NAMED records; the running
        # time is unchanged because the scan already reads the whole file. What changes is
        # which records go onto the heap.
        if suzgec is not None and not suzgec(bas):
            continue
        N += 1
        # THE K-6 CORRECTION (2026-08-03): the counter used to stop at 3, so the score
        # saturated for everything. The real seed set is used now.
        tut = frozenset(i for i in range(nt) if th_l[i] in diz)
        if not tut:
            continue
        for i in tut:
            df[i] += 1             # INVERSE FREQUENCY: counted free in the same pass
        norm = 1.0 - BM25_B + BM25_B * L / ORT
        if gar and any(g in bas.lower() for g in gar):
            zorunlu.append((tut, bas, diz, norm))
            continue
        on = len(tut) / norm       # PRE-FILTER: the raw count normalised for length
        if havuz is None or len(yigin) < havuz:
            heapq.heappush(yigin, (on, _TersBaslik(bas), n, bas, diz, tut, norm))
        elif on > yigin[0][0]:
            heapq.heapreplace(yigin, (on, _TersBaslik(bas), n, bas, diz, tut, norm))
    kisa_liste.son_taranan = n      # HOW MANY RECORDS WERE READ - the output prints this

    # --- STAGE 2: real df -> idf -> re-ranking ---
    idf = [math.log(max(N, 2) / (1.0 + d)) for d in df]

    def _skor(tut, norm):
        return sum(idf[i] for i in tut) / norm

    aday = [(_skor(t, nr), b, d, len(t)) for (_o, _tb, _n, b, d, t, nr) in yigin]
    # THE RANKING CRITERION IS WRITTEN OUT EXPLICITLY: score DESCENDING, header
    # ASCENDING on a tie. Without a fixed tie-break, file order or the internal shape
    # of the heap would decide between equally scored records. This is REQUIRED so that
    # the bulk scanner (all_bin_identities.py) can reproduce the same ranking.
    aday.sort(key=lambda x: (-x[0], x[1]))
    kesme = len(aday) if not ust else ust
    kl = [dict(tohum=int(a[3]), skor=round(a[0], 4), baslik=a[1], dizi=a[2],
               sira=i, kaynak='tohum')
          for i, a in enumerate(aday[:kesme], 1)]

    # The guaranteed records: WHERE would they have landed under normal ranking?
    #
    # THE 2026-08-06 CORRECTION - caught on a clean run. Because records matching the
    # guarantee string were put in a separate bag INDEPENDENTLY of their score, even
    # the ones that would have ENTERED the list on their own merit were stamped
    # 'garanti'. The result: even after the new ranking criterion had solved the
    # problem, the report printed "THE WINNER came in through the patch, the decision
    # DEPENDS on the patch" for 9 of 12 claims, while the same row said "it would have
    # been 1st under seed ranking". The warning refuted itself and hid the real cases
    # of patch dependence.
    # The criterion is now this: if the virtual rank falls INSIDE the cut, the record
    # was going to enter the list anyway, the source is written as 'tohum' and NO
    # warning is printed. Only records coming from OUTSIDE the cut count as 'garanti',
    # and that is where the patch is genuinely needed.
    for tut, b, d, nr in zorunlu[:GARANTI_UST]:
        s = _skor(tut, nr)
        sanal = 1 + sum(1 for a in aday if a[0] > s)
        kendi_gucuyle = bool(kesme) and sanal <= kesme
        kl.append(dict(tohum=int(len(tut)), skor=round(s, 4), baslik=b, dizi=d,
                       sira=sanal,
                       kaynak='tohum' if kendi_gucuyle else 'garanti'))
    return kl


# -------------------------------------------------------------------------
# THE DISCRIMINATING WINDOW - WHY THE CONSERVED REGIONS ARE LEFT OUT
#
# Conserved regions such as 18S, 5.8S and the LSU core are nearly identical across
# every record. Most of an identity measured over the full overlap comes from
# there, and even two completely unrelated organisms give a high percentage. In
# other words a FALSELY HIGH IDENTITY is born exactly in the conserved regions.
#
# This function looks for the window where the best N reference records DIFFER
# FROM ONE ANOTHER: the query's 120 base windows are scored by the DISTRIBUTION of
# reference identities in that window, and the one with the widest spread is
# chosen. In a conserved region the spread is near zero, so it is ruled out by
# itself.
#
# The result puts two numbers side by side: the full overlap identity and the
# identity in the discriminating window. If a claim rests on a high percentage
# coming out of a conserved region, THE GAP BETWEEN THEM makes that visible.
# -------------------------------------------------------------------------
def ayirt_edici_pencere(kayitlar, q):
    """Find the region where the best N reference records DIFFER from one another.

        Conserved regions (18S, 5.8S, the LSU core) are the same in every record and
        are therefore USELESS for discrimination. If a claim rests on a high identity
        coming from there, the full overlap comes out high and the discriminating
        window comes out LOW. That gap is what makes it visible.

        The approach is simple and robust: the query's 120 base windows are scored by
        the DISTRIBUTION of identities among the references themselves, and the most
        divergent window is chosen.

    """
    kayitlar = [k for k in kayitlar if k.get('dizi')][:AYIRT_EDICI_UST]
    if len(kayitlar) < 2:
        return None
    P = 120
    en = None
    for b in range(0, max(1, len(q) - P + 1), 60):
        pen = q[b:b + P]
        if len(pen) < 60 or pen.count('N') > 10:
            continue
        deg = [hizala(pen, k['dizi'])[0] for k in kayitlar]
        yay = max(deg) - min(deg)
        if en is None or yay > en[0]:
            en = (yay, b, pen, deg)
    return en


def kp_yolu(kontrol, etiket, kutu_diz, garanti, kl_ust):
    """The checkpoint file path. A SINGLE SOURCE: both stage I (vtb_tarama) and stage G
        (all_bin_identities.py) use it, so THE SAME BIN IS NEVER SCANNED TWICE; the
        cache is shared between the two stages.

    """
    import hashlib
    imza = hashlib.md5(kutu_diz.encode('utf-8')).hexdigest()[:10]
    g_im = hashlib.md5(('|'.join(sorted(garanti or ()))).encode('utf-8')).hexdigest()[:6]
    return os.path.join(kontrol, 'vtb_%s_%s_%s_%s_kl%d.json'
                        % (re.sub(r'\W+', '_', etiket), imza, g_im, AYAR_IMZASI, kl_ust))


def kl_degerlendir(kl, kutu_diz, kl_ust, taranan=None, t0=None):
    """SHORT LIST -> full alignment -> hits plus self calibration. A SINGLE SOURCE.

        vtb_tarama (one query, one pass per database) and the bulk scanner (every
        query in one pass) call THIS SAME function. The decision logic lives in one
        place; the two routes MUST give the same verdict on the same input.

    """
    t0 = t0 if t0 is not None else time.time()
    t_hiz = time.time()
    isabet = []
    for c in kl:                       # KISA LISTENIN TAMAMI hizalanir
        diz = c['dizi']
        k, d = hizala(kutu_diz if len(kutu_diz) <= len(diz) else diz,
                      diz if len(kutu_diz) <= len(diz) else kutu_diz)
        isabet.append(dict(baslik=c['baslik'][:160], kimlik=k, tohum=c['tohum'],
                           sira=c['sira'], kaynak=c['kaynak'], dizi=diz,
                           hiz_uzunluk=min(len(kutu_diz), len(diz))))
    hiz_sure = round(time.time() - t_hiz, 1)
    # On equal identity the LONGER alignment comes first: 100 per cent over 2900
    # bases is stronger evidence than 100 per cent over 484. Apart from ties the
    # ordering is still by identity.
    isabet.sort(key=lambda x: (-x['kimlik'], -(x.get('hiz_uzunluk') or 0), x['baslik']))

    # --- OZ KALIBRASYON: kazananin tohum sirasi ---
    kazanan_sira = isabet[0]['sira'] if isabet else None
    kazanan_kaynak = isabet[0]['kaynak'] if isabet else None
    uyari = None
    if kazanan_sira is not None and kazanan_sira > SIRA_UYARI_ESIGI:
        uyari = ('THE WINNER CAME FROM RANK %d (the threshold is %d and the '
                 'list holds %d). The cut-off may have started to BIND; raise '
                 '--shortlist and run it again.' % (kazanan_sira, SIRA_UYARI_ESIGI, kl_ust))
    if kazanan_kaynak == 'garanti':
        uyari = ((uyari + u'  ||  ') if uyari else u'') + (
            'THE WINNER entered the short list NOT THROUGH SEEDING but '
            'through the expected taxon patch; on the seed ranking it would '
            'have been number %s. The verdict for this claim DEPENDS ON that '
            'patch, so on a taxon we did not expect the same hit would have '
            'been missed.'
            % kazanan_sira)

    ap = ayirt_edici_pencere(isabet[:AYIRT_EDICI_UST], kutu_diz)
    for i in isabet:                   # diziler JSON'a yazilmaz (devasa olurdu)
        i.pop('dizi', None)
    return dict(durum='TAMAM', isabet=isabet[:10], kayit=len(kl),
                kisa_liste_boyu=kl_ust, hizalanan=len(kl),
                kazanan_sira=kazanan_sira, kazanan_kaynak=kazanan_kaynak,
                sira_uyarisi=uyari, taranan_kayit=taranan,
                sure=round(time.time() - t0, 1), hizalama_suresi=hiz_sure,
                ayirt_edici=(dict(yayilim=round(ap[0], 2), baslangic=ap[1],
                                  kimlikler=[round(x, 2) for x in ap[3]]) if ap else None))


def vtb_tarama(kok, kutu_diz, etiket, dosya, yaz, kontrol, garanti=(), kl_ust=KISA_LISTE,
               suzgec=None, kp_ek=''):
    """In one database: short list -> a full alignment of ALL OF IT -> the best hits.

        EVERY candidate in the short list is aligned (500 candidates take ~10 s). The
        decision is left entirely to alignment; seed ranking is only the criterion
        for COLLECTING CANDIDATES. The winner's seed rank is recorded, and that is
        the evidence for whether the cut off is binding.

    """
    yol = os.path.join(kok, 'REFERENCE_DB', dosya)
    if not os.path.exists(yol):
        return dict(durum='dosya yok')
    # A STABLE key: Python's hash() gives a different value in every PROCESS
    # (PYTHONHASHSEED is random), so the checkpoint never matched. md5 is independent
    # of the process.
    # THE SETTINGS SIGNATURE and THE SHORT LIST SIZE both go into the key, so old
    # vtb_*.json files produced with the 60 item list become invalid and do not come
    # back silently.
    # A FILTERED scan is a different result; if it used the same key the old unfiltered
    # result would come back. The kp_ek key separates them.
    kp = kp_yolu(kontrol, etiket + kp_ek, kutu_diz, garanti, kl_ust)
    if os.path.exists(kp):
        try:
            return json.load(open(kp, encoding='utf-8'))
        # A corrupt checkpoint IS NOT skipped SILENTLY. If it were, a scan that
        # takes hours would run again and nobody would know why; the search for
        # "why is the cache not working" would start in the wrong place.
        except (ValueError, IOError, OSError) as e:
            print(u'     [note] the checkpoint could not be read (%s: %s), '
                  u'this database will be scanned AGAIN'
                  % (type(e).__name__, os.path.basename(kp)))
    t0 = time.time()

    def ilerle(n):
        print(u'     ... %s: %d records scanned (%s)      ' % (etiket, n, sure_metni(time.time() - t0)),
              end='\r', flush=True)
    kl = kisa_liste(yol, kutu_diz, ust=kl_ust, ilerle=ilerle, garanti=garanti,
                    suzgec=suzgec)
    res = kl_degerlendir(kl, kutu_diz, kl_ust,
                         taranan=getattr(kisa_liste, 'son_taranan', None), t0=t0)
    isabet = res['isabet']
    kazanan_sira, kazanan_kaynak = res['kazanan_sira'], res['kazanan_kaynak']
    uyari, hiz_sure = res['sira_uyarisi'], res['hizalama_suresi']
    json.dump(res, open(kp, 'w', encoding='utf-8'), ensure_ascii=False, default=str)
    yaz(u'     %s: %s records SCANNED (all of them), %d short-listed and ALL aligned (%s), best %s%% | WINNER RANK: %s/%d [%s] (%s)'
        % (etiket, '{:,}'.format(res.get('taranan_kayit') or 0).replace(',', ' '),
           len(kl), sure_metni(hiz_sure),
           vir(isabet[0]['kimlik']) if isabet else '-',
           kazanan_sira if kazanan_sira is not None else '-', kl_ust,
           kazanan_kaynak or '-', sure_metni(time.time() - t0)))
    if uyari:
        yaz(u'     >>> WARNING: %s' % uyari)
    return res


# The genus name is taken FROM THE REFERENCE HEADER. The taxonomy tree is
# deliberately NOT USED: Kraken2 already does k-mer plus LCA on the tree, and this
# route exists precisely so as not to repeat it.
def cins_cek(baslik):
    'Pull the genus name out of a reference header. NO taxonomy tree is used.'
    b = baslik
    m = re.search(r'[;|]\s*([A-Z][a-z]+)[ _]([a-z]+)', b)
    if m:
        return m.group(1), '%s %s' % (m.group(1), m.group(2))
    m = re.search(r'\b([A-Z][a-z]{3,})\s+([a-z]{3,})\b', b)
    if m:
        return m.group(1), '%s %s' % (m.group(1), m.group(2))
    m = re.search(r'g__([A-Za-z_]+)', b)
    if m:
        return m.group(1), m.group(1)
    return None, None


# --------------------------------------------------------------- envanter
def envanter_yaz(kok, CIKTI, yaz):
    """Counts EVERY set under REFERENCE_DB and writes down where each one is used.
        For every set that is not used, A REASON IS REQUIRED.

    """
    import glob
    yol = os.path.join(CIKTI, 'VERITABANI_ENVANTERI.md')
    bilinen = {d: (e, kullan, n) for e, d, _t, kullan, n in VTB}
    diskte = sorted(os.path.basename(x) for x in
                    glob.glob(os.path.join(kok, 'REFERENCE_DB', '*.fasta')) +
                    glob.glob(os.path.join(kok, 'REFERENCE_DB', '*.fna')))
    satir = []
    for d in diskte:
        tam = os.path.join(kok, 'REFERENCE_DB', d)
        try:
            boyut = os.path.getsize(tam)
        except OSError:
            boyut = 0
        e, kullan, n = bilinen.get(d, (None, None, None))
        satir.append(dict(dosya=d, mb=round(boyut / 1048576.0), etiket=e or '(tanimsiz)',
                          kimlik=('EVET' if kullan else 'hayir') if e else u'NO, it is not in the list',
                          sebep=n or ('This file IS NOT DEFINED in the '
                                      'database list. To use it, add it to '
                                      'the database list inside '
                                      'identity_verification.py.')))
    with open(yol, 'w', encoding='utf-8') as fh:
        fh.write(u'# REFERENCE_DB inventory: which set is used where\n\n')
        fh.write(u'Generated: %s (rebuilt on every run; do not edit by hand)\n\n'
                 % time.strftime('%Y-%m-%d %H:%M'))
        fh.write(u'| file | MB | label | used at the IDENTITY stage | reason / note |\n')
        fh.write(u'|---|---|---|---|---|\n')
        for r in satir:
            fh.write(u'| `%s` | %d | %s | **%s** | %s |\n'
                     % (r['dosya'], r['mb'], r['etiket'], r['kimlik'], r['sebep']))
        fh.write(u'\n**Rule:** at the identity stage, "there is a cleaner set" is NOT A VALID REASON. Deduplicated sets')
        fh.write(u'**NCBI nt** is not a local file; it is queried as a separate layer, over the network or by manual BLAST. See the report.\n')
    yaz(u'  inventory written: %s' % yol)
    return satir


# --------------------------------------------------------------- NCBI nt katmani
NT_URL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'


def nt_katmani(kutu, dizi, CIKTI, yaz, kip='oto', bekleme=25, tur_ust=40):
    """BLAST against NCBI nt over the URL API (CMD=Put/Get). blastn -remote IS NOT USED.

        On failure it DOES NOT SKIP SILENTLY: it produces a ready input for a manual
        query and the claim is marked "the nt layer did not complete".

    """
    import urllib.request, urllib.parse
    ham = os.path.join(CIKTI, 'nt_ham')
    os.makedirs(ham, exist_ok=True)
    q = dizi[:2500]
    if kip != 'oto':
        return elle_nt_girdi(kutu, q, CIKTI, yaz)
    try:
        # 2026-08-26 (the audit report, finding 5): MEGABLAST WAS TURNED OFF.
        # MEGABLAST is tuned for very similar sequences (a long word size) and
        # misses slightly more distant relatives; and because the result list is
        # cut by bitscore, records that match SHORT but EXACTLY never entered the
        # list at all. MEASURED: in the nt answer for the Petriella bin the number
        # of lines carrying "Petriella" was ZERO, and all of the first 12 hits were
        # long operon records at 92 per cent. Ordinary blastn (word size 11) and a
        # wider list ease both of those. The price: the query takes longer at NCBI.
        p = dict(CMD='Put', PROGRAM='blastn', DATABASE='nt',
                 WORD_SIZE='11', QUERY=q, HITLIST_SIZE='100',
                 FORMAT_TYPE='Text')
        with urllib.request.urlopen(NT_URL, urllib.parse.urlencode(p).encode(),
                                    timeout=90) as f:
            s = f.read().decode('utf-8', 'replace')
        m = re.search(r'RID = (\S+)', s)
        if not m:
            yaz(u'     NCBI nt: could not obtain a RID, falling back to the manual route')
            return elle_nt_girdi(kutu, q, CIKTI, yaz)
        rid = m.group(1)
        yaz(u'     NCBI nt: job submitted (RID %s), waiting...' % rid)
        son = ''
        for i in range(tur_ust):
            time.sleep(bekleme)
            u2 = NT_URL + '?' + urllib.parse.urlencode(
                dict(CMD='Get', RID=rid, FORMAT_TYPE='Text'))
            with urllib.request.urlopen(u2, timeout=90) as f:
                son = f.read().decode('utf-8', 'replace')
            if 'Status=WAITING' not in son:
                break
            print('     ... nt %d. yoklama          ' % (i + 1), end='\r', flush=True)
        open(os.path.join(ham, '%s.txt' % re.sub(r'\W+', '_', kutu)), 'w',
             encoding='utf-8').write(son)
        isabet = []
        for mm in re.finditer(r'^>\s*(\S.*)$', son, re.M):
            isabet.append(dict(baslik=mm.group(1).strip()[:160], kimlik=None))
        if not isabet:
            for satir in son.splitlines():
                mm = re.match(r'^(\S.{20,70}?)\s{2,}[0-9.]+\s+[0-9.]+\s', satir)
                if mm:
                    isabet.append(dict(baslik=mm.group(1).strip(), kimlik=None))
        # --- THE 2026-08-26 FIX (the audit report, findings 3 and 4) ---
        # The old code made two mistakes:
        #   (a) It took the whole number BLAST had ROUNDED ("(99%)"). That number
        #       was then compared against thresholds with two decimals such as
        #       98.70 per cent, so whether it passed the threshold came down to
        #       rounding. MEASURED: all 11 of the 11 nt identities in the output
        #       were whole numbers. Yet the numerator and denominator ("1234/1250")
        #       sit on the same line and give the exact identity.
        #   (b) It assigned an identity only to the FIRST hit (break). The identity
        #       of the other nine stayed None, savunulabilir_duzey filtered them
        #       out, and the "second hit" check NEVER ran in the nt layer.
        # Now each alignment block's own Identities line is matched to that block's
        # header. Blocks start with '>' and Identities sits inside the block.
        bloklar = re.split(r'(?m)^>', son)
        blok_kimlik = []
        for b in bloklar[1:]:
            m2 = re.search(r'Identities\s*=\s*(\d+)/(\d+)', b)
            if m2 and int(m2.group(2)):
                blok_kimlik.append((100.0 * int(m2.group(1)) / int(m2.group(2)),
                                    int(m2.group(2))))
            else:
                blok_kimlik.append((None, None))
        for i, (k, uz) in enumerate(blok_kimlik):
            if i < len(isabet):
                isabet[i]['kimlik'] = k
                isabet[i]['hiz_uzunluk'] = uz
        if isabet and isabet[0].get('kimlik') is None:
            mm = re.search(r'Identities\s*=\s*(\d+)/(\d+)', son)
            if mm and int(mm.group(2)):
                isabet[0]['kimlik'] = 100.0 * int(mm.group(1)) / int(mm.group(2))
        if not isabet:
            yaz(u'     NCBI nt: could not parse the response, falling back to the manual route')
            g = elle_nt_girdi(kutu, q, CIKTI, yaz)
            g['not_'] = 'the automatic reply could not be parsed; the raw reply is under nt_ham/'
            return g
        # --- THE 2026-08-26 EVENING FIX (the other half of finding 5) ---
        # The BLAST text report is in BITSCORE order, NOT in identity order.
        # The code, though, took isabet[0] for "the best". Once the search moved
        # from MEGABLAST to blastn, that silently produced a wrong vote: in bin
        # F1-4_2093780 places one to four went to Lomentospora at 95.21 per cent
        # because its alignment was longer, and Petriella at 96.49 per cent fell
        # to fifth and did not count as "the best hit". We had already corrected
        # the same length bias for the local databases; it had been left in the
        # nt layer. The ordering goes by identity, and on a tie by the longer
        # alignment.
        isabet.sort(key=lambda h: (-(h.get('kimlik') if h.get('kimlik') is not None else -1),
                                   -(h.get('hiz_uzunluk') or 0)))
        yaz(u'     NCBI nt: %d hits (ordered by identity, best %s per cent)'
            % (len(isabet),
               ('%.4f' % isabet[0]['kimlik']) if isabet[0].get('kimlik') is not None else '?'))
        # 25 are kept rather than 10: 100 hits are asked for, 50 came back, and
        # trimming to the first 10 in bitscore order threw away records with a
        # HIGH identity but a short alignment altogether.
        return dict(durum='TAMAM', isabet=isabet[:25], kaynak='NCBI nt (URL API)')
    except Exception as e:
        yaz(u'     NCBI nt FAILED (%s), falling back to the manual route' % type(e).__name__)
        g = elle_nt_girdi(kutu, q, CIKTI, yaz)
        g['not_'] = u'%s: %s' % (type(e).__name__, e)
        return g


def elle_nt_girdi(kutu, q, CIKTI, yaz):
    """If the network cannot be reached: a ready query file for a manual BLAST, plus a
        result template.

    """
    d = os.path.join(CIKTI, 'nt_elle')
    os.makedirs(d, exist_ok=True)
    fa = os.path.join(d, '%s.fasta' % re.sub(r'\W+', '_', kutu))
    with open(fa, 'w', encoding='utf-8') as fh:
        fh.write(u'>%s (identity verification query, first %d bases)\n' % (kutu, len(q)))
        for i in range(0, len(q), 70):
            fh.write(q[i:i + 70] + '\n')
    sab = os.path.join(d, 'NT_SONUC_SABLONU.tsv')
    if not os.path.exists(sab):
        with open(sab, 'w', encoding='utf-8', newline='') as fh:
            fh.write(u'# Write the NCBI nt results HERE, then:\n')
            fh.write(u'#   python3 verification/identity_verification.py --root . --nt-load IDENTITY_RESULT/nt_elle/NT_SONUC_SABLONU.tsv\n')
            fh.write(u'# Address: https://blast.ncbi.nlm.nih.gov/Blast.cgi (Nucleotide BLAST, database = nt)\n')
            fh.write(u'# The query files are in the same directory: <bin>.fasta\n')
            w = csv.writer(fh, delimiter='\t')
            w.writerow(['kutu', 'en_iyi_isabet_basligi', 'kimlik_yuzde', 'notunuz'])
    yaz(u'     NCBI nt: manual query file written -> %s' % os.path.basename(fa))
    return dict(durum='ELLE GEREKIR', isabet=[], kaynak='NCBI nt, by hand',
                girdi=fa, sablon=sab)


# "The best hit" is asked for in one place only. The list order used to be
# trusted; that silently produced a wrong vote once a source whose order is
# bitscore (NCBI nt) was added. Wherever the order comes from, the pick is now
# by identity, and it is done on the consuming side as well so that old cache
# files are read correctly too.
def en_iyi_isabet(v):
    isabetler = [h for h in (v.get('isabet') or []) if h.get('kimlik') is not None]
    if not isabetler:
        liste = v.get('isabet') or []
        return liste[0] if liste else None
    return max(isabetler, key=lambda h: (h['kimlik'], h.get('hiz_uzunluk') or 0))


# Reads the hand filled NT_SONUC_SABLONU.tsv. A row left empty counts as "not
# done"; zero and empty are not the same thing.
def nt_yukle(yol):
    out = {}
    if not yol or not os.path.exists(yol):
        return out
    with open(yol, encoding='utf-8') as fh:
        for r in csv.DictReader((s for s in fh if not s.startswith('#')), delimiter='\t'):
            b = (r.get('en_iyi_isabet_basligi') or '').strip()
            if not b:
                continue
            try:
                k = float((r.get('kimlik_yuzde') or '').replace(',', '.'))
            except ValueError:
                k = None
            out[r['kutu'].strip()] = dict(durum='TAMAM (elle)', kaynak='NCBI nt, by hand',
                                          isabet=[dict(baslik=b, kimlik=k)])
    return out


# --------------------------------------------------------------- NAMING
# HAVING A NAME and CLAIMING AN IDENTITY are two different things. This section
# reports both without confusing them, so that even for a bin we cannot identify we
# can still say "the nearest record is this one, at this percentage".
#
# The species and genus thresholds vary by locus, and using a single number would
# be misleading.
#
# THE 2026-08-26 UPDATE: the values were checked against the literature and the
# FUNGAL thresholds were CORRECTED.
#
#   PROKARYOTIC 16S (SSU) - UNCHANGED, confirmed
#     Kim et al. 2014 (IJSEM), over a million comparisons: the species boundary
#     is 98.65 per cent (which overlaps ANI 95 to 96 per cent). Our 98.70 is the
#     same value in practice (Stackebrandt and Ebers 2006). For the genus figure
#     of 94.5 per cent, Yarza et al. 2014.
#     A WARNING: Hackmann 2025 (IJSEM), 191 million alignments over 19,556 type
#     strains, puts the same species range at 97.2 to 100 and the genus range at
#     90.1 to 99.0, and stresses that THE TWO RANGES OVERLAP. So a single sharp
#     cut point is an approximation in principle; the AYRIM_PAYI rule exists for
#     exactly that reason. Rossi-Tamisier et al. 2015: 57 per cent of valid
#     species sit closer to one another than 98.7 per cent, and only 10.8 per
#     cent of 158 genera obey both thresholds.
#
#   FUNGAL ITS AND LSU - CHANGED (the earlier values were too loose)
#     Vu et al. 2018 (Studies in Mycology), 12,000 ex-type strains of 7,300
#     species:
#       ITS species 99.6 / genus 94.3 ; LSU species 99.8 / genus 98.2
#     Here ITS was on 98.5 / 90.0 and the fungal LSU was using the prokaryotic
#     value of 98.7. MEASURED: with the ITS threshold at 99.6 per cent, 4 of the
#     12 bins that were given a species name in F1 and F2 drop to genus level
#     (the ones between 98.55 and 99.17 per cent). The same study shows that 17
#     to 18 per cent of species cannot be separated AT ALL by ITS or LSU; in
#     fungi "we could not give a species" is usually the limit of the marker.
#     Vu et al. 2022 (dnabarcoder) proposes a local threshold per clade instead
#     of one global threshold; that is the next improvement to make.
#
#   NOTE: a SEPARATE key (LSU_MANTAR) is used for the fungal LSU; the
#   prokaryotic 23S and the fungal 28S cannot share one threshold.
#
#   A KNOWN LIMIT (2026-08-26): the locus is a property of THE DATABASE, not of
#   THE QUERY. SILVA LSU NR99 and Parc carry both prokaryotic and eukaryotic
#   records, so when a fungal query is measured against those sets the
#   prokaryotic threshold (98.7 per cent) is applied where 99.8 per cent would be
#   needed. The sets that are fungal throughout (UNITE ITS, RefSeq fungal ITS and
#   28S) were tied to the right threshold. The correct fix for SILVA LSU is to
#   choose the threshold from the query's domain; that IS NOT DONE at the moment,
#   which means the threshold stays LOOSE for fungal claims. It changes nothing
#   in the present 12 claims (they all win at 100.00 per cent), but it could
#   matter for a new one.
TUR_ESIGI = {'SSU': 98.7, 'LSU': 98.7, 'LSU_MANTAR': 99.8,
             'ITS': 99.6, 'OPERON': 98.7, 'KARISIK': 98.7}
CINS_ESIGI = {'SSU': 94.5, 'LSU': 94.5, 'LSU_MANTAR': 98.2,
              'ITS': 94.3, 'OPERON': 94.5, 'KARISIK': 94.5}
AYRIM_PAYI = 0.5     # if the gap between the best and the second is under this, use "cf."

# ------------------------------------------------- THE MINIMUM ALIGNMENT
# 2026-08-26. The identity percentage ON ITS OWN is not enough: 100 per cent over
# 484 bases and 100 per cent over 2900 bases are not the same evidence. All
# through this module the hits were ordered by identity ALONE, and a short record
# could beat a long one.
#
# MEASURED: claims 1, 2, 3 and 9 take the name "Petriella setifera" from record
# AY882347.1.484, at 100.00 per cent but over ONLY 484 bases, while the query
# consensus is 3,706 bases. The 13 per cent stretch the reference covers matches
# exactly; that is evidence of a GENUS, not of a species.
#
# The answer IS NOT TO THROW THE SHORT HIT AWAY: Petriella's longest record in
# SILVA LSU Parc is 484 bases in the first place, and throwing it away would lose
# the taxon entirely. The answer is to give NO SPECIES NAME from a short
# alignment, to stop at genus level, and to write down the reason. The evidence is
# kept and the claim is not inflated.
#
# The thresholds come from the usual barcode lengths: species level identity in
# 16S generally wants a near full length sequence (>=1200 bp); in the LSU the
# measure is the D1/D2 region (about 600 bp); and because the ITS barcode is about
# 500 to 700 bp, a floor of 400 bp was taken at first.
# THE 2026-08-26 MEASUREMENT: a length sweep trimmed every bin's consensus from
# the centre down to L = 200 to 1400 bases, aligned it again and compared it with
# the name obtained at full length. The result:
#   SSU : from L>=800 upwards the species agreement is 100 per cent and the number
#         of WRONG species is ZERO (a plateau). At shorter lengths the agreement
#         falls as far as 66.7 per cent.
#   ITS : at L=500 three bins were given a WRONG species name; the wrong names
#         only reach zero at L>=1000 and the agreement settles at a plateau of 80
#         per cent.
# That measurement tests THE QUERY length, not the alignment length; the two are
# not the same thing. It still says one thing directly: producing a species name
# from short ITS fragments carries a risk of a wrong name, and that risk WAS
# MEASURED.
# So the ITS floor went from 400 to 600 (the length of a full ITS barcode). It was
# NOT set to 1000: that value comes from the query length, and carrying it over to
# the alignment floor would also remove legitimate full ITS references.
# The SSU floor was LEFT at 1200: the measurement shows 800 is enough, but 1200 is
# the value in the literature (a near full length 16S) and sits well beyond the
# plateau, on the safe side. In this data it costs one bin.
EN_AZ_HIZALAMA = {'SSU': 1200, 'LSU': 600, 'LSU_MANTAR': 600,
                  'ITS': 600, 'OPERON': 1200, 'KARISIK': 1200}


# -------------------------------------------------------------------------
# UNNAMED ENVIRONMENTAL RECORDS - THE 2026-08-21 BUG FIX
#
# NCBI nt and sets like it are full of "unnamed" records:
#     KJ734864.1 Uncultured prokaryote clone D5 16S ribosomal RNA gene
#     KJ957653.1 Uncultured bacterium clone 4B-11 16S ribosomal RNA gene
#     GQ503828.1 Bacterium enrichment culture clone R4-53B 16S ribosomal RNA
#
# ad_coz's second regular expression ( \b[A-Z][a-z]{3,}\s+[a-z]{3,}\b ) was taking
# these for binomials. MEASURED (before this fix, from the real output):
#     'Uncultured prokaryote clone D5...'  -> genus='Uncultured', species='Uncultured prokaryote'
#     'Bacterium enrichment culture...'    -> genus='Bacterium',  species='Bacterium enrichment'
# Because the identity was 99%, savunulabilir_duzey() counted that as a name at
# SPECIES level and the claim was stamped VERIFIED. In other words THE ABSENCE OF
# AN ANSWER was being reported as a confirmed identity. There were four instances
# of it in IDENTITY_RESULT/kimlik_iddialari.tsv ("Uncultured prokaryote", "Uncultured
# bacterium" twice, and "Bacterium enrichment").
#
# This trap was ALREADY KNOWN in the project and was filtered in two other places:
#     screening/exclusion_coverage_check.py:103
#     screening/order_classes.py:200
# and the docstring of verification/ncbi_layer.py separately documents that unnamed
# clones do not get caught by the exclusion filter. The one module actually
# RESPONSIBLE for naming did not have it; it was added.
#
# IMPORTANT: an unnamed hit IS NOT WORTHLESS. "Your bin matches environmental
# clones at 99%" carries information. It simply cannot be A NAME and cannot found a
# species level claim. So the hit is not discarded; its name counts as unresolved
# and the level becomes 'ADLANDIRILAMIYOR'. The module's 'adsiz' claim type
# (around line 200) already carried this idea.
# -------------------------------------------------------------------------
ADSIZ_JETONLARI = (
    'uncultured', 'unclassified', 'unidentified', 'environmental',
    'metagenome', 'enrichment', 'bacterium', 'prokaryote', 'archaeon',
    'eukaryote', 'organism', 'symbiont', 'candidate', 'clone', 'isolate',
    'synthetic', 'construct',
)


def adsiz_mi(ad):
    """Is this string A TAXON NAME, or the description of an unnamed record?

        'Petrimonas sulfuriphila' -> False   (a real name)
        'Uncultured bacterium'    -> True    (a description, not a name)
        'Bacterium enrichment'    -> True

    """
    if not ad:
        return True
    k = [x for x in re.split(r'[^A-Za-z]+', ad) if x]
    if not k:
        return True
    # If the first word (the one sitting where a genus would) is an unnamed token, this
    # is not a name.
    if k[0].lower() in ADSIZ_JETONLARI:
        return True
    # Like 'Bacterium enrichment': the second word is a description, not a species epithet.
    if len(k) > 1 and k[1].lower() in ADSIZ_JETONLARI:
        return True
    return False


def ad_coz(baslik):
    """Extract (genus, species, full_name) from a reference header. The taxonomy tree
        IS NOT USED.

        For an unnamed environmental record it returns (genus, species) = (None,
        None); full_name still comes back, because it has to be shown as "the
        nearest record".

    """
    b = (baslik or '').strip()
    tam = b[:120]
    m = re.search(r'[;|]\s*([A-Z][a-z]{2,})[ _]([a-z]{2,})', b)
    if not m:
        m = re.search(r'\b([A-Z][a-z]{3,})\s+([a-z]{3,})\b', b)
    if m:
        _c, _t = m.group(1), '%s %s' % (m.group(1), m.group(2))
        if adsiz_mi(_t):
            return None, None, tam
        return _c, _t, tam
    m = re.search(r'g__([A-Za-z_]+)', b)
    if m:
        return m.group(1), None, tam
    m = re.search(r'f__([A-Za-z_]+)', b)
    if m:
        return None, None, tam
    return None, None, tam


# -------------------------------------------------------------------------
# HAVING A NAME AND CLAIMING AN IDENTITY ARE TWO DIFFERENT THINGS.
# This function works out what can be DEFENDED from the best three hits: a species
# name, a "cf.", a genus only, or no name at all. The thresholds vary by locus
# (species separation works differently in ITS than in 16S), and one single number
# would be misleading. If the gap between the best and the second is smaller than
# AYRIM_PAYI, a species assignment cannot be defended and "cf." is used.
# -------------------------------------------------------------------------
def erisim_no(baslik):
    """The GenBank or EMBL accession number out of a reference header.

        SILVA   'AY882347.1.484 Eukaryota;...'        -> AY882347
        UNITE   'AY882347|k__Fungi;...'               -> AY882347
        RefSeq  'NR_119564.1 Petriella setifera ...'  -> NR_119564

    """
    b = (baslik or '').strip()
    if not b or b == '-':
        return None
    m = re.match(r'^([A-Z]{1,2}_?[0-9]{5,8})', b)
    return m.group(1) if m else None


def savunulabilir_duzey(isabetler, lokus='SSU'):
    """A wrapper: the decision itself is made by _duzey_karar, and here THE LENGTH
        OF THE ALIGNMENT THAT PRODUCED THE NAME is added.

        WHY: the en_iyi_hiz_uzunluk column was writing the length of THE BEST hit.
        When the name came from a different (named) hit, the column WAS NOT SHOWING
        the alignment that produced the name; a reader saw 600 bases and took the
        name to have come from there, while it had come from another record of 484
        bases. Now, once the decision comes back, the header it used is looked up in
        the hit list and its length is written into a separate key.

    """
    sonuc = _duzey_karar(isabetler, lokus)
    hedef = (sonuc.get('en_iyi') or u'')[:60]
    if hedef:
        for h in (isabetler or []):
            if (h.get('baslik') or u'').startswith(hedef[:40]):
                sonuc['adi_ureten_uzunluk'] = h.get('hiz_uzunluk')
                sonuc['adi_ureten_kimlik'] = h.get('kimlik')
                break
    return sonuc


def _duzey_karar(isabetler, lokus='SSU'):
    """Work out the DEFENSIBLE taxonomic level from the best three hits.

        Returns: dict(duzey, onerilen_ad, gerekce, en_iyi, ikinci, ucuncu)

        The rule, honest naming:
          * identity >= species threshold AND a clear gap to the second -> a SPECIES name
          * identity >= species threshold BUT the second is another species,
            and close                                                   -> "cf." (species uncertain)
          * below the species threshold, above the genus threshold      -> GENUS level
          * identity >= the species threshold BUT the alignment IS SHORT -> GENUS
            level (a short record is not evidence of a species)
          * below the genus threshold                                   -> FAMILY or above;
            NO name is given, only "the nearest record"


    """
    say = [i for i in (isabetler or []) if isinstance(i.get('kimlik'), (int, float))]
    if not say:
        ilk = (isabetler or [{}])[0]
        _c, _t, tam = ad_coz(ilk.get('baslik', ''))
        return dict(duzey='BELIRLENEMEDI', onerilen_ad='-',
                    gerekce='no hit carries a numeric identity',
                    en_iyi=tam or '-', ikinci='-', ucuncu='-')
    say.sort(key=lambda x: -x['kimlik'])
    ilk = say[0]
    c1, t1, tam1 = ad_coz(ilk['baslik'])
    k1 = ilk['kimlik']
    tam2 = tam3 = '-'; t2 = None; k2 = None
    if len(say) > 1:
        c2, t2, tam2 = ad_coz(say[1]['baslik']); k2 = say[1]['kimlik']
    if len(say) > 2:
        _c3, _t3, tam3 = ad_coz(say[2]['baslik'])
    # 2026-08-26: THE LOCUS THRESHOLD IS CORRECTED BY THE QUERY'S DOMAIN.
    # The locus is a property of the database: SILVA LSU carries both prokaryotic
    # and eukaryotic records and arrives under the label 'LSU', so a fungal query
    # would have the prokaryotic threshold (98.7 per cent) applied to it, while
    # the fungal LSU needs 99.8 per cent. The lineage of the best hit tells us the
    # domain (SILVA headers start with "Eukaryota;..." or "Bacteria;/Archaea;...").
    if lokus == 'LSU' and re.search(r'(Eukaryota|Fungi|Dikarya)',
                                    ilk.get('baslik', '') or ''):
        lokus = 'LSU_MANTAR'
    te = TUR_ESIGI.get(lokus, 98.7); ce = CINS_ESIGI.get(lokus, 94.5)

    # 2026-08-21: IF THE BEST HIT IS UNNAMED, NO NAME CAN BE GIVEN WHATEVER THE
    # IDENTITY PERCENTAGE IS. This is a separate branch, because the "no name can be
    # given" branch below states the reason as "the identity is below the genus
    # threshold", and here the situation is the exact opposite: the identity may be
    # 99% while THE MATCHED RECORD ITSELF is unnamed. A wrong reason is more dangerous
    # than a wrong verdict: the reader concludes "so a better reference would produce a
    # name", when the problem is not the reference's CLOSENESS but its NAMELESSNESS.
    # --- 2026-08-26: AN UNNAMED BEST HIT LEVEL WITH A NAMED ONE ---
    # Once the length tie break was added, this came into view: several records
    # sit at 100.00 per cent and the LONGEST of them is unnamed (JX501314.1.600).
    # The old code refused to name outright. But when a NAMED record sits at the
    # SAME closeness (AY882347 at 100.00 per cent), the existence of the unnamed
    # record DOES NOT ERASE the named one; it only says the taxon is also present
    # in unnamed environmental deposits, which is ordinary. So: when the best hit
    # is unnamed, the FIRST NAMED hit within AYRIM_PAYI is looked up and the
    # naming runs through that one, with the unnamed record carried as a note.
    # When there is no named hit nearby, the old behaviour (no name) stands.
    if not t1 and not c1:
        adli = next((h for h in say[1:]
                     if (k1 - h['kimlik']) <= AYRIM_PAYI
                     and any(ad_coz(h.get('baslik', ''))[:2])), None)
        if adli is not None:
            _c2, _t2, _tam2 = ad_coz(adli['baslik'])
            yedek = dict(adli)
            yedek['_adsiz_rakip'] = tam1
            kalan = [adli] + [h for h in say if h is not adli]
            alt = savunulabilir_duzey(kalan, lokus)
            alt['gerekce'] = (u'THE BEST HIT IS UNNAMED (%s, %s per cent) but a named '
                              u'record sits at the SAME closeness (%s per cent); the '
                              u'naming was carried out through that one. || %s'
                              % (tam1[:70], vir(k1), vir(adli['kimlik']),
                                 alt.get('gerekce', '')))
            return alt
        return dict(duzey='ADLANDIRILAMIYOR (referans adsiz)',
                    onerilen_ad='an unnameable lineage - THE NEAREST RECORD: %s (%s per cent)'
                                % (tam1, vir(k1)),
                    gerekce="the best hit's identity is %s per cent (the "
                            'species threshold is %s per cent) BUT the record '
                            'matched is an unnamed environmental sequence: '
                            '"uncultured", "enrichment", "clone" and the '
                            'like. A high identity here IS NOT A NAME, only '
                            'the information that this sequence overlaps '
                            'environmental clones. NO species or genus claim '
                            'can be built on it.'
                            % (vir(k1), vir(te)),
                    en_iyi=tam1, ikinci=tam2, ucuncu=tam3)
    # --- NO SPECIES NAME COMES FROM A SHORT ALIGNMENT (2026-08-26) ---
    # 100 per cent over 484 bases and 100 per cent over 2900 bases are not the same
    # evidence. The short hit IS NOT THROWN AWAY; it is only kept from rising to
    # species level, and the reason is written down.
    enaz = EN_AZ_HIZALAMA.get(lokus, 600)
    uz1 = ilk.get('hiz_uzunluk') or 0
    if k1 >= te and t1 and uz1 and uz1 < enaz:
        return dict(duzey='CINS (hizalama kisa)',
                    onerilen_ad='%s sp.' % (c1 or (t1 or '').split()[0]),
                    gerekce=u'the identity of %s per cent passed the species threshold '
                            u'of %s per cent BUT the alignment is only %d bases, and a '
                            u'species name at the %s locus wants at least %d. A short '
                            u'record is evidence of a GENUS, not of a SPECIES; the '
                            u'nearest record is still shown.%s'
                            % (vir(k1), vir(te), uz1, lokus, enaz,
                               (u'  ON THE EDGE: %d bases short of the threshold, so this '
                                u'decision rests on one threshold choice. Change '
                                u'EN_AZ_HIZALAMA and test the sensitivity.' % (enaz - uz1))
                               if (enaz - uz1) <= max(20, enaz * 0.05) else u''),
                    en_iyi=tam1, ikinci=tam2, ucuncu=tam3)
    if k1 >= te and t1:
        farkli_tur = bool(t2) and t2 != t1
        yakin = (k2 is not None and (k1 - k2) < AYRIM_PAYI)
        # THE 2026-08-26 FIX (the audit report, finding 2).
        # The old rule: if the second hit's NAME COULD NOT BE RESOLVED it did not
        # count as a rival and the species name was handed out freely. But an
        # UNNAMED second hit at the same identity is exactly the situation in which
        # a species name cannot be defended: it says "there is another record
        # standing just as close to this sequence, and we do not know its name".
        # MEASURED: in 4 claims the reason read "there is a clear gap (0.00 per
        # cent)" - we were reporting a gap of zero as a clear gap.
        if (not t2) and yakin and tam2 and tam2 != '-':
            return dict(duzey='CINS (ikinci isabet adsiz, ayni yakinlikta)',
                        onerilen_ad='%s cf. %s' % (c1, t1.split()[-1]) if c1 else t1,
                        gerekce=u'the identity of %s per cent is ABOVE the species '
                                u'threshold of %s per cent, but the second hit is just as '
                                u'close (%s per cent, %s per cent between them) and ITS '
                                u'NAME CANNOT BE RESOLVED. An unnamed rival is the very '
                                u'case in which a species name cannot be defended, so '
                                u'"cf." was used. The second hit: %s'
                                % (vir(k1), vir(te), vir(k2), vir(k1 - k2), tam2[:80]),
                        en_iyi=tam1, ikinci=tam2, ucuncu=tam3)
        if farkli_tur and yakin:
            return dict(duzey='CINS (tur belirsiz)', onerilen_ad='%s cf. %s' % (c1, t1.split()[-1]),
                        gerekce='the identity of %s per cent is ABOVE the '
                                'species threshold of %s per cent, but the '
                                'second hit is another species (%s at %s per '
                                'cent) and only %s per cent separates them, '
                                'so a species assignment is indefensible and '
                                '"cf." was used'
                                % (vir(k1), vir(te), t2, vir(k2), vir(k1 - k2)),
                        en_iyi=tam1, ikinci=tam2, ucuncu=tam3)
        return dict(duzey='TUR', onerilen_ad=t1,
                    gerekce='the identity of %s per cent is above the species '
                            'threshold of %s per cent and there is a clear '
                            'gap to the second hit (%s)'
                            % (vir(k1), vir(te),
                               ('%%%s' % vir(k1 - k2)) if k2 is not None else 'there is no second hit'),
                    en_iyi=tam1, ikinci=tam2, ucuncu=tam3)
    if k1 >= ce and c1:
        return dict(duzey='CINS', onerilen_ad='%s sp.' % c1,
                    gerekce='the identity of %s per cent is BELOW the species '
                            'threshold of %s per cent and above the genus '
                            'threshold of %s per cent, so no species name can '
                            'be given'
                            % (vir(k1), vir(te), vir(ce)),
                    en_iyi=tam1, ikinci=tam2, ucuncu=tam3)
    return dict(duzey='AILE ve USTU (ad VERILEMEZ)',
                onerilen_ad='an unnameable lineage - THE NEAREST RECORD: %s (%s per cent)'
                            % ((t1 or c1 or tam1), vir(k1)),
                gerekce='the identity of %s per cent is BELOW even the genus '
                        'threshold of %s per cent. This IS NOT AN IDENTITY, '
                        'only the nearest reference record, and it has to go '
                        'into the report as "we could not identify it; this '
                        'is the nearest".'
                        % (vir(k1), vir(ce)),
                en_iyi=tam1, ikinci=tam2, ucuncu=tam3)


# -------------------------------------------------------------------------
# THE SINGLE SOURCE OF NAMING (2026-08-21).
#
# This block used to sit inside the main loop. It was pulled into a function of its
# own because it is now called from TWO places: a fresh run, and RE-DERIVATION from
# a checkpoint. Two copies would drift apart over time, and in this code base two
# versions of the same check contradicting one another has already been measured
# twice (cross_check D9 / orientation_code_scan, and _kayit_coz / taxonomy).
# -------------------------------------------------------------------------
def adlandirmayi_turet(bulgular):
    'Derives the naming from the raw hits. Returns: (naming_dict, locus)'
    lokus_tab = {e: t for e, _d, t, _k, _n in VTB}
    havuz = []
    for et, v in (bulgular or {}).items():
        if not str(v.get('durum', '')).startswith('TAMAM'):
            continue
        for i in (v.get('isabet') or [])[:5]:
            havuz.append(dict(i, _vtb=et, _lokus=lokus_tab.get(et, 'SSU')))
    sayisal = [h for h in havuz if isinstance(h.get('kimlik'), (int, float))]
    sayisal.sort(key=lambda x: -x['kimlik'])
    lokus = sayisal[0]['_lokus'] if sayisal else 'SSU'
    adl = savunulabilir_duzey(sayisal or havuz, lokus)

    # THE FIVE NEAREST ORGANISMS (2026-08-21, three before).
    #
    # The deduplication is by ORGANISM, not by record. The reason: the pool holds hits
    # from 13 databases and the same species appears in most sets. Without
    # deduplication the "nearest five" easily becomes five records of THE SAME
    # organism, and the reader learns nothing. Since the point of the question is
    # "what else is close", the distinction has to be at the ORGANISM.
    #
    # Unnamed records (uncultured/enrichment) produce no name, so they are separated by
    # their full header. For each organism the record with the HIGHEST identity is kept
    # (the list is already sorted by identity, so that is the first one seen).
    gorulen = set()
    sira = 0
    for h_ in sayisal:
        c_, t_, tam_ = ad_coz(h_['baslik'])
        anahtar = (t_ or c_ or tam_ or '').strip().lower()
        if not anahtar or anahtar in gorulen:
            continue
        gorulen.add(anahtar)
        sira += 1
        adl['isabet%d' % sira] = dict(
            tam_ad=tam_, cins=c_ or '-', tur=t_ or '-',
            kimlik=h_.get('kimlik'), uzunluk=h_.get('hiz_uzunluk'),
            vtb=h_['_vtb'])
        if sira >= 5:
            break
    # The length of the alignment THAT PRODUCED THE NAME sits in its own key.
    # en_iyi_hiz_uzunluk is the length of THE BEST hit; the name may have come from
    # a different (named) hit, and then the two numbers part company. Both are
    # written so that a reader knows which one to look at.
    adl['adi_ureten_hizalama'] = adl.get('adi_ureten_uzunluk')
    return adl, lokus


def _en_yakin_etiket(isabet):
    """The label to show in the nearest organism list.

        THREE CASES ARE SHOWN SEPARATELY; confusing two of them misleads:

          1) A species or genus name was resolved  -> the name is written
          2) There IS taxonomy but no binomial     -> the DEEPEST taxonomic token is
             written (for example 'Dysgonomonadaceae'). Such a record is NOT unnamed;
             it simply has not been taken down to species level. Calling it 'unnamed'
             would put a classified SILVA record in the same bucket as an
             environmental clone.
          3) There is NO taxonomy either (an environmental clone) -> the record's
             description, with an 'adsiz:' prefix

        Printing a bare '-' is the worst option: the reader cannot see what was
        matched and quietly takes a 99% row seriously.

    """
    if not isabet:
        return '?'
    for alan in ('tur', 'cins'):
        v = (isabet.get(alan) or '').strip()
        if v and v != '-':
            return v
    tam = (isabet.get('tam_ad') or '').strip()
    if not tam:
        return '?'
    # Does it carry taxonomy? Resolved through screening/taxonomy.py, the ONE place
    # that knows all five header formats (SILVA/UNITE/PR2/ROD/RefSeq).
    try:
        import sys as _s, os as _o
        _kok = _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))
        if _kok not in _s.path:
            _s.path.insert(0, _kok)
        from screening import taxonomy as _TX
        _alan, _jet, _org, _tak_var = _TX.coz(tam, isabet.get('vtb') or '')
        if _tak_var and _jet:
            # En derin jeton = taksonominin en ozgul duzeyi.
            derin = _jet[-1]
            ust = _jet[-2] if len(_jet) > 1 else ''
            return u'%s%s' % (derin, (u' (%s)' % ust) if ust and ust != derin else u'')
    # When the screening package is missing or the taxonomy cannot be resolved, we
    # fall onto the simple "adsiz: <description>" path below. The catch is narrow
    # on purpose: a fault in the code INSIDE this module is no longer hidden, it
    # comes back up.
    except (ImportError, AttributeError, IndexError, KeyError, TypeError, ValueError):
        pass
    tanim = tam.split(' ', 1)[1] if ' ' in tam else tam
    return u'adsiz: %s' % tanim[:52]


def _yeniden_turet(kayit, idd, kons):
    """Re-derive the verdict from the RAW measurement in the checkpoint.

        The scan is not repeated (it would take hours); only the derivation, which
        takes seconds, is redone. A fix to the naming logic therefore takes effect
        without re-running the expensive scan.

        Re-derived: the naming (defensible level, suggested name, the five nearest
        organisms). Kept from the cache: the verdict, the literature check, the
        calibration and vtb_detay, which are outputs of the scan itself.

        Re-derived: the naming, THE VERDICT, the evidence and vtb_detay. Kept from
        the cache: the literature check (a network query) and the calibration (a
        measurement).

        THE 2026-08-26 EVENING FIX - THE RE-DERIVATION DID NOT GO DEEP ENOUGH.
          At first only the naming was re-derived; the verdict, the evidence and
          vtb_detay came from the cache. After the nt ordering fix that produced
          rows CONTRADICTING each other: in one and the same claim the vote said
          'Petriella guttulata' while vtb_detay still said 'Lomentospora
          prolificans', because the detail had been frozen at the moment of the
          scan.
          The principle is unchanged: THE MEASUREMENT is stored (it is expensive),
          THE DERIVATION is redone on every run (it is cheap). A verdict is a
          derivation too - counting votes does not take even a second.
    """
    ham = kayit.get('bulgular') or {}
    try:
        adl, _lokus = adlandirmayi_turet(ham)
        yeni = dict(kayit)
        yeni['adlandirma'] = adl
        # The verdict is counted again from the raw measurements. kons IS GENUINELY
        # needed: claims of the "is this the same organism" kind align the bin
        # consensuses against one another. Had an empty dictionary been passed,
        # those claims would silently have become NOT VERIFIED - the fix would
        # have given birth to a new fault.
        try:
            h, kanit, duzeltme = hukum_ver(idd, ham, kons)
            # --- A MISSING nt LAYER MEANS NOT VERIFIED ---
            # This rule stood apart in the scan path; when the re-derivation
            # skipped it, claim 4 silently went from NOT VERIFIED to VERIFIED.
            # Yet that claim's nt layer had ended in "NEEDS TO BE DONE BY HAND":
            # missing evidence cannot count as verification. The rule was put
            # here as well.
            nt = ham.get(NT_ETIKET)
            if (idd.get('tip') != 'gecici' and nt is not None
                    and not str(nt.get('durum', '')).startswith('TAMAM')):
                h = 'DOGRULANAMADI'
                kanit = (kanit + u'  ||  THE NCBI nt LAYER DID NOT COMPLETE (%s), so the '
                         u'claim was counted as not verified rather than silently '
                         u'skipped.' % (nt.get('durum', 'kosulmadi')))
            yeni['hukum'] = h
            yeni['kanit'] = kanit
            yeni['dogru_ifade'] = duzeltme
            yeni['_hukum_turetme'] = 'again'
        except Exception as e:
            yeni['_hukum_turetme'] = ('FAILED (%s), the verdict in the cache was '
                                      'kept' % type(e).__name__)
        # vtb_detay is A REPORT, not a measurement: it is rebuilt from the raw hits
        # on every run so that it points at the same record the vote does.
        detay = {}
        for et, v in ham.items():
            if str(v.get('durum', '')).startswith('TAMAM') and v.get('isabet'):
                en = en_iyi_isabet(v) or {}
                detay[et] = dict(durum=v['durum'], en_iyi=(en.get('baslik') or '')[:120],
                                 kimlik=en.get('kimlik'),
                                 kazanan_sira=v.get('kazanan_sira'),
                                 kazanan_kaynak=v.get('kazanan_kaynak'),
                                 kisa_liste_boyu=v.get('kisa_liste_boyu'),
                                 taranan_kayit=v.get('taranan_kayit'))
            else:
                detay[et] = dict(durum=v.get('durum', '?'), en_iyi='', kimlik=None)
        yeni['vtb_detay'] = detay
        yeni['_turetme'] = 'again, with the raw measurement from the cache'
        return yeni
    except Exception as e:
        # If the re-derivation fails, the OLD record is not returned unchanged; it is
        # flagged, so the report can say "this row was produced by the old logic".
        yeni = dict(kayit)
        yeni['_turetme'] = 'FAILED (%s), so the old naming was used' % type(e).__name__
        return yeni


# --------------------------------------------------------------- hukum
def _say(liste):
    """Keep the numeric ones. The NCBI nt layer can return identity=None (when the
        percentage cannot be parsed out of the URL API's text output), and comparing
        against None raises TypeError. Everything goes through this filter BEFORE any
        numeric comparison.

    """
    return [x for x in liste if isinstance(x, (int, float))]


# -------------------------------------------------------------------------
# THE AT-LEAST-TWO-INDEPENDENT-DATABASES RULE - WHY IT EXISTS
#
# The best hit of a single database IS NOT ENOUGH for an identity claim. Every set
# carries its own bias: dereplicated sets delete rare genera (measured: SILVA
# LSURef NR99 holds 0 Petriella records while the Parc set of the same release
# holds 82), and another set may carry the same record under an outdated name. A
# verdict resting on one source would report that source's mistake as AN IDENTITY.
#
# So the first check is a count: if fewer than two databases returned a result, no
# claim can come out VERIFIED. Then the votes are counted. If >=2 databases support
# the claim it is VERIFIED; if >=2 databases agree on a DIFFERENT answer it is
# CORRECTION NEEDED and the right wording is written; if no two agree it is
# NOT VERIFIED.
#
# NO CONFIRMATION IS INVENTED. Where the evidence is insufficient it says NOT
# VERIFIED; a gap is never rounded up to a positive answer. Independence is
# enforced too: files in the VTB list that are byte for byte twins, or subsets of
# another set, have been taken out of the vote, or the same record would vote twice.
#
# If the NCBI nt layer did not complete, the verdict is pulled down to NOT VERIFIED
# (see calistir): a missing layer is never skipped silently.
# -------------------------------------------------------------------------
def hukum_ver(idd, bulgular, kons):
    """bulgular: {database_label: result}. Returns: (verdict, evidence, correct_wording)"""
    calisan = {k: v for k, v in bulgular.items()
               if str(v.get('durum', '')).startswith('TAMAM') and v.get('isabet')}
    if len(calisan) < 2:
        return ('DOGRULANAMADI',
                'Only %d database returned a result; this round requires at '
                'least TWO independent databases to agree.' % len(calisan), '')

    # her veritabaninin "oyu": en iyi isabetin cinsi/turu ve kimligi
    oylar = {}
    for et, v in calisan.items():
        en = en_iyi_isabet(v)
        cins, tur = cins_cek(en['baslik'])
        oylar[et] = dict(cins=cins, tur=tur, kimlik=en['kimlik'], baslik=en['baslik'])

    tip = idd['tip']
    tol = idd.get('tolerans', UYUM_TOLERANS)

    if tip in ('kimlik', 'cins_duzeyi'):
        bek_c = idd.get('beklenen_cins'); bek_t = idd.get('beklenen_tur')
        bek_y = idd.get('beklenen_yuzde')
        if bek_c:
            uyan = [e for e, o in oylar.items() if o['cins'] and
                    o['cins'].lower().startswith(bek_c.lower()[:6])]
        elif bek_t:
            uyan = [e for e, o in oylar.items() if o['tur'] and
                    bek_t.split()[-1].lower() in (o['tur'] or '').lower()]
        else:
            uyan = list(oylar)
        kanit = '; '.join('%s: %s %%%s' % (e, oylar[e]['tur'] or oylar[e]['cins'] or '?',
                                           vir(oylar[e]['kimlik'])) for e in oylar)
        # --- 2026-08-26 (the audit report, finding 1): INDEPENDENCE AT RECORD LEVEL ---
        # The rule says "at least two INDEPENDENT databases", but the counter was
        # counting databases. One and the same GenBank record can sit in several
        # sets (MEASURED: AY882347 is in both SILVA LSU Parc and UNITE). Two
        # "databases" are then two copies of the same record and do not count as
        # independent testimony. THE VERDICT IS NOT CHANGED - changing it would
        # break comparison with the older verdicts - but how many DIFFERENT records
        # it rests on is measured and written into the evidence text.
        _kayitlar = set()
        for e in uyan:
            en = erisim_no((oylar[e] or {}).get('baslik') or '')
            _kayitlar.add(en or ('vtb:%s' % e))
        if len(uyan) >= 2 and len(_kayitlar) < 2:
            kanit = (u'[AN INDEPENDENCE WARNING: %d databases agree but all of them rest '
                     u'on the SAME record (%s), which is not two independent '
                     u'testimonies] '
                     % (len(uyan), ', '.join(sorted(_kayitlar)))) + kanit
        elif len(uyan) >= 2:
            kanit = (u'[independence: %d databases, %d different records] '
                     % (len(uyan), len(_kayitlar))) + kanit
        if len(uyan) >= 2:
            if bek_y is not None:
                olculen = _say([oylar[e]['kimlik'] for e in uyan])
                if olculen and min(abs(x - bek_y) for x in olculen) > tol:
                    return ('DUZELTILMELI',
                            kanit,
                            'The taxon is right but the PERCENTAGE differs: '
                            'measured %s against a claim of %s per cent. The '
                            'right wording: "%s, with a measured identity of '
                            '%s per cent".'
                            % (', '.join('%%%s' % vir(x) for x in olculen), vir(bek_y),
                               (bek_t or bek_c), vir(sorted(olculen)[len(olculen) // 2])))
            return ('DOGRULANDI', kanit + '  [%d databases agree]' % len(uyan), '')
        # baska bir cevapta birlesiyorlar mi
        say = {}
        for e, o in oylar.items():
            k = (o['cins'] or o['tur'] or '?')
            say.setdefault(k, []).append(e)
        en_cok = max(say.items(), key=lambda kv: len(kv[1]))
        if len(en_cok[1]) >= 2:
            med = sorted(_say([oylar[e]['kimlik'] for e in en_cok[1]])) or [None]
            return ('DUZELTILMELI', kanit,
                    'The claimed taxon is not supported. %d databases agree '
                    'on "%s" at about %s per cent. The right wording: "the '
                    'bin\'s nearest reference is %s, at %s per cent".'
                    % (len(en_cok[1]), en_cok[0], vir(med[len(med) // 2]),
                       en_cok[0], vir(med[len(med) // 2])))
        return ('DOGRULANAMADI', kanit,
                'The databases give DIFFERENT answers; no two of them agree '
                'on the same taxon.')

    if tip == 'kutu2':
        a, b = idd['kutu'][0], idd['karsi'][0]
        if a not in kons or b not in kons:
            return ('DOGRULANAMADI', 'the bin consensus was not found (%s / %s)' % (a, b), '')
        k, _ = hizala(kons[a] if len(kons[a]) <= len(kons[b]) else kons[b],
                      kons[b] if len(kons[a]) <= len(kons[b]) else kons[a])
        ortak = [e for e, o in oylar.items() if o['cins'] and idd.get('beklenen_cins') and
                 o['cins'].lower().startswith(idd['beklenen_cins'].lower()[:6])]
        kanit = ('A direct alignment of %s against %s: %s per cent. The database votes: %s'
                 % (a, b, vir(k),
                    '; '.join('%s=%s' % (e, oylar[e]['tur'] or oylar[e]['cins']) for e in oylar)))
        bek = idd.get('beklenen_yuzde')
        if abs(k - bek) <= tol and len(ortak) >= 2:
            return ('DOGRULANDI', kanit, '')
        if abs(k - bek) > tol:
            return ('DUZELTILMELI', kanit,
                    'The right wording: "the %s bin is %s per cent identical to %s" (the claim said %s per cent).'
                    % (a, b, vir(k), vir(bek)))
        return ('DOGRULANAMADI', kanit,
                'The identity holds, but two databases did not agree on the same genus.')

    if tip == 'dagilim':
        # supheli kutularin konsensusleri capa kutuya ne kadar benziyor
        capa = idd['kutu'][0]
        satir = []
        for k2 in idd.get('supheli', []):
            if capa in kons and k2 in kons:
                v, _ = hizala(kons[k2] if len(kons[k2]) <= len(kons[capa]) else kons[capa],
                              kons[capa] if len(kons[k2]) <= len(kons[capa]) else kons[k2])
                satir.append((k2, v))
        if not satir:
            return ('DOGRULANAMADI', 'the consensus of the suspect bins could not be read', '')
        kanit = '; '.join('%s=%%%s' % (a, vir(v)) for a, v in satir)
        yuksek = [v for _, v in satir if v >= 99.0]
        if len(yuksek) >= max(2, len(satir) - 1):
            return ('DOGRULANDI',
                    kanit + '  [at consensus level %d of %d bins are 99 per cent or more identical]'
                    % (len(yuksek), len(satir)), '')
        return ('DOGRULANAMADI', kanit,
                'At consensus level %d of %d bins are 99 per cent or more '
                'identical. The read fraction claim of 76 to 86 per cent WAS '
                'NOT MEASURED IN THIS ROUND: a read level measurement is not '
                "part of this script's method and needs a round of its own."
                % (len(yuksek), len(satir)))

    if tip == 'ayrilmaz':
        turler = idd['turler']
        bulunan = []
        for et, v in calisan.items():
            for i in v['isabet'][:6]:
                _, t = cins_cek(i['baslik'])
                for ad, bek in turler:
                    if t and ad.split()[-1].lower() in t.lower():
                        bulunan.append((et, ad, i['kimlik']))
        if not bulunan:
            return ('DOGRULANAMADI', 'not one of the claimed species is among the hits', '')
        kanit = '; '.join('%s: %s %%%s' % (e, a, vir(k)) for e, a, k in bulunan[:8])
        adlar = {a for _, a, _ in bulunan}
        if len(adlar) >= 2:
            kk = _say([k for _, _, k in bulunan])
            if not kk:
                return ('DOGRULANAMADI', kanit,
                        'Both species are hit but NO database returned a '
                        'numeric identity, only a header match, so the claim '
                        'that they do not separate could not be tested with a '
                        'number in this round.')
            if max(kk) - min(kk) <= 0.5:
                return ('DOGRULANDI', kanit + '  [the identity difference between the two species is 0.5 per cent or less]', '')
            return ('DUZELTILMELI', kanit,
                    'Both species are hit but the identities diverge over a '
                    'range of %s, so this measurement does not support the '
                    'wording that they do not separate.'
                    % ('%%%s-%%%s' % (vir(min(kk)), vir(max(kk)))))
        return ('DOGRULANAMADI', kanit, 'only one species was hit')

    if tip == 'adsiz':
        en_iyiler = [(e, o['kimlik'], o['tur'] or o['cins']) for e, o in oylar.items()]
        kanit = '; '.join('%s: %s %%%s' % (e, t, vir(k)) for e, k, t in en_iyiler)
        sayisal = [(e, k, t) for e, k, t in en_iyiler if isinstance(k, (int, float))]
        if not sayisal:
            return ('DOGRULANAMADI', kanit,
                    'No database returned a numeric identity, so the claim '
                    'that no name can be given could not be tested with a '
                    'number in this round.')
        if all(k < 90.0 for _, k, _ in sayisal):
            return ('DOGRULANDI',
                    kanit + '  [no database gives a hit above 90 per cent, so '
                            'the claim that no name can be given is '
                            'supported]', '')
        yuksek = [(e, k, t) for e, k, t in sayisal if k >= 97.0]
        if len(yuksek) >= 2:
            return ('DUZELTILMELI', kanit,
                    'A name looks like it CAN be given: %d databases return a '
                    'hit above 97 per cent (%s). The wording has to be looked '
                    'at again.'
                    % (len(yuksek), ', '.join('%s %s' % (t, vir(k)) for _, k, t in yuksek)))
        return ('DOGRULANAMADI', kanit,
                u'Isabetler %90-97 bandinda; ne "ad verilemez" ne de bir ad kesin.')

    if tip == 'gecici':
        return ('DOGRULANAMADI',
                'This is a statement about method and not the identity of a '
                'single bin, so a reference database sweep cannot confirm it.',
                'This row is not a MEASUREMENT claim but a statement of '
                'UNCERTAINTY, and it can be left as it stands. How to test '
                'it: compare the best hits of every class B bin in this '
                'round; when bins carrying the same taxid go to different '
                'taxa, the statement is right.')

    return ('DOGRULANAMADI', u'bilinmeyen iddia tipi', '')


# --------------------------------------------------------------- driver
# -------------------------------------------------------------------------
# THE DRIVER. It works claim by claim, and for each claim the order is fixed:
#   take the consensus -> short list plus full alignment in every database -> the
#   NCBI nt layer -> hukum_ver -> naming (defensible level plus the best three
#   hits) -> literature check -> self calibration -> write to disk.
#
# THE CHECKPOINTS ARE TWO LEVELS DEEP and both mix AYAR_IMZASI and the short list
# size into the key. Without that: even when the database level cache was
# invalidated, the claim level file would still be found, the claim would be
# skipped entirely, and the old verdict produced with the 60 item list would come
# back silently.
#
# The NCBI nt cache is DELIBERATELY outside that seal (it comes over the network,
# it has nothing to do with the short list size, and it should not be re-queried
# for nothing). Against that, A NETWORK ERROR IS NOT A RESULT and is not cached: a
# single Wi-Fi dropout was poisoning every claim permanently.
#
# SELF CALIBRATION: for every query, the seed rank of the winning hit inside the
# short list is recorded. Whether the cut off is binding therefore NEEDS NO
# SEPARATE MEASUREMENT; the evidence comes out of the run itself.
# -------------------------------------------------------------------------
def calistir(kok, yalniz, sifirla, vtb_ust, nt_kip='oto', nt_yukle_yolu=None,
             lit_kip='oto', kl_ust=KISA_LISTE):
    sys.path.insert(0, kok)
    CIKTI = os.path.join(kok, 'IDENTITY_RESULT')
    KONTROL = os.path.join(CIKTI, 'kontrol')
    os.makedirs(KONTROL, exist_ok=True)
    if sifirla:
        for f in os.listdir(KONTROL):
            try:
                os.remove(os.path.join(KONTROL, f))
            except OSError as e:
                print(u'  could not delete: %s (%s)' % (f, e))
    g = open(os.path.join(CIKTI, 'kosu_gunlugu.txt'), 'a', encoding='utf-8')

    def yaz(s=''):
        print(s, flush=True); g.write(s + '\n'); g.flush()

    yaz('=' * 78)
    yaz(u'  IDENTITY VERIFICATION - reported claims are tested independently')
    yaz(u'  version %s   %s' % (VERSIYON, time.strftime('%Y-%m-%d %H:%M')))
    yaz('=' * 78)
    yaz(u'  Method: SEED + ALIGNMENT (NO taxonomy tree, NO primers, NO k-mer LCA).')
    yaz(u'  Rule  : a claim is only VERIFIED when AT LEAST TWO independent')
    yaz(u'          databases agree. No confirmation is ever invented.')
    yaz('')

    from screening import targets as H
    kons = {d['kutu']: d['dizi'] for d in H.konsensusler()}
    var = [(e, d, t) for e, d, t, kullan, _n in VTB
           if kullan and os.path.exists(os.path.join(kok, 'REFERENCE_DB', d))][:vtb_ust]
    envanter_yaz(kok, CIKTI, yaz)
    yaz(u'  usable databases          : %d  (%s)' % (len(var), ', '.join(e for e, _, _ in var)))
    if len(var) < 2:
        yaz(u'  WARNING: fewer than two databases available, so no claim can come out VERIFIED.')

    nt_onceden = nt_yukle(nt_yukle_yolu)
    if nt_onceden:
        yaz(u'  manually loaded NCBI nt results: %d bins' % len(nt_onceden))
    yaz(u'  NCBI nt layer             : %s' % {'oto': u'automatic (URL API)', 'elle': u'manual (a query file is written)', 'yok': u'SKIPPED (on request)'}[nt_kip])

    iddialar = sorted(IDDIALAR, key=lambda x: (x['oncelik'], x['no']))
    if yalniz:
        iddialar = [i for i in iddialar
                    if (yalniz.isdigit() and int(yalniz) == i['no'])
                    or (not yalniz.isdigit() and yalniz.lower() in i['metin'].lower())]
    yaz(u'  claims to test            : %d  (in priority order)' % len(iddialar))
    kutu_say = len({k for i in iddialar for k in (i.get('kutu') or []) + (i.get('karsi') or [])})
    # --- TIME ESTIMATE: the scan and the ALIGNMENT, separately ---
    # The alignment cost was fitted from measurement (with this script's own aligner):
    #     t ~= 6.7e-6 * short + 4.11e-9 * short * long     [seconds]
    # Measured / model: 600x600 0.0050/0.0055 | 1500x1500 0.0187/0.0193 |
    #                   1500x3500 0.0310/0.0316 | 4000x4000 0.0917/0.0926
    _kutular = {k for i in iddialar for k in (i.get('kutu') or []) + (i.get('karsi') or [])}
    _uz = [min(len(kons[k]), 4000) for k in _kutular if k in kons] or [1500]
    _oq = sum(_uz) / float(len(_uz))
    _ref = 2000.0                      # a typical reference record length (SSU and LSU mixed)
    _bir = 6.7e-6 * min(_oq, _ref) + 4.11e-9 * min(_oq, _ref) * max(_oq, _ref)
    _hiz_cift = _bir * kl_ust          # the alignment time of one bin and database pair
    _tara_cift = 420                   # akis taramasi (degismedi)
    _cift = kutu_say * len(var)
    yaz(u'  ESTIMATED TIME: ~%s  (%d bins x %d databases = %d pairs; each pair streams'
        % (sure_metni(_cift * (_tara_cift + _hiz_cift)), kutu_say, len(var), _cift))
    yaz(u'  the whole database). Resumable if interrupted.')
    yaz(u'  SHORT LIST: %d candidates, ALL aligned. Alignment cost per pair ~%s'
        % (kl_ust, sure_metni(_hiz_cift)))
    yaz(u'  (average query %d bp); total alignment cost ~%s.'
        % (int(_oq), sure_metni(_cift * _hiz_cift)))
    yaz(u'  2026-08-05: the list was brought BACK DOWN from 1000 to %d. The ordering criterion is the reverse'
        % kl_ust)
    yaz(u'  frekans + uzunluk normalizasyonuna gecti; olculen en kotu kazanan sirasi')
    yaz(u'  35, so there is %d-fold headroom. The saving is ~%s (measured: with 1000'
        % (kl_ust // 35, sure_metni(_cift * _bir * kl_ust)))
    yaz(u'  candidates alignment was 53% of the total cost, with 500 it is 36%).')
    yaz(u'  Whether the cut-off is binding needs NO SEPARATE MEASUREMENT:')
    yaz(u'  the winning hit\'s short-list rank is recorded for every query.')
    yaz('')

    sonuc = []
    tb = time.time()
    for n, idd in enumerate(iddialar, 1):
        # THE CLAIM LEVEL CHECKPOINT carries the settings signature too. Without it:
        # even with vtb_*.json invalidated, iddia_01.json would be found, the claim
        # would be skipped ENTIRELY, and the old verdict produced with the 60 item
        # list would come back silently. (The NCBI nt cache, nt_*.json, is
        # DELIBERATELY excluded: it comes over the network, has nothing to do with
        # the short list size, and should not be re-queried for nothing.)
        # -------------------------------------------------------------------
        # THE 2026-08-21 ARCHITECTURAL FIX: MEASUREMENT and VERDICT must not be
        # held in the same cache.
        #
        # The bug found: this checkpoint stored the WHOLE claim, 'adlandirma' and
        # 'hukum' included. AYAR_IMZASI, meanwhile, seals only the SCAN
        # parameters. The result: when the naming logic was fixed, a re-run said
        # "taken from the previous run" and handed back THE OLD VERDICT.
        # Measured: in the round run right after the unnamed record fix, 12 of 12
        # claims were unchanged, because the corrected code NEVER RAN AT ALL.
        #
        # The right shape: what is expensive is THE SCAN (hours); naming and the
        # verdict are DERIVATIONS taking seconds. A derivation should be redone on
        # every run, and for that the checkpoint has to carry THE RAW HITS
        # ('bulgular').
        #
        # Old format checkpoints DO NOT carry raw hits; they cannot be re-derived
        # and count as INVALID. The old verdict never comes back silently.
        kp = os.path.join(KONTROL, 'iddia_%02d_%s_kl%d.json'
                          % (idd['no'], AYAR_IMZASI, kl_ust))
        if os.path.exists(kp):
            try:
                _kayit = json.load(open(kp, encoding='utf-8'))
            except Exception:
                _kayit = None
            if _kayit is not None and _kayit.get('_ham_isabet_var'):
                # There is a raw measurement: the SCAN is skipped, the VERDICT is re-derived.
                sonuc.append(_yeniden_turet(_kayit, idd, kons))
                yaz(u'[%2d/%2d] claim %d  (scan from cache, verdict RE-DERIVED)'
                    % (n, len(iddialar), idd['no']))
                continue
            if _kayit is not None:
                yaz(u'[%2d/%2d] claim %d  OLD-FORMAT checkpoint (no raw hits), the verdict cannot be re-derived, RESCANNING'
                    % (n, len(iddialar), idd['no']))
        yaz(u'[%2d/%2d] CLAIM %d (priority %d): %s'
            % (n, len(iddialar), idd['no'], idd['oncelik'], idd['metin']))
        if idd.get('not_'):
            yaz(u'         NOT: %s' % idd['not_'])

        bulgular = {}
        if idd['tip'] == 'gecici':
            kutular = []          # a statement of method; it needs no database scan
        else:
            kutular = (idd.get('kutu') or [])[:1] + (idd.get('karsi') or [])[:1]
        for kutu in kutular:
            if kutu not in kons:
                yaz(u'         no consensus: %s' % kutu)
                continue
            q = kons[kutu]
            if len(q) > 4000:
                q = q[:4000]
            for et, dosya, _tur in var:
                _gar = [x for x in (idd.get('beklenen_cins'),
                                    (idd.get('beklenen_tur') or '').split()[0]
                                    if idd.get('beklenen_tur') else None,
                                    idd.get('beklenen_aile')) if x]
                bulgular[et] = vtb_tarama(kok, q, et, dosya, yaz, KONTROL, _gar,
                                          kl_ust=kl_ust)
            # --- NCBI nt: yerel kumelerin hepsi belirli lokuslara ozeldir, nt en genisi ---
            if nt_kip != 'yok':
                if kutu in nt_onceden:
                    bulgular[NT_ETIKET] = nt_onceden[kutu]
                    yaz(u'     NCBI nt: used the manually loaded result')
                else:
                    ntk = os.path.join(KONTROL, 'nt_%s.json' % re.sub(r'\W+', '_', kutu))
                    onbellek = None
                    if os.path.exists(ntk):
                        try:                      # O-2: korumasiz json.load idi
                            onbellek = json.load(open(ntk, encoding='utf-8'))
                        except Exception as e:
                            yaz(u'     NCBI nt: the cache is corrupt, it will be retried (%s)'
                                % type(e).__name__)
                            onbellek = None
                    if onbellek and str(onbellek.get('durum', '')).startswith('TAMAM'):
                        bulgular[NT_ETIKET] = onbellek
                        yaz(u'     NCBI nt: taken from the previous run')
                    else:
                        # O-3: a network error IS NOT A RESULT. It is not cached and is
                        # retried on every run. Otherwise a single Wi-Fi dropout poisoned
                        # every claim permanently.
                        bulgular[NT_ETIKET] = nt_katmani(kutu, q, CIKTI, yaz, nt_kip)
                        if str(bulgular[NT_ETIKET].get('durum', '')).startswith('TAMAM'):
                            json.dump(bulgular[NT_ETIKET], open(ntk, 'w', encoding='utf-8'),
                                      ensure_ascii=False, default=str)
            break     # the first bin is enough; the rest are compared by alignment directly

        h, kanit, duzeltme = hukum_ver(idd, bulgular, kons)
        nt = bulgular.get(NT_ETIKET)
        nt_eksik = (nt_kip != 'yok' and idd['tip'] != 'gecici'
                    and (nt is None or not str(nt.get('durum', '')).startswith('TAMAM')))
        if nt_eksik:
            h = 'DOGRULANAMADI'
            kanit = (kanit + '  ||  THE NCBI nt LAYER DID NOT FINISH (%s), '
                             'which is why the claim counts as unconfirmed '
                             'rather than being skipped quietly. To finish it '
                             'by hand: BLAST the query file under '
                             'IDENTITY_RESULT/nt_elle/, write the result into '
                             'NT_SONUC_SABLONU.tsv, and give it back with '
                             '--nt-load.'
                     % ((nt or {}).get('durum', 'kosulmadi')))
        # --- NAMING: merge the hits from every database, take the best five
        # and the defensible taxonomic level (the short lists are already in
        # hand, so there is NO re-scan).
        adl, lokus = adlandirmayi_turet(bulgular)
        # --- THE LITERATURE CHECK (A REQUIRED STEP) ---
        try:
            import importlib.util as _lu
            _lp = _lu.spec_from_file_location(
                'lit', os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    'literature_check.py'))
            LIT = _lu.module_from_spec(_lp); _lp.loader.exec_module(LIT)
            _b1 = (adl.get('isabet1') or {}).get('tam_ad', '')
            lit = LIT.kontrol_et(_b1, adl.get('onerilen_ad', ''), lokus,
                                 ag=(lit_kip != 'yok'))
        except Exception as _e:
            lit = dict(durum=u'literatur modulu yuklenemedi: %s' % type(_e).__name__,
                       erisim_no='-', alan='-', vtb_adi='-', ncbi_guncel_ad='-',
                       es_anlamlilar='-', rutbe='-', ad_farkli_mi='-',
                       revizyon_uyarisi='-', revizyon_pmid='-',
                       otorite_kontrolu='GEREKLI', otorite_baglantilari='')
        lit['no'] = idd['no']; lit['onerilen_ad'] = adl.get('onerilen_ad', '-')
        yaz(u'         literature: %s | name differs: %s | revision: %s'
            % (lit.get('durum', '-'), lit.get('ad_farkli_mi', '-'),
               (lit.get('revizyon_uyarisi', '-') or '-')[:48]))

        detay = {}
        for et, v in bulgular.items():
            if str(v.get('durum', '')).startswith('TAMAM') and v.get('isabet'):
                en = en_iyi_isabet(v)
                detay[et] = dict(durum=v['durum'], en_iyi=en.get('baslik', '')[:120],
                                 kimlik=en.get('kimlik'),
                                 kazanan_sira=v.get('kazanan_sira'),
                                 kazanan_kaynak=v.get('kazanan_kaynak'),
                                 kisa_liste_boyu=v.get('kisa_liste_boyu'),
                                 taranan_kayit=v.get('taranan_kayit'))
            else:
                detay[et] = dict(durum=v.get('durum', '?'), en_iyi='', kimlik=None)
        # --- THE SELF CALIBRATION SUMMARY (at claim level) ---
        # The miss rate NEEDS NO SEPARATE MEASUREMENT: the data comes out of the
        # run itself. If the winners always come from the first 100 the cut off is
        # not binding; if any go past 400, even 500 may not be enough.
        _sr = [(e, v.get('kazanan_sira'), v.get('kazanan_kaynak'))
               for e, v in bulgular.items()
               if str(v.get('durum', '')).startswith('TAMAM')
               and isinstance(v.get('kazanan_sira'), int)]
        kal = dict(
            kisa_liste_boyu=kl_ust,
            en_yuksek_kazanan_sira=(max(s for _e, s, _k in _sr) if _sr else None),
            kazanan_sira_vtb=(max(_sr, key=lambda x: x[1])[0] if _sr else '-'),
            guvenli_bolge_disi=len([1 for _e, s, _k in _sr if s > SIRA_GUVENLI_BOLGE]),
            uyari_esigi_ustu=len([1 for _e, s, _k in _sr if s > SIRA_UYARI_ESIGI]),
            garanti_ile_kazanan=len([1 for _e, _s, k in _sr if k == 'garanti']),
            olculen_vtb=len(_sr),
            vtb_siralari={e: s for e, s, _k in _sr})
        yaz(u'         self-calibration: highest winner rank %s/%d (%s) | outside the first %d in %d/%d databases | won via guarantee %d'
            % (kal['en_yuksek_kazanan_sira'] if kal['en_yuksek_kazanan_sira'] is not None else '-',
               kl_ust, kal['kazanan_sira_vtb'], SIRA_GUVENLI_BOLGE,
               kal['guvenli_bolge_disi'], kal['olculen_vtb'], kal['garanti_ile_kazanan']))

        r = dict(no=idd['no'], oncelik=idd['oncelik'], iddia=idd['metin'], hukum=h,
                 adlandirma=adl, literatur=lit, vtb_detay=detay, kalibrasyon=kal,
                 kanit=kanit, dogru_ifade=duzeltme,
                 uyusan_vtb=len([1 for v in bulgular.values() if v.get('durum') == 'TAMAM']),
                 vtb=list(bulgular.keys()),
                 ayirt_edici={e: v.get('ayirt_edici') for e, v in bulgular.items()
                              if v.get('durum') == 'TAMAM'})
        # THE RAW MEASUREMENT IS WRITTEN TO THE CHECKPOINT (2026-08-21).
        # That way, when the naming or verdict logic changes, the derivation can be
        # refreshed WITHOUT RE-RUNNING the scan. Only the fields the derivation needs
        # are stored (the status plus the first 5 hits); the whole short list is not
        # kept, so the file does not swell for nothing.
        r['bulgular'] = {
            e: dict(durum=v.get('durum'),
                    isabet=[{k: i.get(k) for k in ('baslik', 'kimlik', 'hiz_uzunluk')}
                            for i in (v.get('isabet') or [])[:5]])
            for e, v in bulgular.items()}
        r['_ham_isabet_var'] = True
        json.dump(r, open(kp, 'w', encoding='utf-8'), ensure_ascii=False, default=str)
        sonuc.append(r)
        yaz(u'         => %s' % h)
        if duzeltme:
            yaz(u'         CORRECT WORDING: %s' % duzeltme)
        gec = time.time() - tb
        print(u'        elapsed %s | estimated remaining %s'
              % (sure_metni(gec), sure_metni(gec / n * (len(iddialar) - n))), flush=True)

    raporla(CIKTI, sonuc, var, yaz)
    rc = cikti_denetle(yaz, 'I (KIMLIK)', [
        (os.path.join(CIKTI, 'kimlik_iddialari.tsv'), 'kimlik_iddialari.tsv')])
    g.close()
    return rc


# -------------------------------------------------------------------------
# Three outputs: a one line TSV per claim, a manual literature check list and a
# markdown report.
#
# The "short list self calibration" section stands apart in the report: if every
# winner came from inside the first 100 the cut off is NOT binding and that is said
# openly; if any came from past 400 then 500 may not be enough either and raising
# --shortlist is suggested. The adequacy of the list size therefore proves or
# refutes itself on every run.
#
# The CORRECTION NEEDED rows stand at the TOP of the report: those are the ones
# that will reach the write-up.
# -------------------------------------------------------------------------
def raporla(CIKTI, sonuc, var, yaz):
    t = os.path.join(CIKTI, 'kimlik_iddialari.tsv')
    with open(t, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# Every claim was tested INDEPENDENTLY. VERIFIED means AT LEAST TWO databases agree.\n')
        fh.write(u'# Method: seed + alignment. NO taxonomy tree, NO k-mer LCA, NO primers.\n')
        try:
            import sys as _s2, os as _o2
            _k2 = _o2.path.dirname(_o2.path.dirname(_o2.path.abspath(__file__)))
            if _k2 not in _s2.path:
                _s2.path.insert(0, _k2)
            from screening import labels as _L2
            fh.write(_L2.verdict_legend(_L2.KIMLIK, 'VERDICTS'))
            fh.write(_L2.verdict_legend(_L2.DUZEY, 'DEFENSIBLE TAXONOMIC LEVEL'))
            fh.write(_L2.legend(_L2.SUTUN.keys()))
        except Exception:
            pass
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['no', 'oncelik', 'iddia', 'HUKUM',
                    'SAVUNULABILIR_DUZEY', 'ONERILEN_AD', 'adlandirma_gerekcesi',
                    'en_iyi_isabet', 'en_iyi_cins', 'en_iyi_tur', 'en_iyi_kimlik_%',
                    'en_iyi_hiz_uzunluk', 'ADI_URETEN_HIZALAMA', 'en_iyi_vtb',
                    'ikinci_isabet', 'ikinci_tur', 'ikinci_kimlik_%', 'ikinci_vtb',
                    'ucuncu_isabet', 'ucuncu_tur', 'ucuncu_kimlik_%', 'ucuncu_vtb',
                    'dorduncu_isabet', 'dorduncu_tur', 'dorduncu_kimlik_%', 'dorduncu_vtb',
                    'besinci_isabet', 'besinci_tur', 'besinci_kimlik_%', 'besinci_vtb',
                    'EN_YAKIN_5_ORGANIZMA',
                    'LIT_veritabani_adi', 'LIT_ncbi_guncel_ad', 'LIT_es_anlamlilar',
                    'LIT_rutbe', 'LIT_AD_FARKLI_MI', 'LIT_revizyon_uyarisi',
                    'LIT_revizyon_pmid', 'LIT_otorite_kontrolu', 'LIT_baglantilar',
                    'LIT_durum',
                    'kisa_liste_boyu', 'kazanan_sira', 'kazanan_sira_vtb',
                    'kazanan_garanti_ile_mi', 'kesme_baglayici_mi',
                    'sorgulanan_vtb_sayisi', 'sonuc_veren_vtb', 'HER_VTB_NE_DEDI',
                    'kanit', 'DOGRU_IFADE (duzeltilmeliyse)'])
        for s in sonuc:
            d = s.get('vtb_detay') or {}
            hepsi = ' | '.join('%s [%s records scanned]: %s%s%s'
                               % (e, v.get('taranan_kayit', '?'),
                                  v['en_iyi'] or v['durum'],
                                             ('' if v.get('kimlik') is None
                                              else ' (%%%s)' % vir(v['kimlik'])),
                                             ('' if v.get('kazanan_sira') is None
                                              else ' {sira %s/%s%s}'
                                              % (v['kazanan_sira'], v.get('kisa_liste_boyu', '?'),
                                                 ', GARANTI'
                                                 if v.get('kazanan_kaynak') == 'garanti' else '')))
                               for e, v in d.items())
            kal = s.get('kalibrasyon') or {}
            _eys = kal.get('en_yuksek_kazanan_sira')
            if _eys is None:
                _bag = '-'
            elif kal.get('uyari_esigi_ustu'):
                _bag = u'EVET-OLABILIR (>%d)' % SIRA_UYARI_ESIGI
            elif kal.get('guvenli_bolge_disi'):
                _bag = u'HAYIR (ama ilk %d disindan geldi)' % SIRA_GUVENLI_BOLGE
            else:
                _bag = 'NO (all of them are inside the first %d)' % SIRA_GUVENLI_BOLGE
            a = s.get('adlandirma') or {}
            def _i(n, alan, vars_=''):
                return ((a.get('isabet%d' % n) or {}).get(alan) or vars_)
            w.writerow([s['no'], s['oncelik'], s['iddia'], s['hukum'],
                        a.get('duzey', '-'), a.get('onerilen_ad', '-'), a.get('gerekce', '-'),
                        _i(1, 'tam_ad', '-'), _i(1, 'cins', '-'), _i(1, 'tur', '-'),
                        vir(_i(1, 'kimlik', None)), _i(1, 'uzunluk', '-'),
                        # The alignment THAT PRODUCED THE NAME can differ from the
                        # best hit (when the best hit is unnamed the name comes
                        # from the NAMED hit at the same closeness). If the two
                        # numbers part company, a reader has to see it.
                        a.get('adi_ureten_hizalama') or '-',
                        _i(1, 'vtb', '-'),
                        _i(2, 'tam_ad', '-'), _i(2, 'tur', '-'),
                        vir(_i(2, 'kimlik', None)), _i(2, 'vtb', '-'),
                        _i(3, 'tam_ad', '-'), _i(3, 'tur', '-'),
                        vir(_i(3, 'kimlik', None)), _i(3, 'vtb', '-'),
                        _i(4, 'tam_ad', '-'), _i(4, 'tur', '-'),
                        vir(_i(4, 'kimlik', None)), _i(4, 'vtb', '-'),
                        _i(5, 'tam_ad', '-'), _i(5, 'tur', '-'),
                        vir(_i(5, 'kimlik', None)), _i(5, 'vtb', '-'),
                        # A readable summary in one cell: ordered, with the identity percentage.
                        # For an UNNAMED record the species and genus are stored as '-'. If '-' is
                        # not treated as EMPTY the row comes out as "- 99.00%" and tells the reader
                        # nothing, when what they actually want to know is WHAT IT MATCHED. So '-'
                        # is taken as empty, the code falls back to the full header, and the fact
                        # that it is unnamed is written OPENLY.
                        ' | '.join(
                            '%d) %s %%%s [%s]'
                            % (n_, _en_yakin_etiket(a.get('isabet%d' % n_)),
                               vir(_i(n_, 'kimlik', None)), _i(n_, 'vtb', '-'))
                            for n_ in range(1, 6) if a.get('isabet%d' % n_)) or '-',
                        (s.get('literatur') or {}).get('vtb_adi', '-'),
                        (s.get('literatur') or {}).get('ncbi_guncel_ad', '-'),
                        (s.get('literatur') or {}).get('es_anlamlilar', '-'),
                        (s.get('literatur') or {}).get('rutbe', '-'),
                        (s.get('literatur') or {}).get('ad_farkli_mi', '-'),
                        (s.get('literatur') or {}).get('revizyon_uyarisi', '-'),
                        (s.get('literatur') or {}).get('revizyon_pmid', '-'),
                        (s.get('literatur') or {}).get('otorite_kontrolu', 'GEREKLI'),
                        (s.get('literatur') or {}).get('otorite_baglantilari', ''),
                        (s.get('literatur') or {}).get('durum', '-'),
                        kal.get('kisa_liste_boyu', '-'),
                        '-' if _eys is None else _eys,
                        kal.get('kazanan_sira_vtb', '-'),
                        ('EVET' if kal.get('garanti_ile_kazanan') else 'hayir')
                        if kal else '-',
                        _bag,
                        len(d), s['uyusan_vtb'], hepsi, s['kanit'], s['dogru_ifade']])
    yaz(u'  written: %s' % t)
    try:
        import importlib.util as _lu
        _lp = _lu.spec_from_file_location(
            'lit', os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'literature_check.py'))
        LIT = _lu.module_from_spec(_lp); _lp.loader.exec_module(LIT)
        el = LIT.elle_liste_yaz(CIKTI, [dict((s.get('literatur') or {}),
                                             no=s['no'],
                                             onerilen_ad=(s.get('adlandirma') or {})
                                             .get('onerilen_ad', '-')) for s in sonuc])
        yaz(u'  written: %s' % el)
    except Exception as e:
        yaz(u'  could not write the manual check list: %s' % type(e).__name__)

    say = {}
    for s in sonuc:
        say[s['hukum']] = say.get(s['hukum'], 0) + 1
    r = os.path.join(CIKTI, 'KIMLIK_DOGRULAMA_RAPORU.md')
    with open(r, 'w', encoding='utf-8') as fh:
        fh.write(u'# Independent verification of identity claims\n\n')
        fh.write(u'Generated: %s, script %s\n\n' % (time.strftime('%Y-%m-%d %H:%M'), VERSIYON))
        fh.write(u'Local databases used (%d): %s\n\n'
                 % (len(var), ', '.join(e for e, _, _ in var)))
        fh.write(u'Additionally **NCBI nt** is queried as a separate layer. All the local sets are specific to particular loci (S')
        fh.write(u'Full set inventory, and why some sets are unused: `VERITABANI_ENVANTERI.md`\n\n')
        fh.write(u'## Result\n\n')
        for k in ('DOGRULANDI', 'DUZELTILMELI', 'DOGRULANAMADI'):
            if k in say:
                fh.write(u'- **%s**: %d\n' % (k, say[k]))

        # ------------------------------------------- OZ KALIBRASYON OZETI
        kals = [s.get('kalibrasyon') or {} for s in sonuc]
        siralar = [(s['no'], e, v) for s in sonuc
                   for e, v in ((s.get('kalibrasyon') or {}).get('vtb_siralari') or {}).items()
                   if isinstance(v, int)]
        fh.write(u'\n## Short-list self-calibration (is the cut-off binding?)\n\n')
        if not siralar:
            fh.write(u'No measurable winner rank (this run may not have scanned any database).\n\n')
        else:
            ust = max(v for _n, _e, v in siralar)
            n100 = len([1 for _n, _e, v in siralar if v > SIRA_GUVENLI_BOLGE])
            n400 = len([1 for _n, _e, v in siralar if v > SIRA_UYARI_ESIGI])
            gar = sum(k.get('garanti_ile_kazanan', 0) for k in kals)
            boy = next((k.get('kisa_liste_boyu') for k in kals
                        if k.get('kisa_liste_boyu')), KISA_LISTE)
            fh.write(u'The short list is built on **seed count**, but the decision is made by **alignment identity**. Because the two criteria differ')
            fh.write(u'| measure | value |\n|---|---|\n')
            fh.write(u'| short-list size | %d (ALL aligned) |\n' % boy)
            fh.write(u'| queries measured (claim x database) | %d |\n' % len(siralar))
            fh.write(u'| **highest winner rank** | **%d** |\n' % ust)
            fh.write(u'| winner came from outside the first %d | %d / %d queries |\n'
                     % (SIRA_GUVENLI_BOLGE, n100, len(siralar)))
            fh.write(u'| winner came from beyond rank %d | %d / %d queries |\n'
                     % (SIRA_UYARI_ESIGI, n400, len(siralar)))
            fh.write(u'| winner entered via the "expected taxon guarantee" | %d queries |\n\n' % gar)
            if n400:
                fh.write(u'> **WARNING.** In %d queries the winner came from beyond position %d. The cut off may still be binding; raise `--shortlist` (to %d, say) and repeat'
                         % (n400, SIRA_UYARI_ESIGI, boy * 2))
            elif n100:
                fh.write(u'> %d of the winners came from outside the first %d, which '
                         u'means **the old list of 60 WOULD HAVE MISSED those hits**. '
                         u'All of them stayed under %d, so a list of %d is enough.\n\n'
                         % (n100, SIRA_GUVENLI_BOLGE, SIRA_UYARI_ESIGI, boy))
            else:
                fh.write(u'> Every winner came from within the first %d. The cut-off is **not binding**: alignment decides entirely' % SIRA_GUVENLI_BOLGE)
            if gar:
                fh.write(u'> **Caution.** In %d queries the winner entered the short list not by seeding but through the "expected taxon guarantee" patch' % gar)
            enb = sorted(siralar, key=lambda x: -x[2])[:5]
            fh.write(u'Five highest winner ranks: %s\n\n'
                     % '; '.join(u'claim %d / %s: %d' % t for t in enb))

        fh.write(u'\n> **Difference in method.** Kraken2 does k-mer plus LCA on a taxonomy tree; our earlier rounds compared consensus sequences')
        for etiket, baslik in (('DUZELTILMELI', '## What has to be corrected, for the report'),
                               ('DOGRULANDI', u'## Dogrulananlar'),
                               ('DOGRULANAMADI', u'## Dogrulanamayanlar')):
            uy = [s for s in sonuc if s['hukum'] == etiket]
            if not uy:
                continue
            fh.write(baslik + u'\n\n')
            for s in uy:
                fh.write(u'### %d. %s\n\n' % (s['no'], s['iddia']))
                fh.write(u'- **Verdict:** %s  (agreeing databases: %d)\n' % (s['hukum'], s['uyusan_vtb']))
                a = s.get('adlandirma') or {}
                if a:
                    fh.write(u'- **Savunulabilir duzey:** `%s` → **%s**\n'
                             % (a.get('duzey', '-'), a.get('onerilen_ad', '-')))
                    fh.write('  - *The reason:* %s\n' % a.get('gerekce', '-'))
                    fh.write(u'\n  **The five nearest ORGANISMS** (deduplicated by organism, not by record, because the same species appears in several databases')
                    fh.write(u'\n  | # | nearest record | genus | species | identity | database |\n  |---|---|---|---|---|---|\n')
                    for n_ in (1, 2, 3, 4, 5):
                        it = a.get('isabet%d' % n_)
                        if it:
                            fh.write(u'  | %d | %s | %s | %s | %%%s | %s |\n'
                                     % (n_, it['tam_ad'], it['cins'], it['tur'],
                                        vir(it['kimlik']), it['vtb']))
                    fh.write(u'\n  > **Having a name and claiming an identity are different things.** The "nearest record" above is NOT AN IDENTITY')
                fh.write(u'- **Evidence:** %s\n' % s['kanit'])
                d = s.get('vtb_detay') or {}
                if d:
                    fh.write(u'\n  **Databases queried and what each one said (%d sources):**\n\n' % len(d))
                    fh.write(u'  | database | status | best hit | identity | winner rank |\n  |---|---|---|---|---|\n')
                    for e, v in d.items():
                        _ks = ('-' if v.get('kazanan_sira') is None
                               else u'%d / %s%s' % (v['kazanan_sira'],
                                                    v.get('kisa_liste_boyu', '?'),
                                                    u' **(GARANTI)**'
                                                    if v.get('kazanan_kaynak') == 'garanti'
                                                    else ''))
                        fh.write(u'  | %s | %s | %s | %s | %s |\n'
                                 % (e, v['durum'], v['en_iyi'] or '-',
                                    '-' if v.get('kimlik') is None else '%%%s' % vir(v['kimlik']),
                                    _ks))
                    fh.write(u'\n  > *Winning position*: the SEED position, in the short list, of the record that gave the highest identity. Small numbers show that the cut off is not binding.\n\n')
                if s['dogru_ifade']:
                    fh.write(u'- **DOGRU IFADE:** %s\n' % s['dogru_ifade'])
                ae = {k: v for k, v in (s.get('ayirt_edici') or {}).items() if v}
                if ae:
                    fh.write(u'- **Ayirt edici pencere:** %s\n'
                             % '; '.join('%s yayilim %%%s' % (k, vir(v.get('yayilim')))
                                         for k, v in ae.items()))
                fh.write(u'\n')
    yaz(u'  written: %s' % r)
    yaz('')
    yaz('  ' + '   '.join('%s: %d' % kv for kv in say.items()))



# --------------------------------------------------------------- guvenlik agi
def cikti_denetle(yaz, ad, dosyalar, asgari=1):
    """When the stage ends, it audits ITS OWN output.

        If the expected row count is zero, or the file is missing entirely, it DOES
        NOT CARRY ON SILENTLY: it prints a clear error and returns a non-zero code.
        This is so that it cannot produce an empty result overnight and then say
        "nothing was found" in the morning.

    """
    sorun = []
    for yol, etiket in dosyalar:
        if not os.path.exists(yol):
            sorun.append(u'%s WAS NOT PRODUCED (%s)' % (etiket, yol)); continue
        try:
            with open(yol, encoding='utf-8') as fh:
                n = sum(1 for s in fh if s.strip() and not s.startswith('#'))
            n = max(0, n - 1)          # baslik satiri
        except OSError as e:
            sorun.append(u'%s COULD NOT BE READ (%s)' % (etiket, e)); continue
        if n < asgari:
            sorun.append(u'%s IS EMPTY, %d data rows (at least %d were expected)'
                         % (etiket, n, asgari))
    if not sorun:
        return 0
    yaz('')
    yaz('  ' + '!' * 70)
    yaz(u'  STAGE %s PRODUCED EMPTY OUTPUT - THE CHAIN WAS STOPPED HERE' % ad)
    for x in sorun:
        yaz(u'    - %s' % x)
    yaz('')
    yaz(u'  WHY IT STOPPED: the next stage would have read this file as input.')
    yaz(u'  Continuing with empty input does not crash; it produces a MEANINGLESS BUT')
    yaz(u'  CONVINCING summary, which is exactly the silent failure we hunt for.')
    yaz(u'  Read the run log above, fix the cause, then run the same command')
    yaz(u'  again; finished work is skipped from its checkpoints.')
    yaz('  ' + '!' * 70)
    return 4


def girdi_denetle(yaz, ad, dosyalar):
    """Before the stage STARTS: do the files it needs exist, and are they non-empty?"""
    eksik = []
    for yol, etiket, uretici in dosyalar:
        if not os.path.exists(yol):
            eksik.append(u'there is no %s (%s) -> run stage %s first' % (etiket, yol, uretici))
            continue
        with open(yol, encoding='utf-8') as fh:
            n = sum(1 for s in fh if s.strip() and not s.startswith('#'))
        if n <= 1:
            eksik.append(u'%s IS EMPTY (%s) -> stage %s produced no result'
                         % (etiket, yol, uretici))
    if not eksik:
        return 0
    yaz('')
    yaz('  ' + '!' * 70)
    yaz(u'  STAGE %s WAS NOT STARTED - INPUT MISSING' % ad)
    for x in eksik:
        yaz(u'    - %s' % x)
    yaz('  ' + '!' * 70)
    return 5

# The command line: --only a single claim, --db-max how many databases,
# --shortlist how many candidates are fully aligned (changing it invalidates the
# old checkpoints), --nt the NCBI mode, --nt-load a hand filled template,
# --literature, --reset.

# --- CLI value normalisation ------------------------------------------------
# English option values are accepted alongside the original Turkish ones and
# mapped back here. The internal values are unchanged on purpose: they are
# compared in dozens of places and, in some cases, written to output files.
# Translating the interface must not translate the data.
_DEGER = {'auto': 'oto', 'manual': 'elle', 'none': 'yok', 'quick': 'hizli',
          'full': 'tam', 'member': 'uye', 'competitor': 'rakip',
          'exclude': 'disla'}


def _ing_deger(a):
    for _ad in ('nt', 'literatur', 'ncbi', 'karisik', 'moduller', 'mod'):
        _v = getattr(a, _ad, None)
        if isinstance(_v, str) and _v in _DEGER:
            setattr(a, _ad, _DEGER[_v])
    return a

def main():
    p = argparse.ArgumentParser(description='Confirm the identity claims '
                                            'independently')
    p.add_argument('--root', dest='kok', default='.')
    p.add_argument('--only', dest='yalniz', default=None, help='iddia numarasi or metninden parca')
    p.add_argument('--db-max', dest='vtb_ust', type=int, default=len(VTB), help='how many databases to use')
    p.add_argument('--shortlist', type=int, default=KISA_LISTE, dest='kisa_liste',
                   help=('how many candidates from each database to align IN '
                         'FULL (the default is %d). A larger value keeps the '
                         'cut-off from binding; the winner rank column in the '
                         'report says whether it was enough. Changing it '
                         'invalidates the old checkpoints.'
                         % KISA_LISTE))
    p.add_argument('--nt', choices=['auto', 'manual', 'none', 'oto', 'elle', 'yok'], default='oto',
                   help='NCBI nt layer: auto (URL API), manual (write a query file), none')
    p.add_argument('--nt-load', dest='nt_yukle', default=None,
                   help='doldurulmus NT_SONUC_SABLONU.tsv')
    p.add_argument('--literature', dest='literatur', choices=['auto', 'none', 'oto', 'yok'], default='oto',
                   help='NCBI Taxonomy + PubMed literatur kontrolu')
    p.add_argument('--reset', dest='sifirla', action='store_true')
    a = p.parse_args()
    a = _ing_deger(a)
    kok = os.path.abspath(a.kok)
    if not os.path.isdir(os.path.join(kok, 'screening')):
        sys.exit(u'ERROR: there is no screening directory inside %s.' % kok)
    if a.kisa_liste < 1:
        sys.exit(u'ERROR: --shortlist must be at least 1.')
    return calistir(kok, a.yalniz, a.sifirla, a.vtb_ust, a.nt, a.nt_yukle, a.literatur,
                    kl_ust=a.kisa_liste)


if __name__ == '__main__':
    sys.exit(main() or 0)
