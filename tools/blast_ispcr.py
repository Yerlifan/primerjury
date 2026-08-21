#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IN-SILICO PCR INSIDE THE SAMPLE  (version 3)

WHAT IT DOES
It uses the nanopore reads we have as if they were a real PCR template. For every
primer pair it counts whether a product forms in the reads of every organism in
the sample. It is expected to be high on the target and near zero in competitors.
The decision is made on the competitor ratio.

WHY THERE WERE THREE VERSIONS, AND WHAT WAS WRONG IN EACH
  Version 1 counted a BLAST hit as a binding. blastn-short reports a local
    alignment; whether the primer's 3' end is inside it is not asked. So a hit
    covering only the middle of the primer counted as a binding, and the numbers
    came out too high.

  Version 2 added the physical rules but did not account for BLAST's behaviour.
    The rule "the alignment must include the primer's 3' end" was added. But
    blastn-short reports the HIGHEST SCORING alignment, not the whole alignment.
    It gives minus three for a mismatch and plus one for a match. If there is a
    mismatch at the third base from the end, including the last three bases lowers
    the score, so BLAST drops them and reports q1-17.
    The real situation: the last two bases match exactly, there is a single
    difference in the middle, and the primer binds and extends.
    The code said "the 3' end is not in the alignment" and ignored the binding.
    An example, the Podospora pair on its own target:
        read   ACTCGTCGAAGGAGCTTTAC
        primer ACTCGTCGAAGGAGCTTCAC
    BLAST reported this as q1-17 and version 2 said "no binding". Wrong.

VERSION 3'S ANSWER: BRUTE FORCE
Our database is 23.6 million bases. BLAST exists to avoid brute force on databases
of billions of bases; at this scale there is no such necessity. So EVERY position
on every read is now tried one at a time. No seed, no scoring, no truncation, no
misses. Because the comparison is vectorised with numpy, the whole sample is
scanned in about half a minute. blastn is kept as a separate opinion (the --blast
option) but it does not decide.

THE PRODUCT RULE, the physics of PCR
  For each primer
    a. the last two bases must match the target exactly (extension starts there)
    b. at most one mismatch in the last five bases
    c. at most three mismatches over the whole primer
  For the two primers together
    d. one must bind the plus strand and the other the minus strand
    e. their 3' ends must FACE ONE ANOTHER
    f. the distance between them must be inside the product range

THE DECISION
The ratios are found by dividing by the read count, but a raw ratio misleads on a
small sample: in a taxon with seventeen reads, one read is six percent. So the
Wilson lower bound is used in the decision, that is, the lowest value of the ratio
that is statistically defensible.
Two separate axes are reported:
  CROSS REACTION  does a product form in the competitors
  COVERAGE        does a product form in all of the targeted organisms
When coverage is evaluated, each taxon's "reachability" is measured from its own
data: if that taxon gives a high ratio with any primer pair, its reads are sound;
if another pair gives a low ratio in the same taxon, that is not a shortage of
reads but a fault of that pair.

Run:
  python3 blast_ispcr.py --root /path/to/project
  python3 blast_ispcr.py --selftest      (a known-answer test of every rule)
  python3 blast_ispcr.py --root ... --blast   (a second opinion via blastn, slow)

"""
import argparse, csv, glob, math, os, random, re, shutil, statistics, subprocess, sys, tempfile, time
from collections import defaultdict

try:
    import numpy as np
except ImportError:
    sys.exit(u'required: pip install numpy')

RC = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
def rc(s): return "".join(RC.get(b,'N') for b in reversed(s.upper()))
def kod(s): return np.frombuffer(s.encode(), dtype=np.uint8)

LOGF = None
def log(m):
    s = f"[{time.strftime('%H:%M:%S')}] {m}"
    print(s, flush=True)
    if LOGF: LOGF.write(s + "\n"); LOGF.flush()

def arac_var(ad): return shutil.which(ad) is not None
def guvenli(ad): return re.sub(r"[^A-Za-z0-9]+", "_", ad).strip("_")

AYIRAC = 0          # a byte placed between reads that equals no base

# Which taxon each primer pair targets. Where a list is given, that is a group
# primer and a product is expected in EVERY member of the list.
HEDEF_TAXID = {
    "AD01_mazei":              ["2209"],
    "AD01b_Mmazei_tur":        ["2209"],
    "AD02_soehngenii":         ["2223"],
    "AD05_Bthetaiotaomicron":  ["818"],
    "AD06_Alistipes":          ["214856"],
    "AD10_Methanosarcinaceae": ["2209","2208","2210","3078083","1434102"],
    "AD11_Methanomicrobiaceae":["83986","118126","394967","2201"],
    "AD12_Tulosesus_callinus": ["181762"],
    "AD13_Oceanobacillus":     ["2954507"],
    "AD14_asetoklastik":       ["2209","2208","2210","3078083","1434102","2223"],
    "AD15_Methanoculleus":     ["83986","118126","394967"],
    "AD16_Microascaceae":      ["101201"],
    "AD17_Mbarkeri":           ["2208"],
    "AD18_Zoopagomycota":      ["44689"],
    "Methanosarcinaceae_ailesi":  ["2209","2208","2210","3078083","1434102"],
    "Methanomicrobiaceae_ailesi": ["83986","118126","394967","2201"],
    "Dysgonomonadaceae_ailesi":   ["2829812","1642647","1642646","285070"],
    "Bacteroidaceae_ailesi":      ["818","28116"],
    "Asetoklastik_metanojenler_asetattan_metan": ["2209","2208","2210","3078083","1434102","2223"],
    # 224719 (Methanobrevibacter sp. AbM4) was added on 28 July. It has 4723 reads in
    # the sample and is an unarguable H2/CO2 methanogen, but it was not in the member
    # list; the group was therefore answering "full coverage" to an incomplete question.
    # It must match the design stage.
    "Hidrojenotrofik_metanojenler_H2_CO2_den_metan": ["118126","83986","394967","2201","83984","224719"],
    "Metilotrofik_metanojen_metil_bilesiklerinden_metan": ["1406512"],
    "Sakarolitik_bakteriler_seker_parcalayan": ["818","28116","214856"],
    "Podospora_pseudopauciseta": ["2093780"],
    # Added on 28 July. The sample's counterpart to the meeting's "proteolytic /
    # syntrophic" request. Two separate phyla, not expected to be covered by one pair;
    # the group was DEFINED so that the measurement could say so. It must match the
    # design stage.
    "Proteolitik_sintrofik_bakteriler": ["456827","1197717"],
}

ISIMLER = {
    "44689":"Zoopagomycota mantari","2209":"Methanosarcina mazei","40559":"Helotiales askomikot",
    "83984":"Methanocorpusculum labreanum","2208":"Methanosarcina barkeri",
    "2223":"Methanothrix soehngenii","1826872":"Ca. Nitrosocosmicus hydrocola",
    "2740404":"Ca. Sulfurimonas baltica","456827":"Ca. Cloacimonas acidaminovorans",
    "101201":"Microascaceae askomikot","2829812":"Proteiniphilum propionicum",
    "1642646":"Petrimonas mucosa","1760811":"Pseudodifflugia amip","1129264":"Sphaerochaeta associata",
    "1642647":"Proteiniphilum saccharofermentans","500148":"Cyrtohymena siliat",
    "2233851":"Blastochloris tepida","2954507":"Oceanobacillus sp.","285070":"Petrimonas sulfuriphila",
    "224719":"Methanobrevibacter sp. AbM4","80884":"Colletotrichum higginsianum",
    "394967":"Methanoculleus receptaculi","1197717":"Cloacibacillus porcorum",
    "3078083":"Methanosarcina hadiensis","1159556":"Ustilaginoidea virens",
    "2093779":"Podospora pseudocomata","5855":"Vampyrellid amip","214856":"Alistipes finegoldii",
    "181762":"Tulosesus callinus","63577":"Trichoderma atroviride",
    "1406512":"Ca. Methanomassiliicoccus intestinalis","456999":"Rhizoctonia solani",
    "2210":"Methanosarcina thermophila","818":"Bacteroides thetaiotaomicron",
    "2034170":"Trichoderma breve","2093780":"Podospora pseudopauciseta",
    "2545709":"Schizosaccharomyces osmophilus","5811":"Toxoplasma gondii",
    "4896":"Schizosaccharomyces pombe","83986":"Methanoculleus bourgensis",
    "118126":"Methanoculleus chikugoensis","28116":"Bacteroides ovatus",
    "1434102":"Methanosarcina sp. WH1","2201":"Methanofollis liminatans",
}
def isim(tx): return ISIMLER.get(tx, f"taxid {tx}")

# ------------------------------------------------------------------ veri yukleme
def okumalari_yukle(kok, okuma_basina=300, en_az_uzunluk=300):
    """Takes a fixed number of reads from each taxon. Reproducible, because the same seed
        is used.

    """
    dosyalar = defaultdict(list)
    for p in glob.glob(f"{kok}/SONUCLAR/fastq files/*/*reads_*.fastq"):
        m = re.search(r"reads_(\d+)\.fastq$", os.path.basename(p))
        if m: dosyalar[m.group(1)].append(p)
    veri = {}
    for tx, yollar in sorted(dosyalar.items()):
        okumalar = []
        for p in sorted(yollar):
            cap = okuma_basina * 4; tut = []
            with open(p, errors="replace") as fh:
                while len(tut) < cap:
                    h = fh.readline()
                    if not h: break
                    s = fh.readline().strip().upper(); fh.readline(); fh.readline()
                    if len(s) >= en_az_uzunluk: tut.append(s)
            okumalar += tut
        random.Random(7).shuffle(okumalar)
        if okumalar: veri[tx] = okumalar[:okuma_basina]
    return veri

def primerleri_oku(yollar):
    """
        It reads both file formats:
          A) 'Oligo adi' / 'Dizi 5-3' columns, names ending _F and _R  (the order lists)
          B) 'grup' / 'F' / 'R' columns                              (the group and universal list)

    """
    ciftler = {}
    for yol in yollar:
        if not os.path.exists(yol): continue
        with open(yol, encoding="utf-8-sig") as fh:
            oku = csv.DictReader(fh)
            basliklar = [b.strip() for b in (oku.fieldnames or [])]
            b_bicimi = ("F" in basliklar and "R" in basliklar)
            for r in oku:
                r = {(k or "").strip(): (v or "").strip() for k, v in r.items()}
                if b_bicimi:
                    ad = guvenli(r.get("grup") or r.get("tip") or "")
                    F = r.get("F", "").upper(); R = r.get("R", "").upper()
                    if not ad or not F or not R: continue
                    if set(F) - set("ACGT") or set(R) - set("ACGT"): continue
                    ciftler[ad] = {"F": F, "R": R}
                else:
                    ad = r.get("Oligo adi") or r.get("ad") or ""
                    dizi = (r.get("Dizi 5-3") or r.get("dizi") or "").upper()
                    if not ad or not dizi: continue
                    taban, rol = ad.rsplit("_", 1)[0], ad.rsplit("_", 1)[-1].upper()
                    if rol not in ("F", "R"): continue
                    if set(dizi) - set("ACGT"): continue
                    ciftler.setdefault(guvenli(taban), {})[rol] = dizi
    return {k: v for k, v in ciftler.items() if "F" in v and "R" in v}

# ------------------------------------------------------------------ kaba kuvvet tarama
def kalip_kur(okumalar_listesi, en_uzun_primer):
    """
        It appends every read into one byte string, putting separators between them that
        equal no base. That way the whole sample can be scanned in one vectorised pass and
        no window can span two reads.
        Returns: (big_string, start_indices)

    """
    ayirac = np.full(en_uzun_primer, AYIRAC, dtype=np.uint8)
    parcalar = []; baslangic = []; p = 0
    for s in okumalar_listesi:
        a = kod(s)
        baslangic.append(p)
        parcalar.append(a); parcalar.append(ayirac)
        p += len(a) + en_uzun_primer
    if not parcalar: return np.zeros(0, dtype=np.uint8), np.zeros(0, dtype=np.int64)
    return np.concatenate(parcalar), np.array(baslangic, dtype=np.int64)

def baglanma_yerleri(W, primer, ters, en_fazla_uyumsuz=3, uc_pencere=5, uc_pencere_uyumsuz=1):
    """
        Finds every binding site of the primer. There is no seed; every window is tried.
        ters=False : the primer sits on the template directly (the plus strand), and its
                     3' end is at the END of the window
        ters=True  : the primer's reverse complement sits (the minus strand), and its
                     3' end is at the START of the window
        Returns: the window start indices

    """
    k = len(primer)
    pk = kod(rc(primer) if ters else primer)
    if W is None or W.shape[0] == 0 or W.shape[1] < k: return np.zeros(0, dtype=np.int64)
    Wk = W[:, :k]
    # ucuz on eleme: 3' uctaki iki baz birebir uymayan pencereler hemen elenir
    if ters: ank = (Wk[:, 0] == pk[0]) & (Wk[:, 1] == pk[1])
    else:    ank = (Wk[:, k-2] == pk[k-2]) & (Wk[:, k-1] == pk[k-1])
    idx = np.nonzero(ank)[0]
    if idx.size == 0: return idx
    e = (Wk[idx] == pk)
    uc = e[:, :uc_pencere] if ters else e[:, k-uc_pencere:]
    uc_hata = (~uc).sum(1); toplam_hata = (~e).sum(1)
    return idx[(uc_hata <= uc_pencere_uyumsuz) & (toplam_hata <= en_fazla_uyumsuz)]

def tek_primer_orani(W, baslangic, n_okuma, primer):
    """
        Returns in how many reads the primer finds a binding site ON ITS OWN.
        It is needed to tell apart the reasons when no product comes out on a target:
          both primers low        -> that region is not in the reads, and it is not the
                                     primer's fault
          one low and one high    -> that primer does not fit this organism, a real miss

    """
    yer = np.concatenate([baglanma_yerleri(W, primer, False),
                          baglanma_yerleri(W, primer, True)])
    if yer.size == 0: return 0.0
    okuma_no = np.searchsorted(baslangic, yer, side="right") - 1
    return len(set(okuma_no.tolist())) / max(1, n_okuma)

def kapsayan_okuma(W, baslangic, F, R):
    """
        The number of distinct reads AT LEAST ONE of the two primers binds.

        WHY IT IS NEEDED
        The cross reaction ratio was being divided by all of the competitor's reads.
        Because the number of reads covering the locus was computed nowhere, the threshold
        turned in practice into "a product in at most this many reads of the competitor".
        Such a criterion cannot see the difference between "a product in 5 of the 5 reads
        covering the locus" and "a product in 5 of the 130 reads covering the locus"; the
        first means full amplification in the competitor and must not be written CLEAN.

        This number is the denominator: the number of reads in which the product could
        appear. If neither primer binds in a read, that read never reached the locus and
        has no place in the denominator.

    """
    yer = np.concatenate([baglanma_yerleri(W, F, False), baglanma_yerleri(W, F, True),
                          baglanma_yerleri(W, R, False), baglanma_yerleri(W, R, True)])
    if yer.size == 0: return 0
    return len(set((np.searchsorted(baslangic, yer, side="right") - 1).tolist()))

def urun_boyu(arti_5, kp, eksi_3, kq, en_kisa, en_uzun):
    """
        arti_5 : the 5' end of the primer binding the plus strand (the window start); its
                 3' end is arti_5+kp-1
        eksi_3 : the 3' end of the primer binding the minus strand (the window start); its
                 5' end is eksi_3+kq-1
        For a product the two ends must face one another and the distance must be in range.

    """
    if arti_5 + kp - 1 > eksi_3: return None
    boy = (eksi_3 + kq - 1) - arti_5 + 1
    return boy if en_kisa <= boy <= en_uzun else None

def cift_tara(W, baslangic, okuma_boylari, F, R, en_kisa, en_uzun):
    """
        For one primer pair, returns which reads gave a product and the product lengths.
        Returns: (reads_with_product, [lengths])

    """
    kf, kr = len(F), len(R)
    yer = {("F", False): baglanma_yerleri(W, F, False),
           ("F", True):  baglanma_yerleri(W, F, True),
           ("R", False): baglanma_yerleri(W, R, False),
           ("R", True):  baglanma_yerleri(W, R, True)}
    if (yer[("F", False)].size + yer[("F", True)].size == 0): return 0, []
    if (yer[("R", False)].size + yer[("R", True)].size == 0): return 0, []

    def okumaya_dagit(pozlar):
        d = defaultdict(list)
        if pozlar.size == 0: return d
        okuma_no = np.searchsorted(baslangic, pozlar, side="right") - 1
        for o, p in zip(okuma_no.tolist(), pozlar.tolist()):
            d[o].append(p - int(baslangic[o]))
        return d

    dF0, dF1 = okumaya_dagit(yer[("F", False)]), okumaya_dagit(yer[("F", True)])
    dR0, dR1 = okumaya_dagit(yer[("R", False)]), okumaya_dagit(yer[("R", True)])
    sayi = 0; boylar = []
    for o in set(list(dF0) + list(dF1)):
        bulundu = None
        # F arti zincirde, R eksi zincirde
        for i in dF0.get(o, []):
            for j in dR1.get(o, []):
                b = urun_boyu(i, kf, j, kr, en_kisa, en_uzun)
                if b: bulundu = b; break
            if bulundu: break
        # R arti zincirde, F eksi zincirde
        if not bulundu:
            for i in dR0.get(o, []):
                for j in dF1.get(o, []):
                    b = urun_boyu(i, kr, j, kf, en_kisa, en_uzun)
                    if b: bulundu = b; break
                if bulundu: break
        if bulundu: sayi += 1; boylar.append(bulundu)
    return sayi, boylar

# ------------------------------------------------------------------ statistics and the decision
def wilson_alt(k, n, z=1.96):
    """Oranin guvenilen en dusuk degeri. Kucuk orneklemde ham orana guvenilmez."""
    if n <= 0: return 0.0
    p = k / n
    payda = 1 + z*z/n
    orta = p + z*z/(2*n)
    yayilim = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n))
    return max(0.0, (orta - yayilim) / payda)

def capraz_etiketi(rakip_alt, hedef_alt):
    """Tek bir alt sinirdan capraz etiketi. Iki paydayla ayri ayri cagrilir."""
    if rakip_alt <= 0.02: return "TEMIZ"
    if rakip_alt <= 0.10 and hedef_alt >= 10 * rakip_alt: return "SINIRDA"
    return "VAR"

def karar_ver(hedefler, sayimlar, toplamlar, erisilebilir, kapsayanlar=None):
    """
        Returns: dict(capraz, kapsama, hedef_alt, rakip_tx, rakip_alt, eksik, olculemeyen,
                    capraz_ham, capraz_kapsamali, rakip_alt_kapsamali, ayrilik)

        TWO DENOMINATORS, TWO DECISIONS
        capraz_ham        : the denominator is all the taxon's sampled reads
        capraz_kapsamali  : the denominator is the reads covering the locus (measured with
                            kapsayan_okuma)
        If the two diverge, the STRICTER one is taken and the divergence is returned.
        Switching silently to one of them would violate rule number 1. If no covering
        denominator is given only the raw one is computed, and that is visible from
        capraz_kapsamali being empty.

        COVERAGE IS THREE VALUED
        A target with NO reads at all is an unmeasured target; it is neither COMPLETE nor
        MISSING. In the earlier version such targets never entered the missing list, and
        once that list emptied the coverage was written COMPLETE, that is, full coverage
        was claimed for a target that had never been measured.

    """
    hset = set(hedefler)
    # THE DENOMINATOR OF THE TARGET RATIO MUST ALSO BE THE READS COVERING THE LOCUS (29 July).
    #
    # In the earlier version the denominator was ALL of the taxon's reads. That mistake
    # had already been fixed on the competitor side (the covering denominator) but not
    # on the target side; the two sides were measuring ASYMMETRICALLY.
    #
    # The concrete harm: the sample's read library is mixed, holding both short 16S
    # amplicons and full operon reads. Measured: of M. mazei's 102,006 reads only 9,257
    # (>=2000 bp) cover the operon. A pair sitting in the 16S-23S INTERGENIC REGION
    # gives a product in only 9% of all reads even when it works perfectly on its own
    # target, and is DISCARDED as "no product on its target". But a real qPCR runs on
    # genomic DNA, not on a short amplicon library.
    #
    # Methanosarcina species do not separate in 16S; the separation depends on exactly
    # that intergenic region. So this denominator mistake was structurally discarding
    # the project's most critical targets.
    #
    # When no covering denominator is given, the old behaviour continues (backward
    # compatible).
    def _hedef_oran(t):
        n = sayimlar.get(t, 0)
        ham = toplamlar.get(t, 0)
        if not kapsayanlar: return wilson_alt(n, ham)
        # The covering count cannot come out smaller than the count giving a product; if it
        # does, the measurement is inconsistent.
        kap = max(kapsayanlar.get(t, 0), n)
        return wilson_alt(n, kap)
    hedef_alt = min((_hedef_oran(t) for t in hedefler), default=0.0)
    hedef_alt_ham = min((wilson_alt(sayimlar.get(t,0), toplamlar.get(t,0)) for t in hedefler),
                        default=0.0)

    rakip_alt = 0.0; rakip_tx = ""
    rakip_alt_k = 0.0; rakip_tx_k = ""
    for tx, n in sayimlar.items():
        if tx in hset: continue
        a = wilson_alt(n, toplamlar.get(tx, 0))
        if a > rakip_alt: rakip_alt, rakip_tx = a, tx
        if kapsayanlar:
            # the denominator is the number of covering reads. The covering count cannot come
            # out smaller than the count giving a product; if it does, the measurement is
            # inconsistent and the safe side is taken.
            kn = max(kapsayanlar.get(tx, 0), n)
            ak = wilson_alt(n, kn)
            if ak > rakip_alt_k: rakip_alt_k, rakip_tx_k = ak, tx

    capraz_ham = capraz_etiketi(rakip_alt, hedef_alt)
    capraz_k = capraz_etiketi(rakip_alt_k, hedef_alt) if kapsayanlar else ""
    sira = {"TEMIZ": 0, "SINIRDA": 1, "VAR": 2}
    if capraz_k and capraz_k != capraz_ham:
        capraz = capraz_ham if sira[capraz_ham] > sira[capraz_k] else capraz_k
        ayrilik = (f"ham payda {capraz_ham} (rakip alt {rakip_alt:.4f}), "
                   f"kapsayan payda {capraz_k} (rakip alt {rakip_alt_k:.4f}), "
                   f"kati olan alindi: {capraz}")
    else:
        capraz = capraz_ham; ayrilik = ""

    # A target with no reads at all HAS NOT BEEN MEASURED; it is neither missing nor complete.
    olculemeyen = [t for t in hedefler if toplamlar.get(t, 0) == 0]
    eksik = []
    for t in hedefler:
        n = toplamlar.get(t, 0)
        if n == 0: continue
        oran = sayimlar.get(t, 0) / n
        ulasilir = erisilebilir.get(t, 0.0)
        # o takson baska bir ciftle yuksek oran veriyorsa okumalari saglamdir;
        # burada dusuk oran cikiyorsa sucu okumaya atamayiz
        if ulasilir >= 0.30 and oran < 0.25 * ulasilir: eksik.append(t)

    if olculemeyen: kapsama = "OLCULEMEDI"
    elif eksik: kapsama = "EKSIK"
    elif len(hedefler) <= 1: kapsama = "-"
    else: kapsama = "TAM"

    # THE NUMBER BEHIND A DECISION AND THE NAME OF THAT DECISION MUST COME FROM THE
    # SAME CALCULATION (the 29 July fix).
    #
    # In the earlier version rakip_tx came from the COVERING winner while rakip_alt came
    # from the RAW winner. When those were different taxa the sentence came out
    # nonsensical: "the best competitor is Y, the lower bound is 0.071" - the name
    # belonging to Y and the number to X.
    #
    # Worse: the design stage's elimination gate (in the source study, not imported
    # here) used `ra = okuma_karari["rakip_alt"]`, that is, the RAW number. A concrete
    # measurement: only 6 of a competitor's 1000 reads
    # cover the locus, and ALL SIX of those 6 give a product. The raw ratio 6/1000 gives
    # a lower bound of 0.0028; the covering ratio 6/6 gives 0.6097. The decision came out
    # correctly as "VAR", but because the elimination gate compared 0.0028 against a
    # threshold of 0.04, the candidate WAS NOT eliminated, and while the screen said
    # "threshold 0.02" it was printing a number BELOW the threshold. FULL amplification in
    # a competitor was passing as clean.
    #
    # Now the name and the number of a decision both come from THE denominator that led
    # to that decision.
    if capraz_k and capraz_k == capraz and rakip_tx_k:
        karar_tx, karar_alt, karar_payda = rakip_tx_k, rakip_alt_k, "kapsayan"
    else:
        karar_tx, karar_alt, karar_payda = rakip_tx, rakip_alt, "ham"

    return dict(capraz=capraz, kapsama=kapsama, hedef_alt=hedef_alt,
                hedef_alt_ham=hedef_alt_ham,
                rakip_tx=karar_tx,
                # rakip_alt is NOW the one CONSISTENT WITH THE DECISION. The raw value is in a separate field.
                rakip_alt=karar_alt,
                rakip_alt_ham=rakip_alt, rakip_alt_kapsamali=rakip_alt_k,
                rakip_tx_ham=rakip_tx, rakip_tx_kapsamali=rakip_tx_k,
                karar_paydasi=karar_payda,
                eksik=eksik, olculemeyen=olculemeyen,
                capraz_ham=capraz_ham, capraz_kapsamali=capraz_k, ayrilik=ayrilik)

# ------------------------------------------------------------------ birim testleri
def birim_testleri(yazdir=True):
    """
        A known-answer test of every deciding function. No external tool is needed.
        Most of the tests push the "does it correctly find a cross reaction" side, because
        that is the expensive direction to be wrong in: showing a bad primer as clean.

    """
    hata = [0]
    def K(ad, bulunan, beklenen):
        ok = bulunan == beklenen
        if not ok: hata[0] += 1
        if yazdir: print(f"  {u'PASS' if ok else u'FAIL'}  {ad:<50} {bulunan} / {beklenen}")

    def tara(kalip, primer, ters=False):
        b, bas = kalip_kur([kalip], 30)
        W = np.lib.stride_tricks.sliding_window_view(b, 30)
        return [int(x) for x in baglanma_yerleri(W, primer, ters)]
    def urun(kalip, F, R, en_kisa=60, en_uzun=400):
        b, bas = kalip_kur([kalip], 30)
        W = np.lib.stride_tricks.sliding_window_view(b, 30)
        n, boy = cift_tara(W, bas, [len(kalip)], F, R, en_kisa, en_uzun)
        return (boy[0] if boy else None)

    rnd = random.Random(11)
    kalip = "".join(rnd.choice("ACGT") for _ in range(1200))
    F = kalip[300:320]; R = rc(kalip[430:450])
    D = {"A":"C","C":"A","G":"T","T":"G"}
    def yaz(d, p, y): return d[:p] + y + d[p+1:]
    def boz(d, *pozlar):
        for p in pozlar: d = yaz(d, p, D[d[p]])
        return d

    if yazdir: print(u'  --- the 3\' end rule (the last two bases exact, at most one mismatch in the last five)')
    K("a mismatch at the last base does not bind",        tara(boz(kalip,319), F), [])
    K("a mismatch at the second base from the end does not bind",  tara(boz(kalip,318), F), [])
    K("a mismatch AT THE THIRD BASE FROM THE END DOES bind",   tara(boz(kalip,317), F), [300])
    K("a mismatch at the fourth base from the end binds", tara(boz(kalip,316), F), [300])
    K("two mismatches in the last five bases do not bind",tara(boz(kalip,317,316), F), [])
    if yazdir: print(u'  --- the mismatch budget')
    K("an exact match binds",                  tara(kalip, F), [300])
    K("three inner mismatches bind",             tara(boz(kalip,302,305,308), F), [300])
    K("four inner mismatches do not bind",          tara(boz(kalip,302,305,308,311), F), [])
    K("two mismatches at the 5' end bind",       tara(boz(kalip,300,301), F), [300])
    if yazdir: print(u'  --- strand orientation')
    K("the forward primer is found on the plus strand",    tara(kalip, F, False), [300])
    K("the forward primer is not found on the minus strand",   tara(kalip, F, True), [])
    K("the reverse primer is found on the minus strand",     tara(kalip, R, True), [430])
    K("the reverse primer is not found on the plus strand",    tara(kalip, R, False), [])
    if yazdir: print(u'  --- product geometry')
    K("a product at the right orientation and distance",            urun(kalip, F, R), 150)
    K("a product on a flipped template as well",        urun(rc(kalip), F, R), 150)
    K("the primer order does not change the result",      urun(kalip, R, F), 150)
    K("no product on an unrelated template",
      urun("".join(rnd.choice("ACGT") for _ in range(1200)), F, R), None)
    K("no product once the 3' end is spoiled",              urun(boz(kalip,319), F, R), None)
    uzak = kalip[:300] + F + "".join(rnd.choice("ACGT") for _ in range(900)) + rc(R) + kalip[900:]
    K("primers too far apart give no product",        urun(uzak, F, R), None)
    ters_yon = (kalip[:300] + rc(R) + kalip[320:400] + F + kalip[420:])
    K("primers not facing one another give no product", urun(ters_yon, F, R), None)
    if yazdir: print(u'  --- template construction (the reads must not run into one another)')
    b, bas = kalip_kur([kalip[:400], kalip[400:800]], 30)
    W = np.lib.stride_tricks.sliding_window_view(b, 30)
    n, _ = cift_tara(W, bas, [400, 400], F, R, 60, 400)
    K("a region split across two reads gives no product", n, 0)
    n2, _ = cift_tara(W, bas, [400, 400], kalip[300:320], rc(kalip[360:380]), 40, 200)
    K("a region inside a single read gives a product",   n2, 1)
    if yazdir: print(u'  --- statistics')
    K("wilson penalises a small sample",   round(wilson_alt(1, 17), 3), 0.01)
    K("wilson approaches the ratio on a large sample",round(wilson_alt(150, 300), 2), 0.44)
    K("a zero count gives a zero lower bound",           wilson_alt(0, 300), 0.0)

    if yazdir: print(u'  --- the number of covering reads (the denominator of the cross reaction ratio)')
    kl = kalip[:400]
    b3, bas3 = kalip_kur([kl, kl, "".join(rnd.choice("ACGT") for _ in range(400))], 30)
    W3 = np.lib.stride_tricks.sliding_window_view(b3, 30)
    K("the reads the primer bound are counted",
      kapsayan_okuma(W3, bas3, kl[300:320], rc(kl[360:380])), 2)
    K("zero when no primer binds at all",
      kapsayan_okuma(W3, bas3, "ACGTACGTACGTACGTACGT", "TTTTTTTTTTTTTTTTTTTT"), 0)

    if yazdir: print(u'  --- the decision (this function had never been tested before)')
    # 1. A target with NO reads at all. The earlier version did not put this on the
    #    missing list, and once the list emptied the coverage was written as COMPLETE:
    #    a claim of full coverage for a target that was never measured.
    k1 = karar_ver(["83986", "118126", "394967", "2201"], {"2223": 1}, {"2223": 300}, {})
    K("coverage is not written as FULL for a target with no reads", k1["kapsama"], "OLCULEMEDI")
    K("the targets that cannot be measured are listed", len(k1["olculemeyen"]), 4)

    # 2. The competitor covers the locus very little but gives a product in every read
    #    it does cover. CLEAN under the raw denominator, VAR under the covering one.
    #    The stricter one has to be taken.
    hed = {"818": 150}; top = {"818": 300, "28116": 300}
    say = dict(hed); say["28116"] = 5
    k2 = karar_ver(["818"], say, top, {}, {"28116": 6, "818": 300})
    K("the raw denominator alone would have said TEMIZ", k2["capraz_ham"], "TEMIZ")
    K("the covering denominator makes the cross reaction visible", k2["capraz_kapsamali"], "VAR")
    K("the stricter decision is the one taken", k2["capraz"], "VAR")
    K("the disagreement is returned to be written to the log", bool(k2["ayrilik"]), True)
    # 29 JULY: THE NUMBER BEHIND A DECISION MUST COME FROM THAT DECISION'S DENOMINATOR.
    # In the earlier version rakip_alt came from the RAW denominator, and because the
    # design stage's elimination gate compared that number against a threshold of 0.04, a
    # candidate showing FULL amplification in a competitor WAS NOT ELIMINATED. This item
    # locks exactly that gate.
    K("the decision came from the covering denominator", k2["karar_paydasi"], "kapsayan")
    K("rakip_alt agrees with the decision (not the raw one)", k2["rakip_alt"] == k2["rakip_alt_kapsamali"], True)
    K("rakip_alt is above the elimination threshold", k2["rakip_alt"] >= 0.04, True)
    K("the raw value is kept in a separate field", k2["rakip_alt_ham"] < 0.04, True)
    # The name and the number must come from the same taxon
    K("the competitor name and its count come from the same source", k2["rakip_tx"], k2["rakip_tx_kapsamali"])

    # 4. THE TARGET'S DENOMINATOR MUST ALSO BE COVERING READS (29 July).
    # A pair sitting in the intergenic region: only 90 of a target's 1000 reads cover
    # the locus, and 88 of those 90 give a product. The raw denominator gives 88/1000 ->
    # 0.07 (eliminated), the covering denominator 88/90 -> 0.92 (perfect). Because the
    # short 16S amplicons inflated the denominator, intergenic targets such as M. mazei
    # and M. barkeri were being eliminated structurally.
    k4 = karar_ver(["818"], {"818": 88}, {"818": 1000}, {}, {"818": 90})
    K("the target ratio is computed with the covering denominator", k4["hedef_alt"] > 0.80, True)
    K("the raw denominator is kept in a separate field", k4["hedef_alt_ham"] < 0.12, True)
    # Kapsayan payda verilmezse eski davranis surer
    k5 = karar_ver(["818"], {"818": 88}, {"818": 1000}, {})
    K("with no covering denominator the raw behaviour continues", abs(k5["hedef_alt"] - k5["hedef_alt_ham"]) < 1e-9, True)
    # An inconsistent measurement: if the covering count is smaller than the count
    # giving a product, the safe side is taken
    k6 = karar_ver(["818"], {"818": 88}, {"818": 1000}, {}, {"818": 10})
    K("when the covering count is smaller than the product the lower bound stays under 1", k6["hedef_alt"] <= 1.0, True)

    # 3. The same numbers, but the competitor covers the locus widely. Both are CLEAN,
    # there is no divergence.
    k3 = karar_ver(["818"], say, top, {}, {"28116": 250, "818": 300})
    K("on broad coverage both denominators say TEMIZ", (k3["capraz_ham"], k3["capraz_kapsamali"]),
      ("TEMIZ", "TEMIZ"))
    K("empty when there is no disagreement", k3["ayrilik"], "")

    # 4. If the covering count comes out smaller than the count giving a product the
    #    measurement is inconsistent; the denominator is pulled to the count giving a
    #    product, so the ratio does not exceed 1.
    k4 = karar_ver(["818"], say, top, {}, {"28116": 2, "818": 300})
    K("on an inconsistent coverage ratio the lower bound stays under 1", k4["rakip_alt_kapsamali"] <= 1.0, True)
    return hata[0]

# ------------------------------------------------------------------ blastn ikinci gorus
def blast_ikinci_gorus(veri, ciftler, cikti, threads, en_kisa, en_uzun):
    """
        Asks the same question with blastn as well. Because blastn TRUNCATES an alignment
        by its score, it does not report bindings with mismatches near the 3' end; so this
        count generally comes out LOWER than the brute force scan and is not used for the
        decision. The point is only to see that two independent methods point the same way.

    """
    if not (arac_var("blastn") and arac_var("makeblastdb")):
        log(u'blastn was not found, the second opinion was skipped'); return {}
    fa = os.path.join(cikti, "numune_okumalari.fasta")
    with open(fa, "w") as fh:
        for tx, okumalar in veri.items():
            for i, s in enumerate(okumalar): fh.write(f">{tx}|{i}\n{s}\n")
    db = os.path.join(cikti, "numune_db")
    r = subprocess.run(["makeblastdb","-in",fa,"-dbtype","nucl","-out",db],
                       capture_output=True, text=True)
    if r.returncode != 0: log(u'makeblastdb failed, the second opinion was skipped'); return {}
    pf = os.path.join(cikti, "primerler.fasta")
    with open(pf, "w") as fh:
        for taban, d in sorted(ciftler.items()):
            fh.write(f">{taban}_F\n{d['F']}\n>{taban}_R\n{d['R']}\n")
    r = subprocess.run(["blastn","-task","blastn-short","-query",pf,"-db",db,
        "-outfmt","6 qseqid sseqid qstart qend sstart send sstrand qlen",
        "-word_size","7","-evalue","1000","-dust","no","-soft_masking","false",
        "-max_target_seqs","100000","-num_threads",str(threads)],
        capture_output=True, text=True)
    if r.returncode != 0: log(u'blastn failed, the second opinion was skipped'); return {}
    hit = defaultdict(lambda: defaultdict(lambda: {"F": [], "R": []}))
    for satir in r.stdout.splitlines():
        p = satir.split("\t")
        if len(p) < 8: continue
        q, s, qs, qe, ss, se, strand, qlen = p[:8]
        try: qs, qe, ss, se, qlen = int(qs), int(qe), int(ss), int(se), int(qlen)
        except ValueError: continue
        if "_" not in q: continue
        taban, rol = q.rsplit("_", 1); rol = rol.upper()
        if rol not in ("F","R"): continue
        if qe != qlen: continue                      # BLAST truncated it; there is no 3' end
        hit[taban][s][rol].append((ss, se, strand, qlen))
    sonuc = defaultdict(lambda: defaultdict(int))
    for taban, okumalar in hit.items():
        for okuma_id, d in okumalar.items():
            tx = okuma_id.split("|")[0]
            var = False
            for fs, fe, fst, fl in d["F"]:
                for rs, re_, rst, rl in d["R"]:
                    if fst == rst: continue
                    arti = (fs, fe, fl) if fst == "plus" else (rs, re_, rl)
                    eksi = (rs, re_, rl) if fst == "plus" else (fs, fe, fl)
                    if arti[1] > eksi[1]: continue
                    boy = eksi[0] - arti[0] + 1
                    if en_kisa <= boy <= en_uzun: var = True; break
                if var: break
            if var: sonuc[taban][tx] += 1
    return sonuc

# ------------------------------------------------------------------ selftest
def selftest():
    print("=" * 72); print(u'THE KNOWN ANSWER EXAM OF THE RULES'); print("=" * 72)
    h = birim_testleri()
    print("=" * 72)
    print(u'SELFTEST PASSED' if h == 0 else f"SELFTEST FAILED, {h} tests failed")
    print("=" * 72)
    return 0 if h == 0 else 1

# ------------------------------------------------------------------ ana
def main():
    global LOGF
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", "--kok", dest="kok"); ap.add_argument("--primerler", nargs="*", default=[])
    ap.add_argument("--out", default=""); ap.add_argument("--okuma", type=int, default=300)
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--en-kisa", type=int, default=60)
    ap.add_argument("--en-uzun", type=int, default=400)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--blast", action="store_true", help="second opinion via blastn (does not decide)")
    a = ap.parse_args()
    if a.selftest: sys.exit(selftest())
    if not a.kok: ap.error("--root is required")

    zaman = time.strftime("%Y%m%d_%H%M%S")
    cikti = a.out or os.path.join(a.kok, "tools", "0_TESLIM_RAPOR", "ICPCR_" + zaman)
    os.makedirs(cikti, exist_ok=True)
    LOGF = open(os.path.join(cikti, "log.txt"), "w", encoding="utf-8")
    log(f"output: {cikti}")

    log(u'running the known answer exam of the rules')
    if birim_testleri(yazdir=False) != 0:
        log(u'THE EXAM FAILED, stopped. For the detail: --selftest'); sys.exit(2)
    log(u'the exam passed')

    primer_yollari = [p if os.path.isabs(p) else os.path.join(a.kok, p) for p in a.primerler]
    if not primer_yollari:
        primer_yollari = [os.path.join(a.kok, "tools", "0_TESLIM_RAPOR", f) for f in
                          ("OLIGO_SIPARIS_GUNCEL_10ASSAY.csv", "OLIGO_SIPARIS_OPSIYONEL.csv",
                           "GRUP_VE_UNIVERSAL_PRIMERLER.csv")]
    ciftler = primerleri_oku(primer_yollari)
    log(f"{len(ciftler)} primer pairs were read")
    if not ciftler: sys.exit(u'no primer was found')

    log(u'loading the reads')
    veri = okumalari_yukle(a.kok, a.okuma)
    toplam_okuma = sum(len(v) for v in veri.values())
    toplam_baz = sum(len(s) for v in veri.values() for s in v)
    log(f"{len(veri)} taxa, {toplam_okuma} reads, {toplam_baz/1e6:.1f} million bases")

    en_uzun_primer = max(max(len(d["F"]), len(d["R"])) for d in ciftler.values())
    log(f"preparing the brute force scan (window {en_uzun_primer} bases)")
    kalipler = {}
    for tx, okumalar in veri.items():
        b, bas = kalip_kur(okumalar, en_uzun_primer)
        W = np.lib.stride_tricks.sliding_window_view(b, en_uzun_primer)
        kalipler[tx] = (W, bas, [len(s) for s in okumalar])

    log(u'the scan is starting')
    t0 = time.time()
    sayimlar = {}; boy_bilgisi = {}; kapsayanlar = {}
    toplamlar = {tx: len(v) for tx, v in veri.items()}
    # Taxa that appear in the target lists but have no reads in the sample. No
    # measurement can be made for them; "not measured" and "came out clean" must not
    # arrive at the same conclusion.
    hedefte_gecen = {t for hs in HEDEF_TAXID.values() for t in hs}
    okumasiz = sorted(hedefte_gecen - set(toplamlar))
    if okumasiz:
        log("")
        log(f"WARNING: {len(okumasiz)} target taxa have NO reads at all in the sample, "
            f"NO measurement is possible for them:")
        for t in okumasiz: log(f"    {isim(t)} (taxid {t})")
    for i, taban in enumerate(sorted(ciftler), 1):
        F, R = ciftler[taban]["F"], ciftler[taban]["R"]
        s = {}; b = {}
        for tx, (W, bas, boylar) in kalipler.items():
            n, urunler = cift_tara(W, bas, boylar, F, R, a.en_kisa, a.en_uzun)
            if n: s[tx] = n; b[tx] = int(statistics.median(urunler))
        # The number of reads covering the locus, for the denominator of the cross reaction
        # ratio. It is computed only for the taxa and targets that give a product;
        # elsewhere the numerator is zero anyway, so the denominator does not change the
        # decision.
        ilgili = set(s) | set(HEDEF_TAXID.get(taban, []))
        kap = {}
        for tx in ilgili:
            if tx in kalipler:
                W, bas, _ = kalipler[tx]
                kap[tx] = kapsayan_okuma(W, bas, F, R)
        kapsayanlar[taban] = kap
        sayimlar[taban] = s; boy_bilgisi[taban] = b
        log(f"  [{i}/{len(ciftler)}] {taban}  {sum(s.values())} reads giving a product")
    log(f"the scan finished ({int(time.time()-t0)} s)")

    # her taksonun ulasilabilirligi: herhangi bir ciftle elde edilen en yuksek oran
    erisilebilir = {}
    for tx, n in toplamlar.items():
        erisilebilir[tx] = max((sayimlar[c].get(tx, 0) / max(1, n) for c in sayimlar), default=0.0)

    blast = {}
    if a.blast:
        log(u'taking a second opinion with blastn')
        blast = blast_ikinci_gorus(veri, ciftler, cikti, a.threads, a.en_kisa, a.en_uzun)

    satirlar = []; ozet = []
    log("")
    log("SONUCLAR")
    for taban in sorted(ciftler):
        s = sayimlar[taban]; b = boy_bilgisi[taban]
        hedefler = HEDEF_TAXID.get(taban, [])
        K = karar_ver(hedefler, s, toplamlar, erisilebilir, kapsayanlar.get(taban))
        capraz, kapsama, hedef_alt = K["capraz"], K["kapsama"], K["hedef_alt"]
        rakip_tx, rakip_alt, eksik = K["rakip_tx"], K["rakip_alt"], K["eksik"]
        etiket = f"capraz {capraz}" + (f", kapsama {kapsama}" if kapsama != "-" else "")
        log("")
        log(f"  {taban}   [{etiket}]")
        if K["ayrilik"]:
            log(f"    DISAGREEMENT: {K['ayrilik']}")
        if K["olculemeyen"]:
            log(u'    TARGETS THAT CANNOT BE MEASURED (no reads in the sample, no coverage can be claimed): '
                + ", ".join(isim(t) for t in K["olculemeyen"]))
        sebep = {}
        for t in eksik:
            W, bas, boylar = kalipler[t]
            of = tek_primer_orani(W, bas, len(boylar), ciftler[taban]["F"])
            orr = tek_primer_orani(W, bas, len(boylar), ciftler[taban]["R"])
            if max(of, orr) < 0.25:
                sebep[t] = (f"bu bolge okumalarda yok (ileri {of:.2f}, geri {orr:.2f}), "
                            f"primerin sucu degil")
            else:
                zayif = "geri" if orr < of else "ileri"
                sebep[t] = (f"{zayif} primer bu organizmaya UYMUYOR "
                            f"(ileri {of:.2f}, geri {orr:.2f})")
        for t in hedefler:
            n = toplamlar.get(t, 0); k = s.get(t, 0)
            im = f"  <<< {sebep[t]}" if t in eksik else ""
            log(f"    TARGET  {isim(t):<34} {k:>4}/{n:<4} ratio {k/max(1,n):.3f}"
                f"  product {b.get(t,0)} bp{im}")
        rakipler = sorted(((tx, k) for tx, k in s.items() if tx not in set(hedefler)),
                          key=lambda x: -x[1] / max(1, toplamlar.get(x[0], 1)))
        if not rakipler:
            log(u'    no product formed in any of the competitor taxa')
        kap = kapsayanlar.get(taban, {})
        for tx, k in rakipler[:6]:
            n = toplamlar.get(tx, 0); kn = max(kap.get(tx, 0), k)
            # The two ratios are written side by side. The gap between them shows how much of
            # the locus the competitor covers; a large gap is the difference between "a product
            # in few reads" and "a product in every covering read".
            log(f"    competitor  {isim(tx):<34} {k:>4}/{n:<4} ratio {k/max(1,n):.3f}"
                f"  |  the ratio within {kn:>4} covering {k/max(1,kn):.3f}"
                f"  product {b.get(tx,0)} bp")
        for tx in sorted(set(list(s) + hedefler)):
            n = toplamlar.get(tx, 0); k = s.get(tx, 0); kn = max(kap.get(tx, 0), k)
            satirlar.append(dict(cift=taban, taxid=tx, organizma=isim(tx),
                hedef_mi=("hedef" if tx in set(hedefler) else "rakip"),
                urun_veren=k, toplam_okuma=n, oran=round(k/max(1,n), 4),
                oran_alt_sinir=round(wilson_alt(k, n), 4),
                kapsayan_okuma=kn,
                kapsamali_oran=(round(k/kn, 4) if kn else ""),
                kapsamali_alt_sinir=(round(wilson_alt(k, kn), 4) if kn else ""),
                medyan_urun_bp=b.get(tx, 0),
                blastn_ikinci_gorus=(blast.get(taban, {}).get(tx, "") if blast else "")))
        ozet.append(dict(cift=taban, capraz=capraz, kapsama=kapsama,
            hedef_alt_sinir=round(hedef_alt, 4),
            en_iyi_rakip=isim(rakip_tx) if rakip_tx else "",
            en_iyi_rakip_alt_sinir=round(rakip_alt, 4),
            capraz_ham_payda=K["capraz_ham"],
            capraz_kapsayan_payda=K["capraz_kapsamali"],
            rakip_alt_kapsamali=round(K["rakip_alt_kapsamali"], 4),
            olcum_ayriligi=K["ayrilik"],
            olculemeyen_hedefler="; ".join(isim(t) for t in K["olculemeyen"]),
            urun_vermeyen_hedefler="; ".join(isim(t) for t in eksik),
            sebep="; ".join(f"{isim(t)}: {sebep[t]}" for t in eksik)))

    yol = os.path.join(cikti, "icpcr_sonuc.csv")
    with open(yol, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=["cift","taxid","organizma","hedef_mi","urun_veren",
            "toplam_okuma","oran","oran_alt_sinir","kapsayan_okuma","kapsamali_oran",
            "kapsamali_alt_sinir","medyan_urun_bp","blastn_ikinci_gorus"])
        w.writeheader(); w.writerows(satirlar)
    yol2 = os.path.join(cikti, "icpcr_ozet.csv")
    with open(yol2, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=["cift","capraz","kapsama","hedef_alt_sinir",
            "en_iyi_rakip","en_iyi_rakip_alt_sinir","capraz_ham_payda",
            "capraz_kapsayan_payda","rakip_alt_kapsamali","olcum_ayriligi",
            "olculemeyen_hedefler","urun_vermeyen_hedefler","sebep"])
        w.writeheader(); w.writerows(ozet)

    sira = {"TEMIZ": 0, "SINIRDA": 1, "VAR": 2}
    log("")
    log(u'THE COLLECTED VERDICT')
    log(f"  {'capraz':<9}{'kapsama':<9}{'cift':<46}{'hedef':>7}{'rakip':>8}")
    for satir in sorted(ozet, key=lambda x: (sira.get(x["capraz"], 3), x["kapsama"] == "EKSIK")):
        log(f"  {satir['capraz']:<9}{satir['kapsama']:<9}{satir['cift']:<46}"
            f"{satir['hedef_alt_sinir']:>7.3f}{satir['en_iyi_rakip_alt_sinir']:>8.3f}"
            f"  {satir['en_iyi_rakip']}")
    eksikli = [o for o in ozet if o["urun_vermeyen_hedefler"]]
    if eksikli:
        log("")
        log(u'THE ORGANISMS IT TARGETS BUT GIVES NO PRODUCT IN')
        for o in eksikli:
            log(f"  {o['cift']}")
            log(f"     {o['sebep']}")
    log("")
    log(u'The ratios are raw counts; the decision is made on the Wilson lower bound, so')
    log(u'that in a taxon with seventeen reads a single read cannot look like six percent and spoil the decision.')
    log("")
    log(f"written: {yol}")
    log(f"written: {yol2}")
    LOGF.close()

if __name__ == "__main__":
    main()
