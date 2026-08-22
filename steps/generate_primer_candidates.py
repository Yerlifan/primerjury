#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_primer_candidates.py
It produces primer pair candidates from a target consensus by applying the oligo,
thermodynamics, product and region rules of the meeting decisions.

The design principles:
  * Every Tm is measured with two independent libraries (primer3 and Biopython).
    The SYSTEMATIC OFFSET between the two libraries is measured from the data (the
    median of the difference across all candidates), and then the oligos departing
    from that offset by more than the tolerance are eliminated. The offset is not
    written by hand.
  * No threshold is fixed without looking at the data; every one of them can be
    changed from the command line and the value used is written to the log.
  * The product is verified by machine: the head of the fragment cut from the
    template must equal the forward primer exactly and its tail the reverse
    complement of the reverse primer. If they are not equal the candidate is not
    eliminated silently, a counter goes up.
  * No position in the mask file may enter a primer footprint. Masked positions
    inside the product are free, because the rule constrains the primer placement
    only.

WHICH DIRECTORY THE CONSENSUS MUST BE READ FROM (the 2026-08-21 fix)
  Use `konsensus_kanonik/` only. The old examples pointed at the
  `consensus sequences/` directory; that directory is MIXED ORIENTATION (measured:
  71 antisense / 27 sense). On a reversed consensus, in-silico PCR SILENTLY gives 0
  products; the measured loss is 100 percent and the evidence is
  `screening/orientation_impact_test.py`. The prohibition is written in
  `screening/config.py` as well (KONSENSUS_KANONIK).
  The canonical directory is produced with `screening/build_canonical.py`. Which bin
  corresponds to which file is in `konsensus_kanonik/INDEKS.tsv`; do not guess the
  file name, read it from the index.

Usage:
  python3 generate_primer_candidates.py       --consensus "konsensus_kanonik/A1-1_2209.kanonik.fa"       --mask      "N_analizi/maske/A1-1_2209_maske.bed"       --out       "primer_adaylari/A1-1_2209.tsv"

"""
import argparse, csv, itertools, os, statistics, sys

try:
    import primer3
except ImportError:
    sys.exit(u'primer3-py is not installed:  pip install primer3-py')
try:
    from Bio.Seq import Seq
    from Bio.SeqUtils import MeltingTemp as mt
except ImportError:
    sys.exit(u'biopython is not installed:  pip install biopython')

COMP = str.maketrans("ACGTNRYSWKMBDHV", "TGCANYRSWMKVHDB")


def rc(s):
    return s.translate(COMP)[::-1]


# ----------------------------------------------------------------- arguments
def get_args():
    p = argparse.ArgumentParser(description='Produce primer pair candidates')
    p.add_argument("--consensus", required=True)
    p.add_argument("--mask", default=None, help='the BED produced by the '
                                                'ambiguous base analysis')
    p.add_argument("--mask-contig", default=None,
                   help='only BED rows whose first column has this name are '
                        'used; without it every row is taken')
    p.add_argument("--ambig", default=None, help="IUPAC kodlu consensus (bilgi amacli)")
    p.add_argument("--out", required=True)
    p.add_argument("--label", default=None)
    # oligo kurallari
    p.add_argument("--len-min", type=int, default=18)
    p.add_argument("--len-max", type=int, default=25)
    p.add_argument("--gc-min", type=float, default=40.0)
    p.add_argument("--gc-max", type=float, default=60.0)
    p.add_argument("--gc-hard-min", type=float, default=35.0)
    p.add_argument("--gc-hard-max", type=float, default=65.0)
    p.add_argument("--gc-clamp-last", type=int, default=5,
                   help="son kac bazda GC sayilir")
    p.add_argument("--gc-clamp-max", type=int, default=3,
                   help="maximum G or C in the last N bases")
    p.add_argument("--homopolymer-max", type=int, default=4)
    p.add_argument("--require-3p-gc", type=int, default=1,
                   help="1 ise 3' uc G or C with bitmeli")
    p.add_argument("--degeneracy-budget", type=int, default=0,
                   help='NO LONGER HAS ANY EFFECT. Oligos are produced as '
                        'ACGT only and template ambiguity is handled with '
                        '--iupac-max. The flag stays for backward '
                        'compatibility and prints a warning when it is given.')
    p.add_argument("--degeneracy-fold-max", type=int, default=4,
                   help='NO LONGER HAS ANY EFFECT, it stays for backward '
                        'compatibility')
    p.add_argument("--iupac-max", type=int, default=2,
                   help='the number of IUPAC positions allowed in the '
                        'template window; they are resolved to concrete '
                        'bases, so no degenerate base enters the oligo. With '
                        '0, a window holding an IUPAC code is never used.')
    p.add_argument("--iupac-clamp-forbidden", type=int, default=5,
                   help="IUPAC is not accepted in the last this many bases of the oligo")
    # termodinamik
    p.add_argument("--tm-min", type=float, default=58.0)
    p.add_argument("--tm-max", type=float, default=62.0)
    p.add_argument("--tm-hard-min", type=float, default=57.0)
    p.add_argument("--tm-hard-max", type=float, default=63.0)
    p.add_argument("--tm-cross-k", type=float, default=4.0,
                   help="tolerans = k carpi kaymanin standart sapmasi; veriden "
                        "turetilir")
    p.add_argument("--tm-cross-tol", type=float, default=None,
                   help="pins the tolerance manually; k is ignored when given")
    p.add_argument("--hairpin-dg-min", type=float, default=-3000.0)
    p.add_argument("--homodimer-dg-min", type=float, default=-6000.0)
    p.add_argument("--heterodimer-dg-min", type=float, default=-6000.0)
    # the buffer conditions, given identically to both libraries
    p.add_argument("--mv", type=float, default=50.0, help='the monovalent cation '
                                                          'concentration in mM')
    p.add_argument("--dv", type=float, default=1.5, help='the divalent cation '
                                                         'concentration in mM')
    p.add_argument("--dntp", type=float, default=0.6, help="dNTP mM")
    p.add_argument("--dna-conc", type=float, default=50.0, help="oligo nM")
    # urun
    p.add_argument("--prod-min", type=int, default=70)
    p.add_argument("--prod-max", type=int, default=250)
    p.add_argument("--prod-hard-max", type=int, default=300)
    p.add_argument("--prod-best-min", type=int, default=90)
    p.add_argument("--prod-best-max", type=int, default=150)
    p.add_argument("--prod-gc-min", type=float, default=40.0)
    p.add_argument("--prod-gc-max", type=float, default=60.0)
    p.add_argument("--pair-tm-diff-max", type=float, default=1.5)
    p.add_argument("--max-pairs", type=int, default=5000,
                   help="maximum output rows, trimmed by score")
    p.add_argument("--min-locus-spacing", type=int, default=0,
                   help='thins out shifted copies of the same locus: when two '
                        'candidate pairs start within this distance at both '
                        'the forward and the reverse end, the one with the '
                        'worse score is dropped. 0 turns it off.')
    return p.parse_args()


# ------------------------------------------------------------------ yardimcilar
def read_fasta(path):
    """Reads a single record FASTA. It stops on a file with several records, because
    merging the records silently produces an artificial chimeric primer across the
    junction.

    """
    seq, names = [], []
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith(">"):
            names.append(line[1:].strip())
        else:
            seq.append(line.strip())
    if len(names) > 1:
        sys.exit(u'ERROR: %s holds %d records. This script expects a single record consensus; merging the records produces an artificial primer across the junction.' % (path, len(names)))
    return (names[0] if names else None), "".join(seq).upper()


def read_mask(path, seqlen, contig=None, strict_missing=True):
    """BED -> the forbidden position set (0 based), a class counter and a contig report.
    If contig is given, only the rows belonging to that contig are taken; if no row
    matches, a warning is returned so that this does not turn into a silent absence of
    masking.

    """
    bad, classes, seen = set(), {}, {}
    if not path:
        return bad, classes, seen
    if not os.path.exists(path):
        if strict_missing:
            sys.exit(u'ERROR: there is no mask file: %s. Once a path is given the file must exist, otherwise masking is switched off silently.' % path)
        return bad, classes, seen
    for line in open(path, encoding="utf-8", errors="replace"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        try:
            st, e = int(f[1]), int(f[2])
        except ValueError:
            continue
        name = f[0]
        seen[name] = seen.get(name, 0) + 1
        if contig is not None and name != contig:
            continue
        cls = f[3] if len(f) > 3 else "maske"
        for i in range(max(0, st), min(seqlen, e)):
            bad.add(i)
        classes[cls] = classes.get(cls, 0) + (e - st)
    return bad, classes, seen


def gc_pct(s):
    return 100.0 * (s.count("G") + s.count("C")) / len(s) if s else 0.0


def max_run(s):
    best = run = 1
    for i in range(1, len(s)):
        run = run + 1 if s[i] == s[i - 1] else 1
        best = max(best, run)
    return best if s else 0


# N is deliberately absent from the list: an N in the consensus is a coverage gap,
# not a deliberate degeneracy. No oligo holding an N is produced; primer3 raises
# ValueError on an N too, so the risk of a silently wrong Tm disappears as well.
DEGEN_FOLD = {"R": 2, "Y": 2, "S": 2, "W": 2, "K": 2, "M": 2,
              "B": 3, "D": 3, "H": 3, "V": 3}

# The concrete base an IUPAC code is resolved to. Because the consensus is produced
# with -A, the variable positions are written as R, Y, S, W, K or M. The meeting
# decision wants no degenerate base in an oligo; so the code is resolved to a
# concrete base and whether the choice was right is tested experimentally on the raw
# reads in specificity.py. The resolution rule is fixed and deterministic: the
# alphabetically first base of the set.
IUPAC_COZ = {"R": "A", "Y": "C", "S": "C", "W": "A", "K": "G", "M": "A",
             "B": "C", "D": "A", "H": "A", "V": "A"}

# All the bases the code stands for. Reducing an ambiguous position to a single
# base would be throwing information away; instead every alternative is produced,
# all of them pass the same filters, and the raw read scan in specificity.py decides
# which one really binds.
IUPAC_KUME = {"R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT",
              "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG"}


def _iupac_denetle(win, a, uc="her"):
    """Does the window obey the IUPAC rules. (ambiguous_positions, reason).

    uc: which strand's oligo will be produced.
        "F"   -> oligo = win        , the 3' end is the END of the window
        "R"   -> oligo = rc(win)    , the 3' end is the START of the window
        "her" -> forbid both ends (backward compatible, stricter)

    The strand distinction is required: both F and R are produced from the same window
    and R's 3' end falls on the START of the window. A one sided check let the outermost
    base of the R primer be fixed silently to a single allele from an ambiguous
    position. Measured: on real data, in 435 of 2000 rows R's 3' terminal base came from
    a resolved IUPAC position.

    """
    if "N" in win:
        return None, "kalipta_N"
    k = [i for i, c in enumerate(win) if c not in "ACGT"]
    if not k:
        return [], None
    if any(c not in IUPAC_KUME for c in win if c not in "ACGT"):
        return None, "tanimsiz_kod"
    if len(k) > a.iupac_max:
        return None, "iupac_fazla"
    n = a.iupac_clamp_forbidden
    if uc in ("F", "her") and any(i >= len(win) - n for i in k):
        return None, "iupac_3p"
    if uc in ("R", "her") and any(i < n for i in k):
        return None, "iupac_3p"
    return k, None


def iupac_coz(win, a, uc="her"):
    """Geriye donuk uyum: tek, belirlenimci cozum dondurur."""
    k, why = _iupac_denetle(win, a, uc)
    if k is None:
        return None, why
    if not k:
        return win, 0
    return "".join(IUPAC_COZ.get(c, c) for c in win), len(k)


def iupac_varyantlar(win, a, uc="her"):
    """Produces every concrete ACGT counterpart of the window.
    ([(oligo, resolved_position_count), ...], None) or (None, reason).
    uc: the strand of the oligo to be produced; it is passed straight to _iupac_denetle.

    """
    k, why = _iupac_denetle(win, a, uc)
    if k is None:
        return None, why
    if not k:
        return [(win, 0)], None
    secenekler = [IUPAC_KUME[win[i]] for i in k]
    cikti = []
    for kombin in itertools.product(*secenekler):
        L = list(win)
        for i, b in zip(k, kombin):
            L[i] = b
        cikti.append(("".join(L), len(k)))
    return cikti, None


def composition_ok(s, a):
    'The oligo rules from the design decisions. Returns (is_it_suitable, reason).'
    if "N" in s:
        return False, "kalipta_N"
    degen = [c for c in s if c not in "ACGT"]
    if len(degen) > a.degeneracy_budget:
        return False, "dejenere_baz"
    if degen:
        fold = 1
        for c in degen:
            fold *= DEGEN_FOLD.get(c, 4)
        if fold > a.degeneracy_fold_max:
            return False, "dejenere_kat"
        # dejenere baz son bes bazda olamaz, uzama oradan baslar
        nlast = a.gc_clamp_last if a.gc_clamp_last > 0 else 5
        if any(c not in "ACGT" for c in s[-nlast:]):
            return False, "3p_dejenere"
    g = gc_pct(s)
    if not (a.gc_hard_min <= g <= a.gc_hard_max):
        return False, "gc_sert_sinir"
    if a.require_3p_gc and s[-1] not in "GC":
        return False, "3p_gc_degil"
    if a.gc_clamp_last > 0:
        tail = s[-a.gc_clamp_last:]
        if tail.count("G") + tail.count("C") > a.gc_clamp_max:
            return False, "3p_asiri_sabit"
    if max_run(s) > a.homopolymer_max:
        return False, "homopolimer"
    return True, ""


def tm_primer3(s, a):
    return primer3.calc_tm(s, mv_conc=a.mv, dv_conc=a.dv,
                           dntp_conc=a.dntp, dna_conc=a.dna_conc)


def tm_biopython(s, a):
    # The SantaLucia 1998 nearest neighbour table, the Owczarzy 2008 salt correction.
    # The parameter mapping of the two libraries does not have to be identical; the rule
    # rests on measuring the systematic offset and eliminating what departs from it.
    return float(mt.Tm_NN(Seq(s), nn_table=mt.DNA_NN3, Na=a.mv, Mg=a.dv,
                          dNTPs=a.dntp, dnac1=a.dna_conc, dnac2=0,
                          saltcorr=7))


# ------------------------------------------------------------------------ ana
def main():
    a = get_args()
    label = a.label or os.path.basename(a.consensus).split("_consensus")[0]
    name, seq = read_fasta(a.consensus)
    L = len(seq)
    mask_contig = a.mask_contig if a.mask_contig else None
    masked, mclasses, mseen = read_mask(a.mask, L, contig=mask_contig)
    if a.mask and not masked:
        print(u'WARNING: the mask file was read but no position was forbidden. Contig names in the file: %s' % (list(mseen) or "yok"))
    print(u'target           : %s' % label)
    print(u'consensus        : %s (%d bp, header %s)' % (a.consensus, L, name))
    print("maske            : %s" % (a.mask or "yok"))
    if mclasses:
        print("maske siniflari  : %s" % ", ".join("%s=%d" % kv for kv in sorted(mclasses.items())))
    if L < a.len_min * 2 + a.prod_min:
        sys.exit(u'ERROR: the consensus is too short (%d bp), at least %d bp is needed'
                 % (L, a.len_min * 2 + a.prod_min))
    print(u'forbidden positions : %d (%.2f%%)' % (len(masked), 100.0 * len(masked) / L))
    print(u'degeneracy budget   : %d positions, at most %d fold'
          % (a.degeneracy_budget, a.degeneracy_fold_max))

    # --- 1. kompozisyon suzgeci, ucuz olan once ------------------------
    reasons = {}
    raw = []          # (yon, baslangic0, uzunluk, oligo_dizisi, cozulen_iupac)
    for start in range(L):
        for ln in range(a.len_min, a.len_max + 1):
            end = start + ln
            if end > L:
                break
            if any(i in masked for i in range(start, end)):
                reasons["maskeli"] = reasons.get("maskeli", 0) + 2  # F ve R
                continue
            win = seq[start:end]
            for strand in ("F", "R"):
                varyant, why = iupac_varyantlar(win, a, uc=strand)
                if varyant is None:
                    reasons[why] = reasons.get(why, 0) + 1
                    continue
                for coz, kac in varyant:
                    oligo = coz if strand == "F" else rc(coz)
                    ok, why = composition_ok(oligo, a)
                    if not ok:
                        reasons[why] = reasons.get(why, 0) + 1
                        continue
                    raw.append((strand, start, ln, oligo, kac))
    print(u'\noligos after the composition filter: %d' % len(raw))
    for k, v in sorted(reasons.items(), key=lambda x: -x[1]):
        print("   elenen %-16s %d" % (k, v))
    if not raw:
        sys.exit(u'no oligo passed the composition filter')

    # --- 2. iki bagimsiz Tm olcumu ve sistematik kayma ------------------
    tm3, tmb = [], []
    for _, _, _, oligo, _ in raw:
        tm3.append(tm_primer3(oligo, a))
        tmb.append(tm_biopython(oligo, a))
    diffs = [x - y for x, y in zip(tm3, tmb)]
    offset = statistics.median(diffs)
    spread = statistics.pstdev(diffs) if len(diffs) > 1 else 0.0
    if a.tm_cross_tol is not None:
        tol = a.tm_cross_tol
        tol_src = "elle verildi"
    else:
        tol = max(0.10, a.tm_cross_k * spread)
        tol_src = "veriden turetildi: %.1f carpi sd (%.3f), taban 0,10" % (
            a.tm_cross_k, spread)
    print(u'\nTm comparison between the two libraries (%d oligos)' % len(raw))
    print(u'   primer3 minus Biopython, median offset : %+.2f C' % offset)
    print(u'   standard deviation of the offset       : %.3f C' % spread)
    print(u'   tolerance used                         : %.3f C (%s)' % (tol, tol_src))

    kept = []
    n_cross = n_tm = n_hp = n_hd = 0
    for (strand, start, ln, oligo, kac_iupac), t3, tb in zip(raw, tm3, tmb):
        if abs((t3 - tb) - offset) > tol:
            n_cross += 1
            continue
        tm = t3
        if not (a.tm_hard_min <= tm <= a.tm_hard_max):
            n_tm += 1
            continue
        hp = primer3.calc_hairpin(oligo, mv_conc=a.mv, dv_conc=a.dv,
                                  dntp_conc=a.dntp, dna_conc=a.dna_conc).dg
        if hp < a.hairpin_dg_min:
            n_hp += 1
            continue
        hd = primer3.calc_homodimer(oligo, mv_conc=a.mv, dv_conc=a.dv,
                                    dntp_conc=a.dntp, dna_conc=a.dna_conc).dg
        if hd < a.homodimer_dg_min:
            n_hd += 1
            continue
        kept.append(dict(strand=strand, start=start, ln=ln, oligo=oligo,
                         tm3=t3, tmb=tb, hairpin_dg=hp, homodimer_dg=hd,
                         gc=gc_pct(oligo), iupac=kac_iupac))
    print(u'   dropped, the two measurements diverged : %d' % n_cross)
    print(u'   dropped, Tm outside the hard limits    : %d' % n_tm)
    print(u'   dropped, hairpin dG                    : %d' % n_hp)
    print(u'   dropped, self-dimer dG                 : %d' % n_hd)
    print(u'oligos after the thermodynamic filter     : %d' % len(kept))
    if not kept:
        sys.exit(u'no oligo passed the thermodynamics filter')

    F = sorted([k for k in kept if k["strand"] == "F"], key=lambda x: x["start"])
    R = sorted([k for k in kept if k["strand"] == "R"], key=lambda x: x["start"])
    print(u'   forward candidates: %d   reverse candidates: %d' % (len(F), len(R)))

    # --- 3. ciftleme, urun kurallari ve makine dogrulamasi -------------
    pairs = []
    n_prod = n_tmdiff = n_het = n_overlap = 0
    n_verify_fail = 0
    Rby = {}
    for r in R:
        Rby.setdefault(r["start"], []).append(r)
    rstarts = sorted(Rby)
    import bisect
    for f in F:
        lo = f["start"] + a.prod_min - a.len_max
        hi = f["start"] + a.prod_hard_max
        i = bisect.bisect_left(rstarts, lo)
        while i < len(rstarts) and rstarts[i] <= hi:
            for r in Rby[rstarts[i]]:
                pstart, pend = f["start"], r["start"] + r["ln"]
                plen = pend - pstart
                if plen < a.prod_min or plen > a.prod_hard_max:
                    n_prod += 1
                    continue
                if abs(f["tm3"] - r["tm3"]) >= a.pair_tm_diff_max:
                    n_tmdiff += 1
                    continue
                # The F and R footprints cannot overlap; if they do it is not an
                # amplifiable product. It stays hidden at the default thresholds
                # and comes out once --prod-min is lowered.
                if f["start"] + f["ln"] > r["start"]:
                    n_overlap += 1
                    continue
                product = seq[pstart:pend]
                # makine dogrulamasi: urunun basi ileri primere, sonu geri
                # primerin ters tumleyenine birebir esit olmali
                if product[:f["ln"]] != f["oligo"] or product[-r["ln"]:] != rc(r["oligo"]):
                    n_verify_fail += 1
                    continue
                het = primer3.calc_heterodimer(f["oligo"], r["oligo"], mv_conc=a.mv,
                                               dv_conc=a.dv, dntp_conc=a.dntp,
                                               dna_conc=a.dna_conc).dg
                if het < a.heterodimer_dg_min:
                    n_het += 1
                    continue
                pgc = gc_pct(product)
                # the scoring: a smaller penalty sum is better
                pen = 0.0
                pen += 0.0 if a.prod_best_min <= plen <= a.prod_best_max else \
                    min(abs(plen - a.prod_best_min), abs(plen - a.prod_best_max)) * 0.05
                # --prod-max is a soft upper bound: between 250 and 300 is accepted
                # but penalised, so the parameter does not stay dead
                if plen > a.prod_max:
                    pen += (plen - a.prod_max) * 0.10
                for t in (f["tm3"], r["tm3"]):
                    if not (a.tm_min <= t <= a.tm_max):
                        pen += min(abs(t - a.tm_min), abs(t - a.tm_max)) * 2.0
                for g in (f["gc"], r["gc"]):
                    if not (a.gc_min <= g <= a.gc_max):
                        pen += min(abs(g - a.gc_min), abs(g - a.gc_max)) * 0.2
                if not (a.prod_gc_min <= pgc <= a.prod_gc_max):
                    pen += min(abs(pgc - a.prod_gc_min), abs(pgc - a.prod_gc_max)) * 0.2
                pen += abs(f["tm3"] - r["tm3"]) * 1.0
                pen += max(0.0, -het / 1000.0) * 0.5
                pen += max(0.0, -(f["hairpin_dg"] + r["hairpin_dg"]) / 1000.0) * 0.3
                pairs.append(dict(
                    hedef=label,
                    ileri_dizi=f["oligo"], ileri_baslangic=f["start"] + 1,
                    ileri_uzunluk=f["ln"], ileri_tm_primer3=round(f["tm3"], 2),
                    ileri_tm_biopython=round(f["tmb"], 2), ileri_gc=round(f["gc"], 1),
                    ileri_hairpin_dg=round(f["hairpin_dg"], 1),
                    ileri_selfdimer_dg=round(f["homodimer_dg"], 1),
                    geri_dizi=r["oligo"], geri_baslangic=r["start"] + 1,
                    geri_uzunluk=r["ln"], geri_tm_primer3=round(r["tm3"], 2),
                    geri_tm_biopython=round(r["tmb"], 2), geri_gc=round(r["gc"], 1),
                    geri_hairpin_dg=round(r["hairpin_dg"], 1),
                    geri_selfdimer_dg=round(r["homodimer_dg"], 1),
                    cift_heterodimer_dg=round(het, 1),
                    tm_farki=round(abs(f["tm3"] - r["tm3"]), 2),
                    urun_uzunluk=plen, urun_gc=round(pgc, 1),
                    urun_baslangic=pstart + 1, urun_bitis=pend,
                    urun_dogrulandi="evet_ayni_kalip",
                    ceza=round(pen, 3)))
            i += 1

    print(u'\npair construction')
    print(u'   dropped, product length                : %d' % n_prod)
    print(u'   dropped, pair Tm difference            : %d' % n_tmdiff)
    print(u'   dropped, F and R footprints overlap    : %d' % n_overlap)
    print(u'   dropped, machine verification of the product : %d' % n_verify_fail)
    print(u'     (note: at this stage the product is cut from the same copy of the')
    print(u'      template the primers were derived from, so this check passes by')
    print(u'      definition. Its real purpose appears at the specificity stage, where')
    print(u'      it runs against competitor and raw-read templates)')
    print(u'   dropped, hetero-dimer dG               : %d' % n_het)
    print(u'valid pairs                              : %d' % len(pairs))
    if not pairs:
        sys.exit(u'no valid pair was found')

    pairs.sort(key=lambda x: x["ceza"])
    if a.min_locus_spacing > 0:
        kept_pairs, taken = [], []
        for p in pairs:
            if all(abs(p["ileri_baslangic"] - q[0]) >= a.min_locus_spacing or
                   abs(p["geri_baslangic"] - q[1]) >= a.min_locus_spacing
                   for q in taken):
                kept_pairs.append(p)
                taken.append((p["ileri_baslangic"], p["geri_baslangic"]))
        print(u'locus thinning (%d bp): %d -> %d pairs'
              % (a.min_locus_spacing, len(pairs), len(kept_pairs)))
        pairs = kept_pairs
    if len(pairs) > a.max_pairs:
        print(u'output trimmed to %d rows (by score)' % a.max_pairs)
        pairs = pairs[:a.max_pairs]

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pairs[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(pairs)
    print("\nyazildi: %s" % a.out)
    print(u'\nFive best candidates:')
    for p in pairs[:5]:
        print(u'  penalty=%.2f  product=%d bp (GC %.1f)  F=%s (Tm %.1f)  R=%s (Tm %.1f)'
              % (p["ceza"], p["urun_uzunluk"], p["urun_gc"], p["ileri_dizi"],
                 p["ileri_tm_primer3"], p["geri_dizi"], p["geri_tm_primer3"]))
    print(u'\nNext stage: the specificity filter. These candidates will be tested against')
    print(u'competitor consensus sequences and the REFERANS_DB databases with mfeprimer and blastn.')


if __name__ == "__main__":
    main()
