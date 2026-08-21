"""A self-contained in-silico PCR engine (2026-08-02, the final stage).
The engine_gateway.py / kararlar.py code of the earlier session is NOT USED; this
was written from scratch.
The criterion: F and R must bind facing one another, the last 2 bases at the 3' end
must match exactly, the total mismatch count must be <= max_mm, and the product
length must fall inside the given window.

"""
# ---------------------------------------------------------------------------
# ispcr.py - the panel's core in-silico PCR engine. It measures whether the F
#            and R primers bind facing one another and give a product.
#
# INPUT  : fasta files (read_fasta) or sequence text directly; the F and R primer
#          sequences (IUPAC degenerate bases accepted). From the command line:
#          python ispcr.py <F> <R> <fasta...>
# OUTPUT : it WRITES NO FILE, it returns values. find_sites -> [(start, mm)];
#          amplify -> [(start, end, product_bp, mm_F, mm_R)]; scan_file ->
#          (total, with_product, length_distribution, hit_list).
#          Only __main__ prints a summary to the screen.
# CALLED BY: screening/engine_gateway.py finds this file by name and loads it. If
#          it cannot, it says "ispcr.py was not found" and exits, so no
#          measurement runs without engine_gateway.py. engine_gateway.py itself is
#          used by every measuring module in the package (sample.py,
#          global_scan.py, reference.py, generator.py, targets.py,
#          build_consensus.py, panel_measurement.py, membership_check.py,
#          self_test.py) and from outside by verification/recovery_round.py and
#          protocol/single_protocol_measure.py. Beside those, ara.py, hiza.py,
#          mazei*.py, deg_*.py, mmb_*.py and engine/scanner.py, pair.py say
#          "import ispcr" directly.
#          It loads on every measuring key in the menu: P (panel measurement under
#          a single protocol), K (recovery), D (verification), I and G (identity),
#          T (P->K->D->I), U (membership), H (quick test), 1-9 (search).
#          full_chain.py does not call this file directly, but on start-up it
#          checks that engine/ispcr.py exists; without it the tool does not open.
#
# WHY THERE IS NO SEED: this engine scans a full sliding window inside find_sites
# and makes no seed or shortcut assumption, so its losslessness holds by
# definition. The seeded fast paths (screening/read_engine.py, scanner.py Havuz)
# must imitate this file's criterion, and self_test.py compares all three.
# ---------------------------------------------------------------------------
import re, sys, os, glob
from collections import defaultdict

# The IUPAC degenerate base table: which bases in the sequence a primer letter may
# match. An R at a primer position matches both A and G. A degenerate primer is
# not one oligo but a MIXTURE of oligos, and this table is how the search side
# expands it.
IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT',
}
COMP = str.maketrans('ACGTURYSWKMBDHVNacgturyswkmbdhvn',
                     'TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn')


# Reverse complement. In PCR the R primer binds the MINUS strand; since the engine
# scans only the plus strand, it has to convert R to rc(R) before searching for it
# (see amplify).
def rc(s):
    return s.translate(COMP)[::-1]


# Normalises a sequence: every character outside ACGTUN (a space, a degenerate
# base, an alignment dash, a lower case leftover) becomes N, and U becomes T.
# This matters: those Ns become -1 inside encode() and NEVER match in find_sites.
# So an ambiguous base is never silently counted as a match; it counts as a
# mismatch. The measurement is cautious rather than optimistic.
def clean(s):
    return re.sub(r'[^ACGTUN]', 'N', s.upper()).replace('U', 'T')


# The fasta reader. It does not hold the whole file in memory, it yields record by
# record: the panel's reference and consensus files can hold many records.
def read_fasta(path):
    name, buf = None, []
    with open(path, errors='ignore') as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(buf)
                name, buf = line[1:].rstrip('\n'), []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, ''.join(buf)


import numpy as np

# The byte to base code table. A, C, G, T map to 0, 1, 2, 3; EVERY OTHER BYTE
# stays -1. The -1 is deliberate: the "col >= 0" condition in find_sites makes N
# and every undefined character count automatically as a mismatch.
_B = np.zeros(256, dtype=np.int8)
_B[:] = -1
for _i, _c in enumerate('ACGT'):
    _B[ord(_c)] = _i


# Turns the sequence into an int8 array so the scan can use numpy column
# operations instead of a pure Python loop. Speed comes from here, not the criterion.
def encode(seq):
    a = np.frombuffer(seq.encode(), dtype=np.uint8)
    return _B[a]


# ---------------------------------------------------------------------------
# find_sites - every start position where the primer binds in the FORWARD direction.
#
# WHAT IT COMPUTES: for each possible start position, (a) the total mismatch count
# and (b) whether the critical 3' end positions match exactly. The acceptance
# condition: mm <= max_mm and (when need_tail) ALL of tail_pos match exactly.
#
# WHY IT IS COMPUTED THIS WAY - NO SEED:
# This function scans a FULL SLIDING WINDOW over the sequence. The outer loop runs
# over the primer's bases (L steps), and at each step every alignment of that base
# in the sequence is compared in one go, so no start position is skipped without
# being evaluated. That is what makes losslessness hold by definition.
#
# ITS RELATION TO PIGEONHOLE SEEDING: a seeded engine looks first, for speed, at
# the places where a short piece matches EXACTLY. That can only be lossless under
# this guarantee: in a search allowing at most m mismatches, if the sequence is
# split into (m + 1) NON-OVERLAPPING pieces, at most m pieces can be spoiled and
# AT LEAST ONE piece MUST match exactly. That is a guarantee, not a heuristic
# speed-up. If the piece count drops BELOW m + 1 the guarantee collapses and the
# search starts missing binding sites silently. It raises no error either, the
# numbers simply come out small. This file uses no seed at all, so it does not
# carry that risk, and that is exactly why the correctness of the seeded paths
# (read_engine.py, scanner.py) is tested against this function and brute_force.py.
#
# THE 3' END RULE (tail_pos): the polymerase extends a primer from its 3' end. If
# there is a mismatch at the 3' end, extension does not happen in practice, so the
# last 2 bases at the 3' end are required to match EXACTLY. That is the panel's
# criterion. A negative tail_pos (-1, -2) makes the end of the primer critical,
# and (0, 1) makes the start critical; the second form is what the complement of
# the reverse primer needs (see amplify).
# ---------------------------------------------------------------------------
def find_sites(enc, primer, max_mm, need_tail=True, tail_pos=(-1, -2)):
    """The 0-based starts where the primer binds the sequence in the forward direction (numpy).
    tail_pos: the indices, inside the primer, of the critical 3' positions.

    """
    L, n = len(primer), enc.shape[0]
    if n < L:
        return []
    # m = the number of possible start positions, each a separate candidate alignment.
    m = n - L + 1
    # mm[i]  : the total mismatch count of the alignment at position i
    # tail_ok: whether every critical 3' position matched exactly at position i
    mm = np.zeros(m, dtype=np.int16)
    tail_ok = np.ones(m, dtype=bool)
    # Negative indices (-1, -2) are converted to positive against the primer length.
    tpos = {(p % L) for p in tail_pos}
    # The outer loop runs over the PRIMER bases (L steps) and the inner work is a numpy
    # vector, so every start position is evaluated at the same time. No position is
    # skipped and no seed filtering happens.
    for k, p in enumerate(primer):
        allowed = [ 'ACGT'.index(c) for c in IUPAC.get(p, 'ACGT') ]
        # col: the sequence base at position k of each candidate alignment.
        col = enc[k:k + m]
        ok = np.zeros(m, dtype=bool)
        # A degenerate primer base can allow several bases; all of them are tried.
        for a in allowed:
            ok |= (col == a)
        ok &= (col >= 0)          # N or a gap counts as a mismatch
        mm += (~ok)
        # At a critical 3' position there is NO tolerance: one mismatch rules it out.
        if k in tpos:
            tail_ok &= ok
    sel = mm <= max_mm
    if need_tail:
        sel &= tail_ok
    idx = np.nonzero(sel)[0]
    # (start, the mismatch count of that binding) pairs are returned. The caller uses
    # the mismatch count to pick the "best product".
    return list(zip(idx.tolist(), mm[idx].tolist()))


# ---------------------------------------------------------------------------
# amplify - IN-SILICO PCR. It measures whether the F and R primers can really give
#           a product on a sequence.
#
# WHAT COUNTS AS A PRODUCT (the four conditions of real PCR):
#   1) F must bind the sequence in the forward direction (find_sites, the last 2
#      bases at the 3' end exact).
#   2) R must bind the OPPOSITE strand. Since the engine scans only the plus
#      strand, it searches not for R itself but for its reverse complement rc(R).
#      If rc(R) is found on the plus strand, R really does bind the minus strand.
#   3) THE TWO PRIMERS MUST FACE ONE ANOTHER: F extends left to right, R right to
#      left. The condition j >= i + len(fwd) requires R's binding site to lie to
#      the RIGHT of F and NOT OVERLAP it. Without that, two bindings sitting on
#      the same region would produce a false "product"; and a pair in the wrong
#      order (R on the left, F on the right) gives nothing in real PCR, so it
#      gives nothing here either.
#   4) THE PRODUCT LENGTH must fall in the window: lo <= size <= hi. This has a
#      qPCR counterpart. config.py holds URUN_IDEAL (60-150 bp), URUN_KABUL
#      (150-250) and URUN_MUTLAK_UST 400. A short product completes within the
#      30 s extension of the protocol and gives a clean SYBR Green signal; a very
#      long product neither amplifies efficiently nor appears where it should on
#      a gel or a melt curve. With no window the engine would call two bindings
#      5 kb apart a "product".
#
# WHY THE 3' END INDEX FLIPS: the 3' end of the R primer is the 5' END of the
# rc(R) sequence. That is why tail_pos=(-1,-2) is given for F while (0,1) is given
# for rc(R). Miss this and the 3' rule is applied from the wrong end on the R side.
#
# ORIENTATION WARNING: this function scans ONLY THE PLUS STRAND of the sequence it
# is given; it does not try rc(seq) on its own. Nanopore reads come in both
# directions and a consensus may have been built in reverse. On a reversed
# sequence this function loses EVERY product and RAISES NO ERROR, it quietly
# returns 0. So:
#   - on the consensus side the orientation is read FROM THE CANONICAL SOURCE
#     (screening/hedefler.konsensusler -> konsensus_kanonik/INDEKS.tsv, all SENSE).
#     The old "consensus sequences" directory was MIXED orientation
#     (71 antisense / 27 sense) and falling back to it silently is forbidden.
#   - on the raw read side the caller tries every read in both directions
#     ("for seq in (s, rc(s))" inside okuma.kutu_pcr and numune.pcr_kutu), and
#     motor.urun_var tries both directions in the same way.
# ---------------------------------------------------------------------------
def amplify(seq, fwd, rev, max_mm=3, lo=40, hi=600, need_tail=True, enc=None):
    """Return the products: (start, end, product_bp, mm_F, mm_R)."""
    if enc is None:
        enc = encode(clean(seq))
    revrc = rc(rev)
    fs = find_sites(enc, fwd, max_mm, need_tail, tail_pos=(-1, -2))
    # If F does not bind at all there is no point searching for R; there can be no product.
    if not fs:
        return []
    # The 3' criterion for revrc: the 3' end of rev = the 5' end of revrc -> index 0,1
    rs = find_sites(enc, revrc, max_mm, need_tail, tail_pos=(0, 1))
    if not rs:
        return []
    prods = []
    # Every F x R binding pair is tried; one sequence can give more than one product.
    for i, mmf in fs:
        for j, mmr in rs:
            # The product runs from the start of F to the END of rc(R), so both primers are
            # themselves inside the product, as they are in a real amplicon.
            end = j + len(revrc)
            size = end - i
            # The length window, plus the primers facing one another without overlapping.
            if lo <= size <= hi and j >= i + len(fwd):
                prods.append((i, end, size, mmf, mmr))
    return prods


# ---------------------------------------------------------------------------
# scan_file - counts how many sequences in a fasta file give a product.
#
# WHAT IT COMPUTES: (total, with_product, length_distribution, hit_list).
# The DENOMINATOR of the ratio is the number of sequences left AFTER those shorter
# than min_len are dropped. The reason: if a truncated or partial record is not
# even long enough to carry the product, putting it in the denominator pushes the
# ratio down falsely.
#
# A sequence is counted ONCE whether or not it gives a product (this measures
# presence, not depth). If one sequence gives several products, the "best" one is
# the one with the smallest total mismatch count, since that is the likeliest
# binding.
#
# CAUTION: clean(seq) is passed here but rc(seq) is NOT TRIED. A reversed fasta
# quietly gives 0 products in this function; the orientation must come from the
# canonical source (see the orientation warning in the amplify header).
# ---------------------------------------------------------------------------
def scan_file(path, fwd, rev, max_mm=3, lo=40, hi=600, min_len=0):
    """Scan one fasta. The denominator is the sequences longer than min_len."""
    tot = amp = 0
    sizes = defaultdict(int)
    hits = []
    for name, seq in read_fasta(path):
        # Records too short to enter the denominator are skipped entirely.
        if len(seq) < min_len:
            continue
        tot += 1
        p = amplify(clean(seq), fwd, rev, max_mm, lo, hi)
        if p:
            amp += 1
            # The best product = the one with the smallest mm_F + mm_R.
            best = min(p, key=lambda x: x[3] + x[4])
            sizes[best[2]] += 1
            hits.append((name, best[2], best[3], best[4]))
    return tot, amp, dict(sizes), hits


# Manual use: python ispcr.py <F> <R> <fasta...>
# min_len=1200 - command line use targets full length 16S/ITS records, so that
# short fragments do not pollute the denominator.
if __name__ == '__main__':
    fwd, rev = sys.argv[1].upper(), sys.argv[2].upper()
    for path in sorted(sys.argv[3:]):
        t, a, s, _ = scan_file(path, fwd, rev, min_len=1200)
        print(f'{os.path.basename(path):40s} {a:5d}/{t:5d}  {sorted(s.items())[:4]}')
