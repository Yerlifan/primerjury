"""Fast in-silico PCR at raw read level (a seed plus verification).
okuma_pcr.py IS NOT USED; this was written independently. The filter is 200-6000 bp
(corrected).

"""
# -------------------------------------------------------------------------
# reads.py - in-silico PCR on the sample's RAW READS (fastq); it uses ispcr's
#            criterion but looks for candidate positions with a short seed first.
#
# INPUT  : the fastq / fastq.gz bin files (the reads), and the F and R primer
#          sequences. From the command line: python reads.py <F> <R> <fastq...>
# OUTPUT : it WRITES NO FILE. kutu_pcr -> (total_reads, reads_with_product,
#          product_length_counter). __main__ prints a percentage line per bin.
# CALLED BY: screening/engine_gateway.py finds and loads this file by name (unlike
#          ispcr it is NOT REQUIRED; if it cannot be found, okuma stays None), and
#          immediately after loading it OVERWRITES the read length filter with the
#          corrected values from the configuration: okuma.MINL, okuma.MAXL =
#          C.NUMUNE_OKUMA_MIN, C.NUMUNE_OKUMA_MAX. Beside that, sample.py in the
#          same directory and engine/numune_olc.py, kutu_cache.py say "import okuma".
#          In the menu it is loaded through engine_gateway.py on every measuring key:
#          P, K, D, I, G, T, U, H and 1-9. full_chain.py does not call this file
#          directly.
#
# CAUTION - THIS FILE IS THE PANEL'S OLD ENGINE:
# The Sonda class below uses ONE FIXED 13 base seed, and that seed can be lossy (the
# full reasoning is at the head of Sonda.__init__). The panel's AUTHORITATIVE
# measurement route today is screening/read_engine.py (pigeonhole seeding,
# lossless). self_test.py runs the two engines side by side and measures and reports
# how many binding sites the Sonda here misses.
# -------------------------------------------------------------------------
import os, sys, glob, itertools, json, gzip
from collections import Counter
import ispcr

# The read length filter. Nanopore bins hold very short (fragmented) and very long
# (chimeric or joined) reads, and both distort the product ratio. engine_gateway.py
# replaces these two values at load time with the corrected values from the
# configuration.
MINL, MAXL = 200, 6000
# Tohum uzunlugu. SABIT tek tohum - bkz. Sonda.__init__ icindeki uyari.
SEED = 13


# Expands a degenerate (IUPAC) seed into every possible CONCRETE sequence. A seed
# holding two Rs, for instance, becomes 4 separate strings. It is needed because
# str.find can only search for concrete text; the list grows exponentially with the
# number of degenerate bases.
def variants(p):
    return [''.join(x) for x in itertools.product(*[ispcr.IUPAC.get(c, 'ACGT') for c in p])]


# Once the seed matches, it counts mismatches for the FULL primer. It exits early
# with -1 as soon as max_mm is exceeded (there is no point counting the rest).
# Note: like ispcr.find_sites it is IUPAC aware, so a degenerate primer base counts
# every base it allows as a match.
def mm_ok(primer, win, max_mm):
    mm = 0
    for a, b in zip(primer, win):
        if b not in ispcr.IUPAC.get(a, 'ACGT'):
            mm += 1
            if mm > max_mm:
                return -1
    return mm


class Sonda:
    'Finds the binding sites of one primer on a sequence read in a given orientation.'

    # -----------------------------------------------------------------------
    # THE SEEDING - THIS CAN BE LOSSY (there is NO pigeonhole guarantee)
    #
    # This class takes ONE piece of SEED=13 bases from the primer's 3' end (from the
    # start of the reverse complement) and looks only at the places where that piece
    # matches EXACTLY. The problem: with max_mm >= 1 the mismatch may perfectly well
    # fall INSIDE those 13 bases. The seed then never matches, no candidate position is
    # produced, and a real binding site is SILENTLY missed. The code raises no error,
    # the number simply comes out small.
    #
    # WHAT IT WOULD TAKE TO BE LOSSLESS (pigeonhole):
    # In a search allowing at most m mismatches, if the sequence is split into (m + 1)
    # NON-OVERLAPPING pieces, at most m pieces can be spoiled; AT LEAST ONE of the
    # remaining pieces MUST match exactly. So asking "does any of the pieces match
    # exactly" misses no site that brute force would find. That is A GUARANTEE, not a
    # heuristic speed-up. If the piece count drops BELOW m + 1 the guarantee collapses,
    # and using a single piece as here (that is, assuming m = 0) is exactly that
    # collapsed state when it runs with max_mm = 1.
    #
    # The corrected version is screening/read_engine.py: it splits the primer into
    # max_mm + 1 non-overlapping blocks, scans separately for each block, and re-verifies
    # every candidate found under the full rule. That is the panel's measurement route
    # today; this class stands for historical comparison.
    # -----------------------------------------------------------------------
    def __init__(self, primer, uc5=False, max_mm=1):
        # uc5=True -> 3' kritik uc dizinin BASINDA (ters tumleyen primer)
        self.p = primer
        self.max_mm = max_mm
        self.uc5 = uc5
        # The seed is taken from the primer's 3' end: that is the end the polymerase
        # extends, and therefore the best conserved region in a real binding. For rc(R) it
        # is taken from the start, because the 3' end falls at the START of that sequence.
        s = primer[:SEED] if uc5 else primer[-SEED:]
        # off: the seed's start offset inside the primer. When the seed is found at
        # position i in the sequence, the primer starts at i - off.
        self.off = 0 if uc5 else len(primer) - SEED
        self.seeds = variants(s)

    # It tries EVERY place the seed occurs (the find loop continues from i+1, so
    # overlapping repeats are caught too), then verifies each candidate against the full
    # primer. Returns: [(start, mismatches)].
    def bul(self, seq):
        out = []
        L = len(self.p)
        for sd in self.seeds:
            i = seq.find(sd)
            while i != -1:
                st = i - self.off
                # Primer okumanin disina tasiyorsa aday gecersiz.
                if 0 <= st and st + L <= len(seq):
                    mm = mm_ok(self.p, seq[st:st + L], self.max_mm)
                    # mm >= 0 -> tam kural altinda da gecti.
                    if mm >= 0:
                        out.append((st, mm))
                i = seq.find(sd, i + 1)
        return out


# The fastq reader. A fastq is made of 4 line records (@header / SEQUENCE / + /
# quality); k % 4 == 1 is exactly the SEQUENCE line. .gz files are opened
# transparently and malformed bytes are skipped with errors='ignore' (that happens
# in nanopore output). The filter discards reads outside MINL..MAXL entirely, and
# those do not enter the denominator either.
def okumalar(path):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt', errors='ignore') as fh:
        for k, line in enumerate(fh):
            if k % 4 == 1:
                s = line.strip().upper()
                if MINL <= len(s) <= MAXL:
                    yield s


# -------------------------------------------------------------------------
# kutu_pcr - the proportion of reads giving a product among a bin's (fastq) raw reads.
#
# WHAT IT COMPUTES: (total_reads, reads_with_product, product_length_counter). The
# result is A RATIO: pos / tot.
#
# WHY ON THE RAW READS RATHER THAN ON A REFERENCE SEQUENCE:
# A reference or consensus is a single summarised sequence; it does not hold the
# variants actually present in the sample, the within species differences, or the
# nanopore error pattern. A pair can give a perfect product on a reference and not
# hold on the sample's real reads. The panel's deciding number is produced at raw
# read level for that reason.
#
# ORIENTATION - TRYING BOTH DIRECTIONS IS REQUIRED:
# Nanopore reads are BIDIRECTIONAL; the same region is written as the plus strand in
# some reads and the minus strand in others. The engines (ispcr.find_sites and the
# Sonda here) scan only in the forward direction. So every read is tried both as
# itself and as its reverse complement: "for seq in (s, ispcr.rc(s))". Without that
# loop roughly half the reads would be lost silently, with no error raised.
#
# THE COUNTING RULE: the moment a read gives a product the loop "break"s; a read is
# counted AT MOST ONCE. So this is a PRESENCE measurement (in how many reads is
# there a product), not a depth measurement. In the same way, the first matching
# product length is recorded and the inner loops are broken.
# -------------------------------------------------------------------------
def kutu_pcr(path, F, R, lo=40, hi=600, max_mm=1):
    """The number of reads giving a product in a bin, and the total reads."""
    # F is searched in the forward direction; R as its reverse complement, with
    # uc5=True because its 3' end falls at the START of the sequence (the same idea as
    # tail_pos=(0,1) in ispcr.amplify).
    fs = Sonda(F, False, max_mm)
    rs = Sonda(ispcr.rc(R), True, max_mm)
    tot = pos = 0
    sizes = Counter()
    for s in okumalar(path):
        tot += 1
        # Both directions of the read are tried - see the ORIENTATION note above.
        for seq in (s, ispcr.rc(s)):
            a = fs.bul(seq)
            # With no F there can be no product in this direction; no need to look for R.
            if not a:
                continue
            b = rs.bul(seq)
            if not b:
                continue
            got = False
            for i, _ in a:
                for j, _ in b:
                    # Urun boyu: F'nin basindan R'nin (tumleyeninin) sonuna.
                    n = j + len(R) - i
                    # The length window plus the primers facing one another without overlapping
                    # (exactly the same criterion as ispcr.amplify).
                    if lo <= n <= hi and j >= i + len(F):
                        sizes[n] += 1
                        got = True
                        break
                if got:
                    break
            if got:
                pos += 1
                # The read is already counted; there is no need to try the other direction.
                break
    return tot, pos, sizes


if __name__ == '__main__':
    F, R = sys.argv[1], sys.argv[2]
    for p in sys.argv[3:]:
        t, n, sz = kutu_pcr(p, F, R)
        print(f'{os.path.basename(p):34s} {n:7d}/{t:<8d} {100*n/max(t,1):6.2f}%  {sz.most_common(3)}')
