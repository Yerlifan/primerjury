# -*- coding: utf-8 -*-
"""
read_engine.py - THE CORRECTED raw-read in-silico PCR engine.
2026-08-02, the read engine fix.

IT IS SELF-CONTAINED: numpy, ispcr, reads.py, pysam and minimap2 ARE NOT NEEDED.
Pure standard Python. Windows or Linux makes no difference.

---------------------------------------------------------------------------
WHAT IT FIXES
---------------------------------------------------------------------------
The old engine (engine/reads.py -> class `Sonda`, and also engine/scb.py -> class
`S`) searched for the primer in a read using a single 13 base EXACT MATCHING seed:

    s = primer[:13] if uc5 else primer[-13:]
    i = seq.find(sd)          # <-- an EXACT match requirement

Although the criterion is "total mismatches <= max_mm", when the mismatch falls
inside those 13 bases `find` finds nothing and the binding site disappears
SILENTLY. The program raises no error, it says "no product".

The measured effect (bin A1-4_3078083, the first 400 reads, the M. mazei pair,
THE SAME criterion mm<=1):
    the old Sonda   :   2/400  (0.50%)
    brute force     : 174/400  (43.50%)
    -> in 202 of 205 binding sites (98.5%) the mismatch falls inside the seed
       (188 of them at base 6 of the primer, a single recurring variant base).

---------------------------------------------------------------------------
HOW IT WAS FIXED - PIGEONHOLE SEEDING
---------------------------------------------------------------------------
The primer is split into max_mm+1 NON-OVERLAPPING blocks. If there are at most
max_mm mismatches then, by the pigeonhole principle, AT LEAST ONE of those blocks
must match exactly. So searching for an exact match of any block is LOSSLESS: it
misses no site brute force would find. Every candidate found is then verified ONE
BY ONE under the full rule (the total mismatch count plus the last 2 bases at the
3' end).

Also: in the old code the last-two-bases-at-the-3'-end EXACT MATCH rule was
applied nowhere explicitly; it was a side effect of the 13 base seed. Because that
side effect disappears once the seed shortens, the rule is applied EXPLICITLY in
this module (son2=True).

THE LOSSLESSNESS CLAIM IS TESTED: engine_test.py compares this engine against an
independent brute force implementation (brute_force.py) that uses no seed and
tries every position one at a time. Run that test before entering the panel.

---------------------------------------------------------------------------
USAGE
---------------------------------------------------------------------------
As a module:
    import read_engine as om
    okumalar = list(om.okumalar('kutu.fastq'))            # a 200-6000 bp filter
    pos, n   = om.kutu_pcr(okumalar, F, R, max_mm=1)      # reads giving a product
    yerler   = om.Sonda(F, uc5=False, max_mm=1).bul(dizi) # [(start, mismatches), ...]

From the command line:
    python read_engine.py F_PRIMER R_PRIMER [--mm 1] [--lo 40] [--hi 600] \
                           [--nmax 3000] [--seed 3] [--tsv out.tsv] file1.fastq ...
    -> one line per fastq: file, with_product, reads_used, total_reads, percent,
       dominant_lengths

    An example:
      python read_engine.py GCCCTTGGGACCGGCATAA TCGCTGGCTAGTAGGTACATTACA \
             --mm 1 --tsv sonuc.tsv "fastq files/A1-4/*.fastq"

THE CRITERION LABEL: the --mm value is written into the output. The panel's sample
criterion is mm<=1; mm<=3 is the design pipeline's criterion. DO NOT CONFUSE THE
TWO - they had been confused in the panel's old rows (see the "16 Okuma Motoru
Duzeltmesi" sheet).

"""
# -------------------------------------------------------------------------
# read_engine.py - the corrected engine that does in-silico PCR on raw nanopore
#                   reads; it uses pigeonhole seeding and misses no binding site.
#
# INPUT  : fastq / fastq.gz files (through okumalar() and kutu_yukle(), with a
#          200-6000 bp length filter and fixed seed sampling) plus the forward and
#          reverse primer sequences. It is used both as a module and from the
#          command line.
# OUTPUT : it writes no file (on the command line, --tsv writes a one line TSV per
#          file). Sonda.bul() returns a [(start, mismatches)] list and kutu_pcr()
#          returns the triple (with_product, total, length_counter).
# CALLED BY: inside the package, engine_gateway.py imports this file and it is used
#          through numune.KutuOtorite and numune.KutuEski in every measurement stage
#          - that is, full_chain.py stages 1 to 9. Beside that, orientation.py,
#          orientation_audit.py, orientation_impact_test.py, cross_coverage.py,
#          engine_test.py and independent_check.py import it directly.
# -------------------------------------------------------------------------
import sys, os, re, glob, gzip, random, argparse
from collections import Counter

__version__ = '1.0 (2026-08-02)'

MINL, MAXL = 200, 6000          # the read length filter (A2 ~4.2-4.5 kb and F2 ~3.7 kb included)

IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT',
}
_COMP = str.maketrans('ACGTURYSWKMBDHVNacgturyswkmbdhvn',
                      'TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn')


def rc(s):
    """ters tumleyen"""
    return s.translate(_COMP)[::-1]


def temizle(s):
    return re.sub(r'[^ACGTUN]', 'N', s.upper()).replace('U', 'T')


def _blok_regex(alt):
    """bir primer blogunu, IUPAC dejenere bazlari karakter sinifina cevirerek
    ORTUSEN eslesme bulabilen bir regex'e cevirir."""
    govde = ''.join(c if c in 'ACGT' else '[' + IUPAC.get(c, 'ACGT') + ']' for c in alt)
    return re.compile('(?=(' + govde + '))')


class Sonda:
    """Finds ALL the binding sites of a primer in a sequence. Lossless.

        uc5=False -> a forward primer; the critical 3' end is at the END of the primer
                     (index -1, -2)
        uc5=True  -> the complement of a reverse primer; the critical 3' end is at the
                     START of the sequence (index 0, 1)
        son2=True -> the last 2 bases at the 3' end must match EXACTLY (the panel's rule)

    """

    def __init__(self, primer, uc5=False, max_mm=1, son2=True):
        primer = primer.upper()
        self.p = primer
        self.L = len(primer)
        self.max_mm = max_mm
        self.son2 = son2
        self.uc5 = uc5
        self.ok = [set(IUPAC.get(c, 'ACGT')) for c in primer]
        self.kritik = (0, 1) if uc5 else (self.L - 1, self.L - 2)
        # ------------------------------------------------------------------
        # PIGEONHOLE SEEDING - WHY IT IS LOSSLESS
        #
        # If a primer of length k is split into (m + 1) NON-OVERLAPPING pieces in a
        # search allowing at most m mismatches, at most m pieces can be spoiled;
        # AT LEAST ONE of the remaining pieces MUST match exactly.
        # That is a guarantee, not a heuristic speed-up: searching for an exact match
        # of any of the pieces misses no site that brute force would find. Every
        # candidate position found is then tested one by one with _dogrula() under the
        # full rule (the total mismatch count plus the last 2 bases at the 3' end), so
        # it leaves no false positive either.
        #
        # THAT IS WHY THE SEED LENGTH AND THE PIECE COUNT CANNOT BE CHOSEN ARBITRARILY.
        # If the piece count falls BELOW max_mm + 1 the guarantee collapses and the
        # search starts missing sites silently - which was exactly the old engine's bug:
        # it used a single fixed 13 base seed, and when the mismatch fell inside those
        # 13 bases `find` found nothing. The program raises no error, it says "no
        # product". The measured effect: 202 of 205 binding sites (98.5%) were being
        # lost. Raising the piece count is not free either - as the pieces shorten, the
        # number of false candidates and the verification cost both grow.
        # The losslessness claim is verified in engine_test.py against brute_force.py;
        # self_test.py makes the same comparison on every run.
        # ------------------------------------------------------------------
        nb = max_mm + 1                                   # the pigeonhole block count
        kes = [round(i * self.L / nb) for i in range(nb + 1)]
        self.bloklar = [(kes[i], _blok_regex(primer[kes[i]:kes[i + 1]]))
                        for i in range(nb) if kes[i + 1] > kes[i]]

    def _dogrula(self, seq, st):
        """verify a candidate position under the full rule -> the mismatch count, or None"""
        # The last 2 bases at the 3' end EXACT MATCH rule is applied EXPLICITLY HERE. In the
        # old code that rule was written nowhere; it held as a side effect of the long 13
        # base seed. Because that side effect disappears once the seed shortens, the rule was
        # made explicit. Its biological reason: the polymerase starts extension from the 3'
        # end, and a mismatch at that end effectively prevents the product, which makes it a
        # condition independent of the total mismatch count.
        if st < 0 or st + self.L > len(seq):
            return None
        if self.son2:
            for k in self.kritik:
                if seq[st + k] not in self.ok[k]:
                    return None
        mm = 0
        for k in range(self.L):
            if seq[st + k] not in self.ok[k]:
                mm += 1
                if mm > self.max_mm:
                    return None
        return mm

    def bul(self, seq):
        aday = set()
        for off, rx in self.bloklar:
            for m in rx.finditer(seq):
                aday.add(m.start() - off)
        out = []
        for st in sorted(aday):
            mm = self._dogrula(seq, st)
            if mm is not None:
                out.append((st, mm))
        return out


def urun_var(seq, fs, rs, lenF, lenR, lo, hi):
    """tek yonde: F ve R karsilikli baglanip verilen boy penceresinde urun veriyor mu"""
    a = fs.bul(seq)
    if not a:
        return None
    b = rs.bul(seq)
    if not b:
        return None
    for i, _ in a:
        for j, _ in b:
            n = j + lenR - i
            if lo <= n <= hi and j >= i + lenF:
                return n
    return None


def kutu_pcr(okuma_listesi, F, R, lo=40, hi=600, max_mm=1, son2=True):
    """The number of reads giving a product in a bin.
        Returns: (with_product, total, length_counter)

    """
    F = F.upper(); R = R.upper()
    fs = Sonda(F, False, max_mm, son2)
    rs = Sonda(rc(R), True, max_mm, son2)
    pos = 0
    boylar = Counter()
    for s in okuma_listesi:
        for seq in (s, rc(s)):
            n = urun_var(seq, fs, rs, len(F), len(R), lo, hi)
            if n is not None:
                boylar[n] += 1
                pos += 1
                break
    return pos, len(okuma_listesi), boylar


def okumalar(path, minl=MINL, maxl=MAXL):
    """fastq / fastq.gz akisi, uzunluk suzgeciyle"""
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt', errors='ignore') as fh:
        for k, line in enumerate(fh):
            if k % 4 == 1:
                s = line.strip().upper()
                if minl <= len(s) <= maxl:
                    yield s


def kutu_yukle(path, nmax=3000, seed=3, minl=MINL, maxl=MAXL):
    """Bir fastq'tan okumalari al, nmax'i asiyorsa sabit tohumla ornekle.
    Donus: (ornek_okumalar, suzgecten_gecen_toplam)"""
    rs = list(okumalar(path, minl, maxl))
    n0 = len(rs)
    if nmax and len(rs) > nmax:
        random.seed(seed)
        rs = random.sample(rs, nmax)
    return rs, n0


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='Duzeltilmis ham-okuma in-silico PCR (kayipsiz tohumlama).')
    ap.add_argument('F'); ap.add_argument('R')
    ap.add_argument('fastq', nargs='+', help='fastq files or a glob pattern')
    ap.add_argument('--mm', type=int, default=1, help='total mismatches allowed (default 1)')
    ap.add_argument('--lo', type=int, default=40)
    ap.add_argument('--hi', type=int, default=600)
    ap.add_argument('--nmax', type=int, default=3000, help='maximum reads per bin (0 = all)')
    ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--minlen', type=int, default=MINL)
    ap.add_argument('--maxlen', type=int, default=MAXL)
    ap.add_argument('--last-two', '--son2', dest='son2', type=int, default=1, help="1 = require an exact match at the last two 3' bases")
    ap.add_argument('--tsv', default='', help='write the result to this file as TSV')
    a = ap.parse_args(argv)

    yollar = []
    for d in a.fastq:
        g = sorted(glob.glob(d))
        yollar.extend(g if g else [d])

    satirlar = []
    for p in yollar:
        if not os.path.exists(p):
            sys.stderr.write(u'MISSING: %s\n' % p); continue
        rs, n0 = kutu_yukle(p, a.nmax, a.seed, a.minlen, a.maxlen)
        pos, n, boylar = kutu_pcr(rs, a.F, a.R, a.lo, a.hi, a.mm, bool(a.son2))
        pct = 100.0 * pos / n if n else 0.0
        bb = ';'.join('%d:%d' % x for x in boylar.most_common(3))
        satirlar.append((os.path.basename(p), pos, n, n0, round(pct, 2), a.mm, bb))
        print('%-40s %6d/%-6d %6.2f%%  (mm<=%d, suzgecten gecen %d)  %s'
              % (os.path.basename(p), pos, n, pct, a.mm, n0, bb), flush=True)

    if a.tsv:
        with open(a.tsv, 'w', encoding='utf-8') as fh:
            fh.write('dosya\turun_veren\tkullanilan_okuma\tsuzgecten_gecen\tyuzde\tolcut_mm\tbaskin_boylar\n')
            for r in satirlar:
                fh.write('\t'.join(str(x) for x in r) + '\n')
        sys.stderr.write('TSV yazildi: %s\n' % a.tsv)
    return 0


if __name__ == '__main__':
    sys.exit(main())
