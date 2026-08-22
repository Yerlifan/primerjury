# -*- coding: utf-8 -*-
"""
rederive_membership.py
=======================
It re-derives the membership FROM THE MEASURED IDENTITY rather than the Kraken
label, and re-measures every pair in the panel with the corrected membership.

WHY IT IS NEEDED
----------------
A target's "discrimination ratio" depends directly on which bin counts as a
MEMBER and which as a COMPETITOR. The bin labels come from Kraken and, as far as
has been measured, at least 12 of them are WRONG. When a mislabelled bin is
written into the competitor column, the metric compares the target against the
target and even a perfect primer gives a discrimination below 1. A measured
example: Petriella 0.71x -> 8.47x (the same primer, only the membership corrected).

THE BASIC PRINCIPLE
-------------------
AN ABSENCE OF EVIDENCE IS NOT EVIDENCE. A bin changes place only on POSITIVE
measured evidence. If the measurement produced no signal for a bin, that bin's OLD
state is kept. The rule matters: done the other way round, the very bug being
diagnosed repeats in the opposite direction (which happened exactly once while
this script was being developed).

HOW IT WORKS
------------
1) The consensuses within a class are aligned at full length; those at >=99.0%
   identity form a PRE-GROUP.
2) A DISCRIMINATING k-mer set is extracted for each pre-group:
      the group's k-mers MINUS the k-mers of all the other groups
   This step is required. Without it the Trichoderma reads are assigned to
   Petriella over the conserved regions (18S, 5.8S). Measured: 70% misassignment.
3) Each read is assigned to a group by the normalised share of discriminating
   k-mers. The threshold is f=0.30, calibrated against an independent in-silico
   PCR (agreement: within a few points across a 0-72% range).
4) The membership is re-derived (under the principle above).
5) Every panel pair is measured at BOTH mm<=1 AND mm<=3, AT FULL DEPTH.
6) The old and new discrimination ratios are written side by side.

INTERRUPTION TOLERANT
---------------------
It writes to disk as each stage finishes (_ck_*.json). If it is interrupted the
same command continues where it stopped. No stage is recomputed from scratch.

IT DOES NOT WRITE TO THE PANEL
------------------------------
This script WRITES NOTHING into the panel xlsx or tsv files. It only reads, and
produces its own output.

USAGE
-----
    python3 rederive_membership.py [--root PROJECT_DIRECTORY] [--nmax 3000] [--reset]

"""

# -------------------------------------------------------------------------
# rederive_membership.py re-derives which bin is a MEMBER and which a COMPETITOR
# from the MEASURED identity rather than the Kraken label, and re-measures every
# panel pair with the corrected membership.
#
# INPUT  : konsensus_kanonik/*.kanonik.fa (the bin consensuses),
#          "fastq files"/*/*.fastq(.gz) (the raw reads),
#          screening/target_membership.tsv (the current membership definition),
#          primer_final/devir_ciftleri_20260802_sonrotus_TESLIM.tsv (the panel pairs).
# OUTPUT : engine_RESULT/engine_TURETME.md,
#          engine_RESULT/ciftler_yeniden_olcum.tsv,
#          engine_RESULT/kutu_olculen_kimlik.tsv,
#          engine_RESULT/_ck_*.json (interruption checkpoints).
#          It WRITES NOTHING into the panel files.
# CALLED BY: verification/full_chain.py -> key U
#          (python3 engine/rederive_membership.py --root .)
#
# THE TABLE THIS FILE PRODUCES IS THE INPUT OF STAGES P, K AND D: if the membership
# changes, every discrimination ratio in the panel changes. That is why it stands
# first in the chain.
# -------------------------------------------------------------------------
import os, sys, json, glob, random, argparse, time, csv, re

K = 21
NORM_ESIK = 0.30      # ayirt edici k-mer normalize esigi (kalibre edildi)
UYE_ESIK = 50.0       # a MEMBER when >=50% of its reads are in the target group
KARISIK_ESIK = 15.0   # %15-50 arasi KARISIK (ne uye ne rakip)
KONS_ESIK = 99.0      # the consensus identity threshold for pre-grouping
OKUMA_MIN, OKUMA_MAX = 200, 6000
TOHUM = 20260802
ENKOTU_ASGARI = 150
KAPSAM_ESIGI = 0.20
SINIF_ORNEK = 300     # reads per bin in the identity measurement

try:
    import numpy as np
except ImportError:
    sys.exit(u'ERROR: numpy is missing.  Inside WSL:  pip3 install numpy --break-system-packages')

_C = str.maketrans('ACGTURYSWKMBDHVNacgturyswkmbdhvn', 'TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn')
def rc(s): return s.translate(_C)[::-1]
def temizle(s): return re.sub(r'[^ACGTUN]', 'N', s.upper()).replace('U', 'T')
_M = np.full(256, 4, dtype=np.int8)
for _i, _c in enumerate('ACGT'): _M[ord(_c)] = _i
def enc(s): return _M[np.frombuffer(s.encode(), dtype=np.uint8)]

# Basit FASTA okuyucu. Konsensus dosyalari kucuktur, akis yeterlidir.
def fasta(p):
    h = None; b = []
    for line in open(p, 'rt', errors='ignore'):
        if line.startswith('>'):
            if h: yield h, ''.join(b)
            h = line[1:].strip(); b = []
        else: b.append(line.strip())
    if h: yield h, ''.join(b)

# ----------------------------------------------------------------- alignment
# -------------------------------------------------------------------------
# INFIX (HW) LEVENSHTEIN - it returns the last row.
# Why infix: the short query is aligned INSIDE the long target, so the overhang left
# at the start and the end of the target IS NOT PENALISED. Consensus lengths differ
# a great deal (1.5 kb beside 4.5 kb); a global alignment counts that difference as
# mismatch and would make the same organism look different.
#
# The left neighbour dependency (insertion) is vectorised: the relation
# now[j] = min(cand[j], now[j-1]+1) turns into a running minimum by setting
# a[j] = now[j]-j, and drops to one pass with np.minimum.accumulate. With the Python
# inner loop gone, a 1.5 kb x 1.5 kb alignment falls from minutes to seconds.
# -------------------------------------------------------------------------
def _hw_son(q, t):
    n = len(t); ar = np.arange(n + 1, dtype=np.int32)
    prev = np.zeros(n + 1, dtype=np.int32)
    for i in range(len(q)):
        qi = q[i]
        cost = np.where((t == qi) & (t < 4) & (qi < 4), 0, 1).astype(np.int32)
        cur = np.empty(n + 1, dtype=np.int32)
        cur[0] = prev[0] + 1
        cur[1:] = np.minimum(prev[:-1] + cost, prev[1:] + 1)
        cur = np.minimum.accumulate(cur - ar) + ar
        prev = cur
    return prev

# The percent identity. The denominator is ALWAYS the length of the short sequence;
# the long sequence's excess does not distort the ratio. The bin pre-grouping
# (KONS_ESIK = 99%) rests on this number.
def hw_kimlik(a, b):
    """kisa olani sorgu, uzun olanin icine; donus: yuzde kimlik"""
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    d = int(_hw_son(enc(q), enc(t)).min())
    return round(100.0 * (1 - d / max(len(q), 1)), 2)

# ----------------------------------------------------------------- in-silico PCR
# -------------------------------------------------------------------------
# IN-SILICO PCR ON A BIN'S RAW READS.
#
# The criterion: <=max_mm mismatches AND the last two bases at the 3' end matching
# EXACTLY. The last-two condition is not cosmetic; a polymerase does not start
# extending from a primer whose 3' end does not hold.
#
# PIGEONHOLE SEEDING - WHY IT IS LOSSLESS
# The primer is split into max_mm+1 NON-OVERLAPPING blocks. If the sequence holds at
# most max_mm mismatches, those mismatches can spread over at most max_mm separate
# blocks; AT LEAST ONE block is left over and that block MUST match exactly. So a
# seeding scheme that searches for an exact match of any block CANNOT MISS a single
# binding site satisfying the criterion.
#
# This IS NOT A HEURISTIC SPEED-UP, it is A GUARANTEE. The distinction matters: a
# heuristic seeder says "it will probably find it" and never reports what it missed,
# and the seed bug in the panel's old read engine produced exactly that kind of
# silent loss. Here the candidate set is NARROWED by the seed, and the decision is
# made on all candidates by a full mismatch count.
# -------------------------------------------------------------------------
class Kutu:
    """In-silico PCR on the raw reads. The criterion: <=max_mm mismatches plus the last
        2 bases at the 3' end EXACT.
        Pigeonhole seeding: the primer is split into max_mm+1 NON-OVERLAPPING blocks and
        at least one of them must match exactly -> LOSSLESS (for mm<=1 and mm<=3).

    """
    def __init__(self, path, nmax=3000):
        rs = []
        op = open
        with op(path, 'rt', errors='ignore') as fh:
            for k, line in enumerate(fh):
                if k % 4 == 1:
                    s = line.strip().upper()
                    if OKUMA_MIN <= len(s) <= OKUMA_MAX: rs.append(s)
        self.n_suzgec = len(rs)
        if nmax and len(rs) > nmax:
            random.seed(TOHUM); rs = random.sample(rs, nmax)
        self.n_okuma = len(rs)
        parts = []; sid = []
        for i, s in enumerate(rs):
            for ss in (s, rc(s)): parts.append(ss); sid.append((i, len(ss)))
        blob = 'N'.join(parts); self.E = enc(blob)
        off = []; p = 0
        for (rid, L) in sid: off.append((p, p + L, rid)); p += L + 1
        self.starts = np.array([o[0] for o in off], dtype=np.int64)
        self.ends = np.array([o[1] for o in off], dtype=np.int64)
        self.rid = np.array([o[2] for o in off], dtype=np.int32)
        self.idx = {}
        E = self.E.astype(np.int64); n = len(E)
        for kk in (9, 5):
            code = np.zeros(n - kk + 1, dtype=np.int64); bad = np.zeros(n - kk + 1, dtype=bool)
            for j in range(kk):
                seg = E[j:n - kk + 1 + j]
                code = code * 4 + np.where(seg < 4, seg, 0); bad |= seg >= 4
            code[bad] = -1
            ok = np.nonzero(code >= 0)[0]
            pos = ok[np.argsort(code[ok], kind='stable')]
            self.idx[kk] = (pos, code[pos])
            del code, bad
        del E
        self.okumalar = rs
    @staticmethod
    def _kod(s, k):
        c = 0
        for ch in s[:k]:
            v = 'ACGT'.find(ch)
            if v < 0: return -1
            c = c * 4 + v
        return c
    # -----------------------------------------------------------------------
    # Finds every binding site of an oligo (pigeonhole seeding).
    #
    # 1) The oligo is split into max_mm+1 non-overlapping blocks -> at least one must
    #    match exactly.
    # 2) A k is chosen no larger than the block length (9 for mm<=1, 5 for mm<=3); if k
    #    exceeded the block length the seed would spill outside the block and THE
    #    GUARANTEE WOULD BREAK, so k is always bounded by min(block length).
    # 3) Each block's k-mer code is found by binary search in a pre-built sorted index,
    #    and the union of the candidates is taken.
    # 4) The candidate positions are clipped to the read boundaries (the N separator
    #    placed where one read ends and the next begins must not produce a match across
    #    two reads).
    # 5) The decision: a FULL mismatch count <= max_mm AND both of the last two bases
    #    exact. The seed only COLLECTS CANDIDATES; the elimination happens here, on a
    #    full count.
    # -----------------------------------------------------------------------
    def _yerler(self, olig, max_mm):
        L = len(olig); nb = max_mm + 1
        kes = [round(i * L / nb) for i in range(nb + 1)]
        bloklar = [(kes[i], kes[i + 1] - kes[i]) for i in range(nb) if kes[i + 1] > kes[i]]
        k = min(min(b[1] for b in bloklar), 9 if max_mm <= 1 else 5)
        k = max(k, 4)
        if k not in self.idx: k = min(self.idx, key=lambda x: abs(x - k))
        k = min(k, min(b[1] for b in bloklar))
        if k not in self.idx: k = min(self.idx)
        pos, key = self.idx[k]
        cands = []
        for off, _ in bloklar:
            c = self._kod(olig[off:off + k], k)
            if c < 0: continue
            a = np.searchsorted(key, c, 'left'); b = np.searchsorted(key, c, 'right')
            if b > a: cands.append(pos[a:b] - off)
        if not cands: return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
        S = np.unique(np.concatenate(cands))
        si = np.searchsorted(self.starts, S, 'right') - 1
        ok = si >= 0; S = S[ok]; si = si[ok]
        ok = (S >= self.starts[si]) & (S + L <= self.ends[si]); S = S[ok]; si = si[ok]
        if S.size == 0: return S, si
        O = enc(olig)
        blok = self.E[S[:, None] + np.arange(L)[None, :]]
        mm = (blok != O[None, :]).sum(1)
        keep = (mm <= max_mm) & (blok[:, L - 1] == O[L - 1]) & (blok[:, L - 2] == O[L - 2])
        return S[keep], si[keep]
    # -----------------------------------------------------------------------
    # How many reads give a product? The forward primer is searched in the sense
    # direction and the REVERSE COMPLEMENT of the reverse primer in the same read; the
    # two must be in the SAME read, in the right order (the reverse primer after the
    # forward primer ends) and lo-hi bp apart.
    #
    # What is counted is THE NUMBER OF READS GIVING A PRODUCT, not the number of binding
    # sites: even if a read holds several valid pairs, that read is counted once (a
    # giving set). Otherwise repetitive regions would inflate the proportion.
    # -----------------------------------------------------------------------
    def pcr(self, F, R, lo=60, hi=400, max_mm=1):
        Fs, Fi = self._yerler(F, max_mm)
        if Fs.size == 0: return 0, self.n_okuma
        rr = rc(R); Rs, Ri = self._yerler(rr, max_mm)
        if Rs.size == 0: return 0, self.n_okuma
        from collections import defaultdict
        rmap = defaultdict(list)
        for p, i in zip(Rs.tolist(), Ri.tolist()): rmap[i].append(p)
        veren = set(); LF = len(F); LR = len(R)
        for p, i in zip(Fs.tolist(), Fi.tolist()):
            for q in rmap.get(i, ()):
                n = q + LR - p
                if lo <= n <= hi and q >= p + LF:
                    veren.add(int(self.rid[i])); break
        return len(veren), self.n_okuma

# -------------------------------------------------------------------------
# THE WILSON SCORE INTERVAL - WHY THE RAW PROPORTION IS NOT USED
#
# The raw proportion k/n misleads on a small sample: if 3 of 3 reads gave a product
# the raw proportion is 100%, but there is almost no evidence behind that number. In
# the same way, if 0 of 200 reads gave a product the raw proportion is 0% and gives
# the impression of "no cross reaction at all", when the real proportion could be
# 1.5%.
#
# The Wilson interval turns that uncertainty into a number, and THE CONSERVATIVE
# SIDE IS ALWAYS taken:
#   the member side     -> the LOWER bound (the lowest possible estimate of how well
#                          the target is seen)
#   the competitor side -> the UPPER bound (the highest possible estimate of the
#                          cross reaction risk)
# The discrimination ratio is the ratio of those two, so it always measures the
# worst case. Because of that choice, shallow bins give a LOW fold. That is not an
# error but a scarcity of evidence, and it is also why the measurement depth is kept
# the same on every row.
# -------------------------------------------------------------------------
def wilson(k, n, z=1.96):
    import math
    if n == 0: return (0.0, 1.0)
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d
    s = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - s), min(1.0, c + s))

# ----------------------------------------------------------------- helpers
# Writing a checkpoint is ATOMIC: it is written to a .tmp file first and then put in
# place with os.replace. If the run is interrupted while writing, no half JSON is
# left behind; a half checkpoint would silently break the next run.
def ck_yaz(yol, veri):
    tmp = yol + '.tmp'
    json.dump(veri, open(tmp, 'w', encoding='utf-8'))
    os.replace(tmp, yol)

def ck_oku(yol, varsayilan):
    if os.path.exists(yol):
        try: return json.load(open(yol, encoding='utf-8'))
        except Exception: pass
    return varsayilan

def sinifi(kutu): return kutu.split('_')[0].split('-')[0]

# A fastq file name -> a bin name. The files of sample A1-1 are named with
# underscores (A1_1_reads_2223.fastq); without normalising, those bins would not be
# recognised at all and would silently fall outside the measurement.
def kutu_adi(dosya):
    """The bin name from a fastq file name. Sample A1-1 uses underscore naming
        (A1_1_reads_2223.fastq); it is normalised, or those bins are NOT RECOGNISED.

    """
    b = os.path.basename(dosya)
    for uz in ('.fastq.gz', '.fastq'):
        if b.endswith(uz): b = b[:-len(uz)]
    b = b.replace('_reads_', '-reads_')
    m = re.match(r'^([A-Za-z0-9]+)[-_](\d+)-reads_(\d+)$', b)
    if m: return '%s-%s_%s' % (m.group(1), m.group(2), m.group(3))
    m = re.match(r'^(.+)-reads_(\d+)$', b)
    if m: return '%s_%s' % (m.group(1).replace('_', '-'), m.group(2))
    return b

# ----------------------------------------------------------------- main flow
# -------------------------------------------------------------------------
# THE MAIN FLOW - five steps, all of them checkpointed:
#   0) inventory          : match up the consensus and fastq files.
#   1) the identity matrix: align every consensus pair within a class (_ck_kimlik).
#   2) classification     : pre-groups plus DISCRIMINATING k-mers plus read
#                           assignment (_ck_icerik).
#   3) membership         : rebuild the member / mixed / competitor sets.
#   4) re-measurement     : every panel pair at BOTH mm<=1 AND mm<=3 (_ck_olcum).
# -------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', dest='kok', default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument('--nmax', type=int, default=3000)
    ap.add_argument('--reset', dest='sifirla', action='store_true', help='delete checkpoints and start over')
    a = ap.parse_args()
    KOK = a.kok
    CIK = os.path.join(KOK, 'engine_RESULT')
    os.makedirs(CIK, exist_ok=True)
    if a.sifirla:
        for f in glob.glob(os.path.join(CIK, '_ck_*.json')): os.remove(f)
        print(u'checkpoints deleted.')

    KONS = os.path.join(KOK, 'konsensus_kanonik')
    FQ = os.path.join(KOK, 'fastq files')
    for p, adi in ((KONS, 'konsensus_kanonik'), (FQ, 'fastq files')):
        if not os.path.isdir(p): sys.exit(u'ERROR: there is no %s directory: %s' % (adi, p))

    print('=' * 70); print(u'  RE-DERIVING MEMBERSHIP FROM THE MEASURED IDENTITY'); print('=' * 70)
    print('  Proje  : %s' % KOK); print(u'  Output : %s' % CIK); print()

    # --- 0. envanter
    kons = {}
    for p in glob.glob(os.path.join(KONS, '*.kanonik.fa')):
        kons[os.path.basename(p).replace('.kanonik.fa', '')] = temizle(list(dict(fasta(p)).values())[0])
    fq = {}
    for p in glob.glob(os.path.join(FQ, '*', '*.fastq')) + glob.glob(os.path.join(FQ, '*', '*.fastq.gz')):
        fq[kutu_adi(p)] = p
    ortak = sorted(set(kons) & set(fq))
    print(u'  consensus %d, fastq %d, matched bins %d' % (len(kons), len(fq), len(ortak)))
    eksik = sorted(set(fq) - set(kons))
    if eksik: print(u'  WARNING - a fastq with no consensus: %s' % ', '.join(eksik))
    print()

    # STEP 1 - the identity matrix of every consensus pair within a class. Only bins in
    # THE SAME class are compared: different classes are already different loci and the
    # identity between them means nothing. The matrix is written to a checkpoint, so an
    # interrupted run does not recompute it from scratch.
    # --- 1. the identity matrix
    ckm = os.path.join(CIK, '_ck_kimlik.json')
    KM = ck_oku(ckm, {})
    siniflar = {}
    for k in kons: siniflar.setdefault(sinifi(k), []).append(k)
    ciftler = []
    for s, v in siniflar.items():
        v = sorted(v)
        for i in range(len(v)):
            for j in range(i + 1, len(v)): ciftler.append((v[i], v[j]))
    yeni = [c for c in ciftler if '%s|%s' % c not in KM]
    if yeni:
        print(u'  [1/4] Consensus identity matrix: %d pairs to compute' % len(yeni))
        t0 = time.time()
        for n, (x, y) in enumerate(yeni, 1):
            KM['%s|%s' % (x, y)] = hw_kimlik(kons[x], kons[y])
            if n % 25 == 0 or n == len(yeni):
                ck_yaz(ckm, KM)
                print('        %d/%d  (%.0f sn)' % (n, len(yeni), time.time() - t0), flush=True)
        ck_yaz(ckm, KM)
    else:
        print(u'  [1/4] Identity matrix read from checkpoint (%d pairs)' % len(KM))

    # -----------------------------------------------------------------------
    # STEP 2 - THE DISCRIMINATING k-mer set. The k-mers of ALL THE OTHER groups are
    # SUBTRACTED from a pre-group's k-mers; what is left is only what is SPECIFIC to
    # that group.
    #
    # That subtraction IS REQUIRED. Without it the reads match over conserved regions
    # (18S, 5.8S, the LSU core) and the groups run into one another. Measured: 70% of
    # the Trichoderma reads had been assigned to Petriella. A conserved region is THE
    # SAME in every group so it is useless for discrimination, and falsely high
    # similarity comes from exactly there.
    #
    # The assignment threshold NORM_ESIK = 0.30 is a normalised share (divided by the
    # smaller of the group's discriminating set size and the read's k-mer count), and it
    # was calibrated against an independent in-silico PCR.
    # -----------------------------------------------------------------------
    # --- 2. pre-groups plus discriminating k-mers plus read classification
    cki = os.path.join(CIK, '_ck_icerik.json')
    IC = ck_oku(cki, {})
    GRUPLAR = {}
    for s, probs in sorted(siniflar.items()):
        probs = sorted(probs)
        par = {k: k for k in probs}
        def bul(x):
            while par[x] != x: par[x] = par[par[x]]; x = par[x]
            return x
        for (x, y), v in ((tuple(k.split('|')), v) for k, v in KM.items()):
            if x in par and y in par and v >= KONS_ESIK:
                rx, ry = bul(x), bul(y)
                if rx != ry: par[rx] = ry
        g = {}
        for k in probs: g.setdefault(bul(k), []).append(k)
        GRUPLAR[s] = sorted([sorted(v) for v in g.values()], key=lambda v: v[0])
    if any(k not in IC for k in ortak):
        print(u'  [2/4] Classifying reads with discriminating k-mers')
        for s, gruplar in sorted(GRUPLAR.items()):
            hedefler = [k for k in ortak if sinifi(k) == s and k not in IC]
            if not hedefler: continue
            ks = []
            for grp in gruplar:
                st = set()
                for k in grp:
                    for ss in (kons[k], rc(kons[k])):
                        for i in range(len(ss) - K + 1): st.add(ss[i:i + K])
                ks.append(st)
            ayirt = []
            for i, st in enumerate(ks):
                dis = set()
                for j, s2 in enumerate(ks):
                    if i != j: dis |= s2
                ayirt.append(st - dis)
            boy = [len(x) for x in ayirt]
            idx = {}
            for i, st in enumerate(ayirt):
                for x in st: idx[x] = i
            print(u'        class %-3s: %d groups, discriminating set %s' % (s, len(gruplar), boy))
            for kb in hedefler:
                rs = []
                with open(fq[kb], 'rt', errors='ignore') as fh:
                    for n, line in enumerate(fh):
                        if n % 4 == 1:
                            x = line.strip().upper()
                            if OKUMA_MIN <= len(x) <= OKUMA_MAX: rs.append(x)
                random.seed(TOHUM)
                if len(rs) > SINIF_ORNEK: rs = random.sample(rs, SINIF_ORNEK)
                say = [0] * len(gruplar)
                for r in rs:
                    nk = len(r) - K + 1
                    if nk <= 0: continue
                    c = [0] * len(gruplar)
                    for i in range(nk):
                        j = idx.get(r[i:i + K])
                        if j is not None: c[j] += 1
                    for gi in range(len(gruplar)):
                        if c[gi] / max(1, min(boy[gi], nk)) >= NORM_ESIK: say[gi] += 1
                IC[kb] = dict(n=len(rs), sinif=s,
                              pay={gruplar[i][0]: round(100.0 * say[i] / max(len(rs), 1), 1)
                                   for i in range(len(gruplar)) if say[i]})
                ck_yaz(cki, IC)
                print('          %-18s n=%4d  %s' % (kb, len(rs),
                      ', '.join('%s:%.0f%%' % x for x in sorted(IC[kb]['pay'].items(), key=lambda y: -y[1])[:2]) or '-'),
                      flush=True)
    else:
        print(u'  [2/4] Read classification read from checkpoint')

    # -----------------------------------------------------------------------
    # STEP 3 - REBUILDING THE MEMBERSHIP.
    # First the dominant groups of the current member bins are collected (hg = the
    # target groups), then ALL the bins of the class are redistributed against those
    # groups:
    #   share >= UYE_ESIK       -> MEMBER
    #   KARISIK_ESIK .. UYE     -> MIXED (neither a member nor a competitor)
    #   below that              -> COMPETITOR
    # On universal targets the whole class is a member; no separation is made.
    # -----------------------------------------------------------------------
    # --- 3. re-deriving the membership
    print(u'  [3/4] Re-deriving membership')
    uyelik_tsv = os.path.join(KOK, 'screening', 'target_membership.tsv')
    panel_tsv = os.path.join(KOK, 'primer_final', 'devir_ciftleri_20260802_sonrotus_TESLIM.tsv')
    if not os.path.exists(uyelik_tsv): sys.exit(u'ERROR: %s is missing' % uyelik_tsv)
    if not os.path.exists(panel_tsv): sys.exit(u'ERROR: %s is missing' % panel_tsv)
    uyelik = {}
    for line in open(uyelik_tsv, encoding='utf-8'):
        if line.startswith('#') or line.startswith('hedef\t'): continue
        p = line.rstrip('\n').split('\t')
        if len(p) >= 2 and p[0].strip(): uyelik[p[0]] = p[1]
    rows = [r for r in csv.DictReader(open(panel_tsv, encoding='utf-8'), delimiter='\t')
            if (r.get('Hedef') or '').strip() in uyelik]
    UY = {}
    for r in rows:
        hd = r['Hedef'].strip()
        sn = [x.strip() for x in (r.get('Amplikon sinifi') or '').split('/') if x.strip()] or list(siniflar)
        ut = uyelik[hd]; evrensel = ut.strip().startswith('*')
        tum = [k for k in ortak if sinifi(k) in sn]
        if evrensel:
            UY[hd] = dict(sinif=sn, eski=tum, yeni=tum, karisik=[], rakip=[], eklenen=[], cikan=[], evrensel=True)
            continue
        tx = set(x.strip() for x in ut.split(',') if x.strip())
        eski = [k for k in tum if k.split('_')[1] in tx]
        hg = set()
        for k in eski:
            for g, v in (IC.get(k, {}) or {}).get('pay', {}).items():
                if v >= UYE_ESIK: hg.add(g)
        yeni_u = []; kar = []; rak = []
        for k in tum:
            pay = (IC.get(k, {}) or {}).get('pay', {})
            ic = max([pay.get(g, 0.0) for g in hg], default=0.0)
            en = max(pay.items(), key=lambda x: x[1]) if pay else (None, 0.0)
            if k in eski:
                # ---------------------------------------------------------------
                # AN ABSENCE OF EVIDENCE IS NOT EVIDENCE.
                # A bin that is already a member changes place only on POSITIVE measured
                # evidence: the dominant group of its reads must pass UYE_ESIK AND that
                # group must NOT BE one of the target groups. If the measurement produced no
                # signal for this bin (the share table is empty or weak) the bin STAYS A
                # MEMBER.
                #
                # Had it been done the other way round - "remove it when there is no
                # evidence" - the very bug being diagnosed would repeat in the opposite
                # direction: bins with few reads would fall into the competitor column just
                # for being quiet, and would distort the discrimination ratios.
                # ---------------------------------------------------------------
                # with no POSITIVE evidence the old state is kept
                if en[1] >= UYE_ESIK and en[0] not in hg: rak.append(k)
                else: yeni_u.append(k)
            else:
                if ic >= UYE_ESIK: yeni_u.append(k)
                elif ic >= KARISIK_ESIK: kar.append(k)
                else: rak.append(k)
        UY[hd] = dict(sinif=sn, eski=sorted(eski), yeni=sorted(yeni_u), karisik=sorted(kar),
                      rakip=sorted(rak), eklenen=sorted(set(yeni_u) - set(eski)),
                      cikan=sorted(set(eski) - set(yeni_u)), evrensel=False)
    deg = sum(1 for o in UY.values() if o['eklenen'] or o['cikan'])
    print(u'        the membership of %d targets out of %d changed' % (len(UY), deg))

    # STEP 4 - every panel pair, in every bin, measured at BOTH mm<=1 (the primary
    # criterion) AND mm<=3 (the robustness criterion). Each bin's fastq is read once and
    # ALL the pairs in that bin are asked against the same index; a checkpoint is written
    # as each bin finishes.
    # --- 4. re-measure the panel pairs
    CF = []
    for r in rows:
        F = (r.get("Ileri primer (5'->3')") or '').strip().upper()
        R = (r.get("Geri primer (5'->3')") or '').strip().upper()
        if not re.fullmatch(r'[ACGT]+', F or '') or not re.fullmatch(r'[ACGT]+', R or ''): continue
        CF.append(dict(hedef=r['Hedef'].strip(), F=F, R=R,
                       urun=(r.get('Urun (bp)') or '').strip()))
    cko = os.path.join(CIK, '_ck_olcum.json')
    # THE 2026-08-10 SEQUENCE SEAL. The cache was keyed by bin name and the pair results
    # were held under the pair's ORDINAL (str(i)). Two separate bugs:
    #   1) when a pair's sequence changed, the same ordinal was read and the OLD
    #      measurement was written beside the new sequence (a new sequence with an old
    #      number);
    #   2) when a pair was added to or removed from the panel the ordinals shifted and
    #      numbers could be attached to THE WRONG pair.
    # The fix: the pair key is not the ordinal but a digest of the F+R sequence. A seal
    # is also written at the head of the file; if the seal does not match, the cache is
    # rebuilt from scratch.
    import hashlib as _hl

    def _ck_anahtar(c):
        return _hl.md5((c['F'] + '|' + c['R']).encode('utf-8')).hexdigest()[:12]

    _muhur = _hl.md5('|'.join(sorted(_ck_anahtar(c) for c in CF))
                     .encode('utf-8')).hexdigest()[:12]
    OL = ck_oku(cko, {})
    if OL.get('_muhur') != _muhur:
        if OL:
            print(u'  [4/4] the cache SEQUENCE seal does not hold (recorded %s, now %s), measuring from scratch.' % (OL.get('_muhur') or 'yok', _muhur))
        OL = {'_muhur': _muhur}
    kalan = [k for k in ortak if k not in OL]
    print(u'  [4/4] In-silico PCR: %d pairs x %d bins (%d bins remaining), mm<=1 and mm<=3' % (len(CF), len(ortak), len(kalan)))
    t0 = time.time()
    for n, kb in enumerate(kalan, 1):
        Kt = Kutu(fq[kb], nmax=a.nmax)
        r = {}
        for i, c in enumerate(CF):
            p1, nn = Kt.pcr(c['F'], c['R'], 60, 400, 1)
            p3, _ = Kt.pcr(c['F'], c['R'], 60, 400, 3)
            r[_ck_anahtar(c)] = [p1, nn, p3]
        OL[kb] = r; del Kt
        ck_yaz(cko, OL)
        print('        %-18s %d/%d  (%.0f sn)' % (kb, n, len(kalan), time.time() - t0), flush=True)

    # --- the discrimination ratios
    # -----------------------------------------------------------------------
    # THE DISCRIMINATION RATIO = (the LOWEST Wilson LOWER bound of the member bins) /
    #                            (the HIGHEST Wilson UPPER bound of the competitor bins)
    # The worst bin is taken on both sides: in the numerator the member that sees the
    # target least, in the denominator the competitor that cross reacts most. Had an
    # average been used, a single bad bin would dissolve in the crowd and the real risk
    # would be invisible.
    #
    # Competitor bins below ENKOTU_ASGARI (150 reads) do not enter the denominator: at
    # that depth the Wilson upper bound almost always hits the ceiling, and the fold
    # would measure not the competitor's real behaviour but only the scarcity of reads.
    #
    # COVERAGE is a separate axis: how many of the member bins give >=20% product. If
    # the discrimination is high but the coverage low, the pair is specific but does not
    # see the whole target. The two problems must not be confused.
    # -----------------------------------------------------------------------
    def hesap(ci, uye, rakip, mm):
        uy = []; rk = []
        for kb, r in OL.items():
            if kb == '_muhur' or not isinstance(r, dict):
                continue
            v = r.get(_ck_anahtar(CF[ci]))
            if not v: continue
            p, nn = (v[0], v[1]) if mm == 1 else (v[2], v[1])
            if kb in uye: uy.append((kb, p, nn))
            elif kb in rakip: rk.append((kb, p, nn))
        if not uy: return None
        ua = min(wilson(p, nn)[0] for _, p, nn in uy)
        kaps = sum(1 for _, p, nn in uy if nn and p / nn >= KAPSAM_ESIGI)
        enk = None
        for kb, p, nn in rk:
            if nn < ENKOTU_ASGARI: continue
            hi = wilson(p, nn)[1]
            if enk is None or hi > enk[1]: enk = (kb, hi, p, nn)
        return dict(uye_alt=round(100 * ua, 2), kapsam='%d/%d' % (kaps, len(uy)),
                    enkotu=enk[0] if enk else '-', enkotu_ust=round(100 * enk[1], 2) if enk else 0.0,
                    kat=round(ua / enk[1], 2) if enk and enk[1] > 0 else None)
    SON = []
    for ci, c in enumerate(CF):
        o = UY.get(c['hedef'])
        if not o: continue
        tum = [k for k in OL if sinifi(k) in o['sinif']]
        eu = set(o['eski']); er = set(tum) - eu
        yu = set(o['yeni']); kar = set(o['karisik'])
        rA = set(tum) - yu; rB = set(tum) - yu - kar
        d = dict(hedef=c['hedef'], sinif='/'.join(o['sinif']), F=c['F'], R=c['R'], urun=c['urun'],
                 eski_n=len(eu), yeni_n=len(yu), kar_n=len(kar),
                 eklenen=';'.join(o['eklenen']), cikan=';'.join(o['cikan']))
        for mm in (1, 3):
            d['eski_mm%d' % mm] = hesap(ci, eu, er, mm)
            d['A_mm%d' % mm] = hesap(ci, yu, rA, mm)
            d['B_mm%d' % mm] = hesap(ci, yu, rB, mm)
        SON.append(d)

    # --- ciktilar
    def kat(x): return '' if not x or x['kat'] is None else x['kat']
    p1 = os.path.join(CIK, 'ciftler_yeniden_olcum.tsv')
    with open(p1, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'sinif', 'F', 'R', 'urun_bp', 'eski_uye_n', 'yeni_uye_n', 'karisik_n',
                    'eklenen', 'cikan', 'eski_kat_mm1', 'yeniA_kat_mm1', 'yeniB_kat_mm1',
                    'eski_kat_mm3', 'yeniA_kat_mm3', 'yeniB_kat_mm3', 'yeniA_uye_alt',
                    'yeniA_kapsam', 'yeniA_enkotu_kutu', 'yeniA_enkotu_ust'])
        for d in SON:
            A = d['A_mm1'] or {}
            w.writerow([d['hedef'], d['sinif'], d['F'], d['R'], d['urun'], d['eski_n'], d['yeni_n'],
                        d['kar_n'], d['eklenen'], d['cikan'],
                        kat(d['eski_mm1']), kat(d['A_mm1']), kat(d['B_mm1']),
                        kat(d['eski_mm3']), kat(d['A_mm3']), kat(d['B_mm3']),
                        A.get('uye_alt', ''), A.get('kapsam', ''), A.get('enkotu', ''), A.get('enkotu_ust', '')])
    p2 = os.path.join(CIK, 'kutu_olculen_kimlik.tsv')
    with open(p2, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['kutu', 'sinif', 'kraken_taxid', 'baskin_grup', 'baskin_pay', 'ikinci_grup', 'ikinci_pay', 'yorum'])
        for kb in sorted(IC):
            pay = IC[kb]['pay']; s = IC[kb]['sinif']; t = kb.split('_')[1]
            srt = sorted(pay.items(), key=lambda x: -x[1])
            b1 = srt[0] if srt else ('', 0); b2 = srt[1] if len(srt) > 1 else ('', 0)
            yorum = ''
            if b1[0] and b1[1] >= UYE_ESIK and b1[0].split('_')[1] != t:
                yorum = 'KRAKEN ETIKETI YANLIS -> %s' % b1[0]
            elif not srt:
                yorum = 'olcum sinyali yok - eski durum korundu'
            w.writerow([kb, s, t, b1[0], b1[1], b2[0], b2[1], yorum])
    p3 = os.path.join(CIK, 'engine_TURETME.md')
    with open(p3, 'w', encoding='utf-8') as fh:
        fh.write(u'# Re-deriving membership from measured identity\n\n')
        fh.write('Uretim: %s\n\n' % time.strftime('%Y-%m-%d %H:%M:%S'))
        fh.write(u'Criterion: discriminating %d-mer, normalised threshold %.2f, member >=%%%.0f, mixed %%%.0f-%.0f.\n'
                 % (K, NORM_ESIK, UYE_ESIK, KARISIK_ESIK, UYE_ESIK))
        fh.write(u'In-silico PCR: <=1 and <=3 mismatches with an EXACT match at the last two 3\' bases, at most %d reads per bin.\n\n' % a.nmax)
        fh.write(u'Scenario A = mixed bins count as COMPETITORS (the safe one).  B = mixed bins are excluded.\n\n')
        fh.write(u'| target | old members | new members | mixed | old fold (mm1) | new A | new B | new A (mm3) |\n')
        fh.write('|---|---|---|---|---|---|---|---|\n')
        for d in SON:
            fh.write('| %s | %d | %d | %d | %s | %s | %s | %s |\n' % (
                d['hedef'], d['eski_n'], d['yeni_n'], d['kar_n'],
                kat(d['eski_mm1']), kat(d['A_mm1']), kat(d['B_mm1']), kat(d['A_mm3'])))
        fh.write(u'\n## The targets whose membership changed\n\n')
        for d in SON:
            if d['eklenen'] or d['cikan']:
                fh.write('- **%s**: eklenen `%s` / cikan `%s`\n' % (d['hedef'], d['eklenen'] or '-', d['cikan'] or '-'))
        fh.write(u'\n## The bins whose Kraken label is wrong\n\n')
        for kb in sorted(IC):
            pay = IC[kb]['pay']; t = kb.split('_')[1]
            srt = sorted(pay.items(), key=lambda x: -x[1])
            if srt and srt[0][1] >= UYE_ESIK and srt[0][0].split('_')[1] != t:
                fh.write(u'- `%s` -> %%%.0f of its reads belong to the organism `%s`\n' % (kb, srt[0][1], srt[0][0]))
        fh.write(u'\n> This script does NOT write to the panel files. Applying the changes is a separate job.\n')
    print()
    print('=' * 70); print(u'  DONE'); print('=' * 70)
    print('  %s' % p3); print('  %s' % p1); print('  %s' % p2)
    return 0

if __name__ == '__main__':
    sys.exit(main())
