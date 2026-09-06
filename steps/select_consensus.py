# -*- coding: utf-8 -*-
"""WHICH CONSENSUS TO USE: ask the bin's OWN READS, and say when none is trustworthy.

WHY THIS EXISTS
---------------
A bin can have several consensuses (dominant allele, anchored reference, self
consensus, canonical, regenerated). Which one to use was decided twice by the
wrong criterion before this step existed:

1. A FIXED ORDER. "If the dominant one is empty, take the next fallback." The
   next one won even when a better one sat further down the list.
2. FEWEST N. "The least ambiguous wins." Measured on one bacterial bin:
       dominant           1471 bp, N 0.5 %  -> Oceanobacillus chironomi 99.18 %
       anchored reference 1471 bp, N 0.0 %  -> Oceanobacillus chironomi 85.18 %
   The sequence with zero N was fourteen points FURTHER from the reference. N
   says how uncertain a consensus is, not how right it is.

THE CRITERION
-------------
A consensus exists to represent the bin's reads, so the reads are the judge:
what fraction of the consensus's k-mers occur in a sample of the bin's own
reads? This is independent of any reference database; "pick the one closest to
the reference" would be circular. The denominator is the CONSENSUS's k-mer
count, not the read's, on purpose: with the read as denominator a short
consensus is penalised structurally (measured: 4.2 kb reads, 1.47 kb
consensus, ceiling about 35 per cent). Candidates within SUPPORT_BAND (one
point) of the best support are treated as equal, and among them the one with
the fewest ambiguous bases wins, then the longest: the least ambiguous of the
consensuses the reads actually support.

THE FLOOR (2026-09-04)
----------------------
When the chosen consensus has read support below MIN_SUPPORT (60 per cent) the
row says "trusted = NO" out loud. Two bins had been "chosen" at 2.2 and 3.0 per
cent support: the dominant-allele step had aligned 0 of 3,001 reads to its
IUPAC template (WORK_RECORD 13.1). The best candidate is still chosen; the
floor makes such a failure visible instead of hiding it in a table.

RUN
---
    python3 steps/select_consensus.py --root . [--reads 25]
        [--sets referans_konsensus/baskin/konsensus referans_konsensus/konsensus ...]
OUTPUT
    <root>/CONSENSUS_SELECTION.tsv
"""
from __future__ import print_function
import argparse
import glob
import io
import os
import random
import sys

K = 21
MIN_SUPPORT = 60.0
# Candidates within SUPPORT_BAND points of the best read support count as equal,
# and among them the one with the FEWEST N wins, then the longest. Read support
# still comes first: a zero-N sequence at 85 per cent cannot beat one at 99.
SUPPORT_BAND = 1.0
# 'fungal_polish/consensus' (2026-09-05) is the consensus polished from a fungal bin's
# dominant population (verification/fungal_bin_identity.py). It is weighed with the
# SAME read-support criterion as every other candidate; no fungal special case.
DEFAULT_SETS = ('referans_konsensus/baskin/konsensus', 'referans_konsensus/pak_polish/consensus',
                'referans_konsensus/fungal_polish/consensus',
                'referans_konsensus/konsensus',
                'referans_konsensus/self/konsensus', 'consensus sequences')
COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def rc(s):
    return u''.join(COMPLEMENT.get(c, 'N') for c in reversed(s.upper()))


def kmers(s):
    s = s.upper()
    return set(s[i:i + K] for i in range(len(s) - K + 1) if 'N' not in s[i:i + K])


def read_sample(path, how_many, seed=20260827):
    """The first reads of at least 300 bp, then a fixed-seed sample of `how_many`."""
    out = []
    if not (path and os.path.exists(path)):
        return out
    with io.open(path, encoding='utf-8', errors='replace') as fh:
        while True:
            if not fh.readline():
                break
            seq = fh.readline().strip()
            fh.readline()
            fh.readline()
            if len(seq) >= 300:
                out.append(seq)
            if len(out) >= how_many * 8:
                break
    random.Random(seed).shuffle(out)
    return out[:how_many]


def bins(root):
    """(label, reads file) for every <root>/fastq files/<lib>/<lib>[-_]reads[-_]<id>.fastq."""
    out = []
    for path in sorted(glob.glob(os.path.join(root, 'fastq files', '*', '*reads*.fastq'))):
        lib = os.path.basename(os.path.dirname(path))
        base = os.path.basename(path)[:-6]
        ident = base.split('reads')[-1].lstrip('-_')
        if ident:
            out.append((u'%s_%s' % (lib, ident), path))
    return out


def candidates(root, sets, label):
    out = []
    for sub in sets:
        for f in sorted(glob.glob(os.path.join(root, sub, label + '_*'))
                        + glob.glob(os.path.join(root, sub, label + '.*'))):
            seq = u''.join(x.strip() for x in io.open(f, encoding='utf-8', errors='replace')
                           if not x.startswith(u'>')).upper()
            if seq:
                out.append((sub, seq))
    return out


def choose(cands, read_kmers, max_n=50.0):
    """(set, sequence, support %, N %) of the best candidate, or None."""
    measured = []
    for sub, seq in cands:
        n_pct = 100.0 * seq.count('N') / len(seq)
        if n_pct > max_n:
            continue
        kk = kmers(seq)
        support = (100.0 * len(kk & read_kmers) / len(kk)) if kk else 0.0
        measured.append((sub, seq, support, n_pct))
    if not measured:
        return None
    top = max(x[2] for x in measured)
    band = [x for x in measured if x[2] >= top - SUPPORT_BAND]
    return min(band, key=lambda x: (x[3], -len(x[1]), -x[2], x[0]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--root', default='.')
    ap.add_argument('--reads', type=int, default=25, help='reads sampled per bin')
    ap.add_argument('--sets', nargs='*', default=list(DEFAULT_SETS),
                    help='consensus directories, relative to --root')
    ap.add_argument('--out', default='CONSENSUS_SELECTION.tsv')
    a = ap.parse_args()
    root = os.path.abspath(a.root)
    rows = [u'\t'.join([u'bin', u'candidates', u'CHOSEN set', u'support %', u'N %', u'bp',
                        u'fewest-N set', u'its support %', u'criteria disagree', u'trusted'])]
    disagree = untrusted = 0
    for label, fq in bins(root):
        cands = candidates(root, a.sets, label)
        if len(cands) < 2:
            continue
        sample = read_sample(fq, a.reads)
        if len(sample) < 3:
            continue
        rk = set()
        for seq in sample:
            rk |= kmers(seq) | kmers(rc(seq))
        best = choose(cands, rk)
        if not best:
            continue
        sub, seq, support, n_pct = best
        fewest = min(cands, key=lambda x: (100.0 * x[1].count('N') / len(x[1]), -len(x[1])))
        kk = kmers(fewest[1])
        s2 = (100.0 * len(kk & rk) / len(kk)) if kk else 0.0
        diff = fewest[0] != sub
        if diff:
            disagree += 1
            print(u'  %-14s read support chose %-36s (%.1f%%) | fewest N: %s (%.1f%%)'
                  % (label, sub, support, fewest[0], s2))
        if support >= MIN_SUPPORT:
            trusted = u'yes'
        else:
            trusted = u'NO (support %.1f%% < %.0f%%)' % (support, MIN_SUPPORT)
            untrusted += 1
            print(u'  %-14s UNTRUSTED: the best candidate %s has only %.1f%% read support'
                  % (label, sub, support))
        rows.append(u'\t'.join([label, u'%d' % len(cands), sub, u'%.1f' % support,
                                u'%.1f' % n_pct, u'%d' % len(seq), fewest[0], u'%.1f' % s2,
                                u'YES' if diff else u'no', trusted]))
    out = os.path.join(root, a.out)
    io.open(out, 'w', encoding='utf-8', newline='\n').write(u'\n'.join(rows) + u'\n')
    print(u'\n  bins with more than one candidate : %d' % (len(rows) - 1))
    print(u'  the two criteria DISAGREE          : %d' % disagree)
    print(u'  UNTRUSTED (support < %.0f%%)        : %d' % (MIN_SUPPORT, untrusted))
    print(u'  written: %s' % out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
