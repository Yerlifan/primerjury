# -*- coding: utf-8 -*-
"""CHOOSE WHICH BINS ENTER THE IDENTITY CHAIN, largest first, with a floor.

WHY
---
bin_reads.py writes every bin with at least the read floor; one bacterial
barcode gave 1,962 of them (measured). The chain builds a consensus and scans
twelve reference databases per bin, so 99 bins took a night. Thousands of
small bins are mostly chimeras, fragments and single-error reads; scanning
them all would take days and add nothing. The bins that are not chosen are
NOT deleted: they stay where they are and SELECTION_TABLE.tsv records why.

THE CRITERION, kept in one place (here)
---------------------------------------
* the bin's share of the barcode's reads >= MIN_SHARE (0.5 per cent)   OR
* the bin is among the TOP_PER_WINDOW (5) largest of its length window,
  so that small windows are represented as well
* at most MAX_PER_BARCODE (40) bins per barcode
* every bin already carries >= the read floor (bin_reads.py applied it)

RUN
---
    python3 steps/select_bins.py --source "fastq files_bins" --target "fastq files" [--name A2-1 ...]

The chosen bins are COPIED into --target in the layout the chain reads.
"""
from __future__ import print_function
import argparse
import io
import os
import shutil
import sys

MIN_SHARE = 0.5        # per cent of the barcode's reads
TOP_PER_WINDOW = 5     # the largest N bins of every window are always taken
MAX_PER_BARCODE = 40


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--source', required=True, help='bin root written by bin_reads.py')
    ap.add_argument('--target', required=True, help='where the chosen bins are copied ("fastq files")')
    ap.add_argument('--name', nargs='*', default=None, help='only these barcodes (default: all with a BIN_TABLE.tsv)')
    a = ap.parse_args()
    src, dst = os.path.abspath(a.source), os.path.abspath(a.target)
    if not os.path.isdir(dst):
        os.makedirs(dst)
    names = a.name or sorted(d for d in os.listdir(src)
                             if os.path.isfile(os.path.join(src, d, 'BIN_TABLE.tsv')))
    rows = [u'\t'.join([u'barcode', u'bin', u'window', u'reads', u'share %', u'chosen', u'reason'])]
    total = 0
    for name in names:
        t = [s.rstrip(u'\n').split(u'\t') for s in
             io.open(os.path.join(src, name, 'BIN_TABLE.tsv'), encoding='utf-8') if s.strip()]
        head = t[0]
        ir, ib, ish, iw, ifile = (head.index(u'reads'), head.index(u'bin'),
                                  head.index(u'share % (barcode)'), head.index(u'window bp'),
                                  head.index(u'file'))
        bins = [r for r in t[1:] if r[ib].startswith(u'BIN')]
        rank = {}
        for r in sorted(bins, key=lambda r: -int(r[ir])):
            rank.setdefault(r[iw], []).append(r[ib])
        chosen = []
        for r in sorted(bins, key=lambda r: -int(r[ir])):
            share = float(r[ish])
            pos = rank[r[iw]].index(r[ib]) + 1
            if len(chosen) >= MAX_PER_BARCODE:
                reason, take = u'barcode ceiling %d' % MAX_PER_BARCODE, False
            elif share >= MIN_SHARE:
                reason, take = u'share >= %.1f%%' % MIN_SHARE, True
            elif pos <= TOP_PER_WINDOW:
                reason, take = u'among the %d largest of its window (rank %d)' % (TOP_PER_WINDOW, pos), True
            else:
                reason, take = u'share %.2f%% < %.1f%% and rank %d in its window' % (share, MIN_SHARE, pos), False
            if take:
                chosen.append(r)
            rows.append(u'\t'.join([name, r[ib], r[iw], r[ir], r[ish], u'YES' if take else u'no', reason]))
        d = os.path.join(dst, name)
        if not os.path.isdir(d):
            os.makedirs(d)
        for r in chosen:
            s, t2 = os.path.join(src, name, r[ifile]), os.path.join(d, r[ifile])
            if not os.path.exists(t2) or os.path.getsize(t2) != os.path.getsize(s):
                shutil.copyfile(s, t2)
        total += len(chosen)
        print(u'  %-6s bins %4d -> chosen %3d (windows: %s)'
              % (name, len(bins), len(chosen), u', '.join(sorted(rank))))
    io.open(os.path.join(dst, 'SELECTION_TABLE.tsv'), 'w', encoding='utf-8',
            newline='\n').write(u'\n'.join(rows) + u'\n')
    print(u'  chosen bins in total: %d; %s' % (total, dst))
    return 0


if __name__ == '__main__':
    sys.exit(main())
