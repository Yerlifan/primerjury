# -*- coding: utf-8 -*-
"""CHOOSE WHICH BINS ENTER THE IDENTITY CHAIN, largest first, with a floor.

WHY
---
The binning step writes every bin with at least the read floor; one bacterial
barcode gave 1,962 of them (measured). The chain builds a consensus and scans
thirteen reference databases per bin, so 99 bins took a night. Thousands of
small bins are mostly chimeras, fragments and single-error reads; scanning
them all would take days and add nothing. The bins that are not chosen are
NOT deleted: they stay where they are and SELECTION_TABLE.tsv records why.

THE CRITERION, kept in one place (here)
---------------------------------------
For tables written by clean_bins.py (v4, one organism per bin):
* the bin has >= MIN_READS_SELECT (500) reads after the long-read filter, OR
  it is among the TOP_PER_WINDOW (5) largest of the barcode; at most
  MAX_PER_BARCODE (40). Measured 2026-09-07: a share floor of 0.5 per cent
  left 5 to 12 bins per barcode where v3 chose 22 to 38; 500 reads is the
  floor the PAK polish needs (a 150-read sample plus a sub-population of
  >= 60 per cent and >= 5 reads).
* SATELLITE rule. About 3.5 per cent of the reads carry the amplicon plus a
  ~310 bp library tail shared between reads (not 16S). A tailed read covers
  the centre fully but is itself covered at 82 per cent, below the 85 per
  cent assignment floor, so it is left over after round one; round two seeds
  from the leftovers and the tailed reads of an organism that already has a
  bin form a second bin of the same organism (A1-1 BIN7: 31 full reads and
  2,152 segments, 98 per cent of them 98.9 per cent identical to BIN1's
  centre). A candidate whose sampled reads (60, full reads preferred) align
  >= ASSIGN_IDENTITY to an already CHOSEN centre over >= 80 per cent of that
  centre, for >= half of the sample, is that bin's satellite: not chosen,
  the parent written in the reason. Bins whose reads do not align to any
  chosen centre at all are separate organisms (long 16S variants were).
* LONG-READ FILTER. Reads longer than FULL_HIGH x the amplicon are two-copy
  concatemer "segments" that a run before the segment length cap could
  produce (measured 0.6 to 3 per cent; two small bins were almost nothing
  else). They are dropped from the COPY that enters the chain and counted in
  the table; the bin's read count for the criterion is the filtered count.
  The source files are untouched.

For tables written by bin_reads.py (v3, length windows) the old rule holds:
share >= MIN_SHARE (0.5 per cent) or among the TOP_PER_WINDOW largest of its
window, at most MAX_PER_BARCODE.

RUN
---
    python3 steps/select_bins.py --source "fastq files_bins" --target "fastq files" [--name A2-1 ...]

The chosen bins are COPIED into --target in the layout the chain reads.
"""
from __future__ import print_function
import argparse
import collections
import importlib.util
import io
import os
import random
import shutil
import subprocess
import sys
import tempfile

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, KOK)

MIN_SHARE = 0.5          # per cent of the barcode's reads (v3 tables)
MIN_READS_SELECT = 500   # v4 tables: reads (full + segments) after the long-read filter
TOP_PER_WINDOW = 5       # the largest N bins of every window (v4: of the barcode) are always taken
MAX_PER_BARCODE = 40
SATELLITE_SAMPLE = 60
SATELLITE_SHARE = 0.5    # this share of the sample on one chosen centre -> satellite
SATELLITE_COVERAGE = 0.8  # of the chosen centre


def _clean_bins():
    path = os.path.join(KOK, 'steps', 'clean_bins.py')
    spec = importlib.util.spec_from_file_location('clean_bins', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_centres(path):
    """<name>_centres.fasta of clean_bins: '>read_id BINn reads=..' -> {BINn: sequence}; {} when absent."""
    out = {}
    if not os.path.isfile(path):
        return out
    name = None
    for line in io.open(path, encoding='ascii'):
        if line.startswith('>'):
            name = line[1:].split()[1]
            out[name] = []
        elif name:
            out[name].append(line.strip())
    return dict((k, ''.join(v)) for k, v in out.items())


def count_long(src, cap_bp, fastq_stream):
    """(reads, reads longer than cap_bp)."""
    n = long = 0
    for _r, seq, _q in fastq_stream(src):
        n += 1
        if len(seq) > cap_bp:
            long += 1
    return n, long


def copy_filtered(src, dst, cap_bp, fastq_stream):
    """Copy a bin file dropping reads longer than cap_bp; (written, dropped)."""
    written = dropped = 0
    with io.open(dst, 'w', encoding='ascii') as g:
        for rid, seq, qual in fastq_stream(src):
            if len(seq) > cap_bp:
                dropped += 1
                continue
            g.write(u'@%s\n%s\n+\n%s\n' % (rid, seq, qual))
            written += 1
    return written, dropped


def satellite_from_rows(rows, identity, coverage, sample_size):
    """Pure part of the satellite rule: PAF rows of the sample against the chosen centres ->
    (parent, hits) when >= SATELLITE_SHARE of the sample sits on one centre, else None."""
    best = {}
    for rid, _ql, _qs, _qe, centre, tl, ts, te, ident in rows:
        if ident >= identity and (te - ts) >= coverage * tl and ident > best.get(rid, (0.0, None))[0]:
            best[rid] = (ident, centre)
    tally = collections.Counter(c for _i, c in best.values())
    if not tally:
        return None
    parent, hits = tally.most_common(1)[0]
    if hits >= SATELLITE_SHARE * sample_size:
        return parent, hits
    return None


def satellite_of(mm2, bin_fq, chosen_centres, tmp, cb, fastq_stream):
    """(parent, hits, sample) or None; cb is the clean_bins module (PAF rows, identity floor)."""
    if not chosen_centres:
        return None
    full, every = [], []
    for rid, seq, _q in fastq_stream(bin_fq):
        every.append((rid, seq))
        if '_p' not in rid:
            full.append((rid, seq))
    pool = full if len(full) >= 10 else every
    if not pool:
        return None
    sample = random.Random(1).sample(pool, min(SATELLITE_SAMPLE, len(pool)))
    mf, qf = os.path.join(tmp, 'centres.fa'), os.path.join(tmp, 'sample.fa')
    with io.open(mf, 'w', encoding='ascii') as g:
        for name, seq in chosen_centres.items():
            g.write(u'>%s\n%s\n' % (name, seq))
    with io.open(qf, 'w', encoding='ascii') as g:
        for rid, seq in sample:
            g.write(u'>%s\n%s\n' % (rid, seq))
    with io.open(os.devnull, 'w') as null:
        paf = subprocess.check_output([mm2, '-x', 'map-ont', '-c', '--secondary=no', '-t', '2', mf, qf], stderr=null)
    hit = satellite_from_rows(cb.paf_rows(paf.decode('ascii', 'replace')), cb.ASSIGN_IDENTITY,
                              SATELLITE_COVERAGE, len(sample))
    return (hit[0], hit[1], len(sample)) if hit else None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--source', required=True, help='bin root written by clean_bins.py or bin_reads.py')
    ap.add_argument('--target', required=True, help='where the chosen bins are copied ("fastq files")')
    ap.add_argument('--name', nargs='*', default=None, help='only these barcodes (default: all with a BIN_TABLE.tsv)')
    a = ap.parse_args(argv)
    cb = _clean_bins()
    br = cb._bin_reads()
    src, dst = os.path.abspath(a.source), os.path.abspath(a.target)
    if not os.path.isdir(dst):
        os.makedirs(dst)
    names = a.name or sorted(d for d in os.listdir(src)
                             if os.path.isfile(os.path.join(src, d, 'BIN_TABLE.tsv')))
    rows = [u'\t'.join([u'barcode', u'bin', u'window', u'reads', u'share %', u'chosen', u'reason',
                        u'long dropped (>%.2fx amplicon)' % cb.FULL_HIGH])]
    total = satellites = 0
    mm2 = None
    tmp = tempfile.mkdtemp(prefix='select_bins_')
    try:
        for name in names:
            t = [s.rstrip(u'\n').split(u'\t') for s in
                 io.open(os.path.join(src, name, 'BIN_TABLE.tsv'), encoding='utf-8') if s.strip()]
            head = t[0]
            ib, ish, ifile = head.index(u'bin'), head.index(u'share % (barcode)'), head.index(u'file')
            ir = head.index(u'reads')
            iseg = head.index(u'concatemer segments') if u'concatemer segments' in head else None   # v4 table
            iw = head.index(u'window bp') if u'window bp' in head else None
            bins = []
            for r in t[1:]:
                if r[ib].startswith(u'BIN'):
                    r = list(r)
                    r[ir] = str(int(r[ir]) + (int(r[iseg]) if iseg is not None else 0))
                    if iw is None:
                        r.append(u'full')
                    bins.append(r)
            if iw is None:
                iw = len(head)
            v4 = iseg is not None
            centres = read_centres(os.path.join(src, name, u'%s_centres.fasta' % name)) if v4 else {}
            cap = cb.FULL_HIGH * cb.AMPLICON_BP.get(name.split(u'-')[0], 0) if v4 else None
            long_of = {}
            if cap:
                for r in bins:                               # the criterion uses the FILTERED count
                    n, long = count_long(os.path.join(src, name, r[ifile]), cap, br.fastq_stream)
                    long_of[r[ib]] = long
                    r[ir] = str(n - long)
            rank = {}
            for r in sorted(bins, key=lambda r: -int(r[ir])):
                rank.setdefault(r[iw], []).append(r[ib])
            chosen, chosen_centres = [], collections.OrderedDict()
            for r in sorted(bins, key=lambda r: -int(r[ir])):
                share = float(r[ish])
                pos = rank[r[iw]].index(r[ib]) + 1
                if len(chosen) >= MAX_PER_BARCODE:
                    reason, take = u'barcode ceiling %d' % MAX_PER_BARCODE, False
                elif v4 and int(r[ir]) >= MIN_READS_SELECT:
                    reason, take = u'reads >= %d' % MIN_READS_SELECT, True
                elif not v4 and share >= MIN_SHARE:
                    reason, take = u'share >= %.1f%%' % MIN_SHARE, True
                elif pos <= TOP_PER_WINDOW:
                    reason, take = u'among the %d largest of its window (rank %d)' % (TOP_PER_WINDOW, pos), True
                elif v4:
                    reason, take = u'reads %s < %d and rank %d' % (r[ir], MIN_READS_SELECT, pos), False
                else:
                    reason, take = u'share %.2f%% < %.1f%% and rank %d in its window' % (share, MIN_SHARE, pos), False
                if take and centres:
                    if mm2 is None:
                        mm2 = br.tool('minimap2')
                    hit = satellite_of(mm2, os.path.join(src, name, r[ifile]), chosen_centres, tmp, cb, br.fastq_stream)
                    if hit:
                        reason = (u'satellite of %s: %d/%d sampled reads >= %.0f%% identical to its centre '
                                  u'(tailed or truncated reads of the same organism)'
                                  % (hit[0], hit[1], hit[2], 100 * cb.ASSIGN_IDENTITY))
                        take = False
                        satellites += 1
                if take:
                    chosen.append(r)
                    if r[ib] in centres:
                        chosen_centres[r[ib]] = centres[r[ib]]
                rows.append(u'\t'.join([name, r[ib], r[iw], r[ir], r[ish], u'YES' if take else u'no', reason,
                                        str(long_of.get(r[ib], u''))]))
            d = os.path.join(dst, name)
            if not os.path.isdir(d):
                os.makedirs(d)
            for r in chosen:
                s, t2 = os.path.join(src, name, r[ifile]), os.path.join(d, r[ifile])
                if cap:
                    copy_filtered(s, t2, cap, br.fastq_stream)
                elif not os.path.exists(t2) or os.path.getsize(t2) != os.path.getsize(s):
                    shutil.copyfile(s, t2)
            total += len(chosen)
            print(u'  %-6s bins %4d -> chosen %3d (windows: %s)'
                  % (name, len(bins), len(chosen), u', '.join(sorted(rank))))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    io.open(os.path.join(dst, 'SELECTION_TABLE.tsv'), 'w', encoding='utf-8',
            newline='\n').write(u'\n'.join(rows) + u'\n')
    print(u'  chosen bins in total: %d (left out as satellites: %d); %s' % (total, satellites, dst))
    return 0


if __name__ == '__main__':
    sys.exit(main())
