# -*- coding: utf-8 -*-
"""BIN THE READS OF ONE BARCODE WITHOUT A CLASSIFIER.

WHAT IT IS FOR
--------------
Takes the raw reads of one barcode (one library, one time point) and groups
them by amplicon length and by sequence similarity, then writes every group as
a bin in the layout the rest of the chain already reads:

    <out>/<name>/<name>-reads_BIN<n>.fastq
    <out>/<name>/<name>_centres.fasta
    <out>/<name>/BIN_TABLE.tsv

The consensus, identity and design stages then run on the bins unchanged. A
bin is not a species; the name comes from the identity stage, which is the
point: the bin is defined by the data, not by a database label.

WHY
---
The study's bins were "barcode x Kraken2 taxid". That definition needs a
Kraken2 run, and on the machine the study ran on the 196 GB database could
not be held in memory (12 reads per second through memory mapping, measured).
It also inherits every disagreement between the classifier label and the
alignment-based identity that the README opens with. Binning from the reads
themselves needs no database and carries no label into the identity step.

THE METHOD, every setting measured (2026-09-04)
-----------------------------------------------
0. LENGTH PEAKS. A barcode does not carry a single amplicon. Measured on the
   study data: one archaeal barcode had the 4.2-4.4 kb operon at only 6 per
   cent of reads, 1.2-1.4 kb at 12 per cent and 400-600 bp at 10 per cent; a
   bacterial barcode had 1.4-1.6 kb at 54 per cent; a fungal one 3.6-3.8 kb at
   54 per cent. Each peak of the 100 bp length histogram (reads >= 300 bp,
   share >= --min-peak-share, at most --max-peaks) becomes a WINDOW of
   peak +- 15 per cent. Choosing a single window would have ignored most reads.
1. SEEDS. --seeds reads are drawn at random from the window (fixed --seed).
2. SEED CLUSTERING. minimap2 ava-ont over the seeds; two seeds are joined when
   identity (1 - de) >= --identity and the alignment covers >= 80 per cent of
   the shorter one; connected components are the clusters. The centre is the
   member closest to the component's MEDIAN length, not the longest: the
   longest member was a 5.2 kb chimera in the first trial. vsearch was tried
   first and dropped: cluster_fast on 4 kb reads did not finish 1,200 seeds in
   twenty minutes with an exhaustive search, and with its default k-mer
   prefilter it rejected true relatives (none of 1,160 centres merged at 97
   per cent). ava-ont takes seconds.
3. ASSIGNMENT. Every read of the window is mapped to the centres with minimap2
   map-ont; a read joins a bin when the alignment covers >= 80 per cent of the
   read and identity >= 0.90. A read that fits nowhere is COUNTED as
   unassigned (chimera, rare organism), not dropped silently.
4. Bins with >= --min-reads reads are written; the bin number runs over the
   whole barcode and the window is recorded in the table.

WHAT IT DOES NOT PROVE
----------------------
A window is not an amplicon: in one barcode the centre of the 1.3 kb peak hit
no rRNA database at all (an off-target product, most likely). The chain
reports such a bin as "cannot be named" and the bin is still written, so that
it is seen rather than lost.

RUN
---
    python3 steps/bin_reads.py --fastq barcode05.fastq --name A2-1 \
        --out "fastq files" [--identity 0.97] [--threads 2] [--seeds 1200]

Temporary files go beside the output directory, not to /tmp: on WSL /tmp lives
in the virtual disk, which does not shrink when files are deleted.
"""
from __future__ import print_function
import argparse
import collections
import io
import os
import random
import shutil
import subprocess
import sys
import tempfile

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, KOK)
try:
    from verification.identity_verification import EN_AZ_OKUMA as MIN_READS_DEFAULT
except Exception:                       # the module needs the aligner and the databases
    MIN_READS_DEFAULT = 5

ASSIGN_IDENTITY = 0.90
ASSIGN_COVERAGE = 0.80


def tool(name):
    for y in ('/usr/bin/' + name, '/usr/local/bin/' + name):
        if os.path.exists(y):
            return y
    p = shutil.which(name)
    if p:
        return p
    raise SystemExit('STOPPED: %s not found' % name)


def fastq_stream(path):
    with io.open(path, encoding='ascii', errors='replace') as g:
        while True:
            h = g.readline()
            if not h:
                return
            seq = g.readline().rstrip('\n')
            g.readline()
            qual = g.readline().rstrip('\n')
            if h.startswith('@'):
                yield h[1:].split()[0], seq, qual


def run(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p.returncode != 0:
        raise SystemExit('STOPPED: %s -> %d: %s'
                         % (cmd[0], p.returncode, p.stderr.decode('utf-8', 'replace')[:400]))
    return p.stdout.decode('utf-8', 'replace')


def find_length_peaks(lengths, min_share=0.03, max_peaks=4):
    """Read-length peaks -> [(low, high, peak, share_percent)], by share, non-overlapping."""
    n = len(lengths)
    hist = collections.Counter(b // 100 for b in lengths if b >= 300)
    out, used = [], []
    for bucket, _ in hist.most_common():
        peak = bucket * 100 + 50
        low, high = int(peak * 0.85), int(peak * 1.15)
        if any(not (high < a or low > b) for a, b in used):
            continue
        share = sum(1 for b in lengths if low <= b <= high) / float(n)
        if share < min_share:
            continue
        out.append((low, high, peak, 100.0 * share))
        used.append((low, high))
        if len(out) >= max_peaks:
            break
    return out


def paf_pairs(paf, min_identity, min_coverage):
    """(a, b, identity) from PAF lines: >= min_coverage of the shorter one covered, 1-de >= min_identity."""
    for line in paf.splitlines():
        p = line.split('\t')
        if len(p) < 12 or p[0] == p[5]:
            continue
        ql, qs, qe, tl, ts, te = (int(p[1]), int(p[2]), int(p[3]),
                                  int(p[6]), int(p[7]), int(p[8]))
        de = None
        for tag in p[12:]:
            if tag.startswith('de:f:'):
                de = float(tag[5:])
                break
        if de is None:
            continue
        shorter = min(ql, tl)
        cov = min(qe - qs, te - ts) / float(shorter) if shorter else 0
        if cov >= min_coverage and (1.0 - de) >= min_identity:
            yield p[0], p[5], 1.0 - de


def cluster_seeds(seqs, identity, mm2, threads, tmp, label):
    """Connected components over minimap2 ava-ont. Returns {centre: [members]}."""
    fa = os.path.join(tmp, label + '.fa')
    with io.open(fa, 'w', encoding='ascii') as g:
        for name, s in seqs.items():
            g.write(u'>%s\n%s\n' % (name, s))
    paf = run([mm2, '-x', 'ava-ont', '-c', '-t', str(threads), fa, fa])
    parent = {name: name for name in seqs}

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b, _ in paf_pairs(paf, identity, ASSIGN_COVERAGE):
        ra, rb = root(a), root(b)
        if ra != rb:
            parent[rb] = ra
    comp = collections.defaultdict(list)
    for name in seqs:
        comp[root(name)].append(name)
    out = {}
    for members in comp.values():
        lens = sorted(len(seqs[x]) for x in members)
        median = lens[len(lens) // 2]
        centre = min(members, key=lambda x: (abs(len(seqs[x]) - median), x))
        out[centre] = members
    return out


def process_window(a, mm2, tmp, low, high, label):
    """One length window: seeds -> clusters -> assignment. Returns (n, centres, assignment, clusters, window_fq)."""
    window = os.path.join(tmp, 'window_%s.fq' % label)
    n = 0
    with io.open(window, 'w', encoding='ascii') as g:
        for rid, seq, qual in fastq_stream(a.fastq):
            if low <= len(seq) <= high:
                g.write(u'@%s\n%s\n+\n%s\n' % (rid, seq, qual))
                n += 1
    if n == 0:
        return 0, {}, {}, {}, window
    random.seed(a.seed)
    chosen = set(random.sample(range(n), min(a.seeds, n)))
    seeds = collections.OrderedDict()
    for j, (rid, seq, _) in enumerate(fastq_stream(window)):
        if j in chosen:
            seeds[rid] = seq
    clusters = cluster_seeds(seeds, a.identity, mm2, a.threads, tmp, 'seeds_' + label)
    centres = collections.OrderedDict((c, seeds[c]) for c in clusters)
    centre_fa = os.path.join(tmp, 'centres_%s.fa' % label)
    with io.open(centre_fa, 'w', encoding='ascii') as g:
        for c, s in centres.items():
            g.write(u'>%s\n%s\n' % (c, s))
    paf = run([mm2, '-x', 'map-ont', '-c', '--secondary=no', '-t', str(a.threads),
               centre_fa, window])
    assignment = {}
    for line in paf.splitlines():
        p = line.split('\t')
        if len(p) < 12:
            continue
        rid, qlen, qs, qe, target = p[0], int(p[1]), int(p[2]), int(p[3]), p[5]
        de = None
        for tag in p[12:]:
            if tag.startswith('de:f:'):
                de = float(tag[5:])
                break
        if de is None:
            continue
        ident, cov = 1.0 - de, (qe - qs) / float(qlen) if qlen else 0
        if cov >= ASSIGN_COVERAGE and ident >= ASSIGN_IDENTITY:
            if rid not in assignment or ident > assignment[rid][1]:
                assignment[rid] = (target, ident)
    return n, centres, assignment, clusters, window


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--fastq', required=True, help='raw reads of one barcode')
    ap.add_argument('--name', required=True, help='library-timepoint label, e.g. A2-1; becomes the folder and file prefix')
    ap.add_argument('--out', required=True, help='bin root, normally the "fastq files" directory')
    ap.add_argument('--identity', type=float, default=0.97, help='seed clustering identity (1 - de)')
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--seeds', type=int, default=1200, help='seed reads per window')
    ap.add_argument('--seed', type=int, default=7, help='random seed')
    ap.add_argument('--min-reads', type=int, default=MIN_READS_DEFAULT,
                    help='smallest bin that is written (default: the identity step\'s read floor)')
    ap.add_argument('--min-peak-share', type=float, default=0.03)
    ap.add_argument('--max-peaks', type=int, default=4)
    a = ap.parse_args()
    mm2 = tool('minimap2')
    out = os.path.join(os.path.abspath(a.out), a.name)
    if not os.path.isdir(out):
        os.makedirs(out)
    lengths = [len(s) for _, s, _ in fastq_stream(a.fastq)]
    if not lengths:
        raise SystemExit('STOPPED: %s is empty' % a.fastq)
    windows = find_length_peaks(lengths, a.min_peak_share, a.max_peaks)
    if not windows:
        raise SystemExit('STOPPED: no length peak found for %s' % a.name)
    print(u'  %s: %d reads; windows: %s'
          % (a.name, len(lengths), u', '.join(u'%d-%d bp (peak %d, %.1f%%)' % w for w in windows)))
    tmp = tempfile.mkdtemp(prefix='_tmp_bin_reads_', dir=os.path.dirname(out))
    total = len(lengths)
    table = [u'\t'.join([u'bin', u'window bp', u'reads', u'share % (barcode)', u'centre read',
                         u'centre length bp', u'median identity', u'file'])]
    summary = []
    no = 0
    in_windows = 0
    try:
        centre_file = io.open(os.path.join(out, '%s_centres.fasta' % a.name), 'w', encoding='ascii')
        for (low, high, peak, share) in windows:
            label = u'%d-%d' % (low, high)
            n, centres, assignment, clusters, window = process_window(
                a, mm2, tmp, low, high, label.replace('-', '_'))
            in_windows += n
            count = collections.Counter(c for c, _ in assignment.values())
            idents = collections.defaultdict(list)
            for c, ident in assignment.values():
                idents[c].append(ident)
            passing = sorted((c for c in centres if count.get(c, 0) >= a.min_reads),
                             key=lambda c: -count[c])
            local_no = {}
            for c in passing:
                no += 1
                local_no[c] = no
                centre_file.write(u'>%s BIN%d window=%s reads=%d\n%s\n'
                                  % (c, no, label, count[c], centres[c]))
            writers = {c: io.open(os.path.join(out, '%s-reads_BIN%d.fastq' % (a.name, local_no[c])),
                                  'w', encoding='ascii') for c in passing}
            for rid, seq, qual in fastq_stream(window):
                c = assignment.get(rid, (None, 0))[0]
                if c in writers:
                    writers[c].write(u'@%s\n%s\n+\n%s\n' % (rid, seq, qual))
            for g in writers.values():
                g.close()
            for c in passing:
                med = sorted(idents[c])
                table.append(u'\t'.join([u'BIN%d' % local_no[c], label, str(count[c]),
                                         u'%.2f' % (100.0 * count[c] / total), c, str(len(centres[c])),
                                         u'%.4f' % med[len(med) // 2],
                                         u'%s-reads_BIN%d.fastq' % (a.name, local_no[c])]))
            small = sum(v for c, v in count.items() if c not in local_no)
            unassigned = n - len(assignment)
            summary.append(u'\t'.join([u'small bins (<%d)' % a.min_reads, label, str(small),
                                       u'%.2f' % (100.0 * small / total), u'-', u'-', u'-', u'-']))
            summary.append(u'\t'.join([u'unassigned', label, str(unassigned),
                                       u'%.2f' % (100.0 * unassigned / total), u'-', u'-', u'-', u'-']))
            print(u'  window %s bp (%.1f%%): %d reads, %d seed clusters, bins >=%d: %d, '
                  u'small %d, unassigned %d (%.1f%%)'
                  % (label, share, n, len(centres), a.min_reads, len(passing), small, unassigned,
                     100.0 * unassigned / n if n else 0))
            for c in passing[:5]:
                med = sorted(idents[c])
                print(u'    BIN%-3d %7d reads  %5.2f%%  centre %d bp  median identity %.3f'
                      % (local_no[c], count[c], 100.0 * count[c] / total, len(centres[c]),
                         med[len(med) // 2]))
        centre_file.close()
        outside = total - in_windows
        summary.append(u'\t'.join([u'outside the windows (length)', u'-', str(outside),
                                   u'%.2f' % (100.0 * outside / total), u'-', u'-', u'-', u'-']))
        io.open(os.path.join(out, 'BIN_TABLE.tsv'), 'w', encoding='utf-8',
                newline='\n').write(u'\n'.join(table + summary) + u'\n')
        print(u'  bins in total %d; outside the windows %d (%.1f%%); written: %s'
              % (no, outside, 100.0 * outside / total, out))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
