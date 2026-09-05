# -*- coding: utf-8 -*-
"""INDEPENDENT READ-LEVEL IN SILICO PCR. Shares nothing with the panel code.

Every pair is searched in the RAW READS of every bin:
  * mm0     both primers match exactly
  * mm1     at most one mismatch per primer, the last two 3' bases exact
  * product 40-700 bp, both orientations (a read may come from either strand)
At most --max-reads reads per bin, taken from the start of the file.

THE ARMS TRAP (2026-09-02), the reason mm0 is reported next to mm1
------------------------------------------------------------------
Two ordered oligos disagreed with the template at the third base from the 3'
end in 100 per cent of the target reads. Not a typo: the candidate generator
had produced deliberate -2/-3 mismatches (an ARMS-style variant) and the panel
criterion, "at most one mismatch", could not see it, so the ARMS candidate
scored artificially high (70.9x baseline vs 81.5x for the variant in the panel;
no difference when measured independently). Nothing in the delivery said ARMS.

The rule that catches it: when a pair amplifies its members at mm1 but the
mm0 rate over the members is ZERO, the pair carries a SYSTEMATIC mismatch and
the row says so. "mm <= 1" alone is never enough; every primer check reports
mm0 as well.

INPUT
-----
--pairs   TSV with a header; columns: name, source, class, forward, reverse,
          expected bp, members. `members` is either "MEMBERSHIP" (look the pair
          up in --membership) or "bin1;bin2;..." (the rest of the class's bins
          become the competitors).
--membership  TSV with a header; columns: name, ., ., members, mixed, competitors
          (each a ";"-separated list of bin labels). Optional.
The bins are <root>/fastq files/<lib>/<lib>[-_]reads[-_]<id>.fastq; a file whose
prefix does not match its folder is ignored (a stray copy of another bin).

RUN
---
    python3 verification/read_level_primer_check.py --root . --pairs pairs.tsv \
        [--membership membership.tsv] --out CHECK [--max-reads 4000] [--workers 4]
OUTPUT
    <out>_detail.tsv   one row per (pair, bin)
    <out>_summary.tsv  one row per pair with the verdict
"""
from __future__ import print_function
import argparse
import collections
import io
import os
import re
import sys
from multiprocessing import Pool

COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
MAX_READS = 4000


def rc(s):
    return ''.join(COMPLEMENT.get(c, 'N') for c in reversed(s))


def approximate_find(S, p, three_prime_first, max_mm):
    """Positions of p in S with <= max_mm mismatches (pigeonhole: one half exact).

    three_prime_first=True means the FIRST two bases of p (in read coordinates)
    are the primer's 3' end; False means the LAST two. The two 3' bases must be
    exact. Returns [(start, mismatches)].
    """
    L = len(p)
    h = L // 2
    starts = set()
    for part, offset in ((p[:h], 0), (p[h:], h)):
        i = S.find(part)
        while i != -1:
            starts.add(i - offset)
            i = S.find(part, i + 1)
    out = []
    for b in starts:
        if b < 0 or b + L > len(S):
            continue
        w = S[b:b + L]
        if three_prime_first:
            if w[:2] != p[:2]:
                continue
        elif w[-2:] != p[-2:]:
            continue
        mm = 0
        for x, y in zip(w, p):
            if x != y:
                mm += 1
                if mm > max_mm:
                    break
        if mm <= max_mm:
            out.append((b, mm))
    return out


def products(S, F, R, max_mm):
    """Is there an F..rc(R) or R..rc(F) product in the read? (length, total mismatches) or None."""
    rR, rF = rc(R), rc(F)
    best = None
    for fwd, rev in ((F, rR), (R, rF)):
        fs = approximate_find(S, fwd, False, max_mm)
        if not fs:
            continue
        gs = approximate_find(S, rev, True, max_mm)
        for fb, fm in fs:
            for gb, gm in gs:
                length = gb + len(rev) - fb
                if 40 <= length <= 700:
                    if best is None or fm + gm < best[1]:
                        best = (length, fm + gm)
    return best


def bin_files(root):
    base = os.path.join(root, 'fastq files')
    out = {}
    for d in sorted(os.listdir(base)):
        dd = os.path.join(base, d)
        if not os.path.isdir(dd):
            continue
        for f in sorted(os.listdir(dd)):
            m = re.match(r'^(.+?)[-_]reads[-_]([A-Za-z]*\d+)\.fastq$', f)
            if not m:
                continue
            if m.group(1).replace('_', '-') != d.replace('_', '-'):
                continue          # a stray copy of another bin
            out['%s_%s' % (d, m.group(2))] = os.path.join(dd, f)
    return out


def process_bin(arg):
    label, path, pairs, max_reads = arg
    count = collections.Counter()
    lengths = collections.defaultdict(collections.Counter)
    n = 0
    with io.open(path, encoding='ascii', errors='replace') as g:
        while n < max_reads:
            h = g.readline()
            if not h:
                break
            s = g.readline().strip().upper()
            g.readline()
            g.readline()
            if len(s) < 100:
                continue
            n += 1
            for name, F, R in pairs:
                p1 = products(s, F, R, 1)
                if p1:
                    count[(name, 'mm1')] += 1
                    lengths[name][p1[0]] += 1
                    if p1[1] == 0:
                        count[(name, 'mm0')] += 1
    return label, n, count, lengths


def verdict(member_rates, pooled_mm1, pooled_mm0, worst_competitor_rate, pooled_competitor_rate, universal):
    """The rule in one place. Returns (verdict text, systematic mismatch flag)."""
    inf = float('inf')
    ratio_pooled = (pooled_mm1 / pooled_competitor_rate) if pooled_competitor_rate > 0 else inf
    if member_rates and worst_competitor_rate > 0:
        ratio_worst = min(member_rates) / worst_competitor_rate
    elif member_rates and min(member_rates) > 0:
        ratio_worst = inf
    else:
        ratio_worst = 0.0
    systematic = bool(member_rates) and pooled_mm1 >= 10 and pooled_mm0 == 0
    if universal:
        v = u'COVERAGE %d/%d bins >=10%%' % (sum(1 for x in member_rates if x >= 10), len(member_rates))
    elif not member_rates or max(member_rates) < 10:
        v = u'NO / WEAK PRODUCT IN MEMBERS'
    elif ratio_worst >= 8:
        v = u'PASSED (worst bin >=8x)'
    elif ratio_pooled >= 8:
        v = u'PASSED POOLED, FAILED ON THE WORST BIN'
    else:
        v = u'FAILED (<8x)'
    if systematic:
        v += u' | SYSTEMATIC MISMATCH: mm0 = 0 in every member read (ARMS-style variant?)'
    return v, ratio_pooled, ratio_worst, systematic


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--root', default='.')
    ap.add_argument('--pairs', required=True)
    ap.add_argument('--membership', default=None)
    ap.add_argument('--out', required=True, help='output prefix')
    ap.add_argument('--max-reads', type=int, default=MAX_READS)
    ap.add_argument('--workers', type=int, default=4)
    a = ap.parse_args()
    root = os.path.abspath(a.root)
    pairs = []
    for i, s in enumerate(io.open(a.pairs, encoding='utf-8')):
        p = s.rstrip('\n').split('\t')
        if i == 0 or len(p) < 7:
            continue
        pairs.append(p)
    membership = {}
    if a.membership:
        for i, s in enumerate(io.open(a.membership, encoding='utf-8')):
            p = s.rstrip('\n').split('\t')
            if i == 0 or len(p) < 6:
                continue
            membership[p[0]] = tuple(set(x for x in p[k].split(';') if x) for k in (3, 4, 5))
    files = bin_files(root)
    print(u'bins: %d, pairs: %d, at most %d reads per bin, %d workers'
          % (len(files), len(pairs), a.max_reads, a.workers))
    sys.stdout.flush()
    short = [(c[0], c[3], c[4]) for c in pairs]
    jobs = [(k, y, short, a.max_reads) for k, y in sorted(files.items())]
    result = {}
    pool = Pool(a.workers)
    for i, (label, n, count, lengths) in enumerate(pool.imap_unordered(process_bin, jobs)):
        result[label] = (n, count, lengths)
        print(u'  [%3d/%d] %-16s %5d reads' % (i + 1, len(jobs), label, n))
        sys.stdout.flush()
    pool.close()
    pool.join()

    detail = [u'\t'.join([u'pair', u'source', u'bin', u'role', u'reads', u'mm0', u'mm1',
                          u'mm0 %', u'mm1 %', u'expected bp', u'observed bp (mode)'])]
    summary = [u'\t'.join([u'pair', u'source', u'class', u'expected bp', u'observed bp',
                           u'member bins', u'member mm1 % (min-max)', u'member mm1 % (pooled)',
                           u'member mm0 % (pooled)', u'competitor bins', u'worst competitor',
                           u'worst competitor mm1 %', u'competitor mm1 % (pooled)',
                           u'separation (pooled)', u'separation (worst)', u'systematic mismatch',
                           u'VERDICT'])]
    inf = float('inf')
    for c in pairs:
        name, source, klass, F, R, expected, members_field = c[:7]
        classes = set(klass.split('/'))
        if members_field == 'MEMBERSHIP':
            if name not in membership:
                print(u'  WARNING: %s is not in the membership table' % name)
                continue
            members, mixed, competitors = membership[name]
        else:
            members = set(members_field.split(';'))
            mixed = set()
            competitors = set(k for k in result if k.split('-')[0] in classes) - members
        members = set(k for k in members if k in result)
        competitors = set(k for k in competitors if k in result)
        mixed = set(k for k in mixed if k in result)

        def pooled(group, kind):
            n = sum(result[k][0] for k in group)
            h = sum(result[k][1][(name, kind)] for k in group)
            return (100.0 * h / n) if n else 0.0

        rate = {}
        for k in sorted(members | competitors | mixed):
            n, count, lengths = result[k]
            role = 'MEMBER' if k in members else ('MIXED' if k in mixed else 'COMPETITOR')
            m0, m1 = count[(name, 'mm0')], count[(name, 'mm1')]
            r1 = 100.0 * m1 / n if n else 0.0
            rate[k] = r1
            mode = lengths[name].most_common(1)[0][0] if lengths[name] else '-'
            detail.append(u'\t'.join([name, source, k, role, str(n), str(m0), str(m1),
                                      u'%.1f' % (100.0 * m0 / n if n else 0), u'%.1f' % r1,
                                      expected, str(mode)]))
        member_rates = [rate[k] for k in members]
        pooled_mm1, pooled_mm0 = pooled(members, 'mm1'), pooled(members, 'mm0')
        pooled_comp = pooled(competitors, 'mm1')
        worst = max(competitors, key=lambda k: rate[k]) if competitors else '-'
        worst_rate = rate[worst] if competitors else 0.0
        lens = collections.Counter()
        for k in members:
            lens.update(result[k][2][name])
        observed = lens.most_common(1)[0][0] if lens else '-'
        v, ratio_pooled, ratio_worst, systematic = verdict(
            member_rates, pooled_mm1, pooled_mm0, worst_rate, pooled_comp, 'universal' in name)

        def fmt(x):
            return u'inf' if x == inf else u'%.1f' % x
        summary.append(u'\t'.join([name, source, klass, expected, str(observed), str(len(members)),
                                   (u'%.1f-%.1f' % (min(member_rates), max(member_rates))) if member_rates else '-',
                                   u'%.1f' % pooled_mm1, u'%.1f' % pooled_mm0, str(len(competitors)), worst,
                                   u'%.1f' % worst_rate, u'%.2f' % pooled_comp, fmt(ratio_pooled),
                                   fmt(ratio_worst), u'YES' if systematic else u'no', v]))
    io.open(a.out + '_detail.tsv', 'w', encoding='utf-8', newline='\n').write(u'\n'.join(detail) + u'\n')
    io.open(a.out + '_summary.tsv', 'w', encoding='utf-8', newline='\n').write(u'\n'.join(summary) + u'\n')
    print(u'written: %s_summary.tsv, %s_detail.tsv' % (a.out, a.out))
    return 0


if __name__ == '__main__':
    sys.exit(main())
