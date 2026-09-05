# -*- coding: utf-8 -*-
"""CONSENSUS TAXONOMY FROM vsearch HITS, one library against its own reference.

WHY
---
QIIME2's ready-made classifier is trained for 16S. On the fungal libraries and
on the long archaeal amplicon it produced meaningless labels. Here every
library is aligned against the reference that fits it (SILVA SSU, SILVA LSU,
UNITE) and each feature takes the taxonomy the best hits agree on.

THE RULE
--------
Hits within one identity point of the best are pooled; at every rank at least
51 per cent must agree; the chain stops where agreement breaks.

THE ALIGNMENT LENGTH ENTERS THE DECISION (2026-08-26)
-----------------------------------------------------
The hit length used to be read and never used, so a 300-base record at 100 per
cent beat a full-length record at 99.5 over 1,400 bases and its species name
went into the "best hit" column. That was the third instance of the same
length bias in this project. Now:
  * at equal identity the LONGER alignment wins,
  * a SPECIES name is given only when the alignment reaches the locus floor
    (SSU 1200, LSU 600, ITS 600 bases; one source, identity_verification.py);
    otherwise the name stays at genus level,
  * the alignment length and the reason for any demotion are columns.
The species node is appended to the taxonomy path only when the best hit's
name starts with the consensus genus and the identity reaches the species
threshold of the locus; a feature without a hit still gets a row ("Unassigned",
"no hit in any reference"), because missing is not the same as absent.

RUN
---
    python3 steps/qiime2_reference_taxonomy.py out_taxonomy.tsv hits1.tsv:ref1.fasta [hits2.tsv:ref2.fasta ...] [fasta:features.fasta]
"""
from __future__ import print_function
import collections
import io
import os
import re
import sys

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, KOK)
try:
    from verification.identity_verification import EN_AZ_HIZALAMA as _FLOOR, TUR_ESIGI as _SPECIES
    MIN_ALIGNMENT = {'SSU': _FLOOR['SSU'], 'LSU': _FLOOR['LSU'], 'ITS': _FLOOR['ITS']}
    SPECIES_THRESHOLD = dict(_SPECIES)
    THRESHOLD_SOURCE = 'verification/identity_verification.py'
except Exception as _e:                     # the module needs the aligner; keep a copy
    MIN_ALIGNMENT = {'SSU': 1200, 'LSU': 600, 'ITS': 600}
    SPECIES_THRESHOLD = {'SSU': 98.7, 'LSU': 98.7, 'LSU_MANTAR': 99.8, 'ITS': 99.6}
    THRESHOLD_SOURCE = 'FALLBACK COPY (%s)' % type(_e).__name__

UNNAMED = re.compile(r'^(uncultured|unidentified|unclassified|metagenome|environmental|incertae|unknown)', re.I)


def locus_of(path):
    name = os.path.basename(path).upper()
    if 'UNITE' in name or 'ITS' in name:
        return 'ITS'
    if 'LSU' in name or '28S' in name or '23S' in name:
        return 'LSU'
    return 'SSU'


def taxonomy_of_header(header, path):
    name = os.path.basename(path).upper()
    if 'UNITE' in name:
        p = header.split('|')
        return p[1] if len(p) > 1 else ''
    if 'PR2' in name:
        p = header.split('|')
        return ';'.join(p[4:]) if len(p) > 4 else ''
    p = header.split(' ', 1)
    return p[1] if len(p) > 1 else ''


def id_of_header(header, path):
    first = header.split()[0]
    return first.split('|')[0] if 'UNITE' in os.path.basename(path).upper() or first.count('|') else first


def embedded_taxonomy(target):
    """UNITE targets carry the taxonomy in their name: KJ734967|k__Fungi;p__...|SH..."""
    p = target.split('|')
    if len(p) >= 2 and ';' in p[1]:
        return p[1]
    return ''


def read_headers(fasta, ids):
    wanted = set(ids)
    out = {}
    if not wanted:
        return out
    with io.open(fasta, encoding='utf-8', errors='replace') as g:
        for line in g:
            if not line.startswith('>'):
                continue
            header = line[1:].rstrip('\n')
            k = id_of_header(header, fasta)
            if k in wanted:
                out[k] = taxonomy_of_header(header, fasta)
                if len(out) == len(wanted):
                    break
    return out


def ranks(taxonomy):
    """The rank list; the chain ends at the first name that carries no information."""
    out = []
    for x in [x.strip() for x in taxonomy.split(';') if x.strip()]:
        name = re.sub(r'^[kpcofgs]__', '', x).replace('_', ' ').strip()
        if (not name or UNNAMED.match(name) or name.endswith(' sp.') or name.endswith(' sp')
                or 'Incertae sedis' in name or 'metagenome' in name.lower()):
            break
        out.append(name)
    return out


def consensus(rank_lists, threshold=0.51):
    if not rank_lists:
        return ''
    longest = max(len(t) for t in rank_lists)
    shared = []
    for i in range(longest):
        count = collections.Counter(t[i] for t in rank_lists if len(t) > i)
        if not count:
            break
        name, n = count.most_common(1)[0]
        if name and n / float(len(rank_lists)) >= threshold:
            shared.append(name)
        else:
            break
    return ';'.join(shared)


def best_hits(hits):
    """One-reference rule: the hits of the reference that gave the best identity are
    taken (mixing UNITE's k__Fungi with SILVA's Eukaryota breaks agreement at the
    first rank). Within it, hits within one point of the best; equal identity is
    broken by the LONGER alignment, then by the target name for determinism."""
    best_ref = max(hits, key=lambda x: (x[0], x[1]))[3]
    same = [x for x in hits if x[3] == best_ref]
    top = max(x[0] for x in same)
    return sorted([x for x in same if x[0] >= top - 1.0],
                  key=lambda x: (-x[0], -x[1], x[2])), top


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    out_path = sys.argv[1]
    pairs = []
    all_features = []
    for a in sys.argv[2:]:
        if a.startswith('fasta:'):
            f = a.split(':', 1)[1]
            if os.path.exists(f):
                all_features = [x[1:].split()[0] for x in
                                io.open(f, encoding='utf-8', errors='replace') if x.startswith('>')]
            continue
        h, r = a.split(':', 1)
        pairs.append((h, r))
    hits = collections.defaultdict(list)      # feature -> [(identity, length, target, reference)]
    for h, r in pairs:
        if not os.path.exists(h) or os.path.getsize(h) == 0:
            continue
        for line in io.open(h, encoding='utf-8'):
            p = line.rstrip('\n').split('\t')
            if len(p) < 4:
                continue
            try:
                hits[p[0]].append((float(p[2]), int(p[3]), p[1], r))
            except ValueError:
                continue
    if not hits:
        io.open(out_path, 'w', encoding='utf-8').write(u'Feature ID\tTaxon\tIdentity\n')
        print('no hits; an empty taxonomy was written')
        return 0
    needed = collections.defaultdict(set)
    tax = {}
    for q, l in hits.items():
        good, _ = best_hits(l)
        for pid, ln, target, ref in good:
            g = embedded_taxonomy(target)
            if g:
                tax[target] = g
            else:
                needed[ref].add(target)
    for ref, ids in needed.items():
        tax.update(read_headers(ref, ids))
    rows = [u'Feature ID\tTaxon\tIdentity\tBestHit\tAlnLength\tNote']
    for q, l in sorted(hits.items()):
        good, top = best_hits(l)
        rank_lists = [ranks(tax.get(x[2], '')) for x in good]
        k = consensus([t for t in rank_lists if t])
        best = max(good, key=lambda x: (x[0], x[1]))
        er = ranks(tax.get(best[2], ''))
        name = er[-1] if er else u'Unassigned'
        note = u''
        locus = locus_of(best[3])
        floor = MIN_ALIGNMENT.get(locus, 600)
        if name != u'Unassigned' and u' ' in name and best[1] < floor:
            name = name.split()[0]
            note = (u'alignment %d bases, a species name at %s wants at least %d; '
                    u'brought down to genus' % (best[1], locus, floor))
        path = k or u'Unassigned'
        if path != u'Unassigned' and name != u'Unassigned' and u' ' in name and not note:
            last = path.split(u';')[-1].strip()
            te = SPECIES_THRESHOLD.get(locus, 98.7)
            if name.startswith(last + u' '):
                if top >= te:
                    path = path + u';' + name
                else:
                    note = (u'species %s is known but the identity is %.1f, the %s species '
                            u'threshold is %.1f; not added to the path' % (name, top, locus, te))
        rows.append(u'%s\t%s\t%.1f\t%s\t%d\t%s' % (q, path, top, name, best[1], note))
    present = set(x.split(u'\t')[0] for x in rows[1:])
    missing = [x for x in all_features if x not in present]
    for q in missing:
        rows.append(u'%s\tUnassigned\t0.0\tUnassigned\t0\tno hit in any reference '
                    u'(below threshold, or absent from these databases)' % q)
    io.open(out_path, 'w', encoding='utf-8').write(u'\n'.join(rows) + u'\n')
    print('taxonomy written: %s (%d features, %d without a hit; thresholds from %s)'
          % (out_path, len(rows) - 1, len(missing), THRESHOLD_SOURCE))
    return 0


if __name__ == '__main__':
    sys.exit(main())
