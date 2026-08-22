#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
field_audit.py
THE FIELD CONSISTENCY CHECK, in one place.

The problem: a target's taxid list can occur in a locus library the target does
not belong to. An example: the saccharolytic bacteria target is bacterial
(Bacteroides, Proteiniphilum, Sphaerochaeta, Acetobacteroides) and its bins are
in the bacterial class; but taxid 1760811 (Acetobacteroides) also shows up in two
bins of the fungal ITS and 28S library. When no pair can be found in the correct
class, all that is left is a pair designed from a single bin at the wrong locus.
The table then shows that target as covered while it is not. A rule check cannot
catch this, because no rule is broken.

The fix: which fields (A archaea, B bacteria, F fungi) a target's member bins sit
in is derived FROM THE DATA. When a target carries bins in more than one field,
the design in the minority field is flagged.

This file holds no hand written field table. The field information comes from the
consensus file names and the taxid list in targets.tsv alone.
"""
import os
import re

ETIKET = re.compile(r"((?:A1|A2|B|F1|F2))-\d+_(\d+)$")


def taxid_alanlari(kons_klasoru):
    'Derives taxid -> {field: bin_count} from the consensus file names.'
    d = {}
    if not kons_klasoru or not os.path.isdir(kons_klasoru):
        return d
    for ad in sorted(os.listdir(kons_klasoru)):
        if not ad.endswith(".fasta"):
            continue
        et = re.sub(r"_(baskin|ref|self)?_?konsensus\.fasta$", "", ad)
        m = ETIKET.match(et)
        if not m:
            continue
        sinif, taxid = m.group(1), m.group(2)
        d.setdefault(taxid, {})
        d[taxid][sinif[0]] = d[taxid].get(sinif[0], 0) + 1
    return d


def hedef_taxidleri(hedefler_tsv):
    'targets.tsv -> {target_name: [taxid, ...]}'
    d = {}
    if not hedefler_tsv or not os.path.exists(hedefler_tsv):
        return d
    with open(hedefler_tsv, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or p[0] == "karar":
                continue
            d[p[1]] = [t.strip() for t in p[3].split(",") if t.strip()]
    return d


def alan_dagilimi(hedef, sinif, kons_klasoru=None, hedefler_tsv=None,
                  taxid_alan=None, hedef_taxid=None):
    """Returns (is_it_inconsistent, the distribution, the dominant field).

        When the first value is True, this (target, class) pair is not in the
        target's DOMINANT field and must not be delivered. When the target carries
        bins in one field alone, or the data is missing, it returns (False, {}, None):
        where nothing can be decided, nothing is blamed.
"""
    if taxid_alan is None:
        taxid_alan = taxid_alanlari(kons_klasoru)
    if hedef_taxid is None:
        hedef_taxid = hedef_taxidleri(hedefler_tsv)
    tl = hedef_taxid.get(hedef)
    if not tl or not taxid_alan:
        return False, {}, None
    dagilim = {}
    for t in tl:
        for alan, n in taxid_alan.get(t, {}).items():
            dagilim[alan] = dagilim.get(alan, 0) + n
    if len(dagilim) < 2:
        return False, dagilim, (max(dagilim, key=lambda x: dagilim[x])
                                if dagilim else None)
    baskin = max(dagilim, key=lambda x: dagilim[x])
    bu = (sinif or "")[:1]
    return (bool(bu) and bu != baskin), dagilim, baskin


def aciklama(dagilim, sinif, baskin):
    ad = {"A": "archaea", "B": "bacteria", "F": "fungi"}
    return ("the target's member bins are %s, and this pair is in the %s field "
            "while the dominant field is %s. The design was made from the reads "
            "of a locus library the target does not belong to, so it does not "
            "represent the target."
            % (", ".join("%s=%d" % (ad.get(k, k), v)
                         for k, v in sorted(dagilim.items())),
               ad.get((sinif or "")[:1], sinif), ad.get(baskin, baskin)))
