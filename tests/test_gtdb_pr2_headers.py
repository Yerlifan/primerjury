# -*- coding: utf-8 -*-
"""GTDB and PR2 headers reach the name parser (2026-09-06).

GTDB: "RS_GCF_x~contig d__...;g__Proteiniphilum_A;s__Proteiniphilum_A saccharofermentans
[location=..]" keeps the genus suffix and drops the bracket tail; "s__Genus sp002498885"
is a GTDB species cluster and passes as a name; a genus with digits (UBA1234) names
nothing. PR2: "acc|18S_rRNA|...|Genus|Genus_species" gives the species, "Genus_sp."
the genus.

RUN
    python3 tests/test_gtdb_pr2_headers.py
"""
from __future__ import print_function
import os
import sys

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(KOK, 'steps'))
sys.path.insert(0, os.path.join(KOK, 'verification'))
import target_identity as ti   # noqa: E402

CASES = [
    (u'RS_GCF_018344175.1~NZ_JAAMUS010000074.1 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli [location=16..1459] [ssu_len=1444]',
     u'Escherichia coli'),
    (u'GB_GCA_2~y d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Dysgonomonadaceae;g__Proteiniphilum_A;s__Proteiniphilum_A saccharofermentans [location=1..1400]',
     u'Proteiniphilum_A saccharofermentans'),
    (u'GB_GCA_1~x d__Archaea;p__Halobacteriota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina sp002498885 [location=1..1400]',
     u'Methanosarcina sp002498885'),
    (u'GB_GCA_3~z d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae;g__UBA1234;s__UBA1234 sp001234567 [location=1..1400]',
     None),
    (u'AB353770.1.1740_U|18S_rRNA|nucleus||Eukaryota|TSAR|Alveolata|Dinoflagellata|Dinophyceae|Peridiniales|Kryptoperidiniaceae|Unruhdinium|Unruhdinium_kevei',
     u'Unruhdinium kevei'),
    (u'KX1|18S_rRNA|nucleus||Eukaryota|Obazoa|Opisthokonta|Fungi|Ascomycota|Sordariomycetes|Microascales|Microascaceae|Petriella|Petriella_sp.',
     u'Petriella'),
    (u'NR_172285.1 Petriella musispora CBS 745.69 ITS region; from TYPE material', u'Petriella musispora'),
    (u'UDB016649|k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Thelephorales;f__Thelephoraceae;g__Thelephora;s__Thelephora_albomarginata|SH1281904.10FU',
     u'Thelephora albomarginata'),
]


def main():
    failed = 0
    for header, want in CASES:
        got = ti.ad_ayikla(header)
        ok = got == want
        failed += 0 if ok else 1
        print('  %s %-45s -> %r' % ('PASS' if ok else 'FAIL', header[-45:], got))
    g, e = ti.cins_epitet(u'Proteiniphilum_A saccharofermentans')
    ok = (g, e) == (u'Proteiniphilum_A', u'saccharofermentans')
    failed += 0 if ok else 1
    print('  %s cins_epitet keeps the GTDB genus suffix -> %r' % ('PASS' if ok else 'FAIL', (g, e)))
    print('  %d of %d checks passed' % (len(CASES) + 1 - failed, len(CASES) + 1))
    return 0 if not failed else 1


if __name__ == '__main__':
    sys.exit(main())
