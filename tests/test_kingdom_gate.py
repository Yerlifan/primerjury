# -*- coding: utf-8 -*-
"""The kingdom gate, the synonym table and the nt reconciliation helper (2026-09-09).

What must stay true:

1. non_fungal_eukaryote reads SILVA paths (with or without an accession) and PR2 headers; a Fungi path,
   a RefSeq title and a Bacteria path are not "outside"; PR2 placeholders (_X, _sp) are skipped for the name.
2. kingdom_gate compares WITHIN a locus: a ciliate at 99.65 on 18S against fungal 18S at 87.8 opens the gate;
   fungal 18S at 99 (margin under 2 points) keeps it shut; 94 per cent (under the floor) keeps it shut; 120 bp
   (under the evidence floor) keeps it shut; the ITS locus never opens it.
3. current_name maps synonyms and strips 'cf.'.
4. reconcile_with_nt lowers a species to 'cf.' only for a different, non-synonymous species at or above the
   threshold and inside the margin over >= 250 bp; synonyms, genus-rank hits and weaker hits leave it alone.

RUN
    python3 tests/test_kingdom_gate.py
"""
from __future__ import print_function
import os
import sys

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(KOK, 'verification'))
from identity_verification import (non_fungal_eukaryote, current_name, reconcile_with_nt,   # noqa: E402
                                   KINGDOM_GATE_IDENTITY, KINGDOM_GATE_MARGIN)
from locus_decision import kingdom_gate                                                     # noqa: E402

SIL = u'Eukaryota;SAR;Alveolata;Ciliophora;Spirotrichea;Stichotrichia;Parakahliella;Parakahliella macrostoma'
FUN = u'Eukaryota;Obazoa;Opisthokonta;Nucletmycea;Fungi;Ascomycota;Petriella;Petriella musispora'
REF = u'NR_172285.1 Petriella musispora CBS 745.69 ITS region; from TYPE material'
PR2 = u'MH393767.1.1725_U|18S_rRNA|nucleus||Eukaryota|TSAR|Alveolata|Ciliophora|Spirotrichea|Hypotrichia|Hypotrichia_X|Parakahliella|Parakahliella_macrostoma'
PR2F = u'EU984280.1.833_UC|18S_rRNA|nucleus|strain_CBS385.87|Eukaryota|Obazoa|Opisthokonta|Fungi|Ascomycota|Pezizomycotina|Sordariomycetes|Petriella|Petriella_setifera'
PR2X = u'X|18S_rRNA|nucleus||Eukaryota|TSAR|Rhizaria|Cercozoa|Sarcomonadea|Glissomonadida|Sandonidae|Sandonidae_X|Sandonidae_X_sp'
R18 = u'NG_062754.1 Petriella setifera CBS 385.87 18S ribosomal RNA gene, partial sequence'


def main():
    fails = 0

    def ok(cond, label):
        print('  %s %s' % ('ok  ' if cond else 'FAIL', label))
        return 0 if cond else 1

    fails += ok(non_fungal_eukaryote(SIL) == (True, u'Parakahliella macrostoma'), 'SILVA ciliate path -> outside, last named node')
    fails += ok(non_fungal_eukaryote(u'AB1.1.2 ' + SIL)[0] is True, 'accession-prefixed SILVA path')
    fails += ok(non_fungal_eukaryote(FUN) == (False, None), 'Fungi path -> not outside')
    fails += ok(non_fungal_eukaryote(REF) == (False, None), 'RefSeq title -> not outside')
    fails += ok(non_fungal_eukaryote(u'Bacteria;Pseudomonadota;Thiocapsa;Thiocapsa roseopersicina') == (False, None), 'Bacteria path -> not outside')
    fails += ok(non_fungal_eukaryote(PR2) == (True, u'Parakahliella macrostoma'), 'PR2 ciliate header')
    fails += ok(non_fungal_eukaryote(PR2F) == (False, None), 'PR2 fungal header')
    fails += ok(non_fungal_eukaryote(PR2X) == (True, u'Sandonidae'), 'PR2 placeholders skipped')
    L = {u'18S': [(99.65, 1725, 3000, 3700, PR2, 1750, u'PR2'), (87.78, 892, 800, 3700, R18, 1700, u'RefSeq_18S')],
         u'ITS': [(97.34, 488, 800, 3700, REF, 500, u'RefSeq_ITS')]}
    g = kingdom_gate(L)
    fails += ok(g and g[0] == u'outside kingdom: Parakahliella macrostoma' and g[1] == u'99.65', '18S ciliate 99.65 vs fungal 87.8 -> gate opens')
    fails += ok(kingdom_gate({u'18S': [(99.65, 1725, 0, 0, PR2, 1750, u'P'), (99.0, 1700, 0, 0, PR2F, 1700, u'P')]}) is None,
                'fungal 18S at 99 (margin < %.1f) -> shut' % KINGDOM_GATE_MARGIN)
    fails += ok(kingdom_gate({u'18S': [(94.0, 1700, 0, 0, SIL, 1750, u'S')]}) is None, 'ciliate at 94 (< %.0f) -> shut' % KINGDOM_GATE_IDENTITY)
    fails += ok(kingdom_gate({u'18S': [(99.9, 120, 0, 0, SIL, 1750, u'S')]}) is None, '120 bp under the evidence floor -> shut')
    fails += ok(kingdom_gate({u'28S': [(99.9, 1700, 0, 0, FUN, 1750, u'S')]}) is None, 'fungal record -> shut')
    fails += ok(kingdom_gate({u'ITS': [(99.9, 1700, 0, 0, SIL, 1750, u'S')]}) is None, 'ITS never opens the gate')
    fails += ok(current_name(u'Methanosaeta concilii') == u'Methanothrix soehngenii' and current_name(u'Methanothrix cf. soehngenii') == u'Methanothrix soehngenii', 'synonyms and cf. stripped')
    fails += ok(reconcile_with_nt(u'Ruminofilibacter xylanolyticum', 99.80, u'Xiashengella succiniciproducens', 'species', 99.87, 1490, 98.7)
                == (u'Ruminofilibacter cf. xylanolyticum', True), 'different species inside the margin -> cf.')
    fails += ok(reconcile_with_nt(u'Methanothrix soehngenii', 99.93, u'Methanosaeta concilii', 'species', 99.79, 1440, 98.7)
                == (u'Methanothrix soehngenii', False), 'synonym -> unchanged')
    fails += ok(reconcile_with_nt(u'Petriella musispora', 100.0, u'Lomentospora prolificans', 'species', 99.0, 1719, 99.6)
                == (u'Petriella musispora', False), 'nt below the species threshold / margin -> unchanged')
    fails += ok(reconcile_with_nt(u'Petrimonas sulfuriphila', 99.5, u'Petrimonas sp.', 'genus', 99.9, 1400, 98.7)
                == (u'Petrimonas sulfuriphila', False), 'genus-rank nt hit -> unchanged')
    fails += ok(reconcile_with_nt(u'Blastochloris tepida', 99.5, u'Blastochloris sulfoviridis', 'species', 100.0, 200, 98.7)
                == (u'Blastochloris tepida', False), 'nt hit under 250 bp -> unchanged')
    print('  %s: %d failures' % ('PASS' if not fails else 'FAIL', fails))
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
