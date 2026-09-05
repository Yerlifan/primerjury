# -*- coding: utf-8 -*-
"""The consensus taxonomy behind the QIIME2 route (2026-09-05).

1. ranks() stops at the first uninformative name (uncultured, sp., incertae sedis).
2. consensus() keeps a rank only while at least 51 per cent agree.
3. best_hits() takes one reference only (the one with the best identity), pools
   hits within one point, and breaks equal identity by the LONGER alignment.
4. UNITE targets carry their taxonomy in the name.

RUN
    python3 tests/test_qiime2_taxonomy.py
"""
from __future__ import print_function
import importlib.util
import os
import sys

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load():
    spec = importlib.util.spec_from_file_location(
        'qiime2_reference_taxonomy', os.path.join(KOK, 'steps', 'qiime2_reference_taxonomy.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    m = _load()
    ok = m.ranks('k__Fungi;p__Ascomycota;g__Petriella;s__Petriella_musispora') == \
        ['Fungi', 'Ascomycota', 'Petriella', 'Petriella musispora']
    ok = ok and m.ranks('Bacteria;Firmicutes;uncultured bacterium;Bacillus') == ['Bacteria', 'Firmicutes']
    ok = ok and m.ranks('Fungi;Ascomycota;Incertae sedis;X') == ['Fungi', 'Ascomycota']
    print('  ranks: %s' % ('ok' if ok else 'FAIL'))
    c = m.consensus([['A', 'B', 'C'], ['A', 'B', 'D'], ['A', 'B', 'C']])
    ok2 = c == 'A;B;C'
    c2 = m.consensus([['A', 'B', 'C'], ['A', 'B', 'D']])
    ok2 = ok2 and c2 == 'A;B'                       # 50 per cent is not agreement
    print('  consensus: %s (%s | %s)' % ('ok' if ok2 else 'FAIL', c, c2))
    hits = [(99.0, 1400, 't1', 'silva.fasta'), (99.0, 600, 't2', 'silva.fasta'),
            (98.3, 1400, 't3', 'silva.fasta'), (97.0, 1400, 't4', 'silva.fasta'),
            (99.5, 300, 'u1', 'unite.fasta')]
    good, top = m.best_hits(hits)
    ok3 = top == 99.5 and [x[2] for x in good] == ['u1']    # the best reference wins outright
    good2, top2 = m.best_hits(hits[:4])
    ok3 = ok3 and top2 == 99.0 and [x[2] for x in good2] == ['t1', 't2', 't3']  # within one point, longer first
    print('  best hits: %s' % ('ok' if ok3 else 'FAIL'))
    ok4 = m.embedded_taxonomy('KJ734967|k__Fungi;p__Ascomycota|SH123') == 'k__Fungi;p__Ascomycota'
    ok4 = ok4 and m.embedded_taxonomy('NR_043700.1') == ''
    print('  embedded taxonomy: %s' % ('ok' if ok4 else 'FAIL'))
    r = all([ok, ok2, ok3, ok4])
    print('qiime2 taxonomy: %s' % ('PASS' if r else 'FAIL'))
    return 0 if r else 1


if __name__ == '__main__':
    sys.exit(main())
