# -*- coding: utf-8 -*-
"""Two faults that passed every syntax check (2026-09-05).

1. The dominant-allele step aligned reads against the IUPAC-coded consensus.
   minimap2 cannot seed on IUPAC letters, so a mixed bin with a 16 per cent
   IUPAC template aligned 0 of 3,001 reads. The template must be plain ACGT.
2. The naming decision recursed forever when the best hit was unnamed and a
   named hit sat inside the margin: the recursive call re-sorted the list and
   put the unnamed record back on top.

RUN
    python3 tests/test_consensus_template.py
"""
from __future__ import print_function
import importlib.util
import os
import sys

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, KOK)


def _load_dominant():
    path = os.path.join(KOK, 'steps', 'dominant_allele_consensus.py')
    spec = importlib.util.spec_from_file_location('dominant_allele_consensus', path)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except SystemExit:
        return None          # no aligner backend on this machine
    return mod


def test_template_is_plain_acgt():
    mod = _load_dominant()
    if mod is None:
        print('  template test skipped: no aligner backend')
        return True
    core = 'ACGTRYSWKMNACGT'
    ref = 'ACGTAGTCAGTACGT'     # same length: IUPAC/N take the reference base
    t, src = mod.alignment_template(core, ref)
    ok = (src == 'reference' and t == ref
          and all(c in 'ACGT' for c in t) and len(t) == len(core))
    # different length: first base of the code, N stays
    t2, src2 = mod.alignment_template(core, 'ACGT')
    ok = ok and src2 == 'first_base' and t2 == 'ACGTACCAGANACGT' and len(t2) == len(core)
    print('  template: %s' % ('ok' if ok else 'FAIL (%s %s / %s %s)' % (src, t, src2, t2)))
    return ok


def test_unnamed_best_hit_does_not_recurse():
    sys.setrecursionlimit(200)
    from verification import identity_verification as iv
    hits = [
        dict(baslik='JX501314.1.600 uncultured bacterium', kimlik=99.9, hiz_uzunluk=1400),
        dict(baslik='NR_043700.1 Oceanobacillus chironomi strain T3944D', kimlik=99.8, hiz_uzunluk=1400),
        dict(baslik='KC123456.1 uncultured organism clone', kimlik=99.7, hiz_uzunluk=1400),
        dict(baslik='NR_113330.1 Oceanobacillus indicireducens strain A21', kimlik=97.0, hiz_uzunluk=1400),
    ]
    try:
        r = iv.savunulabilir_duzey(hits, 'SSU')
    except RecursionError:
        print('  recursion: FAIL (RecursionError)')
        return False
    ok = 'Oceanobacillus' in (r.get('onerilen_ad') or '')
    print('  recursion: %s (%s | %s)' % ('ok' if ok else 'FAIL', r.get('duzey'), r.get('onerilen_ad')))
    return ok


if __name__ == '__main__':
    results = [test_template_is_plain_acgt(), test_unnamed_best_hit_does_not_recurse()]
    print('consensus template / naming recursion: %s' % ('PASS' if all(results) else 'FAIL'))
    sys.exit(0 if all(results) else 1)
