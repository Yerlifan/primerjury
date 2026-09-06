# -*- coding: utf-8 -*-
"""The fungal bin module: population split, medoid template, polished consensus.

Why it exists (2026-09-05): in a mixed fungal bin the blended consensus sat at
97.25 per cent to the type record while single reads reached 100. The module must
(1) pick the dominant genus population from per-read best hits, (2) choose the
template from the reads themselves (reverse-strand reads included), (3) polish a
consensus that beats the medoid and reaches the species range on synthetic data,
(4) sample deterministically, (5) write rows that match the header.

Pure-python parts always run; the tool-bound part (blastn, minimap2, samtools)
runs when the tools are present and is skipped otherwise, saying so.

RUN
    python3 tests/test_fungal_bin_identity.py
"""
from __future__ import print_function
import collections
import os
import random
import sys

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(KOK, 'verification'))
import fungal_bin_identity as fbi   # noqa: E402
import population_polish as pp      # noqa: E402


def test_medoid_prefers_central_full_length_read():
    rnd = random.Random(3)
    base = u''.join(rnd.choice(u'ACGT') for _ in range(1200))
    reads = [(u'r%d' % i, fbi._mutate(rnd, base, 0.02)) for i in range(12)]
    reads.append((u'short', base[:500]))                       # too short to be the template
    reads.append((u'far', u''.join(rnd.choice(u'ACGT') for _ in range(1200))))   # unrelated
    reads[3] = (reads[3][0], fbi.reverse_complement(reads[3][1]))   # a reverse-strand read
    m = fbi.medoid(reads)
    assert m is not None
    assert m[0] not in (u'short', u'far'), m[0]
    assert m[2] > 0.3, m[2]
    return True


def test_populations_group_by_genus_and_skip_unnamed():
    best = {
        u'r1': (99.0, 500, u'NR_1 Petriella musispora CBS 1 ITS', 'db'),
        u'r2': (98.0, 500, u'NR_2 Petriella guttulata CBS 2 ITS', 'db'),
        u'r3': (97.0, 500, u'NR_3 Scedosporium boydii ITS', 'db'),
        u'r4': (96.0, 500, u'NR_4 uncultured fungus clone X ITS', 'db'),
        u'r5': (96.0, 500, u'NR_5 Petriella sp. CBS 3 ITS', 'db'),
    }
    count, members, scount, spid = fbi.populations(best)
    assert count[u'Petriella'] == 3 and count[u'Scedosporium'] == 1, dict(count)
    assert u'uncultured' not in u' '.join(count), dict(count)
    assert scount[u'Petriella musispora'] == 1 and u'Petriella sp' not in u' '.join(scount), dict(scount)
    assert len(members[u'Petriella']) == 3
    return True


def test_best_hits_keeps_highest_identity_then_length():
    def fake_blast(query, db, threads=1, top=50):
        return {u'q': [(97.0, 600, 900.0, 3000, u'NR_1 Aus alpha ITS', 600),
                       (99.0, 300, 500.0, 3000, u'NR_2 Bus beta ITS', 300),
                       (99.0, 550, 950.0, 3000, u'NR_3 Cus gamma ITS', 550),
                       (100.0, 100, 190.0, 3000, u'NR_4 Dus delta ITS', 100)]}   # too short
    best = fbi.best_hits('x.fa', ['db'], None, 1, blast_fn=fake_blast)
    assert best[u'q'][2].startswith(u'NR_3'), best
    return True


def test_row_matches_header():
    r = dict(bin=u'F1-1_1', group=u'F1', locus=u'ITS', reads=10, assigned=8, genus=u'Petriella', share=75.0,
             genus2=u'Scedosporium', share2=25.0)
    assert len(pp.make_row(r).split(u'\t')) == len(pp.HEADER)
    r = dict(bin=u'F1-1_2', group=u'F1', reads=2, note=u'too few reads (2)')
    assert len(pp.make_row(r).split(u'\t')) == len(pp.HEADER)
    return True


def test_bin_files_ignore_stray_copies(tmp=None):
    import tempfile, shutil, io
    d = tempfile.mkdtemp()
    try:
        for lib, f in ((u'F2-1', u'F2-1-reads_101201.fastq'), (u'F2-1', u'F1-4-reads_101201.fastq'),
                       (u'F1-1', u'F1_1_reads_OBEK3.fastq')):
            os.makedirs(os.path.join(d, 'fastq files', lib), exist_ok=True)
            io.open(os.path.join(d, 'fastq files', lib, f), 'w').write(u'')
        files = fbi.bin_files(d)
        assert sorted(files) == [u'F1-1_OBEK3', u'F2-1_101201'], sorted(files)
    finally:
        shutil.rmtree(d, ignore_errors=True)
    return True


def test_tool_bound_self_test():
    if not fbi.tools_present():
        print('  tool-bound self-test skipped: blastn / minimap2 / samtools not all present')
        return True
    return fbi.self_test() == 0


def main():
    tests = [test_medoid_prefers_central_full_length_read,
             test_populations_group_by_genus_and_skip_unnamed,
             test_best_hits_keeps_highest_identity_then_length,
             test_row_matches_header,
             test_bin_files_ignore_stray_copies,
             test_tool_bound_self_test]
    failed = 0
    for t in tests:
        try:
            ok = t()
        except Exception as e:      # noqa: BLE001 - a test harness reports, it does not hide
            ok = False
            print('  %s raised %s: %s' % (t.__name__, type(e).__name__, e))
        print('  %-52s %s' % (t.__name__, 'PASS' if ok else 'FAIL'))
        failed += 0 if ok else 1
    print('  %d of %d tests passed' % (len(tests) - failed, len(tests)))
    return 0 if not failed else 1


if __name__ == '__main__':
    sys.exit(main())
