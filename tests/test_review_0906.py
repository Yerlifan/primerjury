# -*- coding: utf-8 -*-
"""Six findings of an external code review (2026-09-06, commit 23aad2b), each kept
from coming back.

1. The local specificity layer imported a module name that does not exist in this
   repository (`yapilandirma`; here it is `config`), so the layer stopped with an
   ImportError before scanning anything.
2. The global scan checkpoint was accepted on its format version alone: changed
   primers or a changed database read the OLD result back from the same file.
3. The scan looked for the amplicon on one strand only, (F, rc R); a record that
   carries it on the other strand gave zero products.
4. `./primerjury consensus` wrote CONSENSUS_SELECTION.tsv and nothing consumed it,
   so the polished consensus that won never reached the canonical set.
5. The canonical file-name parser took the digits inside 'BIN123' for a taxid and
   dropped 'BIN1' altogether.
6. A stage stamped 'bitti' was skipped on resume without auditing its output, even
   when the output had been deleted.

RUN
    python3 tests/test_review_0906.py
"""
from __future__ import print_function
import io
import os
import shutil
import subprocess
import sys
import tempfile

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, KOK)
from screening import global_scan as gs, build_canonical as bc, engine_gateway as eg   # noqa: E402
sys.path.insert(0, os.path.join(KOK, 'verification'))
import full_chain as fc                       # noqa: E402
import specificity_round as sp                # noqa: E402

F, R = 'AGCGCAACCCTCATCATTAG', 'GGCTGCTGGCACGGAGTT'
AMPLICON = F + 'ACGT' * 30 + eg.rc(R)          # 158 bp product


def _write(path, seq, name='synthetic'):
    io.open(path, 'w').write(u'>%s\n%s\n' % (name, seq))


def test_specificity_layer_imports():
    d = tempfile.mkdtemp()
    try:
        try:
            sp.katman1_yerel(d, [], lambda *a, **k: None, d)
        except ImportError as e:
            raise AssertionError('katman1_yerel still fails to import: %s' % e)
    finally:
        shutil.rmtree(d, ignore_errors=True)
    return True


def test_scan_counts_both_strands():
    d = tempfile.mkdtemp()
    try:
        db = os.path.join(d, 'ref.fa')
        cand = [dict(ad='a0', F=F, R=R, lo=60, hi=250)]
        _write(db, AMPLICON)
        ileri = gs.tara(cand, db, max_mm=0)['a0']
        _write(db, eg.rc(AMPLICON))
        ters = gs.tara(cand, db, max_mm=0)['a0']
        assert ileri['urun'] == 1 and ters['urun'] == 1, (ileri['urun'], ters['urun'])
        assert ileri['boy'] == ters['boy'] == {158: 1}, (ileri['boy'], ters['boy'])
        # a record is counted ONCE even though both placements are scanned
        assert len(ters['vurus']) == 1
    finally:
        shutil.rmtree(d, ignore_errors=True)
    return True


def test_scan_checkpoint_is_sealed_by_its_inputs():
    d = tempfile.mkdtemp()
    try:
        db = os.path.join(d, 'ref.fa')
        cp = os.path.join(d, 'scan.pkl')
        cand = [dict(ad='a0', F=F, R=R, lo=60, hi=250)]
        _write(db, AMPLICON)
        assert gs.tara(cand, db, cp, max_mm=0)['a0']['urun'] == 1
        assert gs.tara(cand, db, cp, max_mm=0)['a0']['urun'] == 1          # the checkpoint is reused
        changed = [dict(ad='a0', F='T' * 20, R='T' * 18, lo=60, hi=250)]
        assert gs.tara(changed, db, cp, max_mm=0)['a0']['urun'] == 0, 'changed primers read the old result'
        _write(db, 'ACGT' * 100)                                            # a different database, same path
        assert gs.tara(cand, db, cp, max_mm=0)['a0']['urun'] == 0, 'a changed database read the old result'
        assert gs.tara(cand, db, cp, max_mm=1)['a0']['urun'] == 0
    finally:
        shutil.rmtree(d, ignore_errors=True)
    return True


def test_bin_ids_survive_the_canonical_parser():
    cases = {'B-1_123_consensus.fasta': 'B-1_123',
             'B-1_BIN1_consensus.fasta': 'B-1_BIN1',
             'B-1_BIN123_consensus.fasta': 'B-1_BIN123',
             'A1_1_reads_OBEK3_consensus_strict.fasta': 'A1-1_OBEK3',
             'A1-1-reads_2209_consensus_strict.fasta': 'A1-1_2209',
             'notes.fasta': None}
    for name, want in cases.items():
        got = bc.kutu_adi(name)
        assert got == want, (name, got, want)
    return True


def _sense_b_sequence(filler):
    """A synthetic B (bacterial) consensus the orientation module reads as SENSE:
    the panel's universal pair in the forward placement and two SSU motifs."""
    from screening import orientation as yon
    f_b, r_b = yon.PANEL_CIFT['B'][1], yon.PANEL_CIFT['B'][2]
    return (filler * 10 + f_b + filler * 20 + 'GTGCCAGCAGCCGCGGTAA' + filler * 20
            + 'GGATTAGATACCC' + filler * 20 + eg.rc(r_b) + filler * 20)


def test_selection_table_reaches_the_canonical_set():
    d = tempfile.mkdtemp()
    try:
        ozgun = _sense_b_sequence('ACGT')
        cila = _sense_b_sequence('AGCT')
        os.makedirs(os.path.join(d, 'consensus sequences'))
        os.makedirs(os.path.join(d, 'referans_konsensus', 'pak_polish', 'consensus'))
        _write(os.path.join(d, 'consensus sequences', 'B-1_123_consensus_strict.fasta'), ozgun, 'B-1_123')
        _write(os.path.join(d, 'referans_konsensus', 'pak_polish', 'consensus', 'B-1_123.fasta'), cila, 'B-1_123')
        baslik = u'bin\tcandidates\tCHOSEN set\tsupport %\tN %\tbp\tfewest-N set\tits support %\tcriteria disagree\ttrusted\n'
        satir = u'B-1_123\t2\treferans_konsensus/pak_polish/consensus\t97.0\t0.0\t%d\treferans_konsensus/pak_polish/consensus\t97.0\tno\tYES\n' % len(cila)
        io.open(os.path.join(d, 'CONSENSUS_SELECTION.tsv'), 'w').write(baslik + satir)
        r = subprocess.run([sys.executable, os.path.join(KOK, 'screening', 'build_canonical.py'),
                            '--root', d, '--priority', 'selection'], capture_output=True, text=True)
        assert r.returncode == 0, r.stdout[-800:] + r.stderr[-800:]
        idx = io.open(os.path.join(d, 'canonical_consensus', 'INDEX.tsv')).read().splitlines()
        assert len(idx) == 2, idx
        row = dict(zip(idx[0].split('\t'), idx[1].split('\t')))
        assert row['kaynak'] == 'referans_konsensus/pak_polish/consensus', row
        seq = u''.join(x.strip() for x in io.open(os.path.join(d, 'canonical_consensus', 'B-1_123.canonical.fa')) if not x.startswith('>'))
        assert seq == cila, 'the canonical file does not carry the chosen consensus'
        # without the table the old fixed order (ozgun first) still applies
        r = subprocess.run([sys.executable, os.path.join(KOK, 'screening', 'build_canonical.py'),
                            '--root', d, '--priority', 'ozgun', '--output', 'canon2'], capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-800:]
        seq = u''.join(x.strip() for x in io.open(os.path.join(d, 'canon2', 'B-1_123.canonical.fa')) if not x.startswith('>'))
        assert seq == ozgun
    finally:
        shutil.rmtree(d, ignore_errors=True)
    return True


def test_resume_audits_a_finished_stage_again():
    d = tempfile.mkdtemp()
    try:
        calls = []

        def denet(kok, ayar):
            calls.append('denet')
            return (False, 'output missing') if len(calls) == 1 else (True, 'ok')

        def komut_f(kok, ayar):
            calls.append('komut')
            return []
        log = []
        secili = [('S', 'stage', 'g', (1, 2), False, komut_f, denet)]
        durum = {'S': dict(durum='bitti', sure=1)}
        ayar = dict(kraken_var=False, pluspfp='')
        ok = fc.calistir(d, ayar, secili, durum, os.path.join(d, 'durum.json'), io.StringIO(), log.append, True)
        assert ok is True
        assert calls.count('denet') == 2, calls          # once on resume, once after the rerun
        assert durum['S']['durum'] == 'bitti'
        assert any('FAILED' in x for x in log), log
        # a stage whose output still passes is skipped, after exactly one audit
        calls[:] = []
        log[:] = []
        durum = {'S': dict(durum='bitti', sure=1)}
        secili = [('S', 'stage', 'g', (1, 2), False, komut_f, lambda k, a: (calls.append('denet'), (True, 'ok'))[1])]
        fc.calistir(d, ayar, secili, durum, os.path.join(d, 'durum.json'), io.StringIO(), log.append, True)
        assert calls == ['denet'] and any('SKIPPED' in x for x in log), (calls, log)
    finally:
        shutil.rmtree(d, ignore_errors=True)
    return True


def main():
    tests = [test_specificity_layer_imports, test_scan_counts_both_strands,
             test_scan_checkpoint_is_sealed_by_its_inputs, test_bin_ids_survive_the_canonical_parser,
             test_selection_table_reaches_the_canonical_set, test_resume_audits_a_finished_stage_again]
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
