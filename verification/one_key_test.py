# -*- coding: utf-8 -*-
"""THE SCENARIO TESTS of one_key.py. Not one of them touches the real directory.

WHY A SHADOW ROOT
-----------------
The tests have to try a stage FAILING and a stage BEING LEFT HALF DONE. Doing
that in the real directory would corrupt the tables an ordering decision rests
on. So every scenario runs in a SHADOW ROOT set up under the temporary
directory:
  * the heavy, read only directories are SYMLINKED rather than copied
  * the stage scripts are replaced by FAKE ones whose behaviour we set
  * every write stays inside the temporary directory
Not one byte is written into the mounted directory. The shadow root is left
behind after the run.

To run it:  python3 verification/one_key_test.py --root .
"""

import os, sys, io, json, time, shutil, argparse, subprocess, tempfile

BU = os.path.dirname(os.path.abspath(__file__))
TEK_TUS = os.path.join(BU, 'one_key.py')

# The fake script template. The exit code and the files it produces come from
# outside, so "the stage failed", "the output is empty" and "no output was
# produced at all" can each be arranged without touching this file.
SAHTE = u'''# -*- coding: utf-8 -*-
import os, sys, time
KOK = os.path.dirname(os.path.abspath(__file__))
while KOK and not os.path.exists(os.path.join(KOK, '_SINAMA_AYAR.json')):
    ust = os.path.dirname(KOK)
    if ust == KOK:
        break
    KOK = ust
import json
ayar = json.load(open(os.path.join(KOK, '_SINAMA_AYAR.json'), encoding='utf-8'))
ad = %(ad)r
# screening/__main__.py IKI asama tarafindan cagriliyor: 8 (--selftest) ve
# S (--mode ozet). Hangisi oldugunu argv soyler; sahte betik de ayni ayrimi
# yapmali, yoksa S asamasi hicbir zaman ciktisini uretmez.
if ad == 'sina' and '--mode' in sys.argv:
    ad = 'S'
d = ayar.get(ad, {})
print('[sahte %%s] basladi' %% ad)
sys.stdout.flush()
if d.get('bekle'):
    for i in range(int(d['bekle'])):
        print('[sahte %%s] calisiyor %%d' %% (ad, i)); sys.stdout.flush()
        time.sleep(1)
for y, icerik in (d.get('yaz') or {}).items():
    t = os.path.join(KOK, y)
    os.makedirs(os.path.dirname(t), exist_ok=True)
    open(t, 'w', encoding='utf-8').write(icerik)
    print('[sahte %%s] yazdi: %%s' %% (ad, y))
if d.get('selftest_metni'):
    print(d['selftest_metni'])
print('[sahte %%s] bitti, cikis kodu %%s' %% (ad, d.get('rc', 0)))
sys.exit(int(d.get('rc', 0)))
'''

# (the real path, the fake script name). These must be EXACTLY the same as the
# 'betik' fields in one_key.py's stage schedule.
BETIKLER = [
    ('screening/__main__.py', 'sina'),
    ('verification/quick_consistency_test.py', 'H'),
    ('verification/access_check.py', 'E'),
    ('engine/rederive_membership.py', 'U'),
    ('protocol/single_protocol_measure.py', 'P'),
    ('verification/recovery_round.py', 'K'),
    ('verification/specificity_round.py', 'D'),
    ('verification/identity_verification.py', 'I'),
    ('verification/all_bin_identities.py', 'G'),
    # MEASURED: stage N (verification/audit_all.py) was added on 2026-08-10 but was
    # never written into this list. The pre-check treats that script as REQUIRED, so
    # it always came out missing in the shadow root, the run never started, and every
    # scenario test silently returned nothing.
    ('verification/audit_all.py', 'N'),
]

# The directories the pre-check treats as REQUIRED. All of them must EXIST in the
# shadow root even if empty; otherwise the pre-check stops, correctly, and we
# never reach the test we came for. (That the pre-check DOES stop is tested
# separately by scenario S4.)
KLASORLER = ['fastq files', 'consensus sequences', 'final_primers', 'REFERENCE_DB',
             'screening', 'protocol', 'engine', 'steps',
             'canonical_consensus', 'tools', 'verification']

MFE_IX = ['archaea.16S.fna', 'bacteria.16S.fna', 'fungi.ITS.fna',
          'fungi.28SrRNA.fna', 'fungi.18SrRNA.fna', 'SILVA_138.2_SSURef_NR99.fasta']
KUMELER = MFE_IX + ['SILVA_138.2_LSURef_NR99.fasta', 'UNITE_ITS.fasta',
                    'PR2_SSU_taxo_long.fasta', 'ROD_v1.2_operon_variants.fasta',
                    'ref_all2.fna']


def golge_kur(taban):
    u"""Bos ama on kontrolu gecebilen bir golge kok kurar."""
    if os.path.exists(taban):
        shutil.rmtree(taban)
    os.makedirs(taban)
    for k in KLASORLER:
        os.makedirs(os.path.join(taban, k), exist_ok=True)
        # the directories must not be EMPTY (the pre-check prints an item count; empty
        # passes too, but let it look like the real thing)
        with io.open(os.path.join(taban, k, '_yer_tutucu.txt'), 'w',
                     encoding='utf-8') as fh:
            fh.write(u'sinama\n')
    # screening is called AS A PACKAGE (python3 -m screening)
    with io.open(os.path.join(taban, 'screening', '__init__.py'), 'w',
                 encoding='utf-8') as fh:
        fh.write(u'')
    # mfeprimer ikilisi
    mp = os.path.join(taban, 'tools', 'mfeprimer')
    with io.open(mp, 'w', encoding='utf-8') as fh:
        fh.write(u'#!/bin/sh\nexit 0\n')
    os.chmod(mp, 0o755)
    # SILVA: DNA alfabesi kapisini gecen kucuk bir sahte fasta
    with io.open(os.path.join(taban, 'REFERENCE_DB',
                              'SILVA_138.2_SSURef_NR99.fasta'), 'w',
                 encoding='utf-8') as fh:
        for i in range(50):
            fh.write(u'>sahte%d test\nACGTACGTACGTTTTTACGTACGTACGTTTTT\n' % i)
    for f in KUMELER:
        p = os.path.join(taban, 'REFERENCE_DB', f)
        if not os.path.exists(p):
            with io.open(p, 'w', encoding='utf-8') as fh:
                fh.write(u'>sahte\nACGT\n')
    for f in MFE_IX:
        with io.open(os.path.join(taban, 'REFERENCE_DB', f + '.primerqc.bin'),
                     'wb') as fh:
            fh.write(b'0' * 100)
    # Sahte asama betikleri
    for yol, ad in BETIKLER:
        t = os.path.join(taban, yol)
        os.makedirs(os.path.dirname(t), exist_ok=True)
        with io.open(t, 'w', encoding='utf-8') as fh:
            fh.write(SAHTE % dict(ad=ad))
    # The final verdict table (the order table in the summary reads this).
    # The same shape as the real one: comment lines plus a YENI_HUKUM column.
    # Three rows go into the order (SIPARIS EDILEBILIR x2 + KOSULLU x1) and two
    # do not (ESIK ALTI + ONERILMEZ). The expected count is 3.
    with io.open(os.path.join(taban, 'ESIK_VE_OLCUT_2026-08-08.tsv'), 'w',
                 encoding='utf-8') as fh:
        fh.write(u'# test file\n')
        fh.write(u'hedef\tESKI_HUKUM\tYENI_HUKUM\tdCq_olculen\n')
        fh.write(u'A_hedefi\tKOSULLU\tSIPARIS EDILEBILIR (kosullu)\t7,28\n')
        fh.write(u'B_hedefi\tKOSULLU\tSIPARIS EDILEBILIR (kosullu)\t4,01\n')
        fh.write(u'C_universal\tRISKLI\tKOSULLU (control primer)\t-\n')
        fh.write(u'D_hedefi\tKOSULLU\tESIK ALTI (kosullu, silinmedi)\t-0,43\n')
        fh.write(u'E_target\tRISKLI\tONERILMEZ (the target is absent from the sample)\t-\n')
    return taban


def ayar_yaz(taban, ayar):
    with io.open(os.path.join(taban, '_SINAMA_AYAR.json'), 'w',
                 encoding='utf-8') as fh:
        fh.write(json.dumps(ayar, ensure_ascii=False, indent=1))


def basarili_ayar():
    'The default behaviour, which makes every stage succeed.'
    return {
        'sina': dict(rc=0, selftest_metni=u'TUM SINAMALAR GECTI.'),
        # N: the audit stage. Its fake script has to write its own report,
        # otherwise the output check says DUSTU and the chain stops at N.
        'N': dict(rc=0, yaz={'ONE_KEY_RESULT/DENETIM_RAPORU.md':
                             '# the audit report  A fake report, for the '
                             'test.'}),
        'H': dict(rc=0, yaz={'QUICK_TEST/QUICK_TEST_REPORT.md':
                             '# the report  ZINCIR TUTARLI (against its own '
                             'reference)'}),
        'E': dict(rc=0, yaz={'ACCESS_RESULT/erisim_dogrulama.tsv':
                             u'vtb\tsonuc\narchaea\tGECTI\n'}),
        'U': dict(rc=0, yaz={'uyelik_yeniden_turetme_uyelik_29991231.tsv':
                             u'kutu\thedef\nA1\tX\n'}),
        'P': dict(rc=0, yaz={'ONE_PROTOCOL_RESULT/panel_tek_protokol.tsv':
                             u'hedef\tayrim\nX\t9\n',
                             'ONE_PROTOCOL_RESULT/SIPARIS_LISTESI.tsv':
                             u'hedef\tSINIF\tF\tR\nX\tKESIN\tAAA\tTTT\n'}),
        'K': dict(rc=0, yaz={'RECOVERY_RESULT/kurtarma_satirlari.tsv':
                             u'hedef\tgecti\nX\tEVET\n'}),
        'D': dict(rc=0, yaz={'VERIFICATION_RESULT/dogrulama_uc_sutun.tsv':
                             u'hedef\tKARAR\nX\tKOSULLU\n'}),
        'I': dict(rc=0, yaz={'IDENTITY_RESULT/kimlik_iddialari.tsv':
                             u'iddia\tsonuc\n1\tDOGRULANDI\n'}),
        'G': dict(rc=0, yaz={'ALL_IDENTITIES_RESULT/tum_kutu_kimlikleri.tsv':
                             u'kutu\tkimlik\nA1\tX\n'}),
        'S': dict(rc=0, yaz={'SCREENING_RESULT/00_OZET_HEPSI.md':
                             '# the summary  A fake summary file, for the '
                             'test.'}),
    }


def kos(taban, ek=(), zaman_asimi=180, sinyal_sn=None):
    argv = [sys.executable, TEK_TUS, '--root', taban, '--confirm',
            '--ncbi', 'yok', '--liveness', '3'] + list(ek)
    p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if sinyal_sn:
        # KESINTI SENARYOSU: belirlenen saniyede SIGINT gonderilir; bu,
        # Windows'ta Ctrl+C'nin surece ulasmasinin Linux karsiligidir.
        import signal as _s, threading as _t
        def kes():
            time.sleep(sinyal_sn)
            try:
                p.send_signal(_s.SIGINT)
            except Exception:
                pass
        th = _t.Thread(target=kes)
        th.daemon = True
        th.start()
    out = p.communicate(timeout=zaman_asimi)[0].decode('utf-8', 'replace')
    return p.returncode, out


SONUC = []


def sina(ad, kosul, ayrinti=u''):
    SONUC.append((ad, bool(kosul), ayrinti))
    print(u'  [%s] %s%s' % (u'GECTI' if kosul else u'DUSTU', ad,
                            (u'   -> ' + ayrinti) if ayrinti else u''))
    return bool(kosul)


def durum_oku(taban):
    y = os.path.join(taban, 'ONE_KEY_RESULT', 'durum.json')
    if not os.path.exists(y):
        return {}
    return json.load(io.open(y, encoding='utf-8'))


# ===========================================================================
def s1_sifirdan(ana):
    print(u'\n--- S1: A RUN FROM SCRATCH (every stage succeeds) ---')
    t = golge_kur(os.path.join(ana, 's1'))
    ayar_yaz(t, basarili_ayar())
    rc, out = kos(t)
    d = durum_oku(t)
    sina(u'S1 cikis kodu 0', rc == 0, u'rc=%d' % rc)
    sina(u'S1 the pre-check passed', u'PRE-CHECK PASSED' in out)
    for k in ('8', 'N', 'H', 'E', 'U', 'P', 'K', 'D', 'I', 'G', 'S'):
        sina(u'S1 stage %s finished' % k, d.get(k, {}).get('durum') == 'bitti',
             d.get(k, {}).get('durum', 'YOK'))
    sina(u'S1 the Kraken stages were skipped for want of the tool',
         all(d.get(k, {}).get('durum', '').startswith('atlandi') for k in 'WXZ'),
         u'W=%s X=%s Z=%s' % tuple(d.get(k, {}).get('durum', '?') for k in 'WXZ'))
    sina(u'S1 the order is right: U ran before P',
         out.index(u'>> [') >= 0 and
         out.find(u'\n>> [') >= 0 and
         out.find(u'] U ') < out.find(u'] P '),
         u'U konumu %d, P konumu %d' % (out.find(u'] U '), out.find(u'] P ')))
    sina(u'S1 the order is right: P, K, D in that order',
         out.find(u'] P ') < out.find(u'] K ') < out.find(u'] D '))
    sina(u'S1 the order is right: I before G',
         out.find(u'] I ') < out.find(u'] G '))
    ozet = os.path.join(t, 'ONE_KEY_RESULT', '00_SABAH_OZETI.md')
    sina(u'S1 the morning summary was produced', os.path.exists(ozet))
    m = io.open(ozet, encoding='utf-8').read() if os.path.exists(ozet) else u''
    sina(u'S1 the summary says NO STAGE FAILED', u'NO STAGE FAILED' in m)
    sina(u'S1 the summary holds the order table',
         u'ORDERABLE' in m and u'| X |' in m)
    # MEASURED: these two expectations WENT STALE. In one_key.py's order source
    # list the run's own ONE_PROTOCOL_RESULT/SIPARIS_LISTESI.tsv was moved to the
    # front (deliberately: the order before it read a stale file and misled the
    # first table anyone looks at in the morning). The fake P stage writes ONE
    # KESIN row to that file, so the right count is 1 in and 0 out.
    sina(u'S1 the order count is right (1 in, 0 out)',
         u'ORDERABLE: 1 pairs = 2 oligos' in m,
         u'in the summary: %s' % (u' | '.join(
             x.strip() for x in m.splitlines() if u'ORDERABLE' in x) or u'nothing'))
    sina(u'S1 the order table came from the run own SIPARIS_LISTESI',
         u'SIPARIS_LISTESI.tsv' in m and u'SINIF' in m)
    sina(u'S1 the log file is time stamped',
         any(x.startswith('gunluk_') for x in
             os.listdir(os.path.join(t, 'ONE_KEY_RESULT'))))
    return t


def s2_devam(ana, s1_taban):
    print(u'\n--- S2: A SECOND RUN - finished stages must be skipped (0 s) ---')
    t = s1_taban
    t0 = time.time()
    rc, out = kos(t)
    gecen = time.time() - t0
    d = durum_oku(t)
    sina(u'S2 cikis kodu 0', rc == 0, u'rc=%d' % rc)
    sina(u'S2 P was skipped (the checkpoint is valid)',
         u'>> P' in out and u'SKIPPED' in out.split(u'>> P')[1][:400],
         out.split(u'>> P')[1][:120].replace(u'\n', u' ') if u'>> P' in out else u'')
    sina(u'S2 everything finished quickly on the second run (< 40 s)', gecen < 40,
         u'%.1f sn' % gecen)
    sina(u'S2 durum.json korundu', d.get('D', {}).get('durum') == 'bitti')
    return t


def s3_bayat(ana, taban):
    print(u'\n--- S3: A STALE CHECKPOINT - when the input is newer than the output ---')
    # We touch one of D's inputs: target_clades.tsv. This is the D-9 fault of
    # 2026-08-07 itself: the index was refreshed and the checkpoint read the old
    # zeros back.
    hk = os.path.join(taban, 'screening', 'target_clades.tsv')
    with io.open(hk, 'w', encoding='utf-8') as fh:
        fh.write(u'hedef\tklad\nX\tBacteria\n')
    os.utime(hk, (time.time() + 5, time.time() + 5))   # ciktidan KESIN yeni
    rc, out = kos(taban, ek=['--plan'])
    sina(u'S3 D was marked BAYAT',
         u'D  KOSACAK' in out and u'STALE' in out.split(u'D  KOSACAK')[1][:300],
         out.split(u'D  KOSACAK')[1][:150].replace(u'\n', u' ')
         if u'D  KOSACAK' in out else u'there is no D KOSACAK row')
    sina(u'S3 the stages that are not stale are still SKIPPED',
         u'I  SKIPPED' in out and u'G  SKIPPED' in out)
    rc, out = kos(taban)
    d = durum_oku(taban)
    sina(u'S3 it ran again and finished', d.get('D', {}).get('durum') == 'bitti'
         and d.get('D', {}).get('sure', 0) > 0, u'sure=%s' % d.get('D', {}).get('sure'))
    return taban


def s4_eksik_dosya(ana):
    print(u'\n--- S4: A REQUIRED FILE IS MISSING - the pre-check MUST STOP ---')
    t = golge_kur(os.path.join(ana, 's4'))
    ayar_yaz(t, basarili_ayar())
    os.remove(os.path.join(t, 'REFERENCE_DB',
                           'SILVA_138.2_SSURef_NR99.fasta.primerqc.bin'))
    shutil.rmtree(os.path.join(t, 'fastq files'))
    rc, out = kos(t)
    sina(u'S4 cikis kodu 2 (on kontrol kapisi)', rc == 2, u'rc=%d' % rc)
    sina(u'S4 it wrote PRE-CHECK FAILED', u'PRE-CHECK FAILED' in out)
    sina(u'S4 it named the missing SILVA index',
         u'SILVA_138.2_SSURef_NR99.fasta' in out and u'MFE index' in out)
    sina(u'S4 it named the missing fastq directory', u'fastq files' in out)
    sina(u'S4 it wrote the install command', u'mfeprimer index -i' in out)
    sina(u'S4 NO STAGE RAN',
         u'>> [1/' not in out and not os.path.exists(
             os.path.join(t, 'ONE_PROTOCOL_RESULT', 'panel_tek_protokol.tsv')))
    return t


def s5_asama_dustu(ana):
    print(u'\n--- S5: A STAGE FAILED (rc != 0) - the dependents MUST NOT RUN ---')
    t = golge_kur(os.path.join(ana, 's5'))
    a = basarili_ayar()
    a['P'] = dict(rc=3, yaz={  # the output IS PRODUCED but the exit code IS NOT ZERO.
        'ONE_PROTOCOL_RESULT/panel_tek_protokol.tsv': u'hedef\tayrim\nX\t9\n',
        'ONE_PROTOCOL_RESULT/SIPARIS_LISTESI.tsv': u'hedef\tSINIF\tF\tR\nX\tKESIN\tA\tT\n'})
    ayar_yaz(t, a)
    rc, out = kos(t)
    d = durum_oku(t)
    sina(u'S5 the chain exit code is 3', rc == 3, u'rc=%d' % rc)
    sina(u'S5 P FAILED (rc!=0 was not masked even though the output was not empty)',
         d.get('P', {}).get('durum') == 'DUSTU', d.get('P', {}).get('durum', 'YOK'))
    sina(u'S5 dusme sebebi cikis kodunu yaziyor',
         'EXIT CODE 3' in (d.get('P', {}).get('sebep') or ''),
         (d.get('P', {}).get('sebep') or '')[:90])
    sina(u'S5 K was skipped (dependent)', d.get('K', {}).get('durum') == 'atlandi (bagimli)',
         d.get('K', {}).get('durum', 'YOK'))
    sina(u'S5 D was skipped (dependent)', d.get('D', {}).get('durum') == 'atlandi (bagimli)',
         d.get('D', {}).get('durum', 'YOK'))
    sina(u'S5 BAGIMSIZ asamalar (I, G) yine de kostu',
         d.get('I', {}).get('durum') == 'bitti' and d.get('G', {}).get('durum') == 'bitti',
         u'I=%s G=%s' % (d.get('I', {}).get('durum'), d.get('G', {}).get('durum')))
    m = io.open(os.path.join(t, 'ONE_KEY_RESULT', '00_SABAH_OZETI.md'),
                encoding='utf-8').read()
    sina(u'S5 the summary opened a FAILED STAGES section', u'FAILED STAGES' in m)
    sina(u'S5 the summary carries the word DUSTU on the P row', u'DUSTU' in m)
    return t


def s6_bos_cikti(ana):
    print(u'\n--- S6: EXIT CODE 0 BUT EMPTY OUTPUT - must NOT count as "done" ---')
    t = golge_kur(os.path.join(ana, 's6'))
    a = basarili_ayar()
    a['K'] = dict(rc=0, yaz={'RECOVERY_RESULT/kurtarma_satirlari.tsv':
                             u'hedef\tgecti\n'})     # yalniz baslik: 0 veri satiri
    ayar_yaz(t, a)
    rc, out = kos(t)
    d = durum_oku(t)
    sina(u'S6 K FAILED (an empty output did not count as passing)',
         d.get('K', {}).get('durum') == 'DUSTU', d.get('K', {}).get('durum', 'YOK'))
    sina(u'S6 the reason says the output is EMPTY',
         u'EMPTY' in (d.get('K', {}).get('sebep') or u''),
         (d.get('K', {}).get('sebep') or u'')[:90])
    sina(u'S6 D was skipped (dependent)', d.get('D', {}).get('durum') == 'atlandi (bagimli)')
    sina(u'S6 the chain exit code is 3', rc == 3, u'rc=%d' % rc)
    return t


def s7_kesinti(ana):
    print(u'\n--- S7: INTERRUPTED (Ctrl+C) AND RESUMED WITH THE SAME COMMAND ---')
    t = golge_kur(os.path.join(ana, 's7'))
    a = basarili_ayar()
    a['D'] = dict(rc=0, bekle=25, yaz={'VERIFICATION_RESULT/dogrulama_uc_sutun.tsv':
                                       u'hedef\tKARAR\nX\tKOSULLU\n'})
    ayar_yaz(t, a)
    rc, out = kos(t, sinyal_sn=14)
    d = durum_oku(t)
    sina(u'S7 kesme cikis kodu 130', rc == 130, u'rc=%d' % rc)
    sina(u'S7 INTERRUPT RECEIVED was written', u'INTERRUPT RECEIVED' in out)
    sina(u'S7 D was stamped as interrupted', d.get('D', {}).get('durum') == 'kesildi',
         d.get('D', {}).get('durum', 'YOK'))
    sina(u'S7 the stages BEFORE D are stamped finished',
         all(d.get(k, {}).get('durum') == 'bitti' for k in ('H', 'E', 'U', 'P', 'K')),
         u', '.join(u'%s=%s' % (k, d.get(k, {}).get('durum')) for k in 'HEUPK'))
    sina(u'S7 the morning summary was written on an interruption too',
         os.path.exists(os.path.join(t, 'ONE_KEY_RESULT', '00_SABAH_OZETI.md')))
    # RESUME: the same command is run again
    a['D'] = dict(rc=0, bekle=0, yaz={'VERIFICATION_RESULT/dogrulama_uc_sutun.tsv':
                                      u'hedef\tKARAR\nX\tKOSULLU\n'})
    ayar_yaz(t, a)
    t0 = time.time()
    rc2, out2 = kos(t)
    gecen = time.time() - t0
    d2 = durum_oku(t)
    sina(u'S7-devam cikis kodu 0', rc2 == 0, u'rc=%d' % rc2)
    sina(u'S7-continue the finished stages were SKIPPED',
         u'SKIPPED' in out2 and u'>> [1/' in out2)
    sina(u'S7-continue D is finished now', d2.get('D', {}).get('durum') == 'bitti',
         d2.get('D', {}).get('durum', 'YOK'))
    sina(u'S7-continue it finished quickly (< 40 s)', gecen < 40, u'%.1f sn' % gecen)
    return t


def s8_belirlenimci_imza(ana, taban):
    print(u'\n--- S8: IS THE CHECKPOINT KEY DETERMINISTIC (md5) ---')
    # The signature computed in two separate processes from the same inputs must be
    # THE SAME. Had Python's hash() been used it would have differed under
    # PYTHONHASHSEED.
    kod = (u'import sys,os,json;sys.path.insert(0,%r);'
           u'import one_key as T;'
           u'a=[x for x in T.ASAMALAR(dict(ncbi="yok",karac=None)) if x["kod"]=="D"][0];'
           u'print(T.girdi_imzasi(%r,a))' % (os.path.dirname(TEK_TUS), taban))
    im = []
    for tohum in ('0', '1', '12345'):
        cev = dict(os.environ, PYTHONHASHSEED=tohum)
        r = subprocess.run([sys.executable, '-c', kod], capture_output=True, env=cev)
        im.append(r.stdout.decode().strip() or r.stderr.decode()[-200:])
    sina(u'S8 the signature is THE SAME under three different PYTHONHASHSEED values',
         len(set(im)) == 1 and len(im[0]) == 16, u' / '.join(im))
    return taban


def s9_yalniz_atla(ana):
    print(u'\n--- S9: the --only and --skip filters ---')
    t = golge_kur(os.path.join(ana, 's9'))
    ayar_yaz(t, basarili_ayar())
    rc, out = kos(t, ek=['--only', '8S'])
    d = durum_oku(t)
    sina(u'S9 --only 8S: only 8 and S were recorded',
         set(d.keys()) <= {'8', 'S'} and '8' in d, u'anahtarlar: %s' % sorted(d.keys()))
    t2 = golge_kur(os.path.join(ana, 's9b'))
    ayar_yaz(t2, basarili_ayar())
    rc2, out2 = kos(t2, ek=['--skip', 'IGD'])
    d2 = durum_oku(t2)
    sina(u'S9 --skip IGD: I, G, D hic kaydedilmedi',
         not ({'I', 'G', 'D'} & set(d2.keys())), u'anahtarlar: %s' % sorted(d2.keys()))
    sina(u'S9 the exit code after --skip is 0', rc2 == 0, u'rc=%d' % rc2)
    return t


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--base', dest='taban', default=os.path.join(tempfile.gettempdir(),
                                                   'tek_tus_sinama'))
    p.add_argument('--only', dest='yalniz', default='')
    A = p.parse_args()
    ana = A.taban
    if not os.path.exists(TEK_TUS):
        print(u'ERROR: one_key.py was not found: %s' % TEK_TUS)
        return 1
    os.makedirs(ana, exist_ok=True)
    print(u'=' * 78)
    print(u'  one_key.py SCENARIO TESTS')
    print(u'  shadow root: %s   (NOT ONE BYTE is written to the mounted directory)' % ana)
    print(u'=' * 78)

    t0 = time.time()
    hepsi = A.yalniz.split(',') if A.yalniz else None

    def calis(ad):
        return (not hepsi) or ad in hepsi

    taban1 = None
    if calis('s1'):
        taban1 = s1_sifirdan(ana)
    if calis('s2') and taban1:
        s2_devam(ana, taban1)
    if calis('s3') and taban1:
        s3_bayat(ana, taban1)
    if calis('s8') and taban1:
        s8_belirlenimci_imza(ana, taban1)
    if calis('s4'):
        s4_eksik_dosya(ana)
    if calis('s5'):
        s5_asama_dustu(ana)
    if calis('s6'):
        s6_bos_cikti(ana)
    if calis('s7'):
        s7_kesinti(ana)
    if calis('s9'):
        s9_yalniz_atla(ana)

    gecti = [x for x in SONUC if x[1]]
    dustu = [x for x in SONUC if not x[1]]
    print(u'\n' + u'=' * 78)
    # MEASURED: the sentence was reordered into English but its arguments were
    # not, so it printed "total / passed" where "passed / total" belonged.
    print(u'  RESULT: %d of %d tests PASSED, %d FAILED  (%.0f s)'
          % (len(gecti), len(SONUC), len(dustu), time.time() - t0))
    if dustu:
        print(u'  FAILURES:')
        for ad, _, ayr in dustu:
            print(u'    * %s   %s' % (ad, ayr))
    print(u'=' * 78)
    return 0 if not dustu else 1


if __name__ == '__main__':
    sys.exit(main())
