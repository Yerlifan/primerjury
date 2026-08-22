# -*- coding: utf-8 -*-
"""ONE KEY - it runs the whole chain in order, in the right dependency order, from
one command.

CALLED BY : verification/one_key.py  and  verification/full_chain.py -> key B
OUTPUT    : ONE_KEY_RESULT/00_SABAH_OZETI.md   (the ONE file to look at in the morning)
            ONE_KEY_RESULT/durum.json          (the checkpoints)
            ONE_KEY_RESULT/gunluk_<time>.log   (a timestamped copy of the screen)

WHY THIS FILE EXISTS
--------------------
The existing key A (verification/full_chain.py) runs ten stages, but:
  * it does no pre-check, and discovers a missing file, tool or package in the
    middle of a run,
  * it does not know stages U and P at all (P's input is U's output),
  * it asks about checkpoint validity only as "is there a stamp"; when the INPUT
    is newer than the checkpoint it reads the stale stamp and skips the stage.
    That happened exactly on 2026-08-07 (D-9, the poisoned checkpoint).
This file closes those three gaps. full_chain.py WAS NOT TOUCHED; it still works.

THE DESIGN RULES (every one learned from a bug, none of them decoration)
-----------------------------------------------------------------------
1. THE PRE-CHECK IS A GATE. If something is missing, WHAT is missing is written
   out and the run STOPS. There is no half run. (It can be passed deliberately
   with --skip-precheck, and that is printed to the screen.)
2. THE EXIT CODE IS NOT MASKED. If rc != 0 the stage FAILED. In the past
   full_chain.py ignored stage T returning 3 and wrote "BITTI"; the summary came
   out misleading. Here rc AND the output audit are TWO SEPARATE filters and both
   must pass.
3. THE CHECKPOINT KEY IS DETERMINISTIC. md5 is used; Python's hash() function IS
   NOT (it changes between runs because of PYTHONHASHSEED, so every run misses
   the checkpoint).
4. IF THE INPUT IS NEWER THAN THE CHECKPOINT, THE CHECKPOINT IS INVALID. The
   stage's own SCRIPT counts as an input too: if the script changed, the stage
   runs again.
5. DEPENDENCIES ARE DIRECTED. The dependants of a failed stage ARE NOT RUN and
   are written "atlandi (bagimli)". Independent stages carry on.
6. THE TIME ESTIMATES ARE MEASUREMENTS. Beside every number stands the file or
   run it came from. A stage that was never measured gets NO number; it gets
   "olculmedi".

"""

import os, sys, io, csv, json, time, glob, signal, hashlib, argparse
import subprocess, threading

SURUM = u'1.0 (2026-08-08)'
CIKTI_KLASOR = 'ONE_KEY_RESULT'
CANLILIK_SN = 60          # uzun asamalarda seconds between liveness messages
LOG_YAZMA_ARALIGI = 2.0   # bagli klasore en cok bu sikligta yaziyoruz (D-11 kurali)


# ===========================================================================
#  0) KUCUK YARDIMCILAR
# ===========================================================================
def sn_metni(sn):
    sn = int(round(sn or 0))
    if sn < 90:
        return u'%d sn' % sn
    if sn < 5400:
        return u'%d dk %d sn' % (sn // 60, sn % 60)
    return u'%d sa %d dk' % (sn // 3600, (sn % 3600) // 60)


def boyut_metni(b):
    b = float(b)
    for birim in (u'B', u'KB', u'MB', u'GB'):
        if b < 1024 or birim == u'GB':
            return u'%d %s' % (b, birim) if birim == u'B' else u'%.1f %s' % (b, birim)
        b /= 1024.0


def veri_satiri_say(yol):
    """Counts the data rows in a TSV, EXCLUDING THE HEADER. Comment and blank lines do
        not count.

    """
    if not os.path.exists(yol):
        return -1
    n = 0
    with io.open(yol, encoding='utf-8', errors='ignore') as fh:
        for s in fh:
            if s.startswith('#') or not s.strip():
                continue
            n += 1
    return max(n - 1, 0)


def dosya_parmak(kok, yol):
    """A DETERMINISTIC fingerprint of a file: relative path | size | mtime_ns.

        Why not a content md5: a single file inside REFERANS_DB is 1.5 GB. Reading
        30 GB on every run is unacceptable. The size plus mtime_ns pair changes when the
        file changes and is CONSTANT between runs for THE SAME file, and those are the
        two properties we need. For directories: the file count plus the newest mtime.

    """
    try:
        st = os.stat(yol)
    except OSError:
        return u'%s|YOK' % os.path.relpath(yol, kok).replace('\\', '/')
    g = os.path.relpath(yol, kok).replace('\\', '/')
    if os.path.isdir(yol):
        n, enyeni = 0, 0
        for kk, _, dd in os.walk(yol):
            for d in dd:
                n += 1
                try:
                    enyeni = max(enyeni, os.stat(os.path.join(kk, d)).st_mtime_ns)
                except OSError:
                    pass
        return u'%s|DIR|%d|%d' % (g, n, enyeni)
    return u'%s|%d|%d' % (g, st.st_size, st.st_mtime_ns)


def imzala(parcalar):
    "md5, which is deterministic. Python's hash() function IS NOT USED (rule 3)."
    return hashlib.md5(u'\n'.join(parcalar).encode('utf-8')).hexdigest()[:16]


def en_yeni_mtime(yollar):
    en = 0.0
    for y in yollar:
        try:
            if os.path.isdir(y):
                for kk, _, dd in os.walk(y):
                    for d in dd:
                        en = max(en, os.stat(os.path.join(kk, d)).st_mtime)
            else:
                en = max(en, os.stat(y).st_mtime)
        except OSError:
            pass
    return en


def zaman_metni(t):
    return time.strftime('%d.%m %H:%M', time.localtime(t)) if t else u'yok'


def cozumle(kok, yol):
    return yol if os.path.isabs(yol) else os.path.join(kok, yol)


def yollari_ac(kok, listeler):
    """Turns the input and output lists into real paths.

        An entry written 'GLOB:pattern' expands to EVERY file matching the pattern. That
        is needed because the membership output's name carries a date:
        uyelik_yeniden_turetme_uyelik_20260803.tsv. Had a fixed path been written, a file
        produced on a new date would be invisible and the stage would run from scratch on
        every run.

    """
    out = []
    for y in listeler:
        if y.startswith(u'GLOB:'):
            out.extend(sorted(glob.glob(os.path.join(kok, y[5:]))))
        else:
            out.append(cozumle(kok, y))
    return out


# =========================================================================
#  1) THE STAGE OUTPUT AUDITORS
#     AN EXIT CODE OF 0 IS NOT ENOUGH. For every stage, "was the expected output
#     really produced, and is there data in it" is asked separately.
# =========================================================================
def d_tsv_dolu(yollar, en_az=1):
    def f(kok, ayar, cikti_metni):
        eksik, bos, tamam = [], [], []
        for y in yollar:
            t = cozumle(kok, y)
            if not os.path.exists(t):
                eksik.append(y)
            elif y.endswith('.tsv') or y.endswith('.csv'):
                n = veri_satiri_say(t)
                (tamam if n >= en_az else bos).append(u'%s (%d rows)' % (y, n))
            elif os.path.getsize(t) < 40:
                bos.append(u'%s (%d bayt)' % (y, os.path.getsize(t)))
            else:
                tamam.append(u'%s (%s)' % (y, boyut_metni(os.path.getsize(t))))
        if eksik:
            return False, 'the expected output WAS NOT PRODUCED: %s' % u', '.join(eksik)
        if bos:
            return False, 'the output is EMPTY, with no data row: %s' % u', '.join(bos)
        return True, 'the output was confirmed: %s' % u', '.join(tamam)
    return f


def d_selftest(kok, ayar, cikti_metni):
    if u'TUM SINAMALAR GECTI' in (cikti_metni or u''):
        return True, 'every self test passed'
    return False, ('The line saying every self test passed did not appear. The code could not confirm itself, so no measurement is started.')


def d_hizli_test(kok, ayar, cikti_metni):
    y = os.path.join(kok, 'QUICK_TEST', 'QUICK_TEST_REPORT.md')
    if not os.path.exists(y):
        return False, u'QUICK_TEST/QUICK_TEST_REPORT.md uretilmedi'
    m = io.open(y, encoding='utf-8', errors='ignore').read()
    if u'ZINCIR TUTARSIZ' in m:
        return False, ('the report says THE CHAIN IS INCONSISTENT. The long run IS NOT STARTED; look at QUICK_TEST/QUICK_TEST_REPORT.md.')
    if u'ZINCIR TUTARLI' in m:
        return True, 'the report says the chain is consistent'
    return False, 'the report holds neither a consistent nor an inconsistent verdict; the format is not what was expected'


def d_uyelik(kok, ayar, cikti_metni):
    g = sorted(glob.glob(os.path.join(kok, 'uyelik_yeniden_turetme_uyelik_*.tsv')))
    if not g:
        return False, u'uyelik_yeniden_turetme_uyelik_*.tsv uretilmedi'
    n = veri_satiri_say(g[-1])
    if n < 1:
        return False, u'%s BOS' % os.path.basename(g[-1])
    return True, '%s (%d rows)' % (os.path.basename(g[-1]), n)


def d_yok(kok, ayar, cikti_metni):
    return True, 'this stage produces no file, it only writes to the screen'


# =========================================================================
#  2) THE STAGE SCHEDULE  -  DEPENDENCY ORDERED, AND THE ORDER IS REQUIRED
#
#  The fields:
#    kod          the letter shown in the menu and the summary
#    ad           a one line description
#    betik        the file whose presence the pre-check looks for (None for an
#                 external tool)
#    argv(kok,a)  the list of command lists to run
#    girdi[]      the files this stage READS. If one of them is NEWER than the
#                 stage's output, THE CHECKPOINT IS INVALID.
#    cikti[]      the files this stage PRODUCES (for the output audit and freshness)
#    bagimli[]    the stage codes that must finish first
#    sure_sn      the MEASURED time (seconds). None writes "not measured".
#    kaynak       which file or run that number came from
#    denet        the output auditor
#    kraken       True: needs tools/kraken_tool.sh
#    hep_kos      True: fast and side effect free; no checkpoint is kept
#
#  WHY THE ORDER IS WHAT IT IS (read out of the code, not guessed):
#    U -> P : single_protocol_measure.py takes its membership from
#             uyelik_yeniden_turetme_uyelik_*.tsv (verification/full_chain.py line
#             1092 looks for it too).
#    P -> K : recovery_round.py's input is
#             ONE_PROTOCOL_RESULT/panel_tek_protokol.tsv
#             (verification/full_chain.py line 946 verifies it as a precondition).
#    P -> D : dogrulama_turu.siparistekiler() reads
#             ONE_PROTOCOL_RESULT/SIPARIS_LISTESI.tsv (specificity_round.py lines 150-155).
#    K -> D : D also tests the pairs K recovered (kurtarma_satirlari.tsv).
#    I -> G : all_bin_identities.py SHARES its cache with I; if I runs first the
#             same bin is not scanned twice (verification/full_chain.py line 686).
#    W -> X : the threshold scan goes through the environment check first.
#    X -> Z : the table reads the output of the threshold scan, it makes no new
#             measurement.
#    H first: H is A REGRESSION GATE; it tests against the previous reference run.
#             So it runs BEFORE P and does not wait for P's new output.
# =========================================================================
def _py(*a):
    return [sys.executable] + list(a)


KRAKEN_KOMUTLARI = {'W': ('status', 'find-db', 'db-identity', 'selftest'),
                    'X': ('threshold',),
                    'Z': ('table',)}


def ASAMALAR(ayar):
    ncbi = ayar.get('ncbi', 'oto')
    org = ayar.get('organizma', '')
    karac = ayar.get('karac')

    def d_argv(kok, a):
        arg = [os.path.join('verification', 'specificity_round.py'), '--root', '.']
        arg += ['--all'] if a.get('tumu', True) else ['--order']
        if ncbi == 'yok':
            arg += ['--local-only']
        else:
            arg += ['--ncbi', ncbi]
            if org:
                arg += ['--organism', org]
        return [_py(*arg)]

    def kraken_argv(kod):
        return lambda kok, a: [['bash', a['karac'], k] for k in KRAKEN_KOMUTLARI[kod]]

    L = [
        dict(kod='8', ad='SELF TEST: the code confirms itself and measures nothing',
             grup=u'Grup 4', betik='screening/__main__.py',
             argv=lambda kok, a: [_py('-m', 'screening', '--selftest')],
             girdi=['screening'], cikti=[], bagimli=[],
             sure_sn=4.6, kaynak='FULL_CHAIN_RESULT/durum.json, the run of 2026-08-06',
             denet=d_selftest, hep_kos=True),

        # ADDED 2026-08-10. In this project most of the bugs were not in the measurement but
        # in THE TABLE THE MEASUREMENT RESTED ON, and none of them failed a run; each was
        # found only when somebody asked "is there another bug". An audit that depends on
        # being asked is not an audit. This stage asks that question itself ON EVERY RUN.
        # It measures nothing and changes no file, which is why hep_kos=True.
        dict(kod='N', ad='THE AUDIT: the tables, the references and the seals are looked at on every run',
             grup=u'Grup 4', betik='verification/audit_all.py',
             argv=lambda kok, a: [_py(os.path.join('verification', 'audit_all.py'),
                                      '--root', '.')],
             girdi=['verification/audit_all.py',
                    'screening/target_taxids.tsv',
                    'ONE_PROTOCOL_RESULT/SIPARIS_LISTESI.tsv'],
             cikti=['ONE_KEY_RESULT/DENETIM_RAPORU.md'], bagimli=[],
             sure_sn=90.0, kaynak='measured on 2026-08-10: 2 s locally plus about 85 s for the NCBI coverage',
             # MEASURED: when this stage was added the 'denet' key was not written.
             # main() calls a['denet'] for every stage, so the chain died with
             # KeyError the moment N finished. The pre-check had been stopping the
             # run before that, which is why it was never seen. What is verified is
             # that the audit report WAS PRODUCED; findings inside it are not a
             # problem, an advisory stage exists for exactly that.
             denet=d_tsv_dolu(['ONE_KEY_RESULT/DENETIM_RAPORU.md']),
             hep_kos=True, danisma=True),

        dict(kod='H', ad='THE QUICK CONSISTENCY TEST: a regression gate BEFORE the long run',
             grup=u'Grup 4', betik='verification/quick_consistency_test.py',
             argv=lambda kok, a: [_py(os.path.join('verification', 'quick_consistency_test.py'),
                                      '--root', '.')],
             girdi=['verification/quick_consistency_test.py',
                    'ONE_PROTOCOL_RESULT/panel_tek_protokol.tsv'],
             cikti=['QUICK_TEST/QUICK_TEST_REPORT.md'], bagimli=[],
             sure_sn=122.0, kaynak='FULL_CHAIN_RESULT/durum.json, the run of 2026-08-06',
             denet=d_hizli_test),

        dict(kod='E', ad='CONFIRMING DATABASE ACCESS: is every database really readable',
             grup=u'Grup 2', betik='verification/access_check.py',
             argv=lambda kok, a: [_py(os.path.join('verification', 'access_check.py'),
                                      '--root', '.')],
             girdi=['verification/access_check.py'],
             cikti=['ACCESS_RESULT/erisim_dogrulama.tsv'], bagimli=[],
             sure_sn=None,
             kaynak='NOT MEASURED: this stage has never run, since there is no ACCESS_RESULT directory',
             denet=d_tsv_dolu(['ACCESS_RESULT/erisim_dogrulama.tsv'])),

        dict(kod='U', ad='REDERIVE THE MEMBERSHIP FROM THE MEASURED IDENTITY',
             grup=u'Grup 4', betik='engine/rederive_membership.py',
             argv=lambda kok, a: [_py(os.path.join('engine',
                                                   'rederive_membership.py'),
                                      '--root', '.')],
             girdi=['engine/rederive_membership.py', 'consensus sequences'],
             cikti=['GLOB:uyelik_yeniden_turetme_uyelik_*.tsv'], bagimli=[],
             sure_sn=None,
             kaynak='NOT MEASURED: the 1 to 3 hours in the menu is an estimate and not a measurement',
             denet=d_uyelik),

        dict(kod='P', ad=u'TEK PROTOKOL - panelin tamami TEK kuralla olculur',
             grup=u'Grup 1', betik='protocol/single_protocol_measure.py',
             argv=lambda kok, a: [_py(os.path.join('protocol', 'single_protocol_measure.py'),
                                      '--root', '.')],
             # 2026-08-09: pairs.tsv eklendi, gerekce D asamasindaki notta.
             girdi=['protocol/single_protocol_measure.py', 'primer_final',
                    'screening/pairs.tsv',
                    'GLOB:uyelik_yeniden_turetme_uyelik_*.tsv'],
             cikti=['ONE_PROTOCOL_RESULT/panel_tek_protokol.tsv',
                    'ONE_PROTOCOL_RESULT/SIPARIS_LISTESI.tsv'],
             bagimli=['U'],
             sure_sn=36.0,
             kaynak='from the comment in verification/full_chain.py, where P was measured at 36 s; ALL_STAGES_RESULT/durum.json gives 9.2 s on a warm run',
             denet=d_tsv_dolu(['ONE_PROTOCOL_RESULT/panel_tek_protokol.tsv',
                               'ONE_PROTOCOL_RESULT/SIPARIS_LISTESI.tsv'])),

        dict(kod='K', ad='RECOVERY: the rows below the threshold are recovered by four routes',
             grup=u'Grup 1', betik='verification/recovery_round.py',
             argv=lambda kok, a: [_py(os.path.join('verification', 'recovery_round.py'),
                                      '--root', '.')],
             # 2026-08-09: pairs.tsv eklendi, gerekce D asamasindaki notta.
             girdi=['verification/recovery_round.py',
                    'screening/pairs.tsv',
                    'ONE_PROTOCOL_RESULT/panel_tek_protokol.tsv'],
             cikti=['RECOVERY_RESULT/kurtarma_satirlari.tsv'], bagimli=['P'],
             sure_sn=300.0, kaynak=u'verification/full_chain.py yorumu ("K 5 dk")',
             denet=d_tsv_dolu(['RECOVERY_RESULT/kurtarma_satirlari.tsv'])),

        dict(kod='D', ad='VERIFICATION: the pairs in the panel are tested with four evidence layers',
             grup=u'Grup 1', betik='verification/specificity_round.py',
             argv=d_argv,
             # THE 2026-08-09 FIX (input tracking): NONE of the stages counted
             # screening/pairs.tsv as an input. That file holds the panel's PRIMER
             # SEQUENCES. So even if the sequences of the whole panel were changed, the
             # chain would say "the input has not changed" and skip every stage. On
             # 09.08 two pairs were changed and D would have been skipped again; what it
             # verified would have been the old pair. The root SIPARIS_LISTESI.tsv was
             # added too, because D reads the order rows from there.
             girdi=['verification/specificity_round.py', 'verification/mfeprimer_layer.py',
                    'screening/target_clades.tsv',
                    'screening/pairs.tsv',
                    'SIPARIS_LISTESI.tsv',
                    'ONE_PROTOCOL_RESULT/SIPARIS_LISTESI.tsv',
                    'RECOVERY_RESULT/kurtarma_satirlari.tsv',
                    'REFERANS_DB/SILVA_138.2_SSURef_NR99.fasta.primerqc.bin'],
             cikti=['VERIFICATION_RESULT/dogrulama_uc_sutun.tsv'], bagimli=['P', 'K'],
             sure_sn=6900.0,
             kaynak='from the run of 2026-08-07: layer 2 took 81 minutes cold over 22 pairs, layer 3 took 2.5 minutes and layer 4 at NCBI took 31 minutes, about 1 hour 55 in total',
             denet=d_tsv_dolu(['VERIFICATION_RESULT/dogrulama_uc_sutun.tsv'])),

        dict(kod='I', ad='IDENTITY VERIFICATION: the claims that go into the report are tested independently',
             grup=u'Grup 2', betik='verification/identity_verification.py',
             argv=lambda kok, a: [_py(os.path.join('verification', 'identity_verification.py'),
                                      '--root', '.')],
             girdi=['verification/identity_verification.py'],
             cikti=['IDENTITY_RESULT/kimlik_iddialari.tsv'], bagimli=[],
             sure_sn=12007.0,
             kaynak=u'FULL_CHAIN_RESULT/durum.json, 2026-08-06 (3 sa 20 dk)',
             denet=d_tsv_dolu(['IDENTITY_RESULT/kimlik_iddialari.tsv'])),

        dict(kod='G', ad='EVERY BIN IDENTITY: EVERY bin that enters the panel is confirmed',
             grup=u'Grup 2', betik='verification/all_bin_identities.py',
             argv=lambda kok, a: [_py(os.path.join('verification', 'all_bin_identities.py'),
                                      '--root', '.', '--nt', 'yok')],
             girdi=['verification/all_bin_identities.py'],
             cikti=['ALL_IDENTITIES_RESULT/tum_kutu_kimlikleri.tsv'], bagimli=['I'],
             sure_sn=17038.0,
             kaynak=u'FULL_CHAIN_RESULT/durum.json, 2026-08-06 (4 sa 44 dk)',
             denet=d_tsv_dolu(['ALL_IDENTITIES_RESULT/tum_kutu_kimlikleri.tsv'])),

        dict(kod='W', ad='THE KRAKEN2 ENVIRONMENT CHECK: is it installed and against which database',
             grup=u'Grup 3', betik=None, argv=kraken_argv('W'),
             girdi=[], cikti=[], bagimli=[],
             sure_sn=81.0, kaynak=u'FULL_CHAIN_RESULT/durum.json, 2026-08-06',
             denet=d_yok, kraken=True, hep_kos=True),

        dict(kod='X', ad='THE KRAKEN CONFIDENCE THRESHOLD SCAN',
             grup=u'Grup 3', betik=None, argv=kraken_argv('X'),
             girdi=[], cikti=[], bagimli=['W'],
             sure_sn=6916.0,
             kaynak=u'FULL_CHAIN_RESULT/durum.json, 2026-08-06 (1 sa 55 dk)',
             denet=d_yok, kraken=True),

        dict(kod='Z', ad='THE FOUR COLUMN COMPARISON TABLE',
             grup=u'Grup 3', betik=None, argv=kraken_argv('Z'),
             girdi=[], cikti=[], bagimli=['X'],
             sure_sn=4.8, kaynak=u'FULL_CHAIN_RESULT/durum.json, 2026-08-06',
             denet=d_yok, kraken=True),

        dict(kod='S', ad='REFRESH THE SUMMARY: it measures nothing and reads the existing files',
             grup=u'Grup 4', betik='screening/__main__.py',
             argv=lambda kok, a: [_py('-m', 'screening', '--mode', 'ozet')],
             girdi=[], cikti=['SCREENING_RESULT/00_OZET_HEPSI.md'], bagimli=[],
             sure_sn=0.7, kaynak=u'FULL_CHAIN_RESULT/durum.json, 2026-08-06',
             denet=d_tsv_dolu(['SCREENING_RESULT/00_OZET_HEPSI.md']),
             hep_kos=True),
    ]
    for a in L:
        a.setdefault('kraken', False)
        a.setdefault('hep_kos', False)
        a.setdefault('danisma', False)
        a['karac'] = karac
    return L


# ===========================================================================
#  3) ON KONTROL  -  EKSIK VARSA DURUR
# ===========================================================================
def on_kontrol(kok, ayar, yaz):
    """Every condition is measured one by one before the run starts.

        If a single REQUIRED item fails, the run DOES NOT START. That is deliberate: a
        half run produces a misleading summary to be read in the morning.

    """
    zorunlu_dusen, uyari, sayac = [], [], [0]
    yaz(u'')
    yaz(u'=' * 78)
    yaz(u'  PRE-CHECK - every precondition is measured before running')
    yaz(u'=' * 78)

    def satir(ad, ok, ayrinti, zorunlu=True):
        sayac[0] += 1
        yaz(u'  [%s] %-46s %s' % (u' OK  ' if ok else u'MISSING', ad, ayrinti))
        if not ok:
            (zorunlu_dusen if zorunlu else uyari).append(u'%s -> %s' % (ad, ayrinti))

    # --- 1) Python surumu ve paketler -------------------------------------
    pv = '%d.%d.%d' % sys.version_info[:3]
    satir(u'Python 3.8+', sys.version_info[:2] >= (3, 8), pv)
    for paket, kurulum in (('numpy', 'numpy'), ('primer3', 'primer3-py')):
        try:
            m = __import__(paket)
            satir(u'python paketi: %s' % paket, True, getattr(m, '__version__', 'var'))
        except Exception:
            satir(u'python paketi: %s' % paket, False,
                  u'MISSING, to install: pip3 install %s --break-system-packages' % kurulum)

    # --- 2) Betikler -------------------------------------------------------
    for a in ASAMALAR(ayar):
        if not a.get('betik'):
            continue
        t = os.path.join(kok, a['betik'])
        satir(u'script: %s' % a['betik'], os.path.exists(t),
              boyut_metni(os.path.getsize(t)) if os.path.exists(t) else u'NOT FOUND')

    # --- 3) The data directories ------------------------------------------
    # MEASURED: this list counted the 'engine' directory FOUR TIMES. In the
    # rename four separate source directories fell onto one name and the
    # copies were never cleaned up; the pre-check measured the same directory
    # four times and printed four identical lines. The list now names the
    # directories this repository actually has, once each.
    for d, zor in (('fastq files', True), ('consensus sequences', True),
                   ('primer_final', True), ('REFERANS_DB', True),
                   ('screening', True), ('protocol', True),
                   ('engine', True), ('steps', True),
                   ('tools', True), ('verification', True),
                   ('konsensus_kanonik', False)):
        t = os.path.join(kok, d)
        var = os.path.isdir(t)
        satir(u'directory: %s' % d, var,
              (u'%d oge' % len(os.listdir(t))) if var else u'NOT FOUND', zor)

    # --- 4) MFEprimer ikilisi ve indeksleri --------------------------------
    mfe = None
    for aday in (os.path.join(kok, 'tools', 'mfeprimer'), 'mfeprimer'):
        if os.path.exists(aday):
            mfe = aday + (u'' if os.access(aday, os.X_OK) else '  (THERE IS NO EXECUTE PERMISSION)')
            if not os.access(aday, os.X_OK):
                mfe = None
                uyari.append(u'tools/mfeprimer is there but is not executable: chmod +x tools/mfeprimer')
            break
        try:
            from shutil import which
            w = which(aday)
            if w:
                mfe = w
                break
        except Exception:
            pass
    satir(u'MFEprimer ikilisi', bool(mfe),
          mfe or u'tools/mfeprimer IS MISSING or is not executable')

    MFE_IX = ['archaea.16S.fna', 'bacteria.16S.fna', 'fungi.ITS.fna',
              'fungi.28SrRNA.fna', 'fungi.18SrRNA.fna',
              'SILVA_138.2_SSURef_NR99.fasta']
    for f in MFE_IX:
        p = os.path.join(kok, 'REFERANS_DB', f + '.primerqc.bin')
        satir(u'MFE index: %s' % f, os.path.exists(p),
              boyut_metni(os.path.getsize(p)) if os.path.exists(p)
              else 'MISSING, to build it: mfeprimer index -i REFERENCE_DB/%s' % f)

    # --- 5) Katman 2'nin taradigi 11 kume ----------------------------------
    KUMELER = ['SILVA_138.2_SSURef_NR99.fasta', 'SILVA_138.2_LSURef_NR99.fasta',
               'UNITE_ITS.fasta', 'PR2_SSU_taxo_long.fasta',
               'ROD_v1.2_operon_variants.fasta', 'bacteria.16S.fna',
               'archaea.16S.fna', 'fungi.ITS.fna', 'fungi.28SrRNA.fna',
               'fungi.18SrRNA.fna', 'ref_all2.fna']
    eksik_kume = [f for f in KUMELER
                  if not os.path.exists(os.path.join(kok, 'REFERANS_DB', f))]
    satir(u'the layer 2 databases (11 sets)', not eksik_kume,
          u'11/11 yerinde' if not eksik_kume else u'MISSING: %s' % u', '.join(eksik_kume))

    # --- 6) The SILVA SSU RNA/DNA gate -------------------------------------
    # In the past SILVA's RNA alphabet (U) broke the index and every binding came back
    # 0/0. We confirm it is DNA by counting U against T in the first few thousand lines.
    sp = os.path.join(kok, 'REFERANS_DB', 'SILVA_138.2_SSURef_NR99.fasta')
    if os.path.exists(sp):
        try:
            with io.open(sp, encoding='utf-8', errors='ignore') as fh:
                ornek = u''.join(s for s in (fh.readline() for _ in range(4000))
                                 if not s.startswith(u'>'))
            nu, nt = ornek.count(u'U'), ornek.count(u'T')
            satir(u'SILVA SSU alfabesi DNA mi (U=0 olmali)', nu == 0 and nt > 0,
                  u'U=%d  T=%d' % (nu, nt))
        except Exception as e:
            satir(u'SILVA SSU alfabesi', False, u'could not be read: %s' % e, zorunlu=False)

    # --- 7) Disk alani -----------------------------------------------------
    try:
        st = os.statvfs(kok)
        bos_gb = st.f_bavail * st.f_frsize / 1073741824.0
        satir(u'free disk on the mounted directory (>= 5 GB)', bos_gb >= 5, u'%.1f GB bos' % bos_gb)
    except Exception as e:
        satir(u'free disk on the mounted directory', False, u'olculemedi: %s' % e, zorunlu=False)
    try:
        st = os.statvfs('/tmp')
        tmp_gb = st.f_bavail * st.f_frsize / 1073741824.0
        # MFEprimer kaniti YEREL diskte uretilip toplu kopyalaniyor (D-11).
        # Olculen en buyuk kanit 282 MB; 2 GB rahat pay birakir.
        satir(u'free space on the local /tmp (>= 2 GB, for D-11)', tmp_gb >= 2,
              u'%.1f GB bos' % tmp_gb)
    except Exception as e:
        satir(u'yerel /tmp bos alani', False, u'olculemedi: %s' % e, zorunlu=False)

    # --- 8) WSL bellegi ----------------------------------------------------
    try:
        mem = {}
        for s in io.open('/proc/meminfo', encoding='utf-8'):
            p = s.split(':')
            if len(p) == 2:
                mem[p[0]] = int(p[1].strip().split()[0])
        top_gb = mem.get('MemTotal', 0) / 1048576.0
        # Olculen tepe bellek: mfeprimer indeksleme 6,06 GB
        # (REFERANS_DB/SILVA_138.2_SSURef_NR99.fasta.BOZUK_KANIT.txt).
        # Tarama tarafi bundan dusuk; 4 GB alt sinir, 8 GB rahat.
        satir(u'total WSL memory (>= 4 GB)', top_gb >= 3.8, u'%.1f GB' % top_gb)
        if 3.8 <= top_gb < 7.5:
            uyari.append(u'WSL memory %.1f GB. The Kraken steps (X) can choke on a large database. KILITLENME_COZUMU.md, verification/full_chain.py -> key M.' % top_gb)
    except Exception as e:
        satir(u'WSL bellegi', False, u'olculemedi: %s' % e, zorunlu=False)

    # --- 9) Kraken araci (istege bagli) ------------------------------------
    satir(u'Kraken araci (tools/kraken_tool.sh)', bool(ayar.get('karac')),
          ayar.get('karac') or ayar.get('kraken_sebep', u'missing, W/X/Z WILL BE SKIPPED'),
          zorunlu=False)

    # --- 10) The network (for the NCBI layer) ------------------------------
    if ayar.get('ncbi') == 'oto':
        ok, ayr = _ag_dene()
        satir(u'NCBI access (--ncbi oto was chosen)', ok, ayr, zorunlu=False)
        if not ok:
            uyari.append(u'NCBI could not be reached. The 4th layer of stage D stays "BILINMIYOR"; the local layers still run.')

    # --- 11) Write permission -----------------------------------------------
    # The only thing asked is CAN WE WRITE. Permission to delete is A SEPARATE question,
    # and this chain deletes nothing; a test file that cannot be removed does not stop
    # the run.
    try:
        t = os.path.join(kok, CIKTI_KLASOR)
        os.makedirs(t, exist_ok=True)
        p = os.path.join(t, '_yazma_denemesi.txt')
        with io.open(p, 'w', encoding='utf-8') as fh:
            fh.write(u'write attempt %s\n' % time.strftime('%Y-%m-%d %H:%M:%S'))
        try:
            os.remove(p)
            ek = u'yazilabiliyor'
        except Exception:
            ek = 'it can be written to; the trial file could not be deleted because there is no delete permission, and the chain deletes nothing anyway'
        satir(u'write permission on the mounted directory', True, u'%s/ %s' % (CIKTI_KLASOR, ek))
    except Exception as e:
        satir(u'write permission on the mounted directory', False, u'CANNOT BE WRITTEN: %s' % e)

    yaz(u'')
    if uyari:
        yaz(u'  WARNINGS (the run continues, but know what you are passing over):')
        for u_ in uyari:
            yaz(u'    ! %s' % u_)
        yaz(u'')
    if zorunlu_dusen:
        yaz(u'=' * 78)
        yaz(u'  PRE-CHECK FAILED - THE RUN DID NOT START')
        yaz(u'=' * 78)
        yaz(u'  %d items are missing:' % len(zorunlu_dusen))
        for z in zorunlu_dusen:
            yaz(u'    * %s' % z)
        yaz(u'')
        yaz(u'  No stage was run. A half-finished run would produce a misleading')
        yaz(u'  summary, so it was stopped DELIBERATELY.')
        yaz(u'  Fix what is missing and run the same command again.')
        yaz(u'=' * 78)
        return False, zorunlu_dusen, uyari
    yaz(u'  PRE-CHECK PASSED - %d items checked, nothing mandatory is missing.' % sayac[0])
    return True, zorunlu_dusen, uyari


def _ag_dene():
    try:
        import urllib.request
        t0 = time.time()
        r = urllib.request.urlopen(
            'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi', timeout=12)
        r.read(200)
        r.close()
        return True, u'ulasildi (%.1f sn)' % (time.time() - t0)
    except Exception as e:
        return False, u'ulasilamadi: %s' % str(e)[:70]


# ===========================================================================
#  4) KONTROL NOKTASI - GECERLILIK
# ===========================================================================
def girdi_imzasi(kok, a):
    """A deterministic fingerprint of everything the stage READS.

        The stage's OWN SCRIPT is an input too: if the script changed, the old result is
        stale. This closes the same class as the 'poisoned checkpoint' bug of 2026-08-07,
        where the index had been rebuilt but the checkpoint read the old zeros back.

    """
    p = [u'surum=%s' % SURUM, u'kod=%s' % a['kod']]
    if a.get('betik'):
        p.append(dosya_parmak(kok, os.path.join(kok, a['betik'])))
    acilan = yollari_ac(kok, a['girdi'])
    if not acilan and a['girdi']:
        p.append(u'GIRDI_BULUNAMADI')
    for t in acilan:
        p.append(dosya_parmak(kok, t) if os.path.exists(t)
                 else u'%s|YOK' % os.path.relpath(t, kok).replace('\\', '/'))
    return imzala(p)


def _girdi_yollari(kok, a):
    y = yollari_ac(kok, a['girdi'])
    if a.get('betik'):
        y.append(os.path.join(kok, a['betik']))
    return y


def kontrol_noktasi_gecerli(kok, a, durum):
    """(can_it_be_skipped, reason) - WHY it was or was not skipped is written openly.

        Four filters, and all of them must pass:
          1. there is a stamp and it says 'bitti'
          2. the input signature is the same (the input has not changed)
          3. all the expected outputs exist and are non-empty
          4. the newest output is NEWER than the newest input (a stale output is not read)

    """
    if a.get('hep_kos'):
        return False, 'this stage runs again on every run, because it is fast and has no side effects'

    d = durum.get(a['kod'], {})
    imza = girdi_imzasi(kok, a)
    ciktilar = yollari_ac(kok, a['cikti'])
    girdiler = _girdi_yollari(kok, a)
    cy, gy = en_yeni_mtime(ciktilar), en_yeni_mtime(girdiler)

    if d.get('durum') == 'bitti':
        if d.get('imza') != imza:
            return False, ('the checkpoint is STALE: the input signature changed (the stamp says %s and it is now %s). It will run again.'
                           % (str(d.get('imza'))[:10], imza[:10]))
        if ciktilar:
            ok, mesaj = a['denet'](kok, {}, u'')
            if not ok:
                return False, 'there is a stamp but the output check failed: %s' % mesaj
            if gy > cy:
                return False, ('the checkpoint is STALE: the input is NEWER than the output (%s against %s). It will run again.'
                               % (zaman_metni(gy), zaman_metni(cy)))
        return True, ('it finished on the previous run (%s) and the input has not changed'
                      % sn_metni(d.get('sure', 0)))

    # There is no stamp. Is there a READY and FRESH output on disk?
    # This is for recognising results produced before this script was written; re-running
    # them out of caution would waste hours. But "it exists" IS NOT ENOUGH: "newer than
    # the input" is required too, or we accept a stale result as a fresh one.
    if ciktilar:
        ok, mesaj = a['denet'](kok, {}, u'')
        if ok:
            if cy >= gy and cy > 0:
                return True, ('there is no stamp but the output on disk is FRESH (%s against the newest input %s). %s' % (zaman_metni(cy), zaman_metni(gy), mesaj))
            return False, ('there is an output on disk but it is STALE: the input is %s and the output %s. It will run again.' % (zaman_metni(gy), zaman_metni(cy)))
        return False, 'it has not been run before (%s)' % mesaj
    return False, 'it has not been run before'


# ===========================================================================
#  5) TEK ASAMAYI KOS
# ===========================================================================
class Kesildi(Exception):
    pass


KESME = {'var': False}


def asama_kos(kok, a, ayar, yaz):
    """Runs all of a stage's commands in order.

        * The exit code is read and IS NOT MASKED.
        * On the first non-zero code the stage's remaining commands are not run.
        * A sign of life is printed during long silences.

    """
    t0 = time.time()
    rc, son_cikti = 0, u''
    komutlar = a['argv'](kok, ayar)
    cevre = dict(os.environ)
    cevre.setdefault('IPLIK', '3')
    cevre.setdefault('ZORLA_MMAP', '1')
    cevre['PYTHONUNBUFFERED'] = '1'
    if ayar.get('vt_a'):
        cevre['VT_A'] = ayar['vt_a']
    # The Kraken tool has to be called from its own directory (it looks for its
    # companion files relatively).
    calisma = os.path.dirname(a['karac']) if (a.get('kraken') and a.get('karac')) else kok

    for i, argv in enumerate(komutlar, 1):
        if KESME['var']:
            raise Kesildi()
        yaz(u'   $ [%d/%d] %s' % (i, len(komutlar),
                                  u' '.join(os.path.basename(x) if os.sep in x else x
                                            for x in argv)))
        try:
            p = subprocess.Popen(argv, cwd=calisma, env=cevre,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        except Exception as e:
            rc, son_cikti = 127, 'the command could not be started: %s' % e
            yaz(u'   ERROR: %s' % son_cikti)
            break

        satirlar = []
        son = {'t': time.time()}
        dur = threading.Event()

        def canlilik():
            # Uzun asamalarda "asildi mi?" sorusunu ekranda cevaplar.
            while not dur.wait(CANLILIK_SN):
                yaz(u'   ... [alive] %s | in stage %s | since the last output %s'
                    % (a['kod'], sn_metni(time.time() - t0),
                       sn_metni(time.time() - son['t'])))

        th = threading.Thread(target=canlilik)
        th.daemon = True
        th.start()
        try:
            for ham in p.stdout:
                s = ham.decode('utf-8', 'replace').rstrip('\r\n')
                satirlar.append(s)
                son['t'] = time.time()
                yaz(u'      | ' + s)
                if KESME['var']:
                    break
            p.wait()
        finally:
            dur.set()
        rc = p.returncode if p.returncode is not None else -1
        son_cikti = u'\n'.join(satirlar)
        if KESME['var']:
            raise Kesildi()
        if rc != 0:
            yaz(u'   >>> command %d/%d RETURNED A NON-ZERO CODE: %s'
                % (i, len(komutlar), rc))
            break
    return rc, son_cikti, time.time() - t0


# =========================================================================
#  6) THE PLAN, THE SUMMARY, THE MAIN FLOW
# =========================================================================
def kraken_bul(kok):
    # MEASURED: only <root>/../tools/WSL/ was being looked in, that is, a directory
    # that is a SIBLING of the project. In this repository tools/ is INSIDE the root,
    # so the tool was never found and W, X and Z were silently skipped on every run.
    adaylar = [os.path.join(kok, 'tools', 'kraken_tool.sh'),
               os.path.abspath(os.path.join(kok, '..', 'tools', 'WSL', 'kraken_tool.sh'))]
    for aday in adaylar:
        if os.path.exists(aday):
            return aday, u''
    return None, (u'the Kraken tool was not found (looked in: %s). W, X and Z are '
                  u'SKIPPED; the rest of the chain runs.' % u', '.join(adaylar))


def plan_yaz(kok, asamalar, durum, ayar, yaz):
    yaz(u'')
    yaz(u'=' * 78)
    yaz(u'  PLAN - which stages will run and which will be skipped')
    yaz(u'=' * 78)
    toplam, olculmeyen = 0.0, []
    for a in asamalar:
        atla, sebep = kontrol_noktasi_gecerli(kok, a, durum)
        if a.get('kraken') and not ayar.get('karac'):
            atla, sebep = True, ayar.get('kraken_sebep', u'Kraken araci yok')
            a['_kraken_atla'] = True
        a['_atla'], a['_sebep'] = atla, sebep
        if atla:
            sure = u'0 sn'
        elif a['sure_sn'] is None:
            sure = u'OLCULMEDI'
            olculmeyen.append(a['kod'])
        else:
            sure = sn_metni(a['sure_sn'])
            toplam += a['sure_sn']
        yaz(u'  %-2s %-8s %-12s %s' % (a['kod'], u'SKIPPED' if atla else u'KOSACAK',
                                       sure, a['ad'][:44]))
        yaz(u'       %s' % sebep)
    yaz(u'')
    yaz(u'  Sum of the MEASURED stages    : %s' % sn_metni(toplam))
    if olculmeyen:
        yaz(u'  Stages with NO MEASURED duration : %s' % u', '.join(olculmeyen))
        yaz(u'  No number was INVENTED for those stages. The total above is therefore')
        yaz(u'  a LOWER BOUND: not an estimate, but an incomplete measurement.')
    yaz(u'=' * 78)
    return toplam, olculmeyen


# The source of the final verdict table. They are tried in order and the FIRST ONE
# FOUND is used. The NEWEST file stands first. If a new kind of verdict is produced,
# it has to be added at the head of this list; the code does not choose the file
# itself, it follows the order written here.
# THE 2026-08-10 FIX. An ESIK_VE_OLCUT dated 08-08 stood at the head of the list and
# the morning summary said "11 pairs ORDERABLE", while the
# ONE_PROTOCOL_RESULT/SIPARIS_LISTESI.tsv produced by THE SAME RUN said 20 (15 KESIN
# plus 5 EVRENSEL). The difference is that five pairs were changed after 08-08 and
# Microascaceae was brought back. The summary said openly which file it had read,
# so it was not lying, but it was READING THE WRONG FILE, and that is the first
# place anyone looks in the morning. The table THE RUN ITSELF produces was moved to
# the head of the list.
SIPARIS_KAYNAKLARI = [
    ('ONE_PROTOCOL_RESULT/SIPARIS_LISTESI.tsv', 'durum'),
    ('ESIK_VE_OLCUT_2026-08-08.tsv', 'YENI_HUKUM'),
    ('NIHAI_SIPARIS_LISTESI_2026-08-07.tsv', None),
]

# The verdict prefixes that count as ORDERABLE. The rule is written OPENLY rather
# than being an intuition buried in the code:
#   * "SIPARIS EDILEBILIR..." -> group specific, the seven pairs past the threshold
#   * "KOSULLU..."            -> universal and control primers. Their dCq is
#                                undefined (the competitor set approaches empty), so
#                                no threshold verdict can be given; the measure is
#                                COVERAGE.
# "ESIK ALTI...", "ONERILMEZ..." and "UYELIK DOGRULANAMADI" DO NOT enter the order.
#   * "ESIK USTU - SIPARIS EDILEBILIR" -> SIPARIS_LISTESI.tsv's own wording. The
#                                same rule, a different word order.
SIPARIS_ONEKLERI = (u'SIPARIS EDILEBILIR', u'KOSULLU', u'ESIK USTU')


def siparis_tablosu(kok):
    """Reads the final order table FROM THE FILE ON DISK; it DOES NOT RECOMPUTE it."""
    y = hk = None
    for ad, sut_adi in SIPARIS_KAYNAKLARI:
        t = os.path.join(kok, *ad.split('/'))
        if os.path.exists(t):
            y, hk = t, sut_adi
            break
    if not y:
        return (u'The final verdict table was not found. Files looked for: %s\n'
                % u', '.join(u'`%s`' % a for a, _ in SIPARIS_KAYNAKLARI))
    with io.open(y, encoding='utf-8', errors='ignore') as fh:
        r = list(csv.DictReader((s for s in fh if not s.startswith('#')), delimiter='\t'))
    ad = os.path.basename(y)
    if not r:
        return u'`%s` is empty.\n' % ad
    sut = [k for k in r[0].keys() if k]
    if not hk or hk not in sut:
        hk = (next((k for k in sut if 'HUKUM' in k.upper()), None)
              or next((k for k in sut if 'KARAR' in k.upper()), None)
              or next((k for k in sut if 'SINIF' in k.upper()), None))
    hd = next((k for k in sut if 'hedef' in k.lower()), None)
    if not hk or not hd:
        return (u'`%s` was read (%d rows) but the verdict or target column was not '
                u'recognised. Columns: %s\n' % (ad, len(r), u', '.join(sut)))

    def hukum(s):
        return (s.get(hk) or u'?').strip()

    # If there is a CLASS column the verdict is read FROM IT. The reason: on universal
    # primers dCq IS UNDEFINED (the denominator of the competitor set goes to zero) and
    # the table writes "OLCULEMEDI - KARAR YOK" for them. A rule looking only at the
    # text counts those five control primers as not orderable, when they are the panel's
    # controls and they DO go to order. dogrulama_turu.siparistekiler() uses the same
    # rule; if two places do not use the same definition, two different numbers come out.
    sinif_s = next((k for k in sut if k.strip().upper() == 'SINIF'), None)

    def siparise_girer(s):
        if sinif_s:
            return (s.get(sinif_s) or '').strip().upper() in ('KESIN', 'EVRENSEL')
        return hukum(s).upper().startswith(SIPARIS_ONEKLERI)

    girer = [s for s in r if siparise_girer(s)]
    girmez = [s for s in r if not siparise_girer(s)]
    say = {}
    for s in r:
        say[hukum(s)] = say.get(hukum(s), 0) + 1

    out = [u'Source: `%s` (%d rows), column `%s`. '
           u'**This summary DID NOT RECOMPUTE the table**, it read it from the '
           u'file.\n\n' % (ad, len(r), hk),
           u'**ORDERABLE: %d pairs = %d oligos** · the remaining %d pairs are '
           u'not ordered (they were not deleted, they are listed below).\n\n'
           % (len(girer), 2 * len(girer),
                                                        len(girmez)),
           (u'The rule: a row whose `SINIF` column is `KESIN` or `EVRENSEL` goes '
            u'into the order. On universal and control primers dCq IS UNDEFINED '
            u'(the denominator of the competitor set goes to zero) and the table '
            u'writes "OLCULEMEDI - KARAR YOK" for them; the measure is coverage, '
            u'not a threshold.\n\n') if sinif_s else
           (u'The rule: a row whose verdict starts with `SIPARIS EDILEBILIR...` '
            u'or `KOSULLU...` goes into the order.\n\n'),
           u'| %s | how many pairs |\n|---|---|\n' % hk]
    for k in sorted(say, key=lambda x: (-say[x], x)):
        out.append(u'| %s | %d |\n' % (k, say[k]))

    dcq = next((k for k in sut if 'dCq' in k), None)
    out.append(u'\n### The %d pairs to order\n\n| target | verdict | dCq |\n|---|---|---|\n'
               % len(girer))
    for s in girer:
        out.append(u'| %s | %s | %s |\n' % (s.get(hd) or u'?', hukum(s),
                                            (s.get(dcq) or u'-') if dcq else u'-'))
    out.append(u'\n### The %d pairs not ordered (they were not deleted)\n\n'
               u'| target | verdict | dCq |\n|---|---|---|\n' % len(girmez))
    for s in girmez:
        out.append(u'| %s | %s | %s |\n' % (s.get(hd) or u'?', hukum(s),
                                            (s.get(dcq) or u'-') if dcq else u'-'))
    return u''.join(out)


def ozet_yaz(kok, asamalar, durum, ayar, kesildi, on_uyari, baslangic, gunluk_yolu):
    yol = os.path.join(kok, CIKTI_KLASOR, '00_SABAH_OZETI.md')
    basari = [a for a in asamalar if durum.get(a['kod'], {}).get('durum') == 'bitti']
    dusen = [a for a in asamalar if durum.get(a['kod'], {}).get('durum') == 'DUSTU']
    atlanan = [a for a in asamalar
               if (durum.get(a['kod'], {}).get('durum') or '').startswith('atlandi')]
    with io.open(yol, 'w', encoding='utf-8') as fh:
        w = fh.write
        w(u'# Morning summary, one-key run\n\n')
        w(u'Generated: %s, version %s\n\n' % (time.strftime('%Y-%m-%d %H:%M'), SURUM))
        w(u'Run started: %s, elapsed: %s\n\n'
          % (time.strftime('%Y-%m-%d %H:%M', time.localtime(baslangic)),
             sn_metni(time.time() - baslangic)))
        if kesildi:
            w(u'> **RUN INTERRUPTED** (Ctrl+C or the window was closed). Finished stages are saved; running the same command again resumes where it stopped.\n\n')

        w(u'## 1. What happened at each stage\n\n')
        w(u'| stage | status | time | exit code | note |\n|---|---|---|---|---|\n')
        for a in asamalar:
            d = durum.get(a['kod'], {})
            w(u'| **%s** %s | %s | %s | %s | %s |\n'
              % (a['kod'], a['ad'][:40],
                 (u'done (with warnings)' if d.get('uyarili') else d.get('durum', u'not run')),
                 sn_metni(d.get('sure', 0)) if d.get('sure') else '-',
                 d.get('cikis', '-'), (d.get('sebep') or '').replace('|', '/')[:170]))
        w(u'\n**Succeeded %d | Failed %d | Skipped %d**\n\n'
          % (len(basari), len(dusen), len(atlanan)))

        uyarili = [a for a in asamalar if durum.get(a['kod'], {}).get('uyarili')]
        if uyarili:
            # An advisory stage does not stop the chain but it DOES NOT STAY SILENT
            # either. If there is a finding it appears at the top of the summary;
            # otherwise an audit has been made that nobody looked at, and the gate means
            # nothing.
            w(u'## WARNING: the audit found something, but the chain still ran\n\n')
            for a in uyarili:
                d = durum[a['kod']]
                w(u'- **%s** %s — %s\n' % (a['kod'], a['ad'], (d.get('sebep') or '')[:200]))
            w(u'\nDetail: `ONE_KEY_RESULT/DENETIM_RAPORU.md` and `GECE_BULGULARI.md`.\n\n')

        if dusen:
            w(u'## 2. FAILED STAGES: look at these first\n\n')
            for a in dusen:
                d = durum[a['kod']]
                w(u'### %s — %s\n\n' % (a['kod'], a['ad']))
                w(u'- exit code: `%s`\n- reason: %s\n' % (d.get('cikis'), d.get('sebep')))
                if d.get('son_satirlar'):
                    w(u'- the stage\'s last output:\n\n```\n%s\n```' % d['son_satirlar'])
                w(u'\n')
        else:
            w(u'## 2. NO STAGE FAILED\n\nNo stage returned a non-zero code and no output check failed.\n\n')

        w(u'## 3. Final order table\n\n')
        w(siparis_tablosu(kok))

        w(u'\n## 4. Files you should look at\n\n')
        w(u'| question | file |\n|---|---|\n')
        for soru, dosya in (
            (u'What should I order (FINAL)',
             u'`ONE_PROTOCOL_RESULT/SIPARIS_LISTESI.tsv`, the `durum` and '
             u'`siparis_sarti` columns (the table THIS RUN produced itself)'),
            (u'What was asked for at the meeting, and what came of it',
             u'`TOPLANTI_DURUMU.md`'),
            (u'What the NCBI 4th layer said',
             u'`VERIFICATION_RESULT/NCBI_KATMAN4_RAPORU.md`'),
            (u'Are the audits clean in this run',
             u'`ONE_KEY_RESULT/DENETIM_RAPORU.md`'),
            (u'Why the threshold rule changed', u'`ESIK_VE_OLCUT_2026-08-08.md`'),
            (u'Where to copy the sequences from',
             '`PrimerJury_PANEL_*.xlsx`, the order sheet'),
            (u'Which pair is risky, and why', u'`SIPARIS_KARARI_2026-08-07.md`'),
            (u'Did a contradiction come up in this run',
             u'`VERIFICATION_RESULT/CELISKILER.md`'),
            (u'The verification layers side by side',
             u'`VERIFICATION_RESULT/DOGRULAMA_RAPORU.md`'),
            (u'The detail of the off target products',
             u'`HEDEF_DISI_AYRINTI_2026-08-07.tsv`'),
            (u'Are the bin identities right',
             u'`ALL_IDENTITIES_RESULT/TUM_KUTU_KIMLIK_RAPORU.md`'),
            (u'The raw output of this run',
             u'`%s`' % os.path.relpath(gunluk_yolu, kok).replace('\\', '/')),
        ):
            w(u'| %s | %s |\n' % (soru, dosya))

        if on_uyari:
            w(u'\n## 5. Pre-check warnings\n\n')
            for u_ in on_uyari:
                w(u'- %s\n' % u_)

        w(u'\n---\n\n### What this run did not measure\n\n')
        olcumsuz = [a for a in asamalar if a['sure_sn'] is None]
        for a in olcumsuz:
            w(u'- **%s** duration was not measured: %s\n' % (a['kod'], a['kaynak']))
        if not olcumsuz:
            w(u'- Every stage that ran had been timed before.\n')
        if not ayar.get('karac'):
            w(u'- Kraken comparison (W, X, Z) WAS NOT RUN: %s\n'
              % ayar.get('kraken_sebep', ''))
    return yol



# --- CLI value normalisation ------------------------------------------------
# English option values are accepted alongside the original Turkish ones and
# mapped back here. The internal values are unchanged on purpose: they are
# compared in dozens of places and, in some cases, written to output files.
# Translating the interface must not translate the data.
_DEGER = {'auto': 'oto', 'manual': 'elle', 'none': 'yok', 'quick': 'hizli',
          'full': 'tam', 'member': 'uye', 'competitor': 'rakip',
          'exclude': 'disla'}


def _ing_deger(a):
    for _ad in ('nt', 'literatur', 'ncbi', 'karisik', 'moduller', 'mod'):
        _v = getattr(a, _ad, None)
        if isinstance(_v, str) and _v in _DEGER:
            setattr(a, _ad, _DEGER[_v])
    return a

def main():
    global CANLILIK_SN
    p = argparse.ArgumentParser(description='ONE KEY: runs the whole chain in '
                                            'order')
    p.add_argument('--root', dest='kok', default='.')
    p.add_argument('--plan', action='store_true', help=u'only plani yaz, kosma')
    p.add_argument('--confirm', dest='onayla', action='store_true',
                   help='run WITHOUT showing the plan and waiting for a '
                        'confirmation')
    p.add_argument('--rerun', dest='yeniden', action='store_true',
                   help=u're-run finished stages as well')
    p.add_argument('--only', dest='yalniz', default='', help=u'these stages only, e.g. 8HS')
    p.add_argument('--skip', dest='atla', default='', help=u'skip these stages, e.g. IG')
    p.add_argument('--ncbi', choices=['auto', 'manual', 'none', 'oto', 'elle', 'yok'], default='oto')
    p.add_argument('--organism', dest='organizma',
                   default='Bacteria (taxid:2) OR Archaea (taxid:2157) OR Fungi (taxid:4751)')
    p.add_argument('--order-16', dest='siparis_16', action='store_true',
                   help='the ordering stage tests the 16 ordered pairs only, '
                        'instead of all 22')
    p.add_argument('--db-path', dest='vt', default=os.environ.get('VT_A', ''),
                   help='the path of the Kraken2 database')
    p.add_argument('--skip-precheck', dest='on_kontrol_atla', action='store_true',
                   help=u'skip the pre-check - not recommended; it is reported on screen')
    p.add_argument('--liveness', dest='canlilik', type=int, default=CANLILIK_SN,
                   help=u'seconds between liveness messages')
    A = p.parse_args()
    CANLILIK_SN = max(2, A.canlilik)

    kok = os.path.abspath(A.kok)
    os.makedirs(os.path.join(kok, CIKTI_KLASOR), exist_ok=True)
    zaman = time.strftime('%Y%m%d_%H%M%S')
    gunluk_yolu = os.path.join(kok, CIKTI_KLASOR, 'gunluk_%s.log' % zaman)
    gunluk = io.open(gunluk_yolu, 'w', encoding='utf-8')
    tampon = {'t': 0.0}

    def yaz(s=u''):
        # Without a timestamp on the screen, WITH one in the log. Writing to a mounted
        # directory is limited by LOG_YAZMA_ARALIGI (D-11: many small writes to /mnt/c
        # were measured to be slow).
        try:
            print(s, flush=True)
        except UnicodeEncodeError:
            print(s.encode('ascii', 'replace').decode('ascii'), flush=True)
        gunluk.write(u'%s  %s\n' % (time.strftime('%H:%M:%S'), s))
        if time.time() - tampon['t'] > LOG_YAZMA_ARALIGI:
            gunluk.flush()
            tampon['t'] = time.time()

    def sig(signum, frame):
        KESME['var'] = True
        yaz(u'\n  !! INTERRUPT RECEIVED. Writing state and exiting cleanly...')
    signal.signal(signal.SIGINT, sig)
    try:
        signal.signal(signal.SIGTERM, sig)
    except Exception:
        pass

    baslangic = time.time()
    yaz(u'=' * 78)
    yaz(u'  ONE KEY - the full PrimerJury chain   version %s' % SURUM)
    yaz(u'  %s' % time.strftime('%Y-%m-%d %H:%M:%S'))
    yaz(u'  kok    : %s' % kok)
    yaz(u'  gunluk : %s' % os.path.relpath(gunluk_yolu, kok).replace('\\', '/'))
    yaz(u'=' * 78)

    karac, kraken_sebep = kraken_bul(kok)
    ayar = dict(ncbi=A.ncbi, organizma=A.organizma, tumu=not A.siparis_16,
                karac=karac, kraken_sebep=kraken_sebep, vt_a=A.vt)

    on_uyari = []
    if A.on_kontrol_atla:
        yaz(u'\n  !! PRE-CHECK SKIPPED (--skip-precheck). If something is missing it will surface MID-RUN.')
    else:
        ok, _dusen, on_uyari = on_kontrol(kok, ayar, yaz)
        if not ok:
            gunluk.flush()
            gunluk.close()
            return 2

    hepsi = ASAMALAR(ayar)
    sec = [a for a in hepsi
           if (not A.yalniz or a['kod'] in A.yalniz.upper())
           and a['kod'] not in A.atla.upper()]
    if not sec:
        yaz(u'  No stage selected (--only / --skip). Nothing was run.')
        gunluk.close()
        return 0

    dyol = os.path.join(kok, CIKTI_KLASOR, 'durum.json')
    durum = {}
    if os.path.exists(dyol) and not A.yeniden:
        try:
            durum = json.load(io.open(dyol, encoding='utf-8'))
        except Exception as e:
            yaz(u'  WARNING: durum.json could not be read (%s); checkpoints were ignored.' % e)
            durum = {}
    if A.yeniden:
        yaz(u'  --rerun was given: every checkpoint is being ignored.')
        durum = {}

    plan_yaz(kok, sec, durum, ayar, yaz)
    if A.plan:
        gunluk.flush()
        gunluk.close()
        return 0
    if not A.onayla:
        try:
            c = input(u'\n  Baslatilsin mi? (E = evet, baska tus = vazgec): ').strip()
        except EOFError:
            c = u''
        if c.upper() not in (u'E', u'EVET'):
            yaz(u'  Abandoned, nothing was run.')
            gunluk.close()
            return 0

    def kaydet():
        json.dump(durum, io.open(dyol, 'w', encoding='utf-8'),
                  ensure_ascii=False, indent=1)

    kesildi = False
    kalan = [a for a in sec if not a['_atla']]
    yaz(u'\n  Stages to run: %d' % len(kalan))
    sirano = 0
    for a in sec:
        kod = a['kod']
        # THE PLAN IS AN ESTIMATE, THIS IS THE VERDICT. Earlier stages may have
        # changed the disk (for example, running U changes P's input), so the skip
        # decision is recomputed AT THE MOMENT the stage is reached.
        if not a.get('_kraken_atla'):
            yeni_atla, yeni_sebep = kontrol_noktasi_gecerli(kok, a, durum)
            if yeni_atla != a['_atla']:
                yaz(u'\n   NOTE: the decision for %s differed from the PLAN (%s -> %s). Reason: %s' % (kod, u'SKIPPED' if a['_atla'] else u'KOSACAK',
                                    u'SKIPPED' if yeni_atla else u'KOSACAK', yeni_sebep))
            a['_atla'], a['_sebep'] = yeni_atla, yeni_sebep

        if a['_atla']:
            if a.get('_kraken_atla'):
                durum[kod] = dict(durum=u'atlandi (arac yok)', sebep=a['_sebep'], sure=0)
            elif durum.get(kod, {}).get('durum') != 'bitti':
                durum[kod] = dict(durum='bitti', sure=0, cikis=0,
                                  imza=girdi_imzasi(kok, a), sebep=a['_sebep'],
                                  zaman=time.strftime('%Y-%m-%d %H:%M'))
            yaz(u'\n>> %s  %s\n   SKIPPED - %s' % (kod, a['ad'], a['_sebep']))
            kaydet()
            continue

        # Bagimlilik kapisi: bagimli oldugu bir asama dustuyse KOSMA.
        engel = [b for b in a['bagimli']
                 if (durum.get(b, {}).get('durum') or '').startswith(
                     ('DUSTU', 'atlandi (bagimli', 'kesildi'))]
        if engel:
            durum[kod] = dict(durum='atlandi (bagimli)', sure=0,
                              sebep="it did not run because stage %s did not finish. This stage takes that stage's output as its INPUT, and running it on an empty input would have produced a convincing but meaningless result." % u', '.join(engel))
            yaz(u'\n>> %s  %s\n   SKIPPED (it depends on another stage) - %s'
                % (kod, a['ad'], durum[kod]['sebep']))
            kaydet()
            continue

        sirano += 1
        kalan_sn = sum(x['sure_sn'] or 0 for x in kalan[sirano - 1:])
        yaz(u'\n' + u'-' * 78)
        yaz(u'>> [%d/%d] %s  %s' % (sirano, len(kalan), kod, a['ad']))
        yaz(u'   started : %s' % time.strftime('%H:%M:%S'))
        yaz(u'   measured time: %s   (source: %s)'
            % (sn_metni(a['sure_sn']) if a['sure_sn'] is not None else u'NOT MEASURED',
               a['kaynak']))
        yaz(u'   stages left: %d | measured time left: %s'
            % (len(kalan) - sirano, sn_metni(kalan_sn)))
        yaz(u'-' * 78)

        try:
            rc, cikti_metni, sure = asama_kos(kok, a, ayar, yaz)
        except Kesildi:
            kesildi = True
            durum[kod] = dict(durum='kesildi', sure=0,
                              sebep='the user interrupted it. This stage is HALF DONE; its own checkpoints are on record and pressing the same key continues from where it stopped.')
            kaydet()
            yaz(u'   KESILDI - durum kaydedildi.')
            break

        son_satirlar = u'\n'.join((cikti_metni or u'').splitlines()[-15:])
        tamam, mesaj = a['denet'](kok, ayar, cikti_metni)

        # TWO SEPARATE FILTERS, AND BOTH MUST PASS.
        # In the past only the output audit was looked at and a non-zero code was
        # ignored (full_chain.py, stage T, 2026-08-06): T returned exit code 3 and
        # still got a "BITTI" stamp, and the summary came out misleading.
        if rc != 0:
            mesaj = ('EXIT CODE %s, which is not zero. The output check: %s' % (rc, mesaj))
            tamam = False
        # AN ADVISORY STAGE: it reports findings but DOES NOT FAIL the chain.
        # N (THE AUDIT) is such a stage: its job is "to say what is broken right
        # now". Finding something broken does not mean the chain should not run; on
        # the contrary, let the chain run so that what was produced can be seen too.
        # But it does not stay silent either: UYARILI is written into the summary and
        # the findings sit in DENETIM_RAPORU.md. (2026-08-10: when N was added, the
        # tek_tus_sinama S1 scenario started failing; the gate's job is not to fail a
        # test.)
        if a.get('danisma') and rc != 0:
            mesaj = 'WITH A WARNING: %s (an advisory stage, so the chain was not stopped)' % mesaj
            tamam = True
            uyarili = True
        else:
            uyarili = False
        if not tamam:
            durum[kod] = dict(durum='DUSTU', sure=sure, cikis=rc, sebep=mesaj,
                              son_satirlar=son_satirlar,
                              zaman=time.strftime('%Y-%m-%d %H:%M'))
            yaz(u'   << %s DUSTU (%s) - %s' % (kod, sn_metni(sure), mesaj))
            bagimlilar = [x['kod'] for x in sec if kod in x['bagimli']]
            if bagimlilar:
                yaz(u'   The stages that depend on it WILL NOT RUN: %s' % u', '.join(bagimlilar))
            else:
                yaz(u'   Nothing else depends on it; the chain continues.')
        else:
            durum[kod] = dict(durum='bitti', sure=sure, cikis=rc, sebep=mesaj,
                              uyarili=uyarili, imza=girdi_imzasi(kok, a),
                              zaman=time.strftime('%Y-%m-%d %H:%M'))
            yaz(u'   << %s FINISHED (%s) - %s' % (kod, sn_metni(sure), mesaj))
        kaydet()

    ozet = ozet_yaz(kok, sec, durum, ayar, kesildi, on_uyari, baslangic, gunluk_yolu)
    dusen = [a['kod'] for a in sec if durum.get(a['kod'], {}).get('durum') == 'DUSTU']
    atlanan_b = [a['kod'] for a in sec
                 if durum.get(a['kod'], {}).get('durum') == 'atlandi (bagimli)']
    yaz(u'\n' + u'=' * 78)
    if kesildi:
        yaz(u'  RUN INTERRUPTED. Finished stages are saved; run the same command again.')
    elif dusen:
        yaz(u'  RUN FINISHED but %d STAGES FAILED: %s' % (len(dusen), u', '.join(dusen)))
        if atlanan_b:
            yaz(u'   Not run because they depend on it: %s' % u', '.join(atlanan_b))
    else:
        yaz(u'  EVERY STAGE SUCCEEDED.')
    yaz(u'  Total elapsed time: %s' % sn_metni(time.time() - baslangic))
    yaz(u'')
    yaz(u'  LOOK AT THIS FIRST : %s' % os.path.relpath(ozet, kok).replace('\\', '/'))
    yaz(u'  Raw output         : %s' % os.path.relpath(gunluk_yolu, kok).replace('\\', '/'))
    yaz(u'=' * 78)
    gunluk.flush()
    gunluk.close()
    if kesildi:
        return 130
    return 3 if dusen else 0


if __name__ == '__main__':
    sys.exit(main() or 0)
