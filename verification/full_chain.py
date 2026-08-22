# -*- coding: utf-8 -*-
"""
Full-chain driver, runs the whole pipeline from one command.

THIS SCRIPT TAKES NO MEASUREMENTS. Its job is to call ten stages in dependency
order, INSPECT each stage's output, and STOP if anything does not hold. The
measurements live in the scripts it calls.

THREE DESIGN DECISIONS
    1  Resumability is at STAGE level. Each finished stage is stamped in
       durum.json and skipped on re-run. Within-stage resumability belongs to
       the called scripts, which already checkpoint; two layers of it would
       overlap and confuse.

    2  OUTPUT INSPECTION IS MANDATORY; AN EXIT CODE IS NOT ENOUGH. The failure
       that keeps recurring in this domain is a program that produces a WRONG
       or EMPTY answer without erroring. A zero exit code does not mean the
       work was done. So after every stage: does the expected file exist, is it
       non-empty, and, for some stages, what does it actually say? A failed
       inspection STOPS the chain; it is never skipped over quietly.

    3  KRAKEN STAGES DO NOT BREAK THE CHAIN. Those stages need a separate tool
       and a separate database. Refusing an eight-hour run because they are
       absent would be wrong, so they are marked SKIPPED with a reason and the
       chain continues. The summary records which stages were skipped and why.

"""
# -------------------------------------------------------------------------
# full_chain.py - the engine behind key A in the verification/full_chain.py menu.
#
# INPUT  : the stage scripts in the project root (the screening package, the four
#          scripts under verification) and tools/kraken_tool.sh;
#          plus FULL_CHAIN_RESULT/durum.json (where the previous run stopped).
# OUTPUT : FULL_CHAIN_RESULT/00_TAM_ZINCIR_OZET.md (one combined summary),
#          FULL_CHAIN_RESULT/durum.json, FULL_CHAIN_RESULT/kosu_gunlugu.txt.
# CALLED BY: verification/full_chain.py -> key A
#          (python3 verification/full_chain.py --root . --confirm)
#
# THREE DESIGN DECISIONS, with their reasons
#
# 1) INTERRUPTION TOLERANCE IS AT THE STAGE LEVEL.
#    As each stage finishes, a "bitti" stamp and its duration are written into
#    durum.json. When the same option is pressed again, finished stages are
#    SKIPPED. Interruption tolerance WITHIN a stage is not this script's job; the
#    scripts it calls have their own checkpoints and already do that. The two
#    layers do not overlap here.
#
# 2) THE OUTPUT AUDIT IS REQUIRED; AN EXIT CODE IS NOT ENOUGH.
#    The kind of bug that keeps recurring in this project is this: a program
#    produces a WRONG or EMPTY answer without raising an error. An exit code of
#    zero does not mean "the job is done". So after every stage it is asked
#    whether the expected file EXISTS, whether it is EMPTY, and, on some stages,
#    what its CONTENT says. A stage that fails the audit STOPS the chain; there is
#    no silent move to the next one.
#
# 3) THE KRAKEN STEPS DO NOT BREAK THE CHAIN.
#    Steps W, X, Y and Z need a separate tool (kraken2) and a separate database.
#    Refusing an eight hour run outright because those are not installed would be
#    wrong; those steps are marked ATLANDI with their reason and the chain carries
#    on. Which steps were skipped, and why, also stands in the summary.
# -------------------------------------------------------------------------

from __future__ import print_function
import argparse
import csv
import io
import json
import os
import shutil
import subprocess
import sys
import time

VERSIYON = u'1.0 (2026-08-04)'
CIKTI_ADI = 'FULL_CHAIN_RESULT'


# -------------------------------------------------------------------------
# THE STAGE TABLE
#
# The fields:
#   kod      the single letter or digit in the menu - the key in the summary and in
#            durum.json
#   ad       the name shown on screen and in the summary
#   grup     the group in the menu (for grouping in the report)
#   dk       the estimated time in MINUTES (low, high). An estimate only.
#   kraken   when True, the step is SKIPPED if kraken2 is absent; the chain is not
#            broken
#   komut    the command producer to run (kok, ayar) -> [argv, ...]
#   denet    the output auditor (kok, ayar, ciktilar) -> (ok, message)
#
# THE ORDER IS NOT ARBITRARY:
#   8 first, because no measurement is entered before the code has tested itself.
#   H next, because before an eight hour chain is entered, the chain must be shown
#     to be whole on a small subset with a known answer.
#   W before the kraken steps, because starting a threshold scan without knowing
#     which database is in use spends hours on the wrong database.
#   I and G before T: T's membership decisions rest on the bin identities.
#   T in the middle: it runs stages P, K and D in dependency order within itself.
#   X and Y measure and Z reads their output, which is why Z is last.
#   S at the very end, because a summary only means anything once everything has
#     finished.
# -------------------------------------------------------------------------


def _py(kok, *arg):
    return [sys.executable] + list(arg)


def _satir_sayisi(yol):
    'Counts the data rows in a TSV, NOT COUNTING the header. Comment lines are not counted.'
    if not os.path.exists(yol):
        return -1
    n = 0
    with io.open(yol, encoding='utf-8', errors='ignore') as fh:
        for i, s in enumerate(fh):
            if s.startswith('#') or not s.strip():
                continue
            n += 1
    return max(n - 1, 0)


def d_dosya_dolu(yollar, en_az_satir=1):
    """Do ALL the expected files EXIST and are they NON-EMPTY?

        Why a separate audit: in this project a stage finishing with a zero code and
        producing no rows at all is a bug that has actually been seen. 'The file
        exists' is not enough; 'there is data in it' is asked as well.

    """
    def f(kok, ayar):
        eksik, bos = [], []
        for y in yollar:
            t = os.path.join(kok, y)
            if not os.path.exists(t):
                eksik.append(y)
                continue
            if y.endswith('.tsv'):
                if _satir_sayisi(t) < en_az_satir:
                    bos.append(y)
            elif os.path.getsize(t) < 40:
                bos.append(y)
        if eksik:
            return False, 'the expected output IS MISSING: %s' % ', '.join(eksik)
        if bos:
            return False, 'the output is EMPTY or holds no data row: %s' % ', '.join(bos)
        return True, 'the output was confirmed: %s' % ', '.join(yollar)
    return f


def d_sina(kok, ayar):
    'The self test key: its output has to say that every self test passed.'
    g = ayar.get('_son_cikti', '')
    if u'TUM SINAMALAR GECTI' in g:
        return True, 'every self test passed'
    return False, ('the self test text was not found. The code could not '
                   'confirm itself, so no measurement is started.')


def d_hizli_test(kok, ayar):
    """Key H: does the report say TUTARLI. If it says TUTARSIZ the chain STOPS.

        The reason: that is exactly why this step exists. Spending eight hours on an
        inconsistent chain produces a number that will be thrown away at the end.

    """
    y = os.path.join(kok, 'QUICK_TEST', 'QUICK_TEST_REPORT.md')
    if not os.path.exists(y):
        return False, u'QUICK_TEST/QUICK_TEST_REPORT.md uretilmedi'
    m = io.open(y, encoding='utf-8', errors='ignore').read()
    if u'ZINCIR TUTARSIZ' in m:
        return False, ('the report says THE CHAIN IS INCONSISTENT. The full '
                       'run IS NOT STARTED; look at '
                       'QUICK_TEST/QUICK_TEST_REPORT.md.')
    if u'ZINCIR TUTARLI' in m:
        return True, 'the report says the chain is consistent'
    return False, 'the report holds neither a consistent nor an inconsistent verdict; the format is not what was expected'


def d_esik(kok, ayar):
    """X and Y: was the threshold curve produced, and is the AYRILIK column empty?

        A filled AYRILIK means two independent measurements gave different numbers at
        the same threshold. By the project rule the numbers are not trusted in that
        case and the chain stops.

    """
    a = ayar['kraken_is']
    csvy = os.path.join(a, 'esik_egrisi.csv')
    if not os.path.exists(csvy):
        return False, 'the threshold curve file was not produced (%s)' % csvy
    try:
        with io.open(csvy, encoding='utf-8', errors='ignore') as fh:
            r = list(csv.DictReader(fh))
    except Exception as e:
        return False, 'the threshold curve file could not be read: %s' % e
    if not r:
        return False, 'the threshold curve file is EMPTY; no threshold was measured'
    ad = [k for k in (r[0].keys() if r else []) if k and k.strip().lower() == 'ayrilik']
    if ad:
        dolu = [x for x in r if (x.get(ad[0]) or '').strip()]
        if dolu:
            return False, ('the DIVERGENCE column is FILLED IN on %d rows. '
                           'Two independent measurements diverged, so the '
                           'numbers are not to be trusted.' % len(dolu))
    return True, '%d thresholds were measured with no divergence' % len(r)


def d_tablo(kok, ayar):
    """Z: was the table produced, and is it COMPLETE?

        A THREE STATE AUDIT. This stage carries a special case: the table can be
        produced with MISSING columns and that is not an error - once the Kraken
        measurement is made later, the same key completes the table. But an incomplete
        table MUST NOT GET a 'bitti' stamp; if it did, the resume logic would never try
        it again and the table would stay incomplete forever.
        So there is a third state:
            True    -> a complete table, finished
            'eksik' -> the table was produced but a column or columns are empty; the
                       chain DOES NOT STOP and this stage is TRIED AGAIN on the next run
            False   -> the table could not be produced at all; the chain stops

    """
    y = os.path.join(ayar['arac'], '0_TESLIM_RAPOR', 'KRAKEN_KARSILASTIRMA.md')
    if not os.path.exists(y):
        return False, u'KRAKEN_KARSILASTIRMA.md uretilmedi'
    if os.path.getsize(y) < 200:
        return False, u'KRAKEN_KARSILASTIRMA.md neredeyse bos (%d bayt)' % os.path.getsize(y)
    m = io.open(y, encoding='utf-8', errors='ignore').read()
    eksik = []
    for anahtar, ad in (('there is NO PlusPFP run', u'PlusPFP sutunu'),
                        ('there is no threshold scan', u'esik sutunu'),
                        (u'kimlik_sonuc.csv bulunamadi', u'hizalama sutunu')):
        if anahtar in m:
            eksik.append(ad)
    if eksik:
        return 'eksik', ('the table was produced but it is INCOMPLETE: %s. '
                         'Once the data arrives the same stage completes it, '
                         'which is why it does not count as finished.'
                         % ', '.join(eksik))
    return True, 'the comparison table was produced COMPLETE (%d bytes)' % os.path.getsize(y)



def d_ozet(kok, ayar):
    y = os.path.join(kok, 'SCREENING_RESULT', '00_OZET_HEPSI.md')
    if not os.path.exists(y):
        return False, u'00_OZET_HEPSI.md uretilmedi'
    yas = time.time() - os.path.getmtime(y)
    if yas > 3600:
        return False, ('the summary was NOT REFRESHED in this run; it was '
                       'written %d minutes ago' % (yas / 60))
    return True, u'ozet yenilendi'


def d_kraken_ortam(kok, ayar):
    """W: is kraken2 present?

        If kraken2 IS ABSENT that IS NOT AN ERROR; the later kraken steps are skipped
        and the chain carries on. But it does not count as 'bitti' either: when the
        user installs kraken2 later, this stage HAS TO RUN AGAIN. So 'eksik' is
        returned; the resume logic skips only the ones stamped 'bitti'.

    """
    if ayar.get('kraken_var'):
        return True, 'kraken2 was found (%s) and the database check passed' % (
            ayar.get('kraken2_bin') or 'PATH')
    return 'eksik', ('kraken2 or its database is missing, so the Kraken '
                     'stages will be skipped. Once it is installed this stage '
                     'is tried again. %s'
                     % ayar.get('kraken_sebep', ''))


ASAMALAR = [
    ('8', u'SELF-TEST - the code verifies itself; no measurement', u'Group 4',
     (1, 2), False,
     lambda kok, a: [_py(kok, '-m', 'screening', '--selftest')],
     d_sina),

    ('H', u'QUICK CONSISTENCY TEST - is the chain intact on a small subset?', u'Group 4',
     (25, 40), False,
     lambda kok, a: [_py(kok, os.path.join('verification', 'quick_consistency_test.py'),
                         '--root', '.')],
     d_hizli_test),

    # W's kraken flag is deliberately False: this step is a DIAGNOSTIC step. If kraken2
    # is missing it MUST NOT BE SKIPPED, quite the opposite, it MUST RUN, because this
    # is the place that prints what is missing and how to install it. Skipped, the user
    # would never see the installation instructions.
    ('W', u'KRAKEN2 ENVIRONMENT CHECK - installed? which database?', u'Group 3',
     (1, 5), False,
     lambda kok, a: [['bash', a['karac'], 'status'],
                     ['bash', a['karac'], 'find-db'],
                     ['bash', a['karac'], 'db-identity'],
                     ['bash', a['karac'], 'selftest']],
     d_kraken_ortam),

    # MEASURED on the clean run of 2026-08-05: 3 h 20 min (KISA_LISTE 500 plus idf/BM25
    # ranking). The estimate was pulled to the measured value; the old 6-8 hours
    # belonged to the KISA_LISTE=1000 assumption and was twice as wide as reality.
    ('I', u'IDENTITY VERIFICATION - every reported claim is tested', u'Group 2',
     (180, 260), False,
     lambda kok, a: [_py(kok, os.path.join('verification', 'identity_verification.py'),
                         '--root', '.')],
     d_dosya_dolu([os.path.join('IDENTITY_RESULT', 'kimlik_iddialari.tsv')])),

    # MEASURED on the clean run of 2026-08-05: 4 h 43 min. The estimate was pulled to the measured value.
    ('G', u'ALL BIN IDENTITIES - every bin entering the panel is verified', u'Group 2',
     (260, 350), False,
     lambda kok, a: [_py(kok, os.path.join('verification', 'all_bin_identities.py'),
                         '--root', '.')],
     d_dosya_dolu([os.path.join('ALL_IDENTITIES_RESULT', 'tum_kutu_kimlikleri.tsv')])),

    # Measured: P 36 s, K 5 min, D 6 min. The previous estimate (6-16 hours) was far too wide.
    ('T', u'FULL MEASUREMENT - P, K, D and I in dependency order', u'Group 1',
     (30, 90), False,
     lambda kok, a: [_py(kok, os.path.join('verification', 'run_all_stages.py'), '--root', '.')],
     d_dosya_dolu([os.path.join('ALL_STAGES_RESULT', '00_BIRLESIK_OZET.md')])),

    # MEASURED on the clean run of 2026-08-05: 1 h 55 min (6 thresholds x ~19 min). The
    # old estimate was 10-40 min and was out BY A FACTOR OF THREE: the sampling was set
    # to 100 000 reads but the source held 86 780, so the sampling NEVER CAME INTO PLAY
    # and the scan ran on the full data. The estimate was pulled to what actually happened.
    ('X', u'KRAKEN CONFIDENCE THRESHOLD SWEEP (on a sample)', u'Group 3',
     (100, 140), True,
     lambda kok, a: [['bash', a['karac'], 'threshold']],
     d_esik),

    # Y IS KEPT SEPARATE: reclassifying against the full database is a one off and it is
    # heavy. It is left for the night. It runs with mmap and few threads too.
    ('Y', u'RE-RUN WITH PlusPFP (heavy; leave it overnight)', u'Group 3',
     (120, 480), True,
     lambda kok, a: [['bash', a['karac'], 'db-identity'],
                     ['bash', a['karac'], 'threshold'],
                     ['bash', a['karac'], 'table']],
     d_esik),

    ('Z', u'FOUR-COLUMN COMPARISON TABLE', u'Group 3',
     (1, 2), True,
     lambda kok, a: [['bash', a['karac'], 'table']],
     d_tablo),

    ('S', u'REFRESH THE COMBINED SUMMARY - no measurement', u'Group 4',
     (1, 1), False,
     lambda kok, a: [_py(kok, '-m', 'screening', '--mode', 'ozet')],
     d_ozet),
]


# ---------------------------------------------------------------------------
def sure_metni(dk):
    if dk < 60:
        return u'%d min' % dk
    s, k = divmod(int(dk), 60)
    return u'%d h %d min' % (s, k) if k else u'%d h' % s


def sn_metni(sn):
    sn = int(sn)
    if sn < 60:
        return u'%d s' % sn
    if sn < 3600:
        return u'%d min %d s' % (sn // 60, sn % 60)
    return u'%d h %d min' % (sn // 3600, (sn % 3600) // 60)


def kraken2_bul(karac, ortam):
    """Asks THE TOOL ITSELF for kraken2's full path.

        WHY IT IS DONE THIS WAY (the 2026-08-04 fix)
        This used to be shutil.which('kraken2'), that is, only PATH was looked at. But
        in this project kraken2 IS NOT on PATH: in the source study's installation it
        sits in micromamba's "mikro" environment. The result was that, although the
        binary was on disk, every Kraken step was skipped as "not on PATH".

        The fix is NOT to rewrite the search here but to call the tool's own search
        logic. The ortam_ac function inside kraken_tool.sh was inherited from the
        source study's rerun_kraken.sh and it activates the environment, looks in the
        environment directories directly, and as a last resort searches the home
        directory.

    """
    if not os.path.exists(karac):
        return None, 'the Kraken tool was not found (%s)' % karac
    cevre = dict(os.environ)
    if ortam:
        cevre['ORTAM'] = ortam
    try:
        p = subprocess.Popen(['bash', karac, 'kraken-path'],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             env=cevre)
        cikti, hata = p.communicate(timeout=180)
    except Exception as e:
        return None, u'kraken-yol cagrilamadi: %s' % e
    cikti = cikti.decode('utf-8', 'replace')
    hata = hata.decode('utf-8', 'replace')
    alan = {}
    for satir in cikti.splitlines():
        if '=' in satir and satir.split('=', 1)[0].isupper():
            k, v = satir.split('=', 1)
            alan[k.strip()] = v.strip()
    y = alan.get('KRAKEN2_BIN', '')
    if y and os.path.exists(y):
        # the method and the version travel too: the user should see on screen which route
        # it was found by and which version will run, rather than guessing
        return y, u'%s | %s | %s' % (
            y, alan.get('KRAKEN_YONTEM', 'yontem bildirilmedi'),
            alan.get('KRAKEN2_SURUM', 'surum okunamadi'))
    ilk = [s for s in hata.splitlines() if s.strip()]
    return None, (ilk[0] if ilk else u'kraken2 bulunamadi')


def kraken_ortami(kok, pluspfp, vt_a, ortam=''):
    """MEASURES whether the Kraken side is usable; it does not guess.

        Three things are looked for: the tool itself, the kraken2 binary and at least
        one database. If any of the three is missing the kraken steps are skipped and
        the reason is written; the chain DOES NOT STOP.

    """
    # MEASURED: this used to look for <root>/../tools/WSL/kraken_tool.sh, that is, a
    # directory that is a SIBLING of the project. In this repository tools/ is INSIDE
    # the root and the script sits there directly, so the tool was never found and W, X
    # and Z were silently skipped on every run. The right place is tried first and the
    # old layout is kept as a fallback.
    arac = os.path.join(kok, 'tools')
    karac = os.path.join(arac, 'kraken_tool.sh')
    if not os.path.exists(karac):
        eski = os.path.abspath(os.path.join(kok, '..', 'tools'))
        eski_karac = os.path.join(eski, 'WSL', 'kraken_tool.sh')
        if os.path.exists(eski_karac):
            arac, karac = eski, eski_karac
    sebep = []
    kbin, kmesaj = kraken2_bul(karac, ortam)
    if not kbin:
        sebep.append(kmesaj)

    # The database: the given path first, then ~/k2db, then THE TOOL'S own search.
    # The tool's vt-ara key scans the disk; we do not rewrite that here.
    vt = vt_a or os.path.join(os.path.expanduser('~'), 'k2db')
    if not os.path.exists(os.path.join(vt, 'hash.k2d')):
        if kbin:
            try:
                p = subprocess.Popen(['bash', karac, 'find-db'],
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     env=dict(os.environ, KRAKEN2_BIN=kbin))
                c = p.communicate(timeout=300)[0].decode('utf-8', 'replace')
                for s in c.splitlines():
                    s = s.strip()
                    if s.endswith('hash.k2d'):
                        aday = os.path.dirname(s)
                        if os.path.exists(os.path.join(aday, 'hash.k2d')):
                            vt = aday
                            break
            except Exception:
                pass
        if not os.path.exists(os.path.join(vt, 'hash.k2d')):
            sebep.append(u'There is no Kraken2 database (hash.k2d was looked for inside %s; the tool also swept the disk with vt-ara)' % vt)
    return {
        'arac': arac,
        'karac': karac,
        'kraken_is': os.path.join(kok, 'RESULTS', 'kraken_esik_A'),
        'kraken_var': not sebep,
        'kraken_sebep': u'; '.join(sebep),
        'kraken2_bin': kbin,
        'kraken_mesaj': kmesaj,
        'vt_a': vt,
        'pluspfp': pluspfp,
        'ortam': ortam,
    }


def plan_yaz(yaz, secili, ayar, durum):
    yaz(u'')
    yaz(u'=' * 78)
    yaz(u'  FULL CHAIN - PLAN')
    yaz(u'=' * 78)
    alt = ust = 0
    yaz(u'')
    yaz(u'  %-3s %-52s %-12s %s' % (u'#', u'STAGE', u'TIME', u'STATUS'))
    yaz(u'  ' + u'-' * 74)
    for kod, ad, grup, (a, u_), kraken, _k, _d in secili:
        d = durum.get(kod, {})
        if d.get('durum') == 'bitti':
            nd = 'FINISHED, it will be skipped'
        elif d.get('durum') == 'eksik':
            nd = 'INCOMPLETE, it will be tried again'
        elif kraken and not ayar['kraken_var']:
            nd = u'WILL SKIP (no kraken2)'
        elif kod == 'Y' and not ayar['pluspfp']:
            nd = u'WILL SKIP (no PlusPFP path)'
        else:
            nd = u'will run'
            alt += a
            ust += u_
        yaz(u'  %-3s %-52s %-12s %s' % (kod, ad[:52], sure_metni(a) + u'-' + sure_metni(u_), nd))
    yaz(u'')
    yaz(u'  ESTIMATED TOTAL TIME: between %s and %s' % (sure_metni(alt), sure_metni(ust)))
    if ayar.get('kraken2_bin'):
        parca = (ayar.get('kraken_mesaj') or '').split(' | ')
        yaz(u'')
        yaz(u'  kraken2      : %s' % ayar['kraken2_bin'])
        if len(parca) >= 3:
            yaz(u'  found at     : %s' % parca[1])
            yaz(u'  version      : %s' % parca[2])
        yaz(u'  database     : %s' % ayar.get('vt_a'))
        yaz(u'  (which release, PlusPF or PlusPFP, is determined at stage W)')
    if not ayar['kraken_var']:
        yaz(u'')
        yaz(u'  KRAKEN STAGES WILL BE SKIPPED. Reason:')
        for s in ayar['kraken_sebep'].split('; '):
            yaz(u'    * %s' % s)
        yaz(u'  The chain does NOT stop for this; the remaining stages run normally.')
        yaz(u'')
        yaz(u'  YOU CAN COMPLETE THE KRAKEN PART LATER. You do not have to re-run the')
        yaz(u'  whole chain: a skipped stage is NOT stamped "done", so the next run')
        yaz(u'  retries it by itself. For the Kraken part only:')
        yaz(u'      W, then X, then Z')
        yaz(u'  or in a single command:')
        yaz(u'      python3 verification/full_chain.py --root . --only W,X,Z,S --confirm')
        yaz(u'  Different environment name: --env <name>   |   known binary path: KRAKEN2_BIN=/full/path/kraken2')
    yaz(u'')
    yaz(u'  If interrupted, run the same command again: finished stages are skipped.')
    yaz(u'=' * 78)
    return alt, ust


def calistir(kok, ayar, secili, durum, dyol, gunluk, yaz, kuru):
    """Runs the stages in order. A stage that fails the audit STOPS the chain."""
    for kod, ad, grup, dk, kraken, komut_f, denet in secili:
        d = durum.setdefault(kod, {})
        if d.get('durum') == 'bitti':
            yaz(u'\n>> %s  %s\n   SKIPPED - already finished in a previous run (%s)'
                % (kod, ad, sn_metni(d.get('sure', 0))))
            continue

        if kraken and not ayar['kraken_var']:
            d.update(durum='atlandi', sebep=ayar['kraken_sebep'], sure=0)
            yaz(u'\n>> %s  %s\n   SKIPPED - %s' % (kod, ad, ayar['kraken_sebep']))
            json.dump(durum, io.open(dyol, 'w', encoding='utf-8'),
                      ensure_ascii=False, indent=1)
            continue

        if kod == 'Y' and not ayar['pluspfp']:
            d.update(durum='atlandi', sebep='no PlusPFP database path was given', sure=0)
            yaz(u'\n>> %s  %s\n   SKIPPED - no PlusPFP path given (pass it with --pluspfp)' % (kod, ad))
            json.dump(durum, io.open(dyol, 'w', encoding='utf-8'),
                      ensure_ascii=False, indent=1)
            continue

        yaz(u'\n>> %s  %s' % (kod, ad))
        yaz(u'   started: %s' % time.strftime('%H:%M:%S'))
        t0 = time.time()
        cevre = dict(os.environ)
        # kraken2's FULL PATH is passed explicitly to every kraken step. PATH is not
        # trusted: a subprocess opens a new shell and the micromamba environment may not
        # be active in that shell. Giving the full path removes the problem.
        if ayar.get('kraken2_bin'):
            cevre['KRAKEN2_BIN'] = ayar['kraken2_bin']
        if ayar.get('ortam'):
            cevre['ORTAM'] = ayar['ortam']
        # WSL2'de tepe bellegi dusuk tutan iki ayar. Burada acikca gecirilir ki
        # kullanicinin kabuk ortamina bagli kalmasin.
        cevre.setdefault('IPLIK', os.environ.get('IPLIK', '3'))
        cevre.setdefault('ZORLA_MMAP', '1')
        if kod in ('X', 'Y'):
            cevre.setdefault('ORNEK', os.environ.get('ORNEK', '100000'))
        if kod == 'Y':
            cevre['VT_A'] = ayar['pluspfp']
        elif kod in ('X', 'Z', 'W') and ayar.get('vt_a'):
            cevre['VT_A'] = ayar['vt_a']

        rc, son_cikti = 0, ''
        if kuru:
            yaz(u'   [DRY RUN] the command was not executed')
        else:
            for argv in komut_f(kok, ayar):
                yaz(u'   $ %s' % ' '.join(os.path.basename(x) for x in argv))
                try:
                    p = subprocess.Popen(argv, cwd=kok, env=cevre,
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.STDOUT)
                except Exception as e:
                    rc, son_cikti = 127, 'the command could not be started: %s' % e
                    yaz(u'   %s' % son_cikti)
                    break
                bicik = []
                for ham in p.stdout:
                    s = ham.decode('utf-8', 'replace').rstrip('\n')
                    bicik.append(s)
                    gunluk.write(s + u'\n')
                p.wait()
                son_cikti = u'\n'.join(bicik)
                rc = p.returncode
                if rc != 0:
                    break
        sure = time.time() - t0
        ayar['_son_cikti'] = son_cikti

        tamam, mesaj = denet(kok, ayar)
        # A third state: 'eksik'. The stage ran, produced output, but did not complete.
        # The chain DOES NOT STOP (otherwise one incomplete Kraken column would cut the
        # whole run), but it DOES NOT GET a 'bitti' stamp either; so it is tried again on
        # the next run and completes itself once the data arrives.
        if tamam == 'eksik':
            d.update(durum='eksik', sure=sure, cikis=rc, sebep=mesaj)
            json.dump(durum, io.open(dyol, 'w', encoding='utf-8'),
                      ensure_ascii=False, indent=1)
            yaz(u'   INCOMPLETE (%s) - %s' % (sn_metni(sure), mesaj))
            yaz(u'   The chain continues. This stage will be retried on the next run.')
            continue
        # The exit code and the output audit are TWO SEPARATE filters. Both must pass.
        #
        # A NON-ZERO EXIT CODE IS TOLERATED ONLY ON THE KRAKEN STEPS (the ones whose
        # kraken flag is True). The Kraken wrapper can return non-zero even when the job
        # finished; there, if the audit passes, the job counts as done.
        #
        # THE 2026-08-06 FIX - caught on a clean run: that tolerance was being applied to
        # EVERY stage. Stage T returned exit code 3 (P crashed, K/D/I never ran) but,
        # because the expected file had been written as a "P failed" summary, the audit
        # passed and T got a "BITTI" stamp. In the overall summary table a FAILED stage
        # looked SUCCESSFUL, which is contrary to the reason the chain exists.
        # A non-zero exit code on a non-Kraken stage now counts as A FAILURE.
        if tamam and rc != 0 and not kraken:
            tamam = False
            mesaj = ('exit code %s, which is not zero. The output file is '
                     'there but not every substage may have run, which is '
                     'exactly what happened once at the ordering stage. Read '
                     'the stage output above.' % rc)
        if not tamam:
            d.update(durum='DUSTU', sure=sure, cikis=rc, sebep=mesaj)
            json.dump(durum, io.open(dyol, 'w', encoding='utf-8'),
                      ensure_ascii=False, indent=1)
            yaz(u'   OUTPUT CHECK FAILED: %s' % mesaj)
            yaz(u'   exit code: %s | time: %s' % (rc, sn_metni(sure)))
            yaz(u'\n' + u'=' * 78)
            yaz(u'  CHAIN STOPPED at stage %s' % kod)
            yaz(u'  The remaining stages were NOT run. The reason is above.')
            yaz(u'  Fix it and run the same command again; finished stages are skipped.')
            yaz(u'=' * 78)
            return False
        if rc != 0:
            yaz(u'   warning: exit code %s, but the output check passed' % rc)
        d.update(durum='bitti', sure=sure, cikis=rc, sebep=mesaj)
        json.dump(durum, io.open(dyol, 'w', encoding='utf-8'),
                  ensure_ascii=False, indent=1)
        yaz(u'   DONE (%s) - %s' % (sn_metni(sure), mesaj))
    return True


def ozet_yaz(kok, CIKTI, ayar, secili, durum, kesildi):
    """Combines the result of every stage into ONE file.

        It makes no measurement; it reads durum.json and the stages' own outputs. So it
        can be called even when the chain was interrupted, and it shows whatever there
        is.

    """
    yol = os.path.join(CIKTI, '00_TAM_ZINCIR_OZET.md')
    L = []
    A = L.append
    A(u'# The full chain, a combined summary\n')
    A(u'Generated: %s · version %s\n' % (time.strftime('%Y-%m-%d %H:%M'), VERSIYON))
    if kesildi:
        A(u'> **THE CHAIN DID NOT FINISH.** It stopped at the stage marked DUSTU in the table below. The stages after it did not run.\n')

    A(u'## The state of the stages\n')
    A(u'| # | Stage | Group | State | Time | Note |')
    A(u'|---|---|---|---|---|---|')
    for kod, ad, grup, dk, kr, _k, _d in secili:
        d = durum.get(kod, {})
        # The left side is the STORED value and does not change: it is written into
        # durum.json, read back, and the checkpoint decision looks at it. What changes is
        # only the label printed to the screen.
        du = d.get('durum', 'not run')
        etiket = {'bitti': 'FINISHED', 'atlandi': 'SKIPPED',
                  'eksik': 'INCOMPLETE (will be retried)',
                  'DUSTU': '**FAILED**'}.get(du, du)
        A(u'| %s | %s | %s | %s | %s | %s |'
          % (kod, ad, grup, etiket,
             sn_metni(d.get('sure', 0)) if d.get('sure') else '-',
             (d.get('sebep') or '')[:150]))

    bitti = [k for k, v in durum.items() if v.get('durum') == 'bitti']
    atlanan = [k for k, v in durum.items() if v.get('durum') == 'atlandi']
    dusen = [k for k, v in durum.items() if v.get('durum') == 'DUSTU']
    eksikler = [k for k, v in durum.items() if v.get('durum') == 'eksik']
    A(u'\n**Summary:** %d stages finished, %d skipped, %d incomplete, %d failed.\n'
      % (len(bitti), len(atlanan), len(eksikler), len(dusen)))
    if eksikler:
        A(u'The stages marked *incomplete* ran but did not finish; because they never got the `bitti` stamp they are retried by themselves when the same key is pressed again. They finish once the data arrives.')

    if atlanan:
        A(u'## The steps that were skipped\n')
        A(u'These steps did not run. The chain did not stop for it, but their results are not in hand either.\n')
        for k in sorted(atlanan):
            A(u'* **%s**, %s' % (k, durum[k].get('sebep', '')))
        A(u'')
        if not ayar['kraken_var']:
            A(u'### Finishing the Kraken part later\n')
            A(u'You **do not need** to run the chain from the start. Because the skipped stages never got the `bitti` stamp they are retried by themselves on the next run. To run the Kraken part only, press **W**, then **X**, then **Z** in the menu; or in a single command:')
            A(u'```\npython3 verification/full_chain.py --root . --only W,X,Z,S --confirm\n```\n')
            A(u'Stage Z builds the table **from the data in hand, not from scratch**: it leaves the missing columns empty and writes down which ones are missing, and when the data arrives the same key completes the table. That is why the table can be reproduced even if the Kraken measurement is made weeks later.')
            A(u'if kraken2 is not installed: `micromamba create -n mikro -c bioconda -c conda-forge kraken2 bracken`. If it is installed already the environment name may differ; look with `micromamba env list` and pass `--env <name>`, or give the binary\'s full path with `KRAKEN2_BIN=/full/path/kraken2`.')

    A(u'## Where to look\n')
    A(u'| Question | File |')
    A(u'|---|---|')
    bakilacak = [
        (u'Was the chain consistent', 'QUICK_TEST/QUICK_TEST_REPORT.md'),
        (u'What became of the identity claims',
         'IDENTITY_RESULT/KIMLIK_DOGRULAMA_RAPORU.md'),
        (u'The identity of every bin',
         'ALL_IDENTITIES_RESULT/TUM_KUTU_KIMLIK_RAPORU.md'),
        (u'The measurement and recovery result',
         'ALL_STAGES_RESULT/00_BIRLESIK_OZET.md'),
        (u'What should I order', 'SIPARIS_LISTESI.tsv'),
        (u'The Kraken threshold curve',
         '../tools/RESULTS/kraken_esik_A/esik_egrisi.txt'),
        ('The Kraken table that goes into the report', '../tools/0_TESLIM_RAPOR/KRAKEN_KARSILASTIRMA.md'),
    ]
    for soru, dy in bakilacak:
        t = os.path.join(kok, dy)
        A(u'| %s | `%s`%s |' % (soru, dy, u'' if os.path.exists(t) else u' *(yok)*'))

    A(u'\n## What this summary is not\n')
    A(u'This file shows that the stages **ran** and that their outputs are **not empty**. It does not show that the measurements are **right**; correctness is argued in each stage\'s own report. The order decision is made by reading the reports named in the table above, not by reading this summary.')

    io.open(yol, 'w', encoding='utf-8').write(u'\n'.join(L) + u'\n')
    return yol


def main():
    p = argparse.ArgumentParser(description='The full chain: 8 H W I G T X Y '
                                            'Z S')
    p.add_argument('--root', dest='kok', default='.')
    p.add_argument('--confirm', dest='onayla', action='store_true',
                   help='start without asking for a confirmation, which is '
                        'what the menu does')
    p.add_argument('--rerun', dest='yeniden', action='store_true',
                   help=u'durum.json is reset and everything runs from scratch')
    p.add_argument('--from-scratch', dest='sifirdan', action='store_true',
                   help='A CLEAN RUN: it invalidates the checkpoints not only '
                        'of this script but of every stage it CALLS. Nothing '
                        'is deleted; things are moved into time stamped '
                        'directories.')
    p.add_argument('--only', dest='yalniz', default='',
                   help=u'these stages only, comma-separated, e.g. 8,S')
    p.add_argument('--skip', dest='atla', default='',
                   help=u'skip these stages, comma-separated')
    p.add_argument('--pluspfp', default=os.environ.get('PLUSPFP', ''),
                   help='the path of the PlusPFP database; when it is omitted '
                        'the corresponding step is skipped')
    p.add_argument('--db-path', dest='vt', default=os.environ.get('VT_A', ''),
                   help=u'Kraken2 database path (default ~/k2db, then the tool scans the disk)')
    p.add_argument('--env', dest='ortam', default=os.environ.get('ORTAM', ''),
                   help='the micromamba or conda environment name, mikro by '
                        'default. Give it here when kraken2 lives in another '
                        'environment; to see the names: micromamba env list')
    p.add_argument('--dry-run', dest='kuru', action='store_true',
                   help='show the plan and the checks WITHOUT RUNNING the '
                        'commands')
    p.add_argument('--plan', action='store_true',
                   help='print the plan only and exit; nothing is run. The '
                        'driver calls this first, takes the confirmation '
                        'itself, and then starts the real run with --confirm, '
                        'so that the question is asked outside WSL and does '
                        'not depend on how stdin is passed through.')
    a = p.parse_args()

    kok = os.path.abspath(a.kok)
    if not os.path.isdir(os.path.join(kok, 'screening')):
        sys.stderr.write(
            u'ERROR: no screening package in %s. This script runs from the project root, the directory that contains verification/full_chain.py.\n' % kok)
        return 1

    CIKTI = os.path.join(kok, CIKTI_ADI)
    os.makedirs(CIKTI, exist_ok=True)
    dyol = os.path.join(CIKTI, 'durum.json')
    durum = {}
    if os.path.exists(dyol) and not a.yeniden:
        try:
            durum = json.load(io.open(dyol, encoding='utf-8'))
        except Exception:
            durum = {}
    if a.yeniden:
        # CAUTION: the file IS NOT DELETED, it IS OVERWRITTEN.
        # The reason was measured: on mounted directories, and on some Windows
        # installations, there may be no permission to delete (Operation not permitted).
        # Trying to delete made the --rerun option completely unusable. Writing an empty
        # dictionary does the same job and works on every file system.
        durum = {}
        json.dump(durum, io.open(dyol, 'w', encoding='utf-8'),
                  ensure_ascii=False, indent=1)

    if a.sifirdan:
        # A CLEAN RUN. --rerun resets only THIS script's durum.json; but the stages it
        # calls have THEIR OWN checkpoints, and while those stand, a stage skips itself
        # internally as "already finished". To really run from scratch, those have to be
        # invalidated too.
        # NOTHING IS DELETED: there may be no permission to delete on a mounted directory,
        # and we do not want to lose data anyway. The directories are renamed with a
        # timestamp.
        import datetime as _dt
        damga = _dt.datetime.now().strftime('%Y%m%d_%H%M%S')
        hedefler = ['ONE_PROTOCOL_RESULT', 'RECOVERY_RESULT', 'VERIFICATION_RESULT',
                    'IDENTITY_RESULT', 'ALL_IDENTITIES_RESULT', 'ALL_STAGES_RESULT',
                    'QUICK_TEST', 'ACCESS_RESULT', CIKTI_ADI]
        tasinan = []
        for h in hedefler:
            y = os.path.join(kok, h)
            if os.path.isdir(y):
                yeni_ad = os.path.join(kok, '%s_PREVIOUS_%s' % (h, damga))
                try:
                    os.rename(y, yeni_ad)
                    tasinan.append(h)
                except OSError as e:
                    print(u'  WARNING: could not move %s (%s). This stage may resume from its own checkpoint.' % (h, e))
        ali_is = os.path.join(kok, '..', 'tools', 'RESULTS', 'kraken_esik_A')
        if os.path.isdir(ali_is):
            try:
                os.rename(ali_is, ali_is + '_PREVIOUS_' + damga)
                tasinan.append('tools/RESULTS/kraken_esik_A')
            except OSError:
                pass
        print(u'\nCLEAN RUN: %d result directories were moved aside. Nothing was deleted.' % len(tasinan))
        for h in tasinan:
            print(u'    %s -> %s_PREVIOUS_%s' % (h, h, damga))
        print(u'  The old results are still there if you want to look at them.\n')
        durum = {}
        os.makedirs(CIKTI, exist_ok=True)

    ayar = kraken_ortami(kok, a.pluspfp, a.vt, a.ortam)

    yalniz = [x.strip().upper() for x in a.yalniz.split(',') if x.strip()]
    atla = [x.strip().upper() for x in a.atla.split(',') if x.strip()]
    secili = [s for s in ASAMALAR
              if (not yalniz or s[0] in yalniz) and s[0] not in atla]
    if not secili:
        sys.stderr.write(u'ERROR: no stages left after --only / --skip.\n')
        return 1

    gunluk = io.open(os.path.join(CIKTI, 'kosu_gunlugu.txt'), 'a', encoding='utf-8')

    def yaz(s=u''):
        print(s, flush=True)
        gunluk.write(s + u'\n')

    yaz(u'\nFULL CHAIN  version %s  %s' % (VERSIYON, time.strftime('%Y-%m-%d %H:%M')))
    plan_yaz(yaz, secili, ayar, durum)

    if a.plan:
        gunluk.close()
        return 0

    if not a.onayla and not a.kuru:
        try:
            c = input(u'\n  Baslatilsin mi? [E/h] ').strip().lower()
        except EOFError:
            c = 'e'
        if c not in ('', 'e', 'evet', 'y', 'yes'):
            yaz(u'  Cancelled. Nothing was run.')
            return 0

    t0 = time.time()
    tamamlandi = calistir(kok, ayar, secili, durum, dyol, gunluk, yaz, a.kuru)
    yol = ozet_yaz(kok, CIKTI, ayar, secili, durum, not tamamlandi)

    yaz(u'\n' + u'=' * 78)
    yaz(u'  TOTAL TIME: %s' % sn_metni(time.time() - t0))
    yaz(u'  COMBINED SUMMARY: %s' % yol)
    yaz(u'=' * 78)
    gunluk.close()
    return 0 if tamamlandi else 3


if __name__ == '__main__':
    sys.exit(main() or 0)
