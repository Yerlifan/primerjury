# -*- coding: utf-8 -*-
"""OPTION 9 - RUN EVERYTHING IN ORDER.

It runs every stage in the right ORDER, so the user can click once in the evening
and go to bed. The order is not arbitrary but follows the dependencies:

  1. The self test             are the engines right (nothing runs if it fails)
  2. Consensus regeneration    THE SEARCH uses the backbone -> this comes FIRST
  3. Panel re-measurement      full depth, the corrected engine
  4. Membership audit          which number depends on which definition
  5. The full search           the problem targets (with step 2's consensus)
  6. THE COMBINED SUMMARY      what changed, what dropped, where a new candidate is

A checkpoint is written as each stage finishes. If it is interrupted, choosing (9)
again continues from the stage it stopped at; the finished stages ARE NOT RE-RUN.

THE CONSENSUS -> SEARCH LINK
----------------------------
A regenerated consensus is used as the search backbone if it passes THE QUALITY
GATE (a low N proportion and no divergence between the two methods). If it does
not pass, the old consensus is used. Which backbone was used is written into the
report for every target; no change is made silently.

"""
# -------------------------------------------------------------------------
# run_all.py - the batch flow that runs seven stages in dependency order, and the
#            summary producer that combines every stage's output into one file.
#
# INPUT  : the paths in config.py; kontrol/hepsi_durum.json (which stage has
#          finished); panel_yeniden_olcum.tsv, panel_kutu_duzeyi.tsv,
#          uyelik_duyarlilik.tsv, konsensus_yeniden_uretim.tsv, adaylar.tsv and
#          kontrol/hedef_*.json under SCREENING_RESULT.
#          It calls the stages itself: build_canonical.py (as a subprocess),
#          build_consensus.calistir, panel_measurement.calistir,
#          membership_check.calistir, __main__.aramayi_kos.
# OUTPUT : SCREENING_RESULT/00_OZET_HEPSI.md (the path ozet_yaz returns);
#          kontrol/hepsi_durum.json (the stage state). calistir() returns an exit
#          code: 0 fine, 1 the user stopped it, 2 a gate did not pass.
# CALLED BY: verification/full_chain.py key 9 (--mode hepsi -> HP.calistir) and key
#          S (--mode ozet -> HP.yalniz_ozet). Besides that, because the yon_kapisi()
#          function is called by __main__.py, panel_measurement.py,
#          membership_check.py and build_consensus.py, it also runs indirectly on
#          keys 1, 2, 3, 4, 5, 6 and 7.
# -------------------------------------------------------------------------
import os, json, time, csv
from . import config as C
from . import checks

ASAMA_YOLU = None

# the quality gate for a new consensus to be usable as a backbone
KONSENSUS_N_UST = 20.0        # above this N percentage the old consensus is used
KONSENSUS_AYRILMA_UST = 15    # the cap on the number of columns where the two methods disagree


def _durum_yolu():
    return os.path.join(C.KONTROL, 'hepsi_durum.json')


def durum_oku():
    p = _durum_yolu()
    if os.path.exists(p):
        try:
            d = json.load(open(p, encoding='utf-8'))
            if checks.ayar_uyuyor(d):
                return d
        except Exception:
            pass
    return dict(bitmis=[], baslangic=time.strftime('%Y-%m-%d %H:%M:%S'), ciktilar={})


def durum_yaz(d):
    checks.hazirla()
    d['_ayar'] = dict(checks.AYAR)
    p = _durum_yolu()
    with open(p + '.gecici', 'w', encoding='utf-8') as fh:
        json.dump(d, fh, ensure_ascii=False, indent=1, default=str)
    os.replace(p + '.gecici', p)


# ---------------------------------------------------------------- konsensus kapisi
def konsensus_kalite():
    """The quality table of the regenerated consensuses: bin -> (passed, reason).

        NOTE: this function was deleted by accident during an intermediate edit while
        ozet_yaz() was still calling it, so the LAST step of a 5 hour run failed with
        NameError. It was put back, and ozet_yaz can now also be run on its own (a menu
        option and --mode ozet).

    """
    tsv = os.path.join(C.CIKTI, 'konsensus_yeniden_uretim.tsv')
    out = {}
    if not os.path.exists(tsv):
        return out
    try:
        for row in csv.DictReader(open(tsv, encoding='utf-8'), delimiter='\t'):
            try:
                n = float(row.get('N_yuzde') or 0)
                ayr = int(row.get('yontemler_ayrildi') or 0)
            except (ValueError, TypeError):
                continue
            if n > KONSENSUS_N_UST:
                out[row['kutu']] = (False, 'N orani %%%.1f > %%%.0f' % (n, KONSENSUS_N_UST))
            elif ayr > KONSENSUS_AYRILMA_UST:
                out[row['kutu']] = (False, 'iki yontem %d sutunda ayrildi' % ayr)
            else:
                out[row['kutu']] = (True, 'N %%%.1f, ayrilma %d' % (n, ayr))
    except Exception:
        pass
    return out


# ---------------------------------------------------------------- yon kapisi
def yon_kapisi(yaz, asama):
    """Verifies BEFORE a stage starts that its input is canonical.

        An orientation error is silent: on a reversed consensus, in-silico PCR returns 0
        products with no warning at all and a whole night produces the result "no product
        on any target". So every stage tests its own input first; if it does not pass, THE
        STAGE DOES NOT START.

        Returns: (passed, message_lines)

    """
    # WHY THE ORIENTATION IS READ FROM THE CANONICAL SOURCE
    # Nanopore reads are bidirectional; a bin's consensus may have been produced SENSE or
    # ANTISENSE depending on which direction the reads were anchored. The raw "consensus
    # sequences" directory is therefore MIXED orientation (measured: 71 antisense, 27
    # sense). Because the engines scan the sequence only on the plus strand, in-silico
    # PCR on a reversed consensus loses ALL the products, and does so with no warning at
    # all: the output is "no product", not "the orientation is wrong". So the consensuses
    # are read from a single canonical directory, and every stage verifies here, BEFORE
    # it starts, that its own input is canonical.
    # The gate tests three things in order: does orientation.py pass its own self test,
    # can the canonical index be read, and is none of the consensuses in the index
    # ANTISENSE. If any of the three fails, the stage is not started.
    from . import orientation as Y
    from . import targets as H
    m = []
    try:
        h = Y.kendini_sina()
        if h:
            return False, ['orientation.py kendi sinavindan GECEMEDI: ' + '; '.join(h)]
        kons = H.konsensusler()
    except Exception as e:
        return False, ['Kanonik konsensus okunamadi: %s' % e,
                       'To generate it:  python3 screening/build_canonical.py --root . ']
    if not kons:
        return False, ['Kanonik konsensus indeksi BOS.']
    ters = [k['kutu'] for k in kons
            if Y.tespit(k['dizi'], Y.sinifi(k['kutu']))[0] == 'ANTISENSE']
    if ters:
        return False, ['%d konsensus TERS yonde: %s' % (len(ters), ', '.join(ters[:5])),
                       'Ters yonlu konsensuste urunlerin TAMAMI kaybolur.',
                       'To fix it:  python3 screening/build_canonical.py --root . --rerun']
    m.append('yon kapisi [%s]: %d kutu, hepsi SENSE - GECTI' % (asama, len(kons)))
    return True, m


def kanonik_kos(yaz, sure, oncelik='ozgun'):
    """Calls build_canonical.py in the same process (so that no separate command is needed)."""
    import subprocess, sys as _sys
    betik = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'build_canonical.py')
    if not os.path.exists(betik):
        return False, 'build_canonical.py bulunamadi'
    komut = [_sys.executable, betik, '--root', C.KOK, '--priority', oncelik]
    yaz('  > %s' % ' '.join(komut[1:]))
    t = time.time()
    try:
        r = subprocess.run(komut, capture_output=True, text=True, timeout=7200)
    except subprocess.TimeoutExpired:
        return False, 'kanonik uretim 2 saati asti'
    for satir in (r.stdout or '').strip().split('\n')[-6:]:
        if satir.strip():
            yaz('    ' + satir.strip())
    if r.returncode != 0:
        return False, (r.stderr or '')[-400:]
    # onbellekteki konsensus listesi bayatlamasin
    from . import targets as H
    if hasattr(H, '_kons_onbellek'):
        H._kons_onbellek = None
    return True, 'kanonik uretim tamam (%s)' % sure(time.time() - t)


# ---------------------------------------------------------------- asamalar
def calistir(yaz, sure, cizgi, a):
    checks.hazirla()
    d = durum_oku()
    t0 = time.time()

    cizgi('=')
    yaz(u'  RUN EVERYTHING IN ORDER')
    cizgi('=')
    yaz(u'  Stages and estimated times:')
    yaz(u'    1. Self-test                          ~1 min')
    yaz(u'    2. Canonical consensus (orientation)  ~2-5 min   <- ORIENTATION GATE')
    yaz(u'    3. Consensus regeneration             ~30-60 min')
    yaz(u'    4. Move new consensus to canonical    ~2-5 min   <- ORIENTATION GATE')
    yaz(u'    5. Panel re-measurement (full depth)  ~1-2 h')
    yaz(u'    6. Membership audit                   ~20-40 min')
    yaz(u'    7. Comprehensive search               ~1-3 h')
    yaz(u'       + combined summary report          ~1 min')
    yaz(u'  TOTAL: roughly 3-7 hours. Can be left running overnight.')
    yaz('')
    if d['bitmis']:
        yaz(u'  RESUMING - finished stages: %s' % ', '.join(d['bitmis']))
        yaz('')
    yaz(u'  State is written to disk after every stage. If interrupted, choose (9) again;')
    yaz(u'  finished stages are not re-run.')
    cizgi('=')

    from . import self_test

    # ---- 1 kendini sinama (kapi gorevinde)
    yaz(u'\n[STAGE 1/7] Self-test')
    if not getattr(a, 'sinama_atla', False) and not self_test.calistir(yaz):
        yaz(u'\nSELF-TEST FAILED - no stage was run.')
        yaz(u'Fix the line marked *** FAILED *** above first.')
        return 2

    # ---- 2 KANONIK KONSENSUS (referans onceligi) - her seyin girdisi
    yaz(u'\n[STAGE 2/7] Canonical consensus generation (orientation normalisation)')
    yaz(u'  No measurement is trustworthy without this step: the raw consensus directory')
    yaz(u'  is MIXED orientation, and on a reverse-oriented consensus EVERY product is lost.')
    yaz(u'  Source: "consensus sequences" (the strict set the panel was built on).')
    ok, msj = kanonik_kos(yaz, sure, oncelik='ozgun')
    if not ok:
        yaz(u'\n  THE CANONICAL PRODUCTION FAILED: %s' % msj)
        yaz(u'  The following stages were NOT STARTED (so that no wrong result is produced).')
        return 2
    yaz('  ' + msj)
    ok, msj = yon_kapisi(yaz, 'kanonik uretim sonrasi')
    for x in msj:
        yaz('  ' + x)
    if not ok:
        yaz(u'\n  THE ORIENTATION GATE DID NOT PASS - the later stages WERE NOT STARTED.')
        return 2

    # ---- 3 consensus regeneration
    if 'konsensus' in d['bitmis']:
        yaz(u'\n[STAGE 3/7] Consensus regeneration - SKIPPED (already finished)')
    else:
        yaz(u'\n[STAGE 3/7] Consensus regeneration (from the raw reads)')
        from . import build_consensus
        try:
            y = build_consensus.calistir(yaz, sure, yalniz=getattr(a, 'hedef', None))
            d['bitmis'].append('konsensus')
            d['ciktilar']['konsensus'] = y
            durum_yaz(d)
        except KeyboardInterrupt:
            yaz(u'\nSTOPPED - continue with (9).')
            return 1
        except Exception as e:
            yaz(u'\n  STAGE 3 ERROR: %s' % e)
            yaz(u'  No new consensus could be generated; the canonical reference set stays in use.')

    # ---- 4 YENI konsensusleri kanonige al (tek kaynak korunur)
    yaz(u'\n[STAGE 4/7] Moving new consensus sequences to canonical (--priority new)')
    yaz(u'  Newly generated consensus sequences are NOT used directly. They are first')
    yaz(u'  written into the canonical directory through orientation normalisation, so there stays one source.')
    ok, msj = kanonik_kos(yaz, sure, oncelik='yeni')
    if not ok:
        yaz(u'  The canonical update could not be made (%s); the reference set stays valid.' % msj)
    else:
        yaz('  ' + msj)
    ok, msj = yon_kapisi(yaz, 'yeni konsensus sonrasi')
    for x in msj:
        yaz('  ' + x)
    if not ok:
        yaz(u'\n  THE ORIENTATION GATE DID NOT PASS - the later stages WERE NOT STARTED.')
        return 2

    # ---- 5 panel re-measurement
    if 'panel' in d['bitmis']:
        yaz(u'\n[STAGE 5/7] Panel re-measurement - SKIPPED (already finished)')
    else:
        yaz(u'\n[STAGE 5/7] Panel re-measurement (FULL DEPTH)')
        yaz(u'  This is the longest stage: first every fastq file is read (10-25 min),')
        yaz(u'  then 21 pairs are measured against two criteria. One to two hours in total.')
        from . import panel_measurement
        try:
            y = panel_measurement.calistir(yaz, sure, okuma_sayisi=0,
                                     yalniz=getattr(a, 'hedef', None))
            d['bitmis'].append('panel')
            d['ciktilar']['panel'] = y
            durum_yaz(d)
        except KeyboardInterrupt:
            yaz(u'\nSTOPPED - continue with (9).')
            return 1
        except Exception as e:
            yaz(u'\n  STAGE 5 ERROR: %s' % e)

    # ---- 6 uyelik denetimi
    if 'uyelik' in d['bitmis']:
        yaz(u'\n[STAGE 6/7] Membership audit - SKIPPED (already finished)')
    else:
        yaz(u'\n[STAGE 6/7] Membership audit and sensitivity analysis')
        from . import membership_check
        try:
            y = membership_check.calistir(yaz, sure,
                                         okuma_sayisi=C.NUMUNE_OKUMA_SAYISI,
                                         yalniz=getattr(a, 'hedef', None))
            d['bitmis'].append('uyelik')
            d['ciktilar']['uyelik'] = y
            durum_yaz(d)
        except KeyboardInterrupt:
            yaz(u'\nSTOPPED - continue with (9).')
            return 1
        except Exception as e:
            yaz(u'\n  STAGE 6 ERROR: %s' % e)

    # ---- 7 kapsamli arama
    if 'arama' in d['bitmis']:
        yaz(u'\n[STAGE 7/7] Full search - SKIPPED (already finished)')
    else:
        yaz(u'\n[STAGE 7/7] The comprehensive search (problem targets)')
        from . import __main__ as M
        try:
            rc = M.aramayi_kos(a, yaz, sure, cizgi, mod='devam')
            if rc == 0:
                d['bitmis'].append('arama')
                durum_yaz(d)
        except KeyboardInterrupt:
            yaz(u'\nSTOPPED - continue with (9).')
            return 1
        except Exception as e:
            yaz(u'\n  STAGE 7 ERROR: %s' % e)

    # ---- 6 the combined summary
    yaz(u'\n[AFTER STAGE 7/7] Combined summary report')
    yol = ozet_yaz(d, sure, time.time() - t0)
    d['bitmis'] = list(dict.fromkeys(d['bitmis']))
    d['bitis'] = time.strftime('%Y-%m-%d %H:%M:%S')
    durum_yaz(d)

    cizgi('=')
    yaz(u'  ALL DONE (%s)' % sure(time.time() - t0))
    yaz('')
    yaz(u'  READ THIS FIRST:')
    yaz('    %s' % yol)
    yaz('')
    yaz(u'  Detailed reports:')
    for ad, ys in d.get('ciktilar', {}).items():
        for y in (ys or []):
            if str(y).endswith('.md'):
                yaz('    %s' % y)
    cizgi('=')
    return 0


# ---------------------------------------------------------------- the summary report
def _oku_tsv(ad):
    p = os.path.join(C.CIKTI, ad)
    if not os.path.exists(p):
        return []
    try:
        return list(csv.DictReader(open(p, encoding='utf-8'), delimiter='\t'))
    except Exception:
        return []


def _wilson(k, n, z=1.96):
    import math
    if n == 0:
        return (0.0, 1.0)
    pp = k / n
    d = 1 + z * z / n
    c = (pp + z * z / (2 * n)) / d
    sd = z * math.sqrt(pp * (1 - pp) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - sd), min(1.0, c + sd))


def _asama5_yeniden(asgari=None):
    """Recompute stage 5's "worst single bin" ratio from the RECORDED bin counts, with
        THE SAME absolute threshold as stage 6.

        The reason: during the run the threshold was "half the largest bin", and at full
        depth only 1-5 of the 10-33 competitor bins entered the measurement. Measuring the
        same pair with 300 reads let all of them in. Most of the difference between the two
        stages comes from that. Because the raw per bin counts are recorded, it can be
        corrected without repeating the 5 hour run.

    """
    # The threshold MUST be ABSOLUTE: a relative threshold ("half the largest bin", for
    # instance) takes different bins into the measurement as the depth changes, so the
    # numbers of two stages stop measuring the same thing and become incomparable.
    # Here stage 5's recorded bin counts are recomputed with THE SAME absolute threshold
    # as stage 6 (C.ENKOTU_ASGARI_OKUMA); the difference left comes only from the read
    # depth. The Wilson LOWER bound is used on the member side and the UPPER bound on the
    # competitor side - both conservative, so the ratio never comes out larger than it is.
    if asgari is None:
        asgari = getattr(C, 'ENKOTU_ASGARI_OKUMA', 150)
    rows = _oku_tsv('panel_kutu_duzeyi.tsv')
    if not rows:
        return {}
    from collections import defaultdict
    g = defaultdict(lambda: {'uye': [], 'rakip': []})
    for r in rows:
        if not (r.get('olcut', '') or '').startswith('<=1'):
            continue
        try:
            pv = int(r['YENI_motor_urun_veren']); nv = int(r['okuma'])
        except (ValueError, KeyError):
            continue
        g[r['hedef']][r['grup']].append((r['kutu'], pv, nv))
    out = {}
    for h, d in g.items():
        if not d['uye']:
            continue
        uye_alt = min(_wilson(pv, nv)[0] for _, pv, nv in d['uye'])
        rp = sum(pv for _, pv, _ in d['rakip']); rn = sum(nv for _, _, nv in d['rakip'])
        hav = _wilson(rp, rn)[1] if rn else None
        en = None
        giren = 0
        for kadi, pv, nv in d['rakip']:
            if nv < asgari:
                continue
            giren += 1
            hi = _wilson(pv, nv)[1]
            if en is None or hi > en[1]:
                en = (kadi, hi)
        out[h] = dict(
            uye_alt=round(100 * uye_alt, 3),
            kat_havuz=(round(uye_alt / hav, 2) if hav else None),
            kat_enkotu=(round(uye_alt / en[1], 2) if en and en[1] > 0 else None),
            enkotu_kutu=(en[0] if en else ''),
            rakip_kutu=len(d['rakip']), giren_kutu=giren,
            uye_kutu=len(d['uye']),
            derinlik=sum(nv for _, _, nv in d['uye'] + d['rakip']))
    return out


def _arama_kontrolleri():
    """Stage 7's target checkpoint files: the baseline plus the candidates (under THE SAME conditions)."""
    import glob
    out = []
    for f in sorted(glob.glob(os.path.join(C.KONTROL, 'hedef_*.json'))):
        try:
            out.append(json.load(open(f, encoding='utf-8')))
        except Exception:
            pass
    return out


def _oran_metni(v, rakip_kutu=None):
    # THE DISCRIMINATION RATIO IS UNDEFINED ON UNIVERSAL TARGETS. The discrimination
    # ratio divides the member side's Wilson lower bound by the competitor side's Wilson
    # upper bound. On a universal target (Arke_universal, Bakteri_universal,
    # Mantar_universal F1/F2) the competitor set approaches empty by the definition of
    # the membership; the denominator goes to zero and the ratio becomes undefined. On
    # paper, 0.00 and 117 million stand side by side in the same column, and neither
    # number is a measurement. So no ratio is printed here: if there is no competitor bin
    # at all it says "rakip kutusu yok", and if the number exploded because the
    # denominator approached zero it says "olcusuz (rakip ~0)". Those targets are
    # evaluated by coverage (in how many bins do they give a product) and by the outside
    # the domain proportion. THIS IS NOT LOWERING the threshold; it is replacing the
    # measure itself with the right one.
    if rakip_kutu == 0:
        return 'rakip kutusu yok'
    if v in (None, '', 'None'):
        return '-'
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if f > 100000:
        return 'olcusuz (rakip ~0)'
    return '%.2f' % f


def ozet_yaz(durum, sure, gecen):
    os.makedirs(C.CIKTI, exist_ok=True)
    yol = os.path.join(C.CIKTI, '00_OZET_HEPSI.md')
    L = []; A = L.append

    panel = _oku_tsv('panel_yeniden_olcum.tsv')
    uyelik = _oku_tsv('uyelik_duyarlilik.tsv')
    kons = _oku_tsv('konsensus_yeniden_uretim.tsv')
    aday = _oku_tsv('adaylar.tsv')
    aramalar = _arama_kontrolleri()
    yeniden5 = _asama5_yeniden()
    asgari = getattr(C, 'ENKOTU_ASGARI_OKUMA', 150)

    A(u'# The summary, everything that was run, in order')
    A('')
    A('Baslangic: %s · bitis: %s · sure: %s'
      % (durum.get('baslangic', '?'), time.strftime('%Y-%m-%d %H:%M:%S'), sure(gecen)))
    A('')
    A('Kosulan asamalar: %s' % (', '.join(durum.get('bitmis', [])) or 'yok'))
    A('')
    eksik = [ad for ad, v in (('panel_yeniden_olcum.tsv', panel),
                              ('uyelik_duyarlilik.tsv', uyelik),
                              ('konsensus_yeniden_uretim.tsv', kons),
                              ('adaylar.tsv', aday)) if not v]
    if eksik:
        A(u'> **Missing output files:** %s. The sections that use them will look empty.'
          % ', '.join('`%s`' % e for e in eksik))
        A('')

    # ============================================================ 0. EN ONEMLI
    A(u'## 0. READ THIS FIRST: which number was measured under which conditions')
    A('')
    A(u'In this run the same pair was measured **in three separate places** and the numbers **came out different from one another**. The reason was found; both are right but they are **not comparable**:')
    A('')
    A(u'| | Stage 5 (the panel remeasurement) | Stages 6 and 7 (membership plus the search baseline) |')
    A('|---|---|---|')
    A(u'| read depth | **FULL** (a median of about 3 000 per bin, at most about 46 000) | **300 reads per bin** |')
    A(u'| mismatch criterion | `<=1` and `<=3` (two separate rows) | `<=1` |')
    A(u'| membership definition | the same (`hedef_uyelik.tsv`) | the same |')
    A(u'| engine | the same (`read_engine.py`) | the same |')
    A('')
    A('**Farkin iki sebebi var:**')
    A('')
    A(u'1. **The Wilson interval narrows with depth.** The LOWER bound is used for a member and the UPPER bound for a competitor; as the reads grow the lower bound rises, the upper bound falls, and the ratio **grows**. That is why the ratios come out different although the member percentages below them are almost the same: what changes is not the measurement but the **margin of uncertainty**.')
    A(u'2. **The "worst single bin" measure was sliding with depth during the run.** The threshold was "half of the largest bin"; at full depth, with the largest bin at about 46 000 reads, the threshold became about 23 000 and **only 1 to 5 of the 10 to 33 competitor bins** entered the measurement. Measured with 300 reads the threshold became 150 and **all of them** entered. So stage 5\'s "worst bin" number really meant "**the deepest bin**", and the true worst competitor may have been left outside.')
    A('')
    A(u'That second problem **was fixed** (the threshold is absolute now: %d reads). Because the raw per bin counts are recorded, the correct values below were computed **without rerunning stage 5**.' % asgari)
    A('')
    A(u'### The table made comparable')
    A('')
    A(u'The `Stage 5 (corrected)` column recomputes stage 5\'s own full depth numbers **with the same bin threshold** as stage 6. The difference that is left comes from depth alone.')
    A('')
    A(u'| target | S5 pool (full) | S5 worst (as run) | **S5 worst (corrected)** | S6 and S7 pool (300) | S6 and S7 worst (300) | competitor bins entering S5 |')
    A('|---|---|---|---|---|---|---|')
    p1 = {r['hedef']: r for r in panel if (r.get('olcut') or '').startswith('<=1')}
    ubak = {}
    for r in uyelik:
        if (r.get('tanim') or '').startswith('A.'):
            ubak.setdefault(r['hedef'], r)
    for h in sorted(set(list(p1.keys()) + list(yeniden5.keys()))):
        a5 = p1.get(h, {})
        yd = yeniden5.get(h, {})
        u6 = ubak.get(h, {})
        A('| %s | %s | %s | **%s** | %s | %s | %s/%s |' % (
            h[:40], _oran_metni(a5.get('ayrim_havuz_x')),
            _oran_metni(a5.get('ayrim_en_kotu_x')),
            _oran_metni(yd.get('kat_enkotu')),
            _oran_metni(u6.get('ayrim_havuz_x')),
            _oran_metni(u6.get('ayrim_en_kotu_x')),
            yd.get('giren_kutu', '?'), yd.get('rakip_kutu', '?')))
    A('')
    A(u'> **Which one should you look at?** Use the **`S5 worst (corrected)`** column for a decision: it is computed on the deepest data, covering every competitor bin, in the conservative direction. The `S6 and S7` columns are for ranking the candidates **among themselves**; since they were all measured under the same 300 read condition they can be compared with one another, but not with S5.')
    A('')

    # ============================================================ 1. yeni aday
    A(u'## 1. Was a candidate better than the panel\'s current pair found?')
    A('')
    if not aramalar:
        A(u'*The search checkpoint files were not found.*')
    else:
        A(u'For every target the **current pair** in the panel and the **best candidate** were measured **under THE SAME conditions** (300 reads per bin, `<=1` mismatch, the same membership, the same engine), so these columns can be compared directly.')
        A('')
        A(u'| # | target | current pair | best candidate | is it better | gain | 10x threshold | ARMS | grid cell |')
        A('|---|---|---|---|---|---|---|---|---|')
        iyi = kotu = anlamsiz = belirsiz = 0
        sira = 0
        for v in sorted(aramalar, key=lambda x: x.get('hedef', '')):
            h = v.get('hedef', '?')
            t = v.get('panel_olcum') or {}
            ad = v.get('adaylar') or []
            sira += 1
            # rakip kutusu yoksa (evrensel hedefler) ORAN ANLAMSIZDIR
            rakip_yok = not (t.get('rakip') or [])
            tb = t.get('kat_enkotu')
            if tb in (None, ''):
                tb = t.get('kat_havuz')
            en = None
            for c in ad:
                v2 = c.get('numune_kat_enkotu') or c.get('numune_kat_havuz')
                try:
                    v2 = float(v2)
                except (TypeError, ValueError):
                    continue
                if en is None or v2 > en[0]:
                    en = (v2, c)
            hucre = (en[1].get('izgara_hucresi') if en else '') or ''
            hucre = hucre.replace('|', '/')[:44]
            arms = ('EVET' if (en and en[1].get('arms')) else 'hayir')

            if rakip_yok:
                anlamsiz += 1
                A(u'| %d | %s | *no competitor bin* | *no competitor bin* | **the ratio is meaningless** | - | - | %s | %s |'
                  % (sira, h[:38], arms, hucre))
                continue
            try:
                tbf = float(tb)
            except (TypeError, ValueError):
                tbf = None
            if en is None or tbf is None:
                belirsiz += 1
                A('| %d | %s | %s | %s | **olculemedi** | - | - | %s | %s |'
                  % (sira, h[:38], _oran_metni(tb), _oran_metni(en[0] if en else None),
                     arms, hucre))
                continue

            daha = en[0] > tbf * 1.05
            esikte = en[0] >= C.AYRIM_ESIK   # dCq 3 -> 8,00x (tek kaynak)
            if daha:
                iyi += 1
            else:
                kotu += 1
            kazanc = ('%.1fx' % (en[0] / tbf)) if tbf > 0 else ('yeni (taban 0)' if en[0] > 0 else '-')
            A('| %d | %s | %s | %s | %s | %s | %s | %s | %s |' % (
                sira, h[:38], _oran_metni(tb), _oran_metni(en[0]),
                '**EVET**' if daha else 'hayir', kazanc,
                'gecti' if esikte else '**ALTINDA**', arms, hucre))
        A('')
        A(u'**Result: a better candidate EXISTS for %d targets, the current pair should be kept for %d, the ratio is meaningless for %d (universal, with no competitor bin), and %d could not be measured.**'
          % (iyi, kotu, anlamsiz, belirsiz))
        A('')
        A(u'> **Do not skip the `10x threshold` column.** "Better" DOES NOT mean "enough": a candidate can beat the current pair and still stay below the 10x threshold. If the baseline is 0.12x and the candidate 13.5x it is both better and above the threshold; if the baseline is 0.19x and the candidate 1.04x it is better but **still unusable**.')
        A('')
        A(u'> **On the universal targets** (Arke_universal, Bakteri_universal, Mantar_universal F1 and F2) there is no competitor bin by the definition of membership; the discrimination ratio is undefined and those targets have to be judged by **coverage** (in how many bins they give a product). Do not look at the ratio columns.')
        A('')
        A(u'> The cost of every candidate (which rule was relaxed, whether ARMS was needed, how the product length affects the protocol, whether it can be run at 60 C) is written under the target headings inside `KAPSAMLI_ARAMA_RAPORU.md`.')
    A('')

    # ============================================================ 2. panel degisimi
    A(u'## 2. What changed in the panel\'s numbers (stage 5)')
    A('')
    if not panel:
        A(u'*The panel remeasurement was not run, or its output was not found.*')
    else:
        degisen = [r for r in panel if (r.get('DEGISIM') or '').startswith(('YUKARI', 'ASAGI'))]
        A(u'- Rows measured: **%d** (each pair under two criteria)' % len(panel))
        A(u'- Rows whose value **changed by more than %%30**: **%d**' % len(degisen))
        kayip = []
        for r in panel:
            try:
                k = float(r.get('ESKI_motor_kayip_uye_%') or 0)
                if k > 0:
                    kayip.append((k, r['hedef'], r['olcut']))
            except ValueError:
                pass
        if kayip:
            kayip.sort(reverse=True)
            A(u'- The targets the old (faulty) read engine **missed most**:')
            A('')
            A(u'| target | criterion | what the old engine lost |')
            A('|---|---|---|')
            for k, h, o in kayip[:8]:
                A('| %s | %s | %%%.1f |' % (h, o, k))
        A('')
        A('Ayrinti: `PANEL_YENIDEN_OLCUM.md`, `panel_kutu_duzeyi.tsv`')
    A('')

    # ============================================================ 3. below threshold
    A(u'## 3. The pairs that fell BELOW the 10x threshold')
    A('')
    dusen = [r for r in panel if r.get('ESIK_ALTINA_DUSTU')]
    if not panel:
        A(u'*The panel remeasurement was not run.*')
    elif not dusen:
        A(u'In the new measurement **no pair falls** below the 10x threshold.')
    else:
        A(u'**%d rows** fell below the threshold and have to be discussed at the meeting:' % len(dusen))
        A('')
        A(u'| target | criterion | in the panel | new | member coverage |')
        A('|---|---|---|---|---|')
        for r in dusen:
            A('| %s | %s | %s | %s x | %s |' % (
                r['hedef'], r['olcut'], r.get('PANEL_ayrim_sayi', ''),
                r.get('ayrim_en_kotu_x', ''), r.get('uye_kapsam', '')))
    A('')

    # ============================================================ 4. uyelik
    A(u'## 4. The membership definition: which number depends on the definition')
    A('')
    if not uyelik:
        A(u'*The membership audit was not run.*')
    else:
        oyn, tani = {}, {}
        for r in uyelik:
            h = r['hedef']
            try:
                o = float(r.get('oynaklik_x') or 0)
            except ValueError:
                o = 0
            oyn[h] = max(oyn.get(h, 0), o)
            if r.get('TANI') and r['TANI'] != 'SORUN YOK':
                tani[h] = r['TANI']
        riskli = sorted(((v, k) for k, v in oyn.items() if v and v > 1.5), reverse=True)
        A(u'- Targets whose discrimination ratio **moves by more than 1.5 fold** when the definition changes: **%d**'
          % len(riskli))
        if riskli:
            A('')
            A(u'| target | how many fold the discrimination moves when the definition changes |')
            A('|---|---|')
            for v, k in riskli[:10]:
                A('| %s | %.1fx |' % (k, v))
        if tani:
            A('')
            A(u'| target | diagnosis |')
            A('|---|---|')
            for k, v in tani.items():
                A('| %s | %s |' % (k, v))
        A('')
        A(u'Detail: `UYELIK_DENETIMI.md` . Where to correct it: `screening/hedef_uyelik.tsv`')
    A('')

    # ============================================================ 5. konsensus
    A('## 5. Konsensusler')
    A('')
    if not kons:
        A(u'*The consensus reproduction was not run.*')
    else:
        deg = [r for r in kons if (r.get('eski_ile_farkli') or '0') not in ('0', '')]
        kalite = konsensus_kalite()
        gecen_k = sum(1 for v in kalite.values() if v[0])
        A(u'- Bins reproduced: **%d**' % len(kons))
        A(u'- Bins holding a **different base** from the old consensus: **%d**' % len(deg))
        A(u'- Of a quality usable as a search backbone: **%d**' % gecen_k)
        if deg:
            A('')
            A(u'| bin | different bases | N %% | the two methods disagreed |')
            A('|---|---|---|---|')
            for r in deg[:12]:
                A('| %s | %s | %s | %s |' % (r['kutu'], r['eski_ile_farkli'],
                                             r['N_yuzde'], r['yontemler_ayrildi']))
        A('')
        A(u'Detail: `KONSENSUS_YENIDEN_URETIM.md`')
    A('')

    # ============================================================ 6. what is next
    A('## 6. Sirada ne var')
    A('')
    A(u'1. **Section 1**: on the targets where a better candidate was found, read the cost and decide.')
    A(u'2. **The pairs below the threshold in section 3**: they have to be discussed at the meeting.')
    A(u'3. **The unstable targets in section 4**: which membership definition they were reported under has to be written in the panel plainly.')
    A(u'4. Writing into the panel xlsx: `screening/update_panel.py` was written for that but **it is not in the menu**. It is run by hand on purpose, because it changes the delivery file and other sessions write to that file too.')
    A('')
    A(u'> This tool decides nothing by itself and changes no panel file by itself. It measures, it writes down the cost, and it leaves the choice to you.')
    A('')
    with open(yol, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(L))
    return yol


def yalniz_ozet(yaz, sure, cizgi):
    """Regenerate ONLY the summary report - it MAKES NO measurement and takes seconds.

        Separated so that a crash in the last step of a 5 hour run does not mean repeating
        the whole run. It reads the output files at hand and rewrites the summary.

    """
    import time as _t
    t0 = _t.time()
    cizgi('=')
    yaz(u'  REGENERATE THE SUMMARY REPORT ONLY')
    cizgi('=')
    yaz(u'  Nothing is measured; the existing output files are read.')
    yaz(u'  Source directory: %s' % C.CIKTI)
    yaz('')
    var = []
    for ad in ('panel_yeniden_olcum.tsv', 'panel_kutu_duzeyi.tsv',
               'uyelik_duyarlilik.tsv', 'konsensus_yeniden_uretim.tsv',
               'adaylar.tsv'):
        p_ = os.path.join(C.CIKTI, ad)
        v = os.path.exists(p_)
        var.append(v)
        yaz('   %-34s %s' % (ad, 'VAR' if v else 'YOK'))
    kont = len([f for f in os.listdir(C.KONTROL)
                if f.startswith('hedef_') and f.endswith('.json')]) \
        if os.path.isdir(C.KONTROL) else 0
    yaz(u'   %-34s %d files' % (u'kontrol/hedef_*.json (search)', kont))
    yaz('')
    if not any(var) and not kont:
        yaz(u'  No output file at all; run (9) first, or the stages one by one.')
        return 2
    d = durum_oku()
    yol = ozet_yaz(d, sure, _t.time() - t0)
    yaz(u'  Summary written (%s):' % sure(_t.time() - t0))
    yaz('    %s' % yol)
    cizgi('=')
    return 0
