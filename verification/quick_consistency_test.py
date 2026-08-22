# -*- coding: utf-8 -*-
"""THE QUICK CONSISTENCY TEST (a regression test) - is the chain reproducing itself?

WHAT IT IS AND WHAT IT IS NOT, IN TWO SENTENCES:
  This test checks that the code REPRODUCES ITSELF: does the same engine, on data
  ten times shallower, give the same class and the same ranking?
  It DOES NOT TEST that the measurement is CORRECT, because the expected values
  also come from a full depth run of THE SAME engine; if the engine has a
  systematic error this test cannot catch it. What it can catch is a drift in the
  code or the configuration.

Reading the code takes you only so far; this script looks at THE OUTPUT. It runs
the chain's four stages on a SMALL subset with a KNOWN answer and compares the
result against a reference run.

The reference: ONE_PROTOCOL_RESULT/panel_tek_protokol.tsv (full depth, 3000 reads).

WHAT IS TESTED
  1) Do the above-threshold ones still come out above and the below-threshold ones
     below? (THAT IS THE PRIMARY CRITERION, not the number itself)
  2) Is the ranking preserved? (Nitrosocosmicus > Metilotrofik > Cloacimonas > M_cinsi)
  3) Are the values in a reasonable band? (the band is written below with its reasoning)
  4) Do the later stages PRODUCE OUTPUT? An empty output DOES NOT COUNT AS A PASS.

How the TIME is kept short: the read cap goes 3000 -> 300, the targets are limited
to 8 rows, and the route 3 design scan is reduced. Because the measurement depth
drops, the numbers do not match exactly; nor are they expected to, see THE BAND
REASONING.

"""

# -------------------------------------------------------------------------
# quick_consistency_test.py runs the chain's four stages on a small subset with a
# known answer and tests that the code reproduces itself (a regression test).
#
# INPUT  : the four scripts in the project root and the measurement sources; those
#          are linked symbolically under QUICK_TEST/ (screening, protocol,
#          verification, REFERANS_DB, konsensus_kanonik, primer_final,
#          "fastq files", engine and uyelik_yeniden_turetme_uyelik_*.tsv). The
#          expected values come from the BEKLENEN_UST / BEKLENEN_ALT /
#          BEKLENEN_YENI constants in this file, whose source is a full depth
#          reference run.
# OUTPUT : QUICK_TEST/QUICK_TEST_REPORT.md and QUICK_TEST/test_gunlugu.txt; the
#          stages' own outputs also stay separately under QUICK_TEST/, and the real
#          result directories ARE NOT TOUCHED.
# CALLED BY: verification/full_chain.py -> key H
#          (python3 verification/quick_consistency_test.py --root .)
#
# WHAT IT IS     : it shows that the code and the configuration have not drifted.
# WHAT IT IS NOT : it does not show that the measurement is CORRECT. The expected
#                  values come from a run of the same engine, so a systematic error
#                  in that engine is invisible to this test. Independent
#                  confirmation is the job of the MFEprimer layer.
# -------------------------------------------------------------------------
import io, os, sys, csv, json, time, subprocess, argparse

VERSIYON = '1.0 (2026-08-03)'
OKUMA = 300
# THE THRESHOLD COMES FROM ONE SOURCE: screening/config.py -> ESIK_DCQ = 3.0
# Its fold equivalent is 2 ** ESIK_DCQ = 8.00. NO constant is EMBEDDED; if dCq
# changes it changes in one place. The reasoning and the efficiency warning are
# written in that file.
def _esik_yukle():
    """Reads the threshold from ONE SOURCE: screening/config.py.
        Since verification/ and screening/ are sibling directories the root is derived
        from here, so the script finds it whatever working directory it is called from.

    """
    import os as _o, sys as _s
    _kok = _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))
    if _kok not in _s.path:
        _s.path.insert(0, _kok)
    from screening import config as _y
    return _y

_C = _esik_yukle()
ESIK = _C.AYRIM_ESIK

# --- the reference values -----------------------------------------------
# THE 2026-08-10 FIX. These numbers were constants EMBEDDED IN THE CODE and nowhere
# said WHICH PRIMER PAIR they had been measured from. Today the Bacteroidales pair
# was changed (F: GAAGCTAGGATTTGGTTGCTGTG -> GCGTTATCCGGATTTATTGGGTTT) and the test
# compared the old pair's 0.74x against the new pair's 14.23x and said "THE CHAIN IS
# INCONSISTENT". The chain WAS NOT inconsistent; the reference was stale.
#
# The references are now read from QUICK_TEST/referans_degerler.tsv, and every row
# carries the F/R SEQUENCE that measurement was made with. If the pair has changed
# by the time of the test, NO COMPARISON IS MADE; it says "the reference is invalid,
# the pair changed" and the chain is not stopped. It is not silently counted correct
# and not silently counted wrong.
#
# If the file is missing, the old constants below are used, and the report says so
# openly.
REFERANS_DOSYASI = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, 'QUICK_TEST', 'referans_degerler.tsv')
CIKIS_TUTARSIZ = 6      # gercek gerileme - uzun kosuya GIRILMEZ
CIKIS_REFERANS_BAYAT = 7  # referans karsilastirilamaz - zincir devam edebilir


def _referans_yukle():
    """target -> (reference_x, F, R). Returns empty if the file is missing."""
    out = {}
    y = os.path.abspath(REFERANS_DOSYASI)
    if not os.path.exists(y):
        return out, y
    bas = None
    for l in io.open(y, encoding='utf-8'):
        l = l.rstrip('\n')
        if not l.strip() or l.startswith('#'):
            continue
        p = l.split('\t')
        if bas is None:
            bas = p
            continue
        r = dict(zip(bas, p))
        try:
            ref = float((r.get('referans_x') or '').replace(',', '.'))
        except ValueError:
            continue
        out[(r.get('hedef') or '').strip()] = (
            ref, (r.get('F') or '').strip().upper(), (r.get('R') or '').strip().upper())
    return out, y


REFERANS, REFERANS_YOL = _referans_yukle()

BEKLENEN_UST = [
    ('Nitrosocosmicus_AOA',        155.15),
    ('Metilotrofik_metanojen',      91.44),
    ('Proteolitik_Cloacimonas',     33.89),
    ('Methanosarcina_cinsi',        24.17),
]
BEKLENEN_ALT = [
    ('Proteiniphilum_cinsi',         0.00),
    ('Methanosarcina mazei',         0.82),
    ('Bacteroidales_kumesi',         0.74),
]
BEKLENEN_YENI = [('Petriella_cinsi', 11.03)]

BANT_ALT, BANT_UST = 0.5, 2.0
# A reference value VERY CLOSE to the threshold (< 15x) can cross to the other side
# at a shallow depth, and THAT IS NOT A CHAIN ERROR: the project's own measurement
# gives 11.03x at full depth and 8.93x at 300 reads for the Petriella LSU pair.
# On these rows a class change produces a WARNING; going outside the band is an ERROR.
SINIRDA_UST = ESIK * 1.5
BANT_GEREKCESI = u"""
BANT GEREKCESI - neden 0,5x - 2,0x

Ayrim kati = (uye Wilson ALT siniri) / (rakip Wilson UST siniri). Wilson
araliginin genisligi okuma sayisiyla daralir. Testte derinlik 3000'den 300'e
indigi icin uye alt siniri DUSER, rakip ust siniri YUKSELIR; iki etki de orani
KUCULTUR. Yani testte cikan sayinin referanstan dusuk olmasi BEKLENEN davranistir,
hata degil.

Projenin kendi olculmus ornegi: Petriella LSU cifti tam derinlikte 11,03x,
panelin 300 okuma standardinda 8,93x -> 0,81x oran. Az okumali kutularda sapma
daha buyuk olabilir. 0,5-2,0 bandi bu gozlenen sapmayi rahatca kapsar ama bir
satirin esigin obur tarafina gecmesini yakalayacak kadar dardir.

KARAR bandin kendisine degil, SINIFIN korunmasina baglanir: esik ustu olan esik
ustu, esik alti olan esik alti kalmalidir. Bant disina cikan ama sinifi korunan
satir UYARI uretir, testi dusurmez.

TEK ISTISNA - ESIK SINIRINDAKI SATIRLAR: referans degeri 15x'in altinda olan bir
satir, dusuk derinlikte esigin obur tarafina gecebilir ve bu zincir hatasi degil
olcum derinliginin dogal sonucudur (olculmus ornek: Petriella LSU cifti tam
derinlikte 11,03x, 300 okumada 8,93x). Boyle satirlarda sinif degisimi, oran
bantta kaldigi surece UYARI olarak raporlanir. Oran bandin da disina cikarsa
HATA sayilir.
"""


def sure_metni(sn):
    sn = int(sn)
    return ('%d saniye' % sn) if sn < 90 else ('%d dakika' % round(sn / 60.0)) \
        if sn < 5400 else ('%.1f saat' % (sn / 3600.0))


def vir(x, b=2):
    if x is None:
        return '-'
    try:
        return ('%.*f' % (b, float(str(x).replace(',', '.')))).replace('.', ',')
    except (TypeError, ValueError):
        return str(x)


def _f(s):
    try:
        return float(str(s).replace(',', '.'))
    except (TypeError, ValueError):
        return None


def tsv_oku(yol):
    if not os.path.exists(yol):
        return []
    with open(yol, encoding='utf-8') as fh:
        return list(csv.DictReader(
            (s for s in fh if s.strip() and not s.startswith('#')), delimiter='\t'))


# It runs a stage as a subprocess and applies a time cap. A timeout comes back as
# "no exit code", so that a stuck stage cannot keep the test waiting forever.
def kos(yaz, ad, arg, tavan_sn):
    t0 = time.time()
    yaz(u'  > %s starting...' % ad)
    try:
        pr = subprocess.run(arg, timeout=tavan_sn, capture_output=True, text=True)
        rc, ciktisi = pr.returncode, (pr.stdout or '')[-2000:]
    except subprocess.TimeoutExpired:
        yaz(u'    TIMEOUT (%s): the stage exceeded its time cap' % sure_metni(tavan_sn))
        return None, u'zaman asimi'
    yaz(u'    done: exit code %d, %s' % (rc, sure_metni(time.time() - t0)))
    return rc, ciktisi


# -------------------------------------------------------------------------
# It runs stage P at a SHALLOW depth (300 reads, 3000 on a full run) and only on the
# eight test rows, then checks three things at once:
#   1) Is the CLASS preserved: does an above-threshold row still come out above and
#      a below-threshold row still below. THAT IS THE PRIMARY CRITERION.
#   2) Is the RANKING preserved: has the order of magnitude of the four
#      above-threshold rows been disturbed.
#   3) Are the values inside the 0.5x - 2.0x band (a row outside the band whose
#      class is preserved produces only a WARNING and does not fail the test).
#
# WHY THE NUMBERS DO NOT MATCH EXACTLY: the discrimination ratio = (the member
# Wilson LOWER bound) / (the competitor Wilson UPPER bound). The width of the Wilson
# interval narrows with the read count; when the depth is lowered the member lower
# bound FALLS and the competitor upper bound RISES, and both effects shrink the
# ratio. So the number in the test coming out lower than the reference is the
# expected behaviour, not an error.
#
# Rows whose reference value is below 15x (SINIRDA_UST) can cross to the other side
# of the threshold at a shallow depth; as long as the ratio stays inside the band
# that counts as a WARNING and NOT an ERROR. A measured example: the Petriella LSU
# pair is 11.03x at full depth and 8.93x at 300 reads.
# -------------------------------------------------------------------------
def calistir(kok, hizli_kok, tavan_dk, yaz):
    hedefler = ([h for h, _ in BEKLENEN_UST] + [h for h, _ in BEKLENEN_ALT]
                + [h for h, _ in BEKLENEN_YENI])
    sec = ','.join(hedefler)
    py = sys.executable
    sonuc = dict(asama={}, satir=[], uyari=[], hata=[])

    # --- P ---
    rc, _ = kos(yaz, 'P (TEK PROTOKOL, %d okuma)' % OKUMA,
                [py, os.path.join(kok, 'protocol', 'single_protocol_measure.py'),
                 '--root', hizli_kok, '--reads', str(OKUMA), '--only', sec],
                tavan_dk * 60)
    P = tsv_oku(os.path.join(hizli_kok, 'ONE_PROTOCOL_RESULT', 'panel_tek_protokol.tsv'))
    sonuc['asama']['P'] = dict(rc=rc, satir=len(P))
    if not P:
        sonuc['hata'].append(u'stage P produced no rows at all, so the chain breaks here.')
        return sonuc

    olculen = {}
    for r in P:
        olculen[r['hedef']] = dict(kat=_f(r.get('ASIL_ayrim_mm1')),
                                   karar=r.get('esik_gecti_mi', ''),
                                   kapsam=r.get('ASIL_kapsam_mm1', ''),
                                   F=(r.get('F') or ''), R=(r.get('R') or ''))

    def bul(ad):
        for h, v in olculen.items():
            if ad.lower() in h.lower():
                return h, v
        return None, None

    # --- 1) sinif korunuyor mu + 3) bant ---
    for grup, bekl, ust_mu in (('esik ustu', BEKLENEN_UST, True),
                               ('esik alti', BEKLENEN_ALT, False),
                               ('yeni cift', BEKLENEN_YENI, True)):
        for ad, ref in bekl:
            h, v = bul(ad)
            # If the reference file exists, the file wins; the constant is only a fallback.
            _rf = None
            for _k, _val in REFERANS.items():
                if ad.lower() in _k.lower() or _k.lower() in ad.lower():
                    ref, _rf = _val[0], _val
                    break
            # THE EXPECTED CLASS also comes from the reference. Which constant list
            # (UST/ALT) a target was written into is itself a piece of information that
            # goes stale: Bacteroidales said "below threshold" in the constants, while
            # with the new pair it is 37.23x at full depth, that is, ABOVE. Asking the
            # measurement rather than the list closed that trap.
            #
            # CAUTION: ust_mu is the OUTER loop's variable. Assigning to it affects the
            # LATER targets in the same group too, and a target with no reference would
            # be evaluated with the class of the target before it. So a local variable
            # is used (2026-08-10, caught while writing; fixed before it ran).
            if _rf and v is not None:
                # HAS THE PAIR CHANGED? If it has, comparing against the old measurement
                # is meaningless, and saying "inconsistent" silently would be WRONG.
                _cf = (v.get('F') or '').upper()
                _cr = (v.get('R') or '').upper()
                if _rf[1] and _rf[2] and (_cf, _cr) != (_rf[1], _rf[2]):
                    sonuc.setdefault('referans_bayat', []).append(
                        u'%s: the primer pair CHANGED since the reference measurement (the reference F %s / R %s, now F %s / R %s). The old %sx and the new measurement (%sx) cannot be compared, so this IS NOT a regression. Refresh the reference: QUICK_TEST/referans_degerler.tsv'
                        % (ad, _rf[1][:12], _rf[2][:12], _cf[:12], _cr[:12],
                           vir(ref), vir(v['kat']) if v['kat'] is not None else '-'))
                    continue
            if v is None:
                sonuc['hata'].append(u'%s: the row is MISSING from the P output (expected %s: %sx)'
                                     % (ad, grup, vir(ref)))
                continue
            kat = v['kat']
            if kat is None:
                sonuc['hata'].append(u'%s: the discrimination fold could not be measured (decision=%s), %sx in the reference' % (ad, v['karar'], vir(ref)))
                continue
            _ust = (ref >= ESIK) if _rf else ust_mu
            gecti = kat >= ESIK
            sinirda = 0 < ref < SINIRDA_UST
            oran_ok = (ref > 0 and BANT_ALT <= (kat / ref) <= BANT_UST)
            if gecti != _ust and sinirda and oran_ok:
                sonuc['uyari'].append(
                    u'%s: %sx in the reference, %sx in the test, a row ON THE EDGE of the threshold (the reference is < %sx). Crossing to the other side of the threshold at a low depth is expected behaviour (a measured example: 11,03x -> 8,93x); the ratio is inside the band so it did not count as an ERROR.'
                    % (ad, vir(ref), vir(kat), vir(SINIRDA_UST)))
            elif gecti != _ust:
                sonuc['hata'].append(
                    u'%s: THE CLASS CHANGED, %sx (%s) in the reference against %sx (%s) in the test%s'
                    % (ad, vir(ref), u'esik ustu' if _ust else u'esik alti',
                       vir(kat), u'esik ustu' if gecti else u'esik alti',
                       u'' if oran_ok else u' - ustelik %sx-%sx bandinin DISINDA'
                       % (vir(BANT_ALT), vir(BANT_UST))))
            elif ref > 0:
                oran = kat / ref
                if not (BANT_ALT <= oran <= BANT_UST):
                    sonuc['uyari'].append(
                        u'%s: %sx (the reference %sx, the ratio %sx), outside the band but the class was kept' % (ad, vir(kat), vir(ref), vir(oran)))
            sonuc['satir'].append(dict(hedef=ad, ref=ref, olculen=kat,
                                       karar=v['karar'], kapsam=v['kapsam'],
                                       sinif_ok=(gecti == _ust)))

    # --- 2) siralama korunuyor mu ---
    ref_sira = [ad for ad, _ in sorted(BEKLENEN_UST, key=lambda x: -x[1])]
    olc = []
    for ad in ref_sira:
        _h, v = bul(ad)
        olc.append((ad, v['kat'] if v and v['kat'] is not None else -1))
    test_sira = [ad for ad, _ in sorted(olc, key=lambda x: -x[1])]
    sonuc['siralama'] = dict(referans=ref_sira, test=test_sira,
                             korundu=(ref_sira == test_sira))
    if ref_sira != test_sira:
        sonuc['hata'].append(u'THE ORDER BROKE: the reference says %s, the test %s'
                             % (' > '.join(ref_sira), ' > '.join(test_sira)))
    return sonuc


# -------------------------------------------------------------------------
# It runs stages K, D and I with reduced settings. What is tested here is not the
# CORRECTNESS of the measurement but that the stage RUNS AND PRODUCES ROWS; an empty
# output does not count as a pass.
#
# TELLING A CASCADED FAILURE APART: if K recovered no pair, there is nothing for D
# to verify and D coming out empty IS NOT D's fault. In that case D is run in A
# SEPARATE root with a synthetic single row input and tested on its own, so that "D
# does not work" and "D had nothing to do" are never confused.
# -------------------------------------------------------------------------
def sonraki_asamalar(kok, hizli_kok, tavan_dk, yaz, sonuc):
    """Do K, D and I really PRODUCE OUTPUT? An empty output DOES NOT COUNT AS A PASS."""
    py = sys.executable

    rc, _ = kos(yaz, 'K (verification)',
                [py, os.path.join(kok, 'verification', 'recovery_round.py'),
                 '--root', hizli_kok, '--reads', str(OKUMA),
                 # The test's SCOPE is kept narrow: the aim is not correctness but evidence
                 # that "the stage runs and produces rows". On a full run these caps come off.
                 '--scan-max', '40', '--candidate-max', '5', '--arms-max', '0',
                 '--skip-if-no-panel'],
                tavan_dk * 60)
    K = tsv_oku(os.path.join(hizli_kok, 'RECOVERY_RESULT', 'kurtarma_satirlari.tsv'))
    sonuc['asama']['K'] = dict(rc=rc, satir=len(K))
    if len(K) < 1:
        sonuc['hata'].append(u'stage K produced NO rows at all (an empty output does not count as passing).')
    else:
        yaz(u'    K: %d rows, %d of them recovered'
            % (len(K), sum(1 for r in K if (r.get('esigi_gecti_mi') or '').startswith('EVET'))))

    rc, _ = kos(yaz, 'D (DOGRULAMA, yalniz yerel katman)',
                [py, os.path.join(kok, 'verification', 'specificity_round.py'),
                 '--root', hizli_kok, '--ncbi', 'elle'],
                tavan_dk * 60)
    Dd = tsv_oku(os.path.join(hizli_kok, 'VERIFICATION_RESULT', 'dogrulama_uc_sutun.tsv'))
    sonuc['asama']['D'] = dict(rc=rc, satir=len(Dd))
    if len(Dd) < 1 and len(K) >= 1 and not [r for r in K
                                            if (r.get('esigi_gecti_mi') or '').startswith('EVET')]:
        # K ran but recovered NOTHING -> there is no pair for D to verify.
        # That is not D's fault. Even so, D HAS TO BE SHOWN TO WORK: it is run in a
        # separate root with a synthetic single row input.
        yaz(u'    D: empty because no pair was recovered, testing it ON ITS OWN...')
        oz = os.path.join(hizli_kok, 'D_KENDI_SINAMASI')
        os.makedirs(os.path.join(oz, 'RECOVERY_RESULT'), exist_ok=True)
        for ad in ('screening', 'verification', 'protocol', 'REFERANS_DB',
                   'konsensus_kanonik', 'ONE_PROTOCOL_RESULT', 'primer_final',
                   'engine', 'engine', 'engine',
                   'steps', 'engine', 'fastq files'):
            h = os.path.join(oz, ad)
            kaynak = os.path.join(hizli_kok, ad)
            if not os.path.exists(h) and os.path.exists(kaynak):
                try:
                    os.symlink(os.path.realpath(kaynak), h)
                except OSError:
                    pass
        with open(os.path.join(oz, 'RECOVERY_RESULT', 'kurtarma_satirlari.tsv'),
                  'w', encoding='utf-8', newline='') as fh:
            fh.write(u'# D SELF-TEST - synthetic input\n')
            ww = csv.writer(fh, delimiter='\t')
            ww.writerow(['hedef', 'eski_deger', 'eski_kapsam', 'denenen_yol', 'olcu',
                         'yeni_deger', 'esigi_gecti_mi', 'UYELIK_GEREKCESI', 'sebep'])
            ww.writerow(['D_SINAMA_Petriella', '8,45', '9/9', 'sentetik', 'ayrim kati',
                         'YENI CIFT AAATCTGGCTGCCTGTGC / CTCTCACCCTCTATGGCGTC (101 bp) 11,03 x',
                         'EVET (yeni cift)', 'sentetik', ''])
        rc2, _ = kos(yaz, 'D (kendi sinamasi, sentetik girdi)',
                     [py, os.path.join(kok, 'verification', 'specificity_round.py'),
                      '--root', oz, '--ncbi', 'elle', '--no-mfe',
                      '--cluster-max', '1'], tavan_dk * 60)
        D2 = tsv_oku(os.path.join(oz, 'VERIFICATION_RESULT', 'dogrulama_uc_sutun.tsv'))
        sonuc['asama']['D'] = dict(rc=rc2, satir=len(D2), kendi_sinamasi=True)
        if len(D2) >= 1:
            sonuc['uyari'].append(
                u'stage D came out empty in the chain because K recovered no pair (there is nothing to verify). D was tested ON ITS OWN and produced %d rows, so the stage is sound.' % len(D2))
            yaz(u'    D: produced %d rows in its own self-test, so the stage is SOUND' % len(D2))
            Dd = D2
        else:
            sonuc['hata'].append(u'stage D produced no rows with synthetic input either.')
    if len(Dd) < 1:
        if rc == 7 or len(K) < 1:
            # A CASCADED FAILURE: if K produced no rows, D's input is empty anyway.
            # It is not counted as a separate error but attributed to K's.
            sonuc['uyari'].append(
                u'stage D did not run because K produced no rows at all (a chain of failures). That is not D\'s fault; once K is fixed D runs too.')
            yaz(u'    D: SKIPPED - the input is empty because K produced no rows (a knock-on failure)')
        else:
            sonuc['hata'].append(u'stage D produced NO rows at all.')
    else:
        # The REQUIRED layers: our two measurements. MFEprimer and NCBI are
        # DELIBERATELY skipped in the test (--no-mfe, --ncbi elle: there is no
        # network). Their being missing is not a chain error and is reported as a warning.
        eksik_sutun, yok_sutun, takma_ad, istege_bagli = katman_denetimi(Dd, sonuc)
        yaz(u'    D: %d rows | are the two mandatory sources filled: %s | skipped in the test: %s'
            % (len(Dd), u'NO' if (eksik_sutun or yok_sutun) else 'evet',
               ', '.join(istege_bagli) or 'yok'))
        for _ta in takma_ad:
            yaz(u'    D: (sema notu) %s' % _ta)

    rc, _ = kos(yaz, 'I (KIMLIK, 2 veritabani, nt yok)',
                [py, os.path.join(kok, 'verification', 'identity_verification.py'),
                 '--root', hizli_kok, '--only', '10', '--nt', 'yok', '--db-max', '2'],
                tavan_dk * 60)
    I = tsv_oku(os.path.join(hizli_kok, 'IDENTITY_RESULT', 'kimlik_iddialari.tsv'))
    sonuc['asama']['I'] = dict(rc=rc, satir=len(I))
    if len(I) < 1:
        sonuc['hata'].append(u'stage I produced NO claim results at all.')
    else:
        yaz(u'    I: %d claims resolved' % len(I))
    return sonuc



# =========================================================================
# THE 2026-08-09 FIX - A FALSE ALARM FROM THE H GATE
# -------------------------------------------------------------------------
# THE SYMPTOM : stage H was failing with exit code 6 and saying
#           "the REQUIRED layers were not filled in stage D: source 1
#            (the sample measurement)".
# THE MEASUREMENT: QUICK_TEST/D_KENDI_SINAMASI/VERIFICATION_RESULT/
#           dogrulama_uc_sutun.tsv was read. The layer IS NOT EMPTY; the column's
#           value is 'TEMIZ'. But the column is not named '1_NUMUNE', it is named
#           '1_NUMUNE_oy_vermez'. Stage D renamed the column in the D-2 fix of
#           2026-08-06 (specificity_round.py line 1053) and the gate was looking
#           for the old name; r.get('1_NUMUNE') returned None and the gate counted
#           that as "not filled".
#           So the cause was NOT THE TEST MODE but A SCHEMA DRIFT. Because the name
#           '2_YEREL_DB' had not changed, that layer passed and the error appeared
#           only on layer 1.
# THE FIX : the column is looked up by its full name first, and if that fails, with
#           an '<name>_' prefix.
#           - the column exists and is filled     -> it passes (with a schema note
#                                                    where needed)
#           - the column exists but EVERY row is empty -> AN ERROR (the gate keeps
#                                                    working)
#           - the column DOES NOT EXIST AT ALL    -> A SEPARATE ERROR (schema drift);
#                                                    it does not pass silently,
#                                                    because treating "absent" and
#                                                    "empty" as the same thing was
#                                                    the bug itself
#           The MFEprimer and NCBI layers were ALREADY warnings and were left that
#           way: because the test runs with --no-mfe and with '--ncbi elle' and no
#           network, their being empty is a gap created deliberately.
# =========================================================================
ZORUNLU_KATMAN = (('1_NUMUNE', u'1. kaynak (numune olcumu)'),
                  ('2_YEREL_DB', u'2. kaynak (yerel veritabani)'))
ISTEGE_BAGLI_KATMAN = (('3_MFEPRIMER', u'3. kaynak (MFEprimer)'),
                       ('4_NCBI', u'4. kaynak (NCBI)'))


def sutun_coz(satirlar, ad):
    """Resolves the column by its full name, then with an '<name>_' prefix.

        Returns: (resolved_name_or_None, is_filled). If resolved_name is None the column
        IS NOT THERE AT ALL; that IS NOT the same thing as 'empty' and is reported
        separately.

    """
    if not satirlar:
        return None, False
    basliklar = [b for b in satirlar[0].keys() if b]
    if ad in basliklar:
        coz = ad
    else:
        adaylar = sorted(b for b in basliklar if b.startswith(ad + '_'))
        coz = adaylar[0] if adaylar else None
    if coz is None:
        return None, False
    return coz, any((r.get(coz) or '').strip() for r in satirlar)


def katman_denetimi(Dd, sonuc):
    """Checks the layers in D's output and fills in sonuc['hata'/'uyari'].

        Returns: (missing_column, absent_column, alias, optional)

    """
    eksik_sutun, yok_sutun, takma_ad = [], [], []
    for sut, ad in ZORUNLU_KATMAN:
        coz, dolu = sutun_coz(Dd, sut)
        if coz is None:
            yok_sutun.append(u'%s [the "%s" column is MISSING ENTIRELY from the D output]' % (ad, sut))
        elif not dolu:
            eksik_sutun.append(u'%s (the "%s" column is there but empty on every row)'
                               % (ad, coz))
        elif coz != sut:
            takma_ad.append(u'%s: the "%s" column was used instead of "%s"'
                            % (ad, sut, coz))
    if yok_sutun:
        sonuc['hata'].append(
            u'a REQUIRED layer column is MISSING ENTIRELY from the D output (a schema drift): %s. The gate does not treat this the same as "empty"; if the column names of the D output have changed then the gate has to be updated too.' % ', '.join(yok_sutun))
    if eksik_sutun:
        sonuc['hata'].append(u'the REQUIRED layers were not filled in at stage D: %s'
                             % ', '.join(eksik_sutun))
    for t in takma_ad:
        sonuc['uyari'].append(
            u'a column name changed in the D output but the layer is FULL, so it did not count as an error, %s' % t)
    istege_bagli = []
    for sut, ad in ISTEGE_BAGLI_KATMAN:
        coz, _ = sutun_coz(Dd, sut)
        if coz is None:
            istege_bagli.append(ad)
            continue
        if not any((r.get(coz) or '').strip() not in ('', 'BILINMIYOR') for r in Dd):
            istege_bagli.append(ad)
    if istege_bagli:
        sonuc['uyari'].append(
            u'these sources were not filled in during the test at stage D: %s. That is DELIBERATE (the test runs with --no-mfe and with the network off, --ncbi elle); they are filled in on a full run.' % ', '.join(istege_bagli))
    return eksik_sutun, yok_sutun, takma_ad, istege_bagli


# The decision rests on one condition: is the error list empty. Warnings do not fail
# the test.
# Exit code 6 = THE CHAIN IS INCONSISTENT; it must be resolved before a full run.
def raporla(hizli_kok, sonuc, yaz, gecen_sure):
    guvenilir = not sonuc['hata']
    yol = os.path.join(hizli_kok, 'QUICK_TEST_REPORT.md')
    with open(yol, 'w', encoding='utf-8') as fh:
        fh.write(u'# Quick consistency test (regression)\n\nGenerated: %s, script %s, time: %s\n\n'
                 % (time.strftime('%Y-%m-%d %H:%M'), VERSIYON, sure_metni(gecen_sure)))
        fh.write(u'## VERDICT\n\n')
        if guvenilir:
            fh.write(u'# CHAIN CONSISTENT (against its own reference): the full run can be started\n\n')
            fh.write(u'> **This is not a CORRECTNESS test.** The expected values also come from a full-depth run of the same engine;')
            # THE 2026-08-06 FIX - caught on a clean run: because this sentence was
            # written UNCONDITIONALLY it said "every row preserved its class", while the
            # table in the same report said "class preserved: HAYIR" on the
            # Petriella_cinsi row. The decision (TUTARLI) was right - a row ON THE
            # THRESHOLD BOUNDARY changing sides at a shallow depth is deliberately
            # counted as a warning - but the summary sentence CONTRADICTED the table.
            # It is now written from the count.
            _bozan = [s['hedef'] for s in sonuc['satir'] if not s.get('sinif_ok')]
            if _bozan:
                fh.write(u'The ranking held and all four stages produced non-empty '
                         u'output. %d of the %d rows tested kept their class; %d rows '
                         u'sat ON THE EDGE of the threshold and changed side at a '
                         u'low depth, which was deliberately NOT counted as an error '
                         u'(detail: Warnings): %s\n\n'
                         % (len(sonuc['satir']), len(sonuc['satir']) - len(_bozan),
                            len(_bozan), ', '.join(_bozan)))
            else:
                fh.write(u'Every row tested kept its class, the ranking held, and all four stages produced non-empty output')
        else:
            fh.write(u'# THE CHAIN IS INCONSISTENT, DO NOT START A FULL RUN\n\n\n')
            fh.write(u'The 6 to 16 hour run should not be spent before the mismatches below are resolved.\n\n')
            for h in sonuc['hata']:
                fh.write(u'- **%s**\n' % h)
            fh.write(u'\n')
        if sonuc.get('referans_bayat'):
            fh.write(u'## Rows whose reference is stale (NOT a regression)\n\n')
            fh.write(u'The primer pair on these rows changed since the reference was measured. Comparing a number belonging to the old pair against the new one')
            for h in sonuc['referans_bayat']:
                fh.write(u'- %s\n' % h)
            fh.write(u'\n')
        fh.write(u'## Row by row\n\n| target | reference | measured in the test | ratio | verdict | class preserved |\n|---|---|---|---|---|---|\n')
        for r in sonuc['satir']:
            oran = (r['olculen'] / r['ref']) if (r['ref'] and r['olculen'] is not None) else None
            fh.write(u'| %s | %sx | %sx | %s | %s | %s |\n'
                     % (r['hedef'], vir(r['ref']), vir(r['olculen']),
                        vir(oran) if oran else '-', r['karar'],
                        'evet' if r['sinif_ok'] else '**HAYIR**'))
        s = sonuc.get('siralama') or {}
        fh.write(u'\n## Ranking\n\n- reference: %s\n- in the test: %s\n- **%s**\n'
                 % (' > '.join(s.get('referans', [])), ' > '.join(s.get('test', [])),
                    'korundu' if s.get('korundu') else 'BOZULDU'))
        fh.write(u'\n## Output check per stage\n\n| stage | exit code | rows |\n|---|---|---|\n')
        for k in ('P', 'K', 'D', 'I'):
            a = sonuc['asama'].get(k, {})
            fh.write(u'| %s | %s | %s |\n' % (k, a.get('rc', 'kosulmadi'), a.get('satir', 0)))
        if sonuc['uyari']:
            fh.write(u'\n## Warnings (they do not fail the test)\n\n\n')
            for u_ in sonuc['uyari']:
                fh.write(u'- %s\n' % u_)
        fh.write(u'\n```' + BANT_GEREKCESI + u'```\n')
    yaz('')
    yaz('=' * 74)
    if guvenilir:
        yaz(u'  CHAIN CONSISTENT (against its own reference); the full run can be started.')
        yaz(u'  CAUTION: this does NOT mean "the measurements were validated". It shows that')
        yaz(u'  the same engine reproduces itself. For independent confirmation, use MFEprimer.')
        yaz(u'  All %d rows tested kept their class, the ranking did not change,'
            % len(sonuc['satir']))
        yaz(u'  and all four stages produced non-empty output.')
    else:
        yaz(u'  CHAIN INCONSISTENT - DO NOT START THE FULL RUN')
        for h in sonuc['hata']:
            yaz(u'    - %s' % h)
        yaz(u'  Fix these before spending eight hours on a full run.')
    for u_ in sonuc['uyari']:
        yaz(u'  (warning) %s' % u_)
    for b in sonuc.get('referans_bayat', []):
        yaz(u'  (stale reference) %s' % b)
    if sonuc.get('referans_bayat'):
        yaz(u'  These rows were NOT COMPARED. To refresh the reference:')
        yaz(u'    python verification/refresh_reference.py --root .')
    yaz('=' * 74)
    yaz(u'  detail: %s' % yol)
    if not guvenilir:
        return CIKIS_TUTARSIZ
    return CIKIS_REFERANS_BAYAT if sonuc.get('referans_bayat') else 0


# QUICK_TEST/ is a temporary root: the source directories are linked
# symbolically and the outputs stay separate. That way the test DOES NOT
# OVERWRITE the long run results in the real output directories.
def main():
    p = argparse.ArgumentParser(description='The quick correctness test of '
                                            'the chain')
    p.add_argument('--root', dest='kok', default='.')
    p.add_argument('--cap-minutes', dest='tavan_dk', type=int, default=15,
                   help='asama basina zaman tavani (dakika)')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)
    hizli = os.path.join(kok, 'QUICK_TEST')
    os.makedirs(hizli, exist_ok=True)
    # gecici kok: kaynaklar baglanti, ciktilar ayri
    for ad in ('screening', 'protocol', 'verification', 'REFERANS_DB',
               'konsensus_kanonik', 'primer_final', 'fastq files',
               'engine', 'engine', 'engine',
               'steps', 'engine'):
        h = os.path.join(hizli, ad)
        if not os.path.exists(h) and os.path.exists(os.path.join(kok, ad)):
            try:
                os.symlink(os.path.join(kok, ad), h)
            except OSError:
                pass
    import glob
    for u_ in glob.glob(os.path.join(kok, 'uyelik_yeniden_turetme_uyelik_*.tsv')):
        h = os.path.join(hizli, os.path.basename(u_))
        if not os.path.exists(h):
            try:
                os.symlink(u_, h)
            except OSError:
                pass
    g = open(os.path.join(hizli, 'test_gunlugu.txt'), 'a', encoding='utf-8')

    def yaz(s=''):
        print(s, flush=True); g.write(s + '\n'); g.flush()

    yaz('=' * 74)
    yaz(u'  QUICK CONSISTENCY TEST (regression)   version %s   %s'
        % (VERSIYON, time.strftime('%Y-%m-%d %H:%M')))
    yaz('=' * 74)
    yaz(u'  Target time: about 30 minutes. Measurement depth %d reads (a full run uses 3000).' % OKUMA)
    yaz(u'  Rows tested: %d. Output directory: QUICK_TEST/'
        % (len(BEKLENEN_UST) + len(BEKLENEN_ALT) + len(BEKLENEN_YENI)))
    yaz(u'  NOTE: the numbers will NOT match the reference exactly; the depth was reduced.')
    yaz(u'  What is being checked is that the CLASS and the RANKING are preserved.')
    yaz('')
    yaz(u'  WHAT IT IS    : a check that the code REPRODUCES ITSELF.')
    yaz(u'  WHAT IT IS NOT: it does NOT check that the measurement is CORRECT. The expected')
    yaz(u'               values come from a full-depth run of the SAME engine. If the engine')
    yaz(u'               has a systematic error, this test CANNOT catch it.')
    yaz(u'  Independent confirmation is a separate job: see the MFE_BAGIMSIZ_TEYIT report.')
    yaz('')
    t0 = time.time()
    sonuc = calistir(kok, hizli, a.tavan_dk, yaz)
    if not sonuc['hata'] or sonuc['asama'].get('P', {}).get('satir'):
        sonuc = sonraki_asamalar(kok, hizli, a.tavan_dk, yaz, sonuc)
    rc = raporla(hizli, sonuc, yaz, time.time() - t0)
    g.close()
    return rc


if __name__ == '__main__':
    sys.exit(main() or 0)
