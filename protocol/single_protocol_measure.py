# -*- coding: utf-8 -*-
"""Single-protocol measurement, one rule, one depth, for the whole panel.

Every pair in the panel is measured under IDENTICAL settings, so the numbers
are comparable across targets. Mixed protocols were the earlier failure mode:
two pairs measured at different depths cannot be ranked against each other, and
nothing in the output revealed it.

Reports dCq (discrimination power) per pair, and the panel threshold each pair
must clear.

"""

# -------------------------------------------------------------------------
# single_protocol_measure.py re-measures EVERY pair in the panel under one rule and
# at one depth, and produces a single order list. There is no per-row exception.
#
# INPUT  : the panel table under primer_final/ (through
#          screening.targets.panel_oku), protocol/ek_ciftler.tsv (pairs that are
#          not in the panel), uyelik_yeniden_turetme_uyelik_*.tsv (stage U's
#          MEASURED membership), the raw reads under "fastq files".
# OUTPUT : TEK_PROTOKOL_SONUC/panel_tek_protokol.tsv (the full table),
#          TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv (the order decision),
#          TEK_PROTOKOL_SONUC/kutu_bazli_ham_sayilar.tsv (k and n; every verdict
#          can be recomputed from those two columns),
#          TEK_PROTOKOL_SONUC/PROTOKOL_VE_RAPOR.md, kontrol/ .
# CALLED BY: verification/full_chain.py -> key P
#          (python3 protocol/single_protocol_measure.py --root .)
#
# WHY IT EXISTS: in the old panel the rows had been measured under different
# conditions, some at mm<=1 and some at mm<=3, some at a depth of 300 reads and
# some at 46 000. Because the width of the Wilson interval depends on the read
# count, THE SAME real specificity gives a LOWER "x" in a shallow pool. The numbers
# in that column could not be compared with one another and the 10x threshold was
# not measuring the same thing from row to row. This script ends that confusion.
# -------------------------------------------------------------------------
import os, sys, csv, json, time, argparse, math

VERSIYON = '1.0 (2026-08-03)'

# --- THE THRESHOLD FROM ONE SOURCE: screening/config.py -> ESIK_DCQ ---
def _esik_yukle():
    import os as _o, sys as _s
    _kok = _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))
    if _kok not in _s.path:
        _s.path.insert(0, _kok)
    from screening import config as _y
    return _y

_C = _esik_yukle()


def _sinif_yukle():
    from screening import order_classes as _s
    return _s

_S = _sinif_yukle()


# --------------------------------------------------------------- protocol
# --- WHERE THE THRESHOLDS COME FROM ----------------------------------
# 2026-08-06: the threshold was fixed in terms of dCq. dCq >= 3 -> 2**3 = 8.00 fold.
# THAT IS NO LONGER A TOOL THRESHOLD: dCq >= 3 is the accepted floor in the
# literature for a specificity or NTC passing criterion. The earlier 10x really was
# a tool threshold (it first appeared as a code constant) and corresponded to a dCq
# of 3.32.
# The single source: screening/config.py -> ESIK_DCQ.
# The meeting's own criterion is still A DIFFERENT quantity and is reported in a
# SEPARATE column:
#   WORK_RECORD section 1.7 - "a tolerance of 1-2 CROSS REACTING SPECIES; the
#                         measure is the NUMBER OF cross reacting species, not the
#                         number of products formed in them"
#   WORK_RECORD section 1.5 - "no product may form in any competitor" (zero tolerance)
# Neither contains the other, nor is contained by 10x. So both are reported as
# SEPARATE COLUMNS, with a note of who set which.
ESIK_KOKENI = _C.ESIK_KOKENI
ESIK_VERIM_NOTU = _C.ESIK_VERIM_NOTU
# The MIQE and laboratory language: the discrimination ratio -> dCq. At 100%
# efficiency every cycle doubles, so dCq = log2(fold). 10x = 3.32 cycles. In the
# literature the specificity or NTC passing criterion is dCq >= 3 (NEB, high
# efficiency qPCR data analysis), and our 10x threshold falls just above it.
# Efficiency is assumed to be 100%; once the real efficiency is measured this should
# be corrected with dCq = log(fold)/log(1+E).
def dcq(kat, verim=1.0):
    import math
    try:
        k = float(kat)
    except (TypeError, ValueError):
        return None
    if k <= 0:
        return None
    return round(math.log(k) / math.log(1.0 + verim), 2)

TOPLANTI_CAPRAZ_TABAN = 10.0   # a competitor bin counts as "cross reacting" at >=10% product
TOPLANTI_CAPRAZ_HOSGORU = 2    # WORK_RECORD 1.7: a tolerance of 1 to 2 cross reacting species

PROTOKOL = dict(
    olcut_asil=1,
    olcut_yan=3,
    okuma_tavani=3000,
    esik=_C.AYRIM_ESIK,        # dCq 3 -> 8,00 kat (tek kaynak)
    esik_dcq=_C.ESIK_DCQ,
    karisik='rakip',
    kapsam_esigi=0.20,
    enkotu_asgari_okuma=150,
    urun_alt=60, urun_ust=400,
)

GEREKCE = u"""
PROTOKOL VE NEDEN BOYLE SECILDI
===============================

1) DERINLIK: kutu basina EN COK %(okuma_tavani)d okuma, sabit tohum, 200-6000 bp suzgeci.
   Neden tavan var: Wilson araliginin genisligi okuma sayisina baglidir. Tavan
   konmazsa 46 000 okumalik bir kutu ile 300 okumalik bir kutu ayni tabloda yan
   yana gelir ve derin kutunun ayrimi YAPAY olarak yuksek cikar. Tavan butun
   satirlari ayni istatistiksel zemine oturtur.
   Neden 3000: kutularin buyuk cogunlugu zaten bu sayinin altinda, yani veri
   kaybi kucuk; ustelik uyeligin turetildigi olcum de bu tavanla yapildi, boylece
   uyelik ile ayrim ayni zeminde kalir.

2) ASIL OLCUT: <=1 uyumsuzluk + 3' son 2 baz TAM eslesme.
   Panelin tasarim olcutu budur; primerler bu varsayimla secildi. Karar bu
   sutuna gore verilir.

3) YAN OLCUT: <=3 uyumsuzluk + 3' son 2 baz TAM eslesme.
   Karar sutunu DEGILDIR, DAYANIKLILIK gostergesidir. Gercek PCR'de gevsek
   baglanma olur; bir cift mm<=1'de gecip mm<=3'te cokuyorsa o cift kirilgandir
   ve bu gorunmelidir (olculmus ornek: ileri primeri NL1 olan cift 8,47x -> 0,67x).

4) KARISIK KUTULAR TEK KURALA BAGLI: karisik kutu = %(karisik)s sayilir.
   Karisik kutu, hedef organizmayi KISMEN tasidigi olculen kutudur. Uye saymak
   ayrimi yapay olarak yukseltir; tamamen dislamak gercek capraz sinyali gizler.
   RAKIP saymak en kotu durumu olcer - siparis karari icin dogru taraf budur.
   Degistirilebilir (--mixed uye|rakip|disla) ama secim her ciktinin basina yazilir.

5) ESIK: dCq %(esik_dcq).1f = %(esik).2f kat, EN KOTU TEK RAKIP KUTU uzerinden
   (asgari %(enkotu_asgari_okuma)d okuma). Literatur olcutu; arac esigi DEGILDIR.
   VERIM %%100 VARSAYILDI - %%90 verimde ayni dCq 6,86 kat eder.
   Havuz kati da raporlanir ama karar en kotu kutuya gore verilir: havuz kati
   tek bir kotu kutuyu binlerce temiz okumanin icinde eritir.

6) KAPSAM ayri eksendir: bir uye kutu >=%%%(kapsam_yuzde)d urun veriyorsa "kapsandi" sayilir.
   Ayrim yuksek ama kapsam dusukse cift ozguldur fakat hedefin tamamini gormez;
   bu iki sorun birbirine karistirilmamalidir.

7) UYELIK: Kraken etiketinden DEGIL, olculen kimlikten gelir
   (uyelik_yeniden_turetme_uyelik_*.tsv). Yanlis etiketli bir kutu rakip hanesine
   yazilinca metrik hedefi hedefle kiyaslar ve mukemmel bir primer bile 1'in
   altinda cikar - olculmus ornek: 0,71x -> 8,47x, ayni primer.
""" % dict(PROTOKOL, kapsam_yuzde=int(PROTOKOL['kapsam_esigi'] * 100))


def sure_metni(sn):
    sn = int(sn)
    if sn < 90:
        return '%d saniye' % sn
    if sn < 5400:
        return '%d dakika' % round(sn / 60.0)
    return '%.1f saat' % (sn / 3600.0)


def vir(x, basamak=2):
    """Turkce ondalik: 12.96 -> '12,96'"""
    if x is None or x == '':
        return '-'
    return ('%.*f' % (basamak, x)).replace('.', ',')


# --------------------------------------------------------------- inputs
# It verifies the project root. Without the screening directory the measurement
# modules cannot be imported; failing early and clearly beats a half finished run.
def kok_bul(arg):
    kok = os.path.abspath(arg or '.')
    if not os.path.isdir(os.path.join(kok, 'screening')):
        sys.exit(u'ERROR: there is no screening directory inside %s. Give the project directory with --root.' % kok)
    return kok


def uyelik_dosyasi(kok):
    """Finds the newest membership file - BY TIME, not BY NAME.

        The 2026-08-10 fix. The old code concatenated two globs and took a[-1]; that
        DID NOT MEAN "the newest". The ordering was alphabetical and engine_SONUC
        entries always came AFTER the ones in the root. So a file dated 1 August in a
        subdirectory would beat one dated 3 August in the root. There is only one
        candidate at the moment, so the behaviour does not change; but if the next run
        produced a second file, the wrong membership would be picked silently. A wrong
        membership makes the discrimination ratio come out smaller or larger than it is
        (a measured example: the same pair moved between 0.71x and 8.47x).

    """
    import glob
    a = glob.glob(os.path.join(kok, 'uyelik_yeniden_turetme_uyelik_*.tsv'))
    a += glob.glob(os.path.join(kok, 'engine_SONUC', '*uyelik*.tsv'))
    if not a:
        return None
    # time first, then the name on a tie - both criteria written OPENLY
    a.sort(key=lambda p: (os.path.getmtime(p), os.path.basename(p)))
    return a[-1]


def uyelik_oku(yol):
    """target -> dict(uye=[...], karisik=[...], rakip=[...], sinif=...)"""
    out = {}
    with open(yol, encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            bol = lambda s: [x for x in (r.get(s) or '').split(';') if x.strip()]
            out[r['hedef'].strip()] = dict(
                sinif=(r.get('sinif') or '').strip(),
                uye=bol('yeni_uye_kutular') or bol('eski_uye_kutular'),
                karisik=bol('karisik_kutular'),
                rakip=bol('rakip_kutular'))
    return out


def ek_ciftler_oku(kok):
    """protocol/ek_ciftler.tsv - pairs that are NOT IN the panel TSV.
        The user can edit it by hand. The columns: hedef, sinif, F, R, urun_bp, not

    """
    yol = os.path.join(kok, 'protocol', 'ek_ciftler.tsv')
    if not os.path.exists(yol):
        return []
    out = []
    with open(yol, encoding='utf-8') as fh:
        for r in csv.DictReader((s for s in fh if not s.startswith('#')), delimiter='\t'):
            if not (r.get('hedef') or '').strip():
                continue
            out.append(dict(hedef=r['hedef'].strip(), sinif=(r.get('sinif') or '').strip(),
                            F=r['F'].strip().upper(), R=r['R'].strip().upper(),
                            urun_bp=int(float(r.get('urun_bp') or 0)),
                            kaynak='EK', ta=(r.get('ta') or '').strip(),
                            uyelik_hedefi=(r.get('uyelik_hedefi') or '').strip(),
                            duzey=(r.get('duzey') or '').strip(),
                            not_=(r.get('not') or '').strip()))
    return out


def kutu_adi_normalize(kutu):
    """A1_1_2223 and A1-1_2223 are the same bin (the files of sample A1-1 use an
        underscore). This normalises the class and sample separator always to a HYPHEN.

    """
    if '_' not in kutu:
        return kutu
    bas, _, son = kutu.rpartition('_')          # son = taxid
    return bas.replace('_', '-') + '_' + son


# --------------------------------------------------------------- measurement
# -------------------------------------------------------------------------
# THE MAIN MEASUREMENT. The order: collect the pairs -> resolve the membership ->
# build the read pools -> measure every pair BOTH at mm<=1 (primary, deciding) AND
# at mm<=3 (secondary, robustness).
#
# THE MEMBERSHIP DOES NOT COME FROM THE KRAKEN LABEL, it comes from stage U's
# MEASURED identity. When a mislabelled bin is written into the competitor column,
# the metric compares the target against the target and even a perfect primer comes
# out below 1 (a measured example: 0.71x -> 8.47x, the same primer, only the
# membership corrected).
#
# A DEPTH CAP IS REQUIRED: without one, a bin of 46 000 reads and a bin of 300 reads
# stand side by side in the same table and the deep bin's discrimination comes out
# ARTIFICIALLY high. The cap puts every row on the same statistical ground.
#
# THE CHECKPOINT SEAL (_ayar): the read cap, the mixed bin rule, the script version
# and THE NAME OF THE MEMBERSHIP FILE are all part of the seal. When the membership
# table is refreshed the old measurements are not silently reused; the reason this
# script exists is precisely that "numbers produced under different conditions were
# standing side by side".
# -------------------------------------------------------------------------
def calistir(kok, okuma_tavani, karisik_kural, yalniz=None, sifirla=False):
    sys.path.insert(0, kok)
    from screening import config as C, motor, numune as N, hedefler as H

    CIKTI = os.path.join(kok, 'TEK_PROTOKOL_SONUC')
    KONTROL = os.path.join(CIKTI, 'kontrol')
    os.makedirs(KONTROL, exist_ok=True)
    if sifirla:
        for f in os.listdir(KONTROL):
            try:
                os.remove(os.path.join(KONTROL, f))
            except OSError as e:
                print('  silinemedi: %s (%s)' % (f, e))

    gunluk = open(os.path.join(CIKTI, 'kosu_gunlugu.txt'), 'a', encoding='utf-8')

    def yaz(s=''):
        print(s, flush=True)
        gunluk.write(s + '\n'); gunluk.flush()

    yaz('=' * 78)
    yaz(u'  SINGLE PROTOCOL - the whole panel measured again under one rule')
    yaz(u'  version %s   %s' % (VERSIYON, time.strftime('%Y-%m-%d %H:%M')))
    yaz('=' * 78)

    # --- ciftler -------------------------------------------------------
    panel, panel_yolu = H.panel_oku()
    ciftler = []
    for d in panel:
        ciftler.append(dict(hedef=d['hedef'], sinif=d['sinif'], F=d['F'], R=d['R'],
                            urun_bp=d['urun_bp'], ta=d.get('ta', ''),
                            duzey=d.get('duzey', ''), kaynak='PANEL',
                            panel_ayrim=d.get('ayrim', ''), not_=''))
    # The EXTRA pairs are for pairs NOT IN the panel TSV. If a pair is added to the
    # panel and also stays in ek_ciftler.tsv, THE SAME TARGET IS MEASURED TWICE and
    # produces two different numbers. That is exactly what happened on 2026-08-11:
    # Petriella_cinsi came both from the panel (0.88x with its own membership) and from
    # the EXTRA list (11.03x, because its membership target was Petriella_musispora),
    # and two rows appeared in the order list, one "orderable" and the other "below
    # threshold".
    panel_adlari = {c['hedef'].strip() for c in ciftler}
    ek = ek_ciftler_oku(kok)
    for e in ek:
        if e['hedef'].strip() in panel_adlari:
            print(u'  EXTRA skipped (already in the panel): %s' % e['hedef'])
            continue
        ciftler.append(dict(e, panel_ayrim=''))
    if yalniz:
        ciftler = [c for c in ciftler
                   if any(y.strip().lower() in c['hedef'].lower()
                          for y in yalniz.split(','))]
    if not ciftler:
        sys.exit(u'ERROR: no pair to measure was found.')

    # --- uyelik --------------------------------------------------------
    uy_yol = uyelik_dosyasi(kok)
    if not uy_yol:
        sys.exit(u'ERROR: uyelik_yeniden_turetme_uyelik_*.tsv was not found.\n      verification/full_chain.py -> option U has to be run first.')
    uyelik = uyelik_oku(uy_yol)
    yaz(u'  membership source : %s' % os.path.basename(uy_yol))
    yaz(u'  pairs             : %d  (panel %d + extra %d)' % (len(ciftler), len(panel), len(ek)))
    yaz(u'  main criterion    : <=%d mismatches and an EXACT match at the last two 3\' bases' % PROTOKOL['olcut_asil'])
    yaz('  yan olcut      : <=%d uyumsuzluk + 3\' son 2 baz TAM' % PROTOKOL['olcut_yan'])
    yaz(u'  depth             : at most %d reads per bin (NO PER-ROW EXCEPTIONS)' % okuma_tavani)
    yaz(u'  mixed bins        : counted as %s' % karisik_kural.upper())
    yaz(u'  threshold         : %s, judged on the single worst competitor bin'
        % _C.esik_metni())
    yaz(u'  threshold origin  : %s' % _C.ESIK_KOKENI)
    yaz(u'  WARNING           : %s' % _C.ESIK_VERIM_NOTU)
    yaz('')

    # --- kutular -------------------------------------------------------
    kut = {k['kutu']: k for k in H.kutular()}
    eksik_uyari = set()

    def coz(adlar):
        out = []
        for a in adlar:
            a2 = kutu_adi_normalize(a.strip())
            if a2 in kut:
                out.append(kut[a2])
            elif a.strip() in kut:
                out.append(kut[a.strip()])
            else:
                eksik_uyari.add(a.strip())
        return out

    baglam = {}
    for c in ciftler:
        u = uyelik.get(c['hedef'])
        if u is None:
            # an extra pair: if the membership row does not match the target name, use the class default
            u = uyelik.get(c.get('uyelik_hedefi', ''), None)
        if u is None:
            baglam[c['hedef']] = None
            continue
        uye = coz(u['uye'])
        kar = coz(u['karisik'])
        rak = coz(u['rakip'])
        if karisik_kural == 'uye':
            uye, kar_eklenen = uye + kar, []
        elif karisik_kural == 'rakip':
            rak = rak + kar
        # 'disla' -> it is added to neither side
        if not rak:   # uyelik satirinda rakip bos ise sinifin geri kalani rakiptir
            uye_ad = {k['kutu'] for k in uye} | {k['kutu'] for k in kar}
            rak = [k for k in kut.values()
                   if k['sinif'] == (u['sinif'] or c['sinif']) and k['kutu'] not in uye_ad]
        baglam[c['hedef']] = dict(uye=uye, rakip=rak, karisik=kar)

    gerekli = {}
    for b in baglam.values():
        if b:
            for k in b['uye'] + b['rakip']:
                gerekli[k['kutu']] = k
    if eksik_uyari:
        yaz(u'  WARNING: %d bin names from the membership table were not found in the fastq directory: %s'
            % (len(eksik_uyari), ', '.join(sorted(eksik_uyari)[:6])))
        yaz(u'         (the files for the A1-1 sample are named with underscores; the script')
        yaz(u'          corrects that internally, but anything still missing is listed above.)')
    yaz(u'  bins to read      : %d' % len(gerekli))
    yaz('')
    # --- havuz kurulumu -------------------------------------------------
    def ilerK(i, n, ad):
        print(u'   ... read pool %d/%d  %s          ' % (i, n, ad), end='\r', flush=True)

    t0 = time.time()
    yaz(u'Building read pools (%d bins). Only bin names scroll past on screen during' % len(gerekli))
    yaz(u'this step; it is NOT stuck. The real measurement starts after this and each pair is saved separately.')
    nm = N.Numune(list(gerekli.values()), n=okuma_tavani, ilerle=ilerK, otorite=True)
    top_okuma = sum(h.n_okuma for h in nm.havuz.values())
    yaz(u'\nPools ready: %d bins, %d reads  (%s)' % (len(gerekli), top_okuma, sure_metni(time.time() - t0)))
    tahmin = len(ciftler) * 2 * max(1.0, top_okuma / 20000.0)
    yaz(u'ESTIMATED MEASUREMENT TIME: ~%s   (resumable; the same command continues)'
        % sure_metni(tahmin))
    yaz('')

    # --- measurement -----------------------------------------------------
    def kp_yolu(ad):
        t = ''.join(ch if ch.isalnum() else '_' for ch in ad)
        return os.path.join(KONTROL, 'cift_%s.json' % t)

    # O-10: the membership source IS PART of the seal. When the membership table is
    # refreshed the old checkpoints must not be silently reused.
    #
    # THE 2026-08-11 FIX (the membership CONTENT seal). The seal held only the NAME of
    # the membership file. When the file is corrected IN PLACE the name does not change,
    # the seal matches, and the measurement comes back with the OLD membership as "taken
    # from the previous run". That is exactly what happened: four protist bins were
    # removed from the Mantar_universal (F2) membership, the measurement was re-run, and
    # both targets came back from the cache; the change never reached the numbers. The
    # same bug as the sequence seal fixed on 10 August, this time on the membership
    # side. The md5 of the file's CONTENT is now in the seal too: a membership table
    # with a changed row invalidates the checkpoint.
    import hashlib as _hl0
    with open(uy_yol, 'rb') as _fh0:          # okunamazsa PATLASIN: sessiz
        _uy_muhru = _hl0.md5(_fh0.read()).hexdigest()[:12]   # "okunamadi" muhru
    # (there used to be a try/except here, and because the io module was not imported
    #  the seal came out "could not be read" on every run, that is, constant. A constant
    #  seal is not a seal: it would match even when the membership table changed.
    #  Swallowing the error made the check itself invisibly useless.)
    AYAR = dict(okuma=okuma_tavani, karisik=karisik_kural, surum=VERSIYON,
                uyelik=os.path.basename(uy_yol), uyelik_icerik=_uy_muhru)

    # THE 2026-08-10 FIX (the sequence seal). The primer SEQUENCE WAS NOT in the seal.
    # The consequence: when a pair's forward or reverse sequence was changed the
    # checkpoint still counted as valid and the OLD measurement came back as "taken from
    # the previous run". On 10 August two pairs were changed and two separate full runs
    # (5 h 29 min plus 2 h 0 min) measured the old sequences and took them for the new
    # ones. The md5 of each pair's own F+R sequence now goes into its seal; if the
    # sequence changes the checkpoint invalidates itself.
    import hashlib as _hl

    def _ayar_of(c):
        d = dict(AYAR)
        d['dizi'] = _hl.md5(
            ((c.get('F') or '') + '|' + (c.get('R') or '')).encode('utf-8')
        ).hexdigest()[:12]
        return d
    sonuc = []
    tb = time.time()
    for i, c in enumerate(ciftler, 1):
        kp = kp_yolu(c['hedef'] + '|' + c['kaynak'])
        if os.path.exists(kp):
            try:
                v = json.load(open(kp, encoding='utf-8'))
                # Backward compatibility after O-10: old checkpoints have no 'uyelik'
                # key. If the SHARED keys match it is accepted, but a WARNING is printed;
                # it is not silently reused.
                _e = v.get('_ayar') or {}
                _bek = _ayar_of(c)
                # The SEQUENCE seal invalidates the checkpoint not when it is equal but
                # when it EXISTS and DIFFERS. Old checkpoints have no 'dizi' key at all;
                # in that case it is re-measured too. There is NO silent acceptance.
                if _e.get('dizi') != _bek['dizi']:
                    yaz(u'  %s: the checkpoint\'s SEQUENCE seal does not match (recorded %s, now %s); re-measuring.'
                        % (c['hedef'][:40], _e.get('dizi') or 'yok', _bek['dizi']))
                    raise ValueError('dizi muhru tutmadi')
                # The MEMBERSHIP CONTENT seal works like the sequence seal: if it is
                # absent or does not match, the checkpoint is INVALID. "The shared keys
                # matched" is not a reason to accept it; the missing key is precisely the
                # one that sees a change in the membership table. (2026-08-11: because of
                # this backward compatibility path the protist fix never reached the
                # measurement, and 22 of 22 targets came back from the cache.)
                if _e.get('uyelik_icerik') != _bek['uyelik_icerik']:
                    yaz(u'  %s: the checkpoint\'s MEMBERSHIP seal does not match (recorded %s, now %s); re-measuring.'
                        % (c['hedef'][:40], _e.get('uyelik_icerik') or 'yok',
                           _bek['uyelik_icerik']))
                    raise ValueError('uyelik muhru tutmadi')
                _ortak = {k: _e.get(k) for k in _bek if k in _e}
                _uyar = (_e != _bek and _ortak == {k: _bek[k] for k in _ortak})
                if _e == _bek or _uyar:
                    if _uyar:
                        yaz(u'  WARNING: the %s checkpoint was written with an OLD settings seal (missing: %s). It was accepted because the shared keys match'
                            % (c['hedef'][:40], ', '.join(sorted(set(AYAR) - set(_e))) or '-'))
                    sonuc.append(v)
                    yaz(u'[%2d/%2d] %-46s  (taken from the previous run)' % (i, len(ciftler), c['hedef'][:46]))
                    continue
            except Exception:
                pass
        b = baglam.get(c['hedef'])
        r = dict(c); r['_ayar'] = _ayar_of(c); r['olcum'] = {}
        if not b or not b['uye']:
            r['hata'] = 'there is no membership, or no member bin was found'
            yaz(u'[%2d/%2d] %-46s  SKIPPED (%s)' % (i, len(ciftler), c['hedef'][:46], r['hata']))
        else:
            for mm in (PROTOKOL['olcut_asil'], PROTOKOL['olcut_yan']):
                o = nm.olc(c['F'], c['R'], b['uye'], b['rakip'],
                           lo=PROTOKOL['urun_alt'], hi=PROTOKOL['urun_ust'], mm=mm)
                r['olcum'][str(mm)] = o
            r['uye_n'] = len(b['uye']); r['rakip_n'] = len(b['rakip']); r['karisik_n'] = len(b['karisik'])
            o1 = r['olcum'][str(PROTOKOL['olcut_asil'])]
            d1, g1, day1 = karar(o1, c['hedef'], c.get('duzey', ''))
            yaz(u'[%2d/%2d] %-46s  %s x  %-11s | coverage %s'
                % (i, len(ciftler), c['hedef'][:46], vir(g1), d1,
                   o1.get('uye_kapsam_pay') if o1 else '-'))
        json.dump(r, open(kp, 'w', encoding='utf-8'), ensure_ascii=False, indent=1, default=str)
        sonuc.append(r)
        gecen = time.time() - tb
        print(u'        elapsed %s | estimated remaining %s'
              % (sure_metni(gecen), sure_metni(gecen / i * (len(ciftler) - i))), flush=True)

    yaz('')
    yaz(u'Measurement finished (%s). Writing outputs...' % sure_metni(time.time() - tb))
    raporla(CIKTI, sonuc, dict(uyelik=os.path.basename(uy_yol), okuma=okuma_tavani,
                               karisik=karisik_kural, panel=os.path.basename(panel_yolu)), yaz)
    rc = cikti_denetle(yaz, 'P (TEK PROTOKOL)', [
        (os.path.join(CIKTI, 'panel_tek_protokol.tsv'), 'panel_tek_protokol.tsv'),
        (os.path.join(CIKTI, 'SIPARIS_LISTESI.tsv'), 'SIPARIS_LISTESI.tsv')])
    gunluk.close()
    return rc


# --------------------------------------------------------------- ciktilar
def _o(r, mm):
    return (r.get('olcum') or {}).get(str(mm)) or {}


# -------------------------------------------------------------------------
# THE DISCRIMINATION RATIO IS UNDEFINED ON UNIVERSAL TARGETS.
# The discrimination ratio = (the member lower bound) / (the competitor upper
# bound). Bakteri_universal was designed to amplify all bacteria and Arke_universal
# all archaea; on those rows there IS NO set called "competitors". As the competitor
# set approaches empty the denominator goes to zero and the ratio is either 0/0 or a
# huge number. Indeed, the old panel had 0.00 and 117 million standing side by side
# in the same column. Those numbers measure nothing.
#
# So these rows get NO numeric verdict here; they are marked OLCULEMEDI. The right
# measure is COVERAGE plus the OUTSIDE THE DOMAIN proportion, and it is applied at
# stage K. This IS NOT LOWERING the 10x threshold: since the ratio's denominator is
# undefined, that threshold cannot be applied on these rows at all. On every other
# row the 10x stands exactly as it was.
# -------------------------------------------------------------------------
def evrensel_mi(hedef, duzey=''):
    """O-6: evrensel/alan hedeflerinde ayrim katinin PAYDASI tanimsizdir
    (rakip kumesi yoktur). Bu satirlar sayisal verdikt almamalidir; dogru olcu
    KAPSAMA + ALAN DISI'dir ve K asamasinda uygulanir."""
    ad = (hedef or '').lower()
    return ('universal' in ad or 'evrensel' in ad
            or (duzey or '').strip().lower() == 'alan')


# -------------------------------------------------------------------------
# The decision is made on THE WORST SINGLE COMPETITOR BIN, not on the pool: the pool
# ratio dissolves one bad bin among thousands of clean reads, and a pair that really
# does cross react looks clean. The pool ratio is still reported, but it is not the
# deciding column.
#
# THERE ARE THREE SEPARATE STATES and they are never confused:
#   ESIK USTU / ESIK ALTI - measured, decided.
#   OLCULEMEDI            - there is no competitor bin at sufficient depth, or the
#                           target is universal. That IS NOT "below threshold";
#                           writing an absence of decision into the same column as a
#                           failure would be wrong.
# When the worst bin measure cannot be produced it falls back to the pool, but that
# is written OPENLY in the basis column.
# -------------------------------------------------------------------------
def karar(o, hedef='', duzey=''):
    """(state, value, basis). The decision is made on THE WORST SINGLE COMPETITOR BIN;
        if that measure cannot be produced (no competitor bin at sufficient depth) it
        falls back to the pool and that is marked OPENLY. With neither available it is
        OLCULEMEDI, which IS NOT 'below threshold'.

    """
    if not o:
        return ('OLCULEMEDI', None, 'olcum yok')
    if evrensel_mi(hedef, duzey):
        return ('OLCULEMEDI', None,
                'EVRENSEL HEDEF - ayrim katinin paydasi tanimsiz. Dogru olcu '
                'KAPSAMA + ALAN DISI; K asamasinda uygulanir.')
    g = o.get('kat_enkotu')
    if g is not None:
        return ('ESIK USTU' if g >= PROTOKOL['esik'] else 'ESIK ALTI', g, 'en kotu tek kutu')
    h = o.get('kat_havuz')
    if h is not None:
        return ('ESIK USTU' if h >= PROTOKOL['esik'] else 'ESIK ALTI', h,
                'HAVUZ (yeterli derinlikte rakip kutu yok - en kotu kutu olcusu uretilemedi)')
    return ('OLCULEMEDI', None, 'rakip kutu yok')


# -------------------------------------------------------------------------
# It produces three outputs: the full table, the raw counts and the order list.
#
# THE RAW COUNTS (k = reads giving a product, n = reads in the bin) are written to a
# separate file, because until now NO READER COULD RECOMPUTE a verdict. Every
# verdict derives from those two numbers; publishing them makes the decision rules
# auditable.
#
# THE TWO CRITERIA STAND IN SEPARATE COLUMNS and NEITHER STANDS IN for the other:
#   ayrim_mm1_ARAC_OLCUTU        - 10x, this tool's criterion (NOT a meeting decision).
#   TOPLANTI_OLCUTU_capraz_kutu  - the NUMBER OF competitor BINS giving over 10%
#                                  product (WORK_RECORD 1.7, a tolerance of 1-2
#                                  cross reacting species).
# The dCq column is the same number in laboratory language (dCq = log2(fold), on the
# assumption of 100% efficiency); it is not a new criterion but a translation of the
# same measure.
#
# THE FLAGS do not reject a pair, they make it CONDITIONAL: sensitive to the
# criterion (it collapses at mm<=3), a shallow deciding bin, one or two member bins,
# a partial measurement, only the worst bin passed. A pair with no flag at all that
# passes both criteria is marked UNCONDITIONAL.
# -------------------------------------------------------------------------
def raporla(CIKTI, sonuc, meta, yaz):
    E = PROTOKOL['esik']; A = PROTOKOL['olcut_asil']; Y = PROTOKOL['olcut_yan']
    basli = (u'# Bu dosya TEK PROTOKOLLE uretildi - butun satirlar ayni kural ve ayni derinlik.\n'
             u'# uyelik kaynagi : %(uyelik)s   (Kraken etiketi KULLANILMADI)\n'
             u'# derinlik       : kutu basina en cok %(okuma)d okuma, satir bazinda istisna YOK\n'
             u'# asil olcut     : <=1 uyumsuzluk + 3\' son 2 baz TAM  (karar bu sutuna gore)\n'
             u'# yan olcut      : <=3 uyumsuzluk + 3\' son 2 baz TAM  (dayaniklilik gostergesi)\n'
             u'# karisik kutu   : %(karisik)s sayildi\n'
             u'# esik           : %(esik)s, EN KOTU TEK RAKIP KUTU uzerinden\n'
             u'# esik kokeni    : %(koken)s\n'
             u'# VERIM UYARISI  : %(verim)s\n') % dict(
                 meta, esik=_C.esik_metni(E), koken=_C.ESIK_KOKENI,
                 verim=_C.ESIK_VERIM_NOTU)

    # ---------- 1) tam tablo ----------
    yol = os.path.join(CIKTI, 'panel_tek_protokol.tsv')
    with open(yol, 'w', encoding='utf-8', newline='') as fh:
        fh.write(basli)
        w = csv.writer(fh, delimiter='\t')
        _kokp = os.path.dirname(os.path.abspath(CIKTI))
        _kimp = _S.kimlik_tablosu(_kokp)
        _kynp = _S.kaynak_tablosu(_kokp)
        w.writerow(['hedef'] + _S.kimlik_sutun_basliklari() + _S.kaynak_sutun_basliklari()
                   + ['kaynak', 'sinif', 'urun_bp', 'F', 'R',
                    'ASIL_ayrim_mm1', 'ASIL_ayrim_havuz_mm1', 'ASIL_kapsam_mm1',
                    'YAN_ayrim_mm3', 'YAN_ayrim_havuz_mm3', 'YAN_kapsam_mm3',
                    'uye_kutu', 'karisik_kutu', 'rakip_kutu',
                    'uye_alt_%', 'en_kotu_rakip_kutu', 'esik_gecti_mi',
                    'karar_dayanagi', 'olcute_duyarli_mi', 'panelin_eski_degeri',
                    'rakip_olculen', 'rakip_toplam', 'kismi_olcum_mu', 'not'])
        for r in sonuc:
            o1, o3 = _o(r, A), _o(r, Y)
            d1, g1, day1 = karar(o1, r['hedef'], r.get('duzey', ''))
            d3, g3, _ = karar(o3, r['hedef'], r.get('duzey', ''))
            # O-5: the floor of 150 reads was silently discarding some of the competitor bins.
            ro, rt = o1.get('rakip_olculen'), o1.get('rakip_toplam')
            kismi = 'hayir'
            if ro is not None and rt:
                if ro == 0:
                    kismi = 'EVET - hicbir rakip kutu 150 okumayi gecmedi'
                elif ro < rt / 2.0:
                    kismi = 'EVET - rakiplerin %d/%d si olcume girdi' % (ro, rt)
            if kismi.startswith('EVET') and d1 == 'ESIK USTU':
                d1 = 'ESIK USTU (KISMI OLCUM)'
            gecti = d1
            duyarli = 'EVET' if (d1 == 'ESIK USTU' and d3 == 'ESIK ALTI') else 'hayir'
            w.writerow([r['hedef']] + _S.kimlik_sutunlari(_kimp, r['hedef'])
                       + _S.kaynak_sutunlari(_kynp, r['hedef'])
                       + [r['kaynak'], r.get('sinif', ''), r.get('urun_bp', ''),
                        r.get('F', ''), r.get('R', ''),
                        vir(g1), vir(o1.get('kat_havuz')), o1.get('uye_kapsam_pay', ''),
                        vir(g3), vir(o3.get('kat_havuz')), o3.get('uye_kapsam_pay', ''),
                        r.get('uye_n', ''), r.get('karisik_n', ''), r.get('rakip_n', ''),
                        vir(o1.get('uye_alt'), 3), o1.get('enkotu_kutu', ''),
                        gecti, day1, duyarli, r.get('panel_ayrim', ''),
                        o1.get('rakip_olculen', ''), o1.get('rakip_toplam', ''),
                        kismi, r.get('hata', '')])
    yaz(u'  written: %s' % yol)

    # ---------- 1b) THE RAW COUNTS PER BIN (item 7e) ----------
    # Until now NO READER COULD RECOMPUTE a verdict.
    # k = the reads giving a product, n = the total reads in the bin. Every verdict
    # derives from those two numbers; publishing them makes every decision rule auditable.
    yolk = os.path.join(CIKTI, 'kutu_bazli_ham_sayilar.tsv')
    with open(yolk, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# RAW NUMBERS - every verdict is derived from these two columns.\n')
        fh.write(u'# k = reads that gave a product, n = reads in the bin. ratio = k/n.\n')
        fh.write(u'# Wilson: LOWER bound on the member side, UPPER bound on the competitor side (z=1.96).\n')
        fh.write(basli)
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'olcut_mm', 'kutu', 'grup', 'k', 'n', 'oran_%'])
        for r in sonuc:
            for mm in (A, Y):
                o = _o(r, mm)
                for grup in ('uye', 'rakip'):
                    for satir in (o.get(grup) or []):
                        try:
                            ad, k, n, yuzde = satir[0], satir[1], satir[2], satir[3]
                        except (IndexError, TypeError):
                            continue
                        w.writerow([r['hedef'], mm, ad, grup, k, n, vir(yuzde)])
    yaz(u'  written: %s' % yolk)

    # ---------- 2) TEK siparis listesi ----------
    gecen, kalan, olculemeyen = [], [], []
    for r in sonuc:
        d = karar(_o(r, A), r.get('hedef', ''), r.get('duzey', ''))[0]
        (gecen if d == 'ESIK USTU' else kalan if d == 'ESIK ALTI' else olculemeyen).append(r)
    yol2 = os.path.join(CIKTI, 'SIPARIS_LISTESI.tsv')
    with open(yol2, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# SINGLE ORDER LIST - produced with one protocol.\n')
        fh.write(basli)
        # --- 2026-08-06: below-threshold rows ARE NOT REMOVED FROM THE LIST, they are classified.
        _snf = {}
        for _r in sonuc:
            _d, _g, _ = karar(_o(_r, A), _r['hedef'], _r.get('duzey', ''))
            _snf[_r['hedef']] = _S.sinifla(_g, olculemedi=(_d == 'OLCULEMEDI'))
        _kesin = [x for x in _snf.values() if x[0] in (u'KESIN', u'EVRENSEL')]
        _kos = [x for x in _snf.values() if x[0] == u'KOSULLU']
        _oner = [x for x in _snf.values() if x[0] == u'ONERILMEZ']
        fh.write(u'#\n')
        fh.write(u'# ================  THREE NUMBERS, SEPARATELY  ================\n')
        fh.write(u'#   CERTAIN     : %d pairs = %d oligos   (dCq >= %.1f, or universal/coverage)\n'
                 % (len(_kesin), 2 * len(_kesin), _C.ESIK_DCQ))
        fh.write(u'#   CONDITIONAL : %d pairs = %d oligos   (dCq %.1f-%.1f - orderable BUT validation is required)\n'
                 % (len(_kos), 2 * len(_kos), _S.KOSULLU_ALT_DCQ, _C.ESIK_DCQ))
        fh.write(u'#   NOT ADVISED : %d pairs              (dCq < %.1f - stays on the list, with its reason written out)\n'
                 % (len(_oner), _S.KOSULLU_ALT_DCQ))
        fh.write(u'# A row that fails the threshold is NEVER DELETED SILENTLY; the decision is yours.\n')
        fh.write(u'# (tool counts: above threshold %d, below threshold %d, not measurable %d)\n'
                 % (len(gecen), len(kalan), len(olculemeyen)))
        fh.write(u'# Rows marked "OLCUTE DUYARLI" pass at mm<=1 but collapse at mm<=3, so they are fragile.\n')
        fh.write(u'#\n')
        fh.write(u'# THRESHOLD: %s\n' % _C.esik_metni())
        fh.write(u'# ESIGIN KOKENI: %s\n' % ESIK_KOKENI)
        fh.write(u'# EFFICIENCY WARNING: %s\n' % ESIK_VERIM_NOTU)
        fh.write(u'#   The agreed criterion is a DIFFERENT quantity: the NUMBER OF CROSS-REACTING SPECIES\n')
        fh.write(u'#   (WORK_RECORD 1.7, a tolerance of 1 to 2). So it goes in its own column:\n')
        fh.write(u'#   TOPLANTI_OLCUTU_capraz_kutu = how many competitor bins give >=%%%d product.\n'
                 % int(TOPLANTI_CAPRAZ_TABAN))
        fh.write(u'#   The two criteria DO NOT STAND IN for one another.\n')
        fh.write(u'# dCq_karsiligi: laboratuvarin konustugu birim. dCq = log2(kat),\n')
        fh.write(u'#   %%100 verim varsayimiyla. 10x = 3,32 dongu. Literaturde ozgulluk\n')
        fh.write(u'#   the passing criterion is dCq >= 3 (NEB). It must be corrected once the real efficiency is measured.\n')
        w = csv.writer(fh, delimiter='\t')
        _kok0 = os.path.dirname(os.path.abspath(CIKTI))
        _kim = _S.kimlik_tablosu(_kok0)
        _kyn = _S.kaynak_tablosu(_kok0)
        w.writerow(['sira', 'SINIF', 'durum', 'siparis_sarti', 'hedef']
                   + _S.kimlik_sutun_basliklari() + _S.kaynak_sutun_basliklari()
                   + ['oligo_adi_F', 'F', 'oligo_adi_R', 'R', 'urun_bp',
                      'ayrim_mm1', 'dCq_karsiligi', 'esikten_uzaklik_dCq',
                      'LABORATUVARDA_NE_YAPILMALI',
                      'ayrim_mm3', 'havuz_mm1',
                      'TOPLANTI_OLCUTU_capraz_kutu', 'kapsam_mm1',
                      'karar_veren_kutu', 'karar_kutusu_k', 'karar_kutusu_n',
                      'uye_kutu_sayisi', 'damgalar'])
        n = 0
        for etiket, kume in (('ESIK USTU - SIPARIS EDILEBILIR', gecen),
                             ('ESIK ALTI - SIPARIS EDILMEZ', kalan),
                             ('OLCULEMEDI - KARAR YOK', olculemeyen)):
            for r in kume:
                n += 1
                o1, o3 = _o(r, A), _o(r, Y)
                d1, g1, day1 = karar(o1, r['hedef'], r.get('duzey', ''))
                d3, g3, _ = karar(o3, r['hedef'], r.get('duzey', ''))
                kod = ''.join(ch if ch.isalnum() else '_' for ch in r['hedef'])[:24]
                # --- item 2: THE MEETING CRITERION = the number of cross reacting bins (NOT an efficiency ratio)
                capraz = sum(1 for x in (o1.get('rakip') or [])
                             if len(x) > 3 and (x[3] or 0) >= TOPLANTI_CAPRAZ_TABAN)
                # --- item 3: the RAW counts of the deciding bin
                kk, kn = '', ''
                ek = o1.get('enkotu_kutu')
                for x in (o1.get('rakip') or []):
                    if x and x[0] == ek:
                        kk, kn = x[1], x[2]; break
                # --- madde 4 + 3 + 7a: damgalar ve siparis sarti
                damga = []
                if d1 == 'ESIK USTU' and d3 == 'ESIK ALTI':
                    damga.append(u'OLCUTE DUYARLI (mm<=3 te cokuyor)')
                if kn and int(kn) < 300:
                    damga.append(u'SIG KARAR KUTUSU (n=%s)' % kn)
                if (r.get('uye_n') or 0) and int(r['uye_n']) <= 2:
                    damga.append(u'ONE OR TWO MEMBER BINS, the within target variability WAS NOT TESTED')
                if kismi.startswith('EVET'):
                    damga.append(kismi)
                hv = o1.get('kat_havuz')
                if d1 == 'ESIK USTU' and (hv is None or hv < E):
                    damga.append(u'ONLY the worst bin passed, the POOL did not')
                # an unconditional order: both criteria and both floors at once
                if (d1 == 'ESIK USTU' and d3 == 'ESIK USTU'
                        and hv is not None and hv >= E and not damga):
                    sart = 'KOSULSUZ'
                elif d1 == 'ESIK USTU':
                    sart = 'KOSULLU'
                else:
                    sart = '-'
                _sn, _dq, _uz, _lab = _snf[r['hedef']]
                w.writerow([n, _sn, etiket, sart, r['hedef']]
                           + _S.kimlik_sutunlari(_kim, r['hedef'])
                           + _S.kaynak_sutunlari(_kyn, r['hedef'])
                           + [kod + '_F', r.get('F', ''),
                              kod + '_R', r.get('R', ''), r.get('urun_bp', ''),
                              vir(g1), vir(_dq), vir(_uz), _lab,
                              vir(g3), vir(hv), capraz,
                              o1.get('uye_kapsam_pay', ''), ek, kk, kn,
                              r.get('uye_n', ''), '; '.join(damga) or '-'])
    yaz(u'  written: %s' % yol2)

    # ---------- 3) the protocol plus the report ----------
    yol3 = os.path.join(CIKTI, 'PROTOKOL_VE_RAPOR.md')
    with open(yol3, 'w', encoding='utf-8') as fh:
        fh.write(u'# The panel measurement with a single protocol\n\n')
        fh.write(u'Generated: %s, script version %s\n\n' % (time.strftime('%Y-%m-%d %H:%M'), VERSIYON))
        fh.write(u'## Result\n\n')
        fh.write(u'- Pairs PASSING the threshold (%.0fx): **%d** -> **%d oligos**\n' % (E, len(gecen), 2 * len(gecen)))
        fh.write(u'- Pairs BELOW the threshold: **%d**\n' % len(kalan))
        fh.write(u'- NOT MEASURABLE (no verdict, NOT the same as below threshold): **%d**\n' % len(olculemeyen))
        duy = [r for r in gecen if karar(_o(r, Y))[0] == 'ESIK ALTI']
        fh.write(u'- Passing but SENSITIVE TO THE CRITERION (collapses at mm<=3): **%d**\n\n' % len(duy))
        fh.write(u'```\n' + GEREKCE + u'\n```\n\n')
        fh.write(u'## Tablo\n\n')
        fh.write(u'| target | source | mm<=1 (main) | mm<=3 (secondary) | coverage | status |\n|---|---|---|---|---|---|\n')
        for r in sorted(sonuc, key=lambda x: -(karar(_o(x, A))[1] if karar(_o(x, A))[1] is not None else -1)):
            o1, o3 = _o(r, A), _o(r, Y)
            d1, g1, day1 = karar(o1, r['hedef'], r.get('duzey', ''))
            d3, g3, _ = karar(o3, r['hedef'], r.get('duzey', ''))
            fh.write(u'| %s | %s | %s | %s | %s | %s |\n' % (
                r['hedef'], r['kaynak'], vir(g1), vir(g3),
                o1.get('uye_kapsam_pay', '-'),
                d1 + ('' if day1 == u'worst single bin' else ' (%s)' % day1)))
        fh.write(u'\n## Reading order\n\n1. This file first. 2. `SIPARIS_LISTESI.tsv`. 3. `panel_tek_protokol.tsv` for detail.\n')
    yaz(u'  written: %s' % yol3)
    yaz('')
    yaz(u'  PASSED THE THRESHOLD: %d pairs (%d oligos)   BELOW THRESHOLD: %d   NOT MEASURABLE: %d'
        % (len(gecen), 2 * len(gecen), len(kalan), len(olculemeyen)))



# --------------------------------------------------------------- guvenlik agi
def cikti_denetle(yaz, ad, dosyalar, asgari=1):
    """When the stage ends, it audits ITS OWN output.

        If the expected row count is zero, or the file is missing entirely, it DOES NOT
        CARRY ON SILENTLY: it prints a clear error and returns a non-zero code. This is
        so that it cannot produce an empty result overnight and then say "nothing was
        found" in the morning.

    """
    sorun = []
    for yol, etiket in dosyalar:
        if not os.path.exists(yol):
            sorun.append(u'%s WAS NOT PRODUCED (%s)' % (etiket, yol)); continue
        try:
            with open(yol, encoding='utf-8') as fh:
                n = sum(1 for s in fh if s.strip() and not s.startswith('#'))
            n = max(0, n - 1)          # baslik satiri
        except OSError as e:
            sorun.append(u'%s COULD NOT BE READ (%s)' % (etiket, e)); continue
        if n < asgari:
            sorun.append(u'%s IS EMPTY, %d data rows (at least %d were expected)'
                         % (etiket, n, asgari))
    if not sorun:
        return 0
    yaz('')
    yaz('  ' + '!' * 70)
    yaz(u'  STAGE %s PRODUCED EMPTY OUTPUT - THE CHAIN WAS STOPPED HERE' % ad)
    for x in sorun:
        yaz(u'    - %s' % x)
    yaz('')
    yaz(u'  WHY IT STOPPED: the next stage would have read this file as input.')
    yaz(u'  Continuing with empty input does not crash; it produces a MEANINGLESS BUT')
    yaz(u'  CONVINCING summary, which is exactly the silent failure we hunt for.')
    yaz(u'  Read the run log above, fix the cause, then run the same command')
    yaz(u'  again; finished work is skipped from its checkpoints.')
    yaz('  ' + '!' * 70)
    return 4


def girdi_denetle(yaz, ad, dosyalar):
    """Before the stage STARTS: do the files it needs exist, and are they non-empty?"""
    eksik = []
    for yol, etiket, uretici in dosyalar:
        if not os.path.exists(yol):
            eksik.append(u'there is no %s (%s) -> run stage %s first' % (etiket, yol, uretici))
            continue
        with open(yol, encoding='utf-8') as fh:
            n = sum(1 for s in fh if s.strip() and not s.startswith('#'))
        if n <= 1:
            eksik.append(u'%s IS EMPTY (%s) -> stage %s produced no result'
                         % (etiket, yol, uretici))
    if not eksik:
        return 0
    yaz('')
    yaz('  ' + '!' * 70)
    yaz(u'  STAGE %s WAS NOT STARTED - INPUT MISSING' % ad)
    for x in eksik:
        yaz(u'    - %s' % x)
    yaz('  ' + '!' * 70)
    return 5

# The command line: --reads the depth cap, --mixed what to do with mixed bins
# (uye|rakip|disla; the default rakip measures the worst case), --only a subset,
# --reset deletes the checkpoints.

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
    p = argparse.ArgumentParser(description='Measure the panel with a single '
                                            'protocol')
    p.add_argument('--root', dest='kok', default='.')
    p.add_argument('--reads', dest='okuma', type=int, default=PROTOKOL['okuma_tavani'],
                   help='cap on reads per bin (0 = all of them)')
    p.add_argument('--mixed', dest='karisik', choices=['member', 'competitor', 'exclude', 'uye', 'rakip', 'disla'], default=PROTOKOL['karisik'])
    p.add_argument('--only', dest='yalniz', default=None, help='only targets whose name contains this (for testing)')
    p.add_argument('--reset', dest='sifirla', action='store_true')
    a = p.parse_args()
    a = _ing_deger(a)
    kok = kok_bul(a.kok)
    return calistir(kok, a.okuma, a.karisik, a.yalniz, a.sifirla)


if __name__ == '__main__':
    sys.exit(main() or 0)
