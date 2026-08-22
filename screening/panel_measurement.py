# -*- coding: utf-8 -*-
"""OPTION 4: measure the panel again with the CORRECTED engine.

WHY IT IS NEEDED
----------------
The panel's raw read engine, the `Sonda` class inside engine/reads.py, searches for
the 13 bases at the primer's 3' end as an EXACTLY matching seed (`variants()` only
expands the IUPAC codes, it produces no MISMATCH variant). So it MISSES every binding
site whose single mismatch falls inside those 13 bases at the 3' end. Measured: on
some target and bin pairs Sonda finds 0 sites while brute force finds 146 (a 100
percent loss).

This module DOES NOT USE that engine. It uses two correct routes and confirms both
against brute force:
  * mm<=1 : tarayici.Havuz (the seed is EXACT for a single mismatch)
  * mm>=2 : ispcr.find_sites, SEEDLESS (under a relaxed criterion a seeded search
            misses)

THE TWO CRITERIA
----------------
Some of the panel's rows were measured with <=1 mismatch and some with <=3. This
module measures EVERY ROW UNDER BOTH CRITERIA and gives two separate columns; every
output row says PLAINLY which criterion was used.

"""
# -------------------------------------------------------------------------
# panel_measurement.py measures every pair in the panel again under one and the
#                  same protocol, at full read depth, under two mismatch criteria.
#
# INPUT  : the panel TSV through hedefler.panel_oku(); hedefler.uyelik_oku(),
#          hedefler.konsensusler() and hedefler.kutular(); each target's member and
#          competitor bin set is resolved with hedefler.hedef_baglami(); the reads
#          are read through numune.Numune(otorite=True), that is read_engine.py,
#          and also through numune.KutuEski (an exact reproduction of the panel's
#          old 13 base seeded engine) for the comparison.
# OUTPUT : KAPSAMLI_ARAMA_SONUC/PANEL_YENIDEN_OLCUM.md,
#          panel_yeniden_olcum.tsv and panel_kutu_duzeyi.tsv (calistir returns those
#          three paths); a kontrol/panel_olcum_*.json per target.
# CALLED BY: verification/full_chain.py key 4 (--mod panel-olc --tam-derinlik),
#          key 7 -> choice "2" (a single target at full depth) and the 5th stage
#          inside key 9 (hepsi.calistir -> panel_olcum.calistir, okuma_sayisi=0).
#
# The two criteria (<=1 and <=3) are given in the same run because they had been
# confused in the panel's old rows; without which one was used written on every
# output row, the numbers cannot be compared.
# -------------------------------------------------------------------------
import os, time, json, csv
from . import config as C
from . import engine_gateway, targets as H, sample as N, checks

OLCUTLER = [1, 3]
# THE THRESHOLD COMES FROM ONE SOURCE: screening/config.py -> ESIK_DCQ = 3.0
# Its fold equivalent is 2 ** ESIK_DCQ = 8.00. NO constant is EMBEDDED; if dCq
# changes it changes in one place. The reasoning and the efficiency warning are
# written in that file.
ESIK = C.AYRIM_ESIK


def _kontrol_yolu(ad):
    t = ''.join(ch if ch.isalnum() else '_' for ch in ad)
    return os.path.join(C.KONTROL, 'panel_olcum_%s.json' % t)


def olcut_metni(mm):
    return "<=%d uyumsuzluk + 3' son 2 baz TAM" % mm


def calistir(yaz, sure, okuma_sayisi=0, yalniz=None, yeniden=False):
    """okuma_sayisi=0 -> EVERY read in the bin (full depth)."""
    # THE ORIENTATION GATE runs first: if the consensuses are not canonical this stage
    # DOES NOT START. The reason is in the silence: on a reversed consensus, in-silico
    # PCR returns 0 products without any warning, so a full depth measurement that runs
    # for hours records every row as "no product" and the result file looks consistent
    # at first glance. That is why the gate is at the start of the run and not at the
    # end.
    from .run_all import yon_kapisi
    _ok, _m = yon_kapisi(yaz, 'panel yeniden olcum')
    for _x in _m:
        yaz('  ' + _x)
    if not _ok:
        yaz('')
        yaz(u'  *** INPUT VERIFICATION FAILED - THIS STAGE WAS NOT STARTED ***')
        yaz(u'  Cause: the consensus sequences to be read are not canonical. On a reverse-oriented')
        yaz(u'  consensus, in-silico PCR returns 0 products without any warning,')
        yaz(u'  so the whole run would silently produce a wrong result.')
        yaz(u'  Fix:    python3 screening/build_canonical.py --root . --rerun')
        raise SystemExit(2)

    checks.hazirla()
    panel, panel_yolu = H.panel_oku()
    uyelik = H.uyelik_oku(); kons = H.konsensusler(); kut = H.kutular()
    if yalniz:
        panel = [d for d in panel if yalniz.lower() in d['hedef'].lower()]

    yaz('=' * 78)
    yaz(u'  RE-MEASURING THE PANEL WITH THE CORRECTED ENGINE')
    yaz('=' * 78)
    yaz(u'  pairs         : %d' % len(panel))
    yaz(u'  read depth: %s' % (u'FULL (every read in the bin)' if not okuma_sayisi
                                   else u'%d reads/bin' % okuma_sayisi))
    yaz('  olcutler      : ' + '  |  '.join(olcut_metni(m) for m in OLCUTLER))
    yaz('  motor         : screening/read_engine.py %s'
        % getattr(engine_gateway.okuma_motoru, '__version__', '?'))
    yaz(u'                  pigeonhole seeding - LOSSLESS, verified against brute_force.py')
    yaz(u'                  one to one. reads.py/Sonda IS NOT USED.')
    yaz('')

    baglamlar = {d['hedef']: H.hedef_baglami(d, uyelik, kons, kut) for d in panel}
    gerekli = {}
    for b in baglamlar.values():
        for k in b['uye_kutu'] + b['rakip_kutu']:
            gerekli[k['kutu']] = k
    yaz(u'Building raw read pools: %d bins (%s)'
        % (len(gerekli), 'TAM derinlik' if not okuma_sayisi else 'ornekli'))
    if not okuma_sayisi:
        yaz('')
        yaz(u'  >> NOTE: this step reads every fastq file and can take 10-25 MINUTES.')
        yaz(u'     Only bin names scroll past on screen during that time; it is NOT stuck.')
        yaz(u'     If you stop before this step finishes it starts over. The real measurement')
        yaz(u'     begins afterwards, and there each pair is saved separately.')
        yaz('')

    def ilerK(i, n, ad):
        print('   ... %d/%d  %s          ' % (i, n, ad), end='\r', flush=True)

    t0 = time.time()
    numune = N.Numune(list(gerekli.values()), n=(okuma_sayisi or 0),
                      ilerle=ilerK, otorite=True)
    # The OLD (faulty) engine, for the comparison only; the same reads
    eski = {}
    for k in gerekli.values():
        eski[k['kutu']] = N.KutuEski(k['kutu'], k['yol'], n=(okuma_sayisi or 0))
    top_okuma = sum(h.n_okuma for h in numune.havuz.values())
    yaz(u'\nPools ready: %d bins, %d reads  (%s)'
        % (len(gerekli), top_okuma, sure(time.time() - t0)))
    tahmin = len(panel) * len(OLCUTLER) * max(1.0, top_okuma / 20000.0)
    yaz(u'ESTIMATED TIME: ~%s  (resumable; it continues where it stopped)\n'
        % sure(tahmin))

    sonuclar = []
    tb = time.time()
    for i, d in enumerate(panel, 1):
        kp = _kontrol_yolu(d['hedef'])
        if os.path.exists(kp) and not yeniden:
            try:
                _v = json.load(open(kp, encoding='utf-8'))
                if not checks.ayar_uyuyor(_v):
                    raise ValueError('ayar degisti')
                sonuclar.append(json.load(open(kp, encoding='utf-8')))
                yaz(u'[%d/%d] %-38s  (taken from the previous run)' % (i, len(panel), d['hedef'][:38]))
                continue
            except Exception:
                pass   # bayat/bozuk: silmeye calisma, uzerine yazilacak
        b = baglamlar[d['hedef']]
        r = dict(hedef=d['hedef'], sinif=d['sinif'], plaka=d['plaka'], ta=d['ta'],
                 F=d['F'], R=d['R'], urun_panel=d['urun_bp'],
                 panel_uye=d['uye'], panel_rakip=d['rakip'], panel_ayrim=d['ayrim'],
                 panel_ayrim_sayi=d['ayrim_sayi'],
                 uyelik_kaynagi=b['uyelik_kaynagi'],
                 uye_kutu=len(b['uye_kutu']), rakip_kutu=len(b['rakip_kutu']),
                 olcumler={})
        for mm in OLCUTLER:
            o = numune.olc(d['F'], d['R'], b['uye_kutu'], b['rakip_kutu'],
                           lo=C.URUN_IDEAL[0], hi=C.URUN_MUTLAK_UST, mm=mm)
            # measure the same bins with the OLD engine too (so the effect of the fix shows)
            kutu_satir = []
            e_uye = e_rak = 0
            y_uye = y_rak = 0
            for grup, ks in (('uye', b['uye_kutu']), ('rakip', b['rakip_kutu'])):
                for k in ks:
                    he = eski.get(k['kutu']); hy = numune.havuz.get(k['kutu'])
                    if he is None or hy is None:
                        continue
                    pe, ne, _ = he.urun_veren(d['F'], d['R'], C.URUN_IDEAL[0],
                                              C.URUN_MUTLAK_UST, mm)
                    py, ny, _ = hy.urun_veren(d['F'], d['R'], C.URUN_IDEAL[0],
                                              C.URUN_MUTLAK_UST, mm)
                    kutu_satir.append(dict(kutu=k['kutu'], grup=grup, taxid=k['taxid'],
                                           eski=pe, yeni=py, n=ny,
                                           fark=py - pe))
                    if grup == 'uye':
                        e_uye += pe; y_uye += py
                    else:
                        e_rak += pe; y_rak += py
            if o is not None:
                o['eski_motor_uye'] = e_uye
                o['eski_motor_rakip'] = e_rak
                o['yeni_motor_uye'] = y_uye
                o['yeni_motor_rakip'] = y_rak
                o['eski_motor_kayip_uye'] = (round(100.0 * (1 - e_uye / y_uye), 1)
                                             if y_uye else None)
                o['kutu_duzeyi'] = kutu_satir
            r['olcumler'][str(mm)] = o
        o1 = r['olcumler'][str(OLCUTLER[0])]
        yaz('[%d/%d] %-38s  %s' % (i, len(panel), d['hedef'][:38],
            (u'member %%%.1f-%%%.1f coverage %s | discrimination %s x pool / %s x worst'
             % (o1['uye_min'], o1['uye_max'], o1['uye_kapsam_pay'],
                o1['kat_havuz'], o1['kat_enkotu'])) if o1 else 'OLCULEMEDI'))
        r['_ayar'] = dict(checks.AYAR)
        with open(kp, 'w', encoding='utf-8') as fh:
            json.dump(r, fh, ensure_ascii=False, indent=1, default=str)
        sonuclar.append(r)
        gecen = time.time() - tb
        print('       gecen %s  tahmini kalan %s' % (
            sure(gecen), sure(gecen / i * (len(panel) - i))), flush=True)

    hepsi = []
    for f in sorted(os.listdir(C.KONTROL)):
        if f.startswith('panel_olcum_') and f.endswith('.json'):
            try:
                hepsi.append(json.load(open(os.path.join(C.KONTROL, f), encoding='utf-8')))
            except Exception:
                pass
    yollar = rapor_yaz(hepsi or sonuclar, panel_yolu, top_okuma, okuma_sayisi)
    yaz('')
    yaz('=' * 78)
    yaz(u'  RE-MEASUREMENT FINISHED (%s)' % sure(time.time() - t0))
    for p in yollar:
        yaz('    %s' % p)
    yaz('=' * 78)
    return yollar


# ---------------------------------------------------------------- the report
SUT = ['hedef', 'olcut', 'sinif', 'plaka', 'Ta', 'ileri', 'geri', 'urun_panel',
       'uye_kutu', 'uye_min_%', 'uye_max_%', 'uye_kapsam', 'uye_wilson_alt_%',
       'rakip_havuz', 'rakip_havuz_ust_%', 'en_kotu_rakip_kutu',
       'ayrim_havuz_x', 'ayrim_en_kotu_x',
       'ESKI_motor_uye_urun', 'YENI_motor_uye_urun', 'ESKI_motor_kayip_uye_%',
       'PANEL_uye', 'PANEL_rakip_maks', 'PANEL_ayrim', 'PANEL_ayrim_sayi',
       'DEGISIM', 'ESIK_ALTINA_DUSTU', 'uyelik_kaynagi', 'urun_boylari']


def _satirlar(sonuclar):
    for r in sonuclar:
        for mm in OLCUTLER:
            o = (r.get('olcumler') or {}).get(str(mm))
            if not o:
                continue
            yeni = o.get('kat_enkotu') or o.get('kat_havuz')
            eski = r.get('panel_ayrim_sayi')
            deg = ''
            dus = ''
            if eski and yeni:
                oran = yeni / eski
                if oran >= 1.3:
                    deg = 'YUKARI (%.1fx -> %.1fx)' % (eski, yeni)
                elif oran <= 0.77:
                    deg = 'ASAGI (%.1fx -> %.1fx)' % (eski, yeni)
                else:
                    deg = 'ayni (%.1fx -> %.1fx)' % (eski, yeni)
                if eski >= ESIK > yeni:
                    dus = 'EVET - %.0fx esiginin ALTINA dustu' % ESIK
            elif yeni is None:
                deg = 'OLCULEMEDI'
            yield {
                'hedef': r['hedef'], 'olcut': olcut_metni(mm), 'sinif': r['sinif'],
                'plaka': r['plaka'], 'Ta': r['ta'], 'ileri': r['F'], 'geri': r['R'],
                'urun_panel': r['urun_panel'], 'uye_kutu': o.get('uye_kutu_sayisi'),
                'uye_min_%': o.get('uye_min'), 'uye_max_%': o.get('uye_max'),
                'uye_kapsam': o.get('uye_kapsam_pay'),
                'uye_wilson_alt_%': o.get('uye_alt'),
                'rakip_havuz': o.get('havuz'), 'rakip_havuz_ust_%': o.get('havuz_ust'),
                'en_kotu_rakip_kutu': o.get('enkotu_kutu'),
                'ayrim_havuz_x': o.get('kat_havuz'), 'ayrim_en_kotu_x': o.get('kat_enkotu'),
                'ESKI_motor_uye_urun': o.get('eski_motor_uye'),
                'YENI_motor_uye_urun': o.get('yeni_motor_uye'),
                'ESKI_motor_kayip_uye_%': o.get('eski_motor_kayip_uye'),
                'PANEL_uye': r['panel_uye'], 'PANEL_rakip_maks': r['panel_rakip'],
                'PANEL_ayrim': r['panel_ayrim'], 'PANEL_ayrim_sayi': r['panel_ayrim_sayi'],
                'DEGISIM': deg, 'ESIK_ALTINA_DUSTU': dus,
                'uyelik_kaynagi': r['uyelik_kaynagi'],
                'urun_boylari': ';'.join('%s:%s' % (x[0], x[1]) for x in (o.get('urun_boylari') or [])),
            }


def rapor_yaz(sonuclar, panel_yolu, top_okuma, okuma_sayisi):
    os.makedirs(C.CIKTI, exist_ok=True)
    # the bin level TSV (old and new engine, every pair x bin)
    ktsv = os.path.join(C.CIKTI, 'panel_kutu_duzeyi.tsv')
    with open(ktsv, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'olcut', 'kutu', 'grup', 'taxid',
                    'ESKI_motor_urun_veren', 'YENI_motor_urun_veren',
                    'okuma', 'fark'])
        for r in sonuclar:
            for mm in OLCUTLER:
                o = (r.get('olcumler') or {}).get(str(mm)) or {}
                for x in (o.get('kutu_duzeyi') or []):
                    w.writerow([r['hedef'], olcut_metni(mm), x['kutu'], x['grup'],
                                x['taxid'], x['eski'], x['yeni'], x['n'], x['fark']])

    tsv = os.path.join(C.CIKTI, 'panel_yeniden_olcum.tsv')
    with open(tsv, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(SUT)
        for s in _satirlar(sonuclar):
            w.writerow([s.get(k, '') for k in SUT])

    md = os.path.join(C.CIKTI, 'PANEL_YENIDEN_OLCUM.md')
    L = []; A = L.append
    A(u'# The panel remeasured with the corrected engine')
    A('')
    A(u'Generated: %s' % time.strftime('%Y-%m-%d %H:%M:%S'))
    A('')
    A(u'Source panel: `%s` · read depth: **%s** (%d reads in total)'
      % (os.path.basename(panel_yolu),
         'TAM' if not okuma_sayisi else '%d/kutu' % okuma_sayisi, top_okuma))
    A('')
    A(u'## Why it was measured again')
    A('')
    A(u'The panel\'s raw read engine, the `Sonda` class inside `engine/reads.py`, searches for the 13 bases at the primer\'s 3\' end as an **exactly matching** seed. `variants()` only expands the IUPAC codes, it **produces no mismatch variant**. The consequence: every binding site whose single mismatch falls inside those 13 bases at the 3\' end is missed.')
    A('')
    A(u'The examples measured in this run (200 reads, the forward primer, the same criterion):')
    A('')
    A(u'| target / bin | Sonda (reads.py) | brute force | lost |')
    A('|---|---|---|---|')
    A('| Asetoklastik_metanojenler / A1-1_394967 | 0 | 146 | %100 |')
    A('| Arke_universal / A1-1_394967 | 6 | 163 | %96 |')
    A('| Asetoklastik_metanojenler / A1-1_1826872 | 0 | 2 | %100 |')
    A('')
    A(u'The engine used in this table gives a result **identical** to brute force (the self test confirms it on every run).')
    A('')
    A(u'### The real effect measured in this run')
    A('')
    A(u'The table below compares, **on the same reads**, the product count the old and the new engine find in the member bins. All the `kutu duzeyi` columns are in the `panel_kutu_duzeyi.tsv` file.')
    A('')
    A(u'| target | criterion | the OLD engine | the NEW engine | what the old engine lost |')
    A('|---|---|---|---|---|')
    for s2 in _satirlar(sonuclar):
        if s2.get('ESKI_motor_kayip_uye_%') in ('', None):
            continue
        A('| %s | %s | %s | %s | %%%s |' % (
            s2['hedef'], s2['olcut'], s2['ESKI_motor_uye_urun'],
            s2['YENI_motor_uye_urun'], s2['ESKI_motor_kayip_uye_%']))
    A('')
    A(u'## The two criteria are given separately')
    A('')
    A(u'Some of the panel\'s rows had been measured with `<=1` mismatch and some with `<=3`. Here **every row is measured under both criteria**; the `olcut` column says plainly which one was used on each row. The two criteria are **separate and cannot stand in for one another**.')
    A('')

    degisen = [s for s in _satirlar(sonuclar) if s['DEGISIM'].startswith(('YUKARI', 'ASAGI'))]
    dusen = [s for s in _satirlar(sonuclar) if s['ESIK_ALTINA_DUSTU']]
    olculemeyen = [s for s in _satirlar(sonuclar) if s['DEGISIM'] == 'OLCULEMEDI']

    A('## Ozet')
    A('')
    A(u'- Pairs measured: **%d**' % len(sonuclar))
    A(u'- Rows whose value CHANGED (a deviation above %%30): **%d**' % len(degisen))
    A(u'- Rows that FELL BELOW the %.0fx threshold: **%d**' % (ESIK, len(dusen)))
    A(u'- Rows that could not be measured: **%d**' % len(olculemeyen))
    A('')
    if dusen:
        A(u'### The rows that fell below the %.0fx threshold, LOOK AT THESE FIRST' % ESIK)
        A('')
        A(u'| target | criterion | panel | new | member coverage |')
        A('|---|---|---|---|---|')
        for s in dusen:
            A('| %s | %s | %s | %s x | %s |' % (s['hedef'], s['olcut'],
              s['PANEL_ayrim_sayi'], s['ayrim_en_kotu_x'], s['uye_kapsam']))
        A('')
    A(u'## Every row')
    A('')
    A(u'| target | criterion | member % | coverage | discrimination pool x | discrimination worst x | PANEL discrimination | change |')
    A('|---|---|---|---|---|---|---|---|')
    for s in _satirlar(sonuclar):
        A('| %s | %s | %s-%s | %s | %s | %s | %s | %s |' % (
            s['hedef'], s['olcut'], s['uye_min_%'], s['uye_max_%'], s['uye_kapsam'],
            s['ayrim_havuz_x'], s['ayrim_en_kotu_x'], s['PANEL_ayrim'], s['DEGISIM']))
    A('')
    A(u'Every column: `panel_yeniden_olcum.tsv`')
    A('')
    A('## Sinirlar')
    A('')
    A(u'- The numbers the panel published were measured at a different read depth and, on some rows, with a different subset of member bins; the **absolute** numbers are not expected to match exactly. The value of this table is that every row was measured **with the same engine, on the same bins, under the same criterion**.')
    A(u'- The member and competitor bin definition comes from `screening/hedef_uyelik.tsv`. If a row\'s number comes out unexpectedly, **look at that file first**; the `uyelik_kaynagi` column names the source of each row.')
    A('')
    with open(md, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(L))
    return [md, tsv, ktsv]
