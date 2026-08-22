# -*- coding: utf-8 -*-
"""THE FULL SEARCH - the main flow.

Usage (normally verification/full_chain.py calls it):
    python3 -m screening --mode tam
    python3 -m screening --mode sorunlu
    python3 -m screening --mode devam
    python3 -m screening --selftest          (a self test; it measures nothing)

"""
# -------------------------------------------------------------------------
# __main__.py - the package's single entry point; it reads the command line option
#               and runs the matching stage (search, panel measurement, membership
#               audit, consensus generation, summary, all of them) in order.
#
# INPUT  : every path and constant in config.py; the panel TSV through
#          hedefler.panel_oku(); the canonical consensuses through
#          hedefler.konsensusler(); the bin files under "fastq files" through
#          hedefler.kutular(); the per bin raw read pools through numune.Numune;
#          the SILVA/UNITE pools through reference.py; the previous run's
#          checkpoints through checks.py. It takes its arguments with argparse
#          (--mode, --target, --reads, --light, --rerun, --full-depth, --selftest,
#          --skip-tests).
# OUTPUT : the kontrol/hedef_*.json checkpoints under SCREENING_RESULT/;
#          adaylar.tsv, parametre_izgarasi.tsv and KAPSAMLI_ARAMA_RAPORU.md
#          through report.uret(). main() returns the process exit code
#          (0 = success, 2 = a gate or self test failed).
# CALLED BY: the "python3 -m screening --mode <MOD>" line inside
#          verification/full_chain.py. The keys: 1 (--mode tam), 2 (--mode sorunlu),
#          3 (--mode devam), 4 (--mode panel-olc --full-depth), 5 (--mode uyelik),
#          6 (--mode konsensus), 7 (a single target; sorunlu / panel-olc / uyelik
#          plus --target, depending on the choice), 8 (--selftest), 9 (--mode hepsi),
#          S (--mode ozet). The keys P, K, D, T, I, G, E, R, U, H and W/X/Y/Z run
#          separate scripts under verification/, protocol/, engine/ and tools/;
#          they do not call this package directly.
# -------------------------------------------------------------------------
import os, sys, time, argparse, traceback

from . import config as C
from . import engine_gateway, geometry as G, targets as H, generator as U, sample as N
from . import reference as REF, global_scan as KT, checks, report

BASLANGIC = time.time()


def yaz(s=''):
    print(s, flush=True)


def sure(sn):
    sn = int(max(sn, 0))
    if sn < 60:
        return '%d sn' % sn
    if sn < 3600:
        return '%d dk %02d sn' % (sn // 60, sn % 60)
    return '%d sa %02d dk' % (sn // 3600, (sn % 3600) // 60)


def cizgi(ch='-'):
    yaz(ch * 78)


# ---------------------------------------------------------------- one target
# It runs the WHOLE search funnel end to end for a single panel target. The order of
# the stages is not arbitrary but arranged by cost: what is cheap and eliminating
# first (the geometry), then what is expensive (in-silico PCR on the raw reads), and
# last what is most expensive (a global scan of roughly 500 thousand records). That
# way only the few candidates that passed the earlier filters enter the expensive step.
#
# Stage [0] measures A BASELINE for the panel's existing pair. Without that baseline
# the sentence "the candidate is better" cannot be measured: comparing a candidate
# against the existing pair is meaningless unless both are measured with THE SAME
# engine, on THE SAME bins, under THE SAME criterion.
def hedefi_isle(satir, baglam, numune, sira, toplam, hafif=False):
    ad = satir['hedef']
    t0 = time.time()
    cizgi('=')
    yaz(u'[%d/%d] TARGET: %s' % (sira, toplam, ad))
    yaz('       gerekce: ' + '; '.join(g[1] for g in satir['gerekceler'])[:200])
    yaz(u'       pair in the panel: %s / %s  (%d bp)' % (satir['F'], satir['R'], satir['urun_bp']))
    cizgi('=')

    om = baglam['omurga']
    if om is None:
        yaz(u'  SKIPPED: no member consensus found for this target (membership source: %s).'
            % baglam['uyelik_kaynagi'])
        return dict(hedef=ad, durum='ATLANDI - uye konsensusu yok',
                    uyelik_kaynagi=baglam['uyelik_kaynagi'])

    yaz(u'  backbone consensus : %s (%d bp)' % (om['kutu'], len(om['dizi'])))
    yaz(u'  member bins / competitors : %d / %d      member consensus / competitor: %d / %d'
        % (len(baglam['uye_kutu']), len(baglam['rakip_kutu']),
           len(baglam['uye_kons']), len(baglam['rakip_kons'])))

    # ---------------- THE BASELINE: the panel's existing pair, measured with THE SAME engine
    yaz(u'\n  [0] Measuring the existing panel pair with the same engine (comparison baseline)')
    taban = numune.olc(satir['F'], satir['R'], baglam['uye_kutu'], baglam['rakip_kutu'],
                       lo=C.URUN_IDEAL[0], hi=C.URUN_MUTLAK_UST)
    if taban:
        yaz(u'      existing pair: member %%%.1f-%%%.1f | coverage %s bins (>=%d%%) | pool %s'
            % (taban['uye_min'], taban['uye_max'], taban['uye_kapsam_pay'],
               int(100 * C.KAPSAM_ESIGI), taban['havuz']))
        yaz(u'                   discrimination %s x (pool) / %s x (worst bin)   |   covered bins only: %s x / %s x'
            % (taban['kat_havuz'], taban['kat_enkotu'],
               taban['kat_havuz_kapsayan'], taban['kat_enkotu_kapsayan']))
        yaz(u'      -> A candidate only counts as "better" if it beats THESE numbers.')
        uyari = uyelik_uyarisi(satir, taban)
        if uyari:
            yaz('')
            for satir_u in uyari:
                yaz('      ' + satir_u)
    else:
        yaz(u'      the existing pair could not be measured (no member bin).')

    # ---------------- ASAMA A: pencereler + geometri
    yaz(u'\n  [A] Window generation and geometry measurement (18-25 bp, every position, both orientations)')
    ta = time.time()

    def ilerA(n):
        el = time.time() - ta
        print('      ... %d pencere olculdu (%s)' % (n, sure(el)), end='\r', flush=True)

    ad_p = U.aday_primerler(om['dizi'], ilerle=ilerA)
    yaz('      taranan pencere : %d' % ad_p['taranan_pencere'])
    yaz(u'      passing the invariant rules: forward %d, reverse %d   (%s)'
        % (len(ad_p['F']), len(ad_p['R']), sure(time.time() - ta)))

    if not ad_p['F'] or not ad_p['R']:
        return dict(hedef=ad, durum='COZUM YOK - hicbir pencere degismez kurallari gecmedi',
                    pencere=ad_p['taranan_pencere'])

    # ---------------- ASAMA B: cift kurma + izgara
    yaz(u'\n  [B] Compatible primer combinations (product %d-%d bp) and the 144-cell parameter grid'
        % (C.URUN_IDEAL[0], C.URUN_MUTLAK_UST))
    tb = time.time()

    def ilerB(n):
        print(u'      ... %d pairs counted (%s)' % (n, sure(time.time() - tb)), end='\r', flush=True)

    top = U.tara_ve_topla(ad_p, hucre_basina=8, ilerle=ilerB)
    cl = top['temsilciler']
    yaz(u'      pairs matching the rules : %d  (all counted, no upper limit)   (%s)'
        % (top['toplam'], sure(time.time() - tb)))
    if top['toplam'] == 0:
        return dict(hedef=ad, durum='COZUM YOK - urun boyu penceresinde hic cift yok',
                    pencere=ad_p['taranan_pencere'],
                    sayilar=dict(pencere=ad_p['taranan_pencere'], ileri=len(ad_p['F']),
                                 geri=len(ad_p['R']), cift=0))

    izgara = U.izgara_tablosu_sayactan(top['sayac'])
    ornekler = {}
    for c in cl:
        s_, adh = U.hucre_etiketle(c)
        ornekler.setdefault(adh, '%s / %s (%d bp)' % (c['F'], c['R'], c['urun']))
    for x in izgara:
        x['ornek'] = ornekler.get(x['ad'], '')
    yaz(u'      grid: at least one candidate in %d of %d cells'
        % (sum(1 for x in izgara if x['hayatta']), len(izgara)))
    for x in izgara[:4]:
        yaz(u'        %-62s %8d candidates' % (x['ad'], x['hayatta']))
    yaz(u'      representative candidates (one sample per grid cell): %d' % len(cl))

    # ---------------- ASAMA C: numune taramasi
    secili = sec_ornekle(cl, C.HUNI['numuneye_giden'])
    yaz(u'\n  [C] In-silico PCR against the raw sample reads  (%d candidates x %d bins)'
        % (len(secili), len(baglam['uye_kutu']) + len(baglam['rakip_kutu'])))
    tc = time.time()
    olculen = numunede_olc(secili, numune, baglam, tc, yaz)
    yaz(u'\n      candidates with a finished in-sample measurement: %d   (%s)' % (len(olculen), sure(time.time() - tc)))
    olculen.sort(key=puan)

    # ---------------- ASAMA B2: ARMS - EN IYI adaylar uzerinde
    yaz(u'\n  [B2] ARMS variants (discriminating 3\' terminal base plus a deliberate mismatch at -2/-3)')
    yaz(u'       Tried on the candidates that measured best in the sample. The variants go')
    yaz(u'       through the SAME measurement, so the gain is genuinely measured.')
    arms = []
    if baglam['rakip_kons'] and olculen:
        uye_d = [k['dizi'] for k in baglam['uye_kons']][:4]
        rak_d = [k['dizi'] for k in baglam['rakip_kons']][:10]
        tab = time.time()
        for c in olculen[:C.HUNI['arms_taban']]:
            for yon, p in (('F', c['F']), ('R', c['R'])):
                ok, umm, rmm = U.ayirt_edici_mi(p, uye_d, rak_d)
                if not ok:
                    continue
                for v, et in U.arms_varyantlari(p):
                    mv = G.olc(v)
                    if not G.sabit_gecti(mv):
                        continue
                    yeni_c = dict(c)
                    yeni_c.pop('numune', None); yeni_c.pop('pm', None); yeni_c.pop('um', None)
                    if yon == 'F':
                        yeni_c['F'], yeni_c['mF'] = v, mv
                    else:
                        yeni_c['R'], yeni_c['mR'] = v, mv
                    yeni_c['arms'] = '%s %s (uye_mm=%s rakip_mm=%s)' % (yon, et, umm, rmm)
                    yeni_c['arms_taban'] = '%s / %s' % (c['F'], c['R'])
                    arms.append(yeni_c)
        yaz(u'      ARMS variants produced: %d  (%s)' % (len(arms), sure(time.time() - tab)))
        if arms:
            arms = arms[:C.HUNI['arms_ust']]
            ta2 = time.time()
            olculen_arms = numunede_olc(arms, numune, baglam, ta2, yaz)
            yaz(u'\n      ARMS variants measured in the sample: %d   (%s)'
                % (len(olculen_arms), sure(time.time() - ta2)))
            olculen += olculen_arms
            olculen.sort(key=puan)
    else:
        yaz(u'      skipped (light mode, no competitor consensus, or no measured candidate)')
    yaz(u'      NOTE: a deliberate mismatch is NOT a DEGENERATE BASE (one defined base, one oligo);')
    yaz(u'           but it does not match the template exactly, which is a separate agenda item.')

    en_iyi = olculen[:C.HUNI['referansa_giden']]

    # ---------------- STAGE D: the geometry detail plus reference coverage
    yaz(u'\n  [D] Pair geometry (hairpin/dimer dG at 60 C) plus reference coverage and competitor discrimination')
    taxad = H.taxid_adlari()
    cinsler, rakip_cins = REF.uye_ve_rakip_anahtar(baglam, taxad)
    uye_havuz = REF.havuz_cikar(cinsler, baglam['siniflar'][0],
                                'uye_' + ad) if cinsler and not hafif else []
    rak_havuz = REF.havuz_cikar(rakip_cins, baglam['siniflar'][0],
                                'rakip_' + ad) if rakip_cins and not hafif else []
    yaz(u'      reference pool: %d member sequences (%s) | %d competitor sequences (%s)'
        % (len(uye_havuz), ', '.join(cinsler[:4]), len(rak_havuz), ', '.join(rakip_cins[:4])))

    for i, c in enumerate(en_iyi, 1):
        # THE SECOND CRITERION: some rows of the panel were measured at <=3; the best
        # candidates are given under both criteria so that they can be compared.
        try:
            c['numune_mm3'] = numune.olc(c['F'], c['R'], baglam['uye_kutu'],
                                         baglam['rakip_kutu'], lo=C.URUN_IDEAL[0],
                                         hi=C.URUN_MUTLAK_UST, mm=C.KURESEL_MAX_MM
                                         if False else 3)
        except Exception:
            c['numune_mm3'] = {}
        c['cift'] = G.cift_olc(c['mF'], c['mR'], c['urun'])
        c['sikilik'], c['izgara_hucresi'] = U.hucre_etiketle(c)
        if uye_havuz:
            c['ref_uye'] = REF.kapsam(uye_havuz, c['F'], c['R'],
                                      C.URUN_IDEAL[0], C.URUN_MUTLAK_UST)
        if rak_havuz:
            c['ref_rakip'] = REF.kapsam(rak_havuz, c['F'], c['R'],
                                        C.URUN_IDEAL[0], C.URUN_MUTLAK_UST)
        if i % 10 == 0:
            print('      ... %d/%d' % (i, len(en_iyi)), end='\r', flush=True)

    gecen = [c for c in en_iyi if c['cift']['cift_gecti']]
    yaz(u'      passing pair structure (dTm/heterodimer/dG): %d/%d' % (len(gecen), len(en_iyi)))

    # ---------------- ASAMA E: kuresel ozgulluk
    son = gecen[:C.HUNI['kusele_giden']] if gecen else en_iyi[:3]
    if hafif or not son:
        yaz(u'\n  [E] Global specificity SKIPPED (light mode, or no candidate)')
    else:
        yaz(u'\n  [E] GLOBAL SPECIFICITY - the most expensive step, %d candidates, ~500k records' % len(son))
        te = time.time()
        durum_yolu = os.path.join(C.KONTROL, 'kuresel_%s.pkl'
                                  % ''.join(ch if ch.isalnum() else '_' for ch in ad))

        def ilerE(parca, kayit, gecen_sn):
            print(u'      ... chunk %d, %d records, elapsed %s          '
                  % (parca, kayit, sure(gecen_sn)), end='\r', flush=True)

        try:
            kr = KT.tara([dict(ad='a%d' % i, F=c['F'], R=c['R'],
                               lo=C.URUN_IDEAL[0], hi=C.URUN_MUTLAK_UST)
                          for i, c in enumerate(son)],
                         durum_yolu=durum_yolu, ilerle=ilerE)
            for i, c in enumerate(son):
                c['kuresel'] = kr.get('a%d' % i, {})
            yaz(u'\n      global scan finished (%s)' % sure(time.time() - te))
        except Exception as e:
            yaz('\n      GLOBAL SCAN ERROR (the other results are still valid): %s' % e)

    sonuc = dict(
        hedef=ad, durum='TAMAMLANDI',
        gerekceler=[g[1] for g in satir['gerekceler']],
        etiketler=satir['etiketler'],
        panel=dict(F=satir['F'], R=satir['R'], urun=satir['urun_bp'],
                   plaka=satir['plaka'], ta=satir['ta'], ayrim=satir['ayrim'],
                   geo=satir['geo'], jel=satir['jel']),
        panel_olcum=taban,
        uyelik_uyarisi=uyelik_uyarisi(satir, taban),
        omurga=dict(kutu=om['kutu'], uzunluk=len(om['dizi']), yol=om['yol']),
        uyelik_kaynagi=baglam['uyelik_kaynagi'],
        uye_tax=baglam['uye_tax'],
        sayilar=dict(pencere=ad_p['taranan_pencere'], ileri=len(ad_p['F']),
                     geri=len(ad_p['R']), cift=top['toplam'],
                     temsilci=len(cl), arms=len(arms),
                     numune_olculen=len(olculen), cift_yapisi_gecen=len(gecen)),
        izgara=izgara,
        adaylar=[report.aday_ozet(c) for c in en_iyi],
        sure_sn=round(time.time() - t0, 1),
    )
    yaz(u'\n  DONE: %s  (%s)' % (ad, sure(time.time() - t0)))
    return sonuc


def puan(c):
    """The ranking criterion: FIRST the member coverage (how many member bins amplify),
        THEN the discrimination.

        A single member bin coming out empty makes uye_alt 0 and levels every candidate;
        the coverage axis prevents that.

    """
    # WHY COVERAGE COMES FIRST: the discrimination ratio is computed over the smallest
    # of the member bins' Wilson LOWER bounds. If a single member bin gives no product
    # at all, that lower bound is 0 and so is the ratio, so a good and a bad candidate
    # get the same score and the ranking collapses. Coverage (how many member bins give
    # >=20% product) is not affected by that collapse, which is why it is the first key.
    # On universal targets, where the competitor set is empty, the discrimination ratio
    # is undefined anyway; there the only axis carrying the comparison is coverage. The
    # next two keys are the discrimination ratio computed first over "the covered bins
    # only" and then over all the member bins.
    n = c['numune']
    return (-n.get('uye_kapsam', 0),
            -(n.get('kat_enkotu_kapsayan') or n.get('kat_havuz_kapsayan') or 0),
            -(n.get('kat_enkotu') or n.get('kat_havuz') or 0))


def uyelik_uyarisi(satir, taban):
    """Warns when the measured baseline deviates far from the panel's PUBLISHED value.

        The most common cause is the membership definition (the member and competitor bin
        list). Warning loudly is better than measuring the wrong thing silently.

    """
    # THE MEMBERSHIP IS NOT ADOPTED UNCONDITIONALLY. This function only REPORTS the
    # deviation; it does not change the member or competitor bin list on its own and it
    # does not write to target_membership.tsv. The principle: an absence of evidence is not
    # evidence. For a bin to change place, positive measured evidence is needed, and "the
    # number came out different from what was expected" is not such evidence. So the only
    # thing done here is to point the user at which row of the file to look at.
    # The band is 0.34x to 3.0x: if ANY of the measured values falls inside that range
    # around the value the panel published, no warning is printed. The band is kept wide
    # because the panel rows were measured at different read depths, and the width of the
    # Wilson interval depends on depth.
    if not taban:
        return []
    p = satir.get('ayrim_sayi')
    if p is None or 'kapsam' in (satir.get('ayrim') or '').lower() \
            or 'kutu' in (satir.get('ayrim') or '').lower():
        return []
    olculen = [x for x in (taban.get('kat_enkotu'), taban.get('kat_havuz'),
                           taban.get('kat_enkotu_kapsayan'),
                           taban.get('kat_havuz_kapsayan')) if x]
    if not olculen:
        return [' !! WARNING: the discrimination ratio for the current pair could NOT be measured at all, because the member bins',
                '   hicbiri urun vermiyor. Uyelik tanimi yanlis olabilir:',
                '   screening/target_membership.tsv -> satir "%s"' % satir['hedef']]
    if any(0.34 * p <= o <= 3.0 * p for o in olculen):
        return []
    return [
        '!! UYARI: panelin yayimladigi ayrim %.1fx, bu koşuda olculen %s.' % (
            p, ' / '.join('%.1fx' % o for o in olculen)),
        '   The deviation is large. The most likely cause is the MEMBERSHIP DEFINITION (which bin is a member, which is a competitor).',
        '   Once su dosyaya bakin: screening/target_membership.tsv  ->  satir "%s"' % satir['hedef'],
        '   (If the read count is low the Wilson interval widens, and part of the deviation comes from that.)',
    ]


def numunede_olc(adaylar, numune, baglam, t0, yaz):
    """Aday listesini ham okumalarda olc, ilerlemeyi ve tahmini kalan sureyi bas."""
    out = []
    for i, c in enumerate(adaylar, 1):
        r = numune.olc(c['F'], c['R'], baglam['uye_kutu'], baglam['rakip_kutu'],
                       lo=C.URUN_IDEAL[0], hi=C.URUN_MUTLAK_UST)
        if r:
            c = dict(c); c['numune'] = r
            out.append(c)
        if i % 5 == 0 or i == len(adaylar):
            gecen = time.time() - t0
            kalan = gecen / i * (len(adaylar) - i)
            print(u'      ... candidate %d/%d  elapsed %s  estimated remaining %s        '
                  % (i, len(adaylar), sure(gecen), sure(kalan)), end='\r', flush=True)
    return out


def sec_ornekle(adaylar, ust):
    """Take a representative from every grid cell, so that the question 'at which setting
        is there an answer' can be answered for all the settings.

    """
    from collections import defaultdict
    kova = defaultdict(list)
    for c in adaylar:
        s, ad = U.hucre_etiketle(c)
        kova[(s, ad)].append(c)
    for k in kova:
        kova[k].sort(key=lambda c: (abs(c['urun'] - 105), abs(c['mF']['tm'] - c['mR']['tm'])))
    out, i = [], 0
    anahtarlar = sorted(kova, key=lambda k: k[0])
    while len(out) < ust:
        eklendi = False
        for k in anahtarlar:
            if i < len(kova[k]):
                out.append(kova[k][i]); eklendi = True
                if len(out) >= ust:
                    break
        if not eklendi:
            break
        i += 1
    return out


# ---------------------------------------------------------------- ana


def aramayi_kos(a, yaz, sure, cizgi, mod=None):
    """Kapsamli arama akisi. run_all.py de bunu cagirir."""
    # THE ORIENTATION GATE - this stage's first job. Nanopore reads are bidirectional and
    # a consensus may have been produced in reverse. Every in-silico PCR engine in this
    # project (ispcr.amplify, okuma_motoru.Sonda) scans the given sequence ONLY on the
    # plus strand; on a reversed consensus neither the forward primer nor the complement
    # of the reverse primer is found. The engine raises no error, it quietly says "no
    # product", which means a search running all night produces "no answer on any target"
    # with no warning at all. The measured loss is 100% (orientation_impact_test.py: 117
    # products in the correct orientation, 0 in the reverse). So if the consensuses are
    # not canonical the stage DOES NOT START; it stops with SystemExit(2).
    if mod:
        a.mod = mod
    from .run_all import yon_kapisi
    _ok, _m = yon_kapisi(yaz, 'kapsamli arama')
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

    sorunlu, panel, panel_yolu = H.sorunlu_hedefler()
    yaz(u'\nPanel file  : %s' % panel_yolu)
    yaz(u'Pairs in the panel: %d     found problematic: %d' % (len(panel), len(sorunlu)))
    for d in sorunlu:
        yaz('  [%-5s] %-36s %s' % (d['etiketler'], d['hedef'][:36],
                                   '; '.join(g[1] for g in d['gerekceler'])[:90]))

    if a.mod == 'tam':
        liste = panel
        for d in liste:
            d.setdefault('gerekceler', H.sorun_gerekceleri(d) or [('-', 'panelin tamami taraniyor')])
            d.setdefault('etiketler', '-')
    else:
        liste = sorunlu
    if a.hedef:
        liste = [d for d in liste if a.hedef.lower() in d['hedef'].lower()]

    if a.mod == 'devam':
        liste = [d for d in liste if not checks.bitti_mi(d['hedef'])]
        yaz(u'\nCONTINUE MODE: finished targets are skipped. Remaining: %d' % len(liste))

    if not liste:
        yaz(u'\nNo target to process. The report is being built from the existing checkpoint files.')
        report.uret(checks.hepsi(), panel, panel_yolu, yaz)
        return 0
    uyelik = H.uyelik_oku(); kons = H.konsensusler(); kut = H.kutular()
    baglamlar = {d['hedef']: H.hedef_baglami(d, uyelik, kons, kut) for d in liste}

    gerekli = {}
    for b in baglamlar.values():
        for k in b['uye_kutu'] + b['rakip_kutu']:
            gerekli[k['kutu']] = k
    yaz(u'\nBuilding raw read pools: %d bins x %d reads' % (len(gerekli), a.okuma))

    def ilerK(i, n, ad):
        print('   ... %d/%d  %s          ' % (i, n, ad), end='\r', flush=True)

    tk = time.time()
    numune = N.Numune(list(gerekli.values()), n=a.okuma, ilerle=ilerK)
    yaz(u'\nPools ready (%s)' % sure(time.time() - tk))

    if not a.hafif:
        taxad = H.taxid_adlari()
        istekler = []
        for d in liste:
            b = baglamlar[d['hedef']]
            cins, rakip = REF.uye_ve_rakip_anahtar(b, taxad)
            sf = b['siniflar'][0] if b['siniflar'] else 'B'
            istekler.append(('uye_' + d['hedef'], cins, sf))
            istekler.append(('rakip_' + d['hedef'], rakip, sf))
        yaz(u'\nExtracting reference pools (ONE pass per database)...')

        def ilerR(db, n, bulunan):
            print(u'   ... %s: %d records scanned, %d sequences collected        '
                  % (db, n, bulunan), end='\r', flush=True)

        tr = time.time()
        try:
            REF.toplu_cikar(istekler, ilerle=ilerR)
            yaz(u'\nThe reference pools are ready (%s)' % sure(time.time() - tr))
        except Exception as e:
            yaz(u'\nThe reference pool could not be extracted (%s); the reference step will be skipped.' % e)

    yaz(u'\nESTIMATED TIME: ~%s per target; %d targets -> ~%s'
        % (sure(600 if not a.hafif else 90), len(liste),
           sure((600 if not a.hafif else 90) * len(liste))))
    yaz(u'Resumable: you can close the window and reopen it, then continue with (3).\n')

    for i, d in enumerate(liste, 1):
        try:
            s = hedefi_isle(d, baglamlar[d['hedef']], numune, i, len(liste), hafif=a.hafif)
            checks.yaz(d['hedef'], s)
            report.uret(checks.hepsi(), panel, panel_yolu, lambda *x: None)
        except KeyboardInterrupt:
            yaz(u'\n\nSTOPPED BY THE USER. The finished targets are saved; continue with (3).')
            break
        except Exception:
            yaz(u'\n  ERROR (%s) - this target was skipped, the others continue:' % d['hedef'])
            traceback.print_exc()
            checks.yaz(d['hedef'], dict(hedef=d['hedef'], durum='HATA',
                                         hata=traceback.format_exc()[-2000:]))

    yolar = report.uret(checks.hepsi(), panel, panel_yolu, yaz)
    cizgi('=')
    yaz(u'  TOTAL TIME: %s' % sure(time.time() - BASLANGIC))
    yaz('  RESULTS:')
    for p in yolar:
        yaz('    %s' % p)
    cizgi('=')
    return 0




def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', dest='mod', default='sorunlu',
                    choices=['tam', 'sorunlu', 'devam', 'panel-olc', 'konsensus',
                             'uyelik', 'hepsi', 'ozet'])
    ap.add_argument('--selftest', dest='sina', action='store_true')
    ap.add_argument('--list-targets', dest='hedefleri_listele', action='store_true')
    ap.add_argument('--skip-tests', dest='sinama_atla', action='store_true',
                    help='skip the self-test (development and testing only)')
    ap.add_argument('--rerun', dest='yeniden', action='store_true',
                    help='ignore checkpoints and recompute from scratch')
    ap.add_argument('--light', dest='hafif', action='store_true', help='skip the reference and global scan steps')
    ap.add_argument('--full-depth', dest='tam_derinlik', action='store_true',
                    help='use EVERY read in the bin (same as --reads 0)')
    ap.add_argument('--reads', dest='okuma', type=int, default=C.NUMUNE_OKUMA_SAYISI,
                    help='reads per bin; 0 = EVERY read in the bin (slow but exact)')
    ap.add_argument('--target', dest='hedef', default=None, help='this target only (for testing)')
    ap.add_argument('--candidate-max', dest='aday_ust', type=int, default=None)
    a = ap.parse_args(argv)

    if a.tam_derinlik:
        a.okuma = 0
    if a.aday_ust:
        C.HUNI['numuneye_giden'] = a.aday_ust

    if a.hedefleri_listele:
        try:
            sorunlu, panel, _ = H.sorunlu_hedefler()
            sor = {d['hedef'] for d in sorunlu}
            for d in panel:
                print('   %-46s %s' % (d['hedef'][:46],
                                       '<- sorunlu' if d['hedef'] in sor else ''))
        except Exception as e:
            print(u'The target list could not be read:', e)
        return 0

    checks.ayar_kur(okuma=a.okuma, hafif=bool(a.hafif),
                     aday_ust=C.HUNI['numuneye_giden'])
    checks.hazirla()
    cizgi('=')
    yaz(u'  COMPREHENSIVE PRIMER SEARCH - the PrimerJury panel')
    yaz(u'  start: %s     mode: %s' % (time.strftime('%Y-%m-%d %H:%M:%S'), a.mod))
    yaz(u'  output directory: %s' % C.CIKTI)
    cizgi('=')

    from . import self_test
    if a.mod == 'ozet':
        from . import run_all as HP
        rc = HP.yalniz_ozet(yaz, sure, cizgi)
        yaz(u'  TOTAL TIME: %s' % sure(time.time() - BASLANGIC))
        return rc

    if a.mod == 'hepsi':
        from . import run_all as HP
        rc = HP.calistir(yaz, sure, cizgi, a)
        yaz(u'  TOTAL TIME: %s' % sure(time.time() - BASLANGIC))
        return rc

    if a.mod == 'uyelik':
        from . import membership_check
        if not a.sinama_atla and not self_test.calistir(yaz):
            yaz(u'\nTHE SELF TEST FAILED - the audit was not started.')
            return 2
        membership_check.calistir(yaz, sure,
                                 okuma_sayisi=(a.okuma or C.NUMUNE_OKUMA_SAYISI),
                                 yalniz=a.hedef, yeniden=a.yeniden)
        yaz(u'  TOTAL TIME: %s' % sure(time.time() - BASLANGIC))
        return 0

    if a.mod == 'konsensus':
        from . import build_consensus
        build_consensus.calistir(yaz, sure, yalniz=a.hedef, yeniden=a.yeniden)
        yaz(u'  TOTAL TIME: %s' % sure(time.time() - BASLANGIC))
        return 0

    if not a.sinama_atla and not self_test.calistir(yaz):
        yaz(u'\nTHE SELF TEST FAILED, the search was not started.')
        return 2
    if a.sina:
        yaz(u'\nOnly the test was asked for, no search was done.')
        return 0

    if a.mod == 'panel-olc':
        from . import panel_measurement
        panel_measurement.calistir(yaz, sure, okuma_sayisi=a.okuma if a.okuma else 0,
                             yalniz=a.hedef, yeniden=a.yeniden)
        yaz(u'  TOTAL TIME: %s' % sure(time.time() - BASLANGIC))
        return 0

    return aramayi_kos(a, yaz, sure, cizgi)

if __name__ == '__main__':
    sys.exit(main())
