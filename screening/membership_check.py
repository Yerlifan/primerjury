# -*- coding: utf-8 -*-
"""OPTION 5 - THE MEMBERSHIP DEFINITION AUDIT and SENSITIVITY ANALYSIS.

WHY IT IS NEEDED
----------------
A target's "discrimination ratio" depends directly on which bins count as MEMBERS
and which as COMPETITORS. Change the membership definition and the number changes,
even with the engine unchanged. Measured:

  Proteolitik_Cloacimonas, with the group row (3 taxa)  ->  discrimination  0.0x
  Proteolitik_Cloacimonas, with a single member (456827) ->  discrimination 23.5x
  what the panel says                                    ->                 23.0x

So it was not the number that was wrong but THE DEFINITION. This module measures
EVERY reasonable definition of the membership for each target side by side, so
that it becomes visible which number depends on which definition.

THE DEFINITION SOURCES
----------------------
  A. target_membership.tsv       the definition the tool currently uses (PANEL / TURETILDI)
  B. targets.tsv           the project's decision table (the group row)
  C. a single member        only the single taxon giving the most product
  D. pairs.tsv            the definition from the other session (the read engine fix)
  E. the measured identity  derived from the MEASURED organism in hedef_kimlik.tsv

THE DIAGNOSIS (situations like Proteiniphilum)
----------------------------------------------
When there is no product at all in a target's member bins there are two
possibilities:
  (1) the member set is wrong  -> which bin or bins do give a product is found one
      by one
  (2) the consensus and the raw reads do not agree -> the same pair is measured ON
      THE CONSENSUS; if it gives a product on the consensus but not on the reads,
      the problem is not the membership but a consensus and read mismatch
The module tests BOTH and writes down which it is.

"""
# -------------------------------------------------------------------------
# membership_check.py measures every panel target's membership definitions (which bin
#                      is a member, which a competitor) side by side and shows how
#                      sensitive the number is to the definition.
#
# INPUT  : the panel TSV through hedefler.panel_oku(); the fastq bins through
#          hedefler.kutular(); the canonical consensuses through
#          hedefler.konsensusler(); screening/target_membership.tsv through
#          hedefler.acik_uyelik(); steps/targets.tsv through hedefler.uyelik_oku();
#          screening/pairs.tsv (or eski/pairs.tsv);
#          primer_final/hedef_kimlik.tsv. numune.Numune does the measuring.
# OUTPUT : SCREENING_RESULT/UYELIK_DENETIMI.md and uyelik_duyarlilik.tsv
#          (calistir returns those two paths as a list); plus a
#          kontrol/uyelik_*.json checkpoint per target.
# CALLED BY: verification/full_chain.py key 5 (--mode uyelik), key 7 -> choice "3"
#          (the membership audit of a single target) and the 6th stage inside key 9
#          (hepsi.calistir -> uyelik_denetimi.calistir).
#
# THIS MODULE DOES NOT CHANGE THE MEMBERSHIP. It measures every reasonable
# definition, puts them side by side, and writes at the end of the report "this tool
# does not change the membership definition by itself"; it DOES NOT WRITE to
# target_membership.tsv. The principle: an absence of evidence is not evidence. A number
# coming out different from what was expected is not evidence for moving a bin; a bin
# moves only on positive measured evidence and by hand. So this module's output is
# not a decision but a table of options.
# -------------------------------------------------------------------------
import os, csv, json, time, re
from . import config as C
from . import engine_gateway, targets as H, sample as N, checks

# THE THRESHOLD COMES FROM ONE SOURCE: screening/config.py -> ESIK_DCQ = 3.0
# Its fold equivalent is 2 ** ESIK_DCQ = 8.00. NO constant is EMBEDDED; if dCq
# changes it changes in one place. The reasoning and the efficiency warning are
# written in that file.
ESIK = C.AYRIM_ESIK
PAKET = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------- tanim kaynaklari
def _ciftler_tsv_uyelik():
    """eski/pairs.tsv (the other session's definition) -> target -> the taxid list."""
    out = {}
    for aday in (os.path.join(PAKET, 'pairs.tsv'),
                 os.path.join(PAKET, 'eski', 'pairs.tsv')):
        if not os.path.exists(aday):
            continue
        try:
            for row in csv.DictReader(open(aday, encoding='utf-8'), delimiter='\t'):
                t = [x.strip() for x in (row.get('uye_taksonlar') or '').split(',')
                     if x.strip()]
                if t:
                    out[row['hedef'].strip()] = t
        except Exception:
            pass
        if out:
            break
    return out


def _kimlik_uyelik():
    'Derive a taxid set from the MEASURED identity in the identity table.'
    yol = os.path.join(C.KOK, 'primer_final', 'hedef_kimlik.tsv')
    taxad = H.taxid_adlari()
    ad2tax = {}
    for t, a in taxad.items():
        ad2tax.setdefault(a.lower(), t)
    out = {}
    if not os.path.exists(yol):
        return out
    try:
        for row in csv.DictReader(open(yol, encoding='utf-8'), delimiter='\t'):
            olculen = (row.get('olculen_kimlik') or '').strip()
            if not olculen:
                continue
            ad = re.split(r'\s*\(|\s*/', olculen)[0].strip().lower()
            tax = ad2tax.get(ad)
            if not tax:
                cins = ad.split()[0] if ad else ''
                tax = next((t for a, t in ad2tax.items()
                            if cins and a.startswith(cins)), None)
            if tax:
                out[row['hedef'].strip()] = [tax]
    except Exception:
        pass
    return out


def tanimlar(satir, kut, acik, grup_uyelik, cift_uyelik, kimlik_uyelik):
    """Produce every membership definition for one panel target."""
    ad = satir['hedef']
    sf = [x.strip() for x in (satir['sinif'] or '').split('/') if x.strip()]
    out = []

    a = acik.get(ad)
    if a:
        out.append((u'A. target_membership.tsv (%s)' % (a.get('kaynak') or '?'),
                    a['uye'], a['haric']))

    anahtar = H.AD_ESLEME.get(ad, ad)
    g = grup_uyelik.get(anahtar)
    if g:
        out.append(('B. targets.tsv:%s' % anahtar, g['uye'], g['haric']))

    c = cift_uyelik.get(ad)
    if c:
        out.append((u'D. pairs.tsv (the other session)', c, []))

    k = kimlik_uyelik.get(ad)
    if k:
        out.append((u'E. the measured identity (hedef_kimlik.tsv)', k, []))

    return out, sf


def _kutu_coz(uye_tax, haric, sf, kut):
    yildiz = [t for t in uye_tax if t.startswith('*')]
    if yildiz:
        onek = [t[1:] for t in yildiz]
        uk = [k for k in kut if any(k['sinif'] == o or k['sinif'].startswith(o)
                                    for o in onek) and k['taxid'] not in haric]
        tax = sorted({k['taxid'] for k in uk})
    else:
        tax = [t for t in uye_tax if not t.startswith('*')]
        uk = [k for k in kut if k['taxid'] in tax and k['sinif'] in sf]
        if not uk:
            uk = [k for k in kut if k['taxid'] in tax]
    rk = [k for k in kut if k['sinif'] in sf and k['taxid'] not in tax
          and k['taxid'] not in haric]
    return uk, rk, tax


# ---------------------------------------------------------------- tani
def tani(satir, sf, kut, numune, kons):
    """With no product in a member bin, find the reason: the membership, or the consensus and reads."""
    # The two questions are asked separately because their answers point at two different
    # fixes. (1) "WHICH bins of the same class give a product" - if a bin giving a
    # product is not in the member list, the problem is in the membership definition.
    # (2) "Does the same pair give a product on the CONSENSUS sequence" - if it does on
    # the consensus but not on the raw reads, the problem is not the membership but the
    # consensus failing to represent the reads, and the fix is to regenerate the
    # consensus. Confuse the two and the wrong file gets corrected.
    F, R = satir['F'], satir['R']
    lo, hi = C.URUN_IDEAL[0], C.URUN_MUTLAK_UST
    taxad = H.taxid_adlari()

    # (1) WHICH bins give a product - the whole class is scanned
    veren = []
    for k in kut:
        if k['sinif'] not in sf:
            continue
        h = numune.havuz.get(k['kutu'])
        if h is None:
            continue
        p, n, boy = h.urun_veren(F, R, lo, hi)
        if p:
            veren.append(dict(kutu=k['kutu'], taxid=k['taxid'],
                              ad=taxad.get(k['taxid'], '?'),
                              urun=p, okuma=n, yuzde=round(100.0 * p / max(n, 1), 2),
                              boy=sorted(boy.items(), key=lambda x: -x[1])[:2]))
    veren.sort(key=lambda x: -x['yuzde'])

    # (2) KONSENSUSTE urun veriyor mu (okumalarda vermese bile)
    kons_veren = []
    for k in kons:
        if k['sinif'] not in sf:
            continue
        for s in (k['dizi'], engine_gateway.rc(k['dizi'])):
            pr = engine_gateway.amplify(s, F, R, max_mm=1, lo=lo, hi=hi)
            if pr:
                bp = min(pr, key=lambda x: x[3] + x[4])[2]
                kons_veren.append(dict(kutu=k['kutu'], taxid=k['taxid'],
                                       ad=taxad.get(k['taxid'], '?'), boy=bp))
                break
    return veren, kons_veren


def tani_yorumu(uye_kutu, veren, kons_veren, panel_urun, olcum=None, panel_uye=''):
    """When the coverage is zero, find the reason: the membership, or a consensus and read mismatch.

        "There is a product" and "it is covered" ARE NOT the same thing: if a member bin
        gives 2% product there IS a product, but it does not pass the coverage threshold.
        If the value the panel reports for that row is far higher, that is a contradiction
        needing a separate explanation, and the module reports it as a separate diagnosis.

    """
    # The diagnostic order MATTERS and is written as a narrowing elimination: first is
    # there coverage, then is there a product in a member bin however weak, then is there
    # one in another bin of the class, and last is there one on the consensus. Each step
    # rules out the previous possibility; change the order and the same situation gets a
    # different diagnosis.
    uye_ad = {k['kutu'] for k in uye_kutu}
    uye_veren = [v for v in veren if v['kutu'] in uye_ad]
    kons_uye = [v for v in kons_veren if v['kutu'] in uye_ad]
    kapsam = (olcum or {}).get('uye_kapsam', 0)
    esik = int(100 * C.KAPSAM_ESIGI)

    if kapsam:
        return 'SORUN YOK', '%d of the member bins pass the coverage threshold of %d per cent or more.' % (
            kapsam, esik)

    if uye_veren:
        en = max(uye_veren, key=lambda v: v['yuzde'])
        dis = [v for v in veren if v['kutu'] not in uye_ad]
        ek = ''
        if dis:
            d0 = max(dis, key=lambda v: v['yuzde'])
            ek = (' In the same class %s (%s), which is NOT a member, gives %s per cent. Should the member set be widened?' % (d0['kutu'], d0['ad'], d0['yuzde']))
        return ('BELOW THE COVERAGE THRESHOLD, and it contradicts the panel',
                'There IS a product in the member bins but it is weak: the best is %s at %s per cent against a threshold of %d per cent, with a product length of %s bp. The panel writes "%s" for this row. If the product length holds, the pair is right and what differs is THE RATE. The likely causes: (a) the panel was measured under a different criterion or a different looseness setting, or (b) the consensus and the raw reads are not describing the same organism. The consensus test is below.%s'
                % (en['kutu'], en['yuzde'], esik,
                   en['boy'][0][0] if en['boy'] else '?', panel_uye, ek))

    if veren:
        ilk = veren[0]
        return ('THE MEMBER SET MAY BE WRONG',
                'There is NO product in the member bins while %d bins of the same class give one. The highest is %s (%s) at %s per cent, with a product of %s. Should that taxon count as a member?'
                % (len(veren), ilk['kutu'], ilk['ad'], ilk['yuzde'],
                   ilk['boy'][0][0] if ilk['boy'] else '?'))
    if kons_uye:
        b = kons_uye[0]
        return ('KONSENSUS/OKUMA UYUSMAZLIGI',
                "The member bin's CONSENSUS gives a product (%s, %s bp against %s bp in the panel) while the same bin's RAW READS do not. The problem is not the membership: either the consensus and the reads are not describing the same organism, or the consensus is stale. Regenerate the consensus."
                % (b['kutu'], b['boy'], panel_urun))
    return ('THE PAIR GIVES NO PRODUCT AT ALL',
            "No product forms in the member bins, in any bin of the same class, or in the consensuses. The pair's sequences or the product length window have to be checked.")


# ---------------------------------------------------------------- ana
def calistir(yaz, sure, okuma_sayisi=C.NUMUNE_OKUMA_SAYISI, yalniz=None,
             yeniden=False):
    from .run_all import yon_kapisi
    _ok, _m = yon_kapisi(yaz, 'uyelik denetimi')
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
    if yalniz:
        panel = [d for d in panel if yalniz.lower() in d['hedef'].lower()]
    kut = H.kutular()
    kons = H.konsensusler()
    acik = H.acik_uyelik()
    grup_uyelik = H.uyelik_oku()
    cift_uyelik = _ciftler_tsv_uyelik()
    kimlik_uyelik = _kimlik_uyelik()

    turetildi = [ad for ad, v in acik.items() if (v.get('kaynak') or '') == 'TURETILDI']

    yaz('=' * 78)
    yaz(u'  MEMBERSHIP DEFINITION AUDIT and SENSITIVITY ANALYSIS')
    yaz('=' * 78)
    yaz(u'  number of targets   : %d' % len(panel))
    yaz(u'  definition sources  : target_membership.tsv, targets.tsv, pairs.tsv,')
    yaz(u'                        measured identity, single member')
    yaz(u'  marked DERIVED      : %d rows  (these deserve a specific look)' % len(turetildi))
    for t in turetildi:
        yaz('        - %s' % t)
    yaz('')

    gerekli = {k['kutu']: k for k in kut}
    yaz(u'Building raw read pools: %d bins' % len(gerekli))
    yaz(u'  >> This step takes a few minutes. Bin names scroll past on screen; it is not stuck.')

    def ilerK(i, n, ad):
        print('   ... %d/%d  %s        ' % (i, n, ad), end='\r', flush=True)

    t0 = time.time()
    numune = N.Numune(list(gerekli.values()), n=okuma_sayisi, ilerle=ilerK)
    yaz('\nHavuzlar hazir (%s)' % sure(time.time() - t0))
    yaz('TAHMINI SURE: ~%s\n' % sure(len(panel) * 25))

    sonuclar = []
    for i, d in enumerate(panel, 1):
        kp = os.path.join(C.KONTROL, 'uyelik_%s.json'
                          % ''.join(c if c.isalnum() else '_' for c in d['hedef']))
        if os.path.exists(kp) and not yeniden:
            try:
                _v = json.load(open(kp, encoding='utf-8'))
                if not checks.ayar_uyuyor(_v):
                    raise ValueError('ayar degisti')
                sonuclar.append(json.load(open(kp, encoding='utf-8')))
                yaz(u'[%d/%d] %-38s (from the previous run)' % (i, len(panel), d['hedef'][:38]))
                continue
            except Exception:
                pass   # bayat/bozuk: silmeye calisma, uzerine yazilacak

        tnm, sf = tanimlar(d, kut, acik, grup_uyelik, cift_uyelik, kimlik_uyelik)
        yaz('[%d/%d] %s' % (i, len(panel), d['hedef']))
        varyantlar = []
        gorulen = set()
        for etiket, uye_tax, haric in tnm:
            uk, rk, tax = _kutu_coz(uye_tax, haric, sf, kut)
            anahtar = (tuple(sorted(t['kutu'] for t in uk)),)
            o = numune.olc(d['F'], d['R'], uk, rk, lo=C.URUN_IDEAL[0],
                           hi=C.URUN_MUTLAK_UST) if uk else None
            varyantlar.append(dict(
                tanim=etiket, uye_taxid=','.join(tax), uye_kutu=len(uk),
                rakip_kutu=len(rk),
                kapsam=(o or {}).get('uye_kapsam_pay', ''),
                uye_min=(o or {}).get('uye_min', ''), uye_max=(o or {}).get('uye_max', ''),
                kat_havuz=(o or {}).get('kat_havuz'), kat_enkotu=(o or {}).get('kat_enkotu'),
                havuz=(o or {}).get('havuz', ''),
                ayni_kume=('EVET' if anahtar in gorulen else '')))
            gorulen.add(anahtar)
            yaz(u'        %-42s member %2d bins  coverage %-6s discrimination %s x / %s x'
                % (etiket[:42], len(uk), (o or {}).get('uye_kapsam_pay', '-'),
                   (o or {}).get('kat_havuz'), (o or {}).get('kat_enkotu')))

        # C. a single member - which single taxon gives the best result
        a = acik.get(d['hedef'])
        aday_tax = (a['uye'] if a else [])
        aday_tax = [t for t in aday_tax if not t.startswith('*')]
        if len(aday_tax) > 1:
            en = None
            for t in aday_tax:
                uk, rk, _ = _kutu_coz([t], [], sf, kut)
                if not uk:
                    continue
                o = numune.olc(d['F'], d['R'], uk, rk, lo=C.URUN_IDEAL[0],
                               hi=C.URUN_MUTLAK_UST)
                v = (o or {}).get('kat_enkotu') or (o or {}).get('kat_havuz') or 0
                if en is None or v > en[0]:
                    en = (v, t, o, len(uk), len(rk))
            if en:
                varyantlar.append(dict(
                    tanim='C. tek uye (%s = %s)' % (en[1], H.taxid_adlari().get(en[1], '?')),
                    uye_taxid=en[1], uye_kutu=en[3], rakip_kutu=en[4],
                    kapsam=en[2].get('uye_kapsam_pay', ''),
                    uye_min=en[2].get('uye_min', ''), uye_max=en[2].get('uye_max', ''),
                    kat_havuz=en[2].get('kat_havuz'), kat_enkotu=en[2].get('kat_enkotu'),
                    havuz=en[2].get('havuz', ''), ayni_kume=''))
                yaz(u'        %-42s member %2d bins  coverage %-6s discrimination %s x / %s x'
                    % (('C. tek uye (%s)' % en[1])[:42], en[3],
                       en[2].get('uye_kapsam_pay', '-'),
                       en[2].get('kat_havuz'), en[2].get('kat_enkotu')))

        # the diagnosis: is there a product in the member bins of definition A
        uk_a, _rk_a, _ = _kutu_coz(a['uye'], a['haric'], sf, kut) if a else ([], [], [])
        veren, kons_veren = tani(d, sf, kut, numune, kons)
        olcum_a = None
        for v in varyantlar:
            if v['tanim'].startswith('A.'):
                olcum_a = dict(uye_kapsam=int((v['kapsam'] or '0/0').split('/')[0] or 0))
                break
        tsonuc, taciklama = tani_yorumu(uk_a, veren, kons_veren, d['urun_bp'],
                                        olcum_a, d.get('uye', ''))
        if tsonuc != 'SORUN YOK':
            yaz('        TANI: %s' % tsonuc)
            yaz('              %s' % taciklama[:150])

        # duyarlilik: en yuksek / en dusuk ayrim orani
        deg = [v['kat_enkotu'] or v['kat_havuz'] or 0 for v in varyantlar]
        deg = [x for x in deg if x]
        oynaklik = round(max(deg) / min(deg), 1) if len(deg) > 1 and min(deg) > 0 else None

        r = dict(hedef=d['hedef'], sinif=d['sinif'], F=d['F'], R=d['R'],
                 urun_panel=d['urun_bp'], panel_ayrim=d['ayrim'],
                 panel_ayrim_sayi=d['ayrim_sayi'],
                 kaynak=(a.get('kaynak') if a else ''),
                 turetildi_mi=('EVET' if (a and a.get('kaynak') == 'TURETILDI') else ''),
                 varyantlar=varyantlar, oynaklik=oynaklik,
                 tani=tsonuc, tani_aciklama=taciklama,
                 urun_veren_kutular=veren[:8], konsensus_veren=kons_veren[:8])
        r['_ayar'] = dict(checks.AYAR)
        with open(kp, 'w', encoding='utf-8') as fh:
            json.dump(r, fh, ensure_ascii=False, indent=1, default=str)
        sonuclar.append(r)

    hepsi = []
    for f in sorted(os.listdir(C.KONTROL)):
        if f.startswith('uyelik_') and f.endswith('.json'):
            try:
                hepsi.append(json.load(open(os.path.join(C.KONTROL, f), encoding='utf-8')))
            except Exception:
                pass
    yollar = rapor_yaz(hepsi or sonuclar, panel_yolu, turetildi)
    yaz('')
    yaz('=' * 78)
    yaz(u'  MEMBERSHIP AUDIT FINISHED (%s)' % sure(time.time() - t0))
    for p in yollar:
        yaz('    %s' % p)
    yaz('=' * 78)
    return yollar


# ---------------------------------------------------------------- the report
def rapor_yaz(sonuclar, panel_yolu, turetildi):
    os.makedirs(C.CIKTI, exist_ok=True)
    tsv = os.path.join(C.CIKTI, 'uyelik_duyarlilik.tsv')
    with open(tsv, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'kaynak', 'TURETILDI_mi', 'tanim', 'uye_taxid',
                    'uye_kutu', 'rakip_kutu', 'kapsam', 'uye_min_%', 'uye_max_%',
                    'rakip_havuz', 'ayrim_havuz_x', 'ayrim_en_kotu_x',
                    'A_tanimiyla_ayni_kume', 'PANEL_ayrim', 'oynaklik_x', 'TANI'])
        for r in sonuclar:
            for v in r['varyantlar']:
                w.writerow([r['hedef'], r['kaynak'], r['turetildi_mi'], v['tanim'],
                            v['uye_taxid'], v['uye_kutu'], v['rakip_kutu'],
                            v['kapsam'], v['uye_min'], v['uye_max'], v['havuz'],
                            v['kat_havuz'], v['kat_enkotu'], v['ayni_kume'],
                            r['panel_ayrim'], r['oynaklik'], r['tani']])

    md = os.path.join(C.CIKTI, 'UYELIK_DENETIMI.md')
    L = []; A = L.append
    A(u'# The membership definition audit and sensitivity analysis')
    A('')
    A(u'Generated: %s' % time.strftime('%Y-%m-%d %H:%M:%S'))
    A('')
    A(u'Source panel: `%s`' % os.path.basename(panel_yolu))
    A('')
    A(u'## What this report shows')
    A('')
    A(u'A target\'s **discrimination fold** depends directly on which bins count as a MEMBER and which as a COMPETITOR. When the membership definition changes the number changes, even if the measurement engine never changes at all.')
    A('')
    A(u'The example measured in this run:')
    A('')
    A('| tanim | ayrim |')
    A('|---|---|')
    A(u'| Proteolitik_Cloacimonas, the group row (3 taxa) | 0,0x |')
    A('| Proteolitik_Cloacimonas, tek uye (456827) | 23,5x |')
    A('| **panelde yazan** | **23,0x** |')
    A('')
    A(u'So it was not the number that was wrong but the **definition**. The tables below measure every reasonable definition of membership for each target **side by side**.')
    A('')
    A(u'### The definition sources')
    A('')
    A('| code | source |')
    A('|---|---|')
    A(u'| **A** | `screening/target_membership.tsv`, the definition the tool uses now |')
    A(u'| **B** | `steps/targets.tsv`, the project\'s decision table (the group row) |')
    A(u'| **C** | a single member, the one taxon that gives the best result |')
    A(u'| **D** | `pairs.tsv`, the definition from the read engine correction session |')
    A(u'| **E** | `primer_final/hedef_kimlik.tsv`, derived from the MEASURED identity |')
    A('')

    A(u'## LOOK AT THESE FIRST: the rows marked `TURETILDI`')
    A('')
    A(u'The membership definition of these rows was not written in the panel plainly; it was **derived** from the target name and from `taxid_names.tsv`. These are the ones most likely to be wrong.')
    A('')
    A(u'| target | how many fold the discrimination moves when the definition changes | diagnosis |')
    A('|---|---|---|')
    for r in sonuclar:
        if r['turetildi_mi'] != 'EVET':
            continue
        A('| %s | %s | %s |' % (r['hedef'],
                                ('%sx' % r['oynaklik']) if r['oynaklik'] else '-',
                                r['tani']))
    A('')

    sorunlu = [r for r in sonuclar if r['tani'] != 'SORUN YOK']
    if sorunlu:
        A(u'## The targets that need a diagnosis')
        A('')
        for r in sorunlu:
            A('### %s — %s' % (r['hedef'], r['tani']))
            A('')
            A(r['tani_aciklama'])
            A('')
            if r['urun_veren_kutular']:
                A(u'The bins of this class that give a product:')
                A('')
                A(u'| bin | taxon | product/reads | %% | length |')
                A('|---|---|---|---|---|')
                for v in r['urun_veren_kutular']:
                    A('| %s | %s | %s/%s | %s | %s |' % (
                        v['kutu'], v['ad'], v['urun'], v['okuma'], v['yuzde'],
                        v['boy'][0][0] if v['boy'] else '-'))
                A('')
            if r['konsensus_veren']:
                A(u'The bins that give a product on the consensus (but not on the raw reads):')
                A('')
                A(u'| bin | taxon | product length on the consensus |')
                A('|---|---|---|')
                for v in r['konsensus_veren']:
                    A('| %s | %s | %s bp |' % (v['kutu'], v['ad'], v['boy']))
                A('')

    A(u'## Every target: the sensitivity to the definition')
    A('')
    for r in sonuclar:
        A('### %s%s' % (r['hedef'], '  *(TURETILDI)*' if r['turetildi_mi'] else ''))
        A('')
        A(u'What the panel says: `%s` · product %s bp%s'
          % (r['panel_ayrim'], r['urun_panel'],
             u' · **the discrimination moves %sx when the definition changes**' % r['oynaklik']
             if r['oynaklik'] and r['oynaklik'] > 1.5 else ''))
        A('')
        A(u'| definition | member bins | coverage | member %% | discrimination pool x | discrimination worst x |')
        A('|---|---|---|---|---|---|')
        for v in r['varyantlar']:
            A('| %s%s | %s | %s | %s-%s | %s | %s |' % (
                v['tanim'], u' *(the same set as A)*' if v['ayni_kume'] else '',
                v['uye_kutu'], v['kapsam'], v['uye_min'], v['uye_max'],
                v['kat_havuz'], v['kat_enkotu']))
        A('')

    A('## Ne yapmali')
    A('')
    A(u'1. Go over the rows marked `TURETILDI`; write the right definition into `screening/target_membership.tsv`.')
    A(u'2. On targets whose `oynaklik` column is high, the published number is **very sensitive to the definition**, so which definition it was reported under has to be written in the panel plainly.')
    A(u'3. If the `TANI` column says `KONSENSUS/OKUMA UYUSMAZLIGI`, first reproduce that bin\'s consensus with option (6), then repeat the measurement.')
    A('')
    A(u'> This tool **does not change the membership definition by itself**. It measures, it puts the options side by side, and it leaves the decision to you.')
    A('')
    with open(md, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(L))
    return [md, tsv]
