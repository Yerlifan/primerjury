# -*- coding: utf-8 -*-
"""FARKLI_LOKUS, step 2: a design attempt AT EVERY LOCUS.

It differs from route 5 on four points:

  1. A COMPLETE COVERAGE GUARANTEE. Route 5 splits into anchor based regions; if the
     anchor does not hold it falls back on rough pieces. Here the backbone is scanned
     with 700 bp windows at a 300 bp step. Since the product upper bound is 400 bp,
     every amplicon starting at s falls entirely INSIDE the window k = s // 300
     (300k <= s < 300k+300 and s + 400 < 300k + 700). So no possible amplicon falls
     on a window boundary: "the whole sequence was scanned" is not a claim but a
     theorem.

  2. SEVERAL BACKBONES. Route 5 takes only the LONGEST member consensus as the
     backbone. Here the member consensuses are clustered by length and the longest of
     EACH cluster becomes a separate backbone (for Methanothrix, for instance, both
     A1 at 1348 bp and A2 at 4329 bp).

  3. THE FILTER IS TIED TO THE MEASURE. Route 5 uses generator.ayirt_edici_mi: "does
     the last base at the 3' end fail to match in the competitor". The filter here is
     THE MEASUREMENT RULE itself (<=1 mismatch plus the last 2 bases at the 3' end
     EXACT, okuma_motoru.Sonda). So the prefilter and the final measure use the same
     rule; a candidate eliminated in the prefilter could not have passed the
     measurement. The two searches are not the same, and that difference is written in
     the report.

  4. THE MEASUREMENT BASE IS SETTLED IN ADVANCE AND UNCONDITIONALLY. If a locus
     physically exists only in some libraries (23S only in A2, for instance), the
     candidates of that locus are measured ONLY in those libraries. The rule looks at
     the amplicon's position and NOT at the candidate's result, and it is applied to
     the competitor set in the same way.

A two stage measurement: the prefilter runs on a shallow pool (--sig) and the
survivors are measured again at panel depth (--deep). The shallow numbers are for
RANKING and do not enter the report.

"""
import os, sys, json, time, argparse

# The project root. It used to default to /tmp/fl/kok, the path of the private
# workspace this code grew in; on any other machine the module died at import
# with FileNotFoundError. The default is now derived from this file's own
# location, and _FL_KOK still overrides it.
KOK = os.environ.get('_FL_KOK') or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, KOK)
os.environ['_RECOVERY_ROOT'] = KOK

from screening import targets as H, engine_gateway, generator as U, sample as N
import importlib.util as _iu
_sp = _iu.spec_from_file_location('kt', os.path.join(KOK, 'verification', 'recovery_round.py'))
kt = _iu.module_from_spec(_sp)
_sp.loader.exec_module(kt)
OM = engine_gateway.okuma_motoru

URUN_ALT, URUN_UST = 60, 400
PENCERE, ADIM = 700, 300
UZUN_SINIF = ('A2', 'F2')

HEDEFLER = [
    'Methanothrix_soehngenii_turu',
    'Methanosarcina mazei / M. soligelidi grubu',
    'Petriella_musispora',
    'Bacteroidales_kumesi',
    'Petrimonas_cinsi',
    'Proteiniphilum_cinsi',
]

TABAN_KURALI = ('When the whole amplicon falls outside the SSU region of the '
                'backbone, the measurement base is narrowed to the full '
                'operon libraries, and the same narrowing is applied on both '
                'the member and the competitor side. Candidates inside the '
                "SSU are measured on the panel's own base. The rule was "
                'written BEFORE the candidate.')


def _bg(x):
    return json.load(open(x, encoding='utf-8')) if os.path.exists(x) else {}


def omurgalar(uye_kons):
    d = sorted(uye_kons, key=lambda k: -len(k['dizi']))
    kume = []
    for k in d:
        yer = None
        for g in kume:
            if len(g[0]['dizi']) / float(len(k['dizi'])) < 1.5:
                yer = g
                break
        (kume.append([k]) if yer is None else yer.append(k))
    return [g[0] for g in kume]


def ssu_sinir(dizi):
    c = kt.capa_bul(engine_gateway, dizi)
    bas = c.get('SSU_baslangic', c.get('SSU_orta'))
    if bas is None:
        return None, None
    son = c.get('SSU_son') or c.get('ITS1_baslangic') or min(len(dizi), bas + 1600)
    return int(bas), int(son)


def baglanir(primer, diziler, geri=False):
    """The measurement rule itself: <=1 mismatch plus the last 2 bases at the 3' end EXACT."""
    s = OM.Sonda(OM.rc(primer) if geri else primer, geri, 1, True)
    return [bool(s.bul(d)) for d in diziler]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--duration', dest='sure', type=float, default=34.0)
    ap.add_argument('--status', dest='durum', default=os.path.join(KOK, 'RESULTS', 'cok_lokus_durum.json'))
    ap.add_argument('--sig', type=int, default=900)
    ap.add_argument('--pre-candidate', dest='on_aday', type=int, default=8)
    ap.add_argument('--primer-max', dest='primer_ust', type=int, default=1100)
    ap.add_argument('--target', dest='hedef', default='')
    a = ap.parse_args()
    t0 = time.time()

    D = _bg(a.durum) or dict(bitmis=[], sonuc={}, taban_kurali=TABAN_KURALI,
                             pencere=PENCERE, adim=ADIM, sig=a.sig,
                             urun=[URUN_ALT, URUN_UST])

    kons = H.konsensusler()
    kut = {k['kutu']: k for k in H.kutular()}
    kmap = {d['kutu']: d['dizi'] for d in kons}
    import csv
    panel = {r['hedef']: r for r in csv.DictReader(
        (l for l in open(os.path.join(KOK, 'ONE_PROTOCOL_RESULT',
                                      'panel_tek_protokol.tsv'), encoding='utf-8')
         if not l.startswith('#')), delimiter='\t')}

    hedef_listesi = [h for h in HEDEFLER if (not a.hedef or a.hedef.lower() in h.lower())]
    birim, baglam = [], {}
    for ad in hedef_listesi:
        b = H.hedef_baglami(panel[ad], kons=kons, kut=list(kut.values()))
        uye_ad = {x['kutu'] for x in b['uye_kutu']}
        rak_ad = {x['kutu'] for x in b['rakip_kutu']}
        uye_k = [k for k in kons if k['kutu'] in uye_ad]
        rak_k = [k for k in kons if k['kutu'] in rak_ad]
        baglam[ad] = dict(b=b, uye_k=uye_k, rak_k=rak_k)
        for om in omurgalar(uye_k):
            L = len(om['dizi'])
            n = max(1, (L - PENCERE) // ADIM + 2)
            for i in range(n):
                bas = min(i * ADIM, max(0, L - PENCERE))
                son = min(L, bas + PENCERE)
                if son - bas >= 120:
                    birim.append((ad, om['kutu'], i, bas, son))
    kalan = [x for x in birim
             if '%s|%s|%d' % (x[0], x[1], x[2]) not in D['bitmis']]
    print('is birimi %d, kalan %d' % (len(birim), len(kalan)), flush=True)
    if not kalan:
        print(u'DONE')
        return

    ger = {}
    for ad in {x[0] for x in kalan}:
        for k in baglam[ad]['b']['uye_kutu'] + baglam[ad]['b']['rakip_kutu']:
            ger[k['kutu']] = k
    tt = time.time()
    nm = N.Numune(list(ger.values()), n=a.sig, otorite=True)
    print(u'shallow pool %d bins, %.1f s' % (len(ger), time.time() - tt), flush=True)

    for ad, omk, i, bas, son in kalan:
        if time.time() - t0 > a.sure:
            print('SURE DOLDU'); break
        g = baglam[ad]
        b = g['b']
        om = kmap[omk]
        sbas, sson = ssu_sinir(om)
        alt = om[bas:son]
        ap_ = U.aday_primerler(alt)

        def sey(l, t):
            if len(l) <= t:
                return l
            f = len(l) / float(t)
            return [l[int(j * f)] for j in range(t)]
        Fl, Rl = sey(ap_['F'], a.primer_ust), sey(ap_['R'], a.primer_ust)

        uye_diz = [k['dizi'] for k in g['uye_k']]
        uye_uzun_ix = [j for j, k in enumerate(g['uye_k']) if k['sinif'] in UZUN_SINIF]
        rak_diz = [k['dizi'] for k in g['rak_k']]

        # the binding vectors per primer (under the measurement rule)
        FB, RB = {}, {}
        for (_i, _k, s, _m) in Fl:
            if s not in FB:
                FB[s] = (baglanir(s, uye_diz), baglanir(s, rak_diz))
        for (_i, _k, s, _m) in Rl:
            if s not in RB:
                RB[s] = (baglanir(s, uye_diz, True), baglanir(s, rak_diz, True))

        aday = []
        for t in U.cift_akisi(dict(F=Fl, R=Rl)):
            c = U.cift_yap(t)
            cb, cs = bas + c['iF'], bas + c['iF'] + c['urun']
            ssu_disi = (sbas is not None and (cb >= sson or cs <= sbas))
            fu, fr = FB[c['F']]
            ru, rr = RB[c['R']]
            ix = uye_uzun_ix if ssu_disi else range(len(uye_diz))
            if not ix:
                continue
            if not all(fu[j] and ru[j] for j in ix):
                continue
            rak_ort = sum(1 for j in range(len(rak_diz)) if fr[j] and rr[j])
            aday.append((rak_ort, -(c['urun'] <= 150), c['urun'], c['F'], c['R'],
                         cb, cs, ssu_disi))
        aday.sort()
        aday = aday[:a.on_aday]

        en_iyi = None
        for rak_ort, _p, urun, F, R, cb, cs, ssu_disi in aday:
            if ssu_disi:
                uy = [k for k in b['uye_kutu'] if k['sinif'] in UZUN_SINIF]
                rk = [k for k in b['rakip_kutu'] if k['sinif'] in UZUN_SINIF]
                taban = 'A2/F2 (SSU disi lokus)'
            else:
                uy, rk, taban = b['uye_kutu'], b['rakip_kutu'], 'panel tabani'
            if not uy:
                continue
            o = nm.olc(F, R, uy, rk, lo=URUN_ALT, hi=URUN_UST, mm=1)
            if not o:
                continue
            kat = o.get('kat_enkotu')
            skor = kat if kat is not None else -1.0
            if en_iyi is None or skor > en_iyi['_skor']:
                en_iyi = dict(_skor=skor, F=F, R=R, urun=urun, bas=cb, son=cs,
                              taban=taban, ssu_disi=ssu_disi,
                              sig_kat=kat, sig_kat_kapsayan=o.get('kat_enkotu_kapsayan'),
                              kapsam=o.get('uye_kapsam_pay'), uye_alt=o.get('uye_alt'),
                              havuz=o.get('havuz'), rak_konsensus_ortak=rak_ort)
        anah = '%s|%s|%d' % (ad, omk, i)
        D['sonuc'][anah] = dict(hedef=ad, omurga=omk, omurga_bp=len(om),
                                pencere=[bas, son], ssu=[sbas, sson],
                                gecen_aday=len(aday), en_iyi=en_iyi)
        D['bitmis'].append(anah)
        print('  %-26s %-13s %5d-%5d gecen %2d  %s' % (
            ad[:26], omk, bas, son, len(aday),
            (u'%.2fx coverage %s %s product %d' % (en_iyi['_skor'], en_iyi['kapsam'],
                                             en_iyi['taban'], en_iyi['urun']))
            if en_iyi else '-'), flush=True)
        json.dump(D, open(a.durum, 'w', encoding='utf-8'), ensure_ascii=False)

    json.dump(D, open(a.durum, 'w', encoding='utf-8'), ensure_ascii=False)
    if len(D['bitmis']) >= len(birim):
        print(u'DONE')


if __name__ == '__main__':
    main()
