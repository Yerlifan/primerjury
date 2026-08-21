# -*- coding: utf-8 -*-
"""FARKLI_LOKUS - 2. adim: BUTUN LOKUSLARDA tasarim denemesi.

Yol 5'ten farki dort noktada:

  1. TAM KAPSAMA GARANTISI. Yol 5 capa tabanli bolgelere ayirir; capa tutmazsa
     kaba parcalara boler. Burada omurga 700 bp pencerelerle 300 bp adimla
     tanir. Urun ust siniri 400 bp oldugundan, baslangici s olan her amplikon
     k = s // 300 penceresinin ICINE tam duser (300k <= s < 300k+300 ve
     s + 400 < 300k + 700). Yani hicbir olasi amplikon pencere sinirina
     dusmuyor - "butun dizi tarandi" bir iddia degil, bir teorem.

  2. COK OMURGA. Yol 5 yalniz EN UZUN uye konsensusunu omurga alir. Burada uye
     konsensusleri uzunluga gore kumelenir ve HER kumenin en uzunu ayri omurga
     olur (or. Methanothrix icin hem A1 1348 bp hem A2 4329 bp).

  3. SUZGEC OLCUYE BAGLANDI. Yol 5 uretec.ayirt_edici_mi kullanir: "3' son baz
     rakipte uymuyor mu". Buradaki suzgec dogrudan OLCUM KURALIDIR (<=1
     uyumsuzluk + 3' son 2 baz TAM, okuma_motoru.Sonda). Yani on suzgec ile
     nihai olcu ayni kurali kullanir; on suzgecte elenen bir aday olcumde
     gecemez. Iki arama ayni degildir, bu fark raporda yazilidir.

  4. OLCUM TABANI ONCEDEN VE KOSULSUZ BELIRLENIR. Bir lokus fiziksel olarak
     yalniz bazi kitapliklarda varsa (or. 23S yalniz A2'de) o lokusun adaylari
     YALNIZ o kitapliklarda olculur. Kural adayin sonucuna DEGIL amplikonun
     konumuna bakar ve rakip kumesine de ayni sekilde uygulanir.

Iki asamali olcum: on eleme sig havuzda (--sig), hayatta kalanlar panel
derinliginde (--derin) yeniden olculur. Sig sayilar SIRALAMA icindir, rapora
girmez.
"""
import os, sys, json, time, argparse

KOK = os.environ.get('_FL_KOK', '/tmp/fl/kok')
sys.path.insert(0, KOK)
os.environ['_KURTARMA_KOK'] = KOK

from screening import targets as H, motor, uretec as U, numune as N
import importlib.util as _iu
_sp = _iu.spec_from_file_location('kt', os.path.join(KOK, 'verification', 'recovery_round.py'))
kt = _iu.module_from_spec(_sp)
_sp.loader.exec_module(kt)
OM = motor.okuma_motoru

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

TABAN_KURALI = (u'Amplikonun tamami omurgadaki SSU bolgesinin disindaysa olcum '
                u'tabani tam operon kitapliklarina (A2/F2) daraltilir; uye ve '
                u'rakip tarafinda ayni daraltma uygulanir. SSU icindeki adaylar '
                u'panelin kendi tabaninda olculur. Kural adaydan ONCE yazildi.')


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
    c = kt.capa_bul(motor, dizi)
    bas = c.get('SSU_baslangic', c.get('SSU_orta'))
    if bas is None:
        return None, None
    son = c.get('SSU_son') or c.get('ITS1_baslangic') or min(len(dizi), bas + 1600)
    return int(bas), int(son)


def baglanir(primer, diziler, geri=False):
    """Olcum kuralinin ta kendisi: <=1 uyumsuzluk + 3' son 2 baz TAM."""
    s = OM.Sonda(OM.rc(primer) if geri else primer, geri, 1, True)
    return [bool(s.bul(d)) for d in diziler]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--duration', '--sure', dest='sure', type=float, default=34.0)
    ap.add_argument('--status', '--durum', dest='durum', default='/tmp/fl/cl_durum.json')
    ap.add_argument('--sig', type=int, default=900)
    ap.add_argument('--pre-candidate', '--on-aday', dest='on_aday', type=int, default=8)
    ap.add_argument('--primer-max', '--primer-ust', dest='primer_ust', type=int, default=1100)
    ap.add_argument('--target', '--hedef', dest='hedef', default='')
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
        (l for l in open(os.path.join(KOK, 'TEK_PROTOKOL_SONUC',
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
        print('BITTI')
        return

    ger = {}
    for ad in {x[0] for x in kalan}:
        for k in baglam[ad]['b']['uye_kutu'] + baglam[ad]['b']['rakip_kutu']:
            ger[k['kutu']] = k
    tt = time.time()
    nm = N.Numune(list(ger.values()), n=a.sig, otorite=True)
    print('sig havuz %d kutu, %.1f sn' % (len(ger), time.time() - tt), flush=True)

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

        # primer basina baglanma vektorleri (olcum kuraliyla)
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
            ('%.2fx kapsam %s %s urun %d' % (en_iyi['_skor'], en_iyi['kapsam'],
                                             en_iyi['taban'], en_iyi['urun']))
            if en_iyi else '-'), flush=True)
        json.dump(D, open(a.durum, 'w', encoding='utf-8'), ensure_ascii=False)

    json.dump(D, open(a.durum, 'w', encoding='utf-8'), ensure_ascii=False)
    if len(D['bitmis']) >= len(birim):
        print('BITTI')


if __name__ == '__main__':
    main()
