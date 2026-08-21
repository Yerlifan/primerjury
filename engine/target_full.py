# -*- coding: utf-8 -*-
"""FARKLI_LOKUS - tek hedef icin TAM arama + ARMS + derin olcum.

Akis:
  1. Uye konsensusleri uzunluga gore kumelenir, her kumenin en uzunu omurga.
  2. Her omurga 700 bp pencere / 300 bp adimla taranir. Urun ust siniri 400 bp
     oldugu icin baslangici s olan her amplikon k = s//300 penceresinin ICINE
     tam duser - kapsama bir iddia degil, teorem.
  3. On suzgec OLCUM KURALIDIR (<=1 uyumsuzluk, okuma_motoru.Sonda). 3' son iki
     baz sarti ARTIK UYGULANMIYOR (2026-08-08 karari) - son2=False.
  4. En iyi adaylara ARMS varyantlari (uretec.arms_varyantlari) eklenir.
     ARMS dejenere baz DEGILDIR: tek tanimli baz, tek oligo.
  5. Hepsi PANEL DERINLIGINDE (n=3000) yeniden olculur.
  6. Iki esik birden yazilir: sabit dCq 3 ve bollukla agirlikli
     max(log2(R)+4,3 ; 3,32).
"""
import os, sys, json, time, math, argparse

KOK = os.environ.get('_FL_KOK', '/tmp/fl/kok')
sys.path.insert(0, KOK)
os.environ['_KURTARMA_KOK'] = KOK

from screening import targets as H, motor, uretec as U, numune as N
import importlib.util as _iu
_sp = _iu.spec_from_file_location('kt', os.path.join(KOK, 'verification', 'recovery_round.py'))
kt = _iu.module_from_spec(_sp); _sp.loader.exec_module(kt)
OM = motor.okuma_motoru

URUN_ALT, URUN_UST = 60, 400
PENCERE, ADIM = 700, 300
UZUN = ('A2', 'F2')
SON2 = False          # 2026-08-08: 3' son iki baz olcutu gecersiz sayildi


def gerekli_dcq(R):
    return max(math.log(R, 2) + 4.3, 3.32) if R and R > 0 else 3.32


def baglanir(p, diziler, geri=False):
    s = OM.Sonda(OM.rc(p) if geri else p, geri, 1, SON2)
    return [bool(s.bul(d)) for d in diziler]


def omurgalar(uye_kons):
    d = sorted(uye_kons, key=lambda k: -len(k['dizi'])); kume = []
    for k in d:
        yer = next((g for g in kume if len(g[0]['dizi']) / float(len(k['dizi'])) < 1.5), None)
        (kume.append([k]) if yer is None else yer.append(k))
    return [g[0] for g in kume]


def ssu_sinir(dizi):
    c = kt.capa_bul(motor, dizi)
    bas = c.get('SSU_baslangic', c.get('SSU_orta'))
    if bas is None:
        return None, None
    return int(bas), int(c.get('SSU_son') or c.get('ITS1_baslangic')
                          or min(len(dizi), bas + 1600))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target', '--hedef', dest='hedef', required=True)
    ap.add_argument('--R', type=float, default=0.0)
    ap.add_argument('--duration', '--sure', dest='sure', type=float, default=150.0)
    ap.add_argument('--status', '--durum', dest='durum', default='')
    ap.add_argument('--primer-max', '--primer-ust', dest='primer_ust', type=int, default=1100)
    ap.add_argument('--window-candidates', '--pencere-aday', dest='pencere_aday', type=int, default=6)
    ap.add_argument('--arms-max', '--arms-ust', dest='arms_ust', type=int, default=8)
    ap.add_argument('--sig', type=int, default=900)
    ap.add_argument('--deep', '--derin', dest='derin', type=int, default=3000)
    ap.add_argument('--stage', '--asama', dest='asama', default='tara', choices=['tara', 'derin'])
    ap.add_argument('--member-taxid', '--uye-taxid', dest='uye_taxid', default='', help='set membership MANUALLY (comma-separated)')
    ap.add_argument('--competitor-taxid', '--rakip-taxid', dest='rakip_taxid', default='', help='set competitors MANUALLY')
    ap.add_argument('--class', '--sinif', dest='sinif', default='', help='bin class filter, e.g. F2')
    ap.add_argument('--backbone', '--omurga', dest='omurga', default='',
                    help='choose the template bin MANUALLY (a member bin with a solid consensus)')
    a = ap.parse_args()
    t0 = time.time()
    # Varsayilan durum dosyasi /tmp/fl altindaydi. Iki sakincasi vardi: klasor
    # yoksa cokuyordu (asagidaki duzeltme), ve /tmp her yeniden baslatmada
    # siliniyor - saatlerce suren bir tarama iz birakmadan kayboluyordu. Artik
    # proje icinde duruyor, yarim kalan tarama ertesi gun kaldigi yerden devam
    # ediyor ve sonuc dosyasi elle de okunabiliyor.
    dur = a.durum or os.path.join(KOK, 'FARKLI_LOKUS_SONUC', 'durum_%s.json' % ''.join(
        ch if ch.isalnum() else '_' for ch in a.hedef)[:40])
    # 2026-08-11: durum dosyasinin KLASORU yoksa tarama, ilk pencereyi bitirip
    # sonucu yazmaya calistigi anda FileNotFoundError ile cokuyordu - yani en
    # pahali is (bir pencerelik tarama) yapildiktan SONRA. Klasor bastan
    # acilir; acilamiyorsa durum dosyasi TEK_TUS_SONUC altina alinir ve bu
    # ekrana yazilir. Tarama, durum dosyasi yuzunden cokmez.
    # Yazilabilirlik OS'a sorulur, deneme dosyasi ACILMAZ. (Ilk hali deneme
    # dosyasi yazip siliyordu; bagli klasorde silme yasak oldugu icin kendi
    # denemesi hata verip saglam klasoru "yazilamaz" ilan ediyordu.)
    try:
        _d = os.path.dirname(os.path.abspath(dur))
        if _d and not os.path.isdir(_d):
            os.makedirs(_d)
        if not os.access(_d, os.W_OK):
            raise IOError('klasore yazma izni yok')
    except Exception as _e:
        _yedek = os.path.join(KOK, 'TEK_TUS_SONUC', os.path.basename(dur))
        print('  NOT: durum dosyasi %s yazilamadi (%s).\n'
              '       Durum su dosyaya yazilacak: %s' % (dur, _e, _yedek))
        try:
            os.makedirs(os.path.dirname(_yedek))
        except Exception:
            pass
        dur = _yedek
    D = json.load(open(dur, encoding='utf-8')) if os.path.exists(dur) else \
        dict(hedef=a.hedef, R=a.R, bitmis=[], aday=[], derin=[])

    kons = H.konsensusler(); kut = list({k['kutu']: k for k in H.kutular()}.values())
    kmap = {d['kutu']: d['dizi'] for d in kons}
    import csv
    panel = {r['hedef']: r for r in csv.DictReader(
        (l for l in open(os.path.join(KOK, 'TEK_PROTOKOL_SONUC',
         'panel_tek_protokol.tsv'), encoding='utf-8') if not l.startswith('#')),
        delimiter='\t')}
    # YENI HEDEF (panelde satiri yok) DESTEGI.
    # Hocanin listesindeki taksonlarin cogunun panelde karsiligi yok; eskiden
    # betik burada StopIteration ile cokuyordu. --uye-taxid verildiginde panel
    # satirina zaten ihtiyac yok: uyelik elle veriliyor. Panel satiri varsa
    # mevcut cift de olculur ve karsilastirma basilir; yoksa yalnizca yeni
    # tasarim yapilir.
    ad = next((h for h in panel if a.hedef.lower() in h.lower()), None)
    if ad is None:
        if not a.uye_taxid:
            sys.exit('HATA: "%s" panelde yok. Yeni hedef icin --uye-taxid verin.'
                     % a.hedef)
        sahte = {'hedef': a.hedef, 'sinif': (a.sinif or ''), 'F': '', 'R': ''}
        b = dict(hedef=a.hedef, siniflar=[a.sinif] if a.sinif else [],
                 uye_kutu=[], rakip_kutu=[], uye_kons=[], rakip_kons=[],
                 omurga=None, uye_tax=[], haric=[], anahtar='ELLE',
                 uyelik_kaynagi='ELLE (--uye-taxid)', uyelik_notu='')
        panel[a.hedef] = sahte
        ad = a.hedef
    else:
        b = H.hedef_baglami(panel[ad], kons=kons, kut=kut)
    if a.uye_taxid:
        # KUTU: oneki - uyelik KUTU duzeyinde verilebilir (2026-08-11).
        # Gerekce: ayni taxid'i tasiyan kutular ayni organizma OLMAYABILIR.
        # Olculen ornek: F2-1_40559 ile F2-2_40559 + F2-3_40559 dizi
        # benzerligine gore iki AYRI kumeye dusuyor, ikisine de Kraken
        # "Botrytis cinerea" demis. Taxid ile secince ucu tek hedef olur ve
        # hedefin kendi ici birbirinden ayrilmadigi icin hicbir cift bulunmaz.
        # Bu yol screening/targets.py icindeki KUTU: yolunun aynisidir.
        ham_u = [x for x in a.uye_taxid.split(',') if x]
        ham_r = [x for x in a.rakip_taxid.split(',') if x] if a.rakip_taxid else None
        ku_u = set(x[5:] for x in ham_u if x.startswith('KUTU:'))
        ut = set(x for x in ham_u if not x.startswith('KUTU:'))
        ku_r = set(x[5:] for x in ham_r if x.startswith('KUTU:')) if ham_r else None
        rt = set(x for x in ham_r if not x.startswith('KUTU:')) if ham_r else None
        sf = a.sinif.split(',') if a.sinif else None
        sec = lambda k: (not sf) or k['sinif'] in sf
        uye_mi = lambda k: (k['kutu'] in ku_u) if ku_u else (k['taxid'] in ut)
        if ku_r or rt:
            rak_mi = lambda k: (k['kutu'] in ku_r if ku_r else False) or \
                               (k['taxid'] in rt if rt else False)
        else:
            rak_mi = lambda k: not uye_mi(k)
        b['uye_kutu'] = [k for k in kut if sec(k) and uye_mi(k)]
        b['rakip_kutu'] = [k for k in kut if sec(k) and rak_mi(k)
                           and not uye_mi(k)]
        print(u'MANUAL MEMBERSHIP: %d member bins, %d competitor bins'
              % (len(b['uye_kutu']), len(b['rakip_kutu'])), flush=True)
    uye_k = [k for k in kons if k['kutu'] in {x['kutu'] for x in b['uye_kutu']}]
    rak_k = [k for k in kons if k['kutu'] in {x['kutu'] for x in b['rakip_kutu']}]
    uye_diz = [k['dizi'] for k in uye_k]
    uzun_ix = [j for j, k in enumerate(uye_k) if k['sinif'] in UZUN]
    rak_diz = [k['dizi'] for k in rak_k]
    D['mevcut'] = dict(F=panel[ad]['F'], R=panel[ad]['R'],
                       dCq=panel[ad].get('ASIL_ayrim_mm1'))

    if a.asama == 'tara':
        birim = []
        # --omurga: sablonu ELLE sec (2026-08-12).
        # Gerekce: omurgalar() en uzun konsensusu secer, ama uzun olan SAGLAM
        # olmak zorunda degil. Konsensus saglik denetiminde 99 kutunun 30 unda
        # konsensusun kutunun kendi okumalarini temsil etmedigi olculdu
        # (KONSENSUS_SAGLIK_20260812.xlsx). Bozuk bir sablonla taranan bolge
        # numunede hic yoktur; "cift bulunamadi" sonucu o zaman bir sey
        # kanitlamaz. Bu secenek saglam bir uye kutusunu sablon yapmayi saglar.
        om_liste = omurgalar(uye_k)
        if a.omurga:
            zorla = [k for k in kons if k['kutu'] == a.omurga]
            if not zorla:
                sys.exit('HATA: --omurga %s diye bir kutu yok.' % a.omurga)
            om_liste = zorla
            print('OMURGA ELLE: %s (%d bp)' % (a.omurga, len(zorla[0]['dizi'])),
                  flush=True)
        for om in om_liste:
            L = len(om['dizi']); n = max(1, (L - PENCERE) // ADIM + 2)
            for i in range(n):
                x = min(i * ADIM, max(0, L - PENCERE))
                y = min(L, x + PENCERE)
                if y - x >= 120:
                    birim.append((om['kutu'], i, x, y))
        kalan = [u for u in birim if '%s|%d' % (u[0], u[1]) not in D['bitmis']]
        print('pencere %d, kalan %d' % (len(birim), len(kalan)), flush=True)
        if kalan:
            nm = N.Numune(b['uye_kutu'] + b['rakip_kutu'], n=a.sig, otorite=True)
            for omk, i, x, y in kalan:
                if time.time() - t0 > a.sure:
                    print('SURE DOLDU'); break
                om = kmap[omk]; sb, ss = ssu_sinir(om)
                ap_ = U.aday_primerler(om[x:y])

                def sey(l, t):
                    if len(l) <= t: return l
                    f = len(l) / float(t)
                    return [l[int(j * f)] for j in range(t)]
                Fl, Rl = sey(ap_['F'], a.primer_ust), sey(ap_['R'], a.primer_ust)
                FB = {s: (baglanir(s, uye_diz), baglanir(s, rak_diz))
                      for (_a, _b, s, _c) in Fl}
                RB = {s: (baglanir(s, uye_diz, True), baglanir(s, rak_diz, True))
                      for (_a, _b, s, _c) in Rl}
                havuz = []
                for t in U.cift_akisi(dict(F=Fl, R=Rl)):
                    c = U.cift_yap(t)
                    cb, cs = x + c['iF'], x + c['iF'] + c['urun']
                    dis = (sb is not None and (cb >= ss or cs <= sb))
                    fu, fr = FB[c['F']]; ru, rr = RB[c['R']]
                    ix = uzun_ix if dis else range(len(uye_diz))
                    if not ix or not all(fu[j] and ru[j] for j in ix):
                        continue
                    ro = sum(1 for j in range(len(rak_diz)) if fr[j] and rr[j])
                    havuz.append((ro, 0 if c['urun'] <= 150 else 1, c['urun'],
                                  c['F'], c['R'], cb, cs, dis))
                havuz.sort()
                for ro, _p, urun, F, R, cb, cs, dis in havuz[:a.pencere_aday]:
                    uy = ([k for k in b['uye_kutu'] if k['sinif'] in UZUN] if dis
                          else b['uye_kutu'])
                    rk = ([k for k in b['rakip_kutu'] if k['sinif'] in UZUN] if dis
                          else b['rakip_kutu'])
                    if not uy: continue
                    o = nm.olc(F, R, uy, rk, lo=URUN_ALT, hi=URUN_UST, mm=1)
                    if not o: continue
                    D['aday'].append(dict(F=F, R=R, urun=urun, omurga=omk,
                                          bas=cb, son=cs, ssu_disi=dis,
                                          taban=('A2/F2' if dis else 'panel'),
                                          sig_kat=o.get('kat_enkotu'),
                                          kapsam=o.get('uye_kapsam_pay'),
                                          arms=''))
                D['bitmis'].append('%s|%d' % (omk, i))
                en = max([q['sig_kat'] or -1 for q in D['aday']] or [-1])
                print('  %-13s %5d-%5d gecen %3d  simdiye en iyi(sig) %.2fx'
                      % (omk, x, y, len(havuz), en), flush=True)
                json.dump(D, open(dur, 'w', encoding='utf-8'), ensure_ascii=False)
        json.dump(D, open(dur, 'w', encoding='utf-8'), ensure_ascii=False)
        print('TARAMA BITTI' if len(D['bitmis']) >= len(birim) else 'DEVAM')
        return

    # ---- derin asama: en iyiler + ARMS varyantlari, n=3000 ----
    ad_l = sorted(D['aday'], key=lambda q: -(q['sig_kat'] or -1))
    gor, sec = set(), []
    for q in ad_l:
        if (q['F'], q['R']) in gor: continue
        gor.add((q['F'], q['R'])); sec.append(q)
        if len(sec) >= 20: break
    genis = list(sec)
    for q in sec[:a.arms_ust]:
        for v, et in U.arms_varyantlari(q['F']):
            genis.append(dict(q, F=v, arms='F:' + et, sig_kat=None))
        for v, et in U.arms_varyantlari(q['R']):
            genis.append(dict(q, R=v, arms='R:' + et, sig_kat=None))
    # mevcut panel cifti = kontrol
    genis.insert(0, dict(F=D['mevcut']['F'], R=D['mevcut']['R'], urun=0,
                         omurga='-', bas=-1, son=-1, ssu_disi=False,
                         taban='panel', sig_kat=None, kapsam='', arms='MEVCUT'))
    yap = [q for q in genis if [q['F'], q['R'], q['arms']] not in
           [[z['F'], z['R'], z['arms']] for z in D['derin']]]
    print(u'%d to be measured at full depth (of %d)' % (len(yap), len(genis)), flush=True)
    nm = N.Numune(b['uye_kutu'] + b['rakip_kutu'], n=a.derin, otorite=True)
    ger = gerekli_dcq(a.R or D.get('R') or 0)
    for q in yap:
        if time.time() - t0 > a.sure:
            print('SURE DOLDU'); break
        uy = ([k for k in b['uye_kutu'] if k['sinif'] in UZUN] if q['ssu_disi']
              else b['uye_kutu'])
        rk = ([k for k in b['rakip_kutu'] if k['sinif'] in UZUN] if q['ssu_disi']
              else b['rakip_kutu'])
        o = nm.olc(q['F'], q['R'], uy, rk, lo=URUN_ALT, hi=URUN_UST, mm=1)
        if not o: continue
        kat = o.get('kat_enkotu')
        dcq = (math.log(kat, 2) if kat and kat > 0 else None)
        D['derin'].append(dict(q, kat=kat, dCq=dcq,
                               kapsam=o.get('uye_kapsam_pay'),
                               uye_alt=o.get('uye_alt'), havuz=o.get('havuz'),
                               enkotu=o.get('en_kotu_rakip'),
                               gerekli_dCq=ger, gecer_bolluk=(dcq is not None and dcq >= ger),
                               gecer_sabit3=(dcq is not None and dcq >= 3.0)))
        json.dump(D, open(dur, 'w', encoding='utf-8'), ensure_ascii=False)
    D['derin'].sort(key=lambda z: -(z['dCq'] if z['dCq'] is not None else -99))
    json.dump(D, open(dur, 'w', encoding='utf-8'), ensure_ascii=False)
    for z in D['derin'][:12]:
        print(u'  dCq %6s  fold %8s  coverage %-6s product %3s %-14s %s' % (
            ('%.2f' % z['dCq']) if z['dCq'] is not None else '-',
            ('%.2f' % z['kat']) if z['kat'] is not None else '-',
            z['kapsam'], z['urun'], z['arms'] or '-', z['F'] + '/' + z['R']))
    print('gerekli dCq (bolluk) %.2f ; sabit 3,00' % ger)


if __name__ == '__main__':
    main()
