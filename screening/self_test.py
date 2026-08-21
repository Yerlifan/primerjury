# -*- coding: utf-8 -*-
"""Arama baslamadan once calisan kendini sinama.

Amac: kullanicinin makinesinde ilk adimda patlamamak ve daha onemlisi
SESSIZCE YANLIS OLCMEMEK. Sinamalardan biri bile duserse arama baslatilmaz.

1. Gerekli paketler ve dosyalar yerinde mi.
2. Geometri modulu, projede zaten calisan engine/geometry_core.py ile
   BIREBIR ayni sayilari uretiyor mu (42 panel primerinin hepsi).
3. ispcr motoru paneldeki ciftleri kendi konsensuslerinde dogru urun
   boyuyla buluyor mu.
4. UC MOTOR KARSILASTIRMASI (en onemlisi):
      brute_force.py  (tohumsuz referans)
   == read_engine.py (duzeltilmis, guvercin yuvasi tohumu)
   == numpy havuz     (bu aracin hizli yolu)
   Ayrica panelin ESKI motorunun (reads.py/Sonda) ne kadar site kacirdigi
   olculup rapor edilir.
5. Pencere ureteci gercekten HER pozisyon x HER uzunlugu uretiyor mu.
"""
# ---------------------------------------------------------------------------
# self_test.py — olcum baslamadan once motorlarin dogrulugunu sinayan kapi;
#                   bir sinama bile duserse hicbir asama baslatilmaz.
#
# GIRDI  : config.py'deki yollar (panel TSV, hedefler.tsv, konsensus ve
#          fastq klasorleri, hedef_uyelik.tsv, SILVA); engine/geometry_core.py
#          (ayri surecte calistirilip geo.json ciktisi okunur); panel ciftleri
#          ve kanonik konsensusler targets.py uzerinden; motorlar engine_gateway.py
#          uzerinden (ispcr, tarayici, okuma_motoru, kaba_kuvvet, eski reads.py).
# CIKTI  : dosyaya yazmaz. Ekrana satir satir GECTI / *** DUSTU *** basar ve
#          calistir() tek bir True/False dondurur. yon_sinamasi() ayrica tek
#          basina cagrilabilir ve primer3 gibi opsiyonel bagimliliklara takilmaz.
# CAGRAN : __main__.main() her modun basinda (--sinama-atla verilmedikce) ve
#          hepsi.calistir'in 1. asamasi. Yani full_chain.py asamalari 1, 2,
#          3, 4, 5, 6, 7, 9 ve dogrudan tus 8 (--sina, yalniz sinama).
#
# Sinamanin 4. maddesi uc motoru karsilastirir: brute_force.py (tohumsuz
# referans), read_engine.py (guvercin yuvasi tohumu) ve numpy havuz yolu.
# Ucu birebir ayni sayiyi vermek zorundadir; vermezse tohumlamanin kayipsizligi
# artik dogrulanmamis demektir ve kosu baslamaz.
# ---------------------------------------------------------------------------
import os, sys, json, subprocess, tempfile
from . import config as C


def _ok(yaz, ad, iyi, ek=''):
    yaz('   %-58s %s %s' % (ad, 'GECTI' if iyi else '*** DUSTU ***', ek))
    return iyi


def yon_sinamasi(yaz):
    """YON NORMALIZASYONU sinamasi - BAGIMSIZ calistirilabilir.
    primer3 gibi opsiyonel bagimliliklara TAKILMAZ; konsensus adimi bunu
    her kosuda cagirir. Donus: True/False."""
    from . import read_engine as OM
    from . import targets as H
    tum = True
    # Yon hatasi gece boyunca uc ayri yerde ayri ayri yamandi. Bu sinama, kanonik
    # yonun bir daha sessizce kaymamasi icin konur: kanonik kaynak yoksa ya da
    # icinde ters yonlu dosya varsa ANA IS BASLAMAZ.
    try:
        from . import orientation as _Y
        h = _Y.kendini_sina()
        tum &= _ok(yaz, 'orientation.py kendi sinavi (5 madde)', not h, '; '.join(h)[:70])

        kons = H.konsensusler()
        ters = [k['kutu'] for k in kons
                if _Y.tespit(k['dizi'], _Y.sinifi(k['kutu']))[0] == 'ANTISENSE']
        tum &= _ok(yaz, 'kanonik konsensuslerin HEPSI sense yonde',
                   not ters, '%d ters: %s' % (len(ters), ', '.join(ters[:3])))

        # bilinen cevapli: evrensel cift kanonik sette urun vermeli, tersinde vermemeli
        F, R = 'CTGCGGTTTAATTGGATTCAACGC', 'GAACTGACGACGGCCATGC'
        fs = OM.Sonda(F, False, 3); rs = OM.Sonda(OM.rc(R), True, 3)
        A = [k for k in kons if k['sinif'] == 'A']
        d = sum(1 for k in A
                if OM.urun_var(k['dizi'], fs, rs, len(F), len(R), 40, 600) is not None)
        t = sum(1 for k in A
                if OM.urun_var(OM.rc(k['dizi']), fs, rs, len(F), len(R), 40, 600) is not None)
        tum &= _ok(yaz, 'yon etkisi: Arke_universal dogru yonde urun, tersinde yok',
                   d > 0 and t == 0, 'dogru %d/%d, ters %d/%d' % (d, len(A), t, len(A)))
    except Exception as e:
        tum &= _ok(yaz, 'yon normalizasyonu sinamasi', False, str(e)[:90])

    return tum


def calistir(yaz):
    from . import read_engine as OM
    yaz('\nKENDINI SINAMA')
    tum = True

    # ---- 1 paketler
    try:
        import numpy
        tum &= _ok(yaz, 'numpy', True, numpy.__version__)
    except ImportError:
        tum &= _ok(yaz, 'numpy', False, 'pip3 install numpy --break-system-packages')
        return False
    try:
        import primer3
        tum &= _ok(yaz, 'primer3-py', True)
    except ImportError:
        tum &= _ok(yaz, 'primer3-py', False,
                   'pip3 install primer3-py --break-system-packages')
        return False

    from . import engine_gateway, geometri as G, hedefler as H, uretec as U, numune as N

    # ---- 1b dosyalar
    for ad, p in (('panel TSV', C.PANEL_TSV), ('hedefler.tsv', C.HEDEFLER_TSV),
                  ('konsensus klasoru', C.KONSENSUS), ('fastq klasoru', C.FASTQ),
                  ('hedef_uyelik.tsv', H.UYELIK_TSV)):
        tum &= _ok(yaz, ad, os.path.exists(p), '' if os.path.exists(p) else p)
    var_ref = os.path.exists(C.SILVA_SSU)
    _ok(yaz, 'REFERANS_DB/SILVA SSU (kuresel adim icin)', var_ref,
        '' if var_ref else '-> kuresel adim atlanacak')

    bilgi = motor.surum_bilgisi()
    for ad in ('ispcr', 'tarayici', 'cift', 'okuma_motoru', 'kaba_kuvvet'):
        p = bilgi.get(ad)
        tum &= _ok(yaz, 'ice aktarilan motor: %s' % ad, p is not None,
                   os.path.basename(p or '') +
                   (' ' + bilgi.get('okuma_motoru_surum', '')
                    if ad == 'okuma_motoru' else ''))

    # ---- 2 geometri == geometry_core.py
    geo_py = None
    for kok in C.BETIK_YOLLARI:
        q = os.path.join(kok, 'geometry_core.py')
        if os.path.exists(q):
            geo_py = q
            break
    if geo_py:
        try:
            with tempfile.TemporaryDirectory() as td:
                subprocess.run([sys.executable, geo_py], cwd=td,
                               capture_output=True, timeout=600)
                rows = json.load(open(os.path.join(td, 'geo.json'), encoding='utf-8'))
            n = f = 0
            ilk = ''
            for row in rows:
                if row['primer'] == 'CIFT':
                    continue
                p = row['dizi']
                n += 1
                m = G.olc(p)
                for k, ref in (('gc', row['gc']), ('tm', row['tm']),
                               ('hp_tm', row['hp']), ('hd_tm', row['hd']),
                               ('son5', row['gc5']), ('uc', row['uc'])):
                    if str(m[k]) != str(ref):
                        f += 1
                        if not ilk:
                            ilk = '%s %s: bizde %s, geometry_core.py %s' % (p, k, m[k], ref)
                        break
            tum &= _ok(yaz, 'geometri == engine/geometry_core.py (%d primer)' % n,
                       f == 0 and n >= 42, ilk or '%d/%d birebir' % (n - f, n))
        except Exception as e:
            tum &= _ok(yaz, 'geometri == geometry_core.py', False, str(e)[:70])
    else:
        _ok(yaz, 'geometry_core.py bulunamadi (karsilastirma atlandi)', True)

    # ---- 3 ispcr panel urun boylarini dogruluyor mu
    #
    # NEDEN ARALIK KABUL EDIYOR (2026-08-02 duzeltmesi):
    # Bu sinama once "ilk bulunan konsensusteki boy == panelin tek sayisi" diye
    # bakiyordu ve Bakteri_universal'de dustu (olculen 135, panel 130). Teshis:
    # o cift GERCEKTEN tek boy vermiyor - uye konsensuslerinde 130/133/135/140
    # cikiyor. Panel bunu zaten biliyor ve "URUN BOYU ARALIGI (numunede olculen)"
    # sutununda "129-135 bp (iki tepeli: 130 ve 135)" diye kaydetmis; sinama o
    # sutunu HIC OKUMUYORDU. Ayrica hangi konsensusun ilk gelecegi kaynak
    # degisince kayiyordu.
    # Dogru olcut: panelin bildirdigi boy, uye konsensuslerinde olculen boylar
    # kumesinde BULUNMALI (+-2 bp tolerans, konsensus indelleri icin). Sinama
    # artik ilk uyusmazlikta durmuyor, HEPSINI listeliyor.
    try:
        panel, _ = H.panel_oku()
        kons = H.konsensusler()
        uyelik = H.uyelik_oku()
        kut = H.kutular()
        TOLERANS = 2
        tutan = sapan = urunsuz = 0
        ayrinti = []
        for d in panel:
            b = H.hedef_baglami(d, uyelik, kons, kut)
            if not b['uye_kons']:
                continue
            boylar = {}
            for k in b['uye_kons']:
                for s_ in (k['dizi'], motor.rc(k['dizi'])):
                    pr = motor.amplify(s_, d['F'], d['R'], max_mm=1, lo=40, hi=600)
                    if pr:
                        bp = min(pr, key=lambda x: x[3] + x[4])[2]
                        boylar[bp] = boylar.get(bp, 0) + 1
                        break
            if not boylar:
                urunsuz += 1
                ayrinti.append(('urun yok', d['hedef'], ''))
                continue
            if any(abs(bp - d['urun_bp']) <= TOLERANS for bp in boylar):
                tutan += 1
                if len(boylar) > 1:
                    ayrinti.append(('aralik', d['hedef'],
                                    ','.join('%d(%d)' % x for x in sorted(boylar.items()))))
            else:
                sapan += 1
                ayrinti.append(('SAPMA', d['hedef'],
                                'panel %d, olculen %s' % (
                                    d['urun_bp'],
                                    ','.join('%d(%d)' % x for x in sorted(boylar.items())))))
        tum &= _ok(yaz, 'ispcr panel urun boylarini dogruluyor (+-%d bp)' % TOLERANS,
                   sapan == 0 and tutan >= 15,
                   '%d tutan, %d sapan, %d urunsuz' % (tutan, sapan, urunsuz))
        for tip, ad, ek in ayrinti:
            if tip == 'SAPMA':
                yaz('        *** %s: %s' % (ad[:34], ek))
            elif tip == 'aralik':
                yaz('        (aralik) %-32s %s' % (ad[:32], ek))
            else:
                yaz('        (urun yok) %s' % ad[:40])
    except Exception as e:
        tum &= _ok(yaz, 'ispcr panel dogrulamasi', False, str(e)[:70])

    # ---- 4 UC MOTOR KARSILASTIRMASI
    try:
        om = motor.okuma_motoru
        kk = motor.kaba_kuvvet
        panel, _ = H.panel_oku()
        kut = H.kutular()
        farkli_om = farkli_hv = denenen = 0
        en_kotu = (2.0, '', '', 0, 0)
        # panelin BASINDAN degil, YAYILMIS bir ornek: kotu durumlar panelin
        # her yerinde olabilir (Asetoklastik / Arke_universal gibi)
        ornekler = panel[::max(1, len(panel) // 8)][:8]
        for d in ornekler:
            sf = [x.strip() for x in (d['sinif'] or '').split('/') if x.strip()]
            ks = [k for k in kut if k['sinif'] in sf][:1]
            for k in ks:
                rd, _n0 = om.kutu_yukle(k['yol'], nmax=60, seed=C.NUMUNE_TOHUM,
                                        minl=C.NUMUNE_OKUMA_MIN,
                                        maxl=C.NUMUNE_OKUMA_MAX)
                if not rd:
                    continue
                F = d['F']
                n_kk = n_om = n_es = 0
                sonda_yeni = om.Sonda(F, False, 1, True)
                sonda_eski = motor.okuma.Sonda(F, False, 1) if motor.okuma else None
                for s in rd:
                    t = om.temizle(s)
                    for q in (t, om.rc(t)):
                        if kk is not None:
                            n_kk += len(kk.yerler(q, F, 1, son2=True, uc5=False))
                        n_om += len(sonda_yeni.bul(q))
                        if sonda_eski is not None:
                            n_es += len(sonda_eski.bul(q))
                qs = []
                for s in rd:
                    qs.append(motor.clean(s))
                    qs.append(motor.clean(motor.rc(s)))
                n_hv = int(motor.tarayici.Havuz(qs).bul(F, 1).size)
                denenen += 1
                ref = n_kk if kk is not None else n_om
                if kk is not None and n_om != n_kk:
                    farkli_om += 1
                if n_hv != ref:
                    farkli_hv += 1
                if ref > 0 and n_es / ref < en_kotu[0]:
                    en_kotu = (n_es / ref, d['hedef'], k['kutu'], n_es, ref)
        tum &= _ok(yaz, 'read_engine.py == brute_force.py (%d durum)' % denenen,
                   farkli_om == 0 and denenen > 0,
                   'birebir' if farkli_om == 0 else '%d FARK' % farkli_om)
        tum &= _ok(yaz, 'hizli numpy yolu == kaba kuvvet (%d durum)' % denenen,
                   farkli_hv == 0 and denenen > 0,
                   'birebir' if farkli_hv == 0 else '%d FARK' % farkli_hv)
        if en_kotu[1]:
            _ok(yaz, 'panelin ESKI motoru (reads.py/Sonda) site kaciriyor', True,
                '%s: %d yerine %d (%%%.0f kayip)'
                % (en_kotu[2], en_kotu[3], en_kotu[4], 100 * (1 - en_kotu[0])))
            yaz('      -> Bu arac o motoru KULLANMIYOR. Panelin numune tabanli')
            yaz('         sayilari secenek (4) ile yeniden olculmelidir.')
    except Exception as e:
        tum &= _ok(yaz, 'uc motor karsilastirmasi', False, str(e)[:70])

    # ---- 4b gevsek olcut (<=3) yolu
    try:
        k = H.kutular()[0]
        d0 = H.panel_oku()[0][0]
        h = N.KutuHavuzu(k['kutu'], k['yol'], n=60)
        o = N.KutuOtorite(k['kutu'], k['yol'], n=60)
        iyi = True
        ayr = ''
        for mm in (1, 2, 3):
            a = h.urun_veren(d0['F'], d0['R'], 60, 400, mm)[0]
            b = o.urun_veren(d0['F'], d0['R'], 60, 400, mm)[0]
            if a != b:
                iyi = False
                ayr = 'mm<=%d: hizli %d, otorite %d' % (mm, a, b)
        tum &= _ok(yaz, 'cift taramasi hizli == otorite (mm<=1, <=2, <=3)', iyi, ayr)
    except Exception as e:
        tum &= _ok(yaz, 'gevsek olcut karsilastirmasi', False, str(e)[:70])

    # ---- 5 ureteciler
    try:
        omurga = 'ACGT' * 20
        pen = list(U.pencereler(omurga))
        bek = sum(max(0, min(25, len(omurga) - i) - 18 + 1) for i in range(len(omurga)))
        tum &= _ok(yaz, 'pencere ureteci HER pozisyon x HER uzunluk',
                   len(pen) == bek, '%d/%d' % (len(pen), bek))
        v = U.arms_varyantlari('ACGTACGTACGTACGTACGT')
        tum &= _ok(yaz, "ARMS varyanti: -2 ve -3 icin dort baz (16-1)",
                   len(v) == 15, '%d varyant' % len(v))
        tum &= _ok(yaz, 'parametre izgarasi 144 hucre',
                   len(list(U.izgara_hucreleri())) == 144)
    except Exception as e:
        tum &= _ok(yaz, 'uretec sinamasi', False, str(e)[:70])

    tum &= yon_sinamasi(yaz)

    yaz('   ' + ('TUM SINAMALAR GECTI.' if tum else 'EN AZ BIR SINAMA DUSTU.'))
    return tum
