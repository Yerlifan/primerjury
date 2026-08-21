# -*- coding: utf-8 -*-
"""SECENEK 5 - UYELIK TANIMI DENETIMI ve DUYARLILIK ANALIZI.

NEDEN GEREKLI
-------------
Bir hedefin "ayrim kati" sayisi, hangi kutularin UYE hangilerinin RAKIP
sayildigina dogrudan baglidir. Uyelik tanimi degisince sayi degisir - motor
degismese bile. Olculdu:

  Proteolitik_Cloacimonas, grup satiri (3 takson) ile  ->  ayrim  0,0x
  Proteolitik_Cloacimonas, tek uye (456827) ile        ->  ayrim 23,5x
  panelde yazan                                        ->        23,0x

Yani sayinin kendisi degil, TANIM yanlisti. Bu modul her hedef icin uyeligin
BUTUN makul tanimlarini yan yana olcer, boylece hangi sayinin hangi tanima
bagli oldugu gorunur olur.

TANIM KAYNAKLARI
----------------
  A. hedef_uyelik.tsv       aracin su an kullandigi tanim (PANEL / TURETILDI)
  B. hedefler.tsv           projenin karar tablosu (grup satiri)
  C. tek uye                yalniz en cok urun veren tek takson
  D. ciftler.tsv            diger oturumun (okuma motoru duzeltmesi) tanimi
  E. olculen kimlik         hedef_kimlik.tsv'deki OLCULEN organizmadan turetilmis

TANI (Proteiniphilum gibi durumlar)
-----------------------------------
Bir hedefin uye kutularinda hic urun yoksa iki ihtimal vardir:
  (1) uye kumesi yanlis  -> hangi kutu(lar) urun veriyor, tek tek bulunur
  (2) konsensus ile ham okumalar uyusmuyor -> ayni cift KONSENSUSTE olculur;
      konsensuste urun verip okumalarda vermiyorsa sorun uyelik degil,
      konsensus/okuma uyusmazligidir
Modul IKISINI DE sinar ve hangisi oldugunu yazar.
"""
# ---------------------------------------------------------------------------
# uyelik_denetimi.py — her panel hedefinin uyelik tanimlarini (hangi kutu uye,
#                      hangisi rakip) yan yana olcer ve sayinin tanima ne kadar
#                      duyarli oldugunu gosterir.
#
# GIRDI  : hedefler.panel_oku() ile panel TSV'si; hedefler.kutular() ile fastq
#          kutulari; hedefler.konsensusler() ile kanonik konsensusler;
#          hedefler.acik_uyelik() ile KAPSAMLI_ARAMA/hedef_uyelik.tsv;
#          hedefler.uyelik_oku() ile WSL_betikleri/hedefler.tsv;
#          KAPSAMLI_ARAMA/ciftler.tsv (ya da eski/ciftler.tsv);
#          primer_final/hedef_kimlik.tsv; olcumu numune.Numune yapar.
# CIKTI  : KAPSAMLI_ARAMA_SONUC/UYELIK_DENETIMI.md ve uyelik_duyarlilik.tsv
#          (calistir bu iki yolu liste olarak dondurur); ayrica hedef basina
#          kontrol/uyelik_*.json kontrol noktasi.
# CAGRAN : KAPSAMLI_ARAMA.bat tusu 5 (--mod uyelik), tusu 7 -> "3" secimi
#          (tek hedefin uyelik denetimi) ve tusu 9 icindeki 6. asama
#          (hepsi.calistir -> uyelik_denetimi.calistir).
#
# BU MODUL UYELIGI DEGISTIRMEZ. Butun makul tanimlari olcer, yan yana koyar ve
# raporun sonuna "Bu arac uyelik tanimini kendiliginden degistirmez" diye yazar;
# hedef_uyelik.tsv'ye YAZMAZ. Ilke: kanit yoklugu kanit sayilmaz. Bir sayinin
# beklenenden farkli cikmasi, kutunun yerinin degistirilmesi icin kanit degildir;
# yer degistirme ancak pozitif olcum kanitiyla ve elle yapilir. Bu yuzden modulun
# ciktisi bir karar degil, bir secenek tablosudur.
# ---------------------------------------------------------------------------
import os, csv, json, time, re
from . import yapilandirma as C
from . import motor, hedefler as H, numune as N, kontrol

# ESIK TEK KAYNAKTAN GELIR: KAPSAMLI_ARAMA/yapilandirma.py -> ESIK_DCQ = 3.0
# Kat karsiligi 2 ** ESIK_DCQ = 8,00. Sabit sayi GOMULMEZ; dCq degisirse
# tek yerden degisir. Gerekce ve verim uyarisi o dosyada yazili.
ESIK = C.AYRIM_ESIK
PAKET = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------- tanim kaynaklari
def _ciftler_tsv_uyelik():
    """eski/ciftler.tsv (diger oturumun tanimi) -> hedef -> taxid listesi."""
    out = {}
    for aday in (os.path.join(PAKET, 'ciftler.tsv'),
                 os.path.join(PAKET, 'eski', 'ciftler.tsv')):
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
    """hedef_kimlik.tsv'deki OLCULEN kimlikten taxid kumesi turet."""
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
    """Bir panel hedefi icin butun uyelik tanimlarini uret."""
    ad = satir['hedef']
    sf = [x.strip() for x in (satir['sinif'] or '').split('/') if x.strip()]
    out = []

    a = acik.get(ad)
    if a:
        out.append(('A. hedef_uyelik.tsv (%s)' % (a.get('kaynak') or '?'),
                    a['uye'], a['haric']))

    anahtar = H.AD_ESLEME.get(ad, ad)
    g = grup_uyelik.get(anahtar)
    if g:
        out.append(('B. hedefler.tsv:%s' % anahtar, g['uye'], g['haric']))

    c = cift_uyelik.get(ad)
    if c:
        out.append(('D. ciftler.tsv (diger oturum)', c, []))

    k = kimlik_uyelik.get(ad)
    if k:
        out.append(('E. olculen kimlik (hedef_kimlik.tsv)', k, []))

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
    """Uye kutusunda urun yoksa sebebini bul: uyelik mi, konsensus/okuma mi."""
    # Iki soru ayri ayri sorulur cunku cevaplari farkli iki duzeltmeye isaret
    # eder. (1) "Ayni sinifin HANGI kutulari urun veriyor" - urun veren kutu
    # uye listesinde degilse sorun uyelik tanimindadir. (2) "Ayni cift KONSENSUS
    # dizisinde urun veriyor mu" - konsensuste veriyor ama ham okumalarda
    # vermiyorsa sorun uyelik degil, konsensusun okumalari temsil etmemesidir
    # ve cozum konsensusu yeniden uretmektir. Ikisi karistirilirsa yanlis dosya
    # duzeltilir.
    F, R = satir['F'], satir['R']
    lo, hi = C.URUN_IDEAL[0], C.URUN_MUTLAK_UST
    taxad = H.taxid_adlari()

    # (1) HANGI kutular urun veriyor - butun sinif taranir
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
        for s in (k['dizi'], motor.rc(k['dizi'])):
            pr = motor.amplify(s, F, R, max_mm=1, lo=lo, hi=hi)
            if pr:
                bp = min(pr, key=lambda x: x[3] + x[4])[2]
                kons_veren.append(dict(kutu=k['kutu'], taxid=k['taxid'],
                                       ad=taxad.get(k['taxid'], '?'), boy=bp))
                break
    return veren, kons_veren


def tani_yorumu(uye_kutu, veren, kons_veren, panel_urun, olcum=None, panel_uye=''):
    """Kapsam sifirsa sebebini bul: uyelik mi, konsensus/okuma uyusmazligi mi.

    "Urun var" ile "kapsandi" ayni sey DEGILDIR: bir uye kutusu %2 urun
    veriyorsa urun VARDIR ama kapsam esigini gecmez. Panelin o satir icin
    bildirdigi deger cok daha yuksekse bu ayrica aciklanmasi gereken bir
    celiskidir - modul bunu ayri bir tani olarak bildirir.
    """
    # Tani sirasi ONEMLI ve daralan bir eleme olarak yazilmistir: once kapsam
    # var mi, sonra uye kutusunda zayif da olsa urun var mi, sonra sinifin baska
    # kutusunda var mi, en sonda konsensusta var mi. Her basamak bir onceki
    # ihtimali dislar; sirasi degisirse ayni durum farkli teshis alir.
    uye_ad = {k['kutu'] for k in uye_kutu}
    uye_veren = [v for v in veren if v['kutu'] in uye_ad]
    kons_uye = [v for v in kons_veren if v['kutu'] in uye_ad]
    kapsam = (olcum or {}).get('uye_kapsam', 0)
    esik = int(100 * C.KAPSAM_ESIGI)

    if kapsam:
        return 'SORUN YOK', 'Uye kutularinin %d tanesi kapsam esigini (>=%%%d) geciyor.' % (
            kapsam, esik)

    if uye_veren:
        en = max(uye_veren, key=lambda v: v['yuzde'])
        dis = [v for v in veren if v['kutu'] not in uye_ad]
        ek = ''
        if dis:
            d0 = max(dis, key=lambda v: v['yuzde'])
            ek = (' Ayni sinifta uye OLMAYAN %s (%s) %%%s veriyor - uye kumesi '
                  'genisletilmeli mi?' % (d0['kutu'], d0['ad'], d0['yuzde']))
        return ('KAPSAM ESIGININ ALTINDA - PANELLE CELISKILI',
                'Uye kutularinda urun VAR ama zayif: en iyisi %s %%%s (esik >=%%%d), '
                'urun boyu %s bp. Panelde bu satir icin "%s" yaziyor. Urun boyu '
                'tutuyorsa cift dogru; farkli olan ORAN. Muhtemel sebepler: '
                '(a) panel farkli bir olcut ya da farkli gevseklik ayariyla '
                'olculmus, (b) konsensus ile ham okumalar ayni organizmayi '
                'anlatmiyor. Konsensus sinamasi asagida.%s'
                % (en['kutu'], en['yuzde'], esik,
                   en['boy'][0][0] if en['boy'] else '?', panel_uye, ek))

    if veren:
        ilk = veren[0]
        return ('UYE KUMESI YANLIS OLABILIR',
                'Uye kutularinda urun YOK, ama ayni sinifta %d kutu urun veriyor. '
                'En yuksek: %s (%s) %%%s, urun %s. Bu takson uye sayilmali mi?'
                % (len(veren), ilk['kutu'], ilk['ad'], ilk['yuzde'],
                   ilk['boy'][0][0] if ilk['boy'] else '?'))
    if kons_uye:
        b = kons_uye[0]
        return ('KONSENSUS/OKUMA UYUSMAZLIGI',
                'Uye kutusunun KONSENSUSU urun veriyor (%s, %s bp - panelde %s bp) '
                'ama ayni kutunun HAM OKUMALARI vermiyor. Sorun uyelik degil: '
                'konsensus ile okumalar ayni organizmayi anlatmiyor ya da '
                'konsensus bayat. Secenek (6) ile konsensusu yeniden uretin.'
                % (b['kutu'], b['boy'], panel_urun))
    return ('CIFT HIC URUN VERMIYOR',
            'Ne uye kutularinda, ne ayni sinifin herhangi bir kutusunda, ne de '
            'konsensuslerde urun olusuyor. Cift dizisi ya da urun boyu penceresi '
            'kontrol edilmeli.')


# ---------------------------------------------------------------- ana
def calistir(yaz, sure, okuma_sayisi=C.NUMUNE_OKUMA_SAYISI, yalniz=None,
             yeniden=False):
    from .hepsi import yon_kapisi
    _ok, _m = yon_kapisi(yaz, 'uyelik denetimi')
    for _x in _m:
        yaz('  ' + _x)
    if not _ok:
        yaz('')
        yaz('  *** GIRDI DOGRULAMASI BASARISIZ - BU ASAMA BASLATILMADI ***')
        yaz('  Sebep: okunacak konsensusler kanonik degil. Ters yonlu bir')
        yaz('  konsensuste in-silico PCR hicbir uyari vermeden 0 urun dondurur,')
        yaz('  yani butun kosu sessizce yanlis sonuc uretirdi.')
        yaz('  Cozum:  python3 KAPSAMLI_ARAMA/kanonik_uret.py --kok . --yeniden')
        raise SystemExit(2)

    kontrol.hazirla()
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
    yaz('  UYELIK TANIMI DENETIMI ve DUYARLILIK ANALIZI')
    yaz('=' * 78)
    yaz('  hedef sayisi        : %d' % len(panel))
    yaz('  tanim kaynaklari    : hedef_uyelik.tsv, hedefler.tsv, ciftler.tsv,')
    yaz('                        olculen kimlik, tek uye')
    yaz('  TURETILDI isaretli  : %d satir  (ozellikle kontrol edilmeli)' % len(turetildi))
    for t in turetildi:
        yaz('        - %s' % t)
    yaz('')

    gerekli = {k['kutu']: k for k in kut}
    yaz('Ham okuma havuzlari kuruluyor: %d kutu' % len(gerekli))
    yaz('  >> Bu adim birkac dakika surer; ekranda kutu adlari akar, takilmis degildir.')

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
                if not kontrol.ayar_uyuyor(_v):
                    raise ValueError('ayar degisti')
                sonuclar.append(json.load(open(kp, encoding='utf-8')))
                yaz('[%d/%d] %-38s (onceki kosudan)' % (i, len(panel), d['hedef'][:38]))
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
            yaz('        %-42s uye %2d kutu  kapsam %-6s ayrim %s x / %s x'
                % (etiket[:42], len(uk), (o or {}).get('uye_kapsam_pay', '-'),
                   (o or {}).get('kat_havuz'), (o or {}).get('kat_enkotu')))

        # C. tek uye - hangi tek takson en iyi sonucu veriyor
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
                yaz('        %-42s uye %2d kutu  kapsam %-6s ayrim %s x / %s x'
                    % (('C. tek uye (%s)' % en[1])[:42], en[3],
                       en[2].get('uye_kapsam_pay', '-'),
                       en[2].get('kat_havuz'), en[2].get('kat_enkotu')))

        # tani: A tanimindaki uye kutularinda urun var mi
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
        r['_ayar'] = dict(kontrol.AYAR)
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
    yaz('  UYELIK DENETIMI BITTI (%s)' % sure(time.time() - t0))
    for p in yollar:
        yaz('    %s' % p)
    yaz('=' * 78)
    return yollar


# ---------------------------------------------------------------- rapor
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
    A('# Uyelik tanimi denetimi ve duyarlilik analizi')
    A('')
    A('Uretim zamani: %s' % time.strftime('%Y-%m-%d %H:%M:%S'))
    A('')
    A('Kaynak panel: `%s`' % os.path.basename(panel_yolu))
    A('')
    A('## Bu rapor neyi gosteriyor')
    A('')
    A('Bir hedefin **ayrim kati** sayisi, hangi kutularin UYE hangilerinin RAKIP '
      'sayildigina dogrudan baglidir. Uyelik tanimi degisince sayi degisir - '
      'olcum motoru hic degismese bile.')
    A('')
    A('Bu koşuda olculen ornek:')
    A('')
    A('| tanim | ayrim |')
    A('|---|---|')
    A('| Proteolitik_Cloacimonas, grup satiri (3 takson) | 0,0x |')
    A('| Proteolitik_Cloacimonas, tek uye (456827) | 23,5x |')
    A('| **panelde yazan** | **23,0x** |')
    A('')
    A('Yani sayinin kendisi degil, **tanim** yanlisti. Asagidaki tablolar her hedef '
      'icin uyeligin butun makul tanimlarini **yan yana** olcer.')
    A('')
    A('### Tanim kaynaklari')
    A('')
    A('| kod | kaynak |')
    A('|---|---|')
    A('| **A** | `KAPSAMLI_ARAMA/hedef_uyelik.tsv` — aracin su an kullandigi tanim |')
    A('| **B** | `WSL_betikleri/hedefler.tsv` — projenin karar tablosu (grup satiri) |')
    A('| **C** | tek uye — en iyi sonucu veren tek takson |')
    A('| **D** | `ciftler.tsv` — okuma motoru duzeltmesi oturumunun tanimi |')
    A('| **E** | `primer_final/hedef_kimlik.tsv` — OLCULEN kimlikten turetilmis |')
    A('')

    A('## ONCE BUNLARA BAKIN — `TURETILDI` isaretli satirlar')
    A('')
    A('Bu satirlarin uyelik tanimi panelde acikca yazili degildi; hedef adindan ve '
      '`taxid_adlari.tsv`\'den **cikarildi**. Yanlis olma ihtimali en yuksek olanlar '
      'bunlardir.')
    A('')
    A('| hedef | tanim degisince ayrim kac kat oynuyor | tani |')
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
        A('## Tani gereken hedefler')
        A('')
        for r in sorunlu:
            A('### %s — %s' % (r['hedef'], r['tani']))
            A('')
            A(r['tani_aciklama'])
            A('')
            if r['urun_veren_kutular']:
                A('Bu siniftaki kutulardan urun verenler:')
                A('')
                A('| kutu | takson | urun/okuma | %% | boy |')
                A('|---|---|---|---|---|')
                for v in r['urun_veren_kutular']:
                    A('| %s | %s | %s/%s | %s | %s |' % (
                        v['kutu'], v['ad'], v['urun'], v['okuma'], v['yuzde'],
                        v['boy'][0][0] if v['boy'] else '-'))
                A('')
            if r['konsensus_veren']:
                A('Konsensuste urun veren kutular (ham okumalarda vermiyor):')
                A('')
                A('| kutu | takson | konsensuste urun boyu |')
                A('|---|---|---|')
                for v in r['konsensus_veren']:
                    A('| %s | %s | %s bp |' % (v['kutu'], v['ad'], v['boy']))
                A('')

    A('## Butun hedefler — tanim duyarliligi')
    A('')
    for r in sonuclar:
        A('### %s%s' % (r['hedef'], '  *(TURETILDI)*' if r['turetildi_mi'] else ''))
        A('')
        A('Panelde yazan: `%s` · urun %s bp%s'
          % (r['panel_ayrim'], r['urun_panel'],
             ' · **tanim degisince ayrim %sx oynuyor**' % r['oynaklik']
             if r['oynaklik'] and r['oynaklik'] > 1.5 else ''))
        A('')
        A('| tanim | uye kutu | kapsam | uye %% | ayrim havuz x | ayrim en kotu x |')
        A('|---|---|---|---|---|---|')
        for v in r['varyantlar']:
            A('| %s%s | %s | %s | %s-%s | %s | %s |' % (
                v['tanim'], ' *(A ile ayni kume)*' if v['ayni_kume'] else '',
                v['uye_kutu'], v['kapsam'], v['uye_min'], v['uye_max'],
                v['kat_havuz'], v['kat_enkotu']))
        A('')

    A('## Ne yapmali')
    A('')
    A('1. `TURETILDI` isaretli satirlari gozden gecirin; dogru tanimi '
      '`KAPSAMLI_ARAMA/hedef_uyelik.tsv` dosyasina yazin.')
    A('2. `oynaklik` sutunu yuksek olan hedeflerde yayimlanan sayi **tanima cok '
      'duyarlidir** — hangi tanimla bildirildigi panelde acikca yazilmalidir.')
    A('3. `TANI` sutunu `KONSENSUS/OKUMA UYUSMAZLIGI` diyorsa once secenek (6) ile '
      'o kutunun konsensusunu yeniden uretin, sonra olcumu tekrarlayin.')
    A('')
    A('> Bu arac uyelik tanimini **kendiliginden degistirmez**. Olcer, secenekleri '
      'yan yana koyar, karari size birakir.')
    A('')
    with open(md, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(L))
    return [md, tsv]
