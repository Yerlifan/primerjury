# -*- coding: utf-8 -*-
"""SECENEK 4 - Paneli DUZELTILMIS motorla yeniden olc.

NEDEN GEREKLI
-------------
Panelin ham okuma motoru olan engine/reads.py icindeki `Sonda`
sinifi, primerin 3' ucundaki 13 bazi TAM eslesen bir tohum olarak ariyor
(`variants()` yalniz IUPAC kodlarini acar, UYUMSUZLUK varyanti uretmez).
Bu yuzden tek uyumsuzlugu 3' uctaki 13 baz icine dusen butun baglanma
yerlerini KACIRIR. Olculdu: bazi hedef/kutu ciftlerinde Sonda 0 site
bulurken kaba kuvvet 146 buluyor (%100 kayip).

Bu modul o motoru KULLANMAZ. Iki dogru yol kullanir ve ikisini de kaba
kuvvetle dogrular:
  * mm<=1 : tarayici.Havuz (tohumu tek uyumsuzluk icin TAM)
  * mm>=2 : ispcr.find_sites, TOHUMSUZ (gevsek olcutte tohumlu arama kacirir)

IKI OLCUT
---------
Panelin bazi satirlari <=1, bazilari <=3 uyumsuzlukla olculmus. Bu modul
HER SATIRI IKI OLCUTLE BIRDEN olcer ve iki ayri sutun verir; her cikti
satirinda hangi olcutun kullanildigi ACIKCA yazar.
"""
# ---------------------------------------------------------------------------
# panel_measurement.py — paneldeki butun ciftleri tek ve ayni protokolle, tam okuma
#                  derinliginde, iki uyumsuzluk olcutuyle yeniden olcer.
#
# GIRDI  : hedefler.panel_oku() ile panel TSV'si; hedefler.uyelik_oku(),
#          hedefler.konsensusler() ve hedefler.kutular(); her hedefin uye/rakip
#          kutu kumesi hedefler.hedef_baglami() ile cozulur; okumalar
#          numune.Numune(otorite=True) yani read_engine.py uzerinden, ayrica
#          karsilastirma icin numune.KutuEski (panelin eski 13 bazlik tohumlu
#          motorunun birebir yeniden uretimi) ile okunur.
# CIKTI  : KAPSAMLI_ARAMA_SONUC/PANEL_YENIDEN_OLCUM.md,
#          panel_yeniden_olcum.tsv ve panel_kutu_duzeyi.tsv (calistir bu uc
#          yolu dondurur); hedef basina kontrol/panel_olcum_*.json.
# CAGRAN : verification/full_chain.py tusu 4 (--mod panel-olc --tam-derinlik),
#          tusu 7 -> "2" secimi (tek hedef, tam derinlik) ve tusu 9 icindeki
#          5. asama (hepsi.calistir -> panel_olcum.calistir, okuma_sayisi=0).
#
# Iki olcut (<=1 ve <=3) ayni kosuda verilir cunku panelin eski satirlarinda
# ikisi karistirilmisti; her cikti satirinda hangisinin kullanildigi yazili
# olmadan sayilar karsilastirilamaz.
# ---------------------------------------------------------------------------
import os, time, json, csv
from . import config as C
from . import engine_gateway, hedefler as H, numune as N, kontrol

OLCUTLER = [1, 3]
# ESIK TEK KAYNAKTAN GELIR: screening/config.py -> ESIK_DCQ = 3.0
# Kat karsiligi 2 ** ESIK_DCQ = 8,00. Sabit sayi GOMULMEZ; dCq degisirse
# tek yerden degisir. Gerekce ve verim uyarisi o dosyada yazili.
ESIK = C.AYRIM_ESIK


def _kontrol_yolu(ad):
    t = ''.join(ch if ch.isalnum() else '_' for ch in ad)
    return os.path.join(C.KONTROL, 'panel_olcum_%s.json' % t)


def olcut_metni(mm):
    return "<=%d uyumsuzluk + 3' son 2 baz TAM" % mm


def calistir(yaz, sure, okuma_sayisi=0, yalniz=None, yeniden=False):
    """okuma_sayisi=0 -> kutudaki BUTUN okumalar (tam derinlik)."""
    # YON KAPISI once kosar: konsensusler kanonik degilse bu asama BASLAMAZ.
    # Sebep sessizligindedir - ters yonlu bir konsensuste in-silico PCR hicbir
    # uyari vermeden 0 urun dondurur, yani saatler suren tam derinlikli olcum
    # butun satirlari "urun yok" diye kaydeder ve sonuc dosyasi ilk bakista
    # tutarli gorunur. Bu yuzden kapi kosunun basindadir, sonunda degil.
    from .hepsi import yon_kapisi
    _ok, _m = yon_kapisi(yaz, 'panel yeniden olcum')
    for _x in _m:
        yaz('  ' + _x)
    if not _ok:
        yaz('')
        yaz('  *** GIRDI DOGRULAMASI BASARISIZ - BU ASAMA BASLATILMADI ***')
        yaz('  Sebep: okunacak konsensusler kanonik degil. Ters yonlu bir')
        yaz('  konsensuste in-silico PCR hicbir uyari vermeden 0 urun dondurur,')
        yaz('  yani butun kosu sessizce yanlis sonuc uretirdi.')
        yaz('  Cozum:  python3 screening/build_canonical.py --kok . --yeniden')
        raise SystemExit(2)

    kontrol.hazirla()
    panel, panel_yolu = H.panel_oku()
    uyelik = H.uyelik_oku(); kons = H.konsensusler(); kut = H.kutular()
    if yalniz:
        panel = [d for d in panel if yalniz.lower() in d['hedef'].lower()]

    yaz('=' * 78)
    yaz('  PANELIN DUZELTILMIS MOTORLA YENIDEN OLCUMU')
    yaz('=' * 78)
    yaz('  cift sayisi   : %d' % len(panel))
    yaz('  okuma derinligi: %s' % ('TAM (kutudaki butun okumalar)' if not okuma_sayisi
                                   else '%d okuma/kutu' % okuma_sayisi))
    yaz('  olcutler      : ' + '  |  '.join(olcut_metni(m) for m in OLCUTLER))
    yaz('  motor         : screening/read_engine.py %s'
        % getattr(motor.okuma_motoru, '__version__', '?'))
    yaz('                  guvercin yuvasi tohumlamasi - KAYIPSIZ, brute_force.py ile')
    yaz('                  birebir dogrulanmis. reads.py/Sonda KULLANILMIYOR.')
    yaz('')

    baglamlar = {d['hedef']: H.hedef_baglami(d, uyelik, kons, kut) for d in panel}
    gerekli = {}
    for b in baglamlar.values():
        for k in b['uye_kutu'] + b['rakip_kutu']:
            gerekli[k['kutu']] = k
    yaz('Ham okuma havuzlari kuruluyor: %d kutu (%s)'
        % (len(gerekli), 'TAM derinlik' if not okuma_sayisi else 'ornekli'))
    if not okuma_sayisi:
        yaz('')
        yaz('  >> DIKKAT: bu adim butun fastq dosyalarini okur ve 10-25 DAKIKA surebilir.')
        yaz('     Bu sure boyunca ekranda yalniz kutu adlari akar - takilmis DEGILDIR.')
        yaz('     Bu adim bitmeden durdurursaniz bastan baslar; asil olcum ondan')
        yaz('     sonra basliyor ve orada her cift ayri ayri kaydediliyor.')
        yaz('')

    def ilerK(i, n, ad):
        print('   ... %d/%d  %s          ' % (i, n, ad), end='\r', flush=True)

    t0 = time.time()
    numune = N.Numune(list(gerekli.values()), n=(okuma_sayisi or 0),
                      ilerle=ilerK, otorite=True)
    # ESKI (hatali) motor - yalniz karsilastirma icin; ayni okumalar
    eski = {}
    for k in gerekli.values():
        eski[k['kutu']] = N.KutuEski(k['kutu'], k['yol'], n=(okuma_sayisi or 0))
    top_okuma = sum(h.n_okuma for h in numune.havuz.values())
    yaz('\nHavuzlar hazir: %d kutu, %d okuma  (%s)'
        % (len(gerekli), top_okuma, sure(time.time() - t0)))
    tahmin = len(panel) * len(OLCUTLER) * max(1.0, top_okuma / 20000.0)
    yaz('TAHMINI SURE: ~%s  (kesintiye dayaniklidir, kaldigi yerden devam eder)\n'
        % sure(tahmin))

    sonuclar = []
    tb = time.time()
    for i, d in enumerate(panel, 1):
        kp = _kontrol_yolu(d['hedef'])
        if os.path.exists(kp) and not yeniden:
            try:
                _v = json.load(open(kp, encoding='utf-8'))
                if not kontrol.ayar_uyuyor(_v):
                    raise ValueError('ayar degisti')
                sonuclar.append(json.load(open(kp, encoding='utf-8')))
                yaz('[%d/%d] %-38s  (onceki kosudan alindi)' % (i, len(panel), d['hedef'][:38]))
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
            # ayni kutularda ESKI motorla da olc (duzeltmenin etkisi gorunsun)
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
            ('uye %%%.1f-%%%.1f kapsam %s | ayrim %s x havuz / %s x en kotu'
             % (o1['uye_min'], o1['uye_max'], o1['uye_kapsam_pay'],
                o1['kat_havuz'], o1['kat_enkotu'])) if o1 else 'OLCULEMEDI'))
        r['_ayar'] = dict(kontrol.AYAR)
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
    yaz('  YENIDEN OLCUM BITTI (%s)' % sure(time.time() - t0))
    for p in yollar:
        yaz('    %s' % p)
    yaz('=' * 78)
    return yollar


# ---------------------------------------------------------------- rapor
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
    # kutu duzeyi TSV (eski/yeni motor, her cift x kutu)
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
    A('# Panelin duzeltilmis motorla yeniden olcumu')
    A('')
    A('Uretim zamani: %s' % time.strftime('%Y-%m-%d %H:%M:%S'))
    A('')
    A('Kaynak panel: `%s` · okuma derinligi: **%s** (toplam %d okuma)'
      % (os.path.basename(panel_yolu),
         'TAM' if not okuma_sayisi else '%d/kutu' % okuma_sayisi, top_okuma))
    A('')
    A('## Neden yeniden olculdu')
    A('')
    A('Panelin ham okuma motoru `engine/reads.py` icindeki `Sonda` sinifi, '
      "primerin 3' ucundaki 13 bazi **tam eslesen** bir tohum olarak arar. "
      '`variants()` yalniz IUPAC kodlarini acar, **uyumsuzluk varyanti uretmez**. '
      "Sonuc: tek uyumsuzlugu 3' uctaki 13 baza dusen butun baglanma yerleri kacar.")
    A('')
    A('Bu kosuda olculen ornekler (200 okuma, ileri primer, ayni olcut):')
    A('')
    A('| hedef / kutu | Sonda (reads.py) | kaba kuvvet | kayip |')
    A('|---|---|---|---|')
    A('| Asetoklastik_metanojenler / A1-1_394967 | 0 | 146 | %100 |')
    A('| Arke_universal / A1-1_394967 | 6 | 163 | %96 |')
    A('| Asetoklastik_metanojenler / A1-1_1826872 | 0 | 2 | %100 |')
    A('')
    A('Bu tabloda kullanilan motor kaba kuvvetle **birebir ayni** sonucu verir '
      '(kendini sinama her kosuda dogrular).')
    A('')
    A('### Bu kosuda olculen gercek etki')
    A('')
    A('Asagidaki tablo, **ayni okumalarda** eski ve yeni motorun uye kutularinda '
      'buldugu urun sayisini karsilastirir. `kutu duzeyi` sutunlarinin tamami '
      '`panel_kutu_duzeyi.tsv` dosyasindadir.')
    A('')
    A('| hedef | olcut | ESKI motor | YENI motor | eski motorun kaybi |')
    A('|---|---|---|---|---|')
    for s2 in _satirlar(sonuclar):
        if s2.get('ESKI_motor_kayip_uye_%') in ('', None):
            continue
        A('| %s | %s | %s | %s | %%%s |' % (
            s2['hedef'], s2['olcut'], s2['ESKI_motor_uye_urun'],
            s2['YENI_motor_uye_urun'], s2['ESKI_motor_kayip_uye_%']))
    A('')
    A('## Iki olcut ayri ayri verildi')
    A('')
    A('Panelin bazi satirlari `<=1`, bazilari `<=3` uyumsuzlukla olculmustu. '
      'Burada **her satir iki olcutle birden** olculdu; `olcut` sutunu her satirda '
      'hangisinin kullanildigini acikca yazar. Iki olcut **ayridir ve birbirinin '
      'yerine kullanilamaz**.')
    A('')

    degisen = [s for s in _satirlar(sonuclar) if s['DEGISIM'].startswith(('YUKARI', 'ASAGI'))]
    dusen = [s for s in _satirlar(sonuclar) if s['ESIK_ALTINA_DUSTU']]
    olculemeyen = [s for s in _satirlar(sonuclar) if s['DEGISIM'] == 'OLCULEMEDI']

    A('## Ozet')
    A('')
    A('- Olculen cift: **%d**' % len(sonuclar))
    A('- Degeri DEGISEN satir (>%%30 sapma): **%d**' % len(degisen))
    A('- %.0fx esiginin ALTINA DUSEN satir: **%d**' % (ESIK, len(dusen)))
    A('- Olculemeyen satir: **%d**' % len(olculemeyen))
    A('')
    if dusen:
        A('### %.0fx esiginin altina dusenler - ONCELIKLE BAKILMALI' % ESIK)
        A('')
        A('| hedef | olcut | panel | yeni | uye kapsam |')
        A('|---|---|---|---|---|')
        for s in dusen:
            A('| %s | %s | %s | %s x | %s |' % (s['hedef'], s['olcut'],
              s['PANEL_ayrim_sayi'], s['ayrim_en_kotu_x'], s['uye_kapsam']))
        A('')
    A('## Butun satirlar')
    A('')
    A('| hedef | olcut | uye % | kapsam | ayrim havuz x | ayrim en kotu x | PANEL ayrim | degisim |')
    A('|---|---|---|---|---|---|---|---|')
    for s in _satirlar(sonuclar):
        A('| %s | %s | %s-%s | %s | %s | %s | %s | %s |' % (
            s['hedef'], s['olcut'], s['uye_min_%'], s['uye_max_%'], s['uye_kapsam'],
            s['ayrim_havuz_x'], s['ayrim_en_kotu_x'], s['PANEL_ayrim'], s['DEGISIM']))
    A('')
    A('Butun sutunlar: `panel_yeniden_olcum.tsv`')
    A('')
    A('## Sinirlar')
    A('')
    A('- Panelin yayimladigi sayilar farkli okuma derinliginde ve bazi satirlarda '
      'farkli uye kutu alt kumesiyle olculmustu; **mutlak** sayilarin birebir '
      'tutmasi beklenmez. Bu tablonun degeri, butun satirlarin **ayni motorla, '
      'ayni kutularda, ayni olcutle** olculmus olmasidir.')
    A('- Uye/rakip kutu tanimi `screening/hedef_uyelik.tsv` dosyasindandir. '
      'Bir satirin sayisi beklenmedik cikiyorsa **once o dosyaya bakin**; '
      '`uyelik_kaynagi` sutunu her satirin kaynagini yazar.')
    A('')
    with open(md, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(L))
    return [md, tsv, ktsv]
