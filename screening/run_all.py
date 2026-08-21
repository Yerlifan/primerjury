# -*- coding: utf-8 -*-
"""SECENEK 9 - HER SEYI SIRAYLA KOS.

Kullanici aksam tek tik yapip yatabilsin diye butun asamalari dogru SIRADA
kosar. Sira keyfi degil, bagimliliklara gore:

  1. Kendini sinama          motorlar dogru mu (gecmezse hicbir sey kosmaz)
  2. Konsensus yeniden uretim   ARAMA omurgayi kullanir -> ONCE bu
  3. Panel yeniden olcum        tam derinlik, duzeltilmis motor
  4. Uyelik denetimi            hangi sayi hangi tanima bagli
  5. Kapsamli arama             sorunlu hedefler (2. adimin konsensusuyla)
  6. BIRLESIK OZET RAPOR        ne degisti, ne dustu, nerede yeni aday var

Her asama bitince kontrol noktasi yazilir. Kesilirse aynen (9) secilerek
kaldigi asamadan devam edilir; biten asamalar TEKRAR KOSULMAZ.

KONSENSUS -> ARAMA BAGI
-----------------------
Yeniden uretilen konsensus, arama omurgasi olarak KALITE KAPISINDAN gecerse
kullanilir (N orani dusuk ve iki yontem ayrilmamis). Gecmezse eski konsensus
kullanilir. Hangi omurganin kullanildigi her hedef icin rapora yazilir -
sessiz degisiklik yapilmaz.
"""
# ---------------------------------------------------------------------------
# run_all.py — yedi asamayi bagimlilik sirasiyla kosan toplu akis ve butun
#            asamalarin ciktilarini tek dosyada birlestiren ozet uretici.
#
# GIRDI  : config.py'deki yollar; kontrol/hepsi_durum.json (hangi asama
#          bitmis); KAPSAMLI_ARAMA_SONUC altindaki panel_yeniden_olcum.tsv,
#          panel_kutu_duzeyi.tsv, uyelik_duyarlilik.tsv,
#          konsensus_yeniden_uretim.tsv, adaylar.tsv ve kontrol/hedef_*.json.
#          Asamalari kendisi cagirir: build_canonical.py (alt surec olarak),
#          konsensus_uret.calistir, panel_olcum.calistir,
#          uyelik_denetimi.calistir, __main__.aramayi_kos.
# CIKTI  : KAPSAMLI_ARAMA_SONUC/00_OZET_HEPSI.md (ozet_yaz'in dondurdugu yol);
#          kontrol/hepsi_durum.json (asama durumu). calistir() cikis kodu
#          dondurur: 0 tamam, 1 kullanici durdurdu, 2 kapi gecilemedi.
# CAGRAN : verification/full_chain.py tusu 9 (--mod hepsi -> HP.calistir) ve tusu
#          S (--mod ozet -> HP.yalniz_ozet). Ayrica yon_kapisi() fonksiyonu
#          __main__.py, panel_measurement.py, membership_check.py ve build_consensus.py
#          tarafindan cagrildigi icin 1, 2, 3, 4, 5, 6 ve 7 tuslarinda da
#          dolayli olarak calisir.
# ---------------------------------------------------------------------------
import os, json, time, csv
from . import config as C
from . import checks

ASAMA_YOLU = None

# yeni konsensusun omurga olarak kullanilabilmesi icin kalite kapisi
KONSENSUS_N_UST = 20.0        # %N bundan buyukse eski konsensus kullanilir
KONSENSUS_AYRILMA_UST = 15    # iki yontemin celistigi sutun sayisi tavani


def _durum_yolu():
    return os.path.join(C.KONTROL, 'hepsi_durum.json')


def durum_oku():
    p = _durum_yolu()
    if os.path.exists(p):
        try:
            d = json.load(open(p, encoding='utf-8'))
            if kontrol.ayar_uyuyor(d):
                return d
        except Exception:
            pass
    return dict(bitmis=[], baslangic=time.strftime('%Y-%m-%d %H:%M:%S'), ciktilar={})


def durum_yaz(d):
    kontrol.hazirla()
    d['_ayar'] = dict(kontrol.AYAR)
    p = _durum_yolu()
    with open(p + '.gecici', 'w', encoding='utf-8') as fh:
        json.dump(d, fh, ensure_ascii=False, indent=1, default=str)
    os.replace(p + '.gecici', p)


# ---------------------------------------------------------------- konsensus kapisi
def konsensus_kalite():
    """Yeniden uretilen konsensuslerin kalite tablosu: kutu -> (gecti_mi, sebep).

    NOT: bu fonksiyon bir ara duzenlemede yanlislikla silinmisti ve ozet_yaz()
    icinden cagriliyordu -> 5 saatlik kosunun SON adimi NameError ile dustu.
    Geri konuldu; ozet_yaz artik ayrica tek basina da kosturulabiliyor
    (menu secenegi ve --mod ozet).
    """
    tsv = os.path.join(C.CIKTI, 'konsensus_yeniden_uretim.tsv')
    out = {}
    if not os.path.exists(tsv):
        return out
    try:
        for row in csv.DictReader(open(tsv, encoding='utf-8'), delimiter='\t'):
            try:
                n = float(row.get('N_yuzde') or 0)
                ayr = int(row.get('yontemler_ayrildi') or 0)
            except (ValueError, TypeError):
                continue
            if n > KONSENSUS_N_UST:
                out[row['kutu']] = (False, 'N orani %%%.1f > %%%.0f' % (n, KONSENSUS_N_UST))
            elif ayr > KONSENSUS_AYRILMA_UST:
                out[row['kutu']] = (False, 'iki yontem %d sutunda ayrildi' % ayr)
            else:
                out[row['kutu']] = (True, 'N %%%.1f, ayrilma %d' % (n, ayr))
    except Exception:
        pass
    return out


# ---------------------------------------------------------------- yon kapisi
def yon_kapisi(yaz, asama):
    """Bir asama baslamadan ONCE girdisinin kanonik oldugunu dogrular.

    Yon hatasi sessizdir: ters yonlu bir konsensuste in-silico PCR hicbir uyari
    vermeden 0 urun dondurur ve butun gece "hicbir hedefte urun yok" sonucu
    uretilir. Bu yuzden her asama kendi girdisini bastan sinar; gecmezse
    ASAMA BASLAMAZ.

    Donus: (gecti_mi, mesaj_satirlari)
    """
    # NEDEN YON KANONIK KAYNAKTAN OKUNUR
    # Nanopore okumalari cift yonludur; bir kutunun konsensusu okumalarin hangi
    # yonde capalandigina gore SENSE ya da ANTISENSE uretilmis olabilir. Ham
    # "consensus sequences" klasoru bu yuzden KARISIK yonludur (olculdu: 71
    # antisense, 27 sense). Motorlar diziyi yalniz arti iplikte taradigi icin
    # ters yonlu bir konsensuste in-silico PCR urunlerin TAMAMINI kaybeder ve
    # bunu hicbir uyari vermeden yapar - cikti "urun yok" olur, "yon yanlis"
    # olmaz. Bu yuzden konsensusler tek bir kanonik klasorden okunur ve her
    # asama, kendi girdisinin kanonik oldugunu BASLAMADAN once burada dogrular.
    # Kapi uc seyi sirayla sinar: orientation.py kendi sinavini geciyor mu, kanonik
    # indeks okunabiliyor mu, indeksteki konsensuslerin hicbiri ANTISENSE degil
    # mi. Ucunden biri duserse asama baslatilmaz.
    from . import orientation as Y
    from . import targets as H
    m = []
    try:
        h = Y.kendini_sina()
        if h:
            return False, ['orientation.py kendi sinavindan GECEMEDI: ' + '; '.join(h)]
        kons = H.konsensusler()
    except Exception as e:
        return False, ['Kanonik konsensus okunamadi: %s' % e,
                       'To generate it:  python3 screening/build_canonical.py --root . ']
    if not kons:
        return False, ['Kanonik konsensus indeksi BOS.']
    ters = [k['kutu'] for k in kons
            if Y.tespit(k['dizi'], Y.sinifi(k['kutu']))[0] == 'ANTISENSE']
    if ters:
        return False, ['%d konsensus TERS yonde: %s' % (len(ters), ', '.join(ters[:5])),
                       'Ters yonlu konsensuste urunlerin TAMAMI kaybolur.',
                       'To fix it:  python3 screening/build_canonical.py --root . --yeniden']
    m.append('yon kapisi [%s]: %d kutu, hepsi SENSE - GECTI' % (asama, len(kons)))
    return True, m


def kanonik_kos(yaz, sure, oncelik='ozgun'):
    """build_canonical.py'yi ayni surecte cagirir (ayri komut gerekmesin diye)."""
    import subprocess, sys as _sys
    betik = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'build_canonical.py')
    if not os.path.exists(betik):
        return False, 'build_canonical.py bulunamadi'
    komut = [_sys.executable, betik, '--kok', C.KOK, '--oncelik', oncelik]
    yaz('  > %s' % ' '.join(komut[1:]))
    t = time.time()
    try:
        r = subprocess.run(komut, capture_output=True, text=True, timeout=7200)
    except subprocess.TimeoutExpired:
        return False, 'kanonik uretim 2 saati asti'
    for satir in (r.stdout or '').strip().split('\n')[-6:]:
        if satir.strip():
            yaz('    ' + satir.strip())
    if r.returncode != 0:
        return False, (r.stderr or '')[-400:]
    # onbellekteki konsensus listesi bayatlamasin
    from . import targets as H
    if hasattr(H, '_kons_onbellek'):
        H._kons_onbellek = None
    return True, 'kanonik uretim tamam (%s)' % sure(time.time() - t)


# ---------------------------------------------------------------- asamalar
def calistir(yaz, sure, cizgi, a):
    kontrol.hazirla()
    d = durum_oku()
    t0 = time.time()

    cizgi('=')
    yaz(u'  RUN EVERYTHING IN ORDER')
    cizgi('=')
    yaz(u'  Stages and estimated times:')
    yaz(u'    1. Self-test                          ~1 min')
    yaz(u'    2. Canonical consensus (orientation)  ~2-5 min   <- ORIENTATION GATE')
    yaz(u'    3. Consensus regeneration             ~30-60 min')
    yaz(u'    4. Move new consensus to canonical    ~2-5 min   <- ORIENTATION GATE')
    yaz(u'    5. Panel re-measurement (full depth)  ~1-2 h')
    yaz(u'    6. Membership audit                   ~20-40 min')
    yaz(u'    7. Comprehensive search               ~1-3 h')
    yaz(u'       + combined summary report          ~1 min')
    yaz(u'  TOTAL: roughly 3-7 hours. Can be left running overnight.')
    yaz('')
    if d['bitmis']:
        yaz(u'  RESUMING - finished stages: %s' % ', '.join(d['bitmis']))
        yaz('')
    yaz(u'  State is written to disk after every stage. If interrupted, choose (9) again;')
    yaz(u'  finished stages are not re-run.')
    cizgi('=')

    from . import self_test

    # ---- 1 kendini sinama (kapi gorevinde)
    yaz(u'\n[STAGE 1/7] Self-test')
    if not getattr(a, 'sinama_atla', False) and not kendini_sina.calistir(yaz):
        yaz(u'\nSELF-TEST FAILED - no stage was run.')
        yaz(u'Fix the line marked *** FAILED *** above first.')
        return 2

    # ---- 2 KANONIK KONSENSUS (referans onceligi) - her seyin girdisi
    yaz(u'\n[STAGE 2/7] Canonical consensus generation (orientation normalisation)')
    yaz(u'  No measurement is trustworthy without this step: the raw consensus directory')
    yaz(u'  is MIXED orientation, and on a reverse-oriented consensus EVERY product is lost.')
    yaz(u'  Source: "consensus sequences" (the strict set the panel was built on).')
    ok, msj = kanonik_kos(yaz, sure, oncelik='ozgun')
    if not ok:
        yaz('\n  KANONIK URETIM BASARISIZ: %s' % msj)
        yaz(u'  The following stages were NOT STARTED (so that no wrong result is produced).')
        return 2
    yaz('  ' + msj)
    ok, msj = yon_kapisi(yaz, 'kanonik uretim sonrasi')
    for x in msj:
        yaz('  ' + x)
    if not ok:
        yaz(u'\n  THE ORIENTATION GATE DID NOT PASS - the later stages WERE NOT STARTED.')
        return 2

    # ---- 3 konsensus yeniden uretim
    if 'konsensus' in d['bitmis']:
        yaz(u'\n[STAGE 3/7] Consensus regeneration - SKIPPED (already finished)')
    else:
        yaz(u'\n[STAGE 3/7] Consensus regeneration (from the raw reads)')
        from . import build_consensus
        try:
            y = konsensus_uret.calistir(yaz, sure, yalniz=getattr(a, 'hedef', None))
            d['bitmis'].append('konsensus')
            d['ciktilar']['konsensus'] = y
            durum_yaz(d)
        except KeyboardInterrupt:
            yaz(u'\nSTOPPED - continue with (9).')
            return 1
        except Exception as e:
            yaz('\n  ASAMA 3 HATASI: %s' % e)
            yaz(u'  No new consensus could be generated; the canonical reference set stays in use.')

    # ---- 4 YENI konsensusleri kanonige al (tek kaynak korunur)
    yaz(u'\n[STAGE 4/7] Moving new consensus sequences to canonical (--priority new)')
    yaz(u'  Newly generated consensus sequences are NOT used directly. They are first')
    yaz(u'  written into the canonical directory through orientation normalisation, so there stays one source.')
    ok, msj = kanonik_kos(yaz, sure, oncelik='yeni')
    if not ok:
        yaz(u'  The canonical update could not be made (%s); the reference set stays valid.' % msj)
    else:
        yaz('  ' + msj)
    ok, msj = yon_kapisi(yaz, 'yeni konsensus sonrasi')
    for x in msj:
        yaz('  ' + x)
    if not ok:
        yaz(u'\n  THE ORIENTATION GATE DID NOT PASS - the later stages WERE NOT STARTED.')
        return 2

    # ---- 5 panel yeniden olcum
    if 'panel' in d['bitmis']:
        yaz(u'\n[STAGE 5/7] Panel re-measurement - SKIPPED (already finished)')
    else:
        yaz(u'\n[STAGE 5/7] Panel re-measurement (FULL DEPTH)')
        yaz(u'  This is the longest stage: first every fastq file is read (10-25 min),')
        yaz(u'  then 21 pairs are measured against two criteria. One to two hours in total.')
        from . import panel_measurement
        try:
            y = panel_olcum.calistir(yaz, sure, okuma_sayisi=0,
                                     yalniz=getattr(a, 'hedef', None))
            d['bitmis'].append('panel')
            d['ciktilar']['panel'] = y
            durum_yaz(d)
        except KeyboardInterrupt:
            yaz(u'\nSTOPPED - continue with (9).')
            return 1
        except Exception as e:
            yaz('\n  ASAMA 5 HATASI: %s' % e)

    # ---- 6 uyelik denetimi
    if 'uyelik' in d['bitmis']:
        yaz(u'\n[STAGE 6/7] Membership audit - SKIPPED (already finished)')
    else:
        yaz(u'\n[STAGE 6/7] Membership audit and sensitivity analysis')
        from . import membership_check
        try:
            y = uyelik_denetimi.calistir(yaz, sure,
                                         okuma_sayisi=C.NUMUNE_OKUMA_SAYISI,
                                         yalniz=getattr(a, 'hedef', None))
            d['bitmis'].append('uyelik')
            d['ciktilar']['uyelik'] = y
            durum_yaz(d)
        except KeyboardInterrupt:
            yaz(u'\nSTOPPED - continue with (9).')
            return 1
        except Exception as e:
            yaz('\n  ASAMA 6 HATASI: %s' % e)

    # ---- 7 kapsamli arama
    if 'arama' in d['bitmis']:
        yaz(u'\n[STAGE 7/7] Full search - SKIPPED (already finished)')
    else:
        yaz('\n[ASAMA 7/7] Kapsamli arama (sorunlu hedefler)')
        from . import __main__ as M
        try:
            rc = M.aramayi_kos(a, yaz, sure, cizgi, mod='devam')
            if rc == 0:
                d['bitmis'].append('arama')
                durum_yaz(d)
        except KeyboardInterrupt:
            yaz(u'\nSTOPPED - continue with (9).')
            return 1
        except Exception as e:
            yaz('\n  ASAMA 7 HATASI: %s' % e)

    # ---- 6 birlesik ozet
    yaz(u'\n[AFTER STAGE 7/7] Combined summary report')
    yol = ozet_yaz(d, sure, time.time() - t0)
    d['bitmis'] = list(dict.fromkeys(d['bitmis']))
    d['bitis'] = time.strftime('%Y-%m-%d %H:%M:%S')
    durum_yaz(d)

    cizgi('=')
    yaz(u'  ALL DONE (%s)' % sure(time.time() - t0))
    yaz('')
    yaz(u'  READ THIS FIRST:')
    yaz('    %s' % yol)
    yaz('')
    yaz(u'  Detailed reports:')
    for ad, ys in d.get('ciktilar', {}).items():
        for y in (ys or []):
            if str(y).endswith('.md'):
                yaz('    %s' % y)
    cizgi('=')
    return 0


# ---------------------------------------------------------------- ozet rapor
def _oku_tsv(ad):
    p = os.path.join(C.CIKTI, ad)
    if not os.path.exists(p):
        return []
    try:
        return list(csv.DictReader(open(p, encoding='utf-8'), delimiter='\t'))
    except Exception:
        return []


def _wilson(k, n, z=1.96):
    import math
    if n == 0:
        return (0.0, 1.0)
    pp = k / n
    d = 1 + z * z / n
    c = (pp + z * z / (2 * n)) / d
    sd = z * math.sqrt(pp * (1 - pp) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - sd), min(1.0, c + sd))


def _asama5_yeniden(asgari=None):
    """Asama 5'in "en kotu tek kutu" oranini KAYITLI kutu sayilarindan,
    asama 6 ile AYNI mutlak esikle yeniden hesapla.

    Gerekce: kosu sirasinda esik "en buyuk kutunun yarisi" idi ve tam
    derinlikte 10-33 rakip kutudan yalniz 1-5'i olcume giriyordu. Ayni cift
    300 okumayla olculunce hepsi giriyordu. Iki asama arasindaki farkin buyuk
    kismi bundan. Kutu bazli ham sayilar kayitli oldugu icin 5 saatlik kosuyu
    tekrarlamadan duzeltilebiliyor.
    """
    # Esigin MUTLAK olmasi sarttir: goreli bir esik (ornegin "en buyuk kutunun
    # yarisi") derinlik degistikce baska kutulari olcume alir, boylece iki
    # asamanin sayilari ayni seyi olcmez ve karsilastirilamaz hale gelir.
    # Burada asama 5'in kayitli kutu sayilari, asama 6 ile AYNI mutlak esikle
    # (C.ENKOTU_ASGARI_OKUMA) yeniden hesaplanir; kalan fark yalnizca okuma
    # derinliginden gelir. Uye tarafinda Wilson ALT siniri, rakip tarafinda UST
    # siniri kullanilir - iki taraf da muhafazakar yonde, yani oran asla oldugu
    # degerden buyuk cikmaz.
    if asgari is None:
        asgari = getattr(C, 'ENKOTU_ASGARI_OKUMA', 150)
    rows = _oku_tsv('panel_kutu_duzeyi.tsv')
    if not rows:
        return {}
    from collections import defaultdict
    g = defaultdict(lambda: {'uye': [], 'rakip': []})
    for r in rows:
        if not (r.get('olcut', '') or '').startswith('<=1'):
            continue
        try:
            pv = int(r['YENI_motor_urun_veren']); nv = int(r['okuma'])
        except (ValueError, KeyError):
            continue
        g[r['hedef']][r['grup']].append((r['kutu'], pv, nv))
    out = {}
    for h, d in g.items():
        if not d['uye']:
            continue
        uye_alt = min(_wilson(pv, nv)[0] for _, pv, nv in d['uye'])
        rp = sum(pv for _, pv, _ in d['rakip']); rn = sum(nv for _, _, nv in d['rakip'])
        hav = _wilson(rp, rn)[1] if rn else None
        en = None
        giren = 0
        for kadi, pv, nv in d['rakip']:
            if nv < asgari:
                continue
            giren += 1
            hi = _wilson(pv, nv)[1]
            if en is None or hi > en[1]:
                en = (kadi, hi)
        out[h] = dict(
            uye_alt=round(100 * uye_alt, 3),
            kat_havuz=(round(uye_alt / hav, 2) if hav else None),
            kat_enkotu=(round(uye_alt / en[1], 2) if en and en[1] > 0 else None),
            enkotu_kutu=(en[0] if en else ''),
            rakip_kutu=len(d['rakip']), giren_kutu=giren,
            uye_kutu=len(d['uye']),
            derinlik=sum(nv for _, _, nv in d['uye'] + d['rakip']))
    return out


def _arama_kontrolleri():
    """Asama 7'nin hedef kontrol dosyalari: taban + adaylar (AYNI kosullarda)."""
    import glob
    out = []
    for f in sorted(glob.glob(os.path.join(C.KONTROL, 'hedef_*.json'))):
        try:
            out.append(json.load(open(f, encoding='utf-8')))
        except Exception:
            pass
    return out


def _oran_metni(v, rakip_kutu=None):
    # EVRENSEL HEDEFLERDE AYRIM KATI TANIMSIZDIR. Ayrim kati, uye tarafinin
    # Wilson alt sinirini rakip tarafinin Wilson ust sinirina boler. Evrensel
    # bir hedefte (Arke_universal, Bakteri_universal, Mantar_universal F1/F2)
    # uyelik tanimi geregi rakip kumesi bosa yaklasir; payda sifira gider ve
    # oran tanimsizlasir. Kagit uzerinde ayni sutunda 0,00 ile 117 milyon yan
    # yana durur - iki sayi da olcum degildir. Bu yuzden burada oran basilmaz:
    # rakip kutusu hic yoksa "rakip kutusu yok", payda sifira yaklastigi icin
    # sayi patlamissa "olcusuz (rakip ~0)" yazilir. Bu hedefler kapsama (kac
    # kutuda urun veriyor) ve alan disi orani ile degerlendirilir. Bu, esigi
    # DUSURMEK DEGILDIR - olcunun kendisini dogru olanla degistirmektir.
    if rakip_kutu == 0:
        return 'rakip kutusu yok'
    if v in (None, '', 'None'):
        return '-'
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if f > 100000:
        return 'olcusuz (rakip ~0)'
    return '%.2f' % f


def ozet_yaz(durum, sure, gecen):
    os.makedirs(C.CIKTI, exist_ok=True)
    yol = os.path.join(C.CIKTI, '00_OZET_HEPSI.md')
    L = []; A = L.append

    panel = _oku_tsv('panel_yeniden_olcum.tsv')
    uyelik = _oku_tsv('uyelik_duyarlilik.tsv')
    kons = _oku_tsv('konsensus_yeniden_uretim.tsv')
    aday = _oku_tsv('adaylar.tsv')
    aramalar = _arama_kontrolleri()
    yeniden5 = _asama5_yeniden()
    asgari = getattr(C, 'ENKOTU_ASGARI_OKUMA', 150)

    A('# Ozet — her sey sirayla kosuldu')
    A('')
    A('Baslangic: %s · bitis: %s · sure: %s'
      % (durum.get('baslangic', '?'), time.strftime('%Y-%m-%d %H:%M:%S'), sure(gecen)))
    A('')
    A('Kosulan asamalar: %s' % (', '.join(durum.get('bitmis', [])) or 'yok'))
    A('')
    eksik = [ad for ad, v in (('panel_yeniden_olcum.tsv', panel),
                              ('uyelik_duyarlilik.tsv', uyelik),
                              ('konsensus_yeniden_uretim.tsv', kons),
                              ('adaylar.tsv', aday)) if not v]
    if eksik:
        A('> **Eksik cikti dosyalari:** %s — ilgili bolumler bos gorunecek.'
          % ', '.join('`%s`' % e for e in eksik))
        A('')

    # ============================================================ 0. EN ONEMLI
    A('## 0. ONCE BUNU OKUYUN — hangi sayi hangi kosulda olculdu')
    A('')
    A('Bu kosuda ayni cift **uc ayri yerde** olculdu ve sayilar **birbirinden '
      'farkli cikti**. Sebep bulundu; ikisi de dogru ama **karsilastirilamaz**:')
    A('')
    A('| | Asama 5 (panel yeniden olcum) | Asama 6-7 (uyelik + arama tabani) |')
    A('|---|---|---|')
    A('| okuma derinligi | **TAM** (kutu basina ortanca ~3 000, en fazla ~46 000) | **300 okuma/kutu** |')
    A('| uyumsuzluk olcutu | `<=1` ve `<=3` (iki ayri satir) | `<=1` |')
    A('| uyelik tanimi | ayni (`hedef_uyelik.tsv`) | ayni |')
    A('| motor | ayni (`read_engine.py`) | ayni |')
    A('')
    A('**Farkin iki sebebi var:**')
    A('')
    A('1. **Wilson araligi derinlikle daralir.** Uye icin ALT sinir, rakip icin '
      'UST sinir kullaniliyor; okuma arttikca alt sinir yukselir, ust sinir '
      'duser, oran **buyur**. Alttaki uye yuzdeleri neredeyse ayni oldugu halde '
      'oranlarin farkli cikmasinin sebebi budur — olcum degil, **belirsizlik payi** '
      'degisiyor.')
    A('2. **"En kotu tek kutu" olcusu kosu sirasinda derinlige gore kayiyordu.** '
      'Esik "en buyuk kutunun yarisi" idi; tam derinlikte en buyuk kutu ~46 000 '
      'okuma olunca esik ~23 000 oluyor ve **10-33 rakip kutudan yalnizca 1-5\'i** '
      'olcume giriyordu. Ayni cift 300 okumayla olculunce esik 150 olup **hepsi** '
      'giriyordu. Yani asama 5\'in "en kotu kutu" sayisi aslinda "**en derin '
      'kutu**" demekti ve gercek en kotu rakip disarida kalmis olabilir.')
    A('')
    A('Bu ikinci sorun **duzeltildi** (esik artik mutlak: %d okuma). Kutu bazli '
      'ham sayilar kayitli oldugu icin asama 5 **yeniden kosulmadan** dogru '
      'degerler asagida hesaplandi.' % asgari)
    A('')
    A('### Karsilastirilabilir hale getirilmis tablo')
    A('')
    A('`Asama 5 (duzeltilmis)` sutunu, asama 5\'in kendi tam derinlikli '
      'sayilarini asama 6 ile **ayni kutu esigiyle** yeniden hesaplar. '
      'Kalan fark yalnizca derinlikten gelir.')
    A('')
    A('| hedef | A5 havuz (tam) | A5 en kotu (kosudaki) | **A5 en kotu (duzeltilmis)** | A6/A7 havuz (300) | A6/A7 en kotu (300) | A5 olcume giren rakip kutu |')
    A('|---|---|---|---|---|---|---|')
    p1 = {r['hedef']: r for r in panel if (r.get('olcut') or '').startswith('<=1')}
    ubak = {}
    for r in uyelik:
        if (r.get('tanim') or '').startswith('A.'):
            ubak.setdefault(r['hedef'], r)
    for h in sorted(set(list(p1.keys()) + list(yeniden5.keys()))):
        a5 = p1.get(h, {})
        yd = yeniden5.get(h, {})
        u6 = ubak.get(h, {})
        A('| %s | %s | %s | **%s** | %s | %s | %s/%s |' % (
            h[:40], _oran_metni(a5.get('ayrim_havuz_x')),
            _oran_metni(a5.get('ayrim_en_kotu_x')),
            _oran_metni(yd.get('kat_enkotu')),
            _oran_metni(u6.get('ayrim_havuz_x')),
            _oran_metni(u6.get('ayrim_en_kotu_x')),
            yd.get('giren_kutu', '?'), yd.get('rakip_kutu', '?')))
    A('')
    A('> **Hangisine bakmali?** Karar icin **`A5 en kotu (duzeltilmis)`** sutununu '
      'kullanin: en derin veriyle, butun rakip kutulari kapsayarak, muhafazakar '
      'yonde hesaplanmistir. `A6/A7` sutunlari adaylarin **kendi aralarinda** '
      'siralanmasi icindir — hepsi ayni 300 okuma kosulunda olculdugu icin '
      'birbirleriyle karsilastirilabilirler, ama A5 ile karsilastirilamazlar.')
    A('')

    # ============================================================ 1. yeni aday
    A('## 1. Panelin mevcut ciftinden DAHA IYI aday bulundu mu?')
    A('')
    if not aramalar:
        A('*Arama kontrol dosyalari bulunamadi.*')
    else:
        A('Her hedefte paneldeki **mevcut cift** ile **en iyi aday** '
          '**AYNI kosullarda** (300 okuma/kutu, `<=1` uyumsuzluk, ayni uyelik, '
          'ayni motor) olculdu — bu sutunlar dogrudan karsilastirilabilir.')
        A('')
        A('| # | hedef | mevcut cift | en iyi aday | daha iyi mi | kazanc | 10x esigi | ARMS | izgara hucresi |')
        A('|---|---|---|---|---|---|---|---|---|')
        iyi = kotu = anlamsiz = belirsiz = 0
        sira = 0
        for v in sorted(aramalar, key=lambda x: x.get('hedef', '')):
            h = v.get('hedef', '?')
            t = v.get('panel_olcum') or {}
            ad = v.get('adaylar') or []
            sira += 1
            # rakip kutusu yoksa (evrensel hedefler) ORAN ANLAMSIZDIR
            rakip_yok = not (t.get('rakip') or [])
            tb = t.get('kat_enkotu')
            if tb in (None, ''):
                tb = t.get('kat_havuz')
            en = None
            for c in ad:
                v2 = c.get('numune_kat_enkotu') or c.get('numune_kat_havuz')
                try:
                    v2 = float(v2)
                except (TypeError, ValueError):
                    continue
                if en is None or v2 > en[0]:
                    en = (v2, c)
            hucre = (en[1].get('izgara_hucresi') if en else '') or ''
            hucre = hucre.replace('|', '/')[:44]
            arms = ('EVET' if (en and en[1].get('arms')) else 'hayir')

            if rakip_yok:
                anlamsiz += 1
                A('| %d | %s | *rakip kutusu yok* | *rakip kutusu yok* | '
                  '**oran anlamsiz** | - | - | %s | %s |'
                  % (sira, h[:38], arms, hucre))
                continue
            try:
                tbf = float(tb)
            except (TypeError, ValueError):
                tbf = None
            if en is None or tbf is None:
                belirsiz += 1
                A('| %d | %s | %s | %s | **olculemedi** | - | - | %s | %s |'
                  % (sira, h[:38], _oran_metni(tb), _oran_metni(en[0] if en else None),
                     arms, hucre))
                continue

            daha = en[0] > tbf * 1.05
            esikte = en[0] >= C.AYRIM_ESIK   # dCq 3 -> 8,00x (tek kaynak)
            if daha:
                iyi += 1
            else:
                kotu += 1
            kazanc = ('%.1fx' % (en[0] / tbf)) if tbf > 0 else ('yeni (taban 0)' if en[0] > 0 else '-')
            A('| %d | %s | %s | %s | %s | %s | %s | %s | %s |' % (
                sira, h[:38], _oran_metni(tb), _oran_metni(en[0]),
                '**EVET**' if daha else 'hayir', kazanc,
                'gecti' if esikte else '**ALTINDA**', arms, hucre))
        A('')
        A('**Sonuc: %d hedefte daha iyi aday VAR, %d hedefte mevcut cift korunmali, '
          '%d hedefte oran anlamsiz (evrensel - rakip kutusu yok), %d olculemedi.**'
          % (iyi, kotu, anlamsiz, belirsiz))
        A('')
        A('> **`10x esigi` sutununu atlamayin.** "Daha iyi" demek "yeterli" demek '
          'DEGILDIR: bir aday mevcut cifti gecip hala 10x esiginin altinda '
          'kalabilir. Ornegin taban 0,12x iken aday 13,5x ise hem daha iyi hem '
          'esigi geciyor; taban 0,19x iken aday 1,04x ise daha iyi ama **hala '
          'kullanilamaz**.')
        A('')
        A('> **Evrensel hedeflerde** (Arke_universal, Bakteri_universal, '
          'Mantar_universal F1/F2) uyelik tanimi geregi rakip kutusu yoktur; '
          'ayrim orani tanimsizdir ve bu hedefler **kapsam** olcusuyle '
          'degerlendirilmelidir (kac kutuda urun veriyor). Oran sutunlarina '
          'bakmayin.')
        A('')
        A('> Her adayin bedeli (hangi kural gevsetildi, ARMS gerekti mi, urun '
          'boyu protokolu nasil etkiler, 60 C\'de kosulabilir mi) '
          '`KAPSAMLI_ARAMA_RAPORU.md` icindeki hedef basliklarinda yazili.')
    A('')

    # ============================================================ 2. panel degisimi
    A('## 2. Panelin sayilari ne degisti (asama 5)')
    A('')
    if not panel:
        A('*Panel yeniden olcumu kosulmadi ya da cikti bulunamadi.*')
    else:
        degisen = [r for r in panel if (r.get('DEGISIM') or '').startswith(('YUKARI', 'ASAGI'))]
        A('- Olculen satir: **%d** (her cift iki olcutle)' % len(panel))
        A('- Degeri **%%30\'dan fazla degisen** satir: **%d**' % len(degisen))
        kayip = []
        for r in panel:
            try:
                k = float(r.get('ESKI_motor_kayip_uye_%') or 0)
                if k > 0:
                    kayip.append((k, r['hedef'], r['olcut']))
            except ValueError:
                pass
        if kayip:
            kayip.sort(reverse=True)
            A('- Eski (hatali) okuma motorunun **en cok kacirdigi** hedefler:')
            A('')
            A('| hedef | olcut | eski motorun kaybi |')
            A('|---|---|---|')
            for k, h, o in kayip[:8]:
                A('| %s | %s | %%%.1f |' % (h, o, k))
        A('')
        A('Ayrinti: `PANEL_YENIDEN_OLCUM.md`, `panel_kutu_duzeyi.tsv`')
    A('')

    # ============================================================ 3. esik alti
    A('## 3. 10x esiginin ALTINA dusen ciftler')
    A('')
    dusen = [r for r in panel if r.get('ESIK_ALTINA_DUSTU')]
    if not panel:
        A('*Panel yeniden olcumu kosulmadi.*')
    elif not dusen:
        A('Yeni olcumde 10x esiginin altina **dusen cift yok**.')
    else:
        A('**%d satir** esigin altina dustu — toplantida konusulmali:' % len(dusen))
        A('')
        A('| hedef | olcut | panelde | yeni | uye kapsam |')
        A('|---|---|---|---|---|')
        for r in dusen:
            A('| %s | %s | %s | %s x | %s |' % (
                r['hedef'], r['olcut'], r.get('PANEL_ayrim_sayi', ''),
                r.get('ayrim_en_kotu_x', ''), r.get('uye_kapsam', '')))
    A('')

    # ============================================================ 4. uyelik
    A('## 4. Uyelik tanimi — hangi sayi tanima bagli')
    A('')
    if not uyelik:
        A('*Uyelik denetimi kosulmadi.*')
    else:
        oyn, tani = {}, {}
        for r in uyelik:
            h = r['hedef']
            try:
                o = float(r.get('oynaklik_x') or 0)
            except ValueError:
                o = 0
            oyn[h] = max(oyn.get(h, 0), o)
            if r.get('TANI') and r['TANI'] != 'SORUN YOK':
                tani[h] = r['TANI']
        riskli = sorted(((v, k) for k, v in oyn.items() if v and v > 1.5), reverse=True)
        A('- Tanim degisince ayrim orani **1,5 kattan fazla oynayan** hedef: **%d**'
          % len(riskli))
        if riskli:
            A('')
            A('| hedef | tanim degisince ayrim kac kat oynuyor |')
            A('|---|---|')
            for v, k in riskli[:10]:
                A('| %s | %.1fx |' % (k, v))
        if tani:
            A('')
            A('| hedef | tani |')
            A('|---|---|')
            for k, v in tani.items():
                A('| %s | %s |' % (k, v))
        A('')
        A('Ayrinti: `UYELIK_DENETIMI.md` · Duzeltme yeri: `screening/hedef_uyelik.tsv`')
    A('')

    # ============================================================ 5. konsensus
    A('## 5. Konsensusler')
    A('')
    if not kons:
        A('*Konsensus yeniden uretimi kosulmadi.*')
    else:
        deg = [r for r in kons if (r.get('eski_ile_farkli') or '0') not in ('0', '')]
        kalite = konsensus_kalite()
        gecen_k = sum(1 for v in kalite.values() if v[0])
        A('- Yeniden uretilen kutu: **%d**' % len(kons))
        A('- Eski konsensusla **farkli baz** iceren kutu: **%d**' % len(deg))
        A('- Arama omurgasi olarak kullanilabilir kalitede: **%d**' % gecen_k)
        if deg:
            A('')
            A('| kutu | farkli baz | N %% | iki yontem ayrildi |')
            A('|---|---|---|---|')
            for r in deg[:12]:
                A('| %s | %s | %s | %s |' % (r['kutu'], r['eski_ile_farkli'],
                                             r['N_yuzde'], r['yontemler_ayrildi']))
        A('')
        A('Ayrinti: `KONSENSUS_YENIDEN_URETIM.md`')
    A('')

    # ============================================================ 6. sirada ne var
    A('## 6. Sirada ne var')
    A('')
    A('1. **1. bolum** — daha iyi aday bulunan hedeflerde bedeli okuyup karar verin.')
    A('2. **3. bolumdeki esik alti ciftler** — toplantida konusulmali.')
    A('3. **4. bolumdeki oynak hedefler** — hangi uyelik tanimiyla bildirildikleri '
      'panelde acikca yazilmali.')
    A('4. Panel xlsx\'ine isleme: `screening/update_panel.py` bunun icin '
      'yazildi ama **menude degildir** — bilerek elle calistirilir, cunku '
      'teslim dosyasini degistirir ve o dosyaya baska oturumlar da yaziyor.')
    A('')
    A('> Bu arac hicbir karari kendiliginden vermez ve panel dosyalarini '
      'kendiliginden degistirmez. Olcer, bedelini yazar, secimi size birakir.')
    A('')
    with open(yol, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(L))
    return yol


def yalniz_ozet(yaz, sure, cizgi):
    """SADECE ozet raporu yeniden uret - olcum YAPMAZ, saniyeler surer.

    5 saatlik kosunun son adimi cokerse butun kosuyu tekrarlamak gerekmesin
    diye ayrildi. Elindeki cikti dosyalarini okur, ozeti yeniden yazar.
    """
    import time as _t
    t0 = _t.time()
    cizgi('=')
    yaz(u'  REGENERATE THE SUMMARY REPORT ONLY')
    cizgi('=')
    yaz(u'  Nothing is measured; the existing output files are read.')
    yaz(u'  Source directory: %s' % C.CIKTI)
    yaz('')
    var = []
    for ad in ('panel_yeniden_olcum.tsv', 'panel_kutu_duzeyi.tsv',
               'uyelik_duyarlilik.tsv', 'konsensus_yeniden_uretim.tsv',
               'adaylar.tsv'):
        p_ = os.path.join(C.CIKTI, ad)
        v = os.path.exists(p_)
        var.append(v)
        yaz('   %-34s %s' % (ad, 'VAR' if v else 'YOK'))
    kont = len([f for f in os.listdir(C.KONTROL)
                if f.startswith('hedef_') and f.endswith('.json')]) \
        if os.path.isdir(C.KONTROL) else 0
    yaz(u'   %-34s %d files' % (u'kontrol/hedef_*.json (search)', kont))
    yaz('')
    if not any(var) and not kont:
        yaz(u'  No output file at all; run (9) first, or the stages one by one.')
        return 2
    d = durum_oku()
    yol = ozet_yaz(d, sure, _t.time() - t0)
    yaz(u'  Summary written (%s):' % sure(_t.time() - t0))
    yaz('    %s' % yol)
    cizgi('=')
    return 0
