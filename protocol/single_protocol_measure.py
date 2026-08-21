# -*- coding: utf-8 -*-
"""Single-protocol measurement, one rule, one depth, for the whole panel.

Every pair in the panel is measured under IDENTICAL settings, so the numbers
are comparable across targets. Mixed protocols were the earlier failure mode:
two pairs measured at different depths cannot be ranked against each other, and
nothing in the output revealed it.

Reports dCq (discrimination power) per pair, and the panel threshold each pair
must clear.

--- ozgun aciklama ---
TEK PROTOKOL - panelin tamamini AYNI kural ve AYNI derinlikle olcer,
tek bir siparis listesi uretir.

NEDEN VAR
---------
Paneldeki ayrim katlari bugune kadar FARKLI kosullarda uretildi: kimi satir
mm<=1, kimi mm<=3 olcutuyle; kimi 300 okuma, kimi 46 000 okuma derinliginde;
uyelik kimi satirda Kraken etiketinden, kimi satirda olculen kimlikten geldi.
Wilson araliginin genisligi derinlige bagli oldugu icin AYNI gercek ozgulluk
sig havuzda DAHA DUSUK bir "x" degeri verir. Yani o sutundaki sayilar
birbiriyle karsilastirilamaz ve 10x esigi satirdan satira ayni seyi olcmez.

Bu betik o karisikligi bitirir: TEK protokol, satir bazinda ISTISNA YOK.

Panel dosyalarina YAZMAZ. Yalniz okur, TEK_PROTOKOL_SONUC/ altina yazar.
"""

# -------------------------------------------------------------------------
# single_protocol_measure.py, paneldeki BUTUN ciftleri tek kural ve tek derinlikle
# yeniden olcer, tek bir siparis listesi uretir. Satir bazinda istisna yoktur.
#
# GİRDİ  : primer_final/ altindaki panel tablosu (screening.targets.
#          panel_oku ile), protocol/ek_ciftler.tsv (panelde olmayan ciftler),
#          uyelik_yeniden_turetme_uyelik_*.tsv (U asamasinin OLCULEN uyeligi),
#          "fastq files" altindaki ham okumalar.
# ÇIKTI  : TEK_PROTOKOL_SONUC/panel_tek_protokol.tsv (tam tablo),
#          TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv (siparis karari),
#          TEK_PROTOKOL_SONUC/kutu_bazli_ham_sayilar.tsv (k ve n; her verdikt
#          bu iki sutundan yeniden hesaplanabilir),
#          TEK_PROTOKOL_SONUC/PROTOKOL_VE_RAPOR.md, kontrol/ .
# ÇAĞRAN : verification/full_chain.py -> P tusu
#          (bat icinde: wsl -e python3 "protocol/single_protocol_measure.py" --kok .)
#
# NEDEN VAR: eski panelde satirlar farkli kosullarda olculmustu - kimi mm<=1,
# kimi mm<=3; kimi 300 okuma, kimi 46 000 okuma derinliginde. Wilson araliginin
# genisligi okuma sayisina bagli oldugu icin AYNI gercek ozgulluk sig havuzda
# DAHA DUSUK bir "x" verir. O sutundaki sayilar birbiriyle karsilastirilamazdi ve
# 10x esigi satirdan satira ayni seyi olcmuyordu. Bu betik o karisikligi bitirir.
# -------------------------------------------------------------------------
import os, sys, csv, json, time, argparse, math

VERSIYON = '1.0 (2026-08-03)'

# --- ESIK TEK KAYNAKTAN: screening/config.py -> ESIK_DCQ ---
def _esik_yukle():
    import os as _o, sys as _s
    _kok = _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))
    if _kok not in _s.path:
        _s.path.insert(0, _kok)
    from screening import config as _y
    return _y

_C = _esik_yukle()


def _sinif_yukle():
    from screening import order_classes as _s
    return _s

_S = _sinif_yukle()


# --------------------------------------------------------------- protokol
# --- ESIKLERIN KOKENI ------------------------------------------------
# 2026-08-06: esik dCq cinsinden sabitlendi. dCq >= 3 -> 2**3 = 8,00 kat.
# Bu ARTIK BIR ARAC ESIGI DEGILDIR: dCq >= 3 ozgulluk/NTC gecme olcutu olarak
# literaturde kabul gormus tabandir. Onceki 10x gercekten arac esigiydi (ilk kez
# bir kod sabiti olarak ortaya cikmisti) ve dCq 3,32'ye denk geliyordu.
# Tek kaynak: screening/config.py -> ESIK_DCQ.
# Toplantinin kendi olcutu hala FARKLI bir buyukluk ve AYRI sutunda raporlanir:
#   CALISMA_KAYDI §1.7 - "hosgoru 1-2 CAPRAZ TUR; olcu capraz tur SAYISIDIR,
#                         o turlerde olusan urun sayisi degil"
#   CALISMA_KAYDI §1.5 - "rakiplerin hicbirinde urun olusmamali" (sifir hosgoru)
# Ikisi 10x'i ne kapsar ne de onun tarafindan kapsanir. Bu yuzden ikisi de
# AYRI SUTUN olarak raporlanir ve hangisinin kim tarafindan konuldugu yazilir.
ESIK_KOKENI = _C.ESIK_KOKENI
ESIK_VERIM_NOTU = _C.ESIK_VERIM_NOTU
# MIQE/laboratuvar dili: ayrim kati -> dCq. %100 verimde her dongu 2 kat, yani
# dCq = log2(kat). 10x = 3,32 dongu. Literaturde ozgulluk/NTC gecme olcutu
# dCq >= 3 (NEB yuksek verimli qPCR veri analizi) - bizim 10x esigimiz onun
# hemen ustune dusuyor. Verim %100 varsayilir; gercek verim olculunce
# dCq = log(kat)/log(1+E) ile duzeltilmelidir.
def dcq(kat, verim=1.0):
    import math
    try:
        k = float(kat)
    except (TypeError, ValueError):
        return None
    if k <= 0:
        return None
    return round(math.log(k) / math.log(1.0 + verim), 2)

TOPLANTI_CAPRAZ_TABAN = 10.0   # bir rakip kutu "capraz" sayilmak icin en az %10 urun
TOPLANTI_CAPRAZ_HOSGORU = 2    # CALISMA_KAYDI §1.7: 1-2 capraz tur hosgoru

PROTOKOL = dict(
    olcut_asil=1,
    olcut_yan=3,
    okuma_tavani=3000,
    esik=_C.AYRIM_ESIK,        # dCq 3 -> 8,00 kat (tek kaynak)
    esik_dcq=_C.ESIK_DCQ,
    karisik='rakip',
    kapsam_esigi=0.20,
    enkotu_asgari_okuma=150,
    urun_alt=60, urun_ust=400,
)

GEREKCE = u"""
PROTOKOL VE NEDEN BOYLE SECILDI
===============================

1) DERINLIK: kutu basina EN COK %(okuma_tavani)d okuma, sabit tohum, 200-6000 bp suzgeci.
   Neden tavan var: Wilson araliginin genisligi okuma sayisina baglidir. Tavan
   konmazsa 46 000 okumalik bir kutu ile 300 okumalik bir kutu ayni tabloda yan
   yana gelir ve derin kutunun ayrimi YAPAY olarak yuksek cikar. Tavan butun
   satirlari ayni istatistiksel zemine oturtur.
   Neden 3000: kutularin buyuk cogunlugu zaten bu sayinin altinda, yani veri
   kaybi kucuk; ustelik uyeligin turetildigi olcum de bu tavanla yapildi, boylece
   uyelik ile ayrim ayni zeminde kalir.

2) ASIL OLCUT: <=1 uyumsuzluk + 3' son 2 baz TAM eslesme.
   Panelin tasarim olcutu budur; primerler bu varsayimla secildi. Karar bu
   sutuna gore verilir.

3) YAN OLCUT: <=3 uyumsuzluk + 3' son 2 baz TAM eslesme.
   Karar sutunu DEGILDIR, DAYANIKLILIK gostergesidir. Gercek PCR'de gevsek
   baglanma olur; bir cift mm<=1'de gecip mm<=3'te cokuyorsa o cift kirilgandir
   ve bu gorunmelidir (olculmus ornek: ileri primeri NL1 olan cift 8,47x -> 0,67x).

4) KARISIK KUTULAR TEK KURALA BAGLI: karisik kutu = %(karisik)s sayilir.
   Karisik kutu, hedef organizmayi KISMEN tasidigi olculen kutudur. Uye saymak
   ayrimi yapay olarak yukseltir; tamamen dislamak gercek capraz sinyali gizler.
   RAKIP saymak en kotu durumu olcer - siparis karari icin dogru taraf budur.
   Degistirilebilir (--karisik uye|rakip|disla) ama secim her ciktinin basina yazilir.

5) ESIK: dCq %(esik_dcq).1f = %(esik).2f kat, EN KOTU TEK RAKIP KUTU uzerinden
   (asgari %(enkotu_asgari_okuma)d okuma). Literatur olcutu; arac esigi DEGILDIR.
   VERIM %%100 VARSAYILDI - %%90 verimde ayni dCq 6,86 kat eder.
   Havuz kati da raporlanir ama karar en kotu kutuya gore verilir: havuz kati
   tek bir kotu kutuyu binlerce temiz okumanin icinde eritir.

6) KAPSAM ayri eksendir: bir uye kutu >=%%%(kapsam_yuzde)d urun veriyorsa "kapsandi" sayilir.
   Ayrim yuksek ama kapsam dusukse cift ozguldur fakat hedefin tamamini gormez;
   bu iki sorun birbirine karistirilmamalidir.

7) UYELIK: Kraken etiketinden DEGIL, olculen kimlikten gelir
   (uyelik_yeniden_turetme_uyelik_*.tsv). Yanlis etiketli bir kutu rakip hanesine
   yazilinca metrik hedefi hedefle kiyaslar ve mukemmel bir primer bile 1'in
   altinda cikar - olculmus ornek: 0,71x -> 8,47x, ayni primer.
""" % dict(PROTOKOL, kapsam_yuzde=int(PROTOKOL['kapsam_esigi'] * 100))


def sure_metni(sn):
    sn = int(sn)
    if sn < 90:
        return '%d saniye' % sn
    if sn < 5400:
        return '%d dakika' % round(sn / 60.0)
    return '%.1f saat' % (sn / 3600.0)


def vir(x, basamak=2):
    """Turkce ondalik: 12.96 -> '12,96'"""
    if x is None or x == '':
        return '-'
    return ('%.*f' % (basamak, x)).replace('.', ',')


# --------------------------------------------------------------- girdiler
# Proje kokunu dogrular. screening klasoru yoksa olcum modulleri ithal
# edilemez; erken ve acik hata vermek, yarim kosudan iyidir.
def kok_bul(arg):
    kok = os.path.abspath(arg or '.')
    if not os.path.isdir(os.path.join(kok, 'screening')):
        sys.exit('HATA: %s icinde screening klasoru yok. --kok ile proje '
                 'klasorunu verin.' % kok)
    return kok


def uyelik_dosyasi(kok):
    """En yeni uyelik dosyasini bulur - ADA gore degil, ZAMANA gore.

    2026-08-10 duzeltmesi. Eski kod iki globu birlestirip a[-1] aliyordu; bu
    "en yeni" DEMEK DEGILDI. Siralama alfabetikti ve engine_SONUC
    girdileri her zaman kokteki girdilerden SONRA geliyordu. Yani alt
    klasorde 1 Agustos tarihli bir dosya olsa, kokteki 3 Agustos tarihliyi
    yenerdi. Su an tek aday var, o yuzden davranis degismiyor; ama bir
    sonraki kosu ikinci bir dosya uretirse sessizce yanlis uyelik secilirdi.
    Uyelik yanlis olursa ayrim kati oldugundan kucuk ya da buyuk cikar
    (olculmus ornek: ayni cift 0,71x ile 8,47x arasinda oynadi).
    """
    import glob
    a = glob.glob(os.path.join(kok, 'uyelik_yeniden_turetme_uyelik_*.tsv'))
    a += glob.glob(os.path.join(kok, 'engine_SONUC', '*uyelik*.tsv'))
    if not a:
        return None
    # once zaman, esitlikte ad - iki olcut de ACIK yazili
    a.sort(key=lambda p: (os.path.getmtime(p), os.path.basename(p)))
    return a[-1]


def uyelik_oku(yol):
    """hedef -> dict(uye=[...], karisik=[...], rakip=[...], sinif=...)"""
    out = {}
    with open(yol, encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            bol = lambda s: [x for x in (r.get(s) or '').split(';') if x.strip()]
            out[r['hedef'].strip()] = dict(
                sinif=(r.get('sinif') or '').strip(),
                uye=bol('yeni_uye_kutular') or bol('eski_uye_kutular'),
                karisik=bol('karisik_kutular'),
                rakip=bol('rakip_kutular'))
    return out


def ek_ciftler_oku(kok):
    """protocol/ek_ciftler.tsv - panel TSV'sinde OLMAYAN ciftler.
    Kullanici elle duzenleyebilir. Sutunlar: hedef, sinif, F, R, urun_bp, not"""
    yol = os.path.join(kok, 'protocol', 'ek_ciftler.tsv')
    if not os.path.exists(yol):
        return []
    out = []
    with open(yol, encoding='utf-8') as fh:
        for r in csv.DictReader((s for s in fh if not s.startswith('#')), delimiter='\t'):
            if not (r.get('hedef') or '').strip():
                continue
            out.append(dict(hedef=r['hedef'].strip(), sinif=(r.get('sinif') or '').strip(),
                            F=r['F'].strip().upper(), R=r['R'].strip().upper(),
                            urun_bp=int(float(r.get('urun_bp') or 0)),
                            kaynak='EK', ta=(r.get('ta') or '').strip(),
                            uyelik_hedefi=(r.get('uyelik_hedefi') or '').strip(),
                            duzey=(r.get('duzey') or '').strip(),
                            not_=(r.get('not') or '').strip()))
    return out


def kutu_adi_normalize(kutu):
    """A1_1_2223 ve A1-1_2223 ayni kutudur (A1-1 orneginin dosyalari alt cizgili).
    Sinif-ornek ayraci daima TIRE olacak sekilde duzeltir."""
    if '_' not in kutu:
        return kutu
    bas, _, son = kutu.rpartition('_')          # son = taxid
    return bas.replace('_', '-') + '_' + son


# --------------------------------------------------------------- olcum
# ---------------------------------------------------------------------------
# ASIL OLCUM. Sira: ciftleri topla -> uyeligi coz -> okuma havuzlarini kur ->
# her cifti HEM mm<=1 (asil, karar) HEM mm<=3 (yan, dayaniklilik) ile olc.
#
# UYELIK KRAKEN ETIKETINDEN GELMEZ, U asamasinin OLCULEN kimliginden gelir.
# Yanlis etiketli bir kutu rakip hanesine yazilinca metrik hedefi hedefle
# kiyaslar ve mukemmel bir primer bile 1'in altinda cikar (olculmus ornek:
# 0,71x -> 8,47x, ayni primer, yalniz uyelik duzeltildi).
#
# DERINLIK TAVANI SART: tavan konmazsa 46 000 okumalik kutu ile 300 okumalik
# kutu ayni tabloda yan yana gelir ve derin kutunun ayrimi YAPAY olarak yuksek
# cikar. Tavan butun satirlari ayni istatistiksel zemine oturtur.
#
# KONTROL NOKTASI MUHRU (_ayar): okuma tavani, karisik kutu kurali, betik surumu
# ve UYELIK DOSYASININ ADI muhre dahildir. Uyelik tablosu tazelenince eski
# olcumler sessizce yeniden kullanilmaz - bu betigin var olus sebebi tam olarak
# "farkli kosulda uretilmis sayilarin yan yana durmasi" idi.
# ---------------------------------------------------------------------------
def calistir(kok, okuma_tavani, karisik_kural, yalniz=None, sifirla=False):
    sys.path.insert(0, kok)
    from screening import config as C, motor, numune as N, hedefler as H

    CIKTI = os.path.join(kok, 'TEK_PROTOKOL_SONUC')
    KONTROL = os.path.join(CIKTI, 'kontrol')
    os.makedirs(KONTROL, exist_ok=True)
    if sifirla:
        for f in os.listdir(KONTROL):
            try:
                os.remove(os.path.join(KONTROL, f))
            except OSError as e:
                print('  silinemedi: %s (%s)' % (f, e))

    gunluk = open(os.path.join(CIKTI, 'kosu_gunlugu.txt'), 'a', encoding='utf-8')

    def yaz(s=''):
        print(s, flush=True)
        gunluk.write(s + '\n'); gunluk.flush()

    yaz('=' * 78)
    yaz('  TEK PROTOKOL - panelin tamami ayni kuralla yeniden olculuyor')
    yaz('  surum %s   %s' % (VERSIYON, time.strftime('%Y-%m-%d %H:%M')))
    yaz('=' * 78)

    # --- ciftler -------------------------------------------------------
    panel, panel_yolu = H.panel_oku()
    ciftler = []
    for d in panel:
        ciftler.append(dict(hedef=d['hedef'], sinif=d['sinif'], F=d['F'], R=d['R'],
                            urun_bp=d['urun_bp'], ta=d.get('ta', ''),
                            duzey=d.get('duzey', ''), kaynak='PANEL',
                            panel_ayrim=d.get('ayrim', ''), not_=''))
    # EK ciftler panel TSV'sinde OLMAYANLAR icindir. Bir cift panele
    # eklendiginde ek_ciftler.tsv'de de kalirsa AYNI HEDEF IKI KEZ olculur ve
    # iki farkli sayi uretir. 2026-08-11'de tam bu oldu: Petriella_cinsi hem
    # panelden (kendi uyeligiyle 0,88x) hem EK'ten (uyelik_hedefi
    # Petriella_musispora oldugu icin 11,03x) geldi ve siparis listesinde iki
    # satir olustu - biri "siparis edilebilir", oteki "esik alti".
    panel_adlari = {c['hedef'].strip() for c in ciftler}
    ek = ek_ciftler_oku(kok)
    for e in ek:
        if e['hedef'].strip() in panel_adlari:
            print('  EK atlandi (panelde zaten var): %s' % e['hedef'])
            continue
        ciftler.append(dict(e, panel_ayrim=''))
    if yalniz:
        ciftler = [c for c in ciftler
                   if any(y.strip().lower() in c['hedef'].lower()
                          for y in yalniz.split(','))]
    if not ciftler:
        sys.exit('HATA: olculecek cift bulunamadi.')

    # --- uyelik --------------------------------------------------------
    uy_yol = uyelik_dosyasi(kok)
    if not uy_yol:
        sys.exit('HATA: uyelik_yeniden_turetme_uyelik_*.tsv bulunamadi.\n'
                 '      Once verification/full_chain.py -> secenek U kosulmalidir.')
    uyelik = uyelik_oku(uy_yol)
    yaz('  uyelik kaynagi : %s' % os.path.basename(uy_yol))
    yaz('  cift sayisi    : %d  (panel %d + ek %d)' % (len(ciftler), len(panel), len(ek)))
    yaz('  asil olcut     : <=%d uyumsuzluk + 3\' son 2 baz TAM' % PROTOKOL['olcut_asil'])
    yaz('  yan olcut      : <=%d uyumsuzluk + 3\' son 2 baz TAM' % PROTOKOL['olcut_yan'])
    yaz('  derinlik       : kutu basina en cok %d okuma (SATIR BAZINDA ISTISNA YOK)' % okuma_tavani)
    yaz('  karisik kutu   : %s sayilir' % karisik_kural.upper())
    yaz('  esik           : %s, en kotu tek rakip kutu uzerinden'
        % _C.esik_metni())
    yaz('  esik kokeni    : %s' % _C.ESIK_KOKENI)
    yaz('  UYARI          : %s' % _C.ESIK_VERIM_NOTU)
    yaz('')

    # --- kutular -------------------------------------------------------
    kut = {k['kutu']: k for k in H.kutular()}
    eksik_uyari = set()

    def coz(adlar):
        out = []
        for a in adlar:
            a2 = kutu_adi_normalize(a.strip())
            if a2 in kut:
                out.append(kut[a2])
            elif a.strip() in kut:
                out.append(kut[a.strip()])
            else:
                eksik_uyari.add(a.strip())
        return out

    baglam = {}
    for c in ciftler:
        u = uyelik.get(c['hedef'])
        if u is None:
            # ek cift, uyelik satiri hedef adiyla eslesmiyorsa sinif genelini kullan
            u = uyelik.get(c.get('uyelik_hedefi', ''), None)
        if u is None:
            baglam[c['hedef']] = None
            continue
        uye = coz(u['uye'])
        kar = coz(u['karisik'])
        rak = coz(u['rakip'])
        if karisik_kural == 'uye':
            uye, kar_eklenen = uye + kar, []
        elif karisik_kural == 'rakip':
            rak = rak + kar
        # 'disla' -> hicbir tarafa eklenmez
        if not rak:   # uyelik satirinda rakip bos ise sinifin geri kalani rakiptir
            uye_ad = {k['kutu'] for k in uye} | {k['kutu'] for k in kar}
            rak = [k for k in kut.values()
                   if k['sinif'] == (u['sinif'] or c['sinif']) and k['kutu'] not in uye_ad]
        baglam[c['hedef']] = dict(uye=uye, rakip=rak, karisik=kar)

    gerekli = {}
    for b in baglam.values():
        if b:
            for k in b['uye'] + b['rakip']:
                gerekli[k['kutu']] = k
    if eksik_uyari:
        yaz('  UYARI: uyelik tablosundaki %d kutu adi fastq klasorunde bulunamadi: %s'
            % (len(eksik_uyari), ', '.join(sorted(eksik_uyari)[:6])))
        yaz('         (A1-1 orneginin dosyalari alt cizgili adlandirilmis; betik')
        yaz('          bunu kendi icinde duzeltir, yine de bulunamayanlar yukarida.)')
    yaz('  okunacak kutu  : %d' % len(gerekli))
    yaz('')
    # --- havuz kurulumu -------------------------------------------------
    def ilerK(i, n, ad):
        print('   ... okuma havuzu %d/%d  %s          ' % (i, n, ad), end='\r', flush=True)

    t0 = time.time()
    yaz('Okuma havuzlari kuruluyor (%d kutu). Bu adimda ekranda yalniz kutu adlari' % len(gerekli))
    yaz('akar, takilmis DEGILDIR; asil olcum bundan sonra baslar ve her cift ayri kaydedilir.')
    nm = N.Numune(list(gerekli.values()), n=okuma_tavani, ilerle=ilerK, otorite=True)
    top_okuma = sum(h.n_okuma for h in nm.havuz.values())
    yaz('\nHavuzlar hazir: %d kutu, %d okuma  (%s)' % (len(gerekli), top_okuma, sure_metni(time.time() - t0)))
    tahmin = len(ciftler) * 2 * max(1.0, top_okuma / 20000.0)
    yaz('TAHMINI OLCUM SURESI: ~%s   (kesintiye dayaniklidir, ayni secenekle devam eder)'
        % sure_metni(tahmin))
    yaz('')

    # --- olcum ----------------------------------------------------------
    def kp_yolu(ad):
        t = ''.join(ch if ch.isalnum() else '_' for ch in ad)
        return os.path.join(KONTROL, 'cift_%s.json' % t)

    # O-10: uyelik kaynagi muhre DAHIL. Uyelik tablosu tazelenince eski
    # kontrol noktalari sessizce yeniden kullanilmamali.
    #
    # 2026-08-11 DUZELTME (uyelik ICERIK muhru). Muhurde uyelik dosyasinin
    # yalniz ADI vardi. Dosya YERINDE duzeltilince ad degismiyor, muhur tutuyor
    # ve olcum "onceki kosudan alindi" diye ESKI uyelikle geri geliyor. Bugun
    # tam bu oldu: Mantar_universal (F2) uyeliginden dort protist kutu cikarildi,
    # olcum yeniden kosuldu ve iki hedef de onbellekten dondu - degisiklik
    # sayilara hic yansimadi. Dizi muhrunde 10 Agustos'ta duzeltilen hatanin
    # aynisi, bu sefer uyelik tarafinda. Artik dosyanin ICERIGININ md5'i de
    # muhurde: satiri degisen uyelik tablosu kontrol noktasini gecersiz kilar.
    import hashlib as _hl0
    with open(uy_yol, 'rb') as _fh0:          # okunamazsa PATLASIN: sessiz
        _uy_muhru = _hl0.md5(_fh0.read()).hexdigest()[:12]   # "okunamadi" muhru
    # (once try/except vardi ve io modulu import edilmedigi icin muhur her
    #  kosuda "okunamadi" cikiyordu - yani sabit. Sabit muhur, muhur degildir:
    #  uyelik tablosu degisse de tutardi. Hatanin yutulmasi kontrolun kendisini
    #  gorunmez bicimde ise yaramaz hale getiriyordu.)
    AYAR = dict(okuma=okuma_tavani, karisik=karisik_kural, surum=VERSIYON,
                uyelik=os.path.basename(uy_yol), uyelik_icerik=_uy_muhru)

    # 2026-08-10 DUZELTME (dizi muhru). Muhurde primer DIZISI YOKTU. Bunun
    # sonucu: bir ciftin ileri/geri dizisi degistirildiginde kontrol noktasi
    # gecerli sayiliyor ve ESKI olcum "onceki kosudan alindi" diye geri
    # veriliyordu. 10 Agustos'ta iki cift degistirildi ve iki ayri tam kosu
    # (5 sa 29 dk + 2 sa 0 dk) eski dizileri olcup yeni sandi. Artik her
    # ciftin muhrune kendi F+R dizisinin md5'i giriyor; dizi degisirse
    # kontrol noktasi otomatik gecersiz olur.
    import hashlib as _hl

    def _ayar_of(c):
        d = dict(AYAR)
        d['dizi'] = _hl.md5(
            ((c.get('F') or '') + '|' + (c.get('R') or '')).encode('utf-8')
        ).hexdigest()[:12]
        return d
    sonuc = []
    tb = time.time()
    for i, c in enumerate(ciftler, 1):
        kp = kp_yolu(c['hedef'] + '|' + c['kaynak'])
        if os.path.exists(kp):
            try:
                v = json.load(open(kp, encoding='utf-8'))
                # O-10 sonrasi geriye donuk uyum: eski kontrol noktalarinda
                # 'uyelik' anahtari yok. ORTAK anahtarlar tutuyorsa kabul edilir
                # ama UYARI basilir - sessizce yeniden kullanilmaz.
                _e = v.get('_ayar') or {}
                _bek = _ayar_of(c)
                # DIZI muhru esitse degil, VARSA ve FARKLIYSA kontrol noktasi
                # kesin gecersizdir. Eski kontrol noktalarinda 'dizi' anahtari
                # hic yoktur; o durumda da yeniden olculur, sessizce kabul YOK.
                if _e.get('dizi') != _bek['dizi']:
                    yaz(u'  %s: kontrol noktasi DIZI muhru tutmuyor '
                        u'(kayitli %s, simdi %s) - yeniden olculuyor.'
                        % (c['hedef'][:40], _e.get('dizi') or 'yok', _bek['dizi']))
                    raise ValueError('dizi muhru tutmadi')
                # UYELIK ICERIK muhru de dizi muhru gibidir: yoksa ya da
                # tutmuyorsa kontrol noktasi GECERSIZ. "Ortak anahtarlar tuttu"
                # diye kabul edilemez - eksik olan anahtar, tam da uyelik
                # tablosundaki degisikligi goren anahtardir. (2026-08-11: bu
                # geriye donuk uyum yolu yuzunden protist duzeltmesi olcume
                # yansimadi, 22 hedefin 22'si onbellekten dondu.)
                if _e.get('uyelik_icerik') != _bek['uyelik_icerik']:
                    yaz(u'  %s: kontrol noktasi UYELIK muhru tutmuyor '
                        u'(kayitli %s, simdi %s) - yeniden olculuyor.'
                        % (c['hedef'][:40], _e.get('uyelik_icerik') or 'yok',
                           _bek['uyelik_icerik']))
                    raise ValueError('uyelik muhru tutmadi')
                _ortak = {k: _e.get(k) for k in _bek if k in _e}
                _uyar = (_e != _bek and _ortak == {k: _bek[k] for k in _ortak})
                if _e == _bek or _uyar:
                    if _uyar:
                        yaz(u'  UYARI: %s kontrol noktasi ESKI ayar muhruyle yazilmis '
                            u'(eksik: %s). Ortak anahtarlar tuttugu icin kabul edildi; '
                            u'uyelik tablosu degistiyse --sifirla ile yeniden kosun.'
                            % (c['hedef'][:40], ', '.join(sorted(set(AYAR) - set(_e))) or '-'))
                    sonuc.append(v)
                    yaz('[%2d/%2d] %-46s  (onceki kosudan alindi)' % (i, len(ciftler), c['hedef'][:46]))
                    continue
            except Exception:
                pass
        b = baglam.get(c['hedef'])
        r = dict(c); r['_ayar'] = _ayar_of(c); r['olcum'] = {}
        if not b or not b['uye']:
            r['hata'] = 'uyelik yok ya da uye kutu bulunamadi'
            yaz('[%2d/%2d] %-46s  ATLANDI (%s)' % (i, len(ciftler), c['hedef'][:46], r['hata']))
        else:
            for mm in (PROTOKOL['olcut_asil'], PROTOKOL['olcut_yan']):
                o = nm.olc(c['F'], c['R'], b['uye'], b['rakip'],
                           lo=PROTOKOL['urun_alt'], hi=PROTOKOL['urun_ust'], mm=mm)
                r['olcum'][str(mm)] = o
            r['uye_n'] = len(b['uye']); r['rakip_n'] = len(b['rakip']); r['karisik_n'] = len(b['karisik'])
            o1 = r['olcum'][str(PROTOKOL['olcut_asil'])]
            d1, g1, day1 = karar(o1, c['hedef'], c.get('duzey', ''))
            yaz('[%2d/%2d] %-46s  %s x  %-11s | kapsam %s'
                % (i, len(ciftler), c['hedef'][:46], vir(g1), d1,
                   o1.get('uye_kapsam_pay') if o1 else '-'))
        json.dump(r, open(kp, 'w', encoding='utf-8'), ensure_ascii=False, indent=1, default=str)
        sonuc.append(r)
        gecen = time.time() - tb
        print('        gecen %s | tahmini kalan %s'
              % (sure_metni(gecen), sure_metni(gecen / i * (len(ciftler) - i))), flush=True)

    yaz('')
    yaz('Olcum bitti (%s). Ciktilar yaziliyor...' % sure_metni(time.time() - tb))
    raporla(CIKTI, sonuc, dict(uyelik=os.path.basename(uy_yol), okuma=okuma_tavani,
                               karisik=karisik_kural, panel=os.path.basename(panel_yolu)), yaz)
    rc = cikti_denetle(yaz, 'P (TEK PROTOKOL)', [
        (os.path.join(CIKTI, 'panel_tek_protokol.tsv'), 'panel_tek_protokol.tsv'),
        (os.path.join(CIKTI, 'SIPARIS_LISTESI.tsv'), 'SIPARIS_LISTESI.tsv')])
    gunluk.close()
    return rc


# --------------------------------------------------------------- ciktilar
def _o(r, mm):
    return (r.get('olcum') or {}).get(str(mm)) or {}


# ---------------------------------------------------------------------------
# EVRENSEL HEDEFLERDE AYRIM KATI TANIMSIZDIR.
# Ayrim kati = (uye alt siniri) / (rakip ust siniri). Bakteri_universal butun
# bakterileri, Arke_universal butun arkeleri cogaltmak icin tasarlandi; bu
# satirlarda "rakip" diye bir kume YOKTUR. Rakip kumesi bosa yaklastikca payda
# sifira gider ve oran ya 0/0 olur ya da devasa bir sayi - nitekim eski panelde
# ayni sutunda 0,00 ile 117 milyon yan yana duruyordu. O sayilar bir seyi
# olcmuyor.
#
# Bu yuzden bu satirlar burada sayisal verdikt ALMAZ, OLCULEMEDI isaretlenir.
# Dogru olcu KAPSAMA + ALAN DISI oranidir ve K asamasinda uygulanir. Bu, 10x
# esigini DUSURMEK DEGILDIR: oranin paydasi tanimsiz oldugu icin o esik bu
# satirlarda zaten uygulanamaz. Diger butun satirlarda 10x aynen durur.
# ---------------------------------------------------------------------------
def evrensel_mi(hedef, duzey=''):
    """O-6: evrensel/alan hedeflerinde ayrim katinin PAYDASI tanimsizdir
    (rakip kumesi yoktur). Bu satirlar sayisal verdikt almamalidir; dogru olcu
    KAPSAMA + ALAN DISI'dir ve K asamasinda uygulanir."""
    ad = (hedef or '').lower()
    return ('universal' in ad or 'evrensel' in ad
            or (duzey or '').strip().lower() == 'alan')


# ---------------------------------------------------------------------------
# Karar EN KOTU TEK RAKIP KUTU uzerinden verilir, havuz uzerinden degil: havuz
# kati tek bir kotu kutuyu binlerce temiz okumanin icinde eritir ve gercekte
# capraz veren bir cift temiz gorunur. Havuz kati yine de raporlanir, ama karar
# sutunu degildir.
#
# UC AYRI DURUM VARDIR ve birbirine karistirilmaz:
#   ESIK USTU / ESIK ALTI - olculdu, karar verildi.
#   OLCULEMEDI            - yeterli derinlikte rakip kutu yok ya da hedef
#                           evrensel. Bu "esik alti" DEGILDIR; karar yokluguyla
#                           basarisizligi ayni haneye yazmak yanlis olurdu.
# En kotu kutu olcusu uretilemediginde havuza dusulur ama bu ACIKCA dayanak
# sutununa yazilir.
# ---------------------------------------------------------------------------
def karar(o, hedef='', duzey=''):
    """(durum, deger, dayanak). Karar EN KOTU TEK RAKIP KUTU uzerinden verilir;
    o olcu uretilemiyorsa (yeterli derinlikte rakip kutu yok) havuza duser ve
    bu ACIKCA isaretlenir. Hicbiri yoksa OLCULEMEDI - 'esik alti' DEGILDIR."""
    if not o:
        return ('OLCULEMEDI', None, 'olcum yok')
    if evrensel_mi(hedef, duzey):
        return ('OLCULEMEDI', None,
                'EVRENSEL HEDEF - ayrim katinin paydasi tanimsiz. Dogru olcu '
                'KAPSAMA + ALAN DISI; K asamasinda uygulanir.')
    g = o.get('kat_enkotu')
    if g is not None:
        return ('ESIK USTU' if g >= PROTOKOL['esik'] else 'ESIK ALTI', g, 'en kotu tek kutu')
    h = o.get('kat_havuz')
    if h is not None:
        return ('ESIK USTU' if h >= PROTOKOL['esik'] else 'ESIK ALTI', h,
                'HAVUZ (yeterli derinlikte rakip kutu yok - en kotu kutu olcusu uretilemedi)')
    return ('OLCULEMEDI', None, 'rakip kutu yok')


# ---------------------------------------------------------------------------
# Uc cikti uretir: tam tablo, ham sayilar ve siparis listesi.
#
# HAM SAYILAR (k = urun veren okuma, n = kutudaki okuma) ayri bir dosyaya
# yazilir cunku bugune kadar hicbir okuyucu bir verdikti YENIDEN HESAPLAYAMIYORDU.
# Butun verdiktler bu iki sayidan turer; yayimlanmasi karar kurallarini
# denetlenebilir kilar.
#
# IKI OLCUT AYRI SUTUNDA DURUR ve birbirinin yerine GECMEZ:
#   ayrim_mm1_ARAC_OLCUTU        - 10x, bu aracin olcutu (toplanti karari DEGIL).
#   TOPLANTI_OLCUTU_capraz_kutu  - %10 ustu urun veren rakip KUTU SAYISI
#                                  (CALISMA_KAYDI 1.7, hosgoru 1-2 capraz tur).
# dCq sutunu ayni sayinin laboratuvar dilindeki karsiligidir (dCq = log2(kat),
# %100 verim varsayimiyla); yeni bir olcut degil, ayni olcunun cevirisidir.
#
# DAMGALAR bir cifti reddetmez, KOSULLU yapar: olcute duyarli (mm<=3'te cokuyor),
# sig karar kutusu, tek/iki uye kutu, kismi olcum, yalniz en kotu kutu gecti.
# Hicbir damgasi olmayan ve iki olcutu birden gecen cift KOSULSUZ isaretlenir.
# ---------------------------------------------------------------------------
def raporla(CIKTI, sonuc, meta, yaz):
    E = PROTOKOL['esik']; A = PROTOKOL['olcut_asil']; Y = PROTOKOL['olcut_yan']
    basli = (u'# Bu dosya TEK PROTOKOLLE uretildi - butun satirlar ayni kural ve ayni derinlik.\n'
             u'# uyelik kaynagi : %(uyelik)s   (Kraken etiketi KULLANILMADI)\n'
             u'# derinlik       : kutu basina en cok %(okuma)d okuma, satir bazinda istisna YOK\n'
             u'# asil olcut     : <=1 uyumsuzluk + 3\' son 2 baz TAM  (karar bu sutuna gore)\n'
             u'# yan olcut      : <=3 uyumsuzluk + 3\' son 2 baz TAM  (dayaniklilik gostergesi)\n'
             u'# karisik kutu   : %(karisik)s sayildi\n'
             u'# esik           : %(esik)s, EN KOTU TEK RAKIP KUTU uzerinden\n'
             u'# esik kokeni    : %(koken)s\n'
             u'# VERIM UYARISI  : %(verim)s\n') % dict(
                 meta, esik=_C.esik_metni(E), koken=_C.ESIK_KOKENI,
                 verim=_C.ESIK_VERIM_NOTU)

    # ---------- 1) tam tablo ----------
    yol = os.path.join(CIKTI, 'panel_tek_protokol.tsv')
    with open(yol, 'w', encoding='utf-8', newline='') as fh:
        fh.write(basli)
        w = csv.writer(fh, delimiter='\t')
        _kokp = os.path.dirname(os.path.abspath(CIKTI))
        _kimp = _S.kimlik_tablosu(_kokp)
        _kynp = _S.kaynak_tablosu(_kokp)
        w.writerow(['hedef'] + _S.kimlik_sutun_basliklari() + _S.kaynak_sutun_basliklari()
                   + ['kaynak', 'sinif', 'urun_bp', 'F', 'R',
                    'ASIL_ayrim_mm1', 'ASIL_ayrim_havuz_mm1', 'ASIL_kapsam_mm1',
                    'YAN_ayrim_mm3', 'YAN_ayrim_havuz_mm3', 'YAN_kapsam_mm3',
                    'uye_kutu', 'karisik_kutu', 'rakip_kutu',
                    'uye_alt_%', 'en_kotu_rakip_kutu', 'esik_gecti_mi',
                    'karar_dayanagi', 'olcute_duyarli_mi', 'panelin_eski_degeri',
                    'rakip_olculen', 'rakip_toplam', 'kismi_olcum_mu', 'not'])
        for r in sonuc:
            o1, o3 = _o(r, A), _o(r, Y)
            d1, g1, day1 = karar(o1, r['hedef'], r.get('duzey', ''))
            d3, g3, _ = karar(o3, r['hedef'], r.get('duzey', ''))
            # O-5: 150 okuma tabani rakip kutularin bir kismini sessizce eliyordu.
            ro, rt = o1.get('rakip_olculen'), o1.get('rakip_toplam')
            kismi = 'hayir'
            if ro is not None and rt:
                if ro == 0:
                    kismi = 'EVET - hicbir rakip kutu 150 okumayi gecmedi'
                elif ro < rt / 2.0:
                    kismi = 'EVET - rakiplerin %d/%d si olcume girdi' % (ro, rt)
            if kismi.startswith('EVET') and d1 == 'ESIK USTU':
                d1 = 'ESIK USTU (KISMI OLCUM)'
            gecti = d1
            duyarli = 'EVET' if (d1 == 'ESIK USTU' and d3 == 'ESIK ALTI') else 'hayir'
            w.writerow([r['hedef']] + _S.kimlik_sutunlari(_kimp, r['hedef'])
                       + _S.kaynak_sutunlari(_kynp, r['hedef'])
                       + [r['kaynak'], r.get('sinif', ''), r.get('urun_bp', ''),
                        r.get('F', ''), r.get('R', ''),
                        vir(g1), vir(o1.get('kat_havuz')), o1.get('uye_kapsam_pay', ''),
                        vir(g3), vir(o3.get('kat_havuz')), o3.get('uye_kapsam_pay', ''),
                        r.get('uye_n', ''), r.get('karisik_n', ''), r.get('rakip_n', ''),
                        vir(o1.get('uye_alt'), 3), o1.get('enkotu_kutu', ''),
                        gecti, day1, duyarli, r.get('panel_ayrim', ''),
                        o1.get('rakip_olculen', ''), o1.get('rakip_toplam', ''),
                        kismi, r.get('hata', '')])
    yaz('  yazildi: %s' % yol)

    # ---------- 1b) KUTU BAZLI HAM SAYILAR (madde 7e) ----------
    # Hicbir okuyucu bugune kadar bir verdikti YENIDEN HESAPLAYAMIYORDU.
    # k = urun veren okuma, n = kutudaki toplam okuma. Butun verdiktler bu
    # iki sayidan turer; yayimlanmasi butun karar kurallarini denetlenebilir kilar.
    yolk = os.path.join(CIKTI, 'kutu_bazli_ham_sayilar.tsv')
    with open(yolk, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# HAM SAYILAR - her verdikt bu iki sutundan turer.\n')
        fh.write(u'# k = urun veren okuma, n = kutudaki okuma. oran = k/n.\n')
        fh.write(u'# Wilson: uye tarafi ALT sinir, rakip tarafi UST sinir (z=1,96).\n')
        fh.write(basli)
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'olcut_mm', 'kutu', 'grup', 'k', 'n', 'oran_%'])
        for r in sonuc:
            for mm in (A, Y):
                o = _o(r, mm)
                for grup in ('uye', 'rakip'):
                    for satir in (o.get(grup) or []):
                        try:
                            ad, k, n, yuzde = satir[0], satir[1], satir[2], satir[3]
                        except (IndexError, TypeError):
                            continue
                        w.writerow([r['hedef'], mm, ad, grup, k, n, vir(yuzde)])
    yaz('  yazildi: %s' % yolk)

    # ---------- 2) TEK siparis listesi ----------
    gecen, kalan, olculemeyen = [], [], []
    for r in sonuc:
        d = karar(_o(r, A), r.get('hedef', ''), r.get('duzey', ''))[0]
        (gecen if d == 'ESIK USTU' else kalan if d == 'ESIK ALTI' else olculemeyen).append(r)
    yol2 = os.path.join(CIKTI, 'SIPARIS_LISTESI.tsv')
    with open(yol2, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# TEK SIPARIS LISTESI - tek protokolle uretildi.\n')
        fh.write(basli)
        # --- 2026-08-06: esik alti satirlar LISTEDEN CIKARILMAZ, siniflandirilir.
        _snf = {}
        for _r in sonuc:
            _d, _g, _ = karar(_o(_r, A), _r['hedef'], _r.get('duzey', ''))
            _snf[_r['hedef']] = _S.sinifla(_g, olculemedi=(_d == 'OLCULEMEDI'))
        _kesin = [x for x in _snf.values() if x[0] in (u'KESIN', u'EVRENSEL')]
        _kos = [x for x in _snf.values() if x[0] == u'KOSULLU']
        _oner = [x for x in _snf.values() if x[0] == u'ONERILMEZ']
        fh.write(u'#\n')
        fh.write(u'# ================  UC SAYI, AYRI AYRI  ================\n')
        fh.write(u'#   KESIN     : %d cift = %d oligo   (dCq >= %.1f, ya da evrensel/kapsam)\n'
                 % (len(_kesin), 2 * len(_kesin), _C.ESIK_DCQ))
        fh.write(u'#   KOSULLU   : %d cift = %d oligo   (dCq %.1f-%.1f - siparis edilebilir AMA dogrulama sart)\n'
                 % (len(_kos), 2 * len(_kos), _S.KOSULLU_ALT_DCQ, _C.ESIK_DCQ))
        fh.write(u'#   ONERILMEZ : %d cift             (dCq < %.1f - listede kalir, gerekcesi yazili)\n'
                 % (len(_oner), _S.KOSULLU_ALT_DCQ))
        fh.write(u'# Esigi gecemeyen satir SESSIZCE SILINMEZ; karar sizindir.\n')
        fh.write(u'# (arac sayimi: esik ustu %d, esik alti %d, olculemeyen %d)\n'
                 % (len(gecen), len(kalan), len(olculemeyen)))
        fh.write(u'# "OLCUTE DUYARLI" isaretli satirlar mm<=1 de gecip mm<=3 te cokuyor - kirilgandir.\n')
        fh.write(u'#\n')
        fh.write(u'# ESIK: %s\n' % _C.esik_metni())
        fh.write(u'# ESIGIN KOKENI: %s\n' % ESIK_KOKENI)
        fh.write(u'# VERIM UYARISI: %s\n' % ESIK_VERIM_NOTU)
        fh.write(u'#   Toplantinin KENDI olcutu farkli bir buyukluktur: capraz TUR SAYISI\n')
        fh.write(u'#   (CALISMA_KAYDI §1.7, hosgoru 1-2). O yuzden ayri sutunda:\n')
        fh.write(u'#   TOPLANTI_OLCUTU_capraz_kutu = >=%%%d urun veren rakip kutu sayisi.\n'
                 % int(TOPLANTI_CAPRAZ_TABAN))
        fh.write(u'#   Iki olcut birbirinin yerine GECMEZ.\n')
        fh.write(u'# dCq_karsiligi: laboratuvarin konustugu birim. dCq = log2(kat),\n')
        fh.write(u'#   %%100 verim varsayimiyla. 10x = 3,32 dongu. Literaturde ozgulluk\n')
        fh.write(u'#   gecme olcutu dCq >= 3 (NEB). Gercek verim olculunce duzeltilmeli.\n')
        w = csv.writer(fh, delimiter='\t')
        _kok0 = os.path.dirname(os.path.abspath(CIKTI))
        _kim = _S.kimlik_tablosu(_kok0)
        _kyn = _S.kaynak_tablosu(_kok0)
        w.writerow(['sira', 'SINIF', 'durum', 'siparis_sarti', 'hedef']
                   + _S.kimlik_sutun_basliklari() + _S.kaynak_sutun_basliklari()
                   + ['oligo_adi_F', 'F', 'oligo_adi_R', 'R', 'urun_bp',
                      'ayrim_mm1', 'dCq_karsiligi', 'esikten_uzaklik_dCq',
                      'LABORATUVARDA_NE_YAPILMALI',
                      'ayrim_mm3', 'havuz_mm1',
                      'TOPLANTI_OLCUTU_capraz_kutu', 'kapsam_mm1',
                      'karar_veren_kutu', 'karar_kutusu_k', 'karar_kutusu_n',
                      'uye_kutu_sayisi', 'damgalar'])
        n = 0
        for etiket, kume in (('ESIK USTU - SIPARIS EDILEBILIR', gecen),
                             ('ESIK ALTI - SIPARIS EDILMEZ', kalan),
                             ('OLCULEMEDI - KARAR YOK', olculemeyen)):
            for r in kume:
                n += 1
                o1, o3 = _o(r, A), _o(r, Y)
                d1, g1, day1 = karar(o1, r['hedef'], r.get('duzey', ''))
                d3, g3, _ = karar(o3, r['hedef'], r.get('duzey', ''))
                kod = ''.join(ch if ch.isalnum() else '_' for ch in r['hedef'])[:24]
                # --- madde 2: TOPLANTI OLCUTU = capraz kutu sayisi (verim orani DEGIL)
                capraz = sum(1 for x in (o1.get('rakip') or [])
                             if len(x) > 3 and (x[3] or 0) >= TOPLANTI_CAPRAZ_TABAN)
                # --- madde 3: karar veren kutunun HAM sayilari
                kk, kn = '', ''
                ek = o1.get('enkotu_kutu')
                for x in (o1.get('rakip') or []):
                    if x and x[0] == ek:
                        kk, kn = x[1], x[2]; break
                # --- madde 4 + 3 + 7a: damgalar ve siparis sarti
                damga = []
                if d1 == 'ESIK USTU' and d3 == 'ESIK ALTI':
                    damga.append(u'OLCUTE DUYARLI (mm<=3 te cokuyor)')
                if kn and int(kn) < 300:
                    damga.append(u'SIG KARAR KUTUSU (n=%s)' % kn)
                if (r.get('uye_n') or 0) and int(r['uye_n']) <= 2:
                    damga.append(u'TEK/IKI UYE KUTU - hedef ici degiskenlik SINANMADI')
                if kismi.startswith('EVET'):
                    damga.append(kismi)
                hv = o1.get('kat_havuz')
                if d1 == 'ESIK USTU' and (hv is None or hv < E):
                    damga.append(u'YALNIZ en kotu kutu gecti, HAVUZ gecmedi')
                # kosulsuz siparis: iki olcut + iki taban birden
                if (d1 == 'ESIK USTU' and d3 == 'ESIK USTU'
                        and hv is not None and hv >= E and not damga):
                    sart = 'KOSULSUZ'
                elif d1 == 'ESIK USTU':
                    sart = 'KOSULLU'
                else:
                    sart = '-'
                _sn, _dq, _uz, _lab = _snf[r['hedef']]
                w.writerow([n, _sn, etiket, sart, r['hedef']]
                           + _S.kimlik_sutunlari(_kim, r['hedef'])
                           + _S.kaynak_sutunlari(_kyn, r['hedef'])
                           + [kod + '_F', r.get('F', ''),
                              kod + '_R', r.get('R', ''), r.get('urun_bp', ''),
                              vir(g1), vir(_dq), vir(_uz), _lab,
                              vir(g3), vir(hv), capraz,
                              o1.get('uye_kapsam_pay', ''), ek, kk, kn,
                              r.get('uye_n', ''), '; '.join(damga) or '-'])
    yaz('  yazildi: %s' % yol2)

    # ---------- 3) protokol + rapor ----------
    yol3 = os.path.join(CIKTI, 'PROTOKOL_VE_RAPOR.md')
    with open(yol3, 'w', encoding='utf-8') as fh:
        fh.write(u'# Tek protokolle panel olcumu\n\n')
        fh.write(u'Uretim: %s · betik surumu %s\n\n' % (time.strftime('%Y-%m-%d %H:%M'), VERSIYON))
        fh.write(u'## Sonuc\n\n')
        fh.write(u'- Esigi (%.0fx) GECEN cift: **%d** → **%d oligo**\n' % (E, len(gecen), 2 * len(gecen)))
        fh.write(u'- Esigin ALTINDA kalan cift: **%d**\n' % len(kalan))
        fh.write(u'- OLCULEMEYEN (karar yok, esik alti DEGIL): **%d**\n' % len(olculemeyen))
        duy = [r for r in gecen if karar(_o(r, Y))[0] == 'ESIK ALTI']
        fh.write(u'- Esigi gecen ama OLCUTE DUYARLI (mm<=3 te cokuyor): **%d**\n\n' % len(duy))
        fh.write(u'```\n' + GEREKCE + u'\n```\n\n')
        fh.write(u'## Tablo\n\n')
        fh.write(u'| hedef | kaynak | mm<=1 (asil) | mm<=3 (yan) | kapsam | durum |\n|---|---|---|---|---|---|\n')
        for r in sorted(sonuc, key=lambda x: -(karar(_o(x, A))[1] if karar(_o(x, A))[1] is not None else -1)):
            o1, o3 = _o(r, A), _o(r, Y)
            d1, g1, day1 = karar(o1, r['hedef'], r.get('duzey', ''))
            d3, g3, _ = karar(o3, r['hedef'], r.get('duzey', ''))
            fh.write(u'| %s | %s | %s | %s | %s | %s |\n' % (
                r['hedef'], r['kaynak'], vir(g1), vir(g3),
                o1.get('uye_kapsam_pay', '-'),
                d1 + ('' if day1 == 'en kotu tek kutu' else ' (%s)' % day1)))
        fh.write(u'\n## Okuma sirasi\n\n1. Once bu dosya. 2. `SIPARIS_LISTESI.tsv`. '
                 u'3. Ayrinti icin `panel_tek_protokol.tsv`.\n')
    yaz('  yazildi: %s' % yol3)
    yaz('')
    yaz('  ESIGI GECEN: %d cift (%d oligo)   ESIK ALTI: %d   OLCULEMEYEN: %d'
        % (len(gecen), 2 * len(gecen), len(kalan), len(olculemeyen)))



# --------------------------------------------------------------- guvenlik agi
def cikti_denetle(yaz, ad, dosyalar, asgari=1):
    """Asama bittiginde KENDI ciktisini denetler.

    Beklenen satir sayisi sifirsa ya da dosya hic yoksa SESSIZCE DEVAM ETMEZ:
    acik Turkce hata basar ve sifirdan farkli kod dondurur. Gece boyunca bos
    sonuc uretip sabah "hicbir sey bulunamadi" dememesi icin.
    """
    sorun = []
    for yol, etiket in dosyalar:
        if not os.path.exists(yol):
            sorun.append(u'%s URETILMEDI (%s)' % (etiket, yol)); continue
        try:
            with open(yol, encoding='utf-8') as fh:
                n = sum(1 for s in fh if s.strip() and not s.startswith('#'))
            n = max(0, n - 1)          # baslik satiri
        except OSError as e:
            sorun.append(u'%s OKUNAMADI (%s)' % (etiket, e)); continue
        if n < asgari:
            sorun.append(u'%s BOS - %d veri satiri (en az %d bekleniyordu)'
                         % (etiket, n, asgari))
    if not sorun:
        return 0
    yaz('')
    yaz('  ' + '!' * 70)
    yaz(u'  %s ASAMASI BOS CIKTI URETTI - ZINCIR BURADA DURDURULDU' % ad)
    for x in sorun:
        yaz(u'    - %s' % x)
    yaz('')
    yaz(u'  NEDEN DURDURULDU: sonraki asama bu dosyayi girdi olarak okuyacakti.')
    yaz(u'  Bos girdiyle devam etmek cokme uretmez, ANLAMSIZ AMA INANDIRICI bir')
    yaz(u'  ozet uretir - tam da avladigimiz sessiz hata deseni budur.')
    yaz(u'  Yukaridaki kosu gunlugunu okuyup sebebi giderin, sonra ayni secenegi')
    yaz(u'  tekrar secin; bitmis isler kontrol noktalarindan atlanacaktir.')
    yaz('  ' + '!' * 70)
    return 4


def girdi_denetle(yaz, ad, dosyalar):
    """Asama BASLAMADAN once ihtiyac duydugu dosyalar var mi ve dolu mu."""
    eksik = []
    for yol, etiket, uretici in dosyalar:
        if not os.path.exists(yol):
            eksik.append(u'%s yok (%s) -> once %s asamasini kosun' % (etiket, yol, uretici))
            continue
        with open(yol, encoding='utf-8') as fh:
            n = sum(1 for s in fh if s.strip() and not s.startswith('#'))
        if n <= 1:
            eksik.append(u'%s BOS (%s) -> %s asamasi sonuc uretmemis'
                         % (etiket, yol, uretici))
    if not eksik:
        return 0
    yaz('')
    yaz('  ' + '!' * 70)
    yaz(u'  %s ASAMASI BASLATILMADI - GIRDI EKSIK' % ad)
    for x in eksik:
        yaz(u'    - %s' % x)
    yaz('  ' + '!' * 70)
    return 5

# Komut satiri: --okuma derinlik tavani, --karisik karisik kutulara ne yapilacagi
# (uye|rakip|disla; varsayilan rakip = en kotu durumu olcer), --yalniz alt kume,
# --sifirla kontrol noktalarini siler.

# --- CLI value normalisation ------------------------------------------------
# English option values are accepted alongside the original Turkish ones and
# mapped back here. The internal values are unchanged on purpose: they are
# compared in dozens of places and, in some cases, written to output files.
# Translating the interface must not translate the data.
_DEGER = {'auto': 'oto', 'manual': 'elle', 'none': 'yok', 'quick': 'hizli',
          'full': 'tam', 'member': 'uye', 'competitor': 'rakip',
          'exclude': 'disla'}


def _ing_deger(a):
    for _ad in ('nt', 'literatur', 'ncbi', 'karisik', 'moduller', 'mod'):
        _v = getattr(a, _ad, None)
        if isinstance(_v, str) and _v in _DEGER:
            setattr(a, _ad, _DEGER[_v])
    return a

def main():
    p = argparse.ArgumentParser(description='Tek protokolle panel olcumu')
    p.add_argument('--root', '--kok', dest='kok', default='.')
    p.add_argument('--reads', '--okuma', dest='okuma', type=int, default=PROTOKOL['okuma_tavani'],
                   help='cap on reads per bin (0 = all of them)')
    p.add_argument('--mixed', '--karisik', dest='karisik', choices=['member', 'competitor', 'exclude', 'uye', 'rakip', 'disla'], default=PROTOKOL['karisik'])
    p.add_argument('--only', '--yalniz', dest='yalniz', default=None, help='only targets whose name contains this (for testing)')
    p.add_argument('--reset', '--sifirla', dest='sifirla', action='store_true')
    a = p.parse_args()
    a = _ing_deger(a)
    kok = kok_bul(a.kok)
    return calistir(kok, a.okuma, a.karisik, a.yalniz, a.sifirla)


if __name__ == '__main__':
    sys.exit(main() or 0)
