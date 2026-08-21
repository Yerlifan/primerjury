# -*- coding: utf-8 -*-
"""verification TURU - TEK PROTOKOL kosusunda esik altinda kalan satirlari
dort ayri yolla kurtarmayi dener.

TEMEL CERCEVE
-------------
Toplantida secilen hedef aslinda bir KUTUdur; isim o gunku Kraken etiketidir ve
YANLIS OLABILIR. Amac o kutulari cogaltmaktir. Dolayisiyla bir hedefin uye
kumesi = degerlendiricinin isaret ettigi kutu + olcumun AYNI ORGANIZMA oldugunu gosterdigi
butun kutular. Petriella'yi kurtaran sey buydu: organizma dokuz kutuya dagilmisti
ve digerleri "rakip" hanesinde durdugu icin metrik hedefi hedefle kiyasliyordu.

DORT YOL (hepsi bu tek secenekte, sirayla)
  1) Evrensel hedeflerde OLCUYU duzelt   (ayrim kati bu satirlarda tanimsiz)
  2) UYELIK DARALTMA                     (kapsami tam olmayan her hedef)
  3) YENIDEN TASARIM + ARMS              (kil payi kalanlar)
  4) ESLENIGI KALMIS satirlari temizle    (yerine daha iyisi gelmis olanlar)

ESIK INDIRILMEZ. Bir satiri esigi dusurerek gecirmek YASAKTIR. Olcuyu duzeltmek
(yol 1) ile esigi gevsetmek AYRI SEYLERDIR; bu betik yalnizca birincisini yapar
ve hangi olcunun neden kullanildigini her satirda yazar.

Panel dosyalarina YAZMAZ. Yalniz okur, KURTARMA_SONUC/ altina yazar.
"""

# -------------------------------------------------------------------------
# kurtarma_turu.py — TEK PROTOKOL (P) kosusunun esik altinda biraktigi satirlari
# bes ayri yolla kurtarmayi dener. ESIK ASLA INDIRILMEZ.
#
# GİRDİ  : TEK_PROTOKOL_SONUC/panel_tek_protokol.tsv (esik alti satirlar),
#          uyelik_yeniden_turetme_uyelik_*.tsv (U asamasinin olculen uyeligi),
#          protocol/ek_ciftler.tsv (hedef adi takma adlari),
#          konsensus_kanonik/ ve "fastq files"/ (olcum kaynaklari).
# ÇIKTI  : KURTARMA_SONUC/kurtarma_satirlari.tsv (hedef basina TEK satir),
#          KURTARMA_SONUC/yeni_adaylar.tsv (yol 3 ve yol 5 taramasinin adaylari),
#          KURTARMA_SONUC/KURTARMA_RAPORU.md, kontrol/ (hedef basina JSON).
#          Panel dosyalarina YAZMAZ.
# ÇAĞRAN : screening.bat -> K tusu
#          (bat icinde: wsl -e python3 "verification/kurtarma_turu.py" --kok .)
#
# BES YOL, sirayla denenir:
#   yol 1  evrensel hedeflerde OLCUYU duzelt (ayrim kati orada tanimsiz)
#   yol 2  UYELIK DARALTMA (olculen dizi kimligine gore, KOSULSUZ)
#   yol 3  YENIDEN TASARIM + ARMS (kil payi kalan satirlar)
#   yol 4  ESLENIGI KALMIS satirlari dusenlere tasi
#   yol 5  COK LOKUSLU arama (konsensusun tamamini bolgelere ayirip her birinde
#          ayri tasarim; Petriella'yi kurtaran sey buydu - cozum ITS'te degil
#          LSU'daydi)
#
# OLCUYU DUZELTMEK ile ESIGI GEVSETMEK AYRI SEYLERDIR. Bu betik yalnizca
# birincisini yapar ve hangi satirda hangi olcunun kullanildigini "olcu" sutununa
# acikca yazar.
# -------------------------------------------------------------------------
import os, sys, csv, json, time, argparse, math

VERSIYON = '1.0 (2026-08-03)'

# ESIK TEK KAYNAKTAN GELIR: screening/yapilandirma.py -> ESIK_DCQ = 3.0
# Kat karsiligi 2 ** ESIK_DCQ = 8,00. Sabit sayi GOMULMEZ; dCq degisirse
# tek yerden degisir. Gerekce ve verim uyarisi o dosyada yazili.
# ESIGI KOD ICINDEN INDIRMEK YINE YASAKTIR - bkz. modul basligi. Degisiklik
# yalniz ESIK_DCQ uzerinden ve BILEREK yapilir; 2026-08-06'da dCq 3'e sabitlendi.
def _esik_yukle():
    """Esigi TEK KAYNAKTAN okur: screening/yapilandirma.py.
    verification/ ile screening/ kardes klasorler oldugu icin kok buradan
    turetilir; betik hangi calisma dizininden cagrilirsa cagrilsin bulur."""
    import os as _o, sys as _s
    _kok = _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))
    if _kok not in _s.path:
        _s.path.insert(0, _kok)
    from screening import yapilandirma as _y
    return _y

_C = _esik_yukle()
ESIK = _C.AYRIM_ESIK
KAPSAM_ESIGI = 0.20         # bir uye kutu 'kapsandi' sayilmak icin en az %20 urun
OKUMA_TAVANI = 3000         # TEK PROTOKOL ile ayni
ENKOTU_ASGARI = 150
URUN_ALT, URUN_UST = 60, 400
KIMLIK_ESIGI = 99.0         # ayni organizma sayilmak icin konsensus kimligi
KIL_PAYI_ALT = 5.0          # yol 3'e girecek satirin asgari mevcut kati

# --- evrensel hedefler icin AYRI olcu ------------------------------------
EVRENSEL_KAPSAMA_ESIGI = 0.80   # uye kutularin en az %80'i urun vermeli
EVRENSEL_ALANDISI_UST = 5.0     # alan disi kutularda Wilson ust sinir en cok %5

GEREKCE_EVRENSEL = u"""
YOL 1 - EVRENSEL HEDEFLERDE OLCU NEDEN DEGISTI (esik DUSURULMEDI)

Ayrim kati = (uye alt siniri) / (rakip ust siniri). Evrensel bir hedefte
"rakip" diye bir kume YOKTUR: Bakteri_universal butun bakterileri, Arke_universal
butun arkeleri cogaltmak icin tasarlandi. Payda sifira giderse oran ya 0/0
(tanimsiz, betik 0,00 yazar) ya da devasa bir sayi olur - nitekim ayni sutunda
0,00 ile 117 056 685 yan yana duruyor. Bu sayilar bir SEYI OLCMUYOR.

Bu satirlarda iddia da farklidir: "her seyi ayirt ederim" degil, "grubun
tamamini gorurum, grup disina tasmam". Dogru olcu bu iddiayi olcendir:

  KAPSAMA      = uye kutularin kaci >=%%%d urun veriyor
  ALAN DISI    = hedef grubun DISINDAKI kutularda urun veren okuma orani
                 (Wilson UST siniri - muhafazakar taraf)

GECME OLCUTU (ikisi birden):
  KAPSAMA   >= %%%d  ve  ALAN DISI <= %%%.0f

Bu, 10x esiginin gevsetilmesi DEGILDIR: 10x oranini bu satirlarda uygulamak
zaten mumkun degil, cunku oranin paydasi tanimsiz. Diger butun satirlarda 10x
esigi AYNEN durur.
""" % (int(KAPSAM_ESIGI * 100), int(EVRENSEL_KAPSAMA_ESIGI * 100), EVRENSEL_ALANDISI_UST)

# --- daha once OLCULMUS, tekrar denenmeyecek satirlar --------------------
BILINEN = {
    'Proteiniphilum_cinsi': dict(
        sonuc='KURTARILAMAZ',
        sebep=u'Hedef organizma numunede beyan edildigi gibi MEVCUT DEGIL. '
              u'Uye kutularin 2/3\'u olculen kimlikte Fermentimonas caenicola '
              u'(%95,33 ve %97,13, cins FARKLI) ve cift Fermentimonas\'i bilerek '
              u'disliyor (0/137). Uye kutularda urun HIC yok (0/3 kapsam). '
              u'Uyeligi daraltmak burada bir sey kurtarmaz - sorun uyelik degil '
              u'HEDEF TANIMI. Yeniden olculmedi, zaman harcanmadi.',
        yol=u'atlandi (onceden olculdu)'),
    # 2026-08-06: BU KAYIT KALDIRILDI (yorumda birakildi).
    # Elle yazilmis 'atlandi' damgasi yuzunden satir kurtarma merdivenine
    # HIC girmiyordu, dolayisiyla yol 5 (cok lokuslu arama) da kosmuyordu.
    # Gerekcesi 16S'e dayaniyordu ve dogruydu, ama 'baska lokusta da olmaz'
    # SONUCUNU icermiyordu - o olculmemisti. Olculdu: A2 kutulari (4309 bp
    # tam operon) ayni organizma (A1 ile %98,62-99,38 kimlik), operonun
    # UC bolgesi de tarandi, hicbirinde aday cikmadi. Artik satir merdivene
    # girer ve sonucu olcumden gelir.
    #     'Methanosarcina mazei / M. soligelidi grubu': dict(
    #         sonuc='KURTARILAMAZ',
    #         sebep=u'Sinirlayan rakip A1-4_3078083 (M. hadiensis) AYRI BIR ORGANIZMA: '
    #               u'mazei kutularina konsensus kimligi %98,61-98,75, mazei kutulari '
    #               u'kendi aralarinda %99,79-99,93; 16S tur esigi ~%98,7. Okuma duzeyi '
    #               u'prob testi de ayni sonucu verdi (hadiensis okumalari kendi probunda '
    #               u'kaliyor, 171\'e 5). Yani uyelik daraltmasi bu satiri kurtarmaz. '
    #               u'Kurtarma yolu VAR ama primer degisikligi gerektirir: mazei ile '
    #               u'hadiensis 16S\'te ~19 pozisyonda ayriliyor, mevcut cift bunlarin '
    #               u'hicbirini tutmuyor (yol 3\'e girer, ayri tasarim isi).',
    # 
    #     yol=u'atlandi (onceden olculdu)'),

}


def sure_metni(sn):
    sn = int(sn)
    if sn < 90:
        return '%d saniye' % sn
    if sn < 5400:
        return '%d dakika' % round(sn / 60.0)
    return '%.1f saat' % (sn / 3600.0)


def vir(x, b=2):
    if x is None or x == '':
        return '-'
    try:
        return ('%.*f' % (b, float(x))).replace('.', ',')
    except (TypeError, ValueError):
        return str(x)


def kutu_adi_normalize(k):
    if '_' not in k:
        return k
    bas, _, son = k.rpartition('_')
    return bas.replace('_', '-') + '_' + son


# -------------------------------------------------------------------------
# WILSON SKOR ARALIGI - NEDEN HAM ORAN DEGIL
# Ham oran k/n kucuk orneklemde yaniltici olur: 3 okumanin 3'u urun verdiyse ham
# oran %100'dur ama arkasinda kanit yoktur; 200 okumanin 0'i verdiyse ham oran
# %0'dir ama gercek oran %1,5 olabilir. Wilson araligi bu belirsizligi sayiya
# doker ve daima MUHAFAZAKAR taraf secilir: uye tarafinda ALT sinir, rakip
# tarafinda UST sinir. Yol 1'deki "alan disi" orani da bu yuzden UST sinirdir -
# hedef grubun disinda ne kadar tasma OLABILECEGINI olcer.
# -------------------------------------------------------------------------
def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 1.0)
    p = k / float(n); d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    s = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - s), min(1.0, c + s))


# ---------------------------------------------------------------- girdiler
# P asamasinin panel tablosu bu turun TEK giris kapisidir. Dosya yoksa kosu
# baslamaz: yoklugunda "kurtarilacak satir yok" diye sessizce bitmek, gercekte
# hicbir sey denenmedigi halde is bitmis izlenimi verirdi.
def tek_protokol_oku(kok):
    """TEK_PROTOKOL_SONUC/panel_tek_protokol.tsv -> [{hedef, kaynak, karar, ...}]"""
    yol = os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'panel_tek_protokol.tsv')
    if not os.path.exists(yol):
        sys.exit('HATA: %s yok.\n      Once screening.bat -> secenek (P) kosulmalidir.' % yol)
    with open(yol, encoding='utf-8') as fh:
        satirlar = [s for s in fh if not s.startswith('#')]
    return list(csv.DictReader(satirlar, delimiter='\t')), yol


def _f(s):
    """'8,45' -> 8.45 ; '-' -> None"""
    s = (s or '').strip().replace(',', '.')
    try:
        return float(s)
    except ValueError:
        return None


def uyelik_dosyasi(kok):
    import glob
    # 2026-08-10: "a[-1]" EN YENI DEMEK DEGILDI. Iki glob alfabetik siralanip
    # birlestiriliyordu, yani engine_SONUC girdileri tarihine
    # bakilmaksizin kokteki girdileri yeniyordu. tek_protokol_olc.py ayni
    # tuzagi tasiyordu; ikisi ayri dosya secseydi K ile P farkli uyelikle
    # olcer ve dCq'lari karsilastirilamaz olurdu. Ikisi de artik ZAMANA gore
    # seciyor ve ayni dosyayi buluyor.
    a = glob.glob(os.path.join(kok, 'uyelik_yeniden_turetme_uyelik_*.tsv'))
    a += glob.glob(os.path.join(kok, 'engine_SONUC', '*uyelik*.tsv'))
    if not a:
        return None
    a.sort(key=lambda p: (os.path.getmtime(p), os.path.basename(p)))
    return a[-1]


# Uyelik tablosunu okur. yeni_uye_kutular yoksa eski_uye_kutular'a duser -
# tablo hangi surumden gelirse gelsin satir bos kalmaz.
def uyelik_oku(yol):
    out = {}
    with open(yol, encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            bol = lambda s: [kutu_adi_normalize(x.strip())
                             for x in (r.get(s) or '').split(';') if x.strip()]
            out[r['hedef'].strip()] = dict(
                sinif=(r.get('sinif') or '').strip(),
                uye=bol('yeni_uye_kutular') or bol('eski_uye_kutular'),
                karisik=bol('karisik_kutular'), rakip=bol('rakip_kutular'))
    return out


# ek_ciftler.tsv'deki cift adi ile uyelik tablosundaki hedef adi ayni olmayabilir;
# bu tablo ikisini eslestirir. Eslestirilmezse ek ciftler "uyelik yok" diye
# atlanirdi.
def takma_adlar(kok):
    """protocol/ek_ciftler.tsv: hedef -> uyelik_hedefi (uyelik tablosundaki ad)"""
    yol = os.path.join(kok, 'protocol', 'ek_ciftler.tsv')
    out = {}
    if not os.path.exists(yol):
        return out
    with open(yol, encoding='utf-8') as fh:
        for r in csv.DictReader((s for s in fh if not s.startswith('#')), delimiter='\t'):
            h = (r.get('hedef') or '').strip(); u = (r.get('uyelik_hedefi') or '').strip()
            if h and u:
                out[h] = u
    return out


# Bir hedefin evrensel (alan duzeyi) olup olmadigi. Bu satirlarda ayrim kati
# tanimsizdir; yol 1 devreye girer.
def evrensel_mi(hedef, duzey=''):
    ad = hedef.lower()
    return ('universal' in ad or 'evrensel' in ad
            or (duzey or '').strip().lower() == 'alan')


# --- TOPLANTIDA ISTENIP PANELE HIC GIRMEMIS TALEPLER -------------------
# Zincir panelden basliyor; panelde satiri olmayan bir talep zincire HIC
# gorunmez ve "hic denenmeden yapilamadi" diye kalir. Bu talepler icin en
# azindan BIR tasarim denemesi yapilir: omurga olarak kutunun KENDI konsensusu
# kullanilir (M. barkeri hedefinde uygulanan yontemin aynisi).
PANELSIZ_TALEPLER = [
    dict(hedef='Podospora_pseudopauciseta (PANELSIZ TALEP)', karar='Karar 1',
         sinif='F1', uye=['F1-1_2093780', 'F1-4_2093780'],
         not_=u'Toplantida TUR ozgul istendi. Organizmanin KENDISI numunede yok '
               u'(bes referans ciftinden ucu F1 sinifinin 85 804 okumasinin '
               u'tamamina karsi tarandi, sifir urun; bolluk ust siniri %0,011). '
               u'AMA KUTU VAR ve olculen kimligi Petriella (F1-4_2093780, %99,58). '
               u'Bu deneme adlandirilmis turu degil KUTUYU hedefler.'),
    dict(hedef='Dictyostelium_discoideum_44689 (PANELSIZ TALEP)', karar='Karar 1',
         sinif='F1', uye=['F1-1_44689', 'F1-2_44689', 'F1-3_44689', 'F1-4_44689'],
         not_=u'Toplantida TUR ozgul istendi. Kraken etiketi olcumle curutuldu ve '
               u'kutuya YENI bir ad konulamadi; hedef dizisi tanimlanamadigi icin '
               u'panele hic girmedi. Bu deneme kutunun KENDI konsensusunu omurga '
               u'alir - ad bilinmese de kutuyu cogaltan bir cift bulunabilir.'),
]


# ---------------------------------------------------------------- yollar
# -------------------------------------------------------------------------
# YOL 1 - EVRENSEL HEDEFLERDE OLCUNUN DEGISTIRILMESI (esik DUSURULMEZ)
#
# Ayrim kati = (uye alt siniri) / (rakip ust siniri). Evrensel bir hedefte
# "rakip" diye bir kume YOKTUR: Bakteri_universal butun bakterileri,
# Arke_universal butun arkeleri cogaltmak icin tasarlandi. Rakip kumesi bosa
# yaklastikca payda sifira gider ve oran ya 0/0 (tanimsiz) olur ya da devasa bir
# sayi - nitekim eski panelde ayni sutunda 0,00 ile 117 milyon yan yana duruyordu.
# O sayilar bir SEYI OLCMUYOR.
#
# Bu satirlarda iddia da farklidir: "her seyi ayirt ederim" degil, "grubun
# tamamini gorurum, grup disina tasmam". Dogru olcu bu iddiayi olcendir:
#   KAPSAMA   = uye kutularin kaci >=%20 urun veriyor
#   ALAN DISI = hedef grubun DISINDAKI kutularda urun veren okuma orani,
#               Wilson UST siniri (muhafazakar taraf)
# Gecme olcutu IKISI BIRDEN: kapsama >= %80 ve alan disi <= %5.
#
# BU ESIGI DUSURMEK DEGILDIR. 10x oranini bu satirlarda uygulamak zaten mumkun
# degil, cunku paydasi tanimsiz. Diger butun satirlarda 10x AYNEN durur.
#
# Rakip okuma hic yoksa alan disi OLCULEMEZ ve None dondurulur. Eskiden 0,0
# yaziliyordu, yani mumkun olan EN OLUMLU deger uretilip esik kendiliginden
# geciliyordu - olcum yoklugu basari sayiliyordu.
# -------------------------------------------------------------------------
def yol1_evrensel(nm, uye, rakip, F, R):
    """Kapsama + alan disi orani. Ayrim kati KULLANILMAZ (paydasi tanimsiz)."""
    ka = na = 0
    for k in uye:
        h = nm.havuz.get(k['kutu'])
        if h is None:
            continue
        p, n, _ = h.urun_veren(F, R, URUN_ALT, URUN_UST, 1)
        na += 1
        if n and p / float(n) >= KAPSAM_ESIGI:
            ka += 1
    rp = rn = 0
    for k in rakip:
        h = nm.havuz.get(k['kutu'])
        if h is None:
            continue
        p, n, _ = h.urun_veren(F, R, URUN_ALT, URUN_UST, 1)
        rp += p; rn += n
    kapsama = (ka / float(na)) if na else 0.0
    # O-7: rakip okuma yoksa 'alan disi' OLCULEMEZ. Eskiden 0.0 yaziliyordu,
    # yani mumkun olan EN OLUMLU deger uretilip esik kendiliginden geciliyordu.
    alandisi = (100.0 * wilson(rp, rn)[1]) if rn else None
    gecti = (kapsama >= EVRENSEL_KAPSAMA_ESIGI
             and alandisi is not None and alandisi <= EVRENSEL_ALANDISI_UST)
    return dict(kapsama=kapsama, kapsam_pay='%d/%d' % (ka, na), alandisi=alandisi,
                alandisi_pay='%d/%d' % (rp, rn), gecti=gecti)


_KOD = None


def _enc(s):
    """engine/uyelik_yeniden_turet.py ile AYNI kodlama."""
    global _KOD
    import numpy as np
    if _KOD is None:
        m = np.full(256, 4, dtype=np.uint8)
        for i, c in enumerate('ACGT'):
            m[ord(c)] = i
        _KOD = m
    return _KOD[np.frombuffer(s.encode(), dtype=np.uint8)]


# Iki konsensusun yuzde kimligi. Infix (HW) hizalama: kisa dizi uzunun ICINE
# hizalanir, uclardaki fazlalik cezalandirilmaz. Konsensus uzunluklari cok farkli
# oldugu icin (1,5 kb ile 4,5 kb yan yana) global hizalama ayni organizmayi
# farkli gosterirdi.
def hw_kimlik(a, b):
    """Kisa olani sorgu, uzun olanin ICINE (infix/HW) hizala; yuzde kimlik dondur.
    engine/uyelik_yeniden_turet.py'deki hw_kimlik ile ayni tanim; o dosya
    bir betik (paket degil) oldugu icin burada yeniden yazildi, ithal edilmedi."""
    import numpy as np
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    if not q or not t:
        return 0.0
    Q, T = _enc(q), _enc(t)
    onceki = np.zeros(len(T) + 1, dtype=np.int32)   # HW: bastan bosluk bedava
    for i in range(len(Q)):
        simdi = np.empty_like(onceki)
        simdi[0] = i + 1
        esit = (T != Q[i]) | (T == 4) | (Q[i] == 4)
        # Sol komsu bagimliligi (ekleme) VEKTORLESTIRILDI:
        #   simdi[j] = min(aday[j], simdi[j-1]+1)
        # a[j] = simdi[j]-j konursa a[j] = min(aday[j]-j, a[j-1]) olur, yani
        # kosan minimum. np.minimum.accumulate ile tek gecis. Python ic dongusu
        # kaldirildi: 1,5 kb x 1,5 kb hizalama dakikalardan saniyelere iner.
        aday = np.minimum(onceki[:-1] + esit, onceki[1:] + 1)
        aday = np.concatenate(([i + 1], aday))
        idx = np.arange(len(aday))
        simdi = np.minimum.accumulate(aday - idx) + idx
        onceki = simdi
    d = int(onceki.min())
    return round(100.0 * (1 - d / float(max(len(q), 1))), 2)


# -------------------------------------------------------------------------
# YOL 2 - UYELIK DARALTMA
#
# Toplantida secilen hedef aslinda bir KUTUdur; isim o gunku Kraken etiketidir ve
# yanlis olabilir. Bir hedefin uye kumesi = isaret edilen kutu + OLCUMUN ayni
# organizma oldugunu gosterdigi butun kutular. Uye kutular konsensus dizi
# kimligine gore (esik %99) tek baglantili kumelenir; hedefe ait olmayanlar rakip
# hanesine tasinir ve satir yeniden olculur.
#
# KARAR PRIMERIN SONUCUNA BAKILMADAN VERILIR. Daraltma yalnizca olculen dizi
# kimligine dayanir ve KOSULSUZ benimsenir - satiri dusurse bile. "Esigi
# gecirdigi icin benimsendi" gerekcesi bu betikte HICBIR YERDE yoktur.
#
# En buyuk kumenin uye sayilmasi bir VARSAYIMDIR ve oyle isaretlenir: hangi
# kumenin hedef oldugu numune ici dizi kanitiyla degil kume buyuklugu ile
# secilmistir; dis referansla teyit I ve G asamalarinin isidir.
#
# Hicbir iki kutu birbiriyle >=%99 cikmiyorsa daraltma UYGULANMAZ: birini secmek
# keyfi olurdu. O hedef HETEROJEN isaretlenir ve karar dis referans teyidine
# birakilir.
# -------------------------------------------------------------------------
def yol2_uyelik_daralt(kons, uye_adlari, capa=None):
    """Uye kutulari konsensus kimligine gore kumele; capayi iceren kume gercek
    uye kumesidir. Donen: (yeni_uye, cikarilan, kanit_metni)"""
    d = {k: kons[k] for k in uye_adlari if k in kons and len(kons[k]) > 200}
    if len(d) < 2:
        return (list(uye_adlari), [], u'uye kutu sayisi 2\'den az - kumeleme yapilamaz')
    adlar = sorted(d)
    kimlik = {}
    for i, a in enumerate(adlar):
        for b in adlar[i + 1:]:
            kimlik[(a, b)] = kimlik[(b, a)] = hw_kimlik(d[a], d[b])
    # tek baglantili kumeleme, esik KIMLIK_ESIGI
    kume = {a: {a} for a in adlar}
    for (a, b), v in kimlik.items():
        if v >= KIMLIK_ESIGI and kume[a] is not kume[b]:
            yeni = kume[a] | kume[b]
            for x in yeni:
                kume[x] = yeni
    kumeler = []
    for a in adlar:
        if not any(kume[a] is k for k in kumeler):
            kumeler.append(kume[a])
    en_buyuk = max(kumeler, key=len)
    if capa and capa in kume:
        secilen = kume[capa]
    elif len(en_buyuk) >= 2:
        secilen = en_buyuk
    else:
        # hicbir iki kutu ayni organizma cikmadi -> daraltma ANLAMSIZ
        ozet = '; '.join('%d kutu' % len(k) for k in kumeler)
        return (list(uye_adlari), [],
                u'DARALTMA UYGULANMADI: %d uye kutunun hicbiri birbiriyle >=%%%s '
                u'kimlikte degil (kumeler = %s). Kutular ayni organizma degil; '
                u'birini secmek keyfi olurdu. Bu hedef HETEROJEN - once kutu '
                u'kimlikleri referansla dogrulanmali.'
                % (len(adlar), vir(KIMLIK_ESIGI, 1), ozet))
    cikan = [a for a in adlar if a not in secilen]
    ic = [kimlik[(a, b)] for i, a in enumerate(sorted(secilen))
          for b in sorted(secilen)[i + 1:]]
    ozet = '; '.join('%d kutu' % len(k) for k in sorted(kumeler, key=len, reverse=True))
    ick = (u'kume ici kimlik %%%s-%%%s' % (vir(min(ic)), vir(max(ic)))) if ic else \
          u'kumede tek kutu kaldi, kume ici karsilastirma yok'
    kanit = (u'%d uye kutu KONSENSUS DIZI KIMLIGINE gore kumelendi (esik %%%s): '
             u'kumeler = %s. En buyuk kume UYE sayildi (%d kutu, %s) - bu bir '
             u'VARSAYIMDIR: hangi kumenin hedef oldugu numune ici dizi kanitiyla '
             u'degil kume buyuklugu ile secildi; dis referansla teyidi I asamasinin '
             u'isidir. Kume disinda kalan ve RAKIP hanesine tasinan: %s. '
             u'(Bu karar primerin sonucuna HIC bakmadan verildi.)'
             % (len(adlar), vir(KIMLIK_ESIGI, 1), ozet, len(secilen), ick,
                ', '.join(cikan) if cikan else 'yok'))
    return (sorted(secilen) + [a for a in uye_adlari if a not in d], cikan, kanit)


# YOL 4 - ESLENIGI KALMIS SATIRLAR. Ayni uye kumesini (>=%80 ortusme) hedefleyen
# ve esigi GECEN baska bir cift varsa, esik alti satir gereksizdir: panelde
# tutulmasi plaka yeri israfidir. Bu bir kurtarma degil, bir TEMIZLIKtir - satir
# "dusenlere tasindi" olarak isaretlenir.
def yol4_eslenik_bul(satirlar, uyelik, alias=None):
    """Ayni uye kumesini hedefleyen ve esigi GECEN baska bir cift var mi."""
    alias = alias or {}
    U = lambda h: uyelik.get(alias.get(h, h), {})
    gecen = [r for r in satirlar if (r.get('esik_gecti_mi') or '').startswith('ESIK USTU')]
    out = {}
    for r in satirlar:
        if (r.get('esik_gecti_mi') or '').startswith('ESIK USTU'):
            continue
        u1 = set(U(r['hedef']).get('uye', []))
        if not u1:
            continue
        for g in gecen:
            u2 = set(U(g['hedef']).get('uye', []))
            if not u2:
                continue
            ortak = len(u1 & u2) / float(max(len(u1 | u2), 1))
            if ortak >= 0.80:
                out[r['hedef']] = (g['hedef'], _f(g.get('ASIL_ayrim_mm1')), ortak)
                break
    return out


def _ayirt_onbellekli(U, uye_diz, rak_diz):
    """ayirt_edici_mi PRIMER basina onbelleklenir.

    DARBOGAZ DUZELTMESI: cift_akisi N ileri x M geri cift uretir ama ayirt_edici_mi
    PRIMERE bagliddir - ayni primer yuzlerce ciftte tekrar sorulur. Olculdu:
    cagri basina 0,030 sn; 60x60 izgara = 3600 cift x 2 cagri = 216 sn. Onbellekle
    yalniz 120 ayri primer sorulur = 3,6 sn. Tam kosuda etki daha da buyuk.
    """
    bellek = {}

    def sor(primer, geri=False):
        anahtar = (primer, geri)
        if anahtar not in bellek:
            bellek[anahtar] = bool(U.ayirt_edici_mi(primer, uye_diz, rak_diz, geri=geri)[0])
        return bellek[anahtar]
    return sor, bellek


# -------------------------------------------------------------------------
# YOL 3 - YENIDEN TASARIM + ARMS
#
# Kil payi kalan satirlar icin (mevcut kat >= KIL_PAYI_ALT) omurga konsensusu
# uzerinde yeni cift aranir. Aday primerler geometri suzgecinden gecirilir, sonra
# uye/rakip konsensuslerine karsi "ayirt edici mi" diye sorulur.
#
# ARMS = ileri primerin 3' sondan 2. ve 3. bazina KASITLI uyumsuzluk eklenmesi.
# Bu bir DEJENERE BAZ DEGILDIR: tek bir sabit dizidir, oligo sayisini artirmaz ve
# "dejenere baz kullanilmasin" kararini ihlal etmez. Amaci 3' ucundaki tek baz
# farkini uzatma icin belirleyici hale getirmektir.
#
# IKI ASAMALI ELEME: butun adaylar once ASIL olcutle (mm<=1) olculur, YAN olcut
# (mm<=3) yalniz basa gureseyen 25 adayda kosar. Is yariya iner ve karar sutunu
# yine asil olcuttur.
#
# "yalniz_ileri" kipi mevcut geri primeri KORUR ve sadece ileriyi degistirir
# (NL1 sorunu): bu durumda geri primerin omurgadaki yeri once tam, bulunamazsa
# <=1 uyumsuzlukla aranir.
# -------------------------------------------------------------------------
def yol3_yeniden_tasarim(kok, nm, hedef, uye, rakip, kons, mevcut_F, mevcut_R,
                         yalniz_ileri=False, aday_ust=400, tarama_ust=3000,
                         arms_ust=5, yaz=print):
    """Yeni cift ara + ARMS varyantlari. primer3 yoksa duzgunce atlar."""
    import importlib.util
    if importlib.util.find_spec('primer3') is None:
        return dict(durum='ATLANDI', adaylar=[],
                    sebep=u'primer3-py bu makinede kurulu degil - yeniden tasarim '
                          u'taramasi yapilamadi. Kurulum: '
                          u'pip3 install primer3-py --break-system-packages')
    try:
        from screening import yapilandirma as C, motor, uretec as U, geometri as G
        G.tm('ACGTACGTACGTACGTAC')
    except SystemExit:
        return dict(durum='ATLANDI', adaylar=[],
                    sebep=u'geometri modulu primer3 bulamadigi icin durdu', adaylar2=[])
    except Exception as e:
        return dict(durum='ATLANDI', adaylar=[],
                    sebep=u'yeniden tasarim baslatilamadi (%s)' % type(e).__name__)

    capa = None
    for k in uye:
        if k['kutu'] in kons and len(kons[k['kutu']]) > 500:
            capa = k['kutu']; break
    if not capa:
        return dict(durum='ATLANDI', sebep=u'kullanilabilir omurga konsensusu yok', adaylar=[])
    omurga = kons[capa]
    uye_diz = [kons[k['kutu']] for k in uye if k['kutu'] in kons]
    rak_diz = [kons[k['kutu']] for k in rakip if k['kutu'] in kons]

    yaz(u'      omurga: %s (%d bp), uye konsensus %d, rakip konsensus %d'
        % (capa, len(omurga), len(uye_diz), len(rak_diz)))
    ad = U.aday_primerler(omurga)
    yaz(u'      geometriyi gecen aday: %d ileri / %d geri' % (len(ad['F']), len(ad['R'])))

    # HIZ SINIRI: ayirt_edici_mi her aday icin butun konsensuslari tarar. Aday
    # sayisi on binlerce oldugunda bu saatler surer. Omurga BOYUNCA esit araliklarla
    # ornekleyerek tavana indiriyoruz - ilk N'i almak butun adaylari omurganin 5'
    # ucundan secmek olurdu.
    def sey(liste, tavan):
        if len(liste) <= tavan:
            return liste
        adim = len(liste) / float(tavan)
        return [liste[int(i * adim)] for i in range(tavan)]
    once = (len(ad['F']), len(ad['R']))
    ad['F'] = sey(ad['F'], tarama_ust)
    if not yalniz_ileri:
        ad['R'] = sey(ad['R'], tarama_ust)
    if (len(ad['F']), len(ad['R'])) != once:
        yaz(u'      tarama tavani: %d ileri / %d geri aday olculecek (omurga boyunca '
            u'esit aralikli ornekleme)' % (len(ad['F']), len(ad['R'])))

    if yalniz_ileri:
        # mevcut geri primeri KORU, yalniz ileriyi degistir (NL1 sorunu).
        # cift_akisi geri primerin OMURGADAKI yerini ister; rc(R)'yi omurgada ara.
        hedef_diz = motor.rc(mevcut_R)
        iR = omurga.find(hedef_diz)
        if iR < 0:                       # <=1 uyumsuzlukla ara
            L = len(hedef_diz)
            for j in range(len(omurga) - L + 1):
                if sum(1 for a, b in zip(omurga[j:j + L], hedef_diz) if a != b) <= 1:
                    iR = j; break
        if iR < 0:
            return dict(durum='ATLANDI', adaylar=[],
                        sebep=u'mevcut geri primerin omurgadaki baglanma yeri '
                              u'bulunamadi - "yalniz ileri primeri degistir" '
                              u'kipi uygulanamadi')
        ad['R'] = [(iR, len(mevcut_R), mevcut_R, G.olc(mevcut_R))]
        yaz(u'      mevcut geri primer korunuyor, omurgadaki yeri: %d' % iR)

    secilen = []
    bakilan = 0
    t0 = time.time()
    sor, bellek = _ayirt_onbellekli(U, uye_diz, rak_diz)
    BAKILAN_UST = max(20000, aday_ust * 200)   # kabul edilen degil BAKILAN cift tavani
    for t in U.cift_akisi(ad):
        bakilan += 1
        if bakilan > BAKILAN_UST:
            yaz(u'      (bakilan cift tavani %d asildi - tarama durduruldu)' % BAKILAN_UST)
            break
        if bakilan % 5000 == 0:
            print('      ... %d cift tarandi, %d ayirt edici (%s)          '
                  % (bakilan, len(secilen), sure_metni(time.time() - t0)), end='\r', flush=True)
        c = U.cift_yap(t)
        if not sor(c['F']):
            continue
        if not sor(c['R'], True):
            continue
        secilen.append(c)
        if len(secilen) >= aday_ust:
            break
    yaz(u'      %d cift tarandi -> %d ayirt edici aday (%s, %d ayri primer sorgulandi)'
        % (bakilan, len(secilen), sure_metni(time.time() - t0), len(bellek)))

    # ARMS varyantlari: en iyi birkac adayin ileri primerine + MEVCUT cifte.
    # (Kasitli uyumsuzluk dejenere baz DEGILDIR; oligo sayisini artirmaz.)
    arms = []
    for c in secilen[:arms_ust]:
        for v, etiket in U.arms_varyantlari(c['F']):
            arms.append(dict(F=v, R=c['R'], urun=c['urun'], arms='F ' + etiket))
    for v, etiket in U.arms_varyantlari(mevcut_F):
        arms.append(dict(F=v, R=mevcut_R, urun=0, arms='F ' + etiket + ' (mevcut cift)'))

    # IKI ASAMA: once hepsi ASIL olcutle (mm<=1) elenir; YAN olcut (mm<=3)
    # yalniz basa gureseneler icin olculur - is yarisina iner.
    hepsi = secilen + arms
    yaz(u'      olculecek aday: %d (%d cift + %d ARMS varyanti)'
        % (len(hepsi), len(secilen), len(arms)))
    ilk = []
    for j, c in enumerate(hepsi, 1):
        if j % 10 == 0:
            print('      ... aday %d/%d olculuyor          ' % (j, len(hepsi)), end='\r', flush=True)
        o = nm.olc(c['F'], c['R'], uye, rakip, lo=URUN_ALT, hi=URUN_UST, mm=1)
        if not o or o.get('kat_enkotu') is None:
            continue
        ilk.append((o['kat_enkotu'], c, o))
    ilk.sort(key=lambda x: -x[0])
    sonuc = []
    for kat1, c, o in ilk[:25]:
        o3 = nm.olc(c['F'], c['R'], uye, rakip, lo=URUN_ALT, hi=URUN_UST, mm=3)
        sonuc.append(dict(F=c['F'], R=c['R'], urun=c.get('urun', 0),
                          arms=c.get('arms', ''), kat1=kat1,
                          kat3=(o3 or {}).get('kat_enkotu'),
                          kapsam=o.get('uye_kapsam_pay', '')))
    return dict(durum='TARANDI', sebep='', adaylar=sonuc,
                taranan=len(hepsi), omurga=capa)


# --------------------------------------------------------------- YOL 5
# Cok lokuslu arama: konsensusun TAMAMINI bolgelere ayirip her bolgede AYRI
# tasarim dener. Petriella'yi kurtaran sey buydu - ITS'te cozum yoktu, LSU'da
# vardi. Tek bir omurga penceresine sikismayi onler.
#
# BOLGE SINIRLARI NASIL BULUNUYOR
# -------------------------------
# Korunmus CAPA dizileriyle. Capalar, on yillardir evrensel primer olarak
# kullanilan ve tam da korunmus olduklari icin secilmis bolgelerdir; nanopore
# konsensusunde hata payi oldugu icin <=3 uyumsuzlukla ve IUPAC farkindalikli
# aranirlar (screening/motor.find_sites). Bir capa bulunamazsa o sinir
# "bulunamadi" olarak isaretlenir ve bolge YEDEK YOLLA (oransal pencere) kurulur;
# rapor hangi bolgenin capayla hangisinin yedekle kuruldugunu ACIKCA yazar.
CAPALAR = [
    # (ad, dizi, aciklama)  - hepsi 5'->3', sense yonunde aranir
    ('SSU_baslangic', 'AGAGTTTGATCMTGGCTCAG',  u'27F - bakteri/arke 16S basi'),
    ('SSU_orta',      'GTGYCAGCMGCCGCGGTAA',   u'515F - SSU V4 basi (16S ve 18S)'),
    ('SSU_son',       'TACGGYTACCTTGTTACGACTT', u'1492R (sense) - SSU sonu'),
    ('ITS1_baslangic','CTTGGTCATTTAGAGGAAGTAA', u'ITS1F - ITS1 basi'),
    ('58S_baslangic', 'GCATCGATGAAGAACGCAGC',  u'ITS3 - 5.8S sonu / ITS2 basi'),
    ('58S_son',       'GCTGCGTTCTTCATCGATGC',  u'ITS2 (ters tumleyen) - 5.8S basi'),
    ('LSU_baslangic', 'GCATATCAATAAGCGGAGGAAAAG', u'NL1 - LSU D1 basi'),
    ('LSU_D2_son',    'GGTCCGTGTTTCAAGACGG',   u'NL4 - D2 sonu'),
    ('LSU_ic',        'TCCTCCGCTTATTGATATGC',  u'ITS4 - LSU basi (5.8S sonrasi)'),
]


def capa_bul(motor, dizi, max_mm=3):
    """Capalari konsensuste ara. Donen: {ad: konum} (bulunanlar)."""
    out = {}
    try:
        enc = motor.encode(dizi)
    except Exception:
        return out
    for ad, d, _acik in CAPALAR:
        try:
            y = motor.find_sites(enc, d, max_mm, need_tail=False)
        except Exception:
            y = None
        if y:
            out[ad] = int(sorted(y, key=lambda x: x[1])[0][0])
    return out


def bolgeler_kur(motor, dizi, yaz):
    """Konsensusu bolgelere ayir. Donen: [(ad, bas, son, kaynak)]"""
    L = len(dizi)
    c = capa_bul(motor, dizi)
    b = []

    def ekle(ad, bas, son, kaynak):
        bas, son = max(0, int(bas)), min(L, int(son))
        if son - bas >= 200:
            b.append((ad, bas, son, kaynak))

    ssu_bas = c.get('SSU_baslangic', c.get('SSU_orta'))
    ssu_son = c.get('SSU_son')
    its1 = c.get('ITS1_baslangic')
    s58_bas = c.get('58S_son')
    s58_son = c.get('58S_baslangic')
    lsu = c.get('LSU_baslangic', c.get('LSU_ic'))
    d2son = c.get('LSU_D2_son')

    if ssu_bas is not None:
        ekle('SSU (16S/18S)', ssu_bas, ssu_son if ssu_son else (its1 or ssu_bas + 1600), 'capa')
    if its1 is not None:
        ekle('ITS1', its1, s58_bas if s58_bas else its1 + 300, 'capa')
    if s58_bas is not None and s58_son is not None:
        ekle('5.8S', s58_bas, s58_son + 20, 'capa')
    if s58_son is not None:
        ekle('ITS2', s58_son, lsu if lsu else s58_son + 350, 'capa')
    if lsu is not None:
        ekle('LSU D1-D2', lsu, d2son + 20 if d2son else lsu + 650, 'capa')
        if d2son is not None and L - d2son > 300:
            ekle('LSU kalani', d2son, L, 'capa')
    if not b:
        # YEDEK YOL: hicbir capa tutmadi - konsensusu esit parcalara bol
        n = max(2, min(6, L // 600))
        adim = L // n
        for i in range(n):
            ekle('bolge %d/%d (YEDEK - capa bulunamadi)' % (i + 1, n),
                 i * adim, min(L, (i + 1) * adim + 100), 'yedek')
    else:
        # CAPALARIN KAPSAMADIGI YERLER DE TARANIR. Yol 5'in butun amaci
        # konsensusun TAMAMINI denemek; capa tabanli bolgeler operonun bir
        # kismini disarida birakirsa (orn. arke A2 konsensusunda yalniz SSU
        # capasi tutuyor, 4,3 kb'nin 1,6 kb'si kapsaniyor) kalan parcalar
        # "kapsanmayan" adiyla eklenir.
        kapsanan = sorted((x, y) for _a, x, y, _k in b)
        bosluk, imlec = [], 0
        for x, y in kapsanan:
            if x - imlec >= 400:
                bosluk.append((imlec, x))
            imlec = max(imlec, y)
        if L - imlec >= 400:
            bosluk.append((imlec, L))
        for i, (x, y) in enumerate(bosluk, 1):
            ekle('kapsanmayan %d (capa disi)' % i, x, y, 'yedek')
        b.sort(key=lambda t: t[1])
    yaz(u'      bolgeler: %s' % '; '.join('%s %d-%d(%s)' % (a, x, y, k) for a, x, y, k in b))
    yaz(u'      bulunan capa: %s' % (', '.join(sorted(c)) or 'YOK'))
    return b, c


# YOL 5 - COK LOKUSLU ARAMA. Yol 3 tek bir omurga penceresine sikisir; bu yol
# konsensusun TAMAMINI bolgelere ayirip her bolgede AYRI tasarim dener.
# "Yapilamiyor" demeden once butun lokuslar denenmis olur.
def yol5_cok_lokuslu(kok, nm, hedef, uye, rakip, kons, aday_ust=150,
                     tarama_ust=800, yaz=print):
    """Her bolgede AYRI tasarim denemesi. Donen: bolge bolge rapor."""
    import importlib.util
    if importlib.util.find_spec('primer3') is None:
        return dict(durum='ATLANDI', bolge=[],
                    sebep=u'primer3-py kurulu degil - cok lokuslu arama yapilamadi')
    from screening import motor, uretec as U

    # OMURGA = EN UZUN uye konsensusu (2026-08-06 duzeltmesi).
    # Eskiden 800 bp'yi gecen ILK kutu aliniyordu. Uye kumesi hem A1 (16S
    # amplikonu, ~1,4 kb) hem A2 (tam operon, ~4,3 kb) kutulari tasidiginda
    # liste sirasi A1'i one koyuyor ve yol 5 operonun 2,9 kb'sini HIC gormuyordu:
    # "butun lokuslar denendi" derken aslinda yalniz 16S taranmisti. Yol 5'in
    # varlik sebebi en genis diziyi bolgelere ayirmaktir, o yuzden EN UZUNU.
    capa_kutu = None
    _en = 0
    for k in uye:
        _k = k['kutu']
        if _k in kons and len(kons[_k]) > 800 and len(kons[_k]) > _en:
            capa_kutu, _en = _k, len(kons[_k])
    if not capa_kutu:
        return dict(durum='ATLANDI', bolge=[],
                    sebep=u'800 bp ustu konsensus yok - bolgelere ayrilamaz')
    omurga = kons[capa_kutu]
    yaz(u'      omurga: %s (%d bp)' % (capa_kutu, len(omurga)))
    bolge, capalar = bolgeler_kur(motor, omurga, yaz)

    uye_diz = [kons[k['kutu']] for k in uye if k['kutu'] in kons]
    rak_diz = [kons[k['kutu']] for k in rakip if k['kutu'] in kons]
    rapor = []
    for ad, bas, son, kaynak in bolge:
        alt = omurga[bas:son]
        ad_p = U.aday_primerler(alt)
        def sey(l, t):
            if len(l) <= t:
                return l
            a = len(l) / float(t)
            return [l[int(i * a)] for i in range(t)]
        ad_p['F'] = sey(ad_p['F'], tarama_ust)
        ad_p['R'] = sey(ad_p['R'], tarama_ust)
        secilen = []
        sor, _b = _ayirt_onbellekli(U, uye_diz, rak_diz)
        bakilan = 0
        BAKILAN_UST = max(8000, aday_ust * 100)
        for t in U.cift_akisi(ad_p):
            bakilan += 1
            if bakilan > BAKILAN_UST:
                break
            c = U.cift_yap(t)
            if not sor(c['F']):
                continue
            if not sor(c['R'], True):
                continue
            secilen.append(c)
            if len(secilen) >= aday_ust:
                break
        en_iyi = None
        for c in secilen:
            o = nm.olc(c['F'], c['R'], uye, rakip, lo=URUN_ALT, hi=URUN_UST, mm=1)
            if not o or o.get('kat_enkotu') is None:
                continue
            if en_iyi is None or o['kat_enkotu'] > en_iyi['kat1']:
                o3 = nm.olc(c['F'], c['R'], uye, rakip, lo=URUN_ALT, hi=URUN_UST, mm=3)
                en_iyi = dict(F=c['F'], R=c['R'], urun=c['urun'], kat1=o['kat_enkotu'],
                              kat3=(o3 or {}).get('kat_enkotu'),
                              kapsam=o.get('uye_kapsam_pay', ''))
        rapor.append(dict(bolge=ad, bas=bas, son=son, kaynak=kaynak,
                          uzunluk=son - bas, aday=len(secilen), en_iyi=en_iyi))
        yaz(u'      %-32s %4d-%4d  aday %3d  en iyi %s'
            % (ad, bas, son, len(secilen),
               ('%s x' % vir(en_iyi['kat1'])) if en_iyi else '-'))
    return dict(durum='TARANDI', bolge=rapor, omurga=capa_kutu,
                capalar=sorted(capalar), sebep='')


# ---------------------------------------------------------------- surucu
# -------------------------------------------------------------------------
# HAZIRLIK VE SURUCU. Sira: girdi denetimi -> panel/uyelik/konsensus okuma ->
# esik alti satirlarin secimi -> panelsiz taleplerin eklenmesi -> okuma
# havuzlarinin kurulmasi -> _tur().
#
# Okuma tavani (--okuma) TEK PROTOKOL ile AYNI olmalidir: farkli derinlikte
# olculmus kat sayilari birbiriyle karsilastirilamaz ve "kurtarildi" karari
# derinlik farkindan dogabilir.
#
# PANELSIZ TALEPLER: toplantida istenip panele hic girmemis hedefler. Zincir
# panelden basladigi icin bunlar hic gorunmez ve "denenmeden yapilamadi" diye
# kalirdi. Burada omurga olarak kutunun KENDI konsensusu alinir ve en azindan bir
# tasarim denemesi yapilir.
# -------------------------------------------------------------------------
def calistir(kok, aday_ust, yalniz, sifirla, tarama_ust=3000, okuma=OKUMA_TAVANI,
             arms_ust=5, panelsiz_atla=False):
    os.environ['_KURTARMA_KOK'] = kok
    sys.path.insert(0, kok)
    from screening import numune as N, hedefler as H

    CIKTI = os.path.join(kok, 'KURTARMA_SONUC')
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
        print(s, flush=True); gunluk.write(s + '\n'); gunluk.flush()

    yaz('=' * 78)
    yaz('  verification TURU - esik altinda kalan satirlar')
    yaz('  surum %s   %s' % (VERSIYON, time.strftime('%Y-%m-%d %H:%M')))
    yaz('=' * 78)
    yaz('  ESIK %0.0fx DEGISTIRILMEZ. Yol 1 esigi dusurmez, OLCUYU degistirir;' % ESIK)
    yaz('  gerekcesi raporun basinda yazilidir.')
    yaz('')

    rc = girdi_denetle(yaz, 'K (verification)', [
        (os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'panel_tek_protokol.tsv'),
         'P asamasinin panel tablosu', 'P')])
    if rc:
        return rc
    satirlar, tp_yolu = tek_protokol_oku(kok)
    uy_yol = uyelik_dosyasi(kok)
    if not uy_yol:
        sys.exit('HATA: uyelik_yeniden_turetme_uyelik_*.tsv yok. Once secenek (U).')
    uyelik = uyelik_oku(uy_yol)
    kons = {d['kutu']: d['dizi'] for d in H.konsensusler()}
    kut = {k['kutu']: k for k in H.kutular()}

    hedefler = [r for r in satirlar
                if not (r.get('esik_gecti_mi') or '').startswith('ESIK USTU')]
    # panelde satiri OLMAYAN toplanti talepleri de kurtarma kapsamina alinir
    panelsiz = []
    for t in ([] if panelsiz_atla else PANELSIZ_TALEPLER):
        panelsiz.append(dict(hedef=t['hedef'], sinif=t['sinif'], F='', R='',
                             urun_bp='', ASIL_ayrim_mm1='', ASIL_kapsam_mm1='',
                             esik_gecti_mi='PANELDE SATIR YOK', _panelsiz=t))
    hedefler = hedefler + panelsiz
    if yalniz:
        hedefler = [r for r in hedefler
                    if any(y.strip().lower() in r['hedef'].lower()
                           for y in yalniz.split(','))]
    yaz('  TEK PROTOKOL kaynagi : %s' % os.path.basename(tp_yolu))
    yaz('  uyelik kaynagi       : %s' % os.path.basename(uy_yol))
    yaz('  esik alti satir      : %d' % len(hedefler))
    yaz('')

    alias = takma_adlar(kok)
    eslenik = yol4_eslenik_bul(satirlar, uyelik, alias)

    # okunacak kutular
    gerekli = {}
    for r in hedefler:
        u = uyelik.get(alias.get(r['hedef'], r['hedef']))
        if not u:
            continue
        for ad in u['uye'] + u['rakip'] + u['karisik']:
            if ad in kut:
                gerekli[ad] = kut[ad]
        if not u['rakip']:
            for k in kut.values():
                if k['sinif'] == (u['sinif'] or r.get('sinif')):
                    gerekli[k['kutu']] = k
    # panelsiz taleplerin kutulari da havuza girmeli - yoksa nm.olc onlari
    # goremez ve tasarim denemesi sessizce sifir olcumle biter
    for t in PANELSIZ_TALEPLER:
        for ad in t['uye']:
            if ad in kut:
                gerekli[ad] = kut[ad]
        for k in kut.values():
            if k['sinif'] == t['sinif']:
                gerekli[k['kutu']] = k
    yaz('  okunacak kutu        : %d' % len(gerekli))

    def ilerK(i, n, ad):
        print('   ... okuma havuzu %d/%d  %s          ' % (i, n, ad), end='\r', flush=True)

    t0 = time.time()
    yaz('')
    yaz('Okuma havuzlari kuruluyor. Bu adimda ekranda yalniz kutu adlari akar,')
    yaz('takilmis DEGILDIR; asil is bundan sonra baslar ve her hedef ayri kaydedilir.')
    nm = N.Numune(list(gerekli.values()), n=okuma, ilerle=ilerK, otorite=True)
    top = sum(h.n_okuma for h in nm.havuz.values())
    yaz('\nHavuzlar hazir: %d kutu, %d okuma (%s)' % (len(gerekli), top, sure_metni(time.time() - t0)))
    tasarim_n = sum(1 for r in hedefler
                    if (_f(r.get('ASIL_ayrim_mm1')) or 0) >= KIL_PAYI_ALT
                    and r['hedef'] not in BILINEN)
    tahmin = len(hedefler) * 20 + tasarim_n * max(300, top / 60.0)
    yaz('TAHMINI SURE: ~%s  (yol 3 tasarim taramasi %d satirda kosacak)'
        % (sure_metni(tahmin), tasarim_n))
    yaz('')
    return _tur(kok, CIKTI, KONTROL, yaz, nm, hedefler, uyelik, kons, kut,
                eslenik, aday_ust, gunluk, alias, tarama_ust, arms_ust, okuma)


# -------------------------------------------------------------------------
# HEDEF HEDEF verification DONGUSU. Her hedef icin sira sabittir:
#   panelsiz talep -> BILINEN (onceden olculmus cikmaz) -> yol 4 -> yol 1 ->
#   yol 2 -> yol 3 -> yol 5
#
# BILINEN tablosundaki satirlar (Proteiniphilum, M. mazei grubu) YENIDEN
# DENENMEZ: sebepleri daha once olculmustur, satirda gosterilir ve zaman
# harcanmaz.
#
# KONTROL NOKTASI MUHRU (_ayar): okuma derinligi, aday tavani, tarama tavani,
# ARMS tavani ve betik surumu muhre dahildir. Farkli derinlikte olculmus bir
# sonucun sessizce yeniden kullanilmasi, bu zincirin duzeltmeye calistigi hatanin
# ta kendisidir.
# -------------------------------------------------------------------------
def _tur(kok, CIKTI, KONTROL, yaz, nm, hedefler, uyelik, kons, kut, eslenik,
         aday_ust, gunluk, alias=None, tarama_ust=3000, arms_ust=5,
         okuma=OKUMA_TAVANI):
    alias = alias or {}
    def kp(ad):
        t = ''.join(ch if ch.isalnum() else '_' for ch in ad)
        return os.path.join(KONTROL, 'hedef_%s.json' % t)

    AYAR = dict(okuma=okuma, aday_ust=aday_ust, tarama_ust=tarama_ust,
                arms_ust=arms_ust, surum=VERSIYON)
    # 2026-08-10 DIZI MUHRU: kontrol noktasi anahtarina primer DIZISI
    # dahil edilmemisti. Bir ciftin ileri/geri dizisi degistiginde eski
    # sonuc sessizce yeniden kullaniliyordu. Ayni hata P asamasinda iki
    # tam kosuyu bosa harcatti (5 sa 29 dk + 2 sa 0 dk). Artik dizi de
    # anahtara giriyor; dizi degisirse kontrol noktasi gecersiz olur.
    import hashlib as _hl

    def _ayar_of(r):
        d = dict(AYAR)
        d['dizi'] = _hl.md5(((r.get('F') or '') + '|' + (r.get('R') or ''))
                            .encode('utf-8')).hexdigest()[:12]
        return d
    sonuc = []
    tb = time.time()
    for i, r in enumerate(hedefler, 1):
        hedef = r['hedef']
        yol = kp(hedef)
        if os.path.exists(yol):
            try:
                _v = json.load(open(yol, encoding='utf-8'))
                # O-9: ayar muhru - farkli derinlikte olculmus sonuc yeniden
                # kullanilmamali (TEK PROTOKOL modulunun var olus sebebi bu).
                if _v.get('_ayar') == _ayar_of(r):
                    sonuc.append(_v)
                    yaz('[%2d/%2d] %-44s (onceki kosudan alindi)'
                        % (i, len(hedefler), hedef[:44]))
                    continue
                yaz('[%2d/%2d] %-44s (ayar degismis, yeniden olculuyor)'
                    % (i, len(hedefler), hedef[:44]))
            except Exception:
                pass

        eski = _f(r.get('ASIL_ayrim_mm1'))
        s = dict(_ayar=_ayar_of(r), hedef=hedef, kaynak=r.get('kaynak', ''), eski=eski,
                 eski_kapsam=r.get('ASIL_kapsam_mm1', ''), yollar=[], yeni=None,
                 gecti='HAYIR', olcu=u'ayrim kati (%s)' % _C.esik_metni(),
                 sebep='', ayrinti={})
        yaz('[%2d/%2d] %s  (eski %s x, kapsam %s)'
            % (i, len(hedefler), hedef, vir(eski), s['eski_kapsam']))

        pz = r.get('_panelsiz')
        if pz:
            u = dict(uye=pz['uye'], rakip=[], karisik=[], sinif=pz['sinif'])
        else:
            u = uyelik.get(alias.get(hedef, hedef)) or dict(uye=[], rakip=[], karisik=[],
                                                            sinif=r.get('sinif', ''))
        coz = lambda ad: [kut[a] for a in ad if a in kut]
        uye = coz(u['uye']); rakip = coz(u['rakip']) + coz(u['karisik'])
        if not rakip:
            rakip = [k for k in kut.values()
                     if k['sinif'] == (u['sinif'] or r.get('sinif'))
                     and k['kutu'] not in set(u['uye'])]

        # --- PANELSIZ TALEP: dogrudan tasarim denemesi ---
        if pz:
            s.update(olcu=u'panelde satiri yok - kutudan tasarim denemesi',
                     yollar=[u'yol 3 - panelsiz talep, kutu konsensusundan tasarim'])
            yaz(u'      -> PANELSIZ TALEP (%s): %s' % (pz['karar'], pz['not_'][:90]))
            t = yol3_yeniden_tasarim(kok, nm, hedef, uye, rakip, kons, '', '',
                                     False, aday_ust, tarama_ust, arms_ust, yaz)
            s['ayrinti']['yol3'] = t
            en_iyi = t['adaylar'][0] if t.get('adaylar') else None
            if en_iyi and en_iyi['kat1'] >= ESIK:
                s.update(yeni=u'YENI CIFT %s / %s (%d bp) %s x'
                              % (en_iyi['F'], en_iyi['R'], en_iyi['urun'], vir(en_iyi['kat1'])),
                         gecti='EVET (yeni cift)',
                         sebep=u'Panelde satiri yoktu; kutu konsensusundan tasarim '
                               u'denendi ve esigi gecen aday bulundu. %s' % pz['not_'])
                yaz(u'         BULUNDU: %s x' % vir(en_iyi['kat1']))
            else:
                yaz(u'         tek pencerede yok - YOL 5: cok lokuslu arama')
                t5 = yol5_cok_lokuslu(kok, nm, hedef, uye, rakip, kons,
                                      aday_ust=min(aday_ust, 40),
                                      tarama_ust=min(tarama_ust, 200), yaz=yaz)
                s['yollar'].append(u'yol 5 - cok lokuslu arama (%s)' % t5['durum'])
                s['ayrinti']['yol5'] = t5
                iyi = [b for b in t5.get('bolge', [])
                       if b.get('en_iyi') and b['en_iyi']['kat1'] >= ESIK]
                if iyi:
                    en = max(iyi, key=lambda b: b['en_iyi']['kat1']); e = en['en_iyi']
                    s.update(yeni=u'YENI CIFT (%s bolgesi) %s / %s (%d bp) %s x'
                                  % (en['bolge'], e['F'], e['R'], e['urun'], vir(e['kat1'])),
                             gecti='EVET (yeni cift)',
                             sebep=u'%s Panelde satiri yoktu; %s bolgesinde cozum bulundu.'
                                   % (pz['not_'], en['bolge']))
                else:
                    s['sebep'] = (u'%s DENENDI (tek pencere + %d bolgede cok lokuslu '
                                  u'arama): esigi gecen aday yok.'
                                  % (pz['not_'], len(t5.get('bolge', []))))

        # --- bilinen, tekrar denenmeyecek ---
        elif hedef in BILINEN:
            b = BILINEN[hedef]
            s.update(gecti='HAYIR', sebep=b['sebep'], yollar=[b['yol']])
            yaz(u'      -> %s  (%s)' % (b['sonuc'], b['yol']))

        # --- YOL 4: eslenik ---
        elif hedef in eslenik:
            g, kat, ort = eslenik[hedef]
            s.update(gecti='DUSENLERE TASINDI', olcu='eslenik',
                     sebep=u'Ayni uye kumesini hedefleyen ve esigi GECEN baska bir '
                           u'cift var: "%s" (%s x, uye kumesi ortusme %%%d). Bu satir '
                           u'artik gereksiz; panelde tutulmasi plaka yeri israfidir.'
                           % (g, vir(kat), int(100 * ort)),
                     yollar=[u'yol 4 - eslenigi kalmis satir'])
            yaz(u'      -> YOL 4: eslenigi var (%s, %s x) - dusenlere tasinir' % (g, vir(kat)))

        # --- YOL 1: evrensel ---
        elif evrensel_mi(hedef, r.get('duzey', '')):
            o = yol1_evrensel(nm, uye, rakip, r['F'], r['R'])
            s.update(olcu=u'KAPSAMA + ALAN DISI (ayrim kati bu satirda tanimsiz)',
                     yollar=[u'yol 1 - olcu duzeltildi'],
                     yeni=u'kapsama %s (%%%d), alan disi %%%s'
                          % (o['kapsam_pay'], int(100 * o['kapsama']), vir(o['alandisi'])),
                     gecti='EVET' if o['gecti'] else 'HAYIR',
                     ayrinti=o)
            if not o['gecti']:
                s['sebep'] = (u'Kapsama %%%d (olcut %%%d) / alan disi %%%s (olcut en cok %%%.0f).'
                              % (int(100 * o['kapsama']), int(100 * EVRENSEL_KAPSAMA_ESIGI),
                                 vir(o['alandisi']), EVRENSEL_ALANDISI_UST))
            yaz(u'      -> YOL 1: kapsama %s, alan disi %%%s  => %s'
                % (o['kapsam_pay'], vir(o['alandisi']), s['gecti']))

        else:
            # --- YOL 2: uyelik daraltma ---
            kapsam_tam = (r.get('ASIL_kapsam_mm1') or '').split('/')
            tam = (len(kapsam_tam) == 2 and kapsam_tam[0] == kapsam_tam[1])
            # ---------------------------------------------------------------
            # UYELIK KOSULSUZ BENIMSENIR - VE DUSUS BIR KAYIP DEGILDIR.
            # Daraltma sonrasi olculen deger eskisinden DUSUK cikabilir. Bu,
            # primerin kotulesmesi degil OLCUNUN duzelmesidir: eski deger yanlis
            # uyelikten geliyordu (hedefe ait olmayan kutular uye sayilmis ya da
            # ayni organizma olan kutular rakip hanesinde birakilmisti) ve hicbir
            # zaman gecerli degildi. Satirin dusmesine izin verilmesi, kuralin tek
            # yonlu calismadiginin kanitidir.
            # ---------------------------------------------------------------
            # --- YOL 2: UYELIK DARALTMA ---
            # KRITIK DUZELTME (tasarim incelemesi madde 1): uyelik ARTIK primerin
            # sonucuna gore benimsenmiyor. Daraltma yalniz OLCULEN DIZI KIMLIGINE
            # gore yapilir ve KOSULSUZ benimsenir - sonuc hedefi dusurse bile.
            # "Esigi gecirdigi icin benimsendi" gerekcesi hicbir yerde YOKTUR.
            if len(u['uye']) > 1:
                yeni_uye, cikan, kanit = yol2_uyelik_daralt(kons, u['uye'])
                s['ayrinti']['yol2'] = dict(kanit=kanit, cikan=cikan,
                                            uygulandi=bool(cikan))
                if cikan:
                    uye2 = coz(yeni_uye)
                    rakip2 = rakip + coz(cikan)
                    o = nm.olc(r['F'], r['R'], uye2, rakip2, lo=URUN_ALT, hi=URUN_UST, mm=1)
                    o3 = nm.olc(r['F'], r['R'], uye2, rakip2, lo=URUN_ALT, hi=URUN_UST, mm=3)
                    kat = (o or {}).get('kat_enkotu')
                    kat3 = (o3 or {}).get('kat_enkotu')
                    # KOSULSUZ BENIMSEME - yon ne olursa olsun
                    uye, rakip = uye2, rakip2
                    s['yollar'].append(u'yol 2 - uyelik daraltildi (KOSULSUZ benimsendi)')
                    s['ayrinti']['yol2'].update(kat1=kat, kat3=kat3,
                                                kapsam=(o or {}).get('uye_kapsam_pay'))
                    dus = (kat is not None and eski is not None and kat < eski)
                    yon = (u'DEGISMEDI' if (kat is None or eski is None) else
                           (u'YUKSELDI' if kat > eski else
                            u'DUSTU' if kat < eski else u'DEGISMEDI'))
                    yaz(u'      -> YOL 2: %s' % kanit)
                    yaz(u'         daraltilmis uyelikle: %s x (eski %s x, %s) - '
                        u'KOSULSUZ benimsendi' % (vir(kat), vir(eski), yon))
                    if dus:
                        s['dusus_notu'] = (
                            u'DIKKAT - BU BIR KAYIP DEGIL, DUZELTMEDIR. Bu satirin '
                            u'eski %s x degeri YANLIS UYELIKTEN geliyordu: hedefe ait '
                            u'olmayan kutular uye sayilmis, ya da ayni organizma olan '
                            u'kutular rakip hanesinde birakilmisti. Uyelik olculen dizi '
                            u'kimligine gore duzeltilince gercek deger %s x cikti. '
                            u'Dusus, primerin kotulesmesi degil OLCUNUN duzelmesidir; '
                            u'eski deger hicbir zaman gecerli degildi. Uyelik karari '
                            u'primerin sonucuna BAKILMADAN verildi - bu satirin dusmesi '
                            u'kuralın tek yonlu calismadiginin kanitidir.'
                            % (vir(eski), vir(kat)))
                        yaz(u'         NOT: dusus bir KAYIP DEGIL, olcunun duzelmesi. '
                            u'Eski %s x yanlis uyelikten geliyordu.' % vir(eski))
                    s['uyelik_gerekcesi'] = (
                        u'Uye kumesi YALNIZ olculen konsensus kimligine gore '
                        u'belirlendi (esik %%%s), primerin sonucundan BAGIMSIZ olarak '
                        u've kosulsuz benimsendi. Kanit: %s. Yeni deger %s x '
                        u'(eski %s x, %s) - bu deger benimseme kararini ETKILEMEDI.'
                        % (vir(KIMLIK_ESIGI, 1), kanit, vir(kat), vir(eski), yon))
                    s['eski'] = eski = kat      # bundan sonrasi duzeltilmis uyelikle
                    if kat is not None and kat >= ESIK:
                        s.update(yeni='%s x' % vir(kat), gecti='EVET',
                                 sebep=s['uyelik_gerekcesi'])
                    else:
                        s['sebep'] = s['uyelik_gerekcesi']
                else:
                    s['yollar'].append(u'yol 2 - daraltma uygulanmadi')
                    s['sebep'] = kanit
                    s['uyelik_gerekcesi'] = kanit
                    yaz(u'      -> YOL 2: %s' % kanit)

            # --- YOL 3: yeniden tasarim + ARMS ---
            # KIL_PAYI_ALT kapisi YOL 3 icin dogrudur: ayni omurga penceresinde
            # primer oynatmak 0,8x'i 8x'e cikarmaz. YOL 5 icin YANLISTI - bkz.
            # asagidaki not.
            if s['gecti'] != 'EVET' and (eski or 0) >= KIL_PAYI_ALT:
                yalniz_ileri = ('microasca' in hedef.lower())
                yaz(u'      -> YOL 3: yeniden tasarim%s' %
                    (u' (yalniz ILERI primer degistirilecek - NL1)' if yalniz_ileri else ''))
                t = yol3_yeniden_tasarim(kok, nm, hedef, uye, rakip, kons,
                                         r['F'], r['R'], yalniz_ileri, aday_ust,
                                         tarama_ust, arms_ust, yaz)
                s['yollar'].append(u'yol 3 - yeniden tasarim + ARMS (%s)' % t['durum'])
                s['ayrinti']['yol3'] = t
                en_iyi = t['adaylar'][0] if t.get('adaylar') else None
                if en_iyi and en_iyi['kat1'] >= ESIK:
                    s.update(yeni=u'YENI CIFT %s / %s (%d bp) %s x%s'
                                  % (en_iyi['F'], en_iyi['R'], en_iyi['urun'],
                                     vir(en_iyi['kat1']),
                                     (u' [ARMS: %s]' % en_iyi['arms']) if en_iyi['arms'] else ''),
                             gecti='EVET (yeni cift)',
                             sebep=u'Mevcut cift esigi gecmiyor; taramada esigi gecen aday bulundu.')
                    yaz(u'         BULUNDU: %s x  %s / %s'
                        % (vir(en_iyi['kat1']), en_iyi['F'], en_iyi['R']))
                elif t['durum'] == 'ATLANDI':
                    s['sebep'] = s['sebep'] or t['sebep']
                    yaz(u'         %s' % t['sebep'])
                else:
                    s['sebep'] = s['sebep'] or (
                        u'Tarandi, esigi gecen aday yok (en iyi %s x).'
                        % (vir(en_iyi['kat1']) if en_iyi else '-'))
                    yaz(u'         esigi gecen aday yok')

                # --- YOL 5: COK LOKUSLU ARAMA ---
                # Yol 3 tek omurga penceresine dayanir. Petriella'da cozum ITS'te
                # degil LSU'daydi; "yapilamiyor" demeden once BUTUN lokuslar denenir.
                yaz(u'      -> YOL 5: cok lokuslu arama (bolge bolge)')
                t5 = yol5_cok_lokuslu(kok, nm, hedef, uye, rakip, kons,
                                      aday_ust=min(aday_ust, 40),
                                      tarama_ust=min(tarama_ust, 200), yaz=yaz)
                s['yollar'].append(u'yol 5 - cok lokuslu arama (%s)' % t5['durum'])
                s['ayrinti']['yol5'] = t5
                iyi = [b for b in t5.get('bolge', [])
                       if b.get('en_iyi') and b['en_iyi']['kat1'] >= ESIK]
                if iyi:
                    en = max(iyi, key=lambda b: b['en_iyi']['kat1'])
                    e = en['en_iyi']
                    s.update(yeni=u'YENI CIFT (%s bolgesi) %s / %s (%d bp) %s x'
                                  % (en['bolge'], e['F'], e['R'], e['urun'], vir(e['kat1'])),
                             gecti='EVET (yeni cift)',
                             sebep=u'Tek omurga penceresinde cozum yoktu; %s bolgesinde '
                                   u'bulundu. Taranan bolge: %d.'
                                   % (en['bolge'], len(t5.get('bolge', []))))
                    yaz(u'         BULUNDU: %s bolgesi, %s x' % (en['bolge'], vir(e['kat1'])))
                elif t5.get('bolge'):
                    s['sebep'] += (u' COK LOKUSLU ARAMA: %d bolgenin hepsi tarandi '
                                   u'(%s), hicbirinde esigi gecen aday yok.'
                                   % (len(t5['bolge']),
                                      ', '.join(b['bolge'] for b in t5['bolge'])))
            elif s['gecti'] != 'EVET':
                # ---------------------------------------------------------------
                # 2026-08-06 MANTIK HATASI DUZELTMESI.
                # Eskiden bu dalda HICBIR SEY denenmiyordu: taban KIL_PAYI_ALT'in
                # altindaysa satir "yeniden tasarim kosulmadi" notuyla birakiliyordu.
                # YOL 3 icin dogru, YOL 5 icin YANLIS:
                #   - yol 3 AYNI omurga penceresinde primer oynatir; mevcut kat
                #     0,8x ise o pencerede cozum yoktur, kapi yerindedir.
                #   - yol 5 BASKA BIR LOKUSA gecer. Mevcut kat, TERK EDILEN
                #     lokusun olcusudur; yeni lokus hakkinda hicbir sey soylemez.
                # Sonuc: tabani dusuk bes hedefte cok lokuslu arama HIC KOSMADI ve
                # rapor "ayrim yetmedi" diyordu - oysa "butun lokuslarda denendi"
                # diyebilmemiz gerekiyordu. Artik kapi yalniz yol 3'e uygulanir.
                # ---------------------------------------------------------------
                _uzun = max([len(kons[k['kutu']]) for k in uye if k['kutu'] in kons] or [0])
                if _uzun >= 2000:
                    yaz(u'      -> YOL 5: taban dusuk ama konsensus %d bp - '
                        u'BASKA LOKUSLAR denenir' % _uzun)
                    t5 = yol5_cok_lokuslu(kok, nm, hedef, uye, rakip, kons,
                                          aday_ust=min(aday_ust, 40),
                                          tarama_ust=min(tarama_ust, 200), yaz=yaz)
                    s['yollar'].append(u'yol 5 - cok lokuslu arama (%s)' % t5['durum'])
                    s['ayrinti']['yol5'] = t5
                    iyi = [b for b in t5.get('bolge', [])
                           if b.get('en_iyi') and b['en_iyi']['kat1'] >= ESIK]
                    if iyi:
                        en = max(iyi, key=lambda b: b['en_iyi']['kat1']); e = en['en_iyi']
                        s.update(yeni=u'YENI CIFT (%s bolgesi) %s / %s (%d bp) %s x'
                                      % (en['bolge'], e['F'], e['R'], e['urun'], vir(e['kat1'])),
                                 gecti='EVET (yeni cift)',
                                 sebep=u'Mevcut lokusta taban cok dusuktu; %s bolgesinde '
                                       u'esigi gecen aday bulundu. Taranan bolge: %d.'
                                       % (en['bolge'], len(t5.get('bolge', []))))
                        yaz(u'         BULUNDU: %s bolgesi, %s x' % (en['bolge'], vir(e['kat1'])))
                    else:
                        s['sebep'] = ((s['sebep'] + u'  ') if s['sebep'] else u'') + (
                            u'BUTUN LOKUSLAR DENENDI: %d bolge tarandi (%s), hicbirinde '
                            u'esigi gecen aday yok.'
                            % (len(t5.get('bolge', [])),
                               ', '.join(b['bolge'] for b in t5.get('bolge', [])) or '-'))
                elif not s['sebep']:
                    s['sebep'] = (u'Taban cok dusuk (%s x < %s x) ve konsensus %d bp - '
                                  u'tek lokus (16S) var, gidilecek baska bolge YOK.'
                                  % (vir(eski), vir(KIL_PAYI_ALT), _uzun))

        json.dump(s, open(yol, 'w', encoding='utf-8'), ensure_ascii=False, indent=1, default=str)
        sonuc.append(s)
        g = time.time() - tb
        print('        gecen %s | tahmini kalan %s'
              % (sure_metni(g), sure_metni(g / i * (len(hedefler) - i))), flush=True)

    raporla(CIKTI, sonuc, yaz)
    rc = cikti_denetle(yaz, 'K (verification)', [
        (os.path.join(CIKTI, 'kurtarma_satirlari.tsv'), 'kurtarma_satirlari.tsv')])
    gunluk.close()
    return rc


# ---------------------------------------------------------------------------
# Uc cikti: hedef basina tek satirlik tablo, aday primer tablosu ve markdown
# rapor. "olcu" sutunu her satirda HANGI olcunun uygulandigini soyler - evrensel
# satirlarda kapsama + alan disi, digerlerinde ayrim kati (10x). Esigin
# indirilmedigi ve neden indirilmedigi rapor basliginda yazilidir.
# ---------------------------------------------------------------------------
def raporla(CIKTI, sonuc, yaz):
    yol = os.path.join(CIKTI, 'kurtarma_satirlari.tsv')
    with open(yol, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# verification TURU - her hedef icin TEK satir.\n')
        fh.write(u'# UYELIK KURALI: uye kumesi YALNIZ olculen dizi kimligine gore\n')
        fh.write(u'# belirlenir ve KOSULSUZ benimsenir. Primerin sonucu (esigi gecip\n')
        fh.write(u'# gecmedigi) uyelik kararini ETKILEMEZ. "Esigi gecirdigi icin\n')
        fh.write(u'# benimsendi" gibi bir gerekce bu dosyada YOKTUR.\n')
        fh.write(u'# ESIK %0.0fx DEGISTIRILMEDI. "olcu" sutunu hangi olcunun uygulandigini soyler;\n' % ESIK)
        fh.write(u'# evrensel satirlarda ayrim kati TANIMSIZ oldugu icin KAPSAMA + ALAN DISI kullanilir.\n')
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'eski_deger', 'eski_kapsam', 'denenen_yol', 'olcu',
                    'yeni_deger', 'esigi_gecti_mi', 'UYELIK_GEREKCESI',
                    'DUSUS_KAYIP_MI_DUZELTME_MI', 'sebep'])
        for s in sonuc:
            w.writerow([s['hedef'], vir(s['eski']), s['eski_kapsam'],
                        ' + '.join(s['yollar']) or '-', s['olcu'],
                        s['yeni'] or '-', s['gecti'],
                        s.get('uyelik_gerekcesi', '-'), s['sebep']])
    yaz('  yazildi: %s' % yol)

    ay = os.path.join(CIKTI, 'yeni_adaylar.tsv')
    with open(ay, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# Yol 3 taramasinin urettigi adaylar (her hedefin en iyi 25 tanesi).\n')
        fh.write(u'# ARMS = 3\' sondan 2. ve 3. baza KASITLI uyumsuzluk. Dejenere baz DEGILDIR,\n')
        fh.write(u'# oligo sayisini artirmaz, toplanti kararini ihlal etmez - ama ayri bir maddedir.\n')
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'bolge', 'F', 'R', 'urun_bp', 'arms',
                    'ayrim_mm1', 'ayrim_mm3', 'kapsam'])
        for s in sonuc:
            for a in (s.get('ayrinti', {}).get('yol3', {}) or {}).get('adaylar', []):
                w.writerow([s['hedef'], 'tek pencere (yol 3)', a['F'], a['R'],
                            a['urun'], a['arms'], vir(a['kat1']), vir(a['kat3']), a['kapsam']])
            for b in (s.get('ayrinti', {}).get('yol5', {}) or {}).get('bolge', []):
                e = b.get('en_iyi')
                if e:
                    w.writerow([s['hedef'], '%s (%d-%d, %s)' % (b['bolge'], b['bas'],
                                                               b['son'], b['kaynak']),
                                e['F'], e['R'], e['urun'], '', vir(e['kat1']),
                                vir(e['kat3']), e['kapsam']])
                else:
                    w.writerow([s['hedef'], '%s (%d-%d, %s)' % (b['bolge'], b['bas'],
                                                               b['son'], b['kaynak']),
                                '', '', '', '', 'aday yok (%d taranan)' % b['aday'], '', ''])
    yaz('  yazildi: %s' % ay)

    gecen = [s for s in sonuc if s['gecti'].startswith('EVET')]
    tasi = [s for s in sonuc if s['gecti'] == 'DUSENLERE TASINDI']
    kalan = [s for s in sonuc if s not in gecen and s not in tasi]
    rp = os.path.join(CIKTI, 'KURTARMA_RAPORU.md')
    with open(rp, 'w', encoding='utf-8') as fh:
        fh.write(u'# Kurtarma turu\n\nUretim: %s · betik %s\n\n'
                 % (time.strftime('%Y-%m-%d %H:%M'), VERSIYON))
        fh.write(u'## Sonuc\n\n- Kurtarilan: **%d**\n- Dusenlere tasinan (eslenigi var): **%d**\n'
                 u'- Kurtarilamayan: **%d**\n\n' % (len(gecen), len(tasi), len(kalan)))
        fh.write(u'> **Esik indirilmedi.** Yol 1 evrensel satirlarda ayrim katinin '
                 u'yerine kapsama + alan disi olcusunu koyar; bu, esigin gevsetilmesi '
                 u'degildir - o satirlarda oranin paydasi tanimsizdir. Diger butun '
                 u'satirlarda 10x aynen uygulandi.\n\n')
        fh.write(u'```' + GEREKCE_EVRENSEL + u'```\n\n')
        fh.write(u'## Satir satir\n\n')
        fh.write(u'| hedef | eski | yol | yeni | gecti mi |\n|---|---|---|---|---|\n')
        for s in sonuc:
            fh.write(u'| %s | %s | %s | %s | %s |\n'
                     % (s['hedef'], vir(s['eski']), ' + '.join(s['yollar']) or '-',
                        (s['yeni'] or '-')[:70], s['gecti']))
        fh.write(u'\n## Kurtarilamayanlarin sebebi\n\n')
        for s in kalan:
            fh.write(u'**%s** — %s\n\n' % (s['hedef'], s['sebep'] or '-'))
        fh.write(u'\n## Okuma sirasi\n\n1. Bu dosya. 2. `kurtarma_satirlari.tsv` '
                 u'(her hedef tek satir). 3. `yeni_adaylar.tsv` (yol 3 ciktisi).\n')
    yaz('  yazildi: %s' % rp)
    yaz('')
    yaz('  KURTARILAN: %d   DUSENLERE TASINAN: %d   KURTARILAMAYAN: %d'
        % (len(gecen), len(tasi), len(kalan)))



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

# Komut satiri: --aday-ust ve --tarama-ust yol 3 taramasinin buyuklugunu,
# --arms-ust kac adayin ARMS varyantinin uretilecegini, --okuma derinlik tavanini
# (P ile ayni olmali), --panelsiz-atla hizli sinamada panelsiz talepleri atlamayi
# belirler.
def main():
    p = argparse.ArgumentParser(description='Esik alti satirlar icin kurtarma turu')
    p.add_argument('--kok', default='.')
    p.add_argument('--aday-ust', type=int, default=400,
                   help='yol 3 taramasinda olculecek en fazla aday cift')
    p.add_argument('--tarama-ust', type=int, default=3000,
                   help='yol 3 taramasinda omurgadan ornekle nen en fazla primer adayi')
    p.add_argument('--yalniz', default=None, help='yalniz adi bunu iceren hedefler (sinama)')
    p.add_argument('--arms-ust', type=int, default=5,
                   help='kac adayin ARMS varyantlari uretilsin')
    p.add_argument('--okuma', type=int, default=OKUMA_TAVANI,
                   help='kutu basina okuma tavani (TEK PROTOKOL ile ayni olmali)')
    p.add_argument('--panelsiz-atla', action='store_true',
                   help='panelde satiri olmayan talepleri atla (yalniz hizli test icin)')
    p.add_argument('--sifirla', action='store_true')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)
    if not os.path.isdir(os.path.join(kok, 'screening')):
        sys.exit('HATA: %s icinde screening yok. --kok ile proje klasorunu verin.' % kok)
    return calistir(kok, a.aday_ust, a.yalniz, a.sifirla, a.tarama_ust, a.okuma,
                    a.arms_ust, a.panelsiz_atla)


if __name__ == '__main__':
    sys.exit(main() or 0)
