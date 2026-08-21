# -*- coding: utf-8 -*-
"""Butun sabitler ve yollar tek yerde. Kullanici bu dosyayi duzenleyebilir.

qPCR kisitlari QIAGEN Rotor-Gene Q + QuantiNova SYBR Green icin SABIT tutuldu.
"""
# ---------------------------------------------------------------------------
# yapilandirma.py — butun dosya yollarinin ve sayisal sabitlerin tek tanim yeri.
#
# GIRDI  : yalnizca kendi konumunu okur; KOK, bu paketin bir ust dizinidir.
#          Baska dosya okumaz, olcum yapmaz.
# CIKTI  : dosyaya yazmaz. Modul duzeyinde sabitler disari acilir: yollar
#          (PANEL_TSV, HEDEFLER_TSV, KONSENSUS_KANONIK, FASTQ, REFDB, CIKTI,
#          KONTROL, ONBELLEK), qPCR kisitlari, degismez primer kurallari,
#          144 hucrelik parametre izgarasinin eksenleri, numune olcum ayarlari
#          ve huni kapasiteleri.
# CAGRAN : paketteki hemen her modul "from . import yapilandirma as C" ile ice
#          aktarir; disaridan TEK_PROTOKOL/tek_protokol_olc.py (tus P),
#          KURTARMA/kurtarma_turu.py (tus K) ve KURTARMA/dogrulama_turu.py
#          (tus D) de kullanir. Yani butun menu tuslarinda dolayli olarak yuklu.
#
# Bir esigi buradan degistirmek butun asamalari ayni anda etkiler. Ozellikle
# ENKOTU_ASGARI_OKUMA ve KAPSAM_ESIGI degerleri asamalar arasi
# karsilastirilabilirligi tasidigi icin kosu ortasinda degistirilmemelidir.
# ---------------------------------------------------------------------------
import os

# ---------------------------------------------------------------- yollar
# KOK = proje klasoru (bu paketin bir ust dizini)
KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def y(*p):
    return os.path.join(KOK, *p)

PANEL_TSV   = y('primer_final', 'devir_ciftleri_20260802_sonrotus_TESLIM.tsv')
PANEL_TSV_YEDEK = y('primer_final', 'devir_ciftleri_20260802_sonrotus.tsv')
HEDEFLER_TSV= y('WSL_betikleri', 'hedefler.tsv')
TAXID_ADLARI= y('WSL_betikleri', 'taxid_adlari.tsv')
# YON NORMALIZASYONU (2026-08-02): konsensusler artik TEK KANONIK klasorden okunur.
# 'consensus sequences' KARISIK yonludur (71 antisense / 27 sense - olculdu) ve
# dogrudan okunmasi yasaktir; ters yonlu bir konsensuste in-silico PCR SESSIZCE
# 0 urun verir (olculen kayip %100 - KAPSAMLI_ARAMA/yon_etki_testi.py).
# Kanonik klasor KAPSAMLI_ARAMA/kanonik_uret.py ile uretilir. Tanim: yon.py
KONSENSUS_KANONIK = y('konsensus_kanonik')
KONSENSUS_INDEKS  = os.path.join(KONSENSUS_KANONIK, 'INDEKS.tsv')
KONSENSUS_HAM     = y('consensus sequences')      # YALNIZ kanonik uretimi icin
KONSENSUS         = KONSENSUS_KANONIK
FASTQ       = y('fastq files')
REFDB       = y('REFERANS_DB')
SILVA_SSU   = os.path.join(REFDB, 'SILVA_138.2_SSURef_NR99.fasta')
SILVA_LSU   = os.path.join(REFDB, 'SILVA_138.2_LSURef_NR99.fasta')
UNITE_ITS   = os.path.join(REFDB, 'UNITE_ITS.fasta')
PR2         = os.path.join(REFDB, 'PR2_SSU_taxo_long.fasta')

# mevcut olcum kodunun bulundugu klasorler (ICE AKTARILIR, yeniden yazilmaz)
BETIK_YOLLARI = [y('SON_ETAP_betikleri'), y('MADDE123_betikleri'), y('DUZELTME_betikleri')]

CIKTI       = y('KAPSAMLI_ARAMA_SONUC')
KONTROL     = os.path.join(CIKTI, 'kontrol')      # checkpoint
ONBELLEK    = os.path.join(CIKTI, 'onbellek')     # cache

# ---------------------------------------------------------------- qPCR kisitlari (SABIT)
# QuantiNova SYBR Green + Rotor-Gene Q
URUN_IDEAL      = (60, 150)   # tercih edilen
URUN_KABUL      = (150, 250)  # kabul edilebilir, protokolde 30 sn annealing/extension
URUN_MUTLAK_UST = 400         # aramada uretilen en buyuk urun (400 ustu hic denenmez)
URUN_ONERILMEZ  = 250         # bunun ustu "onerilmez" olarak isaretlenir

TA_HEDEF        = 60.0        # butun panelin ayni Ta'da kosmasi hedefi, 60 C oncelikli
TA_KURALI       = 3.0         # Ta = min(Tm) - 3  (panelin kurali)

# SYBR Green -> yapi olcutleri ELEYICI (uyari degil)
HAIRPIN_TM_UST     = 45.0
HOMODIMER_TM_UST   = 45.0
HETERODIMER_TM_UST = 45.0
# delta-G esikleri (kcal/mol, 60 C'de hesaplanir); daha negatif = daha kotu
DG_YAPI_ALT        = -9.0     # herhangi bir yapi
DG_UC_ALT          = -5.0     # 3' uc kaynakli dimer

# primer3 tuz/deriism kosullari - PANELDEKI DEGERLERIN AYNISI (geo.py ile birebir)
P3 = dict(mv_conc=50, dv_conc=1.5, dntp_conc=0.6, dna_conc=50)

# ---------------------------------------------------------------- degismez primer kurallari
UZUNLUK    = (18, 25)   # kullanicinin istegi: her pozisyondan 18-25 arasi her oligo
TEKRAR_UST = 4          # ayni bazin 4'lu tekrari yasak (ara.py ile ayni)
DTM_UST    = 1.5        # cift ici Tm farki

# ---------------------------------------------------------------- PARAMETRE IZGARASI
# Siki ayardan baslar, kademeli gevser. Sira ONEMLI: rapor "hangi ayarda cozum
# cikti" sorusunu bu siraya gore cevaplar.
IZGARA_GC    = [(40, 60), (37, 63), (35, 65)]
IZGARA_TM    = [(58, 62), (57, 63), (56, 64)]
IZGARA_URUN  = [(60, 150), (60, 200), (60, 250), (60, 400)]
IZGARA_UC_GC = [True, False]    # 3' son baz G/C sart mi
IZGARA_SON5  = [True, False]    # 3' son 5 bazda en cok 3 G/C sart mi
# toplam 3*3*4*2*2 = 144 hucre

# ---------------------------------------------------------------- numune olcumu
NUMUNE_OKUMA_MIN, NUMUNE_OKUMA_MAX = 200, 6000   # duzeltilmis filtre (bkz. 10 Olcum Hatalari #2)
NUMUNE_OKUMA_SAYISI = 300      # kutu basina ornek okuma (sabit tohum)
NUMUNE_TOHUM        = 20260802
NUMUNE_MAX_MM       = 1        # <=1 uyumsuzluk + 3' son 2 baz TAM (panel numune olcutu)
KURESEL_MAX_MM      = 5        # kuresel olcut (panel: toplam <=5, F+R ayri yazilir)
REFERANS_MAX_MM     = 1
# Bu esik MUTLAK olmak zorundadir. Goreli bir esik ("en buyuk kutunun yarisi")
# denendiginde, tam derinlikte en buyuk kutu ~46 000 okuma olunca esik ~23 000'e
# cikiyor ve 10-33 rakip kutudan yalniz 1-5'i olcume giriyordu; yani "en kotu
# rakip kutu" fiilen "en derin kutu" anlamina geliyor ve gercek en kotu rakip
# olcum disinda kaliyordu. Ayni cift 300 okumayla olculunce butun kutular
# giriyor ve iki asama arasinda 40 kata varan sahte fark cikiyordu.
ENKOTU_ASGARI_OKUMA = 150      # "en kotu tek rakip kutu" olcusune girmek icin
                               # gereken ASGARI okuma. MUTLAK olmali: derinlige
                               # gore kayarsa asamalar karsilastirilamaz hale gelir.
# KAPSAM ekseni, ayrim kati tanimsizlastiginda (evrensel hedeflerde rakip kumesi
# bosa yaklasir, payda sifira gider) tek anlamli olcu olarak kalir. Ayrica tek
# bir uye kutusunun bos cikmasi ayrim katini sifirlar; kapsam bundan etkilenmez.
KAPSAM_ESIGI        = 0.20     # bir uye kutusu 'kapsandi' sayilmak icin en az %20 urun

# ---------------------------------------------------------------- huni (funnel) kapasiteleri
# Her asama bir sonrakine kac aday gecirir. Buyutmek suresi uzatir.
HUNI = dict(
    cift_ust        = 400000,   # asama B: sayilan en fazla cift (izgara tablosu bunun uzerinden)
    numuneye_giden  = 1200,     # asama C: ham okuma taramasina giden aday cift
    referansa_giden = 120,      # asama D: referans kapsam / rakip ayrimina giden
    arms_taban      = 25,       # ARMS denenen 'en iyi' aday sayisi
    arms_ust        = 400,      # numunede olculen en fazla ARMS varyanti
    kusele_giden    = 12,       # asama E: kuresel taramaya giden (EN PAHALI)
)

# kuresel taramada bir seferde islenen baz sayisi (bellek tavani ~ bunun 6 kati bayt)
KURESEL_PARCA = 40_000_000


# ---------------------------------------------------------------------------
# AYRIM ESIGI - TEK KAYNAK
#
# 2026-08-06: esik artik KAT olarak degil, dCq (delta Cq) olarak tanimlanir.
# Kullanicinin karari: dCq = 3.
#
# NEDEN dCq: laboratuvarin konustugu dil budur. qPCR'da olculen sey dongu
# farkidir; kat farki ondan TURETILIR. Esigi kat cinsinden gomdugumuzde
# ("10x") sayi hem keyfi goruyordu hem de literaturle karsilastirilamiyordu.
#
# DONUSUM: %100 verimde her dongu urunu IKIYE katlar, yani
#     kat = 2 ** dCq        ve      dCq = log2(kat)
# dCq 3 -> 2**3 = 8,00 kat. (Eski arac esigi 10x = dCq 3,32 idi.)
#
# KAYNAK: dCq >= 3, ozgulluk / NTC gecme olcutu olarak literaturde kabul
# gormus tabandir (NEB yuksek verimli qPCR veri analizi kilavuzu ve MIQE
# raporlama dili). Bu sayi ARTIK BIZIM UYDURDUGUMUZ BIR ARAC ESIGI DEGILDIR;
# ciktilarda "literatur olcutu" olarak isaretlenir.
#
# VERIM UYARISI - SAYIYI OKURKEN AKILDA TUTULACAK:
# Yukaridaki donusum verimi %100 varsayar. Gercek verim daha dusukse AYNI dCq
# daha KUCUK bir kat farkina denk gelir:
#     kat = (1 + E) ** dCq        E = verim (0-1)
#     %100 verim -> 3 dongu = 8,00 kat
#     % 90 verim -> 3 dongu = 6,86 kat
#     % 80 verim -> 3 dongu = 5,83 kat
# Laboratuvarda kalibrasyon egrisi cikarilinca gercek verim yerine konmali ve
# esigi gecen satirlar YENIDEN degerlendirilmelidir. Bu not butun ciktilara
# basilir (ESIK_VERIM_NOTU).
# ---------------------------------------------------------------------------
ESIK_DCQ = 3.0                      # <-- DEGISTIRILECEK TEK YER
AYRIM_ESIK = 2.0 ** ESIK_DCQ        # kat karsiligi: 8,00

ESIK_KOKENI = (u'LITERATUR OLCUTU dCq >= %.1f (ozgulluk/NTC gecme tabani; '
               u'NEB qPCR veri analizi, MIQE raporlama dili)' % ESIK_DCQ)

ESIK_VERIM_NOTU = (
    u'VERIM %%100 VARSAYILDI. dCq %.1f = %.2f kat donusumu her dongunun urunu '
    u'ikiye katladigini kabul eder. Gercek verim daha dusukse ayni dCq daha '
    u'kucuk kat farkina denk gelir: %%90 verimde %.2f kat, %%80 verimde %.2f '
    u'kat. Kalibrasyon egrisi cikinca gercek verimle yeniden degerlendirin.'
    % (ESIK_DCQ, 2.0 ** ESIK_DCQ, 1.9 ** ESIK_DCQ, 1.8 ** ESIK_DCQ))


def esik_metni(kat=None):
    """Butun ciktilarda AYNI bicim: 'dCq 3,0 (8,00x)'."""
    k = AYRIM_ESIK if kat is None else kat
    return u'dCq %s (%sx)' % (('%.1f' % ESIK_DCQ).replace('.', ','),
                              ('%.2f' % k).replace('.', ','))


def kat_dcq(kat, verim=1.0):
    """Kat -> dCq. Verim 1.0 = %100. Gercek verim olculunce burasi kullanilir."""
    import math as _m
    try:
        k = float(kat)
    except (TypeError, ValueError):
        return None
    if k <= 0:
        return None
    return round(_m.log(k) / _m.log(1.0 + verim), 2)


def kat_ve_dcq(kat):
    """'8,45x (dCq 3,08)' - kat ve dCq YAN YANA. Olcum satirlarinda kullanilir."""
    d = kat_dcq(kat)
    if d is None:
        return u'-'
    return u'%sx (dCq %s)' % (('%.2f' % float(kat)).replace('.', ','),
                              ('%.2f' % d).replace('.', ','))
