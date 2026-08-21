# -*- coding: utf-8 -*-
"""KIMLIK DOGRULAMA TURU - hocaya gonderilen kimlik iddialarini BAGIMSIZ olarak sinar.

YONTEM NEDEN FARKLI
-------------------
Bu tur, iddialari ureten yontemlerin HICBIRINI tekrar etmez:

  * Kraken2      : k-mer + taksonomi agacinda en kucuk ortak ata (LCA).
                   Bu betik taksonomi agacini HIC kullanmaz.
  * 1. turumuz   : sinif ici konsensus hizalamasi + ayirt edici 21-mer.
                   Bu betik konsensusleri BIRBIRIYLE degil, DIS REFERANS
                   veritabanlarindaki adlandirilmis kayitlarla kiyaslar.
  * 2. turumuz   : in-silico PCR (primer tabanli).
                   Bu betikte primer yok.

Bu turun yontemi: TOHUM + HIZALAMA (BLAST mantigi, taksonomisiz).
  1) Sorgu konsensustan k-mer tohumlari cikarilir, veritabani akis halinde
     taranir ve tohum sayisina gore kisa liste yapilir.
  2) Kisa listedeki her kayitla TAM HIZALAMA yapilir (Levenshtein DP, infix).
  3) Kimlik IKI KEZ olculur:
       - tam ortusme uzerinden
       - AYIRT EDICI PENCERE uzerinden: en iyi N referans kaydin BIRBIRINDEN
         AYRILDIGI kolonlar. Korunmus bolgeler (18S, 5.8S, LSU cekirdegi)
         boylece disarida kalir - iddia tam da o korunmus bolgelerden gelen
         sahte yuksek kimlige dayaniyorsa burada gorunur.
  4) BIRDEN FAZLA veritabaninin UYUSMASI SARTTIR. Tek veritabanindan gelen
     sonuc "DOGRULANDI" sayilmaz, en fazla "DOGRULANAMADI (tek kaynak)" olur.

HUKUM
  DOGRULANDI      : >=2 bagimsiz veritabani iddiayi destekliyor
  DUZELTILMELI    : >=2 bagimsiz veritabani BASKA bir sonucta birlesiyor
                    (dogru ifade yazilir - kullanici hocasina duzeltme atacak)
  DOGRULANAMADI   : kanit yetersiz, celiskili ya da tek kaynakli

Uydurma teyit URETILMEZ. Kanit yetersizse DOGRULANAMADI yazilir.
Panel dosyalarina YAZMAZ; KIMLIK_SONUC/ altina yazar.
"""

# -------------------------------------------------------------------------
# kimlik_dogrulama.py — hocaya gonderilen kimlik iddialarini dis referans
# veritabanlarina karsi TOHUM + HIZALAMA ile bagimsiz olarak sinar.
#
# GİRDİ  : REFERANS_DB/ altindaki yerel FASTA kumeleri (VTB listesi; ikiz ve
#          altkume olanlar oylamadan cikarilmistir),
#          konsensus_kanonik/ (KAPSAMLI_ARAMA.hedefler.konsensusler ile),
#          NCBI nt (ag uzerinden ayri katman) ve elle doldurulmus
#          KIMLIK_SONUC/nt_elle/NT_SONUC_SABLONU.tsv.
# ÇIKTI  : KIMLIK_SONUC/kimlik_iddialari.tsv (asil tablo),
#          KIMLIK_SONUC/KIMLIK_DOGRULAMA_RAPORU.md,
#          KIMLIK_SONUC/VERITABANI_ENVANTERI.md,
#          KIMLIK_SONUC/LITERATUR_ELLE_KONTROL.tsv,
#          KIMLIK_SONUC/nt_ham/, nt_elle/, kontrol/ .
#          Panel dosyalarina YAZMAZ.
# ÇAĞRAN : KAPSAMLI_ARAMA.bat -> I tusu
#          (bat icinde: wsl -e python3 "KURTARMA/kimlik_dogrulama.py" --kok .)
#          Bu dosyanin fonksiyonlari ayrica G asamasi (tum_kutu_kimlikleri.py) ve
#          E asamasi (erisim_dogrulama.py) tarafindan modul olarak cagrilir.
#
# UC KURAL, hepsi asagida kodda isaretli:
#   1) KISA LISTE 500'DUR ve HEPSI hizalanir - kesme noktasi baglayici degildir.
#   2) Bir iddia DOGRULANDI sayilmak icin EN AZ IKI bagimsiz veritabani uyusmali.
#   3) Kimlik ayrica AYIRT EDICI PENCEREDE olculur - korunmus bolgeler disarida.
# -------------------------------------------------------------------------
import os, sys, csv, json, time, re, argparse

VERSIYON = '1.1 (2026-08-04)'

# --------------------------------------------------------------- AYAR IMZASI
# Kontrol noktasi (checkpoint) anahtarina KARISTIRILIR. Kisa liste boyu ya da
# siralama/secim mantigi degistiginde bu satiri artirin: eski kosulardan kalmis
# vtb_*.json dosyalari boylece OTOMATIK gecersiz olur ve sessizce eski (kucuk
# listeyle uretilmis) sonucu geri vermezler.
AYAR_IMZASI = 'idfbm25-kl500-v3'

K_TOHUM = 16            # tohum uzunlugu

# --------------------------------------------------- SIRALAMA OLCUTU (2026-08-05)
# ESKI OLCUT: skor = "sorgunun kac ayri tohumu bu kayitta geciyor" (duz sayim).
# Iki yonden birden bozuluyordu ve ikisi de veritabani buyudukce kotulesiyordu:
#
#   (a) KORUNMUS BOLGE GURULTUSU. 18S/5.8S/LSU cekirdeginden gelen tohumlar
#       veritabaninin neredeyse tamaminda geciyor (olculdu: bir tohum 2 069 188
#       kayitin 1 612 663'unde). Bu tohumlar hicbir sey ayirt etmiyor ama duz
#       sayimda gercek akrabalik kaniti olan nadir tohumla AYNI agirligi aliyor.
#   (b) UZUNLUK YANLILIGI - ASIL SEBEP BUYDU. Skor bir TOPLAM oldugu icin uzun
#       kayit her zaman kazaniyor. Hedef kayit AY882347 UNITE'de 484 bp; sorgunun
#       ancak 19 tohumunu tutabiliyor. Onun yerine ilk siralara oturan kayitlar
#       892-2000 bp uzunlugunda ve 34-59 tohum tutuyor. Kisa ve TAM eslesen kayit,
#       yalnizca KISA oldugu icin gomuluyor.
#
# YENI OLCUT = ters frekans agirligi + BM25 uzunluk normalizasyonu:
#
#       skor(kayit) = SUM_{tutan tohumlar} ln(N / (1 + df(tohum)))
#                     ----------------------------------------------
#                          1 - b + b * uzunluk(kayit) / ORT_UZUNLUK
#
#   df(tohum) : o tohumun KAC KAYITTA gectigi (ayni akista bedavaya sayilir)
#   N         : taranan kayit sayisi,  b = BM25_B = 0,75,  ORT = ortalama uzunluk
#
# OLCULDU (AY882347 kaydinin sirasi, hepsi TAM veritabani taramasi):
#   UNITE ITS (2 069 188 kayit)   F2-1_101201    1869 -> 19      (yalniz idf: 1320)
#                                 F2-1_2034170   2037 -> 18
#                                 F1-4_101201     136 -> 17
#                                 F1-4_2093779    162 -> 18
#                                 F1-4_2093780    153 -> 17
#                                 F1-2_101201     189 -> 30
#                                 F1-1_2093779    225 -> 35
#   SILVA LSU Parc (1 312 521)    F2-1_101201    6794 ->  2
#                                 F2-1_2034170   5749 ->  3
#   RefSeq ref_all2 (65 358)      F2-1_101201    1320 ->  1
#   RefSeq mantar ITS (20 394)    yedi sorguda da 1 -> 1  (BOZULMA YOK)
#
# Yalniz ters frekans YETMEDI (1869 -> 1320). Uzunluk normalizasyonu sart;
# ikisi birlikte calisiyor. Ayrinti: SIRALAMA_COZUMU.md.
BM25_B = 0.75           # BM25 uzunluk normalizasyon katsayisi (0 = normalizasyon yok)
ADAY_HAVUZU = 3000      # akis sirasinda tutulan aday sayisi; idf ile yeniden
                        # siralanip ilk KISA_LISTE tanesi hizalamaya gider.
                        # Olculen en kotu on-siralama sirasi 45 -> 66x pay var.

# KISA LISTE BOYU - 2026-08-04, ikinci duzeltme.
#
# SORUN (olculdu): kisa liste TOHUM SAYISINA gore siralanip kesiliyor, ama
# karari HIZALAMA KIMLIGI veriyor. Iki olcut farkli oldugu icin gercek en iyi
# eslesme listeye hic giremeyebiliyor:
#   * F2-1_101201: liste Parascedosporium putredinis'i (%98,16) en iyi gosterdi;
#     ayni veritabaninda Petriella guttulata %99,65 vardi.
#   * SILVA LSU NR99, %8 bozulmus sorgu: hedef kaydin 47 tohumu SAG oldugu halde
#     ilk 60'a giremedi - ayni korunmus bolgeleri paylasan binlerce kayitla
#     yarisiyordu.
# "Beklenen taksonu garanti et" yamasi bu hatayi ancak NE ARADIGIMIZI BILDIGIMIZDE
# kapatir; kimlik sorusunda bilmiyoruz.
#
# COZUM: kesme noktasini BAGLAYICI OLMAKTAN CIKAR. Liste 60 -> 500 buyutuldu ve
# 500 adayin HEPSI hizalaniyor; karar tamamen hizalamaya birakildi.
# Maliyet: vektorlestirilmis hizalayici 1,5 kb icin ~0,02 sn -> 500 aday ~10 sn
# (kisa liste asamasinin veritabani akisi zaten dakikalar suruyor, yani hizalama
# maliyeti toplamin kucuk bir parcasi).
#
# OZ-KALIBRASYON: her sorguda KAZANAN isabetin kisa listedeki tohum sirasi
# kaydedilir (bkz. kazanan_sira). Kazananlar hep ilk 100'den geliyorsa kesme
# noktasi baglayici degildir ve bunu AYRI BIR OLCUM YAPMADAN kanitlariz.
# 2026-08-05 (sabah): 500 -> 1000. Oz-kalibrasyon 500'un BAGLAYICI oldugunu
# gosterdi (118 sorgunun 13'unde kazanan 400. siranin otesinden geldi; en yuksek
# kazanan sira 4171). 1000 bunun bir kismini kurtardi ama 4171'i KURTARMADI.
#
# 2026-08-05 (aksam): 1000 -> 500 GERI. Listeyi buyutmek sorunun kendisini
# cozmuyordu, yalnizca daha derin kaziyordu. Sorun SIRALAMA OLCUTUNDEYDI (yukari
# bakin: ters frekans + uzunluk normalizasyonu). Yeni olcutle olculen en kotu
# kazanan sirasi 35'tir; 500 bunun 14 katidir ve I/G asamalarinin suresini
# yariya indirir. Buyutmek artik gereksiz.
KISA_LISTE = 500        # her veritabanindan tam hizalanacak kayit sayisi
SIRA_UYARI_ESIGI = 200  # kazanan bunun otesinden geldiyse olcut yine bozuluyor demektir
SIRA_GUVENLI_BOLGE = 50   # kazananlar bunun icinde kaliyorsa kesme baglayici degil
GARANTI_UST = 40        # "beklenen takson garantisi" ile alinacak azami kayit
AYIRT_EDICI_UST = 8     # ayirt edici pencere icin kullanilacak en iyi kayit sayisi
UYUM_TOLERANS = 1.0     # iddia edilen yuzde ile olculen arasindaki kabul edilir fark

# (etiket, dosya, tur, kimlik_asamasinda_kullan, not)
# KIMLIK asamasi icin kural: TEKRARSIZLASTIRILMIS kume YETMEZ. NR99 gibi kumeler
# nadir cinsleri tumden siler - olculdu: SILVA LSURef NR99'da Petriella kaydi
# SIFIR, ayni surumun Parc kumesinde 82. Kimlik sorusunda Parc SART.
VTB = [
    ('SILVA SSU NR99',     'SILVA_138.2_SSURef_NR99.fasta', 'SSU',    True,
     u'510 495 kayit; SSU tekrarsizlastirilmis'),
    ('SILVA LSU NR99',     'SILVA_138.2_LSURef_NR99.fasta', 'LSU',    True,
     u'95 279 kayit; LSU tekrarsizlastirilmis - NADIR CINSLERI SILER'),
    ('SILVA LSU Parc',     'SILVA_138.2_LSUParc.fasta',     'LSU',    True,
     u'1 312 521 kayit; TEKRARSIZLASTIRILMAMIS. Petriella: NR99=0, Parc=82 '
     u'(olculdu). Kimlik asamasinda SART.'),
    ('UNITE ITS',          'UNITE_ITS.fasta',               'ITS',    True,
     u'2 069 189 kayit; mantar ITS. Petriella: 113 kayit (olculdu).'),
    ('PR2 SSU',            'PR2_SSU_taxo_long.fasta',       'SSU',    True,
     u'240 201 kayit; okaryot 18S'),
    ('ROD operon',         'ROD_v1.2_operon_variants.fasta','OPERON', True,
     u'60 320 kayit; rRNA operon varyantlari'),
    ('RefSeq bakteri 16S', 'bacteria.16S.fna',              'SSU',    True,
     u'26 877 kayit; adlandirilmis tip materyali agirlikli'),
    ('RefSeq arke 16S',    'archaea.16S.fna',               'SSU',    True,
     u'1 160 kayit'),
    ('RefSeq mantar ITS',  'fungi.ITS.fna',                 'ITS',    True,
     u'20 394 kayit'),
    ('RefSeq mantar 28S',  'fungi.28SrRNA.fna',             'LSU',    True,
     u'12 890 kayit; Petriella: 2 kayit (olculdu)'),
    ('RefSeq mantar 18S',  'fungi.18SrRNA.fna',             'SSU',    True,
     u'4 037 kayit'),
    ('RefSeq ref_all2',    'ref_all2.fna',                  'KARISIK', True,
     u'65 358 kayit; RefSeq birlesik kume (ref_all\'in ustkumesi)'),
    ('RefSeq ref_all',     'ref_all.fna',                   'KARISIK', False,
     u'48 431 kayit; ref_all2 bunun USTKUMESI - ayni kayitlar iki kez '
     u'sayilmasin diye kimlik oylamasindan cikarildi (bagimsiz kaynak degil)'),
]

# NCBI nt AYRI bir katmandir (yerel dosya degil, ag uzerinden).
NT_ETIKET = 'NCBI nt'

# --------------------------------------------------------------- IDDIALAR
# tip: 'kimlik'  -> kutu ile adlandirilmis takson arasindaki kimlik iddiasi
#      'kutu2'   -> iki kutunun birbirine kimligi
#      'dagilim' -> bir organizmanin birden cok kutuya dagildigi iddiasi
#      'adsiz'   -> "adlandirilamiyor / ad verilemez" iddiasi
#      'ayrilmaz'-> iki turun birbirinden ayrilamadigi iddiasi
#      'gecici'  -> "bu sinifta kimlik ayrimi yapilamiyor" iddiasi
IDDIALAR = [
 dict(no=1, oncelik=2, tip='kimlik', kutu=['F2-1_101201'], sinif='F2',
      beklenen_cins='Petriella', beklenen_yuzde=None,
      metin=u'taxid 101201 kutusu Petriella cinsinden'),
 dict(no=2, oncelik=1, tip='kutu2', kutu=['F1-4_2093780'], karsi=['F2-1_101201'],
      beklenen_yuzde=99.58, beklenen_cins='Petriella',
      metin=u'Podospora pseudopauciseta kutusu (F1-4_2093780) Petriella ile %99,58 ayni',
      not_=u'EN YENI VE EN AZ SINANMIS IDDIA - oncelik yuksek'),
 dict(no=3, oncelik=1, tip='dagilim', beklenen_cins='Petriella',
      kutu=['F2-1_101201','F2-2_101201','F2-3_101201','F2-4_101201',
            'F2-1_2034170','F2-4_2034170','F2-1_500148','F2-2_500148','F2-4_500148'],
      supheli=['F2-1_2034170','F2-4_2034170','F2-1_500148','F2-2_500148','F2-4_500148'],
      beklenen_alt=76.0, beklenen_ust=86.0,
      # 2026-08-06 DUZELTMESI: metin "dokuz kutu" diyordu ve bu sayi hicbir
      # olcumden gelmiyordu - iddianin kendi kutu listesinin uzunluguydu.
      # G asamasi 96 kutunun HEPSINI bagimsiz taradi ve sunu olctu:
      #   * bu iddianin dokuz kutusundan SEKIZI Petriella cikti,
      #   * dokuzuncusu (F2-4_500148) CIKMADI - adlandirilamayan soy, en yakin
      #     kayit Lomentospora prolifica %83,68; iddianin kendi kanit satirinda
      #     da zaten %52,21 ile ayrisiyordu,
      #   * buna karsilik iddianin saymadigi DORT F1 kutusu (F1-2_101201,
      #     F1-4_101201, F1-4_2093779, F1-4_2093780) da Petriella cikti.
      # Toplam 12 kutu. Metin artik olculen sayiyi yaziyor; hocaya giden
      # belgede "dokuz" ile kanit arasindaki celiski kalmasin.
      metin=u'Petriella numunede 12 kutuya dagilmis (bu iddianin dokuz '
            u'kutusundan sekizi + dort F1 kutusu; F2-4_500148 Petriella '
            u'DEGIL); Kraken\'in T. breve ve M. brunneum dedigi kutularin '
            u'okumalarinin %76-86\'si Petriella',
      not_=u'tek tura dayaniyor - oncelikli'),
 dict(no=4, oncelik=3, tip='kimlik', kutu=['A1-4_2208'], sinif='A1',
      beklenen_tur='Methanosarcina vacuolata', beklenen_yuzde=97.6,
      metin=u'M. barkeri kutusu M. vacuolata\'ya %97,6'),
 dict(no=5, oncelik=3, tip='kimlik', kutu=['B-2_818','B-3_818','B-2_214856'], sinif='B',
      beklenen_tur='Alistipes putredinis', beklenen_yuzde=85.0, tolerans=2.0,
      adsiz_bekleniyor=True,
      metin=u'Bacteroides kutulari adlandirilamayan Bacteroidales, en yakin '
            u'Alistipes putredinis %85 civari'),
 dict(no=6, oncelik=3, tip='kimlik', kutu=['B-2_1197717'], sinif='B',
      beklenen_yuzde=99.4, beklenen_aile='Synergistaceae',
      rakip_cins='Cloacibacillus', rakip_yuzde=90.0, adsiz_bekleniyor=True,
      metin=u'Cloacibacillus hedefi adlandirilamayan Synergistaceae %99,4, '
            u'en yakin Cloacibacillus %90'),
 dict(no=7, oncelik=2, tip='adsiz', kutu=['F1-1_44689','F1-2_44689','F1-3_44689','F1-4_44689'],
      metin=u'taxid 44689 etiketi curutuldu ve yerine isim konulamiyor'),
 dict(no=8, oncelik=2, tip='ayrilmaz', kutu=['A1-2_2209'], sinif='A1',
      turler=[('Methanosarcina soligelidi', 99.93), ('Methanosarcina mazei', 99.85)],
      metin=u'taxid 2209 kutusu M. soligelidi %99,93 / M. mazei %99,85, ikisi ayrilmiyor'),
 dict(no=9, oncelik=2, tip='cins_duzeyi', kutu=['F2-1_101201'], beklenen_cins='Petriella',
      metin=u'Petriella cins duzeyinde kalmali, tur adi verilemez, cf. setifera'),
 # --- K asamasinin "HETEROJEN" dedigi hedefler: karar verebilmek icin
 # kutu kimliklerinin DIS REFERANSLA sinanmasi sart. K bunu yapamaz, I yapar.
 dict(no=11, oncelik=1, tip='kimlik', kutu=['B-4_285070'], sinif='B',
      beklenen_cins='Petrimonas',
      metin=u'Petrimonas hedefinin uye kutulari ayni organizma mi? '
            u'(K asamasi HETEROJEN dedi: 3 kutunun hicbiri birbiriyle >=%99 kimlikte)',
      not_=u'K asamasi bu hedef icin daraltma UYGULAYAMADI ve "once kutu kimlikleri '
            u'referansla dogrulanmali" dedi. Bu iddia tam olarak o denetimdir; '
            u'sonucu olmadan Petrimonas hakkinda karar verilemez.'),
 dict(no=12, oncelik=1, tip='kimlik', kutu=['B-2_818'], sinif='B',
      beklenen_tur='Alistipes putredinis', beklenen_yuzde=85.0, tolerans=3.0,
      adsiz_bekleniyor=True,
      metin=u'Bacteroidales_kumesi hedefinin uye kutulari ayni organizma mi? '
            u'(K asamasi HETEROJEN dedi: 12 kutunun hicbiri birbiriyle >=%99 kimlikte)',
      not_=u'Ayni gerekce: K daraltma uygulayamadi, dis referans teyidi sart.'),
 dict(no=10, oncelik=2, tip='gecici', sinif='B',
      metin=u'B sinifi kutu kimlikleri gecici, konsensus kimligi o sinifta '
            u'ayrim yapamiyor'),
]


# --------------------------------------------------------------- temel
_C = str.maketrans('ACGTUNacgtun', 'TGCAANtgcaan')


# Ters tumleyen. Sorgu her iki yonde de aranir: referans kayitlarin yonu kume
# kume degisir, tek yon aranirsa yarisi kacar.
def rc(s):
    return s.translate(_C)[::-1]


# ACGT disindaki her sey N'e cevrilir (U -> T). N'ler hizalamada daima uyumsuz
# sayilir, yani belirsiz baz asla lehte yorumlanmaz.
def temizle(s):
    return re.sub(r'[^ACGT]', 'N', s.upper().replace('U', 'T'))


def sure_metni(sn):
    sn = int(sn)
    return ('%d saniye' % sn) if sn < 90 else ('%d dakika' % round(sn / 60.0)) \
        if sn < 5400 else ('%.1f saat' % (sn / 3600.0))


# Turkce ondalik ayraci. None ve sayiya cevrilemeyen degerler "-" olur; rapor
# icinde "0" ile "olculmedi" birbirine karismasin diye.
def vir(x, b=2):
    if x is None:
        return '-'
    try:
        return ('%.*f' % (b, float(x))).replace('.', ',')
    except (TypeError, ValueError):
        return str(x)


_M = None


def enc(s):
    global _M
    import numpy as np
    if _M is None:
        m = np.full(256, 4, dtype=np.uint8)
        for i, c in enumerate('ACGT'):
            m[ord(c)] = i
        _M = m
    return _M[np.frombuffer(s.encode(), dtype=np.uint8)]


# -------------------------------------------------------------------------
# INFIX (HW) LEVENSHTEIN - KARARI VEREN OLCU.
#
# Kisa sorgu uzun hedefin ICINE hizalanir, yani hedefin basinda ve sonunda kalan
# fazlalik cezalandirilmaz. Referans kayitlar cok farkli uzunluklarda gelir (kimi
# tam operon, kimi tek lokus); global hizalama bu uzunluk farkini uyumsuzluk gibi
# sayar ve dogru kaydi eler.
#
# BU FONKSIYON KARARI VERIR, tohum sayisi vermez. Tohum siralamasi yalnizca ADAY
# TOPLAMA olcutudur; kisa listedeki 500 adayin HEPSI buradan gecirilir ve en iyi
# isabet hizalama kimligine gore secilir. Iki olcutun ayri olmasi 60'lik liste
# doneminde gercek en iyi eslesmelerin kacirilmasina yol acmisti.
# -------------------------------------------------------------------------
def hizala(q, t):
    """Infix (HW) Levenshtein: kisa sorguyu uzun hedefin ICINE hizalar.
    Donen: (yuzde_kimlik, uzaklik). numpy ile satir satir DP."""
    import numpy as np
    if not q or not t:
        return (0.0, len(q or t or ' '))
    Q, T = enc(q), enc(t)
    onceki = np.zeros(len(T) + 1, dtype=np.int32)
    for i in range(len(Q)):
        simdi = np.empty_like(onceki)
        simdi[0] = i + 1
        farkli = (T != Q[i]) | (T == 4) | (Q[i] == 4)
        # Sol komsu bagimliligi (ekleme) VEKTORLESTIRILDI:
        #   simdi[j] = min(aday[j], simdi[j-1]+1)
        # a[j] = simdi[j]-j konursa a[j] = min(aday[j]-j, a[j-1]) olur, yani
        # kosan minimum. np.minimum.accumulate ile tek gecis. Python ic dongusu
        # kaldirildi: 1,5 kb x 1,5 kb hizalama dakikalardan saniyelere iner.
        aday = np.minimum(onceki[:-1] + farkli, onceki[1:] + 1)
        aday = np.concatenate(([i + 1], aday))
        idx = np.arange(len(aday))
        simdi = np.minimum.accumulate(aday - idx) + idx
        onceki = simdi
    d = int(onceki.min())
    return (round(100.0 * (1 - d / float(len(q))), 2), d)


def fasta_akisi(yol):
    """(baslik, dizi) uretir. Bellege sigmayan dosyalar icin akis.

    TAVAN YOKTUR - dosyanin SONUNA kadar okur. (Eskiden kullanilmayan bir
    'parca=200000' parametresi vardi; hicbir yerde ise yaramiyordu ama tavan
    varmis izlenimi veriyordu, o yuzden kaldirildi. Erisim testinde gercek bir
    tavan sorunu yasanmisti: ilk kosu 120 001 kayitta kesiyordu ve SILVA SSU
    NR99 -510 495-, LSU Parc -1 312 521-, UNITE ITS -2 069 189- fiilen budanmis
    halde taraniyordu. Kimlik asamasinda o hata YOK; taranan kayit sayisi her
    satirda basilir ve beklenen kayit sayisiyla karsilastirilir.)"""
    bas, par = None, []
    with open(yol, encoding='utf-8', errors='ignore') as fh:
        for satir in fh:
            if satir.startswith('>'):
                if bas is not None:
                    yield bas, temizle(''.join(par))
                bas = satir[1:].strip(); par = []
            else:
                par.append(satir.strip())
    if bas is not None:
        yield bas, temizle(''.join(par))


# TOHUM ADIMI - olculerek secildi (2026-08-04).
# Kisa liste asamasi butun veritabani akisini tarar ve maliyeti TOHUM SAYISI ile
# dogru orantilidir. Olculdu (SILVA LSU Parc, 2,1 GB / 1,31 M kayit):
#   adim= 7 -> 410 tohum -> ~21 dakika
#   adim=25 -> 116 tohum -> ~ 6 dakika   <- SECILEN
#   adim=60 ->  48 tohum -> ~ 3 dakika
# Kisa listenin bozulup bozulmadigi sinandi: adim=25'te ilk 25'in 23/25 ve 17/25'i
# ayni kaldi ve HER IKI SINAMADA DA EN IYI ISABET AYNI CIKTI. Butun hukumler
# isabet[0]'a dayandigi icin adim=25 guvenli; adim=60 daha cok kayip veriyor.
def tohumlar(q, k=K_TOHUM, adim=25):
    return {q[i:i + k] for i in range(0, max(1, len(q) - k + 1), adim) if 'N' not in q[i:i + k]}


# -------------------------------------------------------------------------
# ORTALAMA KAYIT UZUNLUGU - BM25 normalizasyonunun paydasi
#
# Normalizasyon L/ORT'ye dayandigi icin ORT AKIS BASLAMADAN bilinmelidir:
# yoksa dosyanin ilk kayitlari yanlis puanlanir ve aday havuzundan yanlis elenir
# (kosan ortalama ile denendi, ilk %5'te sira oynuyor). Bir kez hesaplanir ve
# veritabaninin yanina .ortuz dosyasina yazilir; sonraki kosular okur.
# -------------------------------------------------------------------------
def ortalama_uzunluk(yol):
    yan = yol + '.ortuz'
    try:
        if os.path.exists(yan) and os.path.getmtime(yan) >= os.path.getmtime(yol):
            n, t = open(yan).read().split()
            if int(n) > 0:
                return float(t) / int(n)
    except Exception:
        pass
    n = 0
    t = 0
    fai = yol + '.fai'
    if os.path.exists(fai):                  # samtools indeksi varsa bedava
        # A6 DUZELTMESI (2026-08-21): bozuk satir eskiden sessizce atlaniyor ama
        # t/n birikimi suruyordu, yani ortalama KISMI veriden hesaplaniyordu.
        # Bu deger kisa liste puanlamasina girer. Bozuk satir varsa .fai'ye hic
        # guvenilmez: indeks tutarsizsa TAMAMI atilir ve asagidaki tam tarama
        # kolu devreye girer. Sessiz kismi olcum yerine gorunur tam olcum.
        fai_bozuk = 0
        for satir in open(fai, encoding='utf-8', errors='ignore'):
            p = satir.split('\t')
            if len(p) > 1:
                try:
                    t += int(p[1]); n += 1
                except ValueError:
                    fai_bozuk += 1
        if fai_bozuk:
            sys.stderr.write(
                'UYARI: %s icinde %d bozuk satir var; indeks yok sayildi ve '
                'ortalama uzunluk FASTA taranarak olculuyor.\n' % (fai, fai_bozuk))
            n = 0
            t = 0
    if not n:                                # yoksa yalniz uzunluk sayan hizli tarama
        u = 0
        with open(yol, 'rb') as f:
            for ham in f:
                if ham[:1] == b'>':
                    if u:
                        t += u; n += 1
                    u = 0
                else:
                    u += len(ham.strip())
        if u:
            t += u; n += 1
    ort = float(t) / max(1, n)
    try:
        open(yan, 'w').write('%d %d' % (n, t))
    except Exception:
        pass
    return ort if ort > 0 else 1.0


class _TersBaslik(object):
    """Baslik ARTAN sirali olsun diye min-heap'te ters karsilastirilan sarmalayici."""
    __slots__ = ('s',)

    def __init__(self, s):
        self.s = s

    def __lt__(self, o):
        return self.s > o.s

    def __eq__(self, o):
        return self.s == o.s

    def __str__(self):
        return self.s


# -------------------------------------------------------------------------
# KISA LISTE - IKI ASAMALI SECIM
#
# Liste eskiden 60'ti, sonra 500, sonra 1000 yapildi. Buyutmek sorunu cozmedi:
# olculen bir vakada gercek akraba (AY882347) UNITE ITS'te 1869. siradaydi ve
# hicbir makul liste boyu onu yakalayamazdi. Sorun liste BOYUNDA degil,
# SIRALAMA OLCUTUNDEYDI (gerekce icin dosyanin basindaki "SIRALAMA OLCUTU"
# bolumune bakin: korunmus bolge gurultusu + uzunluk yanliligi).
#
# Simdi tek akista iki asama var:
#   1) ON ELEME - uzunluga gore normalize edilmis HAM tohum sayisi ile en iyi
#      ADAY_HAVUZU (3000) kayit tutulur. Bu asamada df henuz bilinmez; ama ayni
#      akista her tohumun kac kayitta gectigi BEDAVAYA sayilir.
#   2) YENIDEN SIRALAMA - akis bitince gercek df'lerden idf hesaplanir, 3000
#      aday idf+BM25 skoruyla yeniden siralanir, ilk 'ust' tanesi hizalamaya
#      gider. Olculen en kotu on-eleme sirasi 45'tir; havuz 3000 oldugu icin
#      on-eleme kesmesi baglayici DEGILDIR.
#
# "Beklenen takson garantisi" (garanti parametresi) yamasi yerinde kalir ama
# ARTIK GEREKMIYOR: olculen yedi sorguda da kazanan tohumla listeye giriyor.
# Yama ile giren bir kaydin normal siralamada kacinci olacagi yine hesaplanir;
# kazanan yamayla geldiyse rapor UYARI basar.
# -------------------------------------------------------------------------
def kisa_liste(yol, q, ust=KISA_LISTE, ilerle=None, garanti=(), havuz=None,
               suzgec=None):
    """En umutlu kayitlari sec (BLAST'in tohumlama adimi), idf + BM25 ile.

    Donen: [dict(tohum, skor, baslik, dizi, sira, kaynak), ...]
      tohum  : kac ayri sorgu tohumu bu kayitta gecti (ham sayi, bilgi amacli)
      skor   : idf agirlikli, uzunluga gore normalize edilmis siralama skoru
      sira   : SON siralamadaki yer (1 = en iyi). Kazananin bu sayisi kaydedilir.
      kaynak : 'tohum' -> olcutun icinden geldi | 'garanti' -> takson yamasiyla
    """
    import heapq, math
    th_l = sorted(tohumlar(q) | tohumlar(rc(q)))
    if not th_l:
        return []
    nt = len(th_l)
    ORT = ortalama_uzunluk(yol)
    # ust=0 ESKI ANLAMINI KORUR: kesme yok, tohum tutan HER kayit hizalanir.
    # O durumda on eleme havuzu da sinirsiz olmali, yoksa "kesme yok" sessizce
    # "3000'de kes"e donerdi.
    havuz = None if not ust else max(int(havuz or ADAY_HAVUZU), ust * 3, 500)
    df = [0] * nt
    N = 0
    yigin = []                     # min-heap: (on_skor, _TersBaslik, sirano, ...)
    zorunlu = []
    gar = [g.lower() for g in (garanti or ()) if g]
    n = 0
    kisa_liste.son_taranan = 0
    for bas, diz in fasta_akisi(yol):
        n += 1
        if ilerle and n % 20000 == 0:
            ilerle(n)
        L = len(diz)
        if L < 100:
            continue
        # SUZGEC (istege bagli, varsayilan YOK - eski davranis birebir korunur).
        # 2026-08-11: kisa listenin 500 slotunun neredeyse tamami adlandirilmamis
        # cevre klonuyla doluyordu ve "en yakin adlandirilmis tur" sutunu bos
        # kaliyordu. Suzgec verilirse liste yalniz ADLI kayitlarla dolar; sure
        # ayni kalir cunku tarama zaten butun dosyayi okuyor, degisen sey
        # hangi kayitlarin yigina alindigi.
        if suzgec is not None and not suzgec(bas):
            continue
        N += 1
        # K-6 DUZELTMESI (2026-08-03): eskiden sayac 3'te duruyordu ve skor
        # herkeste doyuyordu. Gercek tohum kumesi kullaniliyor.
        tut = frozenset(i for i in range(nt) if th_l[i] in diz)
        if not tut:
            continue
        for i in tut:
            df[i] += 1             # TERS FREKANS: ayni akista bedavaya sayilir
        norm = 1.0 - BM25_B + BM25_B * L / ORT
        if gar and any(g in bas.lower() for g in gar):
            zorunlu.append((tut, bas, diz, norm))
            continue
        on = len(tut) / norm       # ON ELEME: uzunluga gore normalize ham sayim
        if havuz is None or len(yigin) < havuz:
            heapq.heappush(yigin, (on, _TersBaslik(bas), n, bas, diz, tut, norm))
        elif on > yigin[0][0]:
            heapq.heapreplace(yigin, (on, _TersBaslik(bas), n, bas, diz, tut, norm))
    kisa_liste.son_taranan = n      # KAC KAYIT OKUNDU - cikti bunu yazar

    # --- ASAMA 2: gercek df -> idf -> yeniden siralama ---
    idf = [math.log(max(N, 2) / (1.0 + d)) for d in df]

    def _skor(tut, norm):
        return sum(idf[i] for i in tut) / norm

    aday = [(_skor(t, nr), b, d, len(t)) for (_o, _tb, _n, b, d, t, nr) in yigin]
    # SIRALAMA OLCUTU ACIKCA YAZILDI: skor AZALAN, esitlikte baslik ARTAN.
    # Tie-break sabitlenmeseydi esit skorlu kayitlarda dosya sirasi / heap ic
    # yapisi karar verirdi. Toplu tarayici (tum_kutu_kimlikleri.py) ayni
    # siralamayi uretebilsin diye SART.
    aday.sort(key=lambda x: (-x[0], x[1]))
    kesme = len(aday) if not ust else ust
    kl = [dict(tohum=int(a[3]), skor=round(a[0], 4), baslik=a[1], dizi=a[2],
               sira=i, kaynak='tohum')
          for i, a in enumerate(aday[:kesme], 1)]

    # Garanti kayitlari: normal siralamada NEREYE dusecekti?
    #
    # 2026-08-06 DUZELTMESI - temiz kosuda yakalandi. Garanti dizgisiyle eslesen
    # kayitlar skordan BAGIMSIZ olarak ayri torbaya alindigi icin, listeye kendi
    # gucuyle GIRECEK olanlar bile 'garanti' damgasi aliyordu. Sonuc: yeni
    # siralama olcutu sorunu cozdukten sonra bile rapor 12 iddianin 9'unda
    # "KAZANAN yamayla girdi, karar yamaya BAGIMLI" uyarisi basiyordu - oysa ayni
    # satirda "tohum siralamasinda 1. olurdu" yaziyordu. Uyari kendi kendini
    # curutuyordu ve gercek yama bagimliligini gorunmez kiliyordu.
    # Artik olcut sudur: sanal sira kesmenin ICINDEYSE kayit zaten listeye
    # girecekti, kaynak 'tohum' yazilir ve uyari BASILMAZ. Yalnizca kesmenin
    # DISINDAN gelenler 'garanti' sayilir - yamanin gercekten gerektigi yer odur.
    for tut, b, d, nr in zorunlu[:GARANTI_UST]:
        s = _skor(tut, nr)
        sanal = 1 + sum(1 for a in aday if a[0] > s)
        kendi_gucuyle = bool(kesme) and sanal <= kesme
        kl.append(dict(tohum=int(len(tut)), skor=round(s, 4), baslik=b, dizi=d,
                       sira=sanal,
                       kaynak='tohum' if kendi_gucuyle else 'garanti'))
    return kl


# -------------------------------------------------------------------------
# AYIRT EDICI PENCERE - KORUNMUS BOLGELER NEDEN DISARIDA BIRAKILIYOR
#
# 18S, 5.8S ve LSU cekirdegi gibi korunmus bolgeler butun kayitlarda neredeyse
# aynidir. Tam ortusme uzerinden olculen kimligin buyuk kismi oradan gelir ve
# birbirinden apayri iki organizma bile yuksek yuzde verir. Yani SAHTE YUKSEK
# KIMLIK tam olarak korunmus bolgelerden dogar.
#
# Bu fonksiyon en iyi N referans kaydin BIRBIRINDEN AYRILDIGI pencereyi arar:
# sorgunun 120 bazlik pencereleri, referanslarin o penceredeki kimlik DAGILIMINA
# gore puanlanir ve yayilimi en buyuk olan secilir. Korunmus bolgede yayilim
# sifira yakindir, o yuzden kendiliginden elenir.
#
# Sonuc iki sayiyi yan yana koyar: tam ortusme kimligi ve ayirt edici penceredeki
# kimlik. Bir iddia korunmus bolgeden gelen yuksek yuzdeye dayaniyorsa ikisi
# ARASINDAKI FARK bunu gorunur kilar.
# -------------------------------------------------------------------------
def ayirt_edici_pencere(kayitlar, q):
    """En iyi N referans kaydin BIRBIRINDEN ayrildigi bolgeyi bul.

    Korunmus bolgeler (18S, 5.8S, LSU cekirdegi) butun kayitlarda ayni oldugu
    icin ayrimda ISE YARAMAZ; iddia oradan gelen yuksek kimlige dayaniyorsa
    tam ortusme yuksek, ayirt edici pencere DUSUK cikar. Fark buradan gorunur.

    Basit ve saglam yaklasim: sorgunun 120 bazlik pencerelerini, referanslarin
    kendi aralarindaki kimlik DAGILIMINA gore puanlar; en cok ayrisan pencere secilir.
    """
    kayitlar = [k for k in kayitlar if k.get('dizi')][:AYIRT_EDICI_UST]
    if len(kayitlar) < 2:
        return None
    P = 120
    en = None
    for b in range(0, max(1, len(q) - P + 1), 60):
        pen = q[b:b + P]
        if len(pen) < 60 or pen.count('N') > 10:
            continue
        deg = [hizala(pen, k['dizi'])[0] for k in kayitlar]
        yay = max(deg) - min(deg)
        if en is None or yay > en[0]:
            en = (yay, b, pen, deg)
    return en


def kp_yolu(kontrol, etiket, kutu_diz, garanti, kl_ust):
    """Kontrol noktasi dosya yolu. TEK KAYNAK: hem I asamasi (vtb_tarama) hem de
    G asamasi (tum_kutu_kimlikleri.py) bunu kullanir, boylece AYNI KUTU IKI KEZ
    TARANMAZ - onbellek iki asama arasinda paylasilir."""
    import hashlib
    imza = hashlib.md5(kutu_diz.encode('utf-8')).hexdigest()[:10]
    g_im = hashlib.md5(('|'.join(sorted(garanti or ()))).encode('utf-8')).hexdigest()[:6]
    return os.path.join(kontrol, 'vtb_%s_%s_%s_%s_kl%d.json'
                        % (re.sub(r'\W+', '_', etiket), imza, g_im, AYAR_IMZASI, kl_ust))


def kl_degerlendir(kl, kutu_diz, kl_ust, taranan=None, t0=None):
    """KISA LISTE -> tam hizalama -> isabetler + oz-kalibrasyon. TEK KAYNAK.

    vtb_tarama (tek sorgu, veritabani basina bir akis) ve toplu tarayici
    (tek akista butun sorgular) BU AYNI fonksiyonu cagirir. Karar mantigi tek
    yerde durur; iki yol ayni girdide ayni hukmu vermek ZORUNDA.
    """
    t0 = t0 if t0 is not None else time.time()
    t_hiz = time.time()
    isabet = []
    for c in kl:                       # KISA LISTENIN TAMAMI hizalanir
        diz = c['dizi']
        k, d = hizala(kutu_diz if len(kutu_diz) <= len(diz) else diz,
                      diz if len(kutu_diz) <= len(diz) else kutu_diz)
        isabet.append(dict(baslik=c['baslik'][:160], kimlik=k, tohum=c['tohum'],
                           sira=c['sira'], kaynak=c['kaynak'], dizi=diz,
                           hiz_uzunluk=min(len(kutu_diz), len(diz))))
    hiz_sure = round(time.time() - t_hiz, 1)
    isabet.sort(key=lambda x: (-x['kimlik'], x['baslik']))

    # --- OZ KALIBRASYON: kazananin tohum sirasi ---
    kazanan_sira = isabet[0]['sira'] if isabet else None
    kazanan_kaynak = isabet[0]['kaynak'] if isabet else None
    uyari = None
    if kazanan_sira is not None and kazanan_sira > SIRA_UYARI_ESIGI:
        uyari = (u'KAZANAN %d. SIRADAN GELDI (esik %d, liste %d). Kesme noktasi '
                 u'BAGLAYICI olmaya baslamis olabilir - --kisa-liste degerini '
                 u'buyutup tekrarlayin.' % (kazanan_sira, SIRA_UYARI_ESIGI, kl_ust))
    if kazanan_kaynak == 'garanti':
        uyari = ((uyari + u'  ||  ') if uyari else u'') + (
            u'KAZANAN kisa listeye TOHUMLA DEGIL, "beklenen takson garantisi" '
            u'yamasiyla girdi (tohum siralamasinda %s. olurdu). Bu iddiada karar '
            u'yamaya BAGIMLI - bilmedigimiz bir taksonda ayni isabet kacardi.'
            % kazanan_sira)

    ap = ayirt_edici_pencere(isabet[:AYIRT_EDICI_UST], kutu_diz)
    for i in isabet:                   # diziler JSON'a yazilmaz (devasa olurdu)
        i.pop('dizi', None)
    return dict(durum='TAMAM', isabet=isabet[:10], kayit=len(kl),
                kisa_liste_boyu=kl_ust, hizalanan=len(kl),
                kazanan_sira=kazanan_sira, kazanan_kaynak=kazanan_kaynak,
                sira_uyarisi=uyari, taranan_kayit=taranan,
                sure=round(time.time() - t0, 1), hizalama_suresi=hiz_sure,
                ayirt_edici=(dict(yayilim=round(ap[0], 2), baslangic=ap[1],
                                  kimlikler=[round(x, 2) for x in ap[3]]) if ap else None))


def vtb_tarama(kok, kutu_diz, etiket, dosya, yaz, kontrol, garanti=(), kl_ust=KISA_LISTE,
               suzgec=None, kp_ek=''):
    """Tek veritabaninda: kisa liste -> TAMAMININ tam hizalanmasi -> en iyi isabetler.

    Kisa listedeki HER aday hizalanir (500 aday ~10 sn). Karar tamamen
    hizalamaya birakilmistir; tohum siralamasi yalnizca ADAY TOPLAMA olcutudur.
    Kazananin tohum sirasi kaydedilir - kesme noktasinin baglayici olup
    olmadiginin kanit tir.
    """
    yol = os.path.join(kok, 'REFERANS_DB', dosya)
    if not os.path.exists(yol):
        return dict(durum='dosya yok')
    # KARARLI anahtar: Python'un hash() fonksiyonu her SURECTE farkli deger
    # uretir (PYTHONHASHSEED rastgele), o yuzden kontrol noktasi hicbir zaman
    # tutmuyordu. md5 sureclerden bagimsizdir.
    # AYAR IMZASI + KISA LISTE BOYU anahtara katilir: 60'lik listeyle uretilmis
    # eski vtb_*.json dosyalari boylece gecersiz olur ve sessizce geri gelmez.
    # SUZGECLI tarama farkli bir sonuctur; ayni anahtari kullanirsa
    # suzgecsiz eski sonuc geri doner. kp_ek anahtari ayirir.
    kp = kp_yolu(kontrol, etiket + kp_ek, kutu_diz, garanti, kl_ust)
    if os.path.exists(kp):
        try:
            return json.load(open(kp, encoding='utf-8'))
        except Exception:
            pass
    t0 = time.time()

    def ilerle(n):
        print('     ... %s: %d kayit tarandi (%s)      ' % (etiket, n, sure_metni(time.time() - t0)),
              end='\r', flush=True)
    kl = kisa_liste(yol, kutu_diz, ust=kl_ust, ilerle=ilerle, garanti=garanti,
                    suzgec=suzgec)
    res = kl_degerlendir(kl, kutu_diz, kl_ust,
                         taranan=getattr(kisa_liste, 'son_taranan', None), t0=t0)
    isabet = res['isabet']
    kazanan_sira, kazanan_kaynak = res['kazanan_sira'], res['kazanan_kaynak']
    uyari, hiz_sure = res['sira_uyarisi'], res['hizalama_suresi']
    json.dump(res, open(kp, 'w', encoding='utf-8'), ensure_ascii=False, default=str)
    yaz('     %s: %s kayit TARANDI (tamami), %d kisa listeye alindi ve HEPSI '
        u'hizalandi (%s), en iyi %s%% | KAZANAN SIRA: %s/%d [%s] (%s)'
        % (etiket, '{:,}'.format(res.get('taranan_kayit') or 0).replace(',', ' '),
           len(kl), sure_metni(hiz_sure),
           vir(isabet[0]['kimlik']) if isabet else '-',
           kazanan_sira if kazanan_sira is not None else '-', kl_ust,
           kazanan_kaynak or '-', sure_metni(time.time() - t0)))
    if uyari:
        yaz(u'     >>> UYARI: %s' % uyari)
    return res


# Cins adi referans BASLIGINDAN cikarilir. Taksonomi agaci bilerek KULLANILMAZ:
# Kraken2 zaten k-mer + agacta LCA yapiyor ve bu tur onu tekrar etmemek icin var.
def cins_cek(baslik):
    """Referans basligindan cins adini cikar (taksonomi agaci KULLANILMAZ)."""
    b = baslik
    m = re.search(r'[;|]\s*([A-Z][a-z]+)[ _]([a-z]+)', b)
    if m:
        return m.group(1), '%s %s' % (m.group(1), m.group(2))
    m = re.search(r'\b([A-Z][a-z]{3,})\s+([a-z]{3,})\b', b)
    if m:
        return m.group(1), '%s %s' % (m.group(1), m.group(2))
    m = re.search(r'g__([A-Za-z_]+)', b)
    if m:
        return m.group(1), m.group(1)
    return None, None


# --------------------------------------------------------------- envanter
def envanter_yaz(kok, CIKTI, yaz):
    """REFERANS_DB altindaki BUTUN kumeleri sayar ve hangisinin nerede
    kullanildigini yazar. Kullanilmayan her kume icin SEBEP zorunludur."""
    import glob
    yol = os.path.join(CIKTI, 'VERITABANI_ENVANTERI.md')
    bilinen = {d: (e, kullan, n) for e, d, _t, kullan, n in VTB}
    diskte = sorted(os.path.basename(x) for x in
                    glob.glob(os.path.join(kok, 'REFERANS_DB', '*.fasta')) +
                    glob.glob(os.path.join(kok, 'REFERANS_DB', '*.fna')))
    satir = []
    for d in diskte:
        tam = os.path.join(kok, 'REFERANS_DB', d)
        try:
            boyut = os.path.getsize(tam)
        except OSError:
            boyut = 0
        e, kullan, n = bilinen.get(d, (None, None, None))
        satir.append(dict(dosya=d, mb=round(boyut / 1048576.0), etiket=e or '(tanimsiz)',
                          kimlik=('EVET' if kullan else 'hayir') if e else 'HAYIR - listede yok',
                          sebep=n or (u'Bu dosya VTB listesinde TANIMLI DEGIL. Kullanilmasi '
                                      u'isteniyorsa kimlik_dogrulama.py icindeki VTB '
                                      u'listesine eklenmelidir.')))
    with open(yol, 'w', encoding='utf-8') as fh:
        fh.write(u'# REFERANS_DB envanteri - hangi kume nerede kullaniliyor\n\n')
        fh.write(u'Uretim: %s (her kosuda yeniden uretilir, elle guncellenmez)\n\n'
                 % time.strftime('%Y-%m-%d %H:%M'))
        fh.write(u'| dosya | MB | etiket | KIMLIK (I) asamasinda | sebep / not |\n')
        fh.write(u'|---|---|---|---|---|\n')
        for r in satir:
            fh.write(u'| `%s` | %d | %s | **%s** | %s |\n'
                     % (r['dosya'], r['mb'], r['etiket'], r['kimlik'], r['sebep']))
        fh.write(u'\n**Kural:** kimlik asamasinda "daha temizi var" GECERLI BIR SEBEP '
                 u'DEGILDIR. Tekrarsizlastirilmis kumeler (NR99) nadir cinsleri siler; '
                 u'olculdu: SILVA LSURef NR99 icinde Petriella kaydi **0**, ayni surumun '
                 u'Parc kumesinde **82**. Bir kume ancak (a) baska bir kumenin bayt bayt '
                 u'ikizi ya da (b) baska bir kumenin altkumesi ise disarida birakilabilir; '
                 u'iki durum da tabloda gerekcesiyle yazilidir.\n\n')
        fh.write(u'**NCBI nt** yerel bir dosya degildir; ayri bir katman olarak '
                 u'sorgulanir (ag ya da elle BLAST). Bkz. rapor.\n')
    yaz(u'  envanter yazildi: %s' % yol)
    return satir


# --------------------------------------------------------------- NCBI nt katmani
NT_URL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'


def nt_katmani(kutu, dizi, CIKTI, yaz, kip='oto', bekleme=25, tur_ust=40):
    """NCBI nt'ye karsi BLAST - URL API (CMD=Put/Get). blastn -remote KULLANILMAZ.

    Basarisiz olursa SESSIZCE ATLAMAZ: elle sorgulama icin hazir girdi uretir
    ve iddia "nt katmani tamamlanmadi" olarak isaretlenir.
    """
    import urllib.request, urllib.parse
    ham = os.path.join(CIKTI, 'nt_ham')
    os.makedirs(ham, exist_ok=True)
    q = dizi[:2500]
    if kip != 'oto':
        return elle_nt_girdi(kutu, q, CIKTI, yaz)
    try:
        p = dict(CMD='Put', PROGRAM='blastn', DATABASE='nt', MEGABLAST='on',
                 QUERY=q, HITLIST_SIZE='20', FORMAT_TYPE='Text')
        with urllib.request.urlopen(NT_URL, urllib.parse.urlencode(p).encode(),
                                    timeout=90) as f:
            s = f.read().decode('utf-8', 'replace')
        m = re.search(r'RID = (\S+)', s)
        if not m:
            yaz(u'     NCBI nt: RID alinamadi - elle yola dusuluyor')
            return elle_nt_girdi(kutu, q, CIKTI, yaz)
        rid = m.group(1)
        yaz(u'     NCBI nt: is gonderildi (RID %s), bekleniyor...' % rid)
        son = ''
        for i in range(tur_ust):
            time.sleep(bekleme)
            u2 = NT_URL + '?' + urllib.parse.urlencode(
                dict(CMD='Get', RID=rid, FORMAT_TYPE='Text'))
            with urllib.request.urlopen(u2, timeout=90) as f:
                son = f.read().decode('utf-8', 'replace')
            if 'Status=WAITING' not in son:
                break
            print('     ... nt %d. yoklama          ' % (i + 1), end='\r', flush=True)
        open(os.path.join(ham, '%s.txt' % re.sub(r'\W+', '_', kutu)), 'w',
             encoding='utf-8').write(son)
        isabet = []
        for mm in re.finditer(r'^>\s*(\S.*)$', son, re.M):
            isabet.append(dict(baslik=mm.group(1).strip()[:160], kimlik=None))
        if not isabet:
            for satir in son.splitlines():
                mm = re.match(r'^(\S.{20,70}?)\s{2,}[0-9.]+\s+[0-9.]+\s', satir)
                if mm:
                    isabet.append(dict(baslik=mm.group(1).strip(), kimlik=None))
        for mm in re.finditer(r'Identities\s*=\s*\d+/\d+\s*\((\d+)%\)', son):
            if isabet:
                isabet[0]['kimlik'] = float(mm.group(1))
                break
        if not isabet:
            yaz(u'     NCBI nt: yanit ayristirilamadi - elle yola dusuluyor')
            g = elle_nt_girdi(kutu, q, CIKTI, yaz)
            g['not_'] = u'otomatik yanit ayristirilamadi; ham yanit nt_ham/ altinda'
            return g
        yaz(u'     NCBI nt: %d isabet' % len(isabet))
        return dict(durum='TAMAM', isabet=isabet[:10], kaynak='NCBI nt (URL API)')
    except Exception as e:
        yaz(u'     NCBI nt BASARISIZ (%s) - elle yola dusuluyor' % type(e).__name__)
        g = elle_nt_girdi(kutu, q, CIKTI, yaz)
        g['not_'] = u'%s: %s' % (type(e).__name__, e)
        return g


def elle_nt_girdi(kutu, q, CIKTI, yaz):
    """Aga cikilamiyorsa: elle BLAST icin hazir sorgu dosyasi + sonuc sablonu."""
    d = os.path.join(CIKTI, 'nt_elle')
    os.makedirs(d, exist_ok=True)
    fa = os.path.join(d, '%s.fasta' % re.sub(r'\W+', '_', kutu))
    with open(fa, 'w', encoding='utf-8') as fh:
        fh.write('>%s (kimlik dogrulama sorgusu, ilk %d baz)\n' % (kutu, len(q)))
        for i in range(0, len(q), 70):
            fh.write(q[i:i + 70] + '\n')
    sab = os.path.join(d, 'NT_SONUC_SABLONU.tsv')
    if not os.path.exists(sab):
        with open(sab, 'w', encoding='utf-8', newline='') as fh:
            fh.write(u'# NCBI nt sonuclarini BURAYA yazin, sonra:\n')
            fh.write(u'#   python3 KURTARMA/kimlik_dogrulama.py --kok . '
                     u'--nt-yukle KIMLIK_SONUC/nt_elle/NT_SONUC_SABLONU.tsv\n')
            fh.write(u'# Adres: https://blast.ncbi.nlm.nih.gov/Blast.cgi '
                     u'(Nucleotide BLAST, database = nt)\n')
            fh.write(u'# Sorgu dosyalari ayni klasorde: <kutu>.fasta\n')
            w = csv.writer(fh, delimiter='\t')
            w.writerow(['kutu', 'en_iyi_isabet_basligi', 'kimlik_yuzde', 'notunuz'])
    yaz(u'     NCBI nt: elle sorgu dosyasi uretildi -> %s' % os.path.basename(fa))
    return dict(durum='ELLE GEREKIR', isabet=[], kaynak='NCBI nt (elle)',
                girdi=fa, sablon=sab)


# Elle doldurulmus NT_SONUC_SABLONU.tsv okur. Bos birakilan satir "yapilmadi"
# sayilir; sifir ile bos ayni sey degildir.
def nt_yukle(yol):
    out = {}
    if not yol or not os.path.exists(yol):
        return out
    with open(yol, encoding='utf-8') as fh:
        for r in csv.DictReader((s for s in fh if not s.startswith('#')), delimiter='\t'):
            b = (r.get('en_iyi_isabet_basligi') or '').strip()
            if not b:
                continue
            try:
                k = float((r.get('kimlik_yuzde') or '').replace(',', '.'))
            except ValueError:
                k = None
            out[r['kutu'].strip()] = dict(durum='TAMAM (elle)', kaynak='NCBI nt (elle)',
                                          isabet=[dict(baslik=b, kimlik=k)])
    return out


# --------------------------------------------------------------- ADLANDIRMA
# ELDE ISIM OLMASI ile KIMLIK IDDIA ETMEK ayri seylerdir. Bu bolum ikisini
# birbirine karistirmadan raporlar: kimliklendiremedigimiz kutuda bile "en yakin
# kayit sudur, su yuzdeyle" diyebilmek icin.
#
# Tur/cins esikleri lokusa gore degisir - tek bir sayi kullanmak yaniltici olur:
TUR_ESIGI = {'SSU': 98.7, 'LSU': 98.7, 'ITS': 98.5, 'OPERON': 98.7, 'KARISIK': 98.7}
CINS_ESIGI = {'SSU': 94.5, 'LSU': 94.5, 'ITS': 90.0, 'OPERON': 94.5, 'KARISIK': 94.5}
AYRIM_PAYI = 0.5     # en iyi ile ikinci arasindaki fark bundan kucukse "cf."


# ---------------------------------------------------------------------------
# ADSIZ CEVRE KAYITLARI - 2026-08-21 HATA DUZELTMESI
#
# NCBI nt ve benzeri kumeler "adlandirilmamis" kayitlarla doludur:
#     KJ734864.1 Uncultured prokaryote clone D5 16S ribosomal RNA gene
#     KJ957653.1 Uncultured bacterium clone 4B-11 16S ribosomal RNA gene
#     GQ503828.1 Bacterium enrichment culture clone R4-53B 16S ribosomal RNA
#
# ad_coz'un ikinci duzenli ifadesi ( \b[A-Z][a-z]{3,}\s+[a-z]{3,}\b ) bunlari
# ikili ad saniyordu. OLCULDU (bu duzeltmeden once, gercek ciktilar):
#     'Uncultured prokaryote clone D5...'  -> cins='Uncultured', tur='Uncultured prokaryote'
#     'Bacterium enrichment culture...'    -> cins='Bacterium',  tur='Bacterium enrichment'
# Kimlik %99 oldugu icin savunulabilir_duzey() bunu TUR duzeyinde bir ad sayiyor
# ve iddia DOGRULANDI damgasi aliyordu. Yani CEVAPSIZLIK, onaylanmis kimlik
# gibi raporlaniyordu. KIMLIK_SONUC/kimlik_iddialari.tsv'de bunun dort ornegi
# vardi ("Uncultured prokaryote", "Uncultured bacterium" x2, "Bacterium enrichment").
#
# Bu tuzak projede ZATEN BILINIYORDU ve iki ayri yerde suzuluyordu:
#     KAPSAMLI_ARAMA/dislama_kapsama_denetimi.py:103
#     KAPSAMLI_ARAMA/siparis_siniflari.py:200
# ve KURTARMA/ncbi_katman4.py docstring'i adsiz klonlarin dislama suzgecine
# takilmadigini ayrica belgeliyor. Adlandirmanin TEK sorumlusu olan bu modulde
# yoktu; eklendi.
#
# ONEMLI: adsiz bir isabet DEGERSIZ DEGILDIR. "Kutunuz cevre klonlariyla %99
# eslesiyor" bilgi tasir. Yalnizca AD olamaz ve TUR duzeyi iddiasi kuramaz.
# Bu yuzden isabet atilmaz; adi cozulemez sayilir ve duzey 'ADLANDIRILAMIYOR'
# olur. Modulun 'adsiz' iddia tipi (satir ~200) zaten bu kavrami tasiyordu.
# ---------------------------------------------------------------------------
ADSIZ_JETONLARI = (
    'uncultured', 'unclassified', 'unidentified', 'environmental',
    'metagenome', 'enrichment', 'bacterium', 'prokaryote', 'archaeon',
    'eukaryote', 'organism', 'symbiont', 'candidate', 'clone', 'isolate',
    'synthetic', 'construct',
)


def adsiz_mi(ad):
    """Bu dizge bir TAKSON ADI mi, yoksa adlandirilmamis bir kaydin tarifi mi?

    'Petrimonas sulfuriphila' -> False   (gercek ad)
    'Uncultured bacterium'    -> True    (ad degil, tarif)
    'Bacterium enrichment'    -> True
    """
    if not ad:
        return True
    k = [x for x in re.split(r'[^A-Za-z]+', ad) if x]
    if not k:
        return True
    # Ilk kelime (cins yerinde duran sozcuk) adsiz jetonuysa bu bir ad degildir.
    if k[0].lower() in ADSIZ_JETONLARI:
        return True
    # 'Bacterium enrichment' gibi: ikinci kelime de tur epiteti degil, tarif.
    if len(k) > 1 and k[1].lower() in ADSIZ_JETONLARI:
        return True
    return False


def ad_coz(baslik):
    """Referans basligindan (cins, tur, tam_ad) cikar. Taksonomi agaci KULLANILMAZ.

    Adlandirilmamis cevre kaydi ise (cins, tur) = (None, None) doner; tam_ad
    yine doner cunku "en yakin kayit" olarak gosterilmesi gerekir.
    """
    b = (baslik or '').strip()
    tam = b[:120]
    m = re.search(r'[;|]\s*([A-Z][a-z]{2,})[ _]([a-z]{2,})', b)
    if not m:
        m = re.search(r'\b([A-Z][a-z]{3,})\s+([a-z]{3,})\b', b)
    if m:
        _c, _t = m.group(1), '%s %s' % (m.group(1), m.group(2))
        if adsiz_mi(_t):
            return None, None, tam
        return _c, _t, tam
    m = re.search(r'g__([A-Za-z_]+)', b)
    if m:
        return m.group(1), None, tam
    m = re.search(r'f__([A-Za-z_]+)', b)
    if m:
        return None, None, tam
    return None, None, tam


# ---------------------------------------------------------------------------
# ELDE ISIM OLMASI ile KIMLIK IDDIA ETMEK AYRI SEYLERDIR.
# Bu fonksiyon en iyi uc isabetten NEYIN savunulabilecegini cikarir: tur adi mi,
# "cf." mi, yalniz cins mi, yoksa ad hic verilemez mi. Esikler lokusa gore
# degisir (ITS'te tur ayrimi 16S'ten farkli calisir), tek bir sayi yaniltici
# olurdu. En iyi ile ikinci arasindaki fark AYRIM_PAYI'ndan kucukse tur atamasi
# savunulamaz ve "cf." kullanilir.
# ---------------------------------------------------------------------------
def savunulabilir_duzey(isabetler, lokus='SSU'):
    """En iyi uc isabetten SAVUNULABILIR taksonomik duzeyi cikar.

    Donen: dict(duzey, onerilen_ad, gerekce, en_iyi, ikinci, ucuncu)

    Kural - dürüst adlandirma:
      * kimlik >= tur esigi VE ikinciyle arada acik fark    -> TUR adi
      * kimlik >= tur esigi AMA ikinci baska bir tur, yakin -> "cf." (tur belirsiz)
      * tur esigi alti, cins esigi ustu                     -> CINS duzeyi
      * cins esigi alti                                     -> AILE/ustu; ad VERILMEZ,
                                                               yalniz "en yakin kayit"
    """
    say = [i for i in (isabetler or []) if isinstance(i.get('kimlik'), (int, float))]
    if not say:
        ilk = (isabetler or [{}])[0]
        _c, _t, tam = ad_coz(ilk.get('baslik', ''))
        return dict(duzey='BELIRLENEMEDI', onerilen_ad='-',
                    gerekce=u'hicbir isabette sayisal kimlik yok',
                    en_iyi=tam or '-', ikinci='-', ucuncu='-')
    say.sort(key=lambda x: -x['kimlik'])
    ilk = say[0]
    c1, t1, tam1 = ad_coz(ilk['baslik'])
    k1 = ilk['kimlik']
    tam2 = tam3 = '-'; t2 = None; k2 = None
    if len(say) > 1:
        c2, t2, tam2 = ad_coz(say[1]['baslik']); k2 = say[1]['kimlik']
    if len(say) > 2:
        _c3, _t3, tam3 = ad_coz(say[2]['baslik'])
    te = TUR_ESIGI.get(lokus, 98.7); ce = CINS_ESIGI.get(lokus, 94.5)

    # 2026-08-21: EN IYI ISABET ADSIZ ISE, KIMLIK YUZDESI NE OLURSA OLSUN AD
    # VERILEMEZ. Bu ayri bir daldir cunku asagidaki "ad verilemez" dali gerekce
    # olarak "kimlik cins esiginin altinda" yaziyor - burada durum tam tersidir:
    # kimlik %99 olabilir ama eslesilen KAYDIN KENDISI adlandirilmamistir.
    # Yanlis gerekce, dogru hukumden daha tehlikelidir: okuyan kisi "demek ki
    # daha iyi bir referans bulunursa ad cikar" diye dusunur, oysa sorun
    # referansin YAKINLIGINDA degil, ADSIZLIGINDADIR.
    if not t1 and not c1:
        return dict(duzey='ADLANDIRILAMIYOR (referans adsiz)',
                    onerilen_ad=u'adlandirilamayan soy - EN YAKIN KAYIT: %s (%%%s)'
                                % (tam1, vir(k1)),
                    gerekce=u'en iyi isabetin kimligi %%%s (tur esigi %%%s) AMA eslesilen '
                            u'kayit adlandirilmamis bir cevre dizisidir ("uncultured", '
                            u'"enrichment", "clone" vb.). Yuksek kimlik burada bir AD '
                            u'DEGIL, yalnizca "bu dizi cevre klonlariyla ortusuyor" '
                            u'bilgisidir. Tur ya da cins iddiasi KURULAMAZ.'
                            % (vir(k1), vir(te)),
                    en_iyi=tam1, ikinci=tam2, ucuncu=tam3)
    if k1 >= te and t1:
        farkli_tur = bool(t2) and t2 != t1
        yakin = (k2 is not None and (k1 - k2) < AYRIM_PAYI)
        if farkli_tur and yakin:
            return dict(duzey='CINS (tur belirsiz)', onerilen_ad='%s cf. %s' % (c1, t1.split()[-1]),
                        gerekce=u'kimlik %%%s tur esiginin (%%%s) USTUNDE ama ikinci isabet '
                                u'baska bir tur (%s, %%%s) ve arada yalnizca %%%s fark var - '
                                u'tur atamasi savunulamaz, "cf." kullanildi'
                                % (vir(k1), vir(te), t2, vir(k2), vir(k1 - k2)),
                        en_iyi=tam1, ikinci=tam2, ucuncu=tam3)
        return dict(duzey='TUR', onerilen_ad=t1,
                    gerekce=u'kimlik %%%s tur esiginin (%%%s) ustunde ve ikinci isabetle '
                            u'arada acik fark var (%s)'
                            % (vir(k1), vir(te),
                               ('%%%s' % vir(k1 - k2)) if k2 is not None else u'ikinci isabet yok'),
                    en_iyi=tam1, ikinci=tam2, ucuncu=tam3)
    if k1 >= ce and c1:
        return dict(duzey='CINS', onerilen_ad='%s sp.' % c1,
                    gerekce=u'kimlik %%%s tur esiginin (%%%s) ALTINDA, cins esiginin '
                            u'(%%%s) ustunde - tur adi verilemez'
                            % (vir(k1), vir(te), vir(ce)),
                    en_iyi=tam1, ikinci=tam2, ucuncu=tam3)
    return dict(duzey='AILE ve USTU (ad VERILEMEZ)',
                onerilen_ad=u'adlandirilamayan soy - EN YAKIN KAYIT: %s (%%%s)'
                            % ((t1 or c1 or tam1), vir(k1)),
                gerekce=u'kimlik %%%s cins esiginin (%%%s) bile ALTINDA. Bu bir KIMLIK '
                        u'DEGILDIR, yalnizca en yakin referans kayittir; hocaya '
                        u'"kimliklendiremedik, en yakini su" diye aktarilmalidir.'
                        % (vir(k1), vir(ce)),
                en_iyi=tam1, ikinci=tam2, ucuncu=tam3)


# ---------------------------------------------------------------------------
# ADLANDIRMANIN TEK KAYNAGI (2026-08-21).
#
# Bu blok eskiden ana dongunun icinde duruyordu. Ayri bir fonksiyona alindi
# cunku artik IKI yerden cagriliyor: taze kosu ve kontrol noktasindan YENIDEN
# TURETME. Iki kopya olsaydi zamanla ayrisirlardi - bu kod tabaninda ayni
# kontrolun iki ayri surumunun celismesi zaten iki kez olculmustu
# (capraz_kontrol D9 / yon_kod_taramasi, ve _kayit_coz / taksonomi).
# ---------------------------------------------------------------------------
def adlandirmayi_turet(bulgular):
    """Ham isabetlerden adlandirmayi uretir. Doner: (adlandirma_dict, lokus)"""
    lokus_tab = {e: t for e, _d, t, _k, _n in VTB}
    havuz = []
    for et, v in (bulgular or {}).items():
        if not str(v.get('durum', '')).startswith('TAMAM'):
            continue
        for i in (v.get('isabet') or [])[:5]:
            havuz.append(dict(i, _vtb=et, _lokus=lokus_tab.get(et, 'SSU')))
    sayisal = [h for h in havuz if isinstance(h.get('kimlik'), (int, float))]
    sayisal.sort(key=lambda x: -x['kimlik'])
    lokus = sayisal[0]['_lokus'] if sayisal else 'SSU'
    adl = savunulabilir_duzey(sayisal or havuz, lokus)

    # EN YAKIN BES ORGANIZMA (2026-08-21, eskiden uctu).
    #
    # ORGANIZMA bazinda TEKILLESTIRILIR, kayit bazinda degil. Gerekce: havuzda
    # 13 veritabanindan isabet var ve ayni tur cogu kumede bulunur.
    # Tekillestirilmezse "en yakin bes" kolayca AYNI organizmanin bes kaydi
    # olur ve okuyan kisi hicbir sey ogrenmez. Sorunun amaci "baska neler
    # yakin" oldugu icin ayrim ORGANIZMADA olmalidir.
    #
    # Adsiz kayitlar (uncultured/enrichment) ad uretmedigi icin tam baslikla
    # ayrilir. Her organizma icin EN YUKSEK kimlikli kayit tutulur (liste
    # zaten kimlige gore sirali oldugundan ilk gorulen odur).
    gorulen = set()
    sira = 0
    for h_ in sayisal:
        c_, t_, tam_ = ad_coz(h_['baslik'])
        anahtar = (t_ or c_ or tam_ or '').strip().lower()
        if not anahtar or anahtar in gorulen:
            continue
        gorulen.add(anahtar)
        sira += 1
        adl['isabet%d' % sira] = dict(
            tam_ad=tam_, cins=c_ or '-', tur=t_ or '-',
            kimlik=h_.get('kimlik'), uzunluk=h_.get('hiz_uzunluk'),
            vtb=h_['_vtb'])
        if sira >= 5:
            break
    return adl, lokus


def _en_yakin_etiket(isabet):
    """En yakin organizma listesinde gosterilecek etiket.

    UC DURUM AYRI AYRI GOSTERILIR - ikisini karistirmak yaniltir:

      1) Tur/cins adi cozuldu            -> adi yazilir
      2) Taksonomi VAR ama tur ikilisi yok -> EN DERIN taksonomik jeton yazilir
         (ornegin 'Dysgonomonadaceae'). Bu kayit adlandirilmamis DEGILDIR;
         yalnizca tur duzeyine inilmemistir. Ona 'adsiz' demek, SILVA'nin
         siniflandirilmis bir kaydini cevre klonuyla ayni kefeye koyar.
      3) Taksonomi de YOK (cevre klonu)  -> kaydin tanimi, 'adsiz:' onekiyle

    Bos bir '-' basmak en kotu secenektir: okuyan kisi neyle eslesildigini
    goremez ve %99'luk bir satiri sessizce ciddiye alir.
    """
    if not isabet:
        return '?'
    for alan in ('tur', 'cins'):
        v = (isabet.get(alan) or '').strip()
        if v and v != '-':
            return v
    tam = (isabet.get('tam_ad') or '').strip()
    if not tam:
        return '?'
    # Taksonomi tasiyor mu? KAPSAMLI_ARAMA/taksonomi.py ile cozulur - bes ayri
    # baslik bicimini (SILVA/UNITE/PR2/ROD/RefSeq) bilen TEK yer orasidir.
    try:
        import sys as _s, os as _o
        _kok = _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))
        if _kok not in _s.path:
            _s.path.insert(0, _kok)
        from KAPSAMLI_ARAMA import taksonomi as _TX
        _alan, _jet, _org, _tak_var = _TX.coz(tam, isabet.get('vtb') or '')
        if _tak_var and _jet:
            # En derin jeton = taksonominin en ozgul duzeyi.
            derin = _jet[-1]
            ust = _jet[-2] if len(_jet) > 1 else ''
            return u'%s%s' % (derin, (u' (%s)' % ust) if ust and ust != derin else u'')
    except Exception:
        pass
    tanim = tam.split(' ', 1)[1] if ' ' in tam else tam
    return u'adsiz: %s' % tanim[:52]


def _yeniden_turet(kayit, idd):
    """Kontrol noktasindaki HAM olcumden hukmu yeniden turetir.

    Tarama tekrarlanmaz (saatler surerdi); yalnizca saniyeler suren turetme
    yeniden yapilir. Boylece adlandirma mantigindaki bir duzeltme, pahali
    taramayi yeniden kosmadan yururluge girer.

    Yeniden turetilen: adlandirma (savunulabilir duzey, onerilen ad, en yakin
    bes organizma). Onbellekten korunan: hukum, literatur, kalibrasyon,
    vtb_detay - bunlar taramanin kendi ciktilaridir.
    """
    ham = kayit.get('bulgular') or {}
    try:
        adl, _lokus = adlandirmayi_turet(ham)
        yeni = dict(kayit)
        yeni['adlandirma'] = adl
        yeni['_turetme'] = 'yeniden (ham olcum onbellekten)'
        return yeni
    except Exception as e:
        # Turetme basarisizsa ESKI kayit AYNEN donmez; isaretlenir ki
        # rapor "bu satir eski mantikla uretildi" diyebilsin.
        yeni = dict(kayit)
        yeni['_turetme'] = 'BASARISIZ (%s) - eski adlandirma kullanildi' % type(e).__name__
        return yeni


# --------------------------------------------------------------- hukum
def _say(liste):
    """Sayisal olanlari ayikla. NCBI nt katmani kimlik=None dondurebiliyor
    (URL API metin ciktisindan yuzde ayristirilamadigi durum); None ile
    karsilastirma TypeError verir. Sayisal karsilastirmaya girmeden once
    HER ZAMAN bu suzgecten gecirilir."""
    return [x for x in liste if isinstance(x, (int, float))]


# -------------------------------------------------------------------------
# EN AZ IKI BAGIMSIZ VERITABANI SARTI - NEDEN VAR
#
# Tek bir veritabaninin en iyi isabeti kimlik iddiasi icin YETMEZ. Her kume kendi
# yanliligini tasir: tekrarsizlastirilmis kumeler nadir cinsleri siler (olculdu:
# SILVA LSURef NR99 icinde Petriella kaydi 0, ayni surumun Parc kumesinde 82),
# bir digeri ayni kaydi eskimis bir adla tasiyor olabilir. Tek kaynaga dayanan bir
# hukum, o kaynagin hatasini KIMLIK diye raporlardi.
#
# Bu yuzden ilk kontrol sayidir: sonuc veren veritabani ikiden azsa hicbir iddia
# DOGRULANDI cikamaz. Sonra oylar sayilir - >=2 veritabani iddiayi destekliyorsa
# DOGRULANDI, >=2 veritabani BASKA bir cevapta birlesiyorsa DUZELTILMELI ve dogru
# ifade yazilir, hicbir iki tanesi birlesmiyorsa DOGRULANAMADI.
#
# UYDURMA TEYIT URETILMEZ. Kanit yetersizse "DOGRULANAMADI" yazilir; bosluk,
# olumlu bir cevaba yuvarlanmaz. Bagimsizlik da denetlenir: VTB listesinde bayt
# bayt ikiz olan ve baska bir kumenin altkumesi olan dosyalar oylamadan
# cikarilmistir, yoksa ayni kayit iki kez oy verirdi.
#
# NCBI nt katmani tamamlanmadiysa hukum DOGRULANAMADI'ya cekilir (bkz. calistir):
# eksik katman sessizce atlanmaz.
# -------------------------------------------------------------------------
def hukum_ver(idd, bulgular, kons):
    """bulgular: {vtb_etiketi: sonuc}. Donen: (hukum, kanit, dogru_ifade)"""
    calisan = {k: v for k, v in bulgular.items()
               if str(v.get('durum', '')).startswith('TAMAM') and v.get('isabet')}
    if len(calisan) < 2:
        return ('DOGRULANAMADI',
                u'Yalnizca %d veritabani sonuc verdi; bu tur en az IKI bagimsiz '
                u'veritabaninin uyusmasini sart kosar.' % len(calisan), '')

    # her veritabaninin "oyu": en iyi isabetin cinsi/turu ve kimligi
    oylar = {}
    for et, v in calisan.items():
        en = v['isabet'][0]
        cins, tur = cins_cek(en['baslik'])
        oylar[et] = dict(cins=cins, tur=tur, kimlik=en['kimlik'], baslik=en['baslik'])

    tip = idd['tip']
    tol = idd.get('tolerans', UYUM_TOLERANS)

    if tip in ('kimlik', 'cins_duzeyi'):
        bek_c = idd.get('beklenen_cins'); bek_t = idd.get('beklenen_tur')
        bek_y = idd.get('beklenen_yuzde')
        if bek_c:
            uyan = [e for e, o in oylar.items() if o['cins'] and
                    o['cins'].lower().startswith(bek_c.lower()[:6])]
        elif bek_t:
            uyan = [e for e, o in oylar.items() if o['tur'] and
                    bek_t.split()[-1].lower() in (o['tur'] or '').lower()]
        else:
            uyan = list(oylar)
        kanit = '; '.join('%s: %s %%%s' % (e, oylar[e]['tur'] or oylar[e]['cins'] or '?',
                                           vir(oylar[e]['kimlik'])) for e in oylar)
        if len(uyan) >= 2:
            if bek_y is not None:
                olculen = _say([oylar[e]['kimlik'] for e in uyan])
                if olculen and min(abs(x - bek_y) for x in olculen) > tol:
                    return ('DUZELTILMELI',
                            kanit,
                            u'Takson dogru ama YUZDE farkli: olculen %s (iddia %%%s). '
                            u'Dogru ifade: "%s, olculen kimlik %%%s".'
                            % (', '.join('%%%s' % vir(x) for x in olculen), vir(bek_y),
                               (bek_t or bek_c), vir(sorted(olculen)[len(olculen) // 2])))
            return ('DOGRULANDI', kanit + u'  [%d veritabani uyusuyor]' % len(uyan), '')
        # baska bir cevapta birlesiyorlar mi
        say = {}
        for e, o in oylar.items():
            k = (o['cins'] or o['tur'] or '?')
            say.setdefault(k, []).append(e)
        en_cok = max(say.items(), key=lambda kv: len(kv[1]))
        if len(en_cok[1]) >= 2:
            med = sorted(_say([oylar[e]['kimlik'] for e in en_cok[1]])) or [None]
            return ('DUZELTILMELI', kanit,
                    u'Iddia edilen takson desteklenmiyor. %d veritabani "%s" uzerinde '
                    u'birlesiyor (kimlik ~%%%s). Dogru ifade: "kutunun en yakin '
                    u'referansi %s, kimlik %%%s".'
                    % (len(en_cok[1]), en_cok[0], vir(med[len(med) // 2]),
                       en_cok[0], vir(med[len(med) // 2])))
        return ('DOGRULANAMADI', kanit,
                u'Veritabanlari FARKLI cevaplar veriyor; hicbir iki tanesi ayni '
                u'taksonda birlesmiyor.')

    if tip == 'kutu2':
        a, b = idd['kutu'][0], idd['karsi'][0]
        if a not in kons or b not in kons:
            return ('DOGRULANAMADI', u'kutu konsensusu bulunamadi (%s / %s)' % (a, b), '')
        k, _ = hizala(kons[a] if len(kons[a]) <= len(kons[b]) else kons[b],
                      kons[b] if len(kons[a]) <= len(kons[b]) else kons[a])
        ortak = [e for e, o in oylar.items() if o['cins'] and idd.get('beklenen_cins') and
                 o['cins'].lower().startswith(idd['beklenen_cins'].lower()[:6])]
        kanit = (u'Dogrudan hizalama %s <-> %s: %%%s. Veritabani oylari: %s'
                 % (a, b, vir(k),
                    '; '.join('%s=%s' % (e, oylar[e]['tur'] or oylar[e]['cins']) for e in oylar)))
        bek = idd.get('beklenen_yuzde')
        if abs(k - bek) <= tol and len(ortak) >= 2:
            return ('DOGRULANDI', kanit, '')
        if abs(k - bek) > tol:
            return ('DUZELTILMELI', kanit,
                    u'Dogru ifade: "%s kutusu %s ile %%%s ayni" (iddia %%%s).'
                    % (a, b, vir(k), vir(bek)))
        return ('DOGRULANAMADI', kanit,
                u'Kimlik tutuyor ama iki veritabani ayni cinste birlesmedi.')

    if tip == 'dagilim':
        # supheli kutularin konsensusleri capa kutuya ne kadar benziyor
        capa = idd['kutu'][0]
        satir = []
        for k2 in idd.get('supheli', []):
            if capa in kons and k2 in kons:
                v, _ = hizala(kons[k2] if len(kons[k2]) <= len(kons[capa]) else kons[capa],
                              kons[capa] if len(kons[k2]) <= len(kons[capa]) else kons[k2])
                satir.append((k2, v))
        if not satir:
            return ('DOGRULANAMADI', u'supheli kutularin konsensusu okunamadi', '')
        kanit = '; '.join('%s=%%%s' % (a, vir(v)) for a, v in satir)
        yuksek = [v for _, v in satir if v >= 99.0]
        if len(yuksek) >= max(2, len(satir) - 1):
            return ('DOGRULANDI',
                    kanit + u'  [konsensus duzeyinde %d/%d kutu >=%%99 ayni]'
                    % (len(yuksek), len(satir)), '')
        return ('DOGRULANAMADI', kanit,
                u'Konsensus duzeyinde %d/%d kutu >=%%99 kimlikte. Okuma orani '
                u'iddiasi (%%76-86) BU TURDE OLCULMEDI - okuma duzeyi olcumu bu '
                u'betigin yontemine girmiyor; ayri bir tur gerekir.'
                % (len(yuksek), len(satir)))

    if tip == 'ayrilmaz':
        turler = idd['turler']
        bulunan = []
        for et, v in calisan.items():
            for i in v['isabet'][:6]:
                _, t = cins_cek(i['baslik'])
                for ad, bek in turler:
                    if t and ad.split()[-1].lower() in t.lower():
                        bulunan.append((et, ad, i['kimlik']))
        if not bulunan:
            return ('DOGRULANAMADI', u'iddia edilen turlerin hicbiri isabetlerde yok', '')
        kanit = '; '.join('%s: %s %%%s' % (e, a, vir(k)) for e, a, k in bulunan[:8])
        adlar = {a for _, a, _ in bulunan}
        if len(adlar) >= 2:
            kk = _say([k for _, _, k in bulunan])
            if not kk:
                return ('DOGRULANAMADI', kanit,
                        u'Iki tur de isabet ediyor ama HICBIR veritabani sayisal '
                        u'kimlik dondurmedi (yalniz baslik esdi); "ayrilmiyor" '
                        u'iddiasi bu turde sayiyla sinanamadi.')
            if max(kk) - min(kk) <= 0.5:
                return ('DOGRULANDI', kanit + u'  [iki tur arasindaki kimlik farki <=%0,5]', '')
            return ('DUZELTILMELI', kanit,
                    u'Iki tur de isabet ediyor ama kimlikler %s araliginda '
                    u'ayrisiyor; "ayrilmiyor" ifadesi bu olcumle desteklenmiyor.'
                    % ('%%%s-%%%s' % (vir(min(kk)), vir(max(kk)))))
        return ('DOGRULANAMADI', kanit, u'yalnizca bir tur isabet etti')

    if tip == 'adsiz':
        en_iyiler = [(e, o['kimlik'], o['tur'] or o['cins']) for e, o in oylar.items()]
        kanit = '; '.join('%s: %s %%%s' % (e, t, vir(k)) for e, k, t in en_iyiler)
        sayisal = [(e, k, t) for e, k, t in en_iyiler if isinstance(k, (int, float))]
        if not sayisal:
            return ('DOGRULANAMADI', kanit,
                    u'Hicbir veritabani sayisal kimlik dondurmedi; "ad verilemez" '
                    u'iddiasi bu turde sayiyla sinanamadi.')
        if all(k < 90.0 for _, k, _ in sayisal):
            return ('DOGRULANDI',
                    kanit + u'  [hicbir veritabaninda %90 ustu isabet yok - '
                    u'ad verilemez iddiasi destekleniyor]', '')
        yuksek = [(e, k, t) for e, k, t in sayisal if k >= 97.0]
        if len(yuksek) >= 2:
            return ('DUZELTILMELI', kanit,
                    u'Ad VERILEBILIR gorunuyor: %d veritabani %%97 ustu isabet '
                    u'veriyor (%s). Dogru ifade gozden gecirilmeli.'
                    % (len(yuksek), ', '.join('%s %s' % (t, vir(k)) for _, k, t in yuksek)))
        return ('DOGRULANAMADI', kanit,
                u'Isabetler %90-97 bandinda; ne "ad verilemez" ne de bir ad kesin.')

    if tip == 'gecici':
        return ('DOGRULANAMADI',
                u'Bu bir yontem beyanidir, tek bir kutunun kimligi degildir; '
                u'referans veritabani taramasiyla dogrulanamaz.',
                u'Bu satir bir OLCUM iddiasi degil, bir BELIRSIZLIK beyanidir - '
                u'oldugu gibi birakilabilir. Sinama yolu: B sinifi kutularinin '
                u'her birinin bu turdeki en iyi isabetleri karsilastirilir; ayni '
                u'taxid\'li kutular farkli taksonlara gidiyorsa beyan dogrudur.')

    return ('DOGRULANAMADI', u'bilinmeyen iddia tipi', '')


# --------------------------------------------------------------- surucu
# -------------------------------------------------------------------------
# SURUCU. Iddia iddia calisir; her iddia icin sira sabittir:
#   konsensus al -> her veritabaninda kisa liste + tam hizalama -> NCBI nt katmani
#   -> hukum_ver -> adlandirma (savunulabilir duzey + en iyi uc isabet) ->
#   literatur kontrolu -> oz-kalibrasyon -> diske yaz.
#
# KONTROL NOKTALARI IKI DUZEYLIDIR ve ikisi de AYAR_IMZASI ile kisa liste boyunu
# anahtara katar. Katmasaydi: veritabani duzeyi onbellek gecersiz kilinsa bile
# iddia duzeyi dosya bulunur, iddia tumden atlanir ve 60'lik listeyle uretilmis
# eski hukum sessizce geri gelirdi.
#
# NCBI nt onbellegi BILEREK bu muhrun disindadir (ag uzerinden gelir, kisa liste
# boyuyla ilgisi yoktur, bosa yeniden sorgulanmasin). Buna karsilik AG HATASI
# SONUC DEGILDIR ve onbelleklenmez: tek bir Wi-Fi kesintisi butun iddialari kalici
# olarak zehirliyordu.
#
# OZ-KALIBRASYON: her sorguda kazanan isabetin kisa listedeki tohum sirasi
# kaydedilir. Kesme noktasinin baglayici olup olmadigi boylece AYRI BIR OLCUM
# GEREKTIRMEZ - kanit kosunun kendisinden cikar.
# -------------------------------------------------------------------------
def calistir(kok, yalniz, sifirla, vtb_ust, nt_kip='oto', nt_yukle_yolu=None,
             lit_kip='oto', kl_ust=KISA_LISTE):
    sys.path.insert(0, kok)
    CIKTI = os.path.join(kok, 'KIMLIK_SONUC')
    KONTROL = os.path.join(CIKTI, 'kontrol')
    os.makedirs(KONTROL, exist_ok=True)
    if sifirla:
        for f in os.listdir(KONTROL):
            try:
                os.remove(os.path.join(KONTROL, f))
            except OSError as e:
                print('  silinemedi: %s (%s)' % (f, e))
    g = open(os.path.join(CIKTI, 'kosu_gunlugu.txt'), 'a', encoding='utf-8')

    def yaz(s=''):
        print(s, flush=True); g.write(s + '\n'); g.flush()

    yaz('=' * 78)
    yaz('  KIMLIK DOGRULAMA TURU - hocaya gonderilen iddialar bagimsiz sinaniyor')
    yaz('  surum %s   %s' % (VERSIYON, time.strftime('%Y-%m-%d %H:%M')))
    yaz('=' * 78)
    yaz(u'  Yontem: TOHUM + HIZALAMA (taksonomi agaci YOK, primer YOK, k-mer LCA YOK).')
    yaz(u'  Sart  : bir iddia DOGRULANDI sayilmak icin EN AZ IKI bagimsiz')
    yaz(u'          veritabaninin uyusmasi gerekir. Uydurma teyit uretilmez.')
    yaz('')

    from KAPSAMLI_ARAMA import hedefler as H
    kons = {d['kutu']: d['dizi'] for d in H.konsensusler()}
    var = [(e, d, t) for e, d, t, kullan, _n in VTB
           if kullan and os.path.exists(os.path.join(kok, 'REFERANS_DB', d))][:vtb_ust]
    envanter_yaz(kok, CIKTI, yaz)
    yaz(u'  kullanilabilir veritabani : %d  (%s)' % (len(var), ', '.join(e for e, _, _ in var)))
    if len(var) < 2:
        yaz(u'  UYARI: ikiden az veritabani var - hicbir iddia DOGRULANDI cikamaz.')

    nt_onceden = nt_yukle(nt_yukle_yolu)
    if nt_onceden:
        yaz(u'  elle yuklenen NCBI nt sonucu: %d kutu' % len(nt_onceden))
    yaz(u'  NCBI nt katmani           : %s' % {'oto': u'otomatik (URL API)', 'elle': u'elle (sorgu dosyasi uretilir)', 'yok': u'ATLANDI (istek uzerine)'}[nt_kip])

    iddialar = sorted(IDDIALAR, key=lambda x: (x['oncelik'], x['no']))
    if yalniz:
        iddialar = [i for i in iddialar
                    if (yalniz.isdigit() and int(yalniz) == i['no'])
                    or (not yalniz.isdigit() and yalniz.lower() in i['metin'].lower())]
    yaz(u'  sinanacak iddia           : %d  (oncelik sirasina gore)' % len(iddialar))
    kutu_say = len({k for i in iddialar for k in (i.get('kutu') or []) + (i.get('karsi') or [])})
    # --- SURE TAHMINI: tarama + HIZALAMA ayri ayri ---
    # Hizalama maliyeti olculerek uyduruldu (bu betigin kendi hizalayicisiyla):
    #     t ~= 6,7e-6 * kisa + 4,11e-9 * kisa * uzun     [saniye]
    # Olculen / model: 600x600 0,0050/0,0055 | 1500x1500 0,0187/0,0193 |
    #                  1500x3500 0,0310/0,0316 | 4000x4000 0,0917/0,0926
    _kutular = {k for i in iddialar for k in (i.get('kutu') or []) + (i.get('karsi') or [])}
    _uz = [min(len(kons[k]), 4000) for k in _kutular if k in kons] or [1500]
    _oq = sum(_uz) / float(len(_uz))
    _ref = 2000.0                      # tipik referans kayit uzunlugu (SSU/LSU karisik)
    _bir = 6.7e-6 * min(_oq, _ref) + 4.11e-9 * min(_oq, _ref) * max(_oq, _ref)
    _hiz_cift = _bir * kl_ust          # bir kutu-veritabani ciftinin hizalama suresi
    _tara_cift = 420                   # akis taramasi (degismedi)
    _cift = kutu_say * len(var)
    yaz(u'  TAHMINI SURE: ~%s  (%d kutu x %d veritabani = %d cift; her ciftte tam'
        % (sure_metni(_cift * (_tara_cift + _hiz_cift)), kutu_say, len(var), _cift))
    yaz(u'  veritabani akisi taranir). Kesintiye dayaniklidir.')
    yaz(u'  KISA LISTE: %d aday - HEPSI hizalanir. Hizalama payi cift basina ~%s'
        % (kl_ust, sure_metni(_hiz_cift)))
    yaz(u'  (ortalama sorgu %d bp); toplam hizalama payi ~%s.'
        % (int(_oq), sure_metni(_cift * _hiz_cift)))
    yaz(u'  2026-08-05: liste 1000 -> %d GERI indirildi. Siralama olcutu ters'
        % kl_ust)
    yaz(u'  frekans + uzunluk normalizasyonuna gecti; olculen en kotu kazanan sirasi')
    yaz(u'  35, yani %d kat pay var. Kazanc ~%s (olculdu: 1000 adayda hizalama'
        % (kl_ust // 35, sure_metni(_cift * _bir * kl_ust)))
    yaz(u'  toplam maliyetin %53\'u, 500 adayda %36\'si).')
    yaz(u'  Kesme noktasinin baglayici olup olmadigi AYRI OLCUM GEREKTIRMEZ:')
    yaz(u'  her sorguda kazanan isabetin kisa liste sirasi kaydedilir (kazanan_sira).')
    yaz('')

    sonuc = []
    tb = time.time()
    for n, idd in enumerate(iddialar, 1):
        # IDDIA DUZEYI KONTROL NOKTASI de ayar imzasini tasir. Tasimasaydi:
        # vtb_*.json gecersiz kilinsa bile iddia_01.json bulunur, iddia
        # TUMDEN atlanir ve 60'lik listeyle uretilmis eski hukum sessizce geri
        # gelirdi. (NCBI nt onbellegi -nt_*.json- BILEREK haric tutuldu: o ag
        # uzerinden gelir, kisa liste boyuyla ilgisi yok, bosa yeniden
        # sorgulanmasin.)
        # -------------------------------------------------------------------
        # 2026-08-21 MIMARI DUZELTMESI: OLCUM ile HUKUM ayni onbellekte
        # tutulmamalidir.
        #
        # Bulunan hata: bu kontrol noktasi iddianin TUMUNU sakliyordu -
        # 'adlandirma' ve 'hukum' dahil. AYAR_IMZASI ise yalnizca TARAMA
        # parametrelerini muhurluyor. Sonuc: adlandirma mantigi duzeltildiginde
        # yeniden kosu "onceki kosudan alindi" deyip ESKI HUKMU geri veriyordu.
        # Olculdu - adsiz kayit duzeltmesi uygulandiktan sonra kosulan turda
        # 12 iddianin 12'si de degismedi, cunku duzeltilmis kod HIC CALISMADI.
        #
        # Dogrusu: pahali olan TARAMADIR (saatler); adlandirma ve hukum
        # saniyeler suren TURETMELERDIR. Turetme her kosuda YENIDEN yapilmali.
        # Bunun icin kontrol noktasi HAM ISABETLERI ('bulgular') tasimalidir.
        #
        # Eski bicimli kontrol noktalari ham isabet TASIMAZ; onlar yeniden
        # turetilemez ve GECERSIZ sayilir. Sessizce eski hukum donmez.
        kp = os.path.join(KONTROL, 'iddia_%02d_%s_kl%d.json'
                          % (idd['no'], AYAR_IMZASI, kl_ust))
        if os.path.exists(kp):
            try:
                _kayit = json.load(open(kp, encoding='utf-8'))
            except Exception:
                _kayit = None
            if _kayit is not None and _kayit.get('_ham_isabet_var'):
                # Ham olcum var: TARAMA atlanir, HUKUM yeniden turetilir.
                sonuc.append(_yeniden_turet(_kayit, idd))
                yaz('[%2d/%2d] iddia %d  (tarama onbellekten, hukum YENIDEN turetildi)'
                    % (n, len(iddialar), idd['no']))
                continue
            if _kayit is not None:
                yaz('[%2d/%2d] iddia %d  ESKI BICIM kontrol noktasi (ham isabet yok) - '
                    'hukum yeniden turetilemez, YENIDEN TARANIYOR'
                    % (n, len(iddialar), idd['no']))
        yaz('[%2d/%2d] IDDIA %d (oncelik %d): %s'
            % (n, len(iddialar), idd['no'], idd['oncelik'], idd['metin']))
        if idd.get('not_'):
            yaz(u'         NOT: %s' % idd['not_'])

        bulgular = {}
        if idd['tip'] == 'gecici':
            kutular = []          # yontem beyani - veritabani taramasi gerektirmez
        else:
            kutular = (idd.get('kutu') or [])[:1] + (idd.get('karsi') or [])[:1]
        for kutu in kutular:
            if kutu not in kons:
                yaz(u'         konsensus yok: %s' % kutu)
                continue
            q = kons[kutu]
            if len(q) > 4000:
                q = q[:4000]
            for et, dosya, _tur in var:
                _gar = [x for x in (idd.get('beklenen_cins'),
                                    (idd.get('beklenen_tur') or '').split()[0]
                                    if idd.get('beklenen_tur') else None,
                                    idd.get('beklenen_aile')) if x]
                bulgular[et] = vtb_tarama(kok, q, et, dosya, yaz, KONTROL, _gar,
                                          kl_ust=kl_ust)
            # --- NCBI nt: yerel kumelerin hepsi belirli lokuslara ozeldir, nt en genisi ---
            if nt_kip != 'yok':
                if kutu in nt_onceden:
                    bulgular[NT_ETIKET] = nt_onceden[kutu]
                    yaz(u'     NCBI nt: elle yuklenen sonuc kullanildi')
                else:
                    ntk = os.path.join(KONTROL, 'nt_%s.json' % re.sub(r'\W+', '_', kutu))
                    onbellek = None
                    if os.path.exists(ntk):
                        try:                      # O-2: korumasiz json.load idi
                            onbellek = json.load(open(ntk, encoding='utf-8'))
                        except Exception as e:
                            yaz(u'     NCBI nt: onbellek bozuk, yeniden denenecek (%s)'
                                % type(e).__name__)
                            onbellek = None
                    if onbellek and str(onbellek.get('durum', '')).startswith('TAMAM'):
                        bulgular[NT_ETIKET] = onbellek
                        yaz(u'     NCBI nt: onceki kosudan alindi')
                    else:
                        # O-3: ag hatasi SONUC DEGILDIR - onbelleklenmez, her
                        # kosuda yeniden denenir. Yoksa tek Wi-Fi kesintisi
                        # butun iddialari kalici olarak zehirliyordu.
                        bulgular[NT_ETIKET] = nt_katmani(kutu, q, CIKTI, yaz, nt_kip)
                        if str(bulgular[NT_ETIKET].get('durum', '')).startswith('TAMAM'):
                            json.dump(bulgular[NT_ETIKET], open(ntk, 'w', encoding='utf-8'),
                                      ensure_ascii=False, default=str)
            break     # ilk kutu yeter; digerleri dogrudan hizalamayla kiyaslanir

        h, kanit, duzeltme = hukum_ver(idd, bulgular, kons)
        nt = bulgular.get(NT_ETIKET)
        nt_eksik = (nt_kip != 'yok' and idd['tip'] != 'gecici'
                    and (nt is None or not str(nt.get('durum', '')).startswith('TAMAM')))
        if nt_eksik:
            h = 'DOGRULANAMADI'
            kanit = (kanit + u'  ||  NCBI nt KATMANI TAMAMLANMADI (%s) - iddia bu '
                     u'yuzden dogrulanamadi sayildi, sessizce atlanmadi. Elle '
                     u'tamamlamak icin: KIMLIK_SONUC/nt_elle/ altindaki sorgu '
                     u'dosyasini BLAST edip NT_SONUC_SABLONU.tsv icine yazin, '
                     u'sonra --nt-yukle ile geri verin.'
                     % ((nt or {}).get('durum', 'kosulmadi')))
        # --- ADLANDIRMA: butun veritabanlarinin isabetlerini birlestir,
        # en iyi besi ve savunulabilir taksonomik duzeyi cikar (kisa listeler
        # zaten elimizde, yeniden tarama YOK).
        adl, lokus = adlandirmayi_turet(bulgular)
        # --- LITERATUR KONTROLU (ZORUNLU ADIM) ---
        try:
            import importlib.util as _lu
            _lp = _lu.spec_from_file_location(
                'lit', os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    'literatur_kontrol.py'))
            LIT = _lu.module_from_spec(_lp); _lp.loader.exec_module(LIT)
            _b1 = (adl.get('isabet1') or {}).get('tam_ad', '')
            lit = LIT.kontrol_et(_b1, adl.get('onerilen_ad', ''), lokus,
                                 ag=(lit_kip != 'yok'))
        except Exception as _e:
            lit = dict(durum=u'literatur modulu yuklenemedi: %s' % type(_e).__name__,
                       erisim_no='-', alan='-', vtb_adi='-', ncbi_guncel_ad='-',
                       es_anlamlilar='-', rutbe='-', ad_farkli_mi='-',
                       revizyon_uyarisi='-', revizyon_pmid='-',
                       otorite_kontrolu='GEREKLI', otorite_baglantilari='')
        lit['no'] = idd['no']; lit['onerilen_ad'] = adl.get('onerilen_ad', '-')
        yaz(u'         literatur: %s | ad farkli mi: %s | revizyon: %s'
            % (lit.get('durum', '-'), lit.get('ad_farkli_mi', '-'),
               (lit.get('revizyon_uyarisi', '-') or '-')[:48]))

        detay = {}
        for et, v in bulgular.items():
            if str(v.get('durum', '')).startswith('TAMAM') and v.get('isabet'):
                en = v['isabet'][0]
                detay[et] = dict(durum=v['durum'], en_iyi=en.get('baslik', '')[:120],
                                 kimlik=en.get('kimlik'),
                                 kazanan_sira=v.get('kazanan_sira'),
                                 kazanan_kaynak=v.get('kazanan_kaynak'),
                                 kisa_liste_boyu=v.get('kisa_liste_boyu'),
                                 taranan_kayit=v.get('taranan_kayit'))
            else:
                detay[et] = dict(durum=v.get('durum', '?'), en_iyi='', kimlik=None)
        # --- OZ KALIBRASYON OZETI (iddia duzeyinde) ---
        # Kacirma oranini AYRICA OLCMEYE GEREK YOK: veri kosunun kendisinden
        # cikiyor. Kazananlar hep ilk 100'den geliyorsa kesme noktasi baglayici
        # degildir; 400'un ustune cikan varsa 500 de yetmiyor olabilir.
        _sr = [(e, v.get('kazanan_sira'), v.get('kazanan_kaynak'))
               for e, v in bulgular.items()
               if str(v.get('durum', '')).startswith('TAMAM')
               and isinstance(v.get('kazanan_sira'), int)]
        kal = dict(
            kisa_liste_boyu=kl_ust,
            en_yuksek_kazanan_sira=(max(s for _e, s, _k in _sr) if _sr else None),
            kazanan_sira_vtb=(max(_sr, key=lambda x: x[1])[0] if _sr else '-'),
            guvenli_bolge_disi=len([1 for _e, s, _k in _sr if s > SIRA_GUVENLI_BOLGE]),
            uyari_esigi_ustu=len([1 for _e, s, _k in _sr if s > SIRA_UYARI_ESIGI]),
            garanti_ile_kazanan=len([1 for _e, _s, k in _sr if k == 'garanti']),
            olculen_vtb=len(_sr),
            vtb_siralari={e: s for e, s, _k in _sr})
        yaz(u'         oz-kalibrasyon: en yuksek kazanan sira %s/%d (%s) | '
            u'ilk %d disinda %d/%d vtb | garanti ile kazanan %d'
            % (kal['en_yuksek_kazanan_sira'] if kal['en_yuksek_kazanan_sira'] is not None else '-',
               kl_ust, kal['kazanan_sira_vtb'], SIRA_GUVENLI_BOLGE,
               kal['guvenli_bolge_disi'], kal['olculen_vtb'], kal['garanti_ile_kazanan']))

        r = dict(no=idd['no'], oncelik=idd['oncelik'], iddia=idd['metin'], hukum=h,
                 adlandirma=adl, literatur=lit, vtb_detay=detay, kalibrasyon=kal,
                 kanit=kanit, dogru_ifade=duzeltme,
                 uyusan_vtb=len([1 for v in bulgular.values() if v.get('durum') == 'TAMAM']),
                 vtb=list(bulgular.keys()),
                 ayirt_edici={e: v.get('ayirt_edici') for e, v in bulgular.items()
                              if v.get('durum') == 'TAMAM'})
        # HAM OLCUM KONTROL NOKTASINA YAZILIR (2026-08-21).
        # Boylece adlandirma/hukum mantigi degistiginde tarama YENIDEN
        # KOSULMADAN turetme tazelenebilir. Yalniz turetme icin gereken alanlar
        # saklanir (durum + ilk 5 isabet); butun kisa liste saklanmaz, dosya
        # gereksiz sismesin.
        r['bulgular'] = {
            e: dict(durum=v.get('durum'),
                    isabet=[{k: i.get(k) for k in ('baslik', 'kimlik', 'hiz_uzunluk')}
                            for i in (v.get('isabet') or [])[:5]])
            for e, v in bulgular.items()}
        r['_ham_isabet_var'] = True
        json.dump(r, open(kp, 'w', encoding='utf-8'), ensure_ascii=False, default=str)
        sonuc.append(r)
        yaz(u'         => %s' % h)
        if duzeltme:
            yaz(u'         DOGRU IFADE: %s' % duzeltme)
        gec = time.time() - tb
        print('        gecen %s | tahmini kalan %s'
              % (sure_metni(gec), sure_metni(gec / n * (len(iddialar) - n))), flush=True)

    raporla(CIKTI, sonuc, var, yaz)
    rc = cikti_denetle(yaz, 'I (KIMLIK)', [
        (os.path.join(CIKTI, 'kimlik_iddialari.tsv'), 'kimlik_iddialari.tsv')])
    g.close()
    return rc


# ---------------------------------------------------------------------------
# Uc cikti: iddia basina tek satirlik TSV, elle literatur kontrol listesi ve
# markdown rapor.
#
# Raporda "kisa liste oz-kalibrasyonu" bolumu ayri durur: kazananlarin hepsi ilk
# 100 icinden geldiyse kesme noktasi baglayici DEGILDIR ve bu acikca yazilir;
# 400'un otesinden gelen varsa 500 de yetmiyor olabilir ve --kisa-liste degerinin
# buyutulmesi onerilir. Boylece liste boyunun yeterliligi her kosuda kendini
# kanitlar ya da kendini curutur.
#
# DUZELTILMELI satirlari raporun BASINDA durur: hocaya gidecek olan onlardir.
# ---------------------------------------------------------------------------
def raporla(CIKTI, sonuc, var, yaz):
    t = os.path.join(CIKTI, 'kimlik_iddialari.tsv')
    with open(t, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# Her iddia BAGIMSIZ olarak sinandi. DOGRULANDI = en az IKI '
                 u'veritabani uyusuyor.\n')
        fh.write(u'# Yontem: tohum + hizalama. Taksonomi agaci, k-mer LCA ve primer '
                 u'KULLANILMADI.\n')
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['no', 'oncelik', 'iddia', 'HUKUM',
                    'SAVUNULABILIR_DUZEY', 'ONERILEN_AD', 'adlandirma_gerekcesi',
                    'en_iyi_isabet', 'en_iyi_cins', 'en_iyi_tur', 'en_iyi_kimlik_%',
                    'en_iyi_hiz_uzunluk', 'en_iyi_vtb',
                    'ikinci_isabet', 'ikinci_tur', 'ikinci_kimlik_%', 'ikinci_vtb',
                    'ucuncu_isabet', 'ucuncu_tur', 'ucuncu_kimlik_%', 'ucuncu_vtb',
                    'dorduncu_isabet', 'dorduncu_tur', 'dorduncu_kimlik_%', 'dorduncu_vtb',
                    'besinci_isabet', 'besinci_tur', 'besinci_kimlik_%', 'besinci_vtb',
                    'EN_YAKIN_5_ORGANIZMA',
                    'LIT_veritabani_adi', 'LIT_ncbi_guncel_ad', 'LIT_es_anlamlilar',
                    'LIT_rutbe', 'LIT_AD_FARKLI_MI', 'LIT_revizyon_uyarisi',
                    'LIT_revizyon_pmid', 'LIT_otorite_kontrolu', 'LIT_baglantilar',
                    'LIT_durum',
                    'kisa_liste_boyu', 'kazanan_sira', 'kazanan_sira_vtb',
                    'kazanan_garanti_ile_mi', 'kesme_baglayici_mi',
                    'sorgulanan_vtb_sayisi', 'sonuc_veren_vtb', 'HER_VTB_NE_DEDI',
                    'kanit', 'DOGRU_IFADE (duzeltilmeliyse)'])
        for s in sonuc:
            d = s.get('vtb_detay') or {}
            hepsi = ' | '.join('%s [%s kayit tarandi]: %s%s%s'
                               % (e, v.get('taranan_kayit', '?'),
                                  v['en_iyi'] or v['durum'],
                                             ('' if v.get('kimlik') is None
                                              else ' (%%%s)' % vir(v['kimlik'])),
                                             ('' if v.get('kazanan_sira') is None
                                              else ' {sira %s/%s%s}'
                                              % (v['kazanan_sira'], v.get('kisa_liste_boyu', '?'),
                                                 ', GARANTI'
                                                 if v.get('kazanan_kaynak') == 'garanti' else '')))
                               for e, v in d.items())
            kal = s.get('kalibrasyon') or {}
            _eys = kal.get('en_yuksek_kazanan_sira')
            if _eys is None:
                _bag = '-'
            elif kal.get('uyari_esigi_ustu'):
                _bag = u'EVET-OLABILIR (>%d)' % SIRA_UYARI_ESIGI
            elif kal.get('guvenli_bolge_disi'):
                _bag = u'HAYIR (ama ilk %d disindan geldi)' % SIRA_GUVENLI_BOLGE
            else:
                _bag = u'HAYIR (hepsi ilk %d icinde)' % SIRA_GUVENLI_BOLGE
            a = s.get('adlandirma') or {}
            def _i(n, alan, vars_=''):
                return ((a.get('isabet%d' % n) or {}).get(alan) or vars_)
            w.writerow([s['no'], s['oncelik'], s['iddia'], s['hukum'],
                        a.get('duzey', '-'), a.get('onerilen_ad', '-'), a.get('gerekce', '-'),
                        _i(1, 'tam_ad', '-'), _i(1, 'cins', '-'), _i(1, 'tur', '-'),
                        vir(_i(1, 'kimlik', None)), _i(1, 'uzunluk', '-'), _i(1, 'vtb', '-'),
                        _i(2, 'tam_ad', '-'), _i(2, 'tur', '-'),
                        vir(_i(2, 'kimlik', None)), _i(2, 'vtb', '-'),
                        _i(3, 'tam_ad', '-'), _i(3, 'tur', '-'),
                        vir(_i(3, 'kimlik', None)), _i(3, 'vtb', '-'),
                        _i(4, 'tam_ad', '-'), _i(4, 'tur', '-'),
                        vir(_i(4, 'kimlik', None)), _i(4, 'vtb', '-'),
                        _i(5, 'tam_ad', '-'), _i(5, 'tur', '-'),
                        vir(_i(5, 'kimlik', None)), _i(5, 'vtb', '-'),
                        # Tek hucrede okunabilir ozet: sirali, kimlik yuzdesiyle.
                        # ADSIZ kayitta tur/cins '-' olarak saklanir. '-' BOS
                        # sayilmazsa satir "- %99,00" diye cikar ve okuyana
                        # hicbir sey soylemez - oysa asil merak edilen NEYLE
                        # eslestigidir. Bu yuzden '-' bos kabul edilip tam
                        # basliga dusulur ve adsiz oldugu ACIKCA yazilir.
                        ' | '.join(
                            '%d) %s %%%s [%s]'
                            % (n_, _en_yakin_etiket(a.get('isabet%d' % n_)),
                               vir(_i(n_, 'kimlik', None)), _i(n_, 'vtb', '-'))
                            for n_ in range(1, 6) if a.get('isabet%d' % n_)) or '-',
                        (s.get('literatur') or {}).get('vtb_adi', '-'),
                        (s.get('literatur') or {}).get('ncbi_guncel_ad', '-'),
                        (s.get('literatur') or {}).get('es_anlamlilar', '-'),
                        (s.get('literatur') or {}).get('rutbe', '-'),
                        (s.get('literatur') or {}).get('ad_farkli_mi', '-'),
                        (s.get('literatur') or {}).get('revizyon_uyarisi', '-'),
                        (s.get('literatur') or {}).get('revizyon_pmid', '-'),
                        (s.get('literatur') or {}).get('otorite_kontrolu', 'GEREKLI'),
                        (s.get('literatur') or {}).get('otorite_baglantilari', ''),
                        (s.get('literatur') or {}).get('durum', '-'),
                        kal.get('kisa_liste_boyu', '-'),
                        '-' if _eys is None else _eys,
                        kal.get('kazanan_sira_vtb', '-'),
                        ('EVET' if kal.get('garanti_ile_kazanan') else 'hayir')
                        if kal else '-',
                        _bag,
                        len(d), s['uyusan_vtb'], hepsi, s['kanit'], s['dogru_ifade']])
    yaz('  yazildi: %s' % t)
    try:
        import importlib.util as _lu
        _lp = _lu.spec_from_file_location(
            'lit', os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'literatur_kontrol.py'))
        LIT = _lu.module_from_spec(_lp); _lp.loader.exec_module(LIT)
        el = LIT.elle_liste_yaz(CIKTI, [dict((s.get('literatur') or {}),
                                             no=s['no'],
                                             onerilen_ad=(s.get('adlandirma') or {})
                                             .get('onerilen_ad', '-')) for s in sonuc])
        yaz('  yazildi: %s' % el)
    except Exception as e:
        yaz('  elle kontrol listesi yazilamadi: %s' % type(e).__name__)

    say = {}
    for s in sonuc:
        say[s['hukum']] = say.get(s['hukum'], 0) + 1
    r = os.path.join(CIKTI, 'KIMLIK_DOGRULAMA_RAPORU.md')
    with open(r, 'w', encoding='utf-8') as fh:
        fh.write(u'# Kimlik iddialarinin bagimsiz dogrulanmasi\n\n')
        fh.write(u'Uretim: %s · betik %s\n\n' % (time.strftime('%Y-%m-%d %H:%M'), VERSIYON))
        fh.write(u'Kullanilan yerel veritabanlari (%d): %s\n\n'
                 % (len(var), ', '.join(e for e, _, _ in var)))
        fh.write(u'Ayrica **NCBI nt** ayri bir katman olarak sorgulanir. Yerel '
                 u'kumelerin hepsi belirli lokuslara ozeldir (SSU, LSU, ITS, operon); '
                 u'nt en genis olanidir ve o yuzden ayrica sorulur. nt katmani '
                 u'tamamlanmazsa iddia **DOGRULANAMADI** sayilir - sessizce atlanmaz.\n\n')
        fh.write(u'Tam kume envanteri ve kullanilmayanlarin sebebi: '
                 u'`VERITABANI_ENVANTERI.md`\n\n')
        fh.write(u'## Sonuc\n\n')
        for k in ('DOGRULANDI', 'DUZELTILMELI', 'DOGRULANAMADI'):
            if k in say:
                fh.write(u'- **%s**: %d\n' % (k, say[k]))

        # ------------------------------------------- OZ KALIBRASYON OZETI
        kals = [s.get('kalibrasyon') or {} for s in sonuc]
        siralar = [(s['no'], e, v) for s in sonuc
                   for e, v in ((s.get('kalibrasyon') or {}).get('vtb_siralari') or {}).items()
                   if isinstance(v, int)]
        fh.write(u'\n## Kisa liste oz-kalibrasyonu (kesme noktasi baglayici mi?)\n\n')
        if not siralar:
            fh.write(u'Olculebilir kazanan sirasi yok (bu kosuda veritabani taramasi '
                     u'yapilmamis olabilir).\n\n')
        else:
            ust = max(v for _n, _e, v in siralar)
            n100 = len([1 for _n, _e, v in siralar if v > SIRA_GUVENLI_BOLGE])
            n400 = len([1 for _n, _e, v in siralar if v > SIRA_UYARI_ESIGI])
            gar = sum(k.get('garanti_ile_kazanan', 0) for k in kals)
            boy = next((k.get('kisa_liste_boyu') for k in kals
                        if k.get('kisa_liste_boyu')), KISA_LISTE)
            fh.write(u'Kisa liste **tohum sayisina** gore kuruluyor ama karari '
                     u'**hizalama kimligi** veriyor. Iki olcut farkli oldugu icin '
                     u'kesme noktasi gercek en iyi eslesmeyi eleyebilir. Bunu ayrica '
                     u'olcmek yerine **her sorguda kazanan isabetin kisa listedeki tohum '
                     u'sirasi** kaydedildi - kanit kosunun kendisinden cikiyor.\n\n')
            fh.write(u'| olcu | deger |\n|---|---|\n')
            fh.write(u'| kisa liste boyu | %d (HEPSI hizalandi) |\n' % boy)
            fh.write(u'| olculen sorgu (iddia x veritabani) | %d |\n' % len(siralar))
            fh.write(u'| **en yuksek kazanan sira** | **%d** |\n' % ust)
            fh.write(u'| kazanan ilk %d disindan geldi | %d / %d sorgu |\n'
                     % (SIRA_GUVENLI_BOLGE, n100, len(siralar)))
            fh.write(u'| kazanan %d. siranin otesinden geldi | %d / %d sorgu |\n'
                     % (SIRA_UYARI_ESIGI, n400, len(siralar)))
            fh.write(u'| kazanan "beklenen takson garantisi" ile geldi | %d sorgu |\n\n' % gar)
            if n400:
                fh.write(u'> **UYARI.** %d sorguda kazanan %d. siranin otesinden geldi. '
                         u'Kesme noktasi hala baglayici olabilir; `--kisa-liste` degerini '
                         u'buyutup (orn. %d) tekrarlayin.\n\n'
                         % (n400, SIRA_UYARI_ESIGI, boy * 2))
            elif n100:
                fh.write(u'> Kazananlarin %d tanesi ilk %d disindan geldi - yani **eski '
                         u'60\'lik liste bu isabetleri KACIRIRDI**. Hepsi %d\'un altinda '
                         u'kaldigi icin %d\'luk liste yeterli.\n\n'
                         % (n100, SIRA_GUVENLI_BOLGE, SIRA_UYARI_ESIGI, boy))
            else:
                fh.write(u'> Butun kazananlar ilk %d icinden geldi. Kesme noktasi '
                         u'**baglayici degil**: karari tamamen hizalama veriyor, kisa '
                         u'liste boyu sonucu etkilemiyor.\n\n' % SIRA_GUVENLI_BOLGE)
            if gar:
                fh.write(u'> **Dikkat.** %d sorguda kazanan kisa listeye tohumla degil, '
                         u'"beklenen takson garantisi" yamasiyla girdi. O sorgularda karar '
                         u'yamaya BAGIMLIDIR - ne aradigimizi bilmedigimiz bir taksonda '
                         u'ayni isabet kacardi.\n\n' % gar)
            enb = sorted(siralar, key=lambda x: -x[2])[:5]
            fh.write(u'En yuksek bes kazanan sira: %s\n\n'
                     % '; '.join(u'iddia %d / %s: %d' % t for t in enb))

        fh.write(u'\n> **Yontem farki.** Kraken2 k-mer + taksonomi agacinda LCA yapar; '
                 u'onceki turlarimiz konsensusleri birbiriyle ve in-silico PCR ile '
                 u'kiyasladi. Bu tur ise DIS referans veritabanlarindaki adlandirilmis '
                 u'kayitlarla tohum+hizalama yapar ve kimligi ayrica AYIRT EDICI '
                 u'PENCEREDE olcer (korunmus bolgeler disarida kalir). Bir iddia '
                 u'DOGRULANDI sayilmak icin en az IKI bagimsiz veritabani uyusmalidir.\n\n')
        for etiket, baslik in (('DUZELTILMELI', u'## Duzeltilmesi gerekenler (hocaya gidecek)'),
                               ('DOGRULANDI', u'## Dogrulananlar'),
                               ('DOGRULANAMADI', u'## Dogrulanamayanlar')):
            uy = [s for s in sonuc if s['hukum'] == etiket]
            if not uy:
                continue
            fh.write(baslik + u'\n\n')
            for s in uy:
                fh.write(u'### %d. %s\n\n' % (s['no'], s['iddia']))
                fh.write(u'- **Hukum:** %s  (uyusan veritabani: %d)\n' % (s['hukum'], s['uyusan_vtb']))
                a = s.get('adlandirma') or {}
                if a:
                    fh.write(u'- **Savunulabilir duzey:** `%s` → **%s**\n'
                             % (a.get('duzey', '-'), a.get('onerilen_ad', '-')))
                    fh.write(u'  - *Gerekce:* %s\n' % a.get('gerekce', '-'))
                    fh.write(u'\n  **En yakin bes ORGANIZMA** (kayit degil organizma '
                             u'bazinda tekillestirildi; ayni tur birden cok '
                             u'veritabaninda bulundugu icin):\n\n')
                    fh.write(u'\n  | # | en yakin kayit | cins | tur | kimlik | veritabani |\n'
                             u'  |---|---|---|---|---|---|\n')
                    for n_ in (1, 2, 3, 4, 5):
                        it = a.get('isabet%d' % n_)
                        if it:
                            fh.write(u'  | %d | %s | %s | %s | %%%s | %s |\n'
                                     % (n_, it['tam_ad'], it['cins'], it['tur'],
                                        vir(it['kimlik']), it['vtb']))
                    fh.write(u'\n  > **Elde isim olmasi ile kimlik iddia etmek ayri seylerdir.** '
                             u'Yukaridaki "en yakin kayit" bir KIMLIK DEGILDIR; savunulabilir '
                             u'duzey sutunu neyin iddia edilebilecegini soyler.\n\n')
                fh.write(u'- **Kanit:** %s\n' % s['kanit'])
                d = s.get('vtb_detay') or {}
                if d:
                    fh.write(u'\n  **Sorgulanan veritabanlari ve her birinin dedigi '
                             u'(%d kaynak):**\n\n' % len(d))
                    fh.write(u'  | veritabani | durum | en iyi isabet | kimlik | kazanan sira |\n'
                             u'  |---|---|---|---|---|\n')
                    for e, v in d.items():
                        _ks = ('-' if v.get('kazanan_sira') is None
                               else u'%d / %s%s' % (v['kazanan_sira'],
                                                    v.get('kisa_liste_boyu', '?'),
                                                    u' **(GARANTI)**'
                                                    if v.get('kazanan_kaynak') == 'garanti'
                                                    else ''))
                        fh.write(u'  | %s | %s | %s | %s | %s |\n'
                                 % (e, v['durum'], v['en_iyi'] or '-',
                                    '-' if v.get('kimlik') is None else '%%%s' % vir(v['kimlik']),
                                    _ks))
                    fh.write(u'\n  > *Kazanan sira*: en yuksek kimligi veren kaydin kisa '
                             u'listedeki TOHUM sirasi. Kucuk sayilar kesme noktasinin '
                             u'baglayici olmadigini gosterir.\n\n')
                if s['dogru_ifade']:
                    fh.write(u'- **DOGRU IFADE:** %s\n' % s['dogru_ifade'])
                ae = {k: v for k, v in (s.get('ayirt_edici') or {}).items() if v}
                if ae:
                    fh.write(u'- **Ayirt edici pencere:** %s\n'
                             % '; '.join('%s yayilim %%%s' % (k, vir(v.get('yayilim')))
                                         for k, v in ae.items()))
                fh.write(u'\n')
    yaz('  yazildi: %s' % r)
    yaz('')
    yaz('  ' + '   '.join('%s: %d' % kv for kv in say.items()))



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

# Komut satiri: --yalniz tek iddia, --vtb-ust kac veritabani, --kisa-liste tam
# hizalanacak aday sayisi (degistirmek eski kontrol noktalarini gecersiz kilar),
# --nt NCBI kipi, --nt-yukle elle doldurulmus sablon, --literatur, --sifirla.
def main():
    p = argparse.ArgumentParser(description='Kimlik iddialarinin bagimsiz dogrulanmasi')
    p.add_argument('--kok', default='.')
    p.add_argument('--yalniz', default=None, help='iddia numarasi ya da metninden parca')
    p.add_argument('--vtb-ust', type=int, default=len(VTB), help='kac veritabani kullanilsin')
    p.add_argument('--kisa-liste', type=int, default=KISA_LISTE, dest='kisa_liste',
                   help=(u'her veritabanindan TAM HIZALANACAK aday sayisi (varsayilan %d). '
                         u'Buyuk deger kesme noktasini baglayici olmaktan cikarir; '
                         u'raporda "kazanan sira" sutunu yeterli olup olmadigini soyler. '
                         u'Degistirmek eski kontrol noktalarini gecersiz kilar.'
                         % KISA_LISTE))
    p.add_argument('--nt', choices=['oto', 'elle', 'yok'], default='oto',
                   help='NCBI nt katmani: oto (URL API), elle (sorgu dosyasi uret), yok')
    p.add_argument('--nt-yukle', default=None,
                   help='doldurulmus NT_SONUC_SABLONU.tsv')
    p.add_argument('--literatur', choices=['oto', 'yok'], default='oto',
                   help='NCBI Taxonomy + PubMed literatur kontrolu')
    p.add_argument('--sifirla', action='store_true')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)
    if not os.path.isdir(os.path.join(kok, 'KAPSAMLI_ARAMA')):
        sys.exit('HATA: %s icinde KAPSAMLI_ARAMA yok.' % kok)
    if a.kisa_liste < 1:
        sys.exit('HATA: --kisa-liste en az 1 olmali.')
    return calistir(kok, a.yalniz, a.sifirla, a.vtb_ust, a.nt, a.nt_yukle, a.literatur,
                    kl_ust=a.kisa_liste)


if __name__ == '__main__':
    sys.exit(main() or 0)
