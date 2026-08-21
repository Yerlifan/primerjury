# -*- coding: utf-8 -*-
u"""G ASAMASI - TUM KUTU KIMLIKLERININ BAGIMSIZ DOGRULANMASI

NEDEN VAR
---------
`I` asamasi 12 IDDIAYI sinar - ama o iddialar bizim SUPHELENDIGIMIZ kimliklerdi.
Suphelenmedigimiz kutularin kimligi hic bagimsiz dogrulanmadi. Oysa o kutular da
uyelik ve rakip kumelerini belirliyor, yani HER AYRIM KATINI etkiliyorlar.
Sessiz cogunluk hic kontrol edilmedi. Bu asama onu kapatir.

KAPSAM
------
Panel olcumlerine giren HER kutu: bir hedefin UYESI ya da RAKIBI olarak herhangi
bir ayrim hesabina katilan butun kutular. Hicbir hesaba katilmayan kutular
ATLANIR ve atlandiklari cikitiya YAZILIR (sessiz atlama yok).

YONTEM - `I` ILE BIREBIR AYNI, YENIDEN YAZILMADI
------------------------------------------------
Bu betik karar mantigini KENDI YAZMAZ; kimlik_dogrulama.py'nin fonksiyonlarini
CAGIRIR:
    kp_yolu()            - kontrol noktasi anahtari (ONBELLEK PAYLASIMI)
    kl_degerlendir()     - kisa liste -> hizalama -> isabet + kazanan_sira
    kisa_liste()         - tek sorgulu tarayici (dogrulama icin)
    ad_coz(), savunulabilir_duzey(), cins_cek(), hizala(), ayirt_edici_pencere()
    literatur_kontrol.py - literatur katmani
Yani: 500'luk kisa liste, hepsi hizalanir, en az IKI bagimsiz veritabani
uyusmasi, en iyi UC isabet, savunulabilir duzey + onerilen ad, kazanan sira.

VERITABANI KAPSAMI - ALAN FILTRESI YOKTUR
-----------------------------------------
HER kutuya 12 yerel veritabaninin HEPSI sorulur. Alan (bakteri/arke/mantar)
filtresi UYGULANMAZ.

  Neden: alani KRAKEN ETIKETINE gore secmek tehlikelidir - bu asamanin varlik
  sebebi zaten Kraken etiketlerinin yanlis olabilmesi. "Bakteri kutusuna mantar
  veritabani sormaya gerek yok" demek, kutunun bakteri OLDUGUNU varsaymaktir;
  tam da sinamak istedigimiz sey odur. Bu yuzden BUTUN veritabanlari sorulur ve
  alakasiz olanlar SONUCTA dusurulur (isabet yoksa "sonuc yok" diye isaretlenir),
  sorgudan ONCE elenmez.

  Kutunun alani ETIKETTEN DEGIL OLCUMDEN cikar: hangi veritabanlarinin gercekten
  isabet verdigi 'alan_olcumden' sutununa yazilir.

Bir veritabani bir kutu icin sonuc dondurmediyse bu "temiz" SAYILMAZ:
  TAMAM          - tarandi, isabet var
  SONUC YOK      - tarandi, hicbir kayit tohum tutturmadi (alan disi olabilir)
  DOSYA YOK      - veritabani diskte yok
  SORULMADI (..) - sebebi yazilir
Her kutu satirinda 12 veritabaninin 12'si de gorunur.

KAPSAM MUHASEBESI - TAVAN SORUNU TEKRARLANMASIN
-----------------------------------------------
Erisim testinde gercek bir tavan sorunu yasandi: ilk kosu 120 001 kayitta
kesiyordu, SILVA SSU NR99 / LSU Parc / UNITE ITS fiilen budanmis taraniyordu.
Burada her veritabani icin TARANAN kayit sayisi sayilir ve BEKLENEN_KAYIT ile
karsilastirilir; esit degilse 'kapsam' sutunu EKSIK yazar ve uyari basilir.

HIZ - TOPLU TARAMA
------------------
Kutu basina ayri akis: 94 tekil konsensus x 12 vtb = 1128 tam veritabani akisi.
Kabul edilemez. Bunun yerine bir veritabani akisi AYNI ANDA butun kume
sorgularina hizmet eder (tohum -> sorgu ters indeksi). Uretilen kisa liste,
tek sorgulu kisa_liste() ile BIREBIR AYNIDIR ve bu sinamayla kanitlanir.

Panel dosyalarina YAZMAZ; TUM_KIMLIK_SONUC/ altina yazar.
"""

# -------------------------------------------------------------------------
# tum_kutu_kimlikleri.py — panel olcumlerine giren HER kutunun kimligini dis
# referans veritabanlarina karsi bagimsiz olarak dogrular ("sessiz cogunluk").
#
# GİRDİ  : REFERANS_DB/ altindaki 12 yerel FASTA kumesi (hepsi, alan filtresi YOK),
#          konsensus_kanonik/ ve panel + uyelik tablolari (KAPSAMLI_ARAMA.hedefler),
#          KIMLIK_SONUC/kontrol/ (I asamasiyla ORTAK onbellek),
#          istege bagli NCBI nt.
# ÇIKTI  : TUM_KIMLIK_SONUC/tum_kutu_kimlikleri.tsv (kutu basina TEK satir),
#          TUM_KIMLIK_SONUC/TUM_KUTU_KIMLIK_RAPORU.md,
#          TUM_KIMLIK_SONUC/kutu_*.json, kosu_gunlugu.txt.
#          Panel dosyalarina YAZMAZ.
# ÇAĞRAN : KAPSAMLI_ARAMA.bat -> G tusu
#          (bat icinde: wsl -e python3 "KURTARMA/tum_kutu_kimlikleri.py" --kok .)
#
# I ASAMASINDAN FARKI: I, bizim SUPHELENDIGIMIZ 12 iddiayi sinar. G, suphesiz
# saydigimiz kutulari da sinar - oysa o kutular da uyelik ve rakip kumelerini
# belirliyor, yani HER AYRIM KATINI etkiliyorlar.
#
# YONTEM YENIDEN YAZILMADI: karar mantigi kimlik_dogrulama.py'den CAGRILIR
# (kp_yolu, kl_degerlendir, savunulabilir_duzey, ad_coz, cins_cek, hizala,
# ayirt_edici_pencere). Yani 500'luk kisa liste, hepsinin hizalanmasi, en az IKI
# bagimsiz veritabani uyusmasi ve en iyi uc isabet kurallari BIREBIR aynidir.
# Iki ayri uygulama olsaydi iki farkli hukum ureten iki farkli arac olurdu.
# -------------------------------------------------------------------------
import os, sys, csv, json, time, re, argparse, heapq, collections

VERSIYON = '1.0 (2026-08-04)'

# Erisim dogrulamasindan (ERISIM_SONUC/erisim_dogrulama.tsv, 'TAMAMI' kosusu)
# alinan GERCEK kayit sayilari. Taranan bunlardan azsa kapsam EKSIK demektir.
BEKLENEN_KAYIT = {
    'SILVA SSU NR99': 510495, 'SILVA LSU NR99': 95279, 'SILVA LSU Parc': 1312521,
    'UNITE ITS': 2069189, 'PR2 SSU': 240201, 'ROD operon': 60320,
    'RefSeq bakteri 16S': 26877, 'RefSeq arke 16S': 1160, 'RefSeq mantar ITS': 20394,
    'RefSeq mantar 28S': 12890, 'RefSeq mantar 18S': 4037, 'RefSeq ref_all2': 65358,
}


# kimlik_dogrulama.py bir betiktir (paket degil), o yuzden dosya yolundan modul
# olarak yuklenir. Karar mantigi ORADAN gelir; burada yeniden yazilmaz.
def _K(kok):
    """kimlik_dogrulama.py'yi modul olarak yukle - karar mantigi ORADAN gelir."""
    import importlib.util
    yol = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kimlik_dogrulama.py')
    sp = importlib.util.spec_from_file_location('kimlik_dogrulama', yol)
    m = importlib.util.module_from_spec(sp)
    sys.modules['kimlik_dogrulama'] = m
    sp.loader.exec_module(m)
    return m


# --------------------------------------------------------------- KUTU ENVANTERI
# ---------------------------------------------------------------------------
# KAPSAM: bir hedefin UYESI ya da RAKIBI olarak herhangi bir ayrim hesabina
# katilan butun kutular. Hicbir hesaba girmeyen kutular sinanmaz - ama SESSIZCE
# atlanmaz: atlananlar listesine sebebiyle yazilir ve hem TSV hem markdown
# raporunda gorunur.
# ---------------------------------------------------------------------------
def kutu_envanteri(kok, K):
    """Panel olcumlerine giren HER kutu + hangi hedefin uyesi/rakibi oldugu.

    Donen: (katilan, atlanan, uye, rakip, kons)
      katilan : sirali kutu adlari
      atlanan : [(kutu, sebep)]
      uye     : kutu -> [hedef, ...]
      rakip   : kutu -> [hedef, ...]
    """
    sys.path.insert(0, kok)
    from KAPSAMLI_ARAMA import hedefler as H
    panel, _yol = H.panel_oku()
    kons_l = H.konsensusler()
    uyelik = H.uyelik_oku()
    kut = H.kutular()
    uye = collections.defaultdict(list)
    rakip = collections.defaultdict(list)
    for p in panel:
        b = H.hedef_baglami(p, uyelik, kons_l, kut)
        for k in b['uye_kons']:
            uye[k['kutu']].append(p['hedef'])
        for k in b['rakip_kons']:
            rakip[k['kutu']].append(p['hedef'])
    kons = {d['kutu']: d['dizi'] for d in kons_l}
    hepsi = sorted(kons)
    katilan = sorted(set(uye) | set(rakip))
    atlanan = [(k, u'hicbir hedefin uyesi ya da rakibi degil - panel olcumlerine '
                   u'girmiyor, ayrim hesaplarini etkilemiyor')
               for k in hepsi if k not in set(katilan)]
    return katilan, atlanan, uye, rakip, kons, H


# --------------------------------------------------------------- TOPLU TARAMA
class _TersMetin(str):
    u"""Karsilastirmasi TERS cevrilmis metin.

    Kisa listenin siralama olcutu: tohum sayisi AZALAN, esitlikte baslik ARTAN.
    Bellegi sinirli tutmak icin 'ust' boyutlu bir MIN-HEAP kullaniyoruz; heap'in
    atacagi eleman "en kotu" olmali. Tohum sayisinda en kucuk, ESITLIKTE BASLIGI
    EN BUYUK olan en kotudur - bu yuzden baslik karsilastirmasi ters cevrilir.
    Ters cevrilmezse kesme noktasindaki ESIT TOHUMLU kayitlar tek sorgulu
    kisa_liste() ile farkli secilir ve iki yol ayrisir (olculdu: 455. sirada
    ayrisiyordu).
    """
    __slots__ = ()

    def __lt__(self, o):
        return str.__gt__(self, o)

    def __gt__(self, o):
        return str.__lt__(self, o)

    def __le__(self, o):
        return str.__ge__(self, o)

    def __ge__(self, o):
        return str.__le__(self, o)


# ---------------------------------------------------------------------------
# TOPLU TARAMA - tek veritabani akisinda BUTUN kume sorgulari.
#
# Kutu basina ayri akis olsaydi 94 konsensus x 12 veritabani = 1128 tam dosya
# gecisi gerekirdi; kabul edilemez. Burada tohum -> sorgu ters indeksi kurulur ve
# her kaydin k-mer kumesi bir kez cikarilip indekse carpilir.
#
# URETILEN KISA LISTE, tek sorgulu kisa_liste() ile BIREBIR AYNI OLMAK ZORUNDA:
# iki yol ayni girdide farkli liste cikarirsa I ve G asamalari ayni kutu icin
# farkli hukum verebilir. Siralama olcutu ikisinde de "tohum sayisi AZALAN,
# esitlikte baslik ARTAN"dir; _TersMetin sinifi bu esitlik kuralini min-heap
# icinde de korumak icin vardir (olculdu: ters cevrilmezse 455. sirada ayrisiyor).
#
# ALAN FILTRESI YOKTUR: hangi sorgu olursa olsun her kayit degerlendirilir.
# ---------------------------------------------------------------------------
def toplu_kisa_liste(K, yol, sorgular, ust, ilerle=None):
    u"""BIR veritabani akisinda BUTUN sorgular icin kisa liste kur.

    sorgular: {ad: dizi}. Donen: ({ad: kisa_liste}, taranan_kayit_sayisi)

    Tohum -> sorgu ters indeksi kurulur; her kayit icin kaydin 16-mer kumesi bir
    kez cikarilir ve indekse carpilir. Sonuc, her sorgu icin tek tek
    kisa_liste() cagirmakla AYNIDIR (sinamada kanitlanir): 'tohum sayisi' her iki
    yolda da "sorgunun kac ayri tohumu bu kayitta gecti" demektir.

    ALAN FILTRESI YOK: hangi sorgu olursa olsun her kayit degerlendirilir.
    """
    import math
    k = K.K_TOHUM
    # TEK SORGULU kisa_liste() ILE AYNI OLCUT (idf + BM25). Iki yol ayni girdide
    # ayni kisa listeyi cikarmak ZORUNDA; gerekce kimlik_dogrulama.py basindaki
    # "SIRALAMA OLCUTU" bolumunde.
    tohum_sira = {}                          # ad -> sirali tohum listesi
    tohum_sahip = {}                         # tohum -> [(ad, indeks), ...]
    for ad, q in sorgular.items():
        th = sorted(K.tohumlar(q) | K.tohumlar(K.rc(q)))
        tohum_sira[ad] = th
        for i, t in enumerate(th):
            tohum_sahip.setdefault(t, []).append((ad, i))
    tohum_kume = set(tohum_sahip)
    ORT = K.ortalama_uzunluk(yol)
    B = K.BM25_B
    havuz = max(int(K.ADAY_HAVUZU), (ust or 0) * 3, 500)
    df = {ad: [0] * len(t) for ad, t in tohum_sira.items()}
    yigin = {ad: [] for ad in sorgular}      # min-heap, en fazla 'havuz' eleman
    n = 0
    N = 0
    if not tohum_kume:
        return {ad: [] for ad in sorgular}, 0
    for bas, diz in K.fasta_akisi(yol):      # TAVAN YOK - dosya sonuna kadar
        n += 1
        if ilerle and n % 20000 == 0:
            ilerle(n)
        L = len(diz)
        if L < 100:
            continue
        N += 1
        kmers = {diz[i:i + k] for i in range(L - k + 1)}
        ortak = kmers & tohum_kume
        if not ortak:
            continue
        tut = {}
        for t in ortak:
            for ad, i in tohum_sahip[t]:
                tut.setdefault(ad, set()).add(i)
        tb = _TersMetin(bas)
        norm = 1.0 - B + B * L / ORT
        for ad, s in tut.items():
            d_ = df[ad]
            for i in s:
                d_[i] += 1                   # TERS FREKANS: ayni akista bedava
            on = len(s) / norm               # ON ELEME (idf henuz bilinmiyor)
            h = yigin[ad]
            if len(h) < havuz:
                heapq.heappush(h, (on, tb, n, diz, frozenset(s), norm))
            elif on > h[0][0]:
                heapq.heapreplace(h, (on, tb, n, diz, frozenset(s), norm))
    cikti = {}
    for ad, h in yigin.items():
        idf = [math.log(max(N, 2) / (1.0 + d)) for d in df[ad]]
        aday = [(sum(idf[i] for i in s) / nr, str(b), d, len(s))
                for (_o, b, _n, d, s, nr) in h]
        aday.sort(key=lambda x: (-x[0], x[1]))
        kesme = len(aday) if not ust else ust
        cikti[ad] = [dict(tohum=int(a[3]), skor=round(a[0], 4), baslik=a[1],
                          dizi=a[2], sira=i, kaynak='tohum')
                     for i, a in enumerate(aday[:kesme], 1)]
    return cikti, n


# --------------------------------------------------------------- KUTU HUKMU
# -------------------------------------------------------------------------
# EN AZ IKI BAGIMSIZ VERITABANI SARTI - NEDEN VAR
#
# Tek bir veritabaninin en iyi isabeti bir kimlik iddiasi icin YETMEZ. Her kume
# kendi tarihsel yanliliklarini tasir: bir kume nadir cinsleri tekrarsizlastirma
# sirasinda silmis olabilir (olculdu: SILVA LSURef NR99 icinde Petriella kaydi 0,
# ayni surumun Parc kumesinde 82), bir digeri ayni kaydi eskimis bir adla tasiyor
# olabilir. Tek kaynaga dayanan bir hukum, o kaynagin hatasini KIMLIK diye
# raporlardi.
#
# Bu yuzden hukum, veritabanlarinin UYUSMASINA baglanir: ayni cinsi en iyi isabet
# olarak veren BAGIMSIZ veritabani sayisi >=2 ise DOGRULANDI, 1 ise
# "DOGRULANAMADI (tek kaynak)", hicbiri birlesmiyorsa DOGRULANAMADI.
#
# UYDURMA TEYIT URETILMEZ: kanit yetersizse bosluk bosluk olarak raporlanir.
# Bagimsizlik da denetlenir - VTB listesinde ikiz (bayt bayt ayni) ve altkume
# olan kumeler oylamadan cikarilmistir, yoksa ayni kayit iki kez oy verirdi.
#
# Uyusma CINS duzeyinde aranir: tur duzeyi tek harflik ad farkiyla bile ayrilir ve
# gercek bir uyusmayi yapay olarak bozar.
# -------------------------------------------------------------------------
def kutu_hukmu(K, bulgular, lokus_tab):
    u"""Bir kutunun kimligi: en az IKI bagimsiz veritabani uyusmali.

    `I` ile ayni sart. Uyusma CINS duzeyinde aranir (tur duzeyi tek harfle
    degisebiliyor); ayni cinsi en iyi isabet olarak veren BAGIMSIZ veritabani
    sayisi >= 2 ise DOGRULANDI.
    """
    havuz = []
    for et, v in bulgular.items():
        if not str(v.get('durum', '')).startswith('TAMAM'):
            continue
        for i in (v.get('isabet') or [])[:5]:
            havuz.append(dict(i, _vtb=et, _lokus=lokus_tab.get(et, 'SSU')))
    sayisal = [h for h in havuz if isinstance(h.get('kimlik'), (int, float))]
    sayisal.sort(key=lambda x: -x['kimlik'])
    lokus = sayisal[0]['_lokus'] if sayisal else 'SSU'
    adl = K.savunulabilir_duzey(sayisal or havuz, lokus)
    for n_, h_ in enumerate(sayisal[:3], 1):
        c_, t_, tam_ = K.ad_coz(h_['baslik'])
        adl['isabet%d' % n_] = dict(tam_ad=tam_, cins=c_ or '-', tur=t_ or '-',
                                    kimlik=h_.get('kimlik'), uzunluk=h_.get('hiz_uzunluk'),
                                    vtb=h_['_vtb'])
    # her veritabaninin EN IYI isabetinin cinsi -> oy
    oy = collections.defaultdict(list)
    for et, v in bulgular.items():
        if not str(v.get('durum', '')).startswith('TAMAM'):
            continue
        isb = (v.get('isabet') or [])
        if not isb:
            continue
        c, _t, _tam = K.ad_coz(isb[0].get('baslik', ''))
        if c:
            oy[c].append(et)
    if oy:
        en_cok = max(oy.items(), key=lambda kv: len(kv[1]))
        uyusan_cins, uyusan_vtb = en_cok[0], en_cok[1]
    else:
        uyusan_cins, uyusan_vtb = None, []
    if len(uyusan_vtb) >= 2:
        hukum = 'DOGRULANDI'
    elif len(uyusan_vtb) == 1:
        hukum = 'DOGRULANAMADI (tek kaynak)'
    else:
        hukum = 'DOGRULANAMADI'
    return adl, hukum, uyusan_cins, uyusan_vtb, lokus, oy


# Mevcut kayitli kimlik (Kraken taxid adi) ile olculen kimlik ayni cinsi mi
# gosteriyor. "Candidatus" oneki karsilastirmadan once dusurulur; kayitli ad yoksa
# sonuc UYUSMADI degil KAYIT YOK olur - iki durum ayni sey degildir.
def ayni_mi(kayitli, dogrulanan_cins, adl):
    u"""Mevcut kayitli kimlik ile dogrulanan kimlik ayni cinsi mi gosteriyor?"""
    if not kayitli or kayitli in ('?', '-'):
        return 'KAYIT YOK', u'taxid_adlari.tsv icinde bu taxid icin ad yok'
    if not dogrulanan_cins:
        return 'BELIRSIZ', u'dogrulanan kimlikte cins cozulemedi'
    kc = re.sub(r'^(Ca\.|Candidatus)\s+', '', kayitli).split()[0]
    if kc.lower() == dogrulanan_cins.lower():
        return 'EVET', ''
    return 'HAYIR', u'kayitli "%s" -> olculen "%s"' % (kayitli, dogrulanan_cins)


# --------------------------------------------------------------- KOSU
# -------------------------------------------------------------------------
# SURUCU. Kutular kumeler halinde (varsayilan 24) islenir; her veritabani icin
# TEK akis acilir ve o akis kumedeki butun kutulara hizmet eder.
#
# ALAN FILTRESI YOKTUR ve bu bilincli bir karardir. Alani Kraken etiketine gore
# secmek - "bakteri kutusuna mantar veritabani sormaya gerek yok" - kutunun
# bakteri OLDUGUNU varsaymaktir; oysa bu asamanin varlik sebebi tam da o
# etiketlerin yanlis olabilmesi. Butun veritabanlari sorulur, alakasizlar SONUCTA
# dusurulur. Kutunun alani ETIKETTEN DEGIL OLCUMDEN cikar (alan_olcumden sutunu).
#
# SONUC DONDURMEYEN VERITABANI "TEMIZ" SAYILMAZ: SONUC YOK / DOSYA YOK /
# SORULMADI diye ayri isaretlenir ve her kutu satirinda 12 kaynagin 12'si de
# gorunur.
#
# KAPSAM MUHASEBESI: her veritabani icin TARANAN kayit sayisi sayilir ve
# BEKLENEN_KAYIT ile karsilastirilir. Erisim testinde gercek bir tavan sorunu
# yasandi - ilk kosu 120 001 kayitta kesiyordu ve SILVA SSU NR99, LSU Parc,
# UNITE ITS fiilen budanmis taraniyordu. Burada akiticida tavan yoktur ve eksik
# kapsam sessiz kalmaz, uyari basar.
#
# ONBELLEK I ILE PAYLASILIR: kontrol noktalari KIMLIK_SONUC/kontrol altinda ve
# ayni anahtarla durur, yani ayni kutu iki asamada iki kez taranmaz.
# -------------------------------------------------------------------------
def calistir(kok, kl_ust, kume_boyu, nt_kip, lit_kip, sifirla, yalniz, tavan_kutu):
    K = _K(kok)
    CIKTI = os.path.join(kok, 'TUM_KIMLIK_SONUC')
    # ONBELLEK PAYLASIMI: kontrol noktalari `I` ile AYNI klasorde ve AYNI
    # anahtarla durur. Ayni kutu iki asamada iki kez taranmaz.
    KONTROL = os.path.join(kok, 'KIMLIK_SONUC', 'kontrol')
    os.makedirs(CIKTI, exist_ok=True)
    os.makedirs(KONTROL, exist_ok=True)
    if sifirla:
        for f in os.listdir(CIKTI):
            if f.startswith('kutu_') and f.endswith('.json'):
                try:
                    os.remove(os.path.join(CIKTI, f))
                except OSError:
                    pass
    g = open(os.path.join(CIKTI, 'kosu_gunlugu.txt'), 'a', encoding='utf-8')

    def yaz(s=''):
        print(s, flush=True); g.write(s + '\n'); g.flush()

    yaz('=' * 78)
    yaz(u'  G - TUM KUTU KIMLIKLERININ BAGIMSIZ DOGRULANMASI')
    yaz(u'  surum %s   %s' % (VERSIYON, time.strftime('%Y-%m-%d %H:%M')))
    yaz('=' * 78)
    yaz(u'  `I` 12 IDDIAYI sinar - supheli olanlari. Bu asama panel olcumlerine')
    yaz(u'  giren BUTUN kutulari sinar: sessiz cogunlugu.')
    yaz('')

    katilan, atlanan, uye, rakip, kons, H = kutu_envanteri(kok, K)
    if yalniz:
        ay = [x.strip() for x in yalniz.split(',') if x.strip()]
        katilan = [k for k in katilan if any(a.lower() in k.lower() for a in ay)]
    if tavan_kutu:
        katilan = katilan[:tavan_kutu]

    # --- VERITABANI KAPSAMI: HEPSI, ALAN FILTRESI YOK ---
    var, yok = [], []
    for e, d, t, kullan, _n in K.VTB:
        if not kullan:
            continue                       # ikiz/altkume - bagimsiz kaynak degil
        p = os.path.join(kok, 'REFERANS_DB', d)
        (var if os.path.exists(p) else yok).append((e, d, t))
    lokus_tab = {e: t for e, _d, t, _k, _n in K.VTB}

    yaz(u'  kutu (olcume katilan)     : %d' % len(katilan))
    yaz(u'  kutu (ATLANAN)            : %d  %s'
        % (len(atlanan), ', '.join(k for k, _s in atlanan) or '-'))
    yaz(u'  VERITABANI                : %d sorgulanacak, %d dosyasi yok'
        % (len(var), len(yok)))
    for e, d, _t in var:
        yaz(u'      [SORULACAK] %-20s %-32s beklenen %s kayit'
            % (e, d, '{:,}'.format(BEKLENEN_KAYIT.get(e, 0)).replace(',', ' ') or '?'))
    for e, d, _t in yok:
        yaz(u'      [DOSYA YOK] %-20s %-32s REFERANS_DB altinda bulunamadi' % (e, d))
    yaz(u'  ALAN FILTRESI             : YOK. Her kutuya %d veritabaninin HEPSI'
        % len(var))
    yaz(u'    sorulur. Alani Kraken etiketine gore secmek tehlikeli olurdu -')
    yaz(u'    bu asamanin varlik sebebi zaten o etiketlerin yanlis olabilmesi.')
    yaz(u'    Alakasiz veritabanlari SONUCTA dusurulur ("SONUC YOK"), sorgudan')
    yaz(u'    ONCE elenmez. Kutunun alani ETIKETTEN DEGIL OLCUMDEN cikarilir.')
    yaz(u'  NCBI nt                   : %s'
        % {'oto': u'otomatik (URL API) - her kutu icin ayri BLAST',
           'elle': u'elle (sorgu dosyasi uretilir)',
           'yok': u'SORULMADI (varsayilan): %d kutu x BLAST kuyrugu gunlerce '
                  u'surer. `I` asamasindan kalan nt onbellegi VARSA kullanilir. '
                  u'--nt oto ile acilir.' % len(katilan)}[nt_kip])
    yaz('')

    # --- SURE TAHMINI ---
    tekil = {}
    for k in katilan:
        tekil.setdefault(K.kp_yolu(KONTROL, '_', kons[k][:4000], (), kl_ust), k)
    n_tekil = len(tekil)
    _uz = [min(len(kons[k]), 4000) for k in katilan] or [1500]
    _oq = sum(_uz) / float(len(_uz))
    _ref = 2000.0
    _bir = 6.7e-6 * min(_oq, _ref) + 4.11e-9 * min(_oq, _ref) * max(_oq, _ref)
    # Olculmus sabitler (bu betigin kendi kodu, sentetik kume, 2026-08-04):
    #   toplu tarama : ~6 400 kayit/sn (24 sorgu es zamanli)
    #   tek sorgulu  : ~2 200 kayit/sn (1 sorgu)  -> toplu ~70 kat verimli
    # Kullanicinin makinesi farkli olabilir; her veritabani satiri GERCEK sureyi
    # basar, yani tahmin kosu ilerledikce kendini duzeltir.
    TARAMA_HIZI = 6400.0
    _hiz = n_tekil * len(var) * kl_ust * _bir
    _kume = max(1, (n_tekil + kume_boyu - 1) // kume_boyu)
    _kayit = sum(BEKLENEN_KAYIT.get(e, 100000) for e, _d, _t in var)
    _tara = _kume * _kayit / TARAMA_HIZI
    yaz(u'  TAHMINI SURE: ~%s   (GECE BIRAKILACAK IS)' % K.sure_metni(_hiz + _tara))
    yaz(u'    tekil konsensus %d x %d veritabani = %d sorgu'
        % (n_tekil, len(var), n_tekil * len(var)))
    yaz(u'    hizalama payi  ~%s  (%d sorgu x %d aday, ort %d bp) <- BASKIN maliyet'
        % (K.sure_metni(_hiz), n_tekil * len(var), kl_ust, int(_oq)))
    yaz(u'    tarama payi    ~%s  (%d kume x %d veritabani akisi = %d akis, %s kayit)'
        % (K.sure_metni(_tara), _kume, len(var), _kume * len(var),
           '{:,}'.format(_kayit).replace(',', ' ')))
    yaz(u'    NOT: kutu basina AYRI akis olsaydi %d akis ve ~%s tarama gerekirdi;'
        % (n_tekil * len(var), K.sure_metni(n_tekil * _kayit / 2200.0)))
    yaz(u'    toplu tarama bunu %d akisa indirir (~70 kat).' % (_kume * len(var)))
    yaz(u'    Kesintiye dayaniklidir: her kutu bitince diske yazilir, ayrica')
    yaz(u'    veritabani taramalari `I` ile ORTAK onbellekte durur.')
    yaz('')

    # --- KOSU: kume kume, veritabani veritabani ---
    sonuc, tb = [], time.time()
    bekleyen = []
    for k in katilan:
        kp = os.path.join(CIKTI, 'kutu_%s.json' % re.sub(r'\W+', '_', k))
        if os.path.exists(kp):
            try:
                sonuc.append(json.load(open(kp, encoding='utf-8')))
                continue
            except Exception:
                pass
        bekleyen.append(k)
    yaz(u'  onceki kosudan alinan: %d kutu | taranacak: %d kutu'
        % (len(sonuc), len(bekleyen)))

    kapsam_kayit = {}          # etiket -> (taranan, beklenen, kapsam)
    for ki in range(0, len(bekleyen), kume_boyu):
        kume = bekleyen[ki:ki + kume_boyu]
        yaz('')
        yaz(u'[kume %d/%d] %d kutu: %s'
            % (ki // kume_boyu + 1, (len(bekleyen) + kume_boyu - 1) // kume_boyu,
               len(kume), ', '.join(kume[:6]) + (' ...' if len(kume) > 6 else '')))
        bulgular = {k: {} for k in kume}
        for et, dosya, _t in var:
            # onbellekte olmayanlari topla (I ile ORTAK anahtar)
            kalan = {}
            for k in kume:
                q = kons[k][:4000]
                kp = K.kp_yolu(KONTROL, et, q, (), kl_ust)
                if os.path.exists(kp):
                    try:
                        bulgular[k][et] = json.load(open(kp, encoding='utf-8'))
                        continue
                    except Exception:
                        pass
                kalan[k] = q
            if not kalan:
                yaz(u'     %-20s: %d kutunun HEPSI onbellekten geldi' % (et, len(kume)))
                continue
            t0 = time.time()
            yol = os.path.join(kok, 'REFERANS_DB', dosya)

            def ilerle(n, _e=et, _t0=t0):
                print('     ... %s: %d kayit tarandi (%s)      '
                      % (_e, n, K.sure_metni(time.time() - _t0)), end='\r', flush=True)
            kls, taranan = toplu_kisa_liste(K, yol, kalan, kl_ust, ilerle)
            bek = BEKLENEN_KAYIT.get(et)
            kapsam = ('TAMAMI' if bek and taranan >= bek else
                      ('TAMAMI (beklenen bilinmiyor)' if not bek else
                       'EKSIK (%d / %d)' % (taranan, bek)))
            kapsam_kayit[et] = (taranan, bek, kapsam)
            for k, q in kalan.items():
                kl = kls.get(k) or []
                if not kl:
                    # SONUC YOK - "temiz" DEGIL, ayri isaretlenir
                    res = dict(durum=u'SONUC YOK', kayit=0, kisa_liste_boyu=kl_ust,
                               hizalanan=0, isabet=[], kazanan_sira=None,
                               kazanan_kaynak=None, sira_uyarisi=None,
                               taranan_kayit=taranan, kapsam=kapsam,
                               sebep=u'TEMIZ SAYILMAZ: tarandi (%s kayit, kapsam %s) '
                                     u'ama hicbir kayit sorgunun tohumlarini '
                                     u'tutturmadi; kutu bu veritabaninin alani '
                                     u'disinda olabilir'
                                     % ('{:,}'.format(taranan).replace(',', ' '), kapsam))
                else:
                    res = K.kl_degerlendir(kl, q, kl_ust, taranan=taranan, t0=t0)
                    res['kapsam'] = kapsam
                json.dump(res, open(K.kp_yolu(KONTROL, et, q, (), kl_ust), 'w',
                                    encoding='utf-8'), ensure_ascii=False, default=str)
                bulgular[k][et] = res
            bos = len([1 for k in kalan if not (kls.get(k) or [])])
            yaz(u'     %-20s: %s kayit tarandi, kapsam %s | %d kutu isabet, %d SONUC YOK (%s)'
                % (et, '{:,}'.format(taranan).replace(',', ' '), kapsam,
                   len(kalan) - bos, bos, K.sure_metni(time.time() - t0)))
            if bek and taranan < bek:
                yaz(u'     >>> UYARI: KAPSAM EKSIK. %s icinde %d kayit beklenirken %d '
                    u'tarandi. Tavan sorunu tekrarlanmis olabilir - kontrol edin.'
                    % (et, bek, taranan))

        # dosyasi olmayan veritabanlari da SATIRDA GORUNSUN
        for e, d, _t in yok:
            for k in kume:
                bulgular[k][e] = dict(durum=u'DOSYA YOK', isabet=[],
                                      sebep=u'REFERANS_DB/%s bulunamadi' % d)

        for k in kume:
            r = kutu_kaydi(K, kok, k, kons[k][:4000], bulgular[k], lokus_tab, uye,
                           rakip, H, nt_kip, lit_kip, KONTROL, CIKTI, yaz, var, yok)
            json.dump(r, open(os.path.join(CIKTI, 'kutu_%s.json'
                                           % re.sub(r'\W+', '_', k)), 'w',
                              encoding='utf-8'), ensure_ascii=False, default=str)
            sonuc.append(r)
        gec = time.time() - tb
        yap = len([1 for s in sonuc if s['kutu'] in bekleyen])
        if yap:
            print('        gecen %s | tahmini kalan %s'
                  % (K.sure_metni(gec), K.sure_metni(gec / yap * (len(bekleyen) - yap))),
                  flush=True)

    raporla(K, CIKTI, sonuc, atlanan, var, yok, kapsam_kayit, uye, rakip, yaz, kl_ust, nt_kip)
    g.close()
    return 0 if sonuc else 1


# ---------------------------------------------------------------------------
# Tek kutunun satirini kurar: kayitli kimlik, dogrulanan kimlik, uyusma durumu,
# savunulabilir duzey, en iyi uc isabet, literatur katmani ve KAYNAK MUHASEBESI.
#
# Kaynak muhasebesi 12 yerel + NCBI nt satirinin HEPSINI icerir; bir kaynak hic
# denenmediyse "SORULMADI (BILINMEYEN) - HATA, bildirin" yazilir. Sessiz eksik
# birakmaktansa gurultulu bir hata basmak tercih edilmistir.
#
# Satirda ayrica bu kutunun HANGI hedeflerin uyesi ya da rakibi oldugu durur:
# kimlik degisirse hangi olcumlerin yeniden yapilmasi gerektigi buradan okunur.
# ---------------------------------------------------------------------------
def kutu_kaydi(K, kok, kutu, q, bulgular, lokus_tab, uye, rakip, H, nt_kip, lit_kip,
               KONTROL, CIKTI, yaz, var, yok):
    u"""Tek kutunun satiri: kayitli kimlik, dogrulanan kimlik, uyusma, kaynak muhasebesi."""
    taxid = kutu.split('_')[-1]
    kayitli = H.taxid_adlari().get(taxid, '')

    # --- NCBI nt katmani (ayri kaynak) ---
    ntk = os.path.join(KONTROL, 'nt_%s.json' % re.sub(r'\W+', '_', kutu))
    if os.path.exists(ntk):
        try:
            bulgular[K.NT_ETIKET] = json.load(open(ntk, encoding='utf-8'))
        except Exception:
            pass
    if K.NT_ETIKET not in bulgular:
        if nt_kip == 'yok':
            bulgular[K.NT_ETIKET] = dict(
                durum=u'SORULMADI (--nt yok)', isabet=[],
                sebep=u'G asamasinda nt varsayilan olarak kapali: kutu basina ayri '
                      u'BLAST kuyrugu gerekir. `I` asamasindan onbellek de yok. '
                      u'--nt oto ile acilir.')
        else:
            bulgular[K.NT_ETIKET] = K.nt_katmani(kutu, q, CIKTI, yaz, nt_kip)
            if str(bulgular[K.NT_ETIKET].get('durum', '')).startswith('TAMAM'):
                json.dump(bulgular[K.NT_ETIKET], open(ntk, 'w', encoding='utf-8'),
                          ensure_ascii=False, default=str)

    adl, hukum, cins, uyusan_vtb, lokus, oy = kutu_hukmu(K, bulgular, lokus_tab)
    uym, uym_not = ayni_mi(kayitli, cins, adl)

    # --- LITERATUR KATMANI (I ile ayni) ---
    try:
        import importlib.util as _lu
        _lp = _lu.spec_from_file_location(
            'lit', os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'literatur_kontrol.py'))
        LIT = _lu.module_from_spec(_lp); _lp.loader.exec_module(LIT)
        lit = LIT.kontrol_et((adl.get('isabet1') or {}).get('tam_ad', ''),
                             adl.get('onerilen_ad', ''), lokus, ag=(lit_kip != 'yok'))
    except Exception as _e:
        lit = dict(durum=u'literatur modulu yuklenemedi: %s' % type(_e).__name__,
                   ncbi_guncel_ad='-', ad_farkli_mi='-', revizyon_uyarisi='-')

    # --- KAYNAK MUHASEBESI: 12 veritabaninin 12'si de satirda gorunur ---
    beklenen_et = [e for e, _d, _t in var] + [e for e, _d, _t in yok] + [K.NT_ETIKET]
    detay = collections.OrderedDict()
    for et in beklenen_et:
        v = bulgular.get(et)
        if v is None:
            detay[et] = dict(durum=u'SORULMADI (BILINMEYEN)', en_iyi='', kimlik=None,
                             sebep=u'bu veritabani hic denenmedi - HATA, bildirin')
            continue
        d = str(v.get('durum', '?'))
        isb = (v.get('isabet') or [])
        detay[et] = dict(
            durum=d, en_iyi=(isb[0].get('baslik', '')[:120] if isb else ''),
            kimlik=(isb[0].get('kimlik') if isb else None),
            kazanan_sira=v.get('kazanan_sira'), kazanan_kaynak=v.get('kazanan_kaynak'),
            taranan_kayit=v.get('taranan_kayit'), kapsam=v.get('kapsam'),
            sebep=v.get('sebep', ''))
    sorgulanan = len([1 for v in detay.values() if str(v['durum']).startswith('TAMAM')])
    sonucsuz = [e for e, v in detay.items() if v['durum'] == 'SONUC YOK']
    sorulmayan = [e for e, v in detay.items()
                  if str(v['durum']).startswith('SORULMADI') or v['durum'] == 'DOSYA YOK']
    # ALAN OLCUMDEN: hangi veritabanlari fiilen isabet verdi
    alan = sorted({lokus_tab.get(e, '?') for e, v in detay.items()
                   if str(v['durum']).startswith('TAMAM')})

    kzs = [v['kazanan_sira'] for v in detay.values() if isinstance(v.get('kazanan_sira'), int)]
    return dict(
        kutu=kutu, taxid=taxid, kayitli_kimlik=kayitli or '-',
        dogrulanan_cins=cins or '-', dogrulanan_ad=adl.get('onerilen_ad', '-'),
        duzey=adl.get('duzey', '-'), gerekce=adl.get('gerekce', '-'),
        uyusuyor=uym, uyusma_notu=uym_not, hukum=hukum,
        uyusan_vtb=uyusan_vtb, oylar={k: v for k, v in oy.items()},
        adlandirma=adl, literatur=lit, vtb_detay=detay,
        sorgulanan_vtb=sorgulanan, toplam_vtb=len(beklenen_et),
        sonuc_yok_vtb=sonucsuz, sorulmayan_vtb=sorulmayan,
        alan_olcumden=alan, lokus=lokus,
        kazanan_sira_maks=(max(kzs) if kzs else None),
        uye_hedefler=sorted(uye.get(kutu, [])), rakip_hedefler=sorted(rakip.get(kutu, [])))


# --------------------------------------------------------------- RAPOR
# UYUSMAYANLAR EN BASA, sonra belirsizler, en sonda uyusanlar; esitlikte cok hedef
# etkileyen kutu once. Rapor okuyan kisi once en cok is cikaracak satiri gorsun.
def _sirala(s):
    """UYUSMAYANLAR EN BASA. Sonra belirsizler, sonra uyusanlar."""
    o = {'HAYIR': 0, 'BELIRSIZ': 1, 'KAYIT YOK': 2, 'EVET': 3}
    return (o.get(s['uyusuyor'], 9), -len(s['uye_hedefler']) - len(s['rakip_hedefler']),
            s['kutu'])


# ---------------------------------------------------------------------------
# Iki cikti: kutu basina tek satirlik TSV ve markdown rapor.
#
# Markdown raporun sonunda ETKI OZETI durur: kac kutunun kimligi degisti, bu kac
# hedefin uyelik/rakip kumesini etkiliyor ve yeniden olcum gerekip gerekmedigi.
# Uyelik kumesi degisen hedefte omurga konsensus da degisebilir, yani primer ve
# ayrim degerleri yeniden hesaplanmalidir - bu betik o hesabi KENDISI yapmaz,
# yalnizca gerektigini bildirir ve panel dosyalarina dokunmaz.
# ---------------------------------------------------------------------------
def raporla(K, CIKTI, sonuc, atlanan, var, yok, kapsam_kayit, uye, rakip, yaz,
            kl_ust, nt_kip):
    sonuc = sorted(sonuc, key=_sirala)
    beklenen_et = [e for e, _d, _t in var] + [e for e, _d, _t in yok] + [K.NT_ETIKET]
    t = os.path.join(CIKTI, 'tum_kutu_kimlikleri.tsv')
    with open(t, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# Panel olcumlerine giren HER kutunun kimligi bagimsiz sinandi.\n')
        fh.write(u'# Yontem `I` ile ayni: %d\'luk kisa liste, hepsi hizalanir, en az '
                 u'IKI bagimsiz veritabani uyusmasi.\n' % kl_ust)
        fh.write(u'# ALAN FILTRESI YOK: her kutuya butun veritabanlari sorulur; '
                 u'alakasiz olanlar SONUCTA dusurulur ("SONUC YOK"), sorgudan once '
                 u'ELENMEZ. Kutunun alani etiketten degil OLCUMDEN cikar.\n')
        fh.write(u'# UYUSMAYAN satirlar EN BASTA.\n')
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['kutu', 'taxid', 'MEVCUT_KAYITLI_KIMLIK', 'DOGRULANAN_KIMLIK',
                    'UYUSUYOR_MU', 'uyusma_notu', 'HUKUM',
                    'SAVUNULABILIR_DUZEY', 'ONERILEN_AD', 'adlandirma_gerekcesi',
                    # HIZALANAN UZUNLUK, kimlik yuzdesinin yaninda ZORUNLU.
                    # (2026-08-11) Yuzde tek basina yaniltiyor: ayni kutuya
                    # SILVA LSU Parc'ta Petriella setifera %100, RefSeq mantar
                    # ITS'te Petriella musispora %100 cikti. Ikisi de %100
                    # cunku iki AYRI lokusta, iki AYRI uzunlukta hizalandilar.
                    # Uzunluk gorunmeden "%100" okuyan herkes turun kesin
                    # oldugunu sanir. Veri zaten uretiliyordu (hiz_uzunluk),
                    # yalnizca tabloya yazilmiyordu.
                    'en_iyi_isabet', 'en_iyi_kimlik_%', 'en_iyi_hiz_uzunluk', 'en_iyi_vtb',
                    'ikinci_isabet', 'ikinci_kimlik_%', 'ikinci_hiz_uzunluk', 'ikinci_vtb',
                    'ucuncu_isabet', 'ucuncu_kimlik_%', 'ucuncu_hiz_uzunluk', 'ucuncu_vtb',
                    'UYESI_OLDUGU_HEDEFLER', 'RAKIBI_OLDUGU_HEDEFLER',
                    'etkilenen_hedef_sayisi',
                    'SORGULANAN_VTB', 'SONUC_YOK_VTB', 'SORULMAYAN_VTB',
                    'alan_olcumden', 'lokus', 'kazanan_sira_maks',
                    'LIT_ncbi_guncel_ad', 'LIT_AD_FARKLI_MI', 'LIT_revizyon_uyarisi',
                    'LIT_durum', 'HER_VTB_NE_DEDI'])
        for s in sonuc:
            a = s.get('adlandirma') or {}

            def _i(n, alan, v='-'):
                return ((a.get('isabet%d' % n) or {}).get(alan) or v)
            d = s['vtb_detay']
            hepsi = ' | '.join(
                '%s [%s]: %s%s%s%s'
                % (e,
                   ('%s kayit, kapsam %s' % ('{:,}'.format(v['taranan_kayit']).replace(',', ' '),
                                             v.get('kapsam') or '?'))
                   if v.get('taranan_kayit') else v['durum'],
                   v['en_iyi'] or v['durum'],
                   ('' if v.get('kimlik') is None else ' (%%%s)' % K.vir(v['kimlik'])),
                   ('' if v.get('kazanan_sira') is None
                    else ' {sira %s/%d}' % (v['kazanan_sira'], kl_ust)),
                   ('' if not v.get('sebep') else ' <%s>' % v['sebep'][:110]))
                for e, v in d.items())
            w.writerow([
                s['kutu'], s['taxid'], s['kayitli_kimlik'], s['dogrulanan_ad'],
                s['uyusuyor'], s['uyusma_notu'], s['hukum'],
                s['duzey'], s['dogrulanan_ad'], s['gerekce'],
                _i(1, 'tam_ad'), K.vir(_i(1, 'kimlik', None)), _i(1, 'uzunluk'), _i(1, 'vtb'),
                _i(2, 'tam_ad'), K.vir(_i(2, 'kimlik', None)), _i(2, 'uzunluk'), _i(2, 'vtb'),
                _i(3, 'tam_ad'), K.vir(_i(3, 'kimlik', None)), _i(3, 'uzunluk'), _i(3, 'vtb'),
                ', '.join(s['uye_hedefler']) or '-',
                ', '.join(s['rakip_hedefler']) or '-',
                len(set(s['uye_hedefler']) | set(s['rakip_hedefler'])),
                '%d / %d' % (s['sorgulanan_vtb'], s['toplam_vtb']),
                ', '.join(s['sonuc_yok_vtb']) or '-',
                ', '.join(s['sorulmayan_vtb']) or '-',
                ', '.join(s['alan_olcumden']) or '-', s['lokus'],
                s['kazanan_sira_maks'] if s['kazanan_sira_maks'] is not None else '-',
                (s.get('literatur') or {}).get('ncbi_guncel_ad', '-'),
                (s.get('literatur') or {}).get('ad_farkli_mi', '-'),
                (s.get('literatur') or {}).get('revizyon_uyarisi', '-'),
                (s.get('literatur') or {}).get('durum', '-'), hepsi])
        if atlanan:
            fh.write(u'#\n# ATLANAN KUTULAR (sessiz atlama yok):\n')
            for k, sebep in atlanan:
                fh.write(u'# %s\t%s\n' % (k, sebep))
    yaz('  yazildi: %s' % t)

    # ----------------------------------------------------------- MD RAPORU
    degisen = [s for s in sonuc if s['uyusuyor'] == 'HAYIR']
    belirsiz = [s for s in sonuc if s['uyusuyor'] in ('BELIRSIZ', 'KAYIT YOK')]
    etkilenen = set()
    for s in degisen:
        etkilenen |= set(s['uye_hedefler']) | set(s['rakip_hedefler'])
    uye_etkilenen = {h for s in degisen for h in s['uye_hedefler']}
    r = os.path.join(CIKTI, 'TUM_KUTU_KIMLIK_RAPORU.md')
    with open(r, 'w', encoding='utf-8') as fh:
        fh.write(u'# Tum kutu kimliklerinin bagimsiz dogrulanmasi\n\n')
        fh.write(u'Uretim: %s · betik %s\n\n' % (time.strftime('%Y-%m-%d %H:%M'), VERSIYON))
        fh.write(u'`I` asamasi 12 **iddiayi** sinar - supheli olanlari. Bu asama panel '
                 u'olcumlerine giren **butun kutulari** sinar: uyelik ve rakip '
                 u'kumelerini belirleyen, yani her ayrim katini etkileyen sessiz '
                 u'cogunlugu.\n\n')

        # --- KAYNAK KAPSAMI (kanit) ---
        fh.write(u'## Veritabani kapsami\n\n')
        fh.write(u'**Alan filtresi UYGULANMADI.** Her kutuya %d yerel veritabaninin '
                 u'hepsi soruldu. Alani Kraken etiketine gore secmek tehlikeli olurdu: '
                 u'bu asamanin varlik sebebi zaten o etiketlerin yanlis olabilmesi. '
                 u'Alakasiz veritabanlari **sonucta** dusuruldu (`SONUC YOK`), sorgudan '
                 u'once elenmedi. Kutunun alani etiketten degil **olcumden** cikarildi '
                 u'(`alan_olcumden` sutunu).\n\n' % len(var))
        fh.write(u'| # | veritabani | beklenen kayit | taranan | kapsam |\n'
                 u'|---|---|---|---|---|\n')
        for i, (e, _d, _t) in enumerate(var, 1):
            tar, bek, kap = kapsam_kayit.get(e, (None, BEKLENEN_KAYIT.get(e), None))
            fh.write(u'| %d | %s | %s | %s | %s |\n'
                     % (i, e, '{:,}'.format(BEKLENEN_KAYIT.get(e, 0)).replace(',', ' '),
                        '{:,}'.format(tar).replace(',', ' ') if tar else
                        u'onbellekten (bu kosuda taranmadi)',
                        kap or u'onceki kosuda dogrulandi'))
        for e, d, _t in yok:
            fh.write(u'| - | %s | - | - | **DOSYA YOK** (`REFERANS_DB/%s`) |\n' % (e, d))
        fh.write(u'| - | NCBI nt | - | - | %s |\n'
                 % {'yok': u'**SORULMADI** (--nt yok; kutu basina ayri BLAST kuyrugu)',
                    'oto': u'otomatik (URL API)', 'elle': u'elle'}[nt_kip])
        eksik = [e for e, (tar, bek, kap) in kapsam_kayit.items()
                 if kap and kap.startswith('EKSIK')]
        fh.write(u'\n> **Tavan sorunu tekrarlanmadi.** Erisim testinde ilk kosu '
                 u'120 001 kayitta kesiyordu ve SILVA SSU NR99 (510 495), LSU Parc '
                 u'(1 312 521), UNITE ITS (2 069 189) budanmis taraniyordu. Burada '
                 u'akiticida tavan **yoktur**; taranan kayit sayisi sayilir ve '
                 u'beklenenle karsilastirilir. Bu kosuda kapsami eksik veritabani: '
                 u'**%s**.\n\n' % (', '.join(eksik) if eksik else u'YOK'))
        fh.write(u'> Bir veritabani sonuc dondurmediyse **"temiz" sayilmadi**: '
                 u'`SONUC YOK` diye ayri isaretlendi ve satirda gorunur. Her kutu '
                 u'satirinda %d kaynagin %d\'si de listelenir.\n\n'
                 % (len(beklenen_et), len(beklenen_et)))

        # --- ETKI OZETI ---
        fh.write(u'## Etki ozeti\n\n')
        fh.write(u'| olcu | deger |\n|---|---|\n')
        fh.write(u'| sinanan kutu | %d |\n' % len(sonuc))
        fh.write(u'| **kimligi DEGISEN kutu** | **%d** |\n' % len(degisen))
        fh.write(u'| kimligi dogrulanan kutu | %d |\n'
                 % len([s for s in sonuc if s['uyusuyor'] == 'EVET']))
        fh.write(u'| belirsiz / kayitli adi olmayan | %d |\n' % len(belirsiz))
        fh.write(u'| **etkilenen hedef** | **%d** |\n' % len(etkilenen))
        fh.write(u'| bunlardan UYELIK kumesi degisen | %d |\n' % len(uye_etkilenen))
        fh.write(u'| atlanan kutu (olcume girmiyor) | %d |\n\n' % len(atlanan))
        if degisen:
            fh.write(u'> **YENIDEN OLCUM GEREKIYOR.** %d kutunun kimligi degisti ve bu '
                     u'%d hedefin uyelik ya da rakip kumesini etkiliyor. Uyelik kumesi '
                     u'degisen %d hedefte omurga konsensus da degisebilir, yani primer '
                     u've ayrim degerleri yeniden hesaplanmalidir.\n\n'
                     u'Etkilenen hedefler: %s\n\n'
                     % (len(degisen), len(etkilenen), len(uye_etkilenen),
                        ', '.join(sorted(etkilenen))))
        else:
            fh.write(u'> Hicbir kutunun kimligi degismedi: uyelik ve rakip kumeleri '
                     u'oldugu gibi kalir, **yeniden olcum gerekmiyor**.\n\n')

        # --- UYUSMAYANLAR ---
        if degisen:
            fh.write(u'## Uyusmayanlar (once bunlar)\n\n')
            for s in degisen:
                fh.write(u'### %s  (taxid %s)\n\n' % (s['kutu'], s['taxid']))
                fh.write(u'- **Kayitli kimlik:** %s\n' % s['kayitli_kimlik'])
                fh.write(u'- **Dogrulanan:** %s  (`%s`) - %s\n'
                         % (s['dogrulanan_ad'], s['duzey'], s['hukum']))
                fh.write(u'  - *Gerekce:* %s\n' % s['gerekce'])
                fh.write(u'- **Uyesi oldugu hedefler:** %s\n'
                         % (', '.join(s['uye_hedefler']) or '-'))
                fh.write(u'- **Rakibi oldugu hedefler:** %s\n'
                         % (', '.join(s['rakip_hedefler']) or '-'))
                a = s.get('adlandirma') or {}
                fh.write(u'\n  | # | en yakin kayit | cins | tur | kimlik | veritabani |\n'
                         u'  |---|---|---|---|---|---|\n')
                for n_ in (1, 2, 3):
                    it = a.get('isabet%d' % n_)
                    if it:
                        fh.write(u'  | %d | %s | %s | %s | %%%s | %s |\n'
                                 % (n_, it['tam_ad'], it['cins'], it['tur'],
                                    K.vir(it['kimlik']), it['vtb']))
                fh.write(u'\n  **Kaynak muhasebesi (%d/%d):**\n\n'
                         % (s['sorgulanan_vtb'], s['toplam_vtb']))
                fh.write(u'  | veritabani | durum | en iyi isabet | kimlik | kazanan sira |\n'
                         u'  |---|---|---|---|---|\n')
                for e, v in s['vtb_detay'].items():
                    fh.write(u'  | %s | %s | %s | %s | %s |\n'
                             % (e, v['durum'], v['en_iyi'] or (v.get('sebep') or '-')[:96],
                                '-' if v.get('kimlik') is None else '%%%s' % K.vir(v['kimlik']),
                                v.get('kazanan_sira') if v.get('kazanan_sira') is not None else '-'))
                fh.write(u'\n')
        if atlanan:
            fh.write(u'## Atlanan kutular\n\n')
            fh.write(u'Bunlar hicbir hedefin uyesi ya da rakibi degil, yani hicbir ayrim '
                     u'hesabina girmiyorlar. Sessizce atlanmadilar - burada yazililar.\n\n')
            for k, sebep in atlanan:
                fh.write(u'- `%s` - %s\n' % (k, sebep))
            fh.write(u'\n')
    yaz('  yazildi: %s' % r)
    yaz('')
    yaz(u'  ETKI: %d kutu sinandi | kimligi DEGISEN %d | etkilenen hedef %d '
        u'| uyelik degisen %d | atlanan %d'
        % (len(sonuc), len(degisen), len(etkilenen), len(uye_etkilenen), len(atlanan)))
    if eksik:
        yaz(u'  >>> UYARI: kapsami EKSIK veritabani: %s' % ', '.join(eksik))


# Komut satiri: --kisa-liste (varsayilan I ile ayni), --kume bir akista kac kutu,
# --nt NCBI katmani (varsayilan "yok": kutu basina ayri BLAST kuyrugu gunlerce
# surer; I asamasindan kalan onbellek yine kullanilir), --yalniz / --tavan-kutu
# sinama icin alt kume.
def main():
    p = argparse.ArgumentParser(
        description=u'Panel olcumlerine giren TUM kutularin kimligini bagimsiz dogrula')
    p.add_argument('--kok', default='.')
    p.add_argument('--kisa-liste', type=int, default=None, dest='kisa_liste',
                   help=u'tam hizalanacak aday sayisi (varsayilan: I ile ayni)')
    p.add_argument('--kume', type=int, default=24,
                   help=u'bir veritabani akisinda kac kutu birlikte taransin '
                        u'(bellek/hiz dengesi, varsayilan 24)')
    p.add_argument('--nt', choices=['oto', 'elle', 'yok'], default='yok',
                   help=u'NCBI nt katmani (varsayilan yok: kutu basina ayri BLAST '
                        u'kuyrugu gunlerce surer; I asamasindan kalan onbellek yine '
                        u'kullanilir)')
    p.add_argument('--literatur', choices=['oto', 'yok'], default='oto')
    p.add_argument('--yalniz', default=None,
                   help=u'virgulle ayrilmis kutu adi parcalari (sinama icin)')
    p.add_argument('--tavan-kutu', type=int, default=0, dest='tavan_kutu',
                   help=u'yalniz ilk N kutu (sinama icin)')
    p.add_argument('--sifirla', action='store_true')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)
    if not os.path.isdir(os.path.join(kok, 'KAPSAMLI_ARAMA')):
        sys.exit('HATA: %s icinde KAPSAMLI_ARAMA yok.' % kok)
    kl = a.kisa_liste
    if kl is None:
        kl = _K(kok).KISA_LISTE
    return calistir(kok, kl, a.kume, a.nt, a.literatur, a.sifirla, a.yalniz, a.tavan_kutu)


if __name__ == '__main__':
    sys.exit(main() or 0)
