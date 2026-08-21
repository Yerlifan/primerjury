# -*- coding: utf-8 -*-
"""Specificity verification — four independent evidence layers before ordering.

WHY A SEPARATE STAGE
    A pair that looks clean in your sample can still fail in the lab: a sample
    is 99 bins, not the world. A pair may face no competitor in the sample and
    still amplify something abundant that the sample never captured.

THE FOUR LAYERS, AND WHY THEY ARE BUILT DIFFERENTLY
    1  in-sample    in-silico PCR against the raw reads.
                    Independent of every reference database.
    2  local DB     scan of the local reference sets. Uses OUR code and OUR
                    engine — so it shares our bugs, and cannot corroborate
                    itself. That is exactly why layers 3 and 4 exist.
    3  MFEprimer    an external binary. Different implementation, different
                    thermodynamic model, written by other people.
    4  NCBI         Primer-BLAST against nt. Independent of our choice of
                    databases, which is layer 2's main blind spot.

    If the layers disagree the row is marked CELISKILI and is NOT orderable.
    Disagreement is treated as information, not noise: a contradiction means at
    least one measurement is wrong, and shipping either would be a gamble.

THREE STATES, NOT TWO
    BILINMIYOR (unknown) is distinct from TEMIZ (clean). A layer that did not
    run never votes in favour. This is the single most common way pipelines of
    this kind mislead: an unmeasured check silently reads as a passed check.

    Related: the in-sample column is DISPLAYED but does not VOTE. It is the
    admission criterion that put the pair on this list, so by construction it
    can never say RISKLI. Letting a constant vote made 16 of 16 rows come out
    CELISKILI mechanically (fixed 2026-08-06).

--- ozgun aciklama ---
DOGRULAMA TURU - kurtarilan ciftleri siparise gondermeden once uc bagimsiz
kanit katmaniyla sinar.

NEDEN AYRI BIR SECENEK
----------------------
Numunede iyi gorunen bir cift laboratuvarda tutmayabilir: numune 99 kutudan
ibarettir, dunya degildir. Bir cift numunede hic rakip gormeyebilir ama
veritabaninda binlerce hedef disi organizmaya baglanabiliyor olabilir.

Bu betik verification turunun (secenek K) esigi gecirdigi HER YENI ya da DEGISMIS
cift icin iki kanit katmani daha kosar ve ucunu YAN YANA koyar:

  1) NUMUNE OLCUMU        - kurtarma turundan gelir (zaten olculdu)
  2) YEREL VERITABANI      - REFERANS_DB altindaki kumeler; mevcut kuresel
                             tarama kodu (screening/kuresel_tarama.py)
                             AYNEN kullanilir, yeniden yazilmaz
  3) NCBI                  - otomatik (URL API) ya da elle (Primer-BLAST)

UC SONUC AYRILIYORSA satir CELISKILI isaretlenir ve SIPARIS EDILEBILIR SAYILMAZ.
Celiskiler bu turun EN DEGERLI ciktisidir; raporun en basinda dururlar.

Panel dosyalarina YAZMAZ. Yalniz okur, DOGRULAMA_SONUC/ altina yazar.
"""

# -------------------------------------------------------------------------
# dogrulama_turu.py — kurtarma turunun esigi gecirdigi YENI/DEGISMIS ciftleri,
# siparise gonderilmeden once dort bagimsiz kanit katmaniyla yan yana sinar.
#
# GİRDİ  : KURTARMA_SONUC/kurtarma_satirlari.tsv (hangi ciftler dogrulanacak),
#          TEK_PROTOKOL_SONUC/panel_tek_protokol.tsv (primer dizileri),
#          REFERANS_DB/ altindaki KUMELER listesi, tools/mfeprimer indeksleri,
#          NCBI Primer-BLAST (ag) ya da elle doldurulmus NCBI_SONUC_SABLONU.tsv.
# ÇIKTI  : DOGRULAMA_SONUC/dogrulama_uc_sutun.tsv (asil tablo),
#          DOGRULAMA_SONUC/CELISKILER.md (once bu okunur),
#          DOGRULAMA_SONUC/yerel_vuruslar.tsv, DOGRULAMA_RAPORU.md,
#          NCBI_PRIMER_BLAST_GIRDI.tsv + NCBI_SONUC_SABLONU.tsv (elle yol),
#          DOGRULAMA_SONUC/kontrol/ (kume basina kontrol noktalari).
# ÇAĞRAN : screening.bat -> D tusu
#          (bat icinde: wsl -e python3 "verification/dogrulama_turu.py" --kok . ...)
#
# NEDEN DORT KATMAN: numunede iyi gorunen bir cift laboratuvarda tutmayabilir,
# cunku numune 99 kutudan ibarettir, dunya degildir. Katman 1 ve 2 BIZIM
# kodumuzdur ve ayni motoru kullanir - o motorda hata varsa ikisi de ayni yonde
# yanilir. Katman 3 (MFEprimer) disaridan gelen bagimsiz bir aractir, katman 4
# (NCBI) ise bizim veritabani secimimizden bagimsiz bir kaynaktir. Katmanlar
# ayrilirsa satir CELISKILI isaretlenir ve SIPARIS EDILEMEZ.
# -------------------------------------------------------------------------
import os, sys, csv, json, time, argparse, re

VERSIYON = '1.0 (2026-08-03)'

URUN_ALT, URUN_UST = 70, 400          # yalanci urun aranan mesafe araligi

# VURUS_ESIGI = 0  -> TEK bir hedef disi urun bile katmani RISKLI yapar.
#
# A3 (2026-08-21): bu deger gerekcesiz duruyordu; projedeki her sabit gerekcesini
# tasidigi icin istisnaydi. Gerekce su ve BILINCLI bir katiliktir:
#   * Panelin kabul kurali "siparise tek bir yanlis primer gitmesin"dir. Yanlis
#     siparis hem pahali hem yavastir; bir yalanci urunu kacirmanin bedeli, bir
#     adayi haksiz yere elemenin bedelinden buyuktur.
#   * Esik TEK BASINA hukum vermez. RISKLI bir katman yalnizca bir OY'dur;
#     karar dort katmanin UYUSMASINA baglanir (bkz. birlestir()). Yani sifir
#     esik, tek basina siparisi engellemez, ikinci bir kaynagin onaylamasini
#     sart kosar.
#   * Sayim zaten bir SUZGECTEN gecmis urunlerdir: URUN_ALT..URUN_UST araliginda
#     ve baglanma kuralini saglayan urunler. Rastgele gurultu buraya girmez.
#
# NE ZAMAN DEGISTIRILIR: kesif kipinde (hangi aday hic sansi var sorusu) 1-2'ye
# cikarilabilir. Siparis oncesi hukumde YUKSELTILMEMELIDIR. Degistirirseniz
# raporda hangi esikle kosuldugu yazilir - sessizce degismez.
VURUS_ESIGI = int(os.environ.get('PT_VURUS_ESIGI', '0'))
NCBI_SONUC_TAVANI = 1000               # Primer-BLAST sonuc sayfasinin urun tavani.
                                       # n==tavan ise deger SAYIM DEGIL. D-3.
# D-12 (2026-08-07): panelin gercek baglanma sicakligi. Bir hedef disi
# amplikonun "olusabilir" sayilip sayilmayacagi bununla karsilastirilir.
# Kaynak: panelin kendi Ta'si (panel karari, 2026-08-07). MFEprimer her amplikon icin
# KENDI Ta'sini yazar; o deger amplikonun GC'sinden turetilir ve bizim
# termosiklerimizde kullanacagimiz sicaklik DEGILDIR.
TA_PANEL = 57.9
BOY_TOL = 10                           # beklenen urun boyuna bu kadar yakin vurus
                                       # HEDEFIN KENDI urunudur (MFEprimer katmani
                                       # ile ayni tolerans). D-1 duzeltmesi.

# Taranacak yerel kumeler: (etiket, dosya adi, aciklama)
KUMELER = [
    ('SILVA SSU NR99', 'SILVA_138.2_SSURef_NR99.fasta', u'510 495 kayit; SSU (16S/18S)'),
    ('SILVA LSU NR99', 'SILVA_138.2_LSURef_NR99.fasta', u'95 279 kayit; LSU (23S/28S)'),
    ('UNITE ITS', 'UNITE_ITS.fasta', u'2 069 189 kayit; mantar ITS'),
    ('PR2 SSU', 'PR2_SSU_taxo_long.fasta', u'240 201 kayit; okaryot 18S'),
    ('ROD operon', 'ROD_v1.2_operon_variants.fasta', u'60 320 kayit; rRNA operon varyantlari'),
    ('RefSeq bakteri 16S', 'bacteria.16S.fna', u'26 877 kayit'),
    ('RefSeq arke 16S', 'archaea.16S.fna', u'1 160 kayit'),
    ('RefSeq mantar ITS', 'fungi.ITS.fna', u'20 394 kayit'),
    ('RefSeq mantar 28S', 'fungi.28SrRNA.fna', u'12 890 kayit'),
    ('RefSeq mantar 18S', 'fungi.18SrRNA.fna', u'4 037 kayit'),
    ('RefSeq ref_all2', 'ref_all2.fna', u'65 358 kayit; RefSeq birlesik'),
]
# OZGULLUK taramasinda SILVA Parc BILEREK kullanilmaz: 1,3 milyon kayitlik
# tekrarsizlastirilmamis kume, yalanci urun riski sorusuna NR99'un uzerine
# yeni bilgi katmaz ama kosuyu saatlerce uzatir. KIMLIK sorusu farklidir -
# orada Parc SARTTIR (bkz. kimlik_dogrulama.py) cunku NR99 nadir cinsleri siler.
# Istenirse acilir:
PARC_ISTEGE_BAGLI = ('SILVA LSU Parc', 'SILVA_138.2_LSUParc.fasta',
                     u'1 312 521 kayit; tekrarsizlastirilmamis')

OLCUT_NOTU = u"""
YEREL TARAMANIN OLCUTU - ACIKCA YAZILIYOR
=========================================
Kullanilan kod: screening/kuresel_tarama.py  (AYNEN, degistirilmedi)

  * iki primer KARSILIKLI YONELIMDE baglanacak
  * aralarindaki mesafe %d-%d bp
  * F ve R uyumsuzluklari TOPLAM en cok %s

3' SON IKI BAZ SARTI BU KATMANDA UYGULANMADI. Sebep: mevcut kuresel tarama
kodu bu sarti tasimiyor ve o kodu yeniden yazmamak icin oldugu gibi kullanildi.
Bu olcut, 3' son iki baz sarti olan olcutten DAHA GEVSEKTIR - yani bulunan
vuruslarin bir kismi gercekte urun VERMEYEBILIR. Risk taramasinda bu GUVENLI
taraftir: gercek riski gozden kacirmaktansa fazladan uyari uretir.
Bir vurus ciddiye alinacaksa 3' ucunun tuttugu ayrica bakilmalidir.
"""


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
        return ('%.*f' % (b, float(str(x).replace(',', '.')))).replace('.', ',')
    except (TypeError, ValueError):
        return str(x)


# ---------------------------------------------------------------- girdi
_ATLANAN = []


# ---------------------------------------------------------------------------
# Dogrulanacak cift kumesini kurar. Iki tur satir vardir:
#   YENI CIFT   - kurtarma turu yeni bir F/R buldu; diziler satirin metninden
#                 dogrudan cikarilir.
#   DEGISMIS    - primerler ayni, degisen sey OLCU ya da UYELIK; diziler
#                 P asamasinin panel tablosundan alinir.
# Primer dizisi hicbir kaynaktan bulunamayan satir SESSIZCE atlanmaz, _ATLANAN
# listesine yazilir ve kosu basliginda gosterilir.
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# SIPARIS LISTESI KIPI (--siparis) - 2026-08-06
#
# NEDEN EKLENDI: bu betik varsayilan olarak yalnizca KURTARILAN ciftleri
# dogruluyordu (kurtarma_satirlari.tsv). Son kosuda o liste tek satirdi, yani
# siparise giden 15 ciftin 14'u Primer-BLAST katmanini HIC gormedi. Siparis
# oncesi sorulan soru "listedeki her cift yalanci urun riski tasiyor mu"
# oldugu icin girdi kumesi siparis listesi olmalidir.
#
# Ayni dort katman kosar (numune / yerel DB / MFEprimer / NCBI Primer-BLAST);
# degisen tek sey HANGI ciftlerin sinandigi.
# ---------------------------------------------------------------------------
def siparistekiler(kok, hepsi=False):
    """SIPARIS_LISTESI.tsv -> siparise giden ciftler (KESIN + EVRENSEL).

    hepsi=True ise KOSULLU ve ONERILMEZ satirlar da alinir - kullanici onlari
    da sinatmak isterse. Varsayilan: yalniz siparise gidecekler.
    """
    yol = os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'SIPARIS_LISTESI.tsv')
    if not os.path.exists(yol):
        yol = os.path.join(kok, 'SIPARIS_LISTESI.tsv')
    if not os.path.exists(yol):
        sys.exit('HATA: SIPARIS_LISTESI.tsv yok.\n'
                 '      Once screening.bat -> secenek (T) kosulmalidir.')
    with open(yol, encoding='utf-8') as fh:
        satirlar = list(csv.DictReader(
            (s for s in fh if not s.startswith('#')), delimiter='\t'))
    kabul = ('KESIN', 'EVRENSEL') if not hepsi else ('KESIN', 'EVRENSEL',
                                                     'KOSULLU', 'ONERILMEZ')
    out = []
    for s in satirlar:
        sn = (s.get('SINIF') or '').strip().upper()
        if sn not in kabul:
            continue
        F, R = (s.get('F') or '').strip(), (s.get('R') or '').strip()
        if not F or not R:
            _ATLANAN.append(s.get('hedef', '?'))
            continue
        out.append(dict(hedef=s['hedef'], F=F, R=R,
                        urun=s.get('urun_bp', ''),
                        sinif=sn,
                        # SIPARISTE_MI: KESIN/EVRENSEL satirlar siparise gider.
                        # --ncbi-yalniz-siparis bu bayragi kullanir (2026-08-07).
                        sipariste=(sn in ('KESIN', 'EVRENSEL')),
                        tur=u'SIPARIS LISTESI (%s)' % sn,
                        numune_deger=s.get('ayrim_mm1', ''),
                        numune_olcu=u'dCq %s' % (s.get('dCq_karsiligi') or '-'),
                        yol=s.get('siparis_sarti', '')))
    return out, yol


def kurtarilanlar(kok):
    """KURTARMA_SONUC/kurtarma_satirlari.tsv -> esigi gecen YENI/DEGISMIS ciftler."""
    yol = os.path.join(kok, 'KURTARMA_SONUC', 'kurtarma_satirlari.tsv')
    if not os.path.exists(yol):
        sys.exit('HATA: %s yok.\n      Once screening.bat -> secenek (K) kosulmalidir.' % yol)
    with open(yol, encoding='utf-8') as fh:
        satirlar = list(csv.DictReader((s for s in fh if not s.startswith('#')), delimiter='\t'))

    tp = os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'panel_tek_protokol.tsv')
    ciftler = {}
    if os.path.exists(tp):
        with open(tp, encoding='utf-8') as fh:
            for r in csv.DictReader((s for s in fh if not s.startswith('#')), delimiter='\t'):
                ciftler[r['hedef']] = (r.get('F', ''), r.get('R', ''), r.get('urun_bp', ''))

    out = []
    for s in satirlar:
        if not (s.get('esigi_gecti_mi') or '').startswith('EVET'):
            continue
        yeni = s.get('yeni_deger') or ''
        m = re.search(r'YENI CIFT\s+([ACGT]+)\s*/\s*([ACGT]+)\s*\((\d+)\s*bp\)', yeni)
        if m:
            F, R, bp = m.group(1), m.group(2), m.group(3)
            tur = 'YENI CIFT'
        else:
            F, R, bp = ciftler.get(s['hedef'], ('', '', ''))
            tur = 'DEGISMIS (ayni primerler, olcu/uyelik degisti)'
        if not F or not R:
            _ATLANAN.append(s['hedef'])
            continue
        out.append(dict(hedef=s['hedef'], F=F, R=R, urun=bp, tur=tur,
                        numune_deger=s.get('yeni_deger') or s.get('eski_deger'),
                        numune_olcu=s.get('olcu', ''), yol=s.get('denenen_yol', '')))
    return out, yol


# ---------------------------------------------------------------- katman 1
# ---------------------------------------------------------------------------
# KATMAN 2 - yerel veritabani taramasi. Aranan sey: hedef disi bir organizmada
# iki primerin karsilikli yonelimde ve 70-400 bp mesafede baglandigi bir yer var
# mi. Mevcut kuresel tarama kodu (screening/kuresel_tarama.py) AYNEN
# kullanilir; ikinci bir tarayici yazmak iki farkli olcut demek olurdu.
#
# Bu katmanin olcutu 3' son iki baz sarti TASIMAZ, yani gercek PCR olcutunden
# DAHA GEVSEKtir: bulunan vuruslarin bir kismi gercekte urun vermeyebilir. Risk
# taramasinda guvenli taraf budur - gercek riski kacirmaktansa fazladan uyari.
#
# SILVA Parc BILEREK disaridadir: 1,3 milyon kayitlik tekrarsizlastirilmamis kume
# yalanci urun sorusuna NR99'un ustune bilgi katmaz ama kosuyu saatlerce uzatir.
# KIMLIK sorusunda durum tersidir (orada Parc SARTTIR) - iki soru ayni degildir.
# ---------------------------------------------------------------------------
def katman1_yerel(kok, ciftler, yaz, kontrol_dizin, parc=False, kume_ust=0):
    """MEVCUT kuresel tarama kodunu kullanir. Her kume icin ayri kontrol noktasi."""
    sys.path.insert(0, kok)
    from screening import kuresel_tarama as KT, yapilandirma as C

    adaylar = [dict(ad=c['hedef'], F=c['F'], R=c['R'], lo=URUN_ALT, hi=URUN_UST)
               for c in ciftler]
    # 2026-08-06 HATA DUZELTMESI (D-1): eskiden yalniz 'urun' (TOPLAM vurus)
    # tutuluyordu ve o sayi '2_hedef_disi_urun' sutununa yaziliyordu. YANLISTI:
    # kuresel_tarama HEDEFIN KENDI UYELERINI de sayar (Methanosarcina cifti
    # SILVA'da 485 Methanosarcina dizisi bulur; bunlar hedef DISI degil, hedefin
    # TA KENDISI). MFEprimer katmani bu ayrimi zaten 'ayni_boyda' / 'hedef_disi'
    # olarak yapiyordu; iki katman AYNI soruyu sormadigi icin birlestir() surekli
    # CELISKILI uretiyordu. Artik yerel katman da beklenen urun boyunu (+-BOY_TOL)
    # ayri sayar ve hukme yalniz 'hedef_disi' girer.
    bek = {c['hedef']: int(c['urun']) for c in ciftler
           if str(c.get('urun', '')).strip().isdigit()}

    # ----------------------------------------------------------------- A2
    # TAKSONOMIK AYRIM (2026-08-21). Bu katman eskiden hedef ICI / hedef DISI
    # ayrimini yalniz BOYA bakarak yapiyordu; katman 3 (MFEprimer, D-12) ayni
    # soruyu TAKSONA bakarak cevapliyordu. D-12'de olculdu: "hedef disi" sayilan
    # 1.605 amplikonun %95,7'si hedef kladin KENDI ICINDENdi, yalniz boyu
    # farkliydi. Iki katman AYNI soruyu sormadikca birlestir() yapisal olarak
    # celiski uretir - D-1 yorumunun uyardigi sey tam da budur.
    #
    # Siniflandirma TARAMA SIRASINDA yapilir ve SAYILIR; kimlikler saklanmaz.
    # Olculdu: Bakteri_universal 483.098 vurus veriyor, kimlikleri saklamak
    # ~100 MB eder. Sayac aday basina sabit bellek tutar ve sayi EKSIKSIZ olur.
    try:
        sys.path.insert(0, kok)
        from screening import taksonomi as TX
        import verification.mfe_katmani as _MK
        _klad = _MK.klad_tablosu(kok)
    except Exception as e:
        TX, _klad = None, {}
        yaz(u'  UYARI: taksonomik ayrim YAPILAMIYOR (%s). Yalniz boy olcutu '
            u'kullanilacak - bu katman katman 3 ile AYNI soruyu sormayacak.' % e)

    _klad_yok = [c['hedef'] for c in ciftler if c['hedef'] not in _klad]
    if _klad and _klad_yok:
        yaz(u'  UYARI: hedef_klad.tsv\'de tanimi olmayan %d hedef - bunlarda '
            u'taksonomik ayrim YAPILMAZ: %s'
            % (len(_klad_yok), ', '.join(_klad_yok)[:160]))

    def _siniflandirici(aday_ad, baslik, db_ad):
        if TX is None or aday_ad not in _klad:
            return 'bilinmiyor'
        alan, jetonlar, _kaynak = _klad[aday_ad]
        return TX.sinifla(baslik, db_ad, jetonlar, alan)

    toplam = {c['hedef']: dict(urun=0, ayni_boyda=0, hedef_disi=0, kume={},
                               boy={}, vurus=[], tarandi=0, atlanan=0,
                               sinif={k: 0 for k in ('a', 'ao', 'b', 'c', 'bilinmiyor')},
                               siniflandirildi=False)
              for c in ciftler}
    kumeler = list(KUMELER) + ([PARC_ISTEGE_BAGLI] if parc else [])
    if kume_ust:
        # yalniz en kucuk N kume (hizli sinama icin - kapsam degil CALISMA kaniti)
        var = [(e, d, a) for e, d, a in kumeler
               if os.path.exists(os.path.join(kok, 'REFERANS_DB', d))]
        var.sort(key=lambda t: os.path.getsize(os.path.join(kok, 'REFERANS_DB', t[1])))
        kumeler = var[:kume_ust]
        yaz(u'  (kume-ust=%d: yalniz %s taraniyor - CALISMA kaniti, kapsam degil)'
            % (kume_ust, ', '.join(e for e, _, _ in kumeler)))
    for etiket, dosya, aciklama in kumeler:
        db = os.path.join(kok, 'REFERANS_DB', dosya)
        if not os.path.exists(db):
            yaz(u'  [%s] ATLANDI - dosya yok: %s' % (etiket, dosya))
            for h in toplam:
                toplam[h]['kume'][etiket] = 'dosya yok'
                toplam[h]['atlanan'] += 1
            continue
        # K-7: pickle'in icinde aday bazli sonuclar var; aday kumesi her kosuda
        # degisiyor. Imza olmadan ikinci gece KeyError ile oluyordu.
        import hashlib
        # 2026-08-10 DIZI MUHRU: imza yalnizca aday ADLARINDAN uretiliyordu,
        # dizi degisince ayni imza cikip eski tarama geri okunuyordu.
        imza = hashlib.md5(
            '|'.join(sorted('%s>%s<%s' % (a['ad'], a.get('F', ''), a.get('R', ''))
                            for a in adaylar)).encode('utf-8')).hexdigest()[:10]
        dy = os.path.join(kontrol_dizin, 'yerel_%s_%s.pkl'
                          % (''.join(ch if ch.isalnum() else '_' for ch in etiket), imza))
        t0 = time.time()
        yaz(u'  [%s] taraniyor (%s)...' % (etiket, aciklama))

        def ilerle(pi, kayit, gecen):
            print('     ... parca %d, %d kayit (%s)          '
                  % (pi, kayit, sure_metni(gecen)), end='\r', flush=True)
        try:
            res = KT.tara(adaylar, db=db, durum_yolu=dy, ilerle=ilerle,
                          siniflandirici=_siniflandirici)
        except TypeError:
            # Eski imzali kuresel_tarama (siniflandirici parametresi yok).
            # Bu SESSIZ bir dususe donusmemeli: taksonomik ayrim yapilmadigi
            # raporda acikca gorunur ('2_klad_ayrimi' sutunu HAYIR olur).
            yaz(u'  UYARI: kuresel_tarama eski imzali - taksonomik ayrim YOK')
            res = KT.tara(adaylar, db=db, durum_yolu=dy)
        for h, r in res.items():
            if r.get('hata'):
                toplam[h]['kume'][etiket] = r['hata']
                continue
            toplam[h]['kume'][etiket] = r.get('urun', 0)
            toplam[h]['urun'] += r.get('urun', 0)
            toplam[h]['tarandi'] += 1
            # D-1: boy histogramini biriktir ve beklenen boy +-BOY_TOL icinde
            # kalanlari 'ayni_boyda' say. Bu, MFEprimer katmaninin kullandigi
            # olcutun AYNISI (mfe_katmani._spec_ayristir, tolerans=10).
            _b = bek.get(h)
            for _sz, _n2 in (r.get('boy') or {}).items():
                try:
                    _sz = int(_sz)
                except (TypeError, ValueError):
                    continue
                toplam[h]['boy'][_sz] = toplam[h]['boy'].get(_sz, 0) + _n2
                if _b is not None and abs(_sz - _b) <= BOY_TOL:
                    toplam[h]['ayni_boyda'] += _n2
                else:
                    toplam[h]['hedef_disi'] += _n2
            if _b is None:
                # beklenen boy bilinmiyorsa ayrim YAPILAMAZ - hukum bilinmiyor olsun
                toplam[h]['boy_ayrimi_yok'] = True
            # A2: taksonomik sayaclari veritabanlari boyunca biriktir
            _s = r.get('sinif') or {}
            for _k, _v in _s.items():
                if _k in toplam[h]['sinif']:
                    toplam[h]['sinif'][_k] += _v
            if r.get('siniflandirildi'):
                toplam[h]['siniflandirildi'] = True
            for v in (r.get('vurus') or [])[:20]:
                toplam[h]['vurus'].append((etiket,) + tuple(v))
        yaz(u'     bitti (%s): %s' % (sure_metni(time.time() - t0),
                                      ', '.join('%s=%s' % (h, toplam[h]['kume'][etiket])
                                                for h in list(toplam)[:4])))
    return toplam


# ---------------------------------------------------------------- katman 2
PB_URL = 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi'

# NCBI'ye ardisik gonderim arasi asgari bekleme (sn). Primer-BLAST'a hizli
# ardisik is gonderimi IP engellemesine yol acar. Eski kod is anahtari
# ALINAMADIGINDA hic beklemeden bir sonraki cifte geciyordu; 16 cift saniyeler
# icinde pespese gonderiliyordu. Artik her gonderim arasi beklenir.
PB_GONDERIM_ARASI = 10

# ---------------------------------------------------------------------------
# D-8 HATA DUZELTMESI (2026-08-07) - IKI AYRI HATA, IKISI DE SIPARIS ONCESI
# NCBI KATMANINI TAMAMEN ISLEVSIZ BIRAKIYORDU.
#
# HATA 1 - "is anahtari alinamadi" (16/16 cift).
#   Gonderilen deger: ORGANISM='Bacteria (taxid:2) OR Archaea (taxid:2157) OR
#   Fungi (taxid:4751)'. NCBI'nin HAM yaniti (tahmin degil, okundu):
#     "Exception error: Invalid organism or taxonomy id input: 2 OR Archaea .
#      Please check the spelling and make sure it is on the suggested organism
#      list in organism input field"
#   Yani ORGANISM alani TEK organizma alir; 'OR' sozdizimi YOKTUR. Primer-BLAST
#   sayfasinin KENDI javascript'i (js/primerInit.js) coklu organizmayi soyle
#   yapar:
#     function AddOneOrgField(e, orgName, orgVal) { ... name=\"ORGANISM\" ... }
#     function AddOrgField(e) { AddOneOrgField(e,"ORGANISM"); ... }
#     function GetOrganismURL(){ jQuery(".multiOrg").each(function(){
#         url += "&ORGANISM=" + $(this).value; }); }
#   Yeni alanlarin name'i de "ORGANISM"dir - yani coklu organizma TEKRARLANAN
#   ORGANISM alanidir. Python'da dict ile bu YAPILAMAZ (tek anahtar); ikili
#   listesi + urlencode(liste) gerekir. Duzeltme budur.
#
# HATA 2 - "BOS SONUC" (onceki kosuda 9/15 cift; organizma kisitindan BAGIMSIZ).
#   Sonuc sayfasi "No target templates were found in selected database" diyordu
#   ve hicbir urun listelenmiyordu. Sebep: istekte Primer3 alanlarinin cogunu
#   HIC gondermiyorduk. Gondermeyince CGI o alanlari BASLATILMAMIS bellekten
#   okuyor. Ham sonuc sayfasindaki "Search Summary" tablosu bunu acikca
#   gosteriyordu:
#         Opt Primer size      1086305756
#         Min Tm               4.94733e-316
#         Opt Tm               2.18186e+243
#         Max Tm               0            <-- OLDURUCU OLAN BU
#         Max Tm difference    6.95299e-310
#   "Max Tm = 0" ile hicbir urun olcute giremez, sayfa bos doner ve bu bizde
#   "temiz" degil ama "veri yok" olarak isaretlenirdi - yani katman hicbir sey
#   olcmezdi. Duzeltme: NCBI'nin KENDI formundaki (primer-blast/ sayfasi,
#   defVal ozniteligi) varsayilan degerlerin TAMAMI gonderilir. Ayni istek
#   duzeltmeden sonra ayni cift icin 1000 urun satiri dondurdu ve Search
#   Summary'de Max Tm=75, Min Tm=45, Opt Primer size=20 yaziyordu.
#
# NOT: primer boyu/Tm sinirlari kasten GENISLETILMISTIR. Bu turda primer
# TASARLAMIYORUZ, ELDEKI oligolari sinatiyoruz; Primer3'un tasarim suzgeci
# bizim sabit oligolarimizi elemesin diye sinirlar genis tutulur.
# ---------------------------------------------------------------------------
PB_VARSAYILAN = {
    # --- NCBI formunun kendi varsayilanlari (defVal ozniteliginden alindi) ---
    'PRIMER_NUM_RETURN': '10', 'PRIMER_MAX_DIFF_TM': '20',
    'PRIMER_ON_SPLICE_SITE': '0',
    'SPLICE_SITE_OVERLAP_5END': '7', 'SPLICE_SITE_OVERLAP_3END': '4',
    'SPLICE_SITE_OVERLAP_3END_MAX': '8',
    'MIN_INTRON_SIZE': '1000', 'MAX_INTRON_SIZE': '1000000',
    'SEARCHMODE': '0', 'MAX_TARGET_SIZE': '4000',
    # D-13 (2026-08-07, OLCULDU): 1000 bir SINIR DEGIL, formun VARSAYILANI.
    # NCBI form kaynagi: <input name="NUM_TARGETS_WITH_PRIMERS" defVal="1000">
    # ("Max targets to show (for pre-designed primers)"). Istemci tarafinda ust
    # sinir dogrulamasi YOK. Ayni cift ayni veritabaniyla kosuldugunda:
    #     gonderilen 3000  -> donen 3000 satir
    #     gonderilen 8000  -> donen 8000 satir
    #     gonderilen 20000 -> donen 11999 satir (sinir baglamadi, hedef tukendi)
    # Ucu birlikte kirpiyor: NUM_TARGETS_WITH_PRIMERS, MAX_TARGET_PER_TEMPLATE,
    # HITSIZE. Ucu de yukseltildi. Sayfa kesildiginde HICBIR uyari basmadigi
    # icin eski 1000 "tavan" gibi gorunuyordu.
    'NUM_TARGETS': '20', 'NUM_TARGETS_WITH_PRIMERS': '20000',
    'MAX_TARGET_PER_TEMPLATE': '1000',
    'TOTAL_MISMATCH_IGNORE': '6', 'HITSIZE': '100000', 'EVALUE': '30000',
    'WORD_SIZE': '7', 'MAX_CANDIDATE_PRIMER': '500',
    'PRIMER_MIN_GC': '10.0', 'PRIMER_MAX_GC': '90.0',
    'GC_CLAMP': '0', 'POLYX': '5',
    'PRIMER_MAX_END_STABILITY': '9', 'PRIMER_MAX_END_GC': '5',
    'PRIMER_MAX_TEMPLATE_MISPRIMING_TH': '40.00',
    'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH': '70.00',
    'PRIMER_MAX_SELF_ANY_TH': '45.0', 'PRIMER_MAX_SELF_END_TH': '35.0',
    'PRIMER_PAIR_MAX_COMPL_ANY_TH': '45.0', 'PRIMER_PAIR_MAX_COMPL_END_TH': '35.0',
    'PRIMER_MAX_HAIRPIN_TH': '24.0',
    'PRIMER_MAX_TEMPLATE_MISPRIMING': '12.00',
    'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING': '24.00',
    'SELF_ANY': '8.00', 'SELF_END': '3.00',
    'PRIMER_PAIR_MAX_COMPL_ANY': '8.00', 'PRIMER_PAIR_MAX_COMPL_END': '3.00',
    'OVERLAP_5END': '7', 'OVERLAP_3END': '4',
    'MONO_CATIONS': '50.0', 'DIVA_CATIONS': '1.5',
    'CON_DNTPS': '0.6', 'CON_ANEAL_OLIGO': '50.0',
    'SALT_FORMULAR': '1', 'TM_METHOD': '1',
    'PRIMER_MISPRIMING_LIBRARY': 'AUTO',
    'ALLOW_NO_ORGANISM': 'NO', 'UNGAPPED_BLAST': 'on',
    'LOW_COMPLEXITY_FILTER': 'on', 'SHOW_SVIEWER': 'on',
    'SEARCH_SPECIFIC_PRIMER': 'on',
    # --- sabit oligo SINAMA modu icin kasten genisletilenler ---
    'PRIMER_MIN_SIZE': '15', 'PRIMER_OPT_SIZE': '20', 'PRIMER_MAX_SIZE': '30',
    'PRIMER_MIN_TM': '45.0', 'PRIMER_OPT_TM': '60.0', 'PRIMER_MAX_TM': '75.0',
    # --- bos birakilanlar (formda da bos) ---
    'PRIMER5_START': '', 'PRIMER5_END': '', 'PRIMER3_START': '', 'PRIMER3_END': '',
    'ENTREZ_QUERY': '', 'PRODUCT_MIN_TM': '', 'PRODUCT_OPT_TM': '', 'PRODUCT_MAX_TM': '',
    'CMD': 'request',
}


def pb_ac(url, veri=None, deneme=4, timeout=90, yaz=None):
    """Primer-BLAST'a istek. Gecici ag/hiz-siniri hatalarinda yeniden dener.

    NCBI yogunlukta baglantiyi yanit vermeden kapatabilir
    (RemoteDisconnected) ya da 429/502/503 dondurebilir. Eski kod ilk
    hatada o cifti BASARISIZ yazip geciyordu; sinamada tam bu oldu.
    Bekleme her denemede iki katina cikar (10, 20, 40 sn).
    """
    import urllib.request, urllib.error
    son = None
    for i in range(deneme):
        try:
            req = urllib.request.Request(
                url, veri, headers={'User-Agent': 'PrimerJury-primer-QC/1.0'})
            with urllib.request.urlopen(req, timeout=timeout) as f:
                return f.read().decode('utf-8', 'replace')
        except Exception as e:
            son = e
            kod = getattr(e, 'code', None)
            # 400/404 gibi KALICI hatalarda yeniden denemenin anlami yok
            if kod is not None and kod not in (429, 500, 502, 503, 504):
                raise
            if i == deneme - 1:
                break
            b = PB_GONDERIM_ARASI * (2 ** i)
            if yaz:
                yaz(u'    ag hatasi (%s) - %d sn sonra yeniden deneniyor (%d/%d)'
                    % (type(e).__name__, b, i + 2, deneme))
            time.sleep(b)
    raise son


def organizma_listesi(s):
    """'--organizma' dizgesini AYRI organizmalara boler.

    Primer-BLAST'in ORGANISM alani TEK organizma alir (bkz. yukaridaki D-8
    notu). Kullanicinin yazdigi 'A OR B OR C' bicimi kabul edilir ama ISTEGE
    o bicimde GONDERILMEZ - burada parcalanip TEKRARLANAN ORGANISM alanlarina
    donusturulur. ';' ve yeni satir da ayirac sayilir.
    """
    if not s or not s.strip():
        return []
    parcalar = re.split(r'\s+OR\s+|;|\n', s, flags=re.I)
    return [x.strip() for x in parcalar if x.strip()]


# KATMAN 4 - NCBI. blastn -remote KULLANILMAZ: 45 saniyelik surec tavanini asiyor
# ve is yarida kesiliyordu. Yerine Primer-BLAST'in URL API'si kullanilir: is
# gonderilir, is anahtari alinir, bitene kadar yoklanir. Her yanit ham olarak
# diske yazilir ki sayilar tartisildiginda kaynak metne bakilabilsin.

# D-13b (2026-08-07): Primer-BLAST sonuc sayfasinda urunler
# <div class="prPairTl">BASLIK</div> ... bloklarina bolunur. Butun sayfayi
# saymak "hedefteki urun" ile "hedef disi urun"u ayni kefeye koyar.
# 2026-08-10 OLCULDU (canli tek cift denemesi, Proteolitik_Cloacimonas):
# Hedefin kendi taksonu ENTREZ_QUERY ile dislandiginda bile Primer-BLAST
# "potentially unintended templates" bolumunu ACMIYOR - sablon dizi
# bildirilmedigi surece bulunan her urunu "target templates" altina koyuyor.
# Sayfadaki "Products on potentially unintended templates" yazisi bir BOLUM
# BASLIGI degil, sayfanin ust kismindaki BAG LISTESIDIR; bos bolum icin de
# basiliyor. Bu yuzden bolum sayimina dayali hukum YAPISAL OLARAK calismaz.
#
# Olculen sayilar (Proteolitik_Cloacimonas, txid112 dislanmis):
#   42 urun -> 40'i "uncultured bacterium clone ...", 2'si "Methanogenic
#   prokaryote enrichment culture ..." -> ADLI takson SIFIR.
# Adsiz cevre klonlari hicbir taksona bagli olmadigi icin ENTREZ_QUERY onlari
# eleyemez; ustelik hedeflerimiz zaten adlandirilmamis soylar oldugundan bu
# klonlarin buyuk kismi HEDEFIN KENDISI olabilir. Etiketten karar verilemez.
#
# Bu yuzden hukum artik boluma degil BASLIGA bakar: kimligi ADLI olan urun
# hedef disi kanitidir, adsiz cevre klonu ise "karar veremez" hanesine yazilir
# ve kimligine dizi karsilastirmasi (katman 2-3) karar verir.
_ADSIZ_IZLERI = (u'uncultured', u'unidentified', u'unclassified', u'metagenome',
                 u'environmental sample', u'enrichment culture', u'clone',
                 u'synthetic construct')


def _ncbi_urunleri(html):
    """Sonuc sayfasindaki her urunu (erisim_no, baslik, urun_boyu) olarak dondurur."""
    d = re.sub(r'<[^>]+>', ' ', html)
    d = re.sub(r'&nbsp;?', ' ', d)
    d = re.sub(r'\s+', ' ', d)
    return [(m.group(1), m.group(2).strip(), int(m.group(3)))
            for m in re.finditer(
                r'>\s*([A-Z]{1,2}[_A-Z]*\d{5,}\.\d)\s+([^>]{5,200}?)\s+'
                r'product length\s*=\s*(\d+)', d)]


def _adsiz_mi(baslik):
    """Kayit ADSIZ mi. TEK TANIM: ncbi_yeniden_siniflandir.adli_mi().

    2026-08-10 gece duzeltmesi. Ilk kural yalniz anahtar kelime ariyordu
    (uncultured, clone, metagenome...) ve su basliklari ADLI sayiyordu:
        "Bacterium LC2012 16S ribosomal RNA gene"
        "Archaeon 2022-TM-MRBT1 gene for 16S rRNA"
        "anaerobic methanogenic archaeon E15-5 16S rRNA gene"
        "Environmental 16s rDNA sequence from Evry wastewater treatment plant"
    Hicbirinin cins adi yok. Bu yuzden hedef disi sayilari sisiyordu:
    Bacteroidales 650 -> 82, Nitrosocosmicus 170 -> 9, Methanothrix cinsi
    22 -> 1. Yeni kural anahtar kelime aramaz, AD arar.
    """
    try:
        from ncbi_yeniden_siniflandir import adli_mi as _adli
    except ImportError:
        b = baslik.lower()
        return any(iz in b for iz in _ADSIZ_IZLERI)
    return not _adli(baslik)


def _ncbi_bolum_say(html, baslik):
    """Verilen bolum basligindan sonraki 'product length' satirlarini sayar."""
    par = re.split(r'<div class="prPairTl">(.*?)</div>', html)
    t = 0
    for i in range(1, len(par), 2):
        b = re.sub(r'<[^>]+>', '', par[i]).strip().lower()
        if baslik.lower() in b:
            t += len(re.findall(r'product length\s*=\s*\d+', par[i + 1], re.I))
    return t


def katman2_oto(ciftler, cikti, yaz, organizma='', bekleme=20, tur_ust=60,
                haric_taxid=''):
    """NCBI Primer-BLAST URL API. blastn -remote KULLANILMAZ (45 sn tavanini asar).

    Gonder -> is anahtari al -> bitene kadar yokla. Her yanit ham olarak diske
    yazilir; agsiz makinede ya da NCBI kuyruklu oldugunda DUZGUNCE vazgecer ve
    elle yola dusulmesini soyler.
    """
    import urllib.request, urllib.parse
    ham = os.path.join(cikti, 'ncbi_ham')
    os.makedirs(ham, exist_ok=True)
    orgs = organizma_listesi(organizma)
    if orgs:
        yaz(u'  organizma kisiti (%d ayri ORGANISM alani olarak gonderilecek): %s'
            % (len(orgs), ' | '.join(orgs)))
    else:
        yaz(u'  organizma kisiti YOK - tum nt taranacak (genis hedeflerde sonuc '
            u'tavanina carpma olasiligi yuksek)')
    # D-13c (2026-08-07, OLCULDU): ciplak 'NOT txidN[Organism]' filtreyi TERSINE
    # CEVIRIR - o taksonu DISLAMAK yerine YALNIZ onu getirir. Olculmus kanit
    # (ayni cift, ayni veritabani, cins dagilimi):
    #   ENTREZ_QUERY bos                      -> Escherichia 45, Bacillus 34,
    #                                            Pseudomonas 32, Staphylococcus 25
    #   'NOT txid1279[Organism]'              -> YALNIZ Staphylococcus 196  (TERS!)
    #   'all[filter] NOT txid1279[Organism]'  -> Escherichia 14, Bacillus 13,
    #                                            ..., Staphylococcus 1  (DOGRU)
    # Bu yuzden onek KODUN ICINDE zorunlu kilinir; kullanici ciplak NOT yazsa da
    # duzeltilir.
    # 2026-08-10 HEDEF BAZLI DISLAMA. Tek ve genel bir taxid butun kosuya
    # uygulaniyordu; oysa her hedefin KENDI taksonu dislanmali. Genel organizma
    # kisiti (Bacteria OR Archaea OR Fungi) ile calisildiginda Primer-BLAST
    # bulunan her urunu "target templates" hanesine koyuyor ve "unintended"
    # bolumu bos kaliyor - 22 sayfanin 22'sinde olculdu. Cozum: hedefin kendi
    # taksonunu ENTREZ_QUERY ile disla, o zaman kalan her urun tanimi geregi
    # hedef disidir. Harita dosyasi: screening/hedef_taxid.tsv
    def _ent_of(_tx):
        # Birden cok taxid virgulle verilebilir; her biri AYRI bir NOT terimi
        # olur. Evrensel primerlerde (ornek: Metanojen_universal) hedef tek bir
        # takson degil, birkac takimin birlesimidir - tek terim yetmez.
        _ler = [x.strip().lstrip('txid') for x in str(_tx).split(',') if x.strip()]
        if not _ler:
            return ''
        return 'all[filter]' + ''.join(' NOT txid%s[Organism]' % x for x in _ler)

    HARITA = {}
    _hy = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'screening', 'hedef_taxid.tsv')
    if os.path.exists(_hy):
        for _l in open(_hy, encoding='utf-8'):
            _l = _l.rstrip('\n')
            if not _l.strip() or _l.startswith('#'):
                continue
            _p = _l.split('\t')
            if len(_p) >= 2 and _p[1].strip():
                HARITA[_p[0].strip()] = _p[1].strip()
        yaz(u'  hedef bazli dislama haritasi okundu: %d hedef (%s)'
            % (len(HARITA), os.path.basename(_hy)))
    ent = ''
    if haric_taxid:
        ent = _ent_of(haric_taxid)
        yaz(u'  ENTREZ_QUERY (GENEL dislama): %s' % ent)
    if HARITA:
        yaz(u'    NOT: "all[filter]" oneki ZORUNLU - onsuz NCBI filtreyi tersine '
            u'cevirip yalniz o taksonu getirir (olculdu).')
    out = {}
    ilk = True
    for c in ciftler:
        ad = c['hedef']
        if not ilk:
            time.sleep(PB_GONDERIM_ARASI)      # NCBI hiz siniri - bkz. PB_GONDERIM_ARASI
        ilk = False
        p = dict(PB_VARSAYILAN)
        p.update(dict(PRIMER_LEFT_INPUT=c['F'], PRIMER_RIGHT_INPUT=c['R'],
                      PRIMER_PRODUCT_MIN=str(URUN_ALT), PRIMER_PRODUCT_MAX=str(URUN_UST),
                      PRIMER_SPECIFICITY_DATABASE='nt',
                      TOTAL_PRIMER_SPECIFICITY_MISMATCH='5',
                      PRIMER_3END_SPECIFICITY_MISMATCH='2',
                      MISMATCH_REGION_LENGTH='5'))
        _ent_c = _ent_of(HARITA[ad]) if ad in HARITA else ent
        _kendi_dislandi = ad in HARITA
        if _ent_c:
            p['ENTREZ_QUERY'] = _ent_c
            if ad in HARITA:
                yaz(u'  [%s] kendi taksonu dislandi: txid%s' % (ad[:40], HARITA[ad]))
        # ORGANISM dict'e KONULMAZ: coklu organizma TEKRARLANAN alan demektir ve
        # bir dict tek anahtar tutar. Ikili listesi olarak eklenir.
        alanlar = list(p.items()) + [('ORGANISM', o) for o in orgs]
        try:
            veri = urllib.parse.urlencode(alanlar).encode()
            s = pb_ac(PB_URL, veri, yaz=yaz)
            m = re.search(r'job_key=([A-Za-z0-9_\-]+)', s)
            if not m:
                # Sebebi TAHMIN ETME - NCBI'nin kendi hata metnini oku ve yaz.
                _h = re.search(r'(?:Exception error|Error)\s*:\s*([^<\n]{5,300})', s, re.I)
                _sebep = _h.group(1).strip() if _h else u'NCBI hata metni bulunamadi'
                open(os.path.join(ham, '%s_ANAHTARSIZ.html' % re.sub(r'\W+', '_', ad)),
                     'w', encoding='utf-8').write(s)
                out[ad] = dict(durum='BASARISIZ',
                               not_=u'is anahtari alinamadi - NCBI yaniti: %s' % _sebep)
                yaz(u'  [%s] NCBI: is anahtari alinamadi - NCBI diyor ki: %s'
                    % (ad, _sebep)); continue
            anahtar = m.group(1)
            yaz(u'  [%s] NCBI isi gonderildi (%s), bekleniyor...' % (ad, anahtar))
            son = ''
            for i in range(tur_ust):
                time.sleep(bekleme)
                u2 = PB_URL + '?job_key=' + anahtar
                son = pb_ac(u2, yaz=yaz)
                # 2026-08-06 HATA DUZELTMESI - siparis oncesi olcumde yakalandi.
                # Eski kosul yalniz 'still running' ve 'please wait' ariyordu.
                # Primer-BLAST'in bekleme sayfasi bu iki dizgeyi ICERMIYOR; sayfa
                # "Status Running" ve "Time since submission" yaziyor. Sonuc:
                # dongu ILK yoklamada kiriliyor, henuz KOSAN bir isin sayfasinda
                # 'product length' bulunamiyor ve satir "hedef disi 0" yani
                # TEMIZ raporlaniyordu. Siparis oncesi yanlis guvence tam olarak
                # kacinmamiz gereken sey. Artik BITTIGI POZITIF olarak dogrulanir.
                _d = son.lower()
                _kosuyor = ('still running' in _d or 'please wait' in _d
                            or 'status</th><td>running' in _d.replace(' ', '')
                            or re.search(r'status[^<]*<[^>]*>\s*running', _d)
                            or 'time since submission' in _d)
                _bitti = bool(re.search(r'primer pair \d|no significant|not find any target'
                                        r'|unintended template', _d))
                if _bitti and not _kosuyor:
                    break
                print('     ... %d. yoklama (%s)          '
                      % (i + 1, 'kosuyor' if _kosuyor else 'bekleniyor'),
                      end='\r', flush=True)
            open(os.path.join(ham, '%s.html' % re.sub(r'\W+', '_', ad)), 'w',
                 encoding='utf-8').write(son)
            _d = son.lower()
            _bitti = bool(re.search(r'primer pair \d|no significant|not find any target'
                                    r'|unintended template', _d))
            if not _bitti:
                # Is bitmeden tavana carptik. TEMIZ demek YASAK - bilinmiyor denir.
                out[ad] = dict(durum='BASARISIZ',
                               not_=u'NCBI isi %d yoklamada bitmedi (hala kuyrukta). '
                                    u'Sonuc BILINMIYOR - temiz sayilmadi.' % tur_ust)
                yaz(u'  [%s] NCBI: is bitmedi (kuyruk). ELLE yola dusun.' % ad)
                continue
            n = len(re.findall(r'product length\s*=\s*\d+', son, re.I))
            hedefsiz = bool(re.search(r'no significant|not find any target', son, re.I))
            # D-13b (2026-08-07, OLCULDU): sayfa urunleri AYRI BOLUMLERE koyar
            # ('Products on intended targets', '... on potentially unintended
            # templates', '... on target templates'). Butun sayfayi saymak bu
            # ayrimi yok eder. 2026-08-07 kosusunun 22 ham HTML'inin 22'sinde
            # "potentially unintended templates" bolumu BOSTU ve butun urunler
            # "target templates" altindaydi - yani max(0,n-1) bir HEDEF DISI
            # SAYIMI DEGIL, nt icinde bulunan TOPLAM urun sayisidir. Ustelik
            # organizma kisiti Bacteria/Archaea/Fungi oldugu icin hedefin kendi
            # uyeleri de o listede. Ornek olcum: Proteiniphilum_cinsi 876 urun,
            # basliklarin 110'u "Proteiniphilum", 760'i adlandirilamayan cevre
            # klonu ("uncultured bacterium clone ...") - yani NCBI baslıklari
            # bu sorunun cevabini TASIMIYOR.
            n_unint = _ncbi_bolum_say(son, 'potentially unintended templates')
            n_target = _ncbi_bolum_say(son, 'target templates')
            # D-3 HATA DUZELTMESI (2026-08-06): 'n' iki ayri sekilde sahte deger
            # uretiyordu ve ikisi de TAMAM sayiliyordu.
            #   (a) TAVAN: Primer-BLAST sonuc sayfasi en cok 1000 urun listeler.
            #       n==1000 ise gercek sayi 1000 ya da DAHA FAZLA'dir; max(0,n-1)
            #       ile 999 yazmak bir SAYIM degil, tavana carpma isaretidir.
            #       Bu kosuda tam bes hedefte (Metanomikrobiyales, Nitrosocosmicus,
            #       Microascaceae, Metanojen_universal, Mantar F2) 999 cikti -
            #       hepsi genis kapsamli hedef, hepsi tavana carpmis.
            #   (b) BOS: 'Products on ...' bolumleri BOS donen sayfada n==0 olur ve
            #       max(0, 0-1)==0 yani TEMIZ raporlanirdi. Hicbir urun (hedefteki
            #       bile) listelenmemis bir sayfa 'temiz' DEGIL, 'veri yok'tur.
            #       Bu kosuda dokuz hedefte boyle oldu.
            _ur = _ncbi_urunleri(son)
            _adli = [(a, b, L) for a, b, L in _ur if not _adsiz_mi(b)]
            _adsiz = [(a, b, L) for a, b, L in _ur if _adsiz_mi(b)]
            if not hedefsiz and n >= NCBI_SONUC_TAVANI:
                # Sayfa kirpilmis. Ama kirpilmis listede ADLI bir hedef disi
                # takson varsa o bir ALT SINIRDIR ve gecerlidir: "en az bu kadar"
                # denebilir. Sifir cikarsa hicbir sey denemez - liste eksik.
                if _kendi_dislandi and _adli:
                    out[ad] = dict(durum='TAMAM (alt sinir)', hedef_disi=len(_adli),
                                   ncbi_toplam_urun=n, ncbi_adsiz_klon=len(_adsiz),
                                   ncbi_ornek=u'; '.join(b[:70] for _a, b, _L in _adli[:3]),
                                   not_=u'Sayfa tavana carpti (%d urun, gercek sayi daha '
                                        u'fazla). Hedefin kendi taksonu (txid%s) dislandi. '
                                        u'Kirpilmis listede ADLI hedef disi takson: %d - bu '
                                        u'bir ALT SINIRDIR, kesin sayi degildir. Adsiz cevre '
                                        u'klonu: %d (etiketten karar verilemez, katman 2-3 '
                                        u'karar verir).'
                                        % (n, HARITA.get(ad, '?'), len(_adli), len(_adsiz)))
                    yaz(u'  [%s] NCBI: tavan (%d) ama ADLI hedef disi >= %d (alt sinir)'
                        % (ad, n, len(_adli)))
                    continue
                out[ad] = dict(durum='BASARISIZ - SONUC TAVANI',
                               ncbi_toplam_urun=n, ncbi_adsiz_klon=len(_adsiz),
                               not_=u'Primer-BLAST %d urun listeledi (sayfa tavani). '
                                    u'Gercek sayi >= %d. Bu bir SAYIM DEGIL; hukum '
                                    u'icin kullanilamaz. Organizma kisiti (--organizma) '
                                    u'ile daraltip yeniden kosun.'
                                    % (n, NCBI_SONUC_TAVANI))
                yaz(u'  [%s] NCBI: SONUC TAVANI (%d) - sayim degil, sinanmadi' % (ad, n))
                continue
            if not hedefsiz and n == 0:
                out[ad] = dict(durum='BASARISIZ - BOS SONUC',
                               not_=u'Sayfa bitti ama hicbir "product length" satiri yok '
                                    u'(hedefteki urun bile listelenmemis). Bu TEMIZ degil, '
                                    u'VERI YOK. Sinanmadi sayilir.')
                yaz(u'  [%s] NCBI: BOS sonuc sayfasi - sinanmadi' % ad)
                continue
            # ORGANIZMA KISITI YOKSA (--organizma bos) Primer-BLAST hedefin KENDI
            # uyelerini de "unintended template" altinda listeler; max(0,n-1) ancak
            # "tam bir tane amaclanan urun var" varsayimiyla dogrudur ve grup/
            # evrensel primerlerde bu varsayim GECERSIZDIR. Bunu acikca isaretle.
            # D-13b: hukme giren deger artik "unintended templates" BOLUMUNUN
            # sayimidir. O bolum bos ve butun urunler "target templates" altinda
            # ise sayfa bu soruya cevap VERMIYOR (hedef sablonu bildirilmedigi
            # icin Primer-BLAST hicbir urunu 'unintended' saymiyor). Bu durumda
            # katman OY VERMEZ - 'temiz' sayilmaz.
            if hedefsiz:
                out[ad] = dict(durum='TAMAM', hedef_disi=0, ncbi_toplam_urun=0,
                               not_=u'Primer-BLAST hic urun bulamadi.')
                yaz(u'  [%s] NCBI: hic urun yok -> hedef disi 0' % ad)
                continue
            if n_unint == 0 and n_target > 0:
                if not _kendi_dislandi:
                    out[ad] = dict(
                        durum='BASARISIZ - DISLAMA HARITASINDA YOK',
                        ncbi_toplam_urun=n_target, ncbi_adsiz_klon=len(_adsiz),
                        not_=u'Sayfa %d urun listeledi ama bu hedef icin '
                             u'screening/hedef_taxid.tsv icinde dislanacak takson '
                             u'yazili degil. Hedefin kendi uyeleri de listede olabilir, '
                             u'ayirt edilemez. SINANMADI.' % n_target)
                    yaz(u'  [%s] NCBI: dislama haritasinda yok - sinanmadi' % ad)
                    continue
                # Kendi taksonu dislandi. Bolum basligi acilmasa da geriye kalan
                # ADLI her urun tanimi geregi hedef disidir.
                out[ad] = dict(
                    durum='TAMAM', hedef_disi=len(_adli),
                    ncbi_toplam_urun=n_target, ncbi_adsiz_klon=len(_adsiz),
                    ncbi_ornek=u'; '.join(b[:70] for _a, b, _L in _adli[:3]),
                    not_=u'Hedefin kendi taksonu (txid%s) ENTREZ_QUERY ile dislandi. '
                         u'%d urunun %d tanesi ADLI takson (hedef disi kaniti), %d '
                         u'tanesi adsiz cevre klonu ("uncultured ...") - adsizlar '
                         u'hicbir taksona bagli olmadigi icin dislama suzgeci onlara '
                         u'islemez ve hedefin KENDISI olabilirler; kimliklerine dizi '
                         u'karsilastirmasi (katman 2-3) karar verir, hukme girmezler. '
                         u'Bolum basligina bakilmadi - sablon dizi bildirilmedigi surece '
                         u'Primer-BLAST "unintended" bolumunu hic acmiyor (olculdu).'
                         % (HARITA.get(ad, '?'), n_target, len(_adli), len(_adsiz)))
                yaz(u'  [%s] NCBI: adli hedef disi %d / adsiz klon %d / toplam %d'
                    % (ad, len(_adli), len(_adsiz), n_target))
                continue
            _kusur = (u'ORGANIZMA KISITI YOK: hedefin kendi uyeleri de "unintended" '
                      u'altinda sayilmis olabilir. ' if not organizma else u'')
            out[ad] = dict(durum='TAMAM', hedef_disi=n_unint,
                           ncbi_toplam_urun=n_target,
                           not_=_kusur + u'"unintended templates" bolumunun sayimi; '
                                u'ham yanit ncbi_ham/ altinda')
            yaz(u'  [%s] NCBI: hedef disi (unintended bolumu) %s / toplam urun %s%s'
                % (ad, n_unint, n_target, u' (organizma kisiti yok)' if _kusur else u''))
        except Exception as e:
            out[ad] = dict(durum='BASARISIZ', not_=u'%s: %s' % (type(e).__name__, e))
            yaz(u'  [%s] NCBI BASARISIZ (%s) - elle yola dusun' % (ad, type(e).__name__))
    return out


def katman2_elle_girdi(ciftler, cikti, yaz, organizma=''):
    """Kullanicinin Chromium'da Primer-BLAST kosmasi icin hazir girdi + sonuc sablonu."""
    g = os.path.join(cikti, 'NCBI_PRIMER_BLAST_GIRDI.tsv')
    with open(g, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# NCBI Primer-BLAST icin HAZIR GIRDI - dogrudan yapistirilabilir.\n')
        fh.write(u'# Adres: https://www.ncbi.nlm.nih.gov/tools/primer-blast/\n')
        fh.write(u'# Her satir bir cift. Sayfada su alanlara yapistirin:\n')
        fh.write(u'#   "Primer Parameters > Forward primer"  <- F sutunu\n')
        fh.write(u'#   "Primer Parameters > Reverse primer"  <- R sutunu\n')
        fh.write(u'#   "Exon/intron selection > PCR product size"  Min/Max <- urun_min/urun_max\n')
        fh.write(u'#   "Primer Pair Specificity Checking Parameters":\n')
        fh.write(u'#       Database = nt ; Organism = organizma_kisiti sutunu\n')
        fh.write(u'#       Total mismatches = 5 ; 3\' end mismatches = 2\n')
        fh.write(u'# Sonuclari NCBI_SONUC_SABLONU.tsv icine yazip --ncbi-yukle ile geri verin.\n')
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'F', 'R', 'urun_min', 'urun_max', 'organizma_kisiti', 'not'])
        for c in ciftler:
            w.writerow([c['hedef'], c['F'], c['R'], URUN_ALT, URUN_UST,
                        organizma or '(bos = tum nt)', c['tur']])
    yaz(u'  yazildi: %s' % g)

    s = os.path.join(cikti, 'NCBI_SONUC_SABLONU.tsv')
    with open(s, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# NCBI sonuclarini BURAYA yazin, sonra:\n')
        fh.write(u'#   screening.bat -> (D) -> "elle sonuc yukle" ya da\n')
        fh.write(u'#   python3 verification/dogrulama_turu.py --kok . --ncbi-yukle '
                 u'DOGRULAMA_SONUC/NCBI_SONUC_SABLONU.tsv\n')
        fh.write(u'# hedef_disi_urun_sayisi: Primer-BLAST\'in "Products on '
                 u'potentially unintended templates" altinda saydigi urun sayisi.\n')
        fh.write(u'# Hic yoksa 0 yazin. Bakmadiysaniz bos birakin (o satir '
                 u'"NCBI yapilmadi" sayilir).\n')
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'hedef_disi_urun_sayisi', 'en_yakin_hedef_disi_organizma', 'notunuz'])
        for c in ciftler:
            w.writerow([c['hedef'], '', '', ''])
    yaz(u'  yazildi: %s' % s)
    return g, s


def ncbi_yukle(yol, yaz=None):
    """Elle doldurulan NCBI sablonunu okur.

    A5 DUZELTMESI (2026-08-21): bozuk sayi alani eskiden sessizce 'pass'
    ediliyordu. Insan sablonu doldururken '~3' ya da '3 (belki)' yazdiginda o
    hedef icin NCBI katmani HIC olusmuyor, hukum tablosunda 'BILINMIYOR'
    goruluyordu. Kullanici degeri girdigini saniyor, katman dusmus oluyordu.
    BOS birakmak ile BOZUK yazmak ayni sey degildir: birincisi 'bakmadim'
    demektir ve mesrudur, ikincisi bir yazim hatasidir ve gorunmelidir.
    """
    out = {}
    bozuk = []
    with open(yol, encoding='utf-8') as fh:
        for r in csv.DictReader((s for s in fh if not s.startswith('#')), delimiter='\t'):
            n = (r.get('hedef_disi_urun_sayisi') or '').strip()
            if n == '':
                continue                      # bilerek bos: "bakmadim", mesru
            try:
                out[r['hedef'].strip()] = dict(durum='TAMAM (elle)', hedef_disi=int(n),
                                               not_=(r.get('en_yakin_hedef_disi_organizma') or ''))
            except ValueError:
                bozuk.append((r.get('hedef', '?').strip(), n))
    if bozuk and yaz:
        yaz(u'  UYARI: %d satirda "hedef_disi_urun_sayisi" sayi DEGIL. Bu hedefler '
            u'icin NCBI katmani OLUSMADI ve hukumde BILINMIYOR gorunecek:' % len(bozuk))
        for h, n in bozuk:
            yaz(u'    %-40s deger: %r  (tam sayi bekleniyor; hic yoksa 0 yazin, '
                u'bakmadiysaniz BOS birakin)' % (h[:40], n))
    return out


# ---------------------------------------------------------------- birlestir
# Bir katmanin ham sayisini uc duruma indirger. "BILINMIYOR" AYRI BIR DURUMDUR:
# olculmemis bir katman "temiz" sayilmaz, yoksa kosulmamis katmanlar sessizce
# olumlu oy verirdi.
# C-1 HATA DUZELTMESI (2026-08-07): ozet satiri okunmuyordu.
#   "INCELEME - gevsek olcut vurusu (1 adet); 3' son iki baz sinanmali: 1
#    KOSULLU ...: 6   RISKLI ...: 3
#    INCELEME - gevsek olcut vurusu (11 adet)...: 1
#    INCELEME - gevsek olcut vurusu (22 adet)...: 1 ..."
# Sebep: sayac HUKUM DIZGESININ TAMAMINI anahtar yapiyordu. Dizgenin icinde
# vurus SAYISI gectigi icin ("(11 adet)", "(22 adet)") ayni kategorideki her
# satir AYRI bir anahtar oluyordu; 16 cift 10+ sozde kategoriye dagiliyor ve
# ozet hicbir sey ozetlemiyordu. Kategori artik SAYIDAN AYRILIR: hukum
# dizgesinin ilk sozcugu kategoridir, gerisi ayrinti olarak AYRI satirda durur.
KATEGORILER = ('KESIN', 'KOSULLU', 'INCELEME', 'RISKLI', 'CELISKILI', 'EKSIK')
# Bu dordu SIFIR olsa da yazilir: "INCELEME gorunmuyor" ile "INCELEME sifir"
# arasindaki farki okuyan kisi gormeli.
ANA_KATEGORILER = ('KESIN', 'KOSULLU', 'INCELEME', 'RISKLI')


def karar_kategorisi(karar):
    """Hukum dizgesinden SAYISIZ kategori anahtarini cikarir."""
    k = (karar or '').strip()
    # D-16 (2026-08-07): maxsplit KONUMSAL argument olarak verilmisti;
    # Python 3.13'ten itibaren DeprecationWarning, sonra hata olacak.
    ilk = re.split(r'[\s\-:(]', k, maxsplit=1)[0].upper()
    return ilk if ilk in KATEGORILER else (k.split(' ')[0].upper() or 'BILINMIYOR')


# --------------------------------------------------------------- D-18 olcut
# "3' son iki baz TAM eslesiyor mu" ikili sorusunun YERINE gecen olcut.
#
# NEDEN DEGISTI. Eski soru bir olcum degildi: MFEprimer 3.0 ve sonrasi 3'
# terminal bazda uyumsuzluga tanim geregi izin vermiyor [K16], dolayisiyla
# 69/69 kayitta terminal bazin uyusmasi algoritmanin zorunlu ciktisiydi.
# Ustelik ikili soru, literaturde VERIMLI cogaldigi olculmus terminal
# uyumsuzluklari (ozellikle terminal T) sistematik olarak eliyordu [K13].
#
# YENI OLCUT. Uyumsuzluk var/yok degil, BEKLENEN DONGU CEZASI. Bir baglanma
# bolgesi ancak beklenen cezasi o hedef icin GEREKLI dCq'dan kucukse gercek
# rakip sayilir. Boylece olcut, bolluga duyarli esikle ayni para biriminde
# konusur (dongu) ve iki olcut carpistirilabilir.
#
# SAYILAR LITERATURDEN, uydurulmadi (LITERATUR_2026-08-07.md, bolum 3):
#   [K13] Kwok 1990 : 3' terminal A:G, G:A, C:C  -> ~100 kat  = 6,6 dongu
#                     3' terminal A:A            -> ~20 kat   = 4,3 dongu
#                     diger 3' terminal uyumsuzluklar VERIMLI cogaldi = 0
#                     terminal T (T:G, T:C, T:T) en az etkili; bitisik bir
#                     uyumsuzlukla birlikte bile anlamli amplifikasyon = 0
#   [K14] Bru 2008  : 3' uctan BIR ONCEKI pozisyon -> ~3 log = 10 dongu
#                     3' uctan 5., 6., 8. pozisyon -> ~1 log = 3,3 dongu
#                     geri primerde 3' uctan 4 bazdan uzagi -> etkisiz = 0
#   [K17] Sozhamannan 2025 : son DORT baz icinde UC uyumsuzluk -> >15 dongu
#
# OLCULMEYEN ARALIK ACIKCA ISARETLI: 3. ve 4. pozisyon icin dogrudan olcum
# YOK; Bru'nun 5-8 degeri ALT SINIR olarak kullanilir ve 'tahmin' dondurulur.
#
# BU FONKSIYON TEK BASINA HUKUM VERMEZ. Ciktisi bir dongu tahminidir; hukum,
# gerekli dCq = log2(R) + 4,3 ile karsilastirilarak verilir.

UC3_TERMINAL_CEZA = {              # (primer_bazi, sablon_bazi) -> dongu
    ('A', 'G'): 6.6, ('G', 'A'): 6.6, ('C', 'C'): 6.6,
    ('A', 'A'): 4.3,
    ('T', 'G'): 0.0, ('T', 'C'): 0.0, ('T', 'T'): 0.0,
}
UC3_TERMINAL_VARSAYILAN = 0.0      # Kwok: geri kalani verimli cogaldi


def uc3_ceza_dongu(uyumsuz_konumlar, terminal_ciftler=None):
    """3' uca yakin uyumsuzluklarin BEKLENEN dongu cezasini dondurur.

    uyumsuz_konumlar : 3' uctan 1 tabanli konum listesi (1 = terminal baz).
    terminal_ciftler : {konum: (primer_bazi, sablon_bazi)} - yalniz konum 1
                       icin kullanilir; verilmezse Kwok'un en kotu degeri
                       (6,6 dongu) alinir, yani TEMKINLI taraf.

    Doner: (ceza_dongu, dayanak_metni, olculdu_mu)
    """
    if not uyumsuz_konumlar:
        return (0.0, 'uyumsuzluk yok', True)
    kon = sorted(set(int(x) for x in uyumsuz_konumlar))
    tahmin = False
    parcalar = []
    toplam = 0.0

    # Sozhamannan: son dort baz icinde uc ve uzeri uyumsuzluk
    son4 = [k for k in kon if k <= 4]
    if len(son4) >= 3:
        return (15.0, 'son 4 baz icinde %d uyumsuzluk -> >15 dongu [K17]'
                % len(son4), True)

    for k in kon:
        if k == 1:
            cift = (terminal_ciftler or {}).get(1)
            if cift:
                c = UC3_TERMINAL_CEZA.get(tuple(cift), UC3_TERMINAL_VARSAYILAN)
                parcalar.append('terminal %s:%s = %.1f [K13]' % (cift[0], cift[1], c))
            else:
                c = 6.6
                parcalar.append('terminal (baz cifti bilinmiyor, en kotu) = 6,6 [K13]')
            toplam += c
        elif k == 2:
            toplam += 10.0
            parcalar.append('sondan ikinci = 10,0 [K14]')
        elif k in (3, 4):
            toplam += 3.3
            tahmin = True
            parcalar.append('konum %d = 3,3 (OLCULMEDI, 5-8 degerinden alt sinir) [K14]' % k)
        elif k <= 8:
            toplam += 3.3
            parcalar.append('konum %d = 3,3 [K14]' % k)
        else:
            parcalar.append('konum %d = 0,0 (3\' uctan uzak) [K13, K14]' % k)
    return (toplam, ' + '.join(parcalar), not tahmin)


def hukum(v):
    """Bir katmanin sonucunu TEMIZ / RISKLI / BILINMIYOR'a indirger."""
    if v is None:
        return 'BILINMIYOR'
    if isinstance(v, str):
        return 'BILINMIYOR'
    return 'TEMIZ' if v <= VURUS_ESIGI else 'RISKLI'


# ---------------------------------------------------------------------------
# DORT KAYNAK YAN YANA. Karar kaynaklarin UYUSMASINA baglanir, bir kaynagin
# kendi sayisina degil:
#   dort kaynak da TEMIZ            -> KESIN
#   uc kaynak TEMIZ, biri eksik     -> KOSULLU
#   kaynaklar AYRILIYOR             -> CELISKILI (siparis edilemez)
#   hicbiri sonuc vermedi           -> EKSIK
# "Sifir veritabani tarandi" ile "hepsi tarandi ve temiz cikti" ayni sey
# degildir; bu yuzden yerel katmanin degeri ancak gercekten tarama yapildiysa
# (tarandi > 0) sayisal kabul edilir, aksi halde BILINMIYOR olur.
# ---------------------------------------------------------------------------
def birlestir(ciftler, yerel, ncbi, mfe=None, klad=None):
    """UC OLCUM katmani yan yana (yerel DB / MFEprimer / NCBI).

    numune katmani KASTEN oy vermez: sabit 'TEMIZ' uretir, dolayisiyla bir
    olcum degil totolojidir (D-2 duzeltmesi, 2026-08-06).
    Uc katman da TEMIZ ve uctu de olctuyse KESIN; ikisi olctuyse KOSULLU;
    olcen katmanlar AYRILIYORSA CELISKILI."""
    out = []
    for c in ciftler:
        h = c['hedef']
        # D-2 HATA DUZELTMESI (2026-08-06): n_ok SABIT 'TEMIZ' idi ve yine de
        # 'bilinen' oy kumesine katiliyordu. Bu bir totoloji: numune katmani bir
        # OLCUM degil, cifti bu listeye SOKAN kabul olcutudur - tanimi geregi
        # asla RISKLI diyemez. Sabit bir TEMIZ oyu, uyusma testine sokuldugunda
        # herhangi bir katmanin TEK bir RISKLI okumasi set'i {TEMIZ,RISKLI}
        # yapiyor ve satir ZORUNLU olarak CELISKILI cikiyordu. 16 ciftin 16'si
        # bu yuzden celiskili isaretlendi. Artik numune SUTUN olarak gosterilir
        # ama OY VERMEZ.
        n_ok = 'TEMIZ'
        yv = yerel.get(h, {}) if yerel else {}
        # O-1: 'sifir veritabani tarandi' ile 'hepsi tarandi, temiz' ayni degildir
        # D-1: hukme TOPLAM vurus degil, beklenen boydan FARKLI olanlar girer.
        y_tum = yv.get('urun') if yv.get('tarandi') else None
        if yv.get('tarandi') and not yv.get('boy_ayrimi_yok'):
            y = yv.get('hedef_disi')
        else:
            y = None            # boy ayrimi yapilamadiysa bu katman hukum vermez
        y_ok = hukum(y)
        nb = ncbi.get(h) if ncbi else None
        nb_v = nb.get('hedef_disi') if nb and nb.get('durum', '').startswith('TAMAM') else None
        nb_ok = hukum(nb_v)
        mv = (mfe or {}).get(h)
        m_ham = mv.get('hedef_disi') if mv else None
        # D-12 (2026-08-07): HUKME GIREN OLCU DEGISTI. MFEprimer'in "hedef disi"
        # sayisi SADECE BOYA dayanir; evrensel/grup primerlerinde hedef kladin
        # kendi uyeleri de farkli boyda amplikon verir. OLCUM (2026-08-07 kosusu,
        # 1605 amplikon, mfe_hedef_disi_kimlikleri.tsv taksonomi dizgeleri):
        #   (a) hedef klad ici, boyu farkli  1536  (b) ayni alan/klad disi  24
        #  (ao) hedef alan ici ama ORGANEL     31  (c) farkli alan          14
        # Yani ham sayinin %95,7'si zararsiz uzunluk varyanti. Hukme artik
        # klad_disi = (b)+(c) giriyor; ham sayi sutun olarak KALIYOR.
        kl = (klad or {}).get(h)
        if kl is not None:
            m_urun = kl['klad_disi']
            m_klad_ayrimi = True
        else:
            m_urun = m_ham
            m_klad_ayrimi = False
        m_ok = hukum(m_urun) if mv else 'BILINMIYOR'
        # D-2: numune BILEREK disarida - sabit deger oy veremez.
        kaynaklar = dict(yerel=y_ok, mfeprimer=m_ok, ncbi=nb_ok)
        bilinen = [x for x in kaynaklar.values() if x != 'BILINMIYOR']
        uyusan = 0
        if bilinen:
            en_cok = max(bilinen.count(x) for x in set(bilinen))
            uyusan = en_cok
        n_kaynak = len(bilinen)
        # D-6 (2026-08-06): 'yerel RISKLI + MFEprimer TEMIZ' bir CELISKI DEGILDIR.
        # Iki olcut IC ICE gecmis: yerel tarama en cok 5 toplam uyumsuzluga izin
        # verir ve 3' son iki baz sartini UYGULAMAZ (kuresel_tarama.py, need_tail=
        # False); MFEprimer ise termodinamik olcut kullanir (Tm kesimi 30 C). Yani
        # yerel olcut MFEprimer olcutunu KAPSAR. Gevsek olcutun sikisindan FAZLA
        # vurus bulmasi BEKLENEN sonuctur, bir kaynak catismasi degil. Gercek
        # celiski TERSIDIR: sikinin bulup gevsegin kacirdigi vurus.
        # Bu satirlar 'temiz' de sayilmaz - siparise girmez, 3' ucu elle sinanir.
        # D-15 (2026-08-07): yerel katmanin hedef disi sayisi TAKSONOMIK olarak
        # suzulemiyor. Sebep olculdu: yerel_vuruslar.tsv kume basina en cok 20
        # vurus tutuyor (raporla(), '_vurus' listesi), yani 4702 farkli-boy
        # vurusun kimligi diskte YOK. Bakteri_universal icin ornekteki 100
        # vurusun yalniz 2'si farkli boyda - ornek bu soruya cevap vermiyor.
        # Bu yuzden yerel katmanin RISKLI oyu, MFEprimer'in KLAD SUZGECINDEN
        # gecmis sayisi 0 iken tek basina RISKLI uretemez -> INCELEME.
        _gevsek_fazla = (y_ok == 'RISKLI' and m_ok == 'TEMIZ'
                         and nb_ok != 'RISKLI')
        # D-17 (2026-08-07, OLCULDU): ORGANEL urunleri gizlenmesin.
        # SILVA kloroplast ve mitokondri kayitlarini "Bacteria;..." ile baslatir
        # (mitokondri Rickettsiales, kloroplast Cyanobacteriota altinda). Alan
        # sinamasi bu yuzden onlari hedef ICI sayar - Bakteri_universal icin
        # taksonomik olarak dogru ama PRATIKTE yanlis: bunlar bitki organeli.
        # OLCUM: Bakteri_universal'in 31 organel amplikonunun 31'inde de F ve R
        # uyumsuzlugu SIFIR ve FpTm 62,97 / RpTm 61,33 - ikisi de Ta 57,9 C'nin
        # USTUNDE. Yani 31 urunun 31'i standart kosulda OLUSUR (91-302 bp,
        # beklenen 130 bp). Konaklar: Azolla, Isoetes, Equisetum, Ipomoea,
        # Welwitschia, Silene... yani bitki besleme yapilan bir curutucude
        # gerceklesebilecek urunler. Bu satir 'temiz' sayilamaz.
        _organel = (kl or {}).get('ao') or 0
        _organel_notu = ''
        if _organel:
            _organel_notu = (u' | ORGANEL UYARISI: %d konak organel (kloroplast/'
                             u'mitokondri) amplikonu; olusabilir %s'
                             % (_organel, (kl or {}).get('olusabilir')))
        if _gevsek_fazla:
            # D-18 (2026-08-09): "3' son iki baz" HUKUM VERMEZ, artik istenmez.
            # Gerekce OLCULDU ve iki katmanlidir:
            #  (a) MFEprimer 3.0 ve sonrasi 3' terminal bazda uyumsuzluga TANIM
            #      GEREGI izin vermiyor; 69/69 kayitta terminal bazin uyusmasi
            #      algoritmanin zorunlu ciktisi, verinin bulgusu degil.
            #  (b) Kutu taramasinda son2 sarti KALDIRILARAK yeniden olculdu
            #      (2026-08-09, 17 hedef): dCq degisimi en cok 0,41 dongu
            #      (Proteolitik_Synergistaceae -0,41; Petrimonas +0,36; geri
            #      kalan 15 hedefte |fark| <= 0,09). Hicbir hukum degismedi.
            # Yerine istenen: uyumsuzlugun 3' uca UZAKLIGI ve TIPI (bkz.
            # uc3_ceza_dongu). Ayrinti: ESIK_VE_OLCUT_2026-08-08.md.
            karar = ('INCELEME - gevsek olcut vurusu (%s adet); 3\' uca yakin '
                     'uyumsuzlugun KONUMU ve TIPI degerlendirilmeli '
                     '(son iki baz sarti hukum vermez)' % y)
        elif len(set(bilinen)) > 1:
            karar = 'CELISKILI'
        elif n_kaynak == 0:
            karar = 'EKSIK - hicbir kaynak sonuc vermedi'
        elif set(bilinen) == {'TEMIZ'}:
            if _organel:
                # D-17: butun katmanlar temiz gorunse bile olusabilir organel
                # urunu varsa satir temiz DEGILDIR - insan karari gerekir.
                karar = ('INCELEME - katmanlar temiz ama %d organel amplikonu var'
                         % _organel)
            elif n_kaynak >= 3:
                karar = 'KESIN - uc olcum katmani da uyusuyor'
            elif n_kaynak == 2:
                karar = 'KOSULLU - iki katman uyusuyor, biri eksik'
            else:
                karar = 'EKSIK - yalnizca %d kaynak sonuc verdi' % n_kaynak
        else:
            karar = 'RISKLI - siparis edilmez'
        if _organel_notu and 'ORGANEL' not in karar:
            karar = karar + _organel_notu
        out.append(dict(c, kategori=karar_kategorisi(karar),
                        numune=n_ok, numune_deger=c['numune_deger'],
                        yerel=y_ok, yerel_urun=y, yerel_tum=y_tum,
                        yerel_ayni_boyda=yv.get('ayni_boyda'),
                        yerel_kume=(yerel.get(h, {}) or {}).get('kume', {}),
                        # A2: taksonomik sayaclar. HUKME HENUZ GIRMIYOR -
                        # boy tabanli 'yerel_urun' hukum kaynagi olarak
                        # duruyor. Ikisi yan yana yazilir ki farki once
                        # OLCULSUN, sonra karar verilsin.
                        yerel_klad_ayrimi=(yv.get('siniflandirildi') or False),
                        yerel_a=(yv.get('sinif') or {}).get('a'),
                        yerel_ao=(yv.get('sinif') or {}).get('ao'),
                        yerel_b=(yv.get('sinif') or {}).get('b'),
                        yerel_c=(yv.get('sinif') or {}).get('c'),
                        yerel_bilinmiyor=(yv.get('sinif') or {}).get('bilinmiyor'),
                        yerel_klad_disi=((yv.get('sinif') or {}).get('b', 0)
                                         + (yv.get('sinif') or {}).get('c', 0))
                        if yv.get('siniflandirildi') else None,
                        mfeprimer=m_ok, mfe_urun=m_urun, mfe_ham=m_ham,
                        mfe_klad_ayrimi=m_klad_ayrimi,
                        mfe_a=(kl or {}).get('a'), mfe_ao=(kl or {}).get('ao'),
                        mfe_b=(kl or {}).get('b'), mfe_c=(kl or {}).get('c'),
                        mfe_olusabilir=(kl or {}).get('olusabilir'),
                        mfe_olusmaz=(kl or {}).get('olusmaz'),
                        mfe_durum=(mv or {}).get('durum', 'YAPILMADI'),
                        ncbi=nb_ok, ncbi_urun=nb_v,
                        ncbi_durum=(nb or {}).get('durum', 'YAPILMADI'),
                        kaynak_sayisi=n_kaynak, uyusan=uyusan, karar=karar))
    return out


# ---------------------------------------------------------------------------
# Dort dosya yazar. CELISKILER.md kasten ayri bir dosyadir: bu turun en degerli
# ciktisi celiskilerdir ve uzun bir tablonun icinde kaybolmamalidir.
#
# Celiski yoksa yazilan cumle "her sey temiz" DEGILDIR: kosulmamis katmanlar
# EKSIK sayilir ve tanimi geregi celiski uretmez. Bu ayrim raporun icinde acikca
# yazilidir, cunku "hic celiski cikmadi" ifadesi yanlis okunmaya cok musait.
# ---------------------------------------------------------------------------
def raporla(cikti, satirlar, yaz):
    t = os.path.join(cikti, 'dogrulama_uc_sutun.tsv')
    with open(t, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# UC OLCUM KATMANI YAN YANA. Ayrilirsa satir CELISKILI - siparis edilemez.\n')
        # A3 (2026-08-21): esik artik PT_VURUS_ESIGI ile degistirilebiliyor.
        # Degistirilebilen bir olcut, hangi degerle kosuldugu YAZILMADIKCA
        # olcumu okunamaz kilar; bu yuzden ciktinin basinda durur.
        fh.write(u'# VURUS_ESIGI = %d%s  (bir katman bu sayidan COK hedef disi urun\n'
                 u'#   gorurse RISKLI oy verir; tek basina hukum vermez, karar\n'
                 u'#   katmanlarin UYUSMASINA baglidir)\n'
                 % (VURUS_ESIGI,
                    u'  [PT_VURUS_ESIGI ile DEGISTIRILMIS - varsayilan 0]'
                    if VURUS_ESIGI != 0 else u'  [varsayilan]'))
        fh.write(u'# 1_NUMUNE OY VERMEZ: sabit TEMIZ uretir, bir olcum degil kabul olcutudur (D-2).\n')
        fh.write(u'# 2_hedef_disi_urun: beklenen urun boyundan +-%d bp FARKLI vuruslar.\n' % BOY_TOL)
        fh.write(u'# A2 (2026-08-21): 2_klad_* sutunlari TAKSONOMIK ayrimdir ve\n'
                 u'#   2_hedef_disi_urun ile AYNI SEYI OLCMEZ. Boy olcutu yaniltici:\n'
                 u'#   D-12\'de olculdu, "hedef disi" sayilan 1.605 amplikonun %95,7\'si\n'
                 u'#   hedef kladin KENDI ICINDENdi (sinif a), yalniz boyu farkliydi.\n'
                 u'#   HUKME HALA 2_hedef_disi_urun giriyor; taksonomik sayilar once\n'
                 u'#   OLCULSUN diye yan yana yaziliyor. Fark buyukse hangi olcutun\n'
                 u'#   hukum vermesi gerektigi AYRI bir karardir.\n'
                 u'#   2_klad_disi_b_c = b + c  (gercek capraz adayi)\n'
                 u'#   2_bilinmiyor: alan cozulemedi - KANIT SAYILMAZ, capraz sayilmaz.\n'
                 u'#   2_klad_ayrimi_yapildi HAYIR ise sinif sayilari SIFIRDIR ama bu\n'
                 u'#   "capraz yok" DEGIL, "olculmedi" demektir.\n')
        fh.write(u'# 3_hedef_disi_amplikon: D-12 sonrasi bu sutun TAKSONOMIK olarak\n'
                 u'#   suzulmustur = (b) ayni alan/klad disi + (c) farkli alan. Boya\n'
                 u'#   dayali HAM sayi 3_HAM_boya_dayali sutununda durur. Ikisinin\n'
                 u'#   farki hedef kladin kendi uzunluk varyantlaridir (3_a) ve bunlar\n'
                 u'#   evrensel/grup primerlerinde TASARIM GEREGIDIR.\n'
                 u'# 3_olusabilir_Tm_yakin: (b)+(c)+(ao) kayitlarindan min(FpTm,RpTm)\n'
                 u'#   panelin Ta\'sinin (%.1f C) 5 C'
                 u' altina DUSMEYENLER - standart\n'
                 u'#   kosulda urun verebilecek olanlar.\n' % TA_PANEL)
        fh.write(u'# 2_ayni_boyda_HEDEFIN_KENDISI: beklenen boydaki vuruslar - grup/evrensel\n'
                 u'#   primerlerde bunlar TASARIM GEREGI olup hedef disi SAYILMAZ (D-1).\n')
        fh.write(u'# 4_NCBI durumu "SONUC TAVANI" ya da "BOS SONUC" ise deger bir SAYIM DEGIL,\n'
                 u'#   o hucre SINANMADI demektir (D-3).\n')
        fh.write('# ' + OLCUT_NOTU.strip().replace('\n', '\n# ') + '\n')
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'cift_turu', 'F', 'R', 'urun_bp',
                    '1_NUMUNE_oy_vermez', '1_numune_deger',
                    '2_YEREL_DB', '2_hedef_disi_urun', '2_ayni_boyda_HEDEFIN_KENDISI',
                    '2_tum_vurus', '2_kume_dagilimi',
                    '2_klad_ayrimi_yapildi', '2_klad_disi_b_c',
                    '2_a_klad_ici', '2_ao_organel', '2_b_ayni_alan_klad_disi',
                    '2_c_farkli_alan', '2_bilinmiyor',
                    '3_MFEPRIMER', '3_hedef_disi_amplikon', '3_durum',
                    '3_HAM_boya_dayali', '3_klad_ayrimi_yapildi',
                    '3_a_klad_ici_uzunluk_varyanti', '3_ao_organel',
                    '3_b_ayni_alan_klad_disi', '3_c_farkli_alan',
                    '3_olusabilir_Tm_yakin', '3_olusmaz_Tm_dusuk',
                    '4_NCBI', '4_hedef_disi_urun', '4_durum',
                    'kaynak_sayisi', 'uyusan_kaynak', 'KARAR'])
        for s in satirlar:
            w.writerow([s['hedef'], s['tur'], s['F'], s['R'], s['urun'],
                        s['numune'], s['numune_deger'],
                        s['yerel'], s['yerel_urun'],
                        s.get('yerel_ayni_boyda', ''), s.get('yerel_tum', ''),
                        '; '.join('%s=%s' % kv for kv in (s['yerel_kume'] or {}).items()),
                        ('EVET' if s.get('yerel_klad_ayrimi') else 'HAYIR - yalniz boy'),
                        s.get('yerel_klad_disi', ''),
                        s.get('yerel_a', ''), s.get('yerel_ao', ''),
                        s.get('yerel_b', ''), s.get('yerel_c', ''),
                        s.get('yerel_bilinmiyor', ''),
                        s.get('mfeprimer', '-'), s.get('mfe_urun', ''), s.get('mfe_durum', ''),
                        s.get('mfe_ham', ''),
                        ('EVET' if s.get('mfe_klad_ayrimi') else 'HAYIR - ham sayi'),
                        s.get('mfe_a', ''), s.get('mfe_ao', ''),
                        s.get('mfe_b', ''), s.get('mfe_c', ''),
                        s.get('mfe_olusabilir', ''), s.get('mfe_olusmaz', ''),
                        s['ncbi'], s['ncbi_urun'], s['ncbi_durum'],
                        s.get('kaynak_sayisi', ''), s.get('uyusan', ''), s['karar']])
    yaz('  yazildi: %s' % t)

    celiskili = [s for s in satirlar if s['karar'] == 'CELISKILI']
    c = os.path.join(cikti, 'CELISKILER.md')
    with open(c, 'w', encoding='utf-8') as fh:
        fh.write(u'# Celiskiler — bu turun en degerli ciktisi\n\n')
        if not celiskili:
            fh.write(u'Bu kosuda uc katman **hicbir satirda ayrilmadi**.\n\n'
                     u'Bu, "her sey temiz" demek DEGILDIR: NCBI ya da yerel katmani '
                     u'kosulmamis satirlar "EKSIK" sayilir ve celiski uretmez.\n')
        for s in celiskili:
            fh.write(u'## %s\n\n' % s['hedef'])
            fh.write(u'| kaynak | sonuc | deger |\n|---|---|---|\n')
            fh.write(u'| 1 numune olcumu | (OY VERMEZ - sabit deger) | %s |\n' % s['numune_deger'])
            fh.write(u'| 2 yerel veritabani (bizim) | %s | %s hedef disi (+ %s tanesi hedefin '
                     u'KENDI boyunda, toplam %s) |\n'
                     % (s['yerel'], s['yerel_urun'], s.get('yerel_ayni_boyda', '-'),
                        s.get('yerel_tum', '-')))
            fh.write(u'| 3 MFEprimer (BAGIMSIZ) | %s | %s amplikon |\n' % (s.get('mfeprimer', '-'), s.get('mfe_urun', '-')))
            fh.write(u'| 4 NCBI (BAGIMSIZ) | %s | %s |\n\n' % (s['ncbi'], s['ncbi_urun']))
            fh.write(u'**Ne anlama gelir:** numune 99 kutudan ibarettir. Numunede temiz '
                     u'gorunup veritabaninda vurus veren bir cift, numunede BULUNMAYAN '
                     u'bir organizmayi cogaltiyor demektir - o organizma laboratuvara '
                     u'baska bir numuneyle girerse yalanci urun verir. Tersi de olur: '
                     u'veritabaninda temiz gorunup numunede rakip veren cift, '
                     u'veritabaninda temsil edilmeyen bir soyu cogaltiyordur.\n\n')
            fh.write(u'**Yapilmasi gereken:** bu satir SIPARIS EDILEMEZ. Vurusun hangi '
                     u'organizmada oldugu `dogrulama_uc_sutun.tsv` ve `yerel_vuruslar.tsv` '
                     u'dosyalarindan okunup, o organizmanin bu matriste bulunup '
                     u'bulunmadigi karara baglanmalidir.\n\n')
    yaz('  yazildi: %s' % c)

    v = os.path.join(cikti, 'yerel_vuruslar.tsv')
    with open(v, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'kume', 'kayit_basligi', 'urun_bp', 'F_uyumsuzluk', 'R_uyumsuzluk'])
        for s in satirlar:
            for vv in (s.get('_vurus') or []):
                w.writerow([s['hedef']] + list(vv))
    yaz('  yazildi: %s' % v)

    r = os.path.join(cikti, 'DOGRULAMA_RAPORU.md')
    # KATEGORI bazinda say (sayidan arindirilmis anahtar - bkz. C-1 notu).
    say = {}
    for s in satirlar:
        kg = s.get('kategori') or karar_kategorisi(s.get('karar'))
        say[kg] = say.get(kg, 0) + 1
    # ayrinti ayri tutulur: kategori -> {tam hukum dizgesi: adet}
    ayrinti = {}
    for s in satirlar:
        kg = s.get('kategori') or karar_kategorisi(s.get('karar'))
        ayrinti.setdefault(kg, {})
        ayrinti[kg][s['karar']] = ayrinti[kg].get(s['karar'], 0) + 1
    with open(r, 'w', encoding='utf-8') as fh:
        fh.write(u'# Dogrulama turu\n\nUretim: %s · betik %s\n\n'
                 % (time.strftime('%Y-%m-%d %H:%M'), VERSIYON))
        if VURUS_ESIGI != 0:
            fh.write(u'> **UYARI: vurus esigi varsayilan degil.** Bu kosu '
                     u'`PT_VURUS_ESIGI=%d` ile yapildi (varsayilan 0). Siparis '
                     u'oncesi hukumde esik YUKSELTILMEMELIDIR; asagidaki sayilar '
                     u'gevsetilmis olcutle uretilmistir.\n\n' % VURUS_ESIGI)
        fh.write(u'## Sonuc\n\n')
        for k in KATEGORILER:
            if k in say or k in ANA_KATEGORILER:
                fh.write(u'- **%s**: %d\n' % (k, say.get(k, 0)))
        for k in sorted(x for x in say if x not in KATEGORILER):
            fh.write(u'- **%s**: %d\n' % (k, say[k]))
        fh.write(u'- _toplam_: %d\n' % sum(say.values()))
        fh.write(u'\n### Ayrinti (kategori icindeki gerekceler)\n\n')
        for k in list(KATEGORILER) + sorted(x for x in ayrinti if x not in KATEGORILER):
            if k not in ayrinti:
                continue
            for gerekce, n in sorted(ayrinti[k].items(), key=lambda t: -t[1]):
                fh.write(u'- %s (%d): %s\n' % (k, n, gerekce))
        fh.write(u'\n> Celiskiler once okunur: `CELISKILER.md`\n\n')
        fh.write(u'```' + OLCUT_NOTU + u'```\n\n## Uc katman yan yana\n\n')
        fh.write(u'| hedef | 1 numune | 2 yerel DB | 3 MFEprimer | 4 NCBI | uyusan | karar |\n'
                 u'|---|---|---|---|---|---|---|\n')
        for s in satirlar:
            fh.write(u'| %s | %s | %s (%s) | %s (%s) | %s (%s) | %s/%s | **%s** |\n'
                     % (s['hedef'], s['numune'], s['yerel'], s['yerel_urun'],
                        s.get('mfeprimer', '-'), s.get('mfe_urun', '-'),
                        s['ncbi'], s['ncbi_urun'],
                        s.get('uyusan', '-'), s.get('kaynak_sayisi', '-'), s['karar']))
        fh.write(u'\n## NCBI nasil tamamlanir\n\n'
                 u'Otomatik yol kosulmadiysa ya da basarisiz olduysa:\n\n'
                 u'1. `NCBI_PRIMER_BLAST_GIRDI.tsv` dosyasindaki satirlari '
                 u'https://www.ncbi.nlm.nih.gov/tools/primer-blast/ sayfasina yapistirin.\n'
                 u'2. Sonuclari `NCBI_SONUC_SABLONU.tsv` icine yazin.\n'
                 u'3. `screening.bat` -> (D) -> "elle sonuclari yukle" secenegini kosun.\n')
    yaz('  yazildi: %s' % r)
    yaz('')
    # OZET: dort ana kategori tek satirda, SAYIDAN arindirilmis anahtarlarla.
    sirali = [k for k in KATEGORILER if k in say or k in ANA_KATEGORILER] + \
             sorted(x for x in say if x not in KATEGORILER)
    yaz(u'  OZET   ' + '   '.join('%s: %d' % (k, say.get(k, 0)) for k in sirali)
        + '   |   toplam: %d' % sum(say.values()))
    yaz(u'  AYRINTI:')
    for k in sirali:
        for gerekce, n in sorted(ayrinti.get(k, {}).items(), key=lambda t: -t[1]):
            yaz(u'    %-10s %2d  %s' % (k, n, gerekce))



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

# ---------------------------------------------------------------------------
# Bu betikte surucu dogrudan main() icindedir. Katmanlar sirayla kosar:
#   2) yerel veritabani  ->  3) MFEprimer  ->  4) NCBI  ->  birlestir/raporla
# (1. katman numune olcumudur ve K turundan hazir gelir, burada kosulmaz.)
#
# CIKIS KODLARI: 7 = ARDISIK COKUS - K asamasi hic satir uretmedigi icin girdi
# bos; bu D'nin hatasi DEGILDIR ve oyle yazilir. 5 = kendi girdisi eksik,
# 4 = kendi ciktisi bos. Ayrim onemlidir: yanlis asamayi ayiklamaya calismak
# saatler kaybettirir.
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description='Kurtarilan ciftlerin uc katmanli dogrulanmasi')
    p.add_argument('--kok', default='.')
    p.add_argument('--ncbi', choices=['oto', 'elle', 'yok'], default='elle',
                   help='oto: NCBI URL API; elle: yapistirilabilir girdi uret; yok: atla')
    p.add_argument('--ncbi-yukle', default=None, help='doldurulmus NCBI_SONUC_SABLONU.tsv')
    p.add_argument('--organizma', default='', help='NCBI organizma kisiti (bos = tum nt)')
    # D-13c (2026-08-07, OLCULDU): ENTREZ_QUERY ile hedefin KENDI taksonu
    # dislanabilir; o zaman sayfada kalan her urun tanimi geregi hedef disidir.
    # DIKKAT - olculmus tuzak: ciplak 'NOT txidN[Organism]' filtreyi TERSINE
    # CEVIRIR (yalniz o taksonu getirir). Dogru bicim 'all[filter] NOT
    # txidN[Organism]'. Kod bu oneki kendisi ekler.
    p.add_argument('--ncbi-haric-taxid', default='',
                   help="NCBI'de DISLANACAK taxid (ornek: 2157). ENTREZ_QUERY'ye "
                        "'all[filter] NOT txid<N>[Organism]' olarak gonderilir.")
    p.add_argument('--yalniz-yerel', action='store_true', help='yalniz katman 2 (yerel DB)')
    p.add_argument('--mfe-yok', action='store_true',
                   help='MFEprimer katmanini atla')
    p.add_argument('--kume-ust', type=int, default=0,
                   help='yalniz en kucuk N veritabani (hizli sinama)')
    p.add_argument('--parc', action='store_true',
                   help='SILVA LSU Parc kumesini de tara (yavas, ozgullukte gerekmez)')
    p.add_argument('--siparis', action='store_true',
                   help='kurtarilanlar yerine SIPARIS LISTESINDEKI ciftleri dogrula '
                        '(siparis oncesi Primer-BLAST kontrolu icin)')
    p.add_argument('--siparis-hepsi', action='store_true',
                   help='--siparis ile: KOSULLU ve ONERILMEZ satirlari da al')
    # -----------------------------------------------------------------------
    # --tumu  (2026-08-07)
    # Kullanici istegi: "hepsini ekstra o veritabaninda da arasin". --siparis
    # yalniz KESIN+EVRENSEL (16 cift) sinar; KOSULLU ve ONERILMEZ satirlar
    # (6 cift) SILVA'ya HIC sokulmadi. --tumu bunlari da katar -> 22 cift,
    # butun indeksli veritabanlarina, SILVA dahil.
    # --siparis kipi DEGISMEDI; --tumu onun uzerine kurulu ayri bir bayraktir.
    # -----------------------------------------------------------------------
    p.add_argument('--tumu', action='store_true',
                   help='paneldeki BUTUN ciftler (KESIN+EVRENSEL+KOSULLU+'
                        'ONERILMEZ) butun indeksli veritabanlarina, SILVA dahil')
    p.add_argument('--ncbi-yalniz-siparis', action='store_true',
                   help='--tumu ile: KATMAN 4 (NCBI) yalniz siparis listesindeki '
                        '(KESIN/EVRENSEL) ciftlere kosar. Listede olmayanlar '
                        'yalniz yerel + MFEprimer katmanlarini gorur. NCBI cift '
                        'basina ~75 sn + 10 sn bekleme oldugu icin sure kalemi.')
    p.add_argument('--sifirla', action='store_true')
    a = p.parse_args()

    # --tumu = --siparis + --siparis-hepsi. Ayri bir yol DEGIL, var olan
    # yolun genis girdi kumesi. Boylece --siparis kipinin davranisi degismez.
    if getattr(a, 'tumu', False):
        a.siparis = True
        a.siparis_hepsi = True
    if getattr(a, 'ncbi_yalniz_siparis', False) and not getattr(a, 'siparis_hepsi', False):
        # Kisitlanacak bir sey yok: girdi kumesi zaten yalniz siparis listesi.
        a.ncbi_yalniz_siparis = False

    kok = os.path.abspath(a.kok)
    if not os.path.isdir(os.path.join(kok, 'screening')):
        sys.exit('HATA: %s icinde screening yok.' % kok)
    CIKTI = os.path.join(kok, 'DOGRULAMA_SONUC')
    KONTROL = os.path.join(CIKTI, 'kontrol')
    os.makedirs(KONTROL, exist_ok=True)
    if a.sifirla:
        for f in os.listdir(KONTROL):
            try:
                os.remove(os.path.join(KONTROL, f))
            except OSError as e:
                print('  silinemedi: %s (%s)' % (f, e))
    g = open(os.path.join(CIKTI, 'kosu_gunlugu.txt'), 'a', encoding='utf-8')

    def yaz(s=''):
        print(s, flush=True); g.write(s + '\n'); g.flush()

    yaz('=' * 78)
    yaz('  DOGRULAMA TURU - kurtarilan ciftler uc kanit katmaniyla sinaniyor')
    yaz('  surum %s   %s' % (VERSIYON, time.strftime('%Y-%m-%d %H:%M')))
    yaz('=' * 78)

    kyol_ = os.path.join(kok, 'KURTARMA_SONUC', 'kurtarma_satirlari.tsv')
    if not os.path.exists(kyol_) or sum(
            1 for x in open(kyol_, encoding='utf-8') if x.strip() and not x.startswith('#')) <= 1:
        yaz('')
        yaz('  ' + '!' * 70)
        yaz(u'  D ASAMASI CALISTIRILMADI - GIRDI BOS')
        yaz(u'  Sebep: K (verification) asamasi hic satir uretmedi, dolayisiyla')
        yaz(u'  dogrulanacak cift YOK. Bu D\'nin hatasi DEGILDIR; K\'nin')
        yaz(u'  duşmesinin ARDISIK sonucudur.')
        yaz(u'  Yapilacak: once K asamasinin neden satir uretmedigine bakin.')
        yaz(u'  D kendi basina saglamdir - elle hazirlanmis girdiyle sinandi.')
        yaz('  ' + '!' * 70)
        g.close()
        return 7           # 7 = ARDISIK COKUS (girdi bos), 5 = kendi girdisi eksik
    if getattr(a, 'siparis', False):
        ciftler, kyol = siparistekiler(kok, hepsi=getattr(a, 'siparis_hepsi', False))
    else:
        ciftler, kyol = kurtarilanlar(kok)
    if not ciftler:
        yaz(u'  Kurtarma turunda esigi gecen YENI/DEGISMIS cift yok - dogrulanacak bir sey yok.')
        return 0
    yaz(u'  kaynak            : %s' % os.path.basename(kyol))
    yaz(u'  kaynak yolu       : %s' % kyol)
    yaz(u'  dogrulanacak cift : %d' % len(ciftler))
    if getattr(a, 'tumu', False):
        _sip = sum(1 for c in ciftler if c.get('sipariste'))
        yaz(u'  KIP: --tumu  (paneldeki BUTUN ciftler; siparis listesinde %d, '
            u'disinda %d)' % (_sip, len(ciftler) - _sip))
    if _ATLANAN:
        yaz(u'  ATLANAN (primer dizisi bulunamadi - protocol ciktisinda yok): %s'
            % ', '.join(_ATLANAN))
    for c in ciftler:
        yaz(u'     - %-42s %s' % (c['hedef'][:42], c['tur']))
    n_kume = sum(1 for _, d, _ in KUMELER if os.path.exists(os.path.join(kok, 'REFERANS_DB', d)))
    yaz('')
    # SURE BEYANI (2026-08-07): eski satir kume basina 240 sn TAHMIN ediyordu ve
    # "11,7 saat" yaziyordu. OLCUM bunu dogrulamadi: kontrol noktalari hazirken
    # yerel katman 11 kume x 16 cift icin 26 sn surdu (2026-08-07 olcumu).
    # Tahmin yerine, HAZIR KONTROL NOKTASI SAYISI bildirilir - okuyan kisi neyin
    # yeniden kosacagini gorsun.
    import hashlib as _h
    _imza = _h.md5('|'.join(sorted('%s>%s<%s' % (a['hedef'], a.get('F', ''),
                   a.get('R', '')) for a in ciftler)).encode()).hexdigest()[:10]
    _hazir = sum(1 for e, d, _a in KUMELER
                 if os.path.exists(os.path.join(KONTROL, 'yerel_%s_%s.pkl'
                                                % (re.sub(r'\W+', '_', e), _imza))))
    yaz(u'  KATMAN 2 (yerel): %d kume, %d tanesi kontrol noktasindan gelecek, '
        u'%d tanesi bastan taranacak.' % (n_kume, _hazir, max(0, n_kume - _hazir)))
    yaz(u'     (hepsi hazirken olculen sure: ~30 sn; bastan tarama kume basina '
        u'dakikalar surer)')
    yaz(u'  KATMAN 3 (MFEprimer): SILVA dahil 6 indeks. SILVA icin 16 ciftin '
        u'olculen toplam spec suresi ~85 sn + kanit kopyalama ~40 sn.')
    # KATMAN 4 kapsami: --ncbi-yalniz-siparis verilmisse NCBI yalniz siparis
    # listesindeki ciftlere kosar. Sure beyani da o sayidan hesaplanir.
    _ncbi_ciftler = ([c for c in ciftler if c.get('sipariste')]
                     if getattr(a, 'ncbi_yalniz_siparis', False) else list(ciftler))
    _n4 = len(_ncbi_ciftler)
    yaz(u'  KATMAN 4 (NCBI): cift basina OLCULEN sure ~75 sn (gonderim + '
        u'yoklama) + %d sn gonderim arasi bekleme -> %d cift icin ~%s.'
        % (PB_GONDERIM_ARASI, _n4,
           sure_metni(_n4 * 75 + max(0, _n4 - 1) * PB_GONDERIM_ARASI)))
    if getattr(a, 'ncbi_yalniz_siparis', False):
        _dis = [c['hedef'] for c in ciftler if not c.get('sipariste')]
        yaz(u'  KATMAN 4 KISITLI (--ncbi-yalniz-siparis): %d cift NCBI GORMEYECEK.'
            % len(_dis))
        yaz(u'     NCBI sutunu bu satirlarda BILINMIYOR kalir - "temiz" DEGIL.')
        for _d in _dis:
            yaz(u'       - %s' % _d)
    yaz(u'  Kesintiye dayaniklidir: her kume bitince kaydeder.')
    yaz('')

    yaz(u'--- KATMAN 2: YEREL VERITABANI TARAMASI ---')
    yerel = katman1_yerel(kok, ciftler, yaz, KONTROL, a.parc, a.kume_ust)

    # --- KATMAN 3: MFEprimer (BAGIMSIZ ARAC) ---
    mfe_sonuc = {}
    if not a.mfe_yok:
        yaz(u'--- KATMAN 3: MFEprimer (BAGIMSIZ ARAC) ---')
        yaz(u'  Ilk iki katman da BIZIM kodumuz ve ayni motoru kullaniyor; o motorda')
        yaz(u'  bir hata varsa ikisi de ayni yonde yanilir. Bu katman disaridan gelen')
        yaz(u'  bagimsiz bir araci ayni sorulara sokar.')
        import importlib.util as _u
        _sp = _u.spec_from_file_location(
            'mfe_katmani', os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        'mfe_katmani.py'))
        MK = _u.module_from_spec(_sp); _sp.loader.exec_module(MK)
        mfe = MK.mfe_bul(kok)
        if not mfe:
            yaz(u'  MFEprimer ikilisi bulunamadi (tools/mfeprimer) - katman ATLANDI.')
        else:
            yaz(u'  ikili: %s' % mfe)
            dur, spec = MK.spec_kos(kok, mfe, ciftler, CIKTI, yaz, KONTROL)
            yapi = MK.yapi_kos(kok, mfe, ciftler, CIKTI, yaz)
            # D-7: 'hedef disi N' sayisinin yanina KIMLIK dosyasi. Sayi tek
            # basina yanlis okunuyor (bkz. mfe_katmani.hedef_disi_kimlikleri).
            try:
                MK.hedef_disi_kimlikleri(CIKTI, ciftler, yaz)
            except Exception as _e:
                yaz(u'  UYARI: hedef disi kimlik dosyasi yazilamadi: %s' % _e)
            if dur.get('durum') == 'TAMAM':
                for c in ciftler:
                    mfe_sonuc[c['hedef']] = dict(
                        durum='TAMAM', hedef_disi=MK.hedef_disi_say(spec, c['hedef'], ciftler),
                        yapi=yapi, kullanilan_db=dur.get('kullanilan', []))
            else:
                yaz(u'  MFEprimer katmani sonuc vermedi: %s' % dur.get('sebep', ''))

    ncbi = {}
    if a.ncbi_yukle:
        yaz(u'--- KATMAN 4: NCBI (elle girilen sonuclar yukleniyor) ---')
        ncbi = ncbi_yukle(a.ncbi_yukle, yaz)
        yaz(u'  %d satir yuklendi' % len(ncbi))
    elif a.yalniz_yerel or a.ncbi == 'yok':
        yaz(u'--- KATMAN 4: NCBI ATLANDI (istek uzerine) ---')
    elif a.ncbi == 'oto':
        yaz(u'--- KATMAN 4: NCBI OTOMATIK (URL API) ---')
        yaz(u'  Not: blastn -remote KULLANILMIYOR (45 sn tavanini asiyor).')
        if _n4 != len(ciftler):
            yaz(u'  KAPSAM: %d/%d cift (siparis listesindekiler). Digerleri '
                u'yalniz katman 2+3 gordu.' % (_n4, len(ciftler)))
        ncbi = katman2_oto(_ncbi_ciftler, CIKTI, yaz, a.organizma,
                          haric_taxid=getattr(a, 'ncbi_haric_taxid', '') or '')
        if not any(v.get('durum', '').startswith('TAMAM') for v in ncbi.values()):
            yaz(u'  Otomatik yol sonuc vermedi - elle yol dosyalari uretiliyor.')
            katman2_elle_girdi(_ncbi_ciftler, CIKTI, yaz, a.organizma)
    else:
        yaz(u'--- KATMAN 4: NCBI ELLE (girdi ve sablon uretiliyor) ---')
        katman2_elle_girdi(_ncbi_ciftler, CIKTI, yaz, a.organizma)

    # D-12: MFEprimer'in ham (boya dayali) sayisi degil, TAKSONOMIK olarak
    # suzulmus klad_disi hukme girsin. Suzgec mfe_hedef_disi_kimlikleri.tsv'yi
    # ve screening/hedef_klad.tsv'yi okur; ikisinden biri yoksa BOS doner
    # ve o zaman ham sayi kullanilir ama rapor bunu acikca yazar.
    klad_sonuc = {}
    try:
        import verification.mfe_katmani as _MK2
    except ImportError:
        _MK2 = None
        try:
            sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
            import mfe_katmani as _MK2
        except ImportError:
            _MK2 = None
    if _MK2 is not None and hasattr(_MK2, 'klad_siniflandir'):
        try:
            klad_sonuc = _MK2.klad_siniflandir(kok, CIKTI, ciftler, TA_PANEL, yaz)
        except Exception as _e:
            yaz(u'  klad suzgeci koselemedi (%s: %s) - ham sayi kullanilacak'
                % (type(_e).__name__, _e))
            klad_sonuc = {}
    satirlar = birlestir(ciftler, yerel, ncbi, mfe_sonuc, klad_sonuc)
    for s in satirlar:
        s['_vurus'] = (yerel.get(s['hedef'], {}) or {}).get('vurus', [])
    raporla(CIKTI, satirlar, yaz)
    rc = cikti_denetle(yaz, 'D (DOGRULAMA)', [
        (os.path.join(CIKTI, 'dogrulama_uc_sutun.tsv'), 'dogrulama_uc_sutun.tsv')])
    g.close()
    return rc


if __name__ == '__main__':
    sys.exit(main() or 0)
