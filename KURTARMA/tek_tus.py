# -*- coding: utf-8 -*-
u"""TEK TUS - butun zinciri sirayla, dogru bagimlilik sirasinda, tek komutla kosar.

CAGIRAN : SABAH_TEK_TUS.bat  ve  KAPSAMLI_ARAMA.bat -> B tusu
CIKTI   : TEK_TUS_SONUC/00_SABAH_OZETI.md   (sabah bakilacak TEK dosya)
          TEK_TUS_SONUC/durum.json          (kontrol noktalari)
          TEK_TUS_SONUC/gunluk_<zaman>.log  (ekranin zaman damgali kopyasi)

BU DOSYA NEDEN VAR
------------------
Var olan A tusu (KURTARMA/tam_zincir.py) on asama kosuyor ama:
  * on kontrol yapmiyor - eksik dosya/arac/paketi kosunun ortasinda kesfediyor,
  * U ve P asamalarini hic tanimiyor (P'nin girdisi U'nun ciktisidir),
  * kontrol noktasi gecerliligini yalnizca "damga var mi" diye soruyor;
    GIRDI kontrol noktasindan yeniyse bayat damgayi okuyup asamayi atliyor.
    Bu tam olarak 2026-08-07'de yasandi (D-9, zehirli kontrol noktasi).
Bu dosya o uc boslugu kapatir. tam_zincir.py'ye DOKUNULMADI; o hala calisiyor.

TASARIM KURALLARI (hepsi bir hatadan ogrenildi, hicbiri suslemek icin degil)
---------------------------------------------------------------------------
1. ON KONTROL BIR KAPIDIR. Eksik varsa NE eksik oldugu yazilir ve DURULUR.
   Yarim kosu yok. (--on-kontrol-atla ile bilerek gecilebilir, ekrana yazilir.)
2. CIKIS KODU MASKELENMEZ. rc != 0 ise asama DUSTU'dur. Gecmiste tam_zincir.py
   T asamasinin 3 dondurdugunu goz ardi edip "BITTI" yazmisti; ozet yaniltici
   cikti. Burada rc VE cikti denetimi AYRI iki suzgectir, ikisi de gecmelidir.
3. KONTROL NOKTASI ANAHTARI BELIRLENIMCIDIR. md5 kullanilir; Python'un
   hash() fonksiyonu KULLANILMAZ (PYTHONHASHSEED yuzunden kosular arasinda
   degisir ve her kosu kontrol noktasini isikalar).
4. GIRDI KONTROL NOKTASINDAN YENIYSE KONTROL NOKTASI GECERSIZDIR. Asamanin
   kendi BETIGI de girdi sayilir: betik degistiyse asama yeniden kosar.
5. BAGIMLILIK YONLUDUR. Dusen bir asamanin bagimlilari KOSULMAZ, "atlandi
   (bagimli)" yazilir. Bagimsiz asamalar devam eder.
6. SURE TAHMINLERI OLCUMDUR. Her sayinin yaninda hangi dosyadan/kosudan
   geldigi yazilidir. Olculmemis asamaya sayi YAZILMAZ, "olculmedi" yazilir.
"""

import os, sys, io, csv, json, time, glob, signal, hashlib, argparse
import subprocess, threading

SURUM = u'1.0 (2026-08-08)'
CIKTI_KLASOR = 'TEK_TUS_SONUC'
CANLILIK_SN = 60          # uzun asamalarda kac saniyede bir canlilik isareti
LOG_YAZMA_ARALIGI = 2.0   # bagli klasore en cok bu sikligta yaziyoruz (D-11 kurali)


# ===========================================================================
#  0) KUCUK YARDIMCILAR
# ===========================================================================
def sn_metni(sn):
    sn = int(round(sn or 0))
    if sn < 90:
        return u'%d sn' % sn
    if sn < 5400:
        return u'%d dk %d sn' % (sn // 60, sn % 60)
    return u'%d sa %d dk' % (sn // 3600, (sn % 3600) // 60)


def boyut_metni(b):
    b = float(b)
    for birim in (u'B', u'KB', u'MB', u'GB'):
        if b < 1024 or birim == u'GB':
            return u'%d %s' % (b, birim) if birim == u'B' else u'%.1f %s' % (b, birim)
        b /= 1024.0


def veri_satiri_say(yol):
    u"""TSV'de BASLIK HARIC veri satiri sayar. Yorum ve bos satir sayilmaz."""
    if not os.path.exists(yol):
        return -1
    n = 0
    with io.open(yol, encoding='utf-8', errors='ignore') as fh:
        for s in fh:
            if s.startswith('#') or not s.strip():
                continue
            n += 1
    return max(n - 1, 0)


def dosya_parmak(kok, yol):
    u"""Bir dosyanin BELIRLENIMCI parmak izi: goreli yol | boy | mtime_ns.

    Neden icerik md5'i degil: REFERANS_DB icindeki tek bir dosya 1,5 GB.
    Her kosuda 30 GB okumak kabul edilemez. Boy+mtime_ns cifti dosya
    degisiminde degisir ve AYNI dosya icin kosular arasinda SABITTIR - bizim
    ihtiyacimiz olan iki ozellik de bu. Klasorler icin: dosya sayisi + en yeni
    mtime.
    """
    try:
        st = os.stat(yol)
    except OSError:
        return u'%s|YOK' % os.path.relpath(yol, kok).replace('\\', '/')
    g = os.path.relpath(yol, kok).replace('\\', '/')
    if os.path.isdir(yol):
        n, enyeni = 0, 0
        for kk, _, dd in os.walk(yol):
            for d in dd:
                n += 1
                try:
                    enyeni = max(enyeni, os.stat(os.path.join(kk, d)).st_mtime_ns)
                except OSError:
                    pass
        return u'%s|DIR|%d|%d' % (g, n, enyeni)
    return u'%s|%d|%d' % (g, st.st_size, st.st_mtime_ns)


def imzala(parcalar):
    u"""md5 - belirlenimci. Python'un hash() fonksiyonu KULLANILMAZ (kural 3)."""
    return hashlib.md5(u'\n'.join(parcalar).encode('utf-8')).hexdigest()[:16]


def en_yeni_mtime(yollar):
    en = 0.0
    for y in yollar:
        try:
            if os.path.isdir(y):
                for kk, _, dd in os.walk(y):
                    for d in dd:
                        en = max(en, os.stat(os.path.join(kk, d)).st_mtime)
            else:
                en = max(en, os.stat(y).st_mtime)
        except OSError:
            pass
    return en


def zaman_metni(t):
    return time.strftime('%d.%m %H:%M', time.localtime(t)) if t else u'yok'


def cozumle(kok, yol):
    return yol if os.path.isabs(yol) else os.path.join(kok, yol)


def yollari_ac(kok, listeler):
    u"""girdi/cikti listesini gercek yollara cevirir.

    'GLOB:desen' yazan bir oge desene uyan BUTUN dosyalara acilir. Buna ihtiyac
    var cunku uyelik ciktisinin adi tarihli: uyelik_yeniden_turetme_uyelik_
    20260803.tsv. Sabit yol yazsaydik yeni bir tarihte uretilen dosya
    goruinmezdi ve asama her kosuda bastan kosardi.
    """
    out = []
    for y in listeler:
        if y.startswith(u'GLOB:'):
            out.extend(sorted(glob.glob(os.path.join(kok, y[5:]))))
        else:
            out.append(cozumle(kok, y))
    return out


# ===========================================================================
#  1) ASAMA CIKTI DENETCILERI
#     Cikis kodu 0 olmasi YETMEZ. Her asama icin "beklenen cikti gercekten
#     uretildi mi ve icinde veri var mi" ayrica sorulur.
# ===========================================================================
def d_tsv_dolu(yollar, en_az=1):
    def f(kok, ayar, cikti_metni):
        eksik, bos, tamam = [], [], []
        for y in yollar:
            t = cozumle(kok, y)
            if not os.path.exists(t):
                eksik.append(y)
            elif y.endswith('.tsv') or y.endswith('.csv'):
                n = veri_satiri_say(t)
                (tamam if n >= en_az else bos).append(u'%s (%d satir)' % (y, n))
            elif os.path.getsize(t) < 40:
                bos.append(u'%s (%d bayt)' % (y, os.path.getsize(t)))
            else:
                tamam.append(u'%s (%s)' % (y, boyut_metni(os.path.getsize(t))))
        if eksik:
            return False, u'beklenen cikti URETILMEDI: %s' % u', '.join(eksik)
        if bos:
            return False, u'cikti BOS / veri satiri yok: %s' % u', '.join(bos)
        return True, u'cikti dogrulandi: %s' % u', '.join(tamam)
    return f


def d_selftest(kok, ayar, cikti_metni):
    if u'TUM SINAMALAR GECTI' in (cikti_metni or u''):
        return True, u'butun selftestler gecti'
    return False, (u'"TUM SINAMALAR GECTI" satiri cikmadi. Kod kendini '
                   u'dogrulayamadi; olcume girilmez.')


def d_hizli_test(kok, ayar, cikti_metni):
    y = os.path.join(kok, 'HIZLI_TEST', 'HIZLI_TEST_RAPORU.md')
    if not os.path.exists(y):
        return False, u'HIZLI_TEST/HIZLI_TEST_RAPORU.md uretilmedi'
    m = io.open(y, encoding='utf-8', errors='ignore').read()
    if u'ZINCIR TUTARSIZ' in m:
        return False, (u'rapor ZINCIR TUTARSIZ diyor. Uzun kosuya GIRILMEZ - '
                       u'HIZLI_TEST/HIZLI_TEST_RAPORU.md dosyasina bakin.')
    if u'ZINCIR TUTARLI' in m:
        return True, u'rapor ZINCIR TUTARLI diyor'
    return False, u'raporda ne TUTARLI ne TUTARSIZ karari var - bicim beklenmedik'


def d_uyelik(kok, ayar, cikti_metni):
    g = sorted(glob.glob(os.path.join(kok, 'uyelik_yeniden_turetme_uyelik_*.tsv')))
    if not g:
        return False, u'uyelik_yeniden_turetme_uyelik_*.tsv uretilmedi'
    n = veri_satiri_say(g[-1])
    if n < 1:
        return False, u'%s BOS' % os.path.basename(g[-1])
    return True, u'%s (%d satir)' % (os.path.basename(g[-1]), n)


def d_yok(kok, ayar, cikti_metni):
    return True, u'bu asama dosya uretmez, yalniz ekrana yazar'


# ===========================================================================
#  2) ASAMA CIZELGESI  -  BAGIMLILIK YONLU VE SIRA ZORUNLUDUR
#
#  Alanlar:
#    kod          menude ve ozette gorunecek harf
#    ad           tek satirlik aciklama
#    betik        varligi on kontrolde aranan dosya (None ise dis arac)
#    argv(kok,a)  calistirilacak komut listelerinin listesi
#    girdi[]      bu asamanin OKUDUGU dosyalar. Bunlardan biri asamanin
#                 ciktisindan YENIYSE kontrol noktasi GECERSIZDIR.
#    cikti[]      bu asamanin URETTIGI dosyalar (cikti denetimi ve tazelik icin)
#    bagimli[]    once bitmesi gereken asama kodlari
#    sure_sn      OLCULMUS sure (saniye). None ise "olculmedi" yazilir.
#    kaynak       o sayinin hangi dosyadan/kosudan geldigi
#    denet        cikti denetcisi
#    kraken       True: ALI klasorundeki araca ihtiyac duyar
#    hep_kos      True: hizli ve yan etkisiz; kontrol noktasi tutulmaz
#
#  SIRA NEDEN BOYLE (kod okunarak cikarildi, tahmin degil):
#    U -> P : tek_protokol_olc.py uyeligi uyelik_yeniden_turetme_uyelik_*.tsv
#             dosyasindan alir (KAPSAMLI_ARAMA.bat satir 1092 de bunu arar).
#    P -> K : kurtarma_turu.py girdisi TEK_PROTOKOL_SONUC/panel_tek_protokol.tsv
#             (KAPSAMLI_ARAMA.bat satir 946 bunu on kosul olarak dogruluyor).
#    P -> D : dogrulama_turu.siparistekiler() TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv
#             okur (dogrulama_turu.py satir 150-155).
#    K -> D : D, K'nin kurtardigi ciftleri de sinar (kurtarma_satirlari.tsv).
#    I -> G : tum_kutu_kimlikleri.py onbellegi I ile PAYLASIR; I once kosarsa
#             ayni kutu iki kez taranmaz (KAPSAMLI_ARAMA.bat satir 686).
#    W -> X : esik taramasi once ortam denetiminden gecer.
#    X -> Z : tablo, esik taramasinin ciktisini okur, yeni olcum yapmaz.
#    H once : H bir GERILEME KAPISIDIR - onceki referans kosuya karsi sinar.
#             Bu yuzden P'den ONCE kosar; P'nin yeni ciktisini beklemez.
# ===========================================================================
def _py(*a):
    return [sys.executable] + list(a)


KRAKEN_KOMUTLARI = {'W': ('durum', 'vt-ara', 'vt-kimlik', 'sinav'),
                    'X': ('esik',),
                    'Z': ('tablo',)}


def ASAMALAR(ayar):
    ncbi = ayar.get('ncbi', 'oto')
    org = ayar.get('organizma', '')
    karac = ayar.get('karac')

    def d_argv(kok, a):
        arg = [os.path.join('KURTARMA', 'dogrulama_turu.py'), '--kok', '.']
        arg += ['--tumu'] if a.get('tumu', True) else ['--siparis']
        if ncbi == 'yok':
            arg += ['--yalniz-yerel']
        else:
            arg += ['--ncbi', ncbi]
            if org:
                arg += ['--organizma', org]
        return [_py(*arg)]

    def kraken_argv(kod):
        return lambda kok, a: [['bash', a['karac'], k] for k in KRAKEN_KOMUTLARI[kod]]

    L = [
        dict(kod='8', ad=u'KENDINI SINA - kod kendini dogrular, olcum yapmaz',
             grup=u'Grup 4', betik='KAPSAMLI_ARAMA/__main__.py',
             argv=lambda kok, a: [_py('-m', 'KAPSAMLI_ARAMA', '--sina')],
             girdi=['KAPSAMLI_ARAMA'], cikti=[], bagimli=[],
             sure_sn=4.6, kaynak=u'TAM_ZINCIR_SONUC/durum.json, 2026-08-06 kosusu',
             denet=d_selftest, hep_kos=True),

        # 2026-08-10 EKLENDI. Bu projede hatalarin cogu olcumde degil olcumun
        # DAYANDIGI tabloda cikti ve hicbiri kosuyu dusurmedi; ancak birisi
        # "baska hata var mi" diye sordugunda bulundu. Sormaya bagli denetim
        # denetim degildir. Bu asama o soruyu HER kosuda kendisi sorar.
        # Olcum yapmaz, dosya degistirmez; bu yuzden hep_kos=True.
        dict(kod='N', ad=u'DENETIM - tablolar, referanslar ve muhurler her kosuda bakilir',
             grup=u'Grup 4', betik='KURTARMA/hepsini_denetle.py',
             argv=lambda kok, a: [_py(os.path.join('KURTARMA', 'hepsini_denetle.py'),
                                      '--kok', '.')],
             girdi=['KURTARMA/hepsini_denetle.py',
                    'KAPSAMLI_ARAMA/hedef_taxid.tsv',
                    'TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv'],
             cikti=['TEK_TUS_SONUC/DENETIM_RAPORU.md'], bagimli=[],
             sure_sn=90.0, kaynak=u'olculdu 2026-08-10 (yerel 2 sn + NCBI kapsama ~85 sn)',
             hep_kos=True, danisma=True),

        dict(kod='H', ad=u'HIZLI TUTARLILIK TESTI - uzun kosudan ONCE gerileme kapisi',
             grup=u'Grup 4', betik='KURTARMA/hizli_tutarlilik_testi.py',
             argv=lambda kok, a: [_py(os.path.join('KURTARMA', 'hizli_tutarlilik_testi.py'),
                                      '--kok', '.')],
             girdi=['KURTARMA/hizli_tutarlilik_testi.py',
                    'TEK_PROTOKOL_SONUC/panel_tek_protokol.tsv'],
             cikti=['HIZLI_TEST/HIZLI_TEST_RAPORU.md'], bagimli=[],
             sure_sn=122.0, kaynak=u'TAM_ZINCIR_SONUC/durum.json, 2026-08-06 kosusu',
             denet=d_hizli_test),

        dict(kod='E', ad=u'VERITABANI ERISIM DOGRULAMASI - her VT gercekten okunuyor mu',
             grup=u'Grup 2', betik='KURTARMA/erisim_dogrulama.py',
             argv=lambda kok, a: [_py(os.path.join('KURTARMA', 'erisim_dogrulama.py'),
                                      '--kok', '.')],
             girdi=['KURTARMA/erisim_dogrulama.py'],
             cikti=['ERISIM_SONUC/erisim_dogrulama.tsv'], bagimli=[],
             sure_sn=None,
             kaynak=u'OLCULMEDI - bu asama hic kosulmamis (ERISIM_SONUC klasoru yok)',
             denet=d_tsv_dolu(['ERISIM_SONUC/erisim_dogrulama.tsv'])),

        dict(kod='U', ad=u'UYELIGI OLCULEN KIMLIKTEN YENIDEN TURET',
             grup=u'Grup 4', betik='UYELIK_YENIDEN/uyelik_yeniden_turet.py',
             argv=lambda kok, a: [_py(os.path.join('UYELIK_YENIDEN',
                                                   'uyelik_yeniden_turet.py'),
                                      '--kok', '.')],
             girdi=['UYELIK_YENIDEN/uyelik_yeniden_turet.py', 'consensus sequences'],
             cikti=['GLOB:uyelik_yeniden_turetme_uyelik_*.tsv'], bagimli=[],
             sure_sn=None,
             kaynak=u'OLCULMEDI - menudeki "1-3 saat" bir tahmindir, olcum degil',
             denet=d_uyelik),

        dict(kod='P', ad=u'TEK PROTOKOL - panelin tamami TEK kuralla olculur',
             grup=u'Grup 1', betik='TEK_PROTOKOL/tek_protokol_olc.py',
             argv=lambda kok, a: [_py(os.path.join('TEK_PROTOKOL', 'tek_protokol_olc.py'),
                                      '--kok', '.')],
             # 2026-08-09: ciftler.tsv eklendi, gerekce D asamasindaki notta.
             girdi=['TEK_PROTOKOL/tek_protokol_olc.py', 'primer_final',
                    'KAPSAMLI_ARAMA/ciftler.tsv',
                    'GLOB:uyelik_yeniden_turetme_uyelik_*.tsv'],
             cikti=['TEK_PROTOKOL_SONUC/panel_tek_protokol.tsv',
                    'TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv'],
             bagimli=['U'],
             sure_sn=36.0,
             kaynak=u'KURTARMA/tam_zincir.py yorumu ("Olculen: P 36 sn"); '
                    u'TUM_KOSU_SONUC/durum.json sicak kosuda 9,2 sn',
             denet=d_tsv_dolu(['TEK_PROTOKOL_SONUC/panel_tek_protokol.tsv',
                               'TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv'])),

        dict(kod='K', ad=u'KURTARMA - esik alti satirlar dort yolla kurtarilir',
             grup=u'Grup 1', betik='KURTARMA/kurtarma_turu.py',
             argv=lambda kok, a: [_py(os.path.join('KURTARMA', 'kurtarma_turu.py'),
                                      '--kok', '.')],
             # 2026-08-09: ciftler.tsv eklendi, gerekce D asamasindaki notta.
             girdi=['KURTARMA/kurtarma_turu.py',
                    'KAPSAMLI_ARAMA/ciftler.tsv',
                    'TEK_PROTOKOL_SONUC/panel_tek_protokol.tsv'],
             cikti=['KURTARMA_SONUC/kurtarma_satirlari.tsv'], bagimli=['P'],
             sure_sn=300.0, kaynak=u'KURTARMA/tam_zincir.py yorumu ("K 5 dk")',
             denet=d_tsv_dolu(['KURTARMA_SONUC/kurtarma_satirlari.tsv'])),

        dict(kod='D', ad=u'DOGRULAMA - paneldeki ciftler dort kanit katmaniyla sinanir',
             grup=u'Grup 1', betik='KURTARMA/dogrulama_turu.py',
             argv=d_argv,
             # 2026-08-09 DUZELTME (girdi takibi): asamalarin HICBIRI
             # KAPSAMLI_ARAMA/ciftler.tsv dosyasini girdi saymiyordu. O dosya
             # panelin PRIMER DIZILERINI tutuyor. Yani butun panelin dizileri
             # degistirilse bile zincir "girdi degismemis" deyip her asamayi
             # atlardi. 09.08'de iki cift degistirildi ve D yine atlanacakti;
             # dogruladigi sey eski cift olurdu. Kok SIPARIS_LISTESI.tsv de
             # eklendi, cunku D siparis satirlarini oradan okuyor.
             girdi=['KURTARMA/dogrulama_turu.py', 'KURTARMA/mfe_katmani.py',
                    'KAPSAMLI_ARAMA/hedef_klad.tsv',
                    'KAPSAMLI_ARAMA/ciftler.tsv',
                    'SIPARIS_LISTESI.tsv',
                    'TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv',
                    'KURTARMA_SONUC/kurtarma_satirlari.tsv',
                    'REFERANS_DB/SILVA_138.2_SSURef_NR99.fasta.primerqc.bin'],
             cikti=['DOGRULAMA_SONUC/dogrulama_uc_sutun.tsv'], bagimli=['P', 'K'],
             sure_sn=6900.0,
             kaynak=u'TUM_CIFTLER_DEVIR_2026-08-07: katman2 81 dk (soguk, 22 cift) '
                    u'+ katman3 2,5 dk + katman4 NCBI 31 dk = ~1 sa 55 dk',
             denet=d_tsv_dolu(['DOGRULAMA_SONUC/dogrulama_uc_sutun.tsv'])),

        dict(kod='I', ad=u'KIMLIK DOGRULAMA - hocaya giden iddialar bagimsiz sinanir',
             grup=u'Grup 2', betik='KURTARMA/kimlik_dogrulama.py',
             argv=lambda kok, a: [_py(os.path.join('KURTARMA', 'kimlik_dogrulama.py'),
                                      '--kok', '.')],
             girdi=['KURTARMA/kimlik_dogrulama.py'],
             cikti=['KIMLIK_SONUC/kimlik_iddialari.tsv'], bagimli=[],
             sure_sn=12007.0,
             kaynak=u'TAM_ZINCIR_SONUC/durum.json, 2026-08-06 (3 sa 20 dk)',
             denet=d_tsv_dolu(['KIMLIK_SONUC/kimlik_iddialari.tsv'])),

        dict(kod='G', ad=u'TUM KUTU KIMLIKLERI - panele giren HER kutu dogrulanir',
             grup=u'Grup 2', betik='KURTARMA/tum_kutu_kimlikleri.py',
             argv=lambda kok, a: [_py(os.path.join('KURTARMA', 'tum_kutu_kimlikleri.py'),
                                      '--kok', '.', '--nt', 'yok')],
             girdi=['KURTARMA/tum_kutu_kimlikleri.py'],
             cikti=['TUM_KIMLIK_SONUC/tum_kutu_kimlikleri.tsv'], bagimli=['I'],
             sure_sn=17038.0,
             kaynak=u'TAM_ZINCIR_SONUC/durum.json, 2026-08-06 (4 sa 44 dk)',
             denet=d_tsv_dolu(['TUM_KIMLIK_SONUC/tum_kutu_kimlikleri.tsv'])),

        dict(kod='W', ad=u'KRAKEN2 ORTAM DENETIMI - kurulu mu, hangi veritabani',
             grup=u'Grup 3', betik=None, argv=kraken_argv('W'),
             girdi=[], cikti=[], bagimli=[],
             sure_sn=81.0, kaynak=u'TAM_ZINCIR_SONUC/durum.json, 2026-08-06',
             denet=d_yok, kraken=True, hep_kos=True),

        dict(kod='X', ad=u'KRAKEN GUVEN ESIGI TARAMASI',
             grup=u'Grup 3', betik=None, argv=kraken_argv('X'),
             girdi=[], cikti=[], bagimli=['W'],
             sure_sn=6916.0,
             kaynak=u'TAM_ZINCIR_SONUC/durum.json, 2026-08-06 (1 sa 55 dk)',
             denet=d_yok, kraken=True),

        dict(kod='Z', ad=u'DORT SUTUNLU KARSILASTIRMA TABLOSU',
             grup=u'Grup 3', betik=None, argv=kraken_argv('Z'),
             girdi=[], cikti=[], bagimli=['X'],
             sure_sn=4.8, kaynak=u'TAM_ZINCIR_SONUC/durum.json, 2026-08-06',
             denet=d_yok, kraken=True),

        dict(kod='S', ad=u'OZETI YENILE - olcum yapmaz, mevcut dosyalari okur',
             grup=u'Grup 4', betik='KAPSAMLI_ARAMA/__main__.py',
             argv=lambda kok, a: [_py('-m', 'KAPSAMLI_ARAMA', '--mod', 'ozet')],
             girdi=[], cikti=['KAPSAMLI_ARAMA_SONUC/00_OZET_HEPSI.md'], bagimli=[],
             sure_sn=0.7, kaynak=u'TAM_ZINCIR_SONUC/durum.json, 2026-08-06',
             denet=d_tsv_dolu(['KAPSAMLI_ARAMA_SONUC/00_OZET_HEPSI.md']),
             hep_kos=True),
    ]
    for a in L:
        a.setdefault('kraken', False)
        a.setdefault('hep_kos', False)
        a.setdefault('danisma', False)
        a['karac'] = karac
    return L


# ===========================================================================
#  3) ON KONTROL  -  EKSIK VARSA DURUR
# ===========================================================================
def on_kontrol(kok, ayar, yaz):
    u"""Kosuya girmeden once her sart tek tek olculur.

    ZORUNLU dusen tek bir madde varsa kosu BASLAMAZ. Bu bilerek boyle:
    yarim kosu, sabah okunacak yaniltici bir ozet uretir.
    """
    zorunlu_dusen, uyari, sayac = [], [], [0]
    yaz(u'')
    yaz(u'=' * 78)
    yaz(u'  ON KONTROL - kosmadan once her sart olculuyor')
    yaz(u'=' * 78)

    def satir(ad, ok, ayrinti, zorunlu=True):
        sayac[0] += 1
        yaz(u'  [%s] %-46s %s' % (u' OK  ' if ok else u'EKSIK', ad, ayrinti))
        if not ok:
            (zorunlu_dusen if zorunlu else uyari).append(u'%s -> %s' % (ad, ayrinti))

    # --- 1) Python surumu ve paketler -------------------------------------
    pv = '%d.%d.%d' % sys.version_info[:3]
    satir(u'Python 3.8+', sys.version_info[:2] >= (3, 8), pv)
    for paket, kurulum in (('numpy', 'numpy'), ('primer3', 'primer3-py')):
        try:
            m = __import__(paket)
            satir(u'python paketi: %s' % paket, True, getattr(m, '__version__', 'var'))
        except Exception:
            satir(u'python paketi: %s' % paket, False,
                  u'YOK - kurulum: pip3 install %s --break-system-packages' % kurulum)

    # --- 2) Betikler -------------------------------------------------------
    for a in ASAMALAR(ayar):
        if not a.get('betik'):
            continue
        t = os.path.join(kok, a['betik'])
        satir(u'betik: %s' % a['betik'], os.path.exists(t),
              boyut_metni(os.path.getsize(t)) if os.path.exists(t) else u'BULUNAMADI')

    # --- 3) Veri klasorleri ------------------------------------------------
    for d, zor in (('fastq files', True), ('consensus sequences', True),
                   ('primer_final', True), ('REFERANS_DB', True),
                   ('KAPSAMLI_ARAMA', True), ('TEK_PROTOKOL', True),
                   ('SON_ETAP_betikleri', True), ('MADDE123_betikleri', True),
                   ('DUZELTME_betikleri', True), ('UYELIK_YENIDEN', True),
                   ('konsensus_kanonik', False)):
        t = os.path.join(kok, d)
        var = os.path.isdir(t)
        satir(u'klasor: %s' % d, var,
              (u'%d oge' % len(os.listdir(t))) if var else u'BULUNAMADI', zor)

    # --- 4) MFEprimer ikilisi ve indeksleri --------------------------------
    mfe = None
    for aday in (os.path.join(kok, 'ARACLAR', 'mfeprimer'), 'mfeprimer'):
        if os.path.exists(aday):
            mfe = aday + (u'' if os.access(aday, os.X_OK) else u'  (CALISTIRMA IZNI YOK)')
            if not os.access(aday, os.X_OK):
                mfe = None
                uyari.append(u'ARACLAR/mfeprimer var ama calistirilabilir degil: '
                             u'chmod +x ARACLAR/mfeprimer')
            break
        try:
            from shutil import which
            w = which(aday)
            if w:
                mfe = w
                break
        except Exception:
            pass
    satir(u'MFEprimer ikilisi', bool(mfe),
          mfe or u'ARACLAR/mfeprimer YOK ya da calistirilabilir degil')

    MFE_IX = ['archaea.16S.fna', 'bacteria.16S.fna', 'fungi.ITS.fna',
              'fungi.28SrRNA.fna', 'fungi.18SrRNA.fna',
              'SILVA_138.2_SSURef_NR99.fasta']
    for f in MFE_IX:
        p = os.path.join(kok, 'REFERANS_DB', f + '.primerqc.bin')
        satir(u'MFE indeksi: %s' % f, os.path.exists(p),
              boyut_metni(os.path.getsize(p)) if os.path.exists(p)
              else u'YOK - kurulum: mfeprimer index -i REFERANS_DB/%s' % f)

    # --- 5) Katman 2'nin taradigi 11 kume ----------------------------------
    KUMELER = ['SILVA_138.2_SSURef_NR99.fasta', 'SILVA_138.2_LSURef_NR99.fasta',
               'UNITE_ITS.fasta', 'PR2_SSU_taxo_long.fasta',
               'ROD_v1.2_operon_variants.fasta', 'bacteria.16S.fna',
               'archaea.16S.fna', 'fungi.ITS.fna', 'fungi.28SrRNA.fna',
               'fungi.18SrRNA.fna', 'ref_all2.fna']
    eksik_kume = [f for f in KUMELER
                  if not os.path.exists(os.path.join(kok, 'REFERANS_DB', f))]
    satir(u'katman 2 veritabanlari (11 kume)', not eksik_kume,
          u'11/11 yerinde' if not eksik_kume else u'EKSIK: %s' % u', '.join(eksik_kume))

    # --- 6) SILVA SSU RNA/DNA kapisi ---------------------------------------
    # Gecmiste SILVA'nin RNA alfabesi (U) indeksi bozmus ve butun baglanmalar
    # 0/0 gelmisti. Ilk birkac bin satirda U/T sayarak DNA oldugunu dogruluyoruz.
    sp = os.path.join(kok, 'REFERANS_DB', 'SILVA_138.2_SSURef_NR99.fasta')
    if os.path.exists(sp):
        try:
            with io.open(sp, encoding='utf-8', errors='ignore') as fh:
                ornek = u''.join(s for s in (fh.readline() for _ in range(4000))
                                 if not s.startswith(u'>'))
            nu, nt = ornek.count(u'U'), ornek.count(u'T')
            satir(u'SILVA SSU alfabesi DNA mi (U=0 olmali)', nu == 0 and nt > 0,
                  u'U=%d  T=%d' % (nu, nt))
        except Exception as e:
            satir(u'SILVA SSU alfabesi', False, u'okunamadi: %s' % e, zorunlu=False)

    # --- 7) Disk alani -----------------------------------------------------
    try:
        st = os.statvfs(kok)
        bos_gb = st.f_bavail * st.f_frsize / 1073741824.0
        satir(u'bagli klasorde bos disk (>= 5 GB)', bos_gb >= 5, u'%.1f GB bos' % bos_gb)
    except Exception as e:
        satir(u'bagli klasorde bos disk', False, u'olculemedi: %s' % e, zorunlu=False)
    try:
        st = os.statvfs('/tmp')
        tmp_gb = st.f_bavail * st.f_frsize / 1073741824.0
        # MFEprimer kaniti YEREL diskte uretilip toplu kopyalaniyor (D-11).
        # Olculen en buyuk kanit 282 MB; 2 GB rahat pay birakir.
        satir(u'yerel /tmp bos alani (>= 2 GB, D-11 icin)', tmp_gb >= 2,
              u'%.1f GB bos' % tmp_gb)
    except Exception as e:
        satir(u'yerel /tmp bos alani', False, u'olculemedi: %s' % e, zorunlu=False)

    # --- 8) WSL bellegi ----------------------------------------------------
    try:
        mem = {}
        for s in io.open('/proc/meminfo', encoding='utf-8'):
            p = s.split(':')
            if len(p) == 2:
                mem[p[0]] = int(p[1].strip().split()[0])
        top_gb = mem.get('MemTotal', 0) / 1048576.0
        # Olculen tepe bellek: mfeprimer indeksleme 6,06 GB
        # (REFERANS_DB/SILVA_138.2_SSURef_NR99.fasta.BOZUK_KANIT.txt).
        # Tarama tarafi bundan dusuk; 4 GB alt sinir, 8 GB rahat.
        satir(u'WSL toplam bellek (>= 4 GB)', top_gb >= 3.8, u'%.1f GB' % top_gb)
        if 3.8 <= top_gb < 7.5:
            uyari.append(u'WSL bellegi %.1f GB - Kraken adimlari (X) buyuk '
                         u'veritabaninda bogulabilir. KILITLENME_COZUMU.md, '
                         u'KAPSAMLI_ARAMA.bat -> M tusu.' % top_gb)
    except Exception as e:
        satir(u'WSL bellegi', False, u'olculemedi: %s' % e, zorunlu=False)

    # --- 9) Kraken araci (istege bagli) ------------------------------------
    satir(u'Kraken araci (../ALI/WSL/130_KRAKEN_ARAC.sh)', bool(ayar.get('karac')),
          ayar.get('karac') or ayar.get('kraken_sebep', u'yok - W/X/Z ATLANACAK'),
          zorunlu=False)

    # --- 10) Ag (NCBI katmani icin) ----------------------------------------
    if ayar.get('ncbi') == 'oto':
        ok, ayr = _ag_dene()
        satir(u'NCBI erisimi (--ncbi oto secildi)', ok, ayr, zorunlu=False)
        if not ok:
            uyari.append(u'NCBI\'ye ulasilamadi. D asamasinin 4. katmani '
                         u'"BILINMIYOR" kalir; yerel katmanlar yine kosar.')

    # --- 11) Yazma izni ----------------------------------------------------
    # Sorulan tek sey YAZABILIYOR MUYUZ. Silme yetkisi AYRI bir sorudur ve bu
    # zincir hicbir sey silmez; silinemeyen deneme dosyasi kosuyu durdurmaz.
    try:
        t = os.path.join(kok, CIKTI_KLASOR)
        os.makedirs(t, exist_ok=True)
        p = os.path.join(t, '_yazma_denemesi.txt')
        with io.open(p, 'w', encoding='utf-8') as fh:
            fh.write(u'yazma denemesi %s\n' % time.strftime('%Y-%m-%d %H:%M:%S'))
        try:
            os.remove(p)
            ek = u'yazilabiliyor'
        except Exception:
            ek = u'yazilabiliyor (deneme dosyasi silinemedi - silme yetkisi yok, ' \
                 u'zincir zaten hicbir sey silmiyor)'
        satir(u'bagli klasore yazma izni', True, u'%s/ %s' % (CIKTI_KLASOR, ek))
    except Exception as e:
        satir(u'bagli klasore yazma izni', False, u'YAZILAMIYOR: %s' % e)

    yaz(u'')
    if uyari:
        yaz(u'  UYARILAR (kosu durmaz, ama bilerek gecmis olun):')
        for u_ in uyari:
            yaz(u'    ! %s' % u_)
        yaz(u'')
    if zorunlu_dusen:
        yaz(u'=' * 78)
        yaz(u'  ON KONTROL DUSTU - KOSU BASLAMADI')
        yaz(u'=' * 78)
        yaz(u'  Eksik olan %d madde:' % len(zorunlu_dusen))
        for z in zorunlu_dusen:
            yaz(u'    * %s' % z)
        yaz(u'')
        yaz(u'  Hicbir asama kosulmadi. Yarim kosu, sabah okunacak yaniltici bir')
        yaz(u'  ozet uretirdi; bu yuzden BILEREK durduruldu.')
        yaz(u'  Eksikleri giderip ayni tusa yeniden basin.')
        yaz(u'=' * 78)
        return False, zorunlu_dusen, uyari
    yaz(u'  ON KONTROL GECTI - %d madde denetlendi, zorunlu eksik yok.' % sayac[0])
    return True, zorunlu_dusen, uyari


def _ag_dene():
    try:
        import urllib.request
        t0 = time.time()
        r = urllib.request.urlopen(
            'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi', timeout=12)
        r.read(200)
        r.close()
        return True, u'ulasildi (%.1f sn)' % (time.time() - t0)
    except Exception as e:
        return False, u'ulasilamadi: %s' % str(e)[:70]


# ===========================================================================
#  4) KONTROL NOKTASI - GECERLILIK
# ===========================================================================
def girdi_imzasi(kok, a):
    u"""Asamanin OKUDUGU her seyin belirlenimci parmak izi.

    Asamanin KENDI BETIGI de girdidir: betik degistiyse eski sonuc bayattir.
    Bu, 2026-08-07'de yasanan 'zehirli kontrol noktasi' hatasinin ayni
    sinifini kapatir - orada indeks yenilenmis ama kontrol noktasi eski
    sifirlari geri okumustu.
    """
    p = [u'surum=%s' % SURUM, u'kod=%s' % a['kod']]
    if a.get('betik'):
        p.append(dosya_parmak(kok, os.path.join(kok, a['betik'])))
    acilan = yollari_ac(kok, a['girdi'])
    if not acilan and a['girdi']:
        p.append(u'GIRDI_BULUNAMADI')
    for t in acilan:
        p.append(dosya_parmak(kok, t) if os.path.exists(t)
                 else u'%s|YOK' % os.path.relpath(t, kok).replace('\\', '/'))
    return imzala(p)


def _girdi_yollari(kok, a):
    y = yollari_ac(kok, a['girdi'])
    if a.get('betik'):
        y.append(os.path.join(kok, a['betik']))
    return y


def kontrol_noktasi_gecerli(kok, a, durum):
    u"""(atlanabilir_mi, sebep) - NEDEN atlandigi/atlanmadigi acikca yazilir.

    Dort suzgec, hepsi gecmeli:
      1. damga var ve 'bitti'
      2. girdi imzasi ayni (girdi degismemis)
      3. beklenen ciktilarin hepsi var ve dolu
      4. en yeni cikti, en yeni girdiden YENI (bayat cikti okunmaz)
    """
    if a.get('hep_kos'):
        return False, u'bu asama her kosuda yeniden kosar (hizli ve yan etkisiz)'

    d = durum.get(a['kod'], {})
    imza = girdi_imzasi(kok, a)
    ciktilar = yollari_ac(kok, a['cikti'])
    girdiler = _girdi_yollari(kok, a)
    cy, gy = en_yeni_mtime(ciktilar), en_yeni_mtime(girdiler)

    if d.get('durum') == 'bitti':
        if d.get('imza') != imza:
            return False, (u'kontrol noktasi BAYAT - girdi imzasi degismis '
                           u'(damga %s, simdi %s). Yeniden kosulacak.'
                           % (str(d.get('imza'))[:10], imza[:10]))
        if ciktilar:
            ok, mesaj = a['denet'](kok, {}, u'')
            if not ok:
                return False, u'damga var ama cikti denetimi dustu: %s' % mesaj
            if gy > cy:
                return False, (u'kontrol noktasi BAYAT - girdi ciktidan YENI '
                               u'(girdi %s > cikti %s). Yeniden kosulacak.'
                               % (zaman_metni(gy), zaman_metni(cy)))
        return True, (u'onceki kosuda bitmisti (%s), girdi degismemis'
                      % sn_metni(d.get('sure', 0)))

    # Damga yok. Diskte HAZIR ve TAZE bir cikti var mi?
    # Bu, bu betik yazilmadan once uretilmis sonuclari tanimak icindir; onlari
    # korkudan yeniden kosmak saatler israf ederdi. Ama "var" YETMEZ: "girdiden
    # yeni" de sart, yoksa bayat sonucu tazeymis gibi kabul ederiz.
    if ciktilar:
        ok, mesaj = a['denet'](kok, {}, u'')
        if ok:
            if cy >= gy and cy > 0:
                return True, (u'damga yok ama cikti diskte TAZE (cikti %s >= en yeni '
                              u'girdi %s). %s' % (zaman_metni(cy), zaman_metni(gy), mesaj))
            return False, (u'diskte cikti var ama BAYAT - girdi %s, cikti %s. '
                           u'Yeniden kosulacak.' % (zaman_metni(gy), zaman_metni(cy)))
        return False, u'daha once kosulmamis (%s)' % mesaj
    return False, u'daha once kosulmamis'


# ===========================================================================
#  5) TEK ASAMAYI KOS
# ===========================================================================
class Kesildi(Exception):
    pass


KESME = {'var': False}


def asama_kos(kok, a, ayar, yaz):
    u"""Bir asamanin butun komutlarini sirayla kosar.

    * Cikis kodu okunur ve MASKELENMEZ.
    * Ilk sifir disi kodda o asamanin kalan komutlari kosulmaz.
    * Uzun sessizliklerde canlilik isareti basilir.
    """
    t0 = time.time()
    rc, son_cikti = 0, u''
    komutlar = a['argv'](kok, ayar)
    cevre = dict(os.environ)
    cevre.setdefault('IPLIK', '3')
    cevre.setdefault('ZORLA_MMAP', '1')
    cevre['PYTHONUNBUFFERED'] = '1'
    if ayar.get('vt_a'):
        cevre['VT_A'] = ayar['vt_a']
    # Kraken araci kendi klasorunden cagrilmali (yan dosyalarini goreli arar).
    calisma = os.path.dirname(a['karac']) if (a.get('kraken') and a.get('karac')) else kok

    for i, argv in enumerate(komutlar, 1):
        if KESME['var']:
            raise Kesildi()
        yaz(u'   $ [%d/%d] %s' % (i, len(komutlar),
                                  u' '.join(os.path.basename(x) if os.sep in x else x
                                            for x in argv)))
        try:
            p = subprocess.Popen(argv, cwd=calisma, env=cevre,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        except Exception as e:
            rc, son_cikti = 127, u'komut baslatilamadi: %s' % e
            yaz(u'   HATA: %s' % son_cikti)
            break

        satirlar = []
        son = {'t': time.time()}
        dur = threading.Event()

        def canlilik():
            # Uzun asamalarda "asildi mi?" sorusunu ekranda cevaplar.
            while not dur.wait(CANLILIK_SN):
                yaz(u'   ... [canlilik] %s | asamada %s | son ciktidan bu yana %s'
                    % (a['kod'], sn_metni(time.time() - t0),
                       sn_metni(time.time() - son['t'])))

        th = threading.Thread(target=canlilik)
        th.daemon = True
        th.start()
        try:
            for ham in p.stdout:
                s = ham.decode('utf-8', 'replace').rstrip('\r\n')
                satirlar.append(s)
                son['t'] = time.time()
                yaz(u'      | ' + s)
                if KESME['var']:
                    break
            p.wait()
        finally:
            dur.set()
        rc = p.returncode if p.returncode is not None else -1
        son_cikti = u'\n'.join(satirlar)
        if KESME['var']:
            raise Kesildi()
        if rc != 0:
            yaz(u'   >>> komut %d/%d SIFIR DISI KOD DONDURDU: %s'
                % (i, len(komutlar), rc))
            break
    return rc, son_cikti, time.time() - t0


# ===========================================================================
#  6) PLAN, OZET, ANA AKIS
# ===========================================================================
def kraken_bul(kok):
    aday = os.path.abspath(os.path.join(kok, '..', 'ALI', 'WSL', '130_KRAKEN_ARAC.sh'))
    if os.path.exists(aday):
        return aday, u''
    return None, (u'ALI klasoru bulunamadi (aranan: %s). W, X ve Z ATLANIR; '
                  u'zincirin geri kalani kosar.' % aday)


def plan_yaz(kok, asamalar, durum, ayar, yaz):
    yaz(u'')
    yaz(u'=' * 78)
    yaz(u'  PLAN - hangi asama kosacak, hangisi atlanacak')
    yaz(u'=' * 78)
    toplam, olculmeyen = 0.0, []
    for a in asamalar:
        atla, sebep = kontrol_noktasi_gecerli(kok, a, durum)
        if a.get('kraken') and not ayar.get('karac'):
            atla, sebep = True, ayar.get('kraken_sebep', u'Kraken araci yok')
            a['_kraken_atla'] = True
        a['_atla'], a['_sebep'] = atla, sebep
        if atla:
            sure = u'0 sn'
        elif a['sure_sn'] is None:
            sure = u'OLCULMEDI'
            olculmeyen.append(a['kod'])
        else:
            sure = sn_metni(a['sure_sn'])
            toplam += a['sure_sn']
        yaz(u'  %-2s %-8s %-12s %s' % (a['kod'], u'ATLANIR' if atla else u'KOSACAK',
                                       sure, a['ad'][:44]))
        yaz(u'       %s' % sebep)
    yaz(u'')
    yaz(u'  OLCULMUS asamalarin toplami : %s' % sn_metni(toplam))
    if olculmeyen:
        yaz(u'  SURESI OLCULMEMIS asamalar : %s' % u', '.join(olculmeyen))
        yaz(u'  Bu asamalara sayi UYDURULMADI. Yukaridaki toplam bu yuzden bir')
        yaz(u'  ALT SINIRDIR - tahmin degil, eksik olcumdur.')
    yaz(u'=' * 78)
    return toplam, olculmeyen


# Nihai hukum tablosunun kaynagi. Sirayla denenir; ilk BULUNAN kullanilir.
# En basta en YENI dosya durur. Yeni bir hukum turu uretilirse bu listenin
# basina eklenmeli - kod dosyayi kendisi secmez, burada yazili sirayi izler.
# 2026-08-10 DUZELTME. Liste basinda 08-08 tarihli ESIK_VE_OLCUT duruyordu ve
# sabah ozeti "SIPARIS EDILEBILIR 11 cift" yaziyordu; oysa AYNI KOSUNUN kendi
# urettigi TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv 20 cift diyor (15 KESIN +
# 5 EVRENSEL). Aradaki fark, 08-08'den sonra bes ciftin degistirilmis ve
# Microascaceae'nin geri alinmis olmasi. Ozet dosyadan okudugunu acikca
# yaziyordu, yani yalan soylemiyordu - ama YANLIS DOSYAYI okuyordu ve sabah
# ilk bakilan yer orasi. Kosunun KENDI urettigi tablo listenin basina alindi.
SIPARIS_KAYNAKLARI = [
    ('TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv', 'durum'),
    ('ESIK_VE_OLCUT_2026-08-08.tsv', 'YENI_HUKUM'),
    ('NIHAI_SIPARIS_LISTESI_2026-08-07.tsv', None),
]

# SIPARIS EDILEBILIR sayilan hukum onekleri. Kural ACIKCA yazili, kodun icine
# gomulu bir sezgi degil:
#   * "SIPARIS EDILEBILIR..." -> grup ozgul, esigi gecen yedi cift
#   * "KOSULLU..."            -> evrensel/kontrol primerleri. Onlarda dCq
#                                tanimsizdir (rakip kumesi bosa yaklasir), bu
#                                yuzden esik hukmu verilemez; olcu KAPSAMDIR.
# "ESIK ALTI...", "ONERILMEZ...", "UYELIK DOGRULANAMADI" siparise GIRMEZ.
#   * "ESIK USTU - SIPARIS EDILEBILIR" -> SIPARIS_LISTESI.tsv'nin kendi
#                                dili. Ayni kural, farkli sozcuk sirasi.
SIPARIS_ONEKLERI = (u'SIPARIS EDILEBILIR', u'KOSULLU', u'ESIK USTU')


def siparis_tablosu(kok):
    u"""Nihai siparis tablosunu DISKTEKI dosyadan okur; YENIDEN HESAPLAMAZ."""
    y = hk = None
    for ad, sut_adi in SIPARIS_KAYNAKLARI:
        t = os.path.join(kok, *ad.split('/'))
        if os.path.exists(t):
            y, hk = t, sut_adi
            break
    if not y:
        return (u'Nihai hüküm tablosu bulunamadı. Aranan dosyalar: %s\n'
                % u', '.join(u'`%s`' % a for a, _ in SIPARIS_KAYNAKLARI))
    with io.open(y, encoding='utf-8', errors='ignore') as fh:
        r = list(csv.DictReader((s for s in fh if not s.startswith('#')), delimiter='\t'))
    ad = os.path.basename(y)
    if not r:
        return u'`%s` boş.\n' % ad
    sut = [k for k in r[0].keys() if k]
    if not hk or hk not in sut:
        hk = (next((k for k in sut if 'HUKUM' in k.upper()), None)
              or next((k for k in sut if 'KARAR' in k.upper()), None)
              or next((k for k in sut if 'SINIF' in k.upper()), None))
    hd = next((k for k in sut if 'hedef' in k.lower()), None)
    if not hk or not hd:
        return (u'`%s` okundu (%d satır) ama hüküm/hedef sütunu tanınamadı. '
                u'Sütunlar: %s\n' % (ad, len(r), u', '.join(sut)))

    def hukum(s):
        return (s.get(hk) or u'?').strip()

    # SINIF sutunu varsa hukum ONDAN okunur. Sebep: evrensel primerlerde dCq
    # TANIMSIZDIR (rakip kumesinin paydasi sifira gider) ve tablo onlara
    # "OLCULEMEDI - KARAR YOK" yazar. Yalniz metne bakan bir kural o bes
    # kontrol primerini siparis disi sayar - oysa panelin kontrolleri onlar ve
    # siparise GIRERLER. dogrulama_turu.siparistekiler() de ayni kurali
    # kullaniyor; iki yer ayni tanimi kullanmazsa iki farkli sayi uretilir.
    sinif_s = next((k for k in sut if k.strip().upper() == 'SINIF'), None)

    def siparise_girer(s):
        if sinif_s:
            return (s.get(sinif_s) or '').strip().upper() in ('KESIN', 'EVRENSEL')
        return hukum(s).upper().startswith(SIPARIS_ONEKLERI)

    girer = [s for s in r if siparise_girer(s)]
    girmez = [s for s in r if not siparise_girer(s)]
    say = {}
    for s in r:
        say[hukum(s)] = say.get(hukum(s), 0) + 1

    out = [u'Kaynak: `%s` (%d satır), sütun `%s`. '
           u'**Bu tabloyu bu koşu yeniden ÜRETMEDİ**, dosyadan okudu.\n\n' % (ad, len(r), hk),
           u'**SİPARİŞ EDİLEBİLİR: %d çift = %d oligo** · kalan %d çift '
           u'sipariş dışı (silinmedi, aşağıda).\n\n' % (len(girer), 2 * len(girer),
                                                        len(girmez)),
           (u'Kural: `SINIF` sütunu `KESIN` ya da `EVRENSEL` olan satır siparişe '
            u'girer. Evrensel/kontrol primerlerinde dCq TANIMSIZDIR (rakip '
            u'kümesinin paydası sıfıra gider) ve tablo onlara "OLCULEMEDI - KARAR '
            u'YOK" yazar; ölçü kapsamdır, eşik değil.\n\n') if sinif_s else
           (u'Kural: hükmü `SIPARIS EDILEBILIR...` ya da `KOSULLU...` ile başlayan '
            u'satır siparişe girer.\n\n'),
           u'| %s | kaç çift |\n|---|---|\n' % hk]
    for k in sorted(say, key=lambda x: (-say[x], x)):
        out.append(u'| %s | %d |\n' % (k, say[k]))

    dcq = next((k for k in sut if 'dCq' in k), None)
    out.append(u'\n### Sipariş edilecek %d çift\n\n| hedef | hüküm | dCq |\n|---|---|---|\n'
               % len(girer))
    for s in girer:
        out.append(u'| %s | %s | %s |\n' % (s.get(hd) or u'?', hukum(s),
                                            (s.get(dcq) or u'-') if dcq else u'-'))
    out.append(u'\n### Sipariş dışı %d çift (silinmedi)\n\n'
               u'| hedef | hüküm | dCq |\n|---|---|---|\n' % len(girmez))
    for s in girmez:
        out.append(u'| %s | %s | %s |\n' % (s.get(hd) or u'?', hukum(s),
                                            (s.get(dcq) or u'-') if dcq else u'-'))
    return u''.join(out)


def ozet_yaz(kok, asamalar, durum, ayar, kesildi, on_uyari, baslangic, gunluk_yolu):
    yol = os.path.join(kok, CIKTI_KLASOR, '00_SABAH_OZETI.md')
    basari = [a for a in asamalar if durum.get(a['kod'], {}).get('durum') == 'bitti']
    dusen = [a for a in asamalar if durum.get(a['kod'], {}).get('durum') == 'DUSTU']
    atlanan = [a for a in asamalar
               if (durum.get(a['kod'], {}).get('durum') or '').startswith('atlandi')]
    with io.open(yol, 'w', encoding='utf-8') as fh:
        w = fh.write
        w(u'# Sabah özeti — tek tuş koşusu\n\n')
        w(u'Üretim: %s · sürüm %s\n\n' % (time.strftime('%Y-%m-%d %H:%M'), SURUM))
        w(u'Koşu başladı: %s · geçen süre: %s\n\n'
          % (time.strftime('%Y-%m-%d %H:%M', time.localtime(baslangic)),
             sn_metni(time.time() - baslangic)))
        if kesildi:
            w(u'> **KOŞU KESİLDİ** (Ctrl+C ya da pencere kapatıldı). Biten aşamalar '
              u'kayıtlı; aynı tuşa yeniden basmak kaldığı yerden devam ettirir.\n\n')

        w(u'## 1. Hangi aşama ne oldu\n\n')
        w(u'| aşama | durum | süre | çıkış kodu | not |\n|---|---|---|---|---|\n')
        for a in asamalar:
            d = durum.get(a['kod'], {})
            w(u'| **%s** %s | %s | %s | %s | %s |\n'
              % (a['kod'], a['ad'][:40],
                 (u'bitti (UYARILI)' if d.get('uyarili') else d.get('durum', 'koşulmadı')),
                 sn_metni(d.get('sure', 0)) if d.get('sure') else '-',
                 d.get('cikis', '-'), (d.get('sebep') or '').replace('|', '/')[:170]))
        w(u'\n**Başarılı %d · Düşen %d · Atlanan %d**\n\n'
          % (len(basari), len(dusen), len(atlanan)))

        uyarili = [a for a in asamalar if durum.get(a['kod'], {}).get('uyarili')]
        if uyarili:
            # Danisma asamasi zinciri durdurmaz ama SESSIZ de kalmaz. Bulgu
            # varsa ozetin basinda gorunur; yoksa denetim yapilmis ama kimse
            # bakmamis olur ve kapinin bir anlami kalmaz.
            w(u'## UYARI — denetim bulgu buldu, zincir yine de koştu\n\n')
            for a in uyarili:
                d = durum[a['kod']]
                w(u'- **%s** %s — %s\n' % (a['kod'], a['ad'], (d.get('sebep') or '')[:200]))
            w(u'\nAyrıntı: `TEK_TUS_SONUC/DENETIM_RAPORU.md` ve '
              u'`GECE_BULGULARI.md`.\n\n')

        if dusen:
            w(u'## 2. DÜŞEN AŞAMALAR — önce bunlara bakın\n\n')
            for a in dusen:
                d = durum[a['kod']]
                w(u'### %s — %s\n\n' % (a['kod'], a['ad']))
                w(u'- çıkış kodu: `%s`\n- sebep: %s\n' % (d.get('cikis'), d.get('sebep')))
                if d.get('son_satirlar'):
                    w(u'- aşamanın son çıktısı:\n\n```\n%s\n```\n' % d['son_satirlar'])
                w(u'\n')
        else:
            w(u'## 2. DÜŞEN AŞAMA YOK\n\nHiçbir aşama sıfır dışı kod döndürmedi ve '
              u'hiçbir çıktı denetimi düşmedi.\n\n')

        w(u'## 3. Nihai sipariş tablosu\n\n')
        w(siparis_tablosu(kok))

        w(u'\n## 4. Burak\'ın bakması gereken dosyalar\n\n')
        w(u'| soru | dosya |\n|---|---|\n')
        for soru, dosya in (
            (u'Ne sipariş edeyim (NİHAİ)',
             u'`TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv` — `durum` ve '
             u'`siparis_sarti` sütunları (bu koşunun KENDİ ürettiği tablo)'),
            (u'Toplantıda ne istendi, hangisi oldu', u'`TOPLANTI_DURUMU.md`'),
            (u'NCBI 4. katman ne dedi', u'`DOGRULAMA_SONUC/NCBI_KATMAN4_RAPORU.md`'),
            (u'Bu koşuda denetimler temiz mi', u'`TEK_TUS_SONUC/DENETIM_RAPORU.md`'),
            (u'Eşik kuralı neden değişti', u'`ESIK_VE_OLCUT_2026-08-08.md`'),
            (u'Diziler nereden kopyalanacak',
             u'`MicRhoBooster_PANEL_*.xlsx` (`1 Siparis` sayfası)'),
            (u'Hangi çift riskli, neden', u'`SIPARIS_KARARI_2026-08-07.md`'),
            (u'Bu koşuda çelişki çıktı mı', u'`DOGRULAMA_SONUC/CELISKILER.md`'),
            (u'Doğrulama katmanları yan yana', u'`DOGRULAMA_SONUC/DOGRULAMA_RAPORU.md`'),
            (u'Hedef dışı ürünlerin ayrıntısı', u'`HEDEF_DISI_AYRINTI_2026-08-07.tsv`'),
            (u'Kutu kimlikleri doğru mu', u'`TUM_KIMLIK_SONUC/TUM_KUTU_KIMLIK_RAPORU.md`'),
            (u'Bu koşunun ham çıktısı',
             u'`%s`' % os.path.relpath(gunluk_yolu, kok).replace('\\', '/')),
        ):
            w(u'| %s | %s |\n' % (soru, dosya))

        if on_uyari:
            w(u'\n## 5. Ön kontrol uyarıları\n\n')
            for u_ in on_uyari:
                w(u'- %s\n' % u_)

        w(u'\n---\n\n### Bu koşunun ölçemedikleri\n\n')
        olcumsuz = [a for a in asamalar if a['sure_sn'] is None]
        for a in olcumsuz:
            w(u'- **%s** süresi ölçülmedi: %s\n' % (a['kod'], a['kaynak']))
        if not olcumsuz:
            w(u'- Koşan bütün aşamaların süresi daha önce ölçülmüştü.\n')
        if not ayar.get('karac'):
            w(u'- Kraken karşılaştırması (W, X, Z) KOŞULMADI: %s\n'
              % ayar.get('kraken_sebep', ''))
    return yol


def main():
    global CANLILIK_SN
    p = argparse.ArgumentParser(description=u'TEK TUS - butun zinciri sirayla kosar')
    p.add_argument('--kok', default='.')
    p.add_argument('--plan', action='store_true', help=u'yalniz plani yaz, kosma')
    p.add_argument('--onayla', action='store_true',
                   help=u'plani gosterip onay BEKLEMEDEN kos (bat bunu verir)')
    p.add_argument('--yeniden', action='store_true',
                   help=u'bitmis asamalari da yeniden kos')
    p.add_argument('--yalniz', default='', help=u'yalniz bu asamalar, orn: 8HS')
    p.add_argument('--atla', default='', help=u'bu asamalari atla, orn: IG')
    p.add_argument('--ncbi', choices=['oto', 'elle', 'yok'], default='oto')
    p.add_argument('--organizma',
                   default='Bacteria (taxid:2) OR Archaea (taxid:2157) OR Fungi (taxid:4751)')
    p.add_argument('--siparis-16', action='store_true',
                   help=u'D asamasi 22 cift yerine yalniz 16 siparis ciftini sinar')
    p.add_argument('--vt', default=os.environ.get('VT_A', ''),
                   help=u'Kraken2 veritabani yolu')
    p.add_argument('--on-kontrol-atla', action='store_true',
                   help=u'ON KONTROLU ATLA - tavsiye edilmez, ekrana yazilir')
    p.add_argument('--canlilik', type=int, default=CANLILIK_SN,
                   help=u'kac saniyede bir canlilik isareti')
    A = p.parse_args()
    CANLILIK_SN = max(2, A.canlilik)

    kok = os.path.abspath(A.kok)
    os.makedirs(os.path.join(kok, CIKTI_KLASOR), exist_ok=True)
    zaman = time.strftime('%Y%m%d_%H%M%S')
    gunluk_yolu = os.path.join(kok, CIKTI_KLASOR, 'gunluk_%s.log' % zaman)
    gunluk = io.open(gunluk_yolu, 'w', encoding='utf-8')
    tampon = {'t': 0.0}

    def yaz(s=u''):
        # Ekrana zamansiz, gunluge ZAMAN DAMGALI. Bagli klasore yazma sikligi
        # LOG_YAZMA_ARALIGI ile sinirli (D-11: /mnt/c'ye cok sayida kucuk
        # yazmanin yavas oldugu olculmustu).
        try:
            print(s, flush=True)
        except UnicodeEncodeError:
            print(s.encode('ascii', 'replace').decode('ascii'), flush=True)
        gunluk.write(u'%s  %s\n' % (time.strftime('%H:%M:%S'), s))
        if time.time() - tampon['t'] > LOG_YAZMA_ARALIGI:
            gunluk.flush()
            tampon['t'] = time.time()

    def sig(signum, frame):
        KESME['var'] = True
        yaz(u'\n  !! KESME ISTEGI ALINDI. Durum yaziliyor, temiz cikiliyor...')
    signal.signal(signal.SIGINT, sig)
    try:
        signal.signal(signal.SIGTERM, sig)
    except Exception:
        pass

    baslangic = time.time()
    yaz(u'=' * 78)
    yaz(u'  TEK TUS - MicRhoBooster tam zinciri   surum %s' % SURUM)
    yaz(u'  %s' % time.strftime('%Y-%m-%d %H:%M:%S'))
    yaz(u'  kok    : %s' % kok)
    yaz(u'  gunluk : %s' % os.path.relpath(gunluk_yolu, kok).replace('\\', '/'))
    yaz(u'=' * 78)

    karac, kraken_sebep = kraken_bul(kok)
    ayar = dict(ncbi=A.ncbi, organizma=A.organizma, tumu=not A.siparis_16,
                karac=karac, kraken_sebep=kraken_sebep, vt_a=A.vt)

    on_uyari = []
    if A.on_kontrol_atla:
        yaz(u'\n  !! ON KONTROL ATLANDI (--on-kontrol-atla). Eksik bir sey varsa '
            u'kosunun ORTASINDA cikacak.')
    else:
        ok, _dusen, on_uyari = on_kontrol(kok, ayar, yaz)
        if not ok:
            gunluk.flush()
            gunluk.close()
            return 2

    hepsi = ASAMALAR(ayar)
    sec = [a for a in hepsi
           if (not A.yalniz or a['kod'] in A.yalniz.upper())
           and a['kod'] not in A.atla.upper()]
    if not sec:
        yaz(u'  Secilen asama yok (--yalniz / --atla). Hicbir sey kosulmadi.')
        gunluk.close()
        return 0

    dyol = os.path.join(kok, CIKTI_KLASOR, 'durum.json')
    durum = {}
    if os.path.exists(dyol) and not A.yeniden:
        try:
            durum = json.load(io.open(dyol, encoding='utf-8'))
        except Exception as e:
            yaz(u'  UYARI: durum.json okunamadi (%s) - kontrol noktalari yok sayildi.' % e)
            durum = {}
    if A.yeniden:
        yaz(u'  --yeniden verildi: butun kontrol noktalari yok sayiliyor.')
        durum = {}

    plan_yaz(kok, sec, durum, ayar, yaz)
    if A.plan:
        gunluk.flush()
        gunluk.close()
        return 0
    if not A.onayla:
        try:
            c = input(u'\n  Baslatilsin mi? (E = evet, baska tus = vazgec): ').strip()
        except EOFError:
            c = u''
        if c.upper() not in (u'E', u'EVET'):
            yaz(u'  Vazgecildi, hicbir sey kosulmadi.')
            gunluk.close()
            return 0

    def kaydet():
        json.dump(durum, io.open(dyol, 'w', encoding='utf-8'),
                  ensure_ascii=False, indent=1)

    kesildi = False
    kalan = [a for a in sec if not a['_atla']]
    yaz(u'\n  Kosacak asama sayisi: %d' % len(kalan))
    sirano = 0
    for a in sec:
        kod = a['kod']
        # PLAN BIR TAHMINDIR, BURASI HUKUMDUR. Onceki asamalar diski
        # degistirmis olabilir (ornek: U kosunca P'nin girdisi degisir), bu
        # yuzden atlama karari asamaya GELINDIGI ANDA yeniden hesaplanir.
        if not a.get('_kraken_atla'):
            yeni_atla, yeni_sebep = kontrol_noktasi_gecerli(kok, a, durum)
            if yeni_atla != a['_atla']:
                yaz(u'\n   NOT: %s icin karar PLANDAN FARKLI cikti (%s -> %s). '
                    u'Sebep: %s' % (kod, u'ATLANIR' if a['_atla'] else u'KOSACAK',
                                    u'ATLANIR' if yeni_atla else u'KOSACAK', yeni_sebep))
            a['_atla'], a['_sebep'] = yeni_atla, yeni_sebep

        if a['_atla']:
            if a.get('_kraken_atla'):
                durum[kod] = dict(durum='atlandi (arac yok)', sebep=a['_sebep'], sure=0)
            elif durum.get(kod, {}).get('durum') != 'bitti':
                durum[kod] = dict(durum='bitti', sure=0, cikis=0,
                                  imza=girdi_imzasi(kok, a), sebep=a['_sebep'],
                                  zaman=time.strftime('%Y-%m-%d %H:%M'))
            yaz(u'\n>> %s  %s\n   ATLANDI - %s' % (kod, a['ad'], a['_sebep']))
            kaydet()
            continue

        # Bagimlilik kapisi: bagimli oldugu bir asama dustuyse KOSMA.
        engel = [b for b in a['bagimli']
                 if (durum.get(b, {}).get('durum') or '').startswith(
                     ('DUSTU', 'atlandi (bagimli', 'kesildi'))]
        if engel:
            durum[kod] = dict(durum='atlandi (bagimli)', sure=0,
                              sebep=u'%s asamasi bitmedigi icin kosulmadi. Bu asama '
                                    u'onun ciktisini GIRDI olarak kullaniyor; bos '
                                    u'girdiyle kosmak inandirici ama anlamsiz bir '
                                    u'sonuc uretirdi.' % u', '.join(engel))
            yaz(u'\n>> %s  %s\n   ATLANDI (bagimli) - %s'
                % (kod, a['ad'], durum[kod]['sebep']))
            kaydet()
            continue

        sirano += 1
        kalan_sn = sum(x['sure_sn'] or 0 for x in kalan[sirano - 1:])
        yaz(u'\n' + u'-' * 78)
        yaz(u'>> [%d/%d] %s  %s' % (sirano, len(kalan), kod, a['ad']))
        yaz(u'   basladi : %s' % time.strftime('%H:%M:%S'))
        yaz(u'   olculmus sure: %s   (kaynak: %s)'
            % (sn_metni(a['sure_sn']) if a['sure_sn'] is not None else u'OLCULMEDI',
               a['kaynak']))
        yaz(u'   kalan asama: %d | kalan olculmus sure: %s'
            % (len(kalan) - sirano, sn_metni(kalan_sn)))
        yaz(u'-' * 78)

        try:
            rc, cikti_metni, sure = asama_kos(kok, a, ayar, yaz)
        except Kesildi:
            kesildi = True
            durum[kod] = dict(durum='kesildi', sure=0,
                              sebep=u'kullanici kesti. Bu asama YARIM; kendi kontrol '
                                    u'noktalari kayitli, ayni tusa basmak kaldigi '
                                    u'yerden devam ettirir.')
            kaydet()
            yaz(u'   KESILDI - durum kaydedildi.')
            break

        son_satirlar = u'\n'.join((cikti_metni or u'').splitlines()[-15:])
        tamam, mesaj = a['denet'](kok, ayar, cikti_metni)

        # IKI AYRI SUZGEC, IKISI DE GECMELI.
        # Gecmiste yalniz cikti denetimine bakilip sifir disi kod goz ardi
        # edilmisti (tam_zincir.py, T asamasi, 2026-08-06): T cikis kodu 3
        # dondurdugu halde "BITTI" damgasi almisti ve ozet yaniltici cikmisti.
        if rc != 0:
            mesaj = (u'CIKIS KODU %s (sifir degil). Cikti denetimi: %s' % (rc, mesaj))
            tamam = False
        # DANISMA ASAMASI: bulgu bildirir ama zinciri DUSURMEZ.
        # N (DENETIM) boyle bir asamadir: isi "su an neyin bozuk oldugunu
        # soylemek". Bozuk bir sey bulmasi zincirin kosmamasi demek degil;
        # tersine, zincir kossun ki ne uretildigi de gorulsun. Ama sessiz de
        # kalmaz - ozete UYARILI yazilir ve bulgular DENETIM_RAPORU.md'de
        # durur. (2026-08-10: N eklenince tek_tus_sinama S1 senaryosu
        # dusuyordu; kapinin isi sinav dusurmek degil.)
        if a.get('danisma') and rc != 0:
            mesaj = u'UYARILI - %s (danisma asamasi, zincir durdurulmadi)' % mesaj
            tamam = True
            uyarili = True
        else:
            uyarili = False
        if not tamam:
            durum[kod] = dict(durum='DUSTU', sure=sure, cikis=rc, sebep=mesaj,
                              son_satirlar=son_satirlar,
                              zaman=time.strftime('%Y-%m-%d %H:%M'))
            yaz(u'   << %s DUSTU (%s) - %s' % (kod, sn_metni(sure), mesaj))
            bagimlilar = [x['kod'] for x in sec if kod in x['bagimli']]
            if bagimlilar:
                yaz(u'   Buna bagimli asamalar KOSULMAYACAK: %s' % u', '.join(bagimlilar))
            else:
                yaz(u'   Buna bagimli baska asama yok; zincir devam ediyor.')
        else:
            durum[kod] = dict(durum='bitti', sure=sure, cikis=rc, sebep=mesaj,
                              uyarili=uyarili, imza=girdi_imzasi(kok, a),
                              zaman=time.strftime('%Y-%m-%d %H:%M'))
            yaz(u'   << %s BITTI (%s) - %s' % (kod, sn_metni(sure), mesaj))
        kaydet()

    ozet = ozet_yaz(kok, sec, durum, ayar, kesildi, on_uyari, baslangic, gunluk_yolu)
    dusen = [a['kod'] for a in sec if durum.get(a['kod'], {}).get('durum') == 'DUSTU']
    atlanan_b = [a['kod'] for a in sec
                 if durum.get(a['kod'], {}).get('durum') == 'atlandi (bagimli)']
    yaz(u'\n' + u'=' * 78)
    if kesildi:
        yaz(u'  KOSU KESILDI. Biten asamalar kayitli; ayni tusa yeniden basin.')
    elif dusen:
        yaz(u'  KOSU BITTI ama %d ASAMA DUSTU: %s' % (len(dusen), u', '.join(dusen)))
        if atlanan_b:
            yaz(u'  Bagimli oldugu icin kosulmayanlar: %s' % u', '.join(atlanan_b))
    else:
        yaz(u'  BUTUN ASAMALAR BASARILI.')
    yaz(u'  Toplam gecen sure: %s' % sn_metni(time.time() - baslangic))
    yaz(u'')
    yaz(u'  SABAH BUNA BAKIN : %s' % os.path.relpath(ozet, kok).replace('\\', '/'))
    yaz(u'  Ham cikti        : %s' % os.path.relpath(gunluk_yolu, kok).replace('\\', '/'))
    yaz(u'=' * 78)
    gunluk.flush()
    gunluk.close()
    if kesildi:
        return 130
    return 3 if dusen else 0


if __name__ == '__main__':
    sys.exit(main() or 0)
