# -*- coding: utf-8 -*-
"""
tam_zincir.py - bastan sona TEK SECENEKLE kosan tam zincir surucusu.

Sira: 8 -> H -> W -> I -> G -> T -> X -> Y -> Z -> S

Bu betik OLCUM YAPMAZ. Yaptigi is, sirasi ve bagimliligi belli on asamayi
dogru sirayla cagirmak, her asamanin ciktisini DENETLEMEK ve bir sey tutmazsa
DURMAKTIR. Olcumun kendisi cagrilan betiklerin icindedir.
"""
# ---------------------------------------------------------------------------
# tam_zincir.py - KAPSAMLI_ARAMA.bat menusundeki A tusunun motoru.
#
# GIRDI  : proje kokundeki asama betikleri (KAPSAMLI_ARAMA paketi, KURTARMA
#          altindaki dort betik) ve ../ALI/WSL/130_KRAKEN_ARAC.sh;
#          ayrica TAM_ZINCIR_SONUC/durum.json (onceki kosunun nerede kaldigi).
# CIKTI  : TAM_ZINCIR_SONUC/00_TAM_ZINCIR_OZET.md (tek birlesik ozet),
#          TAM_ZINCIR_SONUC/durum.json, TAM_ZINCIR_SONUC/kosu_gunlugu.txt.
# CAGRAN : KAPSAMLI_ARAMA.bat -> A tusu
#          (bat icinde: wsl -e python3 "KURTARMA/tam_zincir.py" --kok . --onayla)
#
# UC TASARIM KARARI, gerekceleriyle
#
# 1) KESINTIYE DAYANIKLILIK ASAMA DUZEYINDEDIR.
#    Her asama bitince durum.json'a "bitti" damgasi ve suresi yazilir. Ayni
#    secenege yeniden basildiginda bitmis asamalar ATLANIR. Asama ICINDEKI
#    kesinti dayanikliligi bu betigin isi degildir; cagrilan betiklerin kendi
#    kontrol noktalari o isi zaten yapiyor. Burada iki katman ust uste binmez.
#
# 2) CIKTI DENETIMI ZORUNLUDUR, CIKIS KODU YETMEZ.
#    Bu projede tekrar tekrar cikan hata turu sudur: program hata vermeden
#    YANLIS ya da BOS cevap uretir. Sifir cikis kodu "is bitti" demek degildir.
#    Bu yuzden her asamanin ardindan beklenen dosya VAR MI, BOS MU ve - bazi
#    asamalarda - ICERIGI ne diyor diye bakilir. Denetim dusen asama zinciri
#    DURDURUR; sessizce bir sonrakine gecilmez.
#
# 3) KRAKEN ADIMLARI ZINCIRI KIRMAZ.
#    W, X, Y ve Z adimlari ayri bir arac (kraken2) ve ayri bir veritabani ister.
#    Bunlar kurulu degilse sekiz saatlik bir kosuyu bastan reddetmek yanlis
#    olurdu; o adimlar ATLANDI diye isaretlenir, sebebi yazilir ve zincir
#    devam eder. Ozette hangi adimlarin neden atlandigi ayrica durur.
# ---------------------------------------------------------------------------

from __future__ import print_function
import argparse
import csv
import io
import json
import os
import shutil
import subprocess
import sys
import time

VERSIYON = u'1.0 (2026-08-04)'
CIKTI_ADI = 'TAM_ZINCIR_SONUC'


# ---------------------------------------------------------------------------
# ASAMA TABLOSU
#
# Alanlar:
#   kod      menudeki tek harf/rakam - ozet ve durum.json anahtari
#   ad       ekranda ve ozette gorunen ad
#   grup     menudeki grup (raporda gruplama icin)
#   dk       tahmini sure, DAKIKA cinsinden (alt, ust). Sadece tahmindir.
#   kraken   True ise kraken2 yoksa ATLANIR, zincir kirilmaz
#   komut    calistirilacak komut uretici (kok, ayar) -> [argv, ...]
#   denet    cikti denetleyicisi (kok, ayar, ciktilar) -> (tamam, mesaj)
#
# Sira KEYFI DEGILDIR:
#   8 once, cunku kod kendini sinamadan hicbir olcume girilmez.
#   H sonra, cunku sekiz saatlik zincire girmeden once zincirin butun oldugu
#     kucuk ve cevabi bilinen bir alt kumede gosterilmelidir.
#   W kraken adimlarindan once, cunku hangi veritabaninin kullanildigini
#     bilmeden esik taramasi baslatmak yanlis veritabaninda saatler harcar.
#   I ve G, T'den once: T'nin uyelik kararlari kutu kimliklerine dayanir.
#   T ortada: P, K, D asamalarini kendi icinde bagimlilik sirasiyla kosar.
#   X ve Y olcum, Z onlarin ciktisini okur - bu yuzden Z en sonda.
#   S en sonda, cunku ozet ancak her sey bittikten sonra anlamlidir.
# ---------------------------------------------------------------------------


def _py(kok, *arg):
    return [sys.executable] + list(arg)


def _satir_sayisi(yol):
    """TSV'de BASLIK HARIC veri satiri sayar. Yorum satirlari sayilmaz."""
    if not os.path.exists(yol):
        return -1
    n = 0
    with io.open(yol, encoding='utf-8', errors='ignore') as fh:
        for i, s in enumerate(fh):
            if s.startswith('#') or not s.strip():
                continue
            n += 1
    return max(n - 1, 0)


def d_dosya_dolu(yollar, en_az_satir=1):
    """Beklenen dosyalarin hepsi VAR mi ve BOS DEGIL mi.

    Neden ayri bir denetim: bu projede bir asamanin sifir kodla bitip hicbir
    satir uretmemesi fiilen gorulmus bir hatadir. 'Dosya var' yetmez, 'icinde
    veri var' da sorulur.
    """
    def f(kok, ayar):
        eksik, bos = [], []
        for y in yollar:
            t = os.path.join(kok, y)
            if not os.path.exists(t):
                eksik.append(y)
                continue
            if y.endswith('.tsv'):
                if _satir_sayisi(t) < en_az_satir:
                    bos.append(y)
            elif os.path.getsize(t) < 40:
                bos.append(y)
        if eksik:
            return False, u'beklenen cikti YOK: %s' % ', '.join(eksik)
        if bos:
            return False, u'cikti BOS ya da veri satiri yok: %s' % ', '.join(bos)
        return True, u'cikti dogrulandi: %s' % ', '.join(yollar)
    return f


def d_sina(kok, ayar):
    """8 tusu: selftest ciktisi 'TUM SINAMALAR GECTI' demeli."""
    g = ayar.get('_son_cikti', '')
    if u'TUM SINAMALAR GECTI' in g:
        return True, u'butun selftestler gecti'
    return False, (u'selftest metni bulunamadi. Kod kendini dogrulayamadi; '
                   u'olcume girilmez.')


def d_hizli_test(kok, ayar):
    """H tusu: rapor TUTARLI mi diyor. TUTARSIZ ise zincir DURUR.

    Gerekce: bu adimin varlik sebebi tam olarak budur. Tutarsiz bir zincirde
    sekiz saat harcamak, sonunda cope gidecek bir sayi uretir.
    """
    y = os.path.join(kok, 'HIZLI_TEST', 'HIZLI_TEST_RAPORU.md')
    if not os.path.exists(y):
        return False, u'HIZLI_TEST/HIZLI_TEST_RAPORU.md uretilmedi'
    m = io.open(y, encoding='utf-8', errors='ignore').read()
    if u'ZINCIR TUTARSIZ' in m:
        return False, (u'rapor ZINCIR TUTARSIZ diyor. Tam kosuya GIRILMEZ; '
                       u'HIZLI_TEST/HIZLI_TEST_RAPORU.md dosyasina bakin.')
    if u'ZINCIR TUTARLI' in m:
        return True, u'rapor ZINCIR TUTARLI diyor'
    return False, u'raporda ne TUTARLI ne TUTARSIZ karari var - bicim beklenmedik'


def d_esik(kok, ayar):
    """X ve Y: esik egrisi uretildi mi ve AYRILIK sutunu bos mu.

    AYRILIK dolu demek, iki bagimsiz olcumun ayni esikte farkli sayi vermesi
    demektir. Proje kurali geregi o durumda sayilara guvenilmez ve zincir durur.
    """
    a = ayar['kraken_is']
    csvy = os.path.join(a, 'esik_egrisi.csv')
    if not os.path.exists(csvy):
        return False, u'esik_egrisi.csv uretilmedi (%s)' % csvy
    try:
        with io.open(csvy, encoding='utf-8', errors='ignore') as fh:
            r = list(csv.DictReader(fh))
    except Exception as e:
        return False, u'esik_egrisi.csv okunamadi: %s' % e
    if not r:
        return False, u'esik_egrisi.csv BOS - hicbir esik olculmemis'
    ad = [k for k in (r[0].keys() if r else []) if k and k.strip().lower() == 'ayrilik']
    if ad:
        dolu = [x for x in r if (x.get(ad[0]) or '').strip()]
        if dolu:
            return False, (u'AYRILIK sutunu DOLU (%d satir). Iki bagimsiz olcum '
                           u'ayrildi; sayilara guvenilmez.' % len(dolu))
    return True, u'%d esik olculdu, ayrilik yok' % len(r)


def d_tablo(kok, ayar):
    """Z: tablo uretildi mi ve TAM mi.

    UC DURUMLU DENETIM. Bu asama ozel bir durum tasir: tablo EKSIK sutunlarla da
    uretilebilir ve bu bir hata degildir - Kraken olcumu sonradan yapilinca ayni
    tus tabloyu tamamlar. Ama eksik tablo 'bitti' damgasi ALMAMALIDIR; alsaydi
    devam mantigi onu bir daha hic denemez ve tablo sonsuza kadar eksik kalirdi.
    Bu yuzden ucuncu bir durum var:
        True    -> tam tablo, bitti
        'eksik' -> tablo uretildi ama sutun(lar) bos; zincir DURMAZ, bir sonraki
                   kosuda bu asama YENIDEN denenir
        False   -> tablo hic uretilemedi; zincir durur
    """
    y = os.path.join(ayar['ali'], '0_TESLIM_HOCA', 'KRAKEN_KARSILASTIRMA.md')
    if not os.path.exists(y):
        return False, u'KRAKEN_KARSILASTIRMA.md uretilmedi'
    if os.path.getsize(y) < 200:
        return False, u'KRAKEN_KARSILASTIRMA.md neredeyse bos (%d bayt)' % os.path.getsize(y)
    m = io.open(y, encoding='utf-8', errors='ignore').read()
    eksik = []
    for anahtar, ad in ((u'PlusPFP kosusu YOK', u'PlusPFP sutunu'),
                        (u'esik taramasi yok', u'esik sutunu'),
                        (u'kimlik_sonuc.csv bulunamadi', u'hizalama sutunu')):
        if anahtar in m:
            eksik.append(ad)
    if eksik:
        return 'eksik', (u'tablo uretildi ama EKSIK: %s. Veri gelince ayni asama '
                         u'tabloyu tamamlar; bu yuzden bitti sayilmadi.'
                         % ', '.join(eksik))
    return True, u'karsilastirma tablosu TAM uretildi (%d bayt)' % os.path.getsize(y)



def d_ozet(kok, ayar):
    y = os.path.join(kok, 'KAPSAMLI_ARAMA_SONUC', '00_OZET_HEPSI.md')
    if not os.path.exists(y):
        return False, u'00_OZET_HEPSI.md uretilmedi'
    yas = time.time() - os.path.getmtime(y)
    if yas > 3600:
        return False, (u'00_OZET_HEPSI.md bu kosuda YENILENMEMIS '
                       u'(%d dakika once yazilmis)' % (yas / 60))
    return True, u'ozet yenilendi'


def d_kraken_ortam(kok, ayar):
    """W: kraken2 var mi.

    kraken2 YOKSA bu bir HATA DEGILDIR - sonraki kraken adimlari atlanir ve
    zincir devam eder. Ama 'bitti' de sayilmaz: kullanici kraken2'yi sonradan
    kurdugunda bu asamanin YENIDEN kosmasi gerekir. Bu yuzden 'eksik' donulur;
    devam mantigi yalniz 'bitti' damgalilari atlar.
    """
    if ayar.get('kraken_var'):
        return True, u'kraken2 bulundu (%s), veritabani denetimi gecti' % (
            ayar.get('kraken2_bin') or 'PATH')
    return 'eksik', (u'kraken2 ya da veritabani yok - X, Y ve Z atlanacak. '
                     u'Kurulunca bu asama yeniden denenir. %s'
                     % ayar.get('kraken_sebep', ''))


ASAMALAR = [
    ('8', u'KENDINI SINA - kod kendini dogrular, olcum yapmaz', u'Grup 4',
     (1, 2), False,
     lambda kok, a: [_py(kok, '-m', 'KAPSAMLI_ARAMA', '--sina')],
     d_sina),

    ('H', u'HIZLI TUTARLILIK TESTI - kucuk alt kumede zincir butun mu', u'Grup 4',
     (25, 40), False,
     lambda kok, a: [_py(kok, os.path.join('KURTARMA', 'hizli_tutarlilik_testi.py'),
                         '--kok', '.')],
     d_hizli_test),

    # W'nin kraken bayragi bilerek False'tur: bu adim TANI adimidir. kraken2
    # yoksa ATLANMAMALI, tam tersine KOSMALIDIR - cunku ekranda neyin eksik
    # oldugunu ve nasil kurulacagini yazan yer burasidir. Atlansaydi kullanici
    # hicbir zaman kurulum yonergesini gormezdi.
    ('W', u'KRAKEN2 ORTAM DENETIMI - kurulu mu, hangi veritabani', u'Grup 3',
     (1, 5), False,
     lambda kok, a: [['bash', a['karac'], 'durum'],
                     ['bash', a['karac'], 'vt-ara'],
                     ['bash', a['karac'], 'vt-kimlik'],
                     ['bash', a['karac'], 'sinav']],
     d_kraken_ortam),

    # OLCULDU 2026-08-05 temiz kosu: 3 sa 20 dk (KISA_LISTE 500 + idf/BM25
    # siralama). Tahmin olculen degere cekildi; eski 6-8 saat KISA_LISTE=1000
    # varsayimina aitti ve gercegin iki kati genisti.
    ('I', u'KIMLIK DOGRULAMA - hocaya giden iddialar sinanir', u'Grup 2',
     (180, 260), False,
     lambda kok, a: [_py(kok, os.path.join('KURTARMA', 'kimlik_dogrulama.py'),
                         '--kok', '.')],
     d_dosya_dolu([os.path.join('KIMLIK_SONUC', 'kimlik_iddialari.tsv')])),

    # OLCULDU 2026-08-05 temiz kosu: 4 sa 43 dk. Tahmin olculen degere cekildi.
    ('G', u'TUM KUTU KIMLIKLERI - panele giren her kutu dogrulanir', u'Grup 2',
     (260, 350), False,
     lambda kok, a: [_py(kok, os.path.join('KURTARMA', 'tum_kutu_kimlikleri.py'),
                         '--kok', '.')],
     d_dosya_dolu([os.path.join('TUM_KIMLIK_SONUC', 'tum_kutu_kimlikleri.tsv')])),

    # Olculen: P 36 sn, K 5 dk, D 6 dk. Onceki tahmin (6-16 saat) fazla genisti.
    ('T', u'TAM OLCUM - P, K, D ve I bagimlilik sirasiyla', u'Grup 1',
     (30, 90), False,
     lambda kok, a: [_py(kok, os.path.join('KURTARMA', 'hepsini_kos.py'), '--kok', '.')],
     d_dosya_dolu([os.path.join('TUM_KOSU_SONUC', '00_BIRLESIK_OZET.md')])),

    # OLCULDU 2026-08-05 temiz kosu: 1 sa 55 dk (6 esik x ~19 dk). Eski tahmin
    # 10-40 dk idi ve UC KAT sapiyordu: ornekleme 100 000 okumaya ayarli ama
    # kaynak 86 780 okuma oldugu icin ornekleme HIC DEVREYE GIRMEDI, tarama tam
    # veriyle kostu. Tahmin gerceklesen degere cekildi.
    ('X', u'KRAKEN GUVEN ESIGI TARAMASI (ornek uzerinde)', u'Grup 3',
     (100, 140), True,
     lambda kok, a: [['bash', a['karac'], 'esik']],
     d_esik),

    # Y AYRI TUTULUR: tam veritabaniyla yeniden siniflandirma tek seferlik ve
    # agirdir. Gece birakilir. O da mmap ve dusuk is parcacigi ile kosar.
    ('Y', u'PlusPFP ILE YENIDEN KOSU (agir, gece birakilir)', u'Grup 3',
     (120, 480), True,
     lambda kok, a: [['bash', a['karac'], 'vt-kimlik'],
                     ['bash', a['karac'], 'esik'],
                     ['bash', a['karac'], 'tablo']],
     d_esik),

    ('Z', u'DORT SUTUNLU KARSILASTIRMA TABLOSU', u'Grup 3',
     (1, 2), True,
     lambda kok, a: [['bash', a['karac'], 'tablo']],
     d_tablo),

    ('S', u'BIRLESIK OZETI YENILE - olcum yapmaz', u'Grup 4',
     (1, 1), False,
     lambda kok, a: [_py(kok, '-m', 'KAPSAMLI_ARAMA', '--mod', 'ozet')],
     d_ozet),
]


# ---------------------------------------------------------------------------
def sure_metni(dk):
    if dk < 60:
        return u'%d dk' % dk
    s, k = divmod(int(dk), 60)
    return u'%d sa %d dk' % (s, k) if k else u'%d saat' % s


def sn_metni(sn):
    sn = int(sn)
    if sn < 60:
        return u'%d sn' % sn
    if sn < 3600:
        return u'%d dk %d sn' % (sn // 60, sn % 60)
    return u'%d sa %d dk' % (sn // 3600, (sn % 3600) // 60)


def kraken2_bul(karac, ortam):
    """kraken2'nin tam yolunu ARACIN KENDISINE sordurur.

    NEDEN BOYLE (2026-08-04 duzeltmesi)
    Eskiden burada shutil.which('kraken2') vardi, yani yalniz PATH'e bakiliyordu.
    Ama bu projede kraken2 PATH'te DEGILDIR: Ali'nin kurulumunda micromamba'nin
    "mikro" ortaminda durur. Sonuc su oldu: ikili diskte oldugu halde butun
    Kraken adimlari "PATH uzerinde yok" diye atlandi.

    Duzeltme, aramayi burada yeniden yazmak DEGIL, aracin kendi arama mantigini
    cagirmaktir. 130_KRAKEN_ARAC.sh icindeki ortam_ac fonksiyonu Ali'nin
    86_KRAKEN_YENIDEN.sh betiginden devralinmistir ve ortami etkinlestirmeyi,
    ortam klasorlerine dogrudan bakmayi ve son care olarak ev dizininde aramayi
    sirayla dener. Arama mantigi TEK YERDE durur; iki ayri yerde iki farkli
    cevap uretmesi imkansiz olur.

    Doner: (tam_yol ya da None, aciklama_metni)
    """
    if not os.path.exists(karac):
        return None, u'Kraken araci bulunamadi (%s)' % karac
    cevre = dict(os.environ)
    if ortam:
        cevre['ORTAM'] = ortam
    try:
        p = subprocess.Popen(['bash', karac, 'kraken-yol'],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             env=cevre)
        cikti, hata = p.communicate(timeout=180)
    except Exception as e:
        return None, u'kraken-yol cagrilamadi: %s' % e
    cikti = cikti.decode('utf-8', 'replace')
    hata = hata.decode('utf-8', 'replace')
    alan = {}
    for satir in cikti.splitlines():
        if '=' in satir and satir.split('=', 1)[0].isupper():
            k, v = satir.split('=', 1)
            alan[k.strip()] = v.strip()
    y = alan.get('KRAKEN2_BIN', '')
    if y and os.path.exists(y):
        # yontem ve surum de tasinir: kullanici ekranda hangi yolla bulundugunu
        # ve hangi surumun kosacagini gormeli, tahmin etmemeli
        return y, u'%s | %s | %s' % (
            y, alan.get('KRAKEN_YONTEM', 'yontem bildirilmedi'),
            alan.get('KRAKEN2_SURUM', 'surum okunamadi'))
    ilk = [s for s in hata.splitlines() if s.strip()]
    return None, (ilk[0] if ilk else u'kraken2 bulunamadi')


def kraken_ortami(kok, pluspfp, vt_a, ortam=''):
    """Kraken tarafinin kullanilabilir olup olmadigini OLCER, tahmin etmez.

    Uc sey aranir: aracin kendisi, kraken2 ikilisi ve en az bir veritabani.
    Ucunden biri yoksa kraken adimlari atlanir ve sebebi yazilir; zincir DURMAZ.
    """
    ali = os.path.abspath(os.path.join(kok, '..', 'ALI'))
    karac = os.path.join(ali, 'WSL', '130_KRAKEN_ARAC.sh')
    sebep = []
    kbin, kmesaj = kraken2_bul(karac, ortam)
    if not kbin:
        sebep.append(kmesaj)

    # Veritabani: once verilen yol, sonra ~/k2db, sonra ARACIN kendi aramasi.
    # Aracin vt-ara tusu diski tarar; burada onu tekrar yazmiyoruz.
    vt = vt_a or os.path.join(os.path.expanduser('~'), 'k2db')
    if not os.path.exists(os.path.join(vt, 'hash.k2d')):
        if kbin:
            try:
                p = subprocess.Popen(['bash', karac, 'vt-ara'],
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     env=dict(os.environ, KRAKEN2_BIN=kbin))
                c = p.communicate(timeout=300)[0].decode('utf-8', 'replace')
                for s in c.splitlines():
                    s = s.strip()
                    if s.endswith('hash.k2d'):
                        aday = os.path.dirname(s)
                        if os.path.exists(os.path.join(aday, 'hash.k2d')):
                            vt = aday
                            break
            except Exception:
                pass
        if not os.path.exists(os.path.join(vt, 'hash.k2d')):
            sebep.append(u'Kraken2 veritabani yok (%s icinde hash.k2d aranmisti; '
                         u'arac vt-ara ile diski da taradi)' % vt)
    return {
        'ali': ali,
        'karac': karac,
        'kraken_is': os.path.join(ali, 'SONUCLAR', 'kraken_esik_A'),
        'kraken_var': not sebep,
        'kraken_sebep': u'; '.join(sebep),
        'kraken2_bin': kbin,
        'kraken_mesaj': kmesaj,
        'vt_a': vt,
        'pluspfp': pluspfp,
        'ortam': ortam,
    }


def plan_yaz(yaz, secili, ayar, durum):
    yaz(u'')
    yaz(u'=' * 78)
    yaz(u'  TAM ZINCIR - PLAN')
    yaz(u'=' * 78)
    alt = ust = 0
    yaz(u'')
    yaz(u'  %-3s %-52s %-12s %s' % (u'#', u'ASAMA', u'SURE', u'DURUM'))
    yaz(u'  ' + u'-' * 74)
    for kod, ad, grup, (a, u_), kraken, _k, _d in secili:
        d = durum.get(kod, {})
        if d.get('durum') == 'bitti':
            nd = u'BITTI - atlanacak'
        elif d.get('durum') == 'eksik':
            nd = u'EKSIK - yeniden denenecek'
        elif kraken and not ayar['kraken_var']:
            nd = u'ATLANACAK (kraken2 yok)'
        elif kod == 'Y' and not ayar['pluspfp']:
            nd = u'ATLANACAK (PlusPFP yolu verilmedi)'
        else:
            nd = u'kosulacak'
            alt += a
            ust += u_
        yaz(u'  %-3s %-52s %-12s %s' % (kod, ad[:52], sure_metni(a) + u'-' + sure_metni(u_), nd))
    yaz(u'')
    yaz(u'  TAHMINI TOPLAM SURE: %s ile %s arasi' % (sure_metni(alt), sure_metni(ust)))
    if ayar.get('kraken2_bin'):
        parca = (ayar.get('kraken_mesaj') or '').split(' | ')
        yaz(u'')
        yaz(u'  kraken2      : %s' % ayar['kraken2_bin'])
        if len(parca) >= 3:
            yaz(u'  bulunma yolu : %s' % parca[1])
            yaz(u'  surum        : %s' % parca[2])
        yaz(u'  veritabani   : %s' % ayar.get('vt_a'))
        yaz(u'  (surum tespiti - PlusPF mi PlusPFP mi - W adiminda yapilir)')
    if not ayar['kraken_var']:
        yaz(u'')
        yaz(u'  KRAKEN ADIMLARI ATLANACAK. Sebep:')
        for s in ayar['kraken_sebep'].split('; '):
            yaz(u'    * %s' % s)
        yaz(u'  Zincir bu yuzden DURMAZ; kalan adimlar normal kosar.')
        yaz(u'')
        yaz(u'  KRAKEN KISMINI SONRADAN TAMAMLAYABILIRSINIZ. Zinciri bastan kosmaniz')
        yaz(u'  gerekmez; atlanan adimlar "bitti" damgasi ALMADIGI icin bir sonraki')
        yaz(u'  koşuda kendiliginden yeniden denenir. Yalniz Kraken kismi icin:')
        yaz(u'      menuden W, sonra X, sonra Z')
        yaz(u'  ya da tek komutla:')
        yaz(u'      python3 KURTARMA/tam_zincir.py --kok . --yalniz W,X,Z,S --onayla')
        yaz(u'  Ortam adi farkliysa:  --ortam <ad>   |  ikilinin yolu belliyse: '
            u'KRAKEN2_BIN=/tam/yol/kraken2')
    yaz(u'')
    yaz(u'  Kesilirse ayni secenege yeniden basin: biten asamalar atlanir.')
    yaz(u'=' * 78)
    return alt, ust


def calistir(kok, ayar, secili, durum, dyol, gunluk, yaz, kuru):
    """Asamalari sirayla kosar. Denetim dusen asama zinciri DURDURUR."""
    for kod, ad, grup, dk, kraken, komut_f, denet in secili:
        d = durum.setdefault(kod, {})
        if d.get('durum') == 'bitti':
            yaz(u'\n>> %s  %s\n   ATLANDI - onceki kosuda bitmisti (%s)'
                % (kod, ad, sn_metni(d.get('sure', 0))))
            continue

        if kraken and not ayar['kraken_var']:
            d.update(durum='atlandi', sebep=ayar['kraken_sebep'], sure=0)
            yaz(u'\n>> %s  %s\n   ATLANDI - %s' % (kod, ad, ayar['kraken_sebep']))
            json.dump(durum, io.open(dyol, 'w', encoding='utf-8'),
                      ensure_ascii=False, indent=1)
            continue

        if kod == 'Y' and not ayar['pluspfp']:
            d.update(durum='atlandi', sebep=u'PlusPFP veritabani yolu verilmedi', sure=0)
            yaz(u'\n>> %s  %s\n   ATLANDI - PlusPFP yolu verilmedi '
                u'(--pluspfp ile verilebilir)' % (kod, ad))
            json.dump(durum, io.open(dyol, 'w', encoding='utf-8'),
                      ensure_ascii=False, indent=1)
            continue

        yaz(u'\n>> %s  %s' % (kod, ad))
        yaz(u'   basladi: %s' % time.strftime('%H:%M:%S'))
        t0 = time.time()
        cevre = dict(os.environ)
        # kraken2'nin TAM YOLU her kraken adimina acikca gecirilir. PATH'e
        # guvenilmez: alt surec yeni bir kabuk acar ve micromamba ortami o
        # kabukta etkin olmayabilir. Tam yol verilince bu sorun ortadan kalkar.
        if ayar.get('kraken2_bin'):
            cevre['KRAKEN2_BIN'] = ayar['kraken2_bin']
        if ayar.get('ortam'):
            cevre['ORTAM'] = ayar['ortam']
        # WSL2'de tepe bellegi dusuk tutan iki ayar. Burada acikca gecirilir ki
        # kullanicinin kabuk ortamina bagli kalmasin.
        cevre.setdefault('IPLIK', os.environ.get('IPLIK', '3'))
        cevre.setdefault('ZORLA_MMAP', '1')
        if kod in ('X', 'Y'):
            cevre.setdefault('ORNEK', os.environ.get('ORNEK', '100000'))
        if kod == 'Y':
            cevre['VT_A'] = ayar['pluspfp']
        elif kod in ('X', 'Z', 'W') and ayar.get('vt_a'):
            cevre['VT_A'] = ayar['vt_a']

        rc, son_cikti = 0, ''
        if kuru:
            yaz(u'   [KURU KOSU] komut calistirilmadi')
        else:
            for argv in komut_f(kok, ayar):
                yaz(u'   $ %s' % ' '.join(os.path.basename(x) for x in argv))
                try:
                    p = subprocess.Popen(argv, cwd=kok, env=cevre,
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.STDOUT)
                except Exception as e:
                    rc, son_cikti = 127, u'komut baslatilamadi: %s' % e
                    yaz(u'   %s' % son_cikti)
                    break
                bicik = []
                for ham in p.stdout:
                    s = ham.decode('utf-8', 'replace').rstrip('\n')
                    bicik.append(s)
                    gunluk.write(s + u'\n')
                p.wait()
                son_cikti = u'\n'.join(bicik)
                rc = p.returncode
                if rc != 0:
                    break
        sure = time.time() - t0
        ayar['_son_cikti'] = son_cikti

        tamam, mesaj = denet(kok, ayar)
        # Ucuncu durum: 'eksik'. Asama kostu, cikti uretti, ama tamamlanmadi.
        # Zincir DURMAZ (yoksa eksik bir Kraken sutunu butun koşuyu keserdi),
        # fakat 'bitti' damgasi da ALMAZ; boylece bir sonraki koşuda yeniden
        # denenir ve veri geldiginde kendiliginden tamamlanir.
        if tamam == 'eksik':
            d.update(durum='eksik', sure=sure, cikis=rc, sebep=mesaj)
            json.dump(durum, io.open(dyol, 'w', encoding='utf-8'),
                      ensure_ascii=False, indent=1)
            yaz(u'   EKSIK (%s) - %s' % (sn_metni(sure), mesaj))
            yaz(u'   Zincir devam ediyor. Bu asama bir sonraki koşuda '
                u'yeniden denenecek.')
            continue
        # Cikis kodu ve cikti denetimi AYRI iki suzgectir. Ikisi de gecmeli.
        #
        # SIFIR DISI CIKIS KODUNA MUSAADE YALNIZ KRAKEN ADIMLARINDADIR (kraken
        # bayragi True olanlar). Kraken sarmalayicisi is bitse de sifir disi
        # donebiliyor; orada denetim gecerse is bitmis sayilir.
        #
        # 2026-08-06 DUZELTMESI - temiz kosuda yakalandi: bu musaade BUTUN
        # asamalara uygulaniyordu. T asamasi cikis kodu 3 verdi (P coktu, K/D/I
        # hic kosmadi) ama beklenen dosya "P dustu" ozeti olarak yazildigi icin
        # denetim gecti ve T "BITTI" damgasi aldi. Genel ozet tablosunda dusen
        # bir asama BASARILI gorundu - zincirin varlik sebebine aykiri.
        # Artik kraken disi bir asamada sifir disi cikis kodu DUSME sayilir.
        if tamam and rc != 0 and not kraken:
            tamam = False
            mesaj = (u'cikis kodu %s (sifir degil). Cikti dosyasi var ama alt '
                     u'asamalarin hepsi kosmamis olabilir - T asamasinda tam '
                     u'boyle olmustu. Yukaridaki asama ciktisini okuyun.' % rc)
        if not tamam:
            d.update(durum='DUSTU', sure=sure, cikis=rc, sebep=mesaj)
            json.dump(durum, io.open(dyol, 'w', encoding='utf-8'),
                      ensure_ascii=False, indent=1)
            yaz(u'   CIKTI DENETIMI DUSTU: %s' % mesaj)
            yaz(u'   cikis kodu: %s | sure: %s' % (rc, sn_metni(sure)))
            yaz(u'\n' + u'=' * 78)
            yaz(u'  ZINCIR DURDU - %s asamasinda' % kod)
            yaz(u'  Sonraki asamalar KOSULMADI. Sebep yukarida.')
            yaz(u'  Duzeltip ayni secenege yeniden basin; biten asamalar atlanir.')
            yaz(u'=' * 78)
            return False
        if rc != 0:
            yaz(u'   uyari: cikis kodu %s, ama cikti denetimi gecti' % rc)
        d.update(durum='bitti', sure=sure, cikis=rc, sebep=mesaj)
        json.dump(durum, io.open(dyol, 'w', encoding='utf-8'),
                  ensure_ascii=False, indent=1)
        yaz(u'   BITTI (%s) - %s' % (sn_metni(sure), mesaj))
    return True


def ozet_yaz(kok, CIKTI, ayar, secili, durum, kesildi):
    """Butun asamalarin sonucunu TEK dosyada birlestirir.

    Olcum yapmaz; durum.json ile asamalarin kendi ciktilarini okur. Bu yuzden
    zincir yarida kesilse bile cagrilabilir ve elde ne varsa onu gosterir.
    """
    yol = os.path.join(CIKTI, '00_TAM_ZINCIR_OZET.md')
    L = []
    A = L.append
    A(u'# Tam zincir - birleşik özet\n')
    A(u'Üretim: %s · sürüm %s\n' % (time.strftime('%Y-%m-%d %H:%M'), VERSIYON))
    if kesildi:
        A(u'> **ZİNCİR TAMAMLANMADI.** Aşağıdaki tabloda DÜŞTÜ işaretli aşamada '
          u'durdu. Sonraki aşamalar koşulmadı.\n')

    A(u'## Aşamaların durumu\n')
    A(u'| # | Aşama | Grup | Durum | Süre | Not |')
    A(u'|---|---|---|---|---|---|')
    for kod, ad, grup, dk, kr, _k, _d in secili:
        d = durum.get(kod, {})
        du = d.get('durum', 'koşulmadı')
        etiket = {'bitti': 'BİTTİ', 'atlandi': 'ATLANDI',
                  'eksik': 'EKSİK (yeniden denenecek)',
                  'DUSTU': '**DÜŞTÜ**'}.get(du, du)
        A(u'| %s | %s | %s | %s | %s | %s |'
          % (kod, ad, grup, etiket,
             sn_metni(d.get('sure', 0)) if d.get('sure') else '-',
             (d.get('sebep') or '')[:150]))

    bitti = [k for k, v in durum.items() if v.get('durum') == 'bitti']
    atlanan = [k for k, v in durum.items() if v.get('durum') == 'atlandi']
    dusen = [k for k, v in durum.items() if v.get('durum') == 'DUSTU']
    eksikler = [k for k, v in durum.items() if v.get('durum') == 'eksik']
    A(u'\n**Özet:** %d aşama bitti, %d atlandı, %d eksik, %d düştü.\n'
      % (len(bitti), len(atlanan), len(eksikler), len(dusen)))
    if eksikler:
        A(u'*Eksik* işaretli aşamalar koştu ama tamamlanmadı; `bitti` damgası '
          u'almadıkları için aynı seçeneğe yeniden basıldığında kendiliğinden '
          u'yeniden denenirler. Veri sonradan geldiğinde tamamlanırlar.\n')

    if atlanan:
        A(u'## Atlanan adımlar\n')
        A(u'Bu adımlar koşulmadı. Zincir bu yüzden durmadı, ama sonuçları '
          u'da elde yok.\n')
        for k in sorted(atlanan):
            A(u'* **%s** — %s' % (k, durum[k].get('sebep', '')))
        A(u'')
        if not ayar['kraken_var']:
            A(u'### Kraken kısmını sonradan tamamlamak\n')
            A(u'Zinciri baştan koşmanız **gerekmez.** Atlanan aşamalar `bitti` damgası '
              u'almadığı için bir sonraki koşuda kendiliğinden yeniden denenir. '
              u'Yalnız Kraken kısmını koşmak için menüden **W**, sonra **X**, sonra '
              u'**Z** tuşlarına basın; ya da tek komutla:\n')
            A(u'```\npython3 KURTARMA/tam_zincir.py --kok . --yalniz W,X,Z,S --onayla\n```\n')
            A(u'Z aşaması tabloyu **sıfırdan değil, eldeki veriyle** kurar: eksik '
              u'sütunları boş bırakıp hangilerinin eksik olduğunu yazar, veri '
              u'sonradan geldiğinde aynı tuş tabloyu tamamlar. Bu yüzden Kraken '
              u'ölçümü haftalar sonra yapılsa bile tablo yeniden üretilebilir.\n')
            A(u'kraken2 kurulu değilse: `micromamba create -n mikro -c bioconda '
              u'-c conda-forge kraken2 bracken`. Zaten kuruluysa ortam adı farklı '
              u'olabilir; `micromamba env list` ile bakıp `--ortam <ad>` verin, ya da '
              u'ikilinin tam yolunu `KRAKEN2_BIN=/tam/yol/kraken2` ile geçin.\n')

    A(u'## Nereye bakılacak\n')
    A(u'| Soru | Dosya |')
    A(u'|---|---|')
    bakilacak = [
        (u'Zincir tutarlı mıydı', 'HIZLI_TEST/HIZLI_TEST_RAPORU.md'),
        (u'Kimlik iddiaları ne oldu', 'KIMLIK_SONUC/KIMLIK_DOGRULAMA_RAPORU.md'),
        (u'Bütün kutuların kimliği', 'TUM_KIMLIK_SONUC/TUM_KUTU_KIMLIK_RAPORU.md'),
        (u'Ölçüm ve kurtarma sonucu', 'TUM_KOSU_SONUC/00_BIRLESIK_OZET.md'),
        (u'Ne sipariş edeyim', 'SIPARIS_LISTESI.tsv'),
        (u'Kraken eşik eğrisi', '../ALI/SONUCLAR/kraken_esik_A/esik_egrisi.txt'),
        (u'Hocaya gidecek Kraken tablosu', '../ALI/0_TESLIM_HOCA/KRAKEN_KARSILASTIRMA.md'),
    ]
    for soru, dy in bakilacak:
        t = os.path.join(kok, dy)
        A(u'| %s | `%s`%s |' % (soru, dy, u'' if os.path.exists(t) else u' *(yok)*'))

    A(u'\n## Bu özet ne değildir\n')
    A(u'Bu dosya aşamaların **koştuğunu** ve çıktılarının **boş olmadığını** '
      u'gösterir. Ölçümlerin **doğru** olduğunu göstermez; doğruluk her aşamanın '
      u'kendi raporunda tartışılır. Sipariş kararı bu özete değil, yukarıdaki '
      u'tabloda adı geçen raporlara bakılarak verilir.\n')

    io.open(yol, 'w', encoding='utf-8').write(u'\n'.join(L) + u'\n')
    return yol


def main():
    p = argparse.ArgumentParser(description=u'Tam zincir: 8 H W I G T X Y Z S')
    p.add_argument('--kok', default='.')
    p.add_argument('--onayla', action='store_true',
                   help=u'onay sormadan basla (menuden gelirken kullanilir)')
    p.add_argument('--yeniden', action='store_true',
                   help=u'durum.json sifirlanir, her sey bastan kosar')
    p.add_argument('--sifirdan', action='store_true',
                   help=u'TEMIZ KOSU: yalniz bu betigin degil, CAGRILAN her '
                        u'asamanin kontrol noktalarini da gecersiz kilar. '
                        u'Hicbir sey silinmez, zaman damgali klasorlere tasinir.')
    p.add_argument('--yalniz', default='',
                   help=u'yalniz bu asamalar, virgulle: ornek 8,S')
    p.add_argument('--atla', default='',
                   help=u'bu asamalar atlanir, virgulle')
    p.add_argument('--pluspfp', default=os.environ.get('PLUSPFP', ''),
                   help=u'PlusPFP veritabani yolu (verilmezse Y adimi atlanir)')
    p.add_argument('--vt', default=os.environ.get('VT_A', ''),
                   help=u'Kraken2 veritabani yolu (verilmezse ~/k2db, sonra arac diski tarar)')
    p.add_argument('--ortam', default=os.environ.get('ORTAM', ''),
                   help=u'micromamba/conda ortam adi (varsayilan: mikro). '
                        u'kraken2 baska bir ortamdaysa burada verin; '
                        u'ortam adlarini gormek icin: micromamba env list')
    p.add_argument('--kuru', action='store_true',
                   help=u'komutlari CALISTIRMADAN plani ve denetimi gosterir')
    p.add_argument('--plan', action='store_true',
                   help=u'YALNIZ plani basar ve cikar; hicbir sey kosulmaz. '
                        u'KAPSAMLI_ARAMA.bat once bunu cagirir, onayi kendi alir, '
                        u'sonra --onayla ile asil kosuyu baslatir. Boylece onay '
                        u'sorusu WSL yerine Windows tarafinda sorulur ve stdin '
                        u'aktariminin bicimine bagli kalmaz.')
    a = p.parse_args()

    kok = os.path.abspath(a.kok)
    if not os.path.isdir(os.path.join(kok, 'KAPSAMLI_ARAMA')):
        sys.stderr.write(
            u'HATA: %s icinde KAPSAMLI_ARAMA yok. Bu betik proje kokunden '
            u'calisir; kok, KAPSAMLI_ARAMA.bat ile ayni klasordur.\n' % kok)
        return 1

    CIKTI = os.path.join(kok, CIKTI_ADI)
    os.makedirs(CIKTI, exist_ok=True)
    dyol = os.path.join(CIKTI, 'durum.json')
    durum = {}
    if os.path.exists(dyol) and not a.yeniden:
        try:
            durum = json.load(io.open(dyol, encoding='utf-8'))
        except Exception:
            durum = {}
    if a.yeniden:
        # DIKKAT: dosya SILINMEZ, UZERINE YAZILIR.
        # Sebep olculdu: bagli klasorlerde ve bazi Windows kurulumlarinda silme
        # izni olmayabiliyor (Operation not permitted). Silmeye calismak
        # --yeniden secenegini tumden kullanilamaz hale getiriyordu. Bos sozluk
        # yazmak ayni isi gorur ve her dosya sisteminde calisir.
        durum = {}
        json.dump(durum, io.open(dyol, 'w', encoding='utf-8'),
                  ensure_ascii=False, indent=1)

    if a.sifirdan:
        # TEMIZ KOSU. --yeniden yalniz BU betigin durum.json'unu sifirlar; ama
        # cagrilan asamalarin KENDI kontrol noktalari vardir ve onlar durdukca
        # asama "zaten bitmis" diye kendi icinde atlar. Gercekten bastan kosmak
        # icin onlar da gecersiz kilinmalidir.
        # SILME YOK: bagli klasorde silme izni olmayabiliyor ve zaten veri
        # kaybetmek istemiyoruz. Klasorler zaman damgasiyla yeniden adlandirilir.
        import datetime as _dt
        damga = _dt.datetime.now().strftime('%Y%m%d_%H%M%S')
        hedefler = ['TEK_PROTOKOL_SONUC', 'KURTARMA_SONUC', 'DOGRULAMA_SONUC',
                    'KIMLIK_SONUC', 'TUM_KIMLIK_SONUC', 'TUM_KOSU_SONUC',
                    'HIZLI_TEST', 'ERISIM_SONUC', CIKTI_ADI]
        tasinan = []
        for h in hedefler:
            y = os.path.join(kok, h)
            if os.path.isdir(y):
                yeni_ad = os.path.join(kok, '%s_ONCEKI_%s' % (h, damga))
                try:
                    os.rename(y, yeni_ad)
                    tasinan.append(h)
                except OSError as e:
                    print(u'  UYARI: %s tasinamadi (%s). Bu asama kendi kontrol '
                          u'noktasindan devam edebilir.' % (h, e))
        ali_is = os.path.join(kok, '..', 'ALI', 'SONUCLAR', 'kraken_esik_A')
        if os.path.isdir(ali_is):
            try:
                os.rename(ali_is, ali_is + '_ONCEKI_' + damga)
                tasinan.append('ALI/SONUCLAR/kraken_esik_A')
            except OSError:
                pass
        print(u'\nTEMIZ KOSU: %d sonuc klasoru kenara alindi (silinmedi).' % len(tasinan))
        for h in tasinan:
            print(u'    %s -> %s_ONCEKI_%s' % (h, h, damga))
        print(u'  Eski sonuclara bakmak isterseniz o klasorler duruyor.\n')
        durum = {}
        os.makedirs(CIKTI, exist_ok=True)

    ayar = kraken_ortami(kok, a.pluspfp, a.vt, a.ortam)

    yalniz = [x.strip().upper() for x in a.yalniz.split(',') if x.strip()]
    atla = [x.strip().upper() for x in a.atla.split(',') if x.strip()]
    secili = [s for s in ASAMALAR
              if (not yalniz or s[0] in yalniz) and s[0] not in atla]
    if not secili:
        sys.stderr.write(u'HATA: secili asama kalmadi (--yalniz / --atla).\n')
        return 1

    gunluk = io.open(os.path.join(CIKTI, 'kosu_gunlugu.txt'), 'a', encoding='utf-8')

    def yaz(s=u''):
        print(s, flush=True)
        gunluk.write(s + u'\n')

    yaz(u'\nTAM ZINCIR  surum %s  %s' % (VERSIYON, time.strftime('%Y-%m-%d %H:%M')))
    plan_yaz(yaz, secili, ayar, durum)

    if a.plan:
        gunluk.close()
        return 0

    if not a.onayla and not a.kuru:
        try:
            c = input(u'\n  Baslatilsin mi? [E/h] ').strip().lower()
        except EOFError:
            c = 'e'
        if c not in ('', 'e', 'evet', 'y', 'yes'):
            yaz(u'  Vazgecildi. Hicbir sey kosulmadi.')
            return 0

    t0 = time.time()
    tamamlandi = calistir(kok, ayar, secili, durum, dyol, gunluk, yaz, a.kuru)
    yol = ozet_yaz(kok, CIKTI, ayar, secili, durum, not tamamlandi)

    yaz(u'\n' + u'=' * 78)
    yaz(u'  TOPLAM SURE: %s' % sn_metni(time.time() - t0))
    yaz(u'  BIRLESIK OZET: %s' % yol)
    yaz(u'=' * 78)
    gunluk.close()
    return 0 if tamamlandi else 3


if __name__ == '__main__':
    sys.exit(main() or 0)
