# -*- coding: utf-8 -*-
u"""Cross-check, an INDEPENDENT, READ-ONLY audit of a finished run.

Opens no measurement of its own that the pipeline already made; it re-asks the
pipeline's questions using different code and reports where the answers differ.
It NEVER writes to panel files.

Seven modules: identity, internal consistency, membership, literature, error
patterns, database health, taxon coverage. Findings are graded KRITIK / CIDDI /
UYARI / BILGI, and checks that COULD NOT RUN are reported as ATLANDI, never as
passed. That distinction is the point: a check that did not run is not a check
that succeeded, and the exit code reflects it.

--- ozgun aciklama ---

CAPRAZ KONTROL - PrimerJury paneli icin BAGIMSIZ, SALT OKUNUR denetci.
===========================================================================

NE YAPAR
    Elimizdeki butun iddialari (a) kendi kaynaklarina ve (b) birbirlerine karsi
    sinar. Yedi modul var; her biri ayri ayri da kosulabilir.

      1  KIMLIK              her kutunun tur duzeyi kimligi + Kraken karsilastirmasi
      2  IC TUTARLILIK       ayni sayi farkli dosyalarda ayni mi
      3  UYELIK BUTUNLUGU    her ciftin uye kumesi tanimli mi
      4  LITERATUR KURALLARI 55 kaynaktan cikan sayisal kurallarin panele uygulanmasi
      5  BILINEN HATA DESENI bu projede tekrar tekrar cikan dokuz desen
      6  VERITABANI SAGLIGI  indeksler gercekten calisiyor mu
      7  TAKSON KAPSAMI      toplanti kararlari panelde karsilik buluyor mu

NE YAPMAZ  --  BU ONEMLI
    HICBIR DOSYAYI DEGISTIRMEZ. Yalnizca okur ve rapor yazar. Tek yazdigi yer
    kendi cikti klasorudur (varsayilan KONTROL_SONUC/). Panel dosyalarina,
    betiklere, veritabanlarina, konsensuslere DOKUNMAZ.

    Kosarken degistirilmemesi gereken dosyalar (baska oturumlar kullaniyor):
    verification/full_chain.py, verification/one_key.py, KONSENSUS_YENIDEN/.
    Bu betik onlarin hicbirine yazmaz; KONSENSUS_YENIDEN/ klasorunu okumaz bile.

SESSIZ ATLAMA YOKTUR
    Bir kontrol kosamazsa "ATLANDI" olarak SAYILIR, sebebi yazilir, ozette
    gorunur ve CIKIS KODU SIFIR OLMAZ. Bir kontrolun her zaman gecmesi, o
    kontrolun bir sey olctugu anlamina gelmez; bu yuzden --kendini-sina
    bayragi her module bilerek bozuk bir girdi verip hatayi yakalayip
    yakalamadigini gosterir.

CIKIS KODU  (bit maskesi, toplanir)
    0  temiz
    1  en az bir KRITIK bulgu
    2  en az bir CIDDI bulgu
    4  en az bir kontrol ATLANDI
    8  betigin kendisi cokti (beklenmeyen hata)
  Ornek: 6 = CIDDI bulgu var VE atlanan kontrol var.

KULLANIM
    python cross_check.py --kok .
    python cross_check.py --kok . --moduller 2,3,5      (yalniz secilenler)
    python cross_check.py --kok . --m1-kip tam          (agir kimlik taramasi)
    python cross_check.py --kok . --kendini-sina        (bozuk girdi sinamasi)
    python3 cross_check.py --kok .

Yazan: bu oturum, 2026-08-09.
"""
from __future__ import print_function

import argparse
import ast
import collections
import glob
import hashlib
import io
import json
import math
import os
import re
import sys
import time
import traceback

VERSIYON = u'1.0 (2026-08-09)'

# UTF-8 ciktisi: Windows konsolu varsayilan olarak cp857/cp1254 kullanir ve
# Turkce karakterlerde UnicodeEncodeError atar. Betik ortasinda cokmek yerine
# akisi bastan sarmaliyoruz. (Hata YUTULMAZ, yalnizca kodlama duzeltilir.)
if hasattr(sys.stdout, 'buffer'):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace', line_buffering=True)
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8',
                                  errors='replace', line_buffering=True)


# ===========================================================================
# CIDDIYET DUZEYLERI
# ===========================================================================
# KRITIK : siparisi ya da raporlanacak sayiyi DOGRUDAN yanlis yapar
# CIDDI  : bir iddianin dayanagi cokuyor ama siparis hemen degismeyebilir
# UYARI  : tutarsiz ya da savunmasiz ama olculen sayi degismiyor
# BILGI  : not duselim, hata degil
# ATLANDI: kontrol KOSULAMADI. Bu bir "gecti" DEGILDIR.
KRITIK, CIDDI, UYARI, BILGI, ATLANDI = u'KRITIK', u'CIDDI', u'UYARI', u'BILGI', u'ATLANDI'
_SIRA = {KRITIK: 0, CIDDI: 1, UYARI: 2, BILGI: 3, ATLANDI: 4}


class Bulgu(object):
    u"""Tek bir denetim bulgusu.

    Her bulgu DORT soruya cevap vermek ZORUNDA, yoksa rapor okunmaz olur:
      beklenen : kuralin ne oldugu
      bulunan  : olculen/okunan gercek deger
      dosya    : bulgunun dayandigi kaynak (yol, mumkunse satir)
      ciddiyet : yukaridaki bes duzeyden biri
    """
    __slots__ = ('modul', 'kod', 'ciddiyet', 'beklenen', 'bulunan', 'dosya', 'oneri')

    def __init__(self, modul, kod, ciddiyet, beklenen, bulunan, dosya, oneri=u''):
        self.modul = modul
        self.kod = kod
        self.ciddiyet = ciddiyet
        self.beklenen = beklenen
        self.bulunan = bulunan
        self.dosya = dosya
        self.oneri = oneri


class Rapor(object):
    u"""Butun modullerin bulgularini toplar, sayar ve cikis kodunu uretir."""

    def __init__(self):
        self.bulgular = []
        self.modul_durumu = collections.OrderedDict()   # modul -> dict(durum, sure, not)
        self.olcum = collections.OrderedDict()          # olculen sureler / oranlar
        self.tablolar = collections.OrderedDict()       # rapora gomulecek ek tablolar

    def ekle(self, modul, kod, ciddiyet, beklenen, bulunan, dosya, oneri=u''):
        self.bulgular.append(Bulgu(modul, kod, ciddiyet, beklenen, bulunan, dosya, oneri))

    def atla(self, modul, kod, beklenen, sebep, dosya=u'-'):
        u"""Bir kontrol kosulamadi. SESSIZ ATLAMA YOK: bulgu olarak kaydedilir."""
        self.bulgular.append(Bulgu(modul, kod, ATLANDI, beklenen,
                                   u'KOSULAMADI: %s' % sebep, dosya,
                                   u'Bu kontrol OLCULMEDI, "gecti" sayilmamalidir.'))

    def say(self, ciddiyet=None, modul=None):
        return len([b for b in self.bulgular
                    if (ciddiyet is None or b.ciddiyet == ciddiyet)
                    and (modul is None or b.modul == modul)])

    def cikis_kodu(self):
        k = 0
        if self.say(KRITIK):
            k |= 1
        if self.say(CIDDI):
            k |= 2
        if self.say(ATLANDI):
            k |= 4
        return k


# ===========================================================================
# TEMEL YARDIMCILAR
# ===========================================================================
def unicode_(x):
    if isinstance(x, bytes):
        return x.decode('utf-8', 'replace')
    return u'%s' % (x,)


def yaz(*a):
    u"""Ilerleme satiri. Zaman damgali, ciktisi hemen bosaltilir."""
    print(u'[%s] %s' % (time.strftime('%H:%M:%S'), u' '.join(unicode_(x) for x in a)))
    try:
        sys.stdout.flush()
    except Exception:
        pass


def vir(x, b=2):
    u"""Turkce ondalik. None -> '-'. '0' ile 'olculmedi' birbirine KARISMAZ."""
    if x is None:
        return u'-'
    try:
        return (u'%.*f' % (b, float(x))).replace(u'.', u',')
    except (TypeError, ValueError):
        return unicode_(x)


def sayi(x):
    u"""Turkce ya da Ingilizce ondalikli metni float'a cevir. Cevrilemezse None.

    None donmesi ONEMLI: cevrilemeyen bir hucreyi 0 saymak, "olculmedi" ile
    "sifir olculdu" arasindaki farki yok eder - bu projede tam olarak bu
    yuzden bir ayrim kati 0,00x sanilmisti.
    """
    if x is None:
        return None
    s = unicode_(x).strip().replace(u'%', u'').replace(u'x', u'')
    s = s.replace(u'−', u'-')          # unicode eksi isareti
    if not s or s.lower() in (u'-', u'--', u'yok', u'olculmedi', u'nd', u'nan', u'?'):
        return None
    # "1.234,56" (TR binlik) ile "1,234.56" (EN binlik) ayrimi: SON ayraç ondalik
    if u',' in s and u'.' in s:
        s = s.replace(u'.', u'') if s.rfind(u',') > s.rfind(u'.') else s.replace(u',', u'')
    s = s.replace(u',', u'.')
    try:
        return float(s)
    except ValueError:
        return None


def md5_metin(*parcalar):
    u"""Kontrol noktasi anahtari. hash() KULLANILMAZ.

    Python'un yerlesik hash()'i surec basina TUZLANIR (PYTHONHASHSEED); ayni
    girdi iki kosuda farkli anahtar uretir ve kontrol noktasi ya hic tutmaz ya
    da beteri, cakisir. md5 kararlidir ve kosudan kosuya ayni kalir.
    """
    h = hashlib.md5()
    for p in parcalar:
        h.update(unicode_(p).encode('utf-8', 'replace'))
        h.update(b'\x00')
    return h.hexdigest()


def dosya_imzasi(yol):
    u"""Bir girdinin kimligi: mutlak yol + boyut + degistirilme zamani.

    Kontrol noktasi anahtarina bu girer; girdi degisince anahtar degisir ve
    eski kontrol noktasi kendiliginden gecersizlesir.
    """
    try:
        st = os.stat(yol)
        return u'%s|%d|%d' % (os.path.abspath(yol), st.st_size, int(st.st_mtime))
    except OSError:
        return u'%s|YOK' % os.path.abspath(yol)


class KontrolNoktasi(object):
    u"""md5 anahtarli, girdi tazeligini denetleyen kontrol noktasi deposu.

    IKI KATMANLI GECERSIZLESTIRME:
      1) Anahtar, girdi dosyalarinin imzasindan (yol+boyut+mtime) turetilir.
         Girdi degisirse anahtar degisir, eski kayit bulunamaz.
      2) Ayrica kayitli dosyanin mtime'i girdilerinkiyle karsilastirilir.
         Girdi kontrol noktasindan YENIYSE kayit GECERSIZ sayilir.
    Ikinci katman, ayni anahtarin elle kopyalanmasi gibi durumlara karsi
    emniyet kemeridir. Bayat kontrol noktasi bu projede daha once "temiz"
    sanilan bir kosuya yol acmisti.
    """

    def __init__(self, klasor, etkin=True):
        self.klasor = klasor
        self.etkin = etkin
        self.isabet = 0
        self.iska = 0
        if etkin and not os.path.isdir(klasor):
            os.makedirs(klasor)

    def _yol(self, anahtar):
        return os.path.join(self.klasor, anahtar + '.json')

    def oku(self, anahtar, girdiler=()):
        if not self.etkin:
            return None
        y = self._yol(anahtar)
        if not os.path.exists(y):
            self.iska += 1
            return None
        try:
            kn_zaman = os.path.getmtime(y)
        except OSError:
            self.iska += 1
            return None
        for g in girdiler:                      # GIRDI DAHA YENIYSE -> GECERSIZ
            try:
                if os.path.exists(g) and os.path.getmtime(g) > kn_zaman:
                    self.iska += 1
                    return None
            except OSError:
                self.iska += 1
                return None
        try:
            with io.open(y, encoding='utf-8') as fh:
                v = json.load(fh)
            self.isabet += 1
            return v
        except (ValueError, IOError) as e:
            # Bozuk kontrol noktasi SESSIZCE yok sayilmaz, ekrana basilir.
            yaz(u'  WARNING: could not read checkpoint (%s): %s' % (anahtar, e))
            self.iska += 1
            return None

    def yazdir(self, anahtar, veri):
        if not self.etkin:
            return
        try:
            with io.open(self._yol(anahtar), 'w', encoding='utf-8') as fh:
                fh.write(unicode_(json.dumps(veri, ensure_ascii=False, default=str)))
        except IOError as e:
            yaz(u'  WARNING: could not write checkpoint (%s): %s' % (anahtar, e))


class Canlilik(object):
    u"""Uzun donguler icin duzenli canlilik isareti.

    Ekrana hicbir sey basmayan bir asama, KILITLENMIS bir asamadan ayirt
    edilemez. Bu sinif en az `aralik` saniyede bir satir basar ve olculen
    hizdan kalan sureyi tahmin eder.
    """

    def __init__(self, etiket, toplam=None, aralik=20.0):
        self.etiket = etiket
        self.toplam = toplam
        self.aralik = aralik
        self.t0 = time.time()
        self.son = self.t0

    def vur(self, n, ek=u''):
        t = time.time()
        if t - self.son < self.aralik:
            return
        self.son = t
        gecen = t - self.t0
        if self.toplam:
            oran = n / float(self.toplam)
            kalan = (gecen / oran - gecen) if oran > 0 else None
            yaz(u'  ... %s %d/%d (%%%s) gecen %s, kalan ~%s %s'
                % (self.etiket, n, self.toplam, vir(100 * oran, 1),
                   sure_metni(gecen), sure_metni(kalan), ek))
        else:
            yaz(u'  ... %s %d, gecen %s %s' % (self.etiket, n, sure_metni(gecen), ek))

    def bitti(self, n):
        d = time.time() - self.t0
        yaz(u'  %s done: %d items, %s' % (self.etiket, n, sure_metni(d)))
        return d


def sure_metni(sn):
    if sn is None:
        return u'?'
    sn = float(sn)
    if sn < 90:
        return u'%d saniye' % int(round(sn))
    if sn < 5400:
        return u'%d dakika' % int(round(sn / 60.0))
    return u'%s saat' % vir(sn / 3600.0, 1)


# ===========================================================================
# SALT OKUNUR KAYNAK OKUYUCULARI
# ===========================================================================
def gecersiz_isareti(yol):
    u"""Dosyanin bas kismindaki '# GECERSIZ' isaretini arar.

    Doner: gerekce metni (isaret varsa) ya da None. Yalniz ilk 10 satira
    bakilir; isaret dosyanin EN USTUNDE olmalidir ki gozden kacmasin.
    """
    try:
        with io.open(yol, encoding='utf-8', errors='replace') as fh:
            for i, s in enumerate(fh):
                if i >= 10:
                    break
                if not s.startswith(u'#'):
                    break
                if u'GECERSIZ' in s.upper():
                    return s.lstrip(u'#').strip()
    except (IOError, OSError):
        return None
    return None


def en_yeni_kaynak(kok, desen):
    u"""desen'e uyan dosyalardan GECERSIZ olmayan EN YENISINI dondurur.

    2026-08-09: kaynak yolu tarihli dosya adina sabitlenmisti
    ('ESIK_VE_OLCUT_2026-08-08.tsv'). Yeni bir surum yazildiginda denetim eski
    surumu okumaya devam ediyordu. Artik desenle bulunup ada gore siralanir;
    ad tarih tasidigi icin en buyuk ad en yeni surumdur. Hicbiri yoksa desenin
    kendisi dondurulur ki "dosya yok" bulgusu yine uretilsin.
    """
    adaylar = sorted(glob.glob(os.path.join(kok, desen)))
    for y in reversed(adaylar):
        if not gecersiz_isareti(y):
            return y
    return adaylar[-1] if adaylar else os.path.join(kok, desen)


def tsv_oku(yol, yorum=u'#'):
    u"""Basligi olan TSV -> [dict]. Yorum satirlari atlanir.

    Dosya yoksa None doner (bos liste DEGIL). "Dosya yok" ile "dosya bos" ayri
    seylerdir ve ayri bulgu uretirler; ikisini birlestirmek sessiz atlamadir.
    """
    if not os.path.exists(yol):
        return None
    # 2026-08-09: GECERSIZ isaretli kaynak OKUNMAZ. Gerekce: 2026-08-09 18:09
    # kosusunda 22 M5-D1-SESSIZ-SIFIR bulgusunun 22'si de artik yerini
    # ESIK_VE_OLCUT_2026-08-09.tsv'ye birakmis olan
    # NIHAI_SIPARIS_LISTESI_2026-08-07.tsv'den geldi. Gecersiz bir karar
    # tablosundan bulgu uretmek denetimi gurultuyle dolduruyor. Dosya
    # SILINMEZ; basligina '# GECERSIZ' satiri konur ve denetim onu kaynak
    # saymaz. Atlama SESSIZ DEGILDIR: cagiran taraf None gorur ve rap.atla()
    # ile ATLANDI bulgusu uretir.
    if gecersiz_isareti(yol):
        return None
    satirlar = []
    with io.open(yol, encoding='utf-8', errors='replace') as fh:
        for s in fh:
            if yorum and s.startswith(yorum):
                continue
            if s.strip():
                satirlar.append(s.rstrip(u'\n').rstrip(u'\r'))
    if not satirlar:
        return []
    bas = [b.strip() for b in satirlar[0].split(u'\t')]
    out = []
    for i, s in enumerate(satirlar[1:], 2):
        p = s.split(u'\t')
        p += [u''] * (len(bas) - len(p))
        d = dict(zip(bas, p[:len(bas)]))
        d['_satir'] = i
        out.append(d)
    return out


def metin_oku(yol):
    if not os.path.exists(yol):
        return None
    with io.open(yol, encoding='utf-8', errors='replace') as fh:
        return fh.read()


def _kod_govdesi(metin):
    u"""Python kaynagindan yorum ve docstring'leri atar, KOD GOVDESINI dondurur.

    Neden gerekli: bir klasor ya da API adinin ACIKLAMADA gecmesi, o kod yolunun
    gercekten kullanildigi anlamina gelmez. Duz metin aramasi bu ikisini ayirt
    edemez ve "hala yamalanmamis kod" sorusunu yanlis cevaplar.

    screening/yon_kod_taramasi.kod_govdesi() ile ayni ISI yapar ama olcut
    duzeltilmistir (2026-08-21, olculdu). Oradaki surum dizeyi SATIR SAYISINA
    gore eliyor:  end_lineno - lineno >= 1.  Bu iki yonden de yanlis:
      * TEK SATIRLIK docstring atilmaz -> aciklamada gecen ad kod sanilir
        (yanlis pozitif; sentetik sinamada dogrulandi).
      * Kodda GERCEKTEN kullanilan cok satirli bir dize atilir -> gomulu yol
        gozden kacar (yanlis negatif).
    Dogru ayirt edici satir sayisi degil, dizenin DOCSTRING olup olmadigidir:
    docstring, degeri dize sabiti olan bir ifade deyimidir (ast.Expr). Kodda
    kullanilan dize her zaman baska bir dugumun cocugudur. Bu yuzden yalniz
    ast.Expr altindaki dize sabitleri atilir.

    Ayristirilamayan dosyada docstring kumesi bos kalir, yani en kotu durumda
    eski (asiri genis) davranisa duser - sessizce kontrolu atlamaz.
    """
    try:
        agac = ast.parse(metin)
        ds = set()
        for d in ast.walk(agac):
            if isinstance(d, ast.Expr) and isinstance(d.value, ast.Constant) \
               and isinstance(d.value.value, str):
                b = getattr(d, 'lineno', None)
                s = getattr(d, 'end_lineno', b)
                if b:
                    ds.update(range(b, (s or b) + 1))
    except Exception:
        ds = set()
    # SATIR NUMARASI KORUNUR: docstring satiri ATILMAZ, BOSALTILIR. Atmak
    # numaralari kaydirir; kaydirinca AST'ten gelen satir kumeleriyle (mesaj
    # dizeleri gibi) hizalama bozulur ve suzgec sessizce islemez - olculdu
    # (screening/run_all.py ozgun 222, kirpilmis govdede 183).
    out = []
    for i, l in enumerate(metin.splitlines(), 1):
        out.append(u'' if i in ds else l.split(u'#', 1)[0])
    return u'\n'.join(out)


# Gorevi GEREGI karisik yonlu klasoru okuyan araclar. Bunlari isaretlemek
# yanlistir: yon denetimi ve kanonik uretimi o klasoru okumadan is goremez.
# Liste KISA tutulmali; her giris bir gerekce tasir.
D9_MESRU_OKUYUCULAR = {
    u'orientation_audit.py':    u'yon denetleyicisi - karisik klasoru olcmek gorevidir',
    u'build_canonical.py':    u'kanonik klasoru o klasorden URETIR',
    u'orientation.py':             u'yon tanimlarinin kaynagi',
    u'orientation_code_scan.py': u'ayni riski tarayan kardes arac',
    u'orientation_report.py':     u'yon kararlarini belgeleyen rapor ureteci',
}

# Ekrana/loga metin basan cagrilar. Bunlarin ICINDEKI dize bir KOD YOLU degil,
# kullaniciya gosterilen bir mesajdir (olculdu: screening/run_all.py:183
# yalnizca "Kaynak: ..." satirini basiyor ve bu yuzden RISKLI sayiliyordu).
D9_CIKTI_CAGRILARI = frozenset(
    [u'yaz', u'print', u'write', u'log', u'uyar', u'bilgi', u'hata', u'yazdir'])


def d9_karisik_klasor_yollari(metin, ad=u'', karisik=u'consensus sequences'):
    u"""Karisik yonlu klasoru GERCEKTEN okuyan kod satirlarini dondurur.

    [(satir_no, kaynak_satiri), ...]  - bos liste "temiz" demektir.

    Uc suzgecten geciyor, hepsinin gerekcesi olculdu (2026-08-21):
      1) docstring ve yorumlar atilir  -> aciklamada gecen ad kod sayilmaz
      2) cikti cagrilarinin argumanlari atilir -> ekrana basilan mesaj kod degil
      3) gorevi geregi okuyan araclar muaf -> D9_MESRU_OKUYUCULAR

    Duz metin aramasi bu ucunu de ayirt edemiyordu ve 2026-08-09 kosusunda bes
    yanlis pozitif uretti.
    """
    if os.path.basename(ad) in D9_MESRU_OKUYUCULAR:
        return []
    try:
        agac = ast.parse(metin)
    except Exception:
        # Ayristirilamiyorsa duz metne duseriz: asiri genis ama sessiz degil.
        return [(i, l) for i, l in enumerate(metin.splitlines(), 1)
                if karisik in l and not l.strip().startswith(u'#')]

    mesaj_satirlari = set()
    for d in ast.walk(agac):
        if not isinstance(d, ast.Call):
            continue
        f = d.func
        adi = getattr(f, 'id', None) or getattr(f, 'attr', None)
        if adi in D9_CIKTI_CAGRILARI:
            for arg in list(d.args) + [k.value for k in d.keywords]:
                for s in ast.walk(arg):
                    if isinstance(s, ast.Constant) and isinstance(s.value, str) \
                       and karisik in s.value:
                        b = getattr(s, 'lineno', None)
                        e = getattr(s, 'end_lineno', b)
                        if b:
                            mesaj_satirlari.update(range(b, (e or b) + 1))

    govde = _kod_govdesi(metin).splitlines()
    return [(i, l) for i, l in enumerate(govde, 1)
            if karisik in l and i not in mesaj_satirlari]


def xlsx_sayfalari(yol):
    u"""Excel -> {sayfa: [[hucre,...],...]}. Hata durumunda metin doner.

    SALT OKUNUR acilir (read_only=True): dosyaya yazma ihtimali sifirdir.
    """
    if not os.path.exists(yol):
        return None
    try:
        import openpyxl
    except ImportError:
        return u'OPENPYXL_YOK'
    try:
        wb = openpyxl.load_workbook(yol, read_only=True, data_only=True)
    except Exception as e:
        return u'ACILAMADI: %s: %s' % (type(e).__name__, e)
    out = collections.OrderedDict()
    try:
        for ad in wb.sheetnames:
            sh = wb[ad]
            out[ad] = [[(u'' if c is None else unicode_(c)) for c in r]
                       for r in sh.iter_rows(values_only=True)]
    finally:
        try:
            wb.close()
        except Exception:
            pass
    return out


def modul_yukle(yol, ad):
    u"""Bir .py dosyasini MODUL olarak yukle (paket degil, betik).

    verification/identity_verification.py ve all_bin_identities.py'nin mantigi burada
    YENIDEN YAZILMAZ; oradan ice aktarilir. Yeniden yazmak iki ayri karar
    mantigi yaratirdi ve hangisinin gecerli oldugu belirsizlesirdi.
    """
    import importlib.util
    if not os.path.exists(yol):
        return None, u'dosya yok: %s' % yol
    try:
        sp = importlib.util.spec_from_file_location(ad, yol)
        m = importlib.util.module_from_spec(sp)
        sys.modules[ad] = m
        sp.loader.exec_module(m)
        return m, None
    except Exception as e:
        return None, u'%s: %s' % (type(e).__name__, e)


# ===========================================================================
# PROJE KAYNAKLARI  -  tek yerden tanimlanir
# ===========================================================================
class Kaynaklar(object):
    u"""Denetimin okudugu butun dosyalarin tek listesi.

    Yol degisirse TEK yerden degisir. Her modul buradan okur; boylece
    "hangi dosyaya baktin" sorusunun cevabi rapordan izlenebilir.
    """

    # BU YOLLARA ASLA YAZILMAZ. Baska oturumlar kosuyor olabilir.
    DOKUNULMAZ = (u'install.sh', u'build_index.sh',
                  os.path.join(u'verification', u'one_key.py'), u'KONSENSUS_YENIDEN')

    def __init__(self, kok):
        self.kok = os.path.abspath(kok)
        y = lambda *p: os.path.join(self.kok, *p)
        self.y = y
        # --- karar tablolari
        # 2026-08-09 20:05 DENENDI VE GERI ALINDI: bu yolu
        # ESIK_VE_OLCUT_2026-08-09.tsv'ye cevirmeyi denedim, cunku 08-07
        # dosyasi GECERSIZ isaretlenince dort kontrol ATLANDI'ya dusuyor.
        # Ama iki dosyanin SUTUN YAPISI ayni degil: ESIK_VE_OLCUT'ta ileri/geri
        # primer sutunu yok, M2 tutarlilik kontrolu oradan sayisal alan okuyup
        # 22 SAHTE KRITIK celiski uretti. Yanlis alarm, durust ATLANDI'dan
        # daha kotudur. Eski yol korundu; o dort kontrol ATLANDI kalir ve
        # bu bilerek boyledir. Kalici cozum: kontrolleri yeni tablonun
        # sutunlarina gore yeniden yazmak, dosya yolunu degistirmek degil.
        self.nihai_siparis = y('NIHAI_SIPARIS_LISTESI_2026-08-07.tsv')
        # 2026-08-09: sabit tarihli ad yerine EN YENI gecerli surum secilir.
        self.esik_olcut = en_yeni_kaynak(self.kok, 'ESIK_VE_OLCUT_*.tsv')
        self.siparis_listesi = y('SIPARIS_LISTESI.tsv')
        self.hedef_disi = y('HEDEF_DISI_AYRINTI_2026-08-07.tsv')
        # --- panel tanimi
        self.ciftler = y('screening', 'ciftler.tsv')
        self.hedef_uyelik = y('screening', 'hedef_uyelik.tsv')
        self.hedef_klad = y('screening', 'hedef_klad.tsv')
        self.hedefler_wsl = y('steps', 'hedefler.tsv')
        self.takson_esleme = y('screening', 'target_taxon_mapping.py')
        # --- excel
        # 2026-08-11: teslim xlsx'i arsive tasindi (icindeki alti ciftin dizisi
        # eskiydi). Yerine her kosuda URETILEN tek dosya kondu; adinda tarih
        # oldugu icin en yenisi secilir. Bulunamazsa eski ad denenir ve
        # denetim "dosya yok" der - sessizce eski dosyaya DUSMEZ.
        import glob as _glob
        _ad = sorted(_glob.glob(y('PrimerJury_PANEL_*.xlsx')))
        # Bulunamazsa BOS birakilir; denetim "dosya yok" der. Eski teslim
        # dosyasina DUSMEK, tam da kacinmak istedigimiz sey: o dosyada alti
        # ciftin dizisi eskiydi ve arsive tasindi.
        self.panel_xlsx = _ad[-1] if _ad else None
        # --- olcum ciktilari
        self.tek_protokol = y('TEK_PROTOKOL_SONUC', 'panel_tek_protokol.tsv')
        self.uyelik_turetme = y('uyelik_yeniden_turetme_uyelik_20260803.tsv')
        self.uyelik_ciftler = y('uyelik_yeniden_turetme_ciftler_20260803.tsv')
        # --- literatur ve toplanti
        self.literatur = y('LITERATUR_2026-08-07.md')
        self.toplanti = y('TOPLANTI_KARARLARI_SON_DURUM.md')
        # --- veri
        self.refdb = y('REFERANS_DB')
        self.konsensus_indeks = y('konsensus_kanonik', 'INDEKS.tsv')
        self.konsensus_kok = y('konsensus_kanonik')
        self.fastq = y('fastq files')
        self.kraken = y('kraken results')
        self.bracken = y('bracken results')
        # --- kod (desen taramasi icin). KONSENSUS_YENIDEN BILEREK YOK.
        self.kod_klasorleri = [y('screening'), y('verification'), y('steps')]
        # --- ice aktarilacak mantik
        self.kimlik_dogrulama = y('verification', 'identity_verification.py')
        self.tum_kutu = y('verification', 'all_bin_identities.py')


# ===========================================================================
# MODUL 1 - KIMLIK
# ===========================================================================
# SORU: her kutuda gercekten hangi organizma var, ve Kraken2'nin dedigi bu mu?
#
# NEDEN AYRI BIR OLCUM: Kraken2 k-mer'e dayali bir SINIFLANDIRICIDIR, hizalama
# yapmaz ve kendi veritabaninda olmayan bir organizmayi en yakin akrabasinin
# adiyla etiketler. Bu projede birden fazla Kraken etiketi olcumle curutuldu
# (Trichoderma asperellum -> Petriella, Dictyostelium discoideum -> etiket
# curutuldu). Bu yuzden her kutu, hizalamayla ve OFFLINE veritabanlarina karsi
# bagimsiz olarak yeniden olculur.
#
# ALTI CIKTI URETILIR, her kutu icin:
#   1 tur duzeyi en iyi isabet  (tur, veritabani, kayit no, kimlik %, hizalanan bp)
#   2 TIP KAYDI isareti         (type material - cok daha guclu kanit, AYRI sutun)
#   3 ikinci en iyi isabet + fark (fark kucukse tur atamasi guvenilir DEGILDIR)
#   4 AYIRT EDILEBILIRLIK       (bizim veri kalitemizle bu iki tur ayrilir mi)
#   5 N orani ve karisim gostergesi
#   6 Kraken karsilastirmasi    (etiket, guven, ad degisti mi)
# ---------------------------------------------------------------------------

# --- kimlik esikleri: identity_verification.py'den ALINIR, burada YENIDEN yazilmaz.
#     (K.TUR_ESIGI, K.CINS_ESIGI, K.AYRIM_PAYI)

# N orani bu esigi asarsa konsensus KULLANILMAZ, ham okumalara donulur.
# GEREKCE (olculdu): F2-4_500148 konsensusunun %43,6'si N. N'ler hizalamada
# daima uyumsuz sayildigi icin bu kutu "adlandirilamayan soy" gorunuyordu;
# okuma bazli bakildiginda kutunun aslinda Microascaceae oldugu ve KARISIK
# oldugu ortaya cikti. Esik %20: bu deger, saglam kutularin en kotusu (%0,2)
# ile bozuk kutu (%43,6) arasinda genis bir bosluga dusuyor.
N_ESIGI = 20.0

# Karisim/saflik olcumu icin ornek buyuklukleri. Ham okuma sayisi binlerce;
# hepsini hizalamak gereksiz. 60 okuma, %10'luk bir alt topluluğu %99'un
# uzerinde bir olasilikla yakalar (1 - 0,9^60 = 0,998).
KARISIM_OKUMA = 60
KARISIM_PENCERE = 700      # okuma-okuma hizalamalarinda kullanilan pencere (bp)
SAF_ESIGI = 90.0           # ESKI KURAL - artik TEK BASINA hukum vermez, bkz asagisi

# ---------------------------------------------------------------------------
# 2026-08-09 DUZELTMESI - KARISIK KUTU KURALI
# ---------------------------------------------------------------------------
# BELIRTI: 2026-08-09 18:09 kosusunda 99 kutunun 62'si KARISIK isaretlendi.
#          Buyuk kismi yanlis alarmdi. Ornek A1-2_1826872: 60 okumanin 53'u
#          tek kumede, kalan 7'si BIRER okumalik yedi ayri kume; baskin %88,3
#          ve tek olcut "baskin >= %90" oldugu icin kutu KARISIK sayildi.
#          Oysa tek okumalik kumeler Nanopore hata orani (olculen eps
#          %1,3-2,4) altinda beklenen gurultudur, ayri organizma degildir.
#
# YENI KURAL: bir kutunun KARISIK sayilmasi icin baskin kume DISINDA
#          ANLAMLI buyuklukte en az bir kume bulunmalidir.
#
# ANLAMLI ESIGI OLCULDU, SECILMEDI:
#   Baskin kumesi cok yuksek olan kutular tanim geregi tek organizmadir;
#   oradaki butun baskin-disi kumeler gurultudur. O kumelerin boy dagilimi
#   olculdu (2026-08-09 kosusunun kume verisinden, 98 kutu):
#
#     temiz tanimi   kutu   gurultu kumesi   en buyuk gurultu kumesi
#     baskin >=%96,5    7               10                        2
#     baskin >=%95,0   12               25                        2
#     baskin >=%93,0   19               50                        2
#     baskin >=%91,5   29               96                        4  <- kirilma
#
#   Uc ayri "temiz" tanimi altinda, toplam 50 gurultu kumesinde gozlenen en
#   buyuk gurultu kumesi 2 OKUMADIR; 3 okumalik bir gurultu kumesi HIC
#   gorulmedi. Ucler kurali ile: P(gurultu kumesi >= 3 okuma) < 3/50 = %6.
#   Bagimsiz ikinci turetme, ayni temiz alt kumeden Poisson: lam = 1,040,
#   P(sahte kume >= 3) = %8,8, yani 98 kutuda beklenen 17,9 kume; gozlenen
#   48. Gozlenen fazlalik gurultuyle aciklanamiyor, yani >=3 okumalik kumeler
#   gercek sinyal tasiyor.
#
#   Bu yuzden esik: ANLAMLI KUME = en az 3 okuma VE ornegin en az %5'i.
#   (Ornek 60 okuma oldugunda 3/60 = %5; oran tabani ornek boyu degisirse
#   kuralin kaymamasi icin var.)
#
# TEMSIL TABANI: baskin kume ornegin YARISINDAN azsa kutu, ikincil kumeleri
#   kucuk olsa bile KARISIK sayilir. Bu bir gurultu esigi degil, temsil
#   sartidir: konsensus baskin kumenin medoidinden turetiliyor; azinlikta
#   kalan bir kumeden turetilen konsensus okumalarin cogunlugunu temsil
#   etmez. %50 secilmis bir parametre degil, cogunlugun tanimidir.
#
# ETKI (olculdu, ayni veri): KARISIK 62 -> 37. Yirmi alti kutu SAF'a dondu,
#   bir kutu (F1-4_101201, 55/4/1) SAF iken KARISIK oldu. Gercekten karisik
#   olanlar yakalanmaya devam ediyor; ornegin baskin kumesi %46,4 olan
#   A2-1_1826872 (26/8/3/2/1...) hala KARISIK.
ANLAMLI_KUME_OKUMA = 3     # olculen gurultu tavani 2 okuma; anlamli = >=3
ANLAMLI_KUME_ORAN = 5.0    # ve ornegin en az %5'i
BASKIN_TEMSIL_TABANI = 50.0  # baskin kume cogunlukta degilse konsensus temsil etmez


def saflik_hukmu(boyut, n_okuma):
    """Kume boyutlarindan saflik hukmu. (saflik, anlamli_kumeler, sebep_eki)

    KARISIK olmasi icin: baskin disinda >=3 okumalik ve >=%5'lik bir kume,
    YA DA baskin kumenin cogunlugu kaybetmis olmasi.
    """
    if not boyut or not n_okuma:
        return u'OLCULEMEDI', [], u''
    baskin = 100.0 * boyut[0] / n_okuma
    anlamli = [b for b in boyut[1:]
               if b >= ANLAMLI_KUME_OKUMA and 100.0 * b / n_okuma >= ANLAMLI_KUME_ORAN]
    temsil_yok = baskin < BASKIN_TEMSIL_TABANI
    if anlamli:
        ek = (u'ANLAMLI ikincil kume(ler): %s okuma (esik: >=%d okuma ve >=%%%s)'
              % (u', '.join(str(b) for b in anlamli), ANLAMLI_KUME_OKUMA,
                 vir(ANLAMLI_KUME_ORAN, 0)))
    elif temsil_yok:
        ek = (u'anlamli ikincil kume YOK ama baskin kume cogunlukta degil '
              u'(%s%% < %s%%): konsensus okumalarin cogunlugunu temsil etmiyor'
              % (vir(baskin, 1), vir(BASKIN_TEMSIL_TABANI, 0)))
    else:
        ek = (u'baskin disindaki butun kumeler gurultu boyutunda '
              u'(en buyugu %d okuma < %d): olculen gurultu tavani 2 okuma'
              % (max(boyut[1:]) if len(boyut) > 1 else 0, ANLAMLI_KUME_OKUMA))
    return (u'KARISIK' if (anlamli or temsil_yok) else u'SAF'), anlamli, ek

# Okuma-okuma kumeleme esigi SABIT DEGILDIR, olculen hata oranindan turetilir.
#
# NEDEN: iki okuma AYNI organizmadan gelse bile, her biri kendi okuma hatasini
# tasir; aralarindaki beklenen fark 2*eps'tir. Olculdu (F2-4_500148): eps =
# %2,208, yani ayni organizmanin iki okumasi ortalama %4,4 farkli gorunuyor.
# Sabit %97 esigiyle bu okumalar AYRI kumelere dusuyordu ve 60 okuma 34 kumeye
# bolunup her kutu "karisik" gorunuyordu - olculen sey biyoloji degil, okuma
# hatasiydi. Esik, beklenen ikili farkin KARISIM_KAT katina kadar tolerans
# taniyacak sekilde hesaplanir.
KARISIM_KAT = 2.5
KARISIM_ESIK_TABAN = 88.0   # bundan gevsek olmaz (ayri cinsler birlesmesin)
KARISIM_ESIK_TAVAN = 99.0   # bundan siki olmaz (hatasiz veride bile pay birak)

# Hizalanan uzunluk tabani. Cok kisa bir referans kaydina %99 benzemek, uzun
# bir kayda %98 benzemekten DAHA ZAYIF kanittir; kisa kayit sorgunun yalnizca
# kucuk bir penceresini gorur. Olculdu: 533 bp'lik bir 28S kaydi %99,44 ile
# 1719 bp'lik %99,01'lik kaydin onune geciyordu.
ASGARI_HIZ_UZ = 400

# Kisa liste ve hizalama butcesi.
KL_UST = 120               # her veritabanindan kac aday tam degerlendirilsin
ON_PENCERE = 900           # on eleme hizalamasinda kullanilan sorgu penceresi
KESIN_UST = 24             # on elemeden sonra TAM hizalanacak aday sayisi
ADAY_HAVUZU_KONTROL = 800  # akis sirasinda tutulan aday havuzu (bellek icin)

# AYIRT EDILEBILIRLIK katsayisi: kac sigma. 3 sigma = %99,7 guven.
AYIRT_SIGMA = 3.0

# M1'in "hizli" kipinde taranan kucuk veritabanlari. Buyukler (SILVA, UNITE,
# PR2, ROD) yalniz "tam" kipinde taranir; sebebi olculen suredir (asagida).
HIZLI_VTB = (u'RefSeq bakteri 16S', u'RefSeq arke 16S', u'RefSeq mantar ITS',
             u'RefSeq mantar 28S', u'RefSeq mantar 18S', u'RefSeq ref_all2')


def _tip_kaydi(etiket, baslik):
    u"""Referans basligindan TIP KAYDI (type material) isaretini ayristir.

    UC DURUM VARDIR ve ucu de birbirinden farklidir:
      EVET      : baslik acikca tip materyali diyor
      HAYIR     : veritabani tip bilgisi TASIYOR ama bu kayitta yok
      BILGI_YOK : veritabani basliklarinda tip bilgisi HIC TASIMIYOR

    "BILGI_YOK"u "HAYIR" saymak yanlis olurdu: SILVA/UNITE/ROD basliklarinda
    tip alani hic yoktur, oradaki bir kaydin tip olmadigini iddia edemeyiz.
    Olculen kapsam (2026-08-09, baslik taramasi):
        fungi.ITS.fna      20 271 / 20 394 kayit "from TYPE material"
        fungi.28SrRNA.fna  12 845 / 12 890
        fungi.18SrRNA.fna   4 009 /  4 037
        ref_all2.fna       37 125 / 65 358
        bacteria.16S.fna        0 / 26 877   <- alan YOK (NR_ kayitlari isaretsiz)
        archaea.16S.fna         0 /  1 160   <- alan YOK
        SILVA SSU/LSU, UNITE, ROD   0        <- alan YOK
    """
    b = baslik or u''
    tasiyor = etiket in (u'RefSeq mantar ITS', u'RefSeq mantar 28S',
                         u'RefSeq mantar 18S', u'RefSeq ref_all2')
    if re.search(r'from\s+TYPE\s+material', b, re.I):
        return u'EVET'
    if re.search(r'\btype\s+strain\b', b, re.I):
        return u'EVET'
    if tasiyor:
        return u'HAYIR'
    return u'BILGI_YOK'


def _kayit_no(baslik):
    u"""Baslikttan kayit numarasini cikar (ilk bosluga/boruya kadar olan belirtec)."""
    b = (baslik or u'').strip()
    if not b:
        return u'-'
    return re.split(r'[\s|]', b, 1)[0][:40] or u'-'


def _okuma_ornekle(fastq_yolu, n=KARISIM_OKUMA, en_az_uz=300):
    u"""FASTQ'tan duzenli arayla n okuma ornekle. Donen: (diziler, ort_hata_orani).

    Duzenli aralikla ornekleme, dosyanin basindaki okumalarin (genelde daha
    kisa/erken) baskin cikmasini engeller. Hata orani Phred'den hesaplanir:
    ort_hata = ortalama(10^(-Q/10)). Ortalama Q'dan hata TUREMEZ - ustel
    ortalama ile aritmetik ortalama ayni sey degildir ve hata orani
    oldugundan kucuk cikardi.
    """
    if not os.path.exists(fastq_yolu):
        return [], None
    diziler = []
    hata_top = 0.0
    hata_n = 0
    havuz = []
    with io.open(fastq_yolu, encoding='utf-8', errors='replace') as fh:
        d = None
        for i, satir in enumerate(fh):
            m = i % 4
            if m == 1:
                d = satir.strip().upper()
            elif m == 3 and d is not None:
                q = satir.strip()
                if len(d) >= en_az_uz:
                    havuz.append((d, q))
                d = None
    if not havuz:
        return [], None
    adim = max(1, len(havuz) // n)
    secilen = havuz[::adim][:n]
    for d, q in secilen:
        diziler.append(d)
        for c in q:
            hata_top += 10.0 ** (-(ord(c) - 33) / 10.0)
            hata_n += 1
    return diziler, (hata_top / hata_n if hata_n else None)


def _karisim_esigi(hata_orani):
    u"""Olculen okuma hata oranindan kumeleme esigini turet.

    Ayni organizmadan gelen iki okumanin beklenen farki 2*eps'tir. Esigi buna
    KARISIM_KAT kadar pay birakarak kurariz. Hata orani olculemediyse temkinli
    davranilir (yuksek hata varsayilir), cunku fazla siki bir esik SAF bir
    kutuyu karisik gosterir - bu yanlis yon daha pahalidir.
    """
    eps = hata_orani if hata_orani is not None else 0.02
    beklenen_fark = 2.0 * eps * 100.0
    esik = 100.0 - KARISIM_KAT * beklenen_fark
    return max(KARISIM_ESIK_TABAN, min(KARISIM_ESIK_TAVAN, esik))


def _kumele(K, diziler, esik, pencere=KARISIM_PENCERE):
    u"""Okumalari ikili kimlige gore kumele. Donen: [[indeks,...], ...] buyukten kucuge.

    Basit tek-baglantili kumeleme: iki okuma birbirine `esik`ten yakinsa ayni
    kumeye girer. Amac tur atamasi degil, KAC AYRI ORGANIZMA VAR sorusudur;
    bunun icin tek baglanti yeterlidir ve ucuzdur.

    Maliyet kontrolu: her okumadan yalnizca ortadaki `pencere` bazlik parca
    kullanilir. Tam boy hizalama 60 okuma icin 1770 cift x ~0,1 sn = 3 dakika
    ederdi; pencere ile ~0,02 sn'ye iner.
    """
    n = len(diziler)
    if n == 0:
        return []
    parca = []
    for d in diziler:
        d = K.temizle(d)
        if len(d) > pencere:
            b = (len(d) - pencere) // 2
            d = d[b:b + pencere]
        parca.append(d)
    ebeveyn = list(range(n))

    def bul(a):
        while ebeveyn[a] != a:
            ebeveyn[a] = ebeveyn[ebeveyn[a]]
            a = ebeveyn[a]
        return a

    for i in range(n):
        for j in range(i + 1, n):
            if bul(i) == bul(j):
                continue                      # zaten ayni kumede, hizalama gereksiz
            k1, _d1 = K.hizala(parca[i], parca[j])
            k2, _d2 = K.hizala(parca[i], K.rc(parca[j]))   # yon farki kume bozmasin
            if max(k1, k2) >= esik:
                ebeveyn[bul(i)] = bul(j)
    kume = collections.defaultdict(list)
    for i in range(n):
        kume[bul(i)].append(i)
    return sorted(kume.values(), key=len, reverse=True)


def _medoid(K, diziler, indeksler, pencere=KARISIM_PENCERE):
    u"""Bir kumenin MEDOIDI: kumedeki digerlerine toplam benzerligi en yuksek okuma.

    Ortalama/konsensus almak yerine gercek bir okuma secilir; boylece sorgu
    dizisi uydurma bir kimera olmaz. Kume tek elemanliysa o eleman doner.
    """
    if not indeksler:
        return None
    if len(indeksler) == 1:
        return diziler[indeksler[0]]
    parca = {}
    for i in indeksler:
        d = K.temizle(diziler[i])
        if len(d) > pencere:
            b = (len(d) - pencere) // 2
            d = d[b:b + pencere]
        parca[i] = d
    en_iyi, en_iyi_p = None, -1.0
    for i in indeksler:
        p = 0.0
        for j in indeksler:
            if i == j:
                continue
            p += K.hizala(parca[i], parca[j])[0]
        if p > en_iyi_p:
            en_iyi, en_iyi_p = i, p
    return diziler[en_iyi]


def _pencere_sec(dizi, pencere=ON_PENCERE):
    u"""Sorgudan N'i EN AZ olan `pencere` bazlik parcayi sec.

    On eleme hizalamasinda N'li bolge kullanmak, aday siralamasini gurultuye
    bogar (N daima uyumsuz sayilir). Kayan pencere ile en temiz bolge secilir.
    """
    if len(dizi) <= pencere:
        return dizi
    adim = max(1, (len(dizi) - pencere) // 40)
    en_iyi, en_az = 0, None
    for b in range(0, len(dizi) - pencere + 1, adim):
        n = dizi.count('N', b, b + pencere)
        if en_az is None or n < en_az:
            en_iyi, en_az = b, n
            if n == 0:
                break
    return dizi[en_iyi:en_iyi + pencere]


def _ayirt_edilebilir(K, en_iyi, ikinci, hata_orani, sorgu_uz):
    u"""EN IYI IKI TURU BIZIM VERI KALITEMIZLE AYIRT EDEBILIR MIYIZ?

    Bu kural KODDADIR, elle karar verilmez. Iki bagimsiz sinav uygulanir ve
    IKISI DE gecmek zorundadir:

    (A) REFERANS AYRIMI - iki referans dizi birbirinden yeterince farkli mi?
        Iki referans arasindaki fark D_ref baz. Bizim okuma hatamiz okuma
        basina eps; hizalanan L baz uzerinde beklenen hatali baz sayisi
        E = eps*L, saçilimi ~sqrt(E). Iki turun ayrilabilmesi icin
            D_ref >= AYIRT_SIGMA * sqrt(E + 1)   ve   D_ref >= 3
        olmali. Degilse iki tur bizim veri kalitemizle AYNI gorunur ve tur
        adi vermek olcumun tasiyabileceginden fazlasini iddia etmektir.

    (B) OLCUM AYRIMI - bizim iki isabetimiz arasindaki fark gurultuden buyuk mu?
        En iyi isabette m1, ikincide m2 uyumsuz baz var. Sayim gurultusu
        ~sqrt(m1+m2). Ayrilabilmesi icin
            |m1 - m2| >= AYIRT_SIGMA * sqrt(m1 + m2 + 1)
        olmali. Degilse iki aday olcum hatasi icinde esittir.

    Donen: dict(ayrilir, sebep, d_ref, esik_a, olcum_farki, esik_b)
    """
    out = dict(ayrilir=None, sebep=u'', d_ref=None, esik_a=None,
               olcum_farki=None, esik_b=None)
    if not en_iyi or not ikinci:
        out['ayrilir'] = None
        out['sebep'] = u'ikinci isabet yok - ayirt edilebilirlik OLCULEMEDI'
        return out
    eps = hata_orani if hata_orani is not None else 0.01   # olculemezse temkinli
    L = float(en_iyi.get('hiz_uz') or sorgu_uz or 1)

    # (A) iki referansi BIRBIRINE hizala
    r1 = en_iyi.get('dizi') or u''
    r2 = ikinci.get('dizi') or u''
    if not r1 or not r2:
        out['sebep'] = u'referans dizileri elde yok - ayirt edilebilirlik OLCULEMEDI'
        return out
    kisa, uzun = (r1, r2) if len(r1) <= len(r2) else (r2, r1)
    ref_kimlik, _u = K.hizala(kisa, uzun)
    d_ref = (100.0 - ref_kimlik) / 100.0 * len(kisa)
    esik_a = max(3.0, AYIRT_SIGMA * math.sqrt(eps * L + 1.0))

    # (B) bizim iki olcumumuz arasindaki fark
    m1 = en_iyi.get('uzaklik')
    m2 = ikinci.get('uzaklik')
    if m1 is None or m2 is None:
        out['sebep'] = u'uyumsuz baz sayilari yok - olcum ayrimi OLCULEMEDI'
        out['d_ref'] = d_ref
        out['esik_a'] = esik_a
        return out
    olcum_farki = abs(m1 - m2)
    esik_b = AYIRT_SIGMA * math.sqrt(m1 + m2 + 1.0)

    a_gecti = d_ref >= esik_a
    b_gecti = olcum_farki >= esik_b
    out.update(d_ref=d_ref, esik_a=esik_a, olcum_farki=olcum_farki, esik_b=esik_b,
               ayrilir=bool(a_gecti and b_gecti))
    if a_gecti and b_gecti:
        out['sebep'] = (u'iki referans %s baz farkli (esik %s) ve olcum farki %s baz '
                        u'(esik %s) - AYRILIR'
                        % (vir(d_ref, 1), vir(esik_a, 1), vir(olcum_farki, 1), vir(esik_b, 1)))
    elif not a_gecti:
        out['sebep'] = (u'iki REFERANS birbirinden yalnizca %s baz farkli, bizim hata '
                        u'orani %s ile gereken esik %s - bu iki tur bizim veri '
                        u'kalitemizle AYRILAMAZ'
                        % (vir(d_ref, 1), vir(100 * eps, 3), vir(esik_a, 1)))
    else:
        out['sebep'] = (u'olcum farki %s baz, gurultu esigi %s baz - iki aday olcum '
                        u'hatasi icinde ESIT' % (vir(olcum_farki, 1), vir(esik_b, 1)))
    return out


def _kraken_haritasi(kraken_kok):
    u"""Barkod raporlarini kutu siniflarina esle. Donen: {(sinif,no): rapor_yolu}.

    ESLEME NASIL DOGRULANDI: her sinif klasorundeki raporlar barkod numarasina
    gore siralanir ve sirasiyla o sinifin 1..4 numarali kutularina baglanir
    (A1 -> 01..04, A2 -> 05..08, F2 -> 09..12, F1 -> 13..16, B -> 17..20).
    Esleme dogrulandi (2026-08-09): 20 kutunun 99 taxid'inin TAMAMI kendi
    atanan barkod raporunda mevcut, eslesmeyen 0. Yanlis esleme olsaydi
    taxid'lerin buyuk kismi raporda bulunamazdi.
    """
    harita = {}
    if not os.path.isdir(kraken_kok):
        return harita
    for sinif in sorted(os.listdir(kraken_kok)):
        kls = os.path.join(kraken_kok, sinif)
        if not os.path.isdir(kls):
            continue
        raporlar = sorted(glob.glob(os.path.join(kls, '*_kraken2.report')))
        for i, r in enumerate(raporlar, 1):
            harita[(sinif, i)] = r
    return harita


def _kraken_oku(yol):
    u"""Kraken2 raporu -> {taxid: dict(ad, duzey, yuzde, atanan, altagac)}.

    Kraken2 raporunun sutunlari: yuzde, altagac_okuma, dogrudan_okuma, duzey,
    taxid, ad. "Guven degeri" Kraken2 raporunda AYRI bir sutun DEGILDIR; en
    yakin karsiligi o takson icin dogrudan atanan okuma yuzdesidir, bu yuzden
    "guven" sutununa bu deger yazilir ve rapor bunu boyle acikca soyler.
    """
    out = {}
    if not os.path.exists(yol):
        return out
    with io.open(yol, encoding='utf-8', errors='replace') as fh:
        for satir in fh:
            p = satir.rstrip(u'\n').split(u'\t')
            if len(p) < 6:
                continue
            try:
                out[p[4].strip()] = dict(
                    yuzde=sayi(p[0]), altagac=int(sayi(p[1]) or 0),
                    atanan=int(sayi(p[2]) or 0), duzey=p[3].strip(),
                    ad=p[5].strip())
            except (ValueError, TypeError):
                continue
    return out


def modul_1_kimlik(kay, rap, kn, kip=u'hizli', yalniz=None, tavan=0):
    u"""MODUL 1 - her kutunun kimligini olc ve Kraken2 ile yan yana koy.

    kip: 'yok'   -> hic kosma (bulgu: ATLANDI)
         'hizli' -> yalniz kucuk RefSeq kumeleri (olculen sure asagida)
         'tam'   -> butun offline veritabanlari (SILVA SSU/LSU, UNITE, PR2, ROD dahil)
    """
    M = u'1 KIMLIK'
    t_basla = time.time()

    if kip == u'yok':
        rap.atla(M, u'M1-KAPALI', u'her kutu icin tur duzeyi kimlik olcumu',
                 u'--m1-kip yok verildi, modul bilerek kapatildi', u'-')
        return

    # --- karar mantigini verification'dan ice aktar (yeniden YAZILMAZ)
    K, hata = modul_yukle(kay.kimlik_dogrulama, 'kimlik_dogrulama')
    if K is None:
        rap.atla(M, u'M1-MOTOR', u'identity_verification.py hizalama motoru yuklenebilmeli',
                 hata, kay.kimlik_dogrulama)
        return
    T, hata2 = modul_yukle(kay.tum_kutu, 'tum_kutu_kimlikleri')
    if T is None:
        rap.atla(M, u'M1-MOTOR2', u'all_bin_identities.py toplu tarama yuklenebilmeli',
                 hata2, kay.tum_kutu)
        return
    try:
        import numpy  # noqa: F401  (hizala() numpy'siz cok yavas)
    except ImportError:
        rap.atla(M, u'M1-NUMPY', u'numpy kurulu olmali (hizalama motoru kullaniyor)',
                 u'numpy ice aktarilamadi', u'-')
        return

    # Bellek icin aday havuzunu kucult. Bu yalniz BU SUREC ICINDEKI kopyayi
    # etkiler; verification/identity_verification.py dosyasi DEGISMEZ.
    # Gerekce: dosyanin kendi notuna gore olculen en kotu on-eleme sirasi 45'tir;
    # 800'luk havuz bu sinirin 17 katidir, kesme baglayici degildir.
    K.ADAY_HAVUZU = ADAY_HAVUZU_KONTROL

    # --- kutu envanteri: KANONIK konsensus indeksi (hepsi SENSE yonde)
    ind = tsv_oku(kay.konsensus_indeks, yorum=None)
    if ind is None:
        rap.atla(M, u'M1-ENVANTER', u'konsensus_kanonik/INDEKS.tsv okunabilmeli',
                 u'dosya yok', kay.konsensus_indeks)
        return
    if not ind:
        rap.ekle(M, u'M1-ENVANTER-BOS', KRITIK,
                 u'kanonik konsensus indeksi kutu icermeli',
                 u'indeks BOS - hicbir kutu okunamadi', kay.konsensus_indeks,
                 u'build_canonical.py ile indeksi yeniden uretin')
        return

    kutular = []
    for r in ind:
        yol = os.path.join(kay.konsensus_kok, r.get('dosya', u''))
        if not os.path.exists(yol):
            rap.ekle(M, u'M1-KONSENSUS-YOK', CIDDI,
                     u'indekste yazan konsensus dosyasi diskte bulunmali',
                     u'%s -> dosya yok' % r.get('dosya'), kay.konsensus_indeks,
                     u'Kutu olculemez; indeks ile klasor ayrismis.')
            continue
        try:
            with io.open(yol, encoding='utf-8', errors='replace') as fh:
                diz = u''.join(s.strip() for s in fh if not s.startswith(u'>')).upper()
        except IOError as e:
            rap.ekle(M, u'M1-KONSENSUS-OKUNAMADI', CIDDI,
                     u'konsensus dosyasi okunabilmeli', u'%s: %s' % (r.get('dosya'), e),
                     yol)
            continue
        kutular.append(dict(kutu=r.get('kutu'), sinif=r.get('sinif'),
                            yol=yol, dizi=K.temizle(diz),
                            yon=r.get('eski_yon'), cevrildi=r.get('cevrildi')))

    if yalniz:
        istenen = set(x.strip() for x in yalniz.split(u',') if x.strip())
        kutular = [k for k in kutular if k['kutu'] in istenen]
    if tavan:
        kutular = kutular[:tavan]
    if not kutular:
        rap.atla(M, u'M1-KUTU-YOK', u'en az bir kutu olculmeli',
                 u'suzgecten sonra kutu kalmadi', kay.konsensus_indeks)
        return

    yaz(u'M1: %d bins, mode=%s' % (len(kutular), kip))

    # -----------------------------------------------------------------
    # ASAMA A - N orani, okuma hata orani, karisim/saflik.
    # Bu asama HIC VERITABANI OKUMAZ; tamamen kutunun kendi verisiyle olcer.
    # -----------------------------------------------------------------
    yaz(u'M1/A: N fraction, read error rate and mixture measurement...')
    # Kac kutunun GERCEKTEN hesaplandigini sayariz. Kontrol noktasindan okunan
    # kutular sifir saniye surer; onlari da sayarsak "kutu basina 0 saniye"
    # cikar ve tam tarama tahmini gercek disi bir sekilde kisa gorunur.
    taze_a = [0]
    can = Canlilik(u'M1/A kutu', len(kutular))
    for i, kb in enumerate(kutular, 1):
        can.vur(i, kb['kutu'])
        diz = kb['dizi']
        n_say = diz.count('N')
        kb['n_oran'] = 100.0 * n_say / max(1, len(diz))
        kb['kons_uz'] = len(diz)

        # kutuya ait FASTQ: <sinif-no>/<sinif-no>-reads_<taxid>.fastq
        grup = kb['kutu'].rsplit(u'_', 1)[0]
        taxid = kb['kutu'].rsplit(u'_', 1)[-1]
        kb['grup'] = grup
        kb['taxid'] = taxid
        fq = os.path.join(kay.fastq, grup, u'%s-reads_%s.fastq' % (grup, taxid))
        if not os.path.exists(fq):
            adaylar = glob.glob(os.path.join(kay.fastq, grup, u'*%s*.fastq*' % taxid))
            fq = adaylar[0] if adaylar else None
        kb['fastq'] = fq

        anahtar = md5_metin(u'M1A', VERSIYON, kb['kutu'], dosya_imzasi(kb['yol']),
                            dosya_imzasi(fq or u'-'), KARISIM_OKUMA, KARISIM_KAT,
                            KARISIM_ESIK_TABAN, KARISIM_ESIK_TAVAN, SAF_ESIGI)
        onbellek = kn.oku(anahtar, [kb['yol']] + ([fq] if fq else []))
        if onbellek:
            kb.update(onbellek)
            continue
        taze_a[0] += 1          # bu kutu GERCEKTEN hesaplandi (sure tahmini icin)

        if not fq:
            kb['hata_orani'] = None
            kb['okuma_sayisi'] = None
            kb['baskin_oran'] = None
            kb['saflik'] = u'OLCULEMEDI'
            kb['saflik_sebep'] = u'kutuya ait FASTQ bulunamadi'
            kb['kume_boyutlari'] = []
            kb['okuma_sorgu'] = None
            rap.ekle(M, u'M1-FASTQ-YOK', UYARI,
                     u'her kutunun ham okuma dosyasi bulunmali',
                     u'%s icin FASTQ yok - karisim ve hata orani OLCULEMEDI' % kb['kutu'],
                     os.path.join(kay.fastq, grup))
        else:
            okumalar, eps = _okuma_ornekle(fq)
            kb['hata_orani'] = eps
            kb['okuma_sayisi'] = len(okumalar)
            if len(okumalar) < 5:
                kb['baskin_oran'] = None
                kb['saflik'] = u'OLCULEMEDI'
                kb['saflik_sebep'] = u'ornege yeterli uzunlukta okuma yok (%d)' % len(okumalar)
                kb['kume_boyutlari'] = []
                kb['okuma_sorgu'] = None
            else:
                k_esik = _karisim_esigi(eps)
                kumeler = _kumele(K, okumalar, k_esik)
                boyut = [len(c) for c in kumeler]
                baskin = 100.0 * boyut[0] / len(okumalar)
                kb['kume_boyutlari'] = boyut
                kb['baskin_oran'] = baskin
                kb['karisim_esigi'] = k_esik
                # 2026-08-09: hukum artik yalniz baskin orandan degil, ANLAMLI
                # ikincil kume varliginndan veriliyor. Gerekce saflik_hukmu()
                # ustundeki blokta, olcumle birlikte yazili.
                kb['saflik'], kb['anlamli_kumeler'], _ek = saflik_hukmu(
                    boyut, len(okumalar))
                kb['saflik_sebep'] = (
                    u'ornekteki %d okuma %d kumeye ayrildi (%s); baskin kume %s%%. '
                    u'Kumeleme esigi olculen okuma hata oranindan turetildi: '
                    u'eps=%s%% -> esik %s%%. HUKUM: %s'
                    % (len(okumalar), len(boyut),
                       u'/'.join(str(b) for b in boyut[:8]) +
                       (u'/...' if len(boyut) > 8 else u''),
                       vir(baskin, 1), vir(100 * eps, 3) if eps is not None else u'?',
                       vir(k_esik, 1), _ek))
                kb['okuma_sorgu'] = _medoid(K, okumalar, kumeler[0])
        kn.yazdir(anahtar, dict((a, kb.get(a)) for a in (
            'n_oran', 'kons_uz', 'hata_orani', 'okuma_sayisi', 'baskin_oran',
            'saflik', 'saflik_sebep', 'kume_boyutlari', 'okuma_sorgu',
            'karisim_esigi')))
    a_sure = can.bitti(len(kutular))
    # Kutu basina hazirlik suresi YALNIZ taze hesaplanan kutulardan olculur.
    if taze_a[0] >= 2:
        rap.olcum[u'_m1a_kutu_sn'] = a_sure / float(taze_a[0])
        rap.olcum[u'M1/A hazirlik'] = u'%d kutu (%d taze, %d kontrol noktasindan), ' \
                                      u'%s -> %s sn/kutu' % (
            len(kutular), taze_a[0], len(kutular) - taze_a[0], sure_metni(a_sure),
            vir(a_sure / float(taze_a[0]), 1))
    else:
        rap.olcum[u'M1/A hazirlik'] = (
            u'%d kutunun %d tanesi kontrol noktasindan okundu; kutu basina '
            u'hazirlik suresi BU KOSUDA OLCULMEDI'
            % (len(kutular), len(kutular) - taze_a[0]))

    # --- SORGU SECIMI: N yuksekse konsensus KULLANILMAZ, okumaya donulur
    sorgular = {}
    for kb in kutular:
        if kb['n_oran'] > N_ESIGI and kb.get('okuma_sorgu'):
            kb['sorgu_kaynagi'] = u'HAM OKUMA (baskin kume medoidi)'
            kb['sorgu'] = K.temizle(kb['okuma_sorgu'])
            rap.ekle(M, u'M1-N-YUKSEK', UYARI,
                     u'konsensusun N orani %%%s altinda olmali' % vir(N_ESIGI, 0),
                     u'%s: N orani %%%s - konsensus kullanilmadi, ham okumaya donuldu '
                     u'(saflik: %s)' % (kb['kutu'], vir(kb['n_oran'], 1), kb.get('saflik')),
                     kb['yol'],
                     u'Bu kutu konsensus uzerinden "adlandirilamaz" gorunur; '
                     u'okuma bazli olcum esas alinmalidir.')
        elif kb['n_oran'] > N_ESIGI:
            kb['sorgu_kaynagi'] = u'KONSENSUS (N yuksek ama okuma yok)'
            kb['sorgu'] = kb['dizi']
            rap.ekle(M, u'M1-N-YUKSEK-OKUMASIZ', CIDDI,
                     u'N orani yuksek kutuda ham okumaya donulebilmeli',
                     u'%s: N orani %%%s ama okuma sorgusu uretilemedi'
                     % (kb['kutu'], vir(kb['n_oran'], 1)), kb['yol'])
        else:
            kb['sorgu_kaynagi'] = u'KONSENSUS'
            kb['sorgu'] = kb['dizi']
        sorgular[kb['kutu']] = kb['sorgu']

    # -----------------------------------------------------------------
    # ASAMA B - veritabani taramasi. Her veritabani TEK akista taranir ve
    # butun kutular ayni akistan beslenir (kutu basina ayri gecis 100x12=1200
    # tam dosya gecisi ederdi).
    # -----------------------------------------------------------------
    vtb = [(e, d, t) for e, d, t, kullan, _n in K.VTB if kullan]
    if kip == u'hizli':
        vtb = [v for v in vtb if v[0] in HIZLI_VTB]
    yaz(u'M1/B: %d databases to scan: %s' % (len(vtb), u', '.join(v[0] for v in vtb)))

    bulgular_kutu = collections.defaultdict(list)   # kutu -> [isabet, ...]
    for ei, (etiket, dosya, lokus) in enumerate(vtb, 1):
        yol = os.path.join(kay.refdb, dosya)
        if not os.path.exists(yol):
            rap.ekle(M, u'M1-VTB-YOK', CIDDI,
                     u'VTB listesindeki veritabani dosyasi diskte bulunmali',
                     u'%s (%s) yok' % (etiket, dosya), kay.refdb,
                     u'Bu veritabani hicbir kutu icin oy kullanamadi.')
            continue
        anahtar = md5_metin(u'M1B', VERSIYON, etiket, dosya_imzasi(yol), KL_UST,
                            ON_PENCERE, KESIN_UST, ADAY_HAVUZU_KONTROL,
                            md5_metin(*[u'%s=%s' % (k, sorgular[k]) for k in sorted(sorgular)]))
        onbellek = kn.oku(anahtar, [yol])
        if onbellek is not None:
            for kutu, isabetler in onbellek.items():
                bulgular_kutu[kutu].extend(isabetler)
            yaz(u'  [%d/%d] %s - read from checkpoint' % (ei, len(vtb), etiket))
            continue

        t0 = time.time()
        mb = os.path.getsize(yol) / 1e6
        yaz(u'  [%d/%d] %s (%s MB) scanning...' % (ei, len(vtb), etiket, vir(mb, 1)))
        can2 = Canlilik(u'%s kayit' % etiket, None, aralik=20.0)
        try:
            kls, taranan = T.toplu_kisa_liste(K, yol, sorgular, KL_UST,
                                              ilerle=lambda n: can2.vur(n))
        except MemoryError:
            rap.atla(M, u'M1-VTB-BELLEK', u'%s taranabilmeli' % etiket,
                     u'bellek yetmedi (MemoryError) - kutu sayisini --m1-tavan ile '
                     u'dusurun ya da veritabanini ayri kosun', yol)
            continue
        except Exception as e:
            # HATA YUTULMAZ: hangi veritabani neden taranamadi raporlanir.
            rap.atla(M, u'M1-VTB-HATA', u'%s taranabilmeli' % etiket,
                     u'%s: %s' % (type(e).__name__, e), yol)
            continue
        tarama_sn = time.time() - t0
        rap.olcum[u'M1 tarama %s' % etiket] = (
            u'%s MB, %d kayit, %s  (%s MB/sn)'
            % (vir(mb, 1), taranan, sure_metni(tarama_sn), vir(mb / max(0.001, tarama_sn), 2)))

        # --- beklenen kayit sayisi ile karsilastir (budanmis tarama yakalanir)
        bekl = getattr(T, 'BEKLENEN_KAYIT', {}).get(etiket)
        if bekl and taranan < bekl:
            rap.ekle(M, u'M1-KAPSAM-EKSIK', CIDDI,
                     u'%s taramasi %d kaydin tamamini gormeli' % (etiket, bekl),
                     u'yalnizca %d kayit tarandi (%%%s)'
                     % (taranan, vir(100.0 * taranan / bekl, 1)), yol,
                     u'Budanmis tarama, gercek en iyi isabeti kacirabilir.')

        # --- her kutu icin iki asamali hizalama
        t1 = time.time()
        can3 = Canlilik(u'%s hizalama' % etiket, len(kutular))
        vtb_sonuc = {}
        for ki, kb in enumerate(kutular, 1):
            can3.vur(ki, kb['kutu'])
            kl = kls.get(kb['kutu']) or []
            if not kl:
                continue
            q = kb['sorgu']
            qp = _pencere_sec(q, ON_PENCERE)
            qp_rc = K.rc(qp)
            # ON ELEME: kisa pencereyle butun adaylari sirala (ucuz).
            #
            # IKI SEY BURADA KRITIK:
            # (1) hizala() INFIX'tir: KISA diziyi UZUN olanin icine yerlestirir.
            #     Sorgu hedeften uzunsa hizalama anlamsizlasir ve kimlik
            #     yapay olarak uzunluk oranina duser (olculdu: 3707 bp sorgu
            #     1700 bp'lik bir 18S kaydina karsi %47 veriyordu - bu bir
            #     benzerlik degil, uzunluk farkinin ta kendisiydi). Bu yuzden
            #     her cagride kisa olan sorgu, uzun olan hedef yapilir -
            #     kimlik_dogrulama.kl_degerlendir() de aynen boyle yapar.
            # (2) Referans kayitlarin YONU kume kume degisir. Tek yonde arayan
            #     bir olcum kayitlarin bir kismini kacirir; bu yuzden pencere
            #     asamasinda her iki yon denenir ve KAZANAN YON tam hizalamada
            #     kullanilir (tam boyu iki kez hizalamak iki kat pahali olurdu).
            on = []
            for a in kl:
                t = a['dizi']
                d, u_ = (qp, t) if len(qp) <= len(t) else (t, qp)
                ileri, _x = K.hizala(d, u_)
                d2, u2 = (qp_rc, t) if len(qp_rc) <= len(t) else (t, qp_rc)
                geri, _y = K.hizala(d2, u2)
                on.append((max(ileri, geri), a, u'+' if ileri >= geri else u'-'))
            on.sort(key=lambda x: -x[0])
            # TIP KAYITLARI kesinlikle degerlendirilsin: tip kaydina benzerlik
            # cok daha guclu kanittir, on elemede kaybolmasi kabul edilemez.
            secili = [(a, yon) for _k, a, yon in on[:KESIN_UST]]
            secili_bas = set(id(a) for a, _y in secili)
            for _k, a, yon in on[KESIN_UST:]:
                if _tip_kaydi(etiket, a['baslik']) == u'EVET' and id(a) not in secili_bas:
                    secili.append((a, yon))
                    secili_bas.add(id(a))
                    if len(secili) >= KESIN_UST + 8:
                        break
            # KESIN OLCUM: tam boy hizalama, kazanan yonde, KISA icine UZUN.
            isabetler = []
            q_rc = K.rc(q)
            for a, yon in secili:
                qq = q if yon == u'+' else q_rc
                t = a['dizi']
                d, u_ = (qq, t) if len(qq) <= len(t) else (t, qq)
                kim, uz = K.hizala(d, u_)
                cins, tur, tam = K.ad_coz(a['baslik'])
                isabetler.append(dict(
                    kimlik=kim, uzaklik=int(uz), baslik=a['baslik'][:200],
                    kayit=_kayit_no(a['baslik']), cins=cins, tur=tur,
                    vtb=etiket, lokus=lokus, sira=a.get('sira'),
                    tip=_tip_kaydi(etiket, a['baslik']),
                    hiz_uz=min(len(q), len(a['dizi'])),
                    dizi=a['dizi'][:6000]))
            isabetler.sort(key=lambda x: -x['kimlik'])
            vtb_sonuc[kb['kutu']] = isabetler[:6]
            bulgular_kutu[kb['kutu']].extend(isabetler[:6])
        hiz_sn = time.time() - t1
        rap.olcum[u'M1 hizalama %s' % etiket] = (
            u'%d kutu, %s (%s sn/kutu)'
            % (len(kutular), sure_metni(hiz_sn), vir(hiz_sn / max(1, len(kutular)), 2)))
        kn.yazdir(anahtar, vtb_sonuc)
        yaz(u'  [%d/%d] %s done: scan %s + alignment %s'
            % (ei, len(vtb), etiket, sure_metni(tarama_sn), sure_metni(hiz_sn)))

    # -----------------------------------------------------------------
    # ASAMA C - hukum + Kraken karsilastirmasi
    # -----------------------------------------------------------------
    yaz(u'M1/C: verdict and Kraken comparison...')
    kharita = _kraken_haritasi(kay.kraken)
    if not kharita:
        rap.atla(M, u'M1-KRAKEN-YOK', u'Kraken2 raporlari okunabilmeli',
                 u'"kraken results" altinda rapor bulunamadi', kay.kraken)
    kraken_onbellek = {}

    satirlar = []
    olculemeyen = []
    for kb in kutular:
        isabetler = sorted(bulgular_kutu.get(kb['kutu']) or [],
                           key=lambda x: -x['kimlik'])
        # --- Kraken tarafi
        grup = kb['grup']
        sinif = grup.split(u'-')[0]
        try:
            no = int(grup.split(u'-')[1])
        except (IndexError, ValueError):
            no = None
        krapor = kharita.get((sinif, no))
        k_etiket = k_guven = k_duzey = u'-'
        if krapor:
            if krapor not in kraken_onbellek:
                kraken_onbellek[krapor] = _kraken_oku(krapor)
            kt = kraken_onbellek[krapor].get(kb['taxid'])
            if kt:
                k_etiket = kt['ad']
                k_duzey = kt['duzey']
                # Kraken2 raporunda ayri bir "guven" sutunu YOKTUR; en yakin
                # karsilik, bu taksona DOGRUDAN atanan okumalarin yuzdesidir.
                k_guven = u'%s%% (dogrudan atanan okuma; Kraken2 ayri guven ' \
                          u'sutunu vermez)' % vir(kt['yuzde'], 2)
            else:
                k_etiket = u'RAPORDA YOK'
                rap.ekle(M, u'M1-KRAKEN-TAXID-YOK', UYARI,
                         u'kutunun taxid\'i kendi barkod raporunda bulunmali',
                         u'%s: taxid %s, %s raporunda yok'
                         % (kb['kutu'], kb['taxid'], os.path.basename(krapor)), krapor)
        else:
            k_etiket = u'RAPOR ESLESMEDI'

        if not isabetler:
            olculemeyen.append(kb['kutu'])
            rap.ekle(M, u'M1-OLCULEMEDI', CIDDI,
                     u'her kutu icin en az bir veritabani isabeti olmali',
                     u'%s: hicbir veritabaninda isabet yok - OLCULEMEDI' % kb['kutu'],
                     kb['yol'],
                     u'Sessizce atlanmadi; bu kutu icin kimlik iddiasi YOKTUR.')
            satirlar.append(dict(
                kutu=kb['kutu'], kraken_etiket=k_etiket, kraken_guven=k_guven,
                olculen=u'OLCULEMEDI', duzey=u'OLCULEMEDI', kimlik=None,
                vtb=u'-', kayit=u'-', tip=u'-', ikinci=u'-', fark=None,
                ayrilir=u'-', ayrim_sebep=u'isabet yok', n_oran=kb['n_oran'],
                saflik=kb.get('saflik'), saflik_sebep=kb.get('saflik_sebep'),
                sorgu_kaynagi=kb['sorgu_kaynagi'], uyusan=0, uyusmayan=u'-',
                ad_degisti=u'OLCULEMEDI'))
            continue

        # --- EN IYI ISABET SECIMI: kimlik TEK BASINA yetmez.
        #
        # Cok kisa bir referans kaydina %99,4 benzemek, uzun bir kayda %99,0
        # benzemekten daha zayif kanittir: kisa kayit sorgunun yalnizca kucuk
        # bir penceresini gorur ve o pencere korunmus bir bolgeyse hemen her
        # akraba %99 verir. Olculdu (F2-1_500148): 533 bp'lik bir 28S kaydi
        # %99,44 ile one geciyor, 1719 bp'lik %99,01'lik kaydi eliyordu.
        #
        # KURAL: en yuksek kimligin GURULTU BANDI icinde kalan isabetler
        # istatistiksel olarak esittir; bunlarin arasindan EN UZUN hizalamayi
        # tasiyan secilir. Band, sayim gurultusunden turetilir (3 sigma).
        # Secim IKI ASAMALIDIR ve her iki asamada AYNI kural kullanilir:
        #   1) her veritabani kendi kazananini secer,
        #   2) genel kazanan, veritabani kazananlari arasindan secilir.
        # Iki asamada farkli kural kullanmak, kazanan cinsin kendi
        # veritabanindan oy alamamasina yol aciyordu (olculdu: genel secim
        # Pseudallescheria derken ayni veritabaninin oyu Lomentospora'ya
        # gidiyor ve kazanan "0 veritabani destekli" gorunuyordu).
        def _band_sec(havuz):
            t = max(havuz, key=lambda x: x['kimlik'])
            b = AYIRT_SIGMA * math.sqrt((t.get('uzaklik') or 0) + 1.0) / max(
                1.0, float(t.get('hiz_uz') or 1)) * 100.0
            esit = [h for h in havuz if (t['kimlik'] - h['kimlik']) <= b]
            return max(esit, key=lambda h: (h.get('hiz_uz') or 0, h['kimlik'])), t

        vtb_grup = collections.defaultdict(list)
        for h in isabetler:
            vtb_grup[h['vtb']].append(h)
        vtb_en_iyi = {}
        for v, hh in vtb_grup.items():
            vtb_en_iyi[v] = _band_sec(hh)[0]
        en_iyi, tepe = _band_sec(list(vtb_en_iyi.values()))
        if en_iyi is not tepe:
            rap.ekle(M, u'M1-KISA-HIZALAMA', BILGI,
                     u'en iyi isabet, kimlik ve hizalanan uzunluk birlikte '
                     u'degerlendirilerek secilmeli',
                     u'%s: en yuksek kimlik %s%% ama yalnizca %s bp uzerinden '
                     u'(%s); gurultu bandi icinde kalan %s bp\'lik %s%% '
                     u'isabet secildi (%s)'
                     % (kb['kutu'], vir(tepe['kimlik'], 2), tepe.get('hiz_uz'),
                        tepe['kayit'], en_iyi.get('hiz_uz'),
                        vir(en_iyi['kimlik'], 2), en_iyi['kayit']), kb['yol'])
        if (en_iyi.get('hiz_uz') or 0) < ASGARI_HIZ_UZ:
            rap.ekle(M, u'M1-HIZALAMA-KISA', CIDDI,
                     u'kimlik hukmu en az %d bp uzerinden verilmeli' % ASGARI_HIZ_UZ,
                     u'%s: en iyi isabet yalnizca %s bp uzerinden olculdu (%s%%, %s)'
                     % (kb['kutu'], en_iyi.get('hiz_uz'), vir(en_iyi['kimlik'], 2),
                        en_iyi['kayit']), kb['yol'],
                     u'Kisa pencerede korunmus bolge her akrabaya yuksek kimlik '
                     u'verir; bu bir tur atamasi icin yeterli kanit degildir.')
        # Hukum siralamasi da secilen isabeti basa alsin.
        isabetler = [en_iyi] + [h for h in isabetler if h is not en_iyi]
        # IKINCI EN IYI: en iyiden FARKLI bir tur (ayni turun ikinci kaydi
        # ayrim sorusuna cevap vermez, tautolojik olurdu)
        # Ikinci isabet de EN AZ ASGARI_HIZ_UZ uzerinden olculmus olmali; kisa
        # bir kayit "rakip" diye gosterilirse ayrim sorusu yanlis kurulur.
        def _ikinci_sec(uzunluk_sarti):
            for h in isabetler[1:]:
                if uzunluk_sarti and (h.get('hiz_uz') or 0) < ASGARI_HIZ_UZ:
                    continue
                if (h.get('tur') or h.get('cins')) != (en_iyi.get('tur') or
                                                       en_iyi.get('cins')):
                    return h
            return None
        ikinci = _ikinci_sec(True) or _ikinci_sec(False)
        fark = (en_iyi['kimlik'] - ikinci['kimlik']) if ikinci else None
        if fark is not None and fark < 0:
            # Negatif fark: rakibin ham kimligi daha yuksek ama DAHA KISA bir
            # hizalama uzerinden. Bu bir celiski degil, uzunluk-kimlik takasidir;
            # okuyanin bunu tahmin etmesi beklenmemeli, acikca yazilir.
            rap.ekle(M, u'M1-RAKIP-DAHA-YUKSEK', UYARI,
                     u'secilen en iyi isabetin kimligi rakibinkinden dusuk olmamali',
                     u'%s: secilen isabet %s%% (%s bp, %s), rakip %s%% (%s bp, %s) - '
                     u'rakip daha yuksek kimlikli ama daha KISA hizalama uzerinden'
                     % (kb['kutu'], vir(en_iyi['kimlik'], 2), en_iyi.get('hiz_uz'),
                        en_iyi['kayit'], vir(ikinci['kimlik'], 2),
                        ikinci.get('hiz_uz'), ikinci['kayit']), kb['yol'],
                     u'Tur atamasi bu iki kayit arasinda karara baglanmamalidir.')

        # --- savunulabilir duzey: identity_verification.py'nin kurali (YENIDEN YAZILMAZ)
        sav = K.savunulabilir_duzey(isabetler, en_iyi.get('lokus') or 'SSU')
        duzey = sav['duzey']
        ad = sav['onerilen_ad']
        gerekce = sav['gerekce']

        # --- AYIRT EDILEBILIRLIK kurali (kodda, elle degil)
        ay = _ayirt_edilebilir(K, en_iyi, ikinci, kb.get('hata_orani'), len(kb['sorgu']))
        if ay['ayrilir'] is False and duzey == u'TUR':
            cins = en_iyi.get('cins')
            tur = en_iyi.get('tur') or u''
            if cins and ikinci and ikinci.get('cins') == cins:
                duzey = u'CINS (tur ayrilamiyor)'
                ad = u'%s cf. %s' % (cins, tur.split()[-1]) if tur else u'%s sp.' % cins
            elif cins:
                duzey = u'CINS (tur ayrilamiyor)'
                ad = u'%s sp.' % cins
            else:
                duzey = u'AILE ve USTU (ad VERILEMEZ)'
                ad = u'adlandirilamayan soy'
            gerekce = u'%s | AYIRT EDILEBILIRLIK: %s' % (gerekce, ay['sebep'])
        else:
            gerekce = u'%s | AYIRT EDILEBILIRLIK: %s' % (gerekce, ay['sebep'])

        # --- kac bagimsiz veritabani ayni CINSI en iyi isabet olarak verdi
        # vtb_en_iyi yukarida, en_iyi ile AYNI kuralla hesaplandi.
        cins_oy = collections.Counter(
            (h.get('cins') or u'?') for h in vtb_en_iyi.values())
        kazanan_cins = en_iyi.get('cins') or u'?'
        uyusan = cins_oy.get(kazanan_cins, 0)
        uyusmayan = u'; '.join(u'%s: %s (%%%s)' % (v, h.get('cins') or u'?', vir(h['kimlik'], 1))
                               for v, h in sorted(vtb_en_iyi.items())
                               if (h.get('cins') or u'?') != kazanan_cins) or u'-'
        if uyusan < 2:
            rap.ekle(M, u'M1-TEK-KAYNAK', UYARI,
                     u'bir kimlik iddiasi en az IKI bagimsiz veritabaninda dogrulanmali',
                     u'%s: "%s" yalnizca %d veritabanindan destekli (uyusmayanlar: %s)'
                     % (kb['kutu'], ad, uyusan, uyusmayan), kb['yol'])

        # --- ad degisti mi (raporda gosterilecek asil sutun)
        k_ad_kok = re.sub(r'[^a-z ]', u'', (k_etiket or u'').lower()).strip()
        o_ad_kok = re.sub(r'[^a-z ]', u'', (ad or u'').lower()).strip()
        if k_etiket in (u'-', u'RAPORDA YOK', u'RAPOR ESLESMEDI'):
            ad_degisti = u'KARSILASTIRILAMADI'
        elif not o_ad_kok:
            ad_degisti = u'KARSILASTIRILAMADI'
        elif k_ad_kok.split()[:1] == o_ad_kok.split()[:1]:
            ad_degisti = u'hayir' if k_ad_kok == o_ad_kok else u'kismen (cins ayni, tur farkli)'
        else:
            ad_degisti = u'EVET'
        if ad_degisti == u'EVET':
            rap.ekle(M, u'M1-AD-DEGISTI', CIDDI,
                     u'olculen kimlik ile Kraken2 etiketi ayni cinse isaret etmeli',
                     u'%s: Kraken "%s" -> olculen "%s" (%%%s, %s%s)'
                     % (kb['kutu'], k_etiket, ad, vir(en_iyi['kimlik'], 2),
                        en_iyi['vtb'],
                        u', TIP KAYDI' if en_iyi['tip'] == u'EVET' else u''),
                     krapor or kb['yol'],
                     u'Rapora gosterilecek asil satir budur: Kraken etiketi olcumle degisti.')

        satirlar.append(dict(
            kutu=kb['kutu'], kraken_etiket=k_etiket, kraken_guven=k_guven,
            kraken_duzey=k_duzey, olculen=ad, duzey=duzey,
            kimlik=en_iyi['kimlik'], vtb=en_iyi['vtb'], kayit=en_iyi['kayit'],
            hiz_uz=en_iyi['hiz_uz'], tip=en_iyi['tip'],
            ikinci=(u'%s (%s, %%%s%s)' % (ikinci.get('tur') or ikinci.get('cins') or ikinci['kayit'],
                                          ikinci['vtb'], vir(ikinci['kimlik'], 2),
                                          u', TIP' if ikinci['tip'] == u'EVET' else u'')
                    ) if ikinci else u'yok',
            fark=fark, ayrilir=(u'EVET' if ay['ayrilir'] else
                                (u'HAYIR' if ay['ayrilir'] is False else u'OLCULEMEDI')),
            ayrim_sebep=ay['sebep'], n_oran=kb['n_oran'],
            saflik=kb.get('saflik'), saflik_sebep=kb.get('saflik_sebep'),
            sorgu_kaynagi=kb['sorgu_kaynagi'], uyusan=uyusan, uyusmayan=uyusmayan,
            ad_degisti=ad_degisti, gerekce=gerekce,
            hata_orani=kb.get('hata_orani')))

        if kb.get('saflik') == u'KARISIK':
            rap.ekle(M, u'M1-KARISIK-KUTU', CIDDI,
                     u'bir kutu tek bir organizma icermeli',
                     u'%s: %s' % (kb['kutu'], kb.get('saflik_sebep')),
                     kb.get('fastq') or kb['yol'],
                     u'Karisik kutudan turetilen konsensus ve ondan turetilen '
                     u'butun ayrim hesaplari supheli.')

    rap.tablolar[u'M1 kimlik tablosu'] = satirlar
    rap.olcum[u'M1 kutu sayisi'] = u'%d olculdu, %d OLCULEMEDI' % (
        len(satirlar) - len(olculemeyen), len(olculemeyen))
    rap.olcum[u'M1 toplam sure'] = sure_metni(time.time() - t_basla)

    # --- TAM TARAMA SURESI TAHMINI, BU KOSUDA OLCULEN HIZLARDAN
    # Tahmin, tahmin degil OLCUME dayanir: bu kosuda gercekten olculen
    # (a) kutu basina on hazirlik suresi, (b) MB basina tarama suresi,
    # (c) kutu-veritabani basina hizalama suresi ile olceklenir. Olculmeyen
    # bir bilesen varsa tahmin URETILMEZ, "olculmedi" yazilir.
    try:
        tum_vtb = [(e, d) for e, d, _t, kullan, _n in K.VTB if kullan]
        toplam_mb = sum(os.path.getsize(os.path.join(kay.refdb, d)) / 1e6
                        for _e, d in tum_vtb
                        if os.path.exists(os.path.join(kay.refdb, d)))
        n_kutu_tam = len(ind)
        # olculen bilesenler
        olc_hazirlik = rap.olcum.get(u'_m1a_kutu_sn')
        tarama_mbsn = []
        hiz_kutu_vtb = []
        for k, v in rap.olcum.items():
            m = re.search(r'\(([\d,]+) MB/sn\)', unicode_(v))
            if k.startswith(u'M1 tarama') and m:
                tarama_mbsn.append(sayi(m.group(1)))
            m2 = re.search(r'\(([\d,]+) sn/kutu\)', unicode_(v))
            if k.startswith(u'M1 hizalama') and m2:
                hiz_kutu_vtb.append(sayi(m2.group(1)))
        if olc_hazirlik and tarama_mbsn and hiz_kutu_vtb:
            orta_mbsn = sum(tarama_mbsn) / len(tarama_mbsn)
            orta_hiz = sum(hiz_kutu_vtb) / len(hiz_kutu_vtb)
            t_hazirlik = olc_hazirlik * n_kutu_tam
            t_tarama = toplam_mb / max(0.01, orta_mbsn)
            t_hizalama = orta_hiz * n_kutu_tam * len(tum_vtb)
            rap.olcum[u'M1 TAM TARAMA TAHMINI'] = (
                u'~%s  =  hazirlik %s (%s sn/kutu x %d kutu) + tarama %s '
                u'(%s MB / %s MB/sn) + hizalama %s (%s sn/kutu-vtb x %d kutu x '
                u'%d vtb). Olceklendirme: hazirlik kutu sayisiyla, tarama toplam '
                u'bayt ile, hizalama kutu x veritabani ile dogru orantili.'
                % (sure_metni(t_hazirlik + t_tarama + t_hizalama),
                   sure_metni(t_hazirlik), vir(olc_hazirlik, 1), n_kutu_tam,
                   sure_metni(t_tarama), vir(toplam_mb, 0), vir(orta_mbsn, 2),
                   sure_metni(t_hizalama), vir(orta_hiz, 2), n_kutu_tam,
                   len(tum_vtb))
                + (u'  DIKKAT: tarama hizi bu kosuda %d sorguyla olculdu; tam '
                   u'kosuda %d sorgu olacagi icin gercek tarama bundan YAVAS '
                   u'olur (tarama maliyeti kayit basina sabit, sorgu basina '
                   u'kucuk bir ek yuk tasir).' % (len(kutular), n_kutu_tam)
                   if len(kutular) < n_kutu_tam else u''))
        else:
            rap.olcum[u'M1 TAM TARAMA TAHMINI'] = (
                u'olculmedi - tahmin icin gereken bilesenlerden biri bu kosuda '
                u'olculmedi (hazirlik/tarama/hizalama)')
    except (OSError, TypeError, ZeroDivisionError) as e:
        rap.olcum[u'M1 TAM TARAMA TAHMINI'] = u'olculmedi (%s)' % e
    if kip == u'hizli':
        rap.ekle(M, u'M1-KISMI-KAPSAM', BILGI,
                 u'kimlik hukmu butun offline veritabanlariyla verilmeli',
                 u'"hizli" kipte yalniz %d kucuk RefSeq kumesi tarandi; SILVA SSU/LSU, '
                 u'UNITE, PR2, ROD taranmadi' % len(HIZLI_VTB), kay.refdb,
                 u'Nihai hukum icin --m1-kip tam ile kosun.')
    yaz(u'M1 done: %d rows, %d bins not measurable, %s'
        % (len(satirlar), len(olculemeyen), sure_metni(time.time() - t_basla)))


# ===========================================================================
# MODUL 2 - IC TUTARLILIK
# ===========================================================================
# SORU: ayni sayi birden fazla dosyada geciyorsa hepsi ayni mi?
#
# NEDEN: bu projede karar tablolari elle ve betikle DEFALARCA guncellendi.
# Bir dosyada duzeltilen bir deger digerinde eski haliyle kaldiginda, hangi
# dosyaya bakildigina gore farkli bir siparis karari cikiyor. Gecmiste ayni
# kutu bir tabloda UYE, baska tabloda RAKIP sayilmisti; asagidaki M2-UYE-RAKIP
# kontrolu tam olarak bu deseni arar.
#
# YONTEM: her kaynaktan (hedef, alan) -> deger cikarilir, ayni (hedef, alan)
# icin farkli deger veren kaynaklar CELISKI olarak listelenir. Alan eslemesi
# ACIKCA yazilir, sutun ADINA gore otomatik eslenmez: ornegin "R" sutunu
# ciftler.tsv'de GERI PRIMER, ESIK_VE_OLCUT'ta BOLLUK ORANI demektir ve bu
# ikisini adina bakip eslestirmek sahte bir celiski uretirdi.
# ---------------------------------------------------------------------------

def _primer_norm(s):
    u"""Primer dizisini karsilastirma icin normalle: buyuk harf, bosluksuz."""
    return re.sub(r'[^A-Z]', u'', unicode_(s).upper())


def _ad_norm(s):
    u"""Hedef adini karsilastirma icin normalle (bosluk/alt cizgi/buyuk-kucuk)."""
    return re.sub(r'[^a-z0-9]', u'', unicode_(s).lower())


def modul_2_ic_tutarlilik(kay, rap):
    M = u'2 IC TUTARLILIK'

    # (kaynak_adi, yol, anahtar_sutun, {kanonik_alan: sutun_adi}, tip)
    # tip: 'primer' | 'sayi' | 'metin'
    ESLEME = [
        (u'ciftler.tsv', kay.ciftler, u'hedef', {
            u'ileri_primer': u'F', u'geri_primer': u'R'}, ),
        (u'NIHAI_SIPARIS', kay.nihai_siparis, u'hedef', {
            u'ileri_primer': u'F', u'geri_primer': u'R', u'urun_bp': u'urun_bp',
            u'kraken_etiketi': u'kraken_etiketi', u'olculen_kimlik': u'olculen_kimlik',
            u'dCq': u'dCq'}, ),
        (u'SIPARIS_LISTESI', kay.siparis_listesi, u'hedef', {
            u'ileri_primer': u'F', u'geri_primer': u'R', u'urun_bp': u'urun_bp',
            u'kraken_etiketi': u'kraken_etiketi', u'olculen_kimlik': u'olculen_kimlik',
            u'dCq': u'dCq_karsiligi'}, ),
        (u'ESIK_VE_OLCUT', kay.esik_olcut, u'hedef', {
            u'dCq': u'dCq_olculen'}, ),
    ]
    ALAN_TIPI = {u'ileri_primer': u'primer', u'geri_primer': u'primer',
                 u'urun_bp': u'sayi', u'dCq': u'sayi',
                 u'kraken_etiketi': u'metin', u'olculen_kimlik': u'metin'}

    # deger[(hedef, alan)] = [(kaynak, ham_deger, satir_no)]
    deger = collections.defaultdict(list)
    okunan_kaynak = 0
    for kaynak, yol, anahtar, alanlar in ESLEME:
        satirlar = tsv_oku(yol)
        if satirlar is None:
            _gs = gecersiz_isareti(yol)
            rap.atla(M, u'M2-KAYNAK-YOK', u'%s okunabilmeli' % kaynak,
                     (u'GECERSIZ isaretli, kaynak sayilmadi: %s' % _gs) if _gs
                     else u'dosya yok', yol)
            continue
        if not satirlar:
            rap.ekle(M, u'M2-KAYNAK-BOS', CIDDI,
                     u'%s veri satiri icermeli' % kaynak,
                     u'dosya BOS (yalniz yorum/baslik)', yol)
            continue
        okunan_kaynak += 1
        for r in satirlar:
            h = (r.get(anahtar) or u'').strip()
            if not h:
                continue
            for kanonik, sutun in alanlar.items():
                if sutun not in r:
                    rap.ekle(M, u'M2-SUTUN-YOK', UYARI,
                             u'%s dosyasinda "%s" sutunu bulunmali' % (kaynak, sutun),
                             u'sutun yok - bu alan karsilastirmaya GIRMEDI',
                             yol)
                    continue
                ham = (r.get(sutun) or u'').strip()
                if ham:
                    deger[(_ad_norm(h), kanonik)].append((kaynak, ham, r.get('_satir')))

    if okunan_kaynak < 2:
        rap.atla(M, u'M2-CAPRAZ', u'capraz karsilastirma icin en az iki kaynak gerekli',
                 u'yalnizca %d kaynak okunabildi' % okunan_kaynak, u'-')
        return

    # --- EXCEL: panel sayfasindaki primerleri de karsilastirmaya kat
    sayfalar = xlsx_sayfalari(kay.panel_xlsx)
    if sayfalar is None:
        rap.atla(M, u'M2-XLSX-YOK', u'panel Excel dosyasi okunabilmeli',
                 u'dosya yok', kay.panel_xlsx)
    elif isinstance(sayfalar, type(u'')):
        rap.atla(M, u'M2-XLSX', u'panel Excel dosyasi okunabilmeli', sayfalar,
                 kay.panel_xlsx)
    else:
        bulundu = False
        for ad, satirlar in sayfalar.items():
            if not satirlar:
                continue
            bas = [unicode_(c).strip() for c in satirlar[0]]
            bl = [b.lower() for b in bas]
            if u'hedef' not in bl:
                continue
            ih = next((i for i, b in enumerate(bl) if b.startswith(u'ileri primer')), None)
            gh = next((i for i, b in enumerate(bl) if b.startswith(u'geri primer')), None)
            if ih is None and gh is None:
                continue
            hi = bl.index(u'hedef')
            bulundu = True
            for r in satirlar[1:]:
                if len(r) <= hi or not unicode_(r[hi]).strip():
                    continue
                h = _ad_norm(r[hi])
                if ih is not None and ih < len(r) and unicode_(r[ih]).strip():
                    deger[(h, u'ileri_primer')].append((u'Excel/%s' % ad, r[ih], None))
                if gh is not None and gh < len(r) and unicode_(r[gh]).strip():
                    deger[(h, u'geri_primer')].append((u'Excel/%s' % ad, r[gh], None))
            break
        if not bulundu:
            rap.atla(M, u'M2-XLSX-SAYFA',
                     u'Excel icinde "Hedef" + "Ileri primer" sutunlu bir sayfa olmali',
                     u'boyle bir sayfa bulunamadi', kay.panel_xlsx)

    # --- CELISKI ARAMA
    celiski = 0
    karsilastirilan = 0
    for (h, alan), kayitlar in sorted(deger.items()):
        if len(kayitlar) < 2:
            continue
        karsilastirilan += 1
        tip = ALAN_TIPI.get(alan, u'metin')
        if tip == u'primer':
            norm = [(k, _primer_norm(v), s) for k, v, s in kayitlar]
            ayri = set(n for _k, n, _s in norm)
        elif tip == u'sayi':
            norm = [(k, sayi(v), s) for k, v, s in kayitlar]
            sayilar = [n for _k, n, _s in norm if n is not None]
            # Sayilarda 0,01'lik yuvarlama farki celiski sayilmaz; daha buyugu sayilir.
            ayri = set()
            if sayilar and (max(sayilar) - min(sayilar)) > 0.011:
                ayri = set(sayilar)
            if len(sayilar) != len(norm):
                ayri.add(None)
        else:
            norm = [(k, unicode_(v).strip(), s) for k, v, s in kayitlar]
            ayri = set(n.lower() for _k, n, _s in norm)
        if len(ayri) <= 1:
            continue
        celiski += 1
        ayrinti = u' | '.join(
            u'%s%s = %s' % (k, (u' (satir %s)' % s) if s else u'', v)
            for k, v, s in kayitlar)
        agir = KRITIK if alan in (u'ileri_primer', u'geri_primer', u'urun_bp') else CIDDI
        rap.ekle(M, u'M2-CELISKI', agir,
                 u'"%s" hedefinin "%s" degeri butun dosyalarda AYNI olmali' % (h, alan),
                 ayrinti,
                 u'; '.join(sorted(set(k for k, _v, _s in kayitlar))),
                 u'Siparis hangi dosyaya bakildigina gore degisir.'
                 if agir == KRITIK else u'')

    rap.olcum[u'M2 karsilastirilan alan'] = u'%d (hedef,alan) cifti, %d celiski' % (
        karsilastirilan, celiski)
    if karsilastirilan == 0:
        rap.atla(M, u'M2-BOS', u'en az bir alan capraz karsilastirilmali',
                 u'iki kaynakta birden gecen hicbir (hedef, alan) bulunamadi', u'-')

    # --- AYNI KUTU HEM UYE HEM RAKIP MI (gecmiste yasandi)
    uyelik = tsv_oku(kay.hedef_uyelik)
    if uyelik is None:
        rap.atla(M, u'M2-UYE-RAKIP', u'hedef_uyelik.tsv okunabilmeli', u'dosya yok',
                 kay.hedef_uyelik)
    else:
        for r in uyelik:
            h = (r.get(u'hedef') or u'').strip()
            uye = set(x.strip() for x in (r.get(u'uye_taxid') or u'').split(u',') if x.strip())
            haric = set(x.strip() for x in (r.get(u'haric') or u'').split(u',') if x.strip())
            kesisim = uye & haric
            if kesisim:
                rap.ekle(M, u'M2-UYE-RAKIP', KRITIK,
                         u'bir taxid ayni hedefte hem uye hem haric olamaz',
                         u'%s: %s hem uye_taxid hem haric listesinde'
                         % (h, u', '.join(sorted(kesisim))),
                         u'%s (satir %s)' % (kay.hedef_uyelik, r.get('_satir')),
                         u'Ayrim orani bu taxid iki kez sayildigi icin yanlis.')


# ===========================================================================
# MODUL 3 - UYELIK BUTUNLUGU
# ===========================================================================
# SORU: her ciftin uye kumesi TANIMLI mi, ve o cift icin hesaplanan ΔCq
# gercekten O CIFTIN uyeligiyle mi hesaplanmis?
#
# NEDEN: Petriella_cinsi hedefinde tam bu olmustu - hedefin kendi uyelik satiri
# YOKTU, uyelik baska bir hedeften devralinmisti ve raporlanan ΔCq aslinda o
# baska hedefin olcumuydu. Sayi dolu gorunuyordu, dolayisiyla kimse fark
# etmemisti. Bos uyelikten sessizce sifir donmesi bu projenin en pahali
# desenlerinden biri.
# ---------------------------------------------------------------------------

def modul_3_uyelik(kay, rap):
    M = u'3 UYELIK'

    ciftler = tsv_oku(kay.ciftler)
    uyelik = tsv_oku(kay.hedef_uyelik)
    esik = tsv_oku(kay.esik_olcut)

    if ciftler is None:
        rap.atla(M, u'M3-CIFTLER', u'screening/ciftler.tsv okunabilmeli',
                 u'dosya yok', kay.ciftler)
        return
    if uyelik is None:
        rap.atla(M, u'M3-UYELIK', u'screening/hedef_uyelik.tsv okunabilmeli',
                 u'dosya yok', kay.hedef_uyelik)
        return

    uyelik_indeks = {}
    for r in uyelik:
        h = (r.get(u'hedef') or u'').strip()
        if h:
            uyelik_indeks[_ad_norm(h)] = r

    dcq_indeks = {}
    if esik:
        for r in esik:
            h = (r.get(u'hedef') or u'').strip()
            if h:
                dcq_indeks[_ad_norm(h)] = r

    # --- uye kumelerini imzala: birebir ayni kume iki hedefte tautoloji demek
    imza = collections.defaultdict(list)

    for r in ciftler:
        h = (r.get(u'hedef') or u'').strip()
        if not h:
            continue
        hn = _ad_norm(h)
        uye_ham = (r.get(u'uye_taksonlar') or u'').strip()
        uye = set(x.strip() for x in uye_ham.split(u',') if x.strip())
        olcu = (r.get(u'olcu_tipi') or u'').strip()
        durum = (r.get(u'uye_kumesi_durumu') or u'').strip()
        dsat = dcq_indeks.get(hn)
        dcq = sayi(dsat.get(u'dCq_olculen')) if dsat else None

        # 1) uye kumesi BOS mu
        #    ONEMLI AYRIM: KAPSAM olcen bir hedefte (evrensel primerler) bos uye
        #    kumesi TASARIM GEREGIDIR, hata degil - o hedef bir soyu digerinden
        #    ayirmaya calismaz, alan genelinde kapsam olcer. Bos uyelik ancak
        #    (a) AYRIM olcen bir hedefte, ya da (b) uyelige dayanan bir dCq
        #    raporlanmissa hatadir. Bu ayrimi yapmayan bir denetci, saglam
        #    tasarimi hata diye bagirir ve gercek hatalar arasinda kaybolur.
        if not uye:
            kapsam_olcusu = (olcu.lower().startswith(u'kapsam') or
                             durum.upper() == u'KAPSAM_OLCUSU')
            if dcq is not None:
                rap.ekle(M, u'M3-UYELIK-BOS', KRITIK,
                         u'uyeligi bos bir hedef icin uyelige dayanan dCq '
                         u'raporlanmamali',
                         u'%s: uye_taksonlar BOS ama dCq %s raporlanmis '
                         u'(olcu tipi: %s)' % (h, vir(dcq), olcu or u'?'),
                         u'%s (satir %s)' % (kay.ciftler, r.get('_satir')),
                         u'Bos uyelikle hesaplanan ayrim orani ANLAMSIZDIR.')
            elif kapsam_olcusu:
                rap.ekle(M, u'M3-UYELIK-BOS-KAPSAM', BILGI,
                         u'kapsam olcen hedeflerde uye kumesi bos olabilir',
                         u'%s: uye_taksonlar bos, olcu tipi "%s" - bu TASARIM '
                         u'geregidir, ayrim hesabina girmez' % (h, olcu),
                         u'%s (satir %s)' % (kay.ciftler, r.get('_satir')))
            else:
                rap.ekle(M, u'M3-UYELIK-BOS', KRITIK,
                         u'ayrim olcen her ciftin uye taksonlari tanimli olmali',
                         u'%s: uye_taksonlar BOS, olcu tipi "%s"' % (h, olcu or u'?'),
                         u'%s (satir %s)' % (kay.ciftler, r.get('_satir')),
                         u'Bos uyelikle hesaplanan ayrim orani ANLAMSIZDIR.')
        else:
            imza[u','.join(sorted(uye))].append((h, olcu))

        # 2) hedefin KENDI uyelik satiri var mi
        if hn not in uyelik_indeks:
            agir = KRITIK if dcq is not None else CIDDI
            rap.ekle(M, u'M3-UYELIK-SATIRI-YOK', agir,
                     u'her hedefin hedef_uyelik.tsv icinde KENDI satiri olmali',
                     u'%s: kendi uyelik satiri YOK%s'
                     % (h, u'; buna ragmen dCq %s raporlanmis - bu olcum baska bir '
                            u'hedefin uyeligiyle yapilmis olabilir' % vir(dcq)
                        if dcq is not None else u''),
                     kay.hedef_uyelik,
                     u'Petriella_cinsi desenidir: uyelik devralinmis, sayi dolu '
                     u'gorunuyor ama o ciftin olcumu degil.')

        # 3) uyelik durumu bir sorun bildiriyor mu
        if durum and durum.upper() not in (u'PANELLE_TUTUYOR', u'YENIDEN_KURULDU',
                                           u'KAPSAM_OLCUSU'):
            rap.ekle(M, u'M3-UYELIK-DURUM', CIDDI,
                     u'uye_kumesi_durumu bilinen ve gecerli bir deger olmali',
                     u'%s: uye_kumesi_durumu = "%s"' % (h, durum),
                     u'%s (satir %s)' % (kay.ciftler, r.get('_satir')))

        # 4) ESIK tablosu uyeligi gecersiz diyorsa ama dCq yine de kullanilmissa
        if dsat:
            gecerli = (dsat.get(u'uyelik_gecerli_mi') or u'').strip().upper()
            if gecerli in (u'HAYIR', u'YOK') and dcq is not None:
                rap.ekle(M, u'M3-DCQ-GECERSIZ-UYELIK', KRITIK,
                         u'uyeligi gecersiz bir hedef icin dCq raporlanmamali',
                         u'%s: uyelik_gecerli_mi=%s ama dCq_olculen=%s ve hukum "%s"'
                         % (h, gecerli, vir(dcq), (dsat.get(u'YENI_HUKUM') or u'?')),
                         u'%s (satir %s)' % (kay.esik_olcut, dsat.get('_satir')),
                         u'Bu ΔCq tanimsiz uyelikle hesaplanmistir.')

    # 5) birebir AYNI uye kumesine sahip iki "ayrim" hedefi -> TOTOLOJI
    for anahtar, hedefler in imza.items():
        ayrim = [h for h, o in hedefler if o.lower().startswith(u'ayrim')]
        if len(hedefler) > 1 and ayrim:
            rap.ekle(M, u'M3-TOTOLOJIK-UYELIK', CIDDI,
                     u'ayrim olcen iki hedefin uye kumeleri farkli olmali',
                     u'%s hedeflerinin uye kumeleri BIREBIR AYNI (%d taxid); '
                     u'bunlardan ayrim olcenler: %s'
                     % (u', '.join(h for h, _o in hedefler),
                        len(anahtar.split(u',')), u', '.join(ayrim)),
                     kay.ciftler,
                     u'Ayni kumeyi kendinden ayirmaya calisan bir olcum '
                     u'TOTOLOJIKTIR; cikan oran gercek bir ayrim degildir.')

    rap.olcum[u'M3 denetlenen cift'] = u'%d cift, %d uyelik satiri' % (
        len(ciftler), len(uyelik_indeks))


# ===========================================================================
# MODUL 4 - LITERATUR KURALLARI
# ===========================================================================
# 55 kaynakli literatur degerlendirmesinden cikan SAYISAL kurallari panele
# uygular. Kurallar burada SABIT yazilmaz; once LITERATUR dosyasindan
# dogrulanir - dosya baska bir sayi soyluyorsa kural BAYAT demektir ve bu da
# bir bulgudur.
# ---------------------------------------------------------------------------

# Bollukla agirlikli esik: gerekli_dCq = max(log2(R) + EK, TABAN)
LIT_EK = 4.3
LIT_TABAN = 3.32


def modul_4_literatur(kay, rap):
    M = u'4 LITERATUR'

    lit = metin_oku(kay.literatur)
    if lit is None:
        rap.atla(M, u'M4-LITERATUR-YOK', u'LITERATUR_2026-08-07.md okunabilmeli',
                 u'dosya yok', kay.literatur)
        lit = u''
    else:
        # KURALIN KENDISI BAYAT MI? Sabitlerimiz literaturde geciyor mu.
        for sabit, ad in ((LIT_EK, u'log2(R) eklentisi'), (LIT_TABAN, u'sabit taban')):
            desen = unicode_(sabit).replace(u'.', u'[.,]')
            if not re.search(desen, lit):
                rap.ekle(M, u'M4-KURAL-BAYAT', CIDDI,
                         u'%s (%s) degeri literatur dosyasinda gecmeli' % (ad, vir(sabit)),
                         u'%s degeri LITERATUR_2026-08-07.md icinde bulunamadi'
                         % vir(sabit), kay.literatur,
                         u'Denetcinin kurali literaturden ayrismis olabilir.')

    esik = tsv_oku(kay.esik_olcut)
    siparis = tsv_oku(kay.nihai_siparis)

    # --- KURAL 1: bollukla agirlikli esik
    if esik is None:
        rap.atla(M, u'M4-KURAL1', u'gerekli_dCq = max(log2(R)+%s, %s) dogrulanmali'
                 % (vir(LIT_EK), vir(LIT_TABAN)), u'ESIK_VE_OLCUT dosyasi yok',
                 kay.esik_olcut)
    else:
        denetlenen = 0
        for r in esik:
            h = (r.get(u'hedef') or u'').strip()
            R = sayi(r.get(u'R'))
            yazan = sayi(r.get(u'gerekli_dCq'))
            olculen = sayi(r.get(u'dCq_olculen'))
            durum = (r.get(u'yeni_kural_durum') or u'').strip().upper()
            if R is None or R <= 0:
                rap.ekle(M, u'M4-R-YOK', CIDDI,
                         u'bollukla agirlikli esik icin R (bolluk orani) tanimli olmali',
                         u'%s: R = "%s" (sayiya cevrilemedi ya da <= 0)'
                         % (h, r.get(u'R')),
                         u'%s (satir %s)' % (kay.esik_olcut, r.get('_satir')),
                         u'Bu hedefte gerekli_dCq HESAPLANAMAZ.')
                continue
            hesap = max(math.log(R, 2) + LIT_EK, LIT_TABAN)
            denetlenen += 1
            if yazan is None:
                rap.ekle(M, u'M4-GEREKLI-YOK', CIDDI,
                         u'gerekli_dCq sutunu dolu olmali',
                         u'%s: gerekli_dCq bos, hesaplanan %s' % (h, vir(hesap)),
                         u'%s (satir %s)' % (kay.esik_olcut, r.get('_satir')))
            elif abs(yazan - hesap) > 0.02:
                rap.ekle(M, u'M4-ESIK-YANLIS', KRITIK,
                         u'gerekli_dCq = max(log2(R)+%s ; %s) = %s'
                         % (vir(LIT_EK), vir(LIT_TABAN), vir(hesap)),
                         u'%s: tabloda %s yaziyor (R=%s)' % (h, vir(yazan), vir(R, 3)),
                         u'%s (satir %s)' % (kay.esik_olcut, r.get('_satir')),
                         u'Esik yanlissa GECER/KALIR hukmu de yanlis.')
            # hukum, kendi sayilariyla tutuyor mu
            if olculen is not None and yazan is not None:
                bekl = u'GECER' if olculen >= yazan else u'KALIR'
                if durum and durum != bekl:
                    rap.ekle(M, u'M4-HUKUM-CELISIK', KRITIK,
                             u'dCq %s, gerekli %s ise hukum "%s" olmali'
                             % (vir(olculen), vir(yazan), bekl),
                             u'%s: yeni_kural_durum = "%s"' % (h, durum),
                             u'%s (satir %s)' % (kay.esik_olcut, r.get('_satir')))
        rap.olcum[u'M4 esik kurali'] = u'%d hedefte yeniden hesaplandi' % denetlenen
        if denetlenen == 0:
            rap.atla(M, u'M4-KURAL1-BOS',
                     u'en az bir hedefte gerekli_dCq yeniden hesaplanmali',
                     u'hicbir hedefte gecerli R bulunamadi', kay.esik_olcut)

    # --- KURAL 2: 3' son iki baz olcutu TEK BASINA hukum vermemeli
    if siparis is None:
        rap.atla(M, u'M4-KURAL2', u"3' son iki baz olcutu hukum vermemeli",
                 u'NIHAI_SIPARIS okunmadi: dosya yok ya da GECERSIZ isaretli', kay.nihai_siparis)
    else:
        UC_DESEN = re.compile(r"3'\s*son\s*iki|son\s*iki\s*baz|3'\s*ucundaki?\s*iki", re.I)
        for r in siparis:
            h = (r.get(u'hedef') or u'').strip()
            katman = (r.get(u'hukmu_veren_katman') or u'')
            gerekce = (r.get(u'GEREKCE') or u'')
            if UC_DESEN.search(katman):
                rap.ekle(M, u'M4-UC-BAZ-HUKUM', CIDDI,
                         u"3' son iki baz olcutu TEK BASINA hukum vermemeli "
                         u'(literatur: bolum 3)',
                         u'%s: hukmu veren katman "%s"' % (h, katman.strip()[:160]),
                         u'%s (satir %s)' % (kay.nihai_siparis, r.get('_satir')),
                         u"Literatur, 3' son iki baz olcutunun tek basina hedef disi "
                         u'urun ongoremedigini soyluyor.')
            elif UC_DESEN.search(gerekce) and not re.search(
                    r'dCq|ayrim|kapsam|MFE|erime', gerekce, re.I):
                rap.ekle(M, u'M4-UC-BAZ-TEK-GEREKCE', UYARI,
                         u"3' son iki baz olcutu yaninda olculmus bir kanit da olmali",
                         u'%s: gerekce yalnizca 3\' olcutune dayaniyor: "%s"'
                         % (h, gerekce.strip()[:160]),
                         u'%s (satir %s)' % (kay.nihai_siparis, r.get('_satir')))

    # --- KURAL 3: evrensel primerler UCLU metrikle degerlendirilmeli
    if siparis is not None:
        UCLU = ((u'kapsam', re.compile(r'kapsam', re.I)),
                (u'filum spektrumu', re.compile(r'filum|phylum|spektrum', re.I)),
                (u'organel orani', re.compile(r'organel|mitokondri|kloroplast|plastid', re.I)))
        for r in siparis:
            h = (r.get(u'hedef') or u'').strip()
            if not re.search(r'universal|evrensel', h, re.I):
                continue
            butun = u' '.join(unicode_(v) for k, v in r.items() if k != '_satir')
            eksik = [ad for ad, d in UCLU if not d.search(butun)]
            if eksik:
                rap.ekle(M, u'M4-EVRENSEL-METRIK', CIDDI,
                         u'evrensel primerler UCLU metrikle degerlendirilmeli '
                         u'(kapsam, filum spektrumu, organel orani) - literatur bolum 1-2',
                         u'%s: satirda su olcut(ler) hic gecmiyor: %s'
                         % (h, u', '.join(eksik)),
                         u'%s (satir %s)' % (kay.nihai_siparis, r.get('_satir')),
                         u'Evrensel primeri "klad disi sayisi" ile degerlendirmek '
                         u'literature aykiridir; klad disi olmasi beklenir.')

    # --- KURAL 4: MIQE - in siliko ve deneysel dogrulama AYRI zorunlu maddeler
    #
    # DENEYSEL sart bu projede AYRI bir dosyada tutuluyor (SIPARIS_LISTESI.tsv
    # icindeki LABORATUVARDA_NE_YAPILMALI sutunu). Yalniz NIHAI_SIPARIS'e bakan
    # bir kontrol, sart baska dosyada yazili oldugu halde "eksik" derdi. Bu
    # yuzden iki dosya hedef adiyla BIRLESTIRILIR ve sart ikisinde de yoksa
    # bulgu uretilir.
    if siparis is not None:
        SIL = re.compile(r'in.?siliko|in.?silico|MFE|hizalama|nt |NCBI', re.I)
        DEN = re.compile(r'jel|gel|erime\s*egrisi|dizileme|amplikon\s*dizile|'
                         r'deneysel|laboratuvar|NTC|qPCR\s*kosu', re.I)
        siparis_l = tsv_oku(kay.siparis_listesi) or []
        ek_metin = collections.defaultdict(list)
        for r in siparis_l:
            ek_metin[_ad_norm(r.get(u'hedef') or u'')].append(
                u' '.join(unicode_(v) for k, v in r.items() if k != '_satir'))
        if not siparis_l:
            rap.atla(M, u'M4-MIQE-IKINCI-KAYNAK',
                     u'deneysel dogrulama sarti SIPARIS_LISTESI.tsv ile birlikte '
                     u'denetlenmeli',
                     u'SIPARIS_LISTESI.tsv okunamadi; MIQE kontrolu yalniz tek '
                     u'dosyaya bakti ve eksik gorunebilir', kay.siparis_listesi)
        for r in siparis:
            h = (r.get(u'hedef') or u'').strip()
            hukum = (r.get(u'HUKUM') or u'').strip().upper()
            if hukum in (u'', u'SIPARIS EDILMEZ'):
                continue
            butun = u' '.join(unicode_(v) for k, v in r.items() if k != '_satir')
            butun += u' ' + u' '.join(ek_metin.get(_ad_norm(h), []))
            var_sil = bool(SIL.search(butun))
            var_den = bool(DEN.search(butun))
            if var_sil and not var_den:
                rap.ekle(M, u'M4-MIQE-DENEYSEL-YOK', CIDDI,
                         u'MIQE in siliko ve DENEYSEL dogrulamayi AYRI zorunlu '
                         u'maddeler sayar (literatur bolum 4)',
                         u'%s (hukum: %s): in siliko kanit var, deneysel dogrulama '
                         u'sarti satirda hic gecmiyor' % (h, hukum),
                         u'%s (satir %s)' % (kay.nihai_siparis, r.get('_satir')),
                         u'Siparis edilen her cift icin jel/erime egrisi/dizileme '
                         u'sarti acikca yazilmali.')
            elif not var_sil and not var_den:
                rap.ekle(M, u'M4-MIQE-KANIT-YOK', CIDDI,
                         u'siparise giden her ciftte en az bir dogrulama kanidi olmali',
                         u'%s (hukum: %s): ne in siliko ne deneysel kanit gecmiyor'
                         % (h, hukum),
                         u'%s (satir %s)' % (kay.nihai_siparis, r.get('_satir')))


# ===========================================================================
# MODUL 5 - BILINEN HATA DESENLERI
# ===========================================================================
# Bu projede TEKRAR TEKRAR ciktigi icin ayri bir modul hak eden dokuz desen.
# Her desen icin ayri bir kontrol var ve her biri "ne bekleniyordu / ne
# bulundu" ciftiyle raporlanir. Bir desen hic bulgu uretmiyorsa bu, o desenin
# olculdugu ve temiz cikti anlamina gelir - kontrolun kosmadigi anlamina
# GELMEZ; kosmayanlar ATLANDI olarak ayrica sayilir.
# ---------------------------------------------------------------------------

# Kod taramasinda aranan maskeleme desenleri.
_MASKE_PY = [
    (re.compile(r'except[^\n:]*:\s*\n\s*pass\b'), u'except ... : pass  (hata yutuluyor)'),
    (re.compile(r'except[^\n:]*:\s*\n\s*continue\b'), u'except ... : continue  (hata yutuluyor)'),
    (re.compile(r'\bsys\.exit\(\s*0\s*\)'), u'sys.exit(0)  (kosulsuz basarili cikis)'),
    (re.compile(r'subprocess\.(call|run)\((?![^)]*check\s*=\s*True)[^)]*\)'),
     u'subprocess check=True olmadan  (alt surecin cikis kodu okunmuyor)'),
]
_MASKE_BAT = [
    (re.compile(r'\|\|\s*(true|ver\b)', re.I), u'|| true  (hata bastiriliyor)'),
    (re.compile(r'>\s*nul\s+2>&1', re.I), u'> nul 2>&1  (hata ciktisi yok ediliyor)'),
    (re.compile(r'\bexit\s*/b\s*0\b', re.I), u'exit /b 0  (kosulsuz basarili cikis)'),
]

# "Tavan degeri sayim sanildi" deseninde aranan bilinen tavanlar.
_TAVANLAR = (500, 1000, 3000, 120001, 200000, 100, 50)


def modul_5_desenler(kay, rap):
    M = u'5 HATA DESENLERI'

    siparis = tsv_oku(kay.nihai_siparis)
    ciftler = tsv_oku(kay.ciftler)
    uyelik = tsv_oku(kay.hedef_uyelik)
    hdisi = tsv_oku(kay.hedef_disi)

    # ---- DESEN 1: SESSIZ SIFIR ---------------------------------------
    # Desenin TAM tanimi: bir KATMAN kosmamis ya da hukum verememis, ama o
    # katmanin sayisal alani yine de 0 olarak yazilmis. Boyle bir 0 "hedef disi
    # bulunmadi" gibi okunur; oysa "bakilmadi" demektir.
    #
    # DIKKAT - burada gecmiste bir sinir hatasi yapmak cok kolay: "her 0 supheli"
    # demek 58 sahte bulgu uretiyor ve gercek olani bogar. O yuzden 0, YALNIZCA
    # kendi katmani "kosmadim" diyorsa bulgu sayilir. Alan -> katman eslemesi
    # sutun onekinden okunur.
    KATMAN_ONEKI = ((u'NCBI_', u'NCBI'), (u'MFE_', u'MFEprimer'),
                    (u'yerel_', u'yerel'))
    KOSMADI = re.compile(
        r'sinanmad|hukum\s*veremedi|kosulmad|[cç]al[iı][sş]t[iı]r[iı]lmad|'
        r'atland|BA[SŞ]ARISIZ|yap[iı]lamad|olculemedi', re.I)
    if siparis is None:
        rap.atla(M, u'M5-D1', u'sessiz sifir taramasi',
                 u'NIHAI_SIPARIS okunmadi: dosya yok ya da GECERSIZ isaretli',
                 kay.nihai_siparis)
    else:
        bakilan = 0
        for r in siparis:
            h = (r.get(u'hedef') or u'').strip()
            baglam = u' '.join([r.get(u'GEREKCE') or u'',
                                r.get(u'hukmu_veren_katman') or u'',
                                r.get(u'NCBI_durumu') or u''])
            # CUMLE CUMLE bakilir. Tek bir uzun gerekce icinde "MFEprimer 0
            # buldu. NCBI hukum veremedi." gibi IKI AYRI ifade olabilir; butun
            # metne birden bakan bir kontrol, MFE alanlarini da NCBI'nin
            # kosmamasi yuzunden suclardi (olculdu: 88 bulgunun 80'i boyleydi).
            cumleler = [c for c in re.split(r'[.;|]\s*', baglam) if c.strip()]
            # Katman basina TEK gerekce cumlesi tutulur; ayni katmani iki ayri
            # cumlede sucladigimizda rapor ayni bulguyu iki kez basiyordu.
            kosmayan = collections.OrderedDict()
            for c in cumleler:
                if not KOSMADI.search(c):
                    continue
                for _onek, k in KATMAN_ONEKI:
                    if k.lower() in c.lower() and k not in kosmayan:
                        kosmayan[k] = c.strip()[:140]
            sifir_alan = collections.defaultdict(list)
            for alan, ham in sorted(r.items()):
                if alan == '_satir' or not unicode_(ham).strip():
                    continue
                katman = next((k for onek, k in KATMAN_ONEKI
                               if alan.startswith(onek)), None)
                if katman is None:
                    continue
                bakilan += 1
                v = sayi(ham)
                if v is not None and abs(v) < 1e-9:
                    sifir_alan[katman].append(alan)
            # Katman basina TEK bulgu: alan alan tekrarlamak raporu sisirir.
            for katman, cumle in kosmayan.items():
                alanlar = sifir_alan.get(katman)
                if not alanlar:
                    continue
                rap.ekle(M, u'M5-D1-SESSIZ-SIFIR', CIDDI,
                         u'kosmamis bir katmanin sayisal alani 0 degil '
                         u'"olculmedi" olarak yazilmali',
                         u'%s: "%s" katmani kosmamis ("%s") ama su alan(lar) 0 '
                         u'yaziyor: %s' % (h, katman, cumle, u', '.join(alanlar)),
                         u'%s (satir %s)' % (kay.nihai_siparis, r.get('_satir')),
                         u'0, "hedef disi bulunmadi" gibi okunur; oysa '
                         u'"bakilmadi" demektir.')
        rap.olcum[u'M5 D1 taranan katman alani'] = u'%d' % bakilan
        if bakilan == 0:
            rap.atla(M, u'M5-D1', u'sessiz sifir taramasi',
                     u'katman oneki tasiyan sayisal alan bulunamadi',
                     kay.nihai_siparis)

    # ---- DESEN 2: HEDEFIN KENDI UYELERINI HEDEF DISI SAYMA -----------
    if uyelik is None or hdisi is None:
        rap.atla(M, u'M5-D2', u'hedefin kendi uyelerini hedef disi saymamasi',
                 u'hedef_uyelik.tsv ya da HEDEF_DISI_AYRINTI yok',
                 kay.hedef_disi if uyelik else kay.hedef_uyelik)
    else:
        uye_harita = {}
        for r in uyelik:
            h = _ad_norm(r.get(u'hedef') or u'')
            uye_harita[h] = set(x.strip() for x in
                                (r.get(u'uye_taxid') or u'').split(u',') if x.strip())
        hedef_sut = next((s for s in (hdisi[0].keys() if hdisi else [])
                          if _ad_norm(s) == u'hedef'), None)
        if hedef_sut is None:
            rap.atla(M, u'M5-D2', u'hedefin kendi uyelerini hedef disi saymamasi',
                     u'HEDEF_DISI_AYRINTI dosyasinda "hedef" sutunu yok', kay.hedef_disi)
        else:
            sayac = 0
            for r in hdisi:
                h = _ad_norm(r.get(hedef_sut) or u'')
                uye = uye_harita.get(h)
                if not uye:
                    continue
                butun = u' '.join(unicode_(v) for k, v in r.items() if k != '_satir')
                gecen = [t for t in uye if re.search(r'(?<!\d)%s(?!\d)' % re.escape(t), butun)]
                if gecen:
                    sayac += 1
                    rap.ekle(M, u'M5-D2-UYE-HEDEF-DISI', CIDDI,
                             u'bir hedefin KENDI uyesi o hedefin "hedef disi" '
                             u'listesinde yer almamali',
                             u'%s: kendi uyesi olan %s taxid\'i hedef disi satirinda geciyor'
                             % (r.get(hedef_sut), u', '.join(sorted(gecen))),
                             u'%s (satir %s)' % (kay.hedef_disi, r.get('_satir')),
                             u'Hedef disi sayisi sisirilmis, ozgulluk oldugundan '
                             u'kotu gorunuyor.')
            rap.olcum[u'M5 D2 taranan hedef disi satiri'] = u'%d satir, %d bulgu' % (
                len(hdisi), sayac)

    # ---- DESEN 3: MASKELENMIS CIKIS KODU -----------------------------
    tarandi = 0
    for klasor in kay.kod_klasorleri:
        if not os.path.isdir(klasor):
            continue
        for yol in sorted(glob.glob(os.path.join(klasor, u'*.py'))):
            if u'.orig' in yol:
                continue
            icerik = metin_oku(yol)
            if icerik is None:
                continue
            tarandi += 1
            for desen, ad in _MASKE_PY:
                n = len(desen.findall(icerik))
                if n:
                    rap.ekle(M, u'M5-D3-MASKE-PY', UYARI,
                             u'hata ve cikis kodlari maskelenmemeli',
                             u'%s: %d kez "%s"' % (os.path.basename(yol), n, ad),
                             yol,
                             u'Yutulan hata, kosunun basarili sanilmasina yol acar.')
    for yol in sorted(glob.glob(os.path.join(kay.kok, u'*.bat'))):
        # Kosan .bat dosyalari OKUNUR ama DEGISTIRILMEZ.
        icerik = metin_oku(yol)
        if icerik is None:
            continue
        tarandi += 1
        for desen, ad in _MASKE_BAT:
            n = len(desen.findall(icerik))
            if n:
                rap.ekle(M, u'M5-D3-MASKE-BAT', UYARI,
                         u'bat dosyalari alt surecin cikis kodunu bastirmamali',
                         u'%s: %d kez "%s"' % (os.path.basename(yol), n, ad), yol)
    if tarandi == 0:
        rap.atla(M, u'M5-D3', u'maskelenmis cikis kodu taramasi',
                 u'taranacak kod dosyasi bulunamadi', u'; '.join(kay.kod_klasorleri))
    else:
        rap.olcum[u'M5 D3 taranan kod dosyasi'] = u'%d' % tarandi

    # ---- DESEN 4: BAYAT KONTROL NOKTASI ------------------------------
    kn_klasorleri = [d for d in glob.glob(os.path.join(kay.kok, u'*', u'kontrol'))
                     if os.path.isdir(d)]
    if not kn_klasorleri:
        rap.atla(M, u'M5-D4', u'bayat kontrol noktasi taramasi',
                 u'hicbir */kontrol klasoru bulunamadi', kay.kok)
    else:
        # Girdi olarak konsensus indeksi ve panel tanimi alinir: kontrol
        # noktasi bunlardan ESKIYSE, o kosunun sonucu bayattir.
        girdiler = [kay.konsensus_indeks, kay.ciftler, kay.hedef_uyelik]
        en_yeni = 0
        for g in girdiler:
            if os.path.exists(g):
                en_yeni = max(en_yeni, os.path.getmtime(g))
        bayat = 0
        toplam = 0
        for d in kn_klasorleri:
            dosyalar = glob.glob(os.path.join(d, u'*'))
            toplam += len(dosyalar)
            eski = [f for f in dosyalar
                    if os.path.isfile(f) and os.path.getmtime(f) < en_yeni]
            if eski:
                bayat += len(eski)
                rap.ekle(M, u'M5-D4-BAYAT-KN', CIDDI,
                         u'kontrol noktalari girdilerinden (panel/uyelik/konsensus) '
                         u'YENI olmali',
                         u'%s: %d/%d kontrol noktasi girdilerden eski'
                         % (os.path.relpath(d, kay.kok), len(eski), len(dosyalar)), d,
                         u'Bayat kontrol noktasi, degismis girdiyle eski sonucu '
                         u'yeniden kullandirir.')
        rap.olcum[u'M5 D4 kontrol noktasi'] = u'%d dosya, %d bayat' % (toplam, bayat)

    # ---- DESEN 5: TAVAN DEGERININ SAYIM SANILMASI --------------------
    # Yalniz SAYIM alanlarina bakilir. Urun uzunlugu (urun_bp = 100) bir sayim
    # degildir ve tavanla ilgisi yoktur; butun sutunlara bakan bir kontrol
    # boyle sahte bulgular uretiyordu.
    # Ayrica satir zaten "SONUC TAVANI" diyorsa proje bunu bilerek isaretlemis
    # demektir; bilinen bir sorunu yeniden bagirmak raporu kirletir.
    SAYIM_DESENI = re.compile(
        r'(sayi|sayisi|adet|urun$|_urun|hedef_disi|klad_disi|kayit|hit|isabet)', re.I)
    if siparis is None:
        rap.atla(M, u'M5-D5', u'tavan degeri taramasi',
                 u'NIHAI_SIPARIS okunmadi: dosya yok ya da GECERSIZ isaretli',
                 kay.nihai_siparis)
    else:
        bakilan = 0
        for r in siparis:
            kabul = u' '.join([r.get(u'NCBI_durumu') or u'',
                               r.get(u'GEREKCE') or u''])
            zaten_biliniyor = bool(re.search(r'tavan', kabul, re.I))
            for k, v in sorted(r.items()):
                if k == '_satir' or k.lower() in (u'urun_bp', u'sira'):
                    continue
                if not SAYIM_DESENI.search(k):
                    continue
                n = sayi(v)
                if n is None or int(n) != n or int(n) not in _TAVANLAR:
                    continue
                bakilan += 1
                rap.ekle(M, u'M5-D5-TAVAN',
                         BILGI if zaten_biliniyor else UYARI,
                         u'raporlanan SAYIM bir arac tavaniyla birebir ayni olmamali',
                         u'%s: %s = %d%s' % (
                             r.get(u'hedef'), k, int(n),
                             u' (satir bunu zaten "tavan" olarak isaretlemis)'
                             if zaten_biliniyor else
                             u' - tavan degeri, satirda tavan oldugu SOYLENMEMIS'),
                         u'%s (satir %s)' % (kay.nihai_siparis, r.get('_satir')),
                         u'Tavan degeri "bulunan sayi" degil "kesilen sayi"dir; '
                         u'gercek sayi daha buyuktur.')
        rap.olcum[u'M5 D5 tavana esit sayim'] = u'%d' % bakilan

    # ---- DESEN 6: TOTOLOJIK OLCUM ------------------------------------
    # M3 uyelik tarafini bakiyor; burada KANIT DOSYASI kendisini dogruluyor mu
    # sorusuna bakilir: bir hedefin kaniti olarak gosterilen dosya, o hedefin
    # hukmunun yazildigi dosyanin ta kendisiyse olcum kendini kanitliyordur.
    if siparis is not None:
        for r in siparis:
            kanit = (r.get(u'GEREKCE') or u'') + u' ' + (r.get(u'hukmu_veren_katman') or u'')
            if os.path.basename(kay.nihai_siparis) in kanit:
                rap.ekle(M, u'M5-D6-TOTOLOJI', CIDDI,
                         u'bir hukmun kaniti, hukmun yazildigi dosyadan BASKA bir '
                         u'olcum dosyasi olmali',
                         u'%s: kanit olarak kendi dosyasi (%s) gosteriliyor'
                         % (r.get(u'hedef'), os.path.basename(kay.nihai_siparis)),
                         u'%s (satir %s)' % (kay.nihai_siparis, r.get('_satir')))

    # ---- DESEN 7: AYNI MOTORU KULLANAN IKI KATMAN = BAGIMSIZ KANIT SANMA
    # Iki katman ayni hizalama cekirdegini kullaniyorsa BAGIMSIZ IKI KANIT
    # DEGILDIR: cekirdegin ortak hatasi iki kez oy verir.
    #
    # BAGIMLILIK TAHMIN EDILMEZ, KODDAN OKUNUR. Bir katmanin proje ici olup
    # olmadigi, o katmanin betiginin engine_gateway.py'yi ice aktarip aktarmadigina
    # bakilarak belirlenir. Disaridan gelen ucuncu araclar (MFEprimer ikilisi,
    # NCBI web hizmeti) proje motorunu kullanmaz ve BAGIMSIZDIR - bunlari
    # "ayni motor" saymak, saglam bir capraz dogrulamayi hata diye gostermek
    # olurdu.
    ici_motor = set()
    for klasor in kay.kod_klasorleri:
        for yol in glob.glob(os.path.join(klasor, u'*.py')):
            if u'.orig' in yol:
                continue
            icerik = metin_oku(yol) or u''
            if re.search(r'^\s*(from\s+\S*\s+)?import\s+.*\bmotor\b', icerik, re.M):
                ici_motor.add(os.path.splitext(os.path.basename(yol))[0])
    # Karar tablosunda gecen katman adlarinin hangi betige karsilik geldigi.
    KATMAN_BETIK = {u'yerel': u'motor', u'kapsam': u'motor', u'okuma': u'okuma_motoru',
                    u'in siliko': u'motor', u'uyelik': u'uyelik_denetimi'}
    DIS_ARAC = (u'mfeprimer', u'mfe', u'ncbi', u'blast', u'kraken', u'bracken')
    if siparis is None:
        rap.atla(M, u'M5-D7', u'katman bagimsizligi taramasi',
                 u'NIHAI_SIPARIS okunmadi: dosya yok ya da GECERSIZ isaretli',
                 kay.nihai_siparis)
    elif not ici_motor:
        rap.atla(M, u'M5-D7', u'katman bagimsizligi taramasi',
                 u'engine_gateway.py\'yi ice aktaran betik bulunamadi - bagimlilik '
                 u'haritasi cikarilamadi', u'; '.join(kay.kod_klasorleri))
    else:
        rap.olcum[u'M5 D7 engine_gateway.py kullanan betik'] = u'%d' % len(ici_motor)
        for r in siparis:
            katman = (r.get(u'hukmu_veren_katman') or u'').strip()
            if not katman:
                continue
            paylasan = collections.Counter()
            for ad, betik in KATMAN_BETIK.items():
                if betik not in ici_motor and betik != u'motor':
                    continue
                if re.search(re.escape(ad), katman, re.I):
                    # Dis arac adinin icinde gecen bir eslesmeyi sayma
                    parcalar = re.split(r'[+;,]', katman)
                    for p in parcalar:
                        if re.search(re.escape(ad), p, re.I) and not any(
                                d in p.lower() for d in DIS_ARAC):
                            paylasan[betik] += 1
                            break
            for betik, n in paylasan.items():
                if n >= 2:
                    rap.ekle(M, u'M5-D7-AYNI-MOTOR', CIDDI,
                             u'bir hukmu destekleyen iki katman BAGIMSIZ olmali',
                             u'%s: hukmu veren katmanlarin %d tanesi ayni cekirdegi '
                             u'(%s.py) kullaniyor - "%s"'
                             % (r.get(u'hedef'), n, betik, katman[:140]),
                             u'%s (satir %s)' % (kay.nihai_siparis, r.get('_satir')),
                             u'Ayni motorun iki ciktisi tek kanittir; ortak hatasi '
                             u'iki kez oy verir.')

    # ---- DESEN 8: DEJENERE BAZ KACAGI --------------------------------
    # Dejenere baz tasiyan primerlerin, dejenere bazi ISLEYEN bir olcumle
    # degerlendirildigi gosterilmeli. Aksi halde arac dejenere bazi N sayip
    # sessizce elemis olabilir.
    DEJ = set(u'RYSWKMBDHVN')
    if ciftler is None:
        rap.atla(M, u'M5-D8', u'dejenere baz kacagi taramasi', u'ciftler.tsv yok',
                 kay.ciftler)
    else:
        dej_hedef = []
        for r in ciftler:
            for alan in (u'F', u'R'):
                s = _primer_norm(r.get(alan))
                if s and set(s) & DEJ:
                    dej_hedef.append((r.get(u'hedef'), alan, s,
                                      sorted(set(s) & DEJ), r.get('_satir')))
        rap.olcum[u'M5 D8 dejenere primer'] = u'%d oligo' % len(dej_hedef)
        for h, alan, s, bazlar, sat in dej_hedef:
            ilgili = None
            if siparis:
                ilgili = next((x for x in siparis
                               if _ad_norm(x.get(u'hedef') or u'') == _ad_norm(h or u'')), None)
            butun = u' '.join(unicode_(v) for k, v in (ilgili or {}).items()
                              if k != '_satir')
            if not re.search(r'dejenere|degenerate|IUPAC|acilim|expand', butun, re.I):
                rap.ekle(M, u'M5-D8-DEJENERE', CIDDI,
                         u'dejenere baz tasiyan bir oligonun degerlendirmesi, '
                         u'dejenere bazlarin nasil islendigini belirtmeli',
                         u'%s / %s: %s (dejenere baz: %s) - karar satirinda dejenere '
                         u'baz islemesine dair ifade yok'
                         % (h, alan, s, u', '.join(bazlar)),
                         u'%s (satir %s)' % (kay.ciftler, sat),
                         u'Dejenere bazi N sayan bir arac bu oligoyu sessizce eler '
                         u've kapsam oldugundan dusuk cikar.')

    # ---- DESEN 9: YON HATASI -----------------------------------------
    ind = tsv_oku(kay.konsensus_indeks, yorum=None)
    if ind is None:
        rap.atla(M, u'M5-D9', u'yon hatasi taramasi',
                 u'konsensus_kanonik/INDEKS.tsv yok', kay.konsensus_indeks)
    else:
        yon_yok = [r for r in ind if not (r.get(u'eski_yon') or u'').strip()]
        cevrilen = [r for r in ind if (r.get(u'cevrildi') or u'').strip().lower()
                    in (u'evet', u'yes', u'1')]
        rap.olcum[u'M5 D9 yon'] = u'%d kutu, %d cevrildi, %d yon bilgisi yok' % (
            len(ind), len(cevrilen), len(yon_yok))
        if yon_yok:
            rap.ekle(M, u'M5-D9-YON-BILGISI-YOK', CIDDI,
                     u'kanonik indeksteki her kutunun yon bilgisi kayitli olmali',
                     u'%d kutuda eski_yon bos: %s'
                     % (len(yon_yok), u', '.join(r.get(u'kutu') or u'?'
                                                 for r in yon_yok[:8])),
                     kay.konsensus_indeks,
                     u'Ters yonlu bir konsensuste in siliko PCR SESSIZCE 0 urun verir.')
        # Kanonik olmayan, karisik yonlu klasorun hala okunuyor olmasi risklidir.
        #
        # B DUZELTMESI (2026-08-21): bu kontrol eskiden dosyanin TAMAMINDA duz
        # metin aramasi yapiyordu, yorum ve docstring'ler dahil. Aciklama
        # satirinda klasor adini ANAN dosyalar da RISKLI isaretleniyordu;
        # 2026-08-09 kosusunda bes yanlis pozitif uretti (ornegin
        # steps/generate_primer_candidates.py ve design_group_primers.py, ikisi de
        # yolu CLI argumani olarak alir, gomulu yol tasimaz).
        # Ayni kontrolun DOGRU surumu projede zaten vardi:
        # screening/yon_kod_taramasi.kod_govdesi(). Iki tarayici ayni
        # soruya farkli cevap veriyordu; artik ikisi de govdeye bakiyor.
        for klasor in kay.kod_klasorleri:
            for yol in sorted(glob.glob(os.path.join(klasor, u'*.py'))):
                if u'.orig' in yol:
                    continue
                ham = metin_oku(yol) or u''
                vurus = d9_karisik_klasor_yollari(ham, yol)
                if vurus and u'konsensus_kanonik' not in _kod_govdesi(ham):
                    rap.ekle(M, u'M5-D9-KARISIK-KLASOR', CIDDI,
                             u'konsensus okumalari KANONIK klasorden yapilmali',
                             u'%s: karisik yonlu "consensus sequences" klasorunu '
                             u'KOD ICINDE okuyor (satir %s), kanonik klasore hic '
                             u'deginmiyor'
                             % (os.path.relpath(yol, kay.kok),
                                u', '.join(str(i) for i, _ in vurus[:6])), yol,
                             u'O klasor karisik yonludur; ters yonlu konsensus '
                             u'sessizce 0 urun verir.')


# ===========================================================================
# MODUL 6 - VERITABANI SAGLIGI
# ===========================================================================
# SORU: her veritabani indeksi GERCEKTEN calisiyor mu?
#
# NEDEN: SILVA indeksi aylarca SESSIZCE sifir dondu. Bozuk indeksin kaniti hala
# duruyor (SILVA_138.2_SSURef_NR99.fasta.BOZUK_KANIT.txt): bozuk kurulumda
# "Sorting 19683 kmers" yaziyordu; 19683 = 3^9, yani dort bazdan biri dusmustu.
# Saglam indekste kmer_count = 4^9 = 262144 olmali. Bu yuzden k-mer sayisi
# DOGRUDAN denetlenir.
# ---------------------------------------------------------------------------

KMER_BEKLENEN = 262144      # 4^9 - saglam indeksin k-mer sayisi
KVALUE_BEKLENEN = 9


def modul_6_veritabani(kay, rap, baglanma_sinamasi=True):
    M = u'6 VERITABANI'

    K, hata = modul_yukle(kay.kimlik_dogrulama, 'kimlik_dogrulama_m6')
    if K is None:
        rap.atla(M, u'M6-VTB-LISTESI', u'VTB listesi identity_verification.py\'den okunmali',
                 hata, kay.kimlik_dogrulama)
        return
    T, _h = modul_yukle(kay.tum_kutu, 'tum_kutu_m6')
    beklenen_kayit = getattr(T, 'BEKLENEN_KAYIT', {}) if T else {}

    if not os.path.isdir(kay.refdb):
        rap.atla(M, u'M6-REFDB', u'REFERANS_DB klasoru bulunmali', u'klasor yok',
                 kay.refdb)
        return

    saglikli = 0
    denetlenen = 0
    for etiket, dosya, lokus, kullan, _not in K.VTB:
        if not kullan:
            continue                     # ikiz/altkume kumeler oylamaya girmiyor
        denetlenen += 1
        yol = os.path.join(kay.refdb, dosya)
        if not os.path.exists(yol):
            rap.ekle(M, u'M6-FASTA-YOK', KRITIK,
                     u'%s veritabani dosyasi bulunmali' % etiket,
                     u'%s yok' % dosya, kay.refdb,
                     u'Bu kaynak hicbir kimlik hukmune oy veremez.')
            continue
        if os.path.getsize(yol) == 0:
            rap.ekle(M, u'M6-FASTA-BOS', KRITIK,
                     u'%s veritabani dosyasi bos olmamali' % etiket,
                     u'%s: 0 bayt' % dosya, yol)
            continue

        # --- 1) k-mer indeksinin saglik gunlugu
        log = yol + u'.log'
        icerik = metin_oku(log)
        if icerik is None:
            rap.ekle(M, u'M6-INDEKS-GUNLUK-YOK', CIDDI,
                     u'%s icin indeks kurulum gunlugu bulunmali' % etiket,
                     u'%s.log yok - indeksin saglikli kuruldugu DOGRULANAMIYOR'
                     % dosya, kay.refdb,
                     u'Bu veritabani icin k-mer sayisi denetlenemedi.')
        else:
            m = re.findall(r'kvalue=(\d+),\s*kmer_count=(\d+)', icerik)
            if not m:
                rap.ekle(M, u'M6-KMER-SATIRI-YOK', CIDDI,
                         u'%s gunlugunde "kvalue=..., kmer_count=..." satiri olmali'
                         % etiket,
                         u'gunlukte k-mer satiri yok (indeks eski/farkli bicimde '
                         u'kurulmus olabilir)', log,
                         u'k-mer sayisi dogrulanamadi; SILVA vakasindaki sessiz '
                         u'sifir bu satirla yakalanmisti.')
            else:
                kv, kc = int(m[-1][0]), int(m[-1][1])
                bekl = 4 ** kv
                if kc != bekl or kv != KVALUE_BEKLENEN or kc != KMER_BEKLENEN:
                    rap.ekle(M, u'M6-KMER-BOZUK', KRITIK,
                             u'k-mer sayisi 4^%d = %d olmali'
                             % (KVALUE_BEKLENEN, KMER_BEKLENEN),
                             u'%s: kvalue=%d, kmer_count=%d (4^%d = %d)'
                             % (etiket, kv, kc, kv, bekl), log,
                             u'19683 = 3^9 bozuk indeksin imzasidir: bir baz dusmus '
                             u'demektir ve indeks SESSIZCE sifir doner.')

        # --- 2) indeks dosyasi FASTA'dan yeni mi (bayat indeks)
        indeksler = [yol + u'.primerqc.bin', yol + u'.primerqc']
        var = [i for i in indeksler if os.path.exists(i)]
        if not var:
            rap.ekle(M, u'M6-INDEKS-YOK', CIDDI,
                     u'%s icin k-mer indeksi (.primerqc.bin) bulunmali' % etiket,
                     u'indeks dosyasi yok - bu veritabani yalniz dogrudan FASTA '
                     u'akisiyla taranabilir, indekse dayanan katmanlar bu kaynakta '
                     u'SIFIR doner', kay.refdb,
                     u'SILVA LSU kumelerinde bu durum var; indekse dayanan bir '
                     u'olcum bu kaynagi sessizce atlar.')
        else:
            for i in var:
                if os.path.getmtime(i) < os.path.getmtime(yol):
                    rap.ekle(M, u'M6-INDEKS-BAYAT', CIDDI,
                             u'indeks, FASTA dosyasindan YENI olmali',
                             u'%s: indeks %s, FASTA %s (indeks daha ESKI)'
                             % (etiket,
                                time.strftime('%Y-%m-%d', time.localtime(os.path.getmtime(i))),
                                time.strftime('%Y-%m-%d', time.localtime(os.path.getmtime(yol)))),
                             i, u'FASTA guncellenmis ama indeks yeniden kurulmamis.')

        # --- 3) kayit sayisi beklenen kadar mi
        bekl_n = beklenen_kayit.get(etiket)
        fai = yol + u'.fai'
        gercek_n = None
        if os.path.exists(fai):
            try:
                with io.open(fai, encoding='utf-8', errors='replace') as fh:
                    gercek_n = sum(1 for s in fh if s.strip())
            except IOError:
                gercek_n = None
        if bekl_n and gercek_n is not None and gercek_n != bekl_n:
            rap.ekle(M, u'M6-KAYIT-SAYISI', CIDDI,
                     u'%s icinde %d kayit olmali' % (etiket, bekl_n),
                     u'.fai indeksi %d kayit gosteriyor' % gercek_n, fai)
        elif bekl_n and gercek_n is None:
            rap.ekle(M, u'M6-KAYIT-SAYISI-YOK', UYARI,
                     u'%s icin kayit sayisi dogrulanabilmeli' % etiket,
                     u'.fai yok, kayit sayisi DOGRULANMADI', kay.refdb)

        # --- 4) CANLI BAGLANMA SINAMASI: bilinen bir primer sifirdan buyuk
        #        baglanma donduruyor mu? "Her zaman gecen" bir kontrol olmasin
        #        diye bu sinav gercek dosyayi okur.
        if baglanma_sinamasi:
            sonuc = _baglanma_sinamasi(K, yol, lokus)
            if sonuc is None:
                rap.atla(M, u'M6-BAGLANMA', u'%s icinde bilinen bir korunmus dizi '
                         u'bulunmali' % etiket,
                         u'sinama kosulamadi (dosya okunamadi)', yol)
            elif sonuc == 0:
                rap.ekle(M, u'M6-BAGLANMA-SIFIR', KRITIK,
                         u'%s icinde bilinen korunmus dizi(ler) icin en az bir '
                         u'baglanma bulunmali' % etiket,
                         u'ilk 20 000 kayitta SIFIR baglanma - veritabani '
                         u'okunuyor ama icerigi beklenene uymuyor', yol,
                         u'SILVA vakasindaki sessiz sifir tam olarak boyle '
                         u'gorunuyordu.')
            else:
                saglikli += 1

    rap.olcum[u'M6 veritabani'] = u'%d kaynak denetlendi, %d canli baglanma sinavini gecti' % (
        denetlenen, saglikli)
    if denetlenen == 0:
        rap.atla(M, u'M6-BOS', u'en az bir veritabani denetlenmeli',
                 u'VTB listesinde kullanilan kaynak yok', kay.kimlik_dogrulama)


# Bilinen korunmus bolgeler. Bunlar TASARIM primerlerimiz degildir; genel
# kabul gormus evrensel primer dizileridir ve ilgili lokusta MUTLAKA bulunmalari
# beklenir. Bulunmuyorsa dosya/indeks bozuktur.
_SINAMA_DIZILERI = {
    'SSU': [u'GTGCCAGCMGCCGCGGTAA',    # 515F, SSU'da evrensel
            u'GGATTAGATACCC'],          # 787 bolgesi, cok korunmus
    'LSU': [u'GCATATCAATAAGCGGAGGA',    # LR0R
            u'ACCCGCTGAACTTAAGC'],      # LSU korunmus
    'ITS': [u'GGAAGTAAAAGTCGTAACAAGG',  # ITS1
            u'TCCTCCGCTTATTGATATGC'],   # ITS4
    'OPERON': [u'GTGCCAGCMGCCGCGGTAA', u'GGATTAGATACCC'],
    'KARISIK': [u'GTGCCAGCMGCCGCGGTAA', u'GGAAGTAAAAGTCGTAACAAGG',
                u'GGATTAGATACCC'],
}
_IUPAC = {u'M': u'[AC]', u'R': u'[AG]', u'W': u'[AT]', u'S': u'[GC]', u'Y': u'[CT]',
          u'K': u'[GT]', u'V': u'[ACG]', u'H': u'[ACT]', u'D': u'[AGT]',
          u'B': u'[CGT]', u'N': u'[ACGT]'}


def _dizi_desen(d):
    u"""Dejenere bazlari ACAN duzenli ifade. Dejenere bazi N saymak, bu
    sinamanin yanlislikla sifir donmesine yol acardi (desen 8'in ta kendisi)."""
    return re.compile(u''.join(_IUPAC.get(c, re.escape(c)) for c in d.upper()))


def _baglanma_sinamasi(K, yol, lokus, tavan=20000):
    u"""Bilinen korunmus dizilerden en az biri, ilk `tavan` kayitta geciyor mu?

    Donen: bulunan baglanma sayisi, ya da dosya okunamazsa None.
    Hem duz hem ters tumleyen yonde aranir; tek yon aramak kayitlarin yarisini
    kacirirdi (kume kume yon degisiyor).
    """
    desenler = []
    for d in _SINAMA_DIZILERI.get(lokus, _SINAMA_DIZILERI['KARISIK']):
        desenler.append(_dizi_desen(d))
        desenler.append(_dizi_desen(K.rc(re.sub(r'[^A-Z]', u'', d.upper()))))
    n = 0
    bulunan = 0
    try:
        for _bas, diz in K.fasta_akisi(yol):
            n += 1
            for de in desenler:
                if de.search(diz):
                    bulunan += 1
                    break
            if n >= tavan:
                break
    except (IOError, OSError):
        return None
    return bulunan


# ===========================================================================
# MODUL 7 - TAKSON KAPSAMI
# ===========================================================================
# SORU: toplanti kararlarinda istenen HER hedef panelde bir karsilik buluyor mu?
# Bulmuyorsa sebebi hangi kategoride?
#
# DORT KATEGORI (bunlar birbirinden farkli sonuclar dogurur):
#   ORGANIZMA YOK      - istenen takson numunede hic yok. Primer tasarlanamaz;
#                        bu bir basarisizlik degil, bir bulgudur.
#   AYRIM YOK          - organizma var ama akrabalarindan ayrilamiyor.
#   LOKUS YOK          - organizma var, ayrilabilir, ama elimizdeki lokusta
#                        ayirt edici bolge yok (baska lokus gerekir).
#   UYE KUMESI AYRISIK - hedef tek bir soy degil; uyeler birbirinden uzak.
# Kategorisi belirlenemeyen bir bosluk, KAPATILMAMIS bir bosluktur.
# ---------------------------------------------------------------------------

_KATEGORI = (
    (u'ORGANIZMA YOK', re.compile(
        r'numunede\s*(hic\s*)?yok|t[uü]r\s*numunede\s*yok|cins\s*numunede\s*yok|'
        r'organizma\s*yok|etiket\s*[cç][uü]r[uü]t[uü]ld[uü]', re.I)),
    (u'AYRIM YOK', re.compile(
        r'e[sş]ik\s*alt|ayr[iı]m\s*(kat[iı]\s*)?(yok|tan[iı]ms[iı]z|0[,.]0)|'
        r'ayr[iı]lam[iı]yor|ayr[iı]m\s*e[sş]ik', re.I)),
    (u'LOKUS YOK', re.compile(
        r'lokus|ba[sş]ka\s*b[oö]lge|ay[iı]rt\s*edici\s*b[oö]lge\s*yok|'
        r'ITS.*yetersiz|farkl[iı]\s*lokus', re.I)),
    (u'UYE KUMESI AYRISIK', re.compile(
        r'heterojen|ayr[iı][sş][iı]k|birbirine\s*%?\s*\d+[,.]?\d*\s*[-, ]\s*\d+|'
        r'kutu\s*heterojen|[uü]ye\s*k[uü]mesi\s*ayr', re.I)),
)


def _toplanti_istekleri(metin):
    u"""Toplanti kararlari markdown'indaki tablolardan (istek, durum_metni) cikar.

    Tablolarin ilk sutunu istenen taksonu, son sutunu durumu tasiyor. Basliklar
    ve ayrac satirlari atlanir.
    """
    out = []
    bolum = u''
    for satir in metin.splitlines():
        s = satir.strip()
        if s.startswith(u'## '):
            bolum = s[3:].strip()
            continue
        if not s.startswith(u'|'):
            continue
        hucre = [h.strip() for h in s.strip(u'|').split(u'|')]
        if len(hucre) < 2:
            continue
        if set(u''.join(hucre)) <= set(u'-: '):
            continue                                  # ayrac satiri
        ilk = re.sub(r'[*_`]', u'', hucre[0]).strip()
        if not ilk or ilk.lower().startswith((u'istenen', u'hedef', u'#', u'toplanti')):
            continue                                  # baslik satiri
        out.append((bolum, ilk, u' | '.join(hucre[1:])))
    return out


def modul_7_kapsam(kay, rap):
    M = u'7 TAKSON KAPSAMI'

    toplanti = metin_oku(kay.toplanti)
    if toplanti is None:
        rap.atla(M, u'M7-TOPLANTI-YOK', u'TOPLANTI_KARARLARI_SON_DURUM.md okunabilmeli',
                 u'dosya yok', kay.toplanti)
        return
    istekler = _toplanti_istekleri(toplanti)
    if not istekler:
        rap.atla(M, u'M7-ISTEK-YOK',
                 u'toplanti kararlarindan istenen hedefler cikarilabilmeli',
                 u'markdown icinde tablo satiri bulunamadi', kay.toplanti)
        return

    ciftler = tsv_oku(kay.ciftler)
    panel_hedefleri = set(_ad_norm(r.get(u'hedef') or u'') for r in (ciftler or []))
    if ciftler is None:
        rap.atla(M, u'M7-PANEL-YOK', u'panel cift listesi okunabilmeli',
                 u'ciftler.tsv yok', kay.ciftler)

    # Panel hedeflerinin toplanti kararlarina baglanmasi target_taxon_mapping.py
    # icindeki KARAR sozlugunde tutuluyor. Dosya CALISTIRILMAZ, ast ile
    # ayristirilir - salt okunur denetci bir betigi kosturmaz.
    esleme_karar = {}
    kaynak = metin_oku(kay.takson_esleme)
    if kaynak is None:
        rap.atla(M, u'M7-ESLEME-YOK',
                 u'target_taxon_mapping.py icindeki KARAR sozlugu okunabilmeli',
                 u'dosya yok', kay.takson_esleme)
    else:
        try:
            agac = ast.parse(kaynak)
            for d in agac.body:
                if not isinstance(d, ast.Assign):
                    continue
                adlar = [t.id for t in d.targets if isinstance(t, ast.Name)]
                if u'KARAR' not in adlar:
                    continue
                for anahtar, deger in zip(d.value.keys, d.value.values):
                    try:
                        v = ast.literal_eval(deger)
                        esleme_karar[ast.literal_eval(anahtar)] = v
                    except (ValueError, SyntaxError):
                        continue
        except SyntaxError as e:
            rap.atla(M, u'M7-ESLEME-AYRISTIRILAMADI',
                     u'target_taxon_mapping.py ayristirilabilmeli',
                     u'SyntaxError: %s' % e, kay.takson_esleme)

    karar_hedefleri = set(_ad_norm(v[0]) for v in esleme_karar.values()
                          if isinstance(v, (list, tuple)) and v)

    kapatilmamis = 0
    kategori_sayaci = collections.Counter()
    for bolum, istek, durum in istekler:
        yapilamadi = bool(re.search(
            r'yap[iı]lamad[iı]|sipari[sş]\s*edilmez|panelden\s*[cç][iı]kar', durum, re.I))
        if not yapilamadi:
            continue
        kategoriler = [ad for ad, d in _KATEGORI if d.search(durum)]
        if kategoriler:
            for k in kategoriler:
                kategori_sayaci[k] += 1
        else:
            kapatilmamis += 1
            rap.ekle(M, u'M7-KATEGORISIZ-BOSLUK', CIDDI,
                     u'karsilanmayan her hedefin sebebi dort kategoriden birine '
                     u'girmeli (organizma yok / ayrim yok / lokus yok / uye kumesi ayrisik)',
                     u'"%s" (%s): "%s" - sebep hicbir kategoriye girmiyor'
                     % (istek, bolum, durum.strip()[:200]), kay.toplanti,
                     u'Kategorisiz bir bosluk, raporda "neden olmadi" sorusuna '
                     u'cevap veremez.')

    # --- panelde olup toplanti kararlarina baglanmayan hedefler
    if ciftler is not None and karar_hedefleri:
        baglanmayan = sorted(h for h in panel_hedefleri if h and h not in karar_hedefleri)
        for h in baglanmayan:
            asil = next((r.get(u'hedef') for r in ciftler
                         if _ad_norm(r.get(u'hedef') or u'') == h), h)
            rap.ekle(M, u'M7-PANELDE-FAZLA', UYARI,
                     u'paneldeki her hedef bir toplanti kararina baglanabilmeli',
                     u'"%s" hedefi target_taxon_mapping.py KARAR tablosunda yok'
                     % asil, kay.takson_esleme,
                     u'Rapora "bu neden panelde" diye sorulursa dayanak gosterilemez.')
    elif ciftler is not None and not karar_hedefleri:
        rap.atla(M, u'M7-BAGLANTI',
                 u'paneldeki hedefler toplanti kararlarina baglanmali',
                 u'KARAR tablosu okunamadi, baglanti denetlenemedi',
                 kay.takson_esleme)

    rap.olcum[u'M7 kapsam'] = (
        u'%d toplanti istegi tarandi; kategorilere dagilim: %s; kategorisiz bosluk: %d'
        % (len(istekler),
           u', '.join(u'%s=%d' % (k, v) for k, v in sorted(kategori_sayaci.items()))
           or u'-', kapatilmamis))


# ===========================================================================
# RAPOR YAZIMI
# ===========================================================================
def _tsv_kacis(s):
    return unicode_(s).replace(u'\t', u' ').replace(u'\n', u' ').replace(u'\r', u' ')


def raporla(kay, rap, cikti, kosulan, sureler):
    u"""Tek markdown raporu + makine okunur TSV. Baska hicbir yere yazilmaz."""
    if not os.path.isdir(cikti):
        os.makedirs(cikti)
    damga = time.strftime('%Y-%m-%d_%H%M')
    md_yol = os.path.join(cikti, u'CAPRAZ_KONTROL_%s.md' % damga)
    tsv_yol = os.path.join(cikti, u'CAPRAZ_KONTROL_%s.tsv' % damga)

    sirali = sorted(rap.bulgular,
                    key=lambda b: (_SIRA.get(b.ciddiyet, 9), b.modul, b.kod))

    # ---------------- TSV
    with io.open(tsv_yol, 'w', encoding='utf-8') as fh:
        fh.write(u'\t'.join([u'ciddiyet', u'modul', u'kod', u'ne_bekleniyordu',
                             u'ne_bulundu', u'dosya', u'oneri']) + u'\n')
        for b in sirali:
            fh.write(u'\t'.join(_tsv_kacis(x) for x in
                                [b.ciddiyet, b.modul, b.kod, b.beklenen,
                                 b.bulunan, b.dosya, b.oneri]) + u'\n')

    # ---------------- M1 kimlik tablosu ayri TSV (rapora girecek asil tablo)
    kimlik_yol = None
    satirlar = rap.tablolar.get(u'M1 kimlik tablosu')
    if satirlar:
        kimlik_yol = os.path.join(cikti, u'KIMLIK_VE_KRAKEN_%s.tsv' % damga)
        sut = [(u'kutu', u'kutu'), (u'Kraken_etiketi', u'kraken_etiket'),
               (u'Kraken_guven', u'kraken_guven'), (u'olculen_kimlik', u'olculen'),
               (u'olculen_duzey', u'duzey'), (u'kimlik_yuzde', u'kimlik'),
               (u'hizalanan_bp', u'hiz_uz'), (u'karar_veren_vtb', u'vtb'),
               (u'kayit_no', u'kayit'), (u'TIP_KAYDI', u'tip'),
               (u'ikinci_isabet', u'ikinci'), (u'fark_yuzde_puan', u'fark'),
               (u'ayirt_edilebilir', u'ayrilir'), (u'ayirt_gerekce', u'ayrim_sebep'),
               (u'N_orani', u'n_oran'), (u'saflik', u'saflik'),
               (u'saflik_gerekce', u'saflik_sebep'), (u'sorgu_kaynagi', u'sorgu_kaynagi'),
               (u'uyusan_vtb_sayisi', u'uyusan'), (u'uyusmayanlar', u'uyusmayan'),
               (u'AD_DEGISTI_MI', u'ad_degisti')]
        with io.open(kimlik_yol, 'w', encoding='utf-8') as fh:
            fh.write(u'\t'.join(a for a, _b in sut) + u'\n')
            for r in satirlar:
                fh.write(u'\t'.join(_tsv_kacis(
                    vir(r.get(b), 2) if isinstance(r.get(b), float) else
                    (u'-' if r.get(b) is None else r.get(b))) for _a, b in sut) + u'\n')

    # ---------------- MARKDOWN
    sayim = collections.OrderedDict(
        (c, rap.say(c)) for c in (KRITIK, CIDDI, UYARI, BILGI, ATLANDI))
    kod = rap.cikis_kodu()
    L = []
    L.append(u'# CAPRAZ KONTROL RAPORU')
    L.append(u'')
    L.append(u'PrimerJury qPCR paneli - bagimsiz, salt okunur denetim.')
    L.append(u'')
    L.append(u'| | |')
    L.append(u'|---|---|')
    L.append(u'| Tarih | %s |' % time.strftime('%Y-%m-%d %H:%M'))
    L.append(u'| Betik | cross_check.py %s |' % VERSIYON)
    L.append(u'| Kok klasor | `%s` |' % kay.kok)
    L.append(u'| Kosan moduller | %s |' % (u', '.join(kosulan) or u'-'))
    L.append(u'| Cikis kodu | **%d** |' % kod)
    L.append(u'')
    L.append(u'## Ozet')
    L.append(u'')
    L.append(u'| Ciddiyet | Sayi |')
    L.append(u'|---|---:|')
    for c, n in sayim.items():
        L.append(u'| %s | %d |' % (c, n))
    L.append(u'')
    if sayim[ATLANDI]:
        L.append(u'> **%d kontrol KOSULAMADI.** Bunlar "gecti" DEGILDIR; asagida '
                 u'ATLANDI bolumunde her birinin sebebi yazili. Cikis kodu bu '
                 u'yuzden sifir degil.' % sayim[ATLANDI])
        L.append(u'')

    L.append(u'## Modul durumu')
    L.append(u'')
    L.append(u'| Modul | Durum | Sure | Bulgu (K/C/U/B/A) |')
    L.append(u'|---|---|---|---|')
    for m in kosulan:
        d = rap.modul_durumu.get(m, {})
        L.append(u'| %s | %s | %s | %d/%d/%d/%d/%d |' % (
            m, d.get(u'durum', u'?'),
            sure_metni(sureler.get(m)) if sureler.get(m) is not None else u'olculmedi',
            rap.say(KRITIK, m), rap.say(CIDDI, m), rap.say(UYARI, m),
            rap.say(BILGI, m), rap.say(ATLANDI, m)))
    L.append(u'')

    if rap.olcum:
        L.append(u'## Olculen degerler')
        L.append(u'')
        L.append(u'Bu bolumdeki her sayi bu kosuda OLCULMUSTUR. Olculmeyen hicbir '
                 u'sey buraya yazilmaz.')
        L.append(u'')
        L.append(u'| Olcum | Deger |')
        L.append(u'|---|---|')
        for k, v in rap.olcum.items():
            if k.startswith(u'_'):
                continue          # ic hesaplama degeri, rapora girmez
            L.append(u'| %s | %s |' % (k, v))
        L.append(u'')

    for c in (KRITIK, CIDDI, UYARI, BILGI, ATLANDI):
        grup = [b for b in sirali if b.ciddiyet == c]
        if not grup:
            continue
        L.append(u'## %s (%d)' % (c, len(grup)))
        L.append(u'')
        for b in grup:
            L.append(u'### %s - %s' % (b.kod, b.modul))
            L.append(u'')
            L.append(u'- **Ne bekleniyordu:** %s' % b.beklenen)
            L.append(u'- **Ne bulundu:** %s' % b.bulunan)
            L.append(u'- **Dosya:** `%s`' % b.dosya)
            if b.oneri:
                L.append(u'- **Neden onemli:** %s' % b.oneri)
            L.append(u'')

    if satirlar:
        L.append(u'## Kimlik ve Kraken karsilastirmasi')
        L.append(u'')
        L.append(u'Tam tablo: `%s`' % os.path.basename(kimlik_yol))
        L.append(u'')
        degisen = [r for r in satirlar if r.get(u'ad_degisti') == u'EVET']
        L.append(u'**Adi DEGISEN kutular (%d)** - raporda gosterilecek asil tablo:'
                 % len(degisen))
        L.append(u'')
        if degisen:
            L.append(u'| Kutu | Kraken etiketi | Olculen kimlik | % | Karar veren VTB | '
                     u'Tip kaydi | Ikinci isabet | Fark | Ayirt edilebilir |')
            L.append(u'|---|---|---|---:|---|---|---|---:|---|')
            for r in degisen:
                L.append(u'| %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (
                    r.get(u'kutu'), r.get(u'kraken_etiket'), r.get(u'olculen'),
                    vir(r.get(u'kimlik'), 2), r.get(u'vtb'), r.get(u'tip'),
                    r.get(u'ikinci'), vir(r.get(u'fark'), 2), r.get(u'ayrilir')))
        else:
            L.append(u'_Bu kosuda adi degisen kutu bulunmadi._')
        L.append(u'')

    L.append(u'## Cikis kodu ne demek')
    L.append(u'')
    L.append(u'Bit maskesi, toplanir: 1 = KRITIK bulgu var, 2 = CIDDI bulgu var, '
             u'4 = en az bir kontrol ATLANDI, 8 = betik coktu.')
    L.append(u'Bu kosu: **%d**.' % kod)
    L.append(u'')

    with io.open(md_yol, 'w', encoding='utf-8') as fh:
        fh.write(u'\n'.join(L))
    return md_yol, tsv_yol, kimlik_yol


# ===========================================================================
# KENDINI SINAMA  -  BILEREK BOZUK GIRDI
# ===========================================================================
# "Her zaman gecen bir kontrol aslinda hicbir sey olcmuyordur." Bu bolum her
# module BILEREK BOZULMUS bir girdi verir ve o modulun hatayi GERCEKTEN
# yakalayip yakalamadigini gosterir. Yakalamayan modul, yesil gorunse bile
# ise yaramaz demektir ve sinama BASARISIZ sayilir.
#
# Sinama gecici bir klasorde kurulur; gercek proje dosyalarina DOKUNULMAZ.
# ---------------------------------------------------------------------------

def _sinama_kok_kur(gecici, kay):
    u"""Kucuk, tam ve SAGLAM bir sahte proje kokü kur. Sonra uzerine hata ekilir."""
    os.makedirs(os.path.join(gecici, 'screening'))
    os.makedirs(os.path.join(gecici, 'verification'))
    os.makedirs(os.path.join(gecici, 'REFERANS_DB'))
    os.makedirs(os.path.join(gecici, 'konsensus_kanonik'))
    y = lambda *p: os.path.join(gecici, *p)

    def d(yol, icerik):
        with io.open(yol, 'w', encoding='utf-8') as fh:
            fh.write(icerik)

    # --- panel tanimi: iki hedef, ikisi de saglam
    d(y('screening', 'ciftler.tsv'),
      u'satir\thedef\tsinif\tF\tR\tuye_taksonlar\tolcu_tipi\tuye_kumesi_durumu\n'
      u'2\tHedef_A\tA\tACGTACGTACGTACGTACGT\tTTGCATGCATGCATGCATGC\t111,222\tayrim\tPANELLE_TUTUYOR\n'
      u'3\tHedef_B\tB\tGGGTACGTACGTACGTACGT\tCCGCATGCATGCATGCATGC\t333\tkapsam\tPANELLE_TUTUYOR\n')
    d(y('screening', 'hedef_uyelik.tsv'),
      u'hedef\tuye_taxid\tharic\tkaynak\tnot\n'
      u'Hedef_A\t111,222\t999\tPANEL\t-\n'
      u'Hedef_B\t333\t\tPANEL\t-\n')
    d(y('screening', 'target_taxon_mapping.py'),
      u"KARAR = {\n 2: ('Hedef_A', 'cins', 'Karar 1', 'not'),\n"
      u" 3: ('Hedef_B', 'cins', 'Karar 1', 'not'),\n}\n")
    # --- karar tablolari: tutarli
    d(y('NIHAI_SIPARIS_LISTESI_2026-08-07.tsv'),
      u'sira\tHUKUM\thedef\turun_bp\tF\tR\thukmu_veren_katman\tGEREKCE\t'
      u'yerel_hedef_disi\tkraken_etiketi\tolculen_kimlik\tdCq\n'
      u'1\tSIPARIS\tHedef_A\t150\tACGTACGTACGTACGTACGT\tTTGCATGCATGCATGCATGC\t'
      u'yerel tarama\tolculdu: in siliko hizalama ve jel dogrulamasi gerekli; '
      u'kapsam, filum spektrumu ve organel orani olculdu\t3\tTaxon A\tTaxon A\t5,00\n'
      u'2\tSIPARIS\tHedef_B\t120\tGGGTACGTACGTACGTACGT\tCCGCATGCATGCATGCATGC\t'
      u'MFE katmani\tolculdu: in siliko MFE ve erime egrisi dogrulamasi gerekli; '
      u'kapsam, filum spektrumu ve organel orani olculdu\t2\tTaxon B\tTaxon B\t4,00\n')
    d(y('SIPARIS_LISTESI.tsv'),
      u'sira\thedef\tF\tR\turun_bp\tkraken_etiketi\tolculen_kimlik\tdCq_karsiligi\n'
      u'1\tHedef_A\tACGTACGTACGTACGTACGT\tTTGCATGCATGCATGCATGC\t150\tTaxon A\tTaxon A\t5,00\n'
      u'2\tHedef_B\tGGGTACGTACGTACGTACGT\tCCGCATGCATGCATGCATGC\t120\tTaxon B\tTaxon B\t4,00\n')
    # R=1 -> log2(1)+4,3 = 4,3 ; taban 3,32 -> gerekli 4,30
    d(y('ESIK_VE_OLCUT_2026-08-08.tsv'),
      u'hedef\tESKI_HUKUM\tYENI_HUKUM\tdCq_olculen\tR\tgerekli_dCq\t'
      u'yeni_kural_durum\tuyelik_gecerli_mi\n'
      u'Hedef_A\tGECER\tGECER\t5,00\t1\t4,30\tGECER\tEVET\n'
      u'Hedef_B\tKALIR\tKALIR\t4,00\t1\t4,30\tKALIR\tEVET\n')
    d(y('HEDEF_DISI_AYRINTI_2026-08-07.tsv'),
      u'hedef\ttaxid\tnot\nHedef_A\t777\tbaska organizma\n')
    # --- literatur ve toplanti
    d(y('LITERATUR_2026-08-07.md'),
      u'# Literatur\n\ngerekli dCq >= log2(R) + 4,3 ve taban 3,32 olmali.\n')
    d(y('TOPLANTI_KARARLARI_SON_DURUM.md'),
      u'## Karar 1\n\n| Istenen tur | Durum |\n|---|---|\n'
      u'| *Taxon A* | Var, siparis edilir |\n'
      u'| *Taxon C* | **Yapilamadi.** Tur numunede yok |\n')
    # --- konsensus indeksi + dosyasi
    d(y('konsensus_kanonik', 'A-1_111.kanonik.fa'), u'>x\n%s\n' % (u'ACGT' * 200))
    d(y('konsensus_kanonik', 'INDEKS.tsv'),
      u'kutu\tsinif\tdosya\tkaynak\teski_yon\tcevrildi\tuzunluk\n'
      u'A-1_111\tA\tA-1_111.kanonik.fa\tkons\tSENSE\thayir\t800\n')
    # --- kimlik motorunu kopyala (VTB listesi M6 icin gerekli)
    if os.path.exists(kay.kimlik_dogrulama):
        with io.open(kay.kimlik_dogrulama, encoding='utf-8', errors='replace') as fh:
            d(y('verification', 'identity_verification.py'), fh.read())
    if os.path.exists(kay.tum_kutu):
        with io.open(kay.tum_kutu, encoding='utf-8', errors='replace') as fh:
            d(y('verification', 'all_bin_identities.py'), fh.read())
    return gecici


def _sinama_vtb_dosyalari(gecici, kay, kmer_satiri=u'kvalue=9, kmer_count=262144'):
    u"""VTB listesindeki her dosya icin kucuk ama GERCEKCI bir kopya uret.

    Baglanma sinamasinin gecmesi icin dosyalarda bilinen korunmus diziler
    bulunmali; boylece "saglam" durum gercekten saglam gorunur ve ekilen hata
    ayirt edilebilir.
    """
    K, _h = modul_yukle(kay.kimlik_dogrulama, 'kd_sinama')
    if K is None:
        return []
    govde = (u'GTGCCAGCAGCCGCGGTAA' + u'ACGT' * 60 + u'GGATTAGATACCC' +
             u'ACGT' * 60 + u'GGAAGTAAAAGTCGTAACAAGG' + u'ACGT' * 60 +
             u'GCATATCAATAAGCGGAGGA' + u'ACGT' * 60 + u'TCCTCCGCTTATTGATATGC' +
             u'ACGT' * 60 + u'ACCCGCTGAACTTAAGC' + u'ACGT' * 60)
    yazilan = []
    for etiket, dosya, _lokus, kullan, _n in K.VTB:
        if not kullan:
            continue
        yol = os.path.join(gecici, 'REFERANS_DB', dosya)
        with io.open(yol, 'w', encoding='utf-8') as fh:
            for i in range(3):
                fh.write(u'>RECORD_%d Example organism %d; from TYPE material\n%s\n'
                         % (i, i, govde))
        with io.open(yol + u'.log', 'w', encoding='utf-8') as fh:
            fh.write(u'Binary index v5-single-strand-64 created: 1 MB, %s\n' % kmer_satiri)
        with io.open(yol + u'.primerqc.bin', 'w', encoding='utf-8') as fh:
            fh.write(u'sahte indeks')
        yazilan.append((etiket, yol))
    return yazilan


def kendini_sina(kay, cikti):
    u"""Her module bilerek bozuk girdi ver, hatayi yakalayip yakalamadigini olc.

    Donen: (gecen_modul_sayisi, toplam_modul, ayrinti_listesi)
    """
    import shutil
    import tempfile

    def kur_m1(g):
        # HATA: indekste yazan konsensus dosyasi diskte YOK.
        os.remove(os.path.join(g, 'konsensus_kanonik', 'A-1_111.kanonik.fa'))
        return u'INDEKS.tsv bir konsensus dosyasi gosteriyor ama dosya silindi'

    def kur_m2(g):
        # HATA: ayni hedefin ileri primeri iki dosyada FARKLI.
        yol = os.path.join(g, 'SIPARIS_LISTESI.tsv')
        s = metin_oku(yol).replace(u'ACGTACGTACGTACGTACGT', u'AAAAACGTACGTACGTACGT', 1)
        with io.open(yol, 'w', encoding='utf-8') as fh:
            fh.write(s)
        return u'Hedef_A ileri primeri SIPARIS_LISTESI icinde farkli yazildi'

    def kur_m3(g):
        # HATA: bir ciftin uye kumesi BOSALTILDI ama dCq yerinde duruyor.
        yol = os.path.join(g, 'screening', 'ciftler.tsv')
        s = metin_oku(yol).replace(u'\t111,222\t', u'\t\t', 1)
        with io.open(yol, 'w', encoding='utf-8') as fh:
            fh.write(s)
        return u'Hedef_A uye_taksonlar sutunu bosaltildi, dCq korundu'

    def kur_m4(g):
        # HATA: gerekli_dCq elle yanlis degere cekildi (4,30 -> 2,00).
        yol = os.path.join(g, 'ESIK_VE_OLCUT_2026-08-08.tsv')
        s = metin_oku(yol).replace(u'\t4,30\tGECER', u'\t2,00\tGECER', 1)
        with io.open(yol, 'w', encoding='utf-8') as fh:
            fh.write(s)
        return u'gerekli_dCq 4,30 yerine 2,00 yazildi (log2(R)+4,3 kuralina aykiri)'

    def kur_m5(g):
        # HATA: hedefin KENDI uyesi (111) hedef disi listesine sokuldu.
        yol = os.path.join(g, 'HEDEF_DISI_AYRINTI_2026-08-07.tsv')
        with io.open(yol, 'w', encoding='utf-8') as fh:
            fh.write(u'hedef\ttaxid\tnot\nHedef_A\t111\tits own member was counted as off-target\n')
        return u'Hedef_A kendi uyesi olan 111 taxid\'i hedef disi listesine eklendi'

    def kur_m6(g):
        # HATA: bir veritabaninin k-mer sayisi 3^9 = 19683 (bozuk indeks imzasi).
        K, _h = modul_yukle(kay.kimlik_dogrulama, 'kd_sinama6')
        ilk = next(d for _e, d, _t, k, _n in K.VTB if k)
        yol = os.path.join(g, 'REFERANS_DB', ilk + u'.log')
        with io.open(yol, 'w', encoding='utf-8') as fh:
            fh.write(u'Binary index v5 created: kvalue=9, kmer_count=19683\n')
        return u'%s indeks gunlugu kmer_count=19683 (3^9, bozuk indeks imzasi)' % ilk

    def kur_m7(g):
        # HATA: karsilanmayan bir hedefin sebebi hicbir kategoriye girmiyor.
        yol = os.path.join(g, 'TOPLANTI_KARARLARI_SON_DURUM.md')
        with io.open(yol, 'w', encoding='utf-8') as fh:
            fh.write(u'## Decision 1\n\n| Requested species | Status |\n|---|---|\n| *Taxon C* | **Not achieved.** No reason recorded |\n')
        return u'Karsilanmayan bir hedefin sebebi dort kategoriden hicbirine girmiyor'

    SINAVLAR = [
        (u'1 KIMLIK', kur_m1, u'M1-KONSENSUS-YOK',
         lambda kk, rr: modul_1_kimlik(kk, rr, KontrolNoktasi(u'', etkin=False),
                                       kip=u'hizli')),
        (u'2 IC TUTARLILIK', kur_m2, u'M2-CELISKI', lambda kk, rr: modul_2_ic_tutarlilik(kk, rr)),
        (u'3 UYELIK', kur_m3, u'M3-UYELIK-BOS', lambda kk, rr: modul_3_uyelik(kk, rr)),
        (u'4 LITERATUR', kur_m4, u'M4-ESIK-YANLIS', lambda kk, rr: modul_4_literatur(kk, rr)),
        (u'5 HATA DESENLERI', kur_m5, u'M5-D2-UYE-HEDEF-DISI',
         lambda kk, rr: modul_5_desenler(kk, rr)),
        (u'6 VERITABANI', kur_m6, u'M6-KMER-BOZUK',
         lambda kk, rr: modul_6_veritabani(kk, rr, baglanma_sinamasi=True)),
        (u'7 TAKSON KAPSAMI', kur_m7, u'M7-KATEGORISIZ-BOSLUK',
         lambda kk, rr: modul_7_kapsam(kk, rr)),
    ]

    yaz(u'')
    yaz(u'=== SELF-TEST: each module is given DELIBERATELY BROKEN input ===')
    ayrinti = []
    gecen = 0
    for ad, kur, beklenen_kod, kos in SINAVLAR:
        gecici = tempfile.mkdtemp(prefix='capraz_sinama_')
        try:
            _sinama_kok_kur(gecici, kay)
            _sinama_vtb_dosyalari(gecici, kay)
            ekilen = kur(gecici)
            kk = Kaynaklar(gecici)
            rr = Rapor()
            t0 = time.time()
            try:
                kos(kk, rr)
                cokme = None
            except Exception as e:
                cokme = u'%s: %s' % (type(e).__name__, e)
            kodlar = set(b.kod for b in rr.bulgular)
            yakalandi = beklenen_kod in kodlar
            if yakalandi:
                gecen += 1
            durum = u'YAKALADI' if yakalandi else (
                u'COKTU (%s)' % cokme if cokme else u'YAKALAYAMADI')
            yaz(u'  [%s] %-22s  planted error: %s' % (
                u'OK ' if yakalandi else u'ERROR', ad, ekilen))
            yaz(u'        beklenen bulgu kodu: %s  ->  %s  (%s)'
                % (beklenen_kod, durum, sure_metni(time.time() - t0)))
            ayrinti.append(dict(modul=ad, ekilen=ekilen, beklenen=beklenen_kod,
                                durum=durum, yakalandi=yakalandi,
                                uretilen_kod_sayisi=len(kodlar)))
        finally:
            shutil.rmtree(gecici, ignore_errors=True)

    yaz(u'')
    yaz(u'=== SELF-TEST RESULT: %d/%d modules caught the planted error ==='
        % (gecen, len(SINAVLAR)))
    if not os.path.isdir(cikti):
        os.makedirs(cikti)
    yol = os.path.join(cikti, u'KENDINI_SINAMA_%s.tsv' % time.strftime('%Y-%m-%d_%H%M'))
    with io.open(yol, 'w', encoding='utf-8') as fh:
        fh.write(u'modul\tekilen_hata\tbeklenen_bulgu_kodu\tsonuc\n')
        for a in ayrinti:
            fh.write(u'%s\t%s\t%s\t%s\n' % (a['modul'], _tsv_kacis(a['ekilen']),
                                            a['beklenen'], a['durum']))
    yaz(u'Self-test detail: %s' % yol)
    return gecen, len(SINAVLAR), ayrinti


# ===========================================================================
# SURUCU
# ===========================================================================
MODUL_ADLARI = collections.OrderedDict([
    (u'1', u'1 KIMLIK'), (u'2', u'2 IC TUTARLILIK'), (u'3', u'3 UYELIK'),
    (u'4', u'4 LITERATUR'), (u'5', u'5 HATA DESENLERI'), (u'6', u'6 VERITABANI'),
    (u'7', u'7 TAKSON KAPSAMI'),
])



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
    p = argparse.ArgumentParser(
        description=u'PrimerJury paneli - bagimsiz, salt okunur capraz kontrol')
    p.add_argument('--root', '--kok', dest='kok', default='.', help=u'project root directory')
    p.add_argument('--output', '--cikti', dest='cikti', default=None, help=u'report directory (default KONTROL_SONUC)')
    p.add_argument('--modules', '--moduller', dest='moduller', default='hepsi',
                   help=u'modules to run, comma-separated: 1,2,3 ... or "all"')
    p.add_argument('--m1-mode', '--m1-kip', dest='m1_kip', default='hizli',
                   choices=['none', 'quick', 'full', 'yok', 'hizli', 'tam'],
                   help=u'M1 scope of the identity scan (default: quick)')
    p.add_argument('--m1-only', '--m1-yalniz', dest='m1_yalniz', default=None,
                   help=u'these bins only (comma-separated), example: F2-1_500148,F2-4_500148')
    p.add_argument('--m1-cap', '--m1-tavan', dest='m1_tavan', type=int, default=0,
                   help=u'M1 for maximum number of bins (0 = all)')
    p.add_argument('--reset', '--sifirla', dest='sifirla', action='store_true',
                   help=u'ignore checkpoints and recompute everything')
    p.add_argument('--no-checkpoint', '--kontrol-noktasi-yok', dest='kn_yok', action='store_true',
                   help=u'never use a checkpoint')
    p.add_argument('--self-test', '--kendini-sina', dest='kendini_sina', action='store_true',
                   help=u'feed each module deliberately broken input and show that it catches the error')
    a = p.parse_args()
    a = _ing_deger(a)

    kok = os.path.abspath(a.kok)
    cikti = a.cikti or os.path.join(kok, u'KONTROL_SONUC')
    kay = Kaynaklar(kok)

    yaz(u'CROSS-CHECK %s' % VERSIYON)
    yaz(u'Root directory : %s' % kok)
    yaz(u'Output         : %s' % cikti)
    yaz(u'READ-ONLY: this script writes to nothing outside the output directory.')
    if not os.path.isdir(kok):
        yaz(u'ERROR: root directory does not exist: %s' % kok)
        return 8

    if a.kendini_sina:
        gecen, toplam, _ay = kendini_sina(kay, cikti)
        # Bir modul ekilen hatayi yakalayamazsa sinama BASARISIZDIR.
        return 0 if gecen == toplam else 8

    if a.moduller.strip().lower() in (u'hepsi', u'all', u'*'):
        secili = list(MODUL_ADLARI.keys())
    else:
        secili = [x.strip() for x in a.moduller.split(u',') if x.strip()]
        bilinmeyen = [x for x in secili if x not in MODUL_ADLARI]
        if bilinmeyen:
            yaz(u'ERROR: unknown module: %s (valid: %s)'
                % (u', '.join(bilinmeyen), u', '.join(MODUL_ADLARI)))
            return 8

    kn_klasor = os.path.join(cikti, u'kontrol')
    if a.sifirla and os.path.isdir(kn_klasor):
        import shutil
        shutil.rmtree(kn_klasor, ignore_errors=True)
        yaz(u'Checkpoints deleted (--reset).')
    kn = KontrolNoktasi(kn_klasor, etkin=not a.kn_yok)

    rap = Rapor()
    kosulan = []
    sureler = {}
    t_hepsi = time.time()

    ISLER = {
        u'1': lambda: modul_1_kimlik(kay, rap, kn, kip=a.m1_kip,
                                     yalniz=a.m1_yalniz, tavan=a.m1_tavan),
        u'2': lambda: modul_2_ic_tutarlilik(kay, rap),
        u'3': lambda: modul_3_uyelik(kay, rap),
        u'4': lambda: modul_4_literatur(kay, rap),
        u'5': lambda: modul_5_desenler(kay, rap),
        u'6': lambda: modul_6_veritabani(kay, rap),
        u'7': lambda: modul_7_kapsam(kay, rap),
    }

    for no in secili:
        ad = MODUL_ADLARI[no]
        kosulan.append(ad)
        yaz(u'')
        yaz(u'===== MODUL %s =====' % ad)
        t0 = time.time()
        try:
            ISLER[no]()
            rap.modul_durumu[ad] = dict(durum=u'kostu')
        except KeyboardInterrupt:
            # Kesinti YUTULMAZ: kontrol noktalari yazildigi icin kaldigi yerden
            # devam edilebilir, ama kosu YARIM oldugu acikca soylenir.
            yaz(u'KESILDI (Ctrl+C). Kontrol noktalari korundu, yeniden kosun.')
            rap.modul_durumu[ad] = dict(durum=u'KESILDI')
            rap.atla(ad, u'KESINTI', u'modul bastan sona kosmali',
                     u'kullanici kesti (Ctrl+C)', u'-')
            sureler[ad] = time.time() - t0
            break
        except Exception as e:
            # HATA MASKELENMEZ: modul cokerse bu bir ATLANDI bulgusudur ve
            # cikis kodunda gorunur.
            iz = traceback.format_exc()
            yaz(u'MODUL COKTU: %s: %s' % (type(e).__name__, e))
            rap.modul_durumu[ad] = dict(durum=u'COKTU')
            rap.atla(ad, u'MODUL-COKTU', u'modul hatasiz kosmali',
                     u'%s: %s | %s' % (type(e).__name__, e,
                                       iz.strip().splitlines()[-1]), u'-')
        sureler[ad] = time.time() - t0
        yaz(u'--- %s done: %s, findings C%d/S%d/W%d/I%d/K%d' % (
            ad, sure_metni(sureler[ad]), rap.say(KRITIK, ad), rap.say(CIDDI, ad),
            rap.say(UYARI, ad), rap.say(BILGI, ad), rap.say(ATLANDI, ad)))

    rap.olcum[u'Toplam kosu suresi'] = sure_metni(time.time() - t_hepsi)
    rap.olcum[u'Kontrol noktasi'] = u'%d isabet, %d iska' % (kn.isabet, kn.iska)
    for no in MODUL_ADLARI:
        if MODUL_ADLARI[no] not in kosulan:
            rap.atla(MODUL_ADLARI[no], u'MODUL-KOSULMADI',
                     u'butun moduller kosmali',
                     u'--moduller secimiyle disarida birakildi', u'-')

    md, tsv, kimlik = raporla(kay, rap, cikti, kosulan, sureler)
    kod = rap.cikis_kodu()
    yaz(u'')
    yaz(u'================= SUMMARY =================')
    for c in (KRITIK, CIDDI, UYARI, BILGI, ATLANDI):
        yaz(u'  %-8s : %d' % (c, rap.say(c)))
    yaz(u'Report         : %s' % md)
    yaz(u'Machine TSV    : %s' % tsv)
    if kimlik:
        yaz(u'Identity TSV   : %s' % kimlik)
    yaz(u'Exit code      : %d  (1=CRITICAL, 2=SERIOUS, 4=SKIPPED present; these are summed)' % kod)
    return kod


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        yaz(u'INTERRUPTED.')
        sys.exit(8)
