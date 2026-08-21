# -*- coding: utf-8 -*-
"""
minimap2_aligner.py - kimlik asamasi icin IKINCI, SECILEBILIR hizalayici.

Bu dosya mevcut motoru DEGISTIRMEZ, yanina durur. identity_verification.py icindeki
saf Python + numpy hizalayicisi (hizala) yerinde kalir ve varsayilan olmaya
devam eder. Buradaki minimap2 yolu ancak iki sart birden saglanirsa devreye
girer: mappy kurulu olacak VE kullanici acikca secmis olacak.
"""
# ---------------------------------------------------------------------------
# minimap2_aligner.py
#
# GIRDI  : sorgu dizisi (kutu konsensusu) ve hedef kayitlar (veritabani kayitlari)
# CIKTI  : (yuzde_kimlik, uzaklik) ciftleri - kimlik_dogrulama.hizala ile AYNI
#          bicimde, boylece yer degistirebilirler
# CAGRAN : verification/identity_verification.py (I tusu) ve verification/all_bin_identities.py
#          (G tusu), yalniz HIZALAYICI=minimap2 verildiginde
#
# NEDEN VAR
# En yavas adimimiz kimlik asamasidir: bir kutu konsensusu on iki veritabanina
# karsi taranir ve her veritabaninda 500 aday TAM hizalanir. Saf Python yolu
# her cift icin O(len(q) x len(t)) dinamik programlama yapar; numpy ile
# vektorlestirilmis olsa bile is yine kuadratiktir. minimap2 tam bu is icin
# tasarlandi: uzun ve hatali okumalari buyuk referans kumelerine tohum-zincir-uzat
# ile hizalar, kuadratik DP'yi yalniz zincirin etrafindaki dar bantta yapar.
#
# NEREDE KULLANILMAZ - bu sinir onemlidir
#   * PRIMER BAGLANMA ARAMASI. Primerler 18-25 bazdir. minimap2 varsayilan
#     tohumu (k=15, w=10) ve zincirleme mantigi bu boydaki sorgular icin
#     tasarlanmadi; kisa sorgularda tohum bulamayip baglanma yerini SESSIZCE
#     kacirir. Orada guvercin yuvasi motoru (screening/read_engine.py)
#     kalir ve kalmalidir, cunku o motorun kayipsizligi bir GARANTIDIR.
#   * IN-SILIKO PCR URUN HESABI. Ayni sebep: iki primerin karsilikli ve dogru
#     yonde baglandigi yerleri bulmak, kisa ve TAM eslesmeye yakin arama
#     isidir. minimap2'nin yaklasik cevabi burada kabul edilemez.
#
# KURULUM
#     pip install mappy
#   ya da proje ortamina:
#     micromamba install -n mikro -c bioconda minimap2
#     micromamba run -n mikro pip install mappy
#
# GUVENLIK: mappy yoksa bu dosya HICBIR SEY BOZMAZ. var_mi() False doner,
# cagiran taraf mevcut motorla devam eder ve zincir kirilmaz.
# ---------------------------------------------------------------------------

from __future__ import print_function

_MAPPY = None
_DENENDI = False
_SEBEP = u''


def var_mi():
    """mappy kurulu ve calisir durumda mi. Tek sefer olcer, sonra onbellekten.

    'import edilebiliyor' yetmez: gercekten bir indeks kurulup kurulamadigi da
    denenir, cunku bozuk bir kurulum import asamasini gecip ilk kullanimda
    patlayabilir ve o an saatler suren bir kosunun ortasidir.
    """
    global _MAPPY, _DENENDI, _SEBEP
    if _DENENDI:
        return _MAPPY is not None
    _DENENDI = True
    try:
        import mappy
    except ImportError:
        _SEBEP = u'mappy kurulu degil (pip install mappy)'
        return False
    try:
        deneme = mappy.Aligner(seq='ACGT' * 40, preset='map-ont', n_threads=1)
        if not deneme:
            _SEBEP = u'mappy indeksi bos dondu'
            return False
        list(deneme.map('ACGT' * 40))
    except Exception as e:
        _SEBEP = u'mappy kurulu ama calismadi: %s: %s' % (type(e).__name__, e)
        return False
    _MAPPY = mappy
    return True


def sebep():
    return _SEBEP


def surum():
    if not var_mi():
        return u'yok'
    return getattr(_MAPPY, '__version__', u'bilinmiyor')


# ---------------------------------------------------------------------------
# KIMLIK YUZDESI NASIL HESAPLANIR
#
# Mevcut motor (kimlik_dogrulama.hizala) sunu doner:
#     yuzde = 100 * (1 - duzenleme_uzakligi / len(sorgu))
# yani sorgunun TAMAMI hedefin icine oturtulur (infix/HW hizalama) ve maliyet
# sorgu uzunluguna bolunur.
#
# minimap2 ise yerel (local) hizalama yapar: sorgunun yalniz hizalanan parcasini
# bildirir. Bu iki sayi AYNI SEY DEGILDIR ve dogrudan karsilastirilamaz.
# Karsilastirilabilir kilmak icin yerel sonuc sorgu uzunluguna gore yeniden
# olceklenir: hizalanmayan her baz bir uyumsuzluk sayilir.
#
#     hizalanan_dogru = blen - NM
#     hizalanmayan    = len(q) - (q_en - q_st)
#     uzaklik         = NM + hizalanmayan
#     yuzde           = 100 * (1 - uzaklik / len(q))
#
# Boylece iki motor ayni tanimi olcer ve karsilastirma anlamli olur. Bu
# donusum KARSILASTIRMANIN GECERLILIGI icin sarttir; atlanirsa minimap2
# sistematik olarak daha yuksek yuzde verir ve "uyusuyorlar" sanilir.
# ---------------------------------------------------------------------------
def hizala_mm(q, t):
    """kimlik_dogrulama.hizala ile AYNI imza ve AYNI donus bicimi.

    Doner: (yuzde_kimlik, uzaklik). Hizalama bulunamazsa (0.0, len(q)).
    """
    if not var_mi():
        raise RuntimeError(u'mappy yok: %s' % _SEBEP)
    if not q or not t:
        return (0.0, len(q or t or ' '))
    try:
        ind = _MAPPY.Aligner(seq=t, preset='map-ont', n_threads=1)
        if not ind:
            return (0.0, len(q))
        en_iyi = None
        for h in ind.map(q):
            if en_iyi is None or h.mlen > en_iyi.mlen:
                en_iyi = h
        if en_iyi is None:
            return (0.0, len(q))
        hizalanmayan = len(q) - (en_iyi.q_en - en_iyi.q_st)
        uzaklik = int(en_iyi.NM) + max(hizalanmayan, 0)
        yuzde = round(100.0 * (1 - uzaklik / float(len(q))), 2)
        return (max(yuzde, 0.0), uzaklik)
    except Exception:
        # Tek bir kaydin hizalanamamasi butun kosuyu dusurmemeli.
        return (0.0, len(q))


def toplu_hizala(q, hedefler, iplik=3):
    """ASIL HIZ KAZANCI BURADADIR: TEK indeks, TEK haritalama.

    2026-08-05 DUZELTMESI - SESSIZ YOL DEGISIMI KAPATILDI
    Ilk surum mappy.Aligner(seq=<liste>) cagiriyordu. mappy bunu KABUL ETMIYOR
    ve TypeError atiyor; kod da onu yakalayip sessizce hedef basina tek tek
    indeks kuran YAVAS yola dusuyordu. Olcumler yine dogruydu ama "toplu indeks"
    iddiasi gercek DEGILDI ve yorum bunun tersini soyluyordu. Tam olarak bu
    projenin kovaladigi hata turu: kod hata vermeden BASKA bir sey yapiyor.

    Dogru yol: diziler gecici bir FASTA'ya yazilir ve indeks dosyadan kurulur.
    Boylece h.ctg GERCEK kayit adini tasir, isim eslesmesi tahmine kalmaz.

    hedefler: [(anahtar, dizi), ...]
    Doner   : {anahtar: (yuzde, uzaklik)}   - hizalanmayanlar (0.0, len(q))
    """
    if not var_mi():
        raise RuntimeError(u'mappy yok: %s' % _SEBEP)
    sonuc = dict((a, (0.0, len(q))) for a, _d in hedefler)
    if not q or not hedefler:
        return sonuc

    # Anahtarlar FASTA basligi olacagi icin bosluk ve satir sonu tasiyamaz.
    # Gecici ad -> gercek anahtar eslemesi ayrica tutulur; boylece ayni
    # basligi tasiyan iki kayit birbirine karismaz.
    import os as _os
    import tempfile as _tf
    esleme = {}
    fd, yol = _tf.mkstemp(suffix='.fa', prefix='mm2_')
    try:
        with _os.fdopen(fd, 'w') as fh:
            for i, (a, d) in enumerate(hedefler):
                ad = 's%d' % i
                esleme[ad] = a
                fh.write('>%s\n%s\n' % (ad, d))
        # best_n YUKSEK TUTULUR. Varsayilan minimap2 yalnizca BIRINCIL isabeti
        # bildirir; toplu indekste bu, "en iyi kayit disindaki her sey gorunmez"
        # demek olurdu ve tam da kaldirmaya calistigimiz KESME NOKTASINI geri
        # getirirdi. Aday bulucu olarak kullanacaksak butun makul isabetler
        # gelmelidir; eleme kararini Python puanlayicisi verir, minimap2 degil.
        ind = _MAPPY.Aligner(fn_idx_in=yol, preset='map-ont', n_threads=iplik,
                             best_n=max(len(hedefler), 50))
        if not ind:
            return sonuc
        en_iyi = {}
        for h in ind.map(q):
            ad = h.ctg
            if ad not in en_iyi or h.mlen > en_iyi[ad].mlen:
                en_iyi[ad] = h
        for ad, h in en_iyi.items():
            anahtar = esleme.get(ad)
            if anahtar is None:
                continue
            hizalanmayan = len(q) - (h.q_en - h.q_st)
            uzaklik = int(h.NM) + max(hizalanmayan, 0)
            yuzde = round(100.0 * (1 - uzaklik / float(len(q))), 2)
            sonuc[anahtar] = (max(yuzde, 0.0), uzaklik)
        return sonuc
    finally:
        try:
            _os.unlink(yol)
        except OSError:
            pass



def secili_mi():
    """Kullanici minimap2'yi ACIKCA sectiyse True.

    Iki sart birden: ortam degiskeni HIZALAYICI=minimap2 VE mappy calisiyor.
    Varsayilan olarak KAPALIDIR. Sebep: karsilastirma raporu (bkz.
    MINIMAP2_KARSILASTIRMA.md) iki motorun ayni sonucu verdigini gosterene
    kadar hizli olani dogru saymayiz. Proje kurali 1: hicbir karar tek bir kod
    yoluna birakilmaz ve iki olcum ayrilirsa hizli olan degil, ELLE DOGRULANAN
    kazanir.
    """
    import os
    if os.environ.get('HIZALAYICI', '').strip().lower() != 'minimap2':
        return False
    return var_mi()


if __name__ == '__main__':
    import sys
    print(u'mappy calisiyor mu : %s' % (u'EVET' if var_mi() else u'HAYIR'))
    if not var_mi():
        print(u'sebep              : %s' % _SEBEP)
        print(u'kurulum            : pip install mappy')
        sys.exit(1)
    print(u'mappy surumu       : %s' % surum())
    print(u'HIZALAYICI secili  : %s' % (u'minimap2' if secili_mi() else u'python (varsayilan)'))


# ---------------------------------------------------------------------------
# HIBRIT YOL - asil onerilen kullanim (2026-08-05 olcumu)
#
# NEDEN TAM DEGISIM DEGIL
# Olcum sunu gosterdi: minimap2 ile saf Python motoru AYNI SORUYU SORMUYOR.
#   * Python'un hizala() fonksiyonu infix (HW) hizalamadir: sorguyu hedefin
#     ICINE ZORLA oturtur ve HER hedef icin bir sayi uretir. Alakasiz bir
#     veritabaninda bile yuzde 50-65 arasi degerler dondurur. Bu degerler
#     korunmus bolgelerden ve zorlamali hizalamadan gelen GURULTUDUR.
#   * minimap2 yerel hizalamadir: gercek bir homoloji yoksa HIC hizalama
#     bildirmez. Yani alakasiz veritabaninda 0 aday doner.
#
# Olculen ornek (A1-1_1826872, veritabani basina 150 kayit):
#     SILVA SSU NR99 (dogru lokus) : minimap2 2 aday buldu, hibritin sectigi
#                                    en iyi isabet Python'un tam taramasiyla
#                                    AYNI cikti (%73,2)
#     SILVA LSU NR99 / LSU Parc / UNITE ITS (yanlis lokus):
#                                    minimap2 0 aday, Python %60,4 / %52,1 /
#                                    %51,3 uretti - hicbiri gercek eslesme degil
#
# Yani "ayrilik" minimap2'nin yanilmasi degil, COPU BILDIRMEYI REDDETMESIDIR.
# Zaten bu degerler hukum esiklerinin (tur %98,7, cins %90) cok altinda kalir
# ve hicbir karara girmez.
#
# DOGRU ENTEGRASYON BU YUZDEN HIBRITTIR:
#     minimap2 ADAYLARI BULUR  ->  Python o adaylari PUANLAR
# Kimlik yuzdesi tanimi degismez (mevcut motorun tanimi korunur), ama her
# veritabaninda binlerce kaydi tam hizalamak yerine yalnizca gercekten
# hizalanan bir avuc kayit puanlanir.
#
# EK KAZANC: KISA LISTE SORUNU
# Mevcut yolda adaylar TOHUM SAYISINA gore siralanip ilk 500'u aliniyordu ve
# olcum kesme noktasinin BAGLAYICI oldugunu gosterdi (kazanan bir sorguda
# 4171. siradan geldi, 118 sorgunun 13'unde 400'un otesinden). minimap2
# minimizer zincirlemesiyle arar; "ilk N" diye bir kesme yoktur, hizalanan
# her kayit gelir. Yani hibrit yol kisa liste kesmesini TUMDEN ORTADAN
# KALDIRIR. Bu, hiz kazancindan daha degerli olabilir.
# ---------------------------------------------------------------------------
def hibrit_adaylar(q, hedefler):
    """minimap2 ile GERCEKTEN hizalanan hedefleri secer.

    Doner: [(anahtar, dizi), ...] - yalnizca hizalama bulunanlar.
    mappy yoksa bos liste yerine None doner; cagiran taraf o zaman mevcut
    kisa liste yolunu kullanmalidir. None ile bos liste AYRI SEYLERDIR:
      None     -> minimap2 yok, karar verilemedi, eski yola dus
      []       -> minimap2 var ve bu veritabaninda hicbir sey hizalanmadi
    Bu ayrim onemlidir; karistirilirsa bir veritabani sessizce atlanir.
    """
    if not var_mi():
        return None
    # OLCULDU 2026-08-05: toplu_hizala TEK indeks kurar ve minimap2 o kipte
    # IKINCIL hizalamalari ELER. Dort hedefli sinamada gercek %75,76'lik iki
    # homolog "0" gorundu. Aday bulucu olarak kullanilamaz: kaldirmaya
    # calistigimiz kisa liste kesmesinin yerine daha OPAK bir kesme koyar.
    # Bu yuzden burada hedef basina AYRI indeks kurulur - yavas ama kayipsiz.
    return [(a, d) for a, d in hedefler if hizala_mm(q, d)[0] > 0]
