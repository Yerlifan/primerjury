# -*- coding: utf-8 -*-
"""Referans kapsami ve rakip ayrimi: SILVA/UNITE havuzlarindan cins/aile
duzeyinde dizi kumeleri cikarilir, ispcr.find_sites ile OLCULUR.

Havuz cikarma pahalidir; sonuc CIKTI/onbellek altina yazilir ve bir daha
cikarilmaz. Kesinti olursa onbellekten devam eder.
"""
# ---------------------------------------------------------------------------
# referans.py — uye ve rakip cinsler icin referans veritabanlarindan dizi
#               havuzlari cikarir ve bir cift bu havuzlarda kac kayitta urun
#               veriyor olcer.
#
# GIRDI  : REFERANS_DB altindaki SILVA SSU / LSU, UNITE ITS fasta dosyalari
#          (amplikon sinifina gore secilir); uye ve rakip cins adlari
#          hedefler.taxid_adlari() + hedef_baglami ciktisindan
#          uye_ve_rakip_anahtar() ile turetilir.
# CIKTI  : KAPSAMLI_ARAMA_SONUC/onbellek/havuz_*.pkl onbellek dosyalari.
#          havuz_cikar() (baslik, dizi) listesi, kapsam() veren/toplam/boy/
#          basliklar sozlugu dondurur.
# CAGRAN : __main__.aramayi_kos (toplu havuz cikarimi) ve
#          __main__.hedefi_isle asama D - yani tuslar 1, 2, 3, 7 ve 9'un
#          7. asamasi. --hafif verilirse bu adim atlanir.
#
# toplu_cikar() NEDEN VAR: her hedef icin ayri gecis yapilsaydi ayni dev fasta
# dosyasi hedef sayisi kadar bastan okunurdu. Burada istekler veritabanina gore
# gruplanir ve her veritabani en fazla BIR kez okunur; 21 hedef icin 21 gecis
# yerine veritabani sayisi kadar gecis olur.
# ---------------------------------------------------------------------------
import os, re, pickle
from . import yapilandirma as C
from . import motor

MIN_UZ = 1200          # 16S icin; ITS havuzlarinda 300'e duser


def _uygun_db(sinif):
    """Amplikon sinifina gore hangi referans veritabani."""
    if sinif in ('F1', 'F2'):
        for p in (C.UNITE_ITS, C.SILVA_LSU):
            if os.path.exists(p):
                return p, 300
    return C.SILVA_SSU, MIN_UZ


def havuz_cikar(anahtarlar, sinif, onbellek_adi, ilerle=None):
    """Taksonomi satirinda anahtarlardan biri gecen dizileri topla."""
    os.makedirs(C.ONBELLEK, exist_ok=True)
    cache = _cache_yolu(onbellek_adi)
    if os.path.exists(cache):
        return pickle.load(open(cache, 'rb'))
    db, minuz = _uygun_db(sinif)
    if not os.path.exists(db):
        return []
    anah = [a.lower() for a in anahtarlar if a]
    out = []
    n = 0
    for ad, seq in motor.read_fasta(db):
        n += 1
        if ilerle and n % 200000 == 0:
            ilerle(n, len(out))
        if len(seq) < minuz:
            continue
        low = ad.lower()
        if any(a in low for a in anah):
            out.append((ad[:150], motor.clean(seq.upper())))
    pickle.dump(out, open(cache, 'wb'))
    return out


def _cache_yolu(ad):
    return os.path.join(C.ONBELLEK, 'havuz_%s.pkl' % re.sub(r'\W+', '_', ad))


def toplu_cikar(istekler, ilerle=None):
    """BUTUN havuzlari TEK gecisde cikar.

    istekler: [(onbellek_adi, [anahtar,...], sinif)]
    Her veritabani en fazla BIR kez okunur. 21 hedef icin 21 gecis yerine
    veritabani sayisi kadar gecis olur - saatlerce fark eder.
    """
    os.makedirs(C.ONBELLEK, exist_ok=True)
    kalan = [(ad, [a.lower() for a in anah if a], sinif)
             for ad, anah, sinif in istekler
             if anah and not os.path.exists(_cache_yolu(ad))]
    if not kalan:
        return
    # veritabanina gore grupla
    gruplar = {}
    for ad, anah, sinif in kalan:
        db, minuz = _uygun_db(sinif)
        gruplar.setdefault((db, minuz), []).append((ad, anah))
    for (db, minuz), isler in gruplar.items():
        if not os.path.exists(db):
            for ad, _ in isler:
                pickle.dump([], open(_cache_yolu(ad), 'wb'))
            continue
        kutular = {ad: [] for ad, _ in isler}
        n = 0
        for baslik, seq in motor.read_fasta(db):
            n += 1
            if ilerle and n % 100000 == 0:
                ilerle(os.path.basename(db), n, sum(len(v) for v in kutular.values()))
            if len(seq) < minuz:
                continue
            low = baslik.lower()
            temiz = None
            for ad, anah in isler:
                if any(a in low for a in anah):
                    if temiz is None:
                        temiz = motor.clean(seq.upper())
                    kutular[ad].append((baslik[:150], temiz))
        for ad, _ in isler:
            pickle.dump(kutular[ad], open(_cache_yolu(ad), 'wb'))


def kapsam(havuz, F, R, lo, hi, max_mm=C.REFERANS_MAX_MM):
    """Havuzdaki kac kayit urun veriyor + boy dagilimi + vurus basliklari."""
    veren = 0
    boy = {}
    basliklar = []
    for ad, seq in havuz:
        p = motor.amplify(seq, F, R, max_mm=max_mm, lo=lo, hi=hi)
        if p:
            veren += 1
            b = min(p, key=lambda x: x[3] + x[4])[2]
            boy[b] = boy.get(b, 0) + 1
            if len(basliklar) < 40:
                basliklar.append((ad, b))
    return dict(veren=veren, toplam=len(havuz), boy=boy, basliklar=basliklar)


def uye_ve_rakip_anahtar(baglam, taxad):
    """Uye taxid adlarindan cins anahtarlari; rakip = ayni sinifin geri kalani."""
    uye_adlar = [taxad.get(t, '') for t in baglam['uye_tax']]
    cinsler = sorted({a.split()[0] for a in uye_adlar if a and not a.startswith('Ca.')}
                     | {' '.join(a.split()[:2]) for a in uye_adlar if a.startswith('Ca.')})
    rakip_adlar = [taxad.get(k['taxid'], '') for k in baglam['rakip_kutu']]
    rakip_cins = sorted({a.split()[0] for a in rakip_adlar if a and a.split()[0] not in cinsler})
    return cinsler, rakip_cins
