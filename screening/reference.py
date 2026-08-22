# -*- coding: utf-8 -*-
"""Reference coverage and competitor discrimination: sequence sets are extracted at
genus or family level from the SILVA and UNITE pools and MEASURED with
ispcr.find_sites.

Extracting a pool is expensive; the result is written under CIKTI/onbellek and not
extracted again. After an interruption it continues from the cache.

"""
# -------------------------------------------------------------------------
# reference.py extracts sequence pools from the reference databases for the member
#               and competitor genera, and measures in how many records of those
#               pools a pair gives a product.
#
# INPUT  : the SILVA SSU / LSU and UNITE ITS fasta files under REFERANS_DB (chosen
#          by the amplicon class); the member and competitor genus names are derived
#          from hedefler.taxid_adlari() plus the hedef_baglami output through
#          uye_ve_rakip_anahtar().
# OUTPUT : the KAPSAMLI_ARAMA_SONUC/onbellek/havuz_*.pkl cache files.
#          havuz_cikar() returns a (header, sequence) list and kapsam() a dictionary
#          of giving/total/length/headers.
# CALLED BY: __main__.aramayi_kos (the bulk pool extraction) and stage D of
#          __main__.hedefi_isle, that is keys 1, 2, 3, 7 and the 7th stage of 9.
#          The step is skipped if --light is given.
#
# WHY toplu_cikar() EXISTS: had a separate pass been made for each target, the same
# huge fasta file would have been read from the start once per target. Here the
# requests are grouped by database and each database is read at most ONCE; for 21
# targets that is one pass per database instead of 21.
# -------------------------------------------------------------------------
import os, re, pickle
from . import config as C
from . import engine_gateway

MIN_UZ = 1200          # for 16S; it drops to 300 in the ITS pools


def _uygun_db(sinif):
    """Which reference database, by amplicon class."""
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
    for ad, seq in engine_gateway.read_fasta(db):
        n += 1
        if ilerle and n % 200000 == 0:
            ilerle(n, len(out))
        if len(seq) < minuz:
            continue
        low = ad.lower()
        if any(a in low for a in anah):
            out.append((ad[:150], engine_gateway.clean(seq.upper())))
    pickle.dump(out, open(cache, 'wb'))
    return out


def _cache_yolu(ad):
    return os.path.join(C.ONBELLEK, 'havuz_%s.pkl' % re.sub(r'\W+', '_', ad))


def toplu_cikar(istekler, ilerle=None):
    """Extracts ALL the pools in a SINGLE pass.

    istekler: [(cache_name, [key,...], class)]
    Each database is read at most ONCE. For 21 targets that is one pass per database
    instead of 21, which is a difference of hours.

    """
    os.makedirs(C.ONBELLEK, exist_ok=True)
    kalan = [(ad, [a.lower() for a in anah if a], sinif)
             for ad, anah, sinif in istekler
             if anah and not os.path.exists(_cache_yolu(ad))]
    if not kalan:
        return
    # group by database
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
        for baslik, seq in engine_gateway.read_fasta(db):
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
                        temiz = engine_gateway.clean(seq.upper())
                    kutular[ad].append((baslik[:150], temiz))
        for ad, _ in isler:
            pickle.dump(kutular[ad], open(_cache_yolu(ad), 'wb'))


def kapsam(havuz, F, R, lo, hi, max_mm=C.REFERANS_MAX_MM):
    """Havuzdaki kac kayit urun veriyor + boy dagilimi + vurus basliklari."""
    veren = 0
    boy = {}
    basliklar = []
    for ad, seq in havuz:
        p = engine_gateway.amplify(seq, F, R, max_mm=max_mm, lo=lo, hi=hi)
        if p:
            veren += 1
            b = min(p, key=lambda x: x[3] + x[4])[2]
            boy[b] = boy.get(b, 0) + 1
            if len(basliklar) < 40:
                basliklar.append((ad, b))
    return dict(veren=veren, toplam=len(havuz), boy=boy, basliklar=basliklar)


def uye_ve_rakip_anahtar(baglam, taxad):
    """Genus keys from the member taxid names; the competitor is the rest of the same class."""
    uye_adlar = [taxad.get(t, '') for t in baglam['uye_tax']]
    cinsler = sorted({a.split()[0] for a in uye_adlar if a and not a.startswith('Ca.')}
                     | {' '.join(a.split()[:2]) for a in uye_adlar if a.startswith('Ca.')})
    rakip_adlar = [taxad.get(k['taxid'], '') for k in baglam['rakip_kutu']]
    rakip_cins = sorted({a.split()[0] for a in rakip_adlar if a and a.split()[0] not in cinsler})
    return cinsler, rakip_cins
