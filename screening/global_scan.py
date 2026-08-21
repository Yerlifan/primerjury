# -*- coding: utf-8 -*-
"""Global specificity scan, the most expensive step, run last.

WHAT IT DOES
    Scans a full reference database (SILVA SSURef NR99, ~500k records; UNITE or
    LSU when needed) for places where both primers bind in opposing orientation
    within the product window. The criterion matches the panel's global rule:
    at most 5 mismatches in total, F and R reported separately.

    This is a RAW SCANNER. It counts every product it finds; deciding which of
    those are genuinely off-target is the caller's job (see
    verification/specificity_round.py). Keeping the two apart is deliberate, a
    scanner that also judged would make its own bugs invisible.

EFFICIENCY
    The database is read ONCE and all candidates are measured against each
    chunk together. A per-candidate pass would re-read ~500k records once per
    candidate. Chunk size is bounded because peak memory is roughly six times
    the chunk.

RESUMABILITY
    Every chunk writes its state to disk; restarting continues from the last
    completed chunk.

OPTIONAL TAXONOMIC CLASSIFICATION
    Pass `siniflandirici` to have every hit classified into the four D-12
    classes (inside clade / organelle / same domain outside clade / different
    domain) and COUNTED. Identities are not stored: one universal primer alone
    produced 483,098 hits, and keeping headers for those would cost ~100 MB.
    Counters are constant memory per candidate and the count is complete; the
    `vurus` list stays capped at 300 as an evidence sample.

    'siniflandirildi' False means the classification was NOT RUN. It does not
    mean "no cross-reaction found". Callers must not conflate the two.

--- ozgun aciklama ---
KURESEL OZGULLUK - aramanin EN PAHALI adimi, en sona birakilir.

REFERANS_DB altindaki tam veritabanina (SILVA SSURef NR99 ~500 bin kayit,
gerekirse UNITE/LSU) karsi tarama. Olcut panelin kuresel olcutuyle ayni:
toplam <=5 uyumsuzluk, F ve R ayri bildirilir.

Verimlilik: veritabani BIR KEZ okunur; o parcada BUTUN son adaylar birlikte
olculur. Aday basina ayri gecis yapilmaz.

Kesinti dayanikliligi: her parca bitince ara sonuc diske yazilir; program
yeniden baslatilinca kalinan parcadan devam eder.
"""
# ---------------------------------------------------------------------------
# global_scan.py, hayatta kalan son adaylari REFERANS_DB'deki tam
#                     veritabanina karsi tarar; aramanin en pahali adimidir.
#
# GIRDI  : aday listesi ({'ad','F','R','lo','hi'}); veritabani olarak
#          yapilandirma.SILVA_SSU (ya da cagiranin verdigi baska fasta);
#          varsa kontrol/kuresel_<hedef>.pkl ara durumu. Tarama motor.encode
#          ve motor.find_sites (numpy vektor arama) ile yapilir.
# CIKTI  : durum_yolu verilirse her parca sonunda pickle ara durumu yazar.
#          tara() {aday_ad: {'urun': n, 'boy': {...}, 'vurus': [...]}} sozlugu
#          dondurur.
# CAGRAN : __main__.hedefi_isle icinde asama E - yani verification/full_chain.py
#          tuslari 1, 2, 3, 7 ve 9'un 7. asamasi (--hafif verilirse atlanir).
#          Disaridan verification/specificity_round.py (tus D) AYNEN bu modulu
#          kullanir; dogrulama turunun yerel veritabani katmani budur.
#
# VERIMLILIK GEREKCESI: veritabani BIR KEZ okunur ve her parcada BUTUN adaylar
# birlikte olculur. Aday basina ayri gecis yapilsaydi ~500 bin kayit aday
# sayisi kadar kez taranirdi. Parca boyu KURESEL_PARCA ile sinirlidir cunku
# bellek tavani yaklasik o boyutun alti katidir.
# ---------------------------------------------------------------------------
import os, json, pickle, time
import numpy as np
from . import config as C
from . import engine_gateway

AYIRAC = 'N' * 60


def _parcalar(db, parca_baz=C.KURESEL_PARCA):
    """FASTA'yi ~parca_baz buyuklugunde bloklara bol; her blok:
    (baslik listesi, uzunluk listesi, birlestirilmis dizi)."""
    ad_l, uz_l, buf, tot = [], [], [], 0
    for ad, seq in motor.read_fasta(db):
        s = motor.clean(seq.upper())
        # A2 (2026-08-21): baslik ARTIK KESILMIYOR (eskiden ad[:150] idi).
        # Olculdu: SILVA SSU basliklarinin %16,6'si 150 karakteri asiyor ve
        # kesilen kuyruk tam da CINS ve TUR jetonlaridir - yani hedef kladla
        # eslesme ihtimali en yuksek olanlar. Kesik baslikla siniflandirmak,
        # klad ICI bir kaydi klad DISI saydirir; duzeltmeye calistigimiz
        # hatanin ta kendisini uretir.
        # Maliyet olculdu ve onemsiz: parca 40 MB dizi ~= 28.500 kayit,
        # ortalama baslik 134 karakter -> parca basina ~3,8 MB. ad_l parca
        # yereli, pickle'a YAZILMAZ.
        ad_l.append(ad); uz_l.append(len(s)); buf.append(s)
        tot += len(s) + len(AYIRAC)
        if tot >= parca_baz:
            yield ad_l, uz_l, AYIRAC.join(buf)
            ad_l, uz_l, buf, tot = [], [], [], 0
    if buf:
        yield ad_l, uz_l, AYIRAC.join(buf)


def _kayit_indeksi(uz_l):
    """Birlestirilmis dizide her kaydin baslangici (ayirac dahil)."""
    off = np.zeros(len(uz_l), dtype=np.int64)
    c = 0
    for i, u in enumerate(uz_l):
        off[i] = c
        c += u + len(AYIRAC)
    return off


# Kontrol noktasi bicim surumu. A2 ile 'sinif' sayaclari eklendi; bu alani
# TASIMAYAN eski pickle'lar GECERSIZDIR. Surum kontrolu olmasaydi eski kontrol
# noktasi sessizce geri okunur ve taksonomik sayaclar sifir kalirdi - "olculdu
# ve capraz cikmadi" ile "hic olculmedi" ayirt edilemezdi.
DURUM_SURUMU = 2

SINIFLAR = ('a', 'ao', 'b', 'c', 'bilinmiyor')


def _bos_sonuc():
    return dict(urun=0, boy={}, vurus=[],
                sinif={k: 0 for k in SINIFLAR}, siniflandirildi=False)


def tara(adaylar, db=None, durum_yolu=None, ilerle=None, max_mm=C.KURESEL_MAX_MM,
         siniflandirici=None):
    """adaylar: [{'ad':..,'F':..,'R':..,'lo':..,'hi':..}]  (az sayida olmali)

    siniflandirici : None ya da  f(aday_ad, baslik, db_dosya_adi) -> sinif dizgesi
        Verilirse her vurus D-12'nin siniflarindan birine sokulur ve SAYILIR.
        A2 (2026-08-21): kimlikleri SAKLAMAK yerine SAYMAK bilincli bir karardir.
        Olculdu - Bakteri_universal tek basina 483.098 vurus veriyor; her biri
        icin baslik saklamak ~100 MB eder ve 'vurus' listesi zaten 300'de
        kesiliyor. Sayac ise aday basina SABIT bellek tutar ve sayi EKSIKSIZ
        olur. 'vurus' listesi kanit ornegi olarak 300'de kalir.

    Donen: {aday_ad: {'urun':n, 'boy':{}, 'vurus':[(baslik,boy,mmF,mmR)],
                      'sinif':{'a':n,'ao':n,'b':n,'c':n,'bilinmiyor':n},
                      'siniflandirildi':bool}}
    'siniflandirildi' False ise sinif sayaclari SIFIRDIR ama bu "capraz yok"
    DEMEK DEGILDIR - olcum yapilmamis demektir. Cagiran taraf ikisini
    karistirmamalidir.
    """
    db = db or C.SILVA_SSU
    if not os.path.exists(db):
        return {a['ad']: dict(hata='veritabani yok: %s' % db) for a in adaylar}
    db_ad = os.path.basename(db)

    durum = dict(surum=DURUM_SURUMU, parca=0, toplam_kayit=0,
                 res={a['ad']: _bos_sonuc() for a in adaylar})
    if durum_yolu and os.path.exists(durum_yolu):
        try:
            eski = pickle.load(open(durum_yolu, 'rb'))
            if eski.get('surum') == DURUM_SURUMU:
                durum = eski
            # surum tutmuyorsa BASTAN taranir; sessizce eski sonuc DONMEZ.
        except Exception:
            pass

    baslangic = durum['parca']
    pi = 0
    t0 = time.time()
    for ad_l, uz_l, big in _parcalar(db):
        pi += 1
        if pi <= baslangic:
            continue
        enc = motor.encode(big)
        off = _kayit_indeksi(uz_l)
        durum['toplam_kayit'] += len(ad_l)
        for a in adaylar:
            F, R = a['F'], a['R']
            revrc = motor.rc(R)
            fs = motor.find_sites(enc, F, max_mm, need_tail=False)
            if not fs:
                continue
            rs = motor.find_sites(enc, revrc, max_mm, need_tail=False)
            if not rs:
                continue
            fpos = np.array([x[0] for x in fs]); fmm = np.array([x[1] for x in fs])
            rpos = np.array([x[0] for x in rs]); rmm = np.array([x[1] for x in rs])
            frec = np.searchsorted(off, fpos, 'right') - 1
            rrec = np.searchsorted(off, rpos, 'right') - 1
            # kayit ici konum
            fin = fpos - off[frec]; rin = rpos - off[rrec]
            gecerli_f = fin + len(F) <= np.array(uz_l)[frec]
            gecerli_r = rin + len(revrc) <= np.array(uz_l)[rrec]
            fpos, fmm, frec, fin = fpos[gecerli_f], fmm[gecerli_f], frec[gecerli_f], fin[gecerli_f]
            rpos, rmm, rrec, rin = rpos[gecerli_r], rmm[gecerli_r], rrec[gecerli_r], rin[gecerli_r]
            ort = set(frec.tolist()) & set(rrec.tolist())
            if not ort:
                continue
            r = durum['res'][a['ad']]
            for kid in sorted(ort):
                fi = fin[frec == kid]; fm = fmm[frec == kid]
                ri = rin[rrec == kid]; rm = rmm[rrec == kid]
                en = None
                for x, xm in zip(fi, fm):
                    for yy, ym in zip(ri, rm):
                        bp = int(yy + len(revrc) - x)
                        if a['lo'] <= bp <= a['hi'] and yy >= x + len(F) and xm + ym <= max_mm:
                            if en is None or xm + ym < en[1] + en[2]:
                                en = (bp, int(xm), int(ym))
                if en:
                    r['urun'] += 1
                    r['boy'][en[0]] = r['boy'].get(en[0], 0) + 1
                    # A2: SAYAC tavansiz, KIMLIK listesi tavanli.
                    # Taksonomik hukum sayaclardan uretilir; 'vurus' yalnizca
                    # insana gosterilecek kanit ornegidir.
                    if siniflandirici is not None:
                        try:
                            s = siniflandirici(a['ad'], ad_l[kid], db_ad)
                        except Exception:
                            s = 'bilinmiyor'
                        if s not in r['sinif']:
                            s = 'bilinmiyor'
                        r['sinif'][s] += 1
                        r['siniflandirildi'] = True
                    if len(r['vurus']) < 300:
                        r['vurus'].append((ad_l[kid], en[0], en[1], en[2]))
        durum['parca'] = pi
        if durum_yolu:
            pickle.dump(durum, open(durum_yolu, 'wb'))
        if ilerle:
            ilerle(pi, durum['toplam_kayit'], time.time() - t0)
        del enc, big
    return durum['res']
