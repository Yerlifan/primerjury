# -*- coding: utf-8 -*-
"""Aday uretimi: pencere -> primer -> cift -> ARMS varyanti -> parametre izgarasi.

TAMLIK NOTU
-----------
Omurga tek zincir uzerinde taranir ve bu YETERLIDIR: cift-zincirli sablonda
bir cift (+ zincirinde i konumundaki F, - zincirinde j>i konumundaki R) ile
tam olarak tanimlanir. Ters zincirden uretilen kume ayni amplikonlarin
aynasidir. Konsensuslerin bir kisminin ters tumleyen yonde saklanmis olmasi
(bkz. "11 B-F Yeniden Olcum" duzeltme 1) bu yuzden aramanin kapsamini
etkilemez.
"""
# ---------------------------------------------------------------------------
# generator.py — aday primer ve aday cift ureteci; 144 hucreli parametre
#             izgarasini ve ARMS varyantlarini da bu modul kurar.
#
# GIRDI  : omurga konsensus dizisi (hedefler.hedef_baglami()'nin sectigi en uzun
#          uye konsensusu); config.py'deki UZUNLUK, URUN_IDEAL,
#          URUN_MUTLAK_UST ve IZGARA_* sabitleri; primer olcumleri icin
#          geometri.olc / geometri.hucre_gecti; dizi islemleri icin motor.rc,
#          motor.encode, motor.find_sites.
# CIKTI  : dosyaya yazmaz. aday_primerler() {'F': [...], 'R': [...]} sozlugu,
#          tara_ve_topla() toplam cift sayisi + izgara sayaci + temsilci aday
#          listesi, izgara_tablosu_sayactan() 144 satirlik tablo,
#          arms_varyantlari() (varyant_dizi, etiket) ciftleri dondurur.
# CAGRAN : __main__.hedefi_isle icinden asama A, B ve B2'de; yani
#          full_chain.py asamalari 1, 2, 3, 7 ve 9 (7. asama). Ayrica
#          verification/recovery_round.py (tus K) bu modulu disaridan ice aktarir.
# ---------------------------------------------------------------------------
import itertools
from . import config as C
from . import geometry as G
from . import engine_gateway


# ---------------------------------------------------------------- pencereler
def pencereler(omurga, uz_ar=C.UZUNLUK):
    """Omurgadaki HER pozisyondan baslayan, uzunlugu uz_ar arasinda degisen
    her oligo. Donen: (baslangic, uzunluk, ileri_dizi)."""
    L = len(omurga)
    lo, hi = uz_ar
    for i in range(L):
        for k in range(lo, hi + 1):
            if i + k > L:
                break
            s = omurga[i:i + k]
            if 'N' in s:
                continue
            yield i, k, s


def aday_primerler(omurga, ilerle=None):
    """Her pencereden iki aday: ileri (oldugu gibi) ve geri (ters tumleyen).

    Geri primerin baglanma yeri omurgada [i, i+k); primerin kendisi rc(pencere).
    Donen sozluk: {'F': [(i, uz, dizi, olcum)], 'R': [...]}
    """
    F, R = [], []
    say = 0
    for i, k, s in pencereler(omurga):
        say += 1
        mF = G.olc(s)
        if G.sabit_gecti(mF):
            F.append((i, k, s, mF))
        r = motor.rc(s)
        mR = G.olc(r)
        if G.sabit_gecti(mR):
            R.append((i, k, r, mR))
        if ilerle and say % 2000 == 0:
            ilerle(say)
    return dict(F=F, R=R, taranan_pencere=say)


# ---------------------------------------------------------------- parametre izgarasi
# HIZ NOTU: 144 hucre icin 144 kez arama YAPILMAZ. Her primer, primer duzeyindeki
# 36 alt-kombinasyonun (3 GC x 3 Tm x 2 uc x 2 son5) hangilerini gectigini bir
# bit maskesinde tasir; bir cift bir hucreyi ancak IKI primeri de o alt-kombinasyonu
# gecerse ve urun boyu araligina duserse gecer. Boylece izgara tablosu ciftlerin
# uzerinden TEK gecisle cikar.

PRIMER_KOMBO = [(g, t, u, s)
                for g in C.IZGARA_GC for t in C.IZGARA_TM
                for u in C.IZGARA_UC_GC for s in C.IZGARA_SON5]        # 36
URUN_KOMBO = list(C.IZGARA_URUN)                                       # 4


def primer_maskesi(m):
    """m (geometri.olc ciktisi) -> 36 bitlik maske."""
    mask = 0
    for i, (g, t, u, s) in enumerate(PRIMER_KOMBO):
        if G.hucre_gecti(m, g, t, u, s):
            mask |= (1 << i)
    return mask


def urun_maskesi(bp):
    mask = 0
    for i, (lo, hi) in enumerate(URUN_KOMBO):
        if lo <= bp <= hi:
            mask |= (1 << i)
    return mask


def _hucre(pi, ui):
    g, t, u, s = PRIMER_KOMBO[pi]
    return dict(gc=g, tm=t, urun=URUN_KOMBO[ui], uc_gc=u, son5=s)


def izgara_hucreleri():
    for ui in range(len(URUN_KOMBO)):
        for pi in range(len(PRIMER_KOMBO)):
            yield _hucre(pi, ui)


def hucre_adi(h):
    return 'GC%d-%d|Tm%d-%d|urun%d-%d|3ucGC:%s|son5:%s' % (
        h['gc'][0], h['gc'][1], h['tm'][0], h['tm'][1], h['urun'][0], h['urun'][1],
        'sart' if h['uc_gc'] else 'serbest', '<=3' if h['son5'] else 'serbest')


def hucre_sikilik(h):
    """0 = en siki. Rapor 'en siki hangi ayarda cozum var' sorusunu bununla siralar."""
    return (C.IZGARA_GC.index(h['gc']) + C.IZGARA_TM.index(h['tm'])
            + C.IZGARA_URUN.index(h['urun']) + (0 if h['uc_gc'] else 1)
            + (0 if h['son5'] else 1))


_SIKILIK = [[hucre_sikilik(_hucre(pi, ui)) for ui in range(len(URUN_KOMBO))]
            for pi in range(len(PRIMER_KOMBO))]


def cift_maskesi(c):
    """Ciftin (primer_maske, urun_maske) ikilisi; bir kez hesaplanip saklanir."""
    if 'pm' not in c:
        c['pm'] = primer_maskesi(c['mF']) & primer_maskesi(c['mR'])
        c['um'] = urun_maskesi(c['urun'])
    return c['pm'], c['um']


def izgara_tablosu(cift_listesi):
    """Her izgara hucresi icin kac aday hayatta kaliyor (tek gecis)."""
    say = [[0] * len(URUN_KOMBO) for _ in range(len(PRIMER_KOMBO))]
    ornek = [[None] * len(URUN_KOMBO) for _ in range(len(PRIMER_KOMBO))]
    for c in cift_listesi:
        pm, um = cift_maskesi(c)
        if not pm or not um:
            continue
        for pi in range(len(PRIMER_KOMBO)):
            if not (pm >> pi) & 1:
                continue
            for ui in range(len(URUN_KOMBO)):
                if (um >> ui) & 1:
                    say[pi][ui] += 1
                    if ornek[pi][ui] is None:
                        ornek[pi][ui] = c
    tablo = []
    for pi in range(len(PRIMER_KOMBO)):
        for ui in range(len(URUN_KOMBO)):
            h = _hucre(pi, ui)
            o = ornek[pi][ui]
            tablo.append(dict(hucre=h, ad=hucre_adi(h), sikilik=_SIKILIK[pi][ui],
                              hayatta=say[pi][ui],
                              ornek=(o['F'] + ' / ' + o['R'] + ' (%d bp)' % o['urun'])
                              if o else ''))
    tablo.sort(key=lambda x: (x['sikilik'], -x['hayatta']))
    return tablo


def hucre_etiketle(c):
    """Bir cift hangi izgara hucrelerinde hayatta? En SIKI hucreyi dondurur."""
    pm, um = cift_maskesi(c)
    en = None
    if pm and um:
        for pi in range(len(PRIMER_KOMBO)):
            if not (pm >> pi) & 1:
                continue
            for ui in range(len(URUN_KOMBO)):
                if not (um >> ui) & 1:
                    continue
                s = _SIKILIK[pi][ui]
                if en is None or s < en[0]:
                    en = (s, hucre_adi(_hucre(pi, ui)))
    return en if en else (99, 'hicbir izgara hucresini gecmiyor')


# ---------------------------------------------------------------- cift kurma
def _urun_maske_tablosu(lo, hi):
    t = [0] * (hi + 1)
    for bp in range(lo, hi + 1):
        t[bp] = urun_maskesi(bp)
    return t


def cift_akisi(ad, urun_ar=(C.URUN_IDEAL[0], C.URUN_MUTLAK_UST)):
    """Kurala uyan HER ileri-geri kombinasyonunu AKIS olarak verir.

    Liste kurulmaz: milyonlarca cift bellege sigmaz, ama tek gecisle
    sayilabilir. Donen: (iF, kF, sF, mF, pmF, iR, kR, sR, mR, pmR, bp)
    """
    import bisect
    lo, hi = urun_ar
    Rs = sorted(ad['R'], key=lambda x: x[0] + x[1])
    uc = [x[0] + x[1] for x in Rs]
    pmR = [primer_maskesi(x[3]) for x in Rs]
    for iF, kF, sF, mF in ad['F']:
        pmF = primer_maskesi(mF)
        if not pmF:
            continue
        a = bisect.bisect_left(uc, iF + lo)
        b = bisect.bisect_right(uc, iF + hi)
        for j in range(a, b):
            iR, kR, sR, mR = Rs[j]
            if iR < iF + kF:
                continue
            pm = pmF & pmR[j]
            if not pm:
                continue
            yield (iF, kF, sF, mF, pmF, iR, kR, sR, mR, pmR[j], uc[j] - iF)


def cift_yap(t):
    """Akis demetini sozluge cevir (yalniz saklanacak adaylar icin)."""
    iF, kF, sF, mF, pmF, iR, kR, sR, mR, pmR, bp = t
    return dict(iF=iF, F=sF, mF=mF, iR=iR, R=sR, mR=mR, urun=bp,
                pm=pmF & pmR, um=urun_maskesi(bp))


def tara_ve_topla(ad, hucre_basina=6, urun_ar=(C.URUN_IDEAL[0], C.URUN_MUTLAK_UST),
                  ilerle=None):
    """TEK GECIS: butun ciftleri sayar, izgara tablosunu cikarir ve her
    izgara hucresi icin sinirli sayida temsilci aday saklar.

    Ust sinir yoktur - cift sayisi milyonlarca olsa da tamami sayilir.
    Bellekte yalniz temsilciler tutulur.
    """
    # Sayim ile saklama AYRILMISTIR. Izgara tablosu ciftlerin TAMAMI uzerinden
    # cikar (hicbir kesme yoktur, "kac aday hayatta kaldi" sayisi gercektir),
    # ama bellege yalniz her izgara hucresinden hucre_basina kadar temsilci
    # alinir. Boylece milyonlarca cift sayilabilirken bellek sabit kalir.
    # Temsilci secilirken urun boyu 105 bp'ye en yakin olan tutulur: qPCR icin
    # ideal aralik 60-150 bp'dir ve 105 bp o araligin ortasidir.
    from collections import Counter, defaultdict
    lo, hi = urun_ar
    umt = _urun_maske_tablosu(lo, hi)
    sayac = Counter()
    kova = defaultdict(list)
    toplam = 0
    for t in cift_akisi(ad, urun_ar):
        bp = t[10]
        um = umt[bp]
        if not um:
            continue
        pm = t[4] & t[9]
        toplam += 1
        anahtar = (pm, um)
        sayac[anahtar] += 1
        k = _en_siki_anahtar(anahtar)
        kutu = kova[k]
        if len(kutu) < hucre_basina:
            kutu.append(t)
        else:
            # daha ideal urun boyuna sahip olani tut (105 bp civari)
            en_kotu = max(range(len(kutu)), key=lambda i: abs(kutu[i][10] - 105))
            if abs(bp - 105) < abs(kutu[en_kotu][10] - 105):
                kutu[en_kotu] = t
        if ilerle and toplam % 250000 == 0:
            ilerle(toplam)
    return dict(toplam=toplam, sayac=sayac,
                temsilciler=[cift_yap(t) for kutu in kova.values() for t in kutu])


_ES_ONBELLEK = {}


def _en_siki_anahtar(anahtar):
    """(pm, um) -> en siki hucrenin (pi, ui) indeksi."""
    v = _ES_ONBELLEK.get(anahtar)
    if v is not None:
        return v
    pm, um = anahtar
    en = None
    for pi in range(len(PRIMER_KOMBO)):
        if not (pm >> pi) & 1:
            continue
        for ui in range(len(URUN_KOMBO)):
            if not (um >> ui) & 1:
                continue
            s = _SIKILIK[pi][ui]
            if en is None or s < en[0]:
                en = (s, pi, ui)
    v = (en[1], en[2]) if en else (-1, -1)
    _ES_ONBELLEK[anahtar] = v
    return v


def izgara_tablosu_sayactan(sayac):
    """tara_ve_topla'nin sayacindan 144 hucrelik tabloyu cikar (tam sayim)."""
    say = [[0] * len(URUN_KOMBO) for _ in range(len(PRIMER_KOMBO))]
    for (pm, um), n in sayac.items():
        for pi in range(len(PRIMER_KOMBO)):
            if not (pm >> pi) & 1:
                continue
            for ui in range(len(URUN_KOMBO)):
                if (um >> ui) & 1:
                    say[pi][ui] += n
    tablo = []
    for pi in range(len(PRIMER_KOMBO)):
        for ui in range(len(URUN_KOMBO)):
            h = _hucre(pi, ui)
            tablo.append(dict(hucre=h, ad=hucre_adi(h), sikilik=_SIKILIK[pi][ui],
                              hayatta=say[pi][ui], ornek=''))
    tablo.sort(key=lambda x: (x['sikilik'], -x['hayatta']))
    return tablo


# ---------------------------------------------------------------- ARMS
def ayirt_edici_mi(primer, uye_diziler, rakip_diziler, geri=False):
    """Primerin 3' SON BAZI ayirt edici konumda mi?

    Olcut: primer uye konsensusune 3' son baz TAM oturuyor; en iyi rakip
    baglanma yerinde ise 3' son baz UYMUYOR. Hizalama gerekmez - dogrudan
    olculur (ispcr.find_sites, gevsek uyumsuzluk tavani).
    """
    import numpy as np
    L = len(primer)

    def en_iyi(diziler, tam_uc):
        best = None
        for s in diziler:
            for d in (s, motor.rc(s)):
                enc = motor.encode(d)
                for mm_tavan in (0, 1, 2, 3, 4):
                    h = motor.find_sites(enc, primer, mm_tavan, need_tail=tam_uc,
                                         tail_pos=(-1,))
                    if h:
                        v = min(x[1] for x in h)
                        if best is None or v < best:
                            best = v
                        break
        return best

    uye = en_iyi(uye_diziler, True)
    if uye is None or uye > 1:
        return False, None, None
    # rakipte 3' son baz TAM tutan bir yer var mi?
    rak_tam = en_iyi(rakip_diziler, True)
    rak_gevsek = en_iyi(rakip_diziler, False)
    if rak_gevsek is None:
        return False, uye, None          # rakip zaten hic baglanmiyor, ARMS gereksiz
    if rak_tam is None or rak_tam > rak_gevsek:
        return True, uye, rak_gevsek     # rakipte 3' son baz UYMUYOR -> ayirt edici
    return False, uye, rak_gevsek


def arms_varyantlari(primer):
    """3' sondan 2. ve 3. baza kasitli uyumsuzluk. DORT BAZIN HEPSI denenir.

    NOT: kasitli uyumsuzluk DEJENERE BAZ DEGILDIR - tek tanimli bir bazdir,
    sentezlenen oligo sayisini artirmaz. Yine de sablonla tam eslesmedigi icin
    ayri bir toplanti maddesidir (raporda boyle bildirilir).
    Donen: (varyant_dizi, etiket)
    """
    # NEDEN -2 ve -3: 3' son baz ayirt ediciligi tasir ve DEGISTIRILMEZ; ARMS
    # mantigi, zaten ayirt edici olan son bazin yanina kasitli bir uyumsuzluk
    # koyarak polimerazin uzatmasini rakip sablonda daha da zorlastirmaktir.
    # Son baz degistirilseydi ayirt edicilik kaybolurdu.
    # 4 x 4 = 16 kombinasyondan degisiklik icermeyen tek durum atlandigi icin
    # 15 varyant uretilir; self_test.py bu sayiyi her kosuda dogrular.
    if len(primer) < 4:
        return []
    out = []
    p = list(primer)
    for b2 in 'ACGT':
        for b3 in 'ACGT':
            if b2 == p[-2] and b3 == p[-3]:
                continue                      # degisiklik yok
            q = p[:]
            q[-2], q[-3] = b2, b3
            et = []
            if b3 != p[-3]:
                et.append('-3:%s>%s' % (p[-3], b3))
            if b2 != p[-2]:
                et.append('-2:%s>%s' % (p[-2], b2))
            out.append((''.join(q), ' '.join(et)))
    return out


