# -*- coding: utf-8 -*-
"""Toplanti geometri kurallarinin PARAMETRELESTIRILMIS hali.

Kural metni ve primer3 tuz kosullari DUZELTME_betikleri/geo.py ile birebir
ayni tutulmustur; kendini_sina.py o dosyayi calistirip 42 panel primerinin
her degerini bu modulle karsilastirir (birebir tutmazsa arama baslamaz).

Fark: geo.py kurallari SABIT tutar, bu modul parametre izgarasi icin
gevsetilebilir esikler alir.
"""
# ---------------------------------------------------------------------------
# geometri.py — bir primerin ve bir ciftin fiziksel buyukluklerini (Tm, GC,
#               hairpin, homodimer, heterodimer, delta-G, Ta) olcer ve toplanti
#               kurallarina gore gecip gecmedigini soyler.
#
# GIRDI  : primer dizisi (ve cift icin iki dizi + urun boyu); esikler ve
#          primer3 tuz/derisim kosullari yapilandirma.py'den; hesaplama
#          primer3-py kutuphanesiyle yapilir.
# CIKTI  : dosyaya yazmaz. olc() olcum sozlugu, sabit_gecti() ve hucre_gecti()
#          True/False, ihlaller() ihlal metinleri listesi, cift_olc() cift
#          duzeyi sozluk dondurur.
# CAGRAN : uretec.aday_primerler ve uretec.primer_maskesi uzerinden asama A ve
#          B'de, __main__.hedefi_isle icinde asama D'de; kendini_sina.py
#          degerleri geo.py ile karsilastirir. Yani tuslar 1, 2, 3, 7, 8, 9;
#          disaridan verification/kurtarma_turu.py (tus K).
#
# sabit_gecti() ile hucre_gecti() bilerek AYRILMISTIR: hairpin ve homodimer
# SYBR Green kimyasinda eleyicidir ve hicbir izgara hucresinde gevsetilmez,
# bu yuzden sabit tarafta durur. GC, Tm ve 3' uc kurallari ise gevsetilebilir
# olduklari icin izgara tarafindadir - raporun "hangi kurali gevsetince cozum
# cikti" sorusuna cevap verebilmesi bu ayrima dayanir.
# ---------------------------------------------------------------------------
from . import yapilandirma as C

try:
    import primer3
except ImportError:
    primer3 = None

KW = dict(C.P3)


def _gerek():
    if primer3 is None:
        raise SystemExit(
            'HATA: primer3 modulu yok. WSL icinde su komutu calistirin:\n'
            '  pip3 install primer3-py --break-system-packages')


def tm(p):
    _gerek(); return round(primer3.calc_tm(p, **KW), 2)


def gc(p):
    return round(100.0 * sum(p.count(c) for c in 'GC') / len(p), 1)


def hp_tm(p):
    _gerek(); return round(primer3.calc_hairpin(p, **KW).tm, 1)


def hd_tm(p):
    _gerek(); return round(primer3.calc_homodimer(p, **KW).tm, 1)


def het_tm(a, b):
    _gerek(); return round(primer3.calc_heterodimer(a, b, **KW).tm, 1)


def hp_dg(p, t=C.TA_HEDEF):
    _gerek(); return round(primer3.calc_hairpin(p, temp_c=t, **KW).dg / 1000.0, 2)


def hd_dg(p, t=C.TA_HEDEF):
    _gerek(); return round(primer3.calc_homodimer(p, temp_c=t, **KW).dg / 1000.0, 2)


def het_dg(a, b, t=C.TA_HEDEF):
    _gerek(); return round(primer3.calc_heterodimer(a, b, temp_c=t, **KW).dg / 1000.0, 2)


def uc_dg(a, b, t=C.TA_HEDEF):
    """3' uc kaynakli dimer kararliligi (primer3 end_stability)."""
    _gerek()
    try:
        return round(primer3.calc_end_stability(a, b, temp_c=t, **KW).dg / 1000.0, 2)
    except TypeError:
        return round(primer3.calc_end_stability(a, b, **KW).dg / 1000.0, 2)


def tekrar_var(p, n=C.TEKRAR_UST):
    return any(b * n in p for b in 'ACGT')


def son5_gc(p):
    return sum(1 for c in p[-5:] if c in 'GC')


def dejenere(p):
    return any(c not in 'ACGT' for c in p)


# ---------------------------------------------------------------- olcumler
def olc(p):
    """Bir primerin butun geometri buyukluklerini bir kez hesapla (pahali kisim)."""
    return dict(
        dizi=p, uz=len(p), gc=gc(p), tm=tm(p),
        hp_tm=hp_tm(p), hd_tm=hd_tm(p),
        hp_dg=hp_dg(p), hd_dg=hd_dg(p),
        uc=p[-1], son5=son5_gc(p),
        tekrar=tekrar_var(p), dejenere=dejenere(p),
    )


# ------------------------------------------------- degismez (izgaradan bagimsiz) suzgec
def sabit_gecti(m):
    """Hangi izgara hucresinde olursak olalim gecerli olan kurallar.

    SYBR Green oldugu icin hairpin/homodimer ELEYICIDIR (uyari degil).
    """
    if not (C.UZUNLUK[0] <= m['uz'] <= C.UZUNLUK[1]):
        return False
    if m['dejenere'] or m['tekrar']:
        return False
    if m['hp_tm'] >= C.HAIRPIN_TM_UST:
        return False
    if m['hd_tm'] >= C.HOMODIMER_TM_UST:
        return False
    if m['hp_dg'] < C.DG_YAPI_ALT or m['hd_dg'] < C.DG_YAPI_ALT:
        return False
    return True


def hucre_gecti(m, gc_ar, tm_ar, uc_gc_sart, son5_sart):
    """Bir izgara hucresinin primer maddeleri."""
    if not (gc_ar[0] <= m['gc'] <= gc_ar[1]):
        return False
    if not (tm_ar[0] <= m['tm'] <= tm_ar[1]):
        return False
    if uc_gc_sart and m['uc'] not in 'GC':
        return False
    if son5_sart and m['son5'] > 3:
        return False
    return True


def ihlaller(p, gc_ar=(40, 60), tm_ar=(58, 62), uc_gc_sart=True, son5_sart=True):
    """geo.py'nin viol() fonksiyonuyla ayni cikti (varsayilan esiklerde birebir)."""
    v = []
    if not (C.UZUNLUK[0] <= len(p) <= C.UZUNLUK[1]):
        v.append('uz %d' % len(p))
    g = gc(p)
    if not (gc_ar[0] <= g <= gc_ar[1]):
        v.append('GC %%%.1f' % g)
    t = tm(p)
    if not (tm_ar[0] <= t <= tm_ar[1]):
        v.append('Tm %.2f' % t)
    if hp_tm(p) >= C.HAIRPIN_TM_UST:
        v.append('hairpin Tm %.1f' % hp_tm(p))
    if hd_tm(p) >= C.HOMODIMER_TM_UST:
        v.append('homodimer Tm %.1f' % hd_tm(p))
    if uc_gc_sart and p[-1] not in 'GC':
        v.append("3' uc %s (G/C degil)" % p[-1])
    n5 = son5_gc(p)
    if son5_sart and n5 > 3:
        v.append("3' son 5 bazda %d G/C" % n5)
    if dejenere(p):
        v.append('dejenere baz')
    return v


# ---------------------------------------------------------------- cift olcumleri
def cift_olc(mF, mR, urun_bp):
    """Cift duzeyi buyukluklerin tamami (60 C degerlendirmesi dahil)."""
    F, R = mF['dizi'], mR['dizi']
    dtm = round(abs(mF['tm'] - mR['tm']), 2)
    ht = het_tm(F, R)
    hg = het_dg(F, R)
    ug = uc_dg(F, R)
    tmin = min(mF['tm'], mR['tm'])
    ta_kural = round(tmin - C.TA_KURALI, 2)
    d = dict(
        dTm=dtm, het_tm=ht, het_dg=hg, uc_dg=ug,
        urun=urun_bp,
        Ta_kural=ta_kural,
        Ta60_marj=round(tmin - C.TA_HEDEF, 2),      # <0 ise 60 C'de primerler Tm'in ustunde
        Ta60_uygun=bool(ta_kural >= C.TA_HEDEF - 0.5),
        urun_sinifi=urun_sinifi(urun_bp),
    )
    d['cift_gecti'] = (dtm < C.DTM_UST and ht < C.HETERODIMER_TM_UST
                       and hg >= C.DG_YAPI_ALT and ug >= C.DG_UC_ALT)
    return d


def urun_sinifi(bp):
    if C.URUN_IDEAL[0] <= bp <= C.URUN_IDEAL[1]:
        return 'ideal (60-150)'
    if bp <= C.URUN_ONERILMEZ:
        return 'kabul edilebilir (150-250; protokolde 30 sn ann/ext)'
    return 'ONERILMEZ (>250)'
