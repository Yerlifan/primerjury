# -*- coding: utf-8 -*-
"""The PARAMETERISED form of the meeting's geometry rules.

The rule text and the primer3 salt conditions are kept exactly the same as in
engine/geometry_core.py; self_test.py runs that file and compares every value of the
42 panel primers against this module (if they do not match exactly, the search does
not start).

The difference: geometry_core.py holds the rules FIXED, while this module takes
relaxable thresholds for the parameter grid.

"""
# -------------------------------------------------------------------------
# geometry.py measures a primer's and a pair's physical quantities (Tm, GC,
#               hairpin, homodimer, heterodimer, delta-G, Ta) and says whether they
#               pass the rules from the meeting.
#
# INPUT  : the primer sequence (and, for a pair, two sequences plus the product
#          length); the thresholds and the primer3 salt and concentration conditions
#          come from config.py; the calculation is done with the primer3-py library.
# OUTPUT : it writes no file. olc() returns a measurement dictionary,
#          sabit_gecti() and hucre_gecti() True or False, ihlaller() a list of
#          violation texts, and cift_olc() a pair level dictionary.
# CALLED BY: at stages A and B through uretec.aday_primerler and
#          uretec.primer_maskesi, and at stage D inside __main__.hedefi_isle;
#          self_test.py compares its values against geometry_core.py. That is keys
#          1, 2, 3, 7, 8, 9; from outside, verification/recovery_round.py (key K).
#
# sabit_gecti() and hucre_gecti() ARE DELIBERATELY SEPARATE: hairpin and homodimer
# are eliminating in SYBR Green chemistry and are relaxed in no grid cell, so they
# sit on the fixed side. The GC, Tm and 3' end rules can be relaxed and are
# therefore on the grid side; the report's ability to answer "which rule had to be
# relaxed for a solution to appear" rests on that separation.
# -------------------------------------------------------------------------
from . import config as C

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
    """Compute all of a primer's geometry quantities once (the expensive part)."""
    return dict(
        dizi=p, uz=len(p), gc=gc(p), tm=tm(p),
        hp_tm=hp_tm(p), hd_tm=hd_tm(p),
        hp_dg=hp_dg(p), hd_dg=hd_dg(p),
        uc=p[-1], son5=son5_gc(p),
        tekrar=tekrar_var(p), dejenere=dejenere(p),
    )


# ------------------------------------------------- degismez (izgaradan bagimsiz) suzgec
def sabit_gecti(m):
    """The rules that hold whichever grid cell we are in."""
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
    """geometry_core.py'nin viol() fonksiyonuyla ayni cikti (varsayilan esiklerde birebir)."""
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
