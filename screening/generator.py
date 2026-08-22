# -*- coding: utf-8 -*-
"""Candidate production: window -> primer -> pair -> ARMS variant -> the parameter
grid.

A NOTE ON COMPLETENESS
----------------------
The backbone is scanned on a single strand and that IS ENOUGH: on a double stranded
template a pair is fully defined by (F at position i on the + strand, R at position
j>i on the - strand). The set produced from the reverse strand is the mirror of the
same amplicons. That some of the consensuses are stored reverse complemented (see
"11 B-F Yeniden Olcum" correction 1) therefore does not affect the coverage of the
search.

"""
# -------------------------------------------------------------------------
# generator.py, the candidate primer and candidate pair producer; this module also
#             builds the 144 cell parameter grid and the ARMS variants.
#
# INPUT  : the backbone consensus sequence (the longest member consensus, chosen by
#          hedefler.hedef_baglami()); the UZUNLUK, URUN_IDEAL, URUN_MUTLAK_UST and
#          IZGARA_* constants in config.py; geometri.olc and geometri.hucre_gecti
#          for the primer measurements; engine_gateway.rc, engine_gateway.encode and
#          engine_gateway.find_sites for the sequence work.
# OUTPUT : it writes no file. aday_primerler() returns a {'F': [...], 'R': [...]}
#          dictionary, tara_ve_topla() the total pair count plus the grid counter
#          plus a list of representative candidates, izgara_tablosu_sayactan() a 144
#          row table, and arms_varyantlari() (variant_sequence, label) pairs.
# CALLED BY: inside __main__.hedefi_isle at stages A, B and B2; that is,
#          full_chain.py stages 1, 2, 3, 7 and 9 (the 7th stage). Besides that,
#          verification/recovery_round.py (key K) imports this module from outside.
# -------------------------------------------------------------------------
import itertools
from . import config as C
from . import geometry as G
from . import engine_gateway


# ---------------------------------------------------------------- pencereler
def pencereler(omurga, uz_ar=C.UZUNLUK):
    'Every oligo starting at EVERY position of the backbone, with its length '
    'running     across the allowed range. Returns (start, length, '
    'forward_sequence).'
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
    """Two candidates from each window: forward (as it is) and reverse (the reverse
    complement).

    The reverse primer's binding site on the backbone is [i, i+k); the primer itself is
    rc(window).
    The dictionary returned: {'F': [(i, uz, dizi, olcum)], 'R': [...]}

    """
    F, R = [], []
    say = 0
    for i, k, s in pencereler(omurga):
        say += 1
        mF = G.olc(s)
        if G.sabit_gecti(mF):
            F.append((i, k, s, mF))
        r = engine_gateway.rc(s)
        mR = G.olc(r)
        if G.sabit_gecti(mR):
            R.append((i, k, r, mR))
        if ilerle and say % 2000 == 0:
            ilerle(say)
    return dict(F=F, R=R, taranan_pencere=say)


# ---------------------------------------------------------------- the parameter grid
# A NOTE ON SPEED: the search IS NOT run 144 times for 144 cells. Every primer
# carries, in a bit mask, which of the 36 primer level subcombinations
# (3 GC x 3 Tm x 2 end x 2 last5) it passed; a pair passes a cell only if BOTH its
# primers passed that subcombination and the product length falls in the range. So
# the grid table comes out in a SINGLE pass over the pairs.

PRIMER_KOMBO = [(g, t, u, s)
                for g in C.IZGARA_GC for t in C.IZGARA_TM
                for u in C.IZGARA_UC_GC for s in C.IZGARA_SON5]        # 36
URUN_KOMBO = list(C.IZGARA_URUN)                                       # 4


def primer_maskesi(m):
    'm, the geometry measurement output, becomes a 36 bit mask.'
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
    """0 = the strictest. The report orders the question 'under which strictest setting is there a solution' by this."""
    return (C.IZGARA_GC.index(h['gc']) + C.IZGARA_TM.index(h['tm'])
            + C.IZGARA_URUN.index(h['urun']) + (0 if h['uc_gc'] else 1)
            + (0 if h['son5'] else 1))


_SIKILIK = [[hucre_sikilik(_hucre(pi, ui)) for ui in range(len(URUN_KOMBO))]
            for pi in range(len(PRIMER_KOMBO))]


def cift_maskesi(c):
    "The pair's (primer_mask, product_mask); computed once and cached."
    if 'pm' not in c:
        c['pm'] = primer_maskesi(c['mF']) & primer_maskesi(c['mR'])
        c['um'] = urun_maskesi(c['urun'])
    return c['pm'], c['um']


def izgara_tablosu(cift_listesi):
    """How many candidates survive for each grid cell (a single pass)."""
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
    """In which grid cells does a pair survive? It returns the STRICTEST cell."""
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
    return en if en else (99, 'it passes no cell of the grid')


# ---------------------------------------------------------------- cift kurma
def _urun_maske_tablosu(lo, hi):
    t = [0] * (hi + 1)
    for bp in range(lo, hi + 1):
        t[bp] = urun_maskesi(bp)
    return t


def cift_akisi(ad, urun_ar=(C.URUN_IDEAL[0], C.URUN_MUTLAK_UST)):
    """Yields EVERY forward and reverse combination obeying the rule AS A STREAM.

    No list is built: millions of pairs do not fit in memory, but they can be counted
    in a single pass. Returned: (iF, kF, sF, mF, pmF, iR, kR, sR, mR, pmR, bp)

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
    """Turn a stream tuple into a dictionary (only for the candidates that will be kept)."""
    iF, kF, sF, mF, pmF, iR, kR, sR, mR, pmR, bp = t
    return dict(iF=iF, F=sF, mF=mF, iR=iR, R=sR, mR=mR, urun=bp,
                pm=pmF & pmR, um=urun_maskesi(bp))


def tara_ve_topla(ad, hucre_basina=6, urun_ar=(C.URUN_IDEAL[0], C.URUN_MUTLAK_UST),
                  ilerle=None):
    """A SINGLE PASS: it counts every pair, produces the grid table and keeps a limited
    number of representative candidates for each grid cell.

    There is no upper bound; even if the pair count runs into the millions all of them
    are counted. Only the representatives are held in memory.

    """
    # Counting and storing ARE SEPARATE. The grid table comes out over ALL the pairs
    # (there is no cut anywhere, so the "how many candidates survived" number is real),
    # but only up to hucre_basina representatives per grid cell are kept in memory. That
    # way millions of pairs can be counted while the memory stays constant.
    # When a representative is chosen, the one with the product length closest to 105 bp
    # is kept: the ideal range for qPCR is 60-150 bp and 105 bp is the middle of it.
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
            # keep the one with the more ideal product length (around 105 bp)
            en_kotu = max(range(len(kutu)), key=lambda i: abs(kutu[i][10] - 105))
            if abs(bp - 105) < abs(kutu[en_kotu][10] - 105):
                kutu[en_kotu] = t
        if ilerle and toplam % 250000 == 0:
            ilerle(toplam)
    return dict(toplam=toplam, sayac=sayac,
                temsilciler=[cift_yap(t) for kutu in kova.values() for t in kutu])


_ES_ONBELLEK = {}


def _en_siki_anahtar(anahtar):
    '(pm, um) -> the index (pi, ui) of the strictest cell.'
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
    'Build the 144 cell table out of the scan counter, as a full count.'
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
    """Is the primer's LAST BASE AT THE 3' END in a discriminating position?

    The criterion: the primer sits on the member consensus with the last base at the 3'
    end EXACTLY, while at the best competitor binding site that last base DOES NOT
    match. No alignment is needed, it is measured directly (ispcr.find_sites, with a
    relaxed mismatch ceiling).

    """
    import numpy as np
    L = len(primer)

    def en_iyi(diziler, tam_uc):
        best = None
        for s in diziler:
            for d in (s, engine_gateway.rc(s)):
                enc = engine_gateway.encode(d)
                for mm_tavan in (0, 1, 2, 3, 4):
                    h = engine_gateway.find_sites(enc, primer, mm_tavan, need_tail=tam_uc,
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
    # is there a place in the competitor where the last base at the 3' end matches EXACTLY?
    rak_tam = en_iyi(rakip_diziler, True)
    rak_gevsek = en_iyi(rakip_diziler, False)
    if rak_gevsek is None:
        return False, uye, None          # rakip zaten hic baglanmiyor, ARMS gereksiz
    if rak_tam is None or rak_tam > rak_gevsek:
        return True, uye, rak_gevsek     # rakipte 3' son baz UYMUYOR -> ayirt edici
    return False, uye, rak_gevsek


def arms_varyantlari(primer):
    """A deliberate mismatch at the 2nd and 3rd base from the 3' end. ALL FOUR BASES are
    tried.

    NOTE: a deliberate mismatch IS NOT A DEGENERATE BASE, it is one defined base and it
    does not raise the number of oligos synthesised. All the same, because it does not
    match the template exactly it is a separate meeting item (and it is reported as such
    in the report).
    Returned: (variant_sequence, label)

    """
    # WHY -2 AND -3: the last base at the 3' end carries the discrimination and IS NOT
    # CHANGED; the ARMS idea is to put a deliberate mismatch beside that already
    # discriminating last base, so the polymerase has an even harder time extending on
    # the competitor template. Had the last base been changed, the discrimination would
    # have been lost.
    # Of the 4 x 4 = 16 combinations the one case with no change is skipped, so 15
    # variants are produced; self_test.py confirms that number on every run.
    if len(primer) < 4:
        return []
    out = []
    p = list(primer)
    for b2 in 'ACGT':
        for b3 in 'ACGT':
            if b2 == p[-2] and b3 == p[-3]:
                continue                      # no change
            q = p[:]
            q[-2], q[-3] = b2, b3
            et = []
            if b3 != p[-3]:
                et.append('-3:%s>%s' % (p[-3], b3))
            if b2 != p[-2]:
                et.append('-2:%s>%s' % (p[-2], b2))
            out.append((''.join(q), ' '.join(et)))
    return out


