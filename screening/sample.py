# -*- coding: utf-8 -*-
"""In-silico PCR ON THE SAMPLE'S RAW READS, plus the Wilson discrimination ratio.

The criterion is exactly the panel's sample criterion: mismatches <=1 and the last
2 bases at the 3' end EXACT. The engine is the Havuz class of engine/scanner.py
plus the urunler() function of pair.py; this file IMPORTS them, it does not
rewrite them.

The pool is built ONCE per bin and cached; every later candidate pair then costs
only two index lookups. That is why thousands of candidates can be scanned.

"""
# -------------------------------------------------------------------------
# sample.py builds the raw read pools per bin and measures, for one primer pair, the
#           product ratios in the member and competitor bins with a Wilson interval.
#
# INPUT  : the bin dictionaries hedefler.kutular() returns (the bin name plus the
#          fastq path); NUMUNE_OKUMA_MIN/MAX, NUMUNE_OKUMA_SAYISI, NUMUNE_TOHUM,
#          NUMUNE_MAX_MM, KAPSAM_ESIGI and ENKOTU_ASGARI_OKUMA from config.py; the
#          measurement engines through engine_gateway.py (okuma_motoru,
#          tarayici.Havuz, cift.urunler, ispcr.find_sites).
# OUTPUT : it writes no file. Numune.olc() returns a single dictionary: the member
#          and competitor bin counts, the member percentages, the coverage share, the
#          competitor pool ratio, the discrimination ratios (kat_havuz, kat_enkotu
#          and their counterparts computed over the covered bins only) and the
#          product length distribution.
# CALLED BY: __main__.hedefi_isle, panel_olcum.calistir and
#          uyelik_denetimi.calistir - that is, full_chain.py stages 1, 2, 3, 4, 5, 7
#          and 9. From outside, protocol/single_protocol_measure.py (key P) and
#          verification/recovery_round.py (key K) import it too.
# -------------------------------------------------------------------------
import os, glob, gzip, math, random, pickle
from . import config as C
from . import engine_gateway


def okumalar(yol, n=C.NUMUNE_OKUMA_SAYISI, tohum=C.NUMUNE_TOHUM):
    """The same as okuma_pcr.py: a 200-6000 bp filter, sampling with a fixed seed."""
    ac = gzip.open if yol.endswith('.gz') else open
    hepsi = []
    with ac(yol, 'rt', errors='ignore') as fh:
        for i, satir in enumerate(fh):
            if i % 4 == 1:
                s = satir.strip().upper()
                if C.NUMUNE_OKUMA_MIN <= len(s) <= C.NUMUNE_OKUMA_MAX:
                    hepsi.append(s)
    if n and len(hepsi) > n:
        hepsi = random.Random(tohum).sample(hepsi, n)
    return hepsi


class KutuOtorite:
    """Measurement with the CORRECTED engine (read_engine.py) - THE AUTHORITATIVE route.

        Pigeonhole seeding is lossless and has been verified exactly against
        brute_force.py. It is slower than the fast numpy route; this route is used when
        few pairs are being measured (the panel re-measurement).

    """

    def __init__(self, kutu, yol, n=C.NUMUNE_OKUMA_SAYISI):
        self.kutu = kutu
        rd, n0 = engine_gateway.okuma_motoru.kutu_yukle(
            yol, nmax=(n or 0), seed=C.NUMUNE_TOHUM,
            minl=C.NUMUNE_OKUMA_MIN, maxl=C.NUMUNE_OKUMA_MAX)
        self.reads = rd
        self.n_okuma = len(rd)
        self.suzgecten_gecen = n0

    def urun_veren(self, F, R, lo, hi, mm=C.NUMUNE_MAX_MM):
        pos, n, boy = engine_gateway.okuma_motoru.kutu_pcr(
            self.reads, F, R, lo=lo, hi=hi, max_mm=mm, son2=True)
        return pos, n, dict(boy)


class KutuEski:
    """Measurement with THE PANEL'S OLD (FAULTY) engine - FOR COMPARISON ONLY.

        The 13 base EXACT matching seed of reads.py/Sonda is reproduced here exactly.
        This route is used for no decision at all; its purpose is to make VISIBLE how
        much difference the fix makes on which row.

    """

    def __init__(self, kutu, yol, n=C.NUMUNE_OKUMA_SAYISI):
        self.kutu = kutu
        rd, n0 = engine_gateway.okuma_motoru.kutu_yukle(
            yol, nmax=(n or 0), seed=C.NUMUNE_TOHUM,
            minl=C.NUMUNE_OKUMA_MIN, maxl=C.NUMUNE_OKUMA_MAX)
        self.reads = rd
        self.n_okuma = len(rd)

    @staticmethod
    def _sonda(primer, uc5, max_mm, seed=13):
        import itertools
        om = engine_gateway.okuma_motoru
        sd = primer[:seed] if uc5 else primer[-seed:]
        off = 0 if uc5 else len(primer) - seed
        tohumlar = [''.join(x) for x in itertools.product(
            *[om.IUPAC.get(c, 'ACGT') for c in sd])]
        L = len(primer)

        def bul(seq):
            out = []
            for t in tohumlar:
                i = seq.find(t)
                while i != -1:
                    st = i - off
                    if 0 <= st and st + L <= len(seq):
                        mm = 0; iyi = True
                        for a, b in zip(primer, seq[st:st + L]):
                            if b not in om.IUPAC.get(a, 'ACGT'):
                                mm += 1
                                if mm > max_mm:
                                    iyi = False; break
                        if iyi:
                            out.append((st, mm))
                    i = seq.find(t, i + 1)
            return out
        return bul

    def urun_veren(self, F, R, lo, hi, mm=C.NUMUNE_MAX_MM):
        om = engine_gateway.okuma_motoru
        fs = self._sonda(F, False, mm)
        rs = self._sonda(om.rc(R), True, mm)
        pos = 0
        boy = {}
        for s in self.reads:
            vur = None
            for seq in (s, om.rc(s)):
                a = fs(seq)
                if not a:
                    continue
                b = rs(seq)
                if not b:
                    continue
                for i, _ in a:
                    for j, _ in b:
                        n = j + len(R) - i
                        if lo <= n <= hi and j >= i + len(F):
                            vur = n; break
                    if vur:
                        break
                if vur:
                    break
            if vur:
                pos += 1
                boy[vur] = boy.get(vur, 0) + 1
        return pos, self.n_okuma, boy


class KutuHavuzu:
    """Bir kutunun okumalari + ters tumleyenleri, tarayici.Havuz indeksiyle."""

    def __init__(self, kutu, yol, n=C.NUMUNE_OKUMA_SAYISI):
        self.kutu = kutu
        rd = okumalar(yol, n)
        self.n_okuma = len(rd)
        # every read is tried in both directions -> put both into the pool, then match
        diziler = []
        self.okuma_id = []
        for i, s in enumerate(rd):
            s = engine_gateway.clean(s)
            diziler.append(s); self.okuma_id.append(i)
            diziler.append(engine_gateway.rc(s)); self.okuma_id.append(i)
        self.hv = engine_gateway.tarayici.Havuz(diziler) if diziler else None

    def urun_veren_kaba(self, F, R, lo, hi, mm):
        """A SEEDLESS (brute force) scan - CORRECT for every mismatch cap.

                The seed in tarayici.Havuz is complete only for a SINGLE mismatch (the seed
                variants are produced with one substitution). Under looser criteria such as
                <=3 a seeded search MISSES sites; this route uses ispcr.find_sites directly on
                the pool's concatenated sequence, so nothing is missed.

        """
        import numpy as np
        if self.hv is None or self.n_okuma == 0:
            return 0, 0, {}
        enc, sid = self.hv.enc, self.hv.sid
        revrc = engine_gateway.rc(R)
        fs = engine_gateway.find_sites(enc, F, mm, need_tail=True, tail_pos=(-1, -2))
        if not fs:
            return 0, self.n_okuma, {}
        rs = engine_gateway.find_sites(enc, revrc, mm, need_tail=True, tail_pos=(0, 1))
        if not rs:
            return 0, self.n_okuma, {}
        fpos = np.array([x[0] for x in fs]); rpos = np.array([x[0] for x in rs])
        fid = sid[fpos]; rid = sid[rpos]
        ok = fid >= 0; fpos, fid = fpos[ok], fid[ok]
        ok = rid >= 0; rpos, rid = rpos[ok], rid[ok]
        from collections import defaultdict
        rmap = defaultdict(list)
        for pz, i in zip(rpos.tolist(), rid.tolist()):
            rmap[i].append(pz)
        veren = set(); boy = {}
        for pz, i in zip(fpos.tolist(), fid.tolist()):
            for q in rmap.get(i, ()):
                n = q + len(revrc) - pz
                if lo <= n <= hi and q >= pz + len(F):
                    veren.add(self.okuma_id[i])
                    boy[n] = boy.get(n, 0) + 1
                    break
        return len(veren), self.n_okuma, boy

    def urun_veren(self, F, R, lo, hi, mm=C.NUMUNE_MAX_MM):
        """The number of READS giving a product (both directions count as one read) and the length distribution."""
        if self.hv is None or self.n_okuma == 0:
            return 0, 0, {}
        if mm > 1:
            return self.urun_veren_kaba(F, R, lo, hi, mm)
        m, boy = engine_gateway.cift.urunler(self.hv, F, R, lo=lo, hi=hi, mm=mm)
        veren = set()
        for idx in m.nonzero()[0]:
            veren.add(self.okuma_id[idx])
        return len(veren), self.n_okuma, boy


class Numune:
    """otorite=True -> read_engine.py (lossless, slow; the panel re-measurement)
           otorite=False -> the numpy pool (fast; scanning thousands of candidates)

    """

    def __init__(self, kutular, n=C.NUMUNE_OKUMA_SAYISI, ilerle=None, otorite=False):
        self.havuz = {}
        self.otorite = otorite
        sinif = KutuOtorite if otorite else KutuHavuzu
        for i, k in enumerate(kutular):
            if ilerle:
                ilerle(i + 1, len(kutular), k['kutu'])
            self.havuz[k['kutu']] = sinif(k['kutu'], k['yol'], n)

    def olc(self, F, R, uye_kutu, rakip_kutu, lo, hi, mm=C.NUMUNE_MAX_MM):
        # HOW THE DISCRIMINATION RATIO IS BUILT, AND WHY THIS WAY
        # The numerator: the Wilson LOWER bound of the WORST of the member bins. The worst
        # rather than the average, because the panel's promise is "it amplifies in every
        # member bin"; if a single member bin comes out empty, the pair does not work for
        # that bin.
        # The denominator: the Wilson UPPER bound of the competitor side. Because the
        # conservative bound is chosen on both sides, the ratio never comes out larger than
        # it is.
        # The Wilson interval is used because a raw percentage hides the uncertainty at a
        # low read count; the interval puts the uncertainty into the number. Its side effect
        # is this: the same real specificity gives a LOWER 'x' in a shallow pool, which is
        # why rows measured at different depths cannot be compared with one another.
        uy, rk = [], []
        boylar = {}
        for k in uye_kutu:
            h = self.havuz.get(k['kutu'])
            if h is None:
                continue
            p, n, boy = h.urun_veren(F, R, lo, hi, mm)
            uy.append((k['kutu'], p, n))
            for s, c in boy.items():
                boylar[s] = boylar.get(s, 0) + c
        for k in rakip_kutu:
            h = self.havuz.get(k['kutu'])
            if h is None:
                continue
            p, n, _ = h.urun_veren(F, R, lo, hi, mm)
            rk.append((k['kutu'], p, n))
        if not uy:
            return None
        uye_alt = min(engine_gateway.wilson(p, n)[0] for _, p, n in uy)
        # THE COVERAGE axis: the panel reports some targets as "13/13 bins" or "33/34 bins".
        # When a single member bin gives no product, uye_alt drops to 0 and every candidate
        # becomes indistinguishable; that is why coverage is measured SEPARATELY.
        kapsayan = [(a, p, n) for a, p, n in uy if n and p / n >= C.KAPSAM_ESIGI]
        uye_alt_k = (min(engine_gateway.wilson(p, n)[0] for _, p, n in kapsayan)
                     if kapsayan else 0.0)
        rp = sum(p for _, p, _ in rk)
        rn = sum(n for _, _, n in rk)
        # THE K-4 FIX (2026-08-03): with no competitor, the denominator used to be 1e-9;
        # dividing turned into multiplying by a billion, and universal primers came out
        # ABOVE THRESHOLD with 27 million 'x'. Now it is None -> OLCULEMEDI.
        #
        # THE POINT: on a universal target the competitor set approaches empty by the
        # definition of the membership. As the denominator goes to zero the ratio becomes
        # undefined; indeed 0.00 and 117 million could stand side by side in the same
        # column, and neither was a measurement. Putting a small constant in and carrying on
        # dividing hides that undefinedness, it does not solve it. The right thing is NOT TO
        # PRODUCE the ratio at all (None = not measurable) and to evaluate those targets by
        # coverage and by the outside the domain proportion.
        # THIS IS NOT LOWERING the threshold; the 10x threshold stands, and a different
        # quantity is measured only on the rows where the threshold cannot be applied.
        havuz_ust = engine_gateway.wilson(rp, rn)[1] if rn else None
        enkotu = None
        # The meaningful denominator threshold for the "worst single competitor bin".
        #
        # THE FIX (2026-08-02): the threshold used to be "half the largest bin".
        # That shifted WITH DEPTH: at full depth the largest bin holds 46,472 reads, so the
        # threshold became 23,236 and ONLY 1 of 10 competitor bins entered the measurement -
        # that is, "the worst bin" came to mean "the deepest bin" and the real worst
        # competitor could stay outside. Measuring the same pair with 300 reads made the
        # threshold 150 and let all 18 bins in; a spurious difference of up to 40 fold
        # between the two stages came from that.
        # It is now ABSOLUTE: the same bins enter the measurement whatever the depth.
        esik = C.ENKOTU_ASGARI_OKUMA
        for kadi, p, n in rk:
            if n < esik:
                continue
            hi_ = engine_gateway.wilson(p, n)[1]
            if enkotu is None or hi_ > enkotu[1]:
                enkotu = (kadi, hi_, p, n)
        return dict(
            olcut='<=%d uyumsuzluk + 3\' son 2 baz TAM' % mm,
            max_mm=mm,
            uye=[(a, p, n, round(100.0 * p / max(n, 1), 2)) for a, p, n in uy],
            rakip=sorted([(a, p, n, round(100.0 * p / max(n, 1), 2)) for a, p, n in rk],
                         key=lambda x: -x[3]),
            uye_alt=round(100 * uye_alt, 3),
            uye_min=round(100 * min(p / max(n, 1) for _, p, n in uy), 2),
            uye_max=round(100 * max(p / max(n, 1) for _, p, n in uy), 2),
            uye_kutu_sayisi=len(uy),
            uye_kapsam=len(kapsayan),
            uye_kapsam_pay='%d/%d' % (len(kapsayan), len(uy)),
            uye_alt_kapsayan=round(100 * uye_alt_k, 3),
            havuz='%d/%d' % (rp, rn),
            # A BUG FIX (2026-08-06, caught on a clean run): havuz_ust is deliberately
            # None when rn == 0 ("not measurable", the reasoning above). The kat_havuz /
            # kat_havuz_kapsayan lines below preserved that correctly, THIS LINE DID NOT,
            # and it crashed on the FIRST target with zero competitor reads:
            #   TypeError: unsupported operand type(s) for *: 'int' and 'NoneType'
            # Stage P failed on the 5th of 22 pairs because of it, stage T returned exit
            # code 3, and the K, D and I substages never ran. None is preserved: the field
            # is reported as "not measurable" and NO invented number IS PRODUCED.
            havuz_ust=(round(100 * havuz_ust, 3) if havuz_ust is not None else None),
            # O-8: kat_havuz can only carry a verdict at sufficient depth
            kat_havuz=(round(uye_alt / havuz_ust, 2)
                       if (havuz_ust and rn >= C.ENKOTU_ASGARI_OKUMA) else None),
            havuz_derinligi=rn,
            rakip_olculen=sum(1 for _, _, n in rk if n >= C.ENKOTU_ASGARI_OKUMA),
            rakip_toplam=len(rk),
            enkotu_kutu=enkotu[0] if enkotu else '',
            kat_enkotu=round(uye_alt / enkotu[1], 2) if enkotu and enkotu[1] > 0 else None,
            kat_havuz_kapsayan=(round(uye_alt_k / havuz_ust, 2)
                                if (havuz_ust and rn >= C.ENKOTU_ASGARI_OKUMA
                                    and kapsayan) else None),
            kat_enkotu_kapsayan=(round(uye_alt_k / enkotu[1], 2)
                                 if enkotu and enkotu[1] > 0 and kapsayan else None),
            urun_boylari=sorted(boylar.items(), key=lambda x: -x[1])[:5],
        )
