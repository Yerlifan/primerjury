# -*- coding: utf-8 -*-
"""
brute_force.py, an INDEPENDENT reference implementation. NO seeding, NO shortcuts.

Every start position is tried one at a time. It is slow; its purpose is not to be
fast but to produce a right answer INDEPENDENT of `read_engine.py`, against which
that engine's correctness can be tested. It shares no common code: even the IUPAC
table and the reverse complement are written out again, so that the same fault
cannot sit in both.

Usage: it is called by engine_test.py.

"""
# -------------------------------------------------------------------------
# brute_force.py, a reference implementation written to test read_engine.py for
#                  correctness; no seeding and no shortcuts.
#
# INPUT  : it takes a sequence and a primer directly; it reads no file. The caller
#          supplies the reads (engine_test.py, self_test.py or
#          independent_check.py).
# OUTPUT : it writes no file. yerler() returns a [(start, mismatches)] list,
#          urun_boyu() a product length or None, and kutu_pcr() the pair
#          (with_product, total).
# CALLED BY: engine_test.py and independent_check.py (tests run by hand) and, through
#          engine_gateway.py, self_test.py; that is verification/full_chain.py key 8
#          and the self test step at the head of every measuring key.
#
# WHY IT WAS WRITTEN SEPARATELY: this file is the evidence for the claim that the
# pigeonhole seeding is LOSSLESS. It tries every start position one at a time, so it
# rests on no seeding assumption. Even the IUPAC table and the reverse complement
# are written out again; had common code been shared, THE SAME fault would sit in
# both implementations and the comparison would prove nothing.
# -------------------------------------------------------------------------

IUP = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
       'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC',
       'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}

_C = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def rc(s):
    return ''.join(_C.get(c, 'N') for c in reversed(s))


def yerler(seq, primer, max_mm, son2=True, uc5=False):
    """primer'in seq uzerindeki TUM baglanma baslangiclari - kaba kuvvet.
    uc5=False -> 3' kritik uc primerin sonunda; uc5=True -> basinda."""
    primer = primer.upper()
    L = len(primer)
    ok = [set(IUP.get(c, 'ACGT')) for c in primer]
    kritik = (0, 1) if uc5 else (L - 1, L - 2)
    out = []
    for st in range(len(seq) - L + 1):
        if son2:
            kotu = False
            for k in kritik:
                if seq[st + k] not in ok[k]:
                    kotu = True
                    break
            if kotu:
                continue
        mm = 0
        kotu = False
        for k in range(L):
            if seq[st + k] not in ok[k]:
                mm += 1
                if mm > max_mm:
                    kotu = True
                    break
        if not kotu:
            out.append((st, mm))
    return out


def urun_boyu(seq, F, R, max_mm, lo=40, hi=600, son2=True):
    """Bu dizide urun varsa boyunu dondur, yoksa None."""
    F = F.upper(); R = R.upper()
    a = yerler(seq, F, max_mm, son2, uc5=False)
    if not a:
        return None
    b = yerler(seq, rc(R), max_mm, son2, uc5=True)
    if not b:
        return None
    for i, _ in a:
        for j, _ in b:
            n = j + len(R) - i
            if lo <= n <= hi and j >= i + len(F):
                return n
    return None


def kutu_pcr(okuma_listesi, F, R, lo=40, hi=600, max_mm=1, son2=True):
    """Returns: (with_product, total)"""
    pos = 0
    for s in okuma_listesi:
        for seq in (s, rc(s)):
            if urun_boyu(seq, F, R, max_mm, lo, hi, son2) is not None:
                pos += 1
                break
    return pos, len(okuma_listesi)
