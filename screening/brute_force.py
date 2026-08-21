# -*- coding: utf-8 -*-
"""
brute_force.py - BAGIMSIZ referans uygulama. Tohum YOK, kisayol YOK.

Her baslangic pozisyonu tek tek denenir. Yavastir; amaci hizli olmak degil,
`read_engine.py`'nin dogrulugunu sinayacak, ondan BAGIMSIZ bir dogru cevap
uretmektir. Ortak kod paylasmaz - IUPAC tablosu ve ters tumleyen bile ayrica
yazilmistir ki ikisinde ayni hata olmasin.

Kullanim: engine_test.py tarafindan cagrilir.
"""
# ---------------------------------------------------------------------------
# brute_force.py — read_engine.py'nin dogrulugunu sinamak icin yazilmis,
#                  tohumsuz ve kisayolsuz referans uygulama.
#
# GIRDI  : dogrudan dizi ve primer alir; dosya okumaz. Cagiran taraf okumalari
#          verir (engine_test.py, self_test.py ya da independent_check.py).
# CIKTI  : dosyaya yazmaz. yerler() [(baslangic, uyumsuzluk)] listesi,
#          urun_boyu() urun boyu ya da None, kutu_pcr() (urun_veren, toplam)
#          ikilisi dondurur.
# CAGRAN : engine_test.py ve independent_check.py (elle calistirilan sinamalar)
#          ile engine_gateway.py uzerinden self_test.py - yani verification/full_chain.py
#          tusu 8 ve her olcum tusunun basindaki kendini sinama adimi.
#
# NEDEN AYRI YAZILDI: bu dosya, guvercin yuvasi tohumlamasinin KAYIPSIZ oldugu
# iddiasinin kanitidir. Her baslangic pozisyonunu tek tek dener, yani hicbir
# tohum varsayimina dayanmaz. IUPAC tablosu ve ters tumleyen bile ayrica
# yazilmistir; ortak kod paylasilsaydi iki uygulamada AYNI hata bulunur ve
# karsilastirma hicbir sey kanitlamazdi.
# ---------------------------------------------------------------------------

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
    """Donus: (urun_veren, toplam)"""
    pos = 0
    for s in okuma_listesi:
        for seq in (s, rc(s)):
            if urun_boyu(seq, F, R, max_mm, lo, hi, son2) is not None:
                pos += 1
                break
    return pos, len(okuma_listesi)
