#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
alignment.py
Hizalama katmanı. Önce mappy denenir; kurulu değilse WSL'de zaten bulunan
minimap2 komut satırı aracına düşülür.

Neden gerekli: mappy saf Python değil, minimap2'nin C kaynağını derliyor ve
zlib başlıkları ya da uygun bir tekerlek yoksa kurulamıyor. Boru hattının
tek bir derleyici bağımlılığına takılmaması için aynı ölçüm iki yoldan da
yapılabilir olmalı.

İki arka uç da AYNI kaydı döndürür:
    Hiz(r_st, r_en, q_st, q_en, strand, cigar, mlen, blen, ctg)
        r_st, r_en : referans üzerinde 0 tabanlı yarı açık aralık
        q_st, q_en : sorgu üzerinde 0 tabanlı yarı açık aralık, ÖZGÜN
                     sorgu koordinatında (ters zincirde de öyle)
        strand     : 1 ya da -1
        cigar      : [(uzunluk, işlem)] listesi, işlem 0=eşleşme bloğu,
                     1=ekleme, 2=silme, 4=kırpma
        mlen       : eşleşen baz sayısı
        blen       : hizalama blok uzunluğu

Kullanım:
    import alignment
    A = hizalama.Hizalayici(referans_dizisi, preset="map-ont")
    for h in A.map(okuma):
        ...
    hizalama.ARKA_UC   -> "mappy" ya da "minimap2-cli" ya da None
"""
import os, re, shutil, subprocess, tempfile, collections

Hiz = collections.namedtuple(
    "Hiz", "r_st r_en q_st q_en strand cigar mlen blen ctg")

try:
    import mappy as _mappy
    MAPPY = True
except ImportError:
    _mappy = None
    MAPPY = False

MINIMAP2 = shutil.which("minimap2")
ARKA_UC = "mappy" if MAPPY else ("minimap2-cli" if MINIMAP2 else None)

_CIGAR_KOD = {"M": 0, "I": 1, "D": 2, "S": 4, "H": 4, "=": 0, "X": 0, "N": 3, "P": 6}
_TAMLAYICI = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn",
                           "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


def revcomp(s):
    return s.translate(_TAMLAYICI)[::-1]


def _cigar_coz(c):
    """SAM CIGAR dizesini [(uzunluk, kod)] listesine cevirir."""
    return [(int(n), _CIGAR_KOD.get(op, 0))
            for n, op in re.findall(r"(\d+)([MIDNSHP=X])", c)]


class Hizalayici(object):
    """mappy.Aligner ile ayni arayuz, iki arka ucu da destekler."""

    def __init__(self, seq=None, fn_idx_in=None, preset="map-ont",
                 ekstra=None, yalniz_birincil=True):
        self.preset = preset
        self.ekstra = list(ekstra or [])
        # mappy'de is_primary EK (supplementary, SAM 0x800) hizalamalar icin de
        # True doner, cunku minimap2'de ek hizalamanin parent==id'dir. Komut
        # satiri yolu ise 0x800'i atar. Bu yuzden mappy yolunda okuma basina
        # YALNIZ ILK kayit tutulur; minimap2 birincili once yazar. Olculdu:
        # bu suzgecle iki arka ucun 98 dosyalik ciktisi 145 baz farktan 0'a
        # dustu.
        self.yalniz_birincil = yalniz_birincil
        self._mp = None
        self._ref_yolu = None
        self._gecici = None
        if MAPPY:
            if seq is not None:
                self._mp = _mappy.Aligner(seq=seq, preset=preset)
            else:
                self._mp = _mappy.Aligner(fn_idx_in, preset=preset)
            return
        if not MINIMAP2:
            raise RuntimeError(
                "Ne mappy kurulu ne de minimap2 bulundu. Birini saglayin:\n"
                "  pip install --break-system-packages mappy\n"
                "  ya da: sudo apt-get install -y minimap2")
        if fn_idx_in:
            self._ref_yolu = fn_idx_in
        else:
            fd, yol = tempfile.mkstemp(suffix=".fa", prefix="hizref_")
            with os.fdopen(fd, "w") as fh:
                fh.write(">ref\n")
                for i in range(0, len(seq), 70):
                    fh.write(seq[i:i + 70] + "\n")
            self._ref_yolu = yol
            self._gecici = yol

    def __bool__(self):
        if MAPPY:
            return bool(self._mp)
        return bool(self._ref_yolu)

    __nonzero__ = __bool__

    def __del__(self):
        try:
            if self._gecici and os.path.exists(self._gecici):
                os.unlink(self._gecici)
        except Exception:
            pass

    # ---- tek dizi ----
    def map(self, seq, **kw):
        if MAPPY:
            for h in self._mp.map(seq, **kw):
                if not getattr(h, "is_primary", True):
                    continue
                yield Hiz(h.r_st, h.r_en, h.q_st, h.q_en, h.strand,
                          list(h.cigar), h.mlen, h.blen, h.ctg)
                if self.yalniz_birincil:
                    return          # ek (0x800) hizalamalar atlanir
            return
        for ad, hl in self.map_toplu({"q": seq}):
            for h in hl:
                yield h

    # ---- toplu: CLI arka ucunda cok daha hizli ----
    def map_toplu(self, diziler):
        """diziler: {ad: dizi}. [(ad, [Hiz, ...])] doner.
        CLI arka ucunda tek minimap2 cagrisiyla hepsi hizalanir; okuma basina
        surec baslatmak kabul edilemez derecede yavas olurdu."""
        if MAPPY:
            for ad, s in diziler.items():
                kayit = []
                for h in self._mp.map(s):
                    if not getattr(h, "is_primary", True):
                        continue
                    kayit.append(Hiz(h.r_st, h.r_en, h.q_st, h.q_en, h.strand,
                                     list(h.cigar), h.mlen, h.blen, h.ctg))
                    if self.yalniz_birincil:
                        break       # ek (0x800) hizalamalar atlanir
                yield ad, kayit
            return
        fd, sorgu = tempfile.mkstemp(suffix=".fa", prefix="hizsor_")
        try:
            with os.fdopen(fd, "w") as fh:
                for ad, s in diziler.items():
                    fh.write(">%s\n" % ad)
                    for i in range(0, len(s), 70):
                        fh.write(s[i:i + 70] + "\n")
            cmd = [MINIMAP2, "-a", "--secondary=no", "-x", self.preset,
                   "-t", "1"] + self.ekstra + [self._ref_yolu, sorgu]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError("minimap2 basarisiz: %s"
                                   % r.stderr.strip()[:300])
            uzunluk = {ad: len(s) for ad, s in diziler.items()}
            toplanan = collections.defaultdict(list)
            for line in r.stdout.splitlines():
                if line.startswith("@"):
                    continue
                p = line.split("\t")
                if len(p) < 11:
                    continue
                qad, bayrak, ctg, pos, cig = p[0], int(p[1]), p[2], int(p[3]), p[5]
                if bayrak & 0x4 or ctg == "*" or cig == "*":
                    continue
                if bayrak & 0x100 or bayrak & 0x800:
                    continue                      # ikincil ve ek hizalamalar
                strand = -1 if (bayrak & 0x10) else 1
                cigar = _cigar_coz(cig)
                r_st = pos - 1
                r_uz = sum(n for n, o in cigar if o in (0, 2, 3))
                # sorgu koordinatlari: CIGAR ters tumleyen uzerinde yurur,
                # ozgun sorgu koordinatina cevrilir
                bas_kirp = cigar[0][0] if cigar and cigar[0][1] == 4 else 0
                son_kirp = cigar[-1][0] if len(cigar) > 1 and cigar[-1][1] == 4 else 0
                q_uz = sum(n for n, o in cigar if o in (0, 1))
                if strand == 1:
                    q_st = bas_kirp
                else:
                    q_st = son_kirp
                q_en = q_st + q_uz
                mlen = None
                for alan in p[11:]:
                    if alan.startswith("NM:i:"):
                        nm = int(alan[5:])
                        eslesme_blok = sum(n for n, o in cigar if o == 0)
                        indel = sum(n for n, o in cigar if o in (1, 2))
                        mlen = eslesme_blok + indel - nm
                        break
                blok = sum(n for n, o in cigar if o in (0, 1, 2))
                if mlen is None:
                    mlen = sum(n for n, o in cigar if o == 0)
                toplanan[qad].append(
                    Hiz(r_st, r_st + r_uz, q_st, q_en, strand,
                        [(n, o) for n, o in cigar if o != 4],
                        max(0, mlen), blok, ctg))
            for ad in diziler:
                yield ad, toplanan.get(ad, [])
        finally:
            if os.path.exists(sorgu):
                os.unlink(sorgu)


def durum():
    if ARKA_UC == "mappy":
        return "hizalama arka ucu: mappy (kutuphane)"
    if ARKA_UC == "minimap2-cli":
        return "hizalama arka ucu: minimap2 komut satiri (%s)" % MINIMAP2
    return ("hizalama arka ucu YOK. Kurulum:\n"
            "   pip install --break-system-packages mappy\n"
            "   ya da: sudo apt-get install -y minimap2")


if __name__ == "__main__":
    print(durum())
    if ARKA_UC:
        ref = ("ACGTACGTTTGGCCAATTCCGGAATTCCGGTTAACCGGTTAACCGGAATTCCGG"
               "TTAACCGGTTAACCGGAATTCCGGTTAACCGGTTAACC") * 6
        A = Hizalayici(seq=ref, preset="map-ont")
        sorgu = ref[100:400]
        for h in A.map(sorgu):
            print("   ileri : r=%d-%d q=%d-%d strand=%d mlen=%d blen=%d"
                  % (h.r_st, h.r_en, h.q_st, h.q_en, h.strand, h.mlen, h.blen))
        for h in A.map(revcomp(sorgu)):
            print("   ters  : r=%d-%d q=%d-%d strand=%d mlen=%d blen=%d"
                  % (h.r_st, h.r_en, h.q_st, h.q_en, h.strand, h.mlen, h.blen))
