#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
alignment.py
The alignment layer. mappy is tried first; if it is not installed it falls back on
the minimap2 command line tool, which is present in WSL anyway.

Why it is needed: mappy is not pure Python, it compiles minimap2's C source, and it
cannot be installed without zlib headers or a suitable wheel. So that the pipeline
does not hang on a single compiler dependency, the same measurement has to be
available two ways.

Both backends return THE SAME record:
    Hiz(r_st, r_en, q_st, q_en, strand, cigar, NM, mlen, blen)

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
    'Turns a SAM CIGAR string into a list of (length, code).'
    return [(int(n), _CIGAR_KOD.get(op, 0))
            for n, op in re.findall(r"(\d+)([MIDNSHP=X])", c)]


class Hizalayici(object):
    """The same interface as mappy.Aligner, supporting both backends."""

    def __init__(self, seq=None, fn_idx_in=None, preset="map-ont",
                 ekstra=None, yalniz_birincil=True):
        self.preset = preset
        self.ekstra = list(ekstra or [])
        # In mappy, is_primary returns True for SUPPLEMENTARY alignments (SAM 0x800) as
        # well, because in minimap2 a supplementary alignment's parent equals its id. The
        # command line route, on the other hand, drops 0x800. So on the mappy route ONLY THE
        # FIRST record per read is kept; minimap2 writes the primary one first. Measured:
        # with this filter the difference between the two backends over 98 files fell from
        # 145 bases to 0.
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

    # ---- bulk: much faster on the CLI backend ----
    def map_toplu(self, diziler):
        """diziler: {name: sequence}. Returns [(name, [Hiz, ...])].
        On the CLI backend everything is aligned with a single minimap2 call; starting a
        process per read would be unacceptably slow.

        """
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
        return 'the alignment back end: mappy, as a library'
    if ARKA_UC == "minimap2-cli":
        return 'the alignment back end: the minimap2 command line (%s)' % MINIMAP2
    return ("""THERE IS NO alignment back end. To install one:
               pip install --break-system-packages mappy
               or: sudo apt-get install -y minimap2""")


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
