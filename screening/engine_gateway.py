# -*- coding: utf-8 -*-
"""Engine gateway, the single import point for the low-level sequence code.

Loads ispcr / reads / scanner / pair from engine/ and re-exports the functions
the rest of the package uses (rc, clean, encode, read_fasta, find_sites,
amplify, scan_file, IUPAC).

Everything goes through here on purpose: one import point means one place where
a version or path problem can surface, instead of a dozen modules each loading
their own copy.

WHAT COMES FROM WHERE
  engine/ispcr.py   -> find_sites / amplify / scan_file / rc / clean / encode
  engine/reads.py   -> Sonda / okumalar / kutu_pcr  (the raw read pipeline)
  engine/scanner.py -> Havuz  (the fast 3' seeded pool)
  engine/pair.py    -> urunler (the pair scan over the pool)

Some scripts hold absolute paths from the session they were written in
(sys.path.insert('/tmp/wk2/...')); those paths are made valid before the import.

"""
# -------------------------------------------------------------------------
# engine_gateway.py, the single door that imports the measurement code already
#            running in this project; it rewrites no algorithm, it only exposes
#            them.
#
# INPUT  : the files engine/ispcr.py, engine/reads.py, engine/scanner.py and
#          engine/pair.py under yapilandirma.BETIK_YOLLARI; also read_engine.py and
#          brute_force.py inside this package.
# OUTPUT : it writes no file. It exposes functions and module objects at module
#          level: rc, clean, encode, read_fasta, find_sites, amplify, scan_file,
#          IUPAC, wilson() and urun_bilgisi().
# CALLED BY: every measurement module in the package; from outside,
#          verification/recovery_round.py (key K) and
#          protocol/single_protocol_measure.py (key P). So it is loaded on every
#          measuring key.
#
# WHY IMPORT RATHER THAN REWRITE: the panel's published numbers were produced with
# these scripts. Rewriting the same function would make the old and the new numbers
# incomparable. The one exception is the raw read pipeline: the panel's
# reads.py/Sonda engine misses sites, so in the sample measurements the authority is
# read_engine.py; reads.py stays only to answer the question "how much difference
# did the fix make".
# -------------------------------------------------------------------------
import os, sys, importlib.util, math

from . import config as C


def _yol_hazirla():
    for p in C.BETIK_YOLLARI:
        if p not in sys.path:
            sys.path.insert(0, p)


def _yukle(ad, dosya):
    """Load a .py file as a module (from outside the package)."""
    if ad in sys.modules:
        return sys.modules[ad]
    spec = importlib.util.spec_from_file_location(ad, dosya)
    m = importlib.util.module_from_spec(spec)
    sys.modules[ad] = m
    spec.loader.exec_module(m)
    return m


def _bul(*adaylar):
    for kok in C.BETIK_YOLLARI:
        for a in adaylar:
            p = os.path.join(kok, a)
            if os.path.exists(p):
                return p
    return None


_yol_hazirla()

_p = _bul('ispcr.py')
if _p is None:
    raise SystemExit('HATA: ispcr.py bulunamadi. engine klasoru yerinde mi?')
ispcr = _yukle('ispcr', _p)

_p = _bul('reads.py')
okuma = _yukle('okuma', _p) if _p else None

_p = _bul('scanner.py')
tarayici = _yukle('tarayici', _p) if _p else None

_p = _bul('pair.py')
cift = _yukle('cift', _p) if _p else None

# ---- the functions reused directly ----
rc = ispcr.rc
clean = ispcr.clean
encode = ispcr.encode
read_fasta = ispcr.read_fasta
find_sites = ispcr.find_sites
amplify = ispcr.amplify
scan_file = ispcr.scan_file
IUPAC = ispcr.IUPAC

if okuma is not None:
    # let reads.py's read length filter carry the corrected value
    okuma.MINL, okuma.MAXL = C.NUMUNE_OKUMA_MIN, C.NUMUNE_OKUMA_MAX


# ---- THE CORRECTED raw read engine (inside this package, produced in the
#      "the read engine fault" session). The panel's reads.py/Sonda engine MISSES
#      sites; in the sample measurements this is the authority.
_pk = os.path.dirname(os.path.abspath(__file__))
okuma_motoru = _yukle('okuma_motoru', os.path.join(_pk, 'read_engine.py'))
try:
    kaba_kuvvet = _yukle('kaba_kuvvet', os.path.join(_pk, 'brute_force.py'))
except Exception:
    kaba_kuvvet = None


def wilson(k, n, z=1.96):
    """Exactly the same Wilson interval as numune_olc.py and okuma_pcr.py."""
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    s = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - s), min(1.0, c + s))


def urun_var(seq, F, R, lo, hi, max_mm=1):
    """Tek bir dizide urun olusuyor mu (ispcr.amplify olcutu, iki yon)."""
    for s in (seq, rc(seq)):
        if ispcr.amplify(s, F, R, max_mm=max_mm, lo=lo, hi=hi):
            return True
    return False


def surum_bilgisi():
    return dict(
        ispcr=os.path.abspath(ispcr.__file__),
        okuma=os.path.abspath(okuma.__file__) if okuma else None,
        tarayici=os.path.abspath(tarayici.__file__) if tarayici else None,
        cift=os.path.abspath(cift.__file__) if cift else None,
        okuma_motoru=os.path.abspath(okuma_motoru.__file__),
        okuma_motoru_surum=getattr(okuma_motoru, '__version__', '?'),
        kaba_kuvvet=os.path.abspath(kaba_kuvvet.__file__) if kaba_kuvvet else None,
    )
