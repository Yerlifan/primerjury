# -*- coding: utf-8 -*-
"""Mevcut olcum kodunu ICE AKTARAN tek kapi.

Bu dosya hicbir olcum algoritmasini yeniden yazmaz. Projede zaten calisan ve
panelin sayilarini uretmis olan modulleri bulur ve disari acar:

  engine/ispcr.py   -> find_sites / amplify / scan_file / rc / clean / encode
  engine/okuma.py   -> Sonda / okumalar / kutu_pcr  (ham okuma hatti)
  engine/tarayici.py-> Havuz  (3' tohumlu hizli havuz)
  engine/cift.py    -> urunler (havuz uzerinde cift taramasi)

Betikler oturuma bagli mutlak yollar (sys.path.insert('/tmp/wk2/...')) iceriyor;
ice aktarmadan once o yollar gecerli hale getirilir.
"""
# ---------------------------------------------------------------------------
# motor.py — projede zaten calisan olcum kodunu ice aktaran tek kapi; hicbir
#            algoritmayi yeniden yazmaz, yalnizca disari acar.
#
# GIRDI  : yapilandirma.BETIK_YOLLARI altindaki engine/ispcr.py,
#          engine/okuma.py, engine/tarayici.py ve
#          engine/cift.py dosyalari; ayrica bu paketin icindeki
#          okuma_motoru.py ve kaba_kuvvet.py.
# CIKTI  : dosyaya yazmaz. Modul duzeyinde fonksiyon ve modul nesneleri acar:
#          rc, clean, encode, read_fasta, find_sites, amplify, scan_file,
#          IUPAC, wilson(), urun_var() ve surum_bilgisi().
# CAGRAN : paketteki butun olcum modulleri; disaridan verification/kurtarma_turu.py
#          (tus K) ve protocol/tek_protokol_olc.py (tus P). Yani her olcum
#          tusunda yuklenir.
#
# NEDEN ICE AKTARIM, YENIDEN YAZIM DEGIL: panelin yayimlanmis sayilari bu
# betiklerle uretildi. Ayni islevi yeniden yazmak, eski ve yeni sayilari
# karsilastirilamaz hale getirirdi. Tek istisna ham okuma hattidir: panelin
# okuma.py/Sonda motoru site kaciriyor, bu yuzden numune olcumlerinde otorite
# okuma_motoru.py'dir; okuma.py yalniz "duzeltme ne kadar fark yaratti"
# sorusuna cevap uretmek icin duruyor.
# ---------------------------------------------------------------------------
import os, sys, importlib.util, math

from . import yapilandirma as C


def _yol_hazirla():
    for p in C.BETIK_YOLLARI:
        if p not in sys.path:
            sys.path.insert(0, p)


def _yukle(ad, dosya):
    """Bir .py dosyasini modul olarak yukle (paket disindan)."""
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

_p = _bul('okuma.py')
okuma = _yukle('okuma', _p) if _p else None

_p = _bul('tarayici.py')
tarayici = _yukle('tarayici', _p) if _p else None

_p = _bul('cift.py')
cift = _yukle('cift', _p) if _p else None

# ---- dogrudan yeniden kullanilan fonksiyonlar ----
rc = ispcr.rc
clean = ispcr.clean
encode = ispcr.encode
read_fasta = ispcr.read_fasta
find_sites = ispcr.find_sites
amplify = ispcr.amplify
scan_file = ispcr.scan_file
IUPAC = ispcr.IUPAC

if okuma is not None:
    # okuma.py'nin okuma uzunluk filtresi duzeltilmis degeri tasisin
    okuma.MINL, okuma.MAXL = C.NUMUNE_OKUMA_MIN, C.NUMUNE_OKUMA_MAX


# ---- DUZELTILMIS ham okuma motoru (bu paketin icinde, "Okuma motoru hatasi"
#      oturumunda uretildi). Panelin okuma.py/Sonda motoru site KACIRIYOR;
#      numune olcumlerinde otorite budur.
_pk = os.path.dirname(os.path.abspath(__file__))
okuma_motoru = _yukle('okuma_motoru', os.path.join(_pk, 'okuma_motoru.py'))
try:
    kaba_kuvvet = _yukle('kaba_kuvvet', os.path.join(_pk, 'kaba_kuvvet.py'))
except Exception:
    kaba_kuvvet = None


def wilson(k, n, z=1.96):
    """numune_olc.py / okuma_pcr.py ile birebir ayni Wilson araligi."""
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
