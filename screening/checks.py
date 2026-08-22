# -*- coding: utf-8 -*-
"""Resilience against interruption: each target's result is written to disk as it
finishes.

When the program is closed and opened again the finished targets are SKIPPED and a
half finished target is started over (a half result does not enter the report). The
global scan has a checkpoint of its own as well (per piece).

"""
# -------------------------------------------------------------------------
# checks.py, resilience against interruption: it writes each target's result to
#              disk and skips the finished targets when the run restarts.
#
# INPUT  : the SCREENING_RESULT/kontrol/hedef_*.json files; the run's settings
#          fingerprint is filled in by __main__.main() through ayar_kur() (the read
#          count, the light mode, the candidate ceiling).
# OUTPUT : kontrol/hedef_<target>.json files (written as .gecici first and put in
#          place atomically with os.replace). oku(), hepsi() and bitti_mi() return
#          what is on disk; sifirla() empties the kontrol directory.
# CALLED BY: __main__.py (as each target finishes), run_all.py, panel_measurement.py
#          and membership_check.py manage their own checkpoints through this
#          module's paths and its settings comparison. That is keys 1, 2, 3, 4, 5,
#          6, 7 and 9.
#
# THE SETTINGS FINGERPRINT IS CRITICAL: the same target gives a different result at
# 300 reads and at full depth. Had a checkpoint with no settings record, or one
# that does not match, been reused, the run would silently have reported the result
# of the old settings as the result of the new ones. That is why ayar_uyuyor() does
# not trust old files and returns False.
# -------------------------------------------------------------------------
import os, json, time
from . import config as C

# A checkpoint is reused ONLY if it was produced under the same settings.
# Otherwise (left over from a trial run made with a low read count, for instance)
# it silently returns the wrong result. The AYAR variable is filled in by main() at
# the start of every run.
AYAR = {}


def ayar_kur(**kw):
    AYAR.clear()
    AYAR.update({k: v for k, v in kw.items() if v is not None})


def ayar_uyuyor(veri):
    kayit = (veri or {}).get('_ayar')
    if kayit is None:
        return False          # an old file with no settings record, do not trust it
    return kayit == AYAR


def hazirla():
    os.makedirs(C.CIKTI, exist_ok=True)
    os.makedirs(C.KONTROL, exist_ok=True)
    os.makedirs(C.ONBELLEK, exist_ok=True)


def _yol(hedef):
    ad = ''.join(ch if ch.isalnum() else '_' for ch in hedef)
    return os.path.join(C.KONTROL, 'hedef_%s.json' % ad)


def bitti_mi(hedef):
    v = oku(hedef)
    return v is not None and ayar_uyuyor(v)


def yaz(hedef, veri):
    hazirla()
    veri = dict(veri)
    veri['_zaman'] = time.strftime('%Y-%m-%d %H:%M:%S')
    veri['_ayar'] = dict(AYAR)
    gec = _yol(hedef) + '.gecici'
    with open(gec, 'w', encoding='utf-8') as fh:
        json.dump(veri, fh, ensure_ascii=False, indent=1, default=str)
    os.replace(gec, _yol(hedef))       # atomic: no half written file is left behind


def oku(hedef):
    p = _yol(hedef)
    if not os.path.exists(p):
        return None
    try:
        return json.load(open(p, encoding='utf-8'))
    except Exception:
        return None      # a corrupted file: do not try to delete it (the mount may be read only)


def hepsi():
    hazirla()
    out = []
    for f in sorted(os.listdir(C.KONTROL)):
        if f.startswith('hedef_') and f.endswith('.json'):
            try:
                out.append(json.load(open(os.path.join(C.KONTROL, f), encoding='utf-8')))
            except Exception:
                pass
    return out


def sifirla():
    hazirla()
    for f in os.listdir(C.KONTROL):
        try:
            os.remove(os.path.join(C.KONTROL, f))
        except OSError:
            pass          # if it cannot be deleted that is fine: the settings fingerprint ignores it anyway
