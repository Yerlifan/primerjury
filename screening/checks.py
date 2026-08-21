# -*- coding: utf-8 -*-
"""Kesinti dayanikliligi: her hedefin sonucu bitince diske yazilir.

Program kapatilip yeniden acilinca bitmis hedefler ATLANIR, yarim kalan hedef
bastan alinir (yarim sonuc rapora girmez). Kuresel taramanin kendi ic
checkpoint'i ayrica vardir (parca bazinda).
"""
# ---------------------------------------------------------------------------
# checks.py — kesinti dayanikliligi: her hedefin sonucunu diske yazar ve
#              yeniden baslatildiginda bitmis hedefleri atlatir.
#
# GIRDI  : KAPSAMLI_ARAMA_SONUC/kontrol/hedef_*.json dosyalari; kosunun ayar
#          parmak izi __main__.main() tarafindan ayar_kur() ile doldurulur
#          (okuma sayisi, hafif mod, aday ust siniri).
# CIKTI  : kontrol/hedef_<hedef>.json dosyalari (once .gecici olarak yazilip
#          os.replace ile atomik yerine konur). oku(), hepsi() ve bitti_mi()
#          diskteki veriyi dondurur; sifirla() kontrol klasorunu bosaltir.
# CAGRAN : __main__.py (her hedef bitince), run_all.py, panel_measurement.py ve
#          membership_check.py kendi kontrol noktalarini bu modulun yollari ve
#          ayar karsilastirmasi uzerinden yonetir. Yani tuslar 1, 2, 3, 4, 5,
#          6, 7 ve 9.
#
# AYAR PARMAK IZI KRITIKTIR: ayni hedef 300 okumayla ve tam derinlikle
# olculdugunde farkli sonuc verir. Ayar kaydi olmayan ya da uyusmayan bir
# kontrol noktasi yeniden kullanilsaydi, kosu sessizce eski ayarin sonucunu
# yeni ayarin sonucu gibi rapor ederdi. Bu yuzden ayar_uyuyor() eski dosyalara
# guvenmez ve False dondurur.
# ---------------------------------------------------------------------------
import os, json, time
from . import config as C

# Bir kontrol noktasi, YALNIZ ayni ayarlarla uretilmisse yeniden kullanilir.
# Aksi halde (ornegin dusuk okuma sayisiyla yapilmis bir deneme kosusundan
# kalmissa) sessizce yanlis sonuc dondurur. AYAR degiskeni her kosu basinda
# main() tarafindan doldurulur.
AYAR = {}


def ayar_kur(**kw):
    AYAR.clear()
    AYAR.update({k: v for k, v in kw.items() if v is not None})


def ayar_uyuyor(veri):
    kayit = (veri or {}).get('_ayar')
    if kayit is None:
        return False          # ayar kaydi olmayan eski dosya - guvenme
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
    os.replace(gec, _yol(hedef))       # atomik: yarim dosya kalmaz


def oku(hedef):
    p = _yol(hedef)
    if not os.path.exists(p):
        return None
    try:
        return json.load(open(p, encoding='utf-8'))
    except Exception:
        return None      # bozuk dosya: silmeye calisma (salt-okunur baglama olabilir)


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
            pass          # silinemiyorsa sorun degil: ayar parmak izi zaten yok sayar
