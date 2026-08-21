# -*- coding: utf-8 -*-
"""A2 OLCUMU: boy tabanli ayrim ile taksonomik ayrim yan yana.

NE SORUYOR
----------
Katman 2 (yerel veritabani) bir urunun hedefin ICINDEN mi DISINDAN mi geldigini
eskiden yalniz BOYA bakarak ayiriyordu. A2 ile TAKSONOMIK ayrim eklendi.
Bu betik ikisini AYNI tarama uzerinde olcer ve farki gosterir.

Hicbir panel dosyasina YAZMAZ. Yalniz okur ve tablo basar.

NEDEN ONEMLI
------------
D-12'de (katman 3, MFEprimer) olculmustu: "hedef disi" sayilan 1.605 amplikonun
%95,7'si hedef kladin KENDI ICINDENdi, yalniz boyu farkliydi. Ayni yanilginin
katman 2'de de olup olmadigi bu betikle olculur.

KOSMA
-----
    python3 DENETIM_SINAMALARI/olcum_A2_oncesi_sonrasi.py --kucuk
    python3 DENETIM_SINAMALARI/olcum_A2_oncesi_sonrasi.py --vtb SILVA_138.2_SSURef_NR99.fasta
    python3 DENETIM_SINAMALARI/olcum_A2_oncesi_sonrasi.py --hedef Bakteri_universal

--kucuk : yalniz kucuk RefSeq kumeleri (~75 MB, dakikalar). Kapsam degil,
          YONTEMIN CALISTIGININ kaniti ve ilk buyukluk mertebesi.
Tam kapsam icin SILVA/UNITE gibi buyuk kumeler ayrica verilmelidir; suresi
saatlerdir ve bu betik onu SESSIZCE yapmaz.
"""
from __future__ import print_function
import argparse
import csv
import hashlib
import io
import json
import os
import re
import sys
import time

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, KOK)

from KAPSAMLI_ARAMA import kuresel_tarama as KT          # noqa: E402
from KAPSAMLI_ARAMA import taksonomi as TX               # noqa: E402
import KURTARMA.mfe_katmani as MK                        # noqa: E402

KUCUK = ['archaea.16S.fna', 'fungi.18SrRNA.fna', 'fungi.28SrRNA.fna',
         'fungi.ITS.fna', 'bacteria.16S.fna']

KAYNAK = os.path.join(KOK, 'DOGRULAMA_SONUC', 'dogrulama_uc_sutun.tsv')
URUN_ALT, URUN_UST, BOY_TOL = 70, 400, 10


def ciftleri_oku():
    if not os.path.exists(KAYNAK):
        sys.stderr.write('cift kaynagi yok: %s\n' % KAYNAK)
        sys.exit(2)
    out = []
    with io.open(KAYNAK, encoding='utf-8') as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith('#')), delimiter='\t'):
            h = (r.get('hedef') or '').strip()
            F = (r.get('F') or '').strip().upper()
            R = (r.get('R') or '').strip().upper()
            bp = (r.get('urun_bp') or '').strip()
            if h and F and R:
                out.append(dict(hedef=h, F=F, R=R,
                                bek=int(bp) if bp.isdigit() else None))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--kucuk', action='store_true')
    p.add_argument('--vtb', nargs='*', default=None)
    p.add_argument('--hedef', default=None)
    a = p.parse_args()

    ciftler = ciftleri_oku()
    if a.hedef:
        ciftler = [c for c in ciftler if a.hedef.lower() in c['hedef'].lower()]
    if not ciftler:
        sys.stderr.write('cift bulunamadi\n'); return 2

    klad = MK.klad_tablosu(KOK)
    yok = [c['hedef'] for c in ciftler if c['hedef'] not in klad]
    if yok:
        print('UYARI: hedef_klad.tsv\'de tanimi olmayan %d hedef - bunlarda'
              ' taksonomik ayrim YAPILAMAZ:' % len(yok))
        for h in yok:
            print('   - %s' % h)
        print()

    def sinifla(aday_ad, baslik, db_ad):
        if aday_ad not in klad:
            return 'bilinmiyor'
        alan, jetonlar, _k = klad[aday_ad]
        return TX.sinifla(baslik, db_ad, jetonlar, alan)

    vtbler = a.vtb if a.vtb else (KUCUK if a.kucuk else KUCUK)
    refdb = os.path.join(KOK, 'REFERANS_DB')

    adaylar = [dict(ad=c['hedef'], F=c['F'], R=c['R'], lo=URUN_ALT, hi=URUN_UST)
               for c in ciftler]
    bek = {c['hedef']: c['bek'] for c in ciftler}

    top = {c['hedef']: dict(urun=0, ayni_boyda=0, boy_disi=0,
                            sinif={k: 0 for k in KT.SINIFLAR}, tarandi=0)
           for c in ciftler}

    # ------------------------------------------------------------------
    # KESINTIYE DAYANIKLILIK
    #
    # Ilk surumde tara() 'durum_yolu=None' ile cagriliyordu, yani kontrol
    # noktasi KAPALIYDI. Bu bir hataydi: kuresel_tarama zaten parca parca
    # kontrol noktasi tutar (katman1_yerel de onu kullanir) ve bu projenin
    # kurali uzun kosunun kesintiye dayanikli olmasidir. Kapali kontrol
    # noktasiyla 28 dakikalik SILVA taramasi bir kapanmada tumden kayboluyordu.
    #
    # Iki duzeyli koruma:
    #   1) tara() icin parca duzeyi kontrol noktasi (durum_yolu)
    #   2) HER VERITABANI bitince kismi sonuc diske yazilir; kosu kesilse bile
    #      o ana kadarki olcum durur ve okunabilir.
    #
    # Kontrol noktasi anahtari aday DIZILERINI de tasir: yalniz adlarla
    # muhurlenirse dizi degistiginde eski tarama sessizce geri gelir
    # (dogrulama_turu'nda 2026-08-10'da tam bu hata olculmustu).
    kn_dizin = os.path.join(KOK, 'KAPSAMLI_ARAMA_SONUC', 'kontrol', 'A2_olcum')
    try:
        os.makedirs(kn_dizin)
    except OSError:
        pass
    imza = hashlib.md5('|'.join(sorted(
        '%s>%s<%s' % (x['ad'], x['F'], x['R']) for x in adaylar)
    ).encode('utf-8')).hexdigest()[:10]
    kismi_yol = os.path.join(kn_dizin, 'kismi_sonuc_%s.json' % imza)

    tamamlanan = []
    if os.path.exists(kismi_yol):
        try:
            _k = json.load(io.open(kismi_yol, encoding='utf-8'))
            if _k.get('imza') == imza:
                top = _k['top']
                tamamlanan = _k.get('tamamlanan', [])
                print('onceki kosudan devam: %d veritabani zaten olculmus (%s)'
                      % (len(tamamlanan), ', '.join(tamamlanan)))
        except Exception as e:
            print('kismi sonuc okunamadi (%s) - bastan olculuyor' % type(e).__name__)

    def kismi_yaz():
        gecici = kismi_yol + '.tmp'
        with io.open(gecici, 'w', encoding='utf-8') as fh:
            fh.write(json.dumps(dict(imza=imza, tamamlanan=tamamlanan, top=top),
                                ensure_ascii=False))
        if os.path.exists(kismi_yol):
            os.remove(kismi_yol)
        os.rename(gecici, kismi_yol)

    print('%d cift, %d veritabani' % (len(ciftler), len(vtbler)))
    for d in vtbler:
        if d in tamamlanan:
            print('  ATLANDI (zaten olculmus): %s' % d)
            continue
        yol = os.path.join(refdb, d)
        if not os.path.exists(yol):
            print('  ATLANDI (yok): %s' % d)
            continue
        t0 = time.time()
        print('  taraniyor: %-34s (%s)' % (d, '%.0f MB' % (os.path.getsize(yol) / 1e6)),
              end='', flush=True)
        dy = os.path.join(kn_dizin, 'tarama_%s_%s.pkl'
                          % (re.sub(r'\W+', '_', d), imza))
        res = KT.tara(adaylar, db=yol, durum_yolu=dy, siniflandirici=sinifla)
        print('  %.0f sn' % (time.time() - t0))
        for h, r in res.items():
            if r.get('hata'):
                continue
            top[h]['tarandi'] += 1
            top[h]['urun'] += r.get('urun', 0)
            b = bek.get(h)
            for sz, n in (r.get('boy') or {}).items():
                if b is not None and abs(int(sz) - b) <= BOY_TOL:
                    top[h]['ayni_boyda'] += n
                else:
                    top[h]['boy_disi'] += n
            for k, v in (r.get('sinif') or {}).items():
                top[h]['sinif'][k] += v
        tamamlanan.append(d)
        kismi_yaz()
        print('     kismi sonuc yazildi: %s' % os.path.relpath(kismi_yol, KOK))

    print()
    print('=' * 118)
    print('%-40s %9s | %11s | %9s %9s %9s %9s | %11s' %
          ('hedef', 'tum', 'BOY: disi', 'a klad-ici', 'ao organel',
           'b ayni-alan', 'c fark-alan', 'TAKSON: b+c'))
    print('-' * 118)
    degisen = []
    for c in ciftler:
        h = c['hedef']; t = top[h]
        if not t['tarandi']:
            continue
        s = t['sinif']
        kd = s['b'] + s['c']
        print('%-40s %9d | %11d | %9d %9d %9d %9d | %11d%s' %
              (h[:40], t['urun'], t['boy_disi'], s['a'], s['ao'], s['b'], s['c'], kd,
               '  (bilinmiyor %d)' % s['bilinmiyor'] if s['bilinmiyor'] else ''))
        if t['boy_disi'] != kd:
            degisen.append((h, t['boy_disi'], kd))
    print('-' * 118)

    print()
    print('=== FARK ===')
    if not degisen:
        print('Iki olcut ayni sonucu verdi (bu veritabani alt kumesinde).')
    else:
        print('%-40s %14s %14s %s' % ('hedef', 'BOY olcutu', 'TAKSON olcutu', 'etki'))
        print('-' * 100)
        for h, b, k in degisen:
            if b > 0 and k == 0:
                etki = 'RISKLI -> TEMIZ  (boy olcutu YANLIS ALARM veriyordu)'
            elif b == 0 and k > 0:
                etki = 'TEMIZ -> RISKLI  (boy olcutu GERCEK caprazi KACIRIYORDU)'
            elif k < b:
                etki = 'capraz sayisi %d -> %d (azaldi)' % (b, k)
            else:
                etki = 'capraz sayisi %d -> %d (artti)' % (b, k)
            print('%-40s %14d %14d  %s' % (h[:40], b, k, etki))

    print()
    print('NOT: bu olcum yalniz su veritabanlarinda yapildi: %s' % ', '.join(vtbler))
    print('     Tam kapsam SILVA/UNITE/PR2/ROD gerektirir ve SAATLER surer.')
    print('     Buradaki sayilar KAPSAM DEGIL, yontemin calistiginin kanitidir.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
