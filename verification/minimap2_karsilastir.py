# -*- coding: utf-8 -*-
"""
minimap2_karsilastir.py - iki hizalayiciyi AYNI girdiyle yan yana olcer.

Sorusu tek: minimap2, saf Python motorunun bulduğu kimligi ayni sekilde buluyor
mu ve ne kadar hizli. Cevap "evet" ise minimap2 varsayilan yapilabilir; "hayir"
ise mevcut motorda kalinir. Hizli olan DOGRU SAYILMAZ.
"""
# ---------------------------------------------------------------------------
# minimap2_karsilastir.py
#
# GIRDI  : REFERANS_DB/ altindaki gercek veritabanlari ve gercek kutu
#          konsensusleri (consensus sequences/ ya da konsensus_kanonik/)
# CIKTI  : MINIMAP2_KARSILASTIRMA.md (kok dizinde) ve ayni adin .tsv'si
# CAGRAN : elle. Menuye BAGLANMADI, cunku bu bir karar belgesi uretir,
#          zincirin bir asamasi degildir.
#
# NEYI OLCER
#   1) AYNI EN IYI ISABET MI. Iki motor ayni kaydi mi birinci sirada getiriyor.
#      Bu en onemli olcut: kimlik karari en iyi isabetten cikar.
#   2) KIMLIK YUZDESI SAPMASI. Ayni kayit icin iki motorun verdigi yuzde farki.
#   3) SIRALAMA KORUNUYOR MU. Ilk 5 isabetin kumesi ne kadar ortusuyor.
#   4) HIZLANMA. Ayni is icin gecen sure orani.
#
# NEDEN BU DORT OLCUT
# Yalniz hiza bakmak yanlis olurdu: hizli ama farkli cevap veren bir motor
# bize zaman kazandirmaz, yanlis primer siparis ettirir. Yalniz en iyi isabete
# bakmak da yetmez: siralama bozulursa "en az iki bagimsiz veritabani uyusmali"
# kurali baska kayitlar uzerinden calisir ve hukum degisebilir.
#
# AYRILIK CIKARSA
# Rapor AYRILAN satirlari ayrica listeler. O satirlarda karar minimap2'ye
# birakilmaz; elle hizalama ile hangisinin hakli oldugu belirlenir. Bu betik
# hangisinin hakli oldugunu KENDI BASINA soylemez, yalnizca ayrildiklari yeri
# gosterir.
# ---------------------------------------------------------------------------

from __future__ import print_function
import argparse
import glob
import io
import os
import sys
import time


def _kd_yukle(kok):
    import importlib.util as u
    y = os.path.join(kok, 'verification', 'kimlik_dogrulama.py')
    if not os.path.exists(y):
        sys.stderr.write(u'HATA: %s yok. --kok proje klasorunu gostermeli.\n' % y)
        return None
    sp = u.spec_from_file_location('kd', y)
    m = u.module_from_spec(sp)
    try:
        sp.loader.exec_module(m)
    except SystemExit:
        pass
    return m


def kutu_konsensuslari(kok, en_fazla):
    """Gercek kutu konsensuslerini okur. Kanonik olanlar tercih edilir."""
    adaylar = []
    kn = os.path.join(kok, 'konsensus_kanonik')
    if os.path.isdir(kn):
        adaylar = sorted(glob.glob(os.path.join(kn, '*.kanonik.fa')))
    if not adaylar:
        adaylar = sorted(glob.glob(os.path.join(kok, 'consensus sequences', '*', '*.fasta')))
    cikti = []
    for y in adaylar[:en_fazla]:
        try:
            satir = io.open(y, encoding='utf-8', errors='ignore').read().split('\n')
        except Exception:
            continue
        diz = ''.join(s.strip() for s in satir if s and not s.startswith('>'))
        if len(diz) >= 300:
            cikti.append((os.path.basename(y), diz.upper()))
    return cikti


def main():
    p = argparse.ArgumentParser(description=u'Iki hizalayiciyi yan yana olcer')
    p.add_argument('--kok', default='.')
    p.add_argument('--kutu', type=int, default=4, help=u'kac kutu denenecek')
    p.add_argument('--kayit', type=int, default=120,
                   help=u'veritabani basina kac kayit hizalanacak')
    p.add_argument('--vtb', default='', help=u'yalniz adinda bu gecen veritabanlari')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)

    sys.path.insert(0, os.path.join(kok, 'verification'))
    try:
        import hizalayici_minimap2 as mm
    except ImportError:
        sys.stderr.write(u'HATA: verification/hizalayici_minimap2.py bulunamadi.\n')
        return 1

    if not mm.var_mi():
        sys.stderr.write(
            u'\nKARSILASTIRMA YAPILAMADI: mappy calismiyor.\n'
            u'  Sebep: %s\n'
            u'  Kurulum: pip install mappy\n'
            u'           (ya da: micromamba install -n mikro -c bioconda minimap2\n'
            u'                   micromamba run -n mikro pip install mappy)\n'
            u'\n'
            u'  Mevcut motor calismaya devam eder; zincir bundan ETKILENMEZ.\n'
            u'  Varsayilan hizalayici degismedi: saf Python.\n\n' % mm.sebep())
        return 2

    kd = _kd_yukle(kok)
    if kd is None:
        return 1

    kutular = kutu_konsensuslari(kok, a.kutu)
    if not kutular:
        sys.stderr.write(u'HATA: kutu konsensusu bulunamadi.\n')
        return 1

    vtb = [(e, d) for e, d, _t, kullan, _n in kd.VTB
           if kullan and os.path.exists(os.path.join(kok, 'REFERANS_DB', d))]
    if a.vtb:
        vtb = [v for v in vtb if a.vtb.lower() in v[1].lower()]
    if not vtb:
        sys.stderr.write(u'HATA: REFERANS_DB altinda kullanilabilir veritabani yok.\n')
        return 1

    satirlar = []
    t_py_top = t_mm_top = 0.0
    print(u'%d kutu x %d veritabani, veritabani basina %d kayit'
          % (len(kutular), len(vtb), a.kayit))

    for kutu, q in kutular:
        for etiket, dosya in vtb:
            yol = os.path.join(kok, 'REFERANS_DB', dosya)
            kayitlar = []
            for bas, diz in kd.fasta_akisi(yol):
                kayitlar.append((bas, diz))
                if len(kayitlar) >= a.kayit:
                    break
            if not kayitlar:
                continue

            # --- motor 1: saf Python ---
            t0 = time.time()
            py = {}
            for bas, diz in kayitlar:
                py[bas] = kd.hizala(q, diz)
            t_py = time.time() - t0

            # --- motor 2: minimap2, toplu indeks ---
            t0 = time.time()
            mmr = mm.toplu_hizala(q, [(bas, diz) for bas, diz in kayitlar])
            t_mm = time.time() - t0

            t_py_top += t_py
            t_mm_top += t_mm

            py_s = sorted(py.items(), key=lambda x: -x[1][0])
            mm_s = sorted(mmr.items(), key=lambda x: -x[1][0])
            py1, mm1 = py_s[0][0], mm_s[0][0]
            ust5_py = set(x[0] for x in py_s[:5])
            ust5_mm = set(x[0] for x in mm_s[:5])
            ortusme = len(ust5_py & ust5_mm)
            sapma = abs(py[py1][0] - mmr.get(py1, (0.0, 0))[0])

            satirlar.append(dict(
                kutu=kutu, vtb=etiket, kayit=len(kayitlar),
                ayni_isabet='EVET' if py1 == mm1 else 'HAYIR',
                py_en_iyi=py1[:52], py_yuzde=py[py1][0],
                mm_en_iyi=mm1[:52], mm_yuzde=mmr[mm1][0],
                ayni_kayitta_sapma=round(sapma, 2),
                ust5_ortusme='%d/5' % ortusme,
                py_sn=round(t_py, 2), mm_sn=round(t_mm, 2),
                hizlanma=round(t_py / t_mm, 1) if t_mm > 0.001 else 'olculemedi'))
            print(u'  %-28s %-22s isabet=%s  sapma=%.2f  ust5=%d/5  %.1fx'
                  % (kutu[:28], etiket[:22], satirlar[-1]['ayni_isabet'],
                     sapma, ortusme,
                     (t_py / t_mm) if t_mm > 0.001 else 0))

    if not satirlar:
        sys.stderr.write(u'HATA: hicbir karsilastirma yapilamadi.\n')
        return 1

    ayni = sum(1 for s in satirlar if s['ayni_isabet'] == 'EVET')
    en_buyuk_sapma = max(s['ayni_kayitta_sapma'] for s in satirlar)
    tam_ortusen = sum(1 for s in satirlar if s['ust5_ortusme'] == '5/5')
    hiz = (t_py_top / t_mm_top) if t_mm_top > 0.001 else 0

    # --- KARAR ---
    # Uc sart birden saglanmadikca minimap2 varsayilan YAPILMAZ.
    sartlar = [
        (u'butun satirlarda ayni en iyi isabet', ayni == len(satirlar)),
        (u'kimlik sapmasi her satirda 0,5 puanin altinda', en_buyuk_sapma < 0.5),
        (u'ilk bes isabet her satirda birebir ortusuyor', tam_ortusen == len(satirlar)),
    ]
    gecti = all(x[1] for x in sartlar)

    tsv = os.path.join(kok, 'MINIMAP2_KARSILASTIRMA.tsv')
    bas = ['kutu', 'vtb', 'kayit', 'ayni_isabet', 'py_en_iyi', 'py_yuzde',
           'mm_en_iyi', 'mm_yuzde', 'ayni_kayitta_sapma', 'ust5_ortusme',
           'py_sn', 'mm_sn', 'hizlanma']
    with io.open(tsv, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'\t'.join(bas) + u'\n')
        for s in satirlar:
            fh.write(u'\t'.join(unicode(s[k]) if sys.version_info[0] < 3
                                else str(s[k]) for k in bas) + u'\n')

    md = os.path.join(kok, 'MINIMAP2_KARSILASTIRMA.md')
    L = []
    A = L.append
    A(u'# minimap2 ile saf Python hizalayıcının karşılaştırması\n')
    A(u'Üretim: %s · mappy %s\n' % (time.strftime('%Y-%m-%d %H:%M'), mm.surum()))
    A(u'Ölçüm: %d kutu × %d veritabanı = %d karşılaştırma, veritabanı başına '
      u'%d kayıt.\n' % (len(kutular), len(vtb), len(satirlar), a.kayit))
    A(u'## Karar\n')
    A(u'**%s**\n' % (u'minimap2 VARSAYILAN YAPILABİLİR.' if gecti
                     else u'MEVCUT MOTORDA KALINIR. minimap2 varsayılan yapılmadı.'))
    A(u'| Şart | Sonuç |')
    A(u'|---|---|')
    for ad, ok in sartlar:
        A(u'| %s | %s |' % (ad, u'geçti' if ok else u'**KALDI**'))
    A(u'')
    A(u'Hızlanma: **%.1f kat** (saf Python %.1f sn, minimap2 %.1f sn).\n'
      % (hiz, t_py_top, t_mm_top))
    if not gecti:
        A(u'> Hızlı olan doğru sayılmaz. Yukarıdaki şartlardan biri bile '
          u'kalırsa, ayrılan satırlarda hangisinin haklı olduğu elle hizalamayla '
          u'belirlenmeden minimap2 varsayılan yapılmaz.\n')
    ayrilan = [s for s in satirlar if s['ayni_isabet'] == 'HAYIR'
               or s['ust5_ortusme'] != '5/5' or s['ayni_kayitta_sapma'] >= 0.5]
    if ayrilan:
        A(u'## Ayrılan satırlar — elle doğrulama gerekir\n')
        A(u'| Kutu | Veritabanı | Python en iyi | % | minimap2 en iyi | % | İlk 5 |')
        A(u'|---|---|---|---|---|---|---|')
        for s in ayrilan:
            A(u'| %s | %s | %s | %s | %s | %s | %s |'
              % (s['kutu'], s['vtb'], s['py_en_iyi'], s['py_yuzde'],
                 s['mm_en_iyi'], s['mm_yuzde'], s['ust5_ortusme']))
        A(u'')
    A(u'## Bütün ölçümler\n')
    A(u'Ham tablo: `MINIMAP2_KARSILASTIRMA.tsv`\n')
    A(u'| Kutu | Veritabanı | Aynı isabet | Sapma | İlk 5 | Python sn | minimap2 sn | Hızlanma |')
    A(u'|---|---|---|---|---|---|---|---|')
    for s in satirlar:
        A(u'| %s | %s | %s | %s | %s | %s | %s | %s |'
          % (s['kutu'], s['vtb'], s['ayni_isabet'], s['ayni_kayitta_sapma'],
             s['ust5_ortusme'], s['py_sn'], s['mm_sn'], s['hizlanma']))
    A(u'\n## Nerede kullanılmaz\n')
    A(u'Bu karşılaştırma yalnız **kimlik aşamasının veritabanı taramasını** '
      u'kapsar. minimap2 primer bağlanma aramasında ve in-siliko PCR ürün '
      u'hesabında **kullanılmaz**: primerler 18-25 bazdır, minimap2 bu boydaki '
      u'sorgular için tasarlanmadı ve tohum bulamadığında bağlanma yerini '
      u'sessizce kaçırır. Oralarda güvercin yuvası motoru kalır, çünkü onun '
      u'kayıpsızlığı bir garantidir.\n')
    io.open(md, 'w', encoding='utf-8').write(u'\n'.join(L) + u'\n')

    print(u'\nyazildi: %s' % md)
    print(u'         %s' % tsv)
    print(u'KARAR  : %s' % (u'minimap2 varsayilan YAPILABILIR'
                            if gecti else u'MEVCUT MOTORDA KALINIR'))
    print(u'HIZ    : %.1f kat' % hiz)
    return 0


if __name__ == '__main__':
    sys.exit(main() or 0)
