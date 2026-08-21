# -*- coding: utf-8 -*-
r"""
orientation_code_scan.py - konsensus okuyan HER betigi bulur ve yonu nasil ele aldigini
siniflandirir. Amac: "hangi kod yolu hala yamalanmamis" sorusunu tahminle degil
taramayla cevaplamak.

Siniflar:
  NORMALIZE_KAYNAK   normalize edilmis klasoru okuyor (referans_konsensus/... )
  IKI_YON_DENIYOR    diziyi hem duz hem rc ile deniyor (yon farketmez)
  KENDI_YAMASI       kendi icinde yon duzeltmesi yapiyor (hiza/edlib/en iyi yon secimi)
  HAM_KAYNAK_VARSAYIM  karisik yonlu klasoru okuyup TEK yon varsayiyor  <-- RISKLI
  YON_ILGISIZ        konsensus okuyor ama yon bagimli is yapmiyor (sayim, rapor, ad)
Kullanim: python orientation_code_scan.py --kok .. --out ..\yon_kod_taramasi.tsv
"""
# ---------------------------------------------------------------------------
# orientation_code_scan.py — projedeki her .py dosyasini okuyup konsensus okuyan kod
#                       yollarini bulur ve yonu nasil ele aldiklarina gore
#                       siniflandirir.
#
# GIRDI  : --kok altindaki butun Python kaynak dosyalari; siniflandirma desenleri
#          dosyanin icinde sabittir. Kaynak metin ast ile ayristirilip yorum ve
#          docstring satirlari ELENIR.
# CIKTI  : --out ile verilen TSV: dosya, sinif (NORMALIZE_KAYNAK /
#          IKI_YON_DENIYOR / KENDI_YAMASI / HAM_KAYNAK_VARSAYIM / YON_ILGISIZ)
#          ve gerekce.
# CAGRAN : MENUDE DEGILDIR - elle calistirilan bir denetimdir. Ciktisi
#          orientation_report.py ile panelin yon sayfasina islenir.
#
# YORUM VE DOCSTRING NEDEN ELENIYOR: bir aciklama satirinda gecen klasor adi
# calisan bir kod yolu degildir. Ilk taramada bu yuzden yanlis pozitif cikti -
# targets.py ve config.py "consensus sequences" adini yalniz aciklamada
# aniyor ama tarama onlari RISKLI olarak isaretliyordu. Amac "hangi kod yolu
# hala yamalanmamis" sorusunu tahminle degil olcumle cevaplamaktir.
# ---------------------------------------------------------------------------
import os, re, ast, csv, glob, argparse

KARISIK = ['consensus sequences', 'KAPSAMLI_ARAMA_SONUC/konsensus_yeni', 'konsensus_yeni']
NORMAL = ['referans_konsensus', 'konsN', 'konsensus_kanonik', 'KANONIK']

IKI_YON = [r'for\s+\w+\s+in\s*\(\s*\w+\s*,\s*(?:ispcr\.)?rc\(', r'rc\(seq\)', r'rc\(s\)\s*\)',
           r'both|iki yon|her iki yon', r'\(s,\s*rc\(s\)\)', r'\(seq,\s*rc\(seq\)\)']
YAMA = [r'edlib', r'ters tumleyen', r'ters_tumleyen', r'en iyi yon', r'yon sec', r'yonu duzelt',
        r'align.*rc|rc.*align', r'min\(.*rc\(', r'best.*orient']
YON_BAGIMLI = [r'find_sites', r'amplify', r'Sonda', r'kutu_pcr', r'baglan', r'primer',
               r'pencere', r'tasarim', r'in_silico', r'ispcr']


def kod_govdesi(metin):
    """Yorum ve docstring'leri AT. Aciklama satirinda gecen klasor adi kod yolu
    degildir - ilk taramada bu yuzden yanlis pozitif cikti (targets.py,
    config.py 'consensus sequences'i yalniz aciklamada aniyor)."""
    try:
        agac = ast.parse(metin)
        ds = set()
        for d in ast.walk(agac):
            if isinstance(d, ast.Constant) and isinstance(d.value, str) \
               and hasattr(d, 'lineno') and hasattr(d, 'end_lineno'):
                if d.end_lineno - d.lineno >= 1:
                    ds.update(range(d.lineno, d.end_lineno + 1))
    except Exception:
        ds = set()
    out = []
    for i, l in enumerate(metin.splitlines(), 1):
        if i in ds:
            continue
        out.append(l.split('#', 1)[0])
    return '\n'.join(out)


def sinifla(metin):
    metin = kod_govdesi(metin)
    kaynak_karisik = any(k in metin for k in KARISIK)
    kaynak_normal = any(k in metin for k in NORMAL)
    iki_yon = any(re.search(p, metin, re.I) for p in IKI_YON)
    yama = any(re.search(p, metin, re.I) for p in YAMA)
    yon_bagimli = any(re.search(p, metin, re.I) for p in YON_BAGIMLI)
    if not yon_bagimli:
        return 'YON_ILGISIZ', kaynak_karisik, kaynak_normal, iki_yon, yama
    if iki_yon:
        return 'IKI_YON_DENIYOR', kaynak_karisik, kaynak_normal, iki_yon, yama
    if yama:
        return 'KENDI_YAMASI', kaynak_karisik, kaynak_normal, iki_yon, yama
    if kaynak_normal and not kaynak_karisik:
        return 'NORMALIZE_KAYNAK', kaynak_karisik, kaynak_normal, iki_yon, yama
    if kaynak_karisik:
        return 'HAM_KAYNAK_VARSAYIM', kaynak_karisik, kaynak_normal, iki_yon, yama
    return 'KAYNAK_BELIRSIZ', kaynak_karisik, kaynak_normal, iki_yon, yama


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--kok', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    satirlar = []
    for yol in sorted(glob.glob(os.path.join(a.kok, '**', '*.py'), recursive=True)):
        rel = os.path.relpath(yol, a.kok)
        if rel.startswith('_') or '/_' in rel or '\\_' in rel:
            continue
        try:
            metin = open(yol, encoding='utf-8', errors='ignore').read()
        except Exception:
            continue
        if not re.search(r'konsensus|consensus', metin, re.I):
            continue
        if not re.search(r'konsensus|consensus', kod_govdesi(metin), re.I):
            continue
        s, kk, kn, iy, ym = sinifla(metin)
        satirlar.append(dict(betik=rel, sinif=s,
                             karisik_kaynak='EVET' if kk else '',
                             normalize_kaynak='EVET' if kn else '',
                             iki_yon='EVET' if iy else '',
                             kendi_yamasi='EVET' if ym else '',
                             satir=len(metin.splitlines())))
    with open(a.out, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(satirlar[0].keys()), delimiter='\t')
        w.writeheader()
        for r in satirlar:
            w.writerow(r)
    say = {}
    for r in satirlar:
        say[r['sinif']] = say.get(r['sinif'], 0) + 1
    for k in sorted(say, key=lambda x: -say[x]):
        print('%-22s %3d' % (k, say[k]))
    print('\nRISKLI (HAM_KAYNAK_VARSAYIM):')
    for r in satirlar:
        if r['sinif'] == 'HAM_KAYNAK_VARSAYIM':
            print('  ', r['betik'])
    print('\ntoplam', len(satirlar), '->', a.out)


if __name__ == '__main__':
    main()
