# -*- coding: utf-8 -*-
"""
orientation_code_scan.py finds EVERY script that reads a consensus and classifies
how it handles orientation. The aim: to answer "which code route is still unpatched"
by scanning rather than by guessing.

The classes:
  NORMALIZE_KAYNAK   it reads the normalised directory (referans_konsensus/...)
  IKI_YON_DENIYOR    it tries the sequence both as it is and as rc (orientation does
                     not matter)
  KENDI_YAMASI       it corrects the orientation inside itself (alignment, edlib,
                     picking the better direction)
  HAM_KAYNAK_VARSAYIM  it reads a mixed orientation directory and assumes ONE
                     direction  <-- RISKY
  YON_ILGISIZ        it reads a consensus but does no orientation dependent work
                     (counting, reporting, naming)
Usage: python orientation_code_scan.py --root .. --out ../yon_kod_taramasi.tsv

"""
# -------------------------------------------------------------------------
# orientation_code_scan.py reads every .py file in the project, finds the code
#                       routes that read a consensus and classifies them by how
#                       they handle orientation.
#
# INPUT  : every Python source file under --root; the classification patterns are
#          fixed inside the file. The source text is parsed with ast and the comment
#          and docstring lines are DROPPED.
# OUTPUT : the TSV given with --out: the file, the class (NORMALIZE_KAYNAK /
#          IKI_YON_DENIYOR / KENDI_YAMASI / HAM_KAYNAK_VARSAYIM / YON_ILGISIZ) and
#          the reason.
# CALLED BY: IT IS NOT IN THE MENU, it is an audit run by hand. Its output is
#          written into the panel's orientation sheet by orientation_report.py.
#
# WHY COMMENTS AND DOCSTRINGS ARE DROPPED: a directory name appearing in an
# explanation line is not a running code route. That produced false positives on the
# first scan: targets.py and config.py mention "consensus sequences" only in an
# explanation, and the scan was marking them RISKLI. The aim is to answer "which
# code route is still unpatched" by measurement rather than by guessing.
# -------------------------------------------------------------------------
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
    """Drop the comments and docstrings. A directory name appearing in an explanation
    line is not a code route; that produced false positives on the first scan (targets.py
    and config.py mention 'consensus sequences' only in an explanation).

    """
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
    ap.add_argument('--root', dest='kok', required=True)
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
    print(u'\ntotal', len(satirlar), '->', a.out)


if __name__ == '__main__':
    main()
