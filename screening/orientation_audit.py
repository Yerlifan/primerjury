# -*- coding: utf-8 -*-
"""
orientation_audit.py determines the ORIENTATION EACH consensus file IS STORED IN.

Why it is needed: the orientation fault was found and patched separately in AT
LEAST THREE DIFFERENT PLACES over one night (the strand choice in the source study,
"the consensuses are reversed against SILVA" on the design side, and "39 of 58
consensuses are reversed" in the B and F remeasurement). Three separate patches
means there is no single canonical fix. This script ties "I think it was fixed" to
a measurement.

THE METHOD: two INDEPENDENT criteria, and if they disagree the file counts as
BELIRSIZ (the project rule):

  Criterion 1 (the panel's own universal primers):
      In an SSU consensus stored in the sense (reference, plus) direction, the
      forward primer F is found DIRECTLY and so is rc(R), the complement of the
      reverse primer.
      If it is stored antisense, rc(F) and R are found instead.
  Criterion 2 (universal motifs from the literature, independent of the panel):
      SSU  : 515F  GTGYCAGCMGCCGCGGTAA  and the 806R region in its sense form
             GGATTAGATACCC
      ITS  : ITS1  TCCGTAGGTGAACCTGCGG  and  rc(ITS4)  GCATATCAATAAGCGGAGGA

For each criterion it counts whether the motifs are found in the sequence AS IT IS
or in its REVERSE COMPLEMENT. The mismatch tolerance is <=2 (a nanopore consensus
carries the occasional error) and the engine is lossless.

Usage:
    python orientation_audit.py --root ..  --out ../yon_denetimi_20260802.tsv
    python orientation_audit.py --root .. --dir "consensus sequences"

"""
# -------------------------------------------------------------------------
# orientation_audit.py determines the stored orientation of EVERY file in a
#                   consensus directory under two independent criteria and writes
#                   them into a table.
#
# INPUT  : one or more consensus directories given with --dir under --root (the
#          defaults are defined inside the file). The binding site search is done
#          with okuma_motoru.Sonda at a tolerance of <=2 mismatches.
# OUTPUT : the TSV given with --out: per file, the class, the result of each of the
#          two criteria and the final verdict (SENSE / ANTISENSE / BELIRSIZ).
# CALLED BY: IT IS NOT IN THE MENU, it is an audit run by hand. The table it
#          produces is written into the panel as the "19 Yon Normalizasyonu" sheet
#          by orientation_report.py. The canonical production itself is done in
#          build_canonical.py.
#
# WHY TWO CRITERIA: the decision is not left to a single code route. Criterion 1
# uses the panel's own universal primers and criterion 2 literature motifs that are
# independent of the panel. If the two disagree the file counts as BELIRSIZ; a
# consensus flipped the wrong way does as much silent damage as one never flipped.
# -------------------------------------------------------------------------
import sys, os, glob, csv, argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import read_engine as om

# --- criterion 1: the panel's universal pairs (by class) --------------------
PANEL = {
    'A':  ('Arke_universal',        'CTGCGGTTTAATTGGATTCAACGC', 'GAACTGACGACGGCCATGC'),
    'B':  ('Bakteri_universal',     'ACAAGCGGTGGAGCATGTG',      'ACGACAGCCATGCAGCAC'),
    'F1': ('Mantar_universal (F1)', 'GGTTACCCGCTGAACTTAAGC',    'CGCTTCACTCGCCGTTAC'),
    'F2': ('Mantar_universal (F2)', 'GTGCATGGCCGTTCTTAGTTG',    'CAAACTTCCATCGGCTTGAGC'),
}

# --- criterion 2: literature motifs, expected on the sense (plus) strand -----
MOTIF = {
    'SSU': ['GTGYCAGCMGCCGCGGTAA',      # 515F
            'GGATTAGATACCC',            # 806R bolgesi, sense hali
            'AGTCCCGCAACGAGCGCAACCC'],  # 1100 bolgesi, sense (SSU korunmus)
    'ITS': ['TCCGTAGGTGAACCTGCGG',      # ITS1
            'GCATATCAATAAGCGGAGGA'],    # rc(ITS4)
}


def var_mi(dizi, desen, mm=2):
    """is the pattern in the sequence (lossless engine, <=mm mismatches)"""
    s = om.Sonda(desen, uc5=False, max_mm=mm, son2=False)
    return len(s.bul(dizi)) > 0


def olcut1(dizi, sinif):
    ad, F, R = PANEL[sinif]
    duz = int(var_mi(dizi, F)) + int(var_mi(dizi, om.rc(R)))
    ters = int(var_mi(dizi, om.rc(F))) + int(var_mi(dizi, R))
    if duz > ters:
        return 'SENSE', duz, ters
    if ters > duz:
        return 'ANTISENSE', duz, ters
    return 'BELIRSIZ', duz, ters


def olcut2(dizi, tip):
    duz = sum(1 for m in MOTIF[tip] if var_mi(dizi, m))
    ters = sum(1 for m in MOTIF[tip] if var_mi(om.rc(dizi), m))
    if duz > ters:
        return 'SENSE', duz, ters
    if ters > duz:
        return 'ANTISENSE', duz, ters
    return 'BELIRSIZ', duz, ters


def sinifi(yol):
    b = os.path.basename(yol)
    for s in ('A1-', 'A2-'):
        if s in b or s in yol:
            return 'A'
    if 'F1-' in b or 'F1-' in yol:
        return 'F1'
    if 'F2-' in b or 'F2-' in yol:
        return 'F2'
    if 'B-' in b or 'B-' in yol:
        return 'B'
    return '?'


def oku(yol):
    ad, buf = None, []
    with open(yol, errors='ignore') as fh:
        for l in fh:
            if l.startswith('>'):
                if ad is not None:
                    yield ad, ''.join(buf)
                ad, buf = l[1:].strip(), []
            else:
                buf.append(l.strip())
    if ad is not None:
        yield ad, ''.join(buf)


def tara(kok, klasorler):
    satirlar = []
    for kl in klasorler:
        for yol in sorted(glob.glob(os.path.join(kok, kl, '**', '*.fasta'), recursive=True)):
            sn = sinifi(yol)
            if sn == '?':
                continue
            tip = 'ITS' if sn in ('F1', 'F2') else 'SSU'
            seqs = [om.temizle(s) for _, s in oku(yol)]
            if not seqs:
                satirlar.append(dict(klasor=kl, dosya=os.path.basename(yol), sinif=sn,
                                     uzunluk=0, N_yuzde='', olcut1='DOSYA_BOS', olcut2='DOSYA_BOS',
                                     karar='DOSYA_BOS', ayrinti=''))
                continue
            s = max(seqs, key=len)
            if not s:
                satirlar.append(dict(klasor=kl, dosya=os.path.basename(yol), sinif=sn,
                                     uzunluk=0, N_yuzde='', olcut1='DIZI_BOS', olcut2='DIZI_BOS',
                                     karar='DIZI_BOS', ayrinti=os.path.relpath(yol, kok)))
                continue
            nN = s.count('N')
            o1, d1, t1 = olcut1(s, sn)
            o2, d2, t2 = olcut2(s, tip)
            if o1 == o2 and o1 != 'BELIRSIZ':
                karar = o1
            elif 'BELIRSIZ' in (o1, o2) and o1 != o2:
                karar = (o1 if o2 == 'BELIRSIZ' else o2) + ' (tek olcut)'
            elif o1 != o2:
                karar = 'AYRILIK - MASKELE'
            else:
                karar = 'BELIRLENEMEDI'
            satirlar.append(dict(klasor=kl, dosya=os.path.basename(yol), sinif=sn,
                                 uzunluk=len(s), N_yuzde=round(100.0 * nN / len(s), 1),
                                 olcut1='%s (%d/%d)' % (o1, d1, t1),
                                 olcut2='%s (%d/%d)' % (o2, d2, t2),
                                 karar=karar, ayrinti=os.path.relpath(yol, kok)))
    return satirlar


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', dest='kok', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--dir', dest='klasor', action='append', default=[])
    a = ap.parse_args()
    kl = a.klasor or ['consensus sequences',
                      'referans_konsensus/konsensus',
                      'referans_konsensus/baskin/konsensus',
                      'referans_konsensus/self/konsensus',
                      'SCREENING_RESULT/konsensus_yeni']
    satirlar = tara(a.kok, kl)
    if satirlar:
        with open(a.out, 'w', newline='', encoding='utf-8') as fh:
            w = csv.DictWriter(fh, fieldnames=list(satirlar[0].keys()), delimiter='\t')
            w.writeheader()
            for r in satirlar:
                w.writerow(r)
    ozet = {}
    for r in satirlar:
        ozet.setdefault(r['klasor'], {}).setdefault(r['karar'], 0)
        ozet[r['klasor']][r['karar']] += 1
    print('%-44s %s' % (u'DIRECTORY', u'THE DISTRIBUTION OF VERDICTS'))
    for k in sorted(ozet):
        print('%-44s %s' % (k, ', '.join('%s=%d' % x for x in sorted(ozet[k].items()))))
    print(u'\ntotal files:', len(satirlar), '| TSV:', a.out)


if __name__ == '__main__':
    main()
