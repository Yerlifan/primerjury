# -*- coding: utf-8 -*-
"""cross_coverage.py runs the universal and broad pairs in every bin WITH NO CLASS
BOUNDARY.

The panel's measurements were class based (the A pairs only in the A bins, the F
pairs only in the F bins). Because of that, the question "which taxon falls under no
target" produced artificial gaps that came from the class coverage rather than from
the data. This script measures the five broad targets in ALL 99 bins.

Usage: python cross_coverage.py --fastq "../fastq files" --out capraz.json

"""
# -------------------------------------------------------------------------
# cross_coverage.py measures the five broad or universal pairs in every bin,
#                    without regard to the amplicon class boundary.
#
# INPUT  : the "fastq files" directory given with --fastq (optionally a subset with
#          --dir-path and a per bin read ceiling with --nmax). The five broad pairs are
#          a fixed list inside the file. read_engine.py does the measuring.
# OUTPUT : the json file given with --out; each row carries the product count per
#          bin under the key "<panel_row>|<bin>". If the file exists it is appended
#          to, so an interrupted run continues where it stopped.
# CALLED BY: IT IS NOT IN THE MENU, it is run by hand. Its output becomes the
#          --cross input of target_taxon_mapping.py.
#
# WHY THE CLASS BOUNDARY IS LIFTED: the panel's measurements were class based (the A
# pairs only in the A bins). That produced artificial gaps in the question "which
# taxon falls under no target": the taxon might really be covered while the
# measurement had never tried it. Here the five broad targets are run over every
# bin, so the gaps come out of a real measurement.
# -------------------------------------------------------------------------
import sys, os, glob, json, argparse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import read_engine as om

GENIS = [(2,  'Metanojen_universal',   'GTGGAGCTTGCGGTTTAATTG',   'CAGGATGCTTCACAGTACGAAC'),
         (6,  'Mantar_universal (F1)', 'GGTTACCCGCTGAACTTAAGC',   'CGCTTCACTCGCCGTTAC'),
         (13, 'Mantar_universal (F2)', 'GTGCATGGCCGTTCTTAGTTG',   'CAAACTTCCATCGGCTTGAGC'),
         (17, 'Bakteri_universal',     'ACAAGCGGTGGAGCATGTG',     'ACGACAGCCATGCAGCAC'),
         (20, 'Arke_universal',        'CTGCGGTTTAATTGGATTCAACGC','GAACTGACGACGGCCATGC')]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fastq', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--dir-path', dest='dizin', default='')
    ap.add_argument('--nmax', type=int, default=3000)
    a = ap.parse_args()
    out = json.load(open(a.out)) if os.path.exists(a.out) else {}
    dizinler = a.dizin.split(',') if a.dizin else sorted(
        os.path.basename(d) for d in glob.glob(os.path.join(a.fastq, '*')) if os.path.isdir(d))
    for d in dizinler:
        t0 = time.time()
        for p in sorted(glob.glob(os.path.join(a.fastq, d, '*.fastq'))):
            rs, n0 = om.kutu_yukle(p, a.nmax, 3)
            tax = os.path.basename(p).split('reads_')[1].split('.fastq')[0]
            for s, ad, F, R in GENIS:
                k = '%d|%s_%s' % (s, d, tax)
                if k in out:
                    continue
                y1, n, _ = om.kutu_pcr(rs, F, R, max_mm=1)
                y3, _, _ = om.kutu_pcr(rs, F, R, max_mm=3)
                out[k] = [y1, y3, n, n0]
        json.dump(out, open(a.out, 'w'))
        print(d, round(time.time() - t0, 1), 's', flush=True)


if __name__ == '__main__':
    main()
