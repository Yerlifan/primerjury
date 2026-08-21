# -*- coding: utf-8 -*-
"""cross_coverage.py - evrensel/genis ciftleri SINIF SINIRI OLMADAN her kutuya kosar.

Panelin olcumleri sinif bazliydi (A ciftleri yalniz A kutularinda, F ciftleri yalniz
F kutularinda). Bu yuzden "hangi takson hicbir hedefin kapsamina girmiyor" sorusu
sinif kapsamindan dogan yapay bosluklar uretiyordu. Bu betik bes genis hedefi
99 kutunun HEPSINDE olcer.

Kullanim: python cross_coverage.py --fastq "..\fastq files" --out capraz.json
"""
# ---------------------------------------------------------------------------
# cross_coverage.py — bes genis/evrensel cifti, amplikon sinifi sinirini
#                    gozetmeden butun kutularda olcer.
#
# GIRDI  : --fastq ile "fastq files" klasoru (istege bagli --dizin ile alt
#          kume, --nmax ile kutu basina okuma tavani). Bes genis cift dosyanin
#          icinde sabit listedir. Olcumu read_engine.py yapar.
# CIKTI  : --out ile verilen json dosyasi; her satir "<panel_satiri>|<kutu>"
#          anahtariyla kutu basina urun sayisini tasir. Dosya varsa uzerine
#          eklenir, boylece kosu kesilse de kaldigi yerden devam eder.
# CAGRAN : MENUDE DEGILDIR - elle calistirilir. Ciktisi target_taxon_mapping.py
#          icin --capraz girdisi olur.
#
# NEDEN SINIF SINIRI KALDIRILIYOR: panelin olcumleri sinif bazliydi (A ciftleri
# yalniz A kutularinda). Bu, "hangi takson hicbir hedefin kapsamina girmiyor"
# sorusuna yapay bosluklar uretiyordu - takson aslinda kapsaniyor olabilir ama
# olcum onu hic denememis oluyordu. Burada bes genis hedef butun kutulara
# kosulur, boylece bosluklar gercek olcumden cikar.
# ---------------------------------------------------------------------------
import sys, os, glob, json, argparse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import okuma_motoru as om

GENIS = [(2,  'Metanojen_universal',   'GTGGAGCTTGCGGTTTAATTG',   'CAGGATGCTTCACAGTACGAAC'),
         (6,  'Mantar_universal (F1)', 'GGTTACCCGCTGAACTTAAGC',   'CGCTTCACTCGCCGTTAC'),
         (13, 'Mantar_universal (F2)', 'GTGCATGGCCGTTCTTAGTTG',   'CAAACTTCCATCGGCTTGAGC'),
         (17, 'Bakteri_universal',     'ACAAGCGGTGGAGCATGTG',     'ACGACAGCCATGCAGCAC'),
         (20, 'Arke_universal',        'CTGCGGTTTAATTGGATTCAACGC','GAACTGACGACGGCCATGC')]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fastq', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--dizin', default='')
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
