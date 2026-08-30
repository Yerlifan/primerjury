#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==== MEETING NOTE ====
# WHAT IT IS FOR   : Brings the ambiguous IUPAC letters in a consensus down to the dominant allele, using the raw reads.
# INPUT            : the consensus directory and the fastq directory
# OUTPUT           : the dominant allele consensus and the allele shares per position
# HOW TO RUN IT    : python3 dominant_allele_consensus.py --consensus <consensus> --fastq "fastq files" --out <out>
# WHY IT IS LIKE THIS : Because the binding rule evaluates by set intersection, IUPAC ambiguity was making a primer look as though it matched both alleles. MEASURED: two consensuses are 79.8 per cent identical character by character but 99.9 per cent by the intersection criterion.
# =======================
"""
dominant_allele_consensus.py
For each consensus it produces the DOMINANT ALLELE sequence from that bin's own raw
reads, plus the allele ratios per position.

Why it is needed: the `samtools consensus -A` output writes variable positions with
an IUPAC code. Because the binding rule judges IUPAC by set intersection, the
uncertainty the code carries makes the primer look as though it suits both alleles.
Measured: in class F2, the Trichoderma asperellum and Metarhizium brunneum
consensuses are 79.8 percent identical character by character but 99.9 percent under
the intersection criterion. So there is a real difference and the uncertainty hides
it. This script produces a second, unambiguous measurement from the same data; the
design and the specificity can be run over the two sets separately and the results
compared.

The output:
  <out>/konsensus/<label>_baskin_konsensus.fasta
  <out>/oran/<label>_alel.tsv       position, depth, A, C, G, T, dominant, ratio,
                                    wilson_alt
  <out>/ozet.tsv

Usage:
  python3 dominant_allele_consensus.py       --consensus referans_konsensus/self/konsensus       --fastq "fastq files" --out referans_konsensus/baskin

"""
import argparse, csv, glob, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import alignment

if alignment.ARKA_UC is None:
    sys.exit(alignment.durum())

BAZLAR = "ACGT"


def oku_fasta(f):
    return "".join(l.strip() for l in open(f, encoding="utf-8", errors="replace")
                   if not l.startswith(">")).upper()


def wilson_alt(k, n, z=1.96):
    """The Wilson lower bound of the ratio. It keeps single read noise from looking like
    a high ratio; that is the criterion the meeting decision asks for on competitor
    ratios.

    """
    if n == 0:
        return 0.0
    p = k / n
    d = 1 + z * z / n
    m = p + z * z / (2 * n)
    s = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return max(0.0, (m - s) / d)


def fastq_bul(kok, grp, taxid):
    """An uncompressed fastq is looked for. There is no .gz support; if such a file
    exists it is warned about plainly rather than skipped silently.

    """
    d = os.path.join(kok, grp)
    for p in sorted(glob.glob(os.path.join(d, "*reads[-_]%s.fastq" % taxid))):
        return p
    gz = sorted(glob.glob(os.path.join(d, "*reads[-_]%s.fastq.gz" % taxid)))
    if gz:
        print(u'   WARNING: only a compressed file is present and that is not supported: %s'
              % os.path.basename(gz[0]))
    return None


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--consensus", required=True, help="self consensus directory")
    p.add_argument("--fastq", required=True, help="'fastq files' directory")
    p.add_argument("--out", required=True)
    p.add_argument("--max-reads", type=int, default=3000)
    p.add_argument("--min-depth", type=int, default=20,
                   help="below this depth no base is called and N is written")
    p.add_argument("--min-fraction", type=float, default=0.50,
                   help="write N when the dominant allele is below this fraction")
    p.add_argument("--min-length", type=int, default=400)
    p.add_argument("--write-fractions", type=int, default=1)
    p.add_argument("--batch", type=int, default=2000,
                   help='how many reads to pass in a single call to the '
                        'minimap2 command line backend')
    return p.parse_args()


def main():
    a = get_args()
    os.makedirs(os.path.join(a.out, "konsensus"), exist_ok=True)
    if a.write_fractions:
        os.makedirs(os.path.join(a.out, "oran"), exist_ok=True)
    ozet = []
    dosyalar = sorted(glob.glob(os.path.join(a.consensus, "*_konsensus.fasta")))
    if not dosyalar:
        sys.exit(u'no consensus found: %s' % a.consensus)
    print(alignment.durum())
    print(u'consensus files: %d' % len(dosyalar))
    for f in dosyalar:
        etiket = re.sub(r"_(ref|self)_konsensus\.fasta$", "", os.path.basename(f))
        m = re.match(r"((?:A1|A2|B|F1|F2)-\d+)_(\d+)$", etiket)
        if not m:
            print(u'   SKIPPED, the label could not be resolved: %s' % etiket)
            continue
        grp, taxid = m.group(1), m.group(2)
        ref = oku_fasta(f)
        # The N's ARE NOT DELETED. Deleting them sticks the two sides of an
        # inner coverage gap together and forms a junction that does not exist
        # in nature; because batch_design.py uses this sequence as a design
        # template, a primer can be chosen across that junction and it will
        # never work. The N's at the start and the end are trimmed (the
        # coordinate shift is reported in the output header) and the inner N's
        # are left in place; minimap2 counts an N as a mismatch anyway.
        bas = 0
        while bas < len(ref) and ref[bas] == "N":
            bas += 1
        son = len(ref)
        while son > bas and ref[son - 1] == "N":
            son -= 1
        cekirdek = ref[bas:son]
        ic_n = cekirdek.count("N")
        if len(cekirdek) - ic_n < 200:
            print(u'   SKIPPED, the consensus is too short: %s' % etiket)
            continue
        fq = fastq_bul(a.fastq, grp, taxid)
        if not fq:
            print(u'   SKIPPED, no fastq: %s' % etiket)
            continue
        try:
            A = alignment.Hizalayici(seq=cekirdek, preset="map-ont")
        except RuntimeError as e:
            sys.exit(str(e))
        if not A:
            print(u'   SKIPPED, the index could not be built: %s' % etiket)
            continue
        say = [[0, 0, 0, 0] for _ in cekirdek]
        okumalar = {}
        n = 0
        with open(fq, errors="replace") as fh:
            for i, line in enumerate(fh):
                if i % 4 != 1:
                    continue
                r = line.strip().upper()   # a lower case FASTQ is accepted too
                if len(r) < a.min_length:
                    continue
                n += 1
                if n > a.max_reads:
                    break
                okumalar["o%d" % n] = r
        hiz = 0
        # A bulk alignment: starting a process per read on the minimap2 command
        # line backend would be unacceptably slow.
        for bas in range(0, len(okumalar), a.batch):
            parca = dict(list(okumalar.items())[bas:bas + a.batch])
            for adq, hl in A.map_toplu(parca):
                if not hl:
                    continue
                h = max(hl, key=lambda x: x.blen)
                hiz += 1
                r = parca[adq]
                # Ters zincirde q_st/q_en OZGUN sorgu koordinatindadir,
                # CIGAR ters tumleyen uzerinde yurur.
                if h.strand == 1:
                    q, qp = r, h.q_st
                else:
                    q, qp = alignment.revcomp(r), len(r) - h.q_en
                rp = h.r_st
                for ln, op in h.cigar:
                    if op == 0:
                        for k in range(ln):
                            if rp + k < len(cekirdek) and qp + k < len(q):
                                b = q[qp + k]
                                j = BAZLAR.find(b)
                                if j >= 0:
                                    say[rp + k][j] += 1
                        rp += ln; qp += ln
                    elif op == 1:
                        qp += ln
                    elif op == 2:
                        rp += ln
                    elif op == 4:
                        qp += ln
        cikti = []
        belirsiz = dusuk = 0
        oran_satir = []
        for i, c in enumerate(say):
            d = sum(c)
            if d == 0 or d < a.min_depth:
                cikti.append("N"); dusuk += 1
                if a.write_fractions:
                    oran_satir.append([i + 1, d, c[0], c[1], c[2], c[3], "N", 0.0, 0.0])
                continue
            en = max(c)
            esitler = [x for x in range(4) if c[x] == en]
            if len(esitler) > 1:
                # On an exact tie the alphabetical order must not decide; if two alleles
                # are at an equal ratio there is no dominant allele and an N is written.
                cikti.append("N"); belirsiz += 1
                if a.write_fractions:
                    oran_satir.append([i + 1, d, c[0], c[1], c[2], c[3], "N",
                                       round(en / d, 4), 0.0])
                continue
            j = esitler[0]
            oran = c[j] / d
            w = wilson_alt(c[j], d)
            if oran < a.min_fraction:
                cikti.append("N"); belirsiz += 1
                baz = "N"
            else:
                baz = BAZLAR[j]
                cikti.append(baz)
            if a.write_fractions:
                oran_satir.append([i + 1, d, c[0], c[1], c[2], c[3], baz,
                                   round(oran, 4), round(w, 4)])
        dizi = "".join(cikti)
        yol = os.path.join(a.out, "konsensus", "%s_baskin_konsensus.fasta" % etiket)
        with open(yol, "w", encoding="utf-8") as fh:
            fh.write(u'>%s dominant_allele reads=%d aligned=%d min_depth=%d min_fraction=%.2f trim=%d-%d inner_N_in_input=%d\n'
                     % (etiket, n, hiz, a.min_depth, a.min_fraction,
                        bas + 1, son, ic_n))
            for k in range(0, len(dizi), 70):
                fh.write(dizi[k:k + 70] + "\n")
        if a.write_fractions:
            with open(os.path.join(a.out, "oran", "%s_alel.tsv" % etiket), "w",
                      newline="", encoding="utf-8") as fh:
                w = csv.writer(fh, delimiter="\t", lineterminator="\n")
                w.writerow(["poz", "derinlik", "A", "C", "G", "T", "baskin",
                            "oran", "wilson_alt"])
                w.writerows(oran_satir)
        kaps = len(dizi) - dizi.count("N")
        ozet.append(dict(etiket=etiket, grup=grp, taxid=taxid, okuma=n,
                         hizalanan=hiz, uzunluk=len(dizi), kapsanan=kaps,
                         dusuk_derinlik=dusuk, belirsiz_cogunluk=belirsiz,
                         kirpma_bas=bas + 1, kirpma_son=son,
                         girdideki_ic_N=ic_n))
        print(u'   %-26s reads=%5d aligned=%5d length=%5d covered=%5d low_depth=%4d ambiguous=%4d'
              % (etiket, n, hiz, len(dizi), kaps, dusuk, belirsiz))
    if ozet:
        with open(os.path.join(a.out, "ozet.tsv"), "w", newline="",
                  encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(ozet[0].keys()),
                               delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(ozet)
        print("\nwritten: %s" % os.path.join(a.out, "ozet.tsv"))


if __name__ == "__main__":
    main()
