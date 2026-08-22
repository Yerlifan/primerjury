#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_primer_geometry.py
Audits from scratch whether the forward and reverse primers were designed
correctly, using NOT ONE function of the design code.

Why a separate script: if the design code confirms itself with the coordinates it
produced itself, a coordinate fault runs the same way in both and stays
invisible. This script takes only the template sequence, the primer sequences and
the reported positions; it cuts the product out of the template and tests these
four conditions independently:

  1. The START of the product is the forward primer itself.
  2. The END of the product is the REVERSE COMPLEMENT of the reverse primer.
  3. Read on the minus strand of the template, the 3' end of the reverse primer
     faces into the product, that is, it falls on the 5' end of its counterpart
     on the plus strand.
  4. The reported product length equals the length of the product cut out.

It also checks whether the 3' end of each primer falls on an ambiguous position
in the template, an IUPAC code or an N. That is easy to miss for the reverse
primer, because its 3' end falls at the START of the window.

Usage:
  python3 check_primer_geometry.py --tsv primer_candidates/X__A1.tsv \
      --consensus reference_consensus/dominant/consensus --anchor A1-1_2209
"""
import argparse, csv, glob, os, re, sys

TAM = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn",
                    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T",
         "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def uyar(oligo, kalip_parca):
    """Is the oligo compatible with this piece of the template.

        The template can hold an IUPAC code while the primer is pure ACGT and has
        been resolved to one of the alleles that code stands for. So what is looked
        for is not character equality but SET COMPATIBILITY. An N in the template
        means that position is unknown, and it does not count as compatible."""
    if len(oligo) != len(kalip_parca):
        return False
    for o, k in zip(oligo, kalip_parca):
        if k == "N":
            return False
        if o not in IUPAC.get(k, ""):
            return False
    return True


def rc(s):
    """The reverse complement: complemented first, then reversed."""
    return s.translate(TAM)[::-1]


def oku(f):
    return "".join(l.strip() for l in open(f, encoding="utf-8",
                                           errors="replace")
                   if not l.startswith(">")).upper()


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--tsv", required=True)
    p.add_argument("--consensus", required=True)
    p.add_argument("--anchor", default=None,
                   help="anchor label; if omitted TSV'deki capa column is used")
    p.add_argument("--max", type=int, default=200)
    return p.parse_args()


def main():
    a = get_args()
    rows = list(csv.DictReader(open(a.tsv, encoding="utf-8"), delimiter="\t"))
    if not rows:
        sys.exit("bos TSV")
    capa_ad = a.anchor or rows[0].get("capa", "")
    aday = sorted(glob.glob(os.path.join(a.consensus, "*%s*" % capa_ad))) or \
        sorted(glob.glob(os.path.join(a.consensus, "*")))
    if not aday:
        sys.exit(u'the anchor consensus was not found: %s' % capa_ad)
    kalip = oku(aday[0])
    print('the anchor file : %s' % os.path.basename(aday[0]))
    print("kalip uzunluk: %d" % len(kalip))
    print(u'rows tested: %d\n' % min(len(rows), a.max))

    say = dict(tamam=0, urun_basi=0, urun_sonu=0, yon=0, uzunluk=0,
               f_3p_belirsiz=0, r_3p_belirsiz=0, kalip_disi=0)
    ornek = []
    for r in rows[:a.max]:
        F, R = r["ileri_dizi"], r["geri_dizi"]
        try:
            fb = int(r["ileri_baslangic"]) - 1
            fl = int(r["ileri_uzunluk"])
            gb = int(r["geri_baslangic"]) - 1
            gl = int(r["geri_uzunluk"])
            bildirilen_min = int(r["urun_min"])
            bildirilen_maks = int(r.get("urun_maks") or r["urun_min"])
        except (KeyError, ValueError):
            sys.exit(u'the TSV has no position columns; the current version of design_group_primers.py is needed')
        if fb < 0 or gb + gl > len(kalip):
            say["kalip_disi"] += 1
            continue

        # Urun: ileri primerin 5' ucundan, geri primerin kalip uzerindeki
        # penceresinin sonuna kadar. Geri primer eksi zincirde okundugu icin
        # kalip uzerindeki penceresi [gb, gb+gl) araligidir.
        urun = kalip[fb:gb + gl]
        hata = []

        # 1. urunun basi = ileri primer
        if not uyar(F, urun[:fl]):
            say["urun_basi"] += 1
            hata.append("urun_basi")
        # 2. urunun sonu = geri primerin ters tumleyeni
        if not uyar(rc(R), urun[-gl:]):
            say["urun_sonu"] += 1
            hata.append("urun_sonu")
        # 3. yon: geri primer kalip penceresinin ters tumleyeni olmali.
        #    Kalip kodlarinin ters tumleyeni alinip oligo ona karsi sinanir.
        if not uyar(R, rc(kalip[gb:gb + gl])):
            say["yon"] += 1
            hata.append("yon")
        # 4. uzunluk: bildirilen aralik uyelere gore hesaplanir, capadaki
        #    urun bu araligin icinde olmali
        if not (bildirilen_min <= len(urun) <= bildirilen_maks):
            say["uzunluk"] += 1
            hata.append("uzunluk(%d, bildirilen %d-%d)"
                        % (len(urun), bildirilen_min, bildirilen_maks))
        # 5. 3' uclarin kalipta belirsiz pozisyona denk gelmesi
        #    ileri primerin 3' ucu  -> kalipta fb+fl-1
        #    geri primerin 3' ucu   -> kalipta gb (pencerenin BASI)
        if kalip[fb + fl - 1] not in "ACGT":
            say["f_3p_belirsiz"] += 1
            hata.append("F_3p_belirsiz(%s)" % kalip[fb + fl - 1])
        if kalip[gb] not in "ACGT":
            say["r_3p_belirsiz"] += 1
            hata.append("R_3p_belirsiz(%s)" % kalip[gb])

        if not hata:
            say["tamam"] += 1
        elif len(ornek) < 5:
            ornek.append((r, urun, hata))

    n = min(len(rows), a.max)
    print(u'RESULT')
    print(u'   rows passing all four geometry conditions : %d / %d' % (say["tamam"], n))
    print(u'   the product does not start with the forward primer : %d' % say["urun_basi"])
    print(u'   the product does not end with the rc of the reverse primer: %d' % say["urun_sonu"])
    print(u'   the reverse primer is not the rc of the template   : %d' % say["yon"])
    print(u'   the reported length does not hold                : %d'
          % say["uzunluk"])
    print(u'   the forward 3-prime end is an ambiguous base     : %d'
          % say["f_3p_belirsiz"])
    print(u'   the reverse 3-prime end is an ambiguous base     : %d'
          % say["r_3p_belirsiz"])
    print(u'   rows running off the end of the template          : %d' % say["kalip_disi"])
    for r, urun, hata in ornek:
        print(u'\n   ERROR %s' % ", ".join(hata))
        print("      F=%s  R=%s" % (r["ileri_dizi"], r["geri_dizi"]))
        print(u'      product start=%s ... end=%s' % (urun[:28], urun[-28:]))
        print("      rc(R)   =%s" % rc(r["geri_dizi"]))
    if say["tamam"] == n:
        print(u'\nEvery row passed all four geometry conditions.')


if __name__ == "__main__":
    main()
