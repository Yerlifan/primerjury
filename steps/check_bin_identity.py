#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_bin_identity.py
Kraken2'nin ayirdigi her takson kutusunun okumalarini, referans
veritabanindaki GERCEK dizilere hizalayarak kutunun kimligini bagimsiz
olarak sinar.

Neden gerekli: Kraken2 varsayilan ayarda cekimser kalmaz. Numunedeki
organizma veritabaninda yoksa okumalar kardes turlere dagitilir ve tek bir
populasyon birden cok "tur" kutusuna bolunur. Bu durumda kutularin
konsensusleri birbirinin ayni cikar ve o hedefte ozgullugu saglamak
mantiken imkansiz hale gelir.

Olcum: her kutudan en cok --ornek okuma alinir, minimap2 (mappy) ile
verilen referans FASTA'daki adaylara hizalanir, okuma basina EN IYI eslesen
referans oylanir. Kutunun kimligi bu oylamanin cogunlugudur.

Kullanim:
  python3 check_bin_identity.py \
      --fastq "/.../fastq files" --ref REFERANS_DB/archaea.16S.fna \
      --desen "Methanosarcina" --out kutu_kimlik.tsv
"""
import argparse, csv, glob, os, sys

try:
    import mappy
except ImportError:
    sys.exit("mappy gerekli: pip install mappy")


def fasta(p):
    ad, buf = None, []
    for l in open(p, errors="replace"):
        if l.startswith(">"):
            if ad:
                yield ad, "".join(buf)
            ad, buf = l[1:].strip(), []
        else:
            buf.append(l.strip())
    if ad:
        yield ad, "".join(buf)


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--fastq", required=True, help="'fastq files' klasoru")
    p.add_argument("--ref", required=True, help="referans FASTA (16S/ITS/28S)")
    p.add_argument("--desen", nargs="+", required=True,
                   help="referans basliginda aranacak ad parcalari")
    p.add_argument("--grup", default=None, help="yalniz bu grup (or. A2-4)")
    p.add_argument("--ornek", type=int, default=1000, help="kutu basina okuma")
    p.add_argument("--min-uzunluk", type=int, default=600)
    p.add_argument("--out", default=None)
    return p.parse_args()


def main():
    a = get_args()
    ref = {}
    for ad, s in fasta(a.ref):
        for d in a.desen:
            if d in ad and ad not in ref:
                ref[ad] = s.upper().replace("U", "T")
    if not ref:
        sys.exit("desene uyan referans bulunamadi: %s" % ", ".join(a.desen))
    print("referans kaydi: %d" % len(ref))
    for ad in list(ref)[:10]:
        print("   %-70s %d bp" % (ad[:70], len(ref[ad])))
    alig = {ad: mappy.Aligner(seq=s, preset="map-ont") for ad, s in ref.items()}

    satir = []
    desen = os.path.join(a.fastq, a.grup if a.grup else "*", "*.fastq")
    for fq in sorted(glob.glob(desen)):
        base = os.path.basename(fq)
        grp = os.path.basename(os.path.dirname(fq))
        say, n = {}, 0
        with open(fq, errors="replace") as fh:
            for i, line in enumerate(fh):
                if i % 4 != 1:
                    continue
                r = line.strip()
                if len(r) < a.min_uzunluk:
                    continue
                n += 1
                if n > a.ornek:
                    break
                en, kaz = 0, None
                esit = False
                for ad, A in alig.items():
                    b = 0
                    for h in A.map(r):
                        if h.mlen > b:
                            b = h.mlen
                    if b > en:
                        en, kaz, esit = b, ad, False
                    elif b == en and b > 0:
                        esit = True
                anahtar = "belirsiz" if (en == 0 or esit) else kaz
                say[anahtar] = say.get(anahtar, 0) + 1
        if not n:
            continue
        sirali = sorted(say.items(), key=lambda x: -x[1])
        kaz, kac = sirali[0]
        satir.append(dict(grup=grp, dosya=base, okuma=n,
                          baskin_referans=kaz[:80],
                          baskin_oran=round(100.0 * kac / n, 1),
                          dagilim=";".join("%s=%d" % (k[:40], v) for k, v in sirali[:4])))
        print("%-10s %-34s okuma=%5d  %-52s %%%.1f"
              % (grp, base[:34], n, kaz[:52], 100.0 * kac / n))
    if a.out and satir:
        d = os.path.dirname(os.path.abspath(a.out))
        if d:
            os.makedirs(d, exist_ok=True)
        with open(a.out, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(satir[0].keys()),
                               delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(satir)
        print("\nyazildi: %s" % a.out)


if __name__ == "__main__":
    main()
