#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reference_identity.py
MEASURES WHAT THE REFERENCE BASED PRIMERS ACTUALLY AMPLIFY IN THE SAMPLE.

design_from_reference.py designed primers from reference database sequences for
the targets the sample cannot meet, and measured whether the sample SUPPORTS
them. A pair designed for Methanosarcina barkeri, for example, gives a product in
63 per cent of the reads in the sample.

But "it gives a product" and "it amplifies the target" are not the same thing.
The earlier measurement showed that the bins in the sample do not agree with
their Kraken2 labels: the bin labelled M. barkeri goes to Methanosarcina
vacuolata at sequence level. So which organism the support comes from has to be
measured separately.

The method: for every pair the reads in the sample are scanned, THE PRODUCT
SEQUENCE is cut out of the reads that give one, a dominant allele consensus is
built from those products, and that consensus is asked of the reference database
with blastn. That answers "what does this primer amplify in the sample" with
sequence evidence.

Usage:
  python3 reference_identity.py --reference <the reference primer table>
      --pt . --db REFERENCE_DB --out <the identity table>
"""
import argparse, csv, collections, glob, importlib.util, os, re, shutil
import subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def yukle(ad, dosya):
    s = importlib.util.spec_from_file_location(ad, os.path.join(HERE, dosya))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


G = yukle("G", "design_group_primers.py")

SINIF_DB = {"A1": "archaea.16S.fna", "A2": "archaea.16S.fna",
            "B": "bacteria.16S.fna", "F1": "fungi.ITS.fna",
            "F2": "fungi.ITS.fna"}


class Esik:
    "The same rules as the group engine's binding search."
    tail_len = 12
    tail_max_mm = 1
    exact_last = 2
    total_max_mm = 3
    min_overlap = 12
    prod_min = 50
    prod_hard_max = 400


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--reference", required=True)
    p.add_argument("--pt", required=True, help="the root holding the 'fastq "
                                               "files' directory")
    p.add_argument("--db", required=True)
    p.add_argument("--reference-targets", default=os.path.join(HERE, "reference_targets.tsv"))
    p.add_argument("--out", required=True)
    p.add_argument("--max-reads", type=int, default=20000)
    p.add_argument("--min-product", type=int, default=30,
                   help="minimum number of products needed to build a consensus")
    p.add_argument("--threads", type=int, default=4)
    return p.parse_args()


def urun_kes(oligo_f, oligo_r, dizi, a):
    """Returns the product sequence when the read gives one, otherwise None.
        Because the template is double stranded, both orientations are tried."""
    rc = G.rc
    for s in (dizi, rc(dizi)):
        idx = G.build_index(s, Esik.tail_len)
        bf = G.find_bindings(oligo_f, s, idx, Esik.tail_len, Esik)
        if not bf:
            continue
        br = G.find_bindings(rc(oligo_r), s, idx, Esik.tail_len, Esik)
        # rc(geri primer) arti zincirde arandi; onun 3' ucu urunun SAG
        # tarafinda, 5' ucu ise sagda daha ileride. Geri primerin kendisi
        # eksi zincirde baglanir.
        br2 = G.find_bindings(oligo_r, rc(s), G.build_index(rc(s), Esik.tail_len),
                              Esik.tail_len, Esik)
        L = len(s)
        for fend, _ in bf:
            fstart = fend - len(oligo_f) + 1
            if fstart < 0:
                continue
            for rend, _ in br2:
                r_left = L - 1 - rend
                if fend >= r_left:
                    continue
                r_right = r_left + len(oligo_r) - 1
                boy = r_right - fstart + 1
                if a.min_boy <= boy <= a.max_boy and r_right < L:
                    return s[fstart:r_right + 1]
    return None


def baskin_konsensus(urunler):
    """Without aligning products of equal length, it calls the dominant base per
        position over the products of the most frequent length. Because the length
        distribution is narrow, being a single amplicon, that is enough; on a wide
        distribution a warning is printed."""
    if not urunler:
        return None, 0, 0
    uz = collections.Counter(len(u) for u in urunler)
    en_uz, en_say = uz.most_common(1)[0]
    sec = [u for u in urunler if len(u) == en_uz]
    say = [collections.Counter(u[i] for u in sec) for i in range(en_uz)]
    kons = "".join(c.most_common(1)[0][0] for c in say)
    return kons, len(sec), len(uz)


def main():
    a = get_args()
    if not shutil.which("blastn"):
        sys.exit(u'blastn was not found')
    rows = list(csv.DictReader(open(a.reference, encoding="utf-8"), delimiter="\t"))
    if not rows:
        sys.exit(u'there is no reference pair')
    # hedefin amaclanan organizmasi
    amac = {}
    if os.path.exists(a.reference_targets):
        # Dosya # ile baslayan aciklama satirlariyla basliyor; DictReader
        # bunlari baslik sanip anahtarlari yanlis kurar.
        with open(a.reference_targets, encoding="utf-8") as fh:
            satirlar = [l for l in fh if not l.startswith("#") and l.strip()]
        for r in csv.DictReader(satirlar, delimiter="\t"):
            if r.get("ad"):
                amac[r["ad"]] = (r.get("ic", ""), r.get("taxid", ""))
    if not amac:
        print(u'WARNING: reference_targets.tsv could not be read, the intended organism will stay empty: %s' % a.reference_targets)

    fq = collections.defaultdict(list)
    for p in glob.glob(os.path.join(a.pt, "fastq files", "*", "*.fastq")):
        m = re.search(r"((?:A1|A2|B|F1|F2)-\d+)[-_]reads[-_](\d+)\.fastq$",
                      os.path.basename(p))
        if m:
            fq[m.group(1)[:2].rstrip("-")].append((m.group(1), m.group(2), p))
    # sinif anahtarini duzelt (A1-4 -> A1)
    fq2 = collections.defaultdict(list)
    for p in glob.glob(os.path.join(a.pt, "fastq files", "*", "*.fastq")):
        m = re.search(r"((?:A1|A2|B|F1|F2))-(\d+)[-_]reads[-_](\d+)\.fastq$",
                      os.path.basename(p))
        if m:
            fq2[m.group(1)].append(("%s-%s" % (m.group(1), m.group(2)),
                                    m.group(3), p))
    print("fastq envanteri: %d sinif" % len(fq2))

    calisma = tempfile.mkdtemp(prefix="refkim_")
    sonuc = []
    for i, r in enumerate(rows):
        hedef, sinif = r["hedef"], r["sinif"]
        F, R = r["ileri_dizi"], r["geri_dizi"]
        try:
            a.min_boy = max(30, int(r.get("urun_min", 70)) - 40)
            a.max_boy = int(r.get("urun_maks", 300)) + 40
        except (TypeError, ValueError):
            a.min_boy, a.max_boy = 50, 400
        ic, tx = amac.get(hedef, ("", ""))
        # bu hedefin taxid'lerine ait fastq dosyalari
        tset = set(t.strip() for t in tx.split(",") if t.strip())
        dosyalar = [(g, t, p) for g, t, p in fq2.get(sinif, []) if t in tset]
        if not dosyalar:
            sonuc.append(dict(hedef=hedef, sinif=sinif, ileri_dizi=F, geri_dizi=R,
                              urun_okuma=0, konsensus_uzunluk=0,
                              amaclanan=ic, olculen="fastq bulunamadi",
                              ozdeslik="", hizalanan="", uyum="olculemedi"))
            continue
        urunler = []
        tarandi = 0
        for g, t, p in dosyalar:
            n = 0
            with open(p, errors="replace") as fh:
                for j, line in enumerate(fh):
                    if j % 4 != 1:
                        continue
                    n += 1
                    if n > a.max_reads:
                        break
                    d = line.strip().upper()
                    if len(d) < 100:
                        continue
                    u = urun_kes(F, R, d, a)
                    if u:
                        urunler.append(u)
            tarandi += n
        kons, kac, uzcesit = baskin_konsensus(urunler)
        if not kons or kac < a.min_product:
            sonuc.append(dict(hedef=hedef, sinif=sinif, ileri_dizi=F, geri_dizi=R,
                              urun_okuma=len(urunler), konsensus_uzunluk=0,
                              amaclanan=ic,
                              olculen=('there are few products (%d)' % len(urunler)),
                              ozdeslik="", hizalanan="", uyum="olculemedi"))
            print(u'   [%2d] %-34s product=%-6d the consensus could not be built'
                  % (i + 1, hedef[:33], len(urunler)))
            continue
        # blastn
        sorgu = os.path.join(calisma, "u%d.fa" % i)
        with open(sorgu, "w") as fh:
            fh.write(">u%d\n%s\n" % (i, kons))
        dbad = SINIF_DB.get(sinif)
        fna = os.path.join(a.db, dbad) if dbad else None
        olculen = ozd = hiz = ""
        if fna and os.path.exists(fna):
            db = fna
            if not (os.path.exists(fna + ".nin") or os.path.exists(fna + ".00.nin")):
                db = os.path.join(calisma, dbad)
                if not os.path.exists(db):
                    os.symlink(os.path.abspath(fna), db)
                subprocess.run(["makeblastdb", "-in", db, "-dbtype", "nucl"],
                               capture_output=True, text=True)
            pr = subprocess.run(
                ["blastn", "-query", sorgu, "-db", db, "-outfmt",
                 "6 pident length bitscore stitle", "-max_target_seqs", "5",
                 "-evalue", "1e-5", "-num_threads", str(a.threads)],
                capture_output=True, text=True)
            en = None
            for line in pr.stdout.splitlines():
                q = line.split("\t")
                if len(q) < 4:
                    continue
                if en is None or float(q[2]) > float(en[2]):
                    en = q
            if en:
                m = re.match(r"\S+\s+(\S+\s+\S+)", en[3])
                olculen = m.group(1) if m else en[3][:40]
                ozd, hiz = en[0], en[1]
        # Uyum DUZEYLI raporlanir. "Methanosarcina barkeri" hedefi icin
        # "Methanosarcina thermophila" cikmasi cins duzeyinde uyum, tur
        # duzeyinde uyumsuzluktur; ikisini "uyusuyor" diye birlestirmek,
        # tur ozgul diye tasarlanmis bir cifti dogrulanmis gosterir.
        amac_turler = [x.strip() for x in ic.split(",") if x.strip()]
        amac_cinsler = set(x.split()[0] for x in amac_turler if x)
        olc_cins = olculen.split()[0] if olculen else ""
        if not olculen:
            uyum = "olculemedi"
        elif any(olculen == t for t in amac_turler):
            uyum = "tur_uyusuyor"
        elif olc_cins and olc_cins in amac_cinsler:
            uyum = "cins_uyusuyor_tur_farkli"
        else:
            uyum = "CINS_FARKLI"
        sonuc.append(dict(hedef=hedef, sinif=sinif, ileri_dizi=F, geri_dizi=R,
                          urun_okuma=len(urunler), konsensus_uzunluk=len(kons),
                          amaclanan=ic, olculen=olculen, ozdeslik=ozd,
                          hizalanan=hiz, uyum=uyum))
        print(u'   [%2d] %-30s product=%-6d cons=%-4d  %-28s %%%-6s %s'
              % (i + 1, hedef[:29], len(urunler), len(kons), olculen[:27],
                 ozd, uyum))

    shutil.rmtree(calisma, ignore_errors=True)
    d = os.path.dirname(os.path.abspath(a.out))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sonuc[0].keys()),
                           delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(sonuc)
    print("\nyazildi: %s" % a.out)
    say = collections.Counter(x["uyum"] for x in sonuc)
    for k, v in say.most_common():
        print("   %-14s %d" % (k, v))


if __name__ == "__main__":
    main()
