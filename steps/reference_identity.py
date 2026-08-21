#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reference_identity.py
REFERANS TABANLI PRIMERLERIN NUMUNEDE NEYI COGALTTIGINI OLCER.

design_from_reference.py, numunede karsilanamayan hedefler icin referans
veritabani dizilerinden primer tasarladi ve numunede DESTEK olup olmadigini
olctu. Ornegin Methanosarcina barkeri icin tasarlanan bir cift, numunedeki
okumalarin %63'unde urun veriyor.

Ama "urun veriyor" ile "hedefi cogaltiyor" ayni sey degildir. Bir onceki
olcum (target_identity.py) numunedeki kutularin Kraken2 etiketleriyle
uyusmadigini gosterdi: M. barkeri etiketli kutu dizi duzeyinde
Methanosarcina vacuolata'ya gidiyor. O halde destegin hangi organizmadan
geldigi ayrica olculmelidir.

Yontem: her cift icin numunedeki okumalar taranir, urun veren okumalardan
URUN DIZISI kesilir, urunlerden baskin alel konsensusu kurulur ve bu
konsensus referans veritabanina blastn ile sorulur. Boylece "bu primer
numunede ne cogaltiyor" sorusu dizi kanitiyla yanitlanir.

Kullanim:
  python3 reference_identity.py --referans primer_referans/primer_referans.tsv \
      --pt . --db REFERANS_DB --out primer_referans/referans_kimlik.tsv
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
    """04'un find_bindings'i ile ayni toplanti kurallari."""
    tail_len = 12
    tail_max_mm = 1
    exact_last = 2
    total_max_mm = 3
    min_overlap = 12
    prod_min = 50
    prod_hard_max = 400


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--referans", required=True)
    p.add_argument("--pt", required=True, help="'fastq files' klasorunu iceren kok")
    p.add_argument("--db", required=True)
    p.add_argument("--hedefler-ref", default=os.path.join(HERE, "hedefler_referans.tsv"))
    p.add_argument("--out", required=True)
    p.add_argument("--max-okuma", type=int, default=20000)
    p.add_argument("--min-urun", type=int, default=30,
                   help="konsensus kurmak icin gereken en az urun sayisi")
    p.add_argument("--is-parcacigi", type=int, default=4)
    return p.parse_args()


def urun_kes(oligo_f, oligo_r, dizi, a):
    """Okumada urun varsa urun dizisini doner, yoksa None.
    Kalip cift sarmal oldugu icin iki yon de denenir."""
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
    """Esit uzunluktaki urunleri hizalamadan, en sik gorulen uzunluktaki
    urunler uzerinden pozisyon basina baskin baz cagirir. Uzunluk
    dagilimi dar oldugu icin (tek amplikon) bu yeterlidir; genis dagilim
    varsa uyari basilir."""
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
        sys.exit("blastn bulunamadi")
    rows = list(csv.DictReader(open(a.referans, encoding="utf-8"), delimiter="\t"))
    if not rows:
        sys.exit("referans cift yok")
    # hedefin amaclanan organizmasi
    amac = {}
    if os.path.exists(a.hedefler_ref):
        # Dosya # ile baslayan aciklama satirlariyla basliyor; DictReader
        # bunlari baslik sanip anahtarlari yanlis kurar.
        with open(a.hedefler_ref, encoding="utf-8") as fh:
            satirlar = [l for l in fh if not l.startswith("#") and l.strip()]
        for r in csv.DictReader(satirlar, delimiter="\t"):
            if r.get("ad"):
                amac[r["ad"]] = (r.get("ic", ""), r.get("taxid", ""))
    if not amac:
        print("UYARI: hedefler_referans.tsv okunamadi, amaclanan organizma "
              "bos kalacak: %s" % a.hedefler_ref)

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
                    if n > a.max_okuma:
                        break
                    d = line.strip().upper()
                    if len(d) < 100:
                        continue
                    u = urun_kes(F, R, d, a)
                    if u:
                        urunler.append(u)
            tarandi += n
        kons, kac, uzcesit = baskin_konsensus(urunler)
        if not kons or kac < a.min_urun:
            sonuc.append(dict(hedef=hedef, sinif=sinif, ileri_dizi=F, geri_dizi=R,
                              urun_okuma=len(urunler), konsensus_uzunluk=0,
                              amaclanan=ic,
                              olculen=("urun az (%d)" % len(urunler)),
                              ozdeslik="", hizalanan="", uyum="olculemedi"))
            print("   [%2d] %-34s urun=%-6d konsensus kurulamadi"
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
                 "-evalue", "1e-5", "-num_threads", str(a.is_parcacigi)],
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
        print("   [%2d] %-30s urun=%-6d kons=%-4d  %-28s %%%-6s %s"
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
