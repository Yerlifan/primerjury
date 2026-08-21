#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
target_identity.py
HEDEF ADI ILE VERININ GOSTERDIGI ORGANIZMAYI KARSILASTIRIR.

Toplanti kararindaki hedef adlari Kraken2'nin takson atamalarina dayaniyor.
Kraken2 varsayilan ayarinda cekimser kalmaz: gercek organizma veritabaninda
yoksa okuma en yuksek puanli yaprakla etiketlenir. Bu yuzden bir kutunun
etiketi ile o kutudaki DIZININ gercekte ne oldugu ayrilabilir.

Bu betik kutu konsensuslerini referans veritabanina blastn ile sorgular ve
her hedef icin ucunu yan yana koyar:
  toplanti_adi     kararda yazan hedef adi
  kraken_etiketi   o hedefin taxid'lerinin adlari
  olculen_kimlik   konsensusun referans veritabanindaki en iyi eslesmesi

Ayirt edici bolge tercih edilir: mantarda ITS, bakteri ve arkede 16S.
Korunmus bolgeler (18S, 28S) cins duzeyinde ayirmaz; betik hangi bolgeden
olculdugunu her satirda yazar, boylece guven duzeyi gizlenmez.

Kullanim:
  python3 target_identity.py --kons referans_konsensus/baskin/konsensus \
      --db REFERANS_DB --hedefler hedefler.tsv --adlar taxid_adlari.tsv \
      --out primer_final/hedef_kimlik.tsv
"""
import argparse, csv, collections, glob, os, re, subprocess, sys, tempfile, shutil

# Hangi sinif hangi veritabanina sorulur. Ilk sirada AYIRT EDICI bolge
# bulunur; sonrakiler korunmus bolgelerdir ve yalnizca ilki bos donerse
# kullanilir.
SINIF_DB = {
    "A1": [("archaea.16S.fna", "16S")],
    "A2": [("archaea.16S.fna", "16S")],
    "B":  [("bacteria.16S.fna", "16S")],
    "F1": [("fungi.ITS.fna", "ITS"), ("fungi.28SrRNA.fna", "28S"),
           ("fungi.18SrRNA.fna", "18S")],
    "F2": [("fungi.ITS.fna", "ITS"), ("fungi.28SrRNA.fna", "28S"),
           ("fungi.18SrRNA.fna", "18S")],
}


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--kons", required=True)
    p.add_argument("--db", required=True)
    p.add_argument("--hedefler", default="hedefler.tsv")
    p.add_argument("--adlar", default="taxid_adlari.tsv")
    p.add_argument("--out", required=True)
    p.add_argument("--min-hizalama", type=int, default=250,
                   help="bundan kisa hizalamalar kimlik sayilmaz")
    p.add_argument("--min-ozdeslik", type=float, default=90.0)
    p.add_argument("--is-parcacigi", type=int, default=4)
    return p.parse_args()


def oku_fasta(f):
    return "".join(l.strip() for l in open(f, encoding="utf-8", errors="replace")
                   if not l.startswith(">")).upper()


def db_hazirla(fna, calisma):
    if os.path.exists(fna + ".nin") or os.path.exists(fna + ".00.nin"):
        return fna
    hedef = os.path.join(calisma, os.path.basename(fna))
    if not os.path.exists(hedef):
        os.symlink(os.path.abspath(fna), hedef)
    if not os.path.exists(hedef + ".nin"):
        r = subprocess.run(["makeblastdb", "-in", hedef, "-dbtype", "nucl"],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print("   makeblastdb basarisiz: %s" % r.stderr.strip()[:180])
            return None
    return hedef


def main():
    a = get_args()
    if not shutil.which("blastn"):
        sys.exit("blastn bulunamadi. Kurmak icin: "
                 "sudo apt-get install -y ncbi-blast+")
    ad = {}
    if os.path.exists(a.adlar):
        for l in open(a.adlar, encoding="utf-8"):
            q = l.rstrip("\n").split("\t")
            if len(q) > 1:
                ad[q[0]] = q[1]
    hedef_taxid = {}
    for l in open(a.hedefler, encoding="utf-8"):
        if l.startswith("#") or not l.strip():
            continue
        q = l.rstrip("\n").split("\t")
        if len(q) < 4 or q[0] == "karar":
            continue
        hedef_taxid[q[1]] = [t.strip() for t in q[3].split(",") if t.strip()]

    # kutu envanteri
    kutular = {}
    for p in sorted(glob.glob(os.path.join(a.kons, "*.fasta"))):
        et = re.sub(r"_(baskin|ref|self)?_?konsensus\.fasta$", "",
                    os.path.basename(p))
        m = re.match(r"((?:A1|A2|B|F1|F2))-\d+_(\d+)$", et)
        if m:
            kutular[et] = (m.group(1), m.group(2), oku_fasta(p))
    if not kutular:
        sys.exit("konsensus bulunamadi: %s" % a.kons)
    print(u'bins: %d' % len(kutular))

    calisma = tempfile.mkdtemp(prefix="kimlik_")
    # sinif basina tek blastn cagrisi
    kutu_kimlik, zayif_kimlik = {}, {}
    for sinif in sorted(set(v[0] for v in kutular.values())):
        etler = [e for e, v in kutular.items() if v[0] == sinif]
        sorgu = os.path.join(calisma, "%s.fa" % sinif)
        with open(sorgu, "w", encoding="utf-8") as fh:
            for e in etler:
                fh.write(">%s\n%s\n" % (e, kutular[e][2].replace("N", "")))
        for dbad, bolge in SINIF_DB.get(sinif, []):
            kalan = [e for e in etler if e not in kutu_kimlik]
            if not kalan:
                break
            fna = os.path.join(a.db, dbad)
            if not os.path.exists(fna):
                print(u'no such database: %s' % dbad)
                continue
            db = db_hazirla(fna, calisma)
            if not db:
                continue
            cikti = sorgu + "." + dbad + ".tsv"
            r = subprocess.run(
                ["blastn", "-query", sorgu, "-db", db, "-outfmt",
                 "6 qseqid pident length bitscore stitle",
                 "-max_target_seqs", "5",
                 "-evalue", "1e-20", "-num_threads", str(a.is_parcacigi),
                 "-out", cikti], capture_output=True, text=True)
            print(u'   blastn %-3s x %-20s (%d bins)' % (sinif, dbad, len(kalan)))
            if r.returncode != 0:
                print("      HATA: %s" % r.stderr.strip()[:160])
                continue
            en, zayif_en = {}, {}
            for line in open(cikti, encoding="utf-8"):
                q = line.rstrip("\n").split("\t")
                if len(q) < 5:
                    continue
                et, pid, aln, bit, tit = (q[0], float(q[1]), int(q[2]),
                                          float(q[3]), q[4])
                zayif = aln < a.min_hizalama or pid < a.min_ozdeslik
                if zayif:
                    # Esigi gecemeyen en iyi vurus da tutulur. "eslesme yok"
                    # ile "en yakin akraba %88'de" ayni sey degildir; ilki
                    # veri yoklugu, ikincisi veritabaninda yakin akrabasi
                    # bulunmayan bir organizma demektir.
                    if et not in zayif_en or bit > zayif_en[et][0]:
                        zayif_en[et] = (bit, pid, aln, tit, bolge)
                    continue
                # BITSCORE ile siralanir. Once uzunluga bakmak yanilticiydi:
                # 524 bp'de %94,47 eslesme, 504 bp'de %98,21 eslesmenin
                # onune geciyordu ve kimlik yanlis cikiyordu. Bitscore
                # uzunlukla ozdesligi birlikte tartar, blastn'in kendi
                # olcutudur.
                if et not in en or bit > en[et][0]:
                    en[et] = (bit, pid, aln, tit, bolge)
            for et, (_, pid, aln, tit, bl) in en.items():
                if et not in kutu_kimlik:
                    kutu_kimlik[et] = (pid, aln, tit, bl, True)
            for et, (_, pid, aln, tit, bl) in zayif_en.items():
                if et not in kutu_kimlik and et not in zayif_kimlik:
                    zayif_kimlik[et] = (pid, aln, tit, bl, False)

    shutil.rmtree(calisma, ignore_errors=True)

    # hedef basina topla
    sonuc = []
    for hedef, tl in sorted(hedef_taxid.items()):
        ilgili = [e for e, v in kutular.items() if v[1] in tl]
        if not ilgili:
            continue
        kraken = sorted(set(ad.get(v, v) for v in tl
                            if any(kutular[e][1] == v for e in ilgili)))
        say = collections.Counter()
        detay = {}
        for e in ilgili:
            k = kutu_kimlik.get(e)
            guclu = True
            if not k:
                k = zayif_kimlik.get(e)
                guclu = False
            if not k:
                say["veritabaninda hic vurus yok"] += 1
                continue
            pid, aln, tit, bl = k[0], k[1], k[2], k[3]
            # referans basligindan tur adini cikar: NR_xxxxx.1 Cins tur ...
            m = re.match(r"\S+\s+(\S+\s+\S+)", tit)
            adi = m.group(1) if m else tit[:40]
            if not guclu:
                adi = adi + " (esik alti)"
            say[adi] += 1
            detay.setdefault(adi, []).append("%s %%%.2f/%dbp/%s%s"
                                             % (e, pid, aln, bl,
                                                "" if guclu else " ZAYIF"))
        if not say:
            continue
        baskin, n = say.most_common(1)[0]
        # Uyum ikili degil, DUZEYLI raporlanir: cins uyusup turun farkli
        # olmasi ile hic uyusmamak ayni sey degildir.
        if "vurus yok" in baskin:
            uyum = "vurus_yok"
        elif baskin.endswith("(esik alti)"):
            uyum = "YAKIN_AKRABA_YOK"
        elif any(baskin == k for k in kraken):
            uyum = "tur_uyusuyor"
        elif any(baskin.split()[0] == k.split()[0] for k in kraken if k):
            uyum = "cins_uyusuyor_tur_farkli"
        else:
            uyum = "CINS_FARKLI"
        sonuc.append(dict(
            hedef=hedef,
            kraken_etiketi="; ".join(kraken),
            olculen_kimlik=baskin,
            kutu_sayisi=len(ilgili),
            destekleyen_kutu=n,
            uyum=uyum,
            kanit="; ".join(detay.get(baskin, [])[:4]),
            diger=("; ".join("%s(%d)" % (k, v)
                             for k, v in say.most_common()[1:4]) or "")))

    if not sonuc:
        sys.exit("hicbir hedef icin kimlik olculemedi")
    d = os.path.dirname(os.path.abspath(a.out))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sonuc[0].keys()),
                           delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(sonuc)
    print("\nyazildi: %s" % a.out)
    say_uyum = collections.Counter(x["uyum"] for x in sonuc)
    print(u'targets: %d' % len(sonuc))
    for k in ("tur_uyusuyor", "cins_uyusuyor_tur_farkli", "CINS_FARKLI",
              "YAKIN_AKRABA_YOK", "vurus_yok"):
        if say_uyum.get(k):
            print("   %-26s %d" % (k, say_uyum[k]))
    print("\n%-34s %-30s %-34s %s"
          % ("HEDEF", "KRAKEN2 ETIKETI", "OLCULEN KIMLIK", "UYUM"))
    for x in sonuc:
        print("%-34s %-30s %-34s %s"
              % (x["hedef"][:33], x["kraken_etiketi"][:29],
                 x["olculen_kimlik"][:33], x["uyum"]))


if __name__ == "__main__":
    main()
