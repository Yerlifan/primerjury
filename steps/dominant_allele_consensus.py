#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dominant_allele_consensus.py
Her konsensüs için, o kutunun kendi ham okumalarından BASKIN ALEL dizisini
ve pozisyon başına alel oranlarını üretir.

Neden gerekli: `samtools consensus -A` çıktısı değişken pozisyonları IUPAC
koduyla yazıyor. Bağlanma kuralı IUPAC'ı küme kesişimiyle değerlendirdiği
için, kodun taşıdığı belirsizlik primeri her iki alele de uyar gösteriyor.
Ölçtüm: F2 sınıfında Trichoderma asperellum ile Metarhizium brunneum
konsensüsleri karakter karakter %79,8 özdeş, ama kesişim ölçütüyle %99,9.
Yani gerçek fark var, belirsizlik onu görünmez kılıyor. Bu betik aynı
veriden ikinci ve belirsizliksiz bir ölçüm üretir; tasarım ve özgüllük iki
küme üzerinde ayrı ayrı çalıştırılıp sonuçlar karşılaştırılabilir.

Çıktı:
  <out>/konsensus/<etiket>_baskin_konsensus.fasta
  <out>/oran/<etiket>_alel.tsv       poz, derinlik, A, C, G, T, baskin, oran,
                                     wilson_alt
  <out>/ozet.tsv

Kullanım:
  python3 dominant_allele_consensus.py \
      --kons referans_konsensus/self/konsensus \
      --fastq "fastq files" --out referans_konsensus/baskin
"""
import argparse, csv, glob, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import alignment

if hizalama.ARKA_UC is None:
    sys.exit(hizalama.durum())

BAZLAR = "ACGT"


def oku_fasta(f):
    return "".join(l.strip() for l in open(f, encoding="utf-8", errors="replace")
                   if not l.startswith(">")).upper()


def wilson_alt(k, n, z=1.96):
    """Oranın Wilson alt sınırı. Tek okumalık gürültünün yüksek oran gibi
    görünmesini engeller; toplantı kararında rakip oranları için istenen
    ölçüt budur."""
    if n == 0:
        return 0.0
    p = k / n
    d = 1 + z * z / n
    m = p + z * z / (2 * n)
    s = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return max(0.0, (m - s) / d)


def fastq_bul(kok, grp, taxid):
    """Sikistirilmamis fastq aranir. .gz destegi yoktur; boyle bir dosya
    varsa sessizce atlamak yerine acikca uyarilir."""
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
    p.add_argument("--kons", required=True, help="self consensus directory")
    p.add_argument("--fastq", required=True, help="'fastq files' directory")
    p.add_argument("--out", required=True)
    p.add_argument("--max-okuma", type=int, default=3000)
    p.add_argument("--min-derinlik", type=int, default=20,
                   help="below this depth no base is called and N is written")
    p.add_argument("--min-oran", type=float, default=0.50,
                   help="write N when the dominant allele is below this fraction")
    p.add_argument("--min-uzunluk", type=int, default=400)
    p.add_argument("--oran-yaz", type=int, default=1)
    p.add_argument("--yigin", type=int, default=2000,
                   help="single call on the minimap2 command-line backend "
                        "verilecek okuma sayisi")
    return p.parse_args()


def main():
    a = get_args()
    os.makedirs(os.path.join(a.out, "konsensus"), exist_ok=True)
    if a.oran_yaz:
        os.makedirs(os.path.join(a.out, "oran"), exist_ok=True)
    ozet = []
    dosyalar = sorted(glob.glob(os.path.join(a.kons, "*_konsensus.fasta")))
    if not dosyalar:
        sys.exit("konsensus bulunamadi: %s" % a.kons)
    print(hizalama.durum())
    print(u'consensus files: %d' % len(dosyalar))
    for f in dosyalar:
        etiket = re.sub(r"_(ref|self)_konsensus\.fasta$", "", os.path.basename(f))
        m = re.match(r"((?:A1|A2|B|F1|F2)-\d+)_(\d+)$", etiket)
        if not m:
            print("   ATLANDI, etiket cozulemedi: %s" % etiket)
            continue
        grp, taxid = m.group(1), m.group(2)
        ref = oku_fasta(f)
        # N'ler SILINMEZ. Silinirse ic kapsama bosluklarinin iki yani
        # birbirine yapisir ve dogada olmayan bir kavsak olusur; 08 bu diziyi
        # tasarim kalibi olarak kullandigi icin kavsagin uzerinden primer
        # secilebilir ve o primer hicbir zaman calismaz. Bastaki ve sondaki
        # N'ler kirpilir (koordinat kaymasi cikti basliginda raporlanir),
        # ic N'ler yerinde birakilir; minimap2 N'i zaten uyumsuz sayar.
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
            A = hizalama.Hizalayici(seq=cekirdek, preset="map-ont")
        except RuntimeError as e:
            sys.exit(str(e))
        if not A:
            print("   ATLANDI, indeks kurulamadi: %s" % etiket)
            continue
        say = [[0, 0, 0, 0] for _ in cekirdek]
        okumalar = {}
        n = 0
        with open(fq, errors="replace") as fh:
            for i, line in enumerate(fh):
                if i % 4 != 1:
                    continue
                r = line.strip().upper()   # kucuk harfli FASTQ da kabul edilir
                if len(r) < a.min_uzunluk:
                    continue
                n += 1
                if n > a.max_okuma:
                    break
                okumalar["o%d" % n] = r
        hiz = 0
        # Toplu hizalama: minimap2 komut satiri arka ucunda okuma basina
        # surec baslatmak kabul edilemez derecede yavas olurdu.
        for bas in range(0, len(okumalar), a.yigin):
            parca = dict(list(okumalar.items())[bas:bas + a.yigin])
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
                    q, qp = hizalama.revcomp(r), len(r) - h.q_en
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
            if d == 0 or d < a.min_derinlik:
                cikti.append("N"); dusuk += 1
                if a.oran_yaz:
                    oran_satir.append([i + 1, d, c[0], c[1], c[2], c[3], "N", 0.0, 0.0])
                continue
            en = max(c)
            esitler = [x for x in range(4) if c[x] == en]
            if len(esitler) > 1:
                # Tam esitlikte alfabetik sira belirleyici olmamali; iki alel
                # esit orandaysa baskin alel yoktur, N yazilir.
                cikti.append("N"); belirsiz += 1
                if a.oran_yaz:
                    oran_satir.append([i + 1, d, c[0], c[1], c[2], c[3], "N",
                                       round(en / d, 4), 0.0])
                continue
            j = esitler[0]
            oran = c[j] / d
            w = wilson_alt(c[j], d)
            if oran < a.min_oran:
                cikti.append("N"); belirsiz += 1
                baz = "N"
            else:
                baz = BAZLAR[j]
                cikti.append(baz)
            if a.oran_yaz:
                oran_satir.append([i + 1, d, c[0], c[1], c[2], c[3], baz,
                                   round(oran, 4), round(w, 4)])
        dizi = "".join(cikti)
        yol = os.path.join(a.out, "konsensus", "%s_baskin_konsensus.fasta" % etiket)
        with open(yol, "w", encoding="utf-8") as fh:
            fh.write(">%s baskin_alel okuma=%d hizalanan=%d min_derinlik=%d "
                     "min_oran=%.2f kirpma=%d-%d ic_N_girdide=%d\n"
                     % (etiket, n, hiz, a.min_derinlik, a.min_oran,
                        bas + 1, son, ic_n))
            for k in range(0, len(dizi), 70):
                fh.write(dizi[k:k + 70] + "\n")
        if a.oran_yaz:
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
        print("   %-26s okuma=%5d hizalanan=%5d uzunluk=%5d kapsanan=%5d "
              "dusuk_derinlik=%4d belirsiz=%4d"
              % (etiket, n, hiz, len(dizi), kaps, dusuk, belirsiz))
    if ozet:
        with open(os.path.join(a.out, "ozet.tsv"), "w", newline="",
                  encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(ozet[0].keys()),
                               delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(ozet)
        print("\nyazildi: %s" % os.path.join(a.out, "ozet.tsv"))


if __name__ == "__main__":
    main()
