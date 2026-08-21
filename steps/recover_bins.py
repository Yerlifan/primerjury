#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
recover_bins.py
REFERANSSIZ KONSENSUS: kendi self referansi kullanilamayan kutulari kurtarir.

Sorun: 06/07 her kutu icin bir "self referans" uretiyor ve 12 o referansa
hizalayarak baskin alel konsensusunu kuruyor. Referansin kendisi asiri
belirsizse (cok IUPAC kodu) minimap2 o diziden neredeyse hic k-mer
cikaramiyor. Olculdu: B-1_2233851 referansi 1454 bp ve %19'u IUPAC;
minimap2 yalnizca 2 minimizer buldu ve konsensus sifir uzunlukta cikti.
Kutuda 5914 okuma var, yani veri kaybi degil, kalip kaybi.

Cozum: kalip disaridan gelmesin, okumalarin kendisinden gelsin.
  1. Medyan uzunluga yakin okumalardan aday tohumlar secilir.
  2. Her tohuma butun okumalar (altorneklem) hizalanir; en cok okuma
     toplayan tohum kalip olur. Boylece kalip en temsili okumadir.
  3. O kalip uzerinde baskin alel cagrilir (12 ile ayni kural: derinlik
     esigi, oran esigi, esitlikte N).
  4. Sonuc yeni kalip yapilip 2 ve 3 bir kez daha calistirilir; tohum
     okumasinin kendi hatalari boylece silinir.

Iki bagimsiz olcum: okumalar iki ayrik yariya bolunur ve konsensus iki
kez ayri ayri kurulur. Iki yari birbirinden ayrisiyorsa sonuc SUPHELI
isaretlenir ve fark raporlanir; sessizce tek bir dizi uretilmez.

Kullanim:
  python3 recover_bins.py --fastq "fastq files/B-1/B-1-reads_2233851.fastq" \
      --etiket B-1_2233851 --out referans_konsensus/baskin/konsensus
"""
import argparse, os, random, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import alignment

if hizalama.ARKA_UC is None:
    sys.exit(hizalama.durum())

BAZLAR = "ACGT"


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--fastq", required=True)
    p.add_argument("--etiket", required=True, help="example: B-1_2233851")
    p.add_argument("--out", required=True, help="consensus directory")
    p.add_argument("--tohum-sayisi", type=int, default=25)
    p.add_argument("--max-okuma", type=int, default=3000)
    p.add_argument("--min-uzunluk", type=int, default=400)
    p.add_argument("--min-derinlik", type=int, default=20)
    p.add_argument("--min-oran", type=float, default=0.50)
    p.add_argument("--tur", type=int, default=2, help="kac kez iyilestirilsin")
    p.add_argument("--yigin", type=int, default=2000)
    p.add_argument("--ayrisma-esigi", type=float, default=0.005,
                   help="iki yari arasindaki izin verilen fark orani")
    p.add_argument("--min-kapsam", type=float, default=0.95,
                   help="iki yarinin hizalamasinda beklenen en az kapsama")
    p.add_argument("--indel-uyari", type=float, default=0.005,
                   help="print a note for indel differences above this fraction")
    p.add_argument("--tohum", type=int, default=20260801,
                   help="random seed, so the same input gives the same result")
    return p.parse_args()


def fastq_oku(yol, min_uz, azami):
    okumalar = []
    with open(yol, errors="replace") as fh:
        for i, line in enumerate(fh):
            if i % 4 != 1:
                continue
            r = line.strip().upper()
            if len(r) >= min_uz:
                okumalar.append(r)
                if len(okumalar) >= azami:
                    break
    return okumalar


def kalipta_say(kalip, okumalar, yigin):
    """Okumalari kalibin uzerine yigar. Doner: (sayim_matrisi, hizalanan)."""
    A = hizalama.Hizalayici(seq=kalip, preset="map-ont")
    if not A:
        return None, 0
    say = [[0, 0, 0, 0] for _ in kalip]
    hiz = 0
    ad_okuma = {"o%d" % i: r for i, r in enumerate(okumalar)}
    adlar = list(ad_okuma)
    for bas in range(0, len(adlar), yigin):
        parca = {k: ad_okuma[k] for k in adlar[bas:bas + yigin]}
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
                        if rp + k < len(kalip) and qp + k < len(q):
                            j = BAZLAR.find(q[qp + k])
                            if j >= 0:
                                say[rp + k][j] += 1
                    rp += ln; qp += ln
                elif op == 1:
                    qp += ln
                elif op == 2:
                    rp += ln
                elif op == 4:
                    qp += ln
    return say, hiz


def konsensus_cagir(say, min_derinlik, min_oran):
    out = []
    dusuk = belirsiz = 0
    for c in say:
        d = sum(c)
        if d == 0 or d < min_derinlik:
            out.append("N"); dusuk += 1
            continue
        en = max(c)
        esitler = [x for x in range(4) if c[x] == en]
        if len(esitler) > 1:
            out.append("N"); belirsiz += 1     # esitlikte baskin alel yoktur
            continue
        j = esitler[0]
        if c[j] / d < min_oran:
            out.append("N"); belirsiz += 1
            continue
        out.append(BAZLAR[j])
    return "".join(out), dusuk, belirsiz


def kur(okumalar, a, gunluk):
    """Tohum sec, kalip kur, iyilestir. Doner: (dizi, tohum_hizalanan)."""
    uzunluklar = sorted(len(r) for r in okumalar)
    orta = statistics.median(uzunluklar)
    # medyan uzunluga en yakin okumalar aday tohum
    adaylar = sorted(okumalar, key=lambda r: abs(len(r) - orta))[:a.tohum_sayisi]
    ornek = okumalar if len(okumalar) <= 400 else random.sample(okumalar, 400)
    en_iyi, en_iyi_hiz = None, -1
    for t in adaylar:
        _, hiz = kalipta_say(t, ornek, a.yigin)
        if hiz > en_iyi_hiz:
            en_iyi, en_iyi_hiz = t, hiz
    gunluk("   tohum secildi: %d bp, %d/%d ornek okuma hizalandi"
           % (len(en_iyi), en_iyi_hiz, len(ornek)))
    kalip = en_iyi
    for tur in range(a.tur):
        say, hiz = kalipta_say(kalip, okumalar, a.yigin)
        if say is None:
            return None, 0
        dizi, dusuk, belirsiz = konsensus_cagir(say, a.min_derinlik, a.min_oran)
        # bastaki ve sondaki N'ler kirpilir, ic N'ler yerinde birakilir
        b = 0
        while b < len(dizi) and dizi[b] == "N":
            b += 1
        s = len(dizi)
        while s > b and dizi[s - 1] == "N":
            s -= 1
        dizi = dizi[b:s]
        gunluk("   tur %d: hizalanan=%d uzunluk=%d kapsanan=%d dusuk=%d belirsiz=%d"
               % (tur + 1, hiz, len(dizi), len(dizi) - dizi.count("N"),
                  dusuk, belirsiz))
        if not dizi:
            return None, hiz
        kalip = dizi.replace("N", "A") if "N" in dizi else dizi
        # bir sonraki turda kalip olarak kullanilacak; N'ler hizalamayi
        # bozmasin diye gecici olarak doldurulur, cagrilan diziye dokunulmaz
        son = dizi
    return son, hiz


def fark_olc(x, y):
    """Iki dizinin farkini IKI AYRI SAYI olarak verir: (ikame, indel).

    Ayrimin sebebi olculdu: bu kutunun iki yarisi 1449 bazda yalnizca 3
    ikame ile ayrildi (%0,21), ama 32 insersiyon ve 32 delesyonla. Indel
    farki homopolimer uzunlugundan geliyor; kalip tabanli baskin alel
    yontemi kalip pozisyonlarindaki bazi duzeltir, kalibin uzunlugunu
    duzeltmez, dolayisiyla tohum okumanin homopolimer uzunlugu sonuca
    gecer. Ikisini tek sayida toplamak, biyolojik sinyali ONT'nin bilinen
    homopolimer hatasiyla karistirir.

    Yon farkindan bagimsizdir: tohum okumalardan secildigi icin iki yarinin
    konsensusu birbirinin ters tumleyeni olabilir.
    Doner: (ikame_orani, indel_orani, kapsam, hizalanan_baz)"""
    if not x or not y:
        return 1.0, 1.0, 0.0, 0

    def olc(a, b):
        try:
            A = hizalama.Hizalayici(seq=a, preset="map-ont")
        except RuntimeError:
            return None
        if not A:
            return None
        en_iyi = None
        for h in A.map(b):
            if en_iyi is None or h.blen > en_iyi.blen:
                en_iyi = h
        if en_iyi is None:
            return None
        if en_iyi.strand == 1:
            q, qp = b, en_iyi.q_st
        else:
            q, qp = hizalama.revcomp(b), len(b) - en_iyi.q_en
        rp = en_iyi.r_st
        ortak = ikame = indel = 0
        # mlen'e guvenilmez; karakter karakter sayilir
        for ln, op in en_iyi.cigar:
            if op == 0:
                for k in range(ln):
                    if rp + k < len(a) and qp + k < len(q):
                        u, v = a[rp + k], q[qp + k]
                        if u == "N" or v == "N":
                            continue
                        ortak += 1
                        if u != v:
                            ikame += 1
                rp += ln; qp += ln
            elif op == 1:
                indel += ln; qp += ln
            elif op == 2:
                indel += ln; rp += ln
            elif op == 4:
                qp += ln
        if ortak == 0:
            return None
        kapsam = min(1.0, ortak / min(len(a), len(b)))
        # UC SAYI AYRI VERILIR. Bunlari tek sayida toplamak, kapsama
        # eksigini ikame farki gibi gostermeye yol aciyordu: bu kutuda
        # gercek ikame orani 0,0021 iken birlestirilmis sayi 0,0145
        # cikiyordu, yani karar biyolojik sinyale degil indellerin yedigi
        # pozisyonlara dayaniyordu.
        return (ikame / ortak, indel / ortak, kapsam, ortak)

    d = [z for z in (olc(x, y), olc(y, x)) if z is not None]
    return min(d, key=lambda z: z[0]) if d else (1.0, 1.0, 0.0, 0)


def main():
    a = get_args()
    random.seed(a.tohum)
    if not os.path.exists(a.fastq):
        sys.exit(u'no fastq found: %s' % a.fastq)
    os.makedirs(a.out, exist_ok=True)
    print(hizalama.durum())
    okumalar = fastq_oku(a.fastq, a.min_uzunluk, a.max_okuma)
    print(u'reads: %d (>= %d bp)' % (len(okumalar), a.min_uzunluk))
    if len(okumalar) < a.min_derinlik * 2:
        sys.exit(u'the read count is not enough: %d' % len(okumalar))

    def gun(s):
        print(s)

    print("TAM KUME")
    tam, hiz = kur(okumalar, a, gun)
    if not tam:
        sys.exit(u'the consensus could not be built')

    # iki bagimsiz olcum: ayrik yarilar
    karisik = okumalar[:]
    random.shuffle(karisik)
    yarim = len(karisik) // 2
    print(u'HALF 1 (%d reads)' % yarim)
    y1, _ = kur(karisik[:yarim], a, gun)
    print(u'HALF 2 (%d reads)' % (len(karisik) - yarim))
    y2, _ = kur(karisik[yarim:], a, gun)

    if y1 and y2:
        f, fi, kaps, ortak = fark_olc(y1, y2)
        nedenler = []
        if f > a.ayrisma_esigi:
            nedenler.append("ikame")
        if kaps < a.min_kapsam:
            nedenler.append("kapsama")
        uyum = "iki_olcum_uyustu" if not nedenler else (
            "ayrisan_olcum(%s)" % "+".join(nedenler))
        print("\niki yarinin karsilastirmasi (%d baz ortusuyor)" % ortak)
        print(u'   substitution difference : %.4f   (threshold %.4f)' % (f, a.ayrisma_esigi))
        print(u'   the indel difference : %.4f   (for information, it does not enter the decision)' % fi)
        print(u'   containment             : %.4f   (threshold %.4f)' % (kaps, a.min_kapsam))
        print(u'   result                  : %s' % uyum)
        if fi > a.indel_uyari:
            print(u'   NOTE: the indel difference comes from homopolymer length. A template based method does not correct template length. The design rules already reject a run of more than four identical bases, so no primer sits on those regions, but the consensus length is not exact.')
    else:
        uyum = "yari_kurulamadi"
        f = fi = 1.0; kaps = 0.0
        print(u'\nWARNING: one of the halves produced no consensus, so the split could not be measured')

    kapsanan = len(tam) - tam.count("N")
    yol = os.path.join(a.out, "%s_baskin_konsensus.fasta" % a.etiket)
    with open(yol, "w", encoding="utf-8") as fh:
        fh.write(u'>%s dominant_allele_no_reference reads=%d aligned=%d min_depth=%d min_fraction=%.2f length=%d covered=%d half_substitution_diff=%.4f half_indel_diff=%.4f half_coverage=%.4f %s\n'
                 % (a.etiket, len(okumalar), hiz, a.min_derinlik, a.min_oran,
                    len(tam), kapsanan, f, fi, kaps, uyum))
        for k in range(0, len(tam), 70):
            fh.write(tam[k:k + 70] + "\n")
    print("yazildi: %s" % yol)
    print("uzunluk=%d kapsanan=%d" % (len(tam), kapsanan))
    if uyum.startswith("ayrisan_olcum"):
        print(u'\nCAUTION: the two halves split at the SUBSTITUTION level. This bin may not come from a single organism, and should be inspected before its consensus is used.')
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
