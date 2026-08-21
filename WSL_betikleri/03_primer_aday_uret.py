#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03_primer_aday_uret.py
Toplantı kararlarındaki oligo, termodinamik, ürün ve bölge kurallarını
uygulayarak bir hedef konsensüsünden primer çifti adayları üretir.

Tasarım ilkeleri:
  * Her Tm iki bağımsız kütüphaneyle ölçülür (primer3 ve Biopython). İki
    kütüphane arasındaki SİSTEMATİK KAYMA veriden ölçülür (tüm adayların
    farkının medyanı), sonra bu kaymadan tolerans kadar sapan oligo elenir.
    Kayma elle yazılmaz.
  * Hiçbir eşik veriye bakılmadan sabitlenmez; hepsi komut satırından
    değiştirilebilir ve kullanılan değer log'a yazılır.
  * Ürün doğrulaması makine tarafından yapılır: kalıptan kesilen parçanın
    başı ileri primere, sonu geri primerin ters tümleyenine birebir eşit
    olmalı. Eşit değilse aday sessizce elenmez, sayacı artar.
  * Maske dosyasındaki hiçbir pozisyon primer ayak izine giremez. Ürünün
    iç kısmındaki maskeli pozisyonlar serbesttir, çünkü kural yalnızca
    primer yerleşimini kısıtlar.

KONSENSÜS HANGİ KLASÖRDEN OKUNMALI (2026-08-21 düzeltmesi)
  Yalnız `konsensus_kanonik/` kullanın. Eski örnekler `consensus sequences/`
  klasörünü gösteriyordu; o klasör KARIŞIK YÖNLÜDÜR (ölçülen: 71 antisense /
  27 sense). Ters yönlü bir konsensüste in-silico PCR SESSİZCE 0 ürün verir —
  ölçülen kayıp %100, kanıt `KAPSAMLI_ARAMA/yon_etki_testi.py`. Yasak
  `KAPSAMLI_ARAMA/yapilandirma.py` içinde de yazılıdır (KONSENSUS_KANONIK).
  Kanonik klasör `KAPSAMLI_ARAMA/kanonik_uret.py` ile üretilir. Hangi kutunun
  hangi dosyaya karşılık geldiği `konsensus_kanonik/INDEKS.tsv` içindedir;
  dosya adını tahmin etmeyin, indeksten okuyun.

Kullanım:
  python3 03_primer_aday_uret.py \
      --consensus "konsensus_kanonik/A1-1_2209.kanonik.fa" \
      --mask      "N_analizi/maske/A1-1_2209_maske.bed" \
      --out       "primer_adaylari/A1-1_2209.tsv"
"""
import argparse, csv, itertools, os, statistics, sys

try:
    import primer3
except ImportError:
    sys.exit("primer3-py kurulu degil:  pip install primer3-py")
try:
    from Bio.Seq import Seq
    from Bio.SeqUtils import MeltingTemp as mt
except ImportError:
    sys.exit("biopython kurulu degil:  pip install biopython")

COMP = str.maketrans("ACGTNRYSWKMBDHV", "TGCANYRSWMKVHDB")


def rc(s):
    return s.translate(COMP)[::-1]


# ----------------------------------------------------------------- argümanlar
def get_args():
    p = argparse.ArgumentParser(description="Primer cifti aday uretimi")
    p.add_argument("--consensus", required=True)
    p.add_argument("--mask", default=None, help="02 betiginin urettigi BED")
    p.add_argument("--mask-contig", default=None,
                   help="BED'in 1. sutununda bu ada sahip satirlar kullanilir; "
                        "verilmezse butun satirlar alinir")
    p.add_argument("--ambig", default=None, help="IUPAC kodlu konsensus (bilgi amacli)")
    p.add_argument("--out", required=True)
    p.add_argument("--label", default=None)
    # oligo kurallari
    p.add_argument("--len-min", type=int, default=18)
    p.add_argument("--len-max", type=int, default=25)
    p.add_argument("--gc-min", type=float, default=40.0)
    p.add_argument("--gc-max", type=float, default=60.0)
    p.add_argument("--gc-hard-min", type=float, default=35.0)
    p.add_argument("--gc-hard-max", type=float, default=65.0)
    p.add_argument("--gc-clamp-last", type=int, default=5,
                   help="son kac bazda GC sayilir")
    p.add_argument("--gc-clamp-max", type=int, default=3,
                   help="son N bazda en fazla kac G veya C")
    p.add_argument("--homopolymer-max", type=int, default=4)
    p.add_argument("--require-3p-gc", type=int, default=1,
                   help="1 ise 3' uc G ya da C ile bitmeli")
    p.add_argument("--degeneracy-budget", type=int, default=0,
                   help="ARTIK ETKISIZ. Toplanti karari geregi oligolar salt "
                        "ACGT uretiliyor; kalip belirsizligi --iupac-max ile "
                        "yonetiliyor. Bayrak geriye donuk uyum icin duruyor "
                        "ve verildiginde uyari basilir.")
    p.add_argument("--degeneracy-fold-max", type=int, default=4,
                   help="ARTIK ETKISIZ, geriye donuk uyum icin duruyor")
    p.add_argument("--iupac-max", type=int, default=2,
                   help="kalip penceresinde izin verilen IUPAC pozisyon sayisi; "
                        "bunlar somut baza cozulur, oligoya dejenere baz "
                        "girmez. 0 verilirse IUPAC iceren pencere hic "
                        "kullanilmaz.")
    p.add_argument("--iupac-son-yasak", type=int, default=5,
                   help="oligonun son bu kadar bazinda IUPAC kabul edilmez")
    # termodinamik
    p.add_argument("--tm-min", type=float, default=58.0)
    p.add_argument("--tm-max", type=float, default=62.0)
    p.add_argument("--tm-hard-min", type=float, default=57.0)
    p.add_argument("--tm-hard-max", type=float, default=63.0)
    p.add_argument("--tm-cross-k", type=float, default=4.0,
                   help="tolerans = k carpi kaymanin standart sapmasi; veriden "
                        "turetilir")
    p.add_argument("--tm-cross-tol", type=float, default=None,
                   help="toleransi elle sabitler, verilirse k yok sayilir")
    p.add_argument("--hairpin-dg-min", type=float, default=-3000.0)
    p.add_argument("--homodimer-dg-min", type=float, default=-6000.0)
    p.add_argument("--heterodimer-dg-min", type=float, default=-6000.0)
    # tampon kosullari, iki kutuphaneye de ayni verilir
    p.add_argument("--mv", type=float, default=50.0, help="tek degerlikli katyon mM")
    p.add_argument("--dv", type=float, default=1.5, help="iki degerlikli katyon mM")
    p.add_argument("--dntp", type=float, default=0.6, help="dNTP mM")
    p.add_argument("--dna-conc", type=float, default=50.0, help="oligo nM")
    # urun
    p.add_argument("--prod-min", type=int, default=70)
    p.add_argument("--prod-max", type=int, default=250)
    p.add_argument("--prod-hard-max", type=int, default=300)
    p.add_argument("--prod-best-min", type=int, default=90)
    p.add_argument("--prod-best-max", type=int, default=150)
    p.add_argument("--prod-gc-min", type=float, default=40.0)
    p.add_argument("--prod-gc-max", type=float, default=60.0)
    p.add_argument("--pair-tm-diff-max", type=float, default=1.5)
    p.add_argument("--max-pairs", type=int, default=5000,
                   help="cikti satiri ust siniri, puana gore kirpilir")
    p.add_argument("--min-locus-spacing", type=int, default=0,
                   help="ayni lokusun kaymis kopyalarini seyreltir: iki aday "
                        "cifti hem F hem R baslangici bu kadar yakinsa daha "
                        "kotu puanli olan atilir. 0 kapatir.")
    return p.parse_args()


# ------------------------------------------------------------------ yardimcilar
def read_fasta(path):
    """Tek kayitli FASTA okur. Cok kayitli dosyada durur, cunku kayitlarin
    sessizce birlestirilmesi kavsak uzerinde yapay kimerik primer uretir."""
    seq, names = [], []
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith(">"):
            names.append(line[1:].strip())
        else:
            seq.append(line.strip())
    if len(names) > 1:
        sys.exit("HATA: %s icinde %d kayit var. Bu betik tek kayitli konsensus "
                 "bekler; kayitlarin birlestirilmesi kavsakta yapay primer "
                 "uretir." % (path, len(names)))
    return (names[0] if names else None), "".join(seq).upper()


def read_mask(path, seqlen, contig=None, strict_missing=True):
    """BED -> yasak pozisyon kumesi (0 tabanli), sinif sayaci ve kontig raporu.
    contig verilirse yalnizca o kontige ait satirlar alinir; eslesen satir yoksa
    bu sessiz bir maskesizlige donusmesin diye uyari dondurulur."""
    bad, classes, seen = set(), {}, {}
    if not path:
        return bad, classes, seen
    if not os.path.exists(path):
        if strict_missing:
            sys.exit("HATA: maske dosyasi yok: %s. Yol verildiginde var olmasi "
                     "zorunludur, aksi halde maskeleme sessizce kapanir." % path)
        return bad, classes, seen
    for line in open(path, encoding="utf-8", errors="replace"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        try:
            st, e = int(f[1]), int(f[2])
        except ValueError:
            continue
        name = f[0]
        seen[name] = seen.get(name, 0) + 1
        if contig is not None and name != contig:
            continue
        cls = f[3] if len(f) > 3 else "maske"
        for i in range(max(0, st), min(seqlen, e)):
            bad.add(i)
        classes[cls] = classes.get(cls, 0) + (e - st)
    return bad, classes, seen


def gc_pct(s):
    return 100.0 * (s.count("G") + s.count("C")) / len(s) if s else 0.0


def max_run(s):
    best = run = 1
    for i in range(1, len(s)):
        run = run + 1 if s[i] == s[i - 1] else 1
        best = max(best, run)
    return best if s else 0


# N bilerek listede yok: konsensusteki N bir kapsama boslugudur, kasitli
# dejenerelik degildir. N iceren hicbir oligo uretilmez; primer3 de N'de
# ValueError firlatir, yani sessiz yanlis Tm riski de ortadan kalkar.
DEGEN_FOLD = {"R": 2, "Y": 2, "S": 2, "W": 2, "K": 2, "M": 2,
              "B": 3, "D": 3, "H": 3, "V": 3}

# IUPAC kodunun cozulecegi somut baz. Konsensus -A ile uretildigi icin
# degisken pozisyonlar R, Y, S, W, K, M olarak yaziliyor. Toplanti karari
# oligoda dejenere baz istemiyor; bu yuzden kod somut bir baza cozulur ve
# secimin dogru olup olmadigi 09'da ham okumalarda deneysel olarak sinanir.
# Cozum kurali sabit ve belirlenimci: kumenin alfabetik ilk bazi.
IUPAC_COZ = {"R": "A", "Y": "C", "S": "C", "W": "A", "K": "G", "M": "A",
             "B": "C", "D": "A", "H": "A", "V": "A"}

# Kodun temsil ettigi bazlarin tamami. Belirsiz pozisyonu tek bir baza
# indirgemek bilgi atmak olur; bunun yerine butun alternatifler uretilir,
# hepsi ayni suzgeclerden gecer ve hangisinin gercekten baglandigina
# 09'daki ham okuma taramasi karar verir.
IUPAC_KUME = {"R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT",
              "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG"}


def _iupac_denetle(win, a, uc="her"):
    """Pencere IUPAC kurallarina uyuyor mu. (belirsiz_pozisyonlar, sebep).

    uc: hangi zincirin oligosu uretilecek.
        "F"   -> oligo = win        , 3' uc pencerenin SONU
        "R"   -> oligo = rc(win)    , 3' uc pencerenin BASI
        "her" -> iki ucu da yasakla (geriye donuk, daha sıkı)

    Zincir ayrimi sart: ayni pencereden hem F hem R uretiliyor ve R'nin 3'
    ucu pencerenin BASINA denk geliyor. Tek yonlu denetim, R primerinin en
    uctaki bazinin belirsiz bir pozisyondan tek alele sessizce sabitlenmesine
    izin veriyordu. Olculdu: gercek veride 2000 satirin 435'inde R'nin 3'
    terminal bazi cozulmus bir IUPAC pozisyonundan geliyordu."""
    if "N" in win:
        return None, "kalipta_N"
    k = [i for i, c in enumerate(win) if c not in "ACGT"]
    if not k:
        return [], None
    if any(c not in IUPAC_KUME for c in win if c not in "ACGT"):
        return None, "tanimsiz_kod"
    if len(k) > a.iupac_max:
        return None, "iupac_fazla"
    n = a.iupac_son_yasak
    if uc in ("F", "her") and any(i >= len(win) - n for i in k):
        return None, "iupac_3p"
    if uc in ("R", "her") and any(i < n for i in k):
        return None, "iupac_3p"
    return k, None


def iupac_coz(win, a, uc="her"):
    """Geriye donuk uyum: tek, belirlenimci cozum dondurur."""
    k, why = _iupac_denetle(win, a, uc)
    if k is None:
        return None, why
    if not k:
        return win, 0
    return "".join(IUPAC_COZ.get(c, c) for c in win), len(k)


def iupac_varyantlar(win, a, uc="her"):
    """Pencerenin butun somut ACGT karsiliklarini uretir.
    ([(oligo, cozulen_pozisyon_sayisi), ...], None) veya (None, sebep).
    uc: uretilecek oligonun zinciri; _iupac_denetle'ye aynen gecer."""
    k, why = _iupac_denetle(win, a, uc)
    if k is None:
        return None, why
    if not k:
        return [(win, 0)], None
    secenekler = [IUPAC_KUME[win[i]] for i in k]
    cikti = []
    for kombin in itertools.product(*secenekler):
        L = list(win)
        for i, b in zip(k, kombin):
            L[i] = b
        cikti.append(("".join(L), len(k)))
    return cikti, None


def composition_ok(s, a):
    """Toplanti kararindaki oligo kurallari. (uygun_mu, sebep) doner."""
    if "N" in s:
        return False, "kalipta_N"
    degen = [c for c in s if c not in "ACGT"]
    if len(degen) > a.degeneracy_budget:
        return False, "dejenere_baz"
    if degen:
        fold = 1
        for c in degen:
            fold *= DEGEN_FOLD.get(c, 4)
        if fold > a.degeneracy_fold_max:
            return False, "dejenere_kat"
        # dejenere baz son bes bazda olamaz, uzama oradan baslar
        nlast = a.gc_clamp_last if a.gc_clamp_last > 0 else 5
        if any(c not in "ACGT" for c in s[-nlast:]):
            return False, "3p_dejenere"
    g = gc_pct(s)
    if not (a.gc_hard_min <= g <= a.gc_hard_max):
        return False, "gc_sert_sinir"
    if a.require_3p_gc and s[-1] not in "GC":
        return False, "3p_gc_degil"
    if a.gc_clamp_last > 0:
        tail = s[-a.gc_clamp_last:]
        if tail.count("G") + tail.count("C") > a.gc_clamp_max:
            return False, "3p_asiri_sabit"
    if max_run(s) > a.homopolymer_max:
        return False, "homopolimer"
    return True, ""


def tm_primer3(s, a):
    return primer3.calc_tm(s, mv_conc=a.mv, dv_conc=a.dv,
                           dntp_conc=a.dntp, dna_conc=a.dna_conc)


def tm_biopython(s, a):
    # SantaLucia 1998 komsu-cift tablosu, Owczarzy 2008 tuz duzeltmesi.
    # Iki kutuphanenin parametre eslemesi birebir ayni olmak zorunda degil;
    # kural sistematik kaymayi olcup ondan sapani elemek uzerine kurulu.
    return float(mt.Tm_NN(Seq(s), nn_table=mt.DNA_NN3, Na=a.mv, Mg=a.dv,
                          dNTPs=a.dntp, dnac1=a.dna_conc, dnac2=0,
                          saltcorr=7))


# ------------------------------------------------------------------------ ana
def main():
    a = get_args()
    label = a.label or os.path.basename(a.consensus).split("_consensus")[0]
    name, seq = read_fasta(a.consensus)
    L = len(seq)
    mask_contig = a.mask_contig if a.mask_contig else None
    masked, mclasses, mseen = read_mask(a.mask, L, contig=mask_contig)
    if a.mask and not masked:
        print("UYARI: maske dosyasi okundu ama hicbir pozisyon yasaklanmadi. "
              "Dosyadaki kontig adlari: %s" % (list(mseen) or "yok"))
    print("hedef            : %s" % label)
    print("konsensus        : %s (%d bp, basliği %s)" % (a.consensus, L, name))
    print("maske            : %s" % (a.mask or "yok"))
    if mclasses:
        print("maske siniflari  : %s" % ", ".join("%s=%d" % kv for kv in sorted(mclasses.items())))
    if L < a.len_min * 2 + a.prod_min:
        sys.exit("HATA: konsensus cok kisa (%d bp), en az %d bp gerekiyor"
                 % (L, a.len_min * 2 + a.prod_min))
    print("yasak pozisyon   : %d (%.2f%%)" % (len(masked), 100.0 * len(masked) / L))
    print("dejenere butcesi : %d pozisyon, en fazla %d kat"
          % (a.degeneracy_budget, a.degeneracy_fold_max))

    # --- 1. kompozisyon suzgeci, ucuz olan once ------------------------
    reasons = {}
    raw = []          # (yon, baslangic0, uzunluk, oligo_dizisi, cozulen_iupac)
    for start in range(L):
        for ln in range(a.len_min, a.len_max + 1):
            end = start + ln
            if end > L:
                break
            if any(i in masked for i in range(start, end)):
                reasons["maskeli"] = reasons.get("maskeli", 0) + 2  # F ve R
                continue
            win = seq[start:end]
            for strand in ("F", "R"):
                varyant, why = iupac_varyantlar(win, a, uc=strand)
                if varyant is None:
                    reasons[why] = reasons.get(why, 0) + 1
                    continue
                for coz, kac in varyant:
                    oligo = coz if strand == "F" else rc(coz)
                    ok, why = composition_ok(oligo, a)
                    if not ok:
                        reasons[why] = reasons.get(why, 0) + 1
                        continue
                    raw.append((strand, start, ln, oligo, kac))
    print("\nkompozisyon suzgeci sonrasi oligo: %d" % len(raw))
    for k, v in sorted(reasons.items(), key=lambda x: -x[1]):
        print("   elenen %-16s %d" % (k, v))
    if not raw:
        sys.exit("kompozisyon suzgecinden gecen oligo yok")

    # --- 2. iki bagimsiz Tm olcumu ve sistematik kayma ------------------
    tm3, tmb = [], []
    for _, _, _, oligo, _ in raw:
        tm3.append(tm_primer3(oligo, a))
        tmb.append(tm_biopython(oligo, a))
    diffs = [x - y for x, y in zip(tm3, tmb)]
    offset = statistics.median(diffs)
    spread = statistics.pstdev(diffs) if len(diffs) > 1 else 0.0
    if a.tm_cross_tol is not None:
        tol = a.tm_cross_tol
        tol_src = "elle verildi"
    else:
        tol = max(0.10, a.tm_cross_k * spread)
        tol_src = "veriden turetildi: %.1f carpi sd (%.3f), taban 0,10" % (
            a.tm_cross_k, spread)
    print("\niki kutuphane Tm karsilastirmasi (%d oligo)" % len(raw))
    print("   primer3 eksi Biopython, medyan kayma : %+.2f C" % offset)
    print("   kaymanin standart sapmasi            : %.3f C" % spread)
    print("   kullanilan tolerans                  : %.3f C (%s)" % (tol, tol_src))

    kept = []
    n_cross = n_tm = n_hp = n_hd = 0
    for (strand, start, ln, oligo, kac_iupac), t3, tb in zip(raw, tm3, tmb):
        if abs((t3 - tb) - offset) > tol:
            n_cross += 1
            continue
        tm = t3
        if not (a.tm_hard_min <= tm <= a.tm_hard_max):
            n_tm += 1
            continue
        hp = primer3.calc_hairpin(oligo, mv_conc=a.mv, dv_conc=a.dv,
                                  dntp_conc=a.dntp, dna_conc=a.dna_conc).dg
        if hp < a.hairpin_dg_min:
            n_hp += 1
            continue
        hd = primer3.calc_homodimer(oligo, mv_conc=a.mv, dv_conc=a.dv,
                                    dntp_conc=a.dntp, dna_conc=a.dna_conc).dg
        if hd < a.homodimer_dg_min:
            n_hd += 1
            continue
        kept.append(dict(strand=strand, start=start, ln=ln, oligo=oligo,
                         tm3=t3, tmb=tb, hairpin_dg=hp, homodimer_dg=hd,
                         gc=gc_pct(oligo), iupac=kac_iupac))
    print("   elenen, iki olcum ayrildi            : %d" % n_cross)
    print("   elenen, Tm sert sinir disi           : %d" % n_tm)
    print("   elenen, hairpin dG                   : %d" % n_hp)
    print("   elenen, self-dimer dG                : %d" % n_hd)
    print("termodinamik suzgec sonrasi oligo       : %d" % len(kept))
    if not kept:
        sys.exit("termodinamik suzgecten gecen oligo yok")

    F = sorted([k for k in kept if k["strand"] == "F"], key=lambda x: x["start"])
    R = sorted([k for k in kept if k["strand"] == "R"], key=lambda x: x["start"])
    print("   ileri aday: %d   geri aday: %d" % (len(F), len(R)))

    # --- 3. ciftleme, urun kurallari ve makine dogrulamasi -------------
    pairs = []
    n_prod = n_tmdiff = n_het = n_overlap = 0
    n_verify_fail = 0
    Rby = {}
    for r in R:
        Rby.setdefault(r["start"], []).append(r)
    rstarts = sorted(Rby)
    import bisect
    for f in F:
        lo = f["start"] + a.prod_min - a.len_max
        hi = f["start"] + a.prod_hard_max
        i = bisect.bisect_left(rstarts, lo)
        while i < len(rstarts) and rstarts[i] <= hi:
            for r in Rby[rstarts[i]]:
                pstart, pend = f["start"], r["start"] + r["ln"]
                plen = pend - pstart
                if plen < a.prod_min or plen > a.prod_hard_max:
                    n_prod += 1
                    continue
                if abs(f["tm3"] - r["tm3"]) >= a.pair_tm_diff_max:
                    n_tmdiff += 1
                    continue
                # F ve R ayak izleri cakisamaz; cakisirsa cogaltilabilir bir
                # urun degildir. Varsayilan esiklerde gizli kalir, ancak
                # --prod-min dusurulunce ortaya cikar.
                if f["start"] + f["ln"] > r["start"]:
                    n_overlap += 1
                    continue
                product = seq[pstart:pend]
                # makine dogrulamasi: urunun basi ileri primere, sonu geri
                # primerin ters tumleyenine birebir esit olmali
                if product[:f["ln"]] != f["oligo"] or product[-r["ln"]:] != rc(r["oligo"]):
                    n_verify_fail += 1
                    continue
                het = primer3.calc_heterodimer(f["oligo"], r["oligo"], mv_conc=a.mv,
                                               dv_conc=a.dv, dntp_conc=a.dntp,
                                               dna_conc=a.dna_conc).dg
                if het < a.heterodimer_dg_min:
                    n_het += 1
                    continue
                pgc = gc_pct(product)
                # puanlama: kucuk ceza toplami daha iyi
                pen = 0.0
                pen += 0.0 if a.prod_best_min <= plen <= a.prod_best_max else \
                    min(abs(plen - a.prod_best_min), abs(plen - a.prod_best_max)) * 0.05
                # --prod-max yumusak ust sinir: 250 ile 300 arasi kabul edilir
                # ama cezalandirilir, boylece parametre olu kalmaz
                if plen > a.prod_max:
                    pen += (plen - a.prod_max) * 0.10
                for t in (f["tm3"], r["tm3"]):
                    if not (a.tm_min <= t <= a.tm_max):
                        pen += min(abs(t - a.tm_min), abs(t - a.tm_max)) * 2.0
                for g in (f["gc"], r["gc"]):
                    if not (a.gc_min <= g <= a.gc_max):
                        pen += min(abs(g - a.gc_min), abs(g - a.gc_max)) * 0.2
                if not (a.prod_gc_min <= pgc <= a.prod_gc_max):
                    pen += min(abs(pgc - a.prod_gc_min), abs(pgc - a.prod_gc_max)) * 0.2
                pen += abs(f["tm3"] - r["tm3"]) * 1.0
                pen += max(0.0, -het / 1000.0) * 0.5
                pen += max(0.0, -(f["hairpin_dg"] + r["hairpin_dg"]) / 1000.0) * 0.3
                pairs.append(dict(
                    hedef=label,
                    ileri_dizi=f["oligo"], ileri_baslangic=f["start"] + 1,
                    ileri_uzunluk=f["ln"], ileri_tm_primer3=round(f["tm3"], 2),
                    ileri_tm_biopython=round(f["tmb"], 2), ileri_gc=round(f["gc"], 1),
                    ileri_hairpin_dg=round(f["hairpin_dg"], 1),
                    ileri_selfdimer_dg=round(f["homodimer_dg"], 1),
                    geri_dizi=r["oligo"], geri_baslangic=r["start"] + 1,
                    geri_uzunluk=r["ln"], geri_tm_primer3=round(r["tm3"], 2),
                    geri_tm_biopython=round(r["tmb"], 2), geri_gc=round(r["gc"], 1),
                    geri_hairpin_dg=round(r["hairpin_dg"], 1),
                    geri_selfdimer_dg=round(r["homodimer_dg"], 1),
                    cift_heterodimer_dg=round(het, 1),
                    tm_farki=round(abs(f["tm3"] - r["tm3"]), 2),
                    urun_uzunluk=plen, urun_gc=round(pgc, 1),
                    urun_baslangic=pstart + 1, urun_bitis=pend,
                    urun_dogrulandi="evet_ayni_kalip",
                    ceza=round(pen, 3)))
            i += 1

    print("\ncift olusturma")
    print("   elenen, urun uzunlugu                : %d" % n_prod)
    print("   elenen, cift Tm farki                : %d" % n_tmdiff)
    print("   elenen, F ve R ayak izi cakismasi    : %d" % n_overlap)
    print("   elenen, urun makine dogrulamasi      : %d" % n_verify_fail)
    print("     (not: urun bu asamada primerlerin turetildigi kalibin ayni")
    print("      kopyasindan kesildigi icin bu kontrol tanim geregi gecer;")
    print("      asil islevi 05 ozgulluk asamasinda rakip ve ham okuma")
    print("      kaliplarina karsi calistirildiginda ortaya cikar)")
    print("   elenen, hetero-dimer dG              : %d" % n_het)
    print("gecerli cift sayisi                     : %d" % len(pairs))
    if not pairs:
        sys.exit("gecerli cift bulunamadi")

    pairs.sort(key=lambda x: x["ceza"])
    if a.min_locus_spacing > 0:
        kept_pairs, taken = [], []
        for p in pairs:
            if all(abs(p["ileri_baslangic"] - q[0]) >= a.min_locus_spacing or
                   abs(p["geri_baslangic"] - q[1]) >= a.min_locus_spacing
                   for q in taken):
                kept_pairs.append(p)
                taken.append((p["ileri_baslangic"], p["geri_baslangic"]))
        print("lokus seyreltmesi (%d bp): %d -> %d cift"
              % (a.min_locus_spacing, len(pairs), len(kept_pairs)))
        pairs = kept_pairs
    if len(pairs) > a.max_pairs:
        print("cikti %d satira kirpildi (puan sirasina gore)" % a.max_pairs)
        pairs = pairs[:a.max_pairs]

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pairs[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(pairs)
    print("\nyazildi: %s" % a.out)
    print("\nEn iyi bes aday:")
    for p in pairs[:5]:
        print("  ceza=%.2f  urun=%d bp (GC %.1f)  F=%s (Tm %.1f)  R=%s (Tm %.1f)"
              % (p["ceza"], p["urun_uzunluk"], p["urun_gc"], p["ileri_dizi"],
                 p["ileri_tm_primer3"], p["geri_dizi"], p["geri_tm_primer3"]))
    print("\nSonraki asama: ozgulluk suzgeci. Bu adaylar rakip konsensuslere ve")
    print("REFERANS_DB veritabanlarina karsi mfeprimer ve blastn ile sinanacak.")


if __name__ == "__main__":
    main()
