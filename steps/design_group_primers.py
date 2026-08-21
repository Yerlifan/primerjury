#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
design_group_primers.py
Çok üyeli hedef kümeleri için primer çifti üretir. Karar 2 (cins özgül),
Karar 3 (işlev grubu) ve Karar 4 (universal) aynı motoru kullanır; aradaki
tek fark hedef kümesinin büyüklüğü ve dejenerelik bütçesidir.

Yöntem: hizalama kullanılmaz. Adaylar bir çapa konsensüsten üretilir, sonra
her aday toplantı kararındaki BAĞLANMA KURALIYLA her üyeye ve her rakibe
karşı taranır:
    son iki baz hedefe birebir uymalı (uzama oradan başlar)
    son beş bazda en fazla bir uyumsuzluk
    primerin tamamında en fazla üç uyumsuzluk
    5' tarafta sarkma serbest
    iki primer ters zincirlerde ve 3' uçları birbirine bakacak

Kabul ölçütleri:
    hedeflenen HER üyede ürün oluşmalı, ürün uzunluğu aralıkta kalmalı
    rakiplerin HİÇBİRİNDE ürün oluşmamalı
    ayrım sağlam olmalı: primerlerden en az biri rakiplerde hiç bağlanma
    yeri bulamamalı. İki primerin de zayıf bağlanıp yalnızca birlikte
    yetersiz kalmasıyla oluşan temizlik kabul edilmez.

KONSENSÜS HANGİ KLASÖRDEN OKUNMALI (2026-08-21 düzeltmesi)
    Yalnız `konsensus_kanonik/` kullanın. Eski örnekler `consensus sequences/`
    klasörünü gösteriyordu; o klasör KARIŞIK YÖNLÜDÜR (ölçülen: 71 antisense /
    27 sense). Ters yönlü bir konsensüste in-silico PCR SESSİZCE 0 ürün verir —
    ölçülen kayıp %100, kanıt `screening/orientation_impact_test.py`.
    Bu betikte tehlike daha büyüktür: girdi bir GLOB'dur, yani ters yönlü tek
    bir üye sessizce "bu üyede ürün yok" diye sayılır ve çift haksız yere
    elenir. Kutu -> dosya eşlemesi `konsensus_kanonik/INDEKS.tsv` içindedir.

Kullanım:
  python3 design_group_primers.py \
     --in-group  "konsensus_kanonik/*_2209.kanonik.fa" \
                 "konsensus_kanonik/*_2223.kanonik.fa" \
     --out-group "konsensus_kanonik/*_394967.kanonik.fa" \
     --label Asetoklastik_metanojenler \
     --out primer_adaylari/Asetoklastik.tsv \
     --degeneracy-budget 2
"""
import argparse, csv, glob, importlib.util, os, re, statistics, sys, bisect

HERE = os.path.dirname(os.path.abspath(__file__))


def load_engine():
    """03 betigindeki kural fonksiyonlarini tek kaynak olarak kullanir."""
    p = os.path.join(HERE, "generate_primer_candidates.py")
    if not os.path.exists(p):
        sys.exit("generate_primer_candidates.py ayni klasorde bulunamadi: %s" % HERE)
    spec = importlib.util.spec_from_file_location("engine03", p)
    m = importlib.util.module_from_spec(spec)
    sys.argv_backup, sys.argv = sys.argv, [p]
    try:
        spec.loader.exec_module(m)
    finally:
        sys.argv = sys.argv_backup
    return m


E = load_engine()
rc, gc_pct, composition_ok = E.rc, E.gc_pct, E.composition_ok
read_fasta, read_mask = E.read_fasta, E.read_mask
tm_primer3, tm_biopython = E.tm_primer3, E.tm_biopython
import primer3

IUPAC_SET = {"A": "A", "C": "C", "G": "G", "T": "T",
             "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
             "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def base_match(p, t):
    """Primer bazi p, kalip bazi t ile uyusuyor mu. IUPAC kumesi kesisimi."""
    if t == "N":
        return False          # kalipta N varsa uyum sayilmaz, bilgi yok
    return bool(set(IUPAC_SET.get(p, "")) & set(IUPAC_SET.get(t, "")))


def build_index(seq, k, azami_acilim=64):
    """k-mer -> baslangic pozisyonlari listesi.

    Kalipta IUPAC kodu bulunabilir (konsensus -A ile uretildiginde) ve
    seed_variants yalnizca A, C, G, T, N uretir. Indeks ham k-mer'i anahtar
    yapsaydi, IUPAC kodu iceren pencereler hicbir cekirdek varyantiyla
    eslesmez ve o baglanma yerleri TAMAMEN gozden kacardi. Olculdu: 12
    gercek konsensus uzerinde find_bindings kuralin dogrudan uygulanmasina
    gore 20 baglanma yerini kaciriyordu ve bunlardan bazilari 'rakipte urun
    yok' kararini tersine cevirecek yerlerdi. Bu yuzden IUPAC iceren her
    k-mer, temsil ettigi somut ACGT k-mer'lerinin hepsine kaydedilir;
    boylece indeks anahtarlari salt ACGT (ve N) olur.

    azami_acilim: bir k-mer'in acilim sayisi bunu asarsa acilim yapilmaz ve
    ham k-mer anahtar olarak birakilir; log'a dusmez ama boyle bir pencere
    zaten bes bazin ucunden fazlasi belirsiz demektir ve primer tasarimina
    girmez."""
    idx = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if all(c in "ACGTN" for c in kmer):
            idx.setdefault(kmer, []).append(i)
            continue
        acilim = [""]
        tasti = False
        for c in kmer:
            secenek = IUPAC_SET.get(c, "ACGT") if c != "N" else "N"
            acilim = [p + o for p in acilim for o in secenek]
            if len(acilim) > azami_acilim:
                tasti = True
                break
        if tasti:
            idx.setdefault(kmer, []).append(i)
            continue
        for v in acilim:
            idx.setdefault(v, []).append(i)
    return idx


def seed_variants(tail, max_mm_in_tail, exact_last):
    """Son 'tail' icin izin verilen KALIP k-mer varyantlari.
    Son 'exact_last' baz sabit, geri kalanda en fazla max_mm_in_tail
    uyumsuzluk. Kalipta N bulunabilir; base_match N'i uyumsuzluk sayar,
    dolayisiyla serbest konumlarda N de bir varyanttir. N sadece serbest
    konumlarda uretilir, cunku son iki bazin birebir uymasi gerekir ve
    N ile birebir uyum mumkun degildir."""
    n = len(tail)
    free = n - exact_last
    out = set()

    def rec(i, cur, mm):
        if i == n:
            out.add("".join(cur))
            return
        allowed = IUPAC_SET.get(tail[i], "ACGT")
        fixed = i >= free
        for b in "ACGTN":
            if b == "N":
                if fixed:
                    continue          # son iki bazda N kabul edilemez
                hit = False           # kalipta N her zaman uyumsuzluk
            else:
                hit = b in allowed
                if fixed and not hit:
                    continue
            nm = mm + (0 if hit else 1)
            if nm > max_mm_in_tail:
                continue
            rec(i + 1, cur + [b], nm)

    rec(0, [], 0)
    return out


def find_bindings(oligo, seq, idx, k, a):
    """Oligo'nun seq uzerindeki baglanma yerleri.
    Toplanti karari geregi 5' tarafta sarkma serbesttir: oligonun 5' ucu
    kalibin disina tasabilir ve tasan kisim uyumsuzluk sayilmaz. 3' uc
    kalibin icinde olmak zorundadir, cunku uzama oradan baslar.
    Doner: [(3'_uc_pozisyonu, ortusen_bolgedeki_uyumsuzluk)] , 0 tabanli."""
    L, n = len(seq), len(oligo)
    tail = oligo[-k:]
    hits = []
    seen = set()
    for var in seed_variants(tail, a.tail_max_mm, a.exact_last):
        for pos in idx.get(var, ()):          # pos = tail baslangici
            start = pos - (n - k)
            end = start + n
            if end > L:
                continue                      # 3' uc kalip disinda kalamaz
            j0 = max(0, -start)               # 5' sarkma kadar atla
            if n - j0 < a.min_overlap:
                continue
            mm = 0
            ok = True
            for j in range(j0, n):
                if not base_match(oligo[j], seq[start + j]):
                    mm += 1
                    if mm > a.total_max_mm:
                        ok = False
                        break
            if ok:
                key = (end - 1, mm)
                if key not in seen:
                    seen.add(key)
                    hits.append(key)
    return hits


def load_set(patterns, mask_dir=None):
    """Desenlerden konsensus kumesi yukler. Doner: [(etiket, dizi)]"""
    out, seen = [], set()
    for pat in patterns:
        for p in sorted(glob.glob(pat)):
            if p in seen:
                continue
            seen.add(p)
            name, s = read_fasta(p)
            tag = os.path.basename(p).split("_consensus")[0]
            out.append((tag, s, p))
    return out


def get_args():
    p = argparse.ArgumentParser(description="Cok uyeli hedef icin primer cifti")
    p.add_argument("--in-group", nargs="+", required=True,
                   help="hedeflenen uyelerin consensus files, glob kabul eder")
    p.add_argument("--out-group", nargs="*", default=[],
                   help="competitor consensus files, glob kabul eder")
    p.add_argument("--anchor", default=None,
                   help="consensus the candidates are generated from; if omitted, the one with fewest Ns")
    p.add_argument("--mask-dir", default=None, help="02 betiginin maske directory")
    p.add_argument("--label", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--max-pairs", type=int, default=2000)
    p.add_argument("--max-oligo", type=int, default=400,
                   help="how many to keep per strand after the thermodynamic filter, at most "
                        "fazla oligo. Tek takson hedeflerinde binlerce korunmus "
                        "oligo cikiyor ve F x R carpimi milyonlara ulasiyor. "
                        "Secim pozisyona gore tabakalanir, yani kalip boyunca "
                        "esit dagilir; 0 sinirsiz.")
    p.add_argument("--stop-after", type=int, default=20000,
                   help="stop pairing after this many valid pairs; 0 means unlimited")
    # baglanma kurali
    p.add_argument("--exact-last", type=int, default=2)
    p.add_argument("--tail-len", type=int, default=5)
    p.add_argument("--tail-max-mm", type=int, default=1)
    p.add_argument("--total-max-mm", type=int, default=3)
    p.add_argument("--min-overlap", type=int, default=15,
                   help="5' overhang is allowed, but the part overlapping the template must be at least "
                        "bu kadar baz olmali")
    p.add_argument("--competitor-prod-max", type=int, default=0,
                   help="rakipte urun aranirken kullanilan ust sinir. 0 ise "
                        "sinirsiz, yani dizi boyunca olusan her bant sayilir. "
                        "Toplanti karari 'rakiplerin hicbirinde urun olusmamali' "
                        "dedigi icin varsayilan sinirsizdir.")
    # ozgulluk
    p.add_argument("--rakip-prod-min", type=int, default=1,
                   help="rakipte urun sayilmasi for en kucuk boy; rakipte "
                        "herhangi bir bant istenmedigi icin varsayilan 1")
    p.add_argument("--yetim-min-uyumsuzluk", type=int, default=0,
                   help="0: kati kural, yetim primer rakiplerde HIC baglanmamali. "
                        ">0: rakiplerdeki en iyi yerlesimi bu kadar uyumsuzluk "
                        "tasiyan primer de yetim sayilir (gevsetilmis kademe)")
    p.add_argument("--require-orphan-primer", type=int, default=1,
                   help="1 ise primerlerden biri rakiplerde HIC baglanma yeri "
                        "bulamamali (toplanti karari)")
    # 03 ile ayni oligo, termodinamik ve urun kurallari
    p.add_argument("--len-min", type=int, default=18)
    p.add_argument("--len-max", type=int, default=25)
    p.add_argument("--gc-min", type=float, default=40.0)
    p.add_argument("--gc-max", type=float, default=60.0)
    p.add_argument("--gc-hard-min", type=float, default=35.0)
    p.add_argument("--gc-hard-max", type=float, default=65.0)
    p.add_argument("--gc-clamp-last", type=int, default=5)
    p.add_argument("--gc-clamp-max", type=int, default=3)
    p.add_argument("--homopolymer-max", type=int, default=4)
    p.add_argument("--require-3p-gc", type=int, default=1)
    p.add_argument("--degeneracy-budget", type=int, default=0)
    p.add_argument("--degeneracy-fold-max", type=int, default=4)
    p.add_argument("--varyantlari-tut", action="store_true",
                   help="keep IUPAC sibling variants of the same locus as separate rows")
    p.add_argument("--iupac-max", type=int, default=2)
    p.add_argument("--iupac-son-yasak", type=int, default=5)
    p.add_argument("--tm-min", type=float, default=58.0)
    p.add_argument("--tm-max", type=float, default=62.0)
    p.add_argument("--tm-hard-min", type=float, default=57.0)
    p.add_argument("--tm-hard-max", type=float, default=63.0)
    p.add_argument("--tm-cross-tol", type=float, default=2.0)
    p.add_argument("--hairpin-dg-min", type=float, default=-3000.0)
    p.add_argument("--homodimer-dg-min", type=float, default=-6000.0)
    p.add_argument("--heterodimer-dg-min", type=float, default=-6000.0)
    p.add_argument("--mv", type=float, default=50.0)
    p.add_argument("--dv", type=float, default=1.5)
    p.add_argument("--dntp", type=float, default=0.6)
    p.add_argument("--dna-conc", type=float, default=50.0)
    p.add_argument("--prod-min", type=int, default=70)
    p.add_argument("--prod-max", type=int, default=250)
    p.add_argument("--prod-hard-max", type=int, default=300)
    p.add_argument("--prod-best-min", type=int, default=90)
    p.add_argument("--prod-best-max", type=int, default=150)
    p.add_argument("--pair-tm-diff-max", type=float, default=1.5)
    return p.parse_args()


def main():
    a = get_args()
    ing = load_set(a.in_group)
    outg = load_set(a.out_group) if a.out_group else []
    if not ing:
        sys.exit("hedef kumesi bos, --in-group desenlerini kontrol edin")
    # Etiket cakismasi: seqs sozlugu tag ile anahtarlandigi icin ayni etikete
    # dusen iki uye sessizce tek diziye cokerdi ve cikti yine iki uyeyi
    # kapsadigini iddia ederdi.
    seen_tags = {}
    for tag, _, p_ in ing + outg:
        seen_tags.setdefault(tag, []).append(p_)
    dup = {t: v for t, v in seen_tags.items() if len(v) > 1}
    if dup:
        for t, v in dup.items():
            print("ETIKET CAKISMASI: %s -> %s" % (t, v), file=sys.stderr)
        sys.exit("Ayni etikete dusen dosyalar var. Sessiz uye kaybini onlemek "
                 "icin duruldu; glob desenlerini ayirin.")
    print("etiket            : %s" % a.label)
    print("hedef uye sayisi  : %d" % len(ing))
    for t, s, p in ing:
        print("    hedef  %-28s %5d bp  N=%d" % (t, len(s), s.count("N")))
    print("rakip sayisi      : %d" % len(outg))
    for t, s, p in outg:
        print("    rakip  %-28s %5d bp  N=%d" % (t, len(s), s.count("N")))
    print("baglanma kurali   : son %d baz birebir, son %d bazda en fazla %d "
          "uyumsuzluk, toplamda en fazla %d"
          % (a.exact_last, a.tail_len, a.tail_max_mm, a.total_max_mm))
    print("dejenere butcesi  : %d pozisyon, en fazla %d kat"
          % (a.degeneracy_budget, a.degeneracy_fold_max))

    # --- capa secimi -------------------------------------------------
    if a.anchor:
        anchor = [x for x in ing if a.anchor in x[0] or a.anchor in x[2]]
        if not anchor:
            sys.exit("capa bulunamadi: %s" % a.anchor)
        anchor = anchor[0]
    else:
        anchor = min(ing, key=lambda x: (x[1].count("N"), -len(x[1])))
    print("capa konsensus    : %s" % anchor[0])

    # --- yon normalizasyonu -------------------------------------------
    # Konsensusler yon normalizasyonu yapilmadan uretilmis: consensus2.sh her
    # takson icin bir cekirdek okuma seciyor ve konsensusu o okumanin yonunde
    # kuruyor, dolayisiyla uyelerin bir kismi ters saklanmis olabilir. Ters bir
    # uye, ileri primerin o uyede arti zincir yerine eksi zincire baglanmasina
    # yol acar ve o uyede urun hic olusmaz. Capadan alinan korunmus problarla
    # her dizinin yonu oylanir, ters olanlar ters tumleyene cevrilir.
    K0 = a.tail_len
    probes = [anchor[1][i:i + 20] for i in range(0, len(anchor[1]) - 20, 40)
              if "N" not in anchor[1][i:i + 20]][:40]
    if not probes:
        sys.exit("capada prob uretilemedi, konsensus fazla N iceriyor")

    # Ikinci, bagimsiz yon olcutu: butun canlilarda korunmus SSU motifleri.
    # Capa problari uzak bir uyeye hic baglanmadiginda (arti=eksi=0) oylama
    # karar veremez; o durumda bu motifler devreye girer. Ikisi de karar
    # veremezse uye sessizce gecirilmez, ADI YAZILARAK bildirilir.
    UNIV = ["GTGCCAGCMGCCGCGGTAA", "GGATTAGATACCC", "AAACTCAAAGGAATTGACGG",
            "GTGYCAGCMGCCGCGGTAA", "ATTAGATACCCBDGTAGTCC"]

    def _expand(m):
        alt = [""]
        for ch in m:
            opts = IUPAC_SET.get(ch, "ACGT")
            alt = [p + o for p in alt for o in opts]
            if len(alt) > 64:
                return alt[:64]
        return alt

    UNIV_EXP = []
    for m in UNIV:
        UNIV_EXP.extend(_expand(m))

    def univ_vote(seq_):
        r = rc(seq_)
        f = sum(1 for m in UNIV_EXP if m in seq_)
        b = sum(1 for m in UNIV_EXP if m in r)
        return f, b

    def orient(tag, seq_):
        ip = build_index(seq_, K0)
        r = rc(seq_)
        im = build_index(r, K0)
        np_ = sum(1 for o in probes if find_bindings(o, seq_, ip, K0, a))
        nm = sum(1 for o in probes if find_bindings(o, r, im, K0, a))
        if np_ == 0 and nm == 0:
            uf, ub = univ_vote(seq_)
            if uf == 0 and ub == 0:
                return (seq_, np_, nm, False, "KARARSIZ")
            return ((seq_, np_, nm, False, "motif") if uf >= ub
                    else (r, np_, nm, True, "motif"))
        return ((seq_, np_, nm, False, "prob") if np_ >= nm
                else (r, np_, nm, True, "prob"))

    flipped, undecided = [], []
    ing2, outg2 = [], []
    for grp, src, dst in (("hedef", ing, ing2), ("rakip", outg, outg2)):
        for tag, s_, p_ in src:
            s2, np_, nm, fl, how = orient(tag, s_)
            dst.append((tag, s2, p_))
            if fl:
                flipped.append((tag, np_, nm, how))
            if how == "KARARSIZ":
                undecided.append((grp, tag, np_, nm))
    ing, outg = ing2, outg2
    print("\nyon normalizasyonu: %d capa probu artı %d universal motif varyanti"
          % (len(probes), len(UNIV_EXP)))
    if flipped:
        for t, np_, nm, how in flipped:
            print("   TERS bulundu, cevrildi: %-30s arti=%d eksi=%d (%s)"
                  % (t, np_, nm, how))
    else:
        print("   capayla ayni yonde olmayan dizi yok")
    if undecided:
        print("   UYARI: yonu belirlenemeyen %d dizi var, ikisi de sifir oy aldi:"
              % len(undecided))
        for grp, t, np_, nm in undecided:
            print("      %-6s %-30s arti=%d eksi=%d" % (grp, t, np_, nm))
        print("   Bu diziler capadan cok uzak demektir. Hedef kumesindelerse")
        print("   urun vermeyecek, rakip kumesindelerse ozgulluk denetimi")
        print("   guvenilmez olur. Kumeleri gozden gecirin.")

    masked = set()
    if a.mask_dir:
        # Etiket "A1-1-reads_2209" ya da "A1_1_reads_2223" bicimindedir; grup ve
        # taxid ikisi birden kullanilir, boylece baska taksonlarin koordinatlari
        # capaya bindirilmez.
        m_ = re.search(r"reads[-_](\d+)", anchor[0])
        tid = m_.group(1) if m_ else None
        grp = re.split(r"[-_]reads", anchor[0])[0].replace("_", "-")
        pats = []
        if tid:
            pats = [os.path.join(a.mask_dir, "%s_%s_maske.bed" % (grp, tid)),
                    os.path.join(a.mask_dir, "*%s*%s*maske.bed" % (grp, tid))]
        hits = []
        for pat in pats:
            hits.extend(glob.glob(pat))
        hits = sorted(set(hits))
        if not hits:
            sys.exit("HATA: --mask-dir verildi ama capa '%s' (grup=%s taxid=%s) "
                     "icin maske dosyasi bulunamadi. Sessiz maskesizligi "
                     "onlemek icin duruldu." % (anchor[0], grp, tid))
        if len(hits) > 1:
            print("UYARI: capa icin birden fazla maske dosyasi eslesti: %s" % hits)
        for cand in hits:
            m, _c, _s = read_mask(cand, len(anchor[1]))
            masked |= m
        print("capada yasak poz. : %d  (maske: %s)" % (len(masked), hits[0]))
    else:
        print("capada yasak poz. : 0  (maske verilmedi)")

    # --- 1. capadan aday uretimi, kompozisyon suzgeci -----------------
    seq = anchor[1]
    L = len(seq)
    raw, reasons = [], {}
    for start in range(L):
        for ln in range(a.len_min, a.len_max + 1):
            end = start + ln
            if end > L:
                break
            if any(i in masked for i in range(start, end)):
                reasons["maskeli"] = reasons.get("maskeli", 0) + 1
                continue
            win = seq[start:end]
            # Zincir basina ayri denetim: R'nin 3' ucu pencerenin BASIDIR.
            for strand in ("F", "R"):
                varyant, why = E.iupac_varyantlar(win, a, uc=strand)
                if varyant is None:
                    reasons[why] = reasons.get(why, 0) + 1
                    continue
                for coz, kac in varyant:
                    oligo = coz if strand == "F" else rc(coz)
                    ok, why = composition_ok(oligo, a)
                    if ok:
                        raw.append((strand, start, ln, oligo, kac))
                    else:
                        reasons[why] = reasons.get(why, 0) + 1
    print("\nkompozisyon sonrasi oligo: %d" % len(raw))
    for k, v in sorted(reasons.items(), key=lambda x: -x[1])[:6]:
        print("   elenen %-16s %d" % (k, v))
    if not raw:
        sys.exit("kompozisyon suzgecinden gecen oligo yok")

    # --- 2. iki bagimsiz Tm ve termodinamik ---------------------------
    tm3 = [tm_primer3(o, a) for _, _, _, o, _ in raw]
    tmb = [tm_biopython(o, a) for _, _, _, o, _ in raw]
    offset = statistics.median(x - y for x, y in zip(tm3, tmb))
    print("\nprimer3 eksi Biopython medyan kayma: %+.2f C (tolerans %.2f)"
          % (offset, a.tm_cross_tol))
    kept = []
    for (strand, start, ln, o, kac), t3, tb in zip(raw, tm3, tmb):
        if abs((t3 - tb) - offset) > a.tm_cross_tol:
            continue
        if not (a.tm_hard_min <= t3 <= a.tm_hard_max):
            continue
        hp = primer3.calc_hairpin(o, mv_conc=a.mv, dv_conc=a.dv,
                                  dntp_conc=a.dntp, dna_conc=a.dna_conc).dg
        if hp < a.hairpin_dg_min:
            continue
        hd = primer3.calc_homodimer(o, mv_conc=a.mv, dv_conc=a.dv,
                                    dntp_conc=a.dntp, dna_conc=a.dna_conc).dg
        if hd < a.homodimer_dg_min:
            continue
        kept.append(dict(strand=strand, start=start, ln=ln, oligo=o, tm=t3,
                         tmb=tb, hairpin=hp, homodimer=hd, gc=gc_pct(o),
                         iupac=kac))
    print("termodinamik sonrasi oligo: %d" % len(kept))
    if not kept:
        sys.exit("termodinamik suzgecten gecen oligo yok")

    # --- 3. her uye ve rakip icin baglanma taramasi -------------------
    K = a.tail_len
    seqs = {}
    for tag, s, _ in ing + outg:
        seqs[tag] = dict(plus=s, minus=rc(s),
                         idx_plus=build_index(s, K), idx_minus=build_index(rc(s), K))
    ing_tags = [t for t, _, _ in ing]
    out_tags = [t for t, _, _ in outg]

    print("\n%d oligo x %d dizi taraniyor" % (len(kept), len(seqs)))
    bind = {}     # oligo -> {tag: {"plus":[(3'poz,mm)], "minus":[...]}}
    for k in kept:
        o = k["oligo"]
        if o in bind:
            continue
        rec = {}
        for tag, d in seqs.items():
            rec[tag] = dict(
                L=len(d["plus"]),
                plus=find_bindings(o, d["plus"], d["idx_plus"], K, a),
                minus=find_bindings(o, d["minus"], d["idx_minus"], K, a))
        bind[o] = rec

    # her uyede en az bir baglanma yeri olan oligolar
    universal = [k for k in kept
                 if all(bind[k["oligo"]][t]["plus"] or bind[k["oligo"]][t]["minus"]
                        for t in ing_tags)]
    print("her hedef uyeye baglanan oligo: %d" % len(universal))
    if not universal:
        sys.exit("butun uyelere baglanan oligo yok. Dejenere butcesini artirmayi "
                 "(--degeneracy-budget) ya da hedef kumesini daraltmayi deneyin.")

    # rakiplerde hic baglanmayan oligolar (yetim primer adaylari)
    orphan = set()
    rakip_en_iyi = {}      # oligo -> rakiplerdeki EN IYI (en dusuk) uyumsuzluk
    if out_tags:
        n_tam = 0
        for k in universal:
            o = k["oligo"]
            mm = [m for t in out_tags
                  for _, m in (bind[o][t]["plus"] + bind[o][t]["minus"])]
            en_iyi = min(mm) if mm else None
            rakip_en_iyi[o] = en_iyi
            if en_iyi is None:
                orphan.add(o); n_tam += 1
            elif a.yetim_min_uyumsuzluk and en_iyi >= a.yetim_min_uyumsuzluk:
                # Katı kural sağlanamadığında kullanılan kademe: primerin
                # rakiplerdeki EN IYI yerleşimi bile bu kadar uyumsuzluk
                # taşıyorsa bağlanma tavlama sıcaklığı düştüğünde bile
                # zayıf kalır. Bu kademe çıktıda ayrıca işaretlenir.
                orphan.add(o)
        print("rakiplerde hic baglanmayan oligo: %d" % n_tam)
        if a.yetim_min_uyumsuzluk:
            print("rakiplerde en iyi yerlesimi >=%d uyumsuzluk olan oligo: %d "
                  "(gevsetilmis kademe)" % (a.yetim_min_uyumsuzluk,
                                            len(orphan) - n_tam))
        dag = {}
        for v in rakip_en_iyi.values():
            dag[v] = dag.get(v, 0) + 1
        print("   rakiplerdeki en iyi uyumsuzluk dagilimi: %s"
              % ", ".join("%s=%d" % ("hic" if k is None else k, v)
                          for k, v in sorted(dag.items(),
                                             key=lambda x: (x[0] is not None, x[0]))))
    else:
        print("rakip kumesi bos: yetim primer kurali uygulanamiyor, "
              "ozgulluk guvencesi VERILMEZ")
        if a.require_orphan_primer:
            print("   (--require-orphan-primer 1 ama rakip yok, kural atlaniyor)")

    # --- 4. ciftleme ve her uyede urun dogrulamasi --------------------
    Fs = [k for k in universal if k["strand"] == "F"]
    Rs = [k for k in universal if k["strand"] == "R"]
    print("ileri aday: %d   geri aday: %d" % (len(Fs), len(Rs)))

    def ozgulluk_skoru(k):
        """Oligonun rakiplerdeki EN IYI yerlesiminin uyumsuzluk sayisi.
        Rakiplerde hic baglanmiyorsa 99 sayilir. Buyuk olan daha ozguldur."""
        v = rakip_en_iyi.get(k["oligo"], None)
        return 99 if v is None else v

    def tabakala(lst, n):
        """Pozisyona gore tabakali secim: kalip n dilime bolunur ve her
        dilimden bir oligo alinir. Dilim icindeki secim once OZGULLUGE,
        sonra Tm'in 58-62 bandinin ortasina yakinligina gore yapilir.
        Sadece Tm'e bakmak, korunmus bolgelerde binlerce ayirt edici
        olmayan oligonun ayirt edici olanlari disari itmesine yol aciyordu;
        rDNA'da korunmus omurga ile degisken bolgelerin sayica orani
        buyuk oldugu icin bu secim olcutu sonucu belirliyor."""
        if not n or len(lst) <= n:
            return lst
        lst = sorted(lst, key=lambda k: k["start"])
        L = lst[-1]["start"] - lst[0]["start"] + 1
        mid = (a.tm_min + a.tm_max) / 2.0
        anahtar = lambda k: (-ozgulluk_skoru(k), abs(k["tm"] - mid))
        kova = {}
        for k in lst:
            b = int((k["start"] - lst[0]["start"]) * n / L)
            cur = kova.get(b)
            if cur is None or anahtar(k) < anahtar(cur):
                kova[b] = k
        out = list(kova.values())
        secili = set(id(k) for k in out)
        # kova sayisi n'den azsa kalanlari once ozgulluk sonra Tm ile doldur
        if len(out) < n:
            kalan = [k for k in lst if id(k) not in secili]
            kalan.sort(key=anahtar)
            out.extend(kalan[:n - len(out)])
        return out

    if a.max_oligo:
        f0, r0 = len(Fs), len(Rs)
        Fs, Rs = tabakala(Fs, a.max_oligo), tabakala(Rs, a.max_oligo)
        if len(Fs) < f0 or len(Rs) < r0:
            print("pozisyona gore tabakali secim (--max-oligo %d): ileri %d -> %d, "
                  "geri %d -> %d" % (a.max_oligo, f0, len(Fs), r0, len(Rs)))
    pairs = []
    fail_member = {}
    n_noprod = n_comp = n_orph = n_tmd = n_het = 0
    durdu = False
    for f in Fs:
        if durdu:
            break
        for r in Rs:
            if a.stop_after and len(pairs) >= a.stop_after:
                durdu = True
                break
            if abs(f["tm"] - r["tm"]) >= a.pair_tm_diff_max:
                n_tmd += 1
                continue
            # her hedef uyede urun var mi
            prods = {}
            ok = True
            for t in ing_tags:
                p = product_len(bind[f["oligo"]][t], bind[r["oligo"]][t],
                                f["ln"], r["ln"], a)
                if p is None:
                    ok = False
                    fail_member[t] = fail_member.get(t, 0) + 1
                    break
                prods[t] = p
            if not ok:
                n_noprod += 1
                continue
            # hicbir rakipte urun olmamali
            bad = False
            for t in out_tags:
                # rakipte ANY bant reddedilir, yalnizca hedef uzunluk penceresi
                # degil; 370 bp'lik bir bant da PCR'de olusur
                cpmax = a.competitor_prod_max if a.competitor_prod_max else 0
                # Rakipte HERHANGI bir bant istenmiyor. Alt siniri burada
                # uygulamak, 70 bp altindaki capraz bantlari "urun yok"
                # saymak demekti; jelde gorunen kisa bir bant da capraz
                # cogalmadir.
                if product_len(bind[f["oligo"]][t], bind[r["oligo"]][t],
                               f["ln"], r["ln"], a, pmax=cpmax,
                               pmin=a.rakip_prod_min) is not None:
                    bad = True
                    break
            if bad:
                n_comp += 1
                continue
            # ayrim saglamligi: en az biri rakiplerde hic baglanmamali
            if a.require_orphan_primer and out_tags:
                if f["oligo"] not in orphan and r["oligo"] not in orphan:
                    n_orph += 1
                    continue
            het = primer3.calc_heterodimer(f["oligo"], r["oligo"], mv_conc=a.mv,
                                           dv_conc=a.dv, dntp_conc=a.dntp,
                                           dna_conc=a.dna_conc).dg
            if het < a.heterodimer_dg_min:
                n_het += 1
                continue
            pl = list(prods.values())
            pen = abs(f["tm"] - r["tm"]) * 1.0
            pen += (max(pl) - min(pl)) * 0.02          # uyeler arasi urun tutarliligi
            for t in (f["tm"], r["tm"]):
                if not (a.tm_min <= t <= a.tm_max):
                    pen += min(abs(t - a.tm_min), abs(t - a.tm_max)) * 2.0
            avg = statistics.mean(pl)
            if not (a.prod_best_min <= avg <= a.prod_best_max):
                pen += min(abs(avg - a.prod_best_min), abs(avg - a.prod_best_max)) * 0.05
            # --prod-max yumusak ust sinir: asildiginda ceza, --prod-hard-max
            # ise mutlak sinir. Eski surumde --prod-max hic okunmuyordu.
            if max(pl) > a.prod_max:
                pen += (max(pl) - a.prod_max) * 0.05
            # Her uyede o uyedeki EN IYI baglanmanin uyumsuzlugu alinir
            # (primer orada baglanacaktir), sonra uyeler arasindaki EN KOTUSU
            # raporlanir: yani "en zayif uyede kac uyumsuzlukla bagliyor".
            mmF = max(min(m for _, m in (bind[f["oligo"]][t]["plus"] +
                                         bind[f["oligo"]][t]["minus"]))
                      for t in ing_tags)
            mmR = max(min(m for _, m in (bind[r["oligo"]][t]["plus"] +
                                         bind[r["oligo"]][t]["minus"]))
                      for t in ing_tags)
            pen += (mmF + mmR) * 0.5
            pairs.append(dict(
                hedef_grubu=a.label, capa=anchor[0],
                ileri_baslangic=f["start"] + 1, ileri_uzunluk=f["ln"],
                geri_baslangic=r["start"] + 1, geri_uzunluk=r["ln"],
                ileri_dizi=f["oligo"], ileri_tm=round(f["tm"], 2), ileri_gc=round(f["gc"], 1),
                ileri_iupac_cozulen=f.get("iupac", 0),
                geri_dizi=r["oligo"], geri_tm=round(r["tm"], 2), geri_gc=round(r["gc"], 1),
                geri_iupac_cozulen=r.get("iupac", 0),
                tm_farki=round(abs(f["tm"] - r["tm"]), 2),
                heterodimer_dg=round(het, 1),
                uye_sayisi=len(ing_tags), rakip_sayisi=len(out_tags),
                urun_min=min(pl), urun_maks=max(pl), urun_ortalama=round(avg, 1),
                en_kotu_uyumsuzluk_ileri=mmF, en_kotu_uyumsuzluk_geri=mmR,
                yetim_primer=("rakip_verilmedi" if not out_tags else
                              ("ileri" if f["oligo"] in orphan else
                               ("geri" if r["oligo"] in orphan else "yok"))),
                yetim_kademe=("rakip_verilmedi" if not out_tags else
                              ("kati" if (rakip_en_iyi.get(f["oligo"]) is None or
                                          rakip_en_iyi.get(r["oligo"]) is None)
                               else ("gevsetilmis" if (f["oligo"] in orphan or
                                                      r["oligo"] in orphan)
                                     else "yok"))),
                ileri_rakip_en_iyi_uyumsuzluk=("hic" if rakip_en_iyi.get(f["oligo"]) is None
                                               else rakip_en_iyi.get(f["oligo"], "")),
                geri_rakip_en_iyi_uyumsuzluk=("hic" if rakip_en_iyi.get(r["oligo"]) is None
                                              else rakip_en_iyi.get(r["oligo"], "")),
                uye_urunleri=";".join("%s=%d" % (t, prods[t]) for t in ing_tags),
                ceza=round(pen, 3)))
    print(u'\npair construction')
    print("   elenen, cift Tm farki             : %d" % n_tmd)
    print("   elenen, bir uyede urun yok        : %d" % n_noprod)
    if fail_member:
        worst = sorted(fail_member.items(), key=lambda x: -x[1])[:5]
        print("      en cok engelleyen uyeler: %s"
              % ", ".join("%s=%d" % kv for kv in worst))
    print("   elenen, rakipte urun olusuyor     : %d" % n_comp)
    print("   elenen, yetim primer yok          : %d" % n_orph)
    print("   elenen, hetero-dimer dG           : %d" % n_het)
    print("gecerli cift sayisi                  : %d%s"
          % (len(pairs), "  (--stop-after ile erken durduruldu)" if durdu else ""))
    if not pairs:
        sys.exit("gecerli cift bulunamadi")
    pairs.sort(key=lambda x: x["ceza"])
    # Ayni lokusun IUPAC kardes varyantlari tek satira indirilir. Aksi halde
    # ilk on aday tek bir bolgenin alel varyantlariyla dolar ve tablo
    # cesitlilik yitirir. Tutulan, en dusuk cezali varyanttir; kac kardes
    # oldugu ayri sutunda raporlanir.
    if not a.varyantlari_tut:
        onceki = len(pairs)
        secili, gorulen = [], {}
        for pr in pairs:
            anahtar = (pr["ileri_baslangic"], pr["ileri_uzunluk"],
                       pr["geri_baslangic"], pr["geri_uzunluk"])
            if anahtar in gorulen:
                gorulen[anahtar]["lokus_varyant_sayisi"] += 1
                continue
            pr["lokus_varyant_sayisi"] = 1
            gorulen[anahtar] = pr
            secili.append(pr)
        pairs = secili
        print("   lokus birlestirme: %d cift -> %d lokus" % (onceki, len(pairs)))
    else:
        for pr in pairs:
            pr["lokus_varyant_sayisi"] = 1
    pairs = pairs[:a.max_pairs]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pairs[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(pairs)
    print("\nyazildi: %s" % a.out)
    print(u'\nFive best candidates:')
    for p in pairs[:5]:
        print("  ceza=%.2f  urun %d-%d bp  yetim=%s  F=%s (Tm %.1f)  R=%s (Tm %.1f)"
              % (p["ceza"], p["urun_min"], p["urun_maks"], p["yetim_primer"],
                 p["ileri_dizi"], p["ileri_tm"], p["geri_dizi"], p["geri_tm"]))


def _one_config(b_plus, b_minus, ln_plus, ln_minus, L, pmin, pmax):
    """Bir primer arti zincirde, oteki eksi zincirde. Koordinat cevrimi:
    eksi zincirdeki m pozisyonunun arti zincirdeki karsiligi L-1-m'dir."""
    best = None
    for fend, _ in b_plus:
        fstart = fend - ln_plus + 1
        for rend, _ in b_minus:
            r_left_plus = L - 1 - rend
            if fend >= r_left_plus:             # 3' uclari birbirine bakmiyor
                continue
            r_right_plus = r_left_plus + ln_minus - 1
            plen = r_right_plus - fstart + 1
            if plen < pmin:
                continue
            if pmax and plen > pmax:
                continue
            if best is None or plen < best:
                best = plen
    return best


def product_len(bf, br, lnf, lnr, a, pmax=None, pmin=None):
    """Kalip cift sarmal oldugu icin IKI konfigurasyon da gercek urun verir:
    (1) birinci primer arti zincirde, ikincisi eksi zincirde
    (2) birinci primer eksi zincirde, ikincisi arti zincirde
    Yalnizca birine bakmak, ters saklanmis bir diziyi 'urun yok' diye
    isaretler; rakiplerde bu, ozgulluk denetimini tamamen atlatir.
    Doner: gecerli en kisa urun uzunlugu ya da None."""
    L = bf.get("L") or br.get("L")
    if not L:
        return None
    if pmax is None:
        pmax = a.prod_hard_max
    if pmin is None:
        pmin = a.prod_min
    cands = [
        _one_config(bf["plus"], br["minus"], lnf, lnr, L, pmin, pmax),
        _one_config(br["plus"], bf["minus"], lnr, lnf, L, pmin, pmax),
    ]
    cands = [c for c in cands if c is not None]
    return min(cands) if cands else None


if __name__ == "__main__":
    main()
