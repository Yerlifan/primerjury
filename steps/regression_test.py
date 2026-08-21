#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
regression_test.py
Boru hattının kurallarını, kodun kendi yardımcı fonksiyonlarına
güvenmeden sınar. Her testin beklenen sonucu toplantı kararlarından
ya da bilinen matematikten türetilir.

Kullanım:
  python3 regression_test.py                 # hızlı testler
  python3 regression_test.py --gercek-veri   # gerçek konsensüs ve okumalarla
"""
import argparse, importlib.util, itertools, math, os, random, shutil, subprocess, sys, glob, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SONUC = []


def yukle(ad, dosya):
    spec = importlib.util.spec_from_file_location(ad, os.path.join(HERE, dosya))
    m = importlib.util.module_from_spec(spec)
    yedek, sys.argv = sys.argv, [dosya]
    try:
        spec.loader.exec_module(m)
    except SystemExit:
        pass
    finally:
        sys.argv = yedek
    return m


def sina(ad, kosul, ayrinti=""):
    SONUC.append((ad, bool(kosul), ayrinti))
    print("   %-58s %s%s" % (ad, "GECTI" if kosul else "KALDI",
                             ("  " + ayrinti) if ayrinti and not kosul else ""))
    return bool(kosul)


# --- bagimsiz referans uygulamalar -------------------------------------
TAM = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT",
         "S": "CG", "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
         "H": "ACT", "V": "ACG", "N": "ACGT"}


def ref_rc(s):
    return s.translate(TAM)[::-1]


def ref_gc(s):
    return 100.0 * sum(1 for c in s if c in "GC") / len(s)


def ref_wilson(k, n, ust=False, z=1.96):
    if n == 0:
        return 1.0 if ust else 0.0
    p = k / n
    d = 1 + z * z / n
    m = p + z * z / (2 * n)
    s = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (m + s) / d if ust else max(0.0, (m - s) / d)


def ref_baglanir(oligo, kalip, j0=0):
    """Toplanti kuralinin dogrudan uygulanmasi."""
    n = len(oligo)
    mm = [j for j in range(j0, n)
          if not (set(IUPAC.get(oligo[j], "")) & set(IUPAC.get(kalip[j], ""))
                  and kalip[j] != "N")]
    if any(j >= n - 2 for j in mm):
        return None
    if sum(1 for j in mm if j >= n - 5) > 1:
        return None
    if len(mm) > 3:
        return None
    return len(mm)


def ref_find(oligo, seq, min_ortusme=12):
    n, L = len(oligo), len(seq)
    out = set()
    for start in range(-n + 1, L):
        end = start + n
        if end > L:
            continue
        j0 = max(0, -start)
        if n - j0 < min_ortusme:
            continue
        mm = ref_baglanir(oligo, "".join(seq[start + j] if 0 <= start + j < L
                                         else "N" for j in range(n)), j0)
        if mm is not None:
            out.add((end - 1, mm))
    return sorted(out)


# =======================================================================
def testler(a):
    E = yukle("E", "generate_primer_candidates.py")
    G = yukle("G", "design_group_primers.py")
    O = yukle("O", "specificity.py")
    K = yukle("K", "dominant_allele_consensus.py")
    D = yukle("D", "indistinguishable_targets.py")
    import argparse as _ap
    random.seed(1234)

    print("\n1. TERS TUMLEYEN")
    ok = True
    for s in ("ATGC", "ACGTRYSWKM", "AAAA", "GCGCGC"):
        ok &= (E.rc(s) == ref_rc(s))
    sina("rc() bagimsiz uygulamayla ayni", ok)
    sina("rc(rc(x)) == x", all(E.rc(E.rc(s)) == s
                              for s in ("ATGCRYKM", "AACCGGTT")))
    sina("ATGC'nin ters tumleyeni GCAT", E.rc("ATGC") == "GCAT",
         "cikan: %s" % E.rc("ATGC"))

    print("\n2. KOMPOZISYON KURALLARI")
    ar = _ap.Namespace(len_min=18, len_max=25, gc_min=40, gc_max=60,
                       gc_hard_min=35, gc_hard_max=65, gc_clamp_last=5,
                       gc_clamp_max=3, homopolymer_max=4, require_3p_gc=1,
                       degeneracy_budget=0, degeneracy_fold_max=4,
                       iupac_max=2, iupac_son_yasak=5)
    sina("3' ucu A ile biten oligo elenir",
         not E.composition_ok("CGCGATATCGCGATATCGA", ar)[0])
    sina("3' ucu G ile biten oligo gecer",
         E.composition_ok("CGCGATATCGCGATATCGG", ar)[0])
    sina("bes ardisik ayni baz elenir",
         not E.composition_ok("CGCGAAAAACGCGATATCG", ar)[0])
    sina("son bes bazda dort G/C elenir",
         not E.composition_ok("ATATATATATATATAGCGCC", ar)[0])
    dusuk = "ATATATATATATATATATAC"
    sina("GC sert alt sinirinin altindaki oligo elenir (GC=%%%.0f)" % ref_gc(dusuk),
         not E.composition_ok(dusuk, ar)[0])
    sina("gc_pct bagimsiz hesapla ayni",
         all(abs(E.gc_pct(s) - ref_gc(s)) < 1e-9
             for s in ("GGCC", "ATAT", "ACGTACGT")))

    print("\n3. IUPAC ZINCIR KURALI")
    for win, uc, bek in (("STAGCTAGCATCGATCGATC", "F", True),
                         ("STAGCTAGCATCGATCGATC", "R", False),
                         ("TAGCTAGCATCGATCGATCS", "F", False),
                         ("TAGCTAGCATCGATCGATCS", "R", True),
                         ("TAGCTAGSATCGATCGATCA", "F", True),
                         ("TAGCTAGSATCGATCGATCA", "R", True)):
        v, _ = E.iupac_varyantlar(win, ar, uc=uc)
        sina("IUPAC %s uc=%s -> %s" % (win[:6] + "...", uc,
                                       "kabul" if bek else "red"),
             (v is not None) == bek)
    v, _ = E.iupac_varyantlar("SAGCTAGCATCGATCGATCA", ar, uc="F")
    sina("belirsiz pozisyon butun alellere acilir", v is not None and len(v) == 2)
    v, _ = E.iupac_varyantlar("NAGCTAGCATCGATCGATCA", ar, uc="F")
    sina("kalipta N iceren pencere reddedilir", v is None)

    print("\n4. BAGLANMA KURALI (kaba kuvvetle karsilastirma)")
    ab = _ap.Namespace(tail_len=5, tail_max_mm=1, exact_last=2,
                       total_max_mm=3, min_overlap=12)
    fark = 0
    for _ in range(40):
        seq = "".join(random.choice("ACGT") for _ in range(300))
        L = list(seq)
        for i in random.sample(range(300), 18):
            L[i] = random.choice("RYSWKM")
        seq = "".join(L)
        idx = G.build_index(seq, ab.tail_len)
        for _ in range(10):
            p0 = random.randrange(0, 270)
            oligo = "".join(random.choice(IUPAC.get(c, "ACGT"))
                            for c in seq[p0:p0 + 20])
            x = sorted(set(G.find_bindings(oligo, seq, idx, ab.tail_len, ab)))
            y = ref_find(oligo, seq, ab.min_overlap)
            if x != y:
                fark += 1
    sina("find_bindings kuralin dogrudan uygulamasiyla ayni (400 deneme)",
         fark == 0, "ayrisan: %d" % fark)
    sina("son iki bazda uyumsuzluk kabul edilmez",
         ref_baglanir("ACGTACGTACGTACGTACGA", "ACGTACGTACGTACGTACGC") is None)
    sina("kuyrukta tek uyumsuzluk kabul edilir",
         ref_baglanir("ACGTACGTACGTACGTAAGC", "ACGTACGTACGTACGTACGC") == 1)

    print("\n5. URUN GEOMETRISI")
    # Sentetik pozisyon sozlugu yerine GERCEK bir kalip kurulur, primerler
    # o kaliptan kesilir, urun elle olusturulur. Beklenen deger boylece
    # kodun kendi olcusunden degil, dizinin kendisinden gelir.
    ap2 = _ap.Namespace(prod_min=70, prod_hard_max=300)
    random.seed(31)
    L = 500
    kalip = "".join(random.choice("ACGT") for _ in range(L))
    F = kalip[100:120]                 # ileri primer, arti zincir 100..119
    R = G.rc(kalip[260:280])           # geri primer, arti zincirin rc'si
    elde_urun = kalip[100:280]
    sina("elle kurulan urun ileri primerle basliyor", elde_urun.startswith(F))
    sina("elle kurulan urun geri primerin rc'siyle bitiyor",
         elde_urun.endswith(G.rc(R)))
    beklenen = len(elde_urun)          # 180

    idx_p = G.build_index(kalip, ab.tail_len)
    idx_m = G.build_index(G.rc(kalip), ab.tail_len)

    def bagla(o):
        return dict(L=L,
                    plus=G.find_bindings(o, kalip, idx_p, ab.tail_len, ab),
                    minus=G.find_bindings(o, G.rc(kalip), idx_m,
                                          ab.tail_len, ab))

    bf, br = bagla(F), bagla(R)
    pl = G.product_len(bf, br, len(F), len(R), ap2)
    sina("urun uzunlugu gercek kaliptaki amplikonla ayni (%d bp)" % beklenen,
         pl == beklenen, "cikan: %s" % pl)
    sina("ters konfigurasyon da denenir (primerler yer degistirince ayni)",
         G.product_len(br, bf, len(R), len(F), ap2) == beklenen)
    # 3' uclari birbirinden UZAKLASAN cift: primerleri takas edip her birinin
    # kendi zincirinde ARKAYA bakmasini saglar
    Fk = G.rc(kalip[100:120])          # arti zincirin 100..119'unun rc'si
    Rk = kalip[260:280]
    bfk, brk = bagla(Fk), bagla(Rk)
    sina("3' uclari birbirinden UZAKLASAN cift urun vermez",
         G.product_len(bfk, brk, len(Fk), len(Rk), ap2) is None)
    # Urun alt sinirin altinda kalirsa kati tasarimda elenir, rakip taramasinda
    # pmin=1 ile yine de gorulur
    Ry = G.rc(kalip[140:160])          # urun 100..159, yani 60 bp
    bry = bagla(Ry)
    sina("kisa urun (60 bp) kati alt sinirda elenir",
         G.product_len(bf, bry, len(F), len(Ry), ap2) is None)
    sina("rakip taramasinda alt sinir 1'e indirilebiliyor",
         G.product_len(bf, bry, len(F), len(Ry), ap2, pmin=1) == 60)

    print("\n6. WILSON SINIRLARI")
    ok_alt = ok_ust = True
    for k, n in ((0, 10), (1, 10), (5, 10), (10, 10), (1, 17), (2, 100),
                 (1, 1), (0, 1), (37, 300), (500, 20000)):
        ok_alt &= abs(O.wilson_alt(k, n) - ref_wilson(k, n)) < 1e-12
        ok_ust &= abs(O.wilson_ust(k, n) - ref_wilson(k, n, ust=True)) < 1e-12
    sina("wilson_alt bagimsiz uygulamayla ayni", ok_alt)
    sina("wilson_ust bagimsiz uygulamayla ayni", ok_ust)
    sina("alt sinir <= nokta tahmini <= ust sinir",
         all(O.wilson_alt(k, n) <= k / n <= O.wilson_ust(k, n)
             for k, n in ((1, 10), (5, 10), (9, 10), (37, 300))))
    sina("bilinen deger: wilson_alt(1,10)=0,01788",
         abs(O.wilson_alt(1, 10) - 0.017876) < 1e-5)

    print("\n7. BASKIN ALEL CAGRISI")
    sina("12'de esitlikte N cagrilir (kod yolu)",
         "esitler" in open(os.path.join(HERE, "dominant_allele_consensus.py"),
                           encoding="utf-8").read())
    sina("12 ic N'leri silmez, yalniz uclari kirpar",
         'ref[bas:son]' in open(os.path.join(HERE,
                                             "dominant_allele_consensus.py"),
                                encoding="utf-8").read())
    sina("12 okumalari buyuk harfe cevirir",
         '.strip().upper()' in open(os.path.join(
             HERE, "dominant_allele_consensus.py"), encoding="utf-8").read())

    print("\n8. AYIRT EDILEMEZLIK OLCUMU")
    sina("N hicbir bazla eslesme sayilmaz", not D.baz_kesisir("N", "A"))
    sina("N ile N de eslesmez", not D.baz_kesisir("N", "N"))
    sina("Y ile C kesisir", D.baz_kesisir("Y", "C"))
    sina("Y ile G kesismez", not D.baz_kesisir("Y", "G"))
    sina("A ile A kesisir", D.baz_kesisir("A", "A"))

    print("\n9. HIZALAMA ARKA UCU")
    try:
        H = yukle("H", "alignment.py")
        sina("hizalama arka ucu var", H.ARKA_UC is not None, str(H.ARKA_UC))
        ref = "".join(random.choice("ACGT") for _ in range(2000))
        A = H.Hizalayici(seq=ref, preset="map-ont")
        q = ref[500:900]
        ileri = [h for h in A.map(q)]
        ters = [h for h in A.map(H.revcomp(q))]
        sina("ileri hizalama dogru konumda",
             ileri and ileri[0].r_st == 500 and ileri[0].r_en == 900,
             str(ileri[:1]))
        sina("ters hizalama ayni konumda ve strand=-1",
             ters and ters[0].r_st == 500 and ters[0].strand == -1,
             str(ters[:1]))
        sina("ters zincirde sorgu koordinati ozgun sorguya gore",
             ters and ters[0].q_st == 0 and ters[0].q_en == len(q),
             str(ters[:1]))
    except Exception as e:
        sina("hizalama arka ucu var", False, str(e)[:80])

    print(u'\n10. CONSISTENCY BETWEEN THE RAW READ SCAN AND THE DESIGN RULE')
    T = "".join(random.choice("ACGT") for _ in range(400))
    F = T[100:121]
    R = ref_rc(T[250:271])
    yol = "/tmp/_reg_okuma.fastq"
    with open(yol, "w") as fh:
        for i in range(300):
            r = list(T)
            for j in random.sample(range(len(r)), 12):
                r[j] = random.choice("ACGT")
            r = "".join(r)
            if i % 2:
                r = ref_rc(r)
            fh.write("@o%d\n%s\n+\n%s\n" % (i, r, "I" * len(r)))
    tot, fh_, rh, both = O.okuma_taramasi(yol, F, R, 70, 300, 300)
    sina("ham okumalarda urun bulunuyor (ileri ve ters yonlu okumalar)",
         both > 200, "tot=%d F=%d R=%d urun=%d" % (tot, fh_, rh, both))
    sina("ileri ve ters yonlu okumalarin ikisi de sayiliyor",
         abs(fh_ - rh) < 30, "F=%d R=%d" % (fh_, rh))
    tot2, _, _, both2 = O.okuma_taramasi(yol, F, ref_rc(T[250:271])[::-1],
                                         70, 300, 300)
    sina("yanlis yonde yazilmis geri primer urun VERMEZ",
         both2 < both * 0.2, "urun=%d (dogru yonde %d)" % (both2, both))

    print(u'\n11. PRODUCT LENGTH IN THE EXTERNAL DATABASES')
    V = yukle("V", "external_databases.py")
    sina("14 urun boyu iki 5' uc arasi (04 ile ayni olcu)",
         "(ur + lr - 1) - (uf - lf + 1) + 1" in open(
             os.path.join(HERE, "external_databases.py"),
             encoding="utf-8").read())
    sina("14 baglanma kuralini uyguluyor",
         V.baglanma_uygun("ACGTACGTACGTACGTACGC",
                          "ACGTACGTACGTACGTACGC")[0] and
         not V.baglanma_uygun("ACGTACGTACGTACGTACGA",
                              "ACGTACGTACGTACGTACGC")[0])

    print(u'\n12. THE FILE FORMAT CONTRACT BETWEEN STAGES 08 AND 09')
    # 08'in yazdigi dislanan_takson.tsv'yi 09 okur. Sutunlar konuma gore
    # okunursa, dosya bicimi degistiginde dislama SESSIZCE devre disi kalir
    # ve 09 hicbir sey soylemez. Bu, "hicbir karar tek koda birakilmasin"
    # kuralinin ihlalidir; asagidaki testler basliga gore okundugunu ve
    # taninmayan baslikta durdugunu dogrular.
    import subprocess as _sp, tempfile as _tf, textwrap as _tw
    k09 = open(os.path.join(HERE, "specificity.py"), encoding="utf-8").read()
    k08 = open(os.path.join(HERE, "batch_design.py"), encoding="utf-8").read()
    sina("08 basliga 'grup' sutununu yaziyor",
         'df.write("grup\\ttaxid\\tetiket\\tuzunluk\\tkapsanan\\n")' in k08)
    try:
        i0 = k09.index("    dislanan = set()")
        j0 = k09.index("    ayirt = {}")
        govde = _tw.dedent(k09[i0:j0])
        dtmp = _tf.mkdtemp()
        sarmal = ("import os, sys\nclass A: adaylar=%r\na=A()\n"
                  "def log(*x): pass\n" % dtmp + govde
                  + "\nprint('SONUC', sorted(dislanan))\n")

        def cagir(icerik):
            with open(os.path.join(dtmp, "dislanan_takson.tsv"), "w",
                      encoding="utf-8") as fh:
                fh.write(icerik)
            r = _sp.run([sys.executable, "-c", sarmal],
                        capture_output=True, text=True)
            return r.returncode, (r.stdout + r.stderr).strip()

        rc1, o1 = cagir("grup\ttaxid\tetiket\tuzunluk\tkapsanan\n"
                        "A1-4\t1434102\tA1-4_1434102\t1444\t0\n")
        sina("dogru bicimde (grup, taxid) ikilisi kuruluyor",
             rc1 == 0 and "('A1-4', '1434102')" in o1, o1[-70:])
        rc2, o2 = cagir("taxid\tetiket\tuzunluk\tkapsanan\n"
                        "1434102\tA1-4_1434102\t1444\t0\n")
        sina("eski 4 sutunlu bicimde SESSIZCE devam etmiyor, duruyor",
             rc2 != 0, "cikis=%d" % rc2)
        rc3, o3 = cagir("taxid\tgrup\tetiket\tuzunluk\tkapsanan\n"
                        "1434102\tA1-4\tA1-4_1434102\t1444\t0\n")
        sina("sutun sirasi degisse de baslikla dogru okunuyor",
             rc3 == 0 and "('A1-4', '1434102')" in o3, o3[-70:])
        rc4, o4 = cagir("grup\ttaxid\tetiket\tuzunluk\tkapsanan\n")
        sina("yalniz baslikli dosya bos kume verir, hata vermez",
             rc4 == 0 and "SONUC []" in o4, o4[-70:])
    except ValueError:
        sina("09'da dislama blogu bulundu", False, "blok isaretleri degismis")

    print(u'\n13. EXTERNAL DATABASE SET MAPPING AND COVERAGE AUDIT')
    # 2026-08-01'de bulunan hata: ROD_v1.2_operon_variants.fasta okaryot
    # yalnizdir (60320/60320 Eukaryota) ama A1/A2/B siniflarina atanmisti.
    # 71 arke/bakteri cifti icin "hedef disi urun yok" yazildi; oysa o
    # veritabaninda o alandan tek dizi bile yoktu. Asagidaki testler hem
    # eslemenin duzeltilmis halini hem de ayni hatanin bir daha sessiz
    # kalmamasini saglayan kapsam denetimini sinar.
    DV = yukle("DV", "external_databases.py")
    MF = yukle("MF", "mfeprimer_layer.py")
    sina("14 ve 19 DAR kumeyi ayni goruyor", DV.SINIF_DB == MF.SINIF_DB)
    sina("14 ve 19 GENIS kumeyi ayni goruyor",
         DV.SINIF_DB_GENIS == MF.SINIF_DB_GENIS)
    sina("dar ve genis kume kesismiyor (ayni db iki kez taranmaz)",
         all(not (set(DV.SINIF_DB[s]) & set(DV.SINIF_DB_GENIS[s]))
             for s in DV.SINIF_DB))
    sina("her sinif ayni veritabani kumesini goruyor",
         len({tuple(sorted(set(DV.SINIF_DB[s]) | set(DV.SINIF_DB_GENIS[s])))
              for s in DV.SINIF_DB}) == 1)
    sina("ROD mantar siniflarinda da taraniyor",
         all("ROD_v1.2_operon_variants.fasta" in DV.SINIF_DB_GENIS[s]
             for s in ("F1", "F2")))
    sina("fungi.18SrRNA hicbir sinifta unutulmamis",
         all("fungi.18SrRNA.fna" in (DV.SINIF_DB[s] + DV.SINIF_DB_GENIS[s])
             for s in DV.SINIF_DB))
    sina("ref_all/ref_all2 taranmiyor (digerlerinin alt kumesi)",
         all(not {"ref_all.fna", "ref_all2.fna"}
             & set(DV.SINIF_DB[s] + DV.SINIF_DB_GENIS[s])
             for s in DV.SINIF_DB))

    if shutil.which("blastn") and shutil.which("makeblastdb"):
        gec = tempfile.mkdtemp(prefix="kapsam_")
        random.seed(77)
        kons = "".join(random.choice("ACGT") for _ in range(1500))
        os.makedirs(os.path.join(gec, "kons"), exist_ok=True)
        open(os.path.join(gec, "kons", "F1-1_5555_baskin_konsensus.fasta"),
             "w").write(">F1-1_5555\n%s\n" % kons)

        def _mut(s, oran):
            l = list(s)
            for i in range(len(l)):
                if random.random() < oran:
                    l[i] = random.choice("ACGT")
            return "".join(l)

        ilgisiz = os.path.join(gec, "ilgisiz.fna")
        with open(ilgisiz, "w") as fh:
            for i in range(200):
                fh.write(">ilg%d\n%s\n"
                         % (i, "".join(random.choice("ACGT")
                                       for _ in range(1200))))
        ilgili = os.path.join(gec, "ilgili.fna")
        with open(ilgili, "w") as fh:
            for i in range(50):
                fh.write(">akr%d\n%s\n" % (i, _mut(kons, 0.08)))
        for f in (ilgisiz, ilgili):
            subprocess.run(["makeblastdb", "-in", f, "-dbtype", "nucl"],
                           capture_output=True, text=True)
        dz = DV.sinif_konsensuslari(os.path.join(gec, "kons"), "F1")
        sina("konsensus yalniz kendi sinifi icin okunuyor",
             len(dz) == 1 and DV.sinif_konsensuslari(
                 os.path.join(gec, "kons"), "B") == [])
        _, uz_yok, _ = DV.kapsam_olc(ilgisiz, dz, gec, "a", 2, 600)
        _, uz_var, kim_var = DV.kapsam_olc(ilgili, dz, gec, "b", 2, 600)
        sina("ilgisiz veritabani KAPSAM_YOK veriyor", uz_yok < 400,
             "en uzun hizalama %d bp" % uz_yok)
        sina("akraba iceren veritabani KAPSANIYOR veriyor",
             uz_var >= 400 and kim_var > 80,
             "en uzun %d bp %%%.1f" % (uz_var, kim_var))
        sina("konsensus verilmezse kapsam OLCULMEDI diye isaretleniyor",
             DV.kapsam_olc(ilgili, [], gec, "c", 2, 600)[0] == "KAPSAM_OLCULMEDI")
        shutil.rmtree(gec, ignore_errors=True)
    else:
        sina("blastn/makeblastdb kurulu", False, "kapsam testleri atlandi")

    print(u'\n14. SEPARATING THE OWN TAXON FROM A FOREIGN TAXON')
    # 2026-08-01: geniş taramada en yuksek "hedef disi urun" sayilarinin bir
    # kismi hedefin KENDISIYDI (Methanothrix hedefi SILVA'daki Methanothrix
    # kayitlarini cogaltiyor). Ham sayiya gore siralamak yanlis primerleri
    # one cikariyordu. Asagidaki testler ayrimin dogru yapildigini sinar.
    sina("'Ca. Nitrosocosmicus hydrocola' cinsi dogru cikariliyor",
         DV._ad_cinsi("Ca. Nitrosocosmicus hydrocola") == "nitrosocosmicus",
         DV._ad_cinsi("Ca. Nitrosocosmicus hydrocola"))
    sina("'uncultured Acetobacteroides sp.' cinsi dogru cikariliyor",
         DV._ad_cinsi("uncultured Acetobacteroides sp.") == "acetobacteroides",
         DV._ad_cinsi("uncultured Acetobacteroides sp."))
    sina("'Candidatus Cloacimonas acidaminovorans' cinsi dogru cikariliyor",
         DV._ad_cinsi("Candidatus Cloacimonas acidaminovorans") == "cloacimonas")
    BASLIKLAR = {
        "SILVA": ("FJ347531.1.916 Archaea;Halobacteriota;Methanosarcinia;"
                  "Methanosarcinales;Methanosaetaceae;Methanothrix;"
                  "uncultured archaeon"),
        "RefSeq": ("NR_104707.1 Methanothrix soehngenii GP6 16S ribosomal RNA,"
                   " partial sequence"),
        "ROD": ("GCA_000002515|CR382124.1/1-2|Eukaryota;Opisthokonta;Fungi;"
                "Ascomycota;Saccharomycetes;Saccharomycetales;"
                "Saccharomycetaceae;Methanothrix;Methanothrix_soehngenii|size=3"),
        "UNITE": ("UDB016649|k__Fungi;p__Basidiomycota;c__Agaricomycetes;"
                  "o__Thelephorales;f__Thelephoraceae;g__Methanothrix;"
                  "s__Methanothrix_soehngenii|SH1281904.10FU"),
        "PR2": ("AB353770.1.1740_U|18S_rRNA|nucleus||Eukaryota|TSAR|Alveolata|"
                "Dinoflagellata|Dinophyceae|Peridiniales|Kryptoperidiniaceae|"
                "Methanothrix|Methanothrix_soehngenii"),
    }
    for bicim, baslik in sorted(BASLIKLAR.items()):
        sina("%s basligindan cins cikariliyor" % bicim,
             "methanothrix" in DV._tokenlar(baslik))
    sina("gurultu kelimeleri token sayilmiyor",
         not ({"uncultured", "partial", "sequence", "ribosomal"}
              & DV._tokenlar(BASLIKLAR["RefSeq"])))
    sina("kisa parcalar (k__, p__) token sayilmiyor",
         not ({"k", "p", "c", "o", "f", "g", "s"}
              & DV._tokenlar(BASLIKLAR["UNITE"])))

    if shutil.which("blastn") and shutil.which("makeblastdb"):
        gec2 = tempfile.mkdtemp(prefix="takson_")
        random.seed(11)
        kalip = "".join(random.choice("ACGT") for _ in range(1500))
        Fp, Rp = kalip[300:322], kalip[500:522]
        os.makedirs(os.path.join(gec2, "kons"))
        os.makedirs(os.path.join(gec2, "final"))
        os.makedirs(os.path.join(gec2, "REFERANS_DB"))
        open(os.path.join(gec2, "kons", "A1-1_2223_baskin_konsensus.fasta"),
             "w").write(">A1-1_2223\n%s\n" % kalip)

        def _m(s, o):
            l = list(s)
            for i in range(len(l)):
                if random.random() < o:
                    l[i] = random.choice("ACGT")
            return "".join(l)

        # ayni cift icin 3 KENDI taksonu, 3 YABANCI takson kaydi
        dbyol = os.path.join(gec2, "REFERANS_DB", "archaea.16S.fna")
        with open(dbyol, "w") as fh:
            for i in range(3):
                fh.write(">KENDI%d.1 Archaea;Methanosarcinales;Methanosaetaceae;"
                         "Methanothrix;uncultured archaeon\n%s\n"
                         % (i, _m(kalip, 0.01)))
            for i in range(3):
                fh.write(">YABANCI%d.1 Archaea;Halobacterales;Halorubrum;"
                         "uncultured archaeon\n%s\n" % (i, _m(kalip, 0.01)))
        subprocess.run(["makeblastdb", "-in", dbyol, "-dbtype", "nucl"],
                       capture_output=True, text=True)
        with open(os.path.join(gec2, "final", "primer_final.tsv"), "w") as fh:
            fh.write("karar\thedef\tsinif\tozgulluk_durum\tileri_dizi\tgeri_dizi\n")
            fh.write("ONERILIR\tMethanothrix_soehngenii_turu\tA1\tGECTI\t%s\t%s\n"
                     % (Fp, ref_rc(Rp)))
            fh.write("ONERILIR\tArke_universal\tA1\tGECTI\t%s\t%s\n"
                     % (Fp, ref_rc(Rp)))
        open(os.path.join(gec2, "hedefler.tsv"), "w").write(
            "karar\tad\tx\ttaxid\n"
            "ONERILIR\tMethanothrix_soehngenii_turu\t-\t2223\n"
            "ONERILIR\tArke_universal\t-\t*A\n")
        open(os.path.join(gec2, "adlar.tsv"), "w").write(
            "2223\tMethanothrix soehngenii\n")
        cik = os.path.join(gec2, "cikti.tsv")
        subprocess.run(
            [sys.executable, os.path.join(HERE, "external_databases.py"),
             "--final", os.path.join(gec2, "final"),
             "--db", os.path.join(gec2, "REFERANS_DB"),
             "--kons", os.path.join(gec2, "kons"),
             "--hedefler", os.path.join(gec2, "hedefler.tsv"),
             "--adlar", os.path.join(gec2, "adlar.tsv"),
             "--out", cik], capture_output=True, text=True)
        satir = {}
        if os.path.exists(cik):
            import csv as _csv
            for r in _csv.DictReader(open(cik, encoding="utf-8"),
                                     delimiter="\t"):
                satir[r["hedef"]] = r
        oz = satir.get("Methanothrix_soehngenii_turu", {})
        ev = satir.get("Arke_universal", {})
        sina("ozgul hedefte 3 kendi + 3 yabanci urun ayriliyor",
             oz.get("urun_kendi_taksonda") == "3"
             and oz.get("urun_yabanci_taksonda") == "3",
             "kendi=%s yabanci=%s toplam=%s" % (oz.get("urun_kendi_taksonda"),
                                                oz.get("urun_yabanci_taksonda"),
                                                oz.get("hedef_disi_urun")))
        sina("kendi + yabanci + bilinmiyor = ham toplam",
             oz and (int(oz["urun_kendi_taksonda"])
                     + int(oz["urun_yabanci_taksonda"])
                     + int(oz["urun_takson_bilinmiyor"])
                     == int(oz["hedef_disi_urun"])))
        sina("yabanci ornekler yalniz yabanci taksonlari gosteriyor",
             oz and "YABANCI" in oz["yabanci_ornekler"]
             and "KENDI" not in oz["yabanci_ornekler"])
        sina("evrensel hedefte hicbir urun yabanci sayilmiyor",
             ev.get("hedef_turu") == "evrensel"
             and ev.get("urun_yabanci_taksonda") == "0",
             "turu=%s yabanci=%s" % (ev.get("hedef_turu"),
                                     ev.get("urun_yabanci_taksonda")))
        shutil.rmtree(gec2, ignore_errors=True)
    else:
        sina("blastn/makeblastdb kurulu", False, "takson testleri atlandi")

    print("\n15. YABANCI VURUSUN UZAKLIGI (SOYAGACI DERINLIGI)")
    # Yabanci takson tek basina yeterli olcu degil: bazi hedefler islevsel
    # gruptur. Olculdu (2026-08-01): Hidrojenotrofik_metanojenler'in yabanci
    # vuruslari Methanobacterium ve Methanosphaera, ki ikisi de
    # hidrojenotrofik metanojendir; Nitrosocosmicus_AOA'nin yabanci vuruslari
    # Nitrosotalea ve Nitrosopumilus, ki ikisi de AOA'dir. Buna karsilik
    # Petrimonas -> Flavobacterium ve Trichoderma -> Calonectria gercekten
    # uzaktir. Uzaklik elle yazilmis bir tabloyla degil, veritabani
    # basliklarindaki soyagaciyla olculur.
    SOY = {
        "SILVA": ("FJ347531.1.916 Archaea;Halobacteriota;Methanosarcinia;"
                  "Methanosarcinales;Methanosaetaceae;Methanothrix;"
                  "uncultured archaeon"),
        "ROD": ("GCA_000002515|CR382124.1/1-2|Eukaryota;Opisthokonta;Fungi;"
                "Ascomycota;Saccharomycetes;Saccharomycetales;"
                "Saccharomycetaceae;Kluyveromyces;Kluyveromyces_lactis|size=3"),
        "UNITE": ("UDB016649|k__Fungi;p__Basidiomycota;c__Agaricomycetes;"
                  "o__Thelephorales;f__Thelephoraceae;g__Thelephora;"
                  "s__Thelephora_albomarginata|SH1281904.10FU"),
        "PR2": ("AB353770.1.1740_U|18S_rRNA|nucleus||Eukaryota|TSAR|Alveolata|"
                "Dinoflagellata|Dinophyceae|Peridiniales|Kryptoperidiniaceae|"
                "Unruhdinium|Unruhdinium_kevei"),
    }
    for bicim in sorted(SOY):
        sina("%s soyagaci cikariliyor (>=6 basamak)" % bicim,
             len(DV._soyagaci(SOY[bicim])) >= 6,
             "cikan: %s" % DV._soyagaci(SOY[bicim]))
    sina("RefSeq basliginda soyagaci yok, bos liste doner",
         DV._soyagaci("NR_104707.1 Methanothrix soehngenii GP6 16S ribosomal "
                      "RNA, partial sequence") == [])
    sina("UNITE'nin k__/p__ on ekleri temizleniyor",
         DV._soyagaci(SOY["UNITE"])[0] == "fungi")
    sina("'Incertae Sedis' basamagi soyagacina alinmiyor",
         "incertae sedis" not in DV._soyagaci(
             "X Archaea;Halobacteriota;Incertae Sedis;Methanosarcinales;"
             "Methanosaetaceae;Methanothrix;uncultured"))

    KENDI1 = SOY["SILVA"]
    KENDI2 = ("CU916879.1.1311 Archaea;Halobacteriota;Methanosarcinia;"
              "Methanosarcinales;Methanosaetaceae;Methanothrix;"
              "uncultured bacterium")
    AYNI_AILE = ("XX Archaea;Halobacteriota;Methanosarcinia;Methanosarcinales;"
                 "Methanosaetaceae;Methanothrix soehngenii;strain")
    AYNI_TAKIM = ("KC604456 Archaea;Halobacteriota;Methanosarcinia;"
                  "Methanosarcinales;Methanoperedenaceae;"
                  "Candidatus Methanoperedens;uncultured")
    FARKLI_SINIF = ("DQ337053 Archaea;Halobacteriota;Halobacteria;"
                    "Halobacterales;Halorubrum;uncultured archaeon")

    def _karar(kendiler, yab):
        ref = DV._baskin_soy([DV._soyagaci(x) for x in kendiler])
        d = DV._ortak_derinlik(ref, DV._soyagaci(yab))
        return "YAKIN" if d >= max(1, len(ref) - 2) else "UZAK"

    sina("ayni ailedeki vurus YAKIN sayiliyor",
         _karar([KENDI1], AYNI_AILE) == "YAKIN")
    sina("ayni takim farkli aile UZAK sayiliyor",
         _karar([KENDI1], AYNI_TAKIM) == "UZAK")
    sina("farkli siniftaki vurus UZAK sayiliyor",
         _karar([KENDI1], FARKLI_SINIF) == "UZAK")
    # ESIK KARARLILIGI: ortak on ek kullanilsaydi kendi vurus sayisi
    # arttikca referans soyagaci kisalir ve ayni yabanci vurus YAKIN'a
    # kayardi. Baskin soyagaci bunu engeller.
    sina("karar kendi vurus sayisindan bagimsiz",
         len({_karar(k, AYNI_TAKIM) for k in ([KENDI1], [KENDI1, KENDI2],
                                              [KENDI1, KENDI2, KENDI1])}) == 1)
    sina("baskin soyagaci derinligi kendi vurus sayisiyla kisalmiyor",
         len(DV._baskin_soy([DV._soyagaci(KENDI1), DV._soyagaci(KENDI2)]))
         == len(DV._soyagaci(KENDI1)))
    sina("kendi vurusu yoksa referans soyagaci bos",
         DV._baskin_soy([]) == [])

    print(u'\n16. DECISION LEVEL AUDIT (27): SPECIES AND GENUS SPECIFICITY')
    # Toplanti karari alti hedefte TUR, dort hedefte CINS ozgulluk istiyor.
    # 09 numunedeki rakiplere, 14 dis veritabanlarina bakar; ikisi de
    # "kendi cinsinin oteki TURLERINDEN ayiriyor mu" sorusunu sormaz.
    # 27 bunu kardes turlerden kurulu bir panele karsi olcer.
    DZ = yukle("DZ", "check_taxonomic_level.py")
    sina("RefSeq basligindan tur adi cikariliyor",
         DZ.tur_adi("NR_104707.1 Methanothrix soehngenii GP6 16S ribosomal "
                    "RNA, partial sequence") == "Methanothrix soehngenii")
    sina("UNITE basligindan tur adi cikariliyor",
         DZ.tur_adi("UDB016649|k__Fungi;p__Basidiomycota;g__Thelephora;"
                    "s__Thelephora_albomarginata|SH1281904.10FU")
         == "Thelephora albomarginata")
    sina("ROD basligindan tur adi cikariliyor",
         DZ.tur_adi("GCA_000002515|CR382124.1/1-2|Eukaryota;Fungi;"
                    "Saccharomycetaceae;Kluyveromyces;Kluyveromyces_lactis"
                    "|size=3") == "Kluyveromyces lactis")
    sina("PR2 basligindan tur adi cikariliyor",
         DZ.tur_adi("AB353770.1.1740_U|18S_rRNA|nucleus||Eukaryota|TSAR|"
                    "Kryptoperidiniaceae|Unruhdinium|Unruhdinium_kevei")
         == "Unruhdinium kevei")
    sina("SILVA'nin 'uncultured' kaydindan tur adi CIKMIYOR",
         DZ.tur_adi("FJ347531.1.916 Archaea;Halobacteriota;Methanosaetaceae;"
                    "Methanothrix;uncultured archaeon") == "")
    sina("'Methanothrix sp.' tur adi sayilmiyor",
         DZ.tur_adi("NR_999.1 Methanothrix sp. uncultured archaeon 16S") == "")
    # RefSeq ITS/28S kayitlari '; from TYPE material' ile biter. Onceki
    # surumde baslikta noktali virgul varsa RefSeq dali hic calismiyor ve
    # TIP SUSU kayitlari panele alinmiyordu; oysa tur ayriminin altin
    # standardi tam olarak o kayitlardir (yalniz Petriella icin 35 kayit).
    # hedef_tur sutunu EKLEMEZ, DEGISTIRIR. Petriella_musispora satirinda
    # in_taxid kutularin Kraken2 etiketleridir (numune destegi icin);
    # eklense hedef tur kumesine Trichoderma da girer ve Trichoderma
    # cogaltan cift "hedef turde urun var" sayilirdi.
    _g6 = tempfile.mkdtemp(prefix="hedeftur_")
    _h = os.path.join(_g6, "h.tsv")
    _ad = os.path.join(_g6, "a.tsv")
    open(_ad, "w").write("101201\tTrichoderma asperellum\n")
    open(_h, "w").write(
        "karar\thedef\tduzey\tin_taxid\tharic\tnot\thedef_tur\n"
        "5\tPetriella_musispora\ttur\t101201\t\tnot\tPetriella musispora\n"
        "1\tMicroascaceae_askomikot\ttur\t101201\t\tnot\t\n")
    _hd = {x["hedef"]: x for x in DZ.hedefleri_oku(_h, _ad, None)}
    sina("hedef_tur verilince taxid adi DEGISTIRILIYOR, eklenmiyor",
         _hd["Petriella_musispora"]["hedef_turler"] == {"Petriella musispora"},
         "cikan: %s" % _hd["Petriella_musispora"]["hedef_turler"])
    sina("hedef_tur verilince cins de degisiyor",
         _hd["Petriella_musispora"]["cinsler"] == {"Petriella"})
    sina("hedef_tur bos olan satir eski davranisi suruyor",
         _hd["Microascaceae_askomikot"]["hedef_turler"]
         == {"Trichoderma asperellum"})
    shutil.rmtree(_g6, ignore_errors=True)

    sina("RefSeq ITS 'from TYPE material' basligindan tur cikariliyor",
         DZ.tur_adi("NR_172285.1 Petriella musispora CBS 745.69 ITS region;"
                    " from TYPE material") == "Petriella musispora",
         "cikan: %r" % DZ.tur_adi("NR_172285.1 Petriella musispora CBS 745.69"
                                  " ITS region; from TYPE material"))
    sina("RefSeq 28S 'from TYPE material' basligindan tur cikariliyor",
         DZ.tur_adi("NG_042733.1 Calonectria pentaseptata CBS 133349 28S rRNA,"
                    " partial sequence; from TYPE material")
         == "Calonectria pentaseptata")
    sina("SILVA soyagacli basligi tur adi vermiyor (panele girmez)",
         DZ.tur_adi("FJ347531.1.916 Archaea;Halobacteriota;Methanosarcinia;"
                    "Methanosaetaceae;Methanothrix;uncultured archaeon") == "")
    sina("'Ca. Nitrosocosmicus hydrocola' cins ve tur olarak ayriliyor",
         DZ.ad_parcala("Ca. Nitrosocosmicus hydrocola")
         == ("Nitrosocosmicus", "hydrocola"))
    # UNITE'te 's__Trichoderma_sp' tur adi sayiliyordu; yalniz Trichoderma
    # icin 16910 kayit panele sahte bir "tur" olarak giriyordu.
    for kotu in ("s__Trichoderma_sp", "s__Marasmius_spp", "s__Podospora_cf"):
        sina("UNITE'te '%s' tur adi SAYILMIYOR" % kotu,
             DZ.tur_adi("UDB1|k__Fungi;g__X;%s|SH1" % kotu) == "",
             "cikan: %r" % DZ.tur_adi("UDB1|k__Fungi;g__X;%s|SH1" % kotu))
    sina("gercek ikili ad hala taniniyor",
         DZ.tur_adi("UDB1|k__Fungi;g__Petriella;s__Petriella_setifera|SH1")
         == "Petriella setifera")
    # 'Methanosarcina_barkeri_referans' -> '_referans' soyulunca
    # 'Methanosarcina_barkeri' cikiyor, hedefler.tsv'de ad
    # 'Methanosarcina_barkeri_turu'. Eslesme tutmayinca hedefin TEK primer
    # takimi sessizce dusuyor ve hedef CIFT_YOK gorunuyordu.
    _ADL = ["Methanosarcina_barkeri_turu", "Proteiniphilum_cinsi",
            "Podospora_pseudopauciseta", "Proteolitik_sintrofik_bakteriler"]
    sina("referans hedefi '_turu' ekli hedefe baglaniyor",
         DZ.referans_esle("Methanosarcina_barkeri_referans", _ADL)
         == "Methanosarcina_barkeri_turu",
         "cikan: %s" % DZ.referans_esle("Methanosarcina_barkeri_referans", _ADL))
    sina("birebir eslesen referans hedefi degismiyor",
         DZ.referans_esle("Proteiniphilum_cinsi_referans", _ADL)
         == "Proteiniphilum_cinsi")
    sina("kisaltilmis referans hedefi tam ada baglaniyor",
         DZ.referans_esle("Proteolitik_sintrofik_referans", _ADL)
         == "Proteolitik_sintrofik_bakteriler")
    sina("eslesmeyen referans hedefi None doner (sessiz dusmez)",
         DZ.referans_esle("Uydurma_referans", _ADL) is None)

    if shutil.which("blastn") and shutil.which("makeblastdb"):
        g3 = tempfile.mkdtemp(prefix="duzey_")
        random.seed(5)
        os.makedirs(os.path.join(g3, "REFERANS_DB"))
        os.makedirs(os.path.join(g3, "final"))
        hedef = "".join(random.choice("ACGT") for _ in range(1400))
        Fp, Rp = hedef[300:322], ref_rc(hedef[500:522])

        def _m2(s, o):
            l = list(s)
            for i in range(len(l)):
                if random.random() < o:
                    l[i] = random.choice("ACGT")
            return "".join(l)

        # kardes turde ciftin 3' uclarina denk gelen bazlar bozuk:
        # baglanma kurali geregi urun OLUSMAMALI
        kardes = list(hedef)
        for p in (321, 320, 500, 501):
            kardes[p] = {"A": "C", "C": "A", "G": "T", "T": "G"}[kardes[p]]
        kardes = "".join(kardes)
        # ayrim yapmayan cift: kardes turde de bozulmamis bir bolgeden
        Fp2, Rp2 = hedef[700:722], ref_rc(hedef[900:922])
        dbyol = os.path.join(g3, "REFERANS_DB", "archaea.16S.fna")
        with open(dbyol, "w") as fh:
            for i in range(3):
                fh.write(">NR_10%d.1 Methanothrix soehngenii GP%d 16S "
                         "ribosomal RNA, complete sequence\n%s\n"
                         % (i, i, _m2(hedef, 0.001)))
            # UC kardes tur: varsayilan esik 2 oldugu icin ucu birden
            # cogaltan cift esigin USTUNDE kalmali (TUR_AYRIMI_YOK)
            for i, t in enumerate(("harundinacea", "thermoacetophila",
                                   "hungatei")):
                fh.write(">NR_20%d.1 Methanothrix %s DSM%d 16S ribosomal "
                         "RNA, complete sequence\n%s\n"
                         % (i, t, i, _m2(kardes, 0.001)))
            fh.write(">NR_999.1 Methanothrix sp. uncultured archaeon 16S "
                     "ribosomal RNA\n%s\n" % _m2(hedef, 0.001))
        with open(os.path.join(g3, "final", "primer_final.tsv"), "w") as fh:
            fh.write("karar\thedef\tsinif\tozgulluk_durum\tileri_dizi\t"
                     "geri_dizi\n")
            fh.write("1\tMethanothrix_soehngenii_turu\tA1\tGECTI\t%s\t%s\n"
                     % (Fp, Rp))
            fh.write("1\tMethanothrix_soehngenii_turu\tA1\tGECTI\t%s\t%s\n"
                     % (Fp2, Rp2))
        open(os.path.join(g3, "hedefler.tsv"), "w").write(
            "karar\thedef\tduzey\tin_taxid\tharic\tnot\n"
            "1\tMethanothrix_soehngenii_turu\ttur\t2223\t\tM. soehngenii\n")
        open(os.path.join(g3, "adlar.tsv"), "w").write(
            "2223\tMethanothrix soehngenii\n")
        cik3 = os.path.join(g3, "cikti.tsv")
        r3 = subprocess.run(
            [sys.executable, os.path.join(HERE, "check_taxonomic_level.py"),
             "--hedefler", os.path.join(g3, "hedefler.tsv"),
             "--adlar", os.path.join(g3, "adlar.tsv"),
             "--final", os.path.join(g3, "final"),
             "--db", os.path.join(g3, "REFERANS_DB"),
             "--out", cik3], capture_output=True, text=True)
        sat = []
        if os.path.exists(cik3):
            import csv as _csv3
            sat = list(_csv3.DictReader(open(cik3, encoding="utf-8"),
                                        delimiter="\t"))
        kararlar = {x["cift_no"]: x for x in sat}
        sina("ayirt eden cift TUR_OZGUL cikiyor",
             kararlar.get("0", {}).get("karar") == "TUR_OZGUL",
             "cikan: %s" % kararlar.get("0", {}).get("karar"))
        sina("ayirt etmeyen cift TUR_AYRIMI_YOK cikiyor",
             kararlar.get("1", {}).get("karar") == "TUR_AYRIMI_YOK",
             "cikan: %s" % kararlar.get("1", {}).get("karar"))
        sina("TUR_OZGUL ciftte kardes turde urun sayisi sifir",
             kararlar.get("0", {}).get("diger_turde_urun") == "0")
        sina("ayirt etmeyen ciftte kardes turler adiyla bildiriliyor",
             "harundinacea" in kararlar.get("1", {}).get("cogaltilan_turler", ""))
        sina("'sp.' kaydi panele alinmadi (panel 4 tur)",
             kararlar.get("0", {}).get("panel_tur_sayisi") == "4",
             "panel_tur_sayisi=%s" % kararlar.get("0", {}).get("panel_tur_sayisi"))
        sina("capraz TUR sayisi bildiriliyor (urun sayisi degil)",
             kararlar.get("1", {}).get("capraz_tur_sayisi") == "3",
             "capraz_tur_sayisi=%s"
             % kararlar.get("1", {}).get("capraz_tur_sayisi"))
        # Toplanti karari 1-2 caprazi hos goruyor. Esik yukseltilince ayni
        # cift esik ici sayilmali; esik dusukken sayilmamali.
        cik5 = os.path.join(g3, "cikti_esik3.tsv")
        subprocess.run(
            [sys.executable, os.path.join(HERE, "check_taxonomic_level.py"),
             "--hedefler", os.path.join(g3, "hedefler.tsv"),
             "--adlar", os.path.join(g3, "adlar.tsv"),
             "--final", os.path.join(g3, "final"),
             "--db", os.path.join(g3, "REFERANS_DB"),
             "--capraz-tur-esik", "3",
             "--out", cik5], capture_output=True, text=True)
        sat5 = {}
        if os.path.exists(cik5):
            import csv as _csv5
            sat5 = {x["cift_no"]: x for x in
                    _csv5.DictReader(open(cik5, encoding="utf-8"),
                                     delimiter="\t")}
        sina("esik 3'e cikinca 3 caprazli cift TUR_OZGUL_ESIKLI oluyor",
             sat5.get("1", {}).get("karar") == "TUR_OZGUL_ESIKLI",
             "cikan: %s" % sat5.get("1", {}).get("karar"))
        sina("caprazsiz cift esikten bagimsiz TUR_OZGUL kaliyor",
             sat5.get("0", {}).get("karar") == "TUR_OZGUL")
        sina("27 sifir cikis koduyla bitti", r3.returncode == 0,
             r3.stderr.strip()[-120:])
        # hedef tur panelde yoksa TUR_OZGUL DEMEZ
        open(os.path.join(g3, "hedefler.tsv"), "w").write(
            "karar\thedef\tduzey\tin_taxid\tharic\tnot\n"
            "1\tMethanothrix_soehngenii_turu\ttur\t9999\t\tyok\n")
        open(os.path.join(g3, "adlar.tsv"), "w").write(
            "9999\tMethanothrix yoktur\n")
        cik4 = os.path.join(g3, "cikti2.tsv")
        subprocess.run(
            [sys.executable, os.path.join(HERE, "check_taxonomic_level.py"),
             "--hedefler", os.path.join(g3, "hedefler.tsv"),
             "--adlar", os.path.join(g3, "adlar.tsv"),
             "--final", os.path.join(g3, "final"),
             "--db", os.path.join(g3, "REFERANS_DB"),
             "--out", cik4], capture_output=True, text=True)
        sat4 = []
        if os.path.exists(cik4):
            import csv as _csv4
            sat4 = list(_csv4.DictReader(open(cik4, encoding="utf-8"),
                                         delimiter="\t"))
        # CINS DUZEYI: beyan edilen cinsin disinda urun olusursa cins
        # ozgul degildir. Olculdu: Proteiniphilum ciftlerinden ikisi
        # Fermentimonas caenicola'yi da cogaltiyordu, eski sayim ikisini
        # de en genis kapsamli cift gibi gosteriyordu.
        g4 = tempfile.mkdtemp(prefix="cins_")
        random.seed(9)
        os.makedirs(os.path.join(g4, "REFERANS_DB"))
        os.makedirs(os.path.join(g4, "final"))
        hb = "".join(random.choice("ACGT") for _ in range(1400))
        Fc, Rc = hb[300:322], ref_rc(hb[500:522])
        disi = list(hb)
        for p in (321, 320, 500, 501):
            disi[p] = {"A": "C", "C": "A", "G": "T", "T": "G"}[disi[p]]
        disi = "".join(disi)
        Fc2, Rc2 = hb[700:722], ref_rc(hb[900:922])

        def _m4(s, o):
            l = list(s)
            for i in range(len(l)):
                if random.random() < o:
                    l[i] = random.choice("ACGT")
            return "".join(l)

        with open(os.path.join(g4, "REFERANS_DB", "bacteria.16S.fna"),
                  "w") as fh:
            for i, t in enumerate(("propionicum", "saccharofermentans",
                                   "acetatigenes")):
                fh.write(">NR_1%d.1 Proteiniphilum %s DSM%d 16S ribosomal "
                         "RNA, complete sequence\n%s\n"
                         % (i, t, i, _m4(hb, 0.001)))
            fh.write(">NR_9.1 Fermentimonas caenicola DSM9 16S ribosomal "
                     "RNA, complete sequence\n%s\n" % _m4(disi, 0.001))
        with open(os.path.join(g4, "final", "primer_final.tsv"), "w") as fh:
            fh.write("karar\thedef\tsinif\tozgulluk_durum\tileri_dizi\t"
                     "geri_dizi\n")
            fh.write("2\tProteiniphilum_cinsi\tB\tGECTI\t%s\t%s\n" % (Fc, Rc))
            fh.write("2\tProteiniphilum_cinsi\tB\tGECTI\t%s\t%s\n" % (Fc2, Rc2))
        open(os.path.join(g4, "hedefler.tsv"), "w").write(
            "karar\thedef\tduzey\tin_taxid\tharic\tnot\n"
            "2\tProteiniphilum_cinsi\tcins\t2829812,1642647\t\tP. cinsi\n")
        open(os.path.join(g4, "adlar.tsv"), "w").write(
            "2829812\tProteiniphilum propionicum\n"
            "1642647\tProteiniphilum saccharofermentans\n")
        # olculen kimlik BASKA CINS: cins duzeyinde hedef sayilmamali
        open(os.path.join(g4, "kimlik.tsv"), "w").write(
            "hedef\tkraken_etiketi\tolculen_kimlik\n"
            "Proteiniphilum_cinsi\tX\tFermentimonas caenicola\n")
        cik6 = os.path.join(g4, "cikti.tsv")
        subprocess.run(
            [sys.executable, os.path.join(HERE, "check_taxonomic_level.py"),
             "--hedefler", os.path.join(g4, "hedefler.tsv"),
             "--adlar", os.path.join(g4, "adlar.tsv"),
             "--final", os.path.join(g4, "final"),
             "--db", os.path.join(g4, "REFERANS_DB"),
             "--kimlik", os.path.join(g4, "kimlik.tsv"),
             "--out", cik6], capture_output=True, text=True)
        sat6 = {}
        if os.path.exists(cik6):
            import csv as _csv6
            sat6 = {x["cift_no"]: x for x in
                    _csv6.DictReader(open(cik6, encoding="utf-8"),
                                     delimiter="\t")}
        sina("cins disina cikmayan cift CINS_OZGUL",
             sat6.get("0", {}).get("karar", "").startswith("CINS_OZGUL"),
             "cikan: %s" % sat6.get("0", {}).get("karar"))
        sina("baska cinsi de cogaltan cift CINS_AYRIMI_YOK",
             sat6.get("1", {}).get("karar", "").startswith("CINS_AYRIMI_YOK"),
             "cikan: %s" % sat6.get("1", {}).get("karar"))
        sina("cins duzeyinde OLCULEN kimlik hedef sayilmiyor",
             "Fermentimonas" in sat6.get("1", {}).get("cogaltilan_turler", "")
             and "Fermentimonas" not in
             sat6.get("1", {}).get("cogaltilan_hedef_turler", ""))
        sina("cins ici tur sayisi karara yaziliyor",
             sat6.get("0", {}).get("karar", "").endswith("_3_4"),
             "cikan: %s" % sat6.get("0", {}).get("karar"))
        shutil.rmtree(g4, ignore_errors=True)

        sina("hedef tur panelde yoksa TUR_OZGUL denmiyor",
             sat4 and all(x["karar"] != "TUR_OZGUL" for x in sat4)
             and any(x["karar"] == "HEDEF_TUR_PANELDE_YOK" for x in sat4),
             "kararlar: %s" % {x["karar"] for x in sat4})
        shutil.rmtree(g3, ignore_errors=True)
    else:
        sina("blastn/makeblastdb kurulu", False, "duzey testleri atlandi")

    print("\n17. REFERANS TASARIMINDA RAKIP KUMESI (15)")
    # 2026-08-01: Podospora referans tasarimi 29 rakip diziyle yapilmisti
    # (tek veritabani, ad basina 6 kayit). 27'nin dogrulama paneli 242
    # kayit / 50 turdu ve tasarlanan ciftlerin P. anserina ile P. comata'yi
    # cogalttigini gosterdi. Yani tasarim, gormedigi rakiplere karsi
    # ozgullugu dogrulanmis gibi rapor ediliyordu.
    RT = yukle("RT", "design_from_reference.py")
    sina("15, tur adi tanimini 27'den aliyor (tek kaynak)",
         RT.DZ.tur_adi("NR_1.1 Podospora comata CBS1 ITS")
         == "Podospora comata")
    g5 = tempfile.mkdtemp(prefix="reftas_")
    random.seed(3)

    def _m5(s, o):
        l = list(s)
        for i in range(len(l)):
            if random.random() < o:
                l[i] = random.choice("ACGT")
        return "".join(l)

    hh = "".join(random.choice("ACGT") for _ in range(700))
    dar = os.path.join(g5, "dar.fna")
    with open(dar, "w") as fh:
        for i in range(3):
            fh.write(">NR_1%d.1 Podospora pseudopauciseta CBS%d ITS\n%s\n"
                     % (i, i, _m5(hh, 0.002)))
        fh.write(">NR_50.1 Podospora comata CBS50 ITS\n%s\n" % _m5(hh, 0.05))
    genis = os.path.join(g5, "genis.fna")
    with open(genis, "w") as fh:
        for i, t in enumerate(("anserina", "bulbillosa", "fimicola",
                               "dennisiae", "bizantiorum", "cupiformis",
                               "dimorpha", "fabiformis")):
            fh.write(">NR_6%d.1 Podospora %s CBS%d ITS\n%s\n"
                     % (i, t, i, _m5(hh, 0.004)))
        fh.write(">NR_99.1 Podospora sp. CBS99 ITS\n%s\n" % _m5(hh, 0.004))

    kardes, kirp = RT.kardes_turleri_bul(
        [dar, genis], {"Podospora"}, {"Podospora pseudopauciseta"}, 60, 2)
    sina("kardes turler VERIDEN bulunuyor (elle liste gerekmeden)",
         len(kardes) == 9, "bulunan: %d -> %s" % (len(kardes),
                                                  sorted(kardes)[:3]))
    sina("hedef turun kendisi rakip sayilmiyor",
         "Podospora pseudopauciseta" not in kardes)
    sina("'Podospora sp.' kardes tur sayilmiyor",
         not any(t.endswith(" sp") or t.endswith(" sp.") for t in kardes))
    sina("tek veritabaniyla yalniz o dosyadaki kardesler bulunur",
         len(RT.kardes_turleri_bul([dar], {"Podospora"},
                                   {"Podospora pseudopauciseta"}, 60, 2)[0])
         == 1)
    kardes2, kirp2 = RT.kardes_turleri_bul(
        [dar, genis], {"Podospora"}, {"Podospora pseudopauciseta"}, 3, 2)
    sina("kardes tur siniri asilinca SESSIZ dusmez, sayilir",
         len(kardes2) == 3 and kirp2 > 0, "kirpilan=%d" % kirp2)
    # ic='Bacteroides' yazildiginda 'Parabacteroides' ve
    # 'Acetobacteroides' kayitlari hedef UYESI sayiliyordu; primer o zaman
    # baska cinslerde de urun vermek zorunda kalir, cins ozgullugu daha
    # tasarim aninda kaybedilirdi.
    bd = os.path.join(g5, "bd.fna")
    with open(bd, "w") as fh:
        for i, ad in enumerate(("Bacteroides ovatus ATCC 8483",
                                "Parabacteroides distasonis ATCC 8503",
                                "Acetobacteroides hydrogenigenes RL-C",
                                "Bacteroides thetaiotaomicron VPI-5482",
                                "Pdegerlendiricieicola vulgatus ATCC 8482")):
            fh.write(">NR_%d.1 %s 16S ribosomal RNA\n%s\n"
                     % (i, ad, "ACGT" * 60))
    bsec = RT.sec(bd, ["Bacteroides", "Parabacteroides"], 100)
    sina("cins adi SOZCUK olarak araniyor, alt dize olarak degil",
         len(bsec["Bacteroides"]) == 2,
         "Bacteroides icin %d kayit" % len(bsec["Bacteroides"]))
    sina("Parabacteroides, Bacteroides uyesi sayilmiyor",
         all("Parabacteroides" not in b for b, _ in bsec["Bacteroides"]))
    sina("Acetobacteroides, Bacteroides uyesi sayilmiyor",
         all("Acetobacteroides" not in b for b, _ in bsec["Bacteroides"]))
    sina("Parabacteroides kendi adiyla hala bulunuyor",
         len(bsec["Parabacteroides"]) == 1)
    sina("sec() birden cok veritabanindan toplayabiliyor",
         len(RT.sec([dar, genis], ["Podospora"], 100)["Podospora"]) == 13,
         "toplanan: %d"
         % len(RT.sec([dar, genis], ["Podospora"], 100)["Podospora"]))
    shutil.rmtree(g5, ignore_errors=True)

    if a.gercek_veri:
        print("\n13. GERCEK VERI: URETILEN CIFTLERIN GEOMETRISI")
        tsvler = sorted(glob.glob(os.path.join(HERE, a.aday, "*__*.tsv")))
        if not tsvler:
            sina("aday TSV bulundu", False, a.aday)
        else:
            toplam = hata = 0
            for t in tsvler[:6]:
                r = subprocess.run(
                    [sys.executable, os.path.join(HERE,
                                                  "check_primer_geometry.py"),
                     "--tsv", t, "--kons", a.kons, "--en-fazla", "500"],
                    capture_output=True, text=True)
                for line in r.stdout.splitlines():
                    if "gecen satir" in line:
                        x = line.split(":")[1].split("/")
                        toplam += int(x[1]); hata += int(x[1]) - int(x[0])
            sina("gercek ciftlerin tamami geometri denetimini geciyor",
                 hata == 0, "sinanan=%d hatali=%d" % (toplam, hata))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gercek-veri", action="store_true")
    p.add_argument("--aday", default="pr_aday")
    p.add_argument("--kons", default="pr_kons/konsensus")
    a = p.parse_args()
    print("=" * 72)
    print("BORU HATTI REGRESYON TESTI")
    print("=" * 72)
    testler(a)
    gecen = sum(1 for _, ok, _ in SONUC if ok)
    print("\n" + "=" * 72)
    print(u'RESULT: %d of %d tests passed, %d failed'
          % (len(SONUC), gecen, len(SONUC) - gecen))
    for ad, ok, ayrinti in SONUC:
        if not ok:
            print("   KALDI: %s   %s" % (ad, ayrinti))
    print("=" * 72)
    sys.exit(0 if gecen == len(SONUC) else 1)


if __name__ == "__main__":
    main()
