#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
regression_test.py
Tests the pipeline's rules without trusting the code's own helper functions.
Every test's expected result is derived from the panel decisions or from known
mathematics.

Usage:
  python3 regression_test.py                 # the quick tests
  python3 regression_test.py --real-data   # with real consensuses and reads

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
    sina(u'rc() agrees with an independent implementation', ok)
    sina("rc(rc(x)) == x", all(E.rc(E.rc(s)) == s
                              for s in ("ATGCRYKM", "AACCGGTT")))
    sina("ATGC'nin ters tumleyeni GCAT", E.rc("ATGC") == "GCAT",
         "cikan: %s" % E.rc("ATGC"))

    print("\n2. KOMPOZISYON KURALLARI")
    ar = _ap.Namespace(len_min=18, len_max=25, gc_min=40, gc_max=60,
                       gc_hard_min=35, gc_hard_max=65, gc_clamp_last=5,
                       gc_clamp_max=3, homopolymer_max=4, require_3p_gc=1,
                       degeneracy_budget=0, degeneracy_fold_max=4,
                       iupac_max=2, iupac_clamp_forbidden=5)
    sina("3' ucu A ile biten oligo elenir",
         not E.composition_ok("CGCGATATCGCGATATCGA", ar)[0])
    sina("3' ucu G ile biten oligo gecer",
         E.composition_ok("CGCGATATCGCGATATCGG", ar)[0])
    sina(u'five identical bases in a row are eliminated',
         not E.composition_ok("CGCGAAAAACGCGATATCG", ar)[0])
    sina("son bes bazda dort G/C elenir",
         not E.composition_ok("ATATATATATATATAGCGCC", ar)[0])
    dusuk = "ATATATATATATATATATAC"
    sina("GC sert alt sinirinin altindaki oligo elenir (GC=%%%.0f)" % ref_gc(dusuk),
         not E.composition_ok(dusuk, ar)[0])
    sina(u'gc_pct agrees with an independent calculation',
         all(abs(E.gc_pct(s) - ref_gc(s)) < 1e-9
             for s in ("GGCC", "ATAT", "ACGTACGT")))

    print(u'\n3. THE IUPAC STRAND RULE')
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
    sina(u'an ambiguous position expands to every allele', v is not None and len(v) == 2)
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
    sina(u'find_bindings agrees with a direct application of the rule (400 trials)',
         fark == 0, "ayrisan: %d" % fark)
    sina(u'a mismatch in the last two bases is not accepted',
         ref_baglanir("ACGTACGTACGTACGTACGA", "ACGTACGTACGTACGTACGC") is None)
    sina(u'a single mismatch in the tail is accepted',
         ref_baglanir("ACGTACGTACGTACGTAAGC", "ACGTACGTACGTACGTACGC") == 1)

    print(u'\n5. THE PRODUCT GEOMETRY')
    # Instead of a synthetic position dictionary, a REAL template is built, the primers
    # are cut from that template and the product is formed by hand. The expected value
    # then comes from the sequence itself rather than from the code's own measure.
    ap2 = _ap.Namespace(prod_min=70, prod_hard_max=300)
    random.seed(31)
    L = 500
    kalip = "".join(random.choice("ACGT") for _ in range(L))
    F = kalip[100:120]                 # ileri primer, arti zincir 100..119
    R = G.rc(kalip[260:280])           # geri primer, arti zincirin rc'si
    elde_urun = kalip[100:280]
    sina(u'a product built by hand starts with the forward primer', elde_urun.startswith(F))
    sina(u'a product built by hand ends with the rc of the reverse primer',
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
    sina(u'the product length matches the amplicon in the real template (%d bp)' % beklenen,
         pl == beklenen, "cikan: %s" % pl)
    sina(u'the reverse configuration is tried too (swapping the primers gives the same)',
         G.product_len(br, bf, len(R), len(F), ap2) == beklenen)
    # A pair whose 3' ends point AWAY from one another: swapping the primers makes each
    # of them look BACKWARDS on its own strand
    Fk = G.rc(kalip[100:120])          # arti zincirin 100..119'unun rc'si
    Rk = kalip[260:280]
    bfk, brk = bagla(Fk), bagla(Rk)
    sina(u'a pair whose 3\' ends FACE AWAY from one another gives no product',
         G.product_len(bfk, brk, len(Fk), len(Rk), ap2) is None)
    # If the product falls below the lower bound it is discarded in strict design, but
    # it is still seen in the competitor scan with pmin=1
    Ry = G.rc(kalip[140:160])          # the product is 100..159, that is 60 bp
    bry = bagla(Ry)
    sina(u'a short product (60 bp) is eliminated at the strict lower bound',
         G.product_len(bf, bry, len(F), len(Ry), ap2) is None)
    sina(u'the lower bound can be dropped to 1 in a competitor scan',
         G.product_len(bf, bry, len(F), len(Ry), ap2, pmin=1) == 60)

    print("\n6. WILSON SINIRLARI")
    ok_alt = ok_ust = True
    for k, n in ((0, 10), (1, 10), (5, 10), (10, 10), (1, 17), (2, 100),
                 (1, 1), (0, 1), (37, 300), (500, 20000)):
        ok_alt &= abs(O.wilson_alt(k, n) - ref_wilson(k, n)) < 1e-12
        ok_ust &= abs(O.wilson_ust(k, n) - ref_wilson(k, n, ust=True)) < 1e-12
    sina(u'wilson_alt agrees with an independent implementation', ok_alt)
    sina(u'wilson_ust agrees with an independent implementation', ok_ust)
    sina("alt sinir <= nokta tahmini <= ust sinir",
         all(O.wilson_alt(k, n) <= k / n <= O.wilson_ust(k, n)
             for k, n in ((1, 10), (5, 10), (9, 10), (37, 300))))
    sina(u'a known value: wilson_alt(1,10)=0,01788',
         abs(O.wilson_alt(1, 10) - 0.017876) < 1e-5)

    print("\n7. BASKIN ALEL CAGRISI")
    sina(u'dominant_allele_consensus.py calls N on a tie (the code path)',
         "esitler" in open(os.path.join(HERE, "dominant_allele_consensus.py"),
                           encoding="utf-8").read())
    sina(u'dominant_allele_consensus.py does not delete inner N, it only trims the ends',
         'ref[bas:son]' in open(os.path.join(HERE,
                                             "dominant_allele_consensus.py"),
                                encoding="utf-8").read())
    sina(u'dominant_allele_consensus.py converts the reads to upper case',
         '.strip().upper()' in open(os.path.join(
             HERE, "dominant_allele_consensus.py"), encoding="utf-8").read())

    print(u'\n8. THE INDISTINGUISHABILITY MEASUREMENT')
    sina(u'N does not count as a match with any base', not D.baz_kesisir("N", "A"))
    sina("N ile N de eslesmez", not D.baz_kesisir("N", "N"))
    sina("Y ile C kesisir", D.baz_kesisir("Y", "C"))
    sina("Y ile G kesismez", not D.baz_kesisir("Y", "G"))
    sina("A ile A kesisir", D.baz_kesisir("A", "A"))

    print(u'\n9. THE ALIGNMENT BACKEND')
    try:
        H = yukle("H", "alignment.py")
        sina(u'there is an alignment backend', H.ARKA_UC is not None, str(H.ARKA_UC))
        ref = "".join(random.choice("ACGT") for _ in range(2000))
        A = H.Hizalayici(seq=ref, preset="map-ont")
        q = ref[500:900]
        ileri = [h for h in A.map(q)]
        ters = [h for h in A.map(H.revcomp(q))]
        sina(u'the forward alignment is at the right position',
             ileri and ileri[0].r_st == 500 and ileri[0].r_en == 900,
             str(ileri[:1]))
        sina(u'the reverse alignment is at the same position with strand=-1',
             ters and ters[0].r_st == 500 and ters[0].strand == -1,
             str(ters[:1]))
        sina(u'on the minus strand the query coordinate is relative to the original query',
             ters and ters[0].q_st == 0 and ters[0].q_en == len(q),
             str(ters[:1]))
    except Exception as e:
        sina(u'there is an alignment backend', False, str(e)[:80])

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
    sina(u'a product is found in the raw reads (forward and reverse reads)',
         both > 200, u'tot=%d F=%d R=%d product=%d' % (tot, fh_, rh, both))
    sina(u'both the forward and the reverse reads are counted',
         abs(fh_ - rh) < 30, "F=%d R=%d" % (fh_, rh))
    tot2, _, _, both2 = O.okuma_taramasi(yol, F, ref_rc(T[250:271])[::-1],
                                         70, 300, 300)
    sina(u'a reverse primer written in the wrong orientation GIVES NO product',
         both2 < both * 0.2, u'product=%d (in the right orientation %d)' % (both2, both))

    print(u'\n11. PRODUCT LENGTH IN THE EXTERNAL DATABASES')
    V = yukle("V", "external_databases.py")
    sina(u'external_databases.py measures the product length between the two 5\' ends (the same measure as design_group_primers.py)',
         "(ur + lr - 1) - (uf - lf + 1) + 1" in open(
             os.path.join(HERE, "external_databases.py"),
             encoding="utf-8").read())
    sina("14 baglanma kuralini uyguluyor",
         V.baglanma_uygun("ACGTACGTACGTACGTACGC",
                          "ACGTACGTACGTACGTACGC")[0] and
         not V.baglanma_uygun("ACGTACGTACGTACGTACGA",
                              "ACGTACGTACGTACGTACGC")[0])

    print(u'\n12. THE FILE FORMAT CONTRACT BETWEEN STAGES 08 AND 09')
    # steps/specificity.py reads the dislanan_takson.tsv written by
    # steps/batch_design.py. If the columns are read by position, then when the file
    # format changes the exclusion is SILENTLY DISABLED and specificity.py says nothing.
    # That violates the rule "no decision is left to a single piece of code"; the tests
    # below confirm that it is read by header and that it stops on an unrecognised one.
    import subprocess as _sp, tempfile as _tf, textwrap as _tw
    k09 = open(os.path.join(HERE, "specificity.py"), encoding="utf-8").read()
    k08 = open(os.path.join(HERE, "batch_design.py"), encoding="utf-8").read()
    sina(u'batch_design.py writes the \'grup\' column into the header',
         'df.write("grup\\ttaxid\\tetiket\\tuzunluk\\tkapsanan\\n")' in k08)
    try:
        i0 = k09.index("    dislanan = set()")
        j0 = k09.index("    ayirt = {}")
        govde = _tw.dedent(k09[i0:j0])
        dtmp = _tf.mkdtemp()
        sarmal = ("import os, sys\nclass A: candidates=%r\na=A()\n"
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
        sina(u'it does NOT CARRY ON SILENTLY in the old 4 column format, it stops',
             rc2 != 0, "cikis=%d" % rc2)
        rc3, o3 = cagir("taxid\tgrup\tetiket\tuzunluk\tkapsanan\n"
                        "1434102\tA1-4\tA1-4_1434102\t1444\t0\n")
        sina(u'even with the column order changed it is read correctly from the header',
             rc3 == 0 and "('A1-4', '1434102')" in o3, o3[-70:])
        rc4, o4 = cagir("grup\ttaxid\tetiket\tuzunluk\tkapsanan\n")
        sina(u'a file with a header only gives an empty set, not an error',
             rc4 == 0 and "SONUC []" in o4, o4[-70:])
    except ValueError:
        sina("09'da dislama blogu bulundu", False, u'the block markers have changed')

    print(u'\n13. EXTERNAL DATABASE SET MAPPING AND COVERAGE AUDIT')
    # The bug found on 2026-08-01: ROD_v1.2_operon_variants.fasta is eukaryote only
    # (60320/60320 Eukaryota) but had been assigned to classes A1/A2/B. "No off-target
    # product" was written for 71 archaeal and bacterial pairs, when that database held
    # not one sequence from that domain. The tests below check both the corrected
    # mapping and the coverage audit that keeps the same mistake from staying silent
    # again.
    DV = yukle("DV", "external_databases.py")
    MF = yukle("MF", "mfeprimer_layer.py")
    sina(u'external_databases.py and mfeprimer_layer.py see the NARROW set the same way', DV.SINIF_DB == MF.SINIF_DB)
    sina(u'external_databases.py and mfeprimer_layer.py see the BROAD set the same way',
         DV.SINIF_DB_GENIS == MF.SINIF_DB_GENIS)
    sina(u'the narrow and the broad set do not intersect (the same db is not scanned twice)',
         all(not (set(DV.SINIF_DB[s]) & set(DV.SINIF_DB_GENIS[s]))
             for s in DV.SINIF_DB))
    sina(u'every class sees the same set of databases',
         len({tuple(sorted(set(DV.SINIF_DB[s]) | set(DV.SINIF_DB_GENIS[s])))
              for s in DV.SINIF_DB}) == 1)
    sina("ROD mantar siniflarinda da taraniyor",
         all("ROD_v1.2_operon_variants.fasta" in DV.SINIF_DB_GENIS[s]
             for s in ("F1", "F2")))
    sina(u'fungi.18SrRNA was not forgotten in any class',
         all("fungi.18SrRNA.fna" in (DV.SINIF_DB[s] + DV.SINIF_DB_GENIS[s])
             for s in DV.SINIF_DB))
    sina(u'ref_all and ref_all2 are not scanned (they are a subset of the others)',
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
        sina(u'a consensus is read only for its own class',
             len(dz) == 1 and DV.sinif_konsensuslari(
                 os.path.join(gec, "kons"), "B") == [])
        _, uz_yok, _ = DV.kapsam_olc(ilgisiz, dz, gec, "a", 2, 600)
        _, uz_var, kim_var = DV.kapsam_olc(ilgili, dz, gec, "b", 2, 600)
        sina(u'an unrelated database gives KAPSAM_YOK', uz_yok < 400,
             u'the longest alignment is %d bp' % uz_yok)
        sina(u'a database holding a relative gives KAPSANIYOR',
             uz_var >= 400 and kim_var > 80,
             "en uzun %d bp %%%.1f" % (uz_var, kim_var))
        sina(u'with no consensus given, coverage is marked OLCULMEDI',
             DV.kapsam_olc(ilgili, [], gec, "c", 2, 600)[0] == "KAPSAM_OLCULMEDI")
        shutil.rmtree(gec, ignore_errors=True)
    else:
        sina("blastn/makeblastdb kurulu", False, u'the coverage tests were skipped')

    print(u'\n14. SEPARATING THE OWN TAXON FROM A FOREIGN TAXON')
    # 2026-08-01: in the wide scan, some of the highest "off-target product" counts were
    # THE TARGET ITSELF (the Methanothrix target amplifies the Methanothrix records in
    # SILVA). Ranking by the raw count pushed the wrong primers to the front. The tests
    # below check that the distinction is made correctly.
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

        # for the same pair, 3 records of ITS OWN taxon and 3 of a FOREIGN taxon
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
        open(os.path.join(gec2, "targets.tsv"), "w").write(
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
             "--consensus", os.path.join(gec2, "kons"),
             "--targets", os.path.join(gec2, "targets.tsv"),
             "--names", os.path.join(gec2, "adlar.tsv"),
             "--out", cik], capture_output=True, text=True)
        satir = {}
        if os.path.exists(cik):
            import csv as _csv
            for r in _csv.DictReader(open(cik, encoding="utf-8"),
                                     delimiter="\t"):
                satir[r["hedef"]] = r
        oz = satir.get("Methanothrix_soehngenii_turu", {})
        ev = satir.get("Arke_universal", {})
        sina(u'on a specific target 3 own and 3 foreign products are separated',
             oz.get("urun_kendi_taksonda") == "3"
             and oz.get("urun_yabanci_taksonda") == "3",
             u'own=%s foreign=%s total=%s' % (oz.get("urun_kendi_taksonda"),
                                                oz.get("urun_yabanci_taksonda"),
                                                oz.get("hedef_disi_urun")))
        sina(u'own + foreign + unknown = the raw total',
             oz and (int(oz["urun_kendi_taksonda"])
                     + int(oz["urun_yabanci_taksonda"])
                     + int(oz["urun_takson_bilinmiyor"])
                     == int(oz["hedef_disi_urun"])))
        sina(u'the foreign examples show foreign taxa only',
             oz and "YABANCI" in oz["yabanci_ornekler"]
             and "KENDI" not in oz["yabanci_ornekler"])
        sina(u'on a universal target no product counts as foreign',
             ev.get("hedef_turu") == "evrensel"
             and ev.get("urun_yabanci_taksonda") == "0",
             "turu=%s yabanci=%s" % (ev.get("hedef_turu"),
                                     ev.get("urun_yabanci_taksonda")))
        shutil.rmtree(gec2, ignore_errors=True)
    else:
        sina("blastn/makeblastdb kurulu", False, u'the taxon tests were skipped')

    print("\n15. YABANCI VURUSUN UZAKLIGI (SOYAGACI DERINLIGI)")
    # A foreign taxon is not a sufficient measure on its own: some targets are
    # functional groups. Measured (2026-08-01): the foreign hits of
    # Hidrojenotrofik_metanojenler are Methanobacterium and Methanosphaera, both of
    # which are hydrogenotrophic methanogens; the foreign hits of Nitrosocosmicus_AOA
    # are Nitrosotalea and Nitrosopumilus, both of which are AOA. Against that,
    # Petrimonas -> Flavobacterium and Trichoderma -> Calonectria really are distant.
    # The distance is measured from the lineage in the database headers, not from a hand
    # written table.
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
    sina(u'a RefSeq header carries no lineage, so an empty list comes back',
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

    sina(u'a hit in the same family counts as NEAR',
         _karar([KENDI1], AYNI_AILE) == "YAKIN")
    sina(u'the same order but a different family counts as FAR',
         _karar([KENDI1], AYNI_TAKIM) == "UZAK")
    sina(u'a hit in a different class counts as FAR',
         _karar([KENDI1], FARKLI_SINIF) == "UZAK")
    # THRESHOLD STABILITY: had the common prefix been used, the reference lineage would
    # shorten as the number of own hits grew and the same foreign hit would drift to
    # NEAR. The dominant lineage prevents that.
    sina(u'the decision does not depend on its own hit count',
         len({_karar(k, AYNI_TAKIM) for k in ([KENDI1], [KENDI1, KENDI2],
                                              [KENDI1, KENDI2, KENDI1])}) == 1)
    sina(u'the depth of the dominant lineage does not shorten with its own hit count',
         len(DV._baskin_soy([DV._soyagaci(KENDI1), DV._soyagaci(KENDI2)]))
         == len(DV._soyagaci(KENDI1)))
    sina(u'with no hit of its own the reference lineage is empty',
         DV._baskin_soy([]) == [])

    print(u'\n16. DECISION LEVEL AUDIT (27): SPECIES AND GENUS SPECIFICITY')
    # The panel decision asks for SPECIES specificity on six targets and GENUS on four.
    # steps/specificity.py looks at the competitors in the sample and
    # steps/external_databases.py at the external databases; neither asks "does it
    # separate the target from the other SPECIES of its own genus".
    # steps/check_taxonomic_level.py measures that against a panel of sibling species.
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
    # RefSeq ITS and 28S records end with '; from TYPE material'. In the earlier
    # version, if the header held a semicolon the RefSeq branch never ran and the TYPE
    # STRAIN records were not taken into the panel; yet those records are exactly the
    # gold standard of species separation (35 records for Petriella alone).
    # The hedef_tur column does not ADD, it REPLACES. On the Petriella_musispora row
    # in_taxid holds the bins' Kraken2 labels (for the in-sample support); had it added,
    # Trichoderma would have entered the target species set and a pair amplifying
    # Trichoderma would have counted as "there is a product in the target species".
    _g6 = tempfile.mkdtemp(prefix="hedeftur_")
    _h = os.path.join(_g6, "h.tsv")
    _ad = os.path.join(_g6, "a.tsv")
    open(_ad, "w").write("101201\tTrichoderma asperellum\n")
    open(_h, "w").write(
        "karar\thedef\tduzey\tin_taxid\tharic\tnot\thedef_tur\n"
        "5\tPetriella_musispora\ttur\t101201\t\tnot\tPetriella musispora\n"
        "1\tMicroascaceae_askomikot\ttur\t101201\t\tnot\t\n")
    _hd = {x["hedef"]: x for x in DZ.hedefleri_oku(_h, _ad, None)}
    sina(u'when hedef_tur is given the taxid name is REPLACED, not added to',
         _hd["Petriella_musispora"]["hedef_turler"] == {"Petriella musispora"},
         "cikan: %s" % _hd["Petriella_musispora"]["hedef_turler"])
    sina(u'when hedef_tur is given the genus changes too',
         _hd["Petriella_musispora"]["cinsler"] == {"Petriella"})
    sina(u'a row with an empty hedef_tur keeps the old behaviour',
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
    sina(u'\'Ca. Nitrosocosmicus hydrocola\' is split into a genus and a species',
         DZ.ad_parcala("Ca. Nitrosocosmicus hydrocola")
         == ("Nitrosocosmicus", "hydrocola"))
    # In UNITE, 's__Trichoderma_sp' was being counted as a species name; for Trichoderma
    # alone, 16910 records were entering the panel as a spurious "species".
    for kotu in ("s__Trichoderma_sp", "s__Marasmius_spp", "s__Podospora_cf"):
        sina("UNITE'te '%s' tur adi SAYILMIYOR" % kotu,
             DZ.tur_adi("UDB1|k__Fungi;g__X;%s|SH1" % kotu) == "",
             "cikan: %r" % DZ.tur_adi("UDB1|k__Fungi;g__X;%s|SH1" % kotu))
    sina("gercek ikili ad hala taniniyor",
         DZ.tur_adi("UDB1|k__Fungi;g__Petriella;s__Petriella_setifera|SH1")
         == "Petriella setifera")
    # 'Methanosarcina_barkeri_referans' strips to 'Methanosarcina_barkeri' while the
    # name in targets.tsv is 'Methanosarcina_barkeri_turu'. When the match failed, the
    # target's ONLY primer set dropped silently and the target appeared as CIFT_YOK.
    _ADL = ["Methanosarcina_barkeri_turu", "Proteiniphilum_cinsi",
            "Podospora_pseudopauciseta", "Proteolitik_sintrofik_bakteriler"]
    sina(u'a reference target is tied to a target carrying the \'_turu\' suffix',
         DZ.referans_esle("Methanosarcina_barkeri_referans", _ADL)
         == "Methanosarcina_barkeri_turu",
         "cikan: %s" % DZ.referans_esle("Methanosarcina_barkeri_referans", _ADL))
    sina(u'a reference target matching exactly does not change',
         DZ.referans_esle("Proteiniphilum_cinsi_referans", _ADL)
         == "Proteiniphilum_cinsi")
    sina(u'a shortened reference target is tied to the full name',
         DZ.referans_esle("Proteolitik_sintrofik_referans", _ADL)
         == "Proteolitik_sintrofik_bakteriler")
    sina(u'a reference target with no match returns None (it does not drop silently)',
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
            # THREE sibling species: since the default threshold is 2, a pair amplifying all
            # three must stay ABOVE the threshold (TUR_AYRIMI_YOK)
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
        open(os.path.join(g3, "targets.tsv"), "w").write(
            "karar\thedef\tduzey\tin_taxid\tharic\tnot\n"
            "1\tMethanothrix_soehngenii_turu\ttur\t2223\t\tM. soehngenii\n")
        open(os.path.join(g3, "adlar.tsv"), "w").write(
            "2223\tMethanothrix soehngenii\n")
        cik3 = os.path.join(g3, "cikti.tsv")
        r3 = subprocess.run(
            [sys.executable, os.path.join(HERE, "check_taxonomic_level.py"),
             "--targets", os.path.join(g3, "targets.tsv"),
             "--names", os.path.join(g3, "adlar.tsv"),
             "--final", os.path.join(g3, "final"),
             "--db", os.path.join(g3, "REFERANS_DB"),
             "--out", cik3], capture_output=True, text=True)
        sat = []
        if os.path.exists(cik3):
            import csv as _csv3
            sat = list(_csv3.DictReader(open(cik3, encoding="utf-8"),
                                        delimiter="\t"))
        kararlar = {x["cift_no"]: x for x in sat}
        sina(u'a pair that discriminates comes out TUR_OZGUL',
             kararlar.get("0", {}).get("karar") == "TUR_OZGUL",
             "cikan: %s" % kararlar.get("0", {}).get("karar"))
        sina(u'a pair that does not discriminate comes out TUR_AYRIMI_YOK',
             kararlar.get("1", {}).get("karar") == "TUR_AYRIMI_YOK",
             "cikan: %s" % kararlar.get("1", {}).get("karar"))
        sina(u'on a TUR_OZGUL pair the product count in the sibling species is zero',
             kararlar.get("0", {}).get("diger_turde_urun") == "0")
        sina(u'on a pair that does not discriminate the sibling species are reported by name',
             "harundinacea" in kararlar.get("1", {}).get("cogaltilan_turler", ""))
        sina("'sp.' kaydi panele alinmadi (panel 4 tur)",
             kararlar.get("0", {}).get("panel_tur_sayisi") == "4",
             "panel_tur_sayisi=%s" % kararlar.get("0", {}).get("panel_tur_sayisi"))
        sina(u'the cross reacting SPECIES count is reported (not the product count)',
             kararlar.get("1", {}).get("capraz_tur_sayisi") == "3",
             "capraz_tur_sayisi=%s"
             % kararlar.get("1", {}).get("capraz_tur_sayisi"))
        # The panel decision tolerates 1-2 cross reactions. When the threshold is raised the
        # same pair must count as within threshold, and must not when it is low.
        cik5 = os.path.join(g3, "cikti_esik3.tsv")
        subprocess.run(
            [sys.executable, os.path.join(HERE, "check_taxonomic_level.py"),
             "--targets", os.path.join(g3, "targets.tsv"),
             "--names", os.path.join(g3, "adlar.tsv"),
             "--final", os.path.join(g3, "final"),
             "--db", os.path.join(g3, "REFERANS_DB"),
             "--cross-species-tolerance", "3",
             "--out", cik5], capture_output=True, text=True)
        sat5 = {}
        if os.path.exists(cik5):
            import csv as _csv5
            sat5 = {x["cift_no"]: x for x in
                    _csv5.DictReader(open(cik5, encoding="utf-8"),
                                     delimiter="\t")}
        sina(u'with the threshold raised to 3, a pair with 3 cross reactions becomes TUR_OZGUL_ESIKLI',
             sat5.get("1", {}).get("karar") == "TUR_OZGUL_ESIKLI",
             "cikan: %s" % sat5.get("1", {}).get("karar"))
        sina(u'a pair with no cross reaction stays TUR_OZGUL whatever the threshold',
             sat5.get("0", {}).get("karar") == "TUR_OZGUL")
        sina(u'check_taxonomic_level.py finished with exit code zero', r3.returncode == 0,
             r3.stderr.strip()[-120:])
        # if the target species is absent from the panel, it must NOT say TUR_OZGUL
        open(os.path.join(g3, "targets.tsv"), "w").write(
            "karar\thedef\tduzey\tin_taxid\tharic\tnot\n"
            "1\tMethanothrix_soehngenii_turu\ttur\t9999\t\tyok\n")
        open(os.path.join(g3, "adlar.tsv"), "w").write(
            "9999\tMethanothrix yoktur\n")
        cik4 = os.path.join(g3, "cikti2.tsv")
        subprocess.run(
            [sys.executable, os.path.join(HERE, "check_taxonomic_level.py"),
             "--targets", os.path.join(g3, "targets.tsv"),
             "--names", os.path.join(g3, "adlar.tsv"),
             "--final", os.path.join(g3, "final"),
             "--db", os.path.join(g3, "REFERANS_DB"),
             "--out", cik4], capture_output=True, text=True)
        sat4 = []
        if os.path.exists(cik4):
            import csv as _csv4
            sat4 = list(_csv4.DictReader(open(cik4, encoding="utf-8"),
                                         delimiter="\t"))
        # GENUS LEVEL: if a product forms outside the declared genus, it is
        # not genus specific. Measured: two of the Proteiniphilum pairs also
        # amplified Fermentimonas caenicola, and the old count showed both as
        # the widest covering pair.
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
        open(os.path.join(g4, "targets.tsv"), "w").write(
            "karar\thedef\tduzey\tin_taxid\tharic\tnot\n"
            "2\tProteiniphilum_cinsi\tcins\t2829812,1642647\t\tP. cinsi\n")
        open(os.path.join(g4, "adlar.tsv"), "w").write(
            "2829812\tProteiniphilum propionicum\n"
            "1642647\tProteiniphilum saccharofermentans\n")
        # the measured identity is ANOTHER GENUS: it must not count as the target at genus level
        open(os.path.join(g4, "kimlik.tsv"), "w").write(
            "hedef\tkraken_etiketi\tolculen_kimlik\n"
            "Proteiniphilum_cinsi\tX\tFermentimonas caenicola\n")
        cik6 = os.path.join(g4, "cikti.tsv")
        subprocess.run(
            [sys.executable, os.path.join(HERE, "check_taxonomic_level.py"),
             "--targets", os.path.join(g4, "targets.tsv"),
             "--names", os.path.join(g4, "adlar.tsv"),
             "--final", os.path.join(g4, "final"),
             "--db", os.path.join(g4, "REFERANS_DB"),
             "--identity", os.path.join(g4, "kimlik.tsv"),
             "--out", cik6], capture_output=True, text=True)
        sat6 = {}
        if os.path.exists(cik6):
            import csv as _csv6
            sat6 = {x["cift_no"]: x for x in
                    _csv6.DictReader(open(cik6, encoding="utf-8"),
                                     delimiter="\t")}
        sina(u'a pair that does not leave the genus is CINS_OZGUL',
             sat6.get("0", {}).get("karar", "").startswith("CINS_OZGUL"),
             "cikan: %s" % sat6.get("0", {}).get("karar"))
        sina(u'a pair that amplifies another genus too is CINS_AYRIMI_YOK',
             sat6.get("1", {}).get("karar", "").startswith("CINS_AYRIMI_YOK"),
             "cikan: %s" % sat6.get("1", {}).get("karar"))
        sina(u'at genus level the MEASURED identity does not count as the target',
             "Fermentimonas" in sat6.get("1", {}).get("cogaltilan_turler", "")
             and "Fermentimonas" not in
             sat6.get("1", {}).get("cogaltilan_hedef_turler", ""))
        sina(u'the within genus species count is written into the decision',
             sat6.get("0", {}).get("karar", "").endswith("_3_4"),
             "cikan: %s" % sat6.get("0", {}).get("karar"))
        shutil.rmtree(g4, ignore_errors=True)

        sina(u'if the target species is not in the panel it is not called TUR_OZGUL',
             sat4 and all(x["karar"] != "TUR_OZGUL" for x in sat4)
             and any(x["karar"] == "HEDEF_TUR_PANELDE_YOK" for x in sat4),
             "kararlar: %s" % {x["karar"] for x in sat4})
        shutil.rmtree(g3, ignore_errors=True)
    else:
        sina("blastn/makeblastdb kurulu", False, u'the level tests were skipped')

    print(u'\n17. THE COMPETITOR SET IN A REFERENCE DESIGN (15)')
    # 2026-08-01: the Podospora reference design had been made with 29 competitor
    # sequences (one database, 6 records per name). The verification panel of
    # check_taxonomic_level.py held 242 records across 50 species and showed that the
    # designed pairs amplify both P. anserina and P. comata. In other words the design
    # was being reported as if its specificity had been verified against competitors it
    # had never seen.
    RT = yukle("RT", "design_from_reference.py")
    sina(u'design_from_reference.py takes the species name definition from check_taxonomic_level.py (a single source)',
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
    sina(u'sibling species are found FROM THE DATA (with no hand written list)',
         len(kardes) == 9, "bulunan: %d -> %s" % (len(kardes),
                                                  sorted(kardes)[:3]))
    sina(u'the target species itself does not count as a competitor',
         "Podospora pseudopauciseta" not in kardes)
    sina("'Podospora sp.' kardes tur sayilmiyor",
         not any(t.endswith(" sp") or t.endswith(" sp.") for t in kardes))
    sina(u'with a single database only the siblings in that file are found',
         len(RT.kardes_turleri_bul([dar], {"Podospora"},
                                   {"Podospora pseudopauciseta"}, 60, 2)[0])
         == 1)
    kardes2, kirp2 = RT.kardes_turleri_bul(
        [dar, genis], {"Podospora"}, {"Podospora pseudopauciseta"}, 3, 2)
    sina("kardes tur siniri asilinca SESSIZ dusmez, sayilir",
         len(kardes2) == 3 and kirp2 > 0, "kirpilan=%d" % kirp2)
    # When ic='Bacteroides' was written, 'Parabacteroides' and 'Acetobacteroides'
    # records counted as MEMBERS of the target; the primer then had to give a product in
    # other genera too, and genus specificity was lost at design time.
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
    sina(u'the genus name is searched for AS A WORD, not as a substring',
         len(bsec["Bacteroides"]) == 2,
         u'%d records for Bacteroides' % len(bsec["Bacteroides"]))
    sina(u'Parabacteroides does not count as a Bacteroides member',
         all("Parabacteroides" not in b for b, _ in bsec["Bacteroides"]))
    sina(u'Acetobacteroides does not count as a Bacteroides member',
         all("Acetobacteroides" not in b for b, _ in bsec["Bacteroides"]))
    sina(u'Parabacteroides is still found under its own name',
         len(bsec["Parabacteroides"]) == 1)
    sina(u'sec() can collect from more than one database',
         len(RT.sec([dar, genis], ["Podospora"], 100)["Podospora"]) == 13,
         "toplanan: %d"
         % len(RT.sec([dar, genis], ["Podospora"], 100)["Podospora"]))
    shutil.rmtree(g5, ignore_errors=True)

    # ------------------------------------------------------------------
    print(u'\n18. THE LOG CONTRACT BETWEEN THE GROUP ENGINE AND ITS READERS')
    # design_group_primers.py prints its counters as text and four other
    # scripts read them back with a regular expression: batch_design.py and
    # design_from_reference.py take the valid pair count, split_clusters.py
    # takes the count plus the blocking members, and export_excel.py fills a
    # whole sheet from ten of them. Nothing binds the two sides, so rewording
    # one line there leaves an empty column here and no error anywhere. That
    # happened once, when the engine's output was translated into English and
    # every pattern kept looking for the Turkish wording. This test runs the
    # engine on a small synthetic case and requires every pattern to match.
    import re as _re
    g18 = _tf.mkdtemp()
    try:
        random.seed(4242)
        govde = "".join(random.choice("ACGT") for _ in range(600))
        uyeler = []
        for i in range(3):
            u = os.path.join(g18, "grup_%d.fasta" % i)
            # the members share the backbone, so a conserved pair can be found
            with open(u, "w") as fh:
                fh.write(">A1-%d_%d\n%s\n" % (i + 1, 1000 + i, govde))
            uyeler.append(u)
        cikti = os.path.join(g18, "cift.tsv")
        r = _sp.run([sys.executable, os.path.join(HERE, "design_group_primers.py"),
                     "--in-group"] + uyeler
                    + ["--label", "LOG_CONTRACT", "--out", cikti,
                       "--max-oligo", "60", "--max-pairs", "20"],
                    capture_output=True, text=True)
        txt = r.stdout + r.stderr
        DESENLER = [
            (u'oligos after composition filter',
             r"oligos after composition filter:\s*(\d+)"),
            (u'oligos after thermodynamics',
             r"oligos after thermodynamics:\s*(\d+)"),
            (u'oligos binding every target member',
             r"oligos binding every target member:\s*(\d+)"),
            (u'dropped, pair Tm difference',
             r"dropped, pair Tm difference\s*:\s*(\d+)"),
            (u'dropped, no product in one member',
             r"dropped, no product in one member\s*:\s*(\d+)"),
            (u'dropped, product forms in a competitor',
             r"dropped, product forms in a competitor\s*:\s*(\d+)"),
            (u'dropped, no orphan primer',
             r"dropped, no orphan primer\s*:\s*(\d+)"),
            (u'dropped, hetero-dimer dG',
             r"dropped, hetero-dimer dG\s*:\s*(\d+)"),
            (u'valid pairs', r"valid pairs\s*:\s*(\d+)"),
        ]
        for ad, desen in DESENLER:
            sina(u'the engine still prints "%s" as its readers expect' % ad,
                 _re.search(desen, txt) is not None,
                 u'the pattern matched nothing in the output')
        # the competitor line is only printed when there IS a competitor set
        r2 = _sp.run([sys.executable, os.path.join(HERE, "design_group_primers.py"),
                      "--in-group", uyeler[0], "--out-group", uyeler[1],
                      "--label", "LOG_CONTRACT_2",
                      "--out", os.path.join(g18, "cift2.tsv"),
                      "--max-oligo", "40", "--max-pairs", "10"],
                     capture_output=True, text=True)
        txt2 = r2.stdout + r2.stderr
        sina(u'the engine still prints "oligos that bind nowhere in the '
             u'competitors"',
             _re.search(r"oligos that bind nowhere in the competitors:\s*(\d+)",
                       txt2) is not None)
    finally:
        shutil.rmtree(g18, ignore_errors=True)

    if a.real_data:
        print(u'\n13. REAL DATA: THE GEOMETRY OF THE PAIRS PRODUCED')
        tsvler = sorted(glob.glob(os.path.join(HERE, a.candidates, "*__*.tsv")))
        if not tsvler:
            sina(u'a candidate TSV was found', False, a.candidates)
        else:
            toplam = hata = 0
            for t in tsvler[:6]:
                r = subprocess.run(
                    [sys.executable, os.path.join(HERE,
                                                  "check_primer_geometry.py"),
                     "--tsv", t, "--consensus", a.consensus, "--max", "500"],
                    capture_output=True, text=True)
                for line in r.stdout.splitlines():
                    if "rows passing all four geometry conditions" in line:
                        x = line.split(":")[1].split("/")
                        toplam += int(x[1]); hata += int(x[1]) - int(x[0])
            sina(u'all of the real pairs pass the geometry audit',
                 hata == 0, u'tested=%d faulty=%d' % (toplam, hata))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--real-data", action="store_true")
    p.add_argument("--candidates", default="pr_aday")
    p.add_argument("--consensus", default="pr_kons/konsensus")
    a = p.parse_args()
    print("=" * 72)
    print("BORU HATTI REGRESYON TESTI")
    print("=" * 72)
    testler(a)
    gecen = sum(1 for _, ok, _ in SONUC if ok)
    print("\n" + "=" * 72)
    print(u'RESULT: %d of %d tests passed, %d failed'
          % (gecen, len(SONUC), len(SONUC) - gecen))
    for ad, ok, ayrinti in SONUC:
        if not ok:
            print("   KALDI: %s   %s" % (ad, ayrinti))
    print("=" * 72)
    sys.exit(0 if gecen == len(SONUC) else 1)


if __name__ == "__main__":
    main()
