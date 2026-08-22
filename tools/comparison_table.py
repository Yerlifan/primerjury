#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DORT SUTUNLU KARSILASTIRMA TABLOSU

Kutu kutu, dort yontem yan yana:
  1. ozgun Kraken ciktisi
  2. Ayni veriye yuksek guven esigiyle bakan Kraken
  3. Genis kapsamli veritabaniyla (PlusPFP) kosan Kraken
  4. Bizim hizalama tabanli kimligimiz

1. SUTUN NEDEN KUTUNUN ETIKETIDIR
Kutular (reads_<taxid>.fastq) zaten kaynak calismanin Kraken ciktisina gore ayrildi;
extract_kraken_reads.py o ciktidaki taxid'e gore okumalari cekti. Yani bir
kutunun etiketi, ozgun Kraken'inin o okumalar hakkindaki iddiasinin
kendisidir. Ayri bir dosyadan okumak degil, tanimi geregi budur.

UC OLASI SONUC, UCU DE ACIKCA YAZILIR
  a) PlusPFP bizim kimliklerimizi dogruluyorsa, hem teshis hem sonuc Kraken'in
     kendi diliyle onaylanmis olur. En guclu senaryo.
  b) PlusPFP eski etiketleri tekrarliyorsa, teshisimiz yanlistir. Boyle
     yazilir, yumusatilmaz.
  c) PlusPFP ucuncu bir sey diyorsa, uc yontemin ayrildigi yerler isaretlenir
     ve hangisinin neye dayandigi gosterilir.
Tablo bizim aleyhimize cikarsa da oyle sunulur.

Calistirma:
  python3 comparison_table.py --root <PROJE> --job-a <PlusPFP tarama> --name-a PlusPFP \
                               --job-b <eski VT tarama> --name-b eski --threshold 0.1
  python3 comparison_table.py --selftest
"""
import argparse, ast, csv, glob, os, re, sys
from collections import Counter, defaultdict

# ------------------------------------------------------------------ yardimci
def isimleri_oku(kok):
    """taxid -> ad. ozgun calismanin Kraken ozet betigiyle ayni yol; modul CALISTIRILMAZ,
    kaynak metinden ast ile okunur, numpy'a bagimlilik olmaz."""
    adaylar = [os.path.join(kok, "tools", "blast_ispcr.py") if kok else "",
               os.path.join(os.path.dirname(os.path.abspath(__file__)), "blast_ispcr.py")]
    for yol in adaylar:
        if not yol or not os.path.exists(yol):
            continue
        try:
            agac = ast.parse(open(yol, encoding="utf-8", errors="replace").read())
        except SyntaxError:
            continue
        for d in agac.body:
            if not isinstance(d, ast.Assign):
                continue
            for h in d.targets:
                if isinstance(h, ast.Name) and h.id == "ISIMLER":
                    try:
                        return ast.literal_eval(d.value)
                    except (ValueError, SyntaxError):
                        return {}
    return {}

UST_DUZEY = re.compile(r"(aceae|ales|mycota|mycetes|bacteria|archaea|idae|inae)$", re.I)

def cins(ad):
    """
    Bir addan cins belirtecini cikarir.
    'Methanosarcina mazei (taxid 2209)' -> 'methanosarcina'
    'Microascaceae askomikot'           -> ''   (aile duzeyi, cins degil)
    Cins cikmayan ad karsilastirmaya SOKULMAZ; zorlama esitlik, bu projede
    tekrar tekrar cikan "yanlis cevabi temiz gosterme" hatasinin ta kendisidir.
    """
    if not ad:
        return ""
    a = re.sub(r"\(taxid[^)]*\)", "", ad).strip()
    a = re.sub(r"^(unclassified|Candidatus)\s+", "", a, flags=re.I).strip()
    ilk = a.split()[0] if a.split() else ""
    if not ilk or not ilk[0].isalpha():
        return ""
    if UST_DUZEY.search(ilk):
        return ""
    return ilk.lower()

def ust_taksonlar(ad):
    """Addaki aile/takim gibi ust duzey belirtecleri. Aile duzeyi etiketler
    icin tek karsilastirma zemini budur."""
    if not ad:
        return set()
    a = re.sub(r"\(taxid[^)]*\)", "", ad)
    return {k.lower() for k in re.findall(r"[A-Za-z]+", a) if UST_DUZEY.search(k)}

def karsilastir(a, b):
    """
    'uyusuyor' / 'ayrisiyor' / 'karsilastirilamaz'
    Once cins duzeyinde bakilir. Iki taraftan biri cins tasimiyorsa (aile duzeyi
    etiket gibi) ust duzey belirtec ortakligina bakilir. O da yoksa karar
    verilmez; "karsilastirilamaz" ile "ayrisiyor" ayni cumleye sokulmaz.
    """
    if not a or not b:
        return "karsilastirilamaz"
    ca, cb = cins(a), cins(b)
    if ca and cb:
        return "uyusuyor" if ca == cb else "ayrisiyor"
    ua, ub = ust_taksonlar(a), ust_taksonlar(b)
    if ua and ub:
        return "uyusuyor" if (ua & ub) else "ayrisiyor"
    return "karsilastirilamaz"

# ------------------------------------------------------------------ kraken okuma
def tur_haritasi(rapor):
    """
    tum.report'tan taxid -> (tur_taxid, tur_adi). ozgun calismanin Kraken ozet betigindeki
    haritanin aynisi: kraken2 LCA'yi turun ALTINDA birakabilir (S1, S2 sus
    duzeyleri) ve ayni turun iki susu kutuyu bosuna KARISIK gosterir.
    """
    harita = {}
    if not rapor or not os.path.exists(rapor):
        return harita
    yigin = []
    with open(rapor, errors="replace") as fh:
        for satir in fh:
            a = satir.rstrip("\n").split("\t")
            if len(a) < 6:
                continue
            rank, tx, ham = a[3].strip(), a[4].strip(), a[5]
            derinlik = (len(ham) - len(ham.lstrip(" "))) // 2
            ad = ham.strip()
            while yigin and yigin[-1][0] >= derinlik:
                yigin.pop()
            yigin.append((derinlik, rank, tx, ad))
            tur = next(((t, n) for _, r, t, n in reversed(yigin) if r == "S"), None)
            harita[tx] = tur if tur else (tx, ad)
    return harita

def taxid_cek(k):
    if "(taxid" in k:
        p = k.rsplit("(taxid", 1)[1].strip().rstrip(")").strip()
        return p if p.isdigit() else ""
    return ""

def kutulari_oku(out_yol, rapor_yol):
    """
    esik_<C>.out dosyasini kaynak kutuya boler ve her kutunun hakim kimligini
    tur duzeyinde bulur.
    Doner: {kaynak_taxid: (hakim_ad, oran, toplam, sinifsiz_oran)}
    """
    if not os.path.exists(out_yol):
        return {}
    TUR = tur_haritasi(rapor_yol)
    sayac = defaultdict(Counter)
    adlar = {}
    toplamlar = Counter()
    sinifsiz = Counter()
    with open(out_yol, errors="replace") as fh:
        for satir in fh:
            a = satir.rstrip("\n").split("\t")
            if len(a) < 3:
                continue
            ad = a[1]
            if not (ad.startswith("tx") and "_" in ad):
                continue
            kay = ad[2:ad.index("_")]
            if not kay.isdigit():
                continue
            toplamlar[kay] += 1
            if a[0] != "C":
                sinifsiz[kay] += 1
                continue
            k = a[2].strip()
            tx = taxid_cek(k)
            ttx, tad = TUR.get(tx, (tx, k))
            anahtar = ttx or k
            sayac[kay][anahtar] += 1
            adlar.setdefault(anahtar, tad if ttx else k)
    sonuc = {}
    for kay, top in toplamlar.items():
        c = sayac.get(kay)
        if not c:
            sonuc[kay] = ("", 0.0, top, sinifsiz[kay] / top if top else 0.0)
            continue
        anahtar, n = c.most_common(1)[0]
        sonuc[kay] = (adlar.get(anahtar, anahtar), n / top, top,
                      sinifsiz[kay] / top if top else 0.0)
    return sonuc

def esik_dosyalari(klasor, esik):
    """Istenen esige en yakin esik_<C> dosya ciftini bulur. Bulunamazsa bos."""
    if not klasor or not os.path.isdir(klasor):
        return "", "", None
    en_iyi = None
    for y in glob.glob(os.path.join(klasor, "esik_*.report")):
        m = re.match(r"esik_(.+)\.report$", os.path.basename(y))
        if not m:
            continue
        try:
            c = float(m.group(1))
        except ValueError:
            continue
        if en_iyi is None or abs(c - esik) < abs(en_iyi[0] - esik):
            en_iyi = (c, m.group(1), y)
    if not en_iyi:
        return "", "", None
    c, metin, rap = en_iyi
    return os.path.join(klasor, f"esik_{metin}.out"), rap, c

# ------------------------------------------------------------------ hizalama
def hizalama_oku(kok):
    """
    Bizim hizalama tabanli kimligimiz: 0_TESLIM_RAPOR/kimlik_sonuc.csv
    Sutunlar: taxid, iddia_tur, iddia_cins, eslesen_baz, en_iyi_referans,
              eslesme_cins, SONUC
    Doner: {taxid: (gosterilecek_ad, sonuc_metni, eslesen_baz)}
    """
    adaylar = [os.path.join(kok, "tools", "0_TESLIM_RAPOR", "kimlik_sonuc.csv"),
               os.path.join(kok, "VALIDASYON_v2", "primerler", "PIPELINE_TEMIZ",
                            "cikti", "NIHAI", "kimlik_sonuc.csv")]
    for y in adaylar:
        if not os.path.exists(y):
            continue
        out = {}
        with open(y, encoding="utf-8-sig", errors="replace", newline="") as fh:
            for r in csv.DictReader(fh):
                tx = (r.get("taxid") or "").strip()
                if not tx:
                    continue
                out[tx] = ((r.get("eslesme_cins") or "").strip(),
                           (r.get("SONUC") or "").strip(),
                           (r.get("eslesen_baz") or "").strip(),
                           (r.get("en_iyi_referans") or "").strip())
        if out:
            return out, y
    return {}, ""

# ------------------------------------------------------------------ tablo
def tablo_kur(kok, is_a, ad_a, is_b, ad_b, esik):
    ISIMLER = isimleri_oku(kok)
    hiz, hiz_yol = hizalama_oku(kok)

    a_out0, a_rap0, _ = esik_dosyalari(is_a, 0.0)
    a_oute, a_rape, a_c = esik_dosyalari(is_a, esik)
    b_oute, b_rape, b_c = esik_dosyalari(is_b, esik) if is_b else ("", "", None)

    A0 = kutulari_oku(a_out0, a_rap0) if a_out0 else {}
    AE = kutulari_oku(a_oute, a_rape) if a_oute else {}
    BE = kutulari_oku(b_oute, b_rape) if b_oute else {}

    # Esik sutunu: ikinci veritabani varsa ondan, yoksa birinciden.
    ESIK_KAYNAK = BE if BE else AE
    esik_ad = (ad_b if BE else ad_a) + f" esik={b_c if BE else a_c}"
    genis_ad = f"{ad_a} esik=0"

    kutular = set(A0) | set(AE) | set(BE) | set(hiz) | set(ISIMLER)
    kutular = {k for k in kutular if str(k).isdigit()}

    satirlar = []
    for tx in sorted(kutular, key=lambda x: ISIMLER.get(x, "zzz")):
        kaynak_ad = ISIMLER.get(tx, f"taxid {tx}")
        e_ad, e_oran = (ESIK_KAYNAK.get(tx, ("", 0.0, 0, 0.0))[0],
                        ESIK_KAYNAK.get(tx, ("", 0.0, 0, 0.0))[1])
        g_ad, g_oran = (A0.get(tx, ("", 0.0, 0, 0.0))[0],
                        A0.get(tx, ("", 0.0, 0, 0.0))[1])
        h = hiz.get(tx, ("", "", "", ""))
        h_ad, h_sonuc, h_baz, h_ref = h

        # Zayif atama, atama sayilmaz. %20'nin altinda kalan hakim kimlik
        # "karar yok"tur; ozgun calismanin Kraken ozet betigiyle ayni esik kullanilir.
        e_gos = e_ad if e_oran >= 0.20 else ""
        g_gos = g_ad if g_oran >= 0.20 else ""

        uyum_genis_hiz = karsilastir(g_gos, h_ad)
        uyum_genis_kaynak = karsilastir(g_gos, kaynak_ad)
        uyum_kaynak_hiz = karsilastir(kaynak_ad, h_ad)

        if not g_gos:
            senaryo = "PlusPFP karar vermedi"
        elif uyum_genis_hiz == "uyusuyor" and uyum_kaynak_hiz == "ayrisiyor":
            senaryo = "a) PlusPFP BIZI DOGRULUYOR"
        elif uyum_genis_hiz == "uyusuyor":
            senaryo = "a) PlusPFP bizimle ayni (etiketle de cakisiyor)"
        elif uyum_genis_kaynak == "uyusuyor":
            senaryo = "b) PlusPFP ESKI ETIKETI TEKRARLIYOR"
        elif uyum_genis_hiz == "karsilastirilamaz":
            senaryo = "karsilastirilamaz"
        else:
            senaryo = "c) PlusPFP UCUNCU BIR SEY DIYOR"

        satirlar.append(dict(
            taxid=tx, kaynak=kaynak_ad,
            esik_ad=e_gos, esik_oran=e_oran,
            genis_ad=g_gos, genis_oran=g_oran,
            hiz_ad=h_ad, hiz_sonuc=h_sonuc, hiz_baz=h_baz, hiz_ref=h_ref,
            esikle_coktu=("evet" if (e_oran < 0.20 and tx in ESIK_KAYNAK) else
                          ("hayir" if tx in ESIK_KAYNAK else "olculmedi")),
            senaryo=senaryo))
    return satirlar, dict(esik_ad=esik_ad, genis_ad=genis_ad, hiz_yol=hiz_yol,
                          a_var=bool(A0), e_var=bool(ESIK_KAYNAK), h_var=bool(hiz))

def kisalt(s, n):
    s = s or ""
    return s if len(s) <= n else s[:n - 1] + "…"

def markdown(satirlar, bilgi, esik):
    g = []
    g.append(u'# The Kraken2 rerun, four methods side by side')
    g.append("")
    g.append(u'Every row is a bin. The bins are piles of reads separated')
    g.append(u'according to the original Kraken output (`reads_<taxid>.fastq`).')
    g.append("")
    g.append(u'| # | Column | What it means |')
    g.append("|---|---|---|")
    g.append(u'| 1 | the original Kraken output | The bin\'s label. The reads were already separated by that claim. |')
    g.append(f"| 2 | {bilgi['esik_ad']} | The same data at a high confidence threshold. Weak assignments drop out. |")
    g.append(f"| 3 | {bilgi['genis_ad']} | A database of broad coverage, with no threshold. |")
    g.append(u'| 4 | Alignment based identity | Our own measurement (`kimlik_sonuc.csv`). |')
    g.append("")
    g.append(u'An empty cell: no identity reached %20 of the reads, so there is no decision.')
    g.append("")
    g.append(u'| Bin (the source study\'s Kraken) | Threshold raised | PlusPFP | Alignment | Result |')
    g.append("|---|---|---|---|---|")
    for s in satirlar:
        e = f"{kisalt(s['esik_ad'],32)} ({s['esik_oran']:.0%})" if s["esik_ad"] else "_karar yok_"
        gg = f"{kisalt(s['genis_ad'],32)} ({s['genis_oran']:.0%})" if s["genis_ad"] else "_karar yok_"
        h = kisalt(s["hiz_ad"], 28) or "_yok_"
        isaret = {"a) PlusPFP BIZI DOGRULUYOR": "**dogruluyor**",
                  "b) PlusPFP ESKI ETIKETI TEKRARLIYOR": "**eski etiket**",
                  "c) PlusPFP UCUNCU BIR SEY DIYOR": "**ayrisiyor**"}.get(
                      s["senaryo"], s["senaryo"])
        g.append(f"| {kisalt(s['kaynak'],34)} | {e} | {gg} | {h} | {isaret} |")
    g.append("")
    return "\n".join(g)

def yorum(satirlar, bilgi):
    n = Counter(s["senaryo"] for s in satirlar)
    dogrulayan = n["a) PlusPFP BIZI DOGRULUYOR"]
    ayni = n["a) PlusPFP bizimle ayni (etiketle de cakisiyor)"]
    eski = n["b) PlusPFP ESKI ETIKETI TEKRARLIYOR"]
    ucuncu = n["c) PlusPFP UCUNCU BIR SEY DIYOR"]
    karar_yok = n["PlusPFP karar vermedi"]
    karsilastirilamaz = n["karsilastirilamaz"]
    karar_verilen = dogrulayan + ayni + eski + ucuncu

    g = ["", "## Sonuc", ""]
    if not bilgi["a_var"]:
        g.append(u'There is NO PlusPFP run. The table is incomplete and cannot be interpreted.')
        g.append(u'To run it: `bash kraken_tool.sh threshold-a`')
        return "\n".join(g)
    g.append(f"- Bins decided: {karar_verilen}")
    g.append(f"- PlusPFP confirms our identity: **{dogrulayan}** "
             f"(and in {ayni} bins all three methods already say the same thing)")
    g.append(f"- PlusPFP repeats the old label: **{eski}**")
    g.append(f"- PlusPFP says a third thing: **{ucuncu}**")
    g.append(f"- PlusPFP could not decide: {karar_yok}   not comparable: {karsilastirilamaz}")
    g.append("")
    if karar_verilen == 0:
        g.append(u'No bin was decided. No conclusion can be drawn from this table.')
    elif dogrulayan >= max(eski, ucuncu) and dogrulayan > 0:
        g.append(u'**Reading:** in most bins PlusPFP confirms our alignment based')
        g.append(u'identity. That is the strongest case: both the diagnosis (the problem was the')
        g.append(u'database coverage) and the result are confirmed in Kraken\'s own language.')
    elif eski >= max(dogrulayan, ucuncu) and eski > 0:
        g.append(u'**Reading:** in most bins PlusPFP REPEATS the old labels.')
        g.append(u'That means our diagnosis is WRONG. The problem was not the database')
        g.append(u'coverage; the broad database says the same thing. We write this')
        g.append(u'without softening it: this table is against us.')
        g.append(u'Our alignment based result then has to be defended separately,')
        g.append(u'because two independent Kraken runs stand against it.')
    elif ucuncu > 0:
        g.append(u'**Reading:** in most bins PlusPFP says a THIRD thing, neither')
        g.append("the old label nor our identity. The three methods disagree.")
        g.append(u'What they rest on differs: the source study\'s label is an LCA falling on the nearest')
        g.append("relative in a narrow database, PlusPFP is broad but still")
        g.append(u'whole genome k-mer equality, while our identity is the percent similarity of a')
        g.append(u'marker gene alignment. On the bins that disagree the decision cannot be left to any one of them alone.')
    else:
        g.append(u'**Reading:** there is no clear tendency; it has to be looked at bin by bin.')
    g.append("")
    g.append(u'### The bins that disagree')
    ayr = [s for s in satirlar if s["senaryo"].startswith(("b)", "c)"))]
    if not ayr:
        g.append("Yok.")
    for s in ayr:
        g.append(f"- **{s['kaynak']}**  threshold: {s['esik_ad'] or 'no decision'}  |  "
                 f"PlusPFP: {s['genis_ad'] or 'no decision'}  |  alignment: {s['hiz_ad'] or 'none'}"
                 + (f"  ({s['hiz_baz']} bases, {s['hiz_sonuc']})" if s["hiz_baz"] else ""))
    g.append("")
    g.append(u'### What this table does not say')
    g.append(u'Kraken2 looks at k-mer equality, not at alignment. An organism that is not')
    g.append(u'represented in the database is labelled with its nearest relative, and that')
    g.append(u'happens without an error. Alignment based identity gives a percent similarity')
    g.append(u'but looks only at the marker gene window. The two methods answer the same')
    g.append(u'question from different places; where they overlap that strengthens both, and')
    g.append(u'where they diverge that refutes neither one on its own.')
    return "\n".join(g)

# ------------------------------------------------------------------ selftest
def selftest():
    print("=" * 72)
    print(u'THE COMPARISON TABLE, A TEST WITH KNOWN ANSWERS')
    print("=" * 72)
    hata = 0
    def K(ad, bul, bek):
        nonlocal hata
        ok = bul == bek
        if not ok:
            hata += 1
        print(f"  {'PASS' if ok else 'FAIL'}  {ad:<58} {bul} / {bek}")

    K("cins binomialden cikar", cins("Methanosarcina mazei (taxid 2209)"), "methanosarcina")
    K("aile duzeyi etiket cins vermez", cins("Microascaceae askomikot"), "")
    K("unclassified oneki atilir", cins("unclassified Proteiniphilum"), "proteiniphilum")
    K("bos ad bos doner", cins(""), "")
    K(u'the same genus agrees', karsilastir("Sphaerochaeta associata", "Sphaerochaeta"), "uyusuyor")
    K(u'a different genus disagrees', karsilastir("Trichoderma asperellum", "Petriella"), "ayrisiyor")
    # Aile duzeyi etiket ile cins, ayni cumleye sokulmamali.
    K("cins olmayan taraf karsilastirilamaz sayilir",
      karsilastir("Microascaceae askomikot", "Petriella"), "karsilastirilamaz")
    K("iki aile duzeyi etiket ortak belirtecle uyusur",
      karsilastir("Microascaceae askomikot", "Microascaceae sp."), "uyusuyor")
    K(u'an empty side does not settle it', karsilastir("", "Petriella"), "karsilastirilamaz")

    import tempfile
    with tempfile.TemporaryDirectory() as d:
        rap = ("  50.00\t50\t0\tS\t2208\t    Methanosarcina barkeri\n"
               "  30.00\t30\t30\tS1\t1434107\t      Methanosarcina barkeri 3\n"
               "  20.00\t20\t20\tS\t101201\t    Trichoderma asperellum\n")
        open(os.path.join(d, "esik_0.report"), "w").write(rap)
        satir = (["C\ttx2208_a%d\tMethanosarcina barkeri 3 (taxid 1434107)\t1\t\n" % i for i in range(60)]
                 + ["C\ttx2208_b%d\tMethanosarcina barkeri (taxid 2208)\t1\t\n" % i for i in range(20)]
                 + ["U\ttx2208_u%d\tunclassified (taxid 0)\t1\t\n" % i for i in range(20)]
                 + ["C\ttx101201_c%d\tTrichoderma asperellum (taxid 101201)\t1\t\n" % i for i in range(10)]
                 + ["U\ttx101201_u%d\tunclassified (taxid 0)\t1\t\n" % i for i in range(90)])
        open(os.path.join(d, "esik_0.out"), "w").writelines(satir)
        k = kutulari_oku(os.path.join(d, "esik_0.out"), os.path.join(d, "esik_0.report"))
        K(u'the strains of the same species are gathered under one heading', k["2208"][0], "Methanosarcina barkeri")
        K("hakim oran paydaya sinifsizi katar", round(k["2208"][1], 2), 0.80)
        K(u'a weak bin gives a low ratio', round(k["101201"][1], 2), 0.10)
        o, r, c = esik_dosyalari(d, 0.1)
        K(u'the nearest threshold file is found', c, 0.0)
        o2, r2, c2 = esik_dosyalari(os.path.join(d, "yok"), 0.1)
        K(u'a directory that does not exist does not crash it', c2, None)

    # Senaryo siniflandirmasi. Uc olasilik da ayri ayri sinanir.
    s_a = dict(genis_ad="Petriella musispora", kaynak="Trichoderma asperellum", hiz_ad="Petriella")
    K("senaryo a: PlusPFP bizi dogruluyor",
      _senaryo(s_a), "a) PlusPFP BIZI DOGRULUYOR")
    s_b = dict(genis_ad="Trichoderma asperellum", kaynak="Trichoderma asperellum", hiz_ad="Petriella")
    K("senaryo b: PlusPFP eski etiketi tekrarliyor",
      _senaryo(s_b), "b) PlusPFP ESKI ETIKETI TEKRARLIYOR")
    s_c = dict(genis_ad="Fusarium oxysporum", kaynak="Trichoderma asperellum", hiz_ad="Petriella")
    K("senaryo c: PlusPFP ucuncu bir sey diyor",
      _senaryo(s_c), "c) PlusPFP UCUNCU BIR SEY DIYOR")
    s_d = dict(genis_ad="", kaynak="Trichoderma asperellum", hiz_ad="Petriella")
    K(u'a PlusPFP that does not decide does not count as a disagreement',
      _senaryo(s_d), u'PlusPFP karar vermedi')

    print("=" * 72)
    print("THE TEST PASSED" if hata == 0 else f"THE TEST FAILED, {hata} items")
    print("=" * 72)
    return 0 if hata == 0 else 1

def _senaryo(s):
    """selftest icin senaryo kurali, tablo_kur icindekiyle ayni mantik."""
    g, a, h = s["genis_ad"], s["kaynak"], s["hiz_ad"]
    ugh, uga, uah = karsilastir(g, h), karsilastir(g, a), karsilastir(a, h)
    if not g:
        return "PlusPFP karar vermedi"
    if ugh == "uyusuyor" and uah == "ayrisiyor":
        return "a) PlusPFP BIZI DOGRULUYOR"
    if ugh == "uyusuyor":
        return "a) PlusPFP bizimle ayni (etiketle de cakisiyor)"
    if uga == "uyusuyor":
        return "b) PlusPFP ESKI ETIKETI TEKRARLIYOR"
    if ugh == "karsilastirilamaz":
        return "karsilastirilamaz"
    return "c) PlusPFP UCUNCU BIR SEY DIYOR"

def selftest_sessiz():
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        return selftest()

# ------------------------------------------------------------------ ana
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", dest="kok", default="")
    ap.add_argument("--job-a", dest="is_a", default="")
    ap.add_argument("--name-a", dest="ad_a", default="PlusPFP")
    ap.add_argument("--job-b", dest="is_b", default="")
    ap.add_argument("--name-b", dest="ad_b", default="eski VT")
    ap.add_argument("--threshold", type=float, default=0.1)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(selftest())
    if selftest_sessiz() != 0:
        print("THE TEST FAILED, stopped (project rule 2)")
        sys.exit(2)

    satirlar, bilgi = tablo_kur(a.kok, a.is_a, a.ad_a,
                                a.is_b if a.ad_b else "", a.ad_b, a.threshold)
    if not satirlar:
        print(u'ERROR: there is no data to compare.')
        print(u'  The scan has to be run first:  bash kraken_tool.sh threshold')
        sys.exit(1)

    eksikler = []
    if not bilgi["a_var"]:
        eksikler.append(u'there is no PlusPFP run (column 3 is empty)')
    if not bilgi["e_var"]:
        eksikler.append(u'there is no threshold scan (column 2 is empty)')
    if not bilgi["h_var"]:
        eksikler.append(u'kimlik_sonuc.csv was not found (column 4 is empty)')
    if eksikler:
        print(u'WARNING, the table is produced incomplete:')
        for e in eksikler:
            print("  " + e)
        print(u'  A missing column must not be read as a column that was measured and came out empty.\n')

    md = markdown(satirlar, bilgi, a.threshold) + "\n" + yorum(satirlar, bilgi)
    print(md)

    cikti = os.path.join(a.kok, "tools", "0_TESLIM_RAPOR") if a.kok else "."
    os.makedirs(cikti, exist_ok=True)
    mdy = os.path.join(cikti, "KRAKEN_KARSILASTIRMA.md")
    open(mdy, "w", encoding="utf-8").write(md + "\n")
    csvy = os.path.join(cikti, "KRAKEN_KARSILASTIRMA.csv")
    with open(csvy, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(satirlar[0]))
        w.writeheader()
        w.writerows(satirlar)
    print(f"\nwritten:\n  {mdy}\n  {csvy}")

if __name__ == "__main__":
    main()
