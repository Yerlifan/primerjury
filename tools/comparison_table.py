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
  python3 comparison_table.py --kok <PROJE> --is-a <PlusPFP tarama> --ad-a PlusPFP \
                               --is-b <eski VT tarama> --ad-b eski --esik 0.1
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
    g.append("# Kraken2 yeniden kosu, dort yontem yan yana")
    g.append("")
    g.append("Her satir bir kutudur. Kutular ozgun Kraken ciktisina gore")
    g.append("ayrilmis okuma yiginlaridir (`reads_<taxid>.fastq`).")
    g.append("")
    g.append("| # | Sutun | Ne demek |")
    g.append("|---|---|---|")
    g.append("| 1 | ozgun Kraken ciktisi | Kutunun etiketi. Okumalar zaten bu iddiaya gore ayrildi. |")
    g.append(f"| 2 | {bilgi['esik_ad']} | Ayni veri, yuksek guven esigi. Zayif atamalar duser. |")
    g.append(f"| 3 | {bilgi['genis_ad']} | Genis kapsamli veritabani, esik yok. |")
    g.append("| 4 | Hizalama tabanli kimlik | Bizim olcumumuz (`kimlik_sonuc.csv`). |")
    g.append("")
    g.append("Bos hucre: hicbir kimlik okumalarin %20'sine ulasmadi, yani karar yok.")
    g.append("")
    g.append("| Kutu (kaynak calismanin Kraken'i) | Esik yukseltilmis | PlusPFP | Hizalama | Sonuc |")
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
        g.append("PlusPFP kosusu YOK. Tablo eksik, yorum yapilamaz.")
        g.append("Kosmak icin: `bash kraken_tool.sh esik-a`")
        return "\n".join(g)
    g.append(f"- Karar verilen kutu: {karar_verilen}")
    g.append(f"- PlusPFP bizim kimligimizi dogruluyor: **{dogrulayan}** "
             f"(ayrica {ayni} kutuda zaten uc yontem de ayni seyi soyluyor)")
    g.append(f"- PlusPFP eski etiketi tekrarliyor: **{eski}**")
    g.append(f"- PlusPFP ucuncu bir sey diyor: **{ucuncu}**")
    g.append(f"- PlusPFP karar veremedi: {karar_yok}   karsilastirilamayan: {karsilastirilamaz}")
    g.append("")
    if karar_verilen == 0:
        g.append("Karar verilen kutu yok. Bu tablodan sonuc cikarilamaz.")
    elif dogrulayan >= max(eski, ucuncu) and dogrulayan > 0:
        g.append("**Okuma:** en cok kutuda PlusPFP bizim hizalama tabanli kimligimizi")
        g.append("dogruluyor. Bu, en guclu senaryodur: hem teshis (sorun veritabani")
        g.append("kapsamiydi) hem de sonuc, Kraken'in kendi diliyle onaylanmis olur.")
    elif eski >= max(dogrulayan, ucuncu) and eski > 0:
        g.append("**Okuma:** en cok kutuda PlusPFP eski etiketleri TEKRARLIYOR.")
        g.append("Bu, teshisimizin YANLIS oldugu anlamina gelir. Sorun veritabani")
        g.append("kapsami degildi; genis veritabani da ayni seyi soyluyor. Bunu")
        g.append("yumusatmadan yaziyoruz: bu tablo bizim aleyhimize.")
        g.append("Hizalama tabanli sonucumuz bu durumda ayrica savunulmalidir,")
        g.append("cunku iki bagimsiz Kraken kosusu ona karsi duruyor.")
    elif ucuncu > 0:
        g.append("**Okuma:** en cok kutuda PlusPFP UCUNCU bir sey soyluyor, yani ne")
        g.append("eski etiketi ne bizim kimligimizi. Uc yontem ayrisiyor.")
        g.append("Neye dayandiklari farkli: kaynak calismanin etiketi dar veritabaninda en")
        g.append("yakin akrabaya dusen LCA, PlusPFP genis ama yine tam genom k-mer")
        g.append("esitligi, bizim kimligimiz ise marker gen hizalamasinin yuzde")
        g.append("benzerligi. Ayrisan kutularda karar tek basina hicbirine birakilamaz.")
    else:
        g.append("**Okuma:** belirgin bir egilim yok, kutu kutu bakilmali.")
    g.append("")
    g.append("### Ayrisan kutular")
    ayr = [s for s in satirlar if s["senaryo"].startswith(("b)", "c)"))]
    if not ayr:
        g.append("Yok.")
    for s in ayr:
        g.append(f"- **{s['kaynak']}**  esik: {s['esik_ad'] or 'karar yok'}  |  "
                 f"PlusPFP: {s['genis_ad'] or 'karar yok'}  |  hizalama: {s['hiz_ad'] or 'yok'}"
                 + (f"  ({s['hiz_baz']} baz, {s['hiz_sonuc']})" if s["hiz_baz"] else ""))
    g.append("")
    g.append("### Bu tablonun soylemedigi")
    g.append("Kraken2 k-mer esitligine bakar, hizalamaya bakmaz. Veritabaninda")
    g.append("temsil edilmeyen bir organizma en yakin akrabaya etiketlenir ve bu")
    g.append("hata verilmeden olur. Hizalama tabanli kimlik ise yuzde benzerlik")
    g.append("verir ama yalnizca marker gen penceresine bakar. Iki yontem ayni")
    g.append("soruya farkli yerlerden cevap veriyor; ortusmeleri guclendirici,")
    g.append("ayrismalari ise tek basina hicbirini curutucu degildir.")
    return "\n".join(g)

# ------------------------------------------------------------------ selftest
def selftest():
    print("=" * 72)
    print("KARSILASTIRMA TABLOSU, BILINEN CEVAPLI SINAV")
    print("=" * 72)
    hata = 0
    def K(ad, bul, bek):
        nonlocal hata
        ok = bul == bek
        if not ok:
            hata += 1
        print(f"  {'GECTI' if ok else 'KALDI'}  {ad:<58} {bul} / {bek}")

    K("cins binomialden cikar", cins("Methanosarcina mazei (taxid 2209)"), "methanosarcina")
    K("aile duzeyi etiket cins vermez", cins("Microascaceae askomikot"), "")
    K("unclassified oneki atilir", cins("unclassified Proteiniphilum"), "proteiniphilum")
    K("bos ad bos doner", cins(""), "")
    K("ayni cins uyusur", karsilastir("Sphaerochaeta associata", "Sphaerochaeta"), "uyusuyor")
    K("farkli cins ayrisir", karsilastir("Trichoderma asperellum", "Petriella"), "ayrisiyor")
    # Aile duzeyi etiket ile cins, ayni cumleye sokulmamali.
    K("cins olmayan taraf karsilastirilamaz sayilir",
      karsilastir("Microascaceae askomikot", "Petriella"), "karsilastirilamaz")
    K("iki aile duzeyi etiket ortak belirtecle uyusur",
      karsilastir("Microascaceae askomikot", "Microascaceae sp."), "uyusuyor")
    K("bos taraf karar verdirmez", karsilastir("", "Petriella"), "karsilastirilamaz")

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
        K("ayni turun suslari tek kalemde toplanir", k["2208"][0], "Methanosarcina barkeri")
        K("hakim oran paydaya sinifsizi katar", round(k["2208"][1], 2), 0.80)
        K("zayif kutu dusuk oran verir", round(k["101201"][1], 2), 0.10)
        o, r, c = esik_dosyalari(d, 0.1)
        K("en yakin esik dosyasi bulunur", c, 0.0)
        o2, r2, c2 = esik_dosyalari(os.path.join(d, "yok"), 0.1)
        K("olmayan klasor cokmez", c2, None)

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
    K("karar vermeyen PlusPFP ayrisma sayilmaz",
      _senaryo(s_d), "PlusPFP karar vermedi")

    print("=" * 72)
    print("SINAV GECTI" if hata == 0 else f"SINAV KALDI, {hata} madde")
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
    ap.add_argument("--root", "--kok", dest="kok", default="")
    ap.add_argument("--is-a", dest="is_a", default="")
    ap.add_argument("--ad-a", dest="ad_a", default="PlusPFP")
    ap.add_argument("--is-b", dest="is_b", default="")
    ap.add_argument("--ad-b", dest="ad_b", default="eski VT")
    ap.add_argument("--esik", type=float, default=0.1)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(selftest())
    if selftest_sessiz() != 0:
        print("SINAV BASARISIZ, durduruldu (proje kurali 2)")
        sys.exit(2)

    satirlar, bilgi = tablo_kur(a.kok, a.is_a, a.ad_a,
                                a.is_b if a.ad_b else "", a.ad_b, a.esik)
    if not satirlar:
        print(u'ERROR: there is no data to compare.')
        print(u'  The scan has to be run first:  bash kraken_tool.sh esik')
        sys.exit(1)

    eksikler = []
    if not bilgi["a_var"]:
        eksikler.append("PlusPFP kosusu yok (sutun 3 bos)")
    if not bilgi["e_var"]:
        eksikler.append("esik taramasi yok (sutun 2 bos)")
    if not bilgi["h_var"]:
        eksikler.append("kimlik_sonuc.csv bulunamadi (sutun 4 bos)")
    if eksikler:
        print("UYARI, tablo eksik uretiliyor:")
        for e in eksikler:
            print("  " + e)
        print(u'  A missing column must not be read as a column that was measured and came out empty.\n')

    md = markdown(satirlar, bilgi, a.esik) + "\n" + yorum(satirlar, bilgi)
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
    print(f"\nyazildi:\n  {mdy}\n  {csvy}")

if __name__ == "__main__":
    main()
