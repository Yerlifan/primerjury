#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GUVEN ESIGI TARAMASI OZETI VE IKI VERITABANININ KARSILASTIRMASI

kraken_tool.sh'in urettigi esik_<C>.report ve esik_<C>.out dosyalarini okur.

NE SORULUYOR
Kraken2 bir okumayi ancak k-mer'lerinin belli bir bolumu ayni klada gidiyorsa
atar. Esik 0 iken tek bir k-mer bile atamaya yeter. Esik yukseldikce zayif
atamalar duser. Ne kadar dustugu, o atamalarin bastan ne kadar zayif oldugunun
Kraken'in kendi agziyla olcusudur.

ASIL SORU (iki veritabani verildiginde)
Eski veritabaninda esik yukselince coken atamalar, PlusPFP'de ayakta kaliyor mu?
Kaliyorsa teshis dogrulanir: sorun kapsamdi, veritabani o organizmalari
icermiyordu ve Kraken en yakin akrabaya atiyordu. Kalmiyorsa teshis yanlistir
ve bu acikca yazilir.

IKI BAGIMSIZ OLCUM (proje kurali 1)
Alan yuzdeleri iki ayri dosyadan, iki ayri kod yoluyla hesaplanir:
  olcum 1: esik_<C>.report icindeki klad okuma sayilari
  olcum 2: esik_<C>.out icindeki okuma basina atamalar, agac uzerinden toplanir
Ikisi ayrilirsa satir AYRILIK olarak isaretlenir ve sessizce birine gecilmez.
Bu, ayristirma hatasini yakalamak icindir; bu projede on hatanin hepsi
"program hata vermeden yanlis cevap uretiyor" turundendi.

Calistirma:
  python3 threshold_summary.py --is <klasor> --kok <PROJE>
  python3 threshold_summary.py --is <eski> --ad "eski VT" --is2 <yeni> --ad2 "PlusPFP" --kok <PROJE>
  python3 threshold_summary.py --selftest
"""
import argparse, csv, glob, os, re, sys
from collections import Counter, defaultdict

# Alan tanimlari. taxid NCBI icindir; ad, sentetik taksonomili ozel veritabani
# icin yedektir. Mantar okaryotun ALTINDADIR, ikisi ayri sutun ama ic ice.
ALANLAR = [
    ("arke",    "2157",  {"archaea"}),
    ("bakteri", "2",     {"bacteria"}),
    ("okaryot", "2759",  {"eukaryota", "eukarya"}),
    ("mantar",  "4751",  {"fungi"}),
    ("bitki",   "33090", {"viridiplantae", "plantae"}),
    ("virus",   "10239", {"viruses"}),
]

# ------------------------------------------------------------------ okuma
def rapor_oku(yol):
    """
    Kraken2 raporunu okur.
    Sutunlar: yuzde, klad_okuma, dugum_okuma, rank, taxid, girintili ad.
    Doner: (dugumler, ebeveyn, toplam_okuma, siniflandirilmayan)

    Girinti agaci verir, iki bosluk bir duzey. Ebeveyn haritasi buradan cikar
    ve olcum 2'nin toplamasi bu haritayla yapilir.
    """
    dugumler = []
    ebeveyn = {}
    yigin = []
    siniflandirilmayan = 0
    kok_klad = 0
    if not os.path.exists(yol):
        return [], {}, 0, 0
    with open(yol, errors="replace") as fh:
        for satir in fh:
            a = satir.rstrip("\n").split("\t")
            if len(a) < 6:
                continue
            try:
                klad = int(a[1]); dugum = int(a[2])
            except ValueError:
                continue
            rank, tx, ham = a[3].strip(), a[4].strip(), a[5]
            derinlik = (len(ham) - len(ham.lstrip(" "))) // 2
            ad = ham.strip()
            if rank == "U":
                siniflandirilmayan = klad
                continue
            if tx == "1":
                kok_klad = klad
            while yigin and yigin[-1][0] >= derinlik:
                yigin.pop()
            ebeveyn[tx] = yigin[-1][1] if yigin else ""
            yigin.append((derinlik, tx))
            dugumler.append(dict(klad=klad, dugum=dugum, rank=rank,
                                 taxid=tx, ad=ad, derinlik=derinlik))
    toplam = siniflandirilmayan + kok_klad
    return dugumler, ebeveyn, toplam, siniflandirilmayan

def alan_dugumu(dugumler, taxid, adlar):
    """Alanin dugumunu once taxid ile, bulamazsa adla arar."""
    for d in dugumler:
        if d["taxid"] == taxid:
            return d
    for d in dugumler:
        if d["ad"].strip().lower() in adlar:
            return d
    return None

def alan_raporundan(dugumler, toplam):
    """OLCUM 1: rapordaki klad okuma sayilarindan alan yuzdeleri."""
    sonuc = {}
    for ad, tx, adlar in ALANLAR:
        d = alan_dugumu(dugumler, tx, adlar)
        sonuc[ad] = d["klad"] if d else 0
    return sonuc

def taxid_cek(kimlik):
    """'Methanosarcina mazei (taxid 2209)' -> '2209'. Bulamazsa bos."""
    if "(taxid" in kimlik:
        p = kimlik.rsplit("(taxid", 1)[1].strip().rstrip(")").strip()
        return p if p.isdigit() else ""
    return ""

def out_oku(yol):
    """
    Kraken2 ciktisi: C/U, okuma adi, kimlik, uzunluk, k-mer haritasi.
    Doner: (sayac_taxid, toplam, siniflandirilmayan, kutu_sayac)
    kutu_sayac: kaynak taxid (@tx<taxid>_ oneki) -> Counter(atanan taxid)
    """
    sayac = Counter()
    kutu = defaultdict(Counter)
    toplam = 0
    sinifsiz = 0
    if not os.path.exists(yol):
        return sayac, 0, 0, kutu
    with open(yol, errors="replace") as fh:
        for satir in fh:
            a = satir.rstrip("\n").split("\t")
            if len(a) < 3:
                continue
            toplam += 1
            ad = a[1]
            kaynak = ""
            if ad.startswith("tx") and "_" in ad:
                k = ad[2:ad.index("_")]
                if k.isdigit():
                    kaynak = k
            if a[0] != "C":
                sinifsiz += 1
                if kaynak:
                    kutu[kaynak]["U"] += 1
                continue
            tx = taxid_cek(a[2]) or a[2].strip()
            sayac[tx] += 1
            if kaynak:
                kutu[kaynak][tx] += 1
    return sayac, toplam, sinifsiz, kutu

def atalar(tx, ebeveyn, sinir=200):
    """tx ve butun atalari. Dongu olusursa sinirda durur, sonsuza gitmez."""
    out = []
    g = tx
    n = 0
    while g and n < sinir:
        out.append(g)
        g = ebeveyn.get(g, "")
        n += 1
    return out

def alan_outtan(sayac, ebeveyn, dugumler):
    """
    OLCUM 2: okuma basina atamalari agac uzerinden alanlara toplar.
    Rapordan tamamen bagimsiz bir yol degildir (agac rapordan gelir) ama
    SAYIM bagimsizdir; bu, klad sayisi ayristirma hatalarini yakalar.
    """
    hedef = {}
    for ad, tx, adlar in ALANLAR:
        d = alan_dugumu(dugumler, tx, adlar)
        if d:
            hedef[d["taxid"]] = ad
    sonuc = {ad: 0 for ad, _, _ in ALANLAR}
    onbellek = {}
    for tx, n in sayac.items():
        if tx not in onbellek:
            bulunan = [hedef[a] for a in atalar(tx, ebeveyn) if a in hedef]
            onbellek[tx] = bulunan
        for ad in onbellek[tx]:
            sonuc[ad] += n
    return sonuc

# ------------------------------------------------------------------ tarama
def esik_listesi(klasor):
    """esik_<C>.report dosyalarini sayisal siraya dizer."""
    out = []
    for y in glob.glob(os.path.join(klasor, "esik_*.report")):
        m = re.match(r"esik_(.+)\.report$", os.path.basename(y))
        if not m:
            continue
        try:
            c = float(m.group(1))
        except ValueError:
            continue
        out.append((c, m.group(1), y))
    return sorted(out)

def tarama_oku(klasor):
    """
    Bir tarama klasorunu okur.
    Doner: (satirlar, kutular)
      satirlar: [dict(esik, toplam, sinifsiz, sinifsiz_oran, alanlar, alanlar2, ayrilik)]
      kutular : {esik_metni: {kaynak_taxid: Counter(atanan)}}
    """
    satirlar = []
    kutular = {}
    for c, metin, rap in esik_listesi(klasor):
        dugumler, ebeveyn, toplam, sinifsiz = rapor_oku(rap)
        a1 = alan_raporundan(dugumler, toplam)
        out_yol = os.path.join(klasor, f"esik_{metin}.out")
        sayac, toplam2, sinifsiz2, kutu = out_oku(out_yol)
        a2 = alan_outtan(sayac, ebeveyn, dugumler) if sayac else {}
        if kutu:
            kutular[metin] = kutu
        # Iki olcum karsilastirilir. Payda olarak raporun toplami kullanilir;
        # out dosyasi varsa toplamlarin da esit olmasi beklenir.
        ayrilik = []
        if toplam2 and toplam and toplam2 != toplam:
            ayrilik.append(f"okuma sayisi rapor {toplam} / out {toplam2}")
        if a2:
            for ad, _, _ in ALANLAR:
                if a1.get(ad, 0) != a2.get(ad, 0):
                    ayrilik.append(f"{ad} rapor {a1.get(ad,0)} / out {a2.get(ad,0)}")
        satirlar.append(dict(
            esik=c, esik_metni=metin, toplam=toplam, sinifsiz=sinifsiz,
            alanlar=a1, alanlar2=a2, ayrilik="; ".join(ayrilik)))
    return satirlar, kutular

def yuzde(n, toplam):
    return (100.0 * n / toplam) if toplam else 0.0

# ------------------------------------------------------------------ yazim
def egri_metni(satirlar, baslik):
    g = []
    g.append("=" * 78)
    g.append(baslik)
    g.append("=" * 78)
    g.append(f"{'esik':>6} {'okuma':>8} {'sinifsiz':>9} " +
             " ".join(f"{a:>9}" for a, _, _ in ALANLAR))
    g.append("-" * 78)
    for s in satirlar:
        t = s["toplam"]
        g.append(f"{s['esik']:>6} {t:>8} {yuzde(s['sinifsiz'], t):>8.2f}% " +
                 " ".join(f"{yuzde(s['alanlar'].get(a,0), t):>8.2f}%" for a, _, _ in ALANLAR))
    g.append("")
    g.append("mantar okaryotun ALTINDADIR, iki sutun ic icedir.")
    ayr = [s for s in satirlar if s["ayrilik"]]
    if ayr:
        g.append("")
        g.append("AYRILIK VAR. Iki olcum ayni sayiyi vermedi, sayilara guvenmeyin:")
        for s in ayr:
            g.append(f"  esik {s['esik']}: {s['ayrilik']}")
    else:
        g.append("iki bagimsiz olcum butun esiklerde ayni sonucu verdi.")
    return "\n".join(g)

def kutu_hakim(kutu_sayac):
    """Bir kutunun en sik atamasi ve orani (siniflandirilamayanlar paydaya dahil)."""
    toplam = sum(kutu_sayac.values())
    if not toplam:
        return "", 0.0
    en = [(tx, n) for tx, n in kutu_sayac.most_common() if tx != "U"]
    if not en:
        return "U", 0.0
    return en[0][0], en[0][1] / toplam

def cokme_esigi(kutular, kaynak, alt=0.20):
    """
    Bir kutunun hakim atamasinin ilk defa %20'nin altina dustugu esik.
    Hic dusmuyorsa None doner. "Coken atama" tam olarak budur.
    """
    for metin in sorted(kutular, key=lambda m: float(m)):
        k = kutular[metin].get(kaynak)
        if not k:
            continue
        _, oran = kutu_hakim(k)
        if oran < alt:
            return float(metin)
    return None

def ayakta_kalma(kutular_a, kutular_b, ad_a, ad_b):
    """
    ASIL SORU. Eski veritabaninda coken atamalar yeni veritabaninda ayakta mi?
    Doner: (satirlar, ozet)
    """
    kaynaklar = set()
    for k in (kutular_a, kutular_b):
        for m in k.values():
            kaynaklar |= set(m.keys())
    satirlar = []
    for kay in sorted(kaynaklar, key=lambda x: (len(x), x)):
        ca = cokme_esigi(kutular_a, kay)
        cb = cokme_esigi(kutular_b, kay)
        if ca is None and cb is None:
            durum = "ikisinde de ayakta"
        elif ca is not None and cb is None:
            durum = "ESKIDE COKTU, YENIDE AYAKTA"
        elif ca is None and cb is not None:
            durum = "eskide ayakta, yenide coktu"
        elif cb > ca:
            durum = "yenide daha dayanikli"
        elif cb < ca:
            durum = "yenide daha kirilgan"
        else:
            durum = "ikisi de ayni esikte coktu"
        satirlar.append(dict(kaynak=kay, cokme_a=ca, cokme_b=cb, durum=durum))
    ozet = Counter(s["durum"] for s in satirlar)
    return satirlar, ozet

# ------------------------------------------------------------------ selftest
def selftest():
    print("=" * 72)
    print("ESIK OZETI, BILINEN CEVAPLI SINAV")
    print("=" * 72)
    hata = 0
    def K(ad, bul, bek):
        nonlocal hata
        ok = bul == bek
        if not ok:
            hata += 1
        print(f"  {'GECTI' if ok else 'KALDI'}  {ad:<58} {bul} / {bek}")

    import tempfile
    # Elle kurulmus, cevabi kagitta hesaplanmis bir rapor.
    # 1000 okuma: 100 sinifsiz, 900 root. Arke 400, bakteri 300, okaryot 200,
    # bunun 150'si mantar.
    rap = ("  10.00\t100\t100\tU\t0\tunclassified\n"
           "  90.00\t900\t10\tR\t1\troot\n"
           "  89.00\t890\t0\tR1\t131567\t  cellular organisms\n"
           "  40.00\t400\t50\tD\t2157\t    Archaea\n"
           "  35.00\t350\t350\tS\t2209\t      Methanosarcina mazei\n"
           "  30.00\t300\t100\tD\t2\t    Bacteria\n"
           "  20.00\t200\t200\tS\t1642647\t      Proteiniphilum saccharofermentans\n"
           "  20.00\t200\t50\tD\t2759\t    Eukaryota\n"
           "  15.00\t150\t0\tK\t4751\t      Fungi\n"
           "  15.00\t150\t150\tS\t101201\t        Trichoderma asperellum\n")
    with tempfile.TemporaryDirectory() as d:
        ry = os.path.join(d, "esik_0.report")
        open(ry, "w").write(rap)
        dug, eb, top, sifsiz = rapor_oku(ry)
        K("toplam okuma sinifsiz + kok kladdir", top, 1000)
        K("siniflandirilamayan okunur", sifsiz, 100)
        K("agac girintiden kurulur, mantar okaryotun altinda", eb.get("4751"), "2759")
        K("tur, alanin altinda dogru baglanir", eb.get("2209"), "2157")
        a1 = alan_raporundan(dug, top)
        K("arke klad sayisi", a1["arke"], 400)
        K("bakteri klad sayisi", a1["bakteri"], 300)
        K("okaryot klad sayisi", a1["okaryot"], 200)
        K("mantar klad sayisi", a1["mantar"], 150)
        K("bulunmayan alan sifirdir, cokmez", a1["bitki"], 0)
        K("yuzde hesabi", round(yuzde(a1["arke"], top), 2), 40.0)

        # OLCUM 2 icin ayni sayilari veren bir out dosyasi.
        oy = os.path.join(d, "esik_0.out")
        satir = []
        satir += ["C\ttx2209_r%d\tMethanosarcina mazei (taxid 2209)\t1500\t\n" % i for i in range(350)]
        satir += ["C\ttx2157_r%d\tArchaea (taxid 2157)\t1500\t\n" % i for i in range(50)]
        satir += ["C\ttx1642647_r%d\tP. saccharofermentans (taxid 1642647)\t1500\t\n" % i for i in range(200)]
        satir += ["C\ttx2_r%d\tBacteria (taxid 2)\t1500\t\n" % i for i in range(100)]
        satir += ["C\ttx101201_r%d\tTrichoderma asperellum (taxid 101201)\t1500\t\n" % i for i in range(150)]
        satir += ["C\ttx2759_r%d\tEukaryota (taxid 2759)\t1500\t\n" % i for i in range(50)]
        satir += ["U\ttx101201_u%d\tunclassified (taxid 0)\t1500\t\n" % i for i in range(100)]
        open(oy, "w").writelines(satir)
        sayac, top2, sifsiz2, kutu = out_oku(oy)
        K("out toplam okuma", top2, 1000)
        K("out siniflandirilamayan", sifsiz2, 100)
        a2 = alan_outtan(sayac, eb, dug)
        K("olcum 2 arke, olcum 1 ile ayni", a2["arke"], 400)
        K("olcum 2 bakteri, olcum 1 ile ayni", a2["bakteri"], 300)
        K("olcum 2 mantar, olcum 1 ile ayni", a2["mantar"], 150)
        K("olcum 2 okaryot mantari da icerir", a2["okaryot"], 200)
        K("kutu okuma adindan cikarilir", sorted(kutu.keys())[0:2], ["101201", "1642647"])

        satirlar, kutular = tarama_oku(d)
        K("tarama tek esik okur", len(satirlar), 1)
        K("ayrilik yok, iki olcum uyusuyor", satirlar[0]["ayrilik"], "")

        # AYRILIK GERCEKTEN YAKALANIYOR MU. Bu madde, olcutun kor olmadigini
        # sinar: bilerek bozulmus bir out dosyasi AYRILIK vermeli.
        open(oy, "a").write("C\ttx2_x1\tBacteria (taxid 2)\t1500\t\n")
        s2, _ = tarama_oku(d)
        K("bozulmus dosya AYRILIK olarak yakalanir", s2[0]["ayrilik"] != "", True)

    # Kutu cokmesi ve ayakta kalma.
    A = {"0":   {"K1": Counter({"9": 90, "U": 10})},
         "0.1": {"K1": Counter({"9": 10, "U": 90})}}
    B = {"0":   {"K1": Counter({"7": 95, "U": 5})},
         "0.1": {"K1": Counter({"7": 92, "U": 8})}}
    K("hakim atama bulunur", kutu_hakim(A["0"]["K1"])[0], "9")
    K("hakim oran paydaya sinifsizi katar", round(kutu_hakim(A["0"]["K1"])[1], 2), 0.9)
    K("cokme esigi bulunur", cokme_esigi(A, "K1"), 0.1)
    K("cokmeyen kutu None doner", cokme_esigi(B, "K1"), None)
    sat, ozet = ayakta_kalma(A, B, "eski", "yeni")
    K("eskide coken, yenide ayakta kalan yakalanir",
      sat[0]["durum"], "ESKIDE COKTU, YENIDE AYAKTA")

    K("dongulu agac sonsuza gitmez", len(atalar("a", {"a": "b", "b": "a"})), 200)
    K("bos rapor cokmez", rapor_oku("/yok/olan/dosya")[2], 0)

    print("=" * 72)
    print("SINAV GECTI" if hata == 0 else f"SINAV KALDI, {hata} madde")
    print("=" * 72)
    return 0 if hata == 0 else 1

def selftest_sessiz():
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        return selftest()

# ------------------------------------------------------------------ ana
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--is", dest="is1", default="")
    ap.add_argument("--ad", default="veritabani 1")
    ap.add_argument("--is2", default="")
    ap.add_argument("--ad2", default="veritabani 2")
    ap.add_argument("--ali", default="")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(selftest())
    if not a.is1:
        ap.error("--is gerekli")
    if selftest_sessiz() != 0:
        print("SINAV BASARISIZ, durduruldu (proje kurali 2)")
        sys.exit(2)

    s1, k1 = tarama_oku(a.is1)
    if not s1:
        print(f"HATA: {a.is1} icinde esik_<C>.report dosyasi yok.")
        print("  Once tarama kosulmali:  bash kraken_tool.sh esik-eski")
        sys.exit(1)

    metin = egri_metni(s1, f"GUVEN ESIGI EGRISI, {a.ad}")
    print(metin)
    yaz_csv(os.path.join(a.is1, "esik_egrisi.csv"), s1, a.ad)
    open(os.path.join(a.is1, "esik_egrisi.txt"), "w", encoding="utf-8").write(metin + "\n")
    print(f"\nyazildi: {os.path.join(a.is1, 'esik_egrisi.csv')}")

    if not a.is2:
        return
    s2, k2 = tarama_oku(a.is2)
    if not s2:
        print(f"\nUYARI: {a.is2} icinde tarama yok, iki veritabani karsilastirilmadi.")
        print("  Kosmak icin:  bash kraken_tool.sh esik-yeni")
        return
    m2 = egri_metni(s2, f"GUVEN ESIGI EGRISI, {a.ad2}")
    print("\n" + m2)
    yaz_csv(os.path.join(a.is2, "esik_egrisi.csv"), s2, a.ad2)
    open(os.path.join(a.is2, "esik_egrisi.txt"), "w", encoding="utf-8").write(m2 + "\n")

    birlesik = yan_yana(s1, s2, a.ad, a.ad2)
    print("\n" + birlesik)
    sat, ozet = ayakta_kalma(k1, k2, a.ad, a.ad2)
    kalma = kalma_metni(sat, ozet, a.ad, a.ad2, a.ali)
    print("\n" + kalma)
    hedef = os.path.join(a.is1, "esik_iki_veritabani.txt")
    open(hedef, "w", encoding="utf-8").write(birlesik + "\n\n" + kalma + "\n")
    with open(os.path.join(a.is1, "esik_ayakta_kalma.csv"), "w",
              newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["kaynak_taxid", f"cokme_esigi_{a.ad}", f"cokme_esigi_{a.ad2}", "durum"])
        for s in sat:
            w.writerow([s["kaynak"], s["cokme_a"] if s["cokme_a"] is not None else "cokmedi",
                        s["cokme_b"] if s["cokme_b"] is not None else "cokmedi", s["durum"]])
    print(f"\nyazildi: {hedef}")

def yaz_csv(yol, satirlar, ad):
    with open(yol, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        basliklar = ["veritabani", "esik", "toplam_okuma", "siniflandirilamayan",
                     "siniflandirilamayan_yuzde"]
        for x, _, _ in ALANLAR:
            basliklar += [f"{x}_okuma", f"{x}_yuzde", f"{x}_yuzde_olcum2"]
        basliklar.append("ayrilik")
        w.writerow(basliklar)
        for s in satirlar:
            t = s["toplam"]
            sat = [ad, s["esik"], t, s["sinifsiz"], round(yuzde(s["sinifsiz"], t), 3)]
            for x, _, _ in ALANLAR:
                sat += [s["alanlar"].get(x, 0),
                        round(yuzde(s["alanlar"].get(x, 0), t), 3),
                        round(yuzde(s["alanlar2"].get(x, 0), t), 3) if s["alanlar2"] else ""]
            sat.append(s["ayrilik"])
            w.writerow(sat)

def yan_yana(s1, s2, ad1, ad2):
    g = ["=" * 96,
         f"IKI VERITABANI YAN YANA   ({ad1}  ve  {ad2})",
         "=" * 96,
         f"{'esik':>6} | " + " | ".join(
             f"{x:>10}" for x in ["sinifsiz", "arke", "bakteri", "mantar"]),
         f"{'':>6} | " + " | ".join(
             f"{ad1[:4]:>4}/{ad2[:5]:>5}" for _ in range(4)),
         "-" * 96]
    h2 = {s["esik"]: s for s in s2}
    for s in s1:
        o = h2.get(s["esik"])
        hucre = []
        for anahtar in ["sinifsiz", "arke", "bakteri", "mantar"]:
            v1 = (yuzde(s["sinifsiz"], s["toplam"]) if anahtar == "sinifsiz"
                  else yuzde(s["alanlar"].get(anahtar, 0), s["toplam"]))
            if o:
                v2 = (yuzde(o["sinifsiz"], o["toplam"]) if anahtar == "sinifsiz"
                      else yuzde(o["alanlar"].get(anahtar, 0), o["toplam"]))
                hucre.append(f"{v1:>4.1f}/{v2:>5.1f}")
            else:
                hucre.append(f"{v1:>4.1f}/{'yok':>5}")
        g.append(f"{s['esik']:>6} | " + " | ".join(f"{h:>10}" for h in hucre))
    eksik = [s["esik"] for s in s1 if s["esik"] not in h2]
    if eksik:
        g.append("")
        g.append(f"UYARI: su esikler {ad2} tarafinda YOK, karsilastirilmadi: {eksik}")
        g.append("  Eksik esik, karsilastirilmis esik gibi okunmamalidir.")
    return "\n".join(g)

def kalma_metni(sat, ozet, ad1, ad2, ali):
    isimler = isimleri_oku(ali)
    g = ["=" * 96,
         "ASIL SORU: eskide coken atamalar yenide ayakta kaliyor mu",
         "=" * 96,
         "",
         f"Bir kutunun 'coktugu' esik, hakim atamasinin ilk defa okumalarin %20'sinin",
         f"altina dustugu esiktir. Sol sutun {ad1}, sag sutun {ad2}.",
         ""]
    for d, n in ozet.most_common():
        g.append(f"  {n:>3} kutu   {d}")
    g.append("")
    g.append(f"{'kutu':<40}{ad1[:12]:>12}{ad2[:12]:>12}   durum")
    g.append("-" * 96)
    for s in sat:
        a = "cokmedi" if s["cokme_a"] is None else str(s["cokme_a"])
        b = "cokmedi" if s["cokme_b"] is None else str(s["cokme_b"])
        ad = isimler.get(s["kaynak"], f"taxid {s['kaynak']}")
        g.append(f"{ad[:39]:<40}{a:>12}{b:>12}   {s['durum']}")
    g.append("")
    n_dogrulayan = ozet.get("ESKIDE COKTU, YENIDE AYAKTA", 0) + ozet.get("yenide daha dayanikli", 0)
    n_toplam = sum(ozet.values())
    g.append("YORUM")
    if n_toplam == 0:
        g.append("  Karsilastirilacak kutu yok.")
    elif n_dogrulayan > n_toplam / 2:
        g.append(f"  {n_dogrulayan}/{n_toplam} kutuda atamalar yeni veritabaninda daha dayanikli.")
        g.append("  Bu, sorunun KAPSAM oldugu teshisini destekler: eski veritabani bu")
        g.append("  organizmalari icermiyordu, Kraken en yakin akrabaya atiyordu ve o")
        g.append("  atamalar zayifti. Yeni veritabani organizmayi icerdigi icin atama guclu.")
    else:
        g.append(f"  Yalnizca {n_dogrulayan}/{n_toplam} kutuda atamalar yenide daha dayanikli.")
        g.append("  Bu, KAPSAM teshisini DESTEKLEMIYOR. Zayiflik veritabani kapsamindan")
        g.append("  degil, okumalarin kendisinden geliyor olabilir. Teshis gozden gecirilmeli.")
    return "\n".join(g)

def isimleri_oku(ali):
    """taxid -> ad. 86_KRAKEN_OZET.py ile ayni yol, numpy'a bagimlilik yok."""
    import ast
    adaylar = [os.path.join(ali, "WSL", "90_BLAST_ICPCR.py") if ali else "",
               os.path.join(os.path.dirname(os.path.abspath(__file__)), "90_BLAST_ICPCR.py")]
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

if __name__ == "__main__":
    main()
