#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
community_trends.py
TOPLULUK TREND CALISMA KITABI, ham Bracken ciktisindan yeniden uretilir.

Neden yeniden: onceki kitap tur duzeyindeki sayilari cins duzeyindekilerle
esit agirlikta sunuyordu. Oysa olctugumuz iki bulgu tur duzeyini bazi
taksonlar icin okunamaz kiliyor:
  1. Ayirt edilemez kutular: iki tur kutusunun dizileri birbirinden
     ayrilmiyor (ayirt_edilemez.tsv, olculmus).
  2. Kutu kimligi: bir kutunun ham okumalari kendi atandigi turu tercih
     etmiyor, baska bir referansa gidiyor (kimlik_*.tsv, olculmus).
Her iki dosya da olcumdur; bu betik guvenilirlik isaretini onlardan
TURETIR, elle yazilmis bir liste tutmaz. Dosyalar verilmezse isaret
konmaz ve bu durum kapak sayfasinda acikca yazilir.

Grup ve yil eslemesi barkod NUMARASINDAN kurulur, klasor adindan degil:
kaynak calismanin klasor adlarinda F1 ve F2 yer degistirmisti.

Kullanim:
  python3 community_trends.py --bracken "bracken results" \
      --ayirt primer_adaylari/ayirt_edilemez.tsv \
      --kimlik t_kimlik/kimlik_A.tsv t_kimlik/kimlik_B.tsv \
      --adlar taxid_adlari.tsv --out Topluluk_Trend.xlsx
"""
import argparse, csv, datetime, glob, math, os, re, sys
from collections import defaultdict

from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
from openpyxl.chart import BarChart, LineChart, Reference

YAZI = "Arial"
BASLIK_DOLGU = PatternFill("solid", fgColor="1F4E78")
BASLIK_YAZI = Font(name=YAZI, bold=True, color="FFFFFF", size=10)
NORMAL = Font(name=YAZI, size=10)
KALIN = Font(name=YAZI, bold=True, size=10)
BUYUK = Font(name=YAZI, bold=True, size=14, color="1F4E78")
INCE = Side(style="thin", color="BFBFBF")
KENAR = Border(left=INCE, right=INCE, top=INCE, bottom=INCE)
YESIL = PatternFill("solid", fgColor="C6EFCE")
SARI = PatternFill("solid", fgColor="FFEB9C")
KIRMIZI = PatternFill("solid", fgColor="FFC7CE")

# Otoriter esleme: ornekleme haritasi. Barkod numarasi belirleyicidir.
GRUP_ARALIK = [("A1", 1, "Arke, kisa amplikon"), ("A2", 5, "Arke, uzun amplikon"),
               ("F2", 9, "Mantar, uzun amplikon"), ("F1", 13, "Mantar, kisa amplikon"),
               ("B", 17, "Bakteri")]
YILLAR = (2021, 2023, 2024, 2025)
BARKOD_GRUP, BARKOD_YIL, GRUP_ACIKLAMA = {}, {}, {}
for _g, _b0, _ac in GRUP_ARALIK:
    GRUP_ACIKLAMA[_g] = _ac
    for _i, _y in enumerate(YILLAR):
        BARKOD_GRUP[_b0 + _i] = _g
        BARKOD_YIL[_b0 + _i] = _y
GRUPLAR = [g for g, _, _ in GRUP_ARALIK]
RUTBE_ADI = {"D": "alem üstü", "K": "alem", "P": "şube", "C": "sınıf",
             "O": "takım", "F": "aile", "G": "cins", "S": "tür"}


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bracken", required=True, help="'bracken results' directory")
    p.add_argument("--ayirt", default=None, help="ayirt_edilemez.tsv")
    p.add_argument("--kimlik", nargs="*", default=[], help="kimlik_*.tsv")
    p.add_argument("--adlar", default=None, help="taxid_adlari.tsv")
    p.add_argument("--out", required=True)
    p.add_argument("--ust", type=int, default=10, help="taxon shown on the page")
    p.add_argument("--rutbe", default=None,
                   help="abundance_rank.py'nin output directory; if given "
                        "'Rutbe Kapsamasi' ve 'Guvenilir Bolluk' sayfalari eklenir")
    p.add_argument("--kimlik-esik", type=float, default=50.0,
                   help="bin identity is suspect when the dominant fraction is below this value")
    return p.parse_args()


# ------------------------------------------------------------------ veri
def bracken_oku(kok, duzey):
    """Doner: {barkod: {takson: yuzde}}"""
    out = {}
    desen = os.path.join(kok, duzey, "*", "*_bracken_%s.txt" % duzey)
    for p in sorted(glob.glob(desen)):
        if p.endswith("_report.txt"):
            continue
        m = re.search(r"barcode(\d+)", os.path.basename(p))
        if not m:
            continue
        bc = int(m.group(1))
        if bc not in BARKOD_GRUP:
            continue
        d = {}
        with open(p, encoding="utf-8") as fh:
            basliklar = fh.readline().rstrip("\n").split("\t")
            try:
                i_ad = basliklar.index("name")
                i_or = basliklar.index("fraction_total_reads")
            except ValueError:
                continue
            for line in fh:
                q = line.rstrip("\n").split("\t")
                if len(q) > max(i_ad, i_or):
                    d[q[i_ad]] = float(q[i_or]) * 100.0
        out[bc] = d
    return out


def shannon(v):
    s = sum(v)
    if s <= 0:
        return 0.0
    return -sum((x / s) * math.log(x / s) for x in v if x > 0)


def simpson(v):
    s = sum(v)
    if s <= 0:
        return 0.0
    return 1.0 - sum((x / s) ** 2 for x in v if x > 0)


def bray_curtis(a, b):
    tk = set(a) | set(b)
    pay = sum(abs(a.get(t, 0.0) - b.get(t, 0.0)) for t in tk)
    payda = sum(a.get(t, 0.0) for t in tk) + sum(b.get(t, 0.0) for t in tk)
    return pay / payda if payda else 0.0


def cins_adi(tur_adi):
    """Tur adindan cins adini cikarir. 'Candidatus' ve kisaltmasi 'Ca.'
    cins adi degildir; bunlar atlanir, yoksa butun Candidatus taksonlari
    tek bir 'Ca.' cinsinde toplanir."""
    if not tur_adi:
        return ""
    p = tur_adi.split()
    i = 0
    while i < len(p) and p[i].rstrip(".").lower() in ("candidatus", "ca"):
        i += 1
    return p[i] if i < len(p) else ""


def guvenilirlik_kur(a):
    """Olculmus dosyalardan supheli TUR ve CINS kumelerini turetir.
    Doner: (supheli_tur -> gerekce, supheli_cins -> gerekce)"""
    ad = {}
    if a.adlar and os.path.exists(a.adlar):
        for l in open(a.adlar, encoding="utf-8"):
            q = l.rstrip("\n").split("\t")
            if len(q) > 1:
                ad[q[0]] = q[1]
    tur, cins = {}, {}

    def ekle(sozluk, anahtar, metin):
        if not anahtar:
            return
        sozluk.setdefault(anahtar, [])
        if metin not in sozluk[anahtar]:
            sozluk[anahtar].append(metin)

    if a.ayirt and os.path.exists(a.ayirt):
        for r in csv.DictReader(open(a.ayirt, encoding="utf-8"), delimiter="\t"):
            t1, t2 = r.get("taxid1", ""), r.get("taxid2", "")
            a1, a2 = ad.get(t1, t1), ad.get(t2, t2)
            oz = r.get("kati_ozdeslik", r.get("ozdeslik_yuzde", ""))
            g = ("dizi düzeyinde %s ile ayrılmıyor (katı özdeşlik %%%s, "
                 "sınıf %s)" % (a2, oz, r.get("sinif", "")))
            ekle(tur, a1, g)
            ekle(tur, a2, ("dizi düzeyinde %s ile ayrılmıyor (katı özdeşlik "
                           "%%%s, sınıf %s)" % (a1, oz, r.get("sinif", ""))))
            for x in (a1, a2):
                ekle(cins, cins_adi(x),
                     "bu cinste ayırt edilemeyen tür çifti var")

    for yol in a.kimlik:
        if not os.path.exists(yol):
            continue
        for r in csv.DictReader(open(yol, encoding="utf-8"), delimiter="\t"):
            m = re.search(r"reads[-_](\d+)\.fastq", r.get("dosya", ""))
            if not m:
                continue
            tx = m.group(1)
            adi = ad.get(tx, tx)
            try:
                oran = float(r.get("baskin_oran", "0"))
            except ValueError:
                continue
            bref = r.get("baskin_referans", "")
            kendi = cins_adi(adi)
            # IKI AYRI SINYAL, ayri ayri raporlanir:
            #  (a) okumalar baska bir CINSE gidiyor  -> kutu kimligi yanlis
            #  (b) hicbir referans cogunluk saglamiyor -> kutu belirsiz
            yanlis_cins = bool(kendi) and kendi not in bref
            cogunluk_yok = oran < a.kimlik_esik
            if yanlis_cins:
                ekle(tur, adi,
                     "kutunun ham okumaları başka bir cinse gidiyor "
                     "(%s kutusunda baskın referans: %s, %%%.1f)"
                     % (r.get("grup", ""), bref[:52], oran))
                ekle(cins, kendi,
                     "bu cinsin kutularında okumalar başka cinslere gidiyor")
            elif cogunluk_yok:
                ekle(tur, adi,
                     "hiçbir referans okumaların çoğunluğunu almıyor, kutu "
                     "tür düzeyinde belirsiz (%s kutusu, en iyi referans "
                     "%s, yalnız %%%.1f)"
                     % (r.get("grup", ""), bref[:44], oran))
                ekle(cins, kendi,
                     "bu cinsin kutularında hiçbir referans çoğunluk sağlamıyor")
    return tur, cins


# ------------------------------------------------------------------ yazim
def yaz_baslik(ws, basliklar, satir=1, dondur=True):
    for j, b in enumerate(basliklar, 1):
        c = ws.cell(row=satir, column=j, value=b)
        c.font = BASLIK_YAZI
        c.fill = BASLIK_DOLGU
        c.alignment = Alignment(horizontal="center", vertical="center",
                                wrap_text=True)
        c.border = KENAR
    if dondur:
        ws.freeze_panes = ws.cell(row=satir + 1, column=1)


def genislik(ws, gs):
    for j, g in enumerate(gs, 1):
        ws.column_dimensions[get_column_letter(j)].width = g


def bolum_yaz(ws, satir, veri, ust, supheli, duzey):
    """Bir grup icin takson x yil tablosu yazar. Doner: sonraki satir."""
    basliklar = ["Takson"] + ["%d" % y for y in YILLAR] + ["Dört yıl ortalaması"]
    if duzey == "tur":
        basliklar.append("Güvenilirlik")
    yaz_baslik(ws, basliklar, satir=satir, dondur=False)
    ilk = satir + 1
    for i, (tk, vals) in enumerate(veri[:ust]):
        r = ilk + i
        c = ws.cell(row=r, column=1, value=tk); c.font = NORMAL; c.border = KENAR
        for j, v in enumerate(vals):
            c = ws.cell(row=r, column=2 + j, value=round(v, 4))
            c.font = NORMAL; c.number_format = "0.00"; c.border = KENAR
        c = ws.cell(row=r, column=6,
                    value="=AVERAGE(%s%d:%s%d)" % ("B", r, "E", r))
        c.font = NORMAL; c.number_format = "0.00"; c.border = KENAR
        if duzey == "tur":
            g = supheli.get(tk)
            c = ws.cell(row=r, column=7,
                        value=("CİNS DÜZEYİNDE OKUYUN: " + "; ".join(g)) if g
                        else "ölçüm bu türü şüpheli göstermedi")
            c.font = NORMAL; c.border = KENAR
            c.alignment = Alignment(wrap_text=True, vertical="top")
            c.fill = KIRMIZI if g else YESIL
    son = ilk + min(len(veri), ust) - 1
    # toplam kontrolu: gosterilen taksonlar + digerleri
    r = son + 1
    c = ws.cell(row=r, column=1, value="Gösterilenlerin toplamı")
    c.font = KALIN; c.border = KENAR
    for j in range(4):
        col = get_column_letter(2 + j)
        c = ws.cell(row=r, column=2 + j,
                    value="=SUM(%s%d:%s%d)" % (col, ilk, col, son))
        c.font = KALIN; c.number_format = "0.00"; c.border = KENAR
    return r + 2, ilk, son


def main():
    a = get_args()
    cins = bracken_oku(a.bracken, "genus")
    tur = bracken_oku(a.bracken, "species")
    if not cins:
        sys.exit(u'no bracken file at genus level was found: %s' % a.bracken)
    print("cins duzeyi barkod: %d, tur duzeyi barkod: %d" % (len(cins), len(tur)))
    eksik = [b for b in sorted(BARKOD_GRUP) if b not in cins]
    if eksik:
        print(u'WARNING: a barcode missing at genus level: %s' % eksik)

    s_tur, s_cins = guvenilirlik_kur(a)
    print(u'suspect species from the measured files: %d, suspect genera: %d'
          % (len(s_tur), len(s_cins)))

    wb = Workbook()

    # ---------------- Kapak ----------------
    ws = wb.active; ws.title = "Kapak ve Yöntem"
    genislik(ws, [30, 104])
    ws["A1"] = "PrimerJury topluluk trend analizi"
    ws["A1"].font = BUYUK; ws.merge_cells("A1:B1")
    satirlar = [
        ("Üretim zamanı", datetime.datetime.now().strftime("%d.%m.%Y %H:%M")),
        ("Kaynak", os.path.abspath(a.bracken)),
        ("", ""),
        ("OKUMA SIRASI", ""),
        ("0", "Önce 'Rütbe Kapsaması' sayfası: her örneğin bolluğu hangi "
              "rütbede okunabilir. Ölçüldü: arke örneklerinde okumaların "
              "%86-97'si cins düzeyine iner, mantar örneklerinde %0,1-2,5, "
              "bakteri örneklerinde %8-10. Bolluk sayıları 'Güvenilir Bolluk' "
              "sayfasından alınmalıdır."),
        ("0b", "Bracken ÇALIŞTIRILMADI. Bracken üst rütbede kalan okumaları "
               "veritabanı önceliklerine göre aşağı dağıtır; bu, gerçek "
               "organizmanın veritabanında bulunduğu varsayımına dayanır ve "
               "bu numunede mantar ile bakteri tarafında o varsayım "
               "kurulamıyor. Ayrıntı: bolluk_rutbe_kaniti.md"),
        ("1", "Baskın Cins ve Baskın Tür sayfaları ÖZGÜN Bracken çıktısıdır, "
              "güven eşiği uygulanmamıştır; karşılaştırma için bırakılmıştır."),
        ("2", "Tür düzeyi sayfaları yalnızca Güvenilirlik sütunu yeşil olan "
              "satırlar için tür düzeyinde okunabilir. Kırmızı satırlar cins "
              "düzeyinde toplanarak yorumlanmalıdır."),
        ("", ""),
        ("NEDEN", ""),
        ("Ayırt edilemeyen kutular",
         "Bazı tür kutularının dizileri birbirinden ayrılmıyor; hangi kutuya "
         "kaç okuma düştüğü dizi kanıtına dayanmıyor. Kaynak: ayirt_edilemez.tsv"),
        ("Kutu kimliği",
         "Bazı kutuların ham okumaları kendi atandıkları türü tercih etmiyor, "
         "başka bir referansa gidiyor. Kraken2 varsayılan ayarında çekimser "
         "kalmaz; gerçek tür veritabanında yoksa okuma en yakın kardeş türe "
         "düşer. Kaynak: kimlik_*.tsv"),
        ("Bu işaret elle konmadı",
         "Şüpheli tür ve cins listesi yukarıdaki iki ÖLÇÜM dosyasından "
         "türetilmiştir. Dosya verilmezse işaret konmaz ve bu satırda belirtilir."),
        ("", ""),
        ("Örnekleme eşlemesi", "Grup ve yıl BARKOD NUMARASINDAN kurulur, "
                               "klasör adından değil."),
    ]
    r = 3
    for k, v in satirlar:
        ws.cell(row=r, column=1, value=k).font = KALIN
        c = ws.cell(row=r, column=2, value=v)
        c.font = NORMAL; c.alignment = Alignment(wrap_text=True, vertical="top")
        if len(v) > 90:
            ws.row_dimensions[r].height = 15 * (1 + len(v) // 90)
        r += 1
    r += 1
    ws.cell(row=r, column=1, value="Grup").font = BASLIK_YAZI
    ws.cell(row=r, column=1).fill = BASLIK_DOLGU
    ws.cell(row=r, column=2, value="Barkodlar ve yıllar").font = BASLIK_YAZI
    ws.cell(row=r, column=2).fill = BASLIK_DOLGU
    r += 1
    for g, b0, ac in GRUP_ARALIK:
        ws.cell(row=r, column=1, value="%s  %s" % (g, ac)).font = NORMAL
        ws.cell(row=r, column=2, value=", ".join(
            "barcode%02d=%d" % (b0 + i, YILLAR[i]) for i in range(4))).font = NORMAL
        r += 1
    ws.cell(row=r + 1, column=1, value="Kaynak dosya sayısı").font = KALIN
    ws.cell(row=r + 1, column=2,
            value="cins %d barkod, tür %d barkod" % (len(cins), len(tur))).font = NORMAL
    ws.cell(row=r + 2, column=1, value="Şüpheli işareti").font = KALIN
    ws.cell(row=r + 2, column=2,
            value=("%d tür, %d cins işaretlendi" % (len(s_tur), len(s_cins)))
            if (s_tur or s_cins) else
            "ÖLÇÜM DOSYASI VERİLMEDİ, tür düzeyi işaretlenemedi").font = NORMAL
    if not (s_tur or s_cins):
        ws.cell(row=r + 2, column=2).fill = SARI

    # ---------------- Metrik sozlugu ----------------
    ws = wb.create_sheet("Metrik Sözlüğü")
    genislik(ws, [26, 100])
    yaz_baslik(ws, ["Metrik", "Ne ölçer, nasıl okunur"])
    for i, (k, v) in enumerate([
        ("Bolluk (%)", "Bracken'in fraction_total_reads sütunu, yüze çevrilmiş. "
                       "Okuma sayısı değil, o örnekteki okumaların oranıdır."),
        ("Shannon", "Hem kaç takson olduğunu hem ne kadar dengeli dağıldıklarını "
                    "birlikte ölçer. Büyüdükçe çeşitlilik artar. Birkaç takson "
                    "baskınsa düşük kalır."),
        ("Zenginlik", "Sıfırdan büyük bolluğa sahip takson sayısı. Dengeyi "
                      "hesaba katmaz, yalnızca kaç tane olduğunu sayar."),
        ("Simpson (1-D)", "Rastgele iki okumanın FARKLI taksona ait olma "
                          "olasılığı. Sıfıra yakınsa topluluk tek taksonun "
                          "elinde, bire yakınsa dengeli."),
        ("Bray-Curtis", "İki örnek arasındaki farklılık. 0 aynı, 1 tamamen "
                        "ayrı. Yıllar arası değişimin büyüklüğünü verir."),
        ("Güvenilirlik", "Bu satırın tür düzeyinde okunup okunamayacağı. "
                         "Ölçüm dosyalarından türetilir, elle konmaz."),
    ], start=2):
        ws.cell(row=i, column=1, value=k).font = KALIN
        c = ws.cell(row=i, column=2, value=v)
        c.font = NORMAL; c.alignment = Alignment(wrap_text=True, vertical="top")
        ws.row_dimensions[i].height = 30

    # ---------------- Baskin cins / tur ----------------
    grafik_yeri = {}
    for duzey, veri, ad_sayfa in (("cins", cins, "Baskın Cins"),
                                  ("tur", tur, "Baskın Tür")):
        if not veri:
            continue
        ws = wb.create_sheet(ad_sayfa)
        genislik(ws, [46, 12, 12, 12, 12, 18] + ([74] if duzey == "tur" else []))
        r = 1
        # Bu sayfalar ozgun Bracken ciktisindan gelir; guven esigi
        # UYGULANMAMISTIR. Rutbe kapsamasi olculduyse hangi sayfaya
        # bakilacagi acikca yazilir, yoksa iki farkli cins tablosu yan yana
        # durur ve okuyan hangisinin gecerli oldugunu bilemez.
        if a.rutbe and os.path.exists(os.path.join(a.rutbe, "ozet.tsv")):
            c = ws.cell(row=1, column=1,
                        value="BU SAYFA GÜVEN EŞİĞİ UYGULANMAMIŞ Bracken "
                              "çıktısındandır ve karşılaştırma için "
                              "bırakılmıştır. Bolluk okunacaksa 'Güvenilir "
                              "Bolluk' sayfası kullanılmalıdır; hangi örneğin "
                              "hangi rütbede okunabildiği 'Rütbe Kapsaması' "
                              "sayfasındadır. Bu sayfadaki tür ve cins "
                              "atamalarının çoğu, güven eşiği uygulandığında "
                              "ayakta kalmamaktadır.")
            c.font = Font(name=YAZI, bold=True, size=10, color="9C0006")
            c.fill = KIRMIZI
            c.alignment = Alignment(wrap_text=True, vertical="center")
            ws.merge_cells(start_row=1, start_column=1, end_row=1,
                           end_column=7 if duzey == "tur" else 6)
            ws.row_dimensions[1].height = 48
            r = 3
        if duzey == "tur":
            c = ws.cell(row=r, column=1,
                        value="DİKKAT: Kırmızı işaretli satırlar tür düzeyinde "
                              "okunmamalıdır; o taksonun kutusu ölçümle şüpheli "
                              "bulunmuştur. Gerekçe Güvenilirlik sütunundadır. "
                              "Cins düzeyi sayfası esas alınmalıdır.")
            c.font = Font(name=YAZI, bold=True, size=10, color="9C0006")
            c.fill = KIRMIZI
            c.alignment = Alignment(wrap_text=True, vertical="center")
            ws.merge_cells(start_row=r, start_column=1, end_row=r, end_column=7)
            ws.row_dimensions[r].height = 34
            r += 2
        for g in GRUPLAR:
            bcs = [b for b in sorted(BARKOD_GRUP) if BARKOD_GRUP[b] == g
                   and b in veri]
            if not bcs:
                continue
            c = ws.cell(row=r, column=1,
                        value="%s  %s  (%s)" % (g, GRUP_ACIKLAMA[g],
                                                ", ".join("barcode%02d" % b for b in bcs)))
            c.font = BUYUK
            r += 1
            tumu = set()
            for b in bcs:
                tumu |= set(veri[b])
            sirali = sorted(tumu,
                            key=lambda t: -sum(veri[b].get(t, 0.0) for b in bcs))
            tablo = [(t, [veri[b].get(t, 0.0) for b in bcs]) for t in sirali]
            supheli = s_tur if duzey == "tur" else s_cins
            r, ilk, son = bolum_yaz(ws, r, tablo, a.ust, supheli, duzey)
            grafik_yeri.setdefault(ad_sayfa, []).append((g, ilk, son))
        # grafikler
        for g, ilk, son in grafik_yeri.get(ad_sayfa, []):
            ch = BarChart(); ch.type = "col"; ch.grouping = "clustered"
            ch.title = "%s  baskın %s, yıllara göre" % (g, duzey)
            ch.y_axis.title = "bolluk (%)"; ch.x_axis.title = "takson"
            veri_ref = Reference(ws, min_col=2, max_col=5, min_row=ilk - 1,
                                 max_row=son)
            kat = Reference(ws, min_col=1, min_row=ilk, max_row=son)
            ch.add_data(veri_ref, titles_from_data=True)
            ch.set_categories(kat)
            ch.height, ch.width = 8, 24
            ws.add_chart(ch, "I%d" % ilk)

    # ---------------- Alfa cesitlilik ----------------
    for duzey, veri, ad_sayfa in (("cins", cins, "Alfa Çeşitlilik Cins"),
                                  ("tur", tur, "Alfa Çeşitlilik Tür")):
        if not veri:
            continue
        ws = wb.create_sheet(ad_sayfa)
        genislik(ws, [10, 10, 10, 14, 14, 14])
        yaz_baslik(ws, ["Grup", "Yıl", "Barkod", "Shannon", "Zenginlik",
                        "Simpson (1-D)"])
        r = 2
        for g in GRUPLAR:
            for b in sorted(BARKOD_GRUP):
                if BARKOD_GRUP[b] != g or b not in veri:
                    continue
                v = [x for x in veri[b].values() if x > 0]
                for j, val in enumerate([g, BARKOD_YIL[b], "barcode%02d" % b,
                                         round(shannon(v), 4), len(v),
                                         round(simpson(v), 4)], start=1):
                    c = ws.cell(row=r, column=j, value=val)
                    c.font = NORMAL; c.border = KENAR
                    if j in (4, 6):
                        c.number_format = "0.000"
                r += 1
        ch = LineChart()
        ch.title = "Shannon çeşitliliği, yıllara göre (%s düzeyi)" % duzey
        ch.y_axis.title = "Shannon"; ch.x_axis.title = "örnek"
        ch.add_data(Reference(ws, min_col=4, min_row=1, max_row=r - 1),
                    titles_from_data=True)
        ch.set_categories(Reference(ws, min_col=3, min_row=2, max_row=r - 1))
        ch.height, ch.width = 9, 26
        ws.add_chart(ch, "H2")

    # ---------------- Bray-Curtis ----------------
    for duzey, veri, ad_sayfa in (("cins", cins, "Bray-Curtis Cins"),
                                  ("tur", tur, "Bray-Curtis Tür")):
        if not veri:
            continue
        ws = wb.create_sheet(ad_sayfa)
        genislik(ws, [10, 14, 14, 14])
        c = ws.cell(row=1, column=1,
                    value="Her grup icinde ardışık yılların Bray-Curtis "
                          "farklılığı. 0 aynı, 1 tamamen ayrı.")
        c.font = NORMAL
        ws.merge_cells("A1:D1")
        yaz_baslik(ws, ["Grup", "Karşılaştırma", "Bray-Curtis", "Yorum"], satir=3)
        r = 4
        for g in GRUPLAR:
            bcs = [b for b in sorted(BARKOD_GRUP) if BARKOD_GRUP[b] == g
                   and b in veri]
            for i in range(len(bcs) - 1):
                d = bray_curtis(veri[bcs[i]], veri[bcs[i + 1]])
                yorum = ("küçük değişim" if d < 0.3 else
                         "orta değişim" if d < 0.6 else "büyük değişim")
                for j, val in enumerate([g, "%d - %d" % (BARKOD_YIL[bcs[i]],
                                                         BARKOD_YIL[bcs[i + 1]]),
                                         round(d, 4), yorum], start=1):
                    c = ws.cell(row=r, column=j, value=val)
                    c.font = NORMAL; c.border = KENAR
                    if j == 3:
                        c.number_format = "0.000"
                        c.fill = YESIL if d < 0.3 else (SARI if d < 0.6 else KIRMIZI)
                r += 1

    # ---------------- Rutbe kapsamasi ve guvenilir bolluk ----------------
    # Bracken CALISTIRILMADI. Sebep olculdu: guven duzeltmesinden sonra
    # arke okumalarinin %86-97'si cins duzeyine inebiliyor, mantar
    # okumalarinin %0,1-1,5'i. Bracken ust rutbede kalanlari veritabani
    # onceliklerine gore asagi dagitir; gercek organizma veritabaninda
    # yoksa bu olcum degil sayi uretmek olur.
    rk = os.path.join(a.rutbe, "rutbe_kapsamasi.tsv") if a.rutbe else None
    ro = os.path.join(a.rutbe, "ozet.tsv") if a.rutbe else None
    rb = os.path.join(a.rutbe, "bolluk.tsv") if a.rutbe else None
    if rk and os.path.exists(ro):
        ozetler = list(csv.DictReader(open(ro, encoding="utf-8"), delimiter="\t"))
        ws = wb.create_sheet("Rütbe Kapsaması")
        genislik(ws, [12, 8, 8, 14, 14, 14, 18, 16, 14, 14])
        c = ws.cell(row=1, column=1,
                    value="Bolluk hangi rütbede okunabilir. 'Cins oranı' o "
                          "örnekteki okumaların yüzde kaçının cins ya da daha "
                          "dar bir rütbeye yerleşebildiğini gösterir. Düşük "
                          "oran, o gruptaki organizmaların Kraken2 "
                          "veritabanında temsil edilmediği anlamına gelir. "
                          "Bracken bu yüzden çalıştırılmadı: üst rütbede "
                          "kalanları cinse dağıtmak sayı üretmek olurdu.")
        c.font = NORMAL; c.alignment = Alignment(wrap_text=True, vertical="top")
        ws.merge_cells("A1:J1"); ws.row_dimensions[1].height = 46
        yaz_baslik(ws, ["Barkod", "Grup", "Yıl", "Toplam okuma", "Sınıflanan",
                        "Sınıflanmamış", "Seçilen rütbe", "Bu rütbede yerleşen",
                        "Cins oranı %", "Tür oranı %"], satir=3)
        r = 4
        for x in ozetler:
            vals = [x["barkod"], x["grup"], int(x["yil"]), int(x["toplam_okuma"]),
                    int(x["siniflanan"]), int(x["sinifsiz"]),
                    x["secilen_rutbe_adi"], int(x["bu_rutbede_yerlesen"]),
                    float(x["cins_orani"]), float(x["tur_orani"])]
            for j, v in enumerate(vals, 1):
                cc = ws.cell(row=r, column=j, value=v)
                cc.font = NORMAL; cc.border = KENAR
                if j in (9, 10):
                    cc.number_format = "0.0"
            co = float(x["cins_orani"])
            ws.cell(row=r, column=9).fill = (YESIL if co >= 70 else
                                             SARI if co >= 30 else KIRMIZI)
            ws.cell(row=r, column=7).fill = (YESIL if x["secilen_rutbe"] in ("G", "S")
                                             else SARI if x["secilen_rutbe"] == "F"
                                             else KIRMIZI)
            r += 1
        if rb and os.path.exists(rb):
            ws = wb.create_sheet("Güvenilir Bolluk")
            genislik(ws, [12, 8, 8, 12, 52, 12, 14, 12])
            c = ws.cell(row=1, column=1,
                        value="Her örnek, o örnekte verinin desteklediği en "
                              "dar rütbede verilmiştir. Rütbe elle seçilmedi; "
                              "sınıflanmış okumaların yarısından çoğunun "
                              "yerleşebildiği en dar rütbe seçildi. "
                              "Yerleşemeyen okumalar gizlenmedi, ayrı satırda "
                              "duruyor.")
            c.font = NORMAL; c.alignment = Alignment(wrap_text=True, vertical="top")
            ws.merge_cells("A1:H1"); ws.row_dimensions[1].height = 34
            yaz_baslik(ws, ["Barkod", "Grup", "Yıl", "Rütbe", "Takson",
                            "taxid", "Okuma", "Yüzde"], satir=3)
            r = 4
            for x in csv.DictReader(open(rb, encoding="utf-8"), delimiter="\t"):
                vals = [x["barkod"], x["grup"], int(x["yil"]),
                        RUTBE_ADI.get(x["rutbe"], x["rutbe"]), x["takson"],
                        x["taxid"], int(x["okuma"]), float(x["yuzde"])]
                for j, v in enumerate(vals, 1):
                    cc = ws.cell(row=r, column=j, value=v)
                    cc.font = NORMAL; cc.border = KENAR
                    if j == 8:
                        cc.number_format = "0.00"
                if x["takson"].startswith("["):
                    for j in range(1, 9):
                        ws.cell(row=r, column=j).fill = SARI
                r += 1
    elif a.rutbe:
        print(u'WARNING: --rutbe was given but there is no ozet.tsv: %s' % ro)

    # ---------------- Guvenilirlik kaniti ----------------
    ws = wb.create_sheet("Tür Düzeyi Güvenilirlik")
    genislik(ws, [40, 16, 96])
    c = ws.cell(row=1, column=1,
                value="Aşağıdaki işaretlerin tamamı ölçüm dosyalarından "
                      "türetilmiştir; elle yazılmış bir liste yoktur.")
    c.font = NORMAL; ws.merge_cells("A1:C1")
    yaz_baslik(ws, ["Takson", "Düzey", "Ölçülen gerekçe"], satir=3)
    r = 4
    for tk in sorted(s_tur):
        for g in s_tur[tk]:
            for j, v in enumerate([tk, "tür", g], start=1):
                c = ws.cell(row=r, column=j, value=v)
                c.font = NORMAL; c.border = KENAR
                c.alignment = Alignment(wrap_text=True, vertical="top")
            ws.cell(row=r, column=1).fill = KIRMIZI
            r += 1
    for tk in sorted(s_cins):
        for g in s_cins[tk]:
            for j, v in enumerate([tk, "cins", g], start=1):
                c = ws.cell(row=r, column=j, value=v)
                c.font = NORMAL; c.border = KENAR
                c.alignment = Alignment(wrap_text=True, vertical="top")
            ws.cell(row=r, column=1).fill = SARI
            r += 1
    if r == 4:
        ws.cell(row=4, column=1,
                value="Ölçüm dosyası verilmedi, işaret konamadı.").font = NORMAL
        ws.cell(row=4, column=1).fill = SARI

    d = os.path.dirname(os.path.abspath(a.out))
    if d:
        os.makedirs(d, exist_ok=True)
    wb.save(a.out)
    print("yazildi: %s" % a.out)
    print("sayfa: %s" % ", ".join(w.title for w in wb.worksheets))
    print(u'suspect species marks: %d, suspect genus marks: %d'
          % (len(s_tur), len(s_cins)))


if __name__ == "__main__":
    main()
