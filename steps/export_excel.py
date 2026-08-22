#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
export_excel.py gathers every output of the primer design into one Excel file.

The input directories:
  --candidates  the batch design output (ozet.tsv, ayirt_edilemez.tsv, *.log)
  --final       the specificity output (primer_final.tsv)
  --splits      the cluster split output, optional (*_bolme.json)
  --names       the taxid name table
"""
import argparse, csv, datetime, glob, json, os, re, sys
from collections import OrderedDict

from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
from openpyxl.chart import BarChart, Reference

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import field_audit

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


def tsv(p):
    if not p or not os.path.exists(p):
        return []
    with open(p, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def yaz_baslik(ws, basliklar, satir=1):
    for j, b in enumerate(basliklar, 1):
        c = ws.cell(row=satir, column=j, value=b)
        c.font = BASLIK_YAZI
        c.fill = BASLIK_DOLGU
        c.alignment = Alignment(horizontal="center", vertical="center",
                                wrap_text=True)
        c.border = KENAR
    ws.freeze_panes = ws.cell(row=satir + 1, column=1)


def genislik(ws, genislikler):
    for j, g in enumerate(genislikler, 1):
        ws.column_dimensions[get_column_letter(j)].width = g


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--candidates", required=True)
    p.add_argument("--final", required=True)
    p.add_argument("--splits", default=None)
    p.add_argument("--external", default=None,
                   help="output of the external-databases step (dis_veritabani.tsv)")
    p.add_argument("--reference", default=None,
                   help="output of the design-from-reference step (primer_referans.tsv)")
    p.add_argument("--names", default=None)
    p.add_argument("--targets", default="targets.tsv")
    p.add_argument("--identity", default=None,
                   help='the output of the target identity step; when it is '
                        'given, a measured identity column is added to every '
                        'row')
    p.add_argument("--consensus", default=None,
                   help='the dominant allele consensus directory, which the '
                        'field consistency check needs')
    p.add_argument("--out", required=True)
    p.add_argument("--note", default="")
    return p.parse_args()


def main():
    a = get_args()
    ad = {}
    if a.names and os.path.exists(a.names):
        for l in open(a.names, encoding="utf-8"):
            q = l.rstrip("\n").split("\t")
            if len(q) > 1:
                ad[q[0]] = q[1]
    # Rather than producing an Excel that looks valid but is silently empty when the
    # input is missing, it stops plainly; an empty workbook was reading as "not one
    # candidate passed".
    eksik = [y for y in (os.path.join(a.candidates, "ozet.tsv"),
                         os.path.join(a.final, "primer_final.tsv"))
             if not os.path.exists(y)]
    if eksik:
        sys.exit(u'The input file was not found, so no Excel was produced:\n   %s'
                 % "\n   ".join(eksik))
    ozet = tsv(os.path.join(a.candidates, "ozet.tsv"))
    ayirt = tsv(os.path.join(a.candidates, "ayirt_edilemez.tsv"))
    final = tsv(os.path.join(a.final, "primer_final.tsv"))
    hedefler = tsv(a.targets) if os.path.exists(a.targets) else []
    ref_ham = tsv(a.reference) if a.reference else []
    dis_yolu = a.external or os.path.join(a.final, "dis_veritabani.tsv")
    dis_ham = tsv(dis_yolu)
    # (target, class, forward, reverse) -> total off target products
    dis = {}
    dis_detay = {}
    for x in dis_ham:
        k = (x["hedef"], x["sinif"], x["ileri_dizi"], x["geri_dizi"])
        dis[k] = dis.get(k, 0) + int(x.get("hedef_disi_urun", 0) or 0)
        if x.get("ornekler"):
            dis_detay.setdefault(k, []).append("%s: %s"
                                               % (x["veritabani"], x["ornekler"]))

    # The measured identity: the target name comes from the Kraken2 assignment, while
    # what the sequence really is is a separate measurement. The two can come apart;
    # if they do not sit side by side in the table, the reader takes the name to be
    # confirmed by sequence evidence.
    kimlik = {}
    kyol = a.identity or os.path.join(a.final, "hedef_kimlik.tsv")
    for r in tsv(kyol):
        kimlik[r["hedef"]] = r
    if kimlik:
        print(u'measured identity table read: %s (%d targets)' % (kyol, len(kimlik)))
    else:
        print(u'WARNING: there is no measured identity table (%s), so the identity column will stay empty'
              % kyol)

    UYUM_METIN = {
        "tur_uyusuyor": "confirmed at species level",
        "cins_uyusuyor_tur_farkli": "the genus is right, the species differs",
        "CINS_FARKLI": "A DIFFERENT GENUS",
        "YAKIN_AKRABA_YOK": "no close relative in the database",
        "vurus_yok": "no hit in the database",
    }

    def kimlik_destek(hedef):
        'Returns (identity, support, bins, is_there_a_majority)          On a '
        'group target a single measured identity can mislead: the most '
        'frequent name may be a PLURALITY rather than a MAJORITY. Measured: '
        'of the         19 bins of the acetoclastic methanogen target only 6 '
        'go to Methanothrix         soehngenii and the rest spread across '
        'Methanosarcina species. That is why         the support fraction is '
        'always shown, and why anything below half is         called "no '
        'clear majority".'
        r = kimlik.get(hedef)
        if not r:
            return None, 0, 0, False
        try:
            d = int(r.get("destekleyen_kutu", 0))
            n = int(r.get("kutu_sayisi", 0))
        except (TypeError, ValueError):
            d = n = 0
        return r, d, n, (n > 0 and d * 2 > n)

    def kimlik_hucre(hedef):
        r, d, n, cog = kimlik_destek(hedef)
        if not r:
            return "", "", None
        u = r.get("uyum", "")
        ad_metni = r.get("olculen_kimlik", "")
        if n:
            ad_metni = '%s (%d of %d bins)' % (ad_metni, d, n)
        if not cog:
            ad_metni += "  [no clear majority]"
            if r.get("diger"):
                ad_metni += "  the others: " + r["diger"][:70]
            return ad_metni, "the bins do not agree on one identity", SARI
        dolgu = (YESIL if u == "tur_uyusuyor" else
                 SARI if u == "cins_uyusuyor_tur_farkli" else KIRMIZI)
        return ad_metni, UYUM_METIN.get(u, u), dolgu

    wb = Workbook()

    # ---------------- 1. The cover and the method ----------------
    ws = wb.active
    ws.title = "Cover and method"
    genislik(ws, [42, 96])
    ws["A1"] = "PrimerJury primer design: the decisions, applied"
    ws["A1"].font = BUYUK
    ws.merge_cells("A1:B1")
    r = 3
    gecen = [x for x in final if x.get("ozgulluk_durum") == "GECTI"]
    # DOMAIN CONSISTENCY: pairs designed from a locus library the target does not
    # belong to ARE TAKEN OUT of the delivery. They break no rule, but they do not
    # represent the target; if they stay in the table that target looks covered when
    # it is not.
    _ta = field_audit.taxid_alanlari(a.consensus)
    _ht = field_audit.hedef_taxidleri(a.targets)
    alan_disi = []
    temiz = []
    for x in gecen:
        uyumsuz, dag, baskin = field_audit.alan_dagilimi(
            x.get("hedef", ""), x.get("sinif", ""),
            taxid_alan=_ta, hedef_taxid=_ht)
        if uyumsuz:
            x["alan_gerekce"] = field_audit.aciklama(dag, x.get("sinif"), baskin)
            alan_disi.append(x)
        else:
            temiz.append(x)
    if alan_disi:
        print(u'pairs removed from delivery because of a DOMAIN MIX: %d (%s)'
              % (len(alan_disi),
                 ", ".join(sorted(set("%s/%s" % (z["hedef"], z["sinif"])
                                      for z in alan_disi)))))
    gecen = temiz
    for x in gecen:
        k = (x["hedef"], x["sinif"], x["ileri_dizi"], x["geri_dizi"])
        x["dis_hedef_disi"] = dis.get(k, "")
        x["dis_ornek"] = ";".join(dis_detay.get(k, []))[:180]
    kapsanan = sorted(set(x["hedef"] for x in gecen))
    satirlar = [
        ("Produced at", datetime.datetime.now().strftime("%d.%m.%Y %H:%M")),
        ("The candidate directory", a.candidates),
        ("The verification directory", a.final),
        ("Targets in total", str(len(set(x["hedef"] for x in ozet)))),
        ("Targets that produced a candidate", str(len(set(x["hedef"] for x in ozet
                                          if x.get("durum", "").startswith("TAMAM"))))),
        ("Targets that passed every rule", str(len(kapsanan))),
        ("Primer pairs that passed every rule", str(len(gecen))),
        ("", ""),
        ("THE METHOD", ""),
        ("The oligo rule",
         'A, C, G and T only. A length of 18 to 25 bases. GC between 40 and '
         "60 per cent with a hard bound of 35 to 65. A G or a C at the 3' "
         'end. At most three G or C in the last five bases. At most four of '
         'the same base in a row.'),
        ("Termodinamik",
         'Tm between 58 and 62 C with a hard bound of 57 to 63. A hairpin dG '
         'of at least -3000, and a self dimer and hetero dimer of at least '
         '-6000. Tm is measured with two independent libraries, primer3 and '
         'Biopython; the systematic offset between them is computed from the '
         'data and an oligo that departs from that offset is dropped.'),
        ("The product",
         '70 to 250 bases with an upper bound of 300. 90 to 150 scores '
         'highest. Less than 1.5 C between the Tm of the two primers.'),
        ("The binding rule",
         'The last two bases have to match exactly, at most one mismatch in '
         'the last five bases, and at most three across the whole primer. An '
         "overhang on the 5' side is free. The two primers sit on opposite "
         "strands with their 3' ends facing each other."),
        ("Specificity",
         'A product has to form in EVERY targeted member and in NONE of the '
         'competitors. At least one of the primers has to find no binding '
         'site at all in the competitors; a cleanliness that comes from both '
         'primers binding weakly is not accepted.'),
        ("Confirmation on the raw reads",
         'Every pair that passes on the consensus is tested again on the '
         "taxon's raw reads. The product rate in the competitor reads is "
         'judged with a Wilson lower bound, so that single read noise does '
         'not look like a high rate.'),
        ("Cross contamination between bins",
         'Because Kraken bins can leak into one another, how much of the '
         'reads in a competitor bin really belong to the target is measured. '
         'When the product rate seen in a competitor stays below the Wilson '
         'upper bound of that leak, it does not count as evidence that the '
         'competitor was amplified.'),
        ("The consensus",
         'The design is made on the unambiguous consensus produced by calling '
         'the DOMINANT ALLELE at every position of the raw reads. The IUPAC '
         'coded consensus is kept as a second and more cautious measurement; '
         'judging the codes by set intersection was covering up real '
         'differences.'),
        ("The two measurement principle",
         'No decision is left to a single code path. Tm is measured with two '
         'libraries, taxon separation with both an alignment and a k-mer, and '
         'specificity on both the consensus and the raw reads. When the '
         'measurements diverge the candidate is rejected and the divergence '
         'is logged.'),
    ]
    if a.note:
        satirlar.append(("Not", a.note))
    for k, v in satirlar:
        ws.cell(row=r, column=1, value=k).font = KALIN if k else NORMAL
        c = ws.cell(row=r, column=2, value=v)
        c.font = NORMAL
        c.alignment = Alignment(wrap_text=True, vertical="top")
        ws.row_dimensions[r].height = max(15, 13 * (1 + len(v) // 95))
        r += 1

    # ---------------- 2. Primer tablosu ----------------
    ws = wb.create_sheet("Primer Tablosu")
    sut = ["karar", "hedef", "sinif", "ileri_dizi", "ileri_tm", "ileri_gc",
           "geri_dizi", "geri_tm", "geri_gc", "tm_farki", "urun_min",
           "urun_maks", "uye_dogrulanan", "uye_toplam", "rakip_wilson",
           "ileri_baglanma_min", "geri_baglanma_min", "sus_ici_fark",
           "yetim_primer", "heterodimer_dg", "dis_hedef_disi", "ceza"]
    basliklar = ["Karar", "Hedef", "Class", "Forward primer (5' to 3')", "Forward Tm",
                 "Forward GC", "Reverse primer (5' to 3')", "Reverse Tm", "Reverse GC",
                 "Tm difference", "Shortest product", "Longest product", "Members confirmed",
                 "Members in total", "Rakip Wilson", "Forward binding", "Reverse binding",
                 "Within strain difference", "Yetim primer", "Hetero-dimer dG",
                 "Off target products, outside databases", "Ceza"]
    yaz_baslik(ws, basliklar)
    genislik(ws, [7, 34, 7, 30, 9, 9, 30, 9, 9, 9, 11, 11, 13, 11, 13, 13, 13,
                  11, 13, 15, 20, 8])
    i = 2
    for x in sorted(gecen, key=lambda z: (z.get("karar", ""), z["hedef"],
                                          z["sinif"], float(z.get("ceza", 9e9)))):
        for j, k in enumerate(sut, 1):
            v = x.get(k, "")
            try:
                v = float(v) if k not in ("hedef", "sinif", "ileri_dizi",
                                          "geri_dizi", "yetim_primer") else v
            except (TypeError, ValueError):
                pass
            c = ws.cell(row=i, column=j, value=v)
            c.font = Font(name="Consolas", size=10) if k in ("ileri_dizi", "geri_dizi") else NORMAL
            c.border = KENAR
        ws.cell(row=i, column=2).fill = YESIL
        i += 1
    if i == 2:
        ws.cell(row=2, column=1, value="No candidate passed every rule").font = NORMAL

    # ---------------- 2b. The suggested pair ----------------
    ws = wb.create_sheet("Recommended pairs")
    ws["A1"] = ('One pair per target. The ordering criteria, in order: the '
                'number of off target products in the outside databases, the '
                'competitor Wilson lower bound, and the design penalty. When '
                'the outside database measurement was not made, that column '
                'stays empty and the ordering rests on the other two.')
    ws["A1"].font = NORMAL
    ws["A1"].alignment = Alignment(wrap_text=True, vertical="top")
    ws.merge_cells("A1:L1")
    ws.row_dimensions[1].height = 45
    basliklar = ["Karar", "Target (the requested name)", "Class", "Forward primer (5' to 3')",
                 "Forward Tm", "Reverse primer (5' to 3')", "Reverse Tm", "Product (bp)",
                 "Members confirmed", "Rakip Wilson", "Off target, outside databases",
                 "Yetim primer", "Measured identity", 'Does the name agree with the '
                                                      'data']
    yaz_baslik(ws, basliklar, satir=3)
    genislik(ws, [7, 34, 7, 30, 9, 30, 9, 12, 14, 13, 18, 14, 30, 24])

    def _sira(z):
        d = z.get("dis_hedef_disi", "")
        try:
            d = int(d)
        except (TypeError, ValueError):
            d = 10 ** 9
        try:
            w = float(z.get("rakip_wilson", 1))
        except (TypeError, ValueError):
            w = 1.0
        try:
            c = float(z.get("ceza", 9e9))
        except (TypeError, ValueError):
            c = 9e9
        return (d, w, c)

    en_iyi = {}
    for x in gecen:
        k = (x["hedef"], x["sinif"])
        if k not in en_iyi or _sira(x) < _sira(en_iyi[k]):
            en_iyi[k] = x
    i = 4
    for k in sorted(en_iyi, key=lambda z: (en_iyi[z].get("karar", ""), z[0], z[1])):
        x = en_iyi[k]
        urun = ("%s-%s" % (x.get("urun_min", ""), x.get("urun_maks", ""))
                if x.get("urun_min") != x.get("urun_maks")
                else str(x.get("urun_min", "")))
        satir = [x.get("karar", ""), x["hedef"], x["sinif"], x["ileri_dizi"],
                 x.get("ileri_tm", ""), x["geri_dizi"], x.get("geri_tm", ""),
                 urun,
                 "%s/%s" % (x.get("uye_dogrulanan", ""), x.get("uye_toplam", "")),
                 x.get("rakip_wilson", ""), x.get("dis_hedef_disi", ""),
                 x.get("yetim_primer", "")]
        kad, kuy, kdolgu = kimlik_hucre(x["hedef"])
        satir += [kad, kuy]
        for j, v in enumerate(satir, 1):
            if j in (5, 7, 10, 11):
                try:
                    v = float(v)
                except (TypeError, ValueError):
                    pass
            c = ws.cell(row=i, column=j, value=v)
            c.font = Font(name="Consolas", size=10) if j in (4, 6) else NORMAL
            c.border = KENAR
        ws.cell(row=i, column=2).fill = YESIL
        if kdolgu is not None:
            ws.cell(row=i, column=14).fill = kdolgu
        ws.cell(row=i, column=13).alignment = Alignment(wrap_text=True,
                                                        vertical="top")
        i += 1
    if i == 4:
        ws.cell(row=4, column=1, value="There is no pair to recommend").font = NORMAL
    # More than one target with the same measured identity means they target THE SAME
    # organism at sequence level. That is something a reader of the table cannot see
    # by themselves, so it is written out plainly.
    ayni_kimlik = {}
    for k in en_iyi:
        r, d, n, cog = kimlik_destek(k[0])
        # It is called "the same organism" only if MORE THAN HALF of the bins go to the
        # same identity. Settling for a plurality would make group targets whose members
        # spread across different species look merged as well.
        if r and cog and r.get("olculen_kimlik") \
                and "vurus yok" not in r["olculen_kimlik"] \
                and "esik alti" not in r["olculen_kimlik"]:
            ayni_kimlik.setdefault(r["olculen_kimlik"], set()).add(k[0])
    cakisan = {kk: vv for kk, vv in ayni_kimlik.items() if len(vv) > 1}
    if cakisan:
        i += 1
        c = ws.cell(row=i, column=1,
                    value='CAREFUL: the targets below aim at THE SAME '
                          'organism at sequence level. They appear on '
                          'separate rows because the requested names are '
                          'kept, not because two different organisms are '
                          'amplified.')
        c.font = Font(name=YAZI, bold=True, size=10, color="9C0006")
        c.fill = KIRMIZI
        c.alignment = Alignment(wrap_text=True, vertical="center")
        ws.merge_cells(start_row=i, start_column=1, end_row=i, end_column=14)
        ws.row_dimensions[i].height = 32
        for kk, vv in sorted(cakisan.items()):
            i += 1
            c = ws.cell(row=i, column=1, value="%s  ->  %s" % (", ".join(sorted(vv)), kk))
            c.font = NORMAL
            ws.merge_cells(start_row=i, start_column=1, end_row=i, end_column=14)

    # ---------------- 3. The target state ----------------
    ws = wb.create_sheet("Hedef Durumu")
    basliklar = ["Karar", "Hedef", "Level", "Class", "Members", "Rakip",
                 "Pairs at design", "The state of the design", "Kademe",
                 "Passed verification", "Note"]
    yaz_baslik(ws, basliklar)
    genislik(ws, [7, 34, 9, 7, 7, 8, 14, 16, 13, 17, 60])
    gecen_say = {}
    for x in gecen:
        gecen_say[(x["hedef"], x["sinif"])] = gecen_say.get((x["hedef"], x["sinif"]), 0) + 1
    i = 2
    for x in ozet:
        anahtar = (x["hedef"], x["sinif"])
        g = gecen_say.get(anahtar, 0)
        satir = [x.get("karar", ""), x["hedef"], x.get("duzey", ""), x["sinif"],
                 x.get("uye", ""), x.get("rakip", ""), x.get("cift", ""),
                 x.get("durum", ""), x.get("kademe", ""), g, x.get("note", "")]
        for j, v in enumerate(satir, 1):
            try:
                v = int(v)
            except (TypeError, ValueError):
                pass
            c = ws.cell(row=i, column=j, value=v)
            c.font = NORMAL
            c.border = KENAR
            c.alignment = Alignment(wrap_text=(j == 11), vertical="top")
        ws.cell(row=i, column=10).fill = YESIL if g else (
            SARI if str(x.get("durum", "")).startswith("TAMAM") else KIRMIZI)
        i += 1
    son = i - 1
    ws.cell(row=i + 1, column=1, value="Pairs that passed verification in total").font = KALIN
    ws.cell(row=i + 1, column=10, value="=SUM(J2:J%d)" % son).font = KALIN

    ch = BarChart()
    ch.type = "bar"
    ch.title = "Primer pairs per target that passed verification"
    ch.y_axis.title = "the number of pairs"
    ch.x_axis.title = "target and class"
    ch.x_axis.delete = False
    ch.y_axis.delete = False
    veri = Reference(ws, min_col=10, min_row=1, max_row=son)
    kat = Reference(ws, min_col=2, min_row=2, max_row=son)
    ch.add_data(veri, titles_from_data=True)
    ch.set_categories(kat)
    ch.height, ch.width = 14, 26
    ch.legend = None
    ws.add_chart(ch, "M2")

    # ---------------- 4. The indistinguishable bins ----------------
    ws = wb.create_sheet("Indistinguishable bins")
    ws["A1"] = ('Consensuses assigned to different taxa inside the same '
                'amplicon class that cannot be separated at sequence level. '
                'They are taken out of the competitor list, because if the '
                "target's own sequence counts as a competitor the rule that "
                'no product may form in a competitor cannot be met at all.')
    ws["A1"].font = NORMAL
    ws["A1"].alignment = Alignment(wrap_text=True, vertical="top")
    ws.merge_cells("A1:I1")
    ws.row_dimensions[1].height = 45
    basliklar = ["Class", "Taxid 1", "Takson 1", "Taxid 2", "Takson 2",
                 "Hizalanan bp", "Identity over the overlap, per cent", "Strict identity, per cent",
                 "The reason"]
    yaz_baslik(ws, basliklar, satir=3)
    genislik(ws, [7, 11, 34, 11, 34, 13, 20, 18, 20])
    i = 4
    for x in ayirt:
        t1, t2 = x.get("taxid1", ""), x.get("taxid2", "")
        satir = [x.get("sinif", ""), t1, ad.get(t1, x.get("ad1", "")), t2,
                 ad.get(t2, x.get("ad2", "")), x.get("hizalanan_bp", ""),
                 x.get("kesisimli_ozdeslik", x.get("ozdeslik_yuzde", "")),
                 x.get("kati_ozdeslik", ""), x.get("gerekce", "")]
        for j, v in enumerate(satir, 1):
            try:
                v = float(v)
            except (TypeError, ValueError):
                pass
            c = ws.cell(row=i, column=j, value=v)
            c.font = NORMAL
            c.border = KENAR
        i += 1
    if i == 4:
        ws.cell(row=4, column=1, value="No indistinguishable pair was found").font = NORMAL

    # ---------------- 5. The subset split ----------------
    if a.splits and os.path.isdir(a.splits):
        ws = wb.create_sheet("Subset split")
        ws["A1"] = ('Function groups that one primer pair cannot cover were '
                    'split into subsets by designability. The split rests not '
                    'on sequence similarity but on which member blocks a pair '
                    'from forming.')
        ws["A1"].font = NORMAL
        ws["A1"].alignment = Alignment(wrap_text=True, vertical="top")
        ws.merge_cells("A1:E1")
        ws.row_dimensions[1].height = 45
        yaz_baslik(ws, ["Grup", "Subset", "Members", "Valid pairs",
                        "The members"], satir=3)
        genislik(ws, [34, 11, 12, 14, 80])
        i = 4
        for p in sorted(glob.glob(os.path.join(a.splits, "*_bolme.json"))):
            d = json.load(open(p, encoding="utf-8"))
            etiket = os.path.basename(p).replace("_bolme.json", "")
            for s in d.get("sets", []):
                satir = [etiket, "SET%d" % s["set"], len(s["members"]),
                         s["pairs"], ", ".join(s["members"])]
                for j, v in enumerate(satir, 1):
                    c = ws.cell(row=i, column=j, value=v)
                    c.font = NORMAL
                    c.border = KENAR
                    c.alignment = Alignment(wrap_text=(j == 5), vertical="top")
                i += 1
            if d.get("kalan"):
                c = ws.cell(row=i, column=1, value=etiket)
                c.font = NORMAL
                ws.cell(row=i, column=2, value="kapsanamayan").font = NORMAL
                ws.cell(row=i, column=5, value=", ".join(d["kalan"])).font = NORMAL
                ws.cell(row=i, column=2).fill = SARI
                i += 1

    # ---------------- 5b. The external database ----------------
    if dis_ham:
        ws = wb.create_sheet("Outside databases")
        ws["A1"] = ('Every primer was searched against the reference '
                    'databases with blastn, task blastn-short. When two hits '
                    'on the same reference sequence sit on opposite strands '
                    "with their 3' ends facing each other and meet inside the "
                    'product range, it counts as an off target product. Every '
                    'hit was also put through the binding rule. On a '
                    'universal primer a high count is the expected result.')
        ws["A1"].font = NORMAL
        ws["A1"].alignment = Alignment(wrap_text=True, vertical="top")
        ws.merge_cells("A1:F1")
        ws.row_dimensions[1].height = 60
        yaz_baslik(ws, ["Hedef", "Class", "Database", "Forward primer",
                        "Reverse primer", "Off target products", "Examples"], satir=3)
        genislik(ws, [34, 7, 22, 30, 30, 16, 70])
        i = 4
        for x in sorted(dis_ham, key=lambda z: (z["hedef"], z["sinif"],
                                                -int(z.get("hedef_disi_urun", 0) or 0))):
            satir = [x["hedef"], x["sinif"], x["veritabani"], x["ileri_dizi"],
                     x["geri_dizi"], int(x.get("hedef_disi_urun", 0) or 0),
                     x.get("ornekler", "")]
            for j, v in enumerate(satir, 1):
                c = ws.cell(row=i, column=j, value=v)
                c.font = Font(name="Consolas", size=10) if j in (4, 5) else NORMAL
                c.border = KENAR
            if satir[5] == 0:
                ws.cell(row=i, column=6).fill = YESIL
            i += 1

    # ---------------- 5c. Warnings ----------------
    ws = wb.create_sheet("Warnings")
    ws["A1"] = ('Automatic checks. These rows do not invalidate the result, '
                'but they set how the primers have to be read.')
    ws["A1"].font = NORMAL
    ws["A1"].alignment = Alignment(wrap_text=True, vertical="top")
    ws.merge_cells("A1:D1")
    ws.row_dimensions[1].height = 30
    yaz_baslik(ws, ["Kind", "Konu", "Detail", "Yorum"], satir=3)
    genislik(ws, [26, 42, 70, 70])
    uyarilar = []

    # 1. If the same primer pair appears on more than one target, the targets are not separated
    paylasim = {}
    for x in gecen:
        paylasim.setdefault((x["ileri_dizi"], x["geri_dizi"]), set()).add(
            "%s (%s)" % (x["hedef"], x["sinif"]))
    for k, v in sorted(paylasim.items(), key=lambda z: -len(z[1])):
        if len(v) > 1:
            hedef_adlari = set(h.split(" (")[0] for h in v)
            if len(hedef_adlari) < 2:
                yorum = ('The same pair was chosen for one target in two '
                         'different amplicon classes. That is not a problem; '
                         'it shows the pair works in both.')
                tur = "Bilgi"
            else:
                yorum = ('Different targets are met by the same pair. That '
                         'pair DOES NOT separate them; it means the target '
                         'definitions overlap in the data.')
                tur = "Dikkat"
            uyarilar.append((tur, u'The same primer pair on more than one target',
                             "%s / %s" % k, yorum + "  Hedefler: " + ", ".join(sorted(v))))

    # 1a. The target name and the measured identity disagree
    for x in sorted(set(z["hedef"] for z in gecen)):
        r = kimlik.get(x)
        if not r:
            continue
        u = r.get("uyum", "")
        if u in ("tur_uyusuyor", ""):
            continue
        tur = "Dikkat" if u == "cins_uyusuyor_tur_farkli" else "Eksik"
        uyarilar.append((tur, u'The target name and the measured identity disagree',
                         u'%s: Kraken2 says %s, the sequence shows %s (%s/%s bins)'
                         % (x, r.get("kraken_etiketi", "")[:40],
                            r.get("olculen_kimlik", ""),
                            r.get("destekleyen_kutu", "?"),
                            r.get("kutu_sayisi", "?")),
                         u'Evidence: %s. The target name comes from the meeting decision and was kept; which organism the primer really amplifies is in the Measured identity column.'
                         % (r.get("kanit", "")[:150] or "yok")))

    # 1b. Alan karisimi nedeniyle teslimden cikarilan ciftler
    for x in alan_disi:
        uyarilar.append(("Eksik", u'A mixture of domains, taken out of the delivery',
                         "%s (%s): F=%s R=%s"
                         % (x["hedef"], x["sinif"], x["ileri_dizi"],
                            x["geri_dizi"]),
                         x.get("alan_gerekce", "")))
    if alan_disi:
        for h in sorted(set(x["hedef"] for x in alan_disi)):
            uyarilar.append(("Eksik", u'The target was not met',
                             h,
                             u'The design in this target\'s own domain gave no valid pair; the single pair that passed came from the wrong locus library and was removed. The target has to count as NOT MET.'))

    # 2. Ayirt edilemez kutular
    for x in ayirt:
        t1, t2 = x.get("taxid1", ""), x.get("taxid2", "")
        uyarilar.append(("Dikkat", u'An indistinguishable taxon bin',
                         u'%s: %s ~ %s (strict identity %%%s)'
                         % (x.get("sinif", ""), ad.get(t1, t1), ad.get(t2, t2),
                            x.get("kati_ozdeslik", x.get("ozdeslik_yuzde", ""))),
                         u'These two bins do not separate at sequence level. It was taken out of the competitor list; the specificity of that target is not claimed for this taxon.'))

    # 3. Aday uretemeyen hedefler
    for x in ozet:
        if not str(x.get("durum", "")).startswith("TAMAM"):
            uyarilar.append(("Eksik", u'No candidate could be produced',
                             u'%s (%s), members=%s competitors=%s'
                             % (x["hedef"], x["sinif"], x.get("uye", ""),
                                x.get("rakip", "")),
                             u'The member that blocks it: %s. See the subset split sheet.' % (x.get("engelleyen", "") or "no record")))

    # 4. Gevsetilmis kademeyle gecen hedefler
    for x in ozet:
        if x.get("kademe") and x.get("kademe") != "kati":
            uyarilar.append(("Dikkat", u'A candidate was produced with a relaxed step',
                             "%s (%s), kademe=%s" % (x["hedef"], x["sinif"],
                                                     x["kademe"]),
                             u'The strict rule (the orphan primer must not bind anywhere in the competitors) could not be met; at its best placement in the competitors the orphan primer carries at least three mismatches.'))

    # 5. Specific targets giving many off target products in an external database
    for x in gecen:
        try:
            n = int(x.get("dis_hedef_disi", 0) or 0)
        except (TypeError, ValueError):
            continue
        if n > 100 and not x["hedef"].endswith("universal"):
            uyarilar.append(("Dikkat", u'Many off target products in an external database',
                             u'%s (%s): %d products' % (x["hedef"], x["sinif"], n),
                             u'Clean in the sample but it gives a product in the relatives in the reference database. Candidates of the same target with a lower count should be preferred.'))

    i = 4
    for t, konu, ayrinti, yorum in uyarilar:
        for j, v in enumerate((t, konu, ayrinti, yorum), 1):
            c = ws.cell(row=i, column=j, value=v)
            c.font = NORMAL
            c.border = KENAR
            c.alignment = Alignment(wrap_text=True, vertical="top")
        ws.cell(row=i, column=1).fill = (KIRMIZI if t == "Eksik"
                                         else SARI if t == "Dikkat" else YESIL)
        i += 1
    if i == 4:
        ws.cell(row=4, column=1, value="No warning").font = NORMAL

    # ---------------- 5d. The reference design ----------------
    if ref_ham:
        ws = wb.create_sheet("Reference design")
        ws["A1"] = ('Pairs designed from REFERENCE database sequences, for '
                    'targets the sample cannot meet. These pairs ARE NOT '
                    'confirmed against the sample; the columns on the right '
                    'only show whether such a template is present in the '
                    'sample at all. The specificity claim holds within the '
                    'coverage of the reference database.')
        ws["A1"].font = NORMAL
        ws["A1"].alignment = Alignment(wrap_text=True, vertical="top")
        ws.merge_cells("A1:K1")
        ws.row_dimensions[1].height = 60
        yaz_baslik(ws, ["Hedef", "Class", "Database", "Forward primer",
                        "Forward Tm", "Reverse primer", "Reverse Tm", "Product (bp)",
                        "Reads giving a product in the sample", "Wilson alt", "Durum"],
                   satir=3)
        genislik(ws, [36, 7, 20, 28, 9, 28, 9, 11, 24, 12, 24])
        i = 4
        for x in ref_ham:
            urun = ("%s-%s" % (x.get("urun_min", ""), x.get("urun_maks", ""))
                    if x.get("urun_min") != x.get("urun_maks")
                    else x.get("urun_min", ""))
            satir = [x.get("hedef", ""), x.get("sinif", ""),
                     x.get("veritabani", ""), x.get("ileri_dizi", ""),
                     x.get("ileri_tm", ""), x.get("geri_dizi", ""),
                     x.get("geri_tm", ""), urun,
                     x.get("numune_urun_okuma", ""),
                     x.get("numune_wilson_alt", ""), x.get("durum", "")]
            for j, v in enumerate(satir, 1):
                if j in (5, 7, 10):
                    try:
                        v = float(v)
                    except (TypeError, ValueError):
                        pass
                c = ws.cell(row=i, column=j, value=v)
                c.font = Font(name="Consolas", size=10) if j in (4, 6) else NORMAL
                c.border = KENAR
            ws.cell(row=i, column=11).fill = (
                YESIL if x.get("durum") == "numunede_destekli" else SARI)
            i += 1

    # ---------------- 5e. The laboratory protocol ----------------
    ws = wb.create_sheet("Laboratory protocol")
    ws["A1"] = ('A preparation table for ordering and for the first PCR. The '
                'suggested annealing temperature is 5 C below the lower Tm of '
                'the pair, and the gradient range is 4 C either side of it. '
                'The values rest on the Tm computed with primer3 and have to '
                'be confirmed with a gradient on the first run.')
    ws["A1"].font = NORMAL
    ws["A1"].alignment = Alignment(wrap_text=True, vertical="top")
    ws.merge_cells("A1:L1")
    ws.row_dimensions[1].height = 50
    yaz_baslik(ws, ["Karar", "Hedef", "Class", "Oligo name", "Dizi (5'-3')",
                    "Uzunluk", "GC %", "Tm (C)", "Expected product (bp)",
                    "Suggested annealing (C)", "Gradient range (C)", "Kaynak"],
               satir=3)
    genislik(ws, [7, 34, 7, 26, 30, 9, 8, 9, 18, 20, 20, 20])
    i = 4
    protokol_kayit = []
    for k in sorted(en_iyi, key=lambda z: (en_iyi[z].get("karar", ""), z[0], z[1])):
        protokol_kayit.append((en_iyi[k], "numune"))
    for x in ref_ham:
        if x.get("durum") == "numunede_destekli":
            protokol_kayit.append((x, "referans"))
    for x, kaynak in protokol_kayit:
        try:
            tf, tr = float(x.get("ileri_tm", 0)), float(x.get("geri_tm", 0))
        except (TypeError, ValueError):
            tf = tr = 0.0
        dusuk = min(tf, tr) if tf and tr else max(tf, tr)
        ta = round(dusuk - 5.0, 1) if dusuk else ""
        aralik = ("%.1f - %.1f" % (ta - 4, ta + 4)) if ta else ""
        urun = ("%s-%s" % (x.get("urun_min", ""), x.get("urun_maks", ""))
                if x.get("urun_min") != x.get("urun_maks")
                else x.get("urun_min", ""))
        for yon, dizi_a, tm_a, gc_a in (
                ("F", "ileri_dizi", "ileri_tm", "ileri_gc"),
                ("R", "geri_dizi", "geri_tm", "geri_gc")):
            dizi = x.get(dizi_a, "")
            if not dizi:
                continue
            satir = [x.get("karar", ""), x.get("hedef", ""), x.get("sinif", ""),
                     "%s_%s_%s" % (x.get("hedef", "")[:22], x.get("sinif", ""), yon),
                     dizi, len(dizi), x.get(gc_a, ""), x.get(tm_a, ""),
                     urun, ta, aralik, kaynak]
            for j, v in enumerate(satir, 1):
                if j in (6, 7, 8, 10):
                    try:
                        v = float(v)
                    except (TypeError, ValueError):
                        pass
                c = ws.cell(row=i, column=j, value=v)
                c.font = Font(name="Consolas", size=10) if j == 5 else NORMAL
                c.border = KENAR
            ws.cell(row=i, column=12).fill = (
                YESIL if kaynak == "numune" else SARI)
            i += 1
    son_p = i - 1
    i += 1
    ws.cell(row=i, column=1, value="Toplam oligo").font = KALIN
    ws.cell(row=i, column=6, value="=COUNTA(E4:E%d)" % son_p).font = KALIN
    i += 2
    for baslik, metin in [
        ("Reaksiyon", "25 uL: 12,5 uL 2x master mix, her primer 0,2-0,5 uM, "
                      "1-5 ng kalip DNA, kalan hacim nukleaz iceriksiz su."),
        ("Dongu", '95 C for 3 minutes; 35 cycles of 95 C for 15 s, annealing '
                  'for 20 s and 72 C for 20 s; then 72 C for 5 minutes. The '
                  'extension is kept short because the products are under 300 '
                  'bp.'),
        ("Ilk kosu", 'The annealing temperature has to be tested with a '
                     'gradient. The gradient range for each pair is in the '
                     'column above.'),
        ("Kontroller", 'The positive control: the sample where the target has '
                       'the highest read count. The negative control: no '
                       'template. The specificity control: the sample where a '
                       'competitor taxon dominates.'),
        ("Jel", '2 per cent agarose with a 100 bp ladder. A single band is '
                'expected; a second band is a sign of cross amplification and '
                'that pair has to be dropped.'),
        ("Kaynak sutunu", "sample: the primer was designed from the sample's "
                          'own reads and confirmed against the raw reads. '
                          'reference: the primer was designed from a '
                          'reference database and the sample only carries a '
                          'support measurement.'),
    ]:
        ws.cell(row=i, column=1, value=baslik).font = KALIN
        c = ws.cell(row=i, column=2, value=metin)
        c.font = NORMAL
        c.alignment = Alignment(wrap_text=True, vertical="top")
        ws.merge_cells(start_row=i, start_column=2, end_row=i, end_column=12)
        ws.row_dimensions[i].height = max(15, 13 * (1 + len(metin) // 110))
        i += 1

    # ---------------- 6. The elimination reasons ----------------
    ws = wb.create_sheet("Why candidates were dropped")
    yaz_baslik(ws, ["Target and class", "Oligos after composition",
                    "After the thermodynamics", "Binding every member",
                    "Binding no competitor", "Dropped on the pair Tm difference",
                    "No product in one member", "A product in a competitor", "No orphan primer",
                    "Hetero-dimer", "Valid pairs"])
    genislik(ws, [40] + [18] * 10)
    i = 2
    # These patterns read the log that design_group_primers.py writes. When a line
    # there is reworded, the pattern here has to move with it, or the column
    # silently comes out empty.
    kal = [
        (r"oligos after composition filter:\s*(\d+)", 2),
        (r"oligos after thermodynamics:\s*(\d+)", 3),
        (r"oligos binding every target member:\s*(\d+)", 4),
        (r"oligos that bind nowhere in the competitors:\s*(\d+)", 5),
        (r"dropped, pair Tm difference\s*:\s*(\d+)", 6),
        (r"dropped, no product in one member\s*:\s*(\d+)", 7),
        (r"dropped, product forms in a competitor\s*:\s*(\d+)", 8),
        (r"dropped, no orphan primer\s*:\s*(\d+)", 9),
        (r"dropped, hetero-dimer dG\s*:\s*(\d+)", 10),
        (r"valid pairs\s*:\s*(\d+)", 11),
    ]
    for p in sorted(glob.glob(os.path.join(a.candidates, "*__*.log"))):
        t = open(p, encoding="utf-8", errors="replace").read()
        ws.cell(row=i, column=1,
                value=os.path.basename(p)[:-4]).font = NORMAL
        ws.cell(row=i, column=1).border = KENAR
        for kalip, sutun in kal:
            m = re.search(kalip, t)
            c = ws.cell(row=i, column=sutun, value=int(m.group(1)) if m else "")
            c.font = NORMAL
            c.border = KENAR
        i += 1

    wb.save(a.out)
    print(u'written: %s' % a.out)
    print("sayfa: %s" % ", ".join(wb.sheetnames))
    print(u'primer pairs passing: %d, targets covered: %d' % (len(gecen), len(kapsanan)))


if __name__ == "__main__":
    main()
