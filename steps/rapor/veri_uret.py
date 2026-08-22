#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
veri_uret.py
Rapor verisini (veri.json) primer_final.tsv'den ELLE DEGIL, veriden uretir.
Alan karisimi nedeniyle teslimden cikarilan ciftler burada da cikarilir;
rapor ile Excel ayni kumeyi anlatmak zorundadir.

Kullanim:
  python3 rapor/veri_uret.py --final <FINAL> --consensus <KONS> \
      --targets targets.tsv --names taxid_names.tsv --out rapor/veri.json
"""
import argparse, csv, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import field_audit


def vir(x, n=1):
    try:
        return ("%%.%df" % n % float(x)).replace(".", ",")
    except (TypeError, ValueError):
        return str(x)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--final", required=True)
    p.add_argument("--consensus", required=True)
    p.add_argument("--targets", default="targets.tsv")
    p.add_argument("--external", default=None)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    tsv = os.path.join(a.final, "primer_final.tsv")
    rows = list(csv.DictReader(open(tsv, encoding="utf-8"), delimiter="\t"))
    toplam = len(rows)
    gecen = [r for r in rows if r.get("ozgulluk_durum") == "GECTI"]

    ta = field_audit.taxid_alanlari(a.consensus)
    ht = field_audit.hedef_taxidleri(a.targets)
    cikarilan = []
    temiz = []
    for r in gecen:
        uyumsuz, dag, bas = field_audit.alan_dagilimi(
            r.get("hedef", ""), r.get("sinif", ""), taxid_alan=ta, hedef_taxid=ht)
        (cikarilan if uyumsuz else temiz).append(r)
    gecen = temiz

    dis_yolu = a.external or os.path.join(a.final, "dis_veritabani.tsv")
    dis = {}
    if os.path.exists(dis_yolu):
        for x in csv.DictReader(open(dis_yolu, encoding="utf-8"), delimiter="\t"):
            k = (x["hedef"], x["sinif"], x["ileri_dizi"], x["geri_dizi"])
            dis[k] = dis.get(k, 0) + int(x.get("hedef_disi_urun", 0) or 0)

    # hedef ve sinif basina tek cift: once dis veritabani hedef disi urun,
    # sonra rakip Wilson alt siniri, sonra tasarim cezasi
    en_iyi = {}
    for r in gecen:
        k = (r["hedef"], r["sinif"])
        d = dis.get((r["hedef"], r["sinif"], r["ileri_dizi"], r["geri_dizi"]), 0)
        anahtar = (d, float(r.get("rakip_wilson") or 0), float(r.get("ceza") or 0))
        if k not in en_iyi or anahtar < en_iyi[k][0]:
            en_iyi[k] = (anahtar, r)

    onerilen = []
    for (hedef, sinif), (_, r) in sorted(en_iyi.items()):
        onerilen.append(dict(
            karar=r.get("karar", ""),
            hedef=hedef.replace("_", " ").replace(" turu", "").replace(" cinsi", ""),
            sinif=sinif,
            F=r["ileri_dizi"], Ftm=vir(r["ileri_tm"]),
            R=r["geri_dizi"], Rtm=vir(r["geri_tm"]),
            urun=r.get("urun_min", ""),
            uye="%s/%s" % (r.get("uye_dogrulanan", ""), r.get("uye_toplam", "")),
            w=vir(r.get("rakip_wilson"), 4)))

    veri = dict(
        toplam=toplam,
        gecen=len(gecen),
        hedef=len(set(r["hedef"] for r in gecen)),
        alan_disi=len(cikarilan),
        alan_disi_hedef=sorted(set(r["hedef"] for r in cikarilan)),
        onerilen=onerilen)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(veri, fh, ensure_ascii=False, indent=1)
    print("yazildi: %s" % a.out)
    print(u'  candidates tested=%d  passing=%d  targets covered=%d  recommended=%d'
          % (toplam, len(gecen), veri["hedef"], len(onerilen)))
    if cikarilan:
        print("  alan karisimi nedeniyle cikarilan: %d (%s)"
              % (len(cikarilan), ", ".join(veri["alan_disi_hedef"])))


if __name__ == "__main__":
    main()
