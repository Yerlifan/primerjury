#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
design_from_reference.py
Numunede karşılanamayan hedefler için REFERANS VERİTABANI dizilerinden
primer tasarlar ve tasarladığı çifti numunenin ham okumalarına karşı sınar.

Neden ayrı bir betik: numunedeki kutular birbirinden ayrılamadığı için bazı
hedeflerde özgül çift bulunamıyor. Referans veritabanındaki adlandırılmış
türler ayrılabilir; oradan tasarlanan primer bilimsel olarak doğrudur ama bu
numuneyle DOĞRULANAMAZ. Bu ayrım çıktıda açıkça taşınır: her satır
`referanstan_tasarlandi` etiketi ve numunede ölçülen destek oranıyla gelir.

Girdi tablosu (--hedefler-ref), sekmeyle ayrılmış:
  ad            çıktı etiketi
  sinif         A1, A2, B, F1, F2   (numune desteği bu sınıfın okumalarında ölçülür)
  veritabani    REFERANS_DB içindeki dosya adı
  ic            hedef tür/cins adları, virgülle
  dis           rakip adları, virgülle
  taxid         numune desteği ölçülürken kullanılacak taxid'ler, virgülle

Kullanım:
  python3 design_from_reference.py --db REFERANS_DB --pt . \
      --hedefler-ref hedefler_referans.tsv --out primer_referans
"""
import argparse, csv, glob, importlib.util, os, re, subprocess, sys, tempfile, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
MOTOR = os.path.join(HERE, "design_group_primers.py")


def yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    m = importlib.util.module_from_spec(spec)
    yedek, sys.argv = sys.argv, [yol]
    try:
        spec.loader.exec_module(m)
    except SystemExit:
        pass
    finally:
        sys.argv = yedek
    return m


TS = lambda: datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
_LOG = [None]


def log(msg):
    line = "[%s] %s" % (TS(), msg)
    print(line, flush=True)
    if _LOG[0]:
        with open(_LOG[0], "a", encoding="utf-8") as fh:
            fh.write(line + "\n")


def fasta_oku(p):
    ad, buf = None, []
    for l in open(p, encoding="utf-8", errors="replace"):
        if l.startswith(">"):
            if ad:
                yield ad, "".join(buf)
            ad, buf = l[1:].strip(), []
        else:
            buf.append(l.strip())
    if ad:
        yield ad, "".join(buf)


def _ad_kalibi(ad):
    """Ad SOZCUK SINIRIYLA aranir, alt dize olarak DEGIL.

    Alt dize aramasi sessiz bir kirlenme uretiyordu: ic='Bacteroides'
    yazildiginda 'Parabacteroides' ve 'Acetobacteroides' kayitlari da
    hedef UYESI sayiliyordu. bacteria.16S.fna'da sozcuk sinirli
    'Bacteroides' 86 kayit verirken, alt dize aramasi bunlara 17
    Parabacteroides ile 1 Acetobacteroides kaydini da katiyordu; primer
    o zaman baska cinslerde de urun vermek zorunda kalir ve cins
    ozgullugu daha tasarim aninda kaybedilirdi.
    """
    return re.compile(r"(?<![A-Za-z])%s(?![A-Za-z])" % re.escape(ad),
                      re.IGNORECASE)


def sec(veritabani, adlar, azami_kayit=6):
    """Basliginda verilen adlardan biri SOZCUK OLARAK gecen kayitlari
    toplar. {ad: [(baslik, dizi), ...]} doner.
    veritabani tek yol ya da yol listesi olabilir."""
    yollar = veritabani if isinstance(veritabani, (list, tuple)) \
        else [veritabani]
    bulunan = {a: [] for a in adlar}
    kalip = {a: _ad_kalibi(a) for a in adlar if a}
    for yol in yollar:
        if not os.path.exists(yol):
            continue
        for baslik, dizi in fasta_oku(yol):
            for a in adlar:
                if (a and len(bulunan[a]) < azami_kayit
                        and kalip[a].search(baslik)):
                    bulunan[a].append((baslik,
                                       dizi.upper().replace("U", "T")))
    return bulunan


# TASARIM ANINDAKI RAKIP KUMESI, DOGRULAMA PANELINDEN DAR OLMAMALI.
#   OLCULDU (2026-08-01): Podospora_pseudopauciseta_referans tasarimi 29
#   rakip diziyle yapildi; rakipler yalnizca fungi.ITS.fna'dan ve ad basina
#   en fazla 6 kayit alinarak toplanmisti (o dosyada toplam 14 Podospora
#   kaydi var). 27'nin dogrulama paneli ise UNITE dahil 242 kayit ve 50
#   turdu. Tasarim P. anserina ile P. comata'yi disladigini sanirken,
#   genis panelde ikisini de cogaltiyordu. Yani basarisizligin sebebi
#   tasarim motoru degil, tasarim aninda GORULMEYEN rakiplerdi.
#
# Cozum iki parcali:
#   1) veritabani sutunu virgulle birden cok dosya kabul eder.
#   2) Rakip listesi elle yazilanla sinirli kalmaz: hedefin CINSINDEKI
#      butun oteki turler veriden bulunup rakip kumesine eklenir. Elle
#      yazilan liste, veriden turetilen kumeye EKtir, onun yerine gecmez.
#
# Tur adi tanimi 27'den alinir. Ayri bir kopya yazilsaydi, tasarim ile
# dogrulama farkli tur tanimlariyla calisir ve tam da duzeltilen hata
# baska bicimde geri gelirdi.
_s27 = importlib.util.spec_from_file_location(
    "_dz", os.path.join(HERE, "check_taxonomic_level.py"))
DZ = importlib.util.module_from_spec(_s27)
_yedek27, sys.argv = sys.argv, ["check_taxonomic_level.py"]
try:
    _s27.loader.exec_module(DZ)
except SystemExit:
    pass
finally:
    sys.argv = _yedek27


def kardes_turleri_bul(yollar, cinsler, hedef_turler, azami_tur,
                       azami_kayit_tur):
    """Hedefin cinsindeki OTEKI turleri veriden bulur.

    Doner: ({tur: [(baslik, dizi), ...]}, kirpilan_tur_sayisi)
    """
    bulunan = {}
    kirpilan = 0
    for yol in yollar:
        if not os.path.exists(yol):
            continue
        for baslik, dizi in fasta_oku(yol):
            dusuk = baslik.lower()
            if not any(c.lower() in dusuk for c in cinsler):
                continue
            tur = DZ.tur_adi(baslik)
            if not tur or tur in hedef_turler:
                continue
            if tur.split()[0] not in cinsler:
                continue
            if tur not in bulunan:
                if len(bulunan) >= azami_tur:
                    kirpilan += 1
                    continue
                bulunan[tur] = []
            if len(bulunan[tur]) < azami_kayit_tur:
                bulunan[tur].append((baslik, dizi.upper().replace("U", "T")))
    return bulunan, kirpilan


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True, help="REFERANS_DB klasoru")
    p.add_argument("--pt", required=True, help="PrimerTasarlama kok klasoru")
    p.add_argument("--hedefler-ref", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--max-okuma", type=int, default=3000)
    p.add_argument("--top", type=int, default=5)
    p.add_argument("--extra", default="")
    p.add_argument("--azami-uye", type=int, default=6,
                   help="hedef ad basina alinacak en fazla referans dizi")
    p.add_argument("--kardes-rakip", action="store_true", default=True,
                   help="hedefin cinsindeki oteki turleri VERIDEN bulup "
                        "rakip kumesine ekler (varsayilan acik)")
    p.add_argument("--kardes-rakip-kapali", dest="kardes_rakip",
                   action="store_false",
                   help="kardes tur rakiplerini kapatir (eski davranis)")
    p.add_argument("--azami-kardes-tur", type=int, default=60,
                   help="rakip olarak alinacak en fazla kardes TUR sayisi")
    p.add_argument("--azami-kardes-kayit", type=int, default=2,
                   help="her kardes tur icin en fazla dizi")
    return p.parse_args()


def main():
    a = get_args()
    os.makedirs(a.out, exist_ok=True)
    _LOG[0] = os.path.join(a.out, "referans_tasarim.log")
    G = yukle("g04", MOTOR)
    O = yukle("o09", os.path.join(HERE, "specificity.py"))

    hedefler = []
    with open(a.hedefler_ref, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if p[0] == "ad":
                continue
            hedefler.append(dict(ad=p[0], sinif=p[1], veritabani=p[2],
                                 ic=p[3], dis=p[4] if len(p) > 4 else "",
                                 taxid=p[5] if len(p) > 5 else ""))
    log("islenecek referans hedefi: %d" % len(hedefler))

    fq = {}
    for p in glob.glob(os.path.join(a.pt, "fastq files", "*", "*.fastq")):
        grp = os.path.basename(os.path.dirname(p))
        m = re.search(r"reads[-_](\d+)", os.path.basename(p))
        if m:
            fq[(re.split(r"[-_]", grp)[0], grp, m.group(1))] = p
    log("fastq envanteri: %d" % len(fq))

    calisma = tempfile.mkdtemp(prefix="refdiz_")
    sonuc = []
    for h in hedefler:
        # veritabani sutunu virgulle birden cok dosya alabilir
        vtler = [os.path.join(a.db, x.strip())
                 for x in h["veritabani"].split(",") if x.strip()]
        var_olan = [v for v in vtler if os.path.exists(v)]
        if not var_olan:
            log("ATLANDI %-34s veritabani yok: %s" % (h["ad"], h["veritabani"]))
            continue
        for v in vtler:
            if v not in var_olan:
                log("   %-34s veritabani bulunamadi, atlandi: %s"
                    % (h["ad"], os.path.basename(v)))
        ic_adlar = [x.strip() for x in h["ic"].split(",") if x.strip()]
        dis_adlar = [x.strip() for x in h["dis"].split(",") if x.strip()]
        bul = sec(var_olan, ic_adlar + dis_adlar, a.azami_uye)
        eksik = [x for x in ic_adlar if not bul.get(x)]
        if eksik:
            log("ATLANDI %-34s veritabaninda bulunamayan hedef adi: %s"
                % (h["ad"], ", ".join(eksik)))
            continue
        yok_dis = [x for x in dis_adlar if not bul.get(x)]
        if yok_dis:
            log("   %-34s veritabaninda bulunamayan RAKIP adi: %s"
                % (h["ad"], ", ".join(yok_dis)))

        # VERIDEN TURETILEN KARDES TUR RAKIPLERI
        if a.kardes_rakip:
            cinsler = {ad.split()[0] for ad in ic_adlar if ad.split()}
            hedef_turler = {ad for ad in ic_adlar if len(ad.split()) >= 2}
            kardes, kirpilan_tur = kardes_turleri_bul(
                var_olan, cinsler, hedef_turler, a.azami_kardes_tur,
                a.azami_kardes_kayit)
            # Elle yazilan rakipler KORUNUR, kardes turler UZERINE eklenir.
            for tur, kayitlar in kardes.items():
                if tur not in bul or not bul[tur]:
                    bul[tur] = kayitlar
                    dis_adlar.append(tur)
            log("   %-34s kardes tur rakibi: %d tur, %d dizi%s"
                % (h["ad"], len(kardes),
                   sum(len(v) for v in kardes.values()),
                   (" (KIRPILDI: %d tur alinmadi)" % kirpilan_tur)
                   if kirpilan_tur else ""))

        def yaz(adlar, onek):
            yollar = []
            for ad in adlar:
                for i, (baslik, dizi) in enumerate(bul.get(ad, [])):
                    yol = os.path.join(calisma, "%s_%s_%s_%d_consensus.fasta"
                                       % (h["ad"], onek,
                                          re.sub(r"\W+", "", ad)[:24], i))
                    with open(yol, "w", encoding="utf-8") as fh2:
                        fh2.write(">%s\n" % baslik.replace(" ", "_")[:80])
                        for k in range(0, len(dizi), 70):
                            fh2.write(dizi[k:k + 70] + "\n")
                    yollar.append(yol)
            return yollar

        ing, outg = yaz(ic_adlar, "in"), yaz(dis_adlar, "out")
        tsv = os.path.join(a.out, "%s.tsv" % h["ad"])
        cmd = [sys.executable, MOTOR, "--in-group"] + ing
        if outg:
            cmd += ["--out-group"] + outg
        cmd += ["--label", h["ad"], "--out", tsv]
        if a.extra:
            cmd += a.extra.split()
        r = subprocess.run(cmd, capture_output=True, text=True)
        txt = r.stdout + r.stderr
        with open(os.path.join(a.out, "%s.log" % h["ad"]), "w",
                  encoding="utf-8") as lf:
            lf.write(txt)
        m = re.search(r"gecerli cift sayisi\s*:\s*(\d+)", txt)
        n = int(m.group(1)) if m else 0
        log("%-34s %-3s referans_uye=%d referans_rakip=%d cift=%d"
            % (h["ad"], h["sinif"], len(ing), len(outg), n))
        if not n or not os.path.exists(tsv):
            continue

        # Numune destegi: tasarlanan cift, hedefin numunedeki okumalarinda
        # urun veriyor mu. Bu bir OZGULLUK kaniti degil, yalnizca "numunede
        # boyle bir kalip var mi" sorusunun cevabidir.
        taxidler = [x.strip() for x in h["taxid"].split(",") if x.strip()]
        uye_fq = [v for (s, g, t), v in fq.items()
                  if s == h["sinif"] and t in taxidler]
        rows = list(csv.DictReader(open(tsv, encoding="utf-8"), delimiter="\t"))
        rows.sort(key=lambda x: float(x.get("ceza", 9e9)))
        for r2 in rows[:a.top]:
            F, R = r2["ileri_dizi"], r2["geri_dizi"]
            destek, toplam_ok, urun_ok = 0, 0, 0
            for p in uye_fq:
                tot, fh_, rh, both = O.okuma_taramasi(p, F, R, 50, 400,
                                                      a.max_okuma)
                toplam_ok += tot
                urun_ok += both
                if both > 0:
                    destek += 1
            oran = urun_ok / toplam_ok if toplam_ok else 0.0
            w = O.wilson_alt(urun_ok, toplam_ok) if toplam_ok else 0.0
            sonuc.append(dict(
                hedef=h["ad"], sinif=h["sinif"], kaynak="referans_veritabani",
                veritabani=h["veritabani"],
                ileri_dizi=F, ileri_tm=r2.get("ileri_tm", ""),
                geri_dizi=R, geri_tm=r2.get("geri_tm", ""),
                tm_farki=r2.get("tm_farki", ""),
                urun_min=r2.get("urun_min", ""), urun_maks=r2.get("urun_maks", ""),
                yetim_primer=r2.get("yetim_primer", ""),
                ceza=r2.get("ceza", ""),
                numune_dosya_destegi="%d/%d" % (destek, len(uye_fq)),
                numune_urun_okuma="%d/%d" % (urun_ok, toplam_ok),
                numune_urun_orani=round(oran, 5),
                numune_wilson_alt=round(w, 5),
                durum=("numunede_destekli" if w > 0.01 else
                       "numuneden_dogrulanamadi")))
        log("      numune destegi olculdu: %d aday" % min(len(rows), a.top))

    if sonuc:
        yol = os.path.join(a.out, "primer_referans.tsv")
        with open(yol, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(sonuc[0].keys()),
                               delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(sonuc)
        log("yazildi: %s  (%d satir)" % (yol, len(sonuc)))
        d = sum(1 for x in sonuc if x["durum"] == "numunede_destekli")
        log("numunede destekli: %d, numuneden dogrulanamayan: %d"
            % (d, len(sonuc) - d))
    else:
        log("hicbir referans hedefi icin cift uretilemedi")


if __name__ == "__main__":
    main()
