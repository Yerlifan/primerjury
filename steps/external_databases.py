#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
external_databases.py
Doğrulamadan geçen primer çiftlerini DIŞ referans veritabanlarına karşı
sınar. Numunede olmayan ama ortamda bulunabilecek akrabalarda ürün oluşup
oluşmadığını gösterir.

Yöntem: her primer blastn ile (task blastn-short) veritabanına aranır.
Aynı referans dizide, ters zincirlerde ve 3' uçları birbirine bakan iki
vuruş ürün uzunluğu aralığında buluşuyorsa hedef dışı ürün sayılır.
BLAST tek başına "primer bağlanır mı" sorusunu yanıtlamaz; bu yüzden her
vuruş ayrıca toplantı kararındaki BAĞLANMA KURALIYLA da denetlenir:
son iki baz birebir, son beş bazda en fazla bir uyumsuzluk, toplamda en
fazla üç uyumsuzluk.

Kullanım:
  python3 external_databases.py --final primer_final --db REFERANS_DB \
      --out primer_final/dis_veritabani.tsv
"""
import argparse, csv, os, re, shutil, subprocess, sys, tempfile, collections

TAMLAYICI = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T",
         "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}

# Hangi sinif hangi veritabanina karsi sinanir.
#
# DAR KUME (varsayilan): NCBI RefSeq'in kuratorlu 16S ve ITS setleri.
# Bunlar tip susu ve temsilci dizilerden olusur; kulturlenmemis cevresel
# soylar buyuk olcude YOKTUR. Anaerobik curutucu topluluğunun onemli bir
# kismi tam olarak o gruptadir, dolayisiyla dar kume tek basina "hedef
# disi urun yok" demek icin yeterli degildir.
SINIF_DB = {
    "A1": ["archaea.16S.fna"],
    "A2": ["archaea.16S.fna"],
    "B":  ["bacteria.16S.fna"],
    "F1": ["fungi.ITS.fna", "fungi.28SrRNA.fna"],
    "F2": ["fungi.ITS.fna", "fungi.28SrRNA.fna"],
}

# GENIS KUME (--genis ile): cevresel dizileri de iceren buyuk veritabanlari.
# SILVA SSU/LSU NR99 ve UNITE, RefSeq'te bulunmayan kulturlenmemis soylari
# tasir; ROD tam rRNA operon varyantlarini, PR2 okaryot SSU'sunu kapsar.
# Calisma suresi ciddi olcude uzar, bu yuzden varsayilan degildir.
#
# SINIF BASINA ELLE SECIM YAPILMAZ. Onceki surumde genis kume elle
# sinifa gore yazilmisti ve bu SESSIZ BIR OLCUM KAYBINA yol acti:
#   OLCULDU (2026-08-01): ROD_v1.2_operon_variants.fasta'nin 60320
#   kaydinin 60320'si de Eukaryota'dir; icinde 0 Bacteria, 0 Archaea
#   vardir (basliklardaki soyagacindan sayildi). ROD A1/A2/B siniflarina
#   atanmisti. Sonuc: 71 arke/bakteri cifti icin "ROD'da hedef disi urun
#   yok" yazildi. Bu bir ozgulluk kaniti degildir; veritabaninda o alandan
#   hicbir dizi yoktur. Ayni hatanin tersi de vardi: ROD'un 9753 mantar
#   tam operonu, mantar siniflarinda hic taranmamisti.
#
# Yeni kural: TUM SINIFLAR TUM rDNA VERITABANLARINI GORUR. Gercek PCR'da
# primer yalnizca kendi alanindaki rDNA ile degil, ortamdaki tum DNA ile
# karsilasir; bir arke primerinin bakteri 23S'inde yanlis baglanmasi tam
# olarak aranan hata turudur. Alanla sinirlamak icin bilimsel gerekce yok,
# yalnizca hiz gerekcesi vardi; hiz gerekcesi dar kumede zaten karsilaniyor.
#
# SILVA_138.2_LSUParc.fasta bilerek DISARIDA: Parc kumesi kismi ve dusuk
# kaliteli kayitlari da icerir, LSURef_NR99 ise ayni kapsamin %99
# tekrarsiz kuratorlu temsilcisidir. SSU tarafinda da Parc degil NR99
# kullaniliyor; ikisini ayni bicimde ele almak icin LSU'da da NR99 alinir.
GENIS_ORTAK = [
    "SILVA_138.2_SSURef_NR99.fasta",
    "SILVA_138.2_LSURef_NR99.fasta",
    "ROD_v1.2_operon_variants.fasta",
    "UNITE_ITS.fasta",
    "PR2_SSU_taxo_long.fasta",
    "archaea.16S.fna",
    "bacteria.16S.fna",
    "fungi.ITS.fna",
    "fungi.18SrRNA.fna",
    "fungi.28SrRNA.fna",
]

# ref_all.fna ve ref_all2.fna BILEREK YOK.
#   OLCULDU (2026-08-01, .fai kimlik kumeleri karsilastirilarak):
#   ref_all2.fna = archaea.16S + bacteria.16S + fungi.ITS + fungi.18SrRNA
#                + fungi.28SrRNA, tam olarak 65358 kayit, iki yonde de
#                fark sifir.
#   ref_all.fna  = archaea.16S + bacteria.16S + fungi.ITS, 48431 kayit.
# Yani ikisi de yukaridaki listenin alt kumesidir; taramaya eklemek tek
# bir yeni dizi getirmez, yalnizca ayni kayitlari ikinci kez tarar.
# (Bu yuzden bu iki dosya icin BLAST indeksi kurmaya da gerek yoktur.)

SINIF_DB_GENIS = {
    s: [d for d in GENIS_ORTAK if d not in SINIF_DB[s]] for s in SINIF_DB
}


def rc(s):
    return s.translate(TAMLAYICI)[::-1]


def uyar(p, t):
    return bool(set(IUPAC.get(p, "")) & set(IUPAC.get(t, "")))


def baglanma_uygun(oligo, hedef, son_tam=2, son_pencere=5, son_mm=1, toplam_mm=3):
    """Toplanti karari: uzama 3' ucten basladigi icin son bazlar kritik."""
    if len(oligo) != len(hedef):
        return False, None
    mm = [i for i, (p, t) in enumerate(zip(oligo, hedef)) if not uyar(p, t)]
    n = len(oligo)
    if any(i >= n - son_tam for i in mm):
        return False, None
    if sum(1 for i in mm if i >= n - son_pencere) > son_mm:
        return False, None
    if len(mm) > toplam_mm:
        return False, None
    return True, len(mm)


def db_hazirla(fna, calisma):
    """Var olan BLAST indeksini kullanir, yoksa calisma klasorunde kurar."""
    if os.path.exists(fna + ".nin") or os.path.exists(fna + ".00.nin"):
        return fna
    hedef = os.path.join(calisma, os.path.basename(fna))
    if not os.path.exists(hedef + ".nin"):
        os.symlink(os.path.abspath(fna), hedef) if not os.path.exists(hedef) else None
        r = subprocess.run(["makeblastdb", "-in", hedef, "-dbtype", "nucl"],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print("   makeblastdb basarisiz: %s" % r.stderr.strip()[:200])
            return None
    return hedef


KONS_ETIKET = re.compile(r"((?:A1|A2|B|F1|F2))-\d+_(\d+)")


def sinif_konsensuslari(kons_klasoru, sinif):
    """Verilen sinifin konsensus dosyalarindan dizileri doner."""
    diziler = []
    if not kons_klasoru or not os.path.isdir(kons_klasoru):
        return diziler
    for ad in sorted(os.listdir(kons_klasoru)):
        if not ad.endswith(".fasta"):
            continue
        m = KONS_ETIKET.match(ad)
        if not m or m.group(1) != sinif:
            continue
        ad_, dizi = None, []
        with open(os.path.join(kons_klasoru, ad), encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(">"):
                    if ad_ and dizi:
                        diziler.append((ad_, "".join(dizi)))
                    ad_, dizi = line[1:].strip().split()[0], []
                else:
                    dizi.append(line.strip())
        if ad_ and dizi:
            diziler.append((ad_, "".join(dizi)))
    return diziler


def kapsam_olc(db, diziler, calisma, etiket, threads, zaman_asimi):
    """Bu veritabani bu sinifin organizmalarini HIC iceriyor mu?

    ROD hatasindan sonra eklendi. Bir veritabaninda ilgili alandan tek
    dizi bile yoksa, 'hedef disi urun bulunamadi' sonucu ozgulluk kaniti
    DEGILDIR; olculebilecek bir sey olmadigi icin bos cikmistir. Bu iki
    durum ciktida ayirt edilebilmelidir.

    Olcu: sinifin kendi konsensus dizileri veritabanina megablast ile
    aranir ve EN UZUN hizalama alinir. Esik uydurulmaz, veriden gelir:
    aranan urun en fazla prod_max baz uzunlugunda oldugu icin, veritabanindaki
    en uzun benzer bolge prod_max'tan kisaysa o urunun orada olusmasi
    zaten mekanik olarak imkansizdir; boyle bir veritabani o sinif icin
    KAPSAM_YOK sayilir.

    Doner: (durum, en_uzun_hizalama, en_iyi_kimlik)
    """
    if not diziler:
        return ("KAPSAM_OLCULMEDI", 0, 0.0)
    sorgu = os.path.join(calisma, "kapsam_%s.fa" % etiket)
    with open(sorgu, "w") as fh:
        for ad, dizi in diziler:
            fh.write(">%s\n%s\n" % (ad, dizi))
    cmd = ["blastn", "-query", sorgu, "-db", db,
           "-outfmt", "6 length pident", "-max_target_seqs", "5",
           "-num_threads", str(threads)]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=zaman_asimi)
    except subprocess.TimeoutExpired:
        return ("KAPSAM_OLCULMEDI", 0, 0.0)
    if r.returncode != 0:
        return ("KAPSAM_OLCULMEDI", 0, 0.0)
    en_uzun, en_kimlik = 0, 0.0
    for line in r.stdout.splitlines():
        p = line.split("\t")
        if len(p) < 2:
            continue
        try:
            uz, kim = int(p[0]), float(p[1])
        except ValueError:
            continue
        if uz > en_uzun:
            en_uzun, en_kimlik = uz, kim
    return ("", en_uzun, en_kimlik)


# --- hedefin kendi taksonu ile yabanci takson ayrimi -------------------
#
# 2026-08-01'de olculdu: genis tarama "hedef disi urun" sutununda EN YUKSEK
# sayilari veren kayitlarin bir kismi aslinda HEDEFIN KENDISIDIR.
#   Asetoklastik_metanojenler x archaea.16S = 306 urun; ornek vuruslar
#     NR_104707.1 Methanothrix soehngenii GP6
#     NR_028242.1 Methanothrix soehngenii Opfikon
#   Methanothrix_soehngenii_turu x SILVA = 308 urun; vuruslarin hepsi
#     Methanosaetaceae;Methanothrix
# Bunlar ozgulluk hatasi degil, primerin isini yapmasidir: veritabani hedef
# taksonun dizilerini de icerir. Buna karsilik daha DUSUK sayili bazi
# kayitlar gercek hatadir:
#   Nitrosocosmicus_AOA x SILVA = 1119; vuruslar Halobacteria,
#     Methanoperedenaceae, Thermoplasmata, Cenarchaeum
#   Petrimonas_cinsi x SILVA = 707; vuruslar Clostridium, Bacteroides
# Ham sayiya gore siralamak, kullaniciyi yanlis primerleri duzeltmeye
# gonderiyordu. Bu yuzden her urun, olustugu referansin taksonuna gore
# KENDI TAKSONU ve YABANCI TAKSON diye ikiye ayrilir.
#
# Takson adlari elle yazilmaz; iki veri kaynagindan gelir:
#   1) hedefler.tsv taxid listesi + taxid_adlari.tsv  (BEYAN EDILEN ad)
#   2) hedef_kimlik.tsv olculen_kimlik sutunu         (OLCULEN ad)
# Ikisi ayrisabilir ve ayrisiyor: Trichoderma_cinsi hedefinin olculen
# kimligi Petriella musispora'dir. Her ikisi de "kendi taksonu" sayilir;
# hangisinin tuttugu ciktida gorulur.

KUCUK_TOKEN = 4
ATLA_TOKEN = {"uncultured", "candidatus", "bacterium", "archaeon", "sp",
              "strain", "clone", "isolate", "unidentified", "environmental",
              "samples", "incertae", "sedis", "type", "material", "partial",
              "complete", "sequence", "ribosomal", "gene", "rrna", "genes",
              "fungal", "endophyte", "voucher", "culture", "enrichment"}


def _tokenlar(metin):
    """Herhangi bir baslik bicimini ortak bir kelime kumesine indirger.
    SILVA (;), ROD (|;_), UNITE (k__/p__ ve _), PR2 (|_) ve RefSeq (bosluk)
    bicimlerinin hepsi ayni islemden gecer."""
    return {t for t in re.split(r"[^A-Za-z]+", metin.lower())
            if len(t) >= KUCUK_TOKEN and t not in ATLA_TOKEN}


def _ad_cinsi(ad):
    """'Ca. Nitrosocosmicus hydrocola' -> nitrosocosmicus
       'uncultured Acetobacteroides sp.' -> acetobacteroides"""
    for kelime in re.split(r"[^A-Za-z]+", ad):
        k = kelime.lower()
        if len(k) >= KUCUK_TOKEN and k not in ATLA_TOKEN:
            return k
    return ""


def hedef_taksonlari(hedefler_tsv, adlar_tsv, kimlik_tsv):
    """{hedef: {'beyan': set, 'olculen': set, 'evrensel': bool}}"""
    adlar = {}
    if adlar_tsv and os.path.exists(adlar_tsv):
        for line in open(adlar_tsv, encoding="utf-8"):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0].strip().isdigit():
                adlar[p[0].strip()] = p[1].strip()
    out = {}
    if hedefler_tsv and os.path.exists(hedefler_tsv):
        for line in open(hedefler_tsv, encoding="utf-8"):
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or p[0] == "karar":
                continue
            ad, tidler = p[1], [t.strip() for t in p[3].split(",") if t.strip()]
            # '*A', '*B', '*F' evrensel hedeflerin isaretidir: bu hedeflerde
            # "yabanci takson" kavrami yoktur, cok sayida taksonu birden
            # cogaltmak zaten amactir.
            evrensel = any(t.startswith("*") for t in tidler)
            beyan = set()
            for t in tidler:
                c = _ad_cinsi(adlar.get(t, ""))
                if c:
                    beyan.add(c)
            out[ad] = {"beyan": beyan, "olculen": set(), "evrensel": evrensel}
    if kimlik_tsv and os.path.exists(kimlik_tsv):
        for r in csv.DictReader(open(kimlik_tsv, encoding="utf-8"),
                                delimiter="\t"):
            ad = r.get("hedef")
            if not ad:
                continue
            out.setdefault(ad, {"beyan": set(), "olculen": set(),
                                "evrensel": False})
            for parca in re.split(r"[;,]", r.get("olculen_kimlik") or ""):
                parca = re.sub(r"\(.*?\)", " ", parca)
                c = _ad_cinsi(parca)
                if c:
                    out[ad]["olculen"].add(c)
    return out


# --- yabanci vurusun UZAKLIGI ------------------------------------------
#
# 2026-08-01, ikinci turdan sonra olculdu: "yabanci takson" tek basina da
# yeterli bir olcu degil, cunku hedeflerin bir kismi TEK BIR TAKSON degil
# ISLEVSEL GRUPTUR. Ornekler, gercek ciktidan:
#   Hidrojenotrofik_metanojenler -> yabanci vuruslar Methanobacterium
#     alcaliphilum, Methanosphaera stadtmanae. Bunlar beyan edilen taxid
#     listesinde yok, ama ikisi de hidrojenotrofik metanojendir; hedefin
#     amaci zaten bu guruhu yakalamaktir.
#   Nitrosocosmicus_AOA -> yabanci vuruslar Nitrosotalea, Nitrosopumilus.
#     Ikisi de amonyak oksitleyen arkedir, yani AOA'dir.
#   Zoopagomycota_mantari -> yabanci vurus Piptocephalis moniliformis,
#     ki kendisi Zoopagomycota'dir.
# Buna karsilik:
#   Petrimonas_cinsi -> Flavobacterium, Pdegerlendiricieicola vulgatus
#   Trichoderma_cinsi -> Calonectria, Acremonium, Trichothecium, Aniptodera
# Bunlar gercekten uzak taksonlardir.
#
# Ayrimi elle yazilmis bir "islevsel grup" tablosu ile yapmak, tam da
# kacinilmasi gereken sey olurdu. Bunun yerine uzaklik VERIDEN olculur:
# veritabani basliklari soyagacini tasir (SILVA, ROD, PR2, UNITE). Hedefin
# kendi soyagaci, KENDI TAKSONUNDA urun veren referanslarin ortak on ekinden
# cikarilir; her yabanci vurusun bu soyagaci ile paylastigi derinlik olculur.
# Son iki basamak disinda her seyi paylasan vurus YAKIN, otekiler UZAK
# sayilir. Oncelik siralamasi UZAK sayisina gore yapilir.

def _soyagaci(baslik):
    """Baslktan siralanmis soyagaci alanlarini cikarir. Soyagaci tasimayan
    bicimlerde (RefSeq gibi) bos liste doner."""
    if not baslik:
        return []
    govde = baslik.split(None, 1)
    # SILVA: 'KIMLIK Archaea;Halobacteriota;...'
    if len(govde) > 1 and ";" in govde[1]:
        alanlar = govde[1].split(";")
    elif "|" in baslik:
        # ROD: 'GCA|kaynak|Eukaryota;...;Tur|size=1'
        # PR2: 'KIMLIK|18S_rRNA|nucleus||Eukaryota|TSAR|...'
        # UNITE: 'UDB|k__Fungi;p__...;s__Tur|SH...'
        parcalar = baslik.split("|")
        soy = [p for p in parcalar if ";" in p]
        alanlar = soy[0].split(";") if soy else parcalar[3:]
    else:
        return []
    out = []
    for x in alanlar:
        x = re.sub(r"^[a-z]__", "", x.strip())
        x = re.sub(r"[^A-Za-z ]+", " ", x).strip().lower()
        if x and x not in ("incertae sedis", "unclassified"):
            out.append(x)
    return out


def _ortak_derinlik(a, b):
    n = 0
    for x, y in zip(a, b):
        if x != y:
            break
        n += 1
    return n


def _baskin_soy(soylar):
    """Kendi taksonundaki vuruslarin BASKIN (en sik) soyagaci.

    Ortak on ek ALINMAZ. Alinsaydi esik, kendi vuruslarinin sayisina ve
    cesitligine gore kayardi: tek vurusta on ek tam soyagaci (derinlik 7),
    iki cesitli vurusta yalnizca takim duzeyi (derinlik 4) cikiyor ve ayni
    yabanci vurus bir kayitta UZAK, otekinde YAKIN sayilabiliyordu.
    Baskin soyagaci ise derinligi sabit tutar. Kendi vuruslarinin hepsi
    tanim geregi ayni cinsten oldugu icin bu soyagaclari zaten neredeyse
    ozdestir; en sik olani secmek guvenlidir.
    """
    soylar = [tuple(s) for s in soylar if s]
    if not soylar:
        return []
    return list(collections.Counter(soylar).most_common(1)[0][0])


def basliklari_coz(fna, kimlikler, onbellek):
    """Urun olusturan referanslarin tam basligini fasta'dan tek gecisle alir.
    blastn ciktisi yalniz kisa kimligi verir; takson bilgisi baslikta durur."""
    d = onbellek.setdefault(fna, {})
    eksik = {k for k in kimlikler if k not in d}
    if not eksik or not os.path.exists(fna):
        return d
    with open(fna, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            kimlik = line[1:].split(None, 1)[0]
            if kimlik in eksik:
                d[kimlik] = line[1:].rstrip("\n")
                eksik.discard(kimlik)
                if not eksik:
                    break
    for k in eksik:
        d[k] = ""          # cozulemedi; yabanci sayilmaz, bilinmiyor sayilir
    return d


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--final", required=True, help="09'un output directory")
    p.add_argument("--db", required=True, help="REFERANS_DB directory")
    p.add_argument("--hedefler", default=None, help="targets.tsv")
    p.add_argument("--adlar", default=None, help="taxid_adlari.tsv")
    p.add_argument("--kimlik", default=None,
                   help="hedef_kimlik.tsv from the target-identity step (measured identity)")
    p.add_argument("--kons", default=None,
                   help="consensus directory; if given, every (class, database) "
                        "ikilisi icin KAPSAM DENETIMI yapilir")
    p.add_argument("--out", default=None)
    p.add_argument("--prod-min", type=int, default=50)
    p.add_argument("--prod-max", type=int, default=400)
    p.add_argument("--evalue", type=float, default=1000.0)
    p.add_argument("--max-hedef", type=int, default=5000)
    p.add_argument("--is-parcacigi", type=int, default=4)
    p.add_argument("--genis", action="store_true",
                   help="also the large databases that include environmental sequences "
                        "tarar (SILVA, UNITE, ROD, PR2). Uzun surer.")
    p.add_argument("--yalniz-genis", action="store_true",
                   help="only genis kumeyi tarar")
    p.add_argument("--zaman-asimi", type=int, default=14400,
                   help="veritabani basina saniye siniri")
    return p.parse_args()


def main():
    a = get_args()
    tsv = os.path.join(a.final, "primer_final.tsv")
    if not os.path.exists(tsv):
        sys.exit("bulunamadi: %s" % tsv)
    rows = [r for r in csv.DictReader(open(tsv, encoding="utf-8"), delimiter="\t")
            if r.get("ozgulluk_durum") == "GECTI"]
    if not rows:
        sys.exit("gecen aday yok")
    print(u'pairs to test: %d' % len(rows))

    # sinif -> primer kumesi
    sinif_primer = collections.defaultdict(dict)   # sinif -> dizi -> ad
    for i, r in enumerate(rows):
        sinif_primer[r["sinif"]]["%s_F%d" % (r["hedef"][:24], i)] = r["ileri_dizi"]
        sinif_primer[r["sinif"]]["%s_R%d" % (r["hedef"][:24], i)] = r["geri_dizi"]

    calisma = tempfile.mkdtemp(prefix="blastdb_")
    sonuc = []
    atlanan = []
    gecersiz = []          # (sinif, db, en_uzun_hizalama, kimlik)
    kons_onbellek = {}
    baslik_onbellek = {}
    taksonlar = hedef_taksonlari(a.hedefler, a.adlar, a.kimlik)
    if taksonlar:
        ev = sum(1 for v in taksonlar.values() if v["evrensel"])
        ad = sum(1 for v in taksonlar.values()
                 if not v["evrensel"] and not (v["beyan"] | v["olculen"]))
        print("takson adi cozulen hedef: %d (evrensel: %d, adsiz: %d)"
              % (len(taksonlar), ev, ad))
        if ad:
            print("   UYARI: adsiz hedeflerde kendi/yabanci ayrimi yapilamaz,"
                  " urunler 'takson bilinmiyor' sayilir")
    else:
        print("UYARI: --hedefler/--adlar/--kimlik verilmedi; kendi takson ile"
              " yabanci takson AYRILMAYACAK, ham urun sayisi yaniltici olur")
    for sinif, primerler in sorted(sinif_primer.items()):
        dblist = [] if a.yalniz_genis else list(SINIF_DB.get(sinif, []))
        if a.genis or a.yalniz_genis:
            dblist += SINIF_DB_GENIS.get(sinif, [])
        for dbad in dblist:
            fna = os.path.join(a.db, dbad)
            if not os.path.exists(fna):
                print(u'   no database, skipped: %s' % fna)
                atlanan.append((sinif, dbad, "dosya yok"))
                continue
            db = db_hazirla(fna, calisma)
            if not db:
                atlanan.append((sinif, dbad, "indeks kurulamadi"))
                continue
            # KAPSAM DENETIMI, taramadan once. Bkz. kapsam_olc().
            kap_durum, kap_uz, kap_kim = kapsam_olc(
                db, kons_onbellek.setdefault(
                    sinif, sinif_konsensuslari(a.kons, sinif)),
                calisma, "%s_%s" % (sinif, dbad), a.is_parcacigi, a.zaman_asimi)
            if not kap_durum:
                if kap_uz < a.prod_max:
                    kap_durum = "KAPSAM_YOK"
                    print("      KAPSAM YOK: %s icinde %s sinifina benzer en uzun "
                          "bolge %d bp (%.1f%%), aranan urun en fazla %d bp. "
                          "Bu veritabanindaki 'urun yok' sonucu OZGULLUK KANITI "
                          "DEGILDIR." % (dbad, sinif, kap_uz, kap_kim, a.prod_max))
                    gecersiz.append((sinif, dbad, kap_uz, kap_kim))
                else:
                    kap_durum = "KAPSANIYOR"
            sorgu = os.path.join(calisma, "%s_%s.fa" % (sinif, dbad))
            with open(sorgu, "w") as fh:
                for ad, dizi in primerler.items():
                    fh.write(">%s\n%s\n" % (ad, dizi))
            cikti = sorgu + ".tsv"
            cmd = ["blastn", "-task", "blastn-short", "-query", sorgu, "-db", db,
                   "-outfmt", "6 qseqid sseqid sstart send sstrand length "
                              "qstart qend qlen sseq qseq mismatch",
                   "-evalue", str(a.evalue), "-max_target_seqs", str(a.max_hedef),
                   "-num_threads", str(a.is_parcacigi), "-dust", "no",
                   "-out", cikti]
            print("   blastn %-14s x %-22s (%d primer)"
                  % (sinif, dbad, len(primerler)))
            try:
                r = subprocess.run(cmd, capture_output=True, text=True,
                                   timeout=a.zaman_asimi)
            except subprocess.TimeoutExpired:
                # Sessizce atlamak, o veritabanini "temiz" gostermek olur.
                print("      ZAMAN ASIMI (%d sn): %s OLCULEMEDI"
                      % (a.zaman_asimi, dbad))
                atlanan.append((sinif, dbad, "zaman asimi"))
                continue
            if r.returncode != 0:
                print("      HATA: %s" % r.stderr.strip()[:200])
                atlanan.append((sinif, dbad, "blastn hatasi"))
                continue
            # vuruslari referans basina topla
            vurus = collections.defaultdict(list)
            with open(cikti) as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 12:
                        continue
                    q, s, sst, sen, strand = p[0], p[1], int(p[2]), int(p[3]), p[4]
                    qs, qe, ql = int(p[6]), int(p[7]), int(p[8])
                    sseq, qseq = p[9].upper(), p[10].upper()
                    # yalniz 3' ucu kapsayan vurus uzama yapabilir
                    if qe != ql:
                        continue
                    # BLAST kisa hizalama dondurdugunde primerin 5' tarafi
                    # hizalamanin disinda kalir ve oradaki uyumsuzluklar
                    # sayilmaz. Gorulemeyen kisim bilinmedigi icin en kotu
                    # durum varsayilir: hizalanmayan her 5' bazi olasi bir
                    # uyumsuzluk sayilir ve toplam siniri buna gore denetlenir.
                    gorulmeyen = qs - 1
                    ok, mm = baglanma_uygun(qseq, sseq)
                    if not ok:
                        continue
                    if mm + gorulmeyen > 3:
                        continue
                    mm += gorulmeyen
                    uc = sen        # 3' ucun referanstaki konumu
                    # Primer uzunlugu da tasinir; urun boyu iki 3' ucu
                    # arasindaki mesafe DEGIL, ileri primerin 5' ucundan geri
                    # primerin 5' ucuna kadar olan mesafedir. Eski surumde
                    # aradaki fark (lenF + lenR - 2, tipik 44 baz) yuzunden
                    # gercek 70-94 bp'lik hedef disi urunler alt sinirin
                    # altinda kalip hic sayilmiyordu.
                    vurus[s].append((q, uc, strand, mm, ql))
            # ayni referansta ters zincirlerde ve 3' uclari birbirine bakan cift
            ekstra = collections.Counter()
            detay = {}
            urunler = collections.defaultdict(list)   # cift -> [(ref, boy)]
            for s, v in vurus.items():
                arti = [x for x in v if x[2] == "plus"]
                eksi = [x for x in v if x[2] == "minus"]
                for qf, uf, _, mf, lf in arti:
                    for qr, ur, _, mr, lr in eksi:
                        if ur <= uf:
                            continue
                        # 04'un product_len'i ile ayni olcu: ileri primerin
                        # 5' ucu uf - lf + 1, geri primerin 5' ucu ur + lr - 1
                        boy = (ur + lr - 1) - (uf - lf + 1) + 1
                        if a.prod_min <= boy <= a.prod_max:
                            i1 = re.sub(r"_[FR]\d+$", "", qf)
                            i2 = re.sub(r"_[FR]\d+$", "", qr)
                            n1 = qf.rsplit("_", 1)[1][1:]
                            n2 = qr.rsplit("_", 1)[1][1:]
                            if n1 != n2:
                                continue          # farkli ciftlerin primerleri
                            ekstra[n1] += 1
                            urunler[n1].append((s, boy))
                            detay.setdefault(n1, []).append("%s:%d bp" % (s, boy))

            # urun olusturan referanslarin taksonunu coz ve kendi/yabanci ayir
            tum_ref = {s for lst in urunler.values() for s, _ in lst}
            bas = basliklari_coz(fna, tum_ref, baslik_onbellek) if tum_ref else {}
            ref_token = {s: _tokenlar(bas.get(s, "")) for s in tum_ref}

            for i, r2 in enumerate(rows):
                if r2["sinif"] != sinif:
                    continue
                k = str(i)
                tk = taksonlar.get(r2["hedef"], {})
                evrensel = bool(tk.get("evrensel"))
                kendi_adlar = set(tk.get("beyan", ())) | set(tk.get("olculen", ()))
                kendi = yabanci = bilinmiyor = 0
                yab_ornek = []
                kendi_soylar, yab_kayit = [], []
                for s, boy in urunler.get(k, []):
                    tok = ref_token.get(s, set())
                    if evrensel:
                        kendi += 1
                    elif not tok:
                        bilinmiyor += 1
                    elif kendi_adlar & tok:
                        kendi += 1
                        kendi_soylar.append(_soyagaci(bas.get(s, "")))
                    else:
                        yabanci += 1
                        yab_kayit.append((s, boy))
                        if len(yab_ornek) < 5:
                            yab_ornek.append("%s:%d bp"
                                             % (bas.get(s, s)[:70], boy))
                # yabanci vuruslarin hedefin soyagacina uzakligi
                ref_soy = _baskin_soy(kendi_soylar)
                yakin = uzak = soysuz = 0
                for s, _boy in yab_kayit:
                    soy = _soyagaci(bas.get(s, ""))
                    if not soy or not ref_soy:
                        soysuz += 1
                    elif _ortak_derinlik(ref_soy, soy) >= max(1, len(ref_soy) - 2):
                        yakin += 1
                    else:
                        uzak += 1
                sonuc.append(dict(
                    hedef=r2["hedef"], sinif=sinif, veritabani=dbad,
                    ileri_dizi=r2["ileri_dizi"], geri_dizi=r2["geri_dizi"],
                    hedef_disi_urun=ekstra.get(k, 0),
                    urun_kendi_taksonda=kendi,
                    urun_yabanci_taksonda=yabanci,
                    yabanci_yakin=yakin,
                    yabanci_uzak=uzak,
                    yabanci_soyagacsiz=soysuz,
                    urun_takson_bilinmiyor=bilinmiyor,
                    hedef_soyagaci=";".join(ref_soy),
                    hedef_turu="evrensel" if evrensel else "ozgul",
                    ornekler=";".join(detay.get(k, [])[:5]),
                    yabanci_ornekler=";".join(yab_ornek),
                    kapsam_durumu=kap_durum,
                    kapsam_en_uzun_bp=kap_uz,
                    kapsam_kimlik=("%.1f" % kap_kim) if kap_uz else ""))
    try:
        shutil.rmtree(calisma, ignore_errors=True)
    except Exception:
        pass
    if not sonuc:
        print(u'no database could be scanned; no output was written')
        return
    if a.out:
        d = os.path.dirname(os.path.abspath(a.out))
        if d:
            os.makedirs(d, exist_ok=True)
        with open(a.out, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(sonuc[0].keys()),
                               delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(sonuc)
        print("\nyazildi: %s" % a.out)
    if atlanan:
        print("\nOLCULEMEYEN VERITABANLARI (temiz sayilmazlar):")
        for sn, db, sb in atlanan:
            print("   %-4s %-34s %s" % (sn, db, sb))
    if gecersiz:
        print("\nKAPSAM DISI IKILILER (bu satirlardaki 'urun yok' sonucu"
              " ozgulluk kaniti degildir):")
        for sn, db, uz, kim in gecersiz:
            print("   %-4s %-34s en uzun benzer bolge %5d bp  %%%.1f"
                  % (sn, db, uz, kim))
    temiz = sum(1 for x in sonuc if x["hedef_disi_urun"] == 0)
    temiz_gecerli = sum(1 for x in sonuc if x["hedef_disi_urun"] == 0
                        and x.get("kapsam_durumu") == "KAPSANIYOR")
    print("kayit: %d, hicbir urun vermeyen: %d"
          " (bunlarin kapsami dogrulanmis olani: %d)"
          % (len(sonuc), temiz, temiz_gecerli))

    # SIRALAMA HAM SAYIYA GORE YAPILMAZ. Gerekcesi hedef_taksonlari()'nin
    # ustundeki nottadir: en yuksek ham sayilarin bir kismi hedefin kendi
    # taksonudur ve hata degildir.
    ozgul = [x for x in sonuc
             if x.get("hedef_turu") == "ozgul"
             and x.get("kapsam_durumu") == "KAPSANIYOR"
             and x.get("yabanci_uzak", 0) > 0]
    print("\nOZGULLUK BULGULARI, UZAK TAKSONA gore siralanmis")
    print("(yakin = hedefin soyagacini son iki basamak disinda paylasan,"
          " cogu zaman ayni islevsel grup)")
    if not ozgul:
        print("   uzak taksonda urun veren kayit yok")
    for x in sorted(ozgul, key=lambda x: -x["yabanci_uzak"])[:15]:
        print("   %-28s %-3s %-28s uzak=%5d yakin=%5d kendi=%5d"
              % (x["hedef"][:27], x["sinif"], x["veritabani"][:27],
                 x["yabanci_uzak"], x["yabanci_yakin"], x["urun_kendi_taksonda"]))
        if x["yabanci_ornekler"]:
            print("        %s" % x["yabanci_ornekler"].split(";")[0][:100])

    # cift duzeyinde ozet
    cift = {}
    for x in sonuc:
        if x.get("kapsam_durumu") != "KAPSANIYOR":
            continue
        k = (x["hedef"], x["ileri_dizi"], x["geri_dizi"])
        d = cift.setdefault(k, {"uzak": 0, "yakin": 0, "kendi": 0,
                                "turu": x.get("hedef_turu")})
        d["uzak"] += x.get("yabanci_uzak", 0)
        d["yakin"] += x.get("yabanci_yakin", 0)
        d["kendi"] += x.get("urun_kendi_taksonda", 0)
    oz = [v for v in cift.values() if v["turu"] == "ozgul"]
    # KENDI TAKSONUNDA HIC URUN VERMEYEN cift, "temiz" sayilamaz. Bu, ROD
    # hatasinin cift duzeyindeki karsiligidir: olculebilecek bir sey
    # olmadigi icin bos cikmistir. Olculdu (2026-08-01): Sakarolitik F2'nin
    # uc cifti sekiz kapsanan veritabaninda hicbir urun vermiyor, kendi
    # taksonunda bile; bunlar zaten 18'in alan_karisimi diye isaretledigi
    # ciftlerdir.
    olculebilir = [v for v in oz if v["kendi"] > 0]
    print("\nkapsanan veritabaniyla olculen OZGUL cift: %d" % len(oz))
    print("   bunlarin kendi taksonunda urun vererek OLCULEBILIR olani: %d"
          % len(olculebilir))
    print("   uzak taksonda hic urun vermeyen (gercekten temiz)      : %d"
          % sum(1 for v in olculebilir if v["uzak"] == 0))
    inert = [v for v in oz if v["kendi"] == 0]
    if inert:
        print("   kendi taksonunda bile urun vermeyen (OLCUM GECERSIZ)   : %d"
              % len(inert))


if __name__ == "__main__":
    main()
