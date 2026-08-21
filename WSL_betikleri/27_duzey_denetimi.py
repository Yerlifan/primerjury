#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
27_duzey_denetimi.py
TOPLANTI KARARININ DOĞRUDAN SINANMASI: her hedef, kararın istediği düzeyde
(tür ya da cins) ayrım yapıyor mu?

Karar:
  tür özgül  : Methanosarcina mazei, Methanothrix soehngenii,
               Methanosarcina barkeri, Podospora pseudopauciseta,
               Dictyostelium discoideum, Trichoderma asperellum
  cins özgül : Bacteroides, Alistipes, Proteiniphilum, Petrimonas
Bu liste elle yazılmaz; hedefler.tsv'nin "duzey" sütunundan okunur.

NEDEN AYRI BİR ÖLÇÜM GEREKİYOR
09 numunedeki rakiplere karşı sınıyor, 14 dış veritabanlarında hedef dışı
ürün arıyor. İkisi de "bu çift kendi cinsinin ÖTEKİ TÜRLERİNDEN ayırıyor mu"
sorusunu sormuyor. Tür özgüllüğü tam olarak bu sorudur ve ancak kardeş
türlerden kurulu bir panele karşı yanıtlanabilir.

YÖNTEM
Her hedefin cinsi, beyan edilen taxid adından çıkarılır. Referans
veritabanları taranıp o cinse ait, TÜR ADI BELLİ olan bütün kayıtlar
toplanır ve bir panel oluşturulur. Panele karşı blastn koşulur, ürünler
14'ün ürün kuralıyla (aynı referans, ters zincirler, 3' uçları karşı
karşıya, ürün boyu aralıkta) ve aynı bağlanma kuralıyla sayılır. Ürün
oluşan kayıtlar türe göre ayrılır.

Panelde YALNIZ tür adı belli kayıtlar kullanılır. SILVA gibi "uncultured
archaeon" ağırlıklı çevresel veritabanları tür ayrımı paneline giremez;
tür kimliği taşımayan bir kayıt, tür ayrımını ne doğrular ne çürütür.

KARAR
  duzey=tur : hedef türde ürün var ve aynı cinsin hiçbir başka türünde ürün
              yoksa TUR_OZGUL. Toplantı kararı 1-2 çapraz türü hoş gördüğü
              için, çapraz TÜR sayısı eşiğin altındaysa TUR_OZGUL_ESIKLI
              (eşik --capraz-tur-esik ile değişir, varsayılan 2). Eşiğin
              üstündeyse TUR_AYRIMI_YOK. Ölçü çapraz TÜR sayısıdır, o
              türlerde oluşan ürün sayısı değil.
              Hedef tür panelde hiç yoksa HEDEF_TUR_PANELDE_YOK; bu durumda
              tür özgüllüğü gösterilemez, çürütülemez de.
  duzey=cins: cins içi kapsama bildirilir (kaç tür çoğaltılıyor). Cins
              DIŞINDA ürün olup olmadığı 14'ün işidir, burada tekrarlanmaz.

Kullanım:
  python3 27_duzey_denetimi.py --hedefler hedefler.tsv \
      --adlar taxid_adlari.tsv --final primer_final --db REFERANS_DB \
      --kimlik primer_final/hedef_kimlik.tsv \
      --out primer_final/duzey_denetimi.tsv
"""
import argparse, collections, csv, os, re, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
# Bagimlilik BILEREK: baglanma ve urun kurali 14 ile BIREBIR ayni olmali.
# Ayri bir kopya yazmak, iki olcumun sessizce ayrismasina kapi acar.
import importlib.util as _il
_s = _il.spec_from_file_location("_dv", os.path.join(HERE, "14_dis_veritabani.py"))
DV = _il.module_from_spec(_s)
_yedek, sys.argv = sys.argv, ["14_dis_veritabani.py"]
try:
    _s.loader.exec_module(DV)
except SystemExit:
    pass
finally:
    sys.argv = _yedek

# Tur adi tasiyan veritabanlari. SILVA disarida: NR99 kayitlarinin buyuk
# cogunlugu 'uncultured archaeon/bacterium' ile biter ve tur kimligi
# tasimaz. Panele alinsalardi, tur ayrimi olcumu tur adi olmayan
# kayitlarla seyreltilirdi.
PANEL_DB = [
    "archaea.16S.fna",
    "bacteria.16S.fna",
    "fungi.ITS.fna",
    "fungi.18SrRNA.fna",
    "fungi.28SrRNA.fna",
    "UNITE_ITS.fasta",
    "ROD_v1.2_operon_variants.fasta",
    "PR2_SSU_taxo_long.fasta",
]

# 'sp', 'spp', 'cf', 'aff' TUR ADI DEGILDIR.
#   OLCULDU (2026-08-01): UNITE'te 's__Trichoderma_sp' bicimindeki kayitlar
#   tur adi sayiliyordu ve panelde ayri bir "tur" gibi duruyordu; yalniz
#   Trichoderma icin 16910, Marasmius icin 1712, Podospora icin 1326 kayit.
#   Tur adi belli olmayan bir kayitta urun olusmasi, tur ayriminin
#   basarisiz oldugunu GOSTERMEZ; hangi tur oldugu bilinmiyor. Bunlar
#   panele girerse hem panel tur sayisi sisiyor hem de yanlis
#   TUR_AYRIMI_YOK karari uretilebiliyor.
#   (RefSeq bicimindeki 'Trichoderma sp.' zaten noktali oldugu icin
#   duzenli ifadeye takilmiyordu; kacan yalniz alt cizgili bicimlerdi.)
GURULTU = {"uncultured", "unidentified", "environmental", "sample", "clone",
           "isolate", "strain", "voucher", "candidatus", "bacterium",
           "archaeon", "fungal", "endophyte", "culture", "enrichment",
           "sp", "spp", "cf", "aff", "indet", "incertae", "sedis"}


def tur_adi(baslik):
    """Baslktan 'Cins tur' ikilisini cikarir; yoksa bos doner.

    RefSeq : 'NR_104707.1 Methanothrix soehngenii GP6 16S ...'
    UNITE  : '...;s__Thelephora_albomarginata|SH...'
    ROD    : '...;Drosophila;Drosophila_melanogaster|size=1'
    PR2    : '...|Unruhdinium|Unruhdinium_kevei'
    SILVA  : '... ;Methanothrix;uncultured archaeon'  -> bos
    """
    if not baslik:
        return ""
    # once alt cizgili ikili bicimler (UNITE, ROD, PR2)
    for parca in re.split(r"[|;]", baslik):
        p = parca.strip()
        p = re.sub(r"^[a-z]__", "", p)
        m = re.match(r"^([A-Z][a-z]+)_([a-z][a-z\-]+)$", p)
        if (m and m.group(1).lower() not in GURULTU
                and m.group(2).lower() not in GURULTU):
            return "%s %s" % (m.group(1), m.group(2))
    # sonra bosluklu RefSeq bicimi: kimlikten sonraki iki kelime
    #
    # BASLIKTA NOKTALI VIRGUL OLMASI TEK BASINA ELEME SEBEBI DEGIL.
    #   OLCULDU (2026-08-01): RefSeq ITS kayitlari
    #   'NR_172285.1 Petriella musispora CBS 745.69 ITS region; from TYPE
    #   material' bicimindedir. Onceki surumde baslikta noktali virgul
    #   varsa bu dal hic calismiyordu, dolayisiyla fungi.ITS.fna ve
    #   fungi.28SrRNA.fna'daki TIP MATERYALI kayitlarinin tur adi
    #   cikarilamiyor ve panele alinmiyorlardi. Panel tam da bu kayitlara
    #   en cok ihtiyac duyulan yerdi: tur ayriminin altin standardi tip
    #   susu dizileridir. (Olculdu: yalniz Petriella icin 35 kayit boyle
    #   dusuyordu.)
    # SILVA'nin soyagacli basliklari bu daldan zaten gecemez, cunku
    # ikinci kelime noktali virgul icerir ve asagidaki duzenli ifadeye
    # takilir; ayri bir korumaya gerek yok.
    kelime = baslik.split()
    if len(kelime) >= 3:
        g, t = kelime[1], kelime[2]
        if (re.match(r"^[A-Z][a-z]+$", g) and re.match(r"^[a-z][a-z\-]+$", t)
                and g.lower() not in GURULTU and t.lower() not in GURULTU):
            return "%s %s" % (g, t)
    return ""


def ad_parcala(ad):
    """'Ca. Nitrosocosmicus hydrocola' -> ('Nitrosocosmicus','hydrocola')"""
    kelime = [k for k in re.split(r"[^A-Za-z]+", ad or "") if k]
    kelime = [k for k in kelime if k.lower() not in GURULTU and len(k) > 2]
    if not kelime:
        return ("", "")
    cins = kelime[0]
    tur = kelime[1].lower() if len(kelime) > 1 else ""
    return (cins, tur)


def referans_esle(ref_hedef, hedef_adlari):
    """primer_referans.tsv'deki hedef adini hedefler.tsv adina baglar.

    SADECE '_referans' EKINI SOYMAK YETMIYOR.
      OLCULDU (2026-08-01): 'Methanosarcina_barkeri_referans' -> soyulunca
      'Methanosarcina_barkeri' cikiyor, oysa hedefler.tsv'deki ad
      'Methanosarcina_barkeri_turu'. Eslesme tutmayinca hedefin TEK
      primer takimi (de novo hicbir cifti yok) sessizce dusuyor ve hedef
      CIFT_YOK gorunuyordu. Sessiz dusme, olcumu yapilmamis bir hedefi
      "cift bulunamadi" gibi gosterir.
    Bu yuzden once birebir, sonra on ek eslesmesi denenir; hicbiri
    tutmazsa None doner ve cagiran taraf bunu YUKSEK SESLE bildirir.
    """
    kok = re.sub(r"_referans$", "", ref_hedef)
    if kok in hedef_adlari:
        return kok
    adaylar = [h for h in hedef_adlari if h.startswith(kok + "_")]
    if len(adaylar) == 1:
        return adaylar[0]
    return None


def hedefleri_oku(hedefler_tsv, adlar_tsv, kimlik_tsv):
    adlar = {}
    if adlar_tsv and os.path.exists(adlar_tsv):
        for line in open(adlar_tsv, encoding="utf-8"):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0].strip().isdigit():
                adlar[p[0].strip()] = p[1].strip()
    olculen = {}
    if kimlik_tsv and os.path.exists(kimlik_tsv):
        for r in csv.DictReader(open(kimlik_tsv, encoding="utf-8"),
                                delimiter="\t"):
            olculen[r["hedef"]] = re.sub(r"\(.*?\)", "",
                                         r.get("olculen_kimlik") or "").strip()
    out = []
    for line in open(hedefler_tsv, encoding="utf-8"):
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 4 or p[0] == "karar":
            continue
        duzey = p[2].strip()
        if duzey not in ("tur", "cins"):
            continue
        tidler = [t.strip() for t in p[3].split(",") if t.strip()]
        beyan = [adlar.get(t, "") for t in tidler if adlar.get(t)]
        cinsler, turler = set(), set()
        for ad in beyan:
            c, t = ad_parcala(ad)
            if c:
                cinsler.add(c)
            if c and t:
                turler.add("%s %s" % (c, t))
        # ISTEGE BAGLI 7. SUTUN: hedef_tur
        # Bazi hedeflerde beyan edilen taxid, numunede bulunan organizmaya
        # karsilik gelmiyor ve dogru turun bizim taxid_adlari.tsv'mizde
        # karsiligi olmayabiliyor. Ornek: kutulari 101201 (Trichoderma
        # asperellum) diye etiketli olan hedefin olculen kimligi Petriella
        # musispora'dir. Taxid UYDURULMAZ; tur adi bu sutunda dogrudan
        # yazilir ve hedefin kendi ad kumesine eklenir. Sutun bossa hicbir
        # sey degismez, eski davranis surer.
        # Sutun EKLEMEZ, DEGISTIRIR. Konulma sebebi taxid'in adinin yanlis
        # olmasidir; adi eklemek yanlis adi da tutmak olurdu. Olculdu:
        # Petriella_musispora satirinda in_taxid, kutularin Kraken2
        # etiketleri olan 101201/2034170/63577'dir (numune destegi bunlarla
        # olculur). Ekleme yapilsaydi hedef tur kumesine Trichoderma
        # asperellum, atroviride ve breve de girerdi ve Trichoderma
        # cogaltan bir cift "hedef turde urun var" sayilirdi.
        if len(p) > 6 and p[6].strip():
            cinsler, turler = set(), set()
            for parca in p[6].split(","):
                c2, t2 = ad_parcala(parca)
                if c2:
                    cinsler.add(c2)
                if c2 and t2:
                    turler.add("%s %s" % (c2, t2))
        olc = olculen.get(p[1], "")
        oc, ot = ad_parcala(olc.split(";")[0]) if olc else ("", "")
        out.append({"hedef": p[1], "duzey": duzey, "beyan_ad": "; ".join(beyan),
                    "cinsler": cinsler, "hedef_turler": turler,
                    "olculen_ad": olc, "olculen_cins": oc,
                    "olculen_tur": ("%s %s" % (oc, ot)) if oc and ot else ""})
    return out


def panelleri_topla(dbklasor, cinsler, en_fazla_tur_basina, gunluk):
    """Veritabanlarini tek gecisle tarar; cins basina tur etiketli kayitlar.
    Doner: {cins: {tur: [(etiket, dizi), ...]}}"""
    panel = {c: collections.defaultdict(list) for c in cinsler}
    kirpilan = collections.Counter()
    dusen_tursuz = collections.Counter()
    dusen_baska_cins = collections.Counter()
    kucuk = {c.lower() for c in cinsler}
    for dbad in PANEL_DB:
        yol = os.path.join(dbklasor, dbad)
        if not os.path.exists(yol):
            gunluk.append("panel: %s yok, atlandi" % dbad)
            continue
        alinan = 0
        with open(yol, encoding="utf-8", errors="replace") as fh:
            baslik, dizi, sec = None, [], None
            for line in fh:
                if line.startswith(">"):
                    if sec and dizi:
                        c, t = sec
                        if len(panel[c][t]) < en_fazla_tur_basina:
                            panel[c][t].append(
                                ("%s|%s" % (dbad.split(".")[0], baslik[:60]),
                                 "".join(dizi)))
                            alinan += 1
                        else:
                            kirpilan[(c, t)] += 1
                    baslik = line[1:].rstrip("\n")
                    dizi, sec = [], None
                    dusuk = baslik.lower()
                    for c in cinsler:
                        if c.lower() in dusuk:
                            ta = tur_adi(baslik)
                            if ta and ta.split()[0].lower() == c.lower():
                                sec = (c, ta)
                            elif ta:
                                # Cins adi baslikta geciyor ama ikili ad
                                # baska bir cinse ait (su adi, notu, konak
                                # bilgisi olabilir). Panele alinmaz; SESSIZ
                                # DUSMEZ, sayilir ve bildirilir.
                                dusen_baska_cins[c] += 1
                            else:
                                dusen_tursuz[c] += 1
                            break
                elif sec is not None:
                    dizi.append(line.strip())
            if sec and dizi:
                c, t = sec
                if len(panel[c][t]) < en_fazla_tur_basina:
                    panel[c][t].append(
                        ("%s|%s" % (dbad.split(".")[0], baslik[:60]),
                         "".join(dizi)))
                    alinan += 1
                else:
                    kirpilan[(c, t)] += 1
        gunluk.append("panel: %-34s %d kayit alindi" % (dbad, alinan))
    # SESSIZ KIRPMA YOK: kirpilan ve tur adi olmadigi icin dusen kayitlar
    # bildirilir, yoksa panel eksikligi tam kapsama gibi okunur.
    for (c, t), n in sorted(kirpilan.items(), key=lambda x: -x[1])[:10]:
        gunluk.append("panel KIRPILDI: %s / %s icin %d kayit alinmadi"
                      % (c, t, n))
    for c, n in dusen_tursuz.items():
        if n:
            gunluk.append("panel: %s icin tur adi tasimayan %d kayit "
                          "panele ALINMADI" % (c, n))
    for c, n in dusen_baska_cins.items():
        if n:
            gunluk.append("panel: %s adi geciyor ama ikili ad baska cinse "
                          "ait olan %d kayit panele ALINMADI" % (c, n))
    return panel


def urun_say(primerler, panel_fa, calisma, etiket, a):
    """Panele karsi blastn; 14'un baglanma ve urun kuralini uygular.
    Doner: {cift_no: {referans: urun_sayisi}}"""
    db = os.path.join(calisma, "panel_%s.fa" % etiket)
    shutil.copyfile(panel_fa, db)
    r = subprocess.run(["makeblastdb", "-in", db, "-dbtype", "nucl"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return None
    sorgu = os.path.join(calisma, "sorgu_%s.fa" % etiket)
    with open(sorgu, "w") as fh:
        for ad, dizi in primerler.items():
            fh.write(">%s\n%s\n" % (ad, dizi))
    cikti = sorgu + ".tsv"
    cmd = ["blastn", "-task", "blastn-short", "-query", sorgu, "-db", db,
           "-outfmt", "6 qseqid sseqid sstart send sstrand length "
                      "qstart qend qlen sseq qseq mismatch",
           "-evalue", str(a.evalue), "-max_target_seqs", "100000",
           "-num_threads", str(a.is_parcacigi), "-dust", "no", "-out", cikti]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        return None
    vurus = collections.defaultdict(list)
    with open(cikti) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            q, s, strand = p[0], p[1], p[4]
            qs, qe, ql = int(p[6]), int(p[7]), int(p[8])
            sen = int(p[3])
            sseq, qseq = p[9].upper(), p[10].upper()
            if qe != ql:
                continue
            gorulmeyen = qs - 1
            ok, mm = DV.baglanma_uygun(qseq, sseq)
            if not ok or mm + gorulmeyen > 3:
                continue
            vurus[s].append((q, sen, strand, ql))
    say = collections.defaultdict(collections.Counter)
    for s, v in vurus.items():
        arti = [x for x in v if x[2] == "plus"]
        eksi = [x for x in v if x[2] == "minus"]
        for qf, uf, _, lf in arti:
            for qr, ur, _, lr in eksi:
                if ur <= uf:
                    continue
                boy = (ur + lr - 1) - (uf - lf + 1) + 1
                if not (a.prod_min <= boy <= a.prod_max):
                    continue
                n1 = qf.rsplit("_", 1)[1][1:]
                n2 = qr.rsplit("_", 1)[1][1:]
                if n1 == n2:
                    say[n1][s] += 1
    return say


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--hedefler", required=True)
    p.add_argument("--adlar", required=True)
    p.add_argument("--final", required=True, help="primer_final klasoru")
    p.add_argument("--referans", default=None, help="primer_referans.tsv")
    p.add_argument("--db", required=True, help="REFERANS_DB klasoru")
    p.add_argument("--kimlik", default=None, help="hedef_kimlik.tsv")
    p.add_argument("--out", required=True)
    p.add_argument("--prod-min", type=int, default=50)
    p.add_argument("--prod-max", type=int, default=400)
    p.add_argument("--evalue", type=float, default=1000.0)
    p.add_argument("--is-parcacigi", type=int, default=4)
    p.add_argument("--tur-basina-en-fazla", type=int, default=200)
    # Toplanti karari: "1-2 capraz tur olursa onlarda tur ozgul sayilir".
    # Deger burada sabit yazilmaz, secenek olarak durur; degistirilirse
    # ciktinin basinda hangi esikle calisildigi yaziliyor.
    p.add_argument("--capraz-tur-esik", type=int, default=2,
                   help="tur ozgullugunde hosgorulen capraz TUR sayisi "
                        "(urun sayisi degil); varsayilan 2")
    return p.parse_args()


def main():
    a = get_args()
    gunluk = []
    hedefler = hedefleri_oku(a.hedefler, a.adlar, a.kimlik)
    if not hedefler:
        sys.exit("hedefler.tsv icinde duzey=tur/cins satiri yok")
    print("karar duzeyi belirtilen hedef: %d (tur: %d, cins: %d)"
          % (len(hedefler),
             sum(1 for h in hedefler if h["duzey"] == "tur"),
             sum(1 for h in hedefler if h["duzey"] == "cins")))
    print("tur ozgullugunde hosgorulen capraz TUR sayisi: %d"
          % a.capraz_tur_esik)

    # ciftleri topla
    ciftler = collections.defaultdict(list)
    tsv = os.path.join(a.final, "primer_final.tsv")
    for r in csv.DictReader(open(tsv, encoding="utf-8"), delimiter="\t"):
        if r.get("ozgulluk_durum") == "GECTI":
            ciftler[r["hedef"]].append((r["ileri_dizi"], r["geri_dizi"],
                                        "de novo"))
    if a.referans and os.path.exists(a.referans):
        # hedefler.tsv'deki BUTUN adlar (duzey ayrimi yapmadan), cunku
        # referans setinde duzey=grup hedefleri de var
        tum_ad = []
        for line in open(a.hedefler, encoding="utf-8"):
            if line.startswith("#") or not line.strip():
                continue
            pp = line.rstrip("\n").split("\t")
            if len(pp) >= 2 and pp[0] != "karar":
                tum_ad.append(pp[1])
        eslesmeyen = collections.Counter()
        for r in csv.DictReader(open(a.referans, encoding="utf-8"),
                                delimiter="\t"):
            ham = r.get("hedef", "")
            ad = referans_esle(ham, tum_ad)
            if ad is None:
                eslesmeyen[ham] += 1
                continue
            ciftler[ad].append((r["ileri_dizi"], r["geri_dizi"], "referans"))
        for ham, n in eslesmeyen.items():
            print("   UYARI: referans hedefi '%s' (%d cift) hedefler.tsv'deki"
                  " hicbir adla eslesmedi, OLCUM DISI KALDI" % (ham, n))

    tum_cins = set()
    for h in hedefler:
        tum_cins |= h["cinsler"]
        if h["olculen_cins"]:
            tum_cins.add(h["olculen_cins"])
    print("panel kurulacak cins: %s" % ", ".join(sorted(tum_cins)))
    panel = panelleri_topla(a.db, tum_cins, a.tur_basina_en_fazla, gunluk)
    for g in gunluk:
        print("   " + g)

    calisma = tempfile.mkdtemp(prefix="duzey_")
    sonuc = []
    for h in hedefler:
        hd = h["hedef"]
        cf = ciftler.get(hd, [])
        # panel: beyan edilen cins(ler) + olculen cins
        kendi_cins = set(h["cinsler"])
        if h["olculen_cins"]:
            kendi_cins.add(h["olculen_cins"])
        kayitlar = []
        for c in sorted(kendi_cins):
            for t, lst in sorted(panel.get(c, {}).items()):
                for i, (etiket, dizi) in enumerate(lst):
                    kayitlar.append(("%s__%s__%d" % (c, t.replace(" ", "_"), i),
                                     t, dizi))
        panel_turler = sorted({t for _, t, _ in kayitlar})
        if not cf:
            sonuc.append(dict(hedef=hd, duzey=h["duzey"], cift_no="",
                              kaynak="", beyan_ad=h["beyan_ad"],
                              olculen_ad=h["olculen_ad"],
                              panel_tur_sayisi=len(panel_turler),
                              hedef_turde_urun="", diger_turde_urun="",
                              cogaltilan_hedef_turler="",
                              capraz_tur_sayisi="", cogaltilan_turler="",
                              karar="CIFT_YOK"))
            continue
        if not kayitlar:
            for i, (F, R, kaynak) in enumerate(cf):
                sonuc.append(dict(hedef=hd, duzey=h["duzey"], cift_no=str(i),
                                  kaynak=kaynak, beyan_ad=h["beyan_ad"],
                                  olculen_ad=h["olculen_ad"],
                                  panel_tur_sayisi=0, hedef_turde_urun="",
                                  diger_turde_urun="",
                                  cogaltilan_hedef_turler="",
                                  capraz_tur_sayisi="",
                                  cogaltilan_turler="", karar="PANEL_YOK"))
            continue
        panel_fa = os.path.join(calisma, "%s_panel.fa" % hd)
        with open(panel_fa, "w") as fh:
            for kimlik, _t, dizi in kayitlar:
                fh.write(">%s\n%s\n" % (kimlik, dizi))
        primerler = {}
        for i, (F, R, _k) in enumerate(cf):
            primerler["%s_F%d" % (hd[:20], i)] = F
            primerler["%s_R%d" % (hd[:20], i)] = R
        print("   %-32s %2d cift x %3d panel kaydi (%d tur)"
              % (hd, len(cf), len(kayitlar), len(panel_turler)))
        say = urun_say(primerler, panel_fa, calisma, hd[:16], a)
        if say is None:
            print("      blastn/makeblastdb basarisiz, atlandi")
            continue
        kimlik_tur = {k: t for k, t, _ in kayitlar}
        # hedef tur kumesi: beyan edilen tur(ler) + olculen tur
        hedef_turler = set(h["hedef_turler"])
        if h["olculen_tur"]:
            hedef_turler.add(h["olculen_tur"])
        for i, (F, R, kaynak) in enumerate(cf):
            per = say.get(str(i), collections.Counter())
            tur_urun = collections.Counter()
            for ref, n in per.items():
                tur_urun[kimlik_tur.get(ref, "?")] += n
            # HANGI hedef turde urun olustugu yazilir. Hedef tur kumesi
            # beyan edilen adi VE olculen kimligi birlikte iceriyor; ikisi
            # ayrisabilir. Ornek: Zoopagomycota_mantari'nin beyan edilen
            # turu Dictyostelium discoideum, olculen kimligi ise esik alti
            # bir Marasmius. Hangisinde urun olustugu yazilmazsa,
            # TUR_OZGUL_ESIKLI karari "Dictyostelium'a ozgul" diye
            # okunabilir; oysa Marasmius'a ozgul olabilir.
            cogaltilan_hedef = sorted(t for t in tur_urun if t in hedef_turler)
            hedefte = sum(n for t, n in tur_urun.items() if t in hedef_turler)
            digerde = sum(n for t, n in tur_urun.items()
                          if t not in hedef_turler)
            digerler = sorted({t for t in tur_urun if t not in hedef_turler})
            hedef_panelde = bool(hedef_turler & set(panel_turler))
            if h["duzey"] == "tur":
                # CAPRAZ TUR SAYISI, urun sayisi degil. Toplanti karari
                # "1-2 capraz tur olursa yine tur ozgul sayilir" diyor;
                # olcu KAC FARKLI TUR cogaltildigidir, o turlerde kac urun
                # olustugu degil. Sifir capraz ile esik ici durum ayri
                # kararlar olarak yazilir; ikisini tek etikete katlamak,
                # daha zayif olan cifti daha guclusuyle esitlerdi.
                if not hedef_panelde:
                    karar = "HEDEF_TUR_PANELDE_YOK"
                elif hedefte == 0:
                    karar = "HEDEF_TURDE_URUN_YOK"
                elif len(digerler) == 0:
                    karar = "TUR_OZGUL"
                elif len(digerler) <= a.capraz_tur_esik:
                    karar = "TUR_OZGUL_ESIKLI"
                else:
                    karar = "TUR_AYRIMI_YOK"
            else:
                # CINS DUZEYI. Onceki surum yalnizca "panelin kac turunde
                # urun olusuyor" diye sayiyordu; bu, cins DISINDA urun
                # olusmasini gizliyordu.
                #   OLCULDU (2026-08-01): Proteiniphilum_cinsi'nin bes
                #   ciftinden ikisi Fermentimonas caenicola'yi da
                #   cogaltiyor, ki bu baska bir cinstir. Eski sayim ikisini
                #   de CINS_ICI_3_3 diye yaziyordu, yani cins ozgullugunu
                #   ihlal eden cift, en genis kapsamli cift gibi
                #   gorunuyordu.
                # Kabul olcutu "cins ozgul" oldugu icin BEYAN EDILEN cinsin
                # disinda urun olmamasi gerekir. Olculen kimlik burada
                # hedef sayilmaz: Proteiniphilum istenmisse Fermentimonas
                # capraz cogaltmadir, kutulardaki organizma o olsa bile.
                tum_cogaltilan = sorted(t for t in tur_urun if t)
                ici = [t for t in tum_cogaltilan
                       if t.split()[0] in h["cinsler"]]
                disi = [t for t in tum_cogaltilan
                        if t.split()[0] not in h["cinsler"]]
                if not ici:
                    karar = "CINS_ICINDE_URUN_YOK"
                elif not disi:
                    karar = "CINS_OZGUL"
                else:
                    karar = "CINS_AYRIMI_YOK"
                cogaltilan_hedef = ici
                digerler = disi
                hedefte = sum(n for t, n in tur_urun.items() if t in ici)
                digerde = sum(n for t, n in tur_urun.items() if t in disi)
                karar += "_%d_%d" % (len(ici), len(panel_turler))
            sonuc.append(dict(
                hedef=hd, duzey=h["duzey"], cift_no=str(i), kaynak=kaynak,
                beyan_ad=h["beyan_ad"], olculen_ad=h["olculen_ad"],
                panel_tur_sayisi=len(panel_turler),
                hedef_turde_urun=hedefte, diger_turde_urun=digerde,
                cogaltilan_hedef_turler="; ".join(cogaltilan_hedef),
                capraz_tur_sayisi=len(digerler),
                cogaltilan_turler="; ".join(digerler[:8]), karar=karar))
    shutil.rmtree(calisma, ignore_errors=True)

    d = os.path.dirname(os.path.abspath(a.out))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sonuc[0].keys()),
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(sonuc)
    print("\nyazildi: %s" % a.out)

    print("\n%-32s %-5s %-6s %-26s %s"
          % ("hedef", "duzey", "cift", "karar", "cogaltilan diger turler"))
    print("-" * 118)
    for h in hedefler:
        satir = [x for x in sonuc if x["hedef"] == h["hedef"]]
        # Tur listesi KARAR BASINA toplanir. Hedef genelinde toplanirsa,
        # ayrimi basaran TUR_OZGUL satirinin yanina baska bir ciftin
        # cogalttigi turler yaziliyor ve o cift de ayrim yapmiyor gibi
        # okunuyordu.
        gruplu = collections.defaultdict(list)
        for x in satir:
            gruplu[x["karar"]].append(x)
        for karar, lst in sorted(gruplu.items(), key=lambda x: -len(x[1])):
            digerler = sorted({d for x in lst for d in
                               x["cogaltilan_turler"].split("; ") if d})
            hedefteki = sorted({d for x in lst for d in
                                x["cogaltilan_hedef_turler"].split("; ") if d})
            print("%-32s %-5s %-6d %-22s %s"
                  % (h["hedef"][:31], h["duzey"], len(lst), karar,
                     ("; ".join(digerler))[:40]))
            if hedefteki:
                print("%-32s %s"
                      % ("", "  hedef turde cogaltilan: "
                         + "; ".join(hedefteki)[:60]))
    tur_hedef = {h["hedef"] for h in hedefler if h["duzey"] == "tur"}
    kati = {x["hedef"] for x in sonuc if x["karar"] == "TUR_OZGUL"}
    esikli = {x["hedef"] for x in sonuc
              if x["karar"] in ("TUR_OZGUL", "TUR_OZGUL_ESIKLI")}
    print("\ntur ozgullugu istenen hedef: %d" % len(tur_hedef))
    print("   en az bir cifti CAPRAZSIZ (TUR_OZGUL)          : %d"
          % len(kati & tur_hedef))
    print("   en az bir cifti esik ici (<= %d capraz tur)     : %d"
          % (a.capraz_tur_esik, len(esikli & tur_hedef)))


if __name__ == "__main__":
    main()
