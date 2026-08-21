#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NUMUNE ICI IN SILICO PCR  (surum 3)

NE YAPAR
Elimizdeki nanopore okumalarini gercek bir PCR kalibi gibi kullanir. Her primer cifti
icin, numunedeki her organizmanin okumalarinda urun olusup olusmadigini sayar. Hedefte
yuksek, rakiplerde sifira yakin olmasi beklenir. Karar rakip oranina gore verilir.

NEDEN UC SURUM OLDU, HER SURUMDE NE YANLISTI
Bu dosyanin gecmisi, ayni sorunun uc farkli sekilde yanlis cozulmesidir. Uc hatanin
UCU DE ayni yone calisiyordu ya da en tehlikeli yone: kotu bir primeri temiz gostermek.

  Surum 0, elle yazilmis eslestirme
    Primerin son sekiz bazinin BIREBIR eslesmesi tohum olarak araniyordu. Rakipte o sekiz
    bazin icinde tek bir fark varsa arama o bolgeyi hic gormuyordu. Gercek capraz
    reaksiyonlar kaciriliyordu.

  Surum 1, blastn'e devir
    Eslestirme blastn'e verildi ama karar kurallari eksik yazildi. Uc ucun hedefe uymasi
    HIC kontrol edilmiyordu, benzerlik yuzdesi okunup kullanilmiyordu, iki primerin
    birbirine bakmasi aranmiyordu. Bu sefer ters yone hata yapiliyordu: olmayan capraz
    reaksiyonlar var gorunuyor, iyi primerler eleniyordu.

  Surum 2, fiziksel kurallar eklendi ama BLAST'in davranisi hesaba katilmadi
    Kural olarak "hizalama primerin 3' ucunu icermeli" kondu. Ama blastn-short en yuksek
    PUANLI hizalamayi bildirir, hizalamanin tamamini degil. Uyumsuzluga eksi uc, uyuma
    arti bir puan verir. Primerin sondan ucuncu bazinda bir uyumsuzluk varsa, son uc bazi
    hizalamaya katmak puani dusurur, BLAST da onlari atar ve q1-17 diye bildirir.
    Gercek durum: son iki baz tam uyuyor, ortada tek fark var, primer baglanir ve uzar.
    Kod ise "3' uc hizalamada yok" deyip baglanmayi yok sayiyordu.
    Ornek, Podospora cifti kendi hedefinde:
        okuma  ACTCGTCGAAGGAGCTTTAC
        primer ACTCGTCGAAGGAGCTTCAC
    BLAST bunu q1-17 bildirdi, surum 2 "baglanma yok" dedi. Yanlis.

SURUM 3'UN COZUMU, KABA KUVVET
Veritabanimiz 23.6 milyon baz. BLAST'in var olma sebebi milyarlarca bazlik veritabanlarinda
kaba kuvvetten kacinmaktir; bu olcekte boyle bir zorunluluk yok. Bu yuzden artik her okuma
uzerinde HER pozisyon tek tek denenir. Tohum yok, puanlama yok, kirpma yok, kacirma yok.
Karsilastirma numpy ile vektorel yapildigi icin butun numune yaklasik yarim dakikada taranir.
blastn ayri bir gorusu olarak korunur (--blast secenegi) ama karar vermez.

URUN KURALI, PCR'in fizigi
  Her primer icin
    a. son iki baz hedefe birebir uymali (uzama oradan baslar)
    b. son bes bazda en fazla bir uyumsuzluk
    c. primerin tamaminda en fazla uc uyumsuzluk
  Iki primer birlikte
    d. biri arti zincire, digeri eksi zincire baglanmali
    e. 3' uclari BIRBIRINE bakmali
    f. aradaki mesafe urun araliginda olmali

KARAR
Oranlar okuma sayisina bolunerek bulunur ama ham oran kucuk orneklemde yaniltir: on yedi
okumasi olan bir taksonda tek okuma yuzde alti eder. Bu yuzden kararda Wilson alt siniri
kullanilir, yani oranin istatistiksel olarak guvenilen en dusuk degeri.
Iki ayri eksen raporlanir:
  CAPRAZ   rakiplerde urun olusuyor mu
  KAPSAMA  hedeflenen organizmalarin hepsinde urun olusuyor mu
Kapsama degerlendirilirken her taksonun "ulasilabilirligi" kendi verisinden olculur: o
takson herhangi bir primer ciftiyle yuksek oran veriyorsa, okumalari saglamdir; ayni
taksonda baska bir cift dusuk oran veriyorsa bu okuma eksikligi degil, o ciftin hatasidir.

Calistirma:
  python3 blast_ispcr.py --root /path/to/project
  python3 blast_ispcr.py --selftest      (butun kurallarin bilinen cevapli sinavi)
  python3 blast_ispcr.py --root ... --blast   (a second opinion via blastn, slow)
"""
import argparse, csv, glob, math, os, random, re, shutil, statistics, subprocess, sys, tempfile, time
from collections import defaultdict

try:
    import numpy as np
except ImportError:
    sys.exit("gerekli: pip install numpy")

RC = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
def rc(s): return "".join(RC.get(b,'N') for b in reversed(s.upper()))
def kod(s): return np.frombuffer(s.encode(), dtype=np.uint8)

LOGF = None
def log(m):
    s = f"[{time.strftime('%H:%M:%S')}] {m}"
    print(s, flush=True)
    if LOGF: LOGF.write(s + "\n"); LOGF.flush()

def arac_var(ad): return shutil.which(ad) is not None
def guvenli(ad): return re.sub(r"[^A-Za-z0-9]+", "_", ad).strip("_")

AYIRAC = 0          # okumalar arasina konan, hicbir baza esit olmayan bayt

# Hangi primer ciftinin hedefi hangi taksondur. Liste verilenler grup primeridir,
# listedeki HER uyede urun beklenir.
HEDEF_TAXID = {
    "AD01_mazei":              ["2209"],
    "AD01b_Mmazei_tur":        ["2209"],
    "AD02_soehngenii":         ["2223"],
    "AD05_Bthetaiotaomicron":  ["818"],
    "AD06_Alistipes":          ["214856"],
    "AD10_Methanosarcinaceae": ["2209","2208","2210","3078083","1434102"],
    "AD11_Methanomicrobiaceae":["83986","118126","394967","2201"],
    "AD12_Tulosesus_callinus": ["181762"],
    "AD13_Oceanobacillus":     ["2954507"],
    "AD14_asetoklastik":       ["2209","2208","2210","3078083","1434102","2223"],
    "AD15_Methanoculleus":     ["83986","118126","394967"],
    "AD16_Microascaceae":      ["101201"],
    "AD17_Mbarkeri":           ["2208"],
    "AD18_Zoopagomycota":      ["44689"],
    "Methanosarcinaceae_ailesi":  ["2209","2208","2210","3078083","1434102"],
    "Methanomicrobiaceae_ailesi": ["83986","118126","394967","2201"],
    "Dysgonomonadaceae_ailesi":   ["2829812","1642647","1642646","285070"],
    "Bacteroidaceae_ailesi":      ["818","28116"],
    "Asetoklastik_metanojenler_asetattan_metan": ["2209","2208","2210","3078083","1434102","2223"],
    # 224719 (Methanobrevibacter sp. AbM4) 28 Temmuz'da eklendi. Numunede 4723 okuma,
    # tartismasiz H2/CO2 metanojeni, ama uye listesinde yoktu; grup bu yuzden eksik
    # bir soruya "tam kapsama" cevabi veriyordu. tasarim asamasiyla ayni olmali.
    "Hidrojenotrofik_metanojenler_H2_CO2_den_metan": ["118126","83986","394967","2201","83984","224719"],
    "Metilotrofik_metanojen_metil_bilesiklerinden_metan": ["1406512"],
    "Sakarolitik_bakteriler_seker_parcalayan": ["818","28116","214856"],
    "Podospora_pseudopauciseta": ["2093780"],
    # 28 Temmuz'da eklendi. Toplantinin "proteolitik / sintrofik" talebinin numunedeki
    # karsiligi. Iki ayri filum, tek ciftle kapsanmasi beklenmiyor; olcum bunu
    # soylesin diye grup TANIMLI hale getirildi. tasarim asamasiyla ayni olmali.
    "Proteolitik_sintrofik_bakteriler": ["456827","1197717"],
}

ISIMLER = {
    "44689":"Zoopagomycota mantari","2209":"Methanosarcina mazei","40559":"Helotiales askomikot",
    "83984":"Methanocorpusculum labreanum","2208":"Methanosarcina barkeri",
    "2223":"Methanothrix soehngenii","1826872":"Ca. Nitrosocosmicus hydrocola",
    "2740404":"Ca. Sulfurimonas baltica","456827":"Ca. Cloacimonas acidaminovorans",
    "101201":"Microascaceae askomikot","2829812":"Proteiniphilum propionicum",
    "1642646":"Petrimonas mucosa","1760811":"Pseudodifflugia amip","1129264":"Sphaerochaeta associata",
    "1642647":"Proteiniphilum saccharofermentans","500148":"Cyrtohymena siliat",
    "2233851":"Blastochloris tepida","2954507":"Oceanobacillus sp.","285070":"Petrimonas sulfuriphila",
    "224719":"Methanobrevibacter sp. AbM4","80884":"Colletotrichum higginsianum",
    "394967":"Methanoculleus receptaculi","1197717":"Cloacibacillus porcorum",
    "3078083":"Methanosarcina hadiensis","1159556":"Ustilaginoidea virens",
    "2093779":"Podospora pseudocomata","5855":"Vampyrellid amip","214856":"Alistipes finegoldii",
    "181762":"Tulosesus callinus","63577":"Trichoderma atroviride",
    "1406512":"Ca. Methanomassiliicoccus intestinalis","456999":"Rhizoctonia solani",
    "2210":"Methanosarcina thermophila","818":"Bacteroides thetaiotaomicron",
    "2034170":"Trichoderma breve","2093780":"Podospora pseudopauciseta",
    "2545709":"Schizosaccharomyces osmophilus","5811":"Toxoplasma gondii",
    "4896":"Schizosaccharomyces pombe","83986":"Methanoculleus bourgensis",
    "118126":"Methanoculleus chikugoensis","28116":"Bacteroides ovatus",
    "1434102":"Methanosarcina sp. WH1","2201":"Methanofollis liminatans",
}
def isim(tx): return ISIMLER.get(tx, f"taxid {tx}")

# ------------------------------------------------------------------ veri yukleme
def okumalari_yukle(kok, okuma_basina=300, en_az_uzunluk=300):
    """Her taksondan sabit sayida okuma alir. Ayni tohum kullanildigi icin tekrarlanabilir."""
    dosyalar = defaultdict(list)
    for p in glob.glob(f"{kok}/SONUCLAR/fastq files/*/*reads_*.fastq"):
        m = re.search(r"reads_(\d+)\.fastq$", os.path.basename(p))
        if m: dosyalar[m.group(1)].append(p)
    veri = {}
    for tx, yollar in sorted(dosyalar.items()):
        okumalar = []
        for p in sorted(yollar):
            cap = okuma_basina * 4; tut = []
            with open(p, errors="replace") as fh:
                while len(tut) < cap:
                    h = fh.readline()
                    if not h: break
                    s = fh.readline().strip().upper(); fh.readline(); fh.readline()
                    if len(s) >= en_az_uzunluk: tut.append(s)
            okumalar += tut
        random.Random(7).shuffle(okumalar)
        if okumalar: veri[tx] = okumalar[:okuma_basina]
    return veri

def primerleri_oku(yollar):
    """
    Iki dosya bicimini de okur:
      A) 'Oligo adi' / 'Dizi 5-3' sutunlari, isim sonu _F ve _R  (siparis listeleri)
      B) 'grup' / 'F' / 'R' sutunlari                            (grup ve universal listesi)
    """
    ciftler = {}
    for yol in yollar:
        if not os.path.exists(yol): continue
        with open(yol, encoding="utf-8-sig") as fh:
            oku = csv.DictReader(fh)
            basliklar = [b.strip() for b in (oku.fieldnames or [])]
            b_bicimi = ("F" in basliklar and "R" in basliklar)
            for r in oku:
                r = {(k or "").strip(): (v or "").strip() for k, v in r.items()}
                if b_bicimi:
                    ad = guvenli(r.get("grup") or r.get("tip") or "")
                    F = r.get("F", "").upper(); R = r.get("R", "").upper()
                    if not ad or not F or not R: continue
                    if set(F) - set("ACGT") or set(R) - set("ACGT"): continue
                    ciftler[ad] = {"F": F, "R": R}
                else:
                    ad = r.get("Oligo adi") or r.get("ad") or ""
                    dizi = (r.get("Dizi 5-3") or r.get("dizi") or "").upper()
                    if not ad or not dizi: continue
                    taban, rol = ad.rsplit("_", 1)[0], ad.rsplit("_", 1)[-1].upper()
                    if rol not in ("F", "R"): continue
                    if set(dizi) - set("ACGT"): continue
                    ciftler.setdefault(guvenli(taban), {})[rol] = dizi
    return {k: v for k, v in ciftler.items() if "F" in v and "R" in v}

# ------------------------------------------------------------------ kaba kuvvet tarama
def kalip_kur(okumalar_listesi, en_uzun_primer):
    """
    Butun okumalari tek bir bayt dizisine ekler, aralarina hicbir baza esit olmayan
    ayiraclar koyar. Boylece butun numune tek bir vektorel gecisle taranabilir ve
    hicbir pencere iki okumaya birden yayilmaz.
    Doner: (buyuk_dizi, baslangic_indeksleri)
    """
    ayirac = np.full(en_uzun_primer, AYIRAC, dtype=np.uint8)
    parcalar = []; baslangic = []; p = 0
    for s in okumalar_listesi:
        a = kod(s)
        baslangic.append(p)
        parcalar.append(a); parcalar.append(ayirac)
        p += len(a) + en_uzun_primer
    if not parcalar: return np.zeros(0, dtype=np.uint8), np.zeros(0, dtype=np.int64)
    return np.concatenate(parcalar), np.array(baslangic, dtype=np.int64)

def baglanma_yerleri(W, primer, ters, en_fazla_uyumsuz=3, uc_pencere=5, uc_pencere_uyumsuz=1):
    """
    Primerin butun baglanma yerlerini bulur. Tohum yoktur, her pencere denenir.
    ters=False : primer dogrudan kaliba oturur (arti zincir), 3' uc pencerenin SONUNDA
    ters=True  : primerin ters tumleyeni oturur (eksi zincir), 3' uc pencerenin BASINDA
    Doner: pencere baslangic indeksleri
    """
    k = len(primer)
    pk = kod(rc(primer) if ters else primer)
    if W is None or W.shape[0] == 0 or W.shape[1] < k: return np.zeros(0, dtype=np.int64)
    Wk = W[:, :k]
    # ucuz on eleme: 3' uctaki iki baz birebir uymayan pencereler hemen elenir
    if ters: ank = (Wk[:, 0] == pk[0]) & (Wk[:, 1] == pk[1])
    else:    ank = (Wk[:, k-2] == pk[k-2]) & (Wk[:, k-1] == pk[k-1])
    idx = np.nonzero(ank)[0]
    if idx.size == 0: return idx
    e = (Wk[idx] == pk)
    uc = e[:, :uc_pencere] if ters else e[:, k-uc_pencere:]
    uc_hata = (~uc).sum(1); toplam_hata = (~e).sum(1)
    return idx[(uc_hata <= uc_pencere_uyumsuz) & (toplam_hata <= en_fazla_uyumsuz)]

def tek_primer_orani(W, baslangic, n_okuma, primer):
    """
    Primerin TEK BASINA kac okumada baglanma yeri buldugunu doner.
    Hedefte urun cikmadiginda sebebi ayirt etmek icin gereklidir:
      iki primer de dusukse  -> o bolge okumalarda yok, primerin sucu degil
      biri dusuk digeri yuksekse -> o primer bu organizmaya uymuyor, gercek kacirma
    """
    yer = np.concatenate([baglanma_yerleri(W, primer, False),
                          baglanma_yerleri(W, primer, True)])
    if yer.size == 0: return 0.0
    okuma_no = np.searchsorted(baslangic, yer, side="right") - 1
    return len(set(okuma_no.tolist())) / max(1, n_okuma)

def kapsayan_okuma(W, baslangic, F, R):
    """
    Iki primerden EN AZ BIRININ baglandigi ayri okuma sayisi.

    NEDEN GEREKLI
    Capraz orani, rakibin butun okumalarina bolunuyordu. Lokusu kapsayan okuma sayisi
    hicbir yerde hesaplanmadigi icin esik pratikte "rakipte en fazla su kadar okumada
    urun" demeye donusuyordu. Boyle bir olcut, "lokusu kapsayan 5 okumanin 5'inde urun"
    ile "lokusu kapsayan 130 okumanin 5'inde urun" arasindaki farki goremez; birincisi
    rakipte tam cogalma demektir ve TEMIZ yazilmamalidir.

    Bu sayi paydadir: urunun gorunebilecegi okuma sayisi. Bir okumada primerlerden
    hicbiri baglanmiyorsa o okuma lokusa hic ulasmamistir ve paydada yeri yoktur.
    """
    yer = np.concatenate([baglanma_yerleri(W, F, False), baglanma_yerleri(W, F, True),
                          baglanma_yerleri(W, R, False), baglanma_yerleri(W, R, True)])
    if yer.size == 0: return 0
    return len(set((np.searchsorted(baslangic, yer, side="right") - 1).tolist()))

def urun_boyu(arti_5, kp, eksi_3, kq, en_kisa, en_uzun):
    """
    arti_5 : arti zincire baglanan primerin 5' ucu (pencere basi), 3' ucu arti_5+kp-1
    eksi_3 : eksi zincire baglanan primerin 3' ucu (pencere basi), 5' ucu eksi_3+kq-1
    Urun icin uc uclarin birbirine bakmasi ve mesafenin aralikta olmasi gerekir.
    """
    if arti_5 + kp - 1 > eksi_3: return None
    boy = (eksi_3 + kq - 1) - arti_5 + 1
    return boy if en_kisa <= boy <= en_uzun else None

def cift_tara(W, baslangic, okuma_boylari, F, R, en_kisa, en_uzun):
    """
    Bir primer cifti icin, hangi okumalarda urun olustugunu ve urun boylarini doner.
    Doner: (urun_veren_okuma_sayisi, [boylar])
    """
    kf, kr = len(F), len(R)
    yer = {("F", False): baglanma_yerleri(W, F, False),
           ("F", True):  baglanma_yerleri(W, F, True),
           ("R", False): baglanma_yerleri(W, R, False),
           ("R", True):  baglanma_yerleri(W, R, True)}
    if (yer[("F", False)].size + yer[("F", True)].size == 0): return 0, []
    if (yer[("R", False)].size + yer[("R", True)].size == 0): return 0, []

    def okumaya_dagit(pozlar):
        d = defaultdict(list)
        if pozlar.size == 0: return d
        okuma_no = np.searchsorted(baslangic, pozlar, side="right") - 1
        for o, p in zip(okuma_no.tolist(), pozlar.tolist()):
            d[o].append(p - int(baslangic[o]))
        return d

    dF0, dF1 = okumaya_dagit(yer[("F", False)]), okumaya_dagit(yer[("F", True)])
    dR0, dR1 = okumaya_dagit(yer[("R", False)]), okumaya_dagit(yer[("R", True)])
    sayi = 0; boylar = []
    for o in set(list(dF0) + list(dF1)):
        bulundu = None
        # F arti zincirde, R eksi zincirde
        for i in dF0.get(o, []):
            for j in dR1.get(o, []):
                b = urun_boyu(i, kf, j, kr, en_kisa, en_uzun)
                if b: bulundu = b; break
            if bulundu: break
        # R arti zincirde, F eksi zincirde
        if not bulundu:
            for i in dR0.get(o, []):
                for j in dF1.get(o, []):
                    b = urun_boyu(i, kr, j, kf, en_kisa, en_uzun)
                    if b: bulundu = b; break
                if bulundu: break
        if bulundu: sayi += 1; boylar.append(bulundu)
    return sayi, boylar

# ------------------------------------------------------------------ istatistik ve karar
def wilson_alt(k, n, z=1.96):
    """Oranin guvenilen en dusuk degeri. Kucuk orneklemde ham orana guvenilmez."""
    if n <= 0: return 0.0
    p = k / n
    payda = 1 + z*z/n
    orta = p + z*z/(2*n)
    yayilim = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n))
    return max(0.0, (orta - yayilim) / payda)

def capraz_etiketi(rakip_alt, hedef_alt):
    """Tek bir alt sinirdan capraz etiketi. Iki paydayla ayri ayri cagrilir."""
    if rakip_alt <= 0.02: return "TEMIZ"
    if rakip_alt <= 0.10 and hedef_alt >= 10 * rakip_alt: return "SINIRDA"
    return "VAR"

def karar_ver(hedefler, sayimlar, toplamlar, erisilebilir, kapsayanlar=None):
    """
    Doner: dict(capraz, kapsama, hedef_alt, rakip_tx, rakip_alt, eksik, olculemeyen,
                capraz_ham, capraz_kapsamali, rakip_alt_kapsamali, ayrilik)

    IKI PAYDA, IKI KARAR
    capraz_ham        : payda taksonun butun ornek okumalari
    capraz_kapsamali  : payda lokusu kapsayan okumalar (kapsayan_okuma ile olculur)
    Ikisi ayrilirsa KATI olani alinir ve ayrilik dondurulur. Sessizce birine gecmek
    1 numarali kurali ihlal ederdi. Kapsamali payda verilmezse yalnizca ham hesaplanir
    ve bu, capraz_kapsamali'nin bos olmasindan anlasilir.

    KAPSAMA UC DEGERLI
    Okumasi HIC olmayan hedef, olculmemis hedeftir; TAM da EKSIK de degildir. Onceki
    surumde boyle hedefler eksik listesine hic girmiyordu ve liste bosalinca kapsama
    TAM yaziliyordu, yani hic olculmemis bir hedef icin tam kapsama iddia ediliyordu.
    """
    hset = set(hedefler)
    # HEDEF ORANININ PAYDASI DA LOKUSU KAPSAYAN OKUMALAR OLMALI (29 Temmuz).
    #
    # Onceki surumde payda taksonun BUTUN okumalariydi. Rakip tarafinda bu hata
    # zaten duzeltilmisti (kapsayan payda), hedef tarafinda duzeltilmemisti; iki
    # taraf ASIMETRIK olcuyordu.
    #
    # Somut zarari: numunenin okuma kutuphanesi karisik, hem kisa 16S amplikonlari
    # hem tam operon okumalari var. Olculdu: M. mazei'nin 102.006 okumasinin yalnizca
    # 9.257'si (>=2000 bp) operonu kapsiyor. 16S-23S ARA BOLGESINE oturan bir cift,
    # kendi hedefinde mukemmel calissa bile butun okumalarin ancak %9'unda urun
    # verir ve "hedefinde urun yok" diye ELENIR. Oysa gercek qPCR genomik DNA
    # uzerinde kosar, kisa amplikon kutuphanesi uzerinde degil.
    #
    # Methanosarcina turleri 16S'te ayrilmiyor; ayrim tam da bu ara bolgeye bagli.
    # Yani bu payda hatasi, projenin en kritik hedeflerini yapisal olarak eliyordu.
    #
    # Kapsayan payda verilmediginde eski davranis surer (geriye donuk uyumlu).
    def _hedef_oran(t):
        n = sayimlar.get(t, 0)
        ham = toplamlar.get(t, 0)
        if not kapsayanlar: return wilson_alt(n, ham)
        # Kapsayan sayisi urun verenden kucuk cikamaz; cikiyorsa olcum tutarsizdir.
        kap = max(kapsayanlar.get(t, 0), n)
        return wilson_alt(n, kap)
    hedef_alt = min((_hedef_oran(t) for t in hedefler), default=0.0)
    hedef_alt_ham = min((wilson_alt(sayimlar.get(t,0), toplamlar.get(t,0)) for t in hedefler),
                        default=0.0)

    rakip_alt = 0.0; rakip_tx = ""
    rakip_alt_k = 0.0; rakip_tx_k = ""
    for tx, n in sayimlar.items():
        if tx in hset: continue
        a = wilson_alt(n, toplamlar.get(tx, 0))
        if a > rakip_alt: rakip_alt, rakip_tx = a, tx
        if kapsayanlar:
            # payda kapsayan okuma sayisi. Kapsayan sayisi urun veren sayisindan
            # kucuk cikamaz; cikiyorsa olcum tutarsizdir, guvenli tarafta kalinir.
            kn = max(kapsayanlar.get(tx, 0), n)
            ak = wilson_alt(n, kn)
            if ak > rakip_alt_k: rakip_alt_k, rakip_tx_k = ak, tx

    capraz_ham = capraz_etiketi(rakip_alt, hedef_alt)
    capraz_k = capraz_etiketi(rakip_alt_k, hedef_alt) if kapsayanlar else ""
    sira = {"TEMIZ": 0, "SINIRDA": 1, "VAR": 2}
    if capraz_k and capraz_k != capraz_ham:
        capraz = capraz_ham if sira[capraz_ham] > sira[capraz_k] else capraz_k
        ayrilik = (f"ham payda {capraz_ham} (rakip alt {rakip_alt:.4f}), "
                   f"kapsayan payda {capraz_k} (rakip alt {rakip_alt_k:.4f}), "
                   f"kati olan alindi: {capraz}")
    else:
        capraz = capraz_ham; ayrilik = ""

    # Okumasi hic olmayan hedef OLCULMEMISTIR, eksik de degildir tam da degildir.
    olculemeyen = [t for t in hedefler if toplamlar.get(t, 0) == 0]
    eksik = []
    for t in hedefler:
        n = toplamlar.get(t, 0)
        if n == 0: continue
        oran = sayimlar.get(t, 0) / n
        ulasilir = erisilebilir.get(t, 0.0)
        # o takson baska bir ciftle yuksek oran veriyorsa okumalari saglamdir;
        # burada dusuk oran cikiyorsa sucu okumaya atamayiz
        if ulasilir >= 0.30 and oran < 0.25 * ulasilir: eksik.append(t)

    if olculemeyen: kapsama = "OLCULEMEDI"
    elif eksik: kapsama = "EKSIK"
    elif len(hedefler) <= 1: kapsama = "-"
    else: kapsama = "TAM"

    # KARARIN SAYISI ILE KARARIN ADI AYNI HESAPTAN GELMELI (29 Temmuz duzeltmesi).
    #
    # Onceki surumde rakip_tx KAPSAMALI galipten, rakip_alt ise HAM galipten
    # geliyordu. Ikisi farkli taksonlar oldugunda cikan cumle sacmaydi:
    # "en iyi rakip Y, alt sinir 0.071" -- isim Y'nin, sayi X'in.
    #
    # Daha agiri: 100_TASARIM eleme kapisinda `ra = okuma_karari["rakip_alt"]`
    # kullaniyordu, yani HAM sayiyi. Somut olcum: rakibin 1000 okumasindan yalnizca
    # 6'si lokusu kapsiyor ve 6'sinin ALTISINDA da urun var. Ham oran 6/1000 ->
    # alt sinir 0.0028; kapsamali oran 6/6 -> alt sinir 0.6097. Karar dogru sekilde
    # "VAR" cikiyor ama eleme kapisi 0.0028'i 0.04 esigiyle karsilastirdigi icin
    # aday ELENMIYOR ve ekrana "esik 0.02" yazarken esigin ALTINDA bir sayi
    # basiliyordu. Rakipte TAM cogalma, temiz gibi geciyordu.
    #
    # Artik karara hangi payda yol actiysa, adi ve sayisi O paydadan aliniyor.
    if capraz_k and capraz_k == capraz and rakip_tx_k:
        karar_tx, karar_alt, karar_payda = rakip_tx_k, rakip_alt_k, "kapsayan"
    else:
        karar_tx, karar_alt, karar_payda = rakip_tx, rakip_alt, "ham"

    return dict(capraz=capraz, kapsama=kapsama, hedef_alt=hedef_alt,
                hedef_alt_ham=hedef_alt_ham,
                rakip_tx=karar_tx,
                # rakip_alt ARTIK KARARLA TUTARLI olandir. Ham deger ayri alanda.
                rakip_alt=karar_alt,
                rakip_alt_ham=rakip_alt, rakip_alt_kapsamali=rakip_alt_k,
                rakip_tx_ham=rakip_tx, rakip_tx_kapsamali=rakip_tx_k,
                karar_paydasi=karar_payda,
                eksik=eksik, olculemeyen=olculemeyen,
                capraz_ham=capraz_ham, capraz_kapsamali=capraz_k, ayrilik=ayrilik)

# ------------------------------------------------------------------ birim testleri
def birim_testleri(yazdir=True):
    """
    Karar veren her fonksiyonun bilinen cevapli sinavi. Disaridan hicbir arac gerekmez.
    Testlerin cogu "capraz VAR mi dogru buluyor" tarafini zorlar, cunku yanlis cevabin
    pahali yonu odur: kotu bir primeri temiz gostermek.
    """
    hata = [0]
    def K(ad, bulunan, beklenen):
        ok = bulunan == beklenen
        if not ok: hata[0] += 1
        if yazdir: print(f"  {'GECTI' if ok else 'KALDI'}  {ad:<50} {bulunan} / {beklenen}")

    def tara(kalip, primer, ters=False):
        b, bas = kalip_kur([kalip], 30)
        W = np.lib.stride_tricks.sliding_window_view(b, 30)
        return [int(x) for x in baglanma_yerleri(W, primer, ters)]
    def urun(kalip, F, R, en_kisa=60, en_uzun=400):
        b, bas = kalip_kur([kalip], 30)
        W = np.lib.stride_tricks.sliding_window_view(b, 30)
        n, boy = cift_tara(W, bas, [len(kalip)], F, R, en_kisa, en_uzun)
        return (boy[0] if boy else None)

    rnd = random.Random(11)
    kalip = "".join(rnd.choice("ACGT") for _ in range(1200))
    F = kalip[300:320]; R = rc(kalip[430:450])
    D = {"A":"C","C":"A","G":"T","T":"G"}
    def yaz(d, p, y): return d[:p] + y + d[p+1:]
    def boz(d, *pozlar):
        for p in pozlar: d = yaz(d, p, D[d[p]])
        return d

    if yazdir: print(u'  --- the 3\' end rule (the last two bases exact, at most one mismatch in the last five)')
    K("son bazda uyumsuzluk baglanmaz",        tara(boz(kalip,319), F), [])
    K("sondan ikincide uyumsuzluk baglanmaz",  tara(boz(kalip,318), F), [])
    K("SONDAN UCUNCUDE uyumsuzluk BAGLANIR",   tara(boz(kalip,317), F), [300])
    K("sondan dorduncude uyumsuzluk baglanir", tara(boz(kalip,316), F), [300])
    K("son bes bazda iki uyumsuzluk baglanmaz",tara(boz(kalip,317,316), F), [])
    if yazdir: print("  --- uyumsuzluk butcesi")
    K("tam eslesme baglanir",                  tara(kalip, F), [300])
    K("uc ic uyumsuzluk baglanir",             tara(boz(kalip,302,305,308), F), [300])
    K("dort ic uyumsuzluk baglanmaz",          tara(boz(kalip,302,305,308,311), F), [])
    K("5' ucta iki uyumsuzluk baglanir",       tara(boz(kalip,300,301), F), [300])
    if yazdir: print(u'  --- strand orientation')
    K("ileri primer arti zincirde bulunur",    tara(kalip, F, False), [300])
    K("ileri primer eksi zincirde bulunmaz",   tara(kalip, F, True), [])
    K("geri primer eksi zincirde bulunur",     tara(kalip, R, True), [430])
    K("geri primer arti zincirde bulunmaz",    tara(kalip, R, False), [])
    if yazdir: print(u'  --- product geometry')
    K("dogru yon ve mesafede urun",            urun(kalip, F, R), 150)
    K("ters cevrilmis kalipta da urun",        urun(rc(kalip), F, R), 150)
    K("primer sirasi sonucu degistirmez",      urun(kalip, R, F), 150)
    K("alakasiz kalipta urun yok",
      urun("".join(rnd.choice("ACGT") for _ in range(1200)), F, R), None)
    K("3' uc bozulunca urun yok",              urun(boz(kalip,319), F, R), None)
    uzak = kalip[:300] + F + "".join(rnd.choice("ACGT") for _ in range(900)) + rc(R) + kalip[900:]
    K("cok uzak primerler urun vermez",        urun(uzak, F, R), None)
    ters_yon = (kalip[:300] + rc(R) + kalip[320:400] + F + kalip[420:])
    K("birbirine bakmayan primerler urun vermez", urun(ters_yon, F, R), None)
    if yazdir: print("  --- kalip kurulumu (okumalar birbirine karismamali)")
    b, bas = kalip_kur([kalip[:400], kalip[400:800]], 30)
    W = np.lib.stride_tricks.sliding_window_view(b, 30)
    n, _ = cift_tara(W, bas, [400, 400], F, R, 60, 400)
    K("iki ayri okumaya bolunmus bolge urun vermez", n, 0)
    n2, _ = cift_tara(W, bas, [400, 400], kalip[300:320], rc(kalip[360:380]), 40, 200)
    K("tek okuma icindeki bolge urun verir",   n2, 1)
    if yazdir: print("  --- istatistik")
    K("wilson kucuk orneklemi cezalandirir",   round(wilson_alt(1, 17), 3), 0.01)
    K("wilson buyuk orneklemde orana yaklasir",round(wilson_alt(150, 300), 2), 0.44)
    K("sifir sayim sifir alt sinir",           wilson_alt(0, 300), 0.0)

    if yazdir: print("  --- kapsayan okuma sayisi (capraz oraninin paydasi)")
    kl = kalip[:400]
    b3, bas3 = kalip_kur([kl, kl, "".join(rnd.choice("ACGT") for _ in range(400))], 30)
    W3 = np.lib.stride_tricks.sliding_window_view(b3, 30)
    K("primerin bagladigi okumalar sayilir",
      kapsayan_okuma(W3, bas3, kl[300:320], rc(kl[360:380])), 2)
    K("hicbir primerin baglanmadigi durumda sifir",
      kapsayan_okuma(W3, bas3, "ACGTACGTACGTACGTACGT", "TTTTTTTTTTTTTTTTTTTT"), 0)

    if yazdir: print(u'  --- the decision (this function had never been tested before)')
    # 1. Okumasi HIC olmayan hedef. Onceki surum bunu eksik listesine sokmuyordu,
    #    liste bosalinca kapsama TAM yaziliyordu: hic olculmemis hedef icin tam
    #    kapsama iddiasi.
    k1 = karar_ver(["83986", "118126", "394967", "2201"], {"2223": 1}, {"2223": 300}, {})
    K("okumasi olmayan hedefte kapsama TAM yazilmaz", k1["kapsama"], "OLCULEMEDI")
    K("olculemeyen hedefler listelenir", len(k1["olculemeyen"]), 4)

    # 2. Rakip lokusu cok az kapsiyor ama kapsadigi her okumada urun veriyor.
    #    Ham payda ile TEMIZ, kapsayan payda ile VAR cikar. Kati olan alinmali.
    hed = {"818": 150}; top = {"818": 300, "28116": 300}
    say = dict(hed); say["28116"] = 5
    k2 = karar_ver(["818"], say, top, {}, {"28116": 6, "818": 300})
    K("ham payda tek basina TEMIZ derdi", k2["capraz_ham"], "TEMIZ")
    K("kapsayan payda ile capraz gorunur", k2["capraz_kapsamali"], "VAR")
    K("kati olan karar alinir", k2["capraz"], "VAR")
    K("ayrilik loga yazilmak uzere dondurulur", bool(k2["ayrilik"]), True)
    # 29 TEMMUZ: KARARIN SAYISI KARARIN PAYDASINDAN GELMELI.
    # Onceki surumde rakip_alt HAM paydadan geliyordu; 100_TASARIM eleme kapisi o
    # sayiyi 0.04 esigiyle karsilastirdigi icin, rakipte TAM cogalma gosteren aday
    # ELENMIYORDU. Bu madde tam o kapiyi kilitler.
    K("karar kapsayan paydadan geldi", k2["karar_paydasi"], "kapsayan")
    K("rakip_alt kararla tutarli (ham degil)", k2["rakip_alt"] == k2["rakip_alt_kapsamali"], True)
    K("rakip_alt eleme esigini asiyor", k2["rakip_alt"] >= 0.04, True)
    K("ham deger ayri alanda korunuyor", k2["rakip_alt_ham"] < 0.04, True)
    # Isim ve sayi ayni taksondan gelmeli
    K("rakip adi ile sayisi ayni kaynaktan", k2["rakip_tx"], k2["rakip_tx_kapsamali"])

    # 4. HEDEF PAYDASI DA KAPSAYAN OKUMA OLMALI (29 Temmuz).
    # Ara bolgeye oturan bir cift: hedefin 1000 okumasindan yalnizca 90'i lokusu
    # kapsiyor ve 90'inin 88'inde urun var. Ham payda 88/1000 -> 0.07 (elenir),
    # kapsayan payda 88/90 -> 0.92 (mukemmel). Kisa 16S amplikonlari paydayi
    # sisirdigi icin M. mazei/barkeri gibi ara bolge hedefleri yapisal olarak
    # eleniyordu.
    k4 = karar_ver(["818"], {"818": 88}, {"818": 1000}, {}, {"818": 90})
    K("hedef orani kapsayan paydayla hesaplanir", k4["hedef_alt"] > 0.80, True)
    K("ham payda ayri alanda korunuyor", k4["hedef_alt_ham"] < 0.12, True)
    # Kapsayan payda verilmezse eski davranis surer
    k5 = karar_ver(["818"], {"818": 88}, {"818": 1000}, {})
    K("kapsayan payda yoksa ham davranis surer", abs(k5["hedef_alt"] - k5["hedef_alt_ham"]) < 1e-9, True)
    # Tutarsiz olcum: kapsayan sayisi urun verenden kucukse guvenli tarafa gecilir
    k6 = karar_ver(["818"], {"818": 88}, {"818": 1000}, {}, {"818": 10})
    K("kapsayan sayisi urunden kucukse alt sinir 1'i asmaz", k6["hedef_alt"] <= 1.0, True)

    # 3. Ayni sayilar, ama rakip lokusu genis kapsiyor. Ikisi de TEMIZ, ayrilik yok.
    k3 = karar_ver(["818"], say, top, {}, {"28116": 250, "818": 300})
    K("genis kapsamada iki payda da TEMIZ", (k3["capraz_ham"], k3["capraz_kapsamali"]),
      ("TEMIZ", "TEMIZ"))
    K("ayrilik yoksa bos", k3["ayrilik"], "")

    # 4. Kapsayan sayisi urun veren sayisindan kucuk gelirse olcum tutarsizdir;
    #    payda urun veren sayisina cekilir, oran 1'i asmaz.
    k4 = karar_ver(["818"], say, top, {}, {"28116": 2, "818": 300})
    K("tutarsiz kapsama oraninda alt sinir 1'i asmaz", k4["rakip_alt_kapsamali"] <= 1.0, True)
    return hata[0]

# ------------------------------------------------------------------ blastn ikinci gorus
def blast_ikinci_gorus(veri, ciftler, cikti, threads, en_kisa, en_uzun):
    """
    Ayni soruyu blastn ile de sorar. blastn hizalamayi puanina gore KIRPTIGI icin
    3' uca yakin uyumsuzluklari olan baglanmalari bildirmez; bu yuzden bu sayim
    genellikle kaba kuvvet taramasindan DUSUK cikar ve karar icin kullanilmaz.
    Amac yalnizca iki bagimsiz yontemin ayni yone isaret ettigini gormektir.
    """
    if not (arac_var("blastn") and arac_var("makeblastdb")):
        log("blastn bulunamadi, ikinci gorus atlandi"); return {}
    fa = os.path.join(cikti, "numune_okumalari.fasta")
    with open(fa, "w") as fh:
        for tx, okumalar in veri.items():
            for i, s in enumerate(okumalar): fh.write(f">{tx}|{i}\n{s}\n")
    db = os.path.join(cikti, "numune_db")
    r = subprocess.run(["makeblastdb","-in",fa,"-dbtype","nucl","-out",db],
                       capture_output=True, text=True)
    if r.returncode != 0: log("makeblastdb basarisiz, ikinci gorus atlandi"); return {}
    pf = os.path.join(cikti, "primerler.fasta")
    with open(pf, "w") as fh:
        for taban, d in sorted(ciftler.items()):
            fh.write(f">{taban}_F\n{d['F']}\n>{taban}_R\n{d['R']}\n")
    r = subprocess.run(["blastn","-task","blastn-short","-query",pf,"-db",db,
        "-outfmt","6 qseqid sseqid qstart qend sstart send sstrand qlen",
        "-word_size","7","-evalue","1000","-dust","no","-soft_masking","false",
        "-max_target_seqs","100000","-num_threads",str(threads)],
        capture_output=True, text=True)
    if r.returncode != 0: log("blastn basarisiz, ikinci gorus atlandi"); return {}
    hit = defaultdict(lambda: defaultdict(lambda: {"F": [], "R": []}))
    for satir in r.stdout.splitlines():
        p = satir.split("\t")
        if len(p) < 8: continue
        q, s, qs, qe, ss, se, strand, qlen = p[:8]
        try: qs, qe, ss, se, qlen = int(qs), int(qe), int(ss), int(se), int(qlen)
        except ValueError: continue
        if "_" not in q: continue
        taban, rol = q.rsplit("_", 1); rol = rol.upper()
        if rol not in ("F","R"): continue
        if qe != qlen: continue                      # BLAST kirpmasi, 3' uc yok
        hit[taban][s][rol].append((ss, se, strand, qlen))
    sonuc = defaultdict(lambda: defaultdict(int))
    for taban, okumalar in hit.items():
        for okuma_id, d in okumalar.items():
            tx = okuma_id.split("|")[0]
            var = False
            for fs, fe, fst, fl in d["F"]:
                for rs, re_, rst, rl in d["R"]:
                    if fst == rst: continue
                    arti = (fs, fe, fl) if fst == "plus" else (rs, re_, rl)
                    eksi = (rs, re_, rl) if fst == "plus" else (fs, fe, fl)
                    if arti[1] > eksi[1]: continue
                    boy = eksi[0] - arti[0] + 1
                    if en_kisa <= boy <= en_uzun: var = True; break
                if var: break
            if var: sonuc[taban][tx] += 1
    return sonuc

# ------------------------------------------------------------------ selftest
def selftest():
    print("=" * 72); print("KURALLARIN BILINEN CEVAPLI SINAVI"); print("=" * 72)
    h = birim_testleri()
    print("=" * 72)
    print("SELFTEST GECTI" if h == 0 else f"SELFTEST KALDI, {h} test kaldi")
    print("=" * 72)
    return 0 if h == 0 else 1

# ------------------------------------------------------------------ ana
def main():
    global LOGF
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", "--kok", dest="kok"); ap.add_argument("--primerler", nargs="*", default=[])
    ap.add_argument("--out", default=""); ap.add_argument("--okuma", type=int, default=300)
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--en-kisa", type=int, default=60)
    ap.add_argument("--en-uzun", type=int, default=400)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--blast", action="store_true", help="second opinion via blastn (does not decide)")
    a = ap.parse_args()
    if a.selftest: sys.exit(selftest())
    if not a.kok: ap.error("--root is required")

    zaman = time.strftime("%Y%m%d_%H%M%S")
    cikti = a.out or os.path.join(a.kok, "tools", "0_TESLIM_RAPOR", "ICPCR_" + zaman)
    os.makedirs(cikti, exist_ok=True)
    LOGF = open(os.path.join(cikti, "log.txt"), "w", encoding="utf-8")
    log(f"cikti: {cikti}")

    log("kurallarin bilinen cevapli sinavi calistiriliyor")
    if birim_testleri(yazdir=False) != 0:
        log("SINAV BASARISIZ, durduruldu. Ayrinti icin: --selftest"); sys.exit(2)
    log("sinav gecti")

    primer_yollari = [p if os.path.isabs(p) else os.path.join(a.kok, p) for p in a.primerler]
    if not primer_yollari:
        primer_yollari = [os.path.join(a.kok, "tools", "0_TESLIM_RAPOR", f) for f in
                          ("OLIGO_SIPARIS_GUNCEL_10ASSAY.csv", "OLIGO_SIPARIS_OPSIYONEL.csv",
                           "GRUP_VE_UNIVERSAL_PRIMERLER.csv")]
    ciftler = primerleri_oku(primer_yollari)
    log(f"{len(ciftler)} primer cifti okundu")
    if not ciftler: sys.exit("primer bulunamadi")

    log("okumalar yukleniyor")
    veri = okumalari_yukle(a.kok, a.okuma)
    toplam_okuma = sum(len(v) for v in veri.values())
    toplam_baz = sum(len(s) for v in veri.values() for s in v)
    log(f"{len(veri)} takson, {toplam_okuma} okuma, {toplam_baz/1e6:.1f} milyon baz")

    en_uzun_primer = max(max(len(d["F"]), len(d["R"])) for d in ciftler.values())
    log(f"kaba kuvvet taramasi hazirlaniyor (pencere {en_uzun_primer} baz)")
    kalipler = {}
    for tx, okumalar in veri.items():
        b, bas = kalip_kur(okumalar, en_uzun_primer)
        W = np.lib.stride_tricks.sliding_window_view(b, en_uzun_primer)
        kalipler[tx] = (W, bas, [len(s) for s in okumalar])

    log("tarama basliyor")
    t0 = time.time()
    sayimlar = {}; boy_bilgisi = {}; kapsayanlar = {}
    toplamlar = {tx: len(v) for tx, v in veri.items()}
    # Hedef listelerinde gecip numunede okumasi olmayan taksonlar. Bunlar icin hicbir
    # olcum yapilamaz; "olculmedi" ile "temiz cikti" ayni sonuca varmamali.
    hedefte_gecen = {t for hs in HEDEF_TAXID.values() for t in hs}
    okumasiz = sorted(hedefte_gecen - set(toplamlar))
    if okumasiz:
        log("")
        log(f"UYARI: {len(okumasiz)} hedef taksonun numunede HIC okumasi yok, "
            f"bunlar icin olcum YAPILAMAZ:")
        for t in okumasiz: log(f"    {isim(t)} (taxid {t})")
    for i, taban in enumerate(sorted(ciftler), 1):
        F, R = ciftler[taban]["F"], ciftler[taban]["R"]
        s = {}; b = {}
        for tx, (W, bas, boylar) in kalipler.items():
            n, urunler = cift_tara(W, bas, boylar, F, R, a.en_kisa, a.en_uzun)
            if n: s[tx] = n; b[tx] = int(statistics.median(urunler))
        # Capraz oraninin paydasi icin lokusu kapsayan okuma sayisi. Yalnizca urun
        # veren taksonlar ve hedefler icin hesaplanir; digerlerinde paya zaten sifir
        # oldugu icin payda karari degistirmez.
        ilgili = set(s) | set(HEDEF_TAXID.get(taban, []))
        kap = {}
        for tx in ilgili:
            if tx in kalipler:
                W, bas, _ = kalipler[tx]
                kap[tx] = kapsayan_okuma(W, bas, F, R)
        kapsayanlar[taban] = kap
        sayimlar[taban] = s; boy_bilgisi[taban] = b
        log(f"  [{i}/{len(ciftler)}] {taban}  {sum(s.values())} urun veren okuma")
    log(f"tarama bitti ({int(time.time()-t0)} sn)")

    # her taksonun ulasilabilirligi: herhangi bir ciftle elde edilen en yuksek oran
    erisilebilir = {}
    for tx, n in toplamlar.items():
        erisilebilir[tx] = max((sayimlar[c].get(tx, 0) / max(1, n) for c in sayimlar), default=0.0)

    blast = {}
    if a.blast:
        log("blastn ile ikinci gorus aliniyor")
        blast = blast_ikinci_gorus(veri, ciftler, cikti, a.threads, a.en_kisa, a.en_uzun)

    satirlar = []; ozet = []
    log("")
    log("SONUCLAR")
    for taban in sorted(ciftler):
        s = sayimlar[taban]; b = boy_bilgisi[taban]
        hedefler = HEDEF_TAXID.get(taban, [])
        K = karar_ver(hedefler, s, toplamlar, erisilebilir, kapsayanlar.get(taban))
        capraz, kapsama, hedef_alt = K["capraz"], K["kapsama"], K["hedef_alt"]
        rakip_tx, rakip_alt, eksik = K["rakip_tx"], K["rakip_alt"], K["eksik"]
        etiket = f"capraz {capraz}" + (f", kapsama {kapsama}" if kapsama != "-" else "")
        log("")
        log(f"  {taban}   [{etiket}]")
        if K["ayrilik"]:
            log(f"    AYRILIK: {K['ayrilik']}")
        if K["olculemeyen"]:
            log("    OLCULEMEYEN HEDEFLER (numunede okumasi yok, kapsama iddia edilemez): "
                + ", ".join(isim(t) for t in K["olculemeyen"]))
        sebep = {}
        for t in eksik:
            W, bas, boylar = kalipler[t]
            of = tek_primer_orani(W, bas, len(boylar), ciftler[taban]["F"])
            orr = tek_primer_orani(W, bas, len(boylar), ciftler[taban]["R"])
            if max(of, orr) < 0.25:
                sebep[t] = (f"bu bolge okumalarda yok (ileri {of:.2f}, geri {orr:.2f}), "
                            f"primerin sucu degil")
            else:
                zayif = "geri" if orr < of else "ileri"
                sebep[t] = (f"{zayif} primer bu organizmaya UYMUYOR "
                            f"(ileri {of:.2f}, geri {orr:.2f})")
        for t in hedefler:
            n = toplamlar.get(t, 0); k = s.get(t, 0)
            im = f"  <<< {sebep[t]}" if t in eksik else ""
            log(f"    HEDEF  {isim(t):<34} {k:>4}/{n:<4} oran {k/max(1,n):.3f}"
                f"  urun {b.get(t,0)} bp{im}")
        rakipler = sorted(((tx, k) for tx, k in s.items() if tx not in set(hedefler)),
                          key=lambda x: -x[1] / max(1, toplamlar.get(x[0], 1)))
        if not rakipler:
            log("    rakip taksonlarin hicbirinde urun olusmadi")
        kap = kapsayanlar.get(taban, {})
        for tx, k in rakipler[:6]:
            n = toplamlar.get(tx, 0); kn = max(kap.get(tx, 0), k)
            # Iki oran yan yana yazilir. Aralarindaki fark, rakibin lokusu ne kadar
            # kapsadigini gosterir; buyuk fark "az okumada urun" ile "kapsayan
            # okumalarin hepsinde urun" ayrimidir.
            log(f"    rakip  {isim(tx):<34} {k:>4}/{n:<4} oran {k/max(1,n):.3f}"
                f"  |  kapsayan {kn:>4} icinde oran {k/max(1,kn):.3f}"
                f"  urun {b.get(tx,0)} bp")
        for tx in sorted(set(list(s) + hedefler)):
            n = toplamlar.get(tx, 0); k = s.get(tx, 0); kn = max(kap.get(tx, 0), k)
            satirlar.append(dict(cift=taban, taxid=tx, organizma=isim(tx),
                hedef_mi=("hedef" if tx in set(hedefler) else "rakip"),
                urun_veren=k, toplam_okuma=n, oran=round(k/max(1,n), 4),
                oran_alt_sinir=round(wilson_alt(k, n), 4),
                kapsayan_okuma=kn,
                kapsamali_oran=(round(k/kn, 4) if kn else ""),
                kapsamali_alt_sinir=(round(wilson_alt(k, kn), 4) if kn else ""),
                medyan_urun_bp=b.get(tx, 0),
                blastn_ikinci_gorus=(blast.get(taban, {}).get(tx, "") if blast else "")))
        ozet.append(dict(cift=taban, capraz=capraz, kapsama=kapsama,
            hedef_alt_sinir=round(hedef_alt, 4),
            en_iyi_rakip=isim(rakip_tx) if rakip_tx else "",
            en_iyi_rakip_alt_sinir=round(rakip_alt, 4),
            capraz_ham_payda=K["capraz_ham"],
            capraz_kapsayan_payda=K["capraz_kapsamali"],
            rakip_alt_kapsamali=round(K["rakip_alt_kapsamali"], 4),
            olcum_ayriligi=K["ayrilik"],
            olculemeyen_hedefler="; ".join(isim(t) for t in K["olculemeyen"]),
            urun_vermeyen_hedefler="; ".join(isim(t) for t in eksik),
            sebep="; ".join(f"{isim(t)}: {sebep[t]}" for t in eksik)))

    yol = os.path.join(cikti, "icpcr_sonuc.csv")
    with open(yol, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=["cift","taxid","organizma","hedef_mi","urun_veren",
            "toplam_okuma","oran","oran_alt_sinir","kapsayan_okuma","kapsamali_oran",
            "kapsamali_alt_sinir","medyan_urun_bp","blastn_ikinci_gorus"])
        w.writeheader(); w.writerows(satirlar)
    yol2 = os.path.join(cikti, "icpcr_ozet.csv")
    with open(yol2, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=["cift","capraz","kapsama","hedef_alt_sinir",
            "en_iyi_rakip","en_iyi_rakip_alt_sinir","capraz_ham_payda",
            "capraz_kapsayan_payda","rakip_alt_kapsamali","olcum_ayriligi",
            "olculemeyen_hedefler","urun_vermeyen_hedefler","sebep"])
        w.writeheader(); w.writerows(ozet)

    sira = {"TEMIZ": 0, "SINIRDA": 1, "VAR": 2}
    log("")
    log("TOPLU KARAR")
    log(f"  {'capraz':<9}{'kapsama':<9}{'cift':<46}{'hedef':>7}{'rakip':>8}")
    for satir in sorted(ozet, key=lambda x: (sira.get(x["capraz"], 3), x["kapsama"] == "EKSIK")):
        log(f"  {satir['capraz']:<9}{satir['kapsama']:<9}{satir['cift']:<46}"
            f"{satir['hedef_alt_sinir']:>7.3f}{satir['en_iyi_rakip_alt_sinir']:>8.3f}"
            f"  {satir['en_iyi_rakip']}")
    eksikli = [o for o in ozet if o["urun_vermeyen_hedefler"]]
    if eksikli:
        log("")
        log("HEDEFLEDIGI HALDE URUN VERMEDIGI ORGANIZMALAR")
        for o in eksikli:
            log(f"  {o['cift']}")
            log(f"     {o['sebep']}")
    log("")
    log("Oranlar ham sayimdir; karar Wilson alt sinirina gore verilir, boylece on yedi")
    log("okumali bir taksonda tek okuma yuzde alti gibi gorunup karari bozmaz.")
    log("")
    log(f"yazildi: {yol}")
    log(f"yazildi: {yol2}")
    LOGF.close()

if __name__ == "__main__":
    main()
