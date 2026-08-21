#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09_ozgulluk.py
08'in ürettiği aday çiftleri toplantı kararındaki özgüllük ve doğrulama
kurallarından geçirir.

Uygulanan kurallar:
  1. HAM OKUMADA DOĞRULAMA. Konsensüste geçen çift, hedefin ham okumalarında
     tekrar sınanır. Ürün çıkmazsa aday "numuneden_dogrulanamadi" işaretini
     alır ve elenir.
  2. RAKİP ORANI, WILSON ALT SINIRI. Rakip taksonların okumalarında ürün
     veren okuma oranı ham sayı olarak değil Wilson alt sınırıyla
     değerlendirilir; on yedi okumalı bir taksonda tek okuma yüzde altı
     görünüp kararı bozmasın diye.
  3. SUŞ İÇİ DEĞİŞKENLİK. İki primer ayrı ayrı ölçülür. Biri okumaların
     yüzde sekseninde bağlanırken diğeri yüzde kırkında kalıyorsa ikincisi
     değişken bir pozisyonda oturuyordur ve aday cezalandırılır.
  4. DIŞ VERİTABANI ÖZGÜLLÜĞÜ. mfeprimer ile referans veritabanlarında
     amplikon aranır (--bind-amp-only). İkinci bağımsız ölçüm olarak blastn
     ile bağlanma yerleri sayılır; ikisi ayrışırsa aday elenir ve ayrılık
     log'a yazılır.

Her satır tarih ve saatle log'lanır. Her hedeften sonra checkpoint yazılır,
koşu yarıda kesilirse kaldığı yerden devam eder.

Kullanım:
  python3 09_ozgulluk.py \
      --adaylar "/.../primer_adaylari" \
      --pt      "/.../PrimerTasarlama" \
      --out     "/.../primer_final" \
      [--top 30] [--max-okuma 200000] [--atla-mfe] [--atla-blast]
"""
import hashlib
import argparse, csv, datetime, glob, importlib.util, json, math, os, re
import subprocess, sys, time, collections

HERE = os.path.dirname(os.path.abspath(__file__))
TS = lambda: datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
_LOG = [None]


def log(msg):
    line = "[%s] %s" % (TS(), msg)
    print(line, flush=True)
    if _LOG[0]:
        with open(_LOG[0], "a", encoding="utf-8") as fh:
            fh.write(line + "\n")


def load(p, n):
    sp = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(sp)
    bak, sys.argv = sys.argv, [p]
    try:
        sp.loader.exec_module(m)
    finally:
        sys.argv = bak
    return m


E = load(os.path.join(HERE, "03_primer_aday_uret.py"), "e03")
G = load(os.path.join(HERE, "04_grup_primer.py"), "g04")
rc = E.rc


class Kural:
    """04 ile birebir ayni baglanma kurali."""
    exact_last = 2
    tail_len = 5
    tail_max_mm = 1
    total_max_mm = 3
    min_overlap = 15


KURAL = Kural()


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hizalama

MAPPY = hizalama.ARKA_UC is not None


def _girdi_parmak_izi(a):
    """Kosunun bagli oldugu her girdinin ozeti. Degisirse checkpoint duser.
    Aday dosyalari degistiginde eski dogrulama sonuclarinin sessizce
    kullanilmasi, kosunun bitmis gorunmesine ama sonucun eski girdiye ait
    olmasina yol aciyordu."""
    h = hashlib.sha256()
    h.update(("adaylar=%s\n" % os.path.abspath(a.adaylar)).encode())
    for f in sorted(glob.glob(os.path.join(a.adaylar, "*__*.tsv"))):
        try:
            st = os.stat(f)
            h.update(("%s|%d|%d\n" % (os.path.basename(f), st.st_size,
                                      int(st.st_mtime))).encode())
        except OSError:
            pass
    for ek in ("ayirt_edilemez.tsv",):
        yol = os.path.join(a.adaylar, ek)
        try:
            st = os.stat(yol)
            h.update(("%s|%d|%d\n" % (ek, st.st_size, int(st.st_mtime))).encode())
        except OSError:
            h.update(("%s|yok\n" % ek).encode())
    try:
        st = os.stat(a.hedefler)
        h.update(("hedefler|%d|%d\n" % (st.st_size, int(st.st_mtime))).encode())
    except OSError:
        pass
    # Motor betikleri de parmak ize girer: baglanma kurali ya da tarama
    # mantigi degistiginde eski dogrulama sonuclari gecerli degildir.
    _burada = os.path.dirname(os.path.abspath(__file__))
    for _b in ("04_grup_primer.py", "03_primer_aday_uret.py",
               "09_ozgulluk.py", "hizalama.py"):
        try:
            _st = os.stat(os.path.join(_burada, _b))
            h.update(("%s|%d|%d\n" % (_b, _st.st_size, int(_st.st_mtime))).encode())
        except OSError:
            h.update(("%s|yok\n" % _b).encode())
    h.update(("kons=%s|top=%d|maxokuma=%d|wilson=%.5f|susici=%.5f|"
              "minuye=%s|bulasma=%d|mfe=%d|blast=%d\n"
              % (os.path.abspath(a.kons) if a.kons else "-", a.top,
                 a.max_okuma, a.rakip_wilson_max, a.sus_ici_fark_max,
                 getattr(a, "min_uye_orani", ""), a.bulasma_ornek,
                 int(a.atla_mfe), int(a.atla_blast))).encode())
    h.update(("bulasma_min=%s|sizinti_tavan=%s\n"
              % (getattr(a, "bulasma_min_okuma", ""),
                 getattr(a, "sizinti_tavan", ""))).encode())
    return h.hexdigest()[:16]


def wilson_ust(k, n, z=1.96):
    """Oranin Wilson UST siniri."""
    if n == 0:
        return 1.0
    p = k / n
    d = 1 + z * z / n
    m = p + z * z / (2 * n)
    s_ = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return min(1.0, (m + s_) / d)


def wilson_alt(k, n, z=1.96):
    """Bir orani ham sayiyla degil Wilson alt siniriyla degerlendirir.
    k basari, n deneme. n=0 ise 0 doner."""
    if n <= 0:
        return 0.0
    p = k / n
    d = 1 + z * z / n
    merkez = p + z * z / (2 * n)
    yari = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return max(0.0, (merkez - yari) / d)


def fastq_oku(path, limit):
    """Okuma dizilerini uretir, limit kadar."""
    n = 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                yield line.strip().upper()
                n += 1
                if limit and n >= limit:
                    return


class _PArg:
    """product_len icin gereken en kucuk arguman kumesi."""
    prod_min = 50
    prod_hard_max = 400


_OKUMA_BELLEK = {}


def okumalari_al(path, limit):
    """Ayni fastq bir hedefteki butun adaylar icin tekrar tekrar taraniyor;
    diskten bir kez okunup bellekte tutulur."""
    k = (path, limit)
    v = _OKUMA_BELLEK.get(k)
    if v is None:
        v = [(s, rc(s)) for s in fastq_oku(path, limit)]
        # Onbellek en fazla --onbellek-dosya kadar fastq tutar. Eski surumde
        # her isabetsizlikte tamamen bosaltiliyordu, yani isabet orani
        # pratikte sifirdi ve ayni fastq her aday icin yeniden okunuyordu.
        while len(_OKUMA_BELLEK) >= 3:
            _OKUMA_BELLEK.pop(next(iter(_OKUMA_BELLEK)))          # tek hedefin dosyalari; bellek sismesin
        _OKUMA_BELLEK[k] = v
    return v


_BULASMA = {}


def capraz_bulasma(hedef_kons_yolu, rakip_kons_yolu, rakip_fq, ornek=400,
                   min_uzunluk=400):
    """Rakip kutusundaki okumalarin ne kadari aslinda HEDEFE ait.

    Kraken kutulari birbirine sizabiliyor: bir kutudaki okumalarin bir kismi
    baska taksonun molekulleridir. Bu okumalar dogal olarak hedef primerle
    urun verir. O yuzden 'rakipte urun var' karari, olculen bu sizinti
    oraniyla karsilastirilmadan verilemez. Doner: (k, n) = hedefe daha iyi
    uyan okuma sayisi ve incelenen okuma sayisi."""
    if not MAPPY:
        return (0, 0)
    anahtar = (hedef_kons_yolu, rakip_kons_yolu, rakip_fq)
    if anahtar in _BULASMA:
        return _BULASMA[anahtar]
    def _oku(p):
        """Bastaki ve sondaki N kirpilir, IC N yerinde birakilir. Ic N'leri
        silmek kapsama bosluklarinin iki yanini birlestirip kimerik referans
        uretir; hizalama sonuclari o kavsakta anlamsizlasir."""
        d = "".join(l.strip() for l in open(p, encoding="utf-8",
                                            errors="replace")
                    if not l.startswith(">")).upper()
        return d.strip("N")
    try:
        A_h = hizalama.Hizalayici(seq=_oku(hedef_kons_yolu), preset="map-ont")
        A_r = hizalama.Hizalayici(seq=_oku(rakip_kons_yolu), preset="map-ont")
    except Exception:
        _BULASMA[anahtar] = (0, 0)
        return (0, 0)
    if not A_h or not A_r:
        _BULASMA[anahtar] = (0, 0)
        return (0, 0)
    okumalar = {}
    n = 0
    try:
        with open(rakip_fq, errors="replace") as fh:
            for i, line in enumerate(fh):
                if i % 4 != 1:
                    continue
                r = line.strip()
                if len(r) < min_uzunluk:
                    continue
                n += 1
                if n > ornek:
                    break
                okumalar["o%d" % n] = r
    except OSError:
        pass
    # Toplu hizalama: komut satiri arka ucunda okuma basina surec baslatmak
    # kabul edilemez derecede yavas olurdu.
    skor_h = {ad: max((x.mlen for x in hl), default=0)
              for ad, hl in A_h.map_toplu(okumalar)}
    skor_r = {ad: max((x.mlen for x in hl), default=0)
              for ad, hl in A_r.map_toplu(okumalar)}
    k = sum(1 for ad in okumalar if skor_h.get(ad, 0) > skor_r.get(ad, 0))
    n = len(okumalar)
    _BULASMA[anahtar] = (k, n)
    return (k, n)


def okuma_taramasi(path, F, R, prod_min, prod_max, limit):
    """Ham okumalarda cift dogrulamasi.

    Iki asamali. Once 3' ucun son 12 bazi str.find ile hizlica aranir, bu
    C hizinda ve okumalarin cogunu eler. Sonra yalnizca aday okumalarda
    04'un sinanmis find_bindings ve product_len fonksiyonlari calistirilir,
    yani baglanma kurali ile urun koordinati tasarim asamasindakiyle birebir
    ayni kodla degerlendirilir.

    Doner: (toplam_okuma, F_baglanan, R_baglanan, urun_veren)
    """
    pa = _PArg()
    pa.prod_min, pa.prod_hard_max = prod_min, prod_max
    K = KURAL.tail_len
    tot = f_hit = r_hit = both = 0
    # Eski surumde kapi F[-12:] tam eslesmesini sarti kosuyordu. Bu kural
    # metninden KATIYDI: kural son bes bazda bir, toplamda uc uyumsuzluga
    # izin veriyor, dolayisiyla -12..-3 araligindaki tek bir uyumsuzluk bile
    # gecerli bir baglanmayi tarama disi birakiyordu. Olculdu: ONT benzeri
    # %5 hatada gercek baglanmalarin %15'i, %8 hatada %30'u kayboluyordu ve
    # rakip urun orani sistematik olarak DUSUK cikiyordu, yani ozgul olmayan
    # ciftler GECTI aliyordu. Kapi artik find_noidx'in kullandigi cekirdek
    # kumesinin AYNISIYLA kuruluyor; boylece kapi hicbir zaman kuralin
    # gecirdigi bir okumayi elemez.

    # Okuma basina 5-mer indeksi kurmak darbogazdi: 3700 bazlik bir okumada
    # her cagride ~3700 girdilik sozluk olusuyordu. Okuma tek bir dizi
    # oldugu icin indeks yerine cekirdek varyantlarini str.find ile taramak
    # cok daha ucuz; sonuc find_bindings ile birebir ayni olmali ve bu
    # asagida sinanmistir.
    varyant = {}

    def _var(oligo):
        v = varyant.get(oligo)
        if v is None:
            v = sorted(G.seed_variants(oligo[-K:], KURAL.tail_max_mm,
                                       KURAL.exact_last))
            varyant[oligo] = v
        return v

    def find_noidx(oligo, seq):
        L, n = len(seq), len(oligo)
        saf = all(c in "ACGT" for c in oligo)
        hits, gor = [], set()
        for v in _var(oligo):
            st = 0
            while True:
                pos = seq.find(v, st)
                if pos < 0:
                    break
                st = pos + 1
                start = pos - (n - K)
                end = start + n
                if end > L:
                    continue
                j0 = max(0, -start)
                if n - j0 < KURAL.min_overlap:
                    continue
                # Hizli yol: ham okumalar yalnizca A, C, G, T ve N icerir.
                # Oligo da saf ACGT ise base_match'in kume kesisimi yerine
                # dogrudan karakter karsilastirmasi yeterlidir ve ayni sonucu
                # verir (N her zaman uyumsuz sayilir). Kume islemi okuma
                # basina yuzlerce kez cagrildigi icin bu fark buyuk.
                mm = 0
                ok = True
                if saf:
                    for j in range(j0, n):
                        if oligo[j] != seq[start + j]:
                            mm += 1
                            if mm > KURAL.total_max_mm:
                                ok = False
                                break
                else:
                    for j in range(j0, n):
                        if not G.base_match(oligo[j], seq[start + j]):
                            mm += 1
                            if mm > KURAL.total_max_mm:
                                ok = False
                                break
                if ok:
                    key = (end - 1, mm)
                    if key not in gor:
                        gor.add(key)
                        hits.append(key)
        return hits

    def baglanmalar(oligo, seq, seqrc):
        return dict(L=len(seq), plus=find_noidx(oligo, seq),
                    minus=find_noidx(oligo, seqrc))

    kapiF = _var(F)
    kapiR = _var(R)

    def _kapi_gecer(seq, seqrc):
        """find_noidx ile AYNI cekirdek kumesi. Kapi yalnizca hicbir cekirdek
        varyanti bulunmayan okumalari eler; bu okumalarda find_noidx da
        tanimi geregi bos doner, dolayisiyla kapi sonucu degistirmez."""
        for v in kapiF:
            if v in seq or v in seqrc:
                return True
        for v in kapiR:
            if v in seq or v in seqrc:
                return True
        return False

    for seq, seqrc in okumalari_al(path, limit):
        tot += 1
        if not _kapi_gecer(seq, seqrc):
            continue
        bf = baglanmalar(F, seq, seqrc)
        br = baglanmalar(R, seq, seqrc)
        fb = bool(bf["plus"] or bf["minus"])
        rb = bool(br["plus"] or br["minus"])
        if fb:
            f_hit += 1
        if rb:
            r_hit += 1
        if fb and rb:
            p = G.product_len(bf, br, len(F), len(R), pa, pmax=prod_max)
            if p is not None:
                both += 1
    return tot, f_hit, r_hit, both


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--adaylar", required=True, help="08'in cikti klasoru")
    p.add_argument("--pt", required=True, help="PrimerTasarlama kok klasoru")
    p.add_argument("--out", required=True)
    p.add_argument("--hedefler", default=os.path.join(HERE, "hedefler.tsv"))
    p.add_argument("--top", type=int, default=15,
                   help="hedef basina sinanacak en iyi aday sayisi")
    p.add_argument("--max-okuma", type=int, default=20000,
                   help="takson basina taranacak en fazla okuma; 0 sinirsiz. "
                        "Kirpma yapilirsa log'a yazilir.")
    p.add_argument("--min-uye-orani", type=float, default=0.5,
                   help="hedef uyelerin en az bu orani ham okumada dogrulanmali")
    p.add_argument("--kons", default=None,
                   help="konsensus klasoru. Verilirse rakip kutulardaki "
                        "capraz bulasma olculur ve 'rakipte urun' karari bu "
                        "olculen sizintiyla karsilastirilarak verilir.")
    p.add_argument("--prod-min", type=int, default=70,
                   help="tasarim asamasindakiyle ayni olmali")
    p.add_argument("--prod-hard-max", type=int, default=300,
                   help="tasarim asamasindakiyle ayni olmali")
    p.add_argument("--bulasma-ornek", type=int, default=400)
    p.add_argument("--bulasma-min-okuma", type=int, default=100,
                   help="capraz bulasma olcumu bu kadar okumadan az ise "
                        "sizintinin esigi acmasina izin verilmez")
    p.add_argument("--sizinti-tavan", type=float, default=0.15,
                   help="sizinti esigi en fazla bu degere kadar acabilir")
    p.add_argument("--rakip-wilson-max", type=float, default=0.02,
                   help="rakipte urun veren okuma oraninin Wilson alt siniri "
                        "bu degeri asarsa aday elenir")
    p.add_argument("--sus-ici-fark-max", type=float, default=0.40,
                   help="iki primerin baglanma orani arasindaki en buyuk fark; "
                        "asilirsa aday cezalandirilir")
    p.add_argument("--atla-mfe", action="store_true")
    p.add_argument("--atla-blast", action="store_true")
    p.add_argument("--mfe", default=None, help="mfeprimer ikilisinin yolu")
    p.add_argument("--yeniden", action="store_true")
    return p.parse_args()


def main():
    a = get_args()
    os.makedirs(a.out, exist_ok=True)
    _LOG[0] = os.path.join(a.out, "ozgulluk.log")
    CKPT = os.path.join(a.out, "checkpoint.json")
    parmak = _girdi_parmak_izi(a)
    ckpt = {}
    if os.path.exists(CKPT) and not a.yeniden:
        try:
            ham = json.load(open(CKPT, encoding="utf-8"))
            eski = ham.get("_girdi_parmak_izi") if isinstance(ham, dict) else None
            kayitlar = {k: v for k, v in ham.items()
                        if not k.startswith("_")} if isinstance(ham, dict) else {}
            if eski is None:
                log("checkpoint parmak izi TASIMIYOR (eski surumden kalma), "
                    "yok sayiliyor")
            elif eski != parmak:
                log("checkpoint parmak izi UYUSMUYOR, yok sayiliyor")
                log("   kayitli : %s" % eski)
                log("   simdiki : %s" % parmak)
                log("   Aday dosyalari ya da esikler degismis; eski "
                    "dogrulama sonuclari yeniden kullanilmayacak.")
            else:
                ckpt = kayitlar
                log("checkpoint gecerli: %d hedef zaten islenmis, atlanacak"
                    % len(ckpt))
        except Exception as e:
            log("checkpoint okunamadi (%s)" % e)
    t0 = time.time()
    log("baslangic. adaylar=%s" % a.adaylar)

    mfe = a.mfe or os.path.join(a.pt, "ARACLAR", "mfeprimer")
    kullan_mfe = (not a.atla_mfe) and os.path.exists(mfe) and os.access(mfe, os.X_OK)
    if not a.atla_mfe and not kullan_mfe:
        log("UYARI: mfeprimer bulunamadi ya da calistirilabilir degil (%s), "
            "dis veritabani adimi atlanacak" % mfe)
    # Dis veritabani ozgullugu ARTIK 14_dis_veritabani.py'de yapiliyor.
    # Buradaki bayraklar yalnizca geriye donuk uyum icin duruyor; hangi
    # adimin nerede kostugunu log'a acikca yaziyoruz ki belge ile kod
    # arasinda sessiz bir fark kalmasin.
    kullan_blast = (not a.atla_blast) and bool(
        subprocess.run(["bash", "-c", "command -v blastn"],
                       capture_output=True).stdout.strip())
    if not a.atla_blast and not kullan_blast:
        log("UYARI: blastn bulunamadi, ikinci olcum atlanacak")

    # fastq envanteri: sinif_taxid -> yol
    fq = {}
    for p in glob.glob(os.path.join(a.pt, "fastq files", "*", "*.fastq")):
        grp = os.path.basename(os.path.dirname(p))
        m = re.search(r"reads[-_](\d+)", os.path.basename(p))
        if m:
            fq[(re.split(r"[-_]", grp)[0], grp, m.group(1))] = p
    log("fastq envanteri: %d dosya" % len(fq))

    # konsensus envanteri: (grup, taxid) -> yol. Capraz bulasma olcumu icin.
    kons = {}
    if a.kons:
        for p2 in glob.glob(os.path.join(a.kons, "*_konsensus.fasta")):
            m = re.match(r"((?:A1|A2|B|F1|F2)-\d+)_(\d+)_", os.path.basename(p2))
            if m:
                kons[(m.group(1), m.group(2))] = p2
        log("konsensus envanteri: %d dosya" % len(kons))
        log(hizalama.durum())
        if not MAPPY:
            log("UYARI: hizalama arka ucu yok, capraz bulasma olculemiyor. "
                "Sabit esik (--rakip-wilson-max %.3f) kullanilacak."
                % a.rakip_wilson_max)
    else:
        log("--kons verilmedi: capraz bulasma olculmeyecek, sabit esik "
            "kullanilacak (--rakip-wilson-max %.3f)" % a.rakip_wilson_max)

    def _kons_of(fq_yolu):
        b = os.path.basename(fq_yolu)
        g = os.path.basename(os.path.dirname(fq_yolu))
        m = re.search(r"reads[-_](\d+)", b)
        return kons.get((g, m.group(1))) if m else None

    # hedef tanimlari
    hedefler = {}
    for line in open(a.hedefler, encoding="utf-8"):
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        if p[0] == "karar":
            continue
        hedefler[p[1]] = dict(karar=p[0], duzey=p[2], inn=p[3],
                              haric=p[4] if len(p) > 4 else "")

    # 08'in olctugu ayirt edilemez takson ciftleri. Ayni ayiklama burada da
    # uygulanmazsa 08'in rakip listesinden cikardigi takson 09'da geri gelir
    # ve hedefin kendi dizisi rakip sayilir.
    # 08'in bozuk konsensus yuzunden disladigi taksonlar. 09 fastq
    # envanterinden calistigi icin ayni taksonu rakip (ve alan hedeflerinde
    # uye) olarak geri aliyordu; 08 ile 09 ayni kumeler uzerinde calismali.
    # (grup, taxid) ikilisi. Takson bazinda dislamak yanlis olurdu: bir
    # taksonun konsensusu bir numunede bos, otekinde saglam olabilir.
    # Dislama YALNIZ uye kumesine uygulanir; rakip kumesinde bu kutularin
    # ham okumalari degerlidir ve onlari cikarmak ozgulluk denetimini
    # zayiflatir. Konsensusun bos olmasi okumalarin yok oldugu anlamina
    # gelmez, yalnizca konsensus kurulamadigi anlamina gelir.
    dislanan = set()
    dl = os.path.join(a.adaylar, "dislanan_takson.tsv")
    if os.path.exists(dl):
        # Sutunlar KONUMA gore degil BASLIGA gore okunur. Eski surumlerde bu
        # dosya 4 sutunluydu (taxid, etiket, uzunluk, kapsanan); konuma gore
        # okuyan kod o dosyada (taxid, etiket) ikilisini (grup, taxid) sanip
        # hicbir seyle eslesmeyen bir kume kurar ve dislama SESSIZCE devre
        # disi kalir. Basligi tanimazsak durmak, yanlis kumeyle devam
        # etmekten iyidir.
        with open(dl, encoding="utf-8") as fh:
            satirlar = [l.rstrip("\n") for l in fh
                        if l.strip() and not l.startswith("#")]
        if satirlar:
            basliklar = satirlar[0].split("\t")
            if "grup" not in basliklar or "taxid" not in basliklar:
                sys.exit("HATA: %s basligi taninmadi: %s\n"
                         "Beklenen sutunlar: grup, taxid, etiket, uzunluk, "
                         "kapsanan.\nBu dosya eski bir 08 surumunden kalmis "
                         "olabilir; 08'i yeniden calistirin." % (dl, basliklar))
            ig, it = basliklar.index("grup"), basliklar.index("taxid")
            for satir in satirlar[1:]:
                p2 = satir.split("\t")
                if len(p2) > max(ig, it) and p2[ig].strip() and p2[it].strip():
                    dislanan.add((p2[ig].strip(), p2[it].strip()))
        log("08'in disladigi (grup, takson): %d -> uye kumesinden cikarilir, "
            "rakip kumesinde ham okumalariyla KALIR" % len(dislanan))
        for g, t in sorted(dislanan):
            log("   %s %s" % (g, t))

    ayirt = {}
    ae = os.path.join(a.adaylar, "ayirt_edilemez.tsv")
    if os.path.exists(ae):
        for r in csv.DictReader(open(ae, encoding="utf-8"), delimiter="\t"):
            ayirt.setdefault((r["sinif"], r["taxid1"]), set()).add(r["taxid2"])
            ayirt.setdefault((r["sinif"], r["taxid2"]), set()).add(r["taxid1"])
        log("ayirt edilemez cift tablosu okundu: %s (%d kayit)"
            % (ae, sum(len(v) for v in ayirt.values()) // 2))
    else:
        log("UYARI: %s yok, ayirt edilemez takson ayiklamasi yapilmayacak" % ae)

    sonuc = []
    tsvler = sorted(glob.glob(os.path.join(a.adaylar, "*__*.tsv")))
    log("NOT: dis veritabani ozgullugu (mfeprimer ve blastn) bu betikte "
        "DEGIL, 14_dis_veritabani.py adiminda calisir.")
    log("islenecek aday dosyasi: %d" % len(tsvler))

    for ti, tsv in enumerate(tsvler, 1):
        etiket = os.path.basename(tsv)[:-4]
        hedef, sinif = etiket.rsplit("__", 1)
        if etiket in ckpt and not a.yeniden:
            log("[%d/%d] ATLANDI (checkpoint) %s" % (ti, len(tsvler), etiket))
            sonuc.extend(ckpt[etiket])
            continue
        th = time.time()
        rows = list(csv.DictReader(open(tsv, encoding="utf-8"), delimiter="\t"))
        if not rows:
            log("[%d/%d] %s bos, atlandi" % (ti, len(tsvler), etiket))
            ckpt[etiket] = []
            continue
        h = hedefler.get(hedef, {})
        in_t = set(x for x in h.get("inn", "").split(",") if x and x != "*")
        haric = set(x for x in h.get("haric", "").split(",") if x)
        alan = h.get("inn", "").startswith("*")
        # Uye fastq'lari TAKSON basina gruplanir. Eski surumde dogrulama
        # sayaci DOSYA basina artiyordu; ayni taksonun birden cok yildaki
        # dosyasi varsa, hic dogrulanmayan bir uye takson oranin icinde
        # kaybolabiliyordu. Olculdu: Asetoklastik_metanojenler A2'de bes
        # uye taksondan yalnizca ikisi urun verdiginde dosya bazli oran
        # 7/10=0,70 ile GECTI, takson bazli oran 2/5=0,40 ile elenmeliydi.
        if alan:
            # Alan hedefi (*A, *B, *F): uye = ayni harfle baslayan butun
            # siniflar, rakip = OTEKI alanlarin siniflari. 08 tam olarak
            # boyle kuruyor; 09'un eski surumu rakip listesini bos
            # birakiyordu ve evrensel primerler hicbir rakip sinanmadan
            # GECTI aliyordu. Ayrica 'haric' sutunu alan hedeflerinde de
            # uye kumesine uygulanir.
            # 08 ile BIREBIR ayni kurulum: uye = YALNIZ bu amplikon
            # sinifindaki taksonlar, rakip = baska ALANLARIN siniflari.
            harf = sinif[0]
            uye_ikili = [(t, v) for (s, g, t), v in fq.items()
                         if s == sinif and t not in haric
                         and (g, t) not in dislanan]
            rakip_fq = [v for (s, g, t), v in fq.items() if s[0] != harf]
        else:
            uye_ikili = [(t, v) for (s, g, t), v in fq.items()
                         if s == sinif and t in in_t and (g, t) not in dislanan]
            rakip_fq = [v for (s, g, t), v in fq.items()
                        if s == sinif and t not in in_t and t not in haric]
        uye_fq = [v for _, v in uye_ikili]
        uye_takson = {}
        for t, v in uye_ikili:
            uye_takson.setdefault(t, []).append(v)
        atilan_tx = set()
        if ayirt and in_t:
            temiz = []
            for v in rakip_fq:
                tx = None
                for (s2, g2, t2), v2 in fq.items():
                    if v2 == v:
                        tx = t2
                        break
                if tx and (ayirt.get((sinif, tx), set()) & in_t):
                    atilan_tx.add(tx)
                else:
                    temiz.append(v)
            rakip_fq = temiz
        if atilan_tx:
            log("      rakipten cikarildi (hedefle ayirt edilemez): %s"
                % ", ".join(sorted(atilan_tx)))
        log("[%d/%d] %-46s aday=%d uye_takson=%d uye_fastq=%d rakip_fastq=%d "
            "(fastq basina en fazla %d okuma)"
            % (ti, len(tsvler), etiket, len(rows), len(uye_takson),
               len(uye_fq), len(rakip_fq), a.max_okuma))
        if not uye_fq:
            log("      uye fastq bulunamadi, ham okuma dogrulamasi atlanacak")

        uye_kons = [k for k in (_kons_of(x) for x in uye_fq) if k]

        def _bulasma_getir(rakip_yolu):
            """Rakip kutusundaki okumalarin kacinin aslinda HEDEFE ait
            oldugunu olcer. En yuksek sizinti veren hedef uyesi alinir."""
            rk = _kons_of(rakip_yolu)
            if not (uye_kons and rk and MAPPY):
                return (0, 0)
            en = (0, 0)
            for hk in uye_kons:
                k, n = capraz_bulasma(hk, rk, rakip_yolu, a.bulasma_ornek)
                if n and (not en[1] or k / n > en[0] / max(1, en[1])):
                    en = (k, n)
            return en

        rows.sort(key=lambda r: float(r.get("ceza", 9e9)))
        gecen = []
        for ri, r in enumerate(rows[:a.top], 1):
            F, R = r["ileri_dizi"], r["geri_dizi"]
            # Urun penceresi tasarim asamasindakiyle AYNI olmali; farkli
            # olursa sinirdaki urunler iki asamada farkli degerlendirilir.
            pmin, pmax = a.prod_min, a.prod_hard_max
            uye_dogru = 0
            f_or, r_or = [], []
            dogrulanmayan = []
            for tx, yollar in sorted(uye_takson.items()):
                tx_urun = 0
                for p in yollar:
                    tot, fh_, rh, both = okuma_taramasi(p, F, R, pmin, pmax,
                                                        a.max_okuma)
                    if tot:
                        f_or.append(fh_ / tot)
                        r_or.append(rh / tot)
                        tx_urun += both
                if tx_urun > 0:
                    uye_dogru += 1
                else:
                    dogrulanmayan.append(tx)
            uye_orani = (uye_dogru / len(uye_takson)) if uye_takson else None
            # rakipte urun veren okuma orani, Wilson alt siniri
            rak_w = 0.0
            rak_detay = []
            rakipte_gercek = False
            for p in rakip_fq:
                tot, fh_, rh, both = okuma_taramasi(p, F, R, pmin, pmax, a.max_okuma)
                w = wilson_alt(both, tot)
                # Olculen kutu sizintisinin UST siniri. Urun orani bunun
                # altindaysa, urun veren okumalar yanlis kutuya dusmus HEDEF
                # okumalariyla aciklanabilir; rakibin kendisinin cogaldigina
                # dair kanit degildir.
                bk, bn = _bulasma_getir(p)
                # Wilson UST siniri az okumada 1'e yaklasir. Korumasiz
                # birakilirsa 10 okumalik bir rakip kutusunda esik 0,28'e,
                # 5 okumalikta 0,43'e cikar ve okumalarinin yarisinda urun
                # veren bir rakip "temiz" sayilir. Bu yuzden sizinti ancak
                # yeterli okuma varsa esigi acabilir ve actigi miktar
                # tavanlanir.
                if bn >= a.bulasma_min_okuma:
                    sizinti_ust = min(wilson_ust(bk, bn), a.sizinti_tavan)
                    sizinti_not = "%.4f" % sizinti_ust
                else:
                    sizinti_ust = 0.0
                    sizinti_not = "yetersiz_okuma(%d)" % bn
                esik = max(a.rakip_wilson_max, sizinti_ust)
                if w > esik:
                    rakipte_gercek = True
                rak_detay.append("%s=%d/%d(W%.4f,sizinti<%s)"
                                 % (os.path.basename(p)[:18], both, tot, w,
                                    sizinti_not))
                rak_w = max(rak_w, w)
            fmin = min(f_or) if f_or else 0.0
            rmin = min(r_or) if r_or else 0.0
            sus_fark = abs((sum(f_or) / len(f_or) if f_or else 0)
                           - (sum(r_or) / len(r_or) if r_or else 0))
            durum = []
            # Hicbir ham okuma taranmamis bir aday "gecti" sayilamaz. Eski
            # surumde uye_fq ya da rakip_fq bos oldugunda hicbir etiket
            # eklenmiyor ve aday GECTI olarak teslim ediliyordu.
            if not uye_fq:
                durum.append("uye_okumasi_yok")
            elif uye_orani is None or uye_orani < a.min_uye_orani:
                durum.append("numuneden_dogrulanamadi")
            if not rakip_fq:
                durum.append("rakip_sinanmadi")
            elif rakipte_gercek:
                durum.append("rakipte_urun")
            if sus_fark > a.sus_ici_fark_max:
                durum.append("sus_ici_degisken")
            r2 = dict(r)
            r2.update(hedef=hedef, sinif=sinif, karar=h.get("karar", ""),
                      uye_dogrulanan=uye_dogru, uye_toplam=len(uye_takson),
                      dogrulanmayan_uye=",".join(dogrulanmayan),
                      uye_orani=round(uye_orani, 3) if uye_orani is not None else "",
                      rakip_wilson=round(rak_w, 5),
                      rakip_detay=";".join(rak_detay)[:200],
                      ileri_baglanma_min=round(fmin, 3),
                      geri_baglanma_min=round(rmin, 3),
                      sus_ici_fark=round(sus_fark, 3),
                      ozgulluk_durum=",".join(durum) if durum else "GECTI")
            gecen.append(r2)
            if len(durum) == 0:
                log("      [%d] GECTI  F=%s R=%s  uye=%d/%d takson rakipW=%.4f"
                    % (ri, F, R, uye_dogru, len(uye_takson), rak_w))
            if sum(1 for g in gecen if g["ozgulluk_durum"] == "GECTI") >= 5:
                log("      bes gecerli aday bulundu, bu hedef icin duruldu")
                break
        ckpt[etiket] = gecen
        sonuc.extend(gecen)
        ckpt["_girdi_parmak_izi"] = parmak
        with open(CKPT, "w", encoding="utf-8") as cf:
            json.dump(ckpt, cf, ensure_ascii=False)
        yaz(a.out, sonuc)
        log("      bitti, %.1f sn, gecen=%d/%d"
            % (time.time() - th,
               sum(1 for g in gecen if g["ozgulluk_durum"] == "GECTI"), len(gecen)))

    yaz(a.out, sonuc)
    ok = sum(1 for s in sonuc if s["ozgulluk_durum"] == "GECTI")
    log("TOPLAM: %d aday sinandi, %d tanesi butun kurallardan gecti" % (len(sonuc), ok))
    log("toplam sure: %.1f dakika" % ((time.time() - t0) / 60))
    log("cikti: %s" % os.path.join(a.out, "primer_final.tsv"))


def yaz(out, sonuc):
    """Sonuc bos olsa bile dosya yeniden yazilir. Eski surumde erken
    donuldugu icin onceki kosunun bayat primer_final.tsv'si diskte kaliyor
    ve 13 onu gecerli sanip Excel'e basiyordu."""
    if not sonuc:
        with open(os.path.join(out, "primer_final.tsv"), "w",
                  encoding="utf-8") as fh:
            fh.write("karar\thedef\tsinif\tozgulluk_durum\n")
        return
    cols = ["karar", "hedef", "sinif", "ozgulluk_durum",
            "ileri_dizi", "ileri_tm", "ileri_gc", "ileri_baslangic", "ileri_uzunluk",
            "geri_dizi", "geri_tm", "geri_gc", "geri_baslangic", "geri_uzunluk",
            "tm_farki", "urun_min", "urun_maks",
            "ileri_iupac_cozulen", "geri_iupac_cozulen", "lokus_varyant_sayisi",
            "uye_dogrulanan", "uye_toplam", "uye_orani", "rakip_wilson",
            "ileri_baglanma_min", "geri_baglanma_min", "sus_ici_fark",
            "yetim_primer", "heterodimer_dg", "ceza", "rakip_detay"]
    cols = [c for c in cols if any(c in s for s in sonuc)]
    with open(os.path.join(out, "primer_final.tsv"), "w", newline="",
              encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(sorted(sonuc, key=lambda s: (s["ozgulluk_durum"] != "GECTI",
                                                 s.get("karar", ""), s["hedef"],
                                                 float(s.get("ceza", 9e9)))))


if __name__ == "__main__":
    main()
