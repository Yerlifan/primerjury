# -*- coding: utf-8 -*-
r"""
uyelik_yeniden_turet.py
=======================
UYELIGI KRAKEN ETIKETINDEN DEGIL, OLCULEN KIMLIKTEN yeniden turetir ve
paneldeki butun ciftleri duzeltilmis uyelikle yeniden olcer.

NEDEN GEREKLI
-------------
Bir hedefin "ayrim kati" sayisi, hangi kutunun UYE hangisinin RAKIP sayildigina
dogrudan baglidir. Kutu etiketleri Kraken'den geliyor ve olculdugu kadariyla
en az 12 kutuda YANLIS. Yanlis etiketli bir kutu rakip hanesine yazilinca metrik
hedefi hedefle kiyaslar ve mukemmel bir primer bile 1'in altinda ayrim verir.
Olculen ornek: Petriella 0,71x -> 8,47x (ayni primer, yalniz uyelik duzeltildi).

TEMEL ILKE
----------
KANIT YOKLUGU KANIT SAYILMAZ. Bir kutu ancak POZITIF olcum kaniti varsa yer
degistirir. Olcum bir kutu icin sinyal uretemediyse o kutunun ESKI durumu korunur.
Bu kural onemlidir: tersi yapilirsa, teshis edilen hatanin aynisi ters yonde
tekrarlanir (bu betigin gelistirilmesi sirasinda bir kez tam olarak bu oldu).

NASIL CALISIR
-------------
1) Sinif ici konsensusler tam uzunlukta hizalanir; kimligi >=%99,0 olanlar ON GRUP.
2) Her on grup icin AYIRT EDICI k-mer kumesi cikarilir:
      grubun k-mer'leri EKSI diger butun gruplarin k-mer'leri
   Bu adim sarttir. Yapilmazsa korunmus bolgeler (18S, 5.8S) uzerinden Trichoderma
   okumalari Petriella'ya atanir - olculdu, %70 yanlis atama.
3) Her okuma, ayirt edici k-mer'lerin normalize payina gore gruba atanir.
   Esik f=0,30; bagimsiz in-silico PCR ile kalibre edildi (uyum: %0-72 araliginda
   birkac puan).
4) Uyelik yeniden turetilir (yukaridaki ilke ile).
5) Butun panel ciftleri HEM mm<=1 HEM mm<=3 ile, TAM DERINLIKTE olculur.
6) Eski ve yeni ayrim katlari yan yana yazilir.

KESINTIYE DAYANIKLI
-------------------
Her asama bitince diske yazar (_ck_*.json). Kesilirse ayni komutla kaldigi
yerden devam eder. Hicbir asama bastan hesaplanmaz.

PANELE YAZMAZ
-------------
Bu betik panel xlsx/tsv dosyalarina YAZMAZ. Yalniz okur ve kendi ciktisini uretir.

KULLANIM
--------
    python3 uyelik_yeniden_turet.py [--kok PROJE_KLASORU] [--nmax 3000] [--sifirla]
"""

# -------------------------------------------------------------------------
# uyelik_yeniden_turet.py — hangi kutunun UYE hangisinin RAKIP oldugunu Kraken
# etiketinden degil OLCULEN kimlikten yeniden turetir ve butun panel ciftlerini
# duzeltilmis uyelikle yeniden olcer.
#
# GİRDİ  : konsensus_kanonik/*.kanonik.fa (kutu konsensusleri),
#          "fastq files"/*/*.fastq(.gz) (ham okumalar),
#          screening/hedef_uyelik.tsv (mevcut uyelik tanimi),
#          primer_final/devir_ciftleri_20260802_sonrotus_TESLIM.tsv (panel ciftleri).
# ÇIKTI  : engine_SONUC/engine_TURETME.md,
#          engine_SONUC/ciftler_yeniden_olcum.tsv,
#          engine_SONUC/kutu_olculen_kimlik.tsv,
#          engine_SONUC/_ck_*.json (kesinti kontrol noktalari).
#          Panel dosyalarina YAZMAZ.
# ÇAĞRAN : screening.bat -> U tusu
#          (bat icinde: wsl -e python3 "engine/uyelik_yeniden_turet.py" --kok .)
#
# BU DOSYANIN URETTIGI TABLO P, K ve D ASAMALARININ GIRDISIDIR: uyelik degisirse
# panelin butun ayrim katlari degisir. Bu yuzden zincirde en basta durur.
# -------------------------------------------------------------------------
import os, sys, json, glob, random, argparse, time, csv, re

K = 21
NORM_ESIK = 0.30      # ayirt edici k-mer normalize esigi (kalibre edildi)
UYE_ESIK = 50.0       # okumalarinin >=%50'si hedef gruptaysa UYE
KARISIK_ESIK = 15.0   # %15-50 arasi KARISIK (ne uye ne rakip)
KONS_ESIK = 99.0      # on gruplama icin konsensus kimlik esigi
OKUMA_MIN, OKUMA_MAX = 200, 6000
TOHUM = 20260802
ENKOTU_ASGARI = 150
KAPSAM_ESIGI = 0.20
SINIF_ORNEK = 300     # kimlik olcumunde kutu basina okuma

try:
    import numpy as np
except ImportError:
    sys.exit('HATA: numpy yok.  WSL icinde:  pip3 install numpy --break-system-packages')

_C = str.maketrans('ACGTURYSWKMBDHVNacgturyswkmbdhvn', 'TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn')
def rc(s): return s.translate(_C)[::-1]
def temizle(s): return re.sub(r'[^ACGTUN]', 'N', s.upper()).replace('U', 'T')
_M = np.full(256, 4, dtype=np.int8)
for _i, _c in enumerate('ACGT'): _M[ord(_c)] = _i
def enc(s): return _M[np.frombuffer(s.encode(), dtype=np.uint8)]

# Basit FASTA okuyucu. Konsensus dosyalari kucuktur, akis yeterlidir.
def fasta(p):
    h = None; b = []
    for line in open(p, 'rt', errors='ignore'):
        if line.startswith('>'):
            if h: yield h, ''.join(b)
            h = line[1:].strip(); b = []
        else: b.append(line.strip())
    if h: yield h, ''.join(b)

# ----------------------------------------------------------------- hizalama
# ---------------------------------------------------------------------------
# INFIX (HW) LEVENSHTEIN - son satiri dondurur.
# Neden infix: kisa sorgu uzun hedefin ICINE hizalanir, yani hedefin basinda ve
# sonunda kalan fazlalik CEZALANDIRILMAZ. Konsensuslerin uzunluklari cok farkli
# (1,5 kb ile 4,5 kb yan yana); global hizalama bu farki uyumsuzluk gibi sayar ve
# ayni organizmayi farkli gosterirdi.
#
# Sol komsu bagimliligi (ekleme) vektorlestirildi: simdi[j] = min(aday[j],
# simdi[j-1]+1) bagintisi a[j] = simdi[j]-j konularak kosan minimuma cevrilir ve
# np.minimum.accumulate ile tek gecise iner. Python ic dongusu kalkinca 1,5 kb x
# 1,5 kb hizalama dakikalardan saniyelere duser.
# ---------------------------------------------------------------------------
def _hw_son(q, t):
    n = len(t); ar = np.arange(n + 1, dtype=np.int32)
    prev = np.zeros(n + 1, dtype=np.int32)
    for i in range(len(q)):
        qi = q[i]
        cost = np.where((t == qi) & (t < 4) & (qi < 4), 0, 1).astype(np.int32)
        cur = np.empty(n + 1, dtype=np.int32)
        cur[0] = prev[0] + 1
        cur[1:] = np.minimum(prev[:-1] + cost, prev[1:] + 1)
        cur = np.minimum.accumulate(cur - ar) + ar
        prev = cur
    return prev

# Yuzde kimlik. Payda DAIMA kisa dizinin uzunlugudur; uzun dizinin fazlaligi
# orani bozmaz. Kutu on gruplamasi (KONS_ESIK = %99) bu sayiya dayanir.
def hw_kimlik(a, b):
    """kisa olani sorgu, uzun olanin icine; donus: yuzde kimlik"""
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    d = int(_hw_son(enc(q), enc(t)).min())
    return round(100.0 * (1 - d / max(len(q), 1)), 2)

# ----------------------------------------------------------------- in-silico PCR
# -------------------------------------------------------------------------
# BIR KUTUNUN HAM OKUMALARINDA IN-SILICO PCR.
#
# Olcut: <=max_mm uyumsuzluk VE 3' son iki baz TAM eslesme. Son iki baz sarti
# kozmetik degildir - polimeraz 3' ucu tutmayan primerden uzatma baslatmaz.
#
# GUVERCIN YUVASI (pigeonhole) TOHUMLAMASI - NEDEN KAYIPSIZ
# Primer, max_mm+1 tane ORTUSMEYEN bloga bolunur. Dizide en fazla max_mm
# uyumsuzluk varsa, uyumsuzluklar en fazla max_mm ayri bloga dagilabilir; geriye
# EN AZ BIR blok kalir ve o blok TAM eslesmek ZORUNDADIR. Dolayisiyla bloklardan
# herhangi birinin tam eslesmesini arayan bir tohumlama, olcutu saglayan HICBIR
# baglanma yerini kaciramaz.
#
# Bu bir SEZGISEL HIZLANDIRMA DEGIL, bir GARANTIDIR. Ayrimi onemlidir: sezgisel
# bir tohumlayici "muhtemelen bulur" der ve kacirdiklarini haber vermez - panelin
# eski okuma motorundaki tohum hatasi tam olarak boyle sessiz kayip uretiyordu.
# Burada aday kumesi tohumla DARALTILIR, karar ise adaylarin tamaminda tam
# uyumsuzluk sayimiyla verilir.
# -------------------------------------------------------------------------
class Kutu:
    """Ham okumalarda in-silico PCR. Olcut: <=max_mm uyumsuzluk + 3' son 2 baz TAM.
    Guvercin yuvasi tohumlamasi: primer max_mm+1 ORTUSMEYEN bloga bolunur, en az
    biri tam tutmak zorundadir -> KAYIPSIZ (mm<=1 ve mm<=3 icin)."""
    def __init__(self, path, nmax=3000):
        rs = []
        op = open
        with op(path, 'rt', errors='ignore') as fh:
            for k, line in enumerate(fh):
                if k % 4 == 1:
                    s = line.strip().upper()
                    if OKUMA_MIN <= len(s) <= OKUMA_MAX: rs.append(s)
        self.n_suzgec = len(rs)
        if nmax and len(rs) > nmax:
            random.seed(TOHUM); rs = random.sample(rs, nmax)
        self.n_okuma = len(rs)
        parts = []; sid = []
        for i, s in enumerate(rs):
            for ss in (s, rc(s)): parts.append(ss); sid.append((i, len(ss)))
        blob = 'N'.join(parts); self.E = enc(blob)
        off = []; p = 0
        for (rid, L) in sid: off.append((p, p + L, rid)); p += L + 1
        self.starts = np.array([o[0] for o in off], dtype=np.int64)
        self.ends = np.array([o[1] for o in off], dtype=np.int64)
        self.rid = np.array([o[2] for o in off], dtype=np.int32)
        self.idx = {}
        E = self.E.astype(np.int64); n = len(E)
        for kk in (9, 5):
            code = np.zeros(n - kk + 1, dtype=np.int64); bad = np.zeros(n - kk + 1, dtype=bool)
            for j in range(kk):
                seg = E[j:n - kk + 1 + j]
                code = code * 4 + np.where(seg < 4, seg, 0); bad |= seg >= 4
            code[bad] = -1
            ok = np.nonzero(code >= 0)[0]
            pos = ok[np.argsort(code[ok], kind='stable')]
            self.idx[kk] = (pos, code[pos])
            del code, bad
        del E
        self.okumalar = rs
    @staticmethod
    def _kod(s, k):
        c = 0
        for ch in s[:k]:
            v = 'ACGT'.find(ch)
            if v < 0: return -1
            c = c * 4 + v
        return c
    # -----------------------------------------------------------------------
    # Bir oligonun butun baglanma yerlerini bulur (guvercin yuvasi tohumlamasi).
    #
    # 1) Oligo max_mm+1 ortusmeyen bloga bolunur -> en az biri tam tutmak zorunda.
    # 2) Blok uzunlugundan kucuk esit bir k secilir (mm<=1 icin 9, mm<=3 icin 5);
    #    k blok boyunu asarsa tohum blogun disina tasar ve GARANTI BOZULUR, o
    #    yuzden k daima min(blok boyu) ile sinirlanir.
    # 3) Her blogun k-mer kodu, onceden kurulmus sirali indekste ikili aramayla
    #    bulunur; adaylarin birlesimi alinir.
    # 4) Aday konumlar okuma sinirlarina kirpilir (bir okuma bitip digeri baslarken
    #    araya konan N ayraci yuzunden okumalar arasi eslesme uretilmemeli).
    # 5) Karar: TAM uyumsuzluk sayimi <= max_mm VE son iki bazin ikisi de tam.
    #    Tohum yalnizca ADAY TOPLAR; eleme burada, tam sayimla yapilir.
    # -----------------------------------------------------------------------
    def _yerler(self, olig, max_mm):
        L = len(olig); nb = max_mm + 1
        kes = [round(i * L / nb) for i in range(nb + 1)]
        bloklar = [(kes[i], kes[i + 1] - kes[i]) for i in range(nb) if kes[i + 1] > kes[i]]
        k = min(min(b[1] for b in bloklar), 9 if max_mm <= 1 else 5)
        k = max(k, 4)
        if k not in self.idx: k = min(self.idx, key=lambda x: abs(x - k))
        k = min(k, min(b[1] for b in bloklar))
        if k not in self.idx: k = min(self.idx)
        pos, key = self.idx[k]
        cands = []
        for off, _ in bloklar:
            c = self._kod(olig[off:off + k], k)
            if c < 0: continue
            a = np.searchsorted(key, c, 'left'); b = np.searchsorted(key, c, 'right')
            if b > a: cands.append(pos[a:b] - off)
        if not cands: return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
        S = np.unique(np.concatenate(cands))
        si = np.searchsorted(self.starts, S, 'right') - 1
        ok = si >= 0; S = S[ok]; si = si[ok]
        ok = (S >= self.starts[si]) & (S + L <= self.ends[si]); S = S[ok]; si = si[ok]
        if S.size == 0: return S, si
        O = enc(olig)
        blok = self.E[S[:, None] + np.arange(L)[None, :]]
        mm = (blok != O[None, :]).sum(1)
        keep = (mm <= max_mm) & (blok[:, L - 1] == O[L - 1]) & (blok[:, L - 2] == O[L - 2])
        return S[keep], si[keep]
    # -----------------------------------------------------------------------
    # Kac okuma urun veriyor? Ileri primer sense yonde, geri primerin TERS
    # TUMLEYENI ayni okumada aranir; ikisi AYNI okumada, dogru sirada (geri primer
    # ileri primerin bittigi yerden sonra) ve lo-hi bp mesafede olmalidir.
    #
    # Sayilan sey URUN VEREN OKUMA SAYISIDIR, baglanma yeri sayisi degil: bir
    # okumada birden fazla gecerli cift bulunsa da o okuma bir kez sayilir
    # (veren bir kume). Aksi halde tekrarli bolgeler orani sisirirdi.
    # -----------------------------------------------------------------------
    def pcr(self, F, R, lo=60, hi=400, max_mm=1):
        Fs, Fi = self._yerler(F, max_mm)
        if Fs.size == 0: return 0, self.n_okuma
        rr = rc(R); Rs, Ri = self._yerler(rr, max_mm)
        if Rs.size == 0: return 0, self.n_okuma
        from collections import defaultdict
        rmap = defaultdict(list)
        for p, i in zip(Rs.tolist(), Ri.tolist()): rmap[i].append(p)
        veren = set(); LF = len(F); LR = len(R)
        for p, i in zip(Fs.tolist(), Fi.tolist()):
            for q in rmap.get(i, ()):
                n = q + LR - p
                if lo <= n <= hi and q >= p + LF:
                    veren.add(int(self.rid[i])); break
        return len(veren), self.n_okuma

# -------------------------------------------------------------------------
# WILSON SKOR ARALIGI - NEDEN HAM ORAN KULLANILMIYOR
#
# Ham oran k/n kucuk orneklemde yaniltici olur: 3 okumanin 3'u urun verdiyse ham
# oran %100'dur, ama bu sayinin arkasinda neredeyse hicbir kanit yoktur. Ayni
# sekilde 200 okumanin 0'i urun verdiyse ham oran %0'dir ve "hic capraz yok"
# izlenimi verir - oysa gercek oran %1,5 olabilir.
#
# Wilson araligi bu belirsizligi sayiya doker ve HER ZAMAN MUHAFAZAKAR TARAF
# secilir:
#   uye tarafi  -> ALT sinir (hedefi gorme basarisini olabilecek en dusuk tahmin)
#   rakip tarafi-> UST sinir (capraz riskini olabilecek en yuksek tahmin)
# Ayrim kati bu ikisinin oranidir, yani daima en kotu senaryoyu olcer. Bu
# secim yuzunden sig kutular DUSUK kat verir - bu bir hata degil, kanitin
# azligidir; olcum derinliginin butun satirlarda ayni tutulmasinin sebebi de budur.
# -------------------------------------------------------------------------
def wilson(k, n, z=1.96):
    import math
    if n == 0: return (0.0, 1.0)
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d
    s = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - s), min(1.0, c + s))

# ----------------------------------------------------------------- yardimcilar
# Kontrol noktasi yazimi ATOMIKtir: once .tmp dosyasina yazilir, sonra os.replace
# ile yerine konur. Kosu yazma aninda kesilirse yarim JSON kalmaz - yarim kontrol
# noktasi bir sonraki kosuyu sessizce bozardi.
def ck_yaz(yol, veri):
    tmp = yol + '.tmp'
    json.dump(veri, open(tmp, 'w', encoding='utf-8'))
    os.replace(tmp, yol)

def ck_oku(yol, varsayilan):
    if os.path.exists(yol):
        try: return json.load(open(yol, encoding='utf-8'))
        except Exception: pass
    return varsayilan

def sinifi(kutu): return kutu.split('_')[0].split('-')[0]

# fastq dosya adi -> kutu adi. A1-1 orneginin dosyalari alt cizgili adlandirilmis
# (A1_1_reads_2223.fastq); normalize edilmezse o kutular hic taninmaz ve sessizce
# olcum disi kalirdi.
def kutu_adi(dosya):
    """fastq dosya adindan kutu adi. A1-1 orneginde alt cizgili adlandirma var
    (A1_1_reads_2223.fastq) - normalize edilir, yoksa o kutular TANINMAZ."""
    b = os.path.basename(dosya)
    for uz in ('.fastq.gz', '.fastq'):
        if b.endswith(uz): b = b[:-len(uz)]
    b = b.replace('_reads_', '-reads_')
    m = re.match(r'^([A-Za-z0-9]+)[-_](\d+)-reads_(\d+)$', b)
    if m: return '%s-%s_%s' % (m.group(1), m.group(2), m.group(3))
    m = re.match(r'^(.+)-reads_(\d+)$', b)
    if m: return '%s_%s' % (m.group(1).replace('_', '-'), m.group(2))
    return b

# ----------------------------------------------------------------- ana akis
# -------------------------------------------------------------------------
# ANA AKIS - bes adim, hepsi kontrol noktali:
#   0) envanter        : konsensus ve fastq dosyalarini eslestir.
#   1) kimlik matrisi  : sinif ici butun konsensus ciftlerini hizala (_ck_kimlik).
#   2) siniflama       : on gruplar + AYIRT EDICI k-mer + okuma atamasi (_ck_icerik).
#   3) uyelik turetme  : uye / karisik / rakip kumelerini yeniden kur.
#   4) yeniden olcum   : butun panel ciftleri HEM mm<=1 HEM mm<=3 (_ck_olcum).
# -------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--kok', default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument('--nmax', type=int, default=3000)
    ap.add_argument('--sifirla', action='store_true', help='kontrol noktalarini sil, bastan kos')
    a = ap.parse_args()
    KOK = a.kok
    CIK = os.path.join(KOK, 'engine_SONUC')
    os.makedirs(CIK, exist_ok=True)
    if a.sifirla:
        for f in glob.glob(os.path.join(CIK, '_ck_*.json')): os.remove(f)
        print('kontrol noktalari silindi.')

    KONS = os.path.join(KOK, 'konsensus_kanonik')
    FQ = os.path.join(KOK, 'fastq files')
    for p, adi in ((KONS, 'konsensus_kanonik'), (FQ, 'fastq files')):
        if not os.path.isdir(p): sys.exit('HATA: %s klasoru yok: %s' % (adi, p))

    print('=' * 70); print('  UYELIGI OLCULEN KIMLIKTEN YENIDEN TURETME'); print('=' * 70)
    print('  Proje  : %s' % KOK); print('  Cikti  : %s' % CIK); print()

    # --- 0. envanter
    kons = {}
    for p in glob.glob(os.path.join(KONS, '*.kanonik.fa')):
        kons[os.path.basename(p).replace('.kanonik.fa', '')] = temizle(list(dict(fasta(p)).values())[0])
    fq = {}
    for p in glob.glob(os.path.join(FQ, '*', '*.fastq')) + glob.glob(os.path.join(FQ, '*', '*.fastq.gz')):
        fq[kutu_adi(p)] = p
    ortak = sorted(set(kons) & set(fq))
    print('  Konsensus %d, fastq %d, eslesen kutu %d' % (len(kons), len(fq), len(ortak)))
    eksik = sorted(set(fq) - set(kons))
    if eksik: print('  UYARI - konsensusu olmayan fastq: %s' % ', '.join(eksik))
    print()

    # ADIM 1 - sinif ici butun konsensus ciftlerinin kimlik matrisi. Yalniz AYNI
    # sinif icinde karsilastirilir: farkli siniflar zaten farkli lokuslardir ve
    # aralarindaki kimlik bir sey ifade etmez. Matris kontrol noktasina yazilir,
    # kosu kesilirse bastan hesaplanmaz.
    # --- 1. kimlik matrisi
    ckm = os.path.join(CIK, '_ck_kimlik.json')
    KM = ck_oku(ckm, {})
    siniflar = {}
    for k in kons: siniflar.setdefault(sinifi(k), []).append(k)
    ciftler = []
    for s, v in siniflar.items():
        v = sorted(v)
        for i in range(len(v)):
            for j in range(i + 1, len(v)): ciftler.append((v[i], v[j]))
    yeni = [c for c in ciftler if '%s|%s' % c not in KM]
    if yeni:
        print('  [1/4] Konsensus kimlik matrisi: %d cift hesaplanacak' % len(yeni))
        t0 = time.time()
        for n, (x, y) in enumerate(yeni, 1):
            KM['%s|%s' % (x, y)] = hw_kimlik(kons[x], kons[y])
            if n % 25 == 0 or n == len(yeni):
                ck_yaz(ckm, KM)
                print('        %d/%d  (%.0f sn)' % (n, len(yeni), time.time() - t0), flush=True)
        ck_yaz(ckm, KM)
    else:
        print('  [1/4] Kimlik matrisi kontrol noktasindan okundu (%d cift)' % len(KM))

    # -----------------------------------------------------------------------
    # ADIM 2 - AYIRT EDICI k-mer kumesi. Bir on grubun k-mer'lerinden DIGER butun
    # gruplarin k-mer'leri CIKARILIR; geriye yalnizca o gruba OZGU olanlar kalir.
    #
    # Bu cikarma SARTTIR. Yapilmazsa okumalar korunmus bolgeler (18S, 5.8S, LSU
    # cekirdegi) uzerinden eslesir ve gruplar birbirine karisir - olculdu:
    # Trichoderma okumalarinin %70'i Petriella'ya atanmisti. Korunmus bolge butun
    # gruplarda AYNI oldugu icin ayrimda ise yaramaz; sahte yuksek benzerlik tam da
    # oradan gelir.
    #
    # Atama esigi NORM_ESIK = 0,30 normalize paydir (grubun ayirt edici kume boyu
    # ile okumanin k-mer sayisinin kucugune bolunur), bagimsiz in-silico PCR ile
    # kalibre edilmistir.
    # -----------------------------------------------------------------------
    # --- 2. on gruplar + ayirt edici k-mer + okuma siniflamasi
    cki = os.path.join(CIK, '_ck_icerik.json')
    IC = ck_oku(cki, {})
    GRUPLAR = {}
    for s, probs in sorted(siniflar.items()):
        probs = sorted(probs)
        par = {k: k for k in probs}
        def bul(x):
            while par[x] != x: par[x] = par[par[x]]; x = par[x]
            return x
        for (x, y), v in ((tuple(k.split('|')), v) for k, v in KM.items()):
            if x in par and y in par and v >= KONS_ESIK:
                rx, ry = bul(x), bul(y)
                if rx != ry: par[rx] = ry
        g = {}
        for k in probs: g.setdefault(bul(k), []).append(k)
        GRUPLAR[s] = sorted([sorted(v) for v in g.values()], key=lambda v: v[0])
    if any(k not in IC for k in ortak):
        print('  [2/4] Ayirt edici k-mer ile okuma siniflamasi')
        for s, gruplar in sorted(GRUPLAR.items()):
            hedefler = [k for k in ortak if sinifi(k) == s and k not in IC]
            if not hedefler: continue
            ks = []
            for grp in gruplar:
                st = set()
                for k in grp:
                    for ss in (kons[k], rc(kons[k])):
                        for i in range(len(ss) - K + 1): st.add(ss[i:i + K])
                ks.append(st)
            ayirt = []
            for i, st in enumerate(ks):
                dis = set()
                for j, s2 in enumerate(ks):
                    if i != j: dis |= s2
                ayirt.append(st - dis)
            boy = [len(x) for x in ayirt]
            idx = {}
            for i, st in enumerate(ayirt):
                for x in st: idx[x] = i
            print('        sinif %-3s: %d grup, ayirt edici kume %s' % (s, len(gruplar), boy))
            for kb in hedefler:
                rs = []
                with open(fq[kb], 'rt', errors='ignore') as fh:
                    for n, line in enumerate(fh):
                        if n % 4 == 1:
                            x = line.strip().upper()
                            if OKUMA_MIN <= len(x) <= OKUMA_MAX: rs.append(x)
                random.seed(TOHUM)
                if len(rs) > SINIF_ORNEK: rs = random.sample(rs, SINIF_ORNEK)
                say = [0] * len(gruplar)
                for r in rs:
                    nk = len(r) - K + 1
                    if nk <= 0: continue
                    c = [0] * len(gruplar)
                    for i in range(nk):
                        j = idx.get(r[i:i + K])
                        if j is not None: c[j] += 1
                    for gi in range(len(gruplar)):
                        if c[gi] / max(1, min(boy[gi], nk)) >= NORM_ESIK: say[gi] += 1
                IC[kb] = dict(n=len(rs), sinif=s,
                              pay={gruplar[i][0]: round(100.0 * say[i] / max(len(rs), 1), 1)
                                   for i in range(len(gruplar)) if say[i]})
                ck_yaz(cki, IC)
                print('          %-18s n=%4d  %s' % (kb, len(rs),
                      ', '.join('%s:%.0f%%' % x for x in sorted(IC[kb]['pay'].items(), key=lambda y: -y[1])[:2]) or '-'),
                      flush=True)
    else:
        print('  [2/4] Okuma siniflamasi kontrol noktasindan okundu')

    # -----------------------------------------------------------------------
    # ADIM 3 - UYELIGIN YENIDEN KURULMASI.
    # Once mevcut uye kutularin baskin gruplari toplanir (hg = hedef gruplar),
    # sonra sinifin BUTUN kutulari bu gruplara gore yeniden dagitilir:
    #   pay >= UYE_ESIK      -> UYE
    #   KARISIK_ESIK .. UYE  -> KARISIK (ne uye ne rakip)
    #   altinda              -> RAKIP
    # Evrensel hedeflerde sinifin tamami uyedir; ayirma yapilmaz.
    # -----------------------------------------------------------------------
    # --- 3. uyelik yeniden turetme
    print('  [3/4] Uyelik yeniden turetiliyor')
    uyelik_tsv = os.path.join(KOK, 'screening', 'hedef_uyelik.tsv')
    panel_tsv = os.path.join(KOK, 'primer_final', 'devir_ciftleri_20260802_sonrotus_TESLIM.tsv')
    if not os.path.exists(uyelik_tsv): sys.exit('HATA: %s yok' % uyelik_tsv)
    if not os.path.exists(panel_tsv): sys.exit('HATA: %s yok' % panel_tsv)
    uyelik = {}
    for line in open(uyelik_tsv, encoding='utf-8'):
        if line.startswith('#') or line.startswith('hedef\t'): continue
        p = line.rstrip('\n').split('\t')
        if len(p) >= 2 and p[0].strip(): uyelik[p[0]] = p[1]
    rows = [r for r in csv.DictReader(open(panel_tsv, encoding='utf-8'), delimiter='\t')
            if (r.get('Hedef') or '').strip() in uyelik]
    UY = {}
    for r in rows:
        hd = r['Hedef'].strip()
        sn = [x.strip() for x in (r.get('Amplikon sinifi') or '').split('/') if x.strip()] or list(siniflar)
        ut = uyelik[hd]; evrensel = ut.strip().startswith('*')
        tum = [k for k in ortak if sinifi(k) in sn]
        if evrensel:
            UY[hd] = dict(sinif=sn, eski=tum, yeni=tum, karisik=[], rakip=[], eklenen=[], cikan=[], evrensel=True)
            continue
        tx = set(x.strip() for x in ut.split(',') if x.strip())
        eski = [k for k in tum if k.split('_')[1] in tx]
        hg = set()
        for k in eski:
            for g, v in (IC.get(k, {}) or {}).get('pay', {}).items():
                if v >= UYE_ESIK: hg.add(g)
        yeni_u = []; kar = []; rak = []
        for k in tum:
            pay = (IC.get(k, {}) or {}).get('pay', {})
            ic = max([pay.get(g, 0.0) for g in hg], default=0.0)
            en = max(pay.items(), key=lambda x: x[1]) if pay else (None, 0.0)
            if k in eski:
                # ---------------------------------------------------------------
                # KANIT YOKLUGU KANIT SAYILMAZ.
                # Halihazirda uye olan bir kutu ancak POZITIF olcum kanitiyla yer
                # degistirir: okumalarinin baskin grubu UYE_ESIK'i gecmeli VE o
                # grup hedef gruplarindan biri OLMAMALIDIR. Olcum bu kutu icin
                # sinyal uretemediyse (pay tablosu bos ya da zayif) kutu UYE
                # KALIR.
                #
                # Tersi yapilsaydi - "kanit yoksa cikar" - teshis edilen hatanin
                # aynisi ters yonde tekrarlanirdi: az okumali kutular sirf sessiz
                # olduklari icin rakip hanesine dusup ayrim katlarini bozarlardi.
                # ---------------------------------------------------------------
                # POZITIF kanit yoksa eski durum korunur
                if en[1] >= UYE_ESIK and en[0] not in hg: rak.append(k)
                else: yeni_u.append(k)
            else:
                if ic >= UYE_ESIK: yeni_u.append(k)
                elif ic >= KARISIK_ESIK: kar.append(k)
                else: rak.append(k)
        UY[hd] = dict(sinif=sn, eski=sorted(eski), yeni=sorted(yeni_u), karisik=sorted(kar),
                      rakip=sorted(rak), eklenen=sorted(set(yeni_u) - set(eski)),
                      cikan=sorted(set(eski) - set(yeni_u)), evrensel=False)
    deg = sum(1 for o in UY.values() if o['eklenen'] or o['cikan'])
    print('        %d hedeften %d tanesinin uyeligi degisti' % (len(UY), deg))

    # ADIM 4 - butun panel ciftleri, butun kutularda, HEM mm<=1 (asil olcut) HEM
    # mm<=3 (dayaniklilik olcutu) ile olculur. Kutu basina tek fastq okumasi yapilir
    # ve o kutudaki BUTUN ciftler ayni indeks uzerinden sorulur; her kutu bitince
    # kontrol noktasina yazilir.
    # --- 4. panel ciftlerini yeniden olc
    CF = []
    for r in rows:
        F = (r.get("Ileri primer (5'->3')") or '').strip().upper()
        R = (r.get("Geri primer (5'->3')") or '').strip().upper()
        if not re.fullmatch(r'[ACGT]+', F or '') or not re.fullmatch(r'[ACGT]+', R or ''): continue
        CF.append(dict(hedef=r['Hedef'].strip(), F=F, R=R,
                       urun=(r.get('Urun (bp)') or '').strip()))
    cko = os.path.join(CIK, '_ck_olcum.json')
    # 2026-08-10 DIZI MUHRU. Onbellek kutu adiyla anahtarlaniyor, cift sonuclari
    # ise ciftin SIRA NUMARASI (str(i)) altinda tutuluyordu. Iki ayri hata:
    #   1) bir ciftin dizisi degisince ayni sira numarasi okunup ESKI olcum
    #      yeni dizinin yanina yaziliyordu (dizi yeni, sayi eski);
    #   2) panelde cift eklenip cikarilinca sira kayiyor ve sayilar YANLIS
    #      cifte atanabiliyordu.
    # Cozum: cift anahtari sira degil, F+R dizisinin ozeti. Ayrica dosyanin
    # basina bir muhur yazilir; muhur tutmazsa onbellek bastan kurulur.
    import hashlib as _hl

    def _ck_anahtar(c):
        return _hl.md5((c['F'] + '|' + c['R']).encode('utf-8')).hexdigest()[:12]

    _muhur = _hl.md5('|'.join(sorted(_ck_anahtar(c) for c in CF))
                     .encode('utf-8')).hexdigest()[:12]
    OL = ck_oku(cko, {})
    if OL.get('_muhur') != _muhur:
        if OL:
            print('  [4/4] onbellek DIZI muhru tutmuyor (kayitli %s, simdi %s) - '
                  'bastan olculuyor.' % (OL.get('_muhur') or 'yok', _muhur))
        OL = {'_muhur': _muhur}
    kalan = [k for k in ortak if k not in OL]
    print('  [4/4] In-silico PCR: %d cift x %d kutu (kalan %d kutu), mm<=1 ve mm<=3' % (len(CF), len(ortak), len(kalan)))
    t0 = time.time()
    for n, kb in enumerate(kalan, 1):
        Kt = Kutu(fq[kb], nmax=a.nmax)
        r = {}
        for i, c in enumerate(CF):
            p1, nn = Kt.pcr(c['F'], c['R'], 60, 400, 1)
            p3, _ = Kt.pcr(c['F'], c['R'], 60, 400, 3)
            r[_ck_anahtar(c)] = [p1, nn, p3]
        OL[kb] = r; del Kt
        ck_yaz(cko, OL)
        print('        %-18s %d/%d  (%.0f sn)' % (kb, n, len(kalan), time.time() - t0), flush=True)

    # --- ayrim katlari
    # -----------------------------------------------------------------------
    # AYRIM KATI = (uye kutularin EN DUSUK Wilson ALT siniri) /
    #              (rakip kutularin EN YUKSEK Wilson UST siniri)
    # Iki tarafta da en kotu kutu secilir: pay tarafinda hedefi en az goren uye,
    # payda tarafinda en cok capraz veren rakip. Ortalama alinsaydi tek bir kotu
    # kutu kalabaligin icinde erir ve gercek risk gorunmezdi.
    #
    # ENKOTU_ASGARI (150 okuma) altindaki rakip kutular paydaya girmez: o
    # derinlikte Wilson ust siniri neredeyse her zaman tavana vurur ve kat sayisi
    # rakibin gercek davranisini degil sadece okuma azligini olcerdi.
    #
    # KAPSAM ayri bir eksendir: uye kutularin kaci >=%20 urun veriyor. Ayrim
    # yuksek ama kapsam dusukse cift ozguldur ama hedefin tamamini gormez - iki
    # sorun birbirine karistirilmamalidir.
    # -----------------------------------------------------------------------
    def hesap(ci, uye, rakip, mm):
        uy = []; rk = []
        for kb, r in OL.items():
            if kb == '_muhur' or not isinstance(r, dict):
                continue
            v = r.get(_ck_anahtar(CF[ci]))
            if not v: continue
            p, nn = (v[0], v[1]) if mm == 1 else (v[2], v[1])
            if kb in uye: uy.append((kb, p, nn))
            elif kb in rakip: rk.append((kb, p, nn))
        if not uy: return None
        ua = min(wilson(p, nn)[0] for _, p, nn in uy)
        kaps = sum(1 for _, p, nn in uy if nn and p / nn >= KAPSAM_ESIGI)
        enk = None
        for kb, p, nn in rk:
            if nn < ENKOTU_ASGARI: continue
            hi = wilson(p, nn)[1]
            if enk is None or hi > enk[1]: enk = (kb, hi, p, nn)
        return dict(uye_alt=round(100 * ua, 2), kapsam='%d/%d' % (kaps, len(uy)),
                    enkotu=enk[0] if enk else '-', enkotu_ust=round(100 * enk[1], 2) if enk else 0.0,
                    kat=round(ua / enk[1], 2) if enk and enk[1] > 0 else None)
    SON = []
    for ci, c in enumerate(CF):
        o = UY.get(c['hedef'])
        if not o: continue
        tum = [k for k in OL if sinifi(k) in o['sinif']]
        eu = set(o['eski']); er = set(tum) - eu
        yu = set(o['yeni']); kar = set(o['karisik'])
        rA = set(tum) - yu; rB = set(tum) - yu - kar
        d = dict(hedef=c['hedef'], sinif='/'.join(o['sinif']), F=c['F'], R=c['R'], urun=c['urun'],
                 eski_n=len(eu), yeni_n=len(yu), kar_n=len(kar),
                 eklenen=';'.join(o['eklenen']), cikan=';'.join(o['cikan']))
        for mm in (1, 3):
            d['eski_mm%d' % mm] = hesap(ci, eu, er, mm)
            d['A_mm%d' % mm] = hesap(ci, yu, rA, mm)
            d['B_mm%d' % mm] = hesap(ci, yu, rB, mm)
        SON.append(d)

    # --- ciktilar
    def kat(x): return '' if not x or x['kat'] is None else x['kat']
    p1 = os.path.join(CIK, 'ciftler_yeniden_olcum.tsv')
    with open(p1, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'sinif', 'F', 'R', 'urun_bp', 'eski_uye_n', 'yeni_uye_n', 'karisik_n',
                    'eklenen', 'cikan', 'eski_kat_mm1', 'yeniA_kat_mm1', 'yeniB_kat_mm1',
                    'eski_kat_mm3', 'yeniA_kat_mm3', 'yeniB_kat_mm3', 'yeniA_uye_alt',
                    'yeniA_kapsam', 'yeniA_enkotu_kutu', 'yeniA_enkotu_ust'])
        for d in SON:
            A = d['A_mm1'] or {}
            w.writerow([d['hedef'], d['sinif'], d['F'], d['R'], d['urun'], d['eski_n'], d['yeni_n'],
                        d['kar_n'], d['eklenen'], d['cikan'],
                        kat(d['eski_mm1']), kat(d['A_mm1']), kat(d['B_mm1']),
                        kat(d['eski_mm3']), kat(d['A_mm3']), kat(d['B_mm3']),
                        A.get('uye_alt', ''), A.get('kapsam', ''), A.get('enkotu', ''), A.get('enkotu_ust', '')])
    p2 = os.path.join(CIK, 'kutu_olculen_kimlik.tsv')
    with open(p2, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['kutu', 'sinif', 'kraken_taxid', 'baskin_grup', 'baskin_pay', 'ikinci_grup', 'ikinci_pay', 'yorum'])
        for kb in sorted(IC):
            pay = IC[kb]['pay']; s = IC[kb]['sinif']; t = kb.split('_')[1]
            srt = sorted(pay.items(), key=lambda x: -x[1])
            b1 = srt[0] if srt else ('', 0); b2 = srt[1] if len(srt) > 1 else ('', 0)
            yorum = ''
            if b1[0] and b1[1] >= UYE_ESIK and b1[0].split('_')[1] != t:
                yorum = 'KRAKEN ETIKETI YANLIS -> %s' % b1[0]
            elif not srt:
                yorum = 'olcum sinyali yok - eski durum korundu'
            w.writerow([kb, s, t, b1[0], b1[1], b2[0], b2[1], yorum])
    p3 = os.path.join(CIK, 'engine_TURETME.md')
    with open(p3, 'w', encoding='utf-8') as fh:
        fh.write('# Uyeligin olculen kimlikten yeniden turetilmesi\n\n')
        fh.write('Uretim: %s\n\n' % time.strftime('%Y-%m-%d %H:%M:%S'))
        fh.write('Olcut: ayirt edici %d-mer, normalize esik %.2f, uye >=%%%.0f, karisik %%%.0f-%.0f.\n'
                 % (K, NORM_ESIK, UYE_ESIK, KARISIK_ESIK, UYE_ESIK))
        fh.write('In-silico PCR: <=1 ve <=3 uyumsuzluk + 3\' son 2 baz TAM, kutu basina en cok %d okuma.\n\n' % a.nmax)
        fh.write('A senaryosu = karisik kutular RAKIP sayilir (guvenli).  B = karisik kutular dislanir.\n\n')
        fh.write('| hedef | eski uye | yeni uye | karisik | eski kat (mm1) | yeni A | yeni B | yeni A (mm3) |\n')
        fh.write('|---|---|---|---|---|---|---|---|\n')
        for d in SON:
            fh.write('| %s | %d | %d | %d | %s | %s | %s | %s |\n' % (
                d['hedef'], d['eski_n'], d['yeni_n'], d['kar_n'],
                kat(d['eski_mm1']), kat(d['A_mm1']), kat(d['B_mm1']), kat(d['A_mm3'])))
        fh.write('\n## Uyeligi degisen hedefler\n\n')
        for d in SON:
            if d['eklenen'] or d['cikan']:
                fh.write('- **%s**: eklenen `%s` / cikan `%s`\n' % (d['hedef'], d['eklenen'] or '-', d['cikan'] or '-'))
        fh.write('\n## Kraken etiketi yanlis olan kutular\n\n')
        for kb in sorted(IC):
            pay = IC[kb]['pay']; t = kb.split('_')[1]
            srt = sorted(pay.items(), key=lambda x: -x[1])
            if srt and srt[0][1] >= UYE_ESIK and srt[0][0].split('_')[1] != t:
                fh.write('- `%s` -> okumalarinin %%%.0f i `%s` organizmasina ait\n' % (kb, srt[0][1], srt[0][0]))
        fh.write('\n> Bu betik panel dosyalarina YAZMAZ. Degisiklikleri uygulamak ayri bir istir.\n')
    print()
    print('=' * 70); print('  BITTI'); print('=' * 70)
    print('  %s' % p3); print('  %s' % p1); print('  %s' % p2)
    return 0

if __name__ == '__main__':
    sys.exit(main())
