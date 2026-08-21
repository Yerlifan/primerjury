"""Bagimsiz in-silico PCR motoru (2026-08-02 son etap).
Onceki oturumun engine_gateway.py / kararlar.py kodu KULLANILMAZ - sifirdan yazildi.
Kriter: F ve R primerleri karsilikli yonelimde baglanacak, 3' son 2 baz TAM tutacak,
toplam uyumsuzluk <= max_mm, urun boyu araligi verilen pencerede olacak.
"""
# ---------------------------------------------------------------------------
# ispcr.py — panelin cekirdek in-silico PCR motoru; bir dizide F ve R
#            primerlerinin karsilikli baglanip urun verip vermedigini olcer.
#
# GİRDİ  : fasta dosyalari (read_fasta) ya da dogrudan dizi metni; F ve R primer
#          dizileri (IUPAC dejenere baz kabul edilir). Komut satirindan:
#          python ispcr.py <F> <R> <fasta...>
# ÇIKTI  : dosyaya YAZMAZ, deger dondurur. find_sites -> [(baslangic, mm)];
#          amplify -> [(bas, bit, urun_bp, mm_F, mm_R)]; scan_file ->
#          (toplam_dizi, urun_veren_dizi, boy_dagilimi, vurus_listesi).
#          Yalniz __main__ ekrana ozet basar.
# ÇAĞRAN : screening/engine_gateway.py bu dosyayi ada gore bulup yukler; bulamazsa
#          "HATA: ispcr.py bulunamadi" deyip cikar, yani engine_gateway.py'siz hicbir
#          olcum kosamaz. engine_gateway.py'yi de paketin butun olcum modulleri
#          (sample.py, global_scan.py, reference.py, generator.py, targets.py,
#          build_consensus.py, panel_measurement.py, membership_check.py,
#          self_test.py) ve disaridan verification/recovery_round.py ile
#          protocol/single_protocol_measure.py kullanir. Ayrica ayni klasordeki
#          ara.py, hiza.py, mazei*.py, deg_*.py, mmb_*.py ve
#          engine/scanner.py, pair.py dogrudan "import ispcr" der.
#          Menude her olcum tusunda yuklenir: P (tek protokolle panel olcumu),
#          K (kurtarma), D (dogrulama), I ve G (kimlik), T (P->K->D->I),
#          U (uyelik), H (hizli test), 1-9 (kapsamli arama).
#          full_chain.py bu dosyayi DOGRUDAN cagirmaz ama acilista
#          "engine\ispcr.py" var mi diye bakar; yoksa arac acilmaz.
#
# NEDEN TOHUMSUZ: bu motor find_sites icinde tam kayan pencere tarar, hicbir
# tohum/kisayol varsayimi yoktur; dolayisiyla kayipsizligi tanim geregidir.
# Tohumlu hizli yollar (screening/read_engine.py, scanner.py Havuz)
# bu dosyanin olcutunu taklit etmek zorundadir ve self_test.py ucunu
# birebir karsilastirir.
# ---------------------------------------------------------------------------
import re, sys, os, glob
from collections import defaultdict

# IUPAC dejenere baz tablosu: bir primer harfinin dizide hangi bazlara
# karsilik gelmesine izin verildigi. Ornegin R yazan bir primer konumu hem A
# hem G ile eslesir. Dejenere primerler tek bir oligo degil, oligo KARISIMIDIR;
# arama tarafinda bu tabloyla acilir.
IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT',
}
COMP = str.maketrans('ACGTURYSWKMBDHVNacgturyswkmbdhvn',
                     'TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn')


# Ters tumleyen (reverse complement). PCR'de R primeri sablonun EKSI ipligine
# baglanir; motor yalniz arti iplikte tarama yaptigi icin R'yi aramadan once
# rc(R)'ye cevirmek zorundadir (bkz. amplify).
def rc(s):
    return s.translate(COMP)[::-1]


# Diziyi normalize eder: ACGTUN disindaki her karakter (bosluk, dejenere baz,
# hizalama tiresi, kucuk harf artigi) N'ye cevrilir, U -> T yapilir.
# Onemli: buradaki N'ler encode() icinde -1 olur ve find_sites'ta ASLA
# eslesmez. Yani belirsiz baz sessizce "uyar" sayilmaz, uyumsuzluk sayilir -
# olcum bu yuzden iyimser degil, temkinlidir.
def clean(s):
    return re.sub(r'[^ACGTUN]', 'N', s.upper()).replace('U', 'T')


# Fasta okuyucu. Tum dosyayi bellege almaz, kayit kayit uretir (yield): panelin
# referans/konsensus dosyalari cok kayitli olabiliyor.
def read_fasta(path):
    name, buf = None, []
    with open(path, errors='ignore') as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(buf)
                name, buf = line[1:].rstrip('\n'), []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, ''.join(buf)


import numpy as np

# Bayt -> baz kodu cevrim tablosu. A,C,G,T sirasiyla 0,1,2,3; DIGER HER BAYT
# -1 kalir. -1 olmasi kasitlidir: find_sites'ta "col >= 0" kosulu N'i ve her
# tanimsiz karakteri otomatik olarak uyumsuz sayar.
_B = np.zeros(256, dtype=np.int8)
_B[:] = -1
for _i, _c in enumerate('ACGT'):
    _B[ord(_c)] = _i


# Diziyi int8 dizisine cevirir. Tarama saf Python dongusu yerine numpy sutun
# islemleriyle yapilabilsin diye; hiz buradan gelir, olcut degismez.
def encode(seq):
    a = np.frombuffer(seq.encode(), dtype=np.uint8)
    return _B[a]


# ---------------------------------------------------------------------------
# find_sites — primerin dizide ILERI yonde bagladigi tum baslangic konumlari.
#
# NE HESAPLAR: her olasi baslangic konumu icin (a) toplam uyumsuzluk sayisi ve
# (b) 3' kritik uclarin tam tutup tutmadigi. Kabul kosulu: mm <= max_mm ve
# (need_tail ise) tail_pos konumlarinin HEPSI tam eslesmis olmali.
#
# NEDEN BOYLE HESAPLANIR - TOHUM YOK:
# Bu fonksiyon dizide TAM KAYAN PENCERE tarar. Dis dongu primerin bazlari
# uzerindedir (L adim), her adimda o bazin dizideki butun hizalamalari tek
# hamlede karsilastirilir; yani hicbir baslangic konumu elenmeden atlanmaz.
# Bu, kayipsizligin tanim geregi saglanmasi demektir.
#
# GUVERCIN YUVASI (PIGEONHOLE) ILE ILISKISI: tohumlu motorlar hiz icin once
# kisa bir parcanin TAM eslestigi yerlere bakar. Bunun kayipsiz olmasi ancak
# su garantiyle mumkundur: en fazla m uyumsuzluga izin verilen bir aramada
# dizi (m + 1) tane ORTUSMEYEN parcaya bolunurse, en fazla m parca bozulabilir
# ve geriye kalan EN AZ BIR parca tam eslesmek ZORUNDADIR. Bu bir garantidir,
# sezgisel bir hizlandirma degildir. Parca sayisi m + 1'in ALTINA duserse
# garanti coker ve arama sessizce baglanma yeri kacirmaya baslar - hata da
# atmaz, sadece sayilar kucuk cikar. Bu dosya hic tohum kullanmadigi icin o
# riski tasimaz; tam da bu yuzden tohumlu yollarin dogrulugu (read_engine.py,
# scanner.py) bu fonksiyona ve brute_force.py'ye karsi sinanir.
#
# 3' UC KURALI (tail_pos): polimeraz primeri 3' ucundan uzatir. 3' ucta
# uyumsuzluk varsa uzatma pratikte olmaz; bu yuzden 3' son 2 bazin TAM tutmasi
# sart kosulur. Panelin olcutu budur. tail_pos negatif verilirse (-1, -2)
# primerin sonu, (0, 1) verilirse basi kritik olur - ikincisi ters primerin
# tumleyeni icin gerekir (bkz. amplify).
# ---------------------------------------------------------------------------
def find_sites(enc, primer, max_mm, need_tail=True, tail_pos=(-1, -2)):
    """primer'in dizide ileri yonde baglandigi 0-tabanli baslangiclar (numpy).
    tail_pos: 3' kritik konumlarin primer icindeki indeksleri."""
    L, n = len(primer), enc.shape[0]
    if n < L:
        return []
    # m = olasi baslangic konumu sayisi. Her konum ayri bir aday hizalamadir.
    m = n - L + 1
    # mm[i]  : i konumundaki hizalamanin toplam uyumsuzluk sayisi
    # tail_ok: i konumunda 3' kritik uclarin hepsi tam tuttu mu
    mm = np.zeros(m, dtype=np.int16)
    tail_ok = np.ones(m, dtype=bool)
    # Negatif indeksler (-1, -2) primer uzunluguna gore pozitife cevrilir.
    tpos = {(p % L) for p in tail_pos}
    # Dis dongu PRIMER bazlari uzerinde (L adim), ic islem numpy vektoru.
    # Boylece butun baslangic konumlari es zamanli degerlendirilir; konum
    # atlanmaz, tohum elemesi yapilmaz.
    for k, p in enumerate(primer):
        allowed = [ 'ACGT'.index(c) for c in IUPAC.get(p, 'ACGT') ]
        # col: her aday hizalamanin k. sirasindaki dizi bazi.
        col = enc[k:k + m]
        ok = np.zeros(m, dtype=bool)
        # Dejenere primer bazi birden cok baza izin verebilir; hepsi denenir.
        for a in allowed:
            ok |= (col == a)
        ok &= (col >= 0)          # N/bosluk = uyumsuz
        mm += (~ok)
        # 3' kritik konumsa tolerans YOK: tek uyumsuzluk konumu tumden eler.
        if k in tpos:
            tail_ok &= ok
    sel = mm <= max_mm
    if need_tail:
        sel &= tail_ok
    idx = np.nonzero(sel)[0]
    # (baslangic, o baglanmadaki uyumsuzluk sayisi) ciftleri dondurulur;
    # uyumsuzluk sayisi cagiran tarafta "en iyi urun"u secmek icin kullanilir.
    return list(zip(idx.tolist(), mm[idx].tolist()))


# ---------------------------------------------------------------------------
# amplify — IN-SILICO PCR. Bir dizide F ve R primerlerinin gercekten bir urun
#           verip veremeyecegini olcer.
#
# URUN SAYILMA SARTLARI (gercek PCR'in dort kosulu):
#   1) F, dizide ileri yonde baglanmali (find_sites, 3' son 2 baz tam).
#   2) R, KARSI iplige baglanmali. Motor yalniz arti ipligi taradigi icin R'nin
#      kendisi degil, ters tumleyeni rc(R) aranir. rc(R) arti iplikte
#      bulunuyorsa, R gercekte eksi iplige baglaniyor demektir.
#   3) IKI PRIMER BIRBIRINE BAKMALI: F soldan saga, R sagdan sola uzar.
#      j >= i + len(fwd) kosulu R'nin baglanma yerinin F'nin SAGINDA ve F ile
#      ORTUSMEDEN olmasini sart kosar. Bu olmadan ayni bolgeye binmis iki
#      baglanma sahte bir "urun" uretirdi; ayrica yanlis siradaki (R solda,
#      F sagda) bir cift PCR'de hicbir sey vermez, burada da verilmez.
#   4) URUN BOYU PENCEREYE girmeli: lo <= size <= hi. Bunun qPCR karsiligi
#      var - config.py'de URUN_IDEAL (60-150 bp), URUN_KABUL (150-250),
#      URUN_MUTLAK_UST 400. Kisa urun protokoldeki 30 sn uzatmada tamamlanir
#      ve SYBR Green sinyali temiz cikar; cok uzun urun ne verimli cogalir ne
#      de jelde/erime egrisinde beklenen yerde gorunur. Pencere yoksa motor
#      birbirinden 5 kb uzaktaki iki baglanmayi da "urun" sayardi.
#
# 3' UC INDEKSININ YON DEGISTIRMESI: R primerinin 3' ucu, rc(R) dizisinin 5'
# UCUDUR. Bu yuzden F icin tail_pos=(-1,-2) verilirken rc(R) icin (0,1)
# verilir. Bu ayrinti atlanirsa 3' kurali R tarafinda yanlis uctan uygulanir.
#
# YON UYARISI: bu fonksiyon verilen dizinin YALNIZ ARTI IPLIGINI tarar,
# rc(seq)'i kendiliginden denemez. Nanopore okumalari cift yonludur ve bir
# konsensus ters yonde uretilmis olabilir; ters yonlu bir dizide bu fonksiyon
# urunlerin TAMAMINI kaybeder ve HATA DA ATMAZ, sessizce 0 doner. Bu yuzden:
#   - konsensus tarafinda yon KANONIK KAYNAKTAN okunur
#     (screening/hedefler.konsensusler -> konsensus_kanonik/INDEKS.tsv,
#     hepsi SENSE). Eski "consensus sequences" klasoru KARISIK yonluydu
#     (71 antisense / 27 sense) ve ona sessizce dusulmesi yasaklanmistir.
#   - ham okuma tarafinda cagiran her okumayi iki yonde de dener
#     (okuma.kutu_pcr ve numune.pcr_kutu icindeki "for seq in (s, rc(s))"),
#     motor.urun_var da ayni sekilde iki yon dener.
# ---------------------------------------------------------------------------
def amplify(seq, fwd, rev, max_mm=3, lo=40, hi=600, need_tail=True, enc=None):
    """Urunleri dondur: (baslangic, bitis, urun_bp, mm_F, mm_R)."""
    if enc is None:
        enc = encode(clean(seq))
    revrc = rc(rev)
    fs = find_sites(enc, fwd, max_mm, need_tail, tail_pos=(-1, -2))
    # F hic baglanmiyorsa R'yi aramaya gerek yok - urun zaten olamaz.
    if not fs:
        return []
    # revrc icin 3' kriteri: rev primerin 3' ucu = revrc'nin 5' ucu -> indeks 0,1
    rs = find_sites(enc, revrc, max_mm, need_tail, tail_pos=(0, 1))
    if not rs:
        return []
    prods = []
    # Butun F x R baglanma cifti denenir; ayni dizide birden cok urun cikabilir.
    for i, mmf in fs:
        for j, mmr in rs:
            # Urun, F'nin baslangicindan rc(R)'nin BITISINE kadardir; yani iki
            # primerin kendisi de urunun icindedir (gercek amplikonda oldugu gibi).
            end = j + len(revrc)
            size = end - i
            # Boy penceresi + primerlerin ortusmeden karsilikli durmasi.
            if lo <= size <= hi and j >= i + len(fwd):
                prods.append((i, end, size, mmf, mmr))
    return prods


# ---------------------------------------------------------------------------
# scan_file — bir fasta dosyasindaki dizilerin kacinda urun ciktigini sayar.
#
# NE HESAPLAR: (toplam, urun_veren, boy_dagilimi, vurus_listesi).
# Oran hesabinin PAYDASI min_len'den kisa diziler ELENDIKTEN sonraki dizi
# sayisidir. Neden: kirpik/parcali bir kayit urunu tasiyacak kadar uzun bile
# degilse onu paydaya koymak orani sahte olarak asagi ceker.
#
# Bir dizi urun verse de vermese de BIR kez sayilir (derinlik degil, varlik
# olcumu). Ayni dizide birden cok urun varsa toplam uyumsuzlugu en kucuk olan
# "en iyi" urun secilir - en olasi baglanma odur.
#
# DIKKAT: burada clean(seq) verilir ama rc(seq) DENENMEZ. Ters yonlu bir
# fasta bu fonksiyonda sessizce 0 urun verir; yon kanonik kaynaktan gelmelidir
# (bkz. amplify basligindaki yon uyarisi).
# ---------------------------------------------------------------------------
def scan_file(path, fwd, rev, max_mm=3, lo=40, hi=600, min_len=0):
    """Bir fasta uzerinde tara. Payda = min_len'den uzun diziler."""
    tot = amp = 0
    sizes = defaultdict(int)
    hits = []
    for name, seq in read_fasta(path):
        # Paydaya girmeyecek kadar kisa kayitlar tamamen atlanir.
        if len(seq) < min_len:
            continue
        tot += 1
        p = amplify(clean(seq), fwd, rev, max_mm, lo, hi)
        if p:
            amp += 1
            # En iyi urun = mm_F + mm_R toplami en kucuk olan.
            best = min(p, key=lambda x: x[3] + x[4])
            sizes[best[2]] += 1
            hits.append((name, best[2], best[3], best[4]))
    return tot, amp, dict(sizes), hits


# Elle kullanim: python ispcr.py <F> <R> <fasta...>
# min_len=1200 - komut satiri kullanimi tam boy 16S/ITS kayitlarini hedefler,
# kisa parcalar paydayi kirletmesin diye.
if __name__ == '__main__':
    fwd, rev = sys.argv[1].upper(), sys.argv[2].upper()
    for path in sorted(sys.argv[3:]):
        t, a, s, _ = scan_file(path, fwd, rev, min_len=1200)
        print(f'{os.path.basename(path):40s} {a:5d}/{t:5d}  {sorted(s.items())[:4]}')
