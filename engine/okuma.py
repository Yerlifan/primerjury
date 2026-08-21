"""Ham okuma duzeyinde hizli in-silico PCR (tohum + dogrulama).
okuma_pcr.py KULLANILMAZ; bagimsiz yazildi. Filtre 200-6000 bp (duzeltilmis).
"""
# ---------------------------------------------------------------------------
# okuma.py — numunenin HAM OKUMALARI (fastq) uzerinde in-silico PCR; ispcr'in
#            olcutunu kullanir ama once kisa bir tohumla aday konum arar.
#
# GİRDİ  : fastq / fastq.gz kutu dosyalari (okumalar), F ve R primer dizileri.
#          Komut satirindan: python okuma.py <F> <R> <fastq...>
# ÇIKTI  : dosyaya YAZMAZ. kutu_pcr -> (toplam_okuma, urun_veren_okuma,
#          urun_boyu_sayaci). __main__ kutu basina yuzde satiri basar.
# ÇAĞRAN : KAPSAMLI_ARAMA/motor.py bu dosyayi ada gore bulup yukler (ispcr'in
#          aksine ZORUNLU degil, bulunamazsa okuma=None kalir) ve yukledikten
#          hemen sonra okuma uzunluk filtresini yapilandirmadaki duzeltilmis
#          degerlerle EZER: okuma.MINL, okuma.MAXL = C.NUMUNE_OKUMA_MIN,
#          C.NUMUNE_OKUMA_MAX. Ayrica ayni klasordeki numune.py ve
#          MADDE123_betikleri/numune_olc.py, kutu_cache.py "import okuma" der.
#          Menude motor.py uzerinden her olcum tusunda yuklenir: P, K, D, I, G,
#          T, U, H ve 1-9. KAPSAMLI_ARAMA.bat bu dosyayi dogrudan cagirmaz.
#
# DIKKAT - BU DOSYA PANELIN ESKI MOTORUDUR:
# Asagidaki Sonda sinifi TEK ve SABIT 13 bazlik bir tohum kullanir; bu tohum
# kayipli olabilir (ayrintili gerekce Sonda.__init__ basinda). Panelin bugunku
# OTORITE olcum yolu KAPSAMLI_ARAMA/okuma_motoru.py'dir (guvercin yuvasi
# tohumlamasi, kayipsiz). kendini_sina.py bu iki motoru yan yana kosar ve
# buradaki Sonda'nin ne kadar baglanma yeri kacirdigini olcup rapor eder.
# ---------------------------------------------------------------------------
import os, sys, glob, itertools, json, gzip
from collections import Counter
import ispcr

# Okuma uzunluk filtresi. Nanopore kutularinda cok kisa (parcalanmis) ve cok
# uzun (kimerik/birlesmis) okumalar bulunur; ikisi de urun orani hesabini
# bozar. motor.py bu iki degeri yukleme aninda yapilandirmadan gelen
# duzeltilmis degerlerle degistirir.
MINL, MAXL = 200, 6000
# Tohum uzunlugu. SABIT tek tohum - bkz. Sonda.__init__ icindeki uyari.
SEED = 13


# Dejenere (IUPAC) bir tohumu, olasi tum SOMUT dizilere acar. Ornegin iki R
# iceren bir tohum 4 ayri metne cevrilir. str.find yalniz somut metin arayabildigi
# icin gereklidir; dejenere baz sayisi arttikca liste ussel buyur.
def variants(p):
    return [''.join(x) for x in itertools.product(*[ispcr.IUPAC.get(c, 'ACGT') for c in p])]


# Tohum tuttuktan sonra TAM primer icin uyumsuzluk sayar. max_mm asilir asilmaz
# -1 ile erken cikar (kalan bazlari saymanin anlami yok).
# Not: ispcr.find_sites'taki gibi IUPAC farkindalidir, yani dejenere primer
# bazi izin verdigi her bazi uyumlu sayar.
def mm_ok(primer, win, max_mm):
    mm = 0
    for a, b in zip(primer, win):
        if b not in ispcr.IUPAC.get(a, 'ACGT'):
            mm += 1
            if mm > max_mm:
                return -1
    return mm


class Sonda:
    """Bir primerin belirli yonelimde okunmus dizide baglanma yerlerini bulur."""

    # -----------------------------------------------------------------------
    # TOHUMLAMA - BURASI KAYIPLI OLABILIR (guvercin yuvasi garantisi YOK)
    #
    # Bu sinif primerin 3' ucundaki (ters tumleyende 5' basindaki) TEK bir
    # SEED=13 bazlik parcayi alir ve yalnizca o parcanin TAM eslestigi yerlere
    # bakar. Sorun su: max_mm >= 1 iken uyumsuzluk pekala o 13 bazin ICINE
    # dusebilir. O zaman tohum hic tutmaz, aday konum hic uretilmez ve gercek
    # bir baglanma yeri SESSIZCE kacirilir - kod hata atmaz, sadece sayi kucuk
    # cikar.
    #
    # KAYIPSIZ OLMASI ICIN GEREKEN (guvercin yuvasi / pigeonhole):
    # En fazla m uyumsuzluga izin verilen bir aramada dizi (m + 1) tane
    # ORTUSMEYEN parcaya bolunurse, en fazla m parca bozulabilir; geriye kalan
    # EN AZ BIR parca tam eslesmek ZORUNDADIR. Yani "parcalardan herhangi biri
    # tam eslesiyor mu" diye aramak, kaba kuvvetin bulacagi hicbir yeri
    # kacirmaz. Bu bir GARANTIDIR, sezgisel hizlandirma degildir. Parca sayisi
    # m + 1'in ALTINA duserse garanti coker - buradaki gibi tek parca (yani
    # m = 0 varsayimi) kullanmak, max_mm = 1 ile kosuldugunda tam olarak bu
    # cokmus durumdur.
    #
    # Duzeltilmis surum: KAPSAMLI_ARAMA/okuma_motoru.py, primeri max_mm + 1
    # ortusmeyen bloga boler, her blok icin ayri tarama yapar ve bulunan her
    # adayi tam kural altinda yeniden dogrular. Panelin bugunku olcum yolu
    # odur; bu sinif tarihsel karsilastirma icin durur.
    # -----------------------------------------------------------------------
    def __init__(self, primer, uc5=False, max_mm=1):
        # uc5=True -> 3' kritik uc dizinin BASINDA (ters tumleyen primer)
        self.p = primer
        self.max_mm = max_mm
        self.uc5 = uc5
        # Tohum primerin 3' ucundan alinir: polimerazin uzattigi uc orasidir,
        # dolayisiyla gercek baglanmalarda en iyi korunan bolge de orasidir.
        # rc(R) icin 3' uc dizinin BASINA dustugu icin bas taraftan alinir.
        s = primer[:SEED] if uc5 else primer[-SEED:]
        # off: tohumun primer icindeki baslangic kaymasi. Tohum dizide i
        # konumunda bulununca primerin baslangici i - off olur.
        self.off = 0 if uc5 else len(primer) - SEED
        self.seeds = variants(s)

    # Tohumun gectigi HER yeri dener (find dongusu i+1'den devam eder, yani
    # ortusen tekrarlar da yakalanir), sonra her adayi tam primerle dogrular.
    # Donus: [(baslangic, uyumsuzluk)].
    def bul(self, seq):
        out = []
        L = len(self.p)
        for sd in self.seeds:
            i = seq.find(sd)
            while i != -1:
                st = i - self.off
                # Primer okumanin disina tasiyorsa aday gecersiz.
                if 0 <= st and st + L <= len(seq):
                    mm = mm_ok(self.p, seq[st:st + L], self.max_mm)
                    # mm >= 0 -> tam kural altinda da gecti.
                    if mm >= 0:
                        out.append((st, mm))
                i = seq.find(sd, i + 1)
        return out


# fastq okuyucu. fastq 4 satirlik kayitlardan olusur (@baslik / DIZI / + /
# kalite); k % 4 == 1 tam olarak DIZI satiridir. .gz dosyalar seffaf acilir,
# bozuk baytlar errors='ignore' ile atlanir (nanopore ciktilarinda olur).
# Filtre MINL..MAXL disindaki okumalari tamamen eler - bunlar paydaya da girmez.
def okumalar(path):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt', errors='ignore') as fh:
        for k, line in enumerate(fh):
            if k % 4 == 1:
                s = line.strip().upper()
                if MINL <= len(s) <= MAXL:
                    yield s


# ---------------------------------------------------------------------------
# kutu_pcr — bir kutunun (fastq) ham okumalarinda urun veren okuma orani.
#
# NE HESAPLAR: (toplam_okuma, urun_veren_okuma, urun_boyu_sayaci). Sonuc bir
# ORANDIR: pos / tot.
#
# NEDEN HAM OKUMALAR UZERINDE, REFERANS DIZI UZERINDE DEGIL:
# Referans/konsensus tek bir ozetlenmis dizidir; numunede fiilen bulunan
# varyantlari, tur ici farklari ve nanopore hata desenini icermez. Bir cift
# referansta kusursuz urun verip numunenin gercek okumalarinda tutmayabilir.
# Panelin karar sayisi bu yuzden ham okuma duzeyinde uretilir.
#
# YON - IKI YONUN DE DENENMESI SART:
# Nanopore okumalari CIFT YONLUDUR; ayni bolge kimi okumada arti, kimisinde
# eksi iplik olarak yazilir. Motorlar (ispcr.find_sites ve buradaki Sonda)
# yalniz ileri yonde tarar. Bu yuzden her okuma hem kendisi hem ters tumleyeni
# olarak denenir: "for seq in (s, ispcr.rc(s))". Bu dongu olmasaydi okumalarin
# kabaca yarisi sessizce kaybedilir, hata da atilmazdi.
#
# SAYIM KURALI: bir okuma urun verdigi anda "break" ile cikilir; okuma en fazla
# BIR kez sayilir. Yani bu bir VARLIK olcumudur (kac okumada urun var), derinlik
# olcumu degil. Ayni sekilde ilk uyan urun boyu kaydedilip ic donguler kirilir.
# ---------------------------------------------------------------------------
def kutu_pcr(path, F, R, lo=40, hi=600, max_mm=1):
    """Bir kutudaki okumalarda urun veren okuma sayisi ve toplam okuma."""
    # F ileri yonde aranir; R ise ters tumleyeni olarak, 3' ucu dizinin BASINA
    # dustugu icin uc5=True ile (ispcr.amplify'daki tail_pos=(0,1) ile ayni fikir).
    fs = Sonda(F, False, max_mm)
    rs = Sonda(ispcr.rc(R), True, max_mm)
    tot = pos = 0
    sizes = Counter()
    for s in okumalar(path):
        tot += 1
        # Okumanin iki yonu de denenir - bkz. yukaridaki YON aciklamasi.
        for seq in (s, ispcr.rc(s)):
            a = fs.bul(seq)
            # F yoksa bu yonde urun olamaz, R'yi aramaya gerek yok.
            if not a:
                continue
            b = rs.bul(seq)
            if not b:
                continue
            got = False
            for i, _ in a:
                for j, _ in b:
                    # Urun boyu: F'nin basindan R'nin (tumleyeninin) sonuna.
                    n = j + len(R) - i
                    # Boy penceresi + primerlerin ortusmeden karsilikli durmasi
                    # (ispcr.amplify ile birebir ayni olcut).
                    if lo <= n <= hi and j >= i + len(F):
                        sizes[n] += 1
                        got = True
                        break
                if got:
                    break
            if got:
                pos += 1
                # Okuma zaten sayildi; diger yonu denemeye gerek yok.
                break
    return tot, pos, sizes


if __name__ == '__main__':
    F, R = sys.argv[1], sys.argv[2]
    for p in sys.argv[3:]:
        t, n, sz = kutu_pcr(p, F, R)
        print(f'{os.path.basename(p):34s} {n:7d}/{t:<8d} {100*n/max(t,1):6.2f}%  {sz.most_common(3)}')
