# -*- coding: utf-8 -*-
r"""
read_engine.py - DUZELTILMIS ham-okuma in-silico PCR motoru.
2026-08-02, okuma motoru duzeltmesi.

BAGIMSIZDIR: numpy, ispcr, reads.py, pysam, minimap2 GEREKMEZ. Saf standart Python.
Windows/Linux farketmez. .bat icinden dogrudan cagrilabilir.

---------------------------------------------------------------------------
NEYI DUZELTIYOR
---------------------------------------------------------------------------
Eski motor (engine/reads.py -> sinif `Sonda`, ayrica engine/
scb.py -> sinif `S`) primeri okumada ararken 13 bazlik TAM ESLESEN tek bir tohum
kullaniyordu:

    s = primer[:13] if uc5 else primer[-13:]
    i = seq.find(sd)          # <-- TAM eslesme sarti

Olcut "toplam uyumsuzluk <= max_mm" oldugu halde, uyumsuzluk 13 bazlik tohumun
icine dustugunde `find` hicbir sey bulamiyor ve baglanma yeri SESSIZCE kayboluyor.
Program hata vermiyor, "urun yok" diyor.

Olculen etki (A1-4_3078083 kutusu, ilk 400 okuma, M. mazei cifti, AYNI olcut mm<=1):
    eski Sonda        :   2/400  (%0,50)
    kaba kuvvet       : 174/400  (%43,50)
    -> 205 baglanma yerinin 202'sinde (%98,5) uyumsuzluk tohumun icine dusuyor
       (188'i primerin 6. bazinda - tek ve tekrar eden bir varyant baz).

---------------------------------------------------------------------------
NASIL DUZELTILDI - GUVERCIN YUVASI (PIGEONHOLE) TOHUMLAMASI
---------------------------------------------------------------------------
Primer, max_mm+1 tane ORTUSMEYEN bloga bolunur. En fazla max_mm uyumsuzluk varsa
guvercin yuvasi ilkesi geregi bu bloklardan EN AZ BIRI tam tutmak zorundadir.
Dolayisiyla bloklardan herhangi birinin tam eslesmesini aramak KAYIPSIZDIR:
kaba kuvvetin bulacagi hicbir yeri kacirmaz. Bulunan her aday yer sonra tam
kuralla (toplam uyumsuzluk + 3' son 2 baz) TEK TEK dogrulanir.

Ayrica: eski kodda 3' son 2 baz TAM ESLESME kurali hicbir yerde acikca
uygulanmiyordu - 13 bazlik tohumun yan etkisiydi. Tohum kisalinca o yan etki
kaybolacagi icin kural bu modulde ACIKCA uygulanir (son2=True).

KAYIPSIZLIK IDDIASI TEST EDILMISTIR: engine_test.py bu motoru, tohum kullanmayan
ve her pozisyonu tek tek deneyen bagimsiz bir kaba kuvvet uygulamasiyla
(brute_force.py) karsilastirir. Panele girmeden once o testi kosun.

---------------------------------------------------------------------------
KULLANIM
---------------------------------------------------------------------------
Modul olarak:
    import read_engine as om
    okumalar = list(om.okumalar('kutu.fastq'))            # 200-6000 bp suzgeci
    pos, n   = om.kutu_pcr(okumalar, F, R, max_mm=1)      # urun veren okuma sayisi
    yerler   = om.Sonda(F, uc5=False, max_mm=1).bul(dizi) # [(baslangic, uyumsuzluk), ...]

Komut satirindan (.bat icin):
    python read_engine.py F_PRIMER R_PRIMER [--mm 1] [--lo 40] [--hi 600]
                           [--nmax 3000] [--seed 3] [--tsv cikti.tsv] dosya1.fastq ...
    -> her fastq icin satir: dosya, urun_veren, kullanilan_okuma, toplam_okuma, yuzde, baskin_boylar

    Ornek:
      python read_engine.py GCCCTTGGGACCGGCATAA TCGCTGGCTAGTAGGTACATTACA ^
             --mm 1 --tsv sonuc.tsv "fastq files\A1-4\*.fastq"

OLCUT ETIKETI: --mm degeri ciktiya yazilir. Panelin numune olcutu mm<=1'dir;
mm<=3 tasarim boru hattinin olcutudur. IKISINI KARISTIRMAYIN - panelin eski
satirlarinda karistirilmisti (bkz. "16 Okuma Motoru Duzeltmesi" sayfasi).
"""
# ---------------------------------------------------------------------------
# read_engine.py — ham nanopore okumalari uzerinde in-silico PCR yapan
#                   duzeltilmis motor; guvercin yuvasi tohumlamasi kullanir ve
#                   hicbir baglanma yerini kacirmaz.
#
# GIRDI  : fastq / fastq.gz dosyalari (okumalar() ve kutu_yukle() ile,
#          200-6000 bp uzunluk suzgeci ve sabit tohumlu ornekleme) ile ileri ve
#          geri primer dizileri. Modul olarak da, komut satirindan da kullanilir.
# CIKTI  : dosyaya yazmaz (komut satirinda --tsv verilirse dosya basina bir
#          satirlik TSV yazar). Sonda.bul() [(baslangic, uyumsuzluk)] listesi,
#          kutu_pcr() (urun_veren, toplam, boy_sayaci) uclusu dondurur.
# CAGRAN : Paket ici: engine_gateway.py bu dosyayi ice aktarir, numune.KutuOtorite ve
#          numune.KutuEski uzerinden butun olcum asamalarinda kullanilir - yani
#          full_chain.py asamalari 1, 2, 3, 4, 5, 6, 7, 8, 9. Ayrica orientation.py,
#          orientation_audit.py, orientation_impact_test.py, cross_coverage.py, engine_test.py
#          ve independent_check.py dogrudan ice aktarir.
# ---------------------------------------------------------------------------
import sys, os, re, glob, gzip, random, argparse
from collections import Counter

__version__ = '1.0 (2026-08-02)'

MINL, MAXL = 200, 6000          # okuma uzunluk suzgeci (A2 ~4,2-4,5 kb, F2 ~3,7 kb dahil)

IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT',
}
_COMP = str.maketrans('ACGTURYSWKMBDHVNacgturyswkmbdhvn',
                      'TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn')


def rc(s):
    """ters tumleyen"""
    return s.translate(_COMP)[::-1]


def temizle(s):
    return re.sub(r'[^ACGTUN]', 'N', s.upper()).replace('U', 'T')


def _blok_regex(alt):
    """bir primer blogunu, IUPAC dejenere bazlari karakter sinifina cevirerek
    ORTUSEN eslesme bulabilen bir regex'e cevirir."""
    govde = ''.join(c if c in 'ACGT' else '[' + IUPAC.get(c, 'ACGT') + ']' for c in alt)
    return re.compile('(?=(' + govde + '))')


class Sonda:
    """Bir primerin bir dizideki TUM baglanma yerlerini bulur. Kayipsiz.

    uc5=False -> ileri primer; 3' kritik uc primerin SONUNDA (indeks -1, -2)
    uc5=True  -> ters primerin tumleyeni; 3' kritik uc dizinin BASINDA (indeks 0, 1)
    son2=True -> 3' son 2 baz TAM eslesmek zorunda (panelin kurali)
    """

    def __init__(self, primer, uc5=False, max_mm=1, son2=True):
        primer = primer.upper()
        self.p = primer
        self.L = len(primer)
        self.max_mm = max_mm
        self.son2 = son2
        self.uc5 = uc5
        self.ok = [set(IUPAC.get(c, 'ACGT')) for c in primer]
        self.kritik = (0, 1) if uc5 else (self.L - 1, self.L - 2)
        # ------------------------------------------------------------------
        # GUVERCIN YUVASI (PIGEONHOLE) TOHUMLAMASI - NEDEN KAYIPSIZ
        #
        # k uzunlugundaki bir primer, en fazla m uyumsuzluga izin verilen bir
        # aramada (m + 1) tane ORTUSMEYEN parcaya bolunurse, en fazla m parca
        # bozulabilir; geriye kalan EN AZ BIR parca tam eslesmek ZORUNDADIR.
        # Bu bir garantidir, sezgisel bir hizlandirma degildir: parcalardan
        # herhangi birinin tam eslesmesini aramak, kaba kuvvetin bulacagi
        # hicbir yeri kacirmaz. Bulunan her aday konum sonra _dogrula() ile
        # tam kural altinda (toplam uyumsuzluk + 3' son 2 baz) tek tek
        # sinanir, yani yanlis pozitif de birakmaz.
        #
        # BU YUZDEN TOHUM UZUNLUGU VE PARCA SAYISI KEYFI SECILEMEZ. Parca
        # sayisi max_mm + 1'den KUCUK olursa garanti coker ve arama sessizce
        # site kacirmaya baslar - eski motorun hatasi tam olarak buydu: tek,
        # sabit 13 bazlik bir tohum kullaniyordu ve uyumsuzluk o 13 bazin
        # icine dustugunde `find` hicbir sey bulamiyordu. Program hata vermez,
        # "urun yok" der. Olculen etki: 205 baglanma yerinin 202'si (%98,5)
        # kayboluyordu. Parca sayisini artirmak da bedava degildir - parcalar
        # kisaldikca yanlis aday sayisi ve dogrulama maliyeti buyur.
        # Kayipsizlik iddiasi engine_test.py'de brute_force.py'ye karsi
        # dogrulanir; self_test.py her kosuda ayni karsilastirmayi yapar.
        # ------------------------------------------------------------------
        nb = max_mm + 1                                   # guvercin yuvasi blok sayisi
        kes = [round(i * self.L / nb) for i in range(nb + 1)]
        self.bloklar = [(kes[i], _blok_regex(primer[kes[i]:kes[i + 1]]))
                        for i in range(nb) if kes[i + 1] > kes[i]]

    def _dogrula(self, seq, st):
        """aday yeri tam kuralla dogrula -> uyumsuzluk sayisi ya da None"""
        # 3' son 2 baz TAM ESLESME kurali BURADA ACIKCA uygulanir. Eski kodda bu
        # kural hicbir yerde yaziyla yer almiyordu; 13 bazlik uzun tohumun yan
        # etkisi olarak kendiliginden saglaniyordu. Tohum kisalinca o yan etki
        # kaybolacagi icin kural acik hale getirildi. Biyolojik gerekcesi:
        # polimeraz uzatmayi 3' uctan baslatir, o uctaki uyumsuzluk urunu fiilen
        # engeller - dolayisiyla toplam uyumsuzluk sayisindan bagimsiz bir sart.
        if st < 0 or st + self.L > len(seq):
            return None
        if self.son2:
            for k in self.kritik:
                if seq[st + k] not in self.ok[k]:
                    return None
        mm = 0
        for k in range(self.L):
            if seq[st + k] not in self.ok[k]:
                mm += 1
                if mm > self.max_mm:
                    return None
        return mm

    def bul(self, seq):
        aday = set()
        for off, rx in self.bloklar:
            for m in rx.finditer(seq):
                aday.add(m.start() - off)
        out = []
        for st in sorted(aday):
            mm = self._dogrula(seq, st)
            if mm is not None:
                out.append((st, mm))
        return out


def urun_var(seq, fs, rs, lenF, lenR, lo, hi):
    """tek yonde: F ve R karsilikli baglanip verilen boy penceresinde urun veriyor mu"""
    a = fs.bul(seq)
    if not a:
        return None
    b = rs.bul(seq)
    if not b:
        return None
    for i, _ in a:
        for j, _ in b:
            n = j + lenR - i
            if lo <= n <= hi and j >= i + lenF:
                return n
    return None


def kutu_pcr(okuma_listesi, F, R, lo=40, hi=600, max_mm=1, son2=True):
    """Bir kutudaki okumalarda urun veren okuma sayisi.
    Donus: (urun_veren, toplam, boy_sayaci)"""
    F = F.upper(); R = R.upper()
    fs = Sonda(F, False, max_mm, son2)
    rs = Sonda(rc(R), True, max_mm, son2)
    pos = 0
    boylar = Counter()
    for s in okuma_listesi:
        for seq in (s, rc(s)):
            n = urun_var(seq, fs, rs, len(F), len(R), lo, hi)
            if n is not None:
                boylar[n] += 1
                pos += 1
                break
    return pos, len(okuma_listesi), boylar


def okumalar(path, minl=MINL, maxl=MAXL):
    """fastq / fastq.gz akisi, uzunluk suzgeciyle"""
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt', errors='ignore') as fh:
        for k, line in enumerate(fh):
            if k % 4 == 1:
                s = line.strip().upper()
                if minl <= len(s) <= maxl:
                    yield s


def kutu_yukle(path, nmax=3000, seed=3, minl=MINL, maxl=MAXL):
    """Bir fastq'tan okumalari al, nmax'i asiyorsa sabit tohumla ornekle.
    Donus: (ornek_okumalar, suzgecten_gecen_toplam)"""
    rs = list(okumalar(path, minl, maxl))
    n0 = len(rs)
    if nmax and len(rs) > nmax:
        random.seed(seed)
        rs = random.sample(rs, nmax)
    return rs, n0


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='Duzeltilmis ham-okuma in-silico PCR (kayipsiz tohumlama).')
    ap.add_argument('F'); ap.add_argument('R')
    ap.add_argument('fastq', nargs='+', help='fastq files or a glob pattern')
    ap.add_argument('--mm', type=int, default=1, help='total mismatches allowed (default 1)')
    ap.add_argument('--lo', type=int, default=40)
    ap.add_argument('--hi', type=int, default=600)
    ap.add_argument('--nmax', type=int, default=3000, help='maximum reads per bin (0 = all)')
    ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--minlen', type=int, default=MINL)
    ap.add_argument('--maxlen', type=int, default=MAXL)
    ap.add_argument('--last-two', '--son2', dest='son2', type=int, default=1, help="1 = require an exact match at the last two 3' bases")
    ap.add_argument('--tsv', default='', help='write the result to this file as TSV')
    a = ap.parse_args(argv)

    yollar = []
    for d in a.fastq:
        g = sorted(glob.glob(d))
        yollar.extend(g if g else [d])

    satirlar = []
    for p in yollar:
        if not os.path.exists(p):
            sys.stderr.write('YOK: %s\n' % p); continue
        rs, n0 = kutu_yukle(p, a.nmax, a.seed, a.minlen, a.maxlen)
        pos, n, boylar = kutu_pcr(rs, a.F, a.R, a.lo, a.hi, a.mm, bool(a.son2))
        pct = 100.0 * pos / n if n else 0.0
        bb = ';'.join('%d:%d' % x for x in boylar.most_common(3))
        satirlar.append((os.path.basename(p), pos, n, n0, round(pct, 2), a.mm, bb))
        print('%-40s %6d/%-6d %6.2f%%  (mm<=%d, suzgecten gecen %d)  %s'
              % (os.path.basename(p), pos, n, pct, a.mm, n0, bb), flush=True)

    if a.tsv:
        with open(a.tsv, 'w', encoding='utf-8') as fh:
            fh.write('dosya\turun_veren\tkullanilan_okuma\tsuzgecten_gecen\tyuzde\tolcut_mm\tbaskin_boylar\n')
            for r in satirlar:
                fh.write('\t'.join(str(x) for x in r) + '\n')
        sys.stderr.write('TSV yazildi: %s\n' % a.tsv)
    return 0


if __name__ == '__main__':
    sys.exit(main())
