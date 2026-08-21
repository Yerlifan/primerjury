# -*- coding: utf-8 -*-
"""Numunenin HAM OKUMALARINDA in-silico PCR + Wilson ayrim orani.

Olcut panelin numune olcutuyle birebir: uyumsuzluk <=1 ve 3' son 2 baz TAM.
Motor engine/scanner.py'nin Havuz sinifi + pair.py'nin urunler()
fonksiyonudur - bu dosya onlari ICE AKTARIR, yeniden yazmaz.

Havuz kutu basina BIR KEZ kurulur ve onbellege alinir; sonraki her aday cift
yalniz iki indeks sorgusu maliyetindedir. Binlerce adayi bu yuzden
tarayabiliyoruz.
"""
# ---------------------------------------------------------------------------
# sample.py — kutu basina ham okuma havuzlarini kurar ve bir primer cifti icin
#             uye/rakip kutulardaki urun oranlarini Wilson araligiyla olcer.
#
# GIRDI  : hedefler.kutular()'in verdigi kutu sozlukleri (kutu adi + fastq
#          yolu); config.py'deki NUMUNE_OKUMA_MIN/MAX, NUMUNE_OKUMA_SAYISI,
#          NUMUNE_TOHUM, NUMUNE_MAX_MM, KAPSAM_ESIGI, ENKOTU_ASGARI_OKUMA;
#          olcum motorlari engine_gateway.py uzerinden (okuma_motoru, tarayici.Havuz,
#          cift.urunler, ispcr.find_sites).
# CIKTI  : dosyaya yazmaz. Numune.olc() tek bir sozluk dondurur: uye ve rakip
#          kutu sayilari, uye yuzdeleri, kapsam payi, rakip havuz orani, ayrim
#          katlari (kat_havuz, kat_enkotu ve bunlarin yalniz kapsanan kutular
#          uzerinden hesaplanan esleri) ve urun boyu dagilimi.
# CAGRAN : __main__.hedefi_isle, panel_olcum.calistir ve
#          uyelik_denetimi.calistir icinden - yani full_chain.py asamalari
#          1, 2, 3, 4, 5, 7 ve 9. Disaridan protocol/single_protocol_measure.py
#          (tus P) ve verification/recovery_round.py (tus K) de ice aktarir.
# ---------------------------------------------------------------------------
import os, glob, gzip, math, random, pickle
from . import config as C
from . import engine_gateway


def okumalar(yol, n=C.NUMUNE_OKUMA_SAYISI, tohum=C.NUMUNE_TOHUM):
    """okuma_pcr.py ile ayni: 200-6000 bp filtre, sabit tohumla ornekleme."""
    ac = gzip.open if yol.endswith('.gz') else open
    hepsi = []
    with ac(yol, 'rt', errors='ignore') as fh:
        for i, satir in enumerate(fh):
            if i % 4 == 1:
                s = satir.strip().upper()
                if C.NUMUNE_OKUMA_MIN <= len(s) <= C.NUMUNE_OKUMA_MAX:
                    hepsi.append(s)
    if n and len(hepsi) > n:
        hepsi = random.Random(tohum).sample(hepsi, n)
    return hepsi


class KutuOtorite:
    """DUZELTILMIS motorla (read_engine.py) olcum - OTORITE yol.

    Guvercin yuvasi tohumlamasi kayipsizdir ve brute_force.py ile birebir
    dogrulanmistir. Hizli numpy yolundan yavastir; az sayida cift olculurken
    (panel yeniden olcumu) bu yol kullanilir.
    """

    def __init__(self, kutu, yol, n=C.NUMUNE_OKUMA_SAYISI):
        self.kutu = kutu
        rd, n0 = motor.okuma_motoru.kutu_yukle(
            yol, nmax=(n or 0), seed=C.NUMUNE_TOHUM,
            minl=C.NUMUNE_OKUMA_MIN, maxl=C.NUMUNE_OKUMA_MAX)
        self.reads = rd
        self.n_okuma = len(rd)
        self.suzgecten_gecen = n0

    def urun_veren(self, F, R, lo, hi, mm=C.NUMUNE_MAX_MM):
        pos, n, boy = motor.okuma_motoru.kutu_pcr(
            self.reads, F, R, lo=lo, hi=hi, max_mm=mm, son2=True)
        return pos, n, dict(boy)


class KutuEski:
    """PANELIN ESKI (HATALI) motoruyla olcum - YALNIZ KARSILASTIRMA ICIN.

    reads.py/Sonda'nin 13 bazlik TAM eslesen tohumu burada birebir yeniden
    uretilir. Bu yol hicbir karar icin kullanilmaz; amaci duzeltmenin hangi
    satirda ne kadar fark yarattigini GORUNUR kilmaktir.
    """

    def __init__(self, kutu, yol, n=C.NUMUNE_OKUMA_SAYISI):
        self.kutu = kutu
        rd, n0 = motor.okuma_motoru.kutu_yukle(
            yol, nmax=(n or 0), seed=C.NUMUNE_TOHUM,
            minl=C.NUMUNE_OKUMA_MIN, maxl=C.NUMUNE_OKUMA_MAX)
        self.reads = rd
        self.n_okuma = len(rd)

    @staticmethod
    def _sonda(primer, uc5, max_mm, seed=13):
        import itertools
        om = motor.okuma_motoru
        sd = primer[:seed] if uc5 else primer[-seed:]
        off = 0 if uc5 else len(primer) - seed
        tohumlar = [''.join(x) for x in itertools.product(
            *[om.IUPAC.get(c, 'ACGT') for c in sd])]
        L = len(primer)

        def bul(seq):
            out = []
            for t in tohumlar:
                i = seq.find(t)
                while i != -1:
                    st = i - off
                    if 0 <= st and st + L <= len(seq):
                        mm = 0; iyi = True
                        for a, b in zip(primer, seq[st:st + L]):
                            if b not in om.IUPAC.get(a, 'ACGT'):
                                mm += 1
                                if mm > max_mm:
                                    iyi = False; break
                        if iyi:
                            out.append((st, mm))
                    i = seq.find(t, i + 1)
            return out
        return bul

    def urun_veren(self, F, R, lo, hi, mm=C.NUMUNE_MAX_MM):
        om = motor.okuma_motoru
        fs = self._sonda(F, False, mm)
        rs = self._sonda(om.rc(R), True, mm)
        pos = 0
        boy = {}
        for s in self.reads:
            vur = None
            for seq in (s, om.rc(s)):
                a = fs(seq)
                if not a:
                    continue
                b = rs(seq)
                if not b:
                    continue
                for i, _ in a:
                    for j, _ in b:
                        n = j + len(R) - i
                        if lo <= n <= hi and j >= i + len(F):
                            vur = n; break
                    if vur:
                        break
                if vur:
                    break
            if vur:
                pos += 1
                boy[vur] = boy.get(vur, 0) + 1
        return pos, self.n_okuma, boy


class KutuHavuzu:
    """Bir kutunun okumalari + ters tumleyenleri, tarayici.Havuz indeksiyle."""

    def __init__(self, kutu, yol, n=C.NUMUNE_OKUMA_SAYISI):
        self.kutu = kutu
        rd = okumalar(yol, n)
        self.n_okuma = len(rd)
        # her okuma iki yonde denenir -> havuza ikisini de koy, sonra esle
        diziler = []
        self.okuma_id = []
        for i, s in enumerate(rd):
            s = motor.clean(s)
            diziler.append(s); self.okuma_id.append(i)
            diziler.append(motor.rc(s)); self.okuma_id.append(i)
        self.hv = motor.tarayici.Havuz(diziler) if diziler else None

    def urun_veren_kaba(self, F, R, lo, hi, mm):
        """TOHUMSUZ (kaba kuvvet) tarama - her uyumsuzluk tavani icin DOGRU.

        tarayici.Havuz'un tohumu yalniz TEK uyumsuzluk icin tamdir (tohum
        varyantlari bir ikame ile uretilir). <=3 gibi daha gevsek olcutlerde
        tohumlu arama site KACIRIR; bu yol ispcr.find_sites'i havuzun
        birlestirilmis dizisi uzerinde dogrudan kullanir - kacirma yoktur.
        """
        import numpy as np
        if self.hv is None or self.n_okuma == 0:
            return 0, 0, {}
        enc, sid = self.hv.enc, self.hv.sid
        revrc = motor.rc(R)
        fs = motor.find_sites(enc, F, mm, need_tail=True, tail_pos=(-1, -2))
        if not fs:
            return 0, self.n_okuma, {}
        rs = motor.find_sites(enc, revrc, mm, need_tail=True, tail_pos=(0, 1))
        if not rs:
            return 0, self.n_okuma, {}
        fpos = np.array([x[0] for x in fs]); rpos = np.array([x[0] for x in rs])
        fid = sid[fpos]; rid = sid[rpos]
        ok = fid >= 0; fpos, fid = fpos[ok], fid[ok]
        ok = rid >= 0; rpos, rid = rpos[ok], rid[ok]
        from collections import defaultdict
        rmap = defaultdict(list)
        for pz, i in zip(rpos.tolist(), rid.tolist()):
            rmap[i].append(pz)
        veren = set(); boy = {}
        for pz, i in zip(fpos.tolist(), fid.tolist()):
            for q in rmap.get(i, ()):
                n = q + len(revrc) - pz
                if lo <= n <= hi and q >= pz + len(F):
                    veren.add(self.okuma_id[i])
                    boy[n] = boy.get(n, 0) + 1
                    break
        return len(veren), self.n_okuma, boy

    def urun_veren(self, F, R, lo, hi, mm=C.NUMUNE_MAX_MM):
        """Urun veren OKUMA sayisi (iki yon tek okuma sayilir) ve boy dagilimi."""
        if self.hv is None or self.n_okuma == 0:
            return 0, 0, {}
        if mm > 1:
            return self.urun_veren_kaba(F, R, lo, hi, mm)
        m, boy = motor.cift.urunler(self.hv, F, R, lo=lo, hi=hi, mm=mm)
        veren = set()
        for idx in m.nonzero()[0]:
            veren.add(self.okuma_id[idx])
        return len(veren), self.n_okuma, boy


class Numune:
    """otorite=True -> read_engine.py (kayipsiz, yavas; panel yeniden olcumu)
       otorite=False -> numpy havuz (hizli; binlerce aday taramasi)"""

    def __init__(self, kutular, n=C.NUMUNE_OKUMA_SAYISI, ilerle=None, otorite=False):
        self.havuz = {}
        self.otorite = otorite
        sinif = KutuOtorite if otorite else KutuHavuzu
        for i, k in enumerate(kutular):
            if ilerle:
                ilerle(i + 1, len(kutular), k['kutu'])
            self.havuz[k['kutu']] = sinif(k['kutu'], k['yol'], n)

    def olc(self, F, R, uye_kutu, rakip_kutu, lo, hi, mm=C.NUMUNE_MAX_MM):
        # AYRIM KATI NASIL KURULUR VE NEDEN BOYLE
        # Pay: uye kutulari icinde EN KOTU olanin Wilson ALT siniri. Ortalama
        # degil en kotu, cunku panelin vaadi "butun uye kutularinda cogaltir";
        # tek bir uye kutusu bos kaliyorsa cift o kutu icin calismiyordur.
        # Payda: rakip tarafinin Wilson UST siniri. Iki tarafta da muhafazakar
        # sinir secildigi icin oran asla oldugundan buyuk cikmaz.
        # Wilson araligi kullanilir cunku ham yuzde, dusuk okuma sayisinda
        # belirsizligi gizler; aralik belirsizligi sayinin icine katar. Bunun
        # yan etkisi sudur: ayni gercek ozgulluk sig bir havuzda DAHA DUSUK bir
        # 'x' verir, bu yuzden farkli derinliklerde olculmus satirlar birbiriyle
        # karsilastirilamaz.
        uy, rk = [], []
        boylar = {}
        for k in uye_kutu:
            h = self.havuz.get(k['kutu'])
            if h is None:
                continue
            p, n, boy = h.urun_veren(F, R, lo, hi, mm)
            uy.append((k['kutu'], p, n))
            for s, c in boy.items():
                boylar[s] = boylar.get(s, 0) + c
        for k in rakip_kutu:
            h = self.havuz.get(k['kutu'])
            if h is None:
                continue
            p, n, _ = h.urun_veren(F, R, lo, hi, mm)
            rk.append((k['kutu'], p, n))
        if not uy:
            return None
        uye_alt = min(motor.wilson(p, n)[0] for _, p, n in uy)
        # KAPSAM ekseni: panel bazi hedefleri "13/13 kutu", "33/34 kutu" diye bildirir.
        # Tek bir uye kutusu urun vermeyince uye_alt 0'a duser ve butun adaylar
        # ayirt edilemez hale gelir; kapsam bu yuzden AYRI olculur.
        kapsayan = [(a, p, n) for a, p, n in uy if n and p / n >= C.KAPSAM_ESIGI]
        uye_alt_k = (min(motor.wilson(p, n)[0] for _, p, n in kapsayan)
                     if kapsayan else 0.0)
        rp = sum(p for _, p, _ in rk)
        rn = sum(n for _, _, n in rk)
        # K-4 DUZELTMESI (2026-08-03): eskiden rakip yokken payda 1e-9 idi;
        # bolme bir milyarla carpmaya donuyordu ve evrensel primerler
        # 27 milyon 'x' ile ESIK USTU cikiyordu. Artik None -> OLCULEMEDI.
        #
        # ISIN OZU: evrensel bir hedefte uyelik tanimi geregi rakip kumesi bosa
        # yaklasir. Payda sifira giderken oran tanimsizlasir - nitekim ayni
        # sutunda 0,00 ile 117 milyon yan yana durabiliyordu ve ikisi de olcum
        # degildi. Kucuk bir sabit koyup bolmeye devam etmek bu tanimsizligi
        # gizler, cozmez. Dogrusu orani HIC URETMEMEK (None = olculemedi) ve o
        # hedefleri kapsama ile alan disi orani uzerinden degerlendirmektir.
        # Bu esigi DUSURMEK DEGILDIR; 10x esigi yerinde durur, yalnizca esigin
        # uygulanamayacagi satirlarda baska bir buyukluk olculur.
        havuz_ust = motor.wilson(rp, rn)[1] if rn else None
        enkotu = None
        # "en kotu tek rakip kutu" icin anlamli payda esigi.
        #
        # DUZELTME (2026-08-02): esik eskiden "en buyuk kutunun yarisi" idi.
        # Bu, DERINLIGE gore kayiyordu: tam derinlikte en buyuk kutu 46 472
        # okuma olunca esik 23 236 oluyor ve 10 rakip kutudan YALNIZ 1'i olcume
        # giriyordu - yani "en kotu kutu" aslinda "en derin kutu" demek oluyordu
        # ve gercek en kotu rakip disarida kalabiliyordu. Ayni cift 300 okumayla
        # olculunce esik 150 olup 18 kutunun hepsi giriyordu; iki asama arasinda
        # 40 kata varan sahte fark bu yuzden cikiyordu.
        # Artik MUTLAK: derinlik degisse de ayni kutular olcume girer.
        esik = C.ENKOTU_ASGARI_OKUMA
        for kadi, p, n in rk:
            if n < esik:
                continue
            hi_ = motor.wilson(p, n)[1]
            if enkotu is None or hi_ > enkotu[1]:
                enkotu = (kadi, hi_, p, n)
        return dict(
            olcut='<=%d uyumsuzluk + 3\' son 2 baz TAM' % mm,
            max_mm=mm,
            uye=[(a, p, n, round(100.0 * p / max(n, 1), 2)) for a, p, n in uy],
            rakip=sorted([(a, p, n, round(100.0 * p / max(n, 1), 2)) for a, p, n in rk],
                         key=lambda x: -x[3]),
            uye_alt=round(100 * uye_alt, 3),
            uye_min=round(100 * min(p / max(n, 1) for _, p, n in uy), 2),
            uye_max=round(100 * max(p / max(n, 1) for _, p, n in uy), 2),
            uye_kutu_sayisi=len(uy),
            uye_kapsam=len(kapsayan),
            uye_kapsam_pay='%d/%d' % (len(kapsayan), len(uy)),
            uye_alt_kapsayan=round(100 * uye_alt_k, 3),
            havuz='%d/%d' % (rp, rn),
            # HATA DUZELTMESI (2026-08-06, temiz kosuda yakalandi): havuz_ust
            # rn == 0 iken bilerek None'dur ("olculemedi" - yukaridaki gerekce).
            # Asagidaki kat_havuz / kat_havuz_kapsayan satirlari bunu dogru
            # sekilde koruyordu, BU SATIR KORUMUYORDU ve rakip okumasi sifir olan
            # ILK hedefte cokuyordu:
            #   TypeError: unsupported operand type(s) for *: 'int' and 'NoneType'
            # P asamasi 22 ciftin 5.'sinde bu yuzden dustu, T asamasi cikis kodu 3
            # verdi ve K, D, I alt asamalari hic kosmadi. None korunur: alan
            # "olculemedi" olarak raporlanir, uydurma bir sayi URETILMEZ.
            havuz_ust=(round(100 * havuz_ust, 3) if havuz_ust is not None else None),
            # O-8: kat_havuz ancak yeterli derinlikte verdikt tasiyabilir
            kat_havuz=(round(uye_alt / havuz_ust, 2)
                       if (havuz_ust and rn >= C.ENKOTU_ASGARI_OKUMA) else None),
            havuz_derinligi=rn,
            rakip_olculen=sum(1 for _, _, n in rk if n >= C.ENKOTU_ASGARI_OKUMA),
            rakip_toplam=len(rk),
            enkotu_kutu=enkotu[0] if enkotu else '',
            kat_enkotu=round(uye_alt / enkotu[1], 2) if enkotu and enkotu[1] > 0 else None,
            kat_havuz_kapsayan=(round(uye_alt_k / havuz_ust, 2)
                                if (havuz_ust and rn >= C.ENKOTU_ASGARI_OKUMA
                                    and kapsayan) else None),
            kat_enkotu_kapsayan=(round(uye_alt_k / enkotu[1], 2)
                                 if enkotu and enkotu[1] > 0 and kapsayan else None),
            urun_boylari=sorted(boylar.items(), key=lambda x: -x[1])[:5],
        )
