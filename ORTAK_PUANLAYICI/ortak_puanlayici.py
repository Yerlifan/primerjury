# -*- coding: utf-8 -*-
"""ORTAK PUANLAYICI - panelin KENDI olcusuyle puanlar, hicbir sayiyi yeniden yazmaz.

NEDEN VAR
---------
2026-08-10: ayni primer cifti icin farkli betikler farkli ayrim sayisi uretti.
Sebep aritmetik degil, UYELIK KAYNAGI idi:

  tek_protokol_olc.py  -> uyelik_yeniden_turetme_uyelik_*.tsv (OLCULEN uyelik),
                          karisik kutu RAKIP sayilir, kutu adlari normalize edilir.
  cok_lokus.py,        -> hedefler.hedef_baglami() : WSL_betikleri/hedefler.tsv
  referans_tasarim.py     ve KAPSAMLI_ARAMA/hedef_uyelik.tsv uzerinden TAXID
                          kumesi + sinif suzgeci. Karisik kutu kavrami YOK.

Iki yol farkli uye/rakip kumesi kurar. Her ikisi de kendi icinde dogru hesap
yapar ama sonuclar KARSILASTIRILAMAZ.

BU MODUL panelin yolunu TEK dogru yol kabul eder ve onu ICE AKTARIR:
  - uyelik okuma      : tek_protokol_olc.uyelik_oku
  - kutu adi duzeltme : tek_protokol_olc.kutu_adi_normalize
  - verdikt           : tek_protokol_olc.karar
  - esik/protokol     : tek_protokol_olc.PROTOKOL
  - ASIL OLCUM        : KAPSAMLI_ARAMA.numune.Numune.olc  (yeniden yazilmadi)

Bu dosyada HICBIR ayrim sayisi hesaplanmaz; yalnizca panelin hesapladigi sayi
tasinir.

ONBELLEK ANAHTARI: anahtara primer DIZISI girer (F|R). Dizi degisince onbellek
kendiliginden gecersiz olur.
"""
import os, sys, json, hashlib


def _kur(kok):
    kok = os.path.abspath(kok)
    if not os.path.isdir(os.path.join(kok, 'KAPSAMLI_ARAMA')):
        raise SystemExit('HATA: %s icinde KAPSAMLI_ARAMA yok.' % kok)
    if kok not in sys.path:
        sys.path.insert(0, kok)
    tp_dir = os.path.join(kok, 'TEK_PROTOKOL')
    if tp_dir not in sys.path:
        sys.path.insert(0, tp_dir)
    import importlib.util
    ad = 'tek_protokol_olc'
    if ad in sys.modules:
        TP = sys.modules[ad]
    else:
        spec = importlib.util.spec_from_file_location(
            ad, os.path.join(tp_dir, 'tek_protokol_olc.py'))
        TP = importlib.util.module_from_spec(spec)
        sys.modules[ad] = TP
        spec.loader.exec_module(TP)
    from KAPSAMLI_ARAMA import numune as N, hedefler as H, yapilandirma as C, motor
    return TP, N, H, C, motor


class Puanlayici(object):
    """Panel derinliginde puanlayici.

    otorite=True  -> okuma_motoru (KutuOtorite). PANELIN KULLANDIGI YOL.
    otorite=False -> numpy havuz (KutuHavuzu). YALNIZ ON ELEME.
    Ikisi ayni sayiyi vermek ZORUNDA DEGILDIR; panel derinligindeki her karar
    otorite=True ile uretilir.
    """

    def __init__(self, kok, derinlik=None, karisik_kural=None, otorite=True,
                 onbellek_yolu=None):
        self.kok = os.path.abspath(kok)
        self.TP, self.N, self.H, self.C, self.motor = _kur(self.kok)
        P = self.TP.PROTOKOL
        self.derinlik = P['okuma_tavani'] if derinlik is None else derinlik
        self.karisik_kural = P['karisik'] if karisik_kural is None else karisik_kural
        self.otorite = otorite
        self.lo, self.hi = P['urun_alt'], P['urun_ust']
        self.mm_asil, self.mm_yan = P['olcut_asil'], P['olcut_yan']
        self.esik = P['esik']

        self.uy_yol = self.TP.uyelik_dosyasi(self.kok)
        if not self.uy_yol:
            raise SystemExit('HATA: uyelik tablosu bulunamadi.')
        self.uyelik = self.TP.uyelik_oku(self.uy_yol)
        # EK CIFTLER: panelde olmayan ciftler uyeligi kendi adiyla degil
        # 'uyelik_hedefi' sutunuyla cozer (tek_protokol_olc.calistir ile ayni).
        # Bu satir olmadan Petriella_cinsi uyeliksiz kaliyor ve OLCULEMEDI
        # dusuyordu - kalibrasyonda yakalandi.
        self.ek_uyelik = {}
        for e in self.TP.ek_ciftler_oku(self.kok):
            if e.get('uyelik_hedefi'):
                self.ek_uyelik[e['hedef']] = e['uyelik_hedefi']
        self.kut = {k['kutu']: k for k in self.H.kutular()}
        self._baglam = {}
        self._nm = None
        self._yuklu = set()
        self.eksik_kutu = set()

        self.onbellek_yolu = onbellek_yolu
        self._ob = {}
        if onbellek_yolu and os.path.exists(onbellek_yolu):
            try:
                self._ob = json.load(open(onbellek_yolu, encoding='utf-8'))
            except Exception:
                self._ob = {}

    def _coz(self, adlar):
        """tek_protokol_olc.calistir icindeki coz() ile BIREBIR ayni."""
        out = []
        for a in adlar:
            a2 = self.TP.kutu_adi_normalize(a.strip())
            if a2 in self.kut:
                out.append(self.kut[a2])
            elif a.strip() in self.kut:
                out.append(self.kut[a.strip()])
            else:
                self.eksik_kutu.add(a.strip())
        return out

    def baglam(self, hedef, sinif_yedek=''):
        """Panelin uye/rakip kumesi. tek_protokol_olc.calistir ile BIREBIR ayni."""
        if hedef in self._baglam:
            return self._baglam[hedef]
        u = self.uyelik.get(hedef)
        if u is None:
            u = self.uyelik.get(self.ek_uyelik.get(hedef, ''), None)
        if u is None:
            self._baglam[hedef] = None
            return None
        uye = self._coz(u['uye'])
        kar = self._coz(u['karisik'])
        rak = self._coz(u['rakip'])
        if self.karisik_kural == 'uye':
            uye = uye + kar
        elif self.karisik_kural == 'rakip':
            rak = rak + kar
        if not rak:
            uye_ad = {k['kutu'] for k in uye} | {k['kutu'] for k in kar}
            rak = [k for k in self.kut.values()
                   if k['sinif'] == (u['sinif'] or sinif_yedek) and k['kutu'] not in uye_ad]
        b = dict(uye=uye, rakip=rak, karisik=kar)
        self._baglam[hedef] = b
        return b

    def havuz_hazirla(self, hedefler, ilerle=None):
        ger = {}
        for h in hedefler:
            b = self.baglam(h)
            if not b:
                continue
            for k in b['uye'] + b['rakip']:
                ger[k['kutu']] = k
        yeni = [k for a, k in sorted(ger.items()) if a not in self._yuklu]
        if not yeni:
            return 0
        if self._nm is None:
            self._nm = self.N.Numune([], n=self.derinlik, otorite=self.otorite)
        sinif = self.N.KutuOtorite if self.otorite else self.N.KutuHavuzu
        for i, k in enumerate(yeni):
            if ilerle:
                ilerle(i + 1, len(yeni), k['kutu'])
            self._nm.havuz[k['kutu']] = sinif(k['kutu'], k['yol'], self.derinlik)
            self._yuklu.add(k['kutu'])
        return len(yeni)

    def anahtar(self, hedef, F, R, mm):
        """ONBELLEK ANAHTARI - primer DIZISI dahildir."""
        ham = u'%s|%s|%s|%d|%d|%s|%s|%d|%d|%s' % (
            hedef, F.upper(), R.upper(), mm, self.derinlik, self.karisik_kural,
            os.path.basename(self.uy_yol), self.lo, self.hi,
            'otorite' if self.otorite else 'havuz')
        return hashlib.md5(ham.encode('utf-8')).hexdigest()

    def olc_ham(self, hedef, F, R, mm=None):
        """Panelin Numune.olc CIKTISI - oldugu gibi."""
        mm = self.mm_asil if mm is None else mm
        ana = self.anahtar(hedef, F, R, mm)
        if ana in self._ob:
            return self._ob[ana]
        b = self.baglam(hedef)
        if not b or not b['uye']:
            return None
        self.havuz_hazirla([hedef])
        o = self._nm.olc(F.upper(), R.upper(), b['uye'], b['rakip'],
                         lo=self.lo, hi=self.hi, mm=mm)
        self._ob[ana] = o
        return o

    def puanla(self, hedef, F, R, duzey='', yan=False):
        """PANEL DERINLIGINDE puan. 'kat' ve 'durum' panelin karar() fonksiyonundan
        gelir; burada esik karsilastirmasi YAPILMAZ."""
        o1 = self.olc_ham(hedef, F, R, self.mm_asil)
        d1, g1, day1 = self.TP.karar(o1, hedef, duzey)
        s = dict(hedef=hedef, F=F.upper(), R=R.upper(),
                 kat=g1, durum=d1, dayanak=day1,
                 dcq=self.C.kat_dcq(g1) if g1 is not None else None,
                 kapsam=(o1 or {}).get('uye_kapsam_pay', ''),
                 uye_alt=(o1 or {}).get('uye_alt'),
                 enkotu_kutu=(o1 or {}).get('enkotu_kutu', ''),
                 kat_havuz=(o1 or {}).get('kat_havuz'),
                 rakip_olculen=(o1 or {}).get('rakip_olculen'),
                 rakip_toplam=(o1 or {}).get('rakip_toplam'),
                 uye_n=len((o1 or {}).get('uye') or []),
                 rakip_n=len((o1 or {}).get('rakip') or []),
                 urun_boylari=(o1 or {}).get('urun_boylari'))
        if yan:
            o3 = self.olc_ham(hedef, F, R, self.mm_yan)
            d3, g3, _ = self.TP.karar(o3, hedef, duzey)
            s['kat_mm3'], s['durum_mm3'] = g3, d3
        return s

    def onbellek_yaz(self):
        if not self.onbellek_yolu:
            return
        os.makedirs(os.path.dirname(self.onbellek_yolu), exist_ok=True)
        tmp = self.onbellek_yolu + '.tmp'
        json.dump(self._ob, open(tmp, 'w', encoding='utf-8'), default=str)
        os.replace(tmp, self.onbellek_yolu)
