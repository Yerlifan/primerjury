# -*- coding: utf-8 -*-
"""FARKLI_LOKUS - 1. adim: hedef basina FIZIKSEL LOKUS ENVANTERI.

Bu betik hicbir sey degistirmez, yalnizca olcer. Sorular:
  1. Hedefin uye kutulari hangileri, her birinin konsensusu kac bp?
  2. Konsensusta rRNA operonunun hangi parcalari var (16S/18S, ITS, 23S/28S)?
     Olcu: kurtarma_turu.capa_bul + bolgeler_kur (ayni capa listesi).
  3. Ham okumalarin boy dagilimi ne? Konsensus kisa ama okumalar uzunsa
     "baska lokus yok" demek YANLIS olur; bu yuzden okuma boyu ayri olculur.
  4. Uye konsensusleri birbirine ne kadar benziyor (kume heterojen mi)?

Cikti: JSON, /tmp altina. Bagli klasore yalniz betik yazilir.
"""
import os, sys, json, gzip, random, itertools

# The project root. It used to default to /tmp/fl/kok, the path of the private
# workspace this code grew in; on any other machine the module died at import
# with FileNotFoundError. The default is now derived from this file's own
# location, and _FL_KOK still overrides it.
KOK = os.environ.get('_FL_KOK') or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, KOK)
os.environ['_KURTARMA_KOK'] = KOK

from screening import targets as H, engine_gateway
import importlib.util as _iu
_sp = _iu.spec_from_file_location('kt', os.path.join(KOK, 'verification', 'recovery_round.py'))
kt = _iu.module_from_spec(_sp)
_sp.loader.exec_module(kt)

HEDEFLER = [
    'Methanothrix_soehngenii_turu',
    'Methanosarcina mazei / M. soligelidi grubu',
    'Bacteroidales_kumesi',
    'Petrimonas_cinsi',
    'Petriella_musispora',
    'Proteiniphilum_cinsi',
]


def okuma_boylari(yol, n=4000):
    """Ham fastq okuma boyu dagilimi - SUZGECSIZ (200-6000 filtresi yok).
    Filtre uygulanirsa 'uzun okuma var mi' sorusu kendi kendini yanitlar."""
    boy = []
    ac = gzip.open if yol.endswith('.gz') else open
    with ac(yol, 'rt', errors='ignore') as fh:
        for i, s in enumerate(fh):
            if i % 4 == 1:
                boy.append(len(s.strip()))
            if len(boy) >= n:
                break
    if not boy:
        return None
    boy.sort()
    q = lambda p: boy[min(len(boy) - 1, int(len(boy) * p))]
    return dict(n=len(boy), min=boy[0], q25=q(.25), medyan=q(.5), q75=q(.75),
                q95=q(.95), q99=q(.99), maks=boy[-1],
                ust_1500=sum(1 for b in boy if b > 1500),
                ust_2500=sum(1 for b in boy if b > 2500),
                ust_4000=sum(1 for b in boy if b > 4000))


def benzerlik(a, b, ornek=1200):
    """Kaba yerel benzerlik: k-mer Jaccard + en iyi kayma ile ozdeslik.
    Hizalayici yok; amac 'kume heterojen mi' sorusunu tek sayiyla yanitlamak."""
    k = 12
    A = set(a[i:i + k] for i in range(len(a) - k))
    B = set(b[i:i + k] for i in range(len(b) - k))
    jac = len(A & B) / float(len(A | B)) if (A | B) else 0.0
    # kayma taramasi ile en iyi ozdeslik
    kisa, uzun = (a, b) if len(a) <= len(b) else (b, a)
    if len(kisa) > ornek:
        kisa = kisa[:ornek]
    en = 0.0
    adim = max(1, (len(uzun) - len(kisa)) // 200 or 1)
    for off in range(0, max(1, len(uzun) - len(kisa) + 1), adim):
        w = uzun[off:off + len(kisa)]
        e = sum(1 for x, y in zip(kisa, w) if x == y) / float(len(kisa))
        if e > en:
            en = e
    return round(jac * 100, 2), round(en * 100, 2)


def main():
    kons = H.konsensusler()
    kut = {k['kutu']: k for k in H.kutular()}
    kmap = {d['kutu']: d['dizi'] for d in kons}
    panel = {}
    import csv
    p = os.path.join(KOK, 'TEK_PROTOKOL_SONUC', 'panel_tek_protokol.tsv')
    for r in csv.DictReader((l for l in open(p, encoding='utf-8')
                             if not l.startswith('#')), delimiter='\t'):
        panel[r['hedef']] = r

    cikti = {}
    for ad in HEDEFLER:
        r = panel[ad]
        b = H.hedef_baglami(r, kons=kons, kut=list(kut.values()))
        uye = b['uye_kutu']
        rak = b['rakip_kutu']
        d = dict(hedef=ad, sinif=r['sinif'], uyelik_kaynagi=b['anahtar'],
                 uye_tax=b['uye_tax'], haric=b['haric'],
                 uye_kutu=[], rakip_kutu=[k['kutu'] for k in rak])
        for k in uye:
            kd = kmap.get(k['kutu'], '')
            capa = kt.capa_bul(engine_gateway, kd) if kd else {}
            bol = []
            if kd:
                bol = [(a, x, y, s) for a, x, y, s in
                       kt.bolgeler_kur(engine_gateway, kd, lambda *_a, **_k: None)[0]]
            d['uye_kutu'].append(dict(
                kutu=k['kutu'], sinif=k['sinif'], taxid=k['taxid'],
                kons_bp=len(kd), capa=capa,
                bolge=[dict(ad=a, bas=x, son=y, bp=y - x, kaynak=s) for a, x, y, s in bol],
                okuma=okuma_boylari(k['yol'])))
        # uye konsensusleri arasi benzerlik
        cift = []
        us = [(k['kutu'], kmap.get(k['kutu'], '')) for k in uye if kmap.get(k['kutu'])]
        for (n1, s1), (n2, s2) in itertools.combinations(us, 2):
            j, e = benzerlik(s1, s2)
            cift.append(dict(a=n1, b=n2, jaccard12=j, ozdeslik=e,
                             bp_a=len(s1), bp_b=len(s2)))
        d['uye_benzerlik'] = cift
        cikti[ad] = d
        print(u'%-46s member %2d bins, cons %s' % (
            ad[:46], len(uye),
            '/'.join(str(x['kons_bp']) for x in d['uye_kutu'])), flush=True)

    yol = os.environ.get('_FL_CIKTI') or os.path.join(KOK, 'SONUCLAR', 'envanter.json')
    json.dump(cikti, open(yol, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    print('yazildi:', yol)


if __name__ == '__main__':
    main()
