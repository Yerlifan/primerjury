# -*- coding: utf-8 -*-
"""ADAY URETIMI + KONSENSUS ON ELEMESI (asama 1).  [duzeltilmis surum]

IKI DUZELTME (ilk surum yanlisti, olcumle yakalandi):

 1) GERI PRIMER YONU. Numune.olc geri primeri rc(R) olarak arar. Yani
    konsensus uzerinde bulunan pencere rc(R)'dir ve GERI PRIMER = rc(pencere).
    Ilk surum R=pencere yaziyordu; 196 adayin hepsi butun uye kutularda 0,00
    urun verdi. Sifir, biyolojik sonuc degil YON HATASI imzasidir.

 2) "BUTUN uye konsensuslara baglanmali" sarti FAZLA KATIYDI. Panelin KENDI
    Petrimonas ileri primeri bu sarti gecmiyor (kapsami 1/3). Sart kaldirildi;
    baglanma SAYISI artik bir eleme degil, SIRALAMA olcutudur.

Pencereler BUTUN uye konsensuslardan uretilir (tek omurga degil). Eslestirme
en uzun uye konsensusu (referans) uzerinde yapilir - urunun olusabilmesi icin
iki primerin ayni sablona baglanmasi gerekir.
"""
import os, sys, json, argparse, time
sys.path.insert(0, '/tmp/mrb')
from shared_scorer import Puanlayici

KOK = '/sessions/dreamy-elegant-wozniak/mnt/PrimerTasarlama'
UZ_ALT, UZ_UST = 18, 25
GC_ALT, GC_UST = 35.0, 65.0


def konsensus_tablosu(kok):
    import csv
    ix = os.path.join(kok, 'konsensus_kanonik', 'INDEKS.tsv')
    out = {}
    for r in csv.DictReader(open(ix, encoding='utf-8'), delimiter='\t'):
        y = os.path.join(os.path.dirname(ix), r['dosya'])
        if not os.path.exists(y):
            continue
        out[r['kutu']] = ''.join(l.strip() for l in open(y, encoding='utf-8', errors='ignore')
                                 if not l.startswith('>')).upper()
    return out


def ucuz_gecti(s):
    if any(c not in 'ACGT' for c in s):
        return False
    g = 100.0 * sum(s.count(c) for c in 'GC') / len(s)
    if not (GC_ALT <= g <= GC_UST):
        return False
    return not any(b * 4 in s for b in 'ACGT')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target', '--hedef', dest='hedef', required=True)
    ap.add_argument('--output', '--cikti', dest='cikti', required=True)
    ap.add_argument('--sinama_F', default='')
    ap.add_argument('--sinama_R', default='')
    g = ap.parse_args()

    t0 = time.time()
    P = Puanlayici(KOK)
    motor = P.motor
    b = P.baglam(g.hedef)
    if not b:
        sys.exit(u'ERROR: there is no membership for %s' % g.hedef)
    kons = konsensus_tablosu(KOK)
    uye_ad = [k['kutu'] for k in b['uye'] if k['kutu'] in kons]
    rak_ad = [k['kutu'] for k in b['rakip'] if k['kutu'] in kons]
    uye_sq = [kons[a] for a in uye_ad]
    rak_sq = [kons[a] for a in rak_ad]
    eksik_uye = [k['kutu'] for k in b['uye'] if k['kutu'] not in kons]
    if not uye_sq:
        sys.exit(u'ERROR: there is no member consensus')
    ref = max(range(len(uye_sq)), key=lambda i: len(uye_sq[i]))
    print(u'%s | member cons %d (ref=%s %d bp) | competitor cons %d | members with no consensus %d'
          % (g.hedef, len(uye_sq), uye_ad[ref], len(uye_sq[ref]), len(rak_sq), len(eksik_uye)))

    uye_hv = motor.tarayici.Havuz(uye_sq)
    rak_hv = motor.tarayici.Havuz(rak_sq) if rak_sq else None

    pencere = set()
    for s in uye_sq:
        for i in range(len(s)):
            for L in range(UZ_ALT, UZ_UST + 1):
                if i + L > len(s):
                    break
                w = s[i:i + L]
                if w not in pencere and ucuz_gecti(w):
                    pencere.add(w)
    print('  tekil pencere (ucuz suzgec): %d   [%.1f sn]' % (len(pencere), time.time() - t0))

    # --- SINAMA: bilinen calisan bir cift bu pencere kumesinde var mi?
    if g.sinama_F:
        rcR = motor.rc(g.sinama_R)
        print('  SINAMA  panel F pencerede: %s | panel rc(R) pencerede: %s'
              % (g.sinama_F in pencere, rcR in pencere))

    ileri, geri, rakb = {}, {}, {}
    for w in pencere:
        fs = motor.find_sites(uye_hv.enc, w, 1, True, (-1, -2))
        if fs:
            yer = {}
            for p, _ in fs:
                k = int(uye_hv.sid[p])
                if k >= 0:
                    yer.setdefault(k, int(p) - int(uye_hv.sid[:p].__len__() * 0))
            # referans uzerindeki yerel pozisyon
            loc = {}
            for p, _ in fs:
                k = int(uye_hv.sid[p])
                if k >= 0 and k not in loc:
                    loc[k] = int(p)
            ileri[w] = loc
        rs = motor.find_sites(uye_hv.enc, w, 1, True, (0, 1))
        if rs:
            loc = {}
            for p, _ in rs:
                k = int(uye_hv.sid[p])
                if k >= 0 and k not in loc:
                    loc[k] = int(p)
            geri[w] = loc
    print('  ileri baglanan %d | geri baglanan %d   [%.1f sn]'
          % (len(ileri), len(geri), time.time() - t0))

    if rak_hv is not None:
        for w in set(ileri) | set(geri):
            n = 0
            for tp in ((-1, -2), (0, 1)):
                fs = motor.find_sites(rak_hv.enc, w, 1, True, tp)
                n = max(n, len({int(rak_hv.sid[p]) for p, _ in fs if rak_hv.sid[p] >= 0}))
            rakb[w] = n
    print(u'  the competitor binding profile   [%.1f s]' % (time.time() - t0))

    if g.sinama_F:
        print('  SINAMA  panel F ileri-baglandi: %s | panel rc(R) geri-baglandi: %s'
              % (g.sinama_F in ileri, motor.rc(g.sinama_R) in geri))

    json.dump(dict(hedef=g.hedef, uye_ad=uye_ad, rak_ad=rak_ad, ref=ref,
                   uye_uz=[len(x) for x in uye_sq], eksik_uye=eksik_uye,
                   ileri={k: {str(a): c for a, c in v.items()} for k, v in ileri.items()},
                   geri={k: {str(a): c for a, c in v.items()} for k, v in geri.items()},
                   rakb=rakb, pencere_n=len(pencere)),
              open(g.cikti, 'w'), default=str)
    print(u'  written: %s   [%.1f s]' % (g.cikti, time.time() - t0))
    return 0


if __name__ == '__main__':
    sys.exit(main())
