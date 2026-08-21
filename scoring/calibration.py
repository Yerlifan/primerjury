# -*- coding: utf-8 -*-
"""KALIBRASYON SINAMASI - ortak puanlayici panelin sayilarini yeniden uretiyor mu?

Panelin KENDI kontrol noktalari (TEK_PROTOKOL_SONUC/kontrol/*.json) beklenen
deger kaynagidir: her cift icin uye/rakip kutu listesi, kutu basina k ve n,
ve kat_enkotu yazilidir. Bu dosyalar YALNIZ OKUNUR.

Karsilastirma UC KADEMELIDIR - yalnizca son sayiya bakmak yetmez:
  1) uye kutu kumesi     (sira dahil)
  2) rakip kutu kumesi   (sira dahil)
  3) kutu basina k/n     (ham sayilar)
  4) kat_enkotu          (karar sayisi)
Ilk uc kademe tutmadan dorduncunun tutmasi TESADUFTUR ve guvenilmez.

Parti parti kosar (45 sn bash tavani). Durum /tmp/mrb/kontrol/kalib.json.
"""
import os, sys, json, glob, argparse
sys.path.insert(0, '/tmp/mrb')
from ortak_puanlayici import Puanlayici

KOK = '/sessions/dreamy-elegant-wozniak/mnt/PrimerTasarlama'
KP = '/tmp/mrb/kontrol/kalib.json'


def beklenenler():
    """Panelin kontrol noktalarindan beklenen degerleri cikarir."""
    out = []
    for y in sorted(glob.glob(os.path.join(KOK, 'TEK_PROTOKOL_SONUC', 'kontrol', '*.json'))):
        v = json.load(open(y, encoding='utf-8'))
        o1 = (v.get('olcum') or {}).get('1')
        out.append(dict(hedef=v['hedef'], F=v.get('F', ''), R=v.get('R', ''),
                        duzey=v.get('duzey', ''), kaynak=v.get('kaynak', ''),
                        dosya=os.path.basename(y),
                        bek_kat=(o1 or {}).get('kat_enkotu'),
                        bek_uye=[x[0] for x in (o1 or {}).get('uye') or []],
                        bek_rakip=[x[0] for x in (o1 or {}).get('rakip') or []],
                        bek_uye_kn={x[0]: (x[1], x[2]) for x in (o1 or {}).get('uye') or []},
                        bek_rakip_kn={x[0]: (x[1], x[2]) for x in (o1 or {}).get('rakip') or []},
                        bek_enkotu=(o1 or {}).get('enkotu_kutu'),
                        bek_uye_alt=(o1 or {}).get('uye_alt'),
                        var_olcum=o1 is not None))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--bas', type=int, default=0)
    p.add_argument('--adet', type=int, default=4)
    a = p.parse_args()

    bek = beklenenler()
    os.makedirs('/tmp/mrb/kontrol', exist_ok=True)
    durum = json.load(open(KP, encoding='utf-8')) if os.path.exists(KP) else {}

    P = Puanlayici(KOK, onbellek_yolu='/tmp/mrb/onbellek/kalib_olcum.json')
    dilim = bek[a.bas:a.bas + a.adet]
    if not dilim:
        print('parti bos - bitti'); return 0
    P.havuz_hazirla([d['hedef'] for d in dilim])

    for d in dilim:
        h = d['hedef']
        o = P.olc_ham(h, d['F'], d['R'], P.mm_asil)
        s = P.puanla(h, d['F'], d['R'], d['duzey'])
        got_uye = [x[0] for x in (o or {}).get('uye') or []]
        got_rak = [x[0] for x in (o or {}).get('rakip') or []]
        got_uye_kn = {x[0]: (x[1], x[2]) for x in (o or {}).get('uye') or []}
        got_rak_kn = {x[0]: (x[1], x[2]) for x in (o or {}).get('rakip') or []}
        kn_fark = []
        for k, v in d['bek_uye_kn'].items():
            if tuple(got_uye_kn.get(k, ())) != tuple(v):
                kn_fark.append('uye:%s bek%s got%s' % (k, v, got_uye_kn.get(k)))
        for k, v in d['bek_rakip_kn'].items():
            if tuple(got_rak_kn.get(k, ())) != tuple(v):
                kn_fark.append('rakip:%s bek%s got%s' % (k, v, got_rak_kn.get(k)))
        r = dict(
            hedef=h, kaynak=d['kaynak'],
            uye_kume_ayni=(sorted(got_uye) == sorted(d['bek_uye'])),
            rakip_kume_ayni=(sorted(got_rak) == sorted(d['bek_rakip'])),
            kn_ayni=(not kn_fark), kn_fark=kn_fark[:6],
            bek_kat=d['bek_kat'], got_kat=(o or {}).get('kat_enkotu'),
            kat_ayni=((o or {}).get('kat_enkotu') == d['bek_kat']),
            bek_enkotu=d['bek_enkotu'], got_enkotu=(o or {}).get('enkotu_kutu'),
            bek_uye_alt=d['bek_uye_alt'], got_uye_alt=(o or {}).get('uye_alt'),
            durum=s['durum'], olculebilir=d['var_olcum'] and d['bek_kat'] is not None,
            bek_uye_n=len(d['bek_uye']), got_uye_n=len(got_uye),
            bek_rakip_n=len(d['bek_rakip']), got_rakip_n=len(got_rak))
        durum[h] = r
        bayrak = 'TAM' if (r['uye_kume_ayni'] and r['rakip_kume_ayni']
                           and r['kn_ayni'] and r['kat_ayni']) else 'FARK'
        print('%-52s bek=%-8s got=%-8s  %s' % (h[:52], d['bek_kat'], r['got_kat'], bayrak))
        if bayrak == 'FARK':
            print('     uye_kume=%s rakip_kume=%s kn=%s kat=%s' % (
                r['uye_kume_ayni'], r['rakip_kume_ayni'], r['kn_ayni'], r['kat_ayni']))
            for f in r['kn_fark']:
                print('     ', f)

    P.onbellek_yaz()
    json.dump(durum, open(KP, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    print('--- toplam kayit: %d / %d' % (len(durum), len(bek)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
