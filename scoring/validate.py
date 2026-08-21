# -*- coding: utf-8 -*-
"""SON DOGRULAMA - PANELIN KENDI MOTORUYLA, PANEL DERINLIGINDE.

Toplu tarama hiz/bellek icin numpy motorunu (otorite=False) ve bazi siniflarda
sig derinligi kullandi. RAPORLANACAK her aday burada YENIDEN olculur:
    otorite=True (KutuOtorite / okuma_motoru)  +  derinlik 3000  +  mm<=1
yani single_protocol_measure.py'nin kullandigi yolun BIREBIR aynisi.

Tarama degeri ile dogrulama degeri FARKLI cikarsa RAPORA DOGRULAMA DEGERI girer
ve fark acikca yazilir.
"""
import os, sys, json, time, argparse
sys.path.insert(0, '/tmp/mrb')
from ortak_puanlayici import Puanlayici
KOK = '/sessions/dreamy-elegant-wozniak/mnt/PrimerTasarlama'


def kapsam_orani(s):
    k = (s.get('kapsam') or '')
    if '/' in k:
        a, b = k.split('/')
        try:
            return int(a) / max(int(b), 1)
        except ValueError:
            return 0.0
    return 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--threads', '--is', dest='isler', action='append', required=True,
                    help='hedef::kp_dosya::N')
    ap.add_argument('--duration', '--sure', dest='sure', type=float, default=36.0)
    ap.add_argument('--output', '--cikti', dest='cikti', default='/tmp/mrb/kontrol/dogrulama.json')
    g = ap.parse_args()
    t0 = time.time()

    isler = []
    for x in g.isler:
        p = x.split('::'); isler.append((p[0], p[1], int(p[2])))

    son = json.load(open(g.cikti, encoding='utf-8')) if os.path.exists(g.cikti) else {}
    P = Puanlayici(KOK, otorite=True, derinlik=3000,
                   onbellek_yolu='/tmp/mrb/onbellek/dogrulama_olcum.json')
    P.havuz_hazirla([h for h, _, _ in isler])
    print('havuz hazir (%d kutu, otorite=True, derinlik=3000) [%.1f sn]'
          % (len(P._yuklu), time.time() - t0), flush=True)

    for hedef, kpy, N in isler:
        if not os.path.exists(kpy):
            print('%s: %s YOK - atlandi' % (hedef, kpy)); continue
        v = list(json.load(open(kpy, encoding='utf-8'))['olcum'].values())
        evrensel = all(s.get('kat') is None for s in v)
        v.sort(key=(lambda s: -kapsam_orani(s)) if evrensel
               else (lambda s: -(s['kat'] if s.get('kat') is not None else -1)))
        son.setdefault(hedef, {})
        for s in v[:N]:
            ad = '%s|%s' % (s['F'], s['R'])
            if ad in son[hedef]:
                continue
            if time.time() - t0 > g.sure:
                break
            d = P.puanla(hedef, s['F'], s['R'], yan=False)
            d['tarama_kat'] = s.get('kat')
            d['tarama_kapsam'] = s.get('kapsam')
            d['grup'] = s.get('grup'); d['arms'] = s.get('arms', '')
            d['urun'] = s.get('urun'); d['tmF'] = s.get('tmF'); d['tmR'] = s.get('tmR')
            d['dtm'] = s.get('dtm')
            d.pop('urun_boylari', None)
            son[hedef][ad] = d
        json.dump(son, open(g.cikti, 'w', encoding='utf-8'), default=str)
        P.onbellek_yaz()
        print('--- %s  (dogrulanan %d)' % (hedef, len(son[hedef])), flush=True)
        for d in sorted(son[hedef].values(),
                        key=(lambda s: -kapsam_orani(s)) if evrensel
                        else (lambda s: -(s['kat'] if s.get('kat') is not None else -1)))[:5]:
            print('   DOGRU kat=%-8s dCq=%-6s %-11s kapsam=%-6s (tarama %s) %s %s/%s'
                  % (d['kat'], d['dcq'], d['durum'], d['kapsam'], d['tarama_kat'],
                     d['arms'] or '-', d['F'], d['R']), flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
