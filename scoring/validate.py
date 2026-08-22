# -*- coding: utf-8 -*-
"""THE FINAL CONFIRMATION, WITH THE PANEL'S OWN ENGINE, AT PANEL DEPTH.

The bulk scan used the numpy engine (otorite=False) for speed and memory, and a
shallow depth in some classes. EVERY candidate TO BE REPORTED is measured again here:
    otorite=True (KutuOtorite / okuma_motoru)  +  a depth of 3000  +  mm<=1
that is, exactly the route single_protocol_measure.py uses.

If the scan value and the confirmation value COME OUT DIFFERENT, THE CONFIRMATION
VALUE goes into the report and the difference is written out plainly.

"""
import os, sys, json, time, argparse
sys.path.insert(0, '/tmp/mrb')
from shared_scorer import Puanlayici
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
    ap.add_argument('--threads', dest='isler', action='append', required=True,
                    help='target::checkpoint_file::N')
    ap.add_argument('--duration', dest='sure', type=float, default=36.0)
    ap.add_argument('--output', dest='cikti', default='/tmp/mrb/kontrol/dogrulama.json')
    g = ap.parse_args()
    t0 = time.time()

    isler = []
    for x in g.isler:
        p = x.split('::'); isler.append((p[0], p[1], int(p[2])))

    son = json.load(open(g.cikti, encoding='utf-8')) if os.path.exists(g.cikti) else {}
    P = Puanlayici(KOK, otorite=True, derinlik=3000,
                   onbellek_yolu='/tmp/mrb/onbellek/dogrulama_olcum.json')
    P.havuz_hazirla([h for h, _, _ in isler])
    print(u'pool ready (%d bins, authority=True, depth=3000) [%.1f s]'
          % (len(P._yuklu), time.time() - t0), flush=True)

    for hedef, kpy, N in isler:
        if not os.path.exists(kpy):
            print(u'%s: %s ABSENT - skipped' % (hedef, kpy)); continue
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
            print(u'   CORRECT fold=%-8s dCq=%-6s %-11s coverage=%-6s (scan %s) %s %s/%s'
                  % (d['kat'], d['dcq'], d['durum'], d['kapsam'], d['tarama_kat'],
                     d['arms'] or '-', d['F'], d['R']), flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
