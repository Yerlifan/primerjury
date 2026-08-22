# -*- coding: utf-8 -*-
"""A BULK MEASUREMENT OVER MANY TARGETS, with the pool SHARED.

The targets of the same class use the same bins (class B: Petrimonas,
Proteiniphilum, Bacteroidales). Because building the pool takes about 20 s per call,
several targets are measured in a single process.

UNIVERSAL TARGETS: karar() returns OLCULEMEDI on those rows (the denominator of the
discrimination fold is undefined, there is no competitor set). For them the measure
is COVERAGE (uye_kapsam_pay) and the ranking is made by that. A discrimination fold
IS NOT INVENTED.

"""
import os, sys, json, time, argparse
sys.path.insert(0, '/tmp/mrb')
from shared_scorer import Puanlayici
KOK = '/sessions/dreamy-elegant-wozniak/mnt/PrimerTasarlama'


def kapsam_orani(s):
    k = s.get('kapsam') or ''
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
                    help='target::shortlist_file::checkpoint_file')
    ap.add_argument('--duration', dest='sure', type=float, default=36.0)
    # MEMORY: KutuHavuzu (otorite=False) holds about 160 MB per bin (1.5 kb reads) up
    # to about 400 MB (3.7 kb reads); the dominant item is the seed index, because of
    # its int64 key array. 20 bins x 3000 reads goes over 3.9 GB of RAM in the F2, F1
    # and A2 classes and the process is SIGKILLed by the OOM killer (measured: 3.77 GB
    # RSS, 23 s). The fix: the scan runs at a SHALLOWER depth (the memory falls
    # linearly) while the DECISION is always made at panel depth with the panel's own
    # engine (--authority 1).
    ap.add_argument('--depth', dest='derinlik', type=int, default=3000)
    ap.add_argument('--authority', dest='otorite', type=int, default=0)
    g = ap.parse_args()
    t0 = time.time()

    isler = []
    for x in g.isler:
        h, kd, kp = x.split('::')
        isler.append((h, kd, kp))

    P = Puanlayici(KOK, otorite=bool(g.otorite), derinlik=g.derinlik)
    P.havuz_hazirla([h for h, _, _ in isler])
    print(u'pool ready (%d bins, depth=%d, authority=%s) [%.1f s]'
          % (len(P._yuklu), g.derinlik, bool(g.otorite), time.time() - t0), flush=True)

    for h, kd, kpy in isler:
        D = json.load(open(kd))
        aday = []
        for grup, anahtar in (('siralama', 'kisa'), ('ARMS', 'arms'),
                              ('kontrol_rastgele', 'kontrol')):
            for i, c in enumerate(D.get(anahtar) or []):
                c = dict(c); c['grup'] = grup
                c['sira1'] = i if grup == 'siralama' else -1
                c['ad'] = '%s|%s' % (c['F'], c['R'])
                aday.append(c)
        kpy = kpy if (g.derinlik == 3000 and not g.otorite) else (
            kpy.replace('.json', '_d%d%s.json' % (g.derinlik, 'o' if g.otorite else '')))
        kp = json.load(open(kpy, encoding='utf-8')) if os.path.exists(kpy) else {}
        kp.setdefault('olcum', {})
        n = 0
        for c in aday:
            if c['ad'] in kp['olcum']:
                continue
            if time.time() - t0 > g.sure:
                break
            s = P.puanla(h, c['F'], c['R'])
            s.pop('urun_boylari', None)
            s.update(grup=c['grup'], sira1=c['sira1'], urun=c['urun'],
                     arms=c.get('arms', ''), rak_min=c.get('rak_min'),
                     ortak_kons=c.get('ortak_kons'), tmF=c.get('tmF'),
                     tmR=c.get('tmR'), dtm=c.get('dtm'))
            kp['olcum'][c['ad']] = s
            n += 1
            if n % 100 == 0:      # ARA KAYIT: surec oldurulurse ilerleme kaybolmasin
                json.dump(kp, open(kpy, 'w', encoding='utf-8'), default=str)
        json.dump(kp, open(kpy, 'w', encoding='utf-8'), default=str)
        v = list(kp['olcum'].values())
        kalan = len(aday) - len(kp['olcum'])
        evrensel = all(s['kat'] is None for s in v) and v
        print(u'%-44s +%-4d total %4d/%-4d remaining %4d  [%.0f s]'
              % (h[:44], n, len(kp['olcum']), len(aday), kalan, time.time() - t0))
        if evrensel:
            en = sorted(v, key=lambda s: -kapsam_orani(s))[:3]
            for s in en:
                print(u'    COVERAGE=%-7s %-11s %s %s/%s'
                      % (s['kapsam'], s['durum'], s['arms'] or '-', s['F'], s['R']))
        else:
            en = sorted(v, key=lambda s: -(s['kat'] if s['kat'] is not None else -1))[:3]
            for s in en:
                print(u'    layer=%-8s dCq=%-6s %-10s coverage=%-6s %-16s rank1=%-4s %s %s/%s'
                      % (s['kat'], s['dcq'], s['durum'], s['kapsam'], s['grup'],
                         s['sira1'], s['arms'] or '-', s['F'], s['R']))
        if kalan == 0:
            gec = [s for s in v if s['durum'] == 'ESIK USTU']
            r = [s for s in v if s['grup'] == 'siralama' and s['kat'] is not None]
            k = [s for s in v if s['grup'] == 'kontrol_rastgele' and s['kat'] is not None]
            print(u'    DONE. passing the threshold=%d' % len(gec))
            if r and k:
                br, bk = max(s['kat'] for s in r), max(s['kat'] for s in k)
                iyi = max(r, key=lambda s: s['kat'])
                print(u'    PRE-FILTER: best ranked position=%d/%d, ranked=%.2f control=%.2f -> %s'
                      % (iyi['sira1'], len(r), br, bk,
                         'BAGLAYICI' if bk > br else u'not binding'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
