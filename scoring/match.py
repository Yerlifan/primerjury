# -*- coding: utf-8 -*-
"""ASAMA 1b: eslestirme + geometri + KISA LISTE.  [duzeltilmis surum]

GERI PRIMER = rc(pencere).  Konsensus uzerinde bulunan pencere rc(R)'dir,
cunku Numune.olc geri primeri rc(R) olarak arar. Geometri (hairpin/homodimer)
GERCEK primer dizisi rc(pencere) uzerinde olculur - Tm ve GC ters tumleyende
ayni kalir ama yapi olcumleri KALMAZ.

ARMS: 3' son iki baz TAM eslesmek zorunda oldugu icin kasitli ikinci uyumsuzluk
yalniz -3 ve -4 konumlarina konur (-1/-2'ye konursa uye de cogalmaz).
"""
import os, sys, json, argparse, bisect, time
sys.path.insert(0, '/tmp/mrb')
from shared_scorer import Puanlayici

KOK = '/sessions/dreamy-elegant-wozniak/mnt/PrimerTasarlama'
TM_ALT, TM_UST = 58.0, 62.5
DTM_UST, URUN_ALT, URUN_UST = 1.5, 60, 400
IDEAL = (60, 150)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', '--girdi', dest='girdi', required=True)
    ap.add_argument('--output', '--cikti', dest='cikti', required=True)
    ap.add_argument('--short', '--kisa', dest='kisa', type=int, default=70)
    ap.add_argument('--arms', type=int, default=16)
    ap.add_argument('--check', '--kontrol', dest='kontrol', type=int, default=150)
    g = ap.parse_args()

    t0 = time.time()
    D = json.load(open(g.girdi))
    P = Puanlayici(KOK)
    rc = P.motor.rc
    from screening import geometry as G

    REF = str(D['ref'])
    UYE_N = len(D['uye_ad'])
    ileri = {w: {int(k): int(v) for k, v in d.items()} for w, d in D['ileri'].items()}
    geri = {w: {int(k): int(v) for k, v in d.items()} for w, d in D['geri'].items()}
    rakb = {k: int(v) for k, v in D['rakb'].items()}
    refi = int(D['ref'])

    # --- geometri: ileri rolu pencerenin KENDISI, geri rolu rc(pencere)
    gF, gR = {}, {}
    for w, loc in ileri.items():
        if refi not in loc:
            continue
        t = G.tm(w)
        if not (TM_ALT <= t <= TM_UST):
            continue
        m = G.olc(w)
        if G.sabit_gecti(m):
            gF[w] = m
    for w, loc in geri.items():
        if refi not in loc:
            continue
        r = rc(w)
        t = G.tm(r)
        if not (TM_ALT <= t <= TM_UST):
            continue
        m = G.olc(r)
        if G.sabit_gecti(m):
            gR[w] = m
    print('  geometri gecen: ileri %d, geri %d   [%.1f sn]' % (len(gF), len(gR), time.time() - t0))

    fl = sorted((ileri[w][refi], w) for w in gF)
    rl = sorted((geri[w][refi], w) for w in gR)
    rpos = [p for p, _ in rl]

    ciftler = {}
    for pf, wf in fl:
        lo = bisect.bisect_left(rpos, pf + len(wf))
        for j in range(lo, len(rl)):
            pr, wr = rl[j]
            urun = pr + len(wr) - pf
            if urun < URUN_ALT:
                continue
            if urun > URUN_UST:
                break
            mF, mR = gF[wf], gR[wr]
            if abs(mF['tm'] - mR['tm']) > DTM_UST:
                continue
            Rseq = rc(wr)
            ortak = len(set(ileri[wf]) & set(geri[wr]))
            rmin = min(rakb.get(wf, 0), rakb.get(wr, 0))
            ciftler[(wf, Rseq)] = dict(
                F=wf, R=Rseq, urun=urun, pozF=pf, pozR=pr, ortak_kons=ortak,
                tmF=mF['tm'], tmR=mR['tm'], dtm=round(abs(mF['tm'] - mR['tm']), 2),
                gcF=mF['gc'], gcR=mR['gc'], rak_min=rmin,
                rak_top=rakb.get(wf, 0) + rakb.get(wr, 0),
                ideal_ceza=0 if IDEAL[0] <= urun <= IDEAL[1] else (1 if urun <= 250 else 2))
    print(u'  candidate pairs: %d   [%.1f s]' % (len(ciftler), time.time() - t0))
    if not ciftler:
        json.dump(dict(hedef=D['hedef'], uye_ad=D['uye_ad'], rak_ad=D['rak_ad'],
                       pencere_n=D['pencere_n'], aday_cift_n=0, kisa=[], arms=[]),
                  open(g.cikti, 'w'), default=str)
        print(u'  WARNING: no pair could be produced'); return 0

    sirali = sorted(ciftler.values(),
                    key=lambda x: (-x['ortak_kons'], x['rak_min'], x['rak_top'],
                                   x['ideal_ceza'],
                                   round(abs((x['tmF'] + x['tmR']) / 2 - 60.5), 2)))
    kisa, kul = [], []
    for c in sirali:
        if len(kisa) < 5 or all(abs(c['pozF'] - u) >= 25 for u in kul):
            kisa.append(c); kul.append(c['pozF'])
        if len(kisa) >= g.kisa:
            break
    for c in sirali:
        if len(kisa) >= g.kisa:
            break
        if c not in kisa:
            kisa.append(c)

    arms = []
    for c in kisa[:g.arms]:
        for rol in ('F', 'R'):
            s = c[rol]
            for konum in (-3, -4):
                for yeni in 'ACGT':
                    if yeni == s[konum]:
                        continue
                    v = s[:len(s) + konum] + yeni + s[len(s) + konum + 1:]
                    if any(bb * 4 in v for bb in 'ACGT'):
                        continue
                    t = G.tm(v)
                    if not (TM_ALT - 2.0 <= t <= TM_UST):
                        continue
                    d = dict(c); d[rol] = v; d['tm' + rol] = t
                    d['arms'] = '%s %d. konum %s->%s' % (rol, konum, s[konum], yeni)
                    d['dtm'] = round(abs(d['tmF'] - d['tmR']), 2)
                    if d['dtm'] <= DTM_UST + 0.5:
                        arms.append(d)
    print('  ARMS: %d   [%.1f sn]' % (len(arms), time.time() - t0))

    # ON ELEME BAGLAYICILIK KONTROLU: siralamanin GERI KALANINDAN rastgele ornek.
    # Kisa listede olmayan adaylardan alinir; ayni derinlikte olculur.
    import random as _r
    kalan = [c for c in sirali if c not in kisa]
    _r.Random(20260810).shuffle(kalan)
    kontrol = kalan[:g.kontrol]
    json.dump(dict(hedef=D['hedef'], uye_ad=D['uye_ad'], rak_ad=D['rak_ad'],
                   pencere_n=D['pencere_n'], aday_cift_n=len(ciftler),
                   kisa=kisa, arms=arms, kontrol=kontrol), open(g.cikti, 'w'), default=str)
    print(u'  control sample (drawn at random from the rest of the ranking): %d' % len(kontrol))
    print('  kisa %d + ARMS %d -> %s' % (len(kisa), len(arms), g.cikti))
    for c in kisa[:4]:
        print(u'    shared_cons=%d/%d comp_min=%d product=%d  %s / %s'
              % (c['ortak_kons'], UYE_N, c['rak_min'], c['urun'], c['F'], c['R']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
