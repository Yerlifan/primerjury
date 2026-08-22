# -*- coding: utf-8 -*-
"""CANDIDATE PRODUCTION plus A CONSENSUS PREFILTER (stage 1).  [the corrected version]

TWO CORRECTIONS (the first version was wrong and it was caught by measurement):

 1) THE REVERSE PRIMER'S ORIENTATION. Numune.olc searches for the reverse primer as
    rc(R). So the window found on the consensus IS rc(R) and THE REVERSE PRIMER =
    rc(window). The first version wrote R=window; all 196 candidates gave 0.00
    product in every member bin. A zero is not a biological result, it is the
    signature of AN ORIENTATION FAULT.

 2) The condition "it must bind ALL the member consensuses" WAS TOO STRICT. The
    panel's OWN Petrimonas forward primer does not pass it (its coverage is 1 of 3).
    The condition was removed; the binding COUNT is now a RANKING criterion rather
    than an elimination.

The windows are produced from ALL the member consensuses (not from a single
backbone). The matching is done on the longest member consensus (the reference),
because for a product to form the two primers have to bind the same template.

"""
import os, sys, json, argparse, time
sys.path.insert(0, '/tmp/mrb')
from shared_scorer import Puanlayici

KOK = '/sessions/dreamy-elegant-wozniak/mnt/PrimerTasarlama'
UZ_ALT, UZ_UST = 18, 25
GC_ALT, GC_UST = 35.0, 65.0


def konsensus_tablosu(kok):
    import csv
    ix = os.path.join(kok, 'canonical_consensus', 'INDEX.tsv')
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
    ap.add_argument('--target', dest='hedef', required=True)
    ap.add_argument('--output', dest='cikti', required=True)
    ap.add_argument('--test-forward', default='')
    ap.add_argument('--test-reverse', default='')
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

    # --- A TEST: is a pair known to work present in this window set?
    if g.test_forward:
        rcR = motor.rc(g.test_reverse)
        print('  SINAMA  panel F pencerede: %s | panel rc(R) pencerede: %s'
              % (g.test_forward in pencere, rcR in pencere))

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

    if g.test_forward:
        print('  SINAMA  panel F ileri-baglandi: %s | panel rc(R) geri-baglandi: %s'
              % (g.test_forward in ileri, motor.rc(g.test_reverse) in geri))

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
