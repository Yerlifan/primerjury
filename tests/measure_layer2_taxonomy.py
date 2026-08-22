# -*- coding: utf-8 -*-
"""THE LAYER 2 MEASUREMENT: length based separation beside taxonomic separation.

WHAT IT ASKS
------------
Layer 2, the local databases, used to decide whether a product came from INSIDE
or OUTSIDE the target by LENGTH alone. A taxonomic separation was added beside
it. This script measures both ON THE SAME SCAN and shows the difference.

It WRITES to no panel file. It only reads and prints a table.

WHY IT MATTERS
--------------
It was measured in layer 3, MFEprimer: of the 1,605 amplicons counted off target,
95.7 per cent came from INSIDE the target's own clade and differed only in
length. This script measures whether the same mistake is present in layer 2.

TO RUN IT
---------
    python3 tests/measure_layer2_taxonomy.py --small
    python3 tests/measure_layer2_taxonomy.py --db SILVA_138.2_SSURef_NR99.fasta
    python3 tests/measure_layer2_taxonomy.py --target Bakteri_universal

--small uses the small RefSeq sets alone, about 75 MB and a few minutes. That is
not coverage but evidence THAT THE METHOD WORKS, plus a first order of magnitude.
For full coverage the large sets such as SILVA and UNITE have to be given
separately; that takes hours and this script DOES NOT do it SILENTLY.
"""
from __future__ import print_function
import argparse
import csv
import hashlib
import io
import json
import os
import re
import sys
import time

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, KOK)

from screening import global_scan as KT          # noqa: E402
from screening import taxonomy as TX               # noqa: E402
import verification.mfeprimer_layer as MK                        # noqa: E402

KUCUK = ['archaea.16S.fna', 'fungi.18SrRNA.fna', 'fungi.28SrRNA.fna',
         'fungi.ITS.fna', 'bacteria.16S.fna']

KAYNAK = os.path.join(KOK, 'VERIFICATION_RESULT', 'dogrulama_uc_sutun.tsv')
URUN_ALT, URUN_UST, BOY_TOL = 70, 400, 10


def ciftleri_oku():
    if not os.path.exists(KAYNAK):
        sys.stderr.write(u'no source of pairs: %s\n' % KAYNAK)
        sys.exit(2)
    out = []
    with io.open(KAYNAK, encoding='utf-8') as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith('#')), delimiter='\t'):
            h = (r.get('hedef') or '').strip()
            F = (r.get('F') or '').strip().upper()
            R = (r.get('R') or '').strip().upper()
            bp = (r.get('urun_bp') or '').strip()
            if h and F and R:
                out.append(dict(hedef=h, F=F, R=R,
                                bek=int(bp) if bp.isdigit() else None))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--small', dest='kucuk', action='store_true')
    p.add_argument('--db', dest='vtb', nargs='*', default=None)
    p.add_argument('--target', dest='hedef', default=None)
    a = p.parse_args()

    ciftler = ciftleri_oku()
    if a.hedef:
        ciftler = [c for c in ciftler if a.hedef.lower() in c['hedef'].lower()]
    if not ciftler:
        sys.stderr.write(u'no pair found\n'); return 2

    klad = MK.klad_tablosu(KOK)
    yok = [c['hedef'] for c in ciftler if c['hedef'] not in klad]
    if yok:
        print(u'WARNING: %d targets have no definition in target_clades.tsv, so taxonomic separation CANNOT be made for them:' % len(yok))
        for h in yok:
            print('   - %s' % h)
        print()

    def sinifla(aday_ad, baslik, db_ad):
        if aday_ad not in klad:
            return 'bilinmiyor'
        alan, jetonlar, _k = klad[aday_ad]
        return TX.sinifla(baslik, db_ad, jetonlar, alan)

    vtbler = a.vtb if a.vtb else (KUCUK if a.kucuk else KUCUK)
    refdb = os.path.join(KOK, 'REFERENCE_DB')

    adaylar = [dict(ad=c['hedef'], F=c['F'], R=c['R'], lo=URUN_ALT, hi=URUN_UST)
               for c in ciftler]
    bek = {c['hedef']: c['bek'] for c in ciftler}

    top = {c['hedef']: dict(urun=0, ayni_boyda=0, boy_disi=0,
                            sinif={k: 0 for k in KT.SINIFLAR}, tarandi=0)
           for c in ciftler}

    # ------------------------------------------------------------------
    # RESILIENCE AGAINST INTERRUPTION
    #
    # In the first version tara() was called with durum_yolu=None, that is, with the
    # checkpoint OFF. That was a fault: kuresel_tarama keeps a checkpoint chunk by chunk
    # already (katman1_yerel uses it), and this project's rule is that a long run must
    # survive interruption. With the checkpoint off, a 28 minute SILVA scan was lost
    # entirely on a single shutdown.
    #
    # Two levels of protection:
    #   1) a chunk level checkpoint for tara() (durum_yolu)
    #   2) a partial result is written to disk as EVERY DATABASE finishes; even if the
    #      run is cut, the measurement up to that point stands and can be read.
    #
    # The checkpoint key carries the candidate SEQUENCES too: sealed by name alone, an
    # old scan comes back silently when a sequence changes (that fault was measured in
    # dogrulama_turu on 2026-08-10).
    kn_dizin = os.path.join(KOK, 'SCREENING_RESULT', 'kontrol', 'A2_olcum')
    try:
        os.makedirs(kn_dizin)
    except OSError:
        pass
    imza = hashlib.md5('|'.join(sorted(
        '%s>%s<%s' % (x['ad'], x['F'], x['R']) for x in adaylar)
    ).encode('utf-8')).hexdigest()[:10]
    kismi_yol = os.path.join(kn_dizin, 'kismi_sonuc_%s.json' % imza)

    tamamlanan = []
    if os.path.exists(kismi_yol):
        try:
            _k = json.load(io.open(kismi_yol, encoding='utf-8'))
            if _k.get('imza') == imza:
                top = _k['top']
                tamamlanan = _k.get('tamamlanan', [])
                print(u'resuming from an earlier run: %d databases already measured (%s)'
                      % (len(tamamlanan), ', '.join(tamamlanan)))
        except Exception as e:
            print(u'the partial result could not be read (%s), so it is being measured from scratch' % type(e).__name__)

    def kismi_yaz():
        gecici = kismi_yol + '.tmp'
        with io.open(gecici, 'w', encoding='utf-8') as fh:
            fh.write(json.dumps(dict(imza=imza, tamamlanan=tamamlanan, top=top),
                                ensure_ascii=False))
        if os.path.exists(kismi_yol):
            os.remove(kismi_yol)
        os.rename(gecici, kismi_yol)

    print(u'%d pairs, %d databases' % (len(ciftler), len(vtbler)))
    for d in vtbler:
        if d in tamamlanan:
            print(u'  SKIPPED (already measured): %s' % d)
            continue
        yol = os.path.join(refdb, d)
        if not os.path.exists(yol):
            print(u'  SKIPPED (absent): %s' % d)
            continue
        t0 = time.time()
        print('  taraniyor: %-34s (%s)' % (d, '%.0f MB' % (os.path.getsize(yol) / 1e6)),
              end='', flush=True)
        dy = os.path.join(kn_dizin, 'tarama_%s_%s.pkl'
                          % (re.sub(r'\W+', '_', d), imza))
        res = KT.tara(adaylar, db=yol, durum_yolu=dy, siniflandirici=sinifla)
        print('  %.0f sn' % (time.time() - t0))
        for h, r in res.items():
            if r.get('hata'):
                continue
            top[h]['tarandi'] += 1
            top[h]['urun'] += r.get('urun', 0)
            b = bek.get(h)
            for sz, n in (r.get('boy') or {}).items():
                if b is not None and abs(int(sz) - b) <= BOY_TOL:
                    top[h]['ayni_boyda'] += n
                else:
                    top[h]['boy_disi'] += n
            for k, v in (r.get('sinif') or {}).items():
                top[h]['sinif'][k] += v
        tamamlanan.append(d)
        kismi_yaz()
        print(u'     partial result written: %s' % os.path.relpath(kismi_yol, KOK))

    print()
    print('=' * 118)
    print('%-40s %9s | %11s | %9s %9s %9s %9s | %11s' %
          ('hedef', 'tum', 'BOY: disi', 'a klad-ici', 'ao organel',
           'b ayni-alan', 'c fark-alan', 'TAKSON: b+c'))
    print('-' * 118)
    degisen = []
    for c in ciftler:
        h = c['hedef']; t = top[h]
        if not t['tarandi']:
            continue
        s = t['sinif']
        kd = s['b'] + s['c']
        print('%-40s %9d | %11d | %9d %9d %9d %9d | %11d%s' %
              (h[:40], t['urun'], t['boy_disi'], s['a'], s['ao'], s['b'], s['c'], kd,
               '  (bilinmiyor %d)' % s['bilinmiyor'] if s['bilinmiyor'] else ''))
        if t['boy_disi'] != kd:
            degisen.append((h, t['boy_disi'], kd))
    print('-' * 118)

    print()
    print('=== FARK ===')
    if not degisen:
        print(u'The two criteria gave the same answer (on this subset of databases).')
    else:
        print('%-40s %14s %14s %s' % ('hedef', 'BOY olcutu', 'TAKSON olcutu', 'etki'))
        print('-' * 100)
        for h, b, k in degisen:
            if b > 0 and k == 0:
                etki = 'risky becomes clean (the length criterion was raising A FALSE ALARM)'
            elif b == 0 and k > 0:
                etki = 'clean becomes risky (the length criterion was MISSING a REAL cross reaction)'
            elif k < b:
                etki = 'the cross reaction count goes from %d to %d, a fall' % (b, k)
            else:
                etki = 'the cross reaction count goes from %d to %d, a rise' % (b, k)
            print('%-40s %14d %14d  %s' % (h[:40], b, k, etki))

    print()
    print(u'NOTE: this measurement covered only these databases: %s' % ', '.join(vtbler))
    print(u'     Full coverage needs SILVA/UNITE/PR2/ROD and takes HOURS.')
    print(u'     These numbers are NOT COVERAGE, they are evidence that the method runs.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
