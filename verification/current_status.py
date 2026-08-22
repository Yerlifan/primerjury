# -*- coding: utf-8 -*-
"""GUNCEL_DURUM.md  -  panelin BUGUNKU sayilari, TEK yerden ve URETILEREK.

NEDEN VAR
---------
Ayni soruya uc belge uc farkli cevap veriyordu: OKU_ONCE "KESIN 16 cift",
NASIL_DEVAM_EDILIR "KESIN 16 cift", CALISTIRMA_KILAVUZU "siparis edilebilir
11 cift". Panelin bugunku hali 20 cift (15 hedef ozgul + 5 evrensel). Ucu de
yazildiklari gunde dogruydu; sorun sayinin belgeye ELLE yazilmis olmasi.

Bu dosya elle yazilmaz. Her kosuda panelin kendi ciktisindan uretilir. Belgeler
sayiyi tekrar etmek yerine buraya isaret eder.

Kosum:
    python verification/current_status.py --root .
"""
from __future__ import print_function

import argparse
import csv
import io
import os
import sys
import time


def _tsv(yol):
    if not os.path.exists(yol):
        return []
    with io.open(yol, encoding='utf-8') as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith('#')),
                                   delimiter='\t'))


def sayilar(kok):
    sl = _tsv(os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'SIPARIS_LISTESI.tsv'))
    kesin = [r for r in sl if (r.get('SINIF') or '').strip().upper() == 'KESIN']
    evr = [r for r in sl if (r.get('SINIF') or '').strip().upper() == 'EVRENSEL']
    disi = [r for r in sl if (r.get('SINIF') or '').strip().upper()
            not in ('KESIN', 'EVRENSEL')]
    return sl, kesin, evr, disi


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', dest='kok', default='.')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)
    sl, kesin, evr, disi = sayilar(kok)
    if not sl:
        sys.exit(u'ERROR: TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv could not be read.')

    y = os.path.join(kok, 'GUNCEL_DURUM.md')
    with io.open(y, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# The panel as it stands today\n\n')
        fh.write(u'Generated: %s, source `TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv`\n\n'
                 % time.strftime('%Y-%m-%d %H:%M'))
        fh.write(u'> This file is **never written by hand**. It is regenerated from the panel\'s own output on every run. If another document contradicts it')
        fh.write(u'| | pairs | oligos |\n|---|---|---|\n')
        fh.write(u'| **To order** | **%d** | **%d** |\n'
                 % (len(kesin) + len(evr), 2 * (len(kesin) + len(evr))))
        fh.write(u'| of which target-specific | %d | %d |\n' % (len(kesin), 2 * len(kesin)))
        fh.write(u'| of which universal/control | %d | %d |\n' % (len(evr), 2 * len(evr)))
        fh.write(u'| Not ordered (not deleted) | %d | - |\n\n' % len(disi))
        fh.write(u'## Pairs to order\n\n')
        fh.write(u'| target | class | condition | product (bp) | dCq |\n|---|---|---|---|---|\n')
        for r in kesin + evr:
            fh.write(u'| %s | %s | %s | %s | %s |\n'
                     % (r.get('hedef', ''), (r.get('SINIF') or '').strip(),
                        (r.get('siparis_sarti') or '—').strip(),
                        (r.get('urun_bp') or '—').strip(),
                        (r.get('dCq_karsiligi') or '—').strip()))
        if disi:
            fh.write(u'\n## Not ordered\n\n| target | status | dCq |\n|---|---|---|\n')
            for r in disi:
                fh.write(u'| %s | %s | %s |\n'
                         % (r.get('hedef', ''), (r.get('durum') or '—').strip(),
                            (r.get('dCq_karsiligi') or '—').strip()))
        fh.write(u'\n## What this number does not know\n\n')
        fh.write(u'This table is a verdict on **threshold and coverage**. A pair appearing here does not mean it passed the geometry gate, that in-plate gel separation is clean')
    print('yazildi: %s' % y)
    print(u'  to order: %d pairs (%d target-specific + %d universal), %d oligos'
          % (len(kesin) + len(evr), len(kesin), len(evr), 2 * (len(kesin) + len(evr))))
    print(u'  not ordered: %d' % len(disi))
    return 0


if __name__ == '__main__':
    sys.exit(main())
