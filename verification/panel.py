# -*- coding: utf-8 -*-
"""panel.py, THE ONE ENTRY POINT.

WHY THEY WERE MERGED
--------------------
Eight separate scripts piled up in one night: the audit, the geometry, the plate,
the thresholds, the current state, the state of the requests, the NCBI layer and
bin generation. Eight separate commands means running one and forgetting another,
and every fault in this project came out exactly that way: something changed and
what depended on it was not run again.

Here they are all subcommands of ONE command. The code WAS NOT COPIED: every
subcommand calls its own script, so two copies cannot drift apart over time.

Usage:
    python verification/panel.py hepsi           # all of them in order (measures, writes nothing)
    python verification/panel.py denetle         # the audit gate
    python verification/panel.py geometri        # Tm and rule checks with primer3
    python verification/panel.py geometri --write  # correct the Tm and dTm columns
    python verification/panel.py plaka           # a plate suggestion against gel overlap
    python verification/panel.py esik            # the two threshold rules side by side
    python verification/panel.py guncel          # produce the current state document
    python verification/panel.py toplanti        # what was asked for and what came of it
    python verification/panel.py ncbi4           # the NCBI fourth layer, about 30 minutes
    python verification/panel.py kutu --plan     # a plan for the taxa never measured
    python verification/panel.py referans        # refresh the quick test references

"hepsi" measures and changes NO file: the steps with a side effect, such as
geometri --write and ncbi4, ARE NOT INCLUDED in it. A person asks for a side
effect.
"""
from __future__ import print_function

import os
import subprocess
import sys

BURA = os.path.dirname(os.path.abspath(__file__))
KOK = os.path.dirname(BURA)

# alt komut -> (betik yolu, aciklama, "hepsi"ye dahil mi)
KOMUT = {
    'denetle':  ('verification/audit_all.py',
                 'the eighteen item audit gate: the tables, the references, '
                 'the seals, the geometry, the plate and the document counts', True),
    'guncel':   ('verification/current_status.py',
                 u'GUNCEL_DURUM.md - panelin bugunku sayilari', True),
    'toplanti': ('verification/decision_status.py',
                 'what was asked for and which of it came about', True),
    'esik':     ('verification/recompute_thresholds.py',
                 'the flat threshold beside the abundance weighted one', True),
    'geometri': ('verification/refresh_geometry.py',
                 'Tm, GC, length, hairpin and dimer with primer3, plus the '
                 'product length', True),
    'plaka':    ('verification/assign_plate.py',
                 u'plaka ici jel cakismasini azaltan dagilim onerisi', True),
    'referans': ('verification/refresh_reference.py',
                 'refresh the quick test references from a full run', False),
    'ncbi4':    ('verification/ncbi_layer.py',
                 'the NCBI Primer-BLAST fourth layer, about 30 minutes and it '
                 'uses the network', False),
    'kutu':     ('verification/build_bins.py',
                 'turn the taxa never measured into bins, behind the '
                 'calibration gate', False),
    'siparis':  ('verification/order_form.py',
                 'the ONE correct oligo list for the supplier, produced '
                 'rather than written', True),
    'excel':    ('verification/build_excel.py',
                 'write every current value into ONE Excel file, produced '
                 'rather than written', True),
    'lokus':    ('engine/target_full.py',
                 'look for a pair at EVERY locus of one target; primer3 is '
                 'needed', False),
    'arsiv':    ('verification/archive.py',
                 'move the old files aside, printing the PLAN first', False),
    'sinif':    ('verification/ncbi_reclassify.py',
                 'count the NCBI results again under the strict name rule, '
                 'without the network', True),
    'kapsama':  ('screening/exclusion_coverage_check.py',
                 'does the exclusion taxon cover its members, against NCBI '
                 'Taxonomy', False),
}

# "hepsi" sirasi: once durum uretilir, sonra denetim onu da gorur.
# "sabah" = hepsi + NCBI yeniden sayimi. Ikisi de OLCER, DEGISTIRMEZ.
SIRA = ['guncel', 'siparis', 'toplanti', 'esik', 'geometri', 'plaka',
        'sinif', 'excel', 'denetle']


def kos(ad, ek):
    yol = os.path.join(KOK, *KOMUT[ad][0].split('/'))
    if not os.path.exists(yol):
        print(u'  SKIPPED: %s is missing (%s)' % (ad, yol))
        return 127
    # target_full.py --root almaz; kokU _FL_KOK ortam degiskeninden okur.
    if ad == 'lokus':
        os.environ['_FL_KOK'] = KOK
        komut = [sys.executable, yol] + list(ek)
    else:
        komut = [sys.executable, yol, '--root', KOK] + list(ek)
    print()
    print('=' * 78)
    print('  >> %s   (%s)' % (ad, KOMUT[ad][1]))
    print('  $ %s' % ' '.join(komut[1:]))
    print('=' * 78)
    sys.stdout.flush()
    return subprocess.call(komut)


def yardim():
    print(__doc__)
    print('  alt komutlar:')
    for k in sorted(KOMUT):
        print('    %-10s %s%s' % (k, KOMUT[k][1],
                                  '' if KOMUT[k][2] else u'   [not included in \'all\']'))


def main(argv):
    if not argv or argv[0] in ('-h', '--help', 'yardim'):
        yardim()
        return 0
    ad, ek = argv[0], argv[1:]
    if ad in ('hepsi', 'sabah'):
        kodlar = {}
        for k in SIRA:
            kodlar[k] = kos(k, [])
        print()
        print('=' * 78)
        print(u'  ALL DONE')
        for k in SIRA:
            print('    %-10s cikis kodu %s' % (k, kodlar[k]))
        print()
        print('  The steps with a side effect were DELIBERATELY not run. To '
              'run them separately:')
        for k in sorted(KOMUT):
            if not KOMUT[k][2]:
                print('    python verification/panel.py %s' % k)
        print('=' * 78)
        return 1 if any(v not in (0, 2) for v in kodlar.values()) else 0
    if ad not in KOMUT:
        print('Unknown subcommand: %s' % ad)
        yardim()
        return 2
    return kos(ad, ek)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
