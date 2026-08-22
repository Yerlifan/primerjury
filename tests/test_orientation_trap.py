# -*- coding: utf-8 -*-
"THE EVIDENCE for the mixed orientation directory check.\n\nWHAT IS TESTED\n--------------\nThe cross-check's mixed directory scan has to answer this correctly: does this\nfile read the mixed orientation consensus directory IN ITS CODE?\n\nThe question matters because that directory is mixed in orientation (measured:\n71 antisense against 27 sense), and on a reversed consensus an in-silico PCR\nSILENTLY gives 0 products, a measured loss of 100 per cent.\n\nWHY A SEPARATE TEST IS NEEDED\n-----------------------------\nThe check used to search the WHOLE file as plain text. In one run it produced\nfive false positives: files that merely MENTIONED the directory name in their\nprose were marked risky too. A false positive from an auditor is not harmless; it\nburies the real findings in noise and costs the report its credibility.\n\nThe fix added three filters:\n  1) docstrings and comments are stripped, so a name in prose is not code\n  2) the arguments of print calls are stripped, so a message on the screen is not\n     code\n  3) the readers whose job it is are exempt\n\nTO RUN IT\n---------\n    python3 tests/test_orientation_trap.py\n    python3 tests/test_orientation_trap.py --root /another/path\n\nExit code 0 means all seven tests came out as expected.\nTHE PATH IS NOT EMBEDDED: by default the directory above this file is the root.\n"
from __future__ import print_function
import argparse
import importlib.util
import io
import os
import sys

VARSAYILAN_KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

KARISIK = u'consensus sequences'
KANONIK = u'konsensus_kanonik'

# --- sentetik ornekler: her biri TEK bir ayrimi sinar --------------------
GERCEK_RISK = u'''
import os
def yukle(kok):
    return open(os.path.join(kok, "consensus sequences", "A1-1.fasta")).read()
'''

SADECE_ACIKLAMA = '\n"""This script used to read the "consensus sequences" directory."""\n# consensus sequences is no longer used\nimport os\ndef yukle(kok):\n    return open(os.path.join(kok, "kanonik", "A1-1.fasta")).read()\n'

SADECE_MESAJ = '\nimport os\ndef kos(yaz, kok):\n    yaz("  The source: \\"consensus sequences\\" (the old set).")\n    return open(os.path.join(kok, "kanonik", "A1-1.fasta")).read()\n'

# --- gercek dosyalar: 2026-08-09'da yanlis pozitif verenler --------------
GERCEK_DOSYALAR = [
    (os.path.join('screening', 'run_all.py'), False,
     'a message string that is only printed to the screen'),
    (os.path.join('screening', 'orientation_audit.py'), False,
     'it reads the directory because that is its job, so it is exempt'),
    (os.path.join('steps', 'generate_primer_candidates.py'), False,
     u'yol CLI argumani, docstring ornegi'),
    (os.path.join('steps', 'design_group_primers.py'), False,
     u'yol CLI argumani, docstring ornegi'),
]


def ck_yukle(kok):
    yol = os.path.join(kok, 'cross_check.py')
    if not os.path.exists(yol):
        sys.stderr.write(u'cross_check.py was not found: %s\n' % yol)
        sys.exit(2)
    spec = importlib.util.spec_from_file_location('ck', yol)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', dest='kok', default=VARSAYILAN_KOK)
    a = p.parse_args()
    ck = ck_yukle(a.kok)

    def riskli(metin, ad=u''):
        v = ck.d9_karisik_klasor_yollari(metin, ad)
        return bool(v) and KANONIK not in ck._kod_govdesi(metin)

    def duz_metin(metin):          # duzeltme ONCESI mantik, karsilastirma icin
        return KARISIK in metin and KANONIK not in metin

    gecti = True
    print(u'%-42s | %-8s | %-8s | %-8s | %s'
          % (u'sinama', u'ESKI', u'YENI', u'beklenen', u'sonuc'))
    print(u'-' * 100)

    for rel, beklenen, gerekce in GERCEK_DOSYALAR:
        yol = os.path.join(a.kok, rel)
        if not os.path.exists(yol):
            print(u'%-42s | NO SUCH FILE - the test could not be run' % rel)
            gecti = False
            continue
        m = io.open(yol, encoding='utf-8', errors='replace').read()
        y = riskli(m, yol)
        ok = (y == beklenen)
        gecti = gecti and ok
        print(u'%-42s | %-8s | %-8s | %-8s | %-6s (%s)'
              % (rel.replace(os.sep, '/'),
                 u'RISKLI' if duz_metin(m) else u'temiz',
                 u'RISKLI' if y else u'temiz',
                 u'RISKLI' if beklenen else u'temiz',
                 u'DOGRU' if ok else u'YANLIS', gerekce))

    print()
    for ad, metin, beklenen in [
            (u'sentetik: KODDA gercekten okuyor', GERCEK_RISK, True),
            ('synthetic: it only mentions it in prose', SADECE_ACIKLAMA, False),
            ('synthetic: it only prints it to the screen', SADECE_MESAJ, False)]:
        y = riskli(metin, u'<synthetic>')
        ok = (y == beklenen)
        gecti = gecti and ok
        print(u'%-42s | %-8s | %-8s | %-8s | %s'
              % (ad, u'-', u'risky' if y else u'clean',
                 u'risky' if beklenen else u'clean',
                 u'RIGHT' if ok else u'WRONG'))

    print()
    print(u'RESULT: ' + (u'ALL SEVEN TESTS CAME OUT AS EXPECTED'
                        if gecti else u'FAILED'))
    return 0 if gecti else 1


if __name__ == '__main__':
    sys.exit(main())
