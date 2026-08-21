# -*- coding: utf-8 -*-
"""Repository health check, catches the gaps that syntax checking cannot.

WHY THIS EXISTS
---------------
While reorganising this repository, `screening/motor.py` passed every syntax
check and still could not be imported: it loads its sequence engine from disk
at import time, and those modules had not been copied. A parse-only check said
"110 files, 0 errors" about a package that did not work.

That is the failure mode this file guards against. Each check below answers a
question that `ast.parse` cannot:

  1  IMPORT       do the packages actually import, side effects and all?
  2  ENTRY POINTS do scripts referenced by other scripts exist?
  3  STALE NAMES  is anything still pointing at a pre-rename path?
  4  PRIVACY      did any personal or study-specific name survive?
  5  CONFIG       do the paths in the config point at things that can exist?
  6  PACKAGING    are the files git would ship actually sufficient?

RUN
---
    python3 tests/test_repo_health.py
    python3 tests/test_repo_health.py --ayrinti      # list every finding

Exit code 0 only if every check passes.
"""
from __future__ import print_function

import argparse
import ast
import importlib
import io
import os
import re
import subprocess
import sys

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, KOK)

ATLA_DIZIN = {'.git', '__pycache__', 'docs', 'sequences', 'examples'}

# Packages that must import cleanly, side effects included.
PAKETLER = [
    'screening.config', 'screening.taxonomy', 'screening.engine_gateway',
    'screening.targets', 'screening.global_scan', 'screening.geometry',
    'screening.order_classes', 'screening.sample', 'screening.reference',
    'screening.checks', 'screening.report', 'screening.orientation',
    'verification.mfeprimer_layer',
]

# Names that must not survive anywhere in shipped code.
#
# WARNING TO FUTURE MAINTAINERS: this list contains the very strings a bulk
# rename would target. A search-and-replace run over the repository WILL
# corrupt it, that happened once already: the pattern for the old project
# name was rewritten to the new one, which would have made this check reject
# the correct name. Exclude this file from any bulk rewrite, and re-read the
# list afterwards.
YASAK = [
    (r'\bMicRho' + r'Booster\b', 'old project name'),
    (r'\bKAPSAMLI_ARAMA\b', 'old package name (now screening)'),
    (r'\bKURTARMA\b(?!_SONUC)', 'old package name (now verification)'),
    (r'\bWSL_betikleri\b', 'old directory name (now steps)'),
    (r'\bTEK_PROTOKOL\b(?!_SONUC)', 'old directory name (now protocol)'),
    (r'\bORTAK_PUANLAYICI\b', 'old directory name (now scoring)'),
    (r'\bDENETIM_SINAMALARI\b', 'old directory name (now tests)'),
    (r'\bARACLAR\b', 'old directory name (now tools)'),
    (r'\bsekanslar\b', 'old input directory (now sequences)'),
    (r'\bornek_veri\b', 'old directory name (now examples)'),
    (r'\bKUR\.sh\b', 'old entry point (now install.sh)'),
    (r'\bINDEKS_KUR\.sh\b', 'old entry point (now build_index.sh)'),
    (r'\bcapraz_kontrol\.py\b', 'old entry point (now cross_check.py)'),
    # Personal names, split so a bulk rename cannot rewrite them (see warning).
    (r'\bBur' + r'ak\b', 'personal name'),
    # OLCULDU: bu desen once yalniz "Al" + "i'nin" idi. Kaynakta kesme
    # isareti KACISLI durdugu icin (Ali'nin) hicbir zaman eslesmedi. Ustelik
    # ciplak ALI ve Ali gecisleri, bir CLI bayragi (--ali) ve bir ortam
    # degiskeni olarak depoda duruyordu. Desen artik ada bakiyor.
    (r'\bAl' + r'i\b', 'personal name'),
    (r'\bAL' + r'I\b', 'personal name'),
    (r'--al' + r'i\b', 'personal name as a CLI flag'),
    (r'\bhoc' + r'aya\b', 'study-specific role'),
    (r'\bhoc' + r'anin\b', 'study-specific role'),
    (r'\bSON_ETAP_betikleri\b', 'directory that is not in this repository'),
    (r'\bMADDE123_betikleri\b', 'directory that is not in this repository'),
    (r'\bDUZELTME_betikleri\b', 'directory that is not in this repository'),
    (r'^\s*\d\d_[a-z]', 'numbered script name'),
]

# Old numbered script names must be gone from references too.
NUMARALI = re.compile(r'\b\d{2,3}_[A-Za-z_]+\.(py|sh)\b')


def metin_dosyalari():
    for kok, dizinler, dosyalar in os.walk(KOK):
        dizinler[:] = [d for d in dizinler if d not in ATLA_DIZIN]
        for d in dosyalar:
            if d.endswith(('.py', '.sh', '.md', '.txt', '.cff', '.tsv')):
                yield os.path.join(kok, d)


def oku(y):
    try:
        return io.open(y, encoding='utf-8', errors='replace').read()
    except Exception:
        return ''


def rel(y):
    return os.path.relpath(y, KOK).replace(os.sep, '/')


# ---------------------------------------------------------------- 1 IMPORT
def kontrol_import(bulgu):
    for m in PAKETLER:
        try:
            importlib.import_module(m)
        except Exception as e:
            bulgu.append(('IMPORT', m, '%s: %s' % (type(e).__name__, str(e)[:110])))


# ----------------------------------------------------------- 2 ENTRY POINTS
def kontrol_giris_noktalari(bulgu):
    """Scripts referenced by other scripts must exist.

    Only literal, repo-relative references are checked; a false alarm here is
    worse than a miss, because it trains people to ignore the report.
    """
    var = set()
    for kok, dizinler, dosyalar in os.walk(KOK):
        dizinler[:] = [d for d in dizinler if d not in ATLA_DIZIN]
        for d in dosyalar:
            var.add(d)

    desen = re.compile(r'[\'"]([A-Za-z0-9_./-]+\.(?:py|sh))[\'"]')
    for y in metin_dosyalari():
        s = oku(y)
        for m in desen.finditer(s):
            ref = m.group(1)
            ad = os.path.basename(ref)
            # skip obvious non-repo references
            if ref.startswith(('/', 'http')) or ad in ('setup.py',):
                continue
            if '/' in ref and not ref.startswith('.'):
                tam = os.path.join(KOK, ref)
                if os.path.exists(tam):
                    continue
            if ad in var:
                continue
            bulgu.append(('ENTRY', rel(y), 'references missing file: %s' % ref))


# ------------------------------------------------------------ 3+4 STALE/PRIVACY
# Attribution files: the author's name BELONGS here. The rule being enforced
# is "no personal names in the CODE", not "no author anywhere" -- a project
# with no author is not more neutral, only less citable.
ATIF_DOSYALARI = {'CITATION.cff', 'LICENSE', 'README.md'}


def kontrol_yasak_adlar(bulgu):
    for y in metin_dosyalari():
        if rel(y) == 'tests/test_repo_health.py':
            continue                      # this file lists them on purpose
        if os.path.basename(y) in ATIF_DOSYALARI:
            continue                      # author attribution is expected here
        s = oku(y)
        for desen, sebep in YASAK:
            for m in re.finditer(desen, s, re.M):
                satir = s[:m.start()].count('\n') + 1
                bulgu.append(('NAME', '%s:%d' % (rel(y), satir),
                              '%s -> %s' % (m.group(0).strip(), sebep)))
        for m in NUMARALI.finditer(s):
            satir = s[:m.start()].count('\n') + 1
            bulgu.append(('NAME', '%s:%d' % (rel(y), satir),
                          'numbered script reference: %s' % m.group(0)))


# ---------------------------------------------------------------- 5 CONFIG
def kontrol_yapilandirma(bulgu):
    """Config paths must be inside the repo or be documented run outputs."""
    try:
        from screening import config as C
    except Exception as e:
        bulgu.append(('CONFIG', 'screening/config.py',
                      'cannot import: %s' % type(e).__name__))
        return
    # Directories the config points at that the code loads AT IMPORT time
    for ad in ('BETIK_YOLLARI',):
        for p in getattr(C, ad, []):
            if not os.path.isdir(p):
                bulgu.append(('CONFIG', ad,
                              'directory does not exist: %s' % rel(p)))
    # Sanity: the engine modules motor.py needs
    gerekli = ['ispcr.py']
    for g in gerekli:
        if not any(os.path.exists(os.path.join(p, g))
                   for p in getattr(C, 'BETIK_YOLLARI', [])):
            bulgu.append(('CONFIG', 'BETIK_YOLLARI',
                          'required engine module not found: %s' % g))


# ------------------------------------------------------------- 6 BARE IMPORTS
def kontrol_ciplak_import(bulgu):
    """Bare `import X` must name a module that exists.

    Several scripts add their own directory to sys.path and import siblings
    by bare name. Those imports carry no package prefix and no quotes, so a
    rename pass that rewrites package imports and filename strings misses
    them entirely -- measured: 24 bare imports still named modules that had
    been renamed, and every one would have failed at runtime with
    ModuleNotFoundError while every syntax check passed.
    """
    var = set()
    for kok, dizinler, dosyalar in os.walk(KOK):
        dizinler[:] = [d for d in dizinler if d not in ATLA_DIZIN]
        for d in dosyalar:
            if d.endswith('.py'):
                var.add(d[:-3])
    import sys as _sys
    std = getattr(_sys, "stdlib_module_names", set())
    DIS = {'numpy', 'Bio', 'primer3', 'openpyxl', 'pysam', 'matplotlib',
           'mappy', 'requests', 'yaml', 'scipy', 'pandas', 'playwright',
           'setuptools', 'pkg_resources'}
    # MEASURED TWICE: a regex over the raw text reads prose as an import.
    # An English sentence wrapped so that a line begins with "from" ("...different
    # enough / from one another?") matched `^\s*from\s+(\w+)\s` and was reported as
    # a bare import of the module "one". The imports are now read from the AST, so
    # only real import statements are seen and comments and docstrings cannot
    # produce a finding.
    PAKETLER_KOK = ('screening', 'verification', 'steps', 'protocol',
                    'scoring', 'engine', 'tools', 'tests')
    for y in metin_dosyalari():
        if not y.endswith('.py'):
            continue
        try:
            agac = ast.parse(oku(y))
        except SyntaxError:
            continue                  # the PARSE check reports this separately
        for n in ast.walk(agac):
            if isinstance(n, ast.Import):
                adlar = [(a.name.split('.')[0], n.lineno) for a in n.names]
            elif isinstance(n, ast.ImportFrom):
                if n.level:            # a relative import names no top level module
                    continue
                adlar = [((n.module or '').split('.')[0], n.lineno)]
            else:
                continue
            for ad, satir in adlar:
                if not ad or ad in var or ad in std or ad in DIS:
                    continue
                if ad in PAKETLER_KOK:
                    continue
                bulgu.append(('IMPORT', '%s:%d' % (rel(y), satir),
                              'bare import of unknown module: %s' % ad))


# ------------------------------------------------------------ 6 LINE ENDINGS
def kontrol_satir_sonu(bulgu):
    """Shell scripts and Python files must use LF, not CRLF.

    This project runs on Linux and WSL2. A CRLF shell script fails there in a
    way that names characters instead of causes: bash reports a carriage
    return as an unknown command, or a syntax error near an unexpected token.

    Measured once on this repository, after a Windows checkout with git's
    default autocrlf: every shell script broke in WSL, while the same files
    read fine on Windows. .gitattributes now forces LF; this check makes sure
    it stays that way.
    """
    CRLF = bytes((13, 10))
    for y in metin_dosyalari():
        try:
            ham = open(y, 'rb').read()
        except Exception:
            continue
        if CRLF in ham:
            bulgu.append(('EOL', rel(y), 'CRLF line endings; must be LF'))
    # the executable entry point is not caught by metin_dosyalari()
    g = os.path.join(KOK, 'primerjury')
    if os.path.exists(g) and CRLF in open(g, 'rb').read():
        bulgu.append(('EOL', 'primerjury', 'CRLF line endings; must be LF'))


# ------------------------------------------------------------- 7 PACKAGING
def kontrol_paketleme(bulgu):
    """What git ships must be enough to import the packages."""
    try:
        cikti = subprocess.check_output(['git', 'ls-files'], cwd=KOK)
    except Exception:
        return                                    # not a git repo; skip quietly
    izlenen = set(cikti.decode('utf-8', 'replace').split())
    # every .py inside a package directory must be tracked
    for kok, dizinler, dosyalar in os.walk(KOK):
        dizinler[:] = [d for d in dizinler if d not in ATLA_DIZIN]
        for d in dosyalar:
            if not d.endswith('.py'):
                continue
            r = rel(os.path.join(kok, d))
            if r not in izlenen:
                bulgu.append(('PACKAGING', r, 'python file not tracked by git'))
    # required top-level files
    for zorunlu in ('README.md', 'LICENSE', 'requirements.txt', '.gitignore',
                    'install.sh', 'docs/GUIDE.md'):
        if zorunlu not in izlenen:
            bulgu.append(('PACKAGING', zorunlu, 'required file missing from git'))


# ------------------------------------------------------------------ 7 PARSE
def kontrol_parse(bulgu):
    for y in metin_dosyalari():
        if not y.endswith('.py'):
            continue
        try:
            ast.parse(oku(y))
        except SyntaxError as e:
            bulgu.append(('PARSE', '%s:%s' % (rel(y), e.lineno), e.msg))


KONTROLLER = [
    ('PARSE',     'every Python file parses',            kontrol_parse),
    ('IMPORT',    'core packages import (side effects)', kontrol_import),
    ('ENTRY',     'referenced scripts exist',            kontrol_giris_noktalari),
    ('NAME',      'no stale or personal names',          kontrol_yasak_adlar),
    ('CONFIG',    'config paths resolve',                kontrol_yapilandirma),
    ('IMPORT',    'bare imports resolve',                kontrol_ciplak_import),
    ('EOL',       'text files use LF, not CRLF',         kontrol_satir_sonu),
    ('PACKAGING', 'git ships enough to run',             kontrol_paketleme),
]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--detail', '--ayrinti', dest='ayrinti', action='store_true', help='list every finding')
    a = p.parse_args()

    bulgu = []
    print('%-12s %-38s %s' % ('CHECK', 'QUESTION', 'RESULT'))
    print('-' * 78)
    for etiket, soru, fn in KONTROLLER:
        onceki = len(bulgu)
        fn(bulgu)
        yeni = len(bulgu) - onceki
        print('%-12s %-38s %s' % (etiket, soru,
                                  'PASS' if yeni == 0 else 'FAIL (%d)' % yeni))

    if bulgu:
        print()
        gruplu = {}
        for tur, nerede, ne in bulgu:
            gruplu.setdefault(tur, []).append((nerede, ne))
        for tur in sorted(gruplu):
            liste = gruplu[tur]
            print('=== %s (%d) ===' % (tur, len(liste)))
            for nerede, ne in (liste if a.ayrinti else liste[:12]):
                print('   %-46s %s' % (nerede[:46], ne[:80]))
            if not a.ayrinti and len(liste) > 12:
                print('   ... %d more (use --ayrinti)' % (len(liste) - 12))
            print()

    print('-' * 78)
    print('RESULT: %s' % ('ALL CHECKS PASSED' if not bulgu
                          else '%d FINDINGS' % len(bulgu)))
    return 0 if not bulgu else 1


if __name__ == '__main__':
    sys.exit(main())
