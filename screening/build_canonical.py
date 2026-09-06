# -*- coding: utf-8 -*-
"""
build_canonical.py - it produces THE ONE CANONICAL SOURCE: canonical_consensus/

It reads every consensus (from the mixed orientation directories), turns them to the
SENSE direction with orientation.py and writes them into a SINGLE directory. From
then on every script reads HERE; no script has an orientation patch of its own.

The input directories (in order of precedence; if the same bin is in more than one
directory the one with precedence wins and the other is written to the manifest as
atlandi):
    1. SCREENING_RESULT/konsensus_yeni   (if present, the newest production)
    2. referans_konsensus/konsensus          (the set normalised overnight)
    3. consensus sequences                   (the original output)

The output:
    canonical_consensus/<bin>_kanonik.fasta
    canonical_consensus/MANIFEST.tsv           per file: source, old orientation, flipped
    canonical_consensus/UNDECIDED.tsv           files whose orientation could not be settled

Usage:
    python build_canonical.py --root ..
    python build_canonical.py --root .. --rerun        (overwrite if present)

"""
# -------------------------------------------------------------------------
# build_canonical.py - it scans every consensus directory, turns each bin's
#                   sequence to the SENSE direction with orientation.py and writes
#                   it into a single canonical directory; that way the orientation
#                   question is settled in one place.
#
# INPUT  : the three source directories under --root, in the order chosen with
#          --priority: consensus sequences (the original set the panel was built
#          on), SCREENING_RESULT/konsensus_yeni (the new production) and
#          referans_konsensus/konsensus. The orientation decision is made with
#          orientation.dosya_kanonik().
# OUTPUT : <bin>.canonical.fa per bin under canonical_consensus/; besides that
#          INDEX.tsv (the one list consumers must read), MANIFEST.tsv (the source,
#          the old orientation, whether it was flipped), UNDECIDED.tsv and
#          README.txt. Exit code 0 = no file needing a flip is left in the
#          canonical directory.
# CALLED BY: hepsi.kanonik_kos() runs it as a separate process, as the 2nd stage of
#          key 9 (--priority ozgun) and the 4th stage (--priority yeni). It is also
#          the command suggested for running by hand in the error message of every
#          stage whose orientation gate fails.
#
# WHY INDEX.tsv EXISTS: leftover old files on the mounted directory CANNOT BE
# DELETED. A consumer reading with glob collects those leftovers too and takes old
# mixed orientation files for canonical ones. That is why the list of valid files is
# kept separately.
# -------------------------------------------------------------------------
import os, sys, re, csv, glob, argparse, shutil

BURA = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BURA)
import orientation

# The precedence: if the same bin is in more than one directory the one first in
# the order wins.
#
# THE DEFAULT IS ozgun (the 2026-08-02 fix). The previous default was referans and
# it made A SILENT CHANGE OF SOURCE: the referans_konsensus/ directory is a
# DIFFERENT consensus rebuild from the consensus sequences/*_consensus_strict.fasta
# files the panel was built on (even the lengths differ, 1503 against 1534).
# Measured: with referans, Bakteri_universal (the UNIVERSAL bacterial pair) gave a
# product in only 2 of the 20 B bins and the length came out 135 instead of 130;
# with ozgun it gives a product in 7 bins and the panel's 130 bp value is
# reproduced. All the panel's numbers were measured on the ozgun set, so that must
# be the baseline. ORIENTATION normalisation is a separate job and is applied on
# both sources.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'steps'))
from select_consensus import DEFAULT_SETS as _SETS   # noqa: E402

ONCELIK = {
 'referans': [('referans_konsensus', 'referans_konsensus/konsensus'),
              ('konsensus_yeni', 'SCREENING_RESULT/konsensus_yeni'),
              ('ozgun', 'consensus sequences')],
 # yeni: used once the night's production is finished. THE FALLBACK ORDER MATTERS:
 # if konsensus_yeni could not produce a bin, ozgun (the panel's baseline) comes
 # FIRST; referans_konsensus was put last because it is a different rebuild (see the
 # note above).
 'yeni':     [('konsensus_yeni', 'SCREENING_RESULT/konsensus_yeni'),
              ('ozgun', 'consensus sequences'),
              ('referans_konsensus', 'referans_konsensus/konsensus')],
 'ozgun':    [('ozgun', 'consensus sequences'),
              ('konsensus_yeni', 'SCREENING_RESULT/konsensus_yeni'),
              ('referans_konsensus', 'referans_konsensus/konsensus')],
 # 'selection' (2026-09-06, external review): `./primerjury consensus` weighs every candidate
 # set by read support and writes CONSENSUS_SELECTION.tsv; until now nothing consumed that
 # table, so the winning polished consensus never reached identification or design. The
 # canonical set now follows the table: for a bin the table names, only the chosen set is
 # taken (its 'trusted' verdict is carried into the manifest, never hidden); a bin the table
 # does not name falls through the sets in this order. The set list is the one
 # select_consensus measures (single source).
 'selection': [(s, s) for s in _SETS],
}


_KUTU = re.compile(r'(?:^|[^A-Za-z0-9])(A1|A2|F1|F2|B)[-_](\d)(?!\d)')
_TAX = re.compile(r'(?<![0-9])(\d{3,7})(?![0-9])')
_BIN = re.compile(r'(BIN\d+|OBEK\d+)')     # a Kraken-free bin id (bin_reads); kept whole


def kutu_adi(yol):
    """extract the <class>-<barcode>_<taxid> bin from the file name.
    The naming of the source directories is INCONSISTENT
    (A1-1-reads_2209_consensus_strict, A1_1_reads_1826872_consensus_strict,
    A1-1_2209_yeniden_konsensus and so on), so pattern matching is used rather than
    splitting the name apart.

    """
    b = os.path.basename(yol).replace('.fasta', '')
    m = _KUTU.search(b)
    if not m:
        return None
    kutu = '%s-%s' % (m.group(1), m.group(2))
    kalan = b[m.end():]
    # BIN FIRST (2026-09-06, external review): the digits inside 'BIN123' were taken for a
    # taxid and the bin became 'B-1_123'; 'BIN1' was dropped altogether.
    adaylar = _BIN.findall(kalan) or _TAX.findall(kalan)
    if not adaylar:
        return None
    return '%s_%s' % (kutu, adaylar[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', dest='kok', required=True)
    ap.add_argument('--output', dest='cikti', default='canonical_consensus')
    ap.add_argument('--rerun', dest='yeniden', action='store_true')
    ap.add_argument('--priority', dest='oncelik', default='ozgun', choices=sorted(ONCELIK))
    a = ap.parse_args()

    h = orientation.kendini_sina()
    if h:
        sys.exit(u'orientation.py DID NOT PASS its self test: %s' % h)
    print(u'orientation.py self test: PASSED. The canonical orientation =', orientation.KANONIK_YON)
    print(u'source precedence :', ' > '.join(e for e, _ in ONCELIK[a.oncelik]))

    cik = os.path.join(a.kok, a.cikti)
    os.makedirs(cik, exist_ok=True)
    # NOTE: a file CANNOT BE DELETED on the mounted directory (Operation not
    # permitted). That is why valid files are written to the *.canonical.fa pattern and
    # recorded in INDEX.tsv. Consumers must read THE INDEX, NOT a glob; the old
    # *_kanonik.fasta leftovers are inert.

    manifest, belirsiz, gorulen = [], [], {}
    secim, guven = {}, {}
    if a.oncelik == 'selection':
        sy = os.path.join(a.kok, 'CONSENSUS_SELECTION.tsv')
        if not os.path.exists(sy):
            sys.exit('--priority selection needs %s (run ./primerjury consensus first)' % sy)
        with open(sy, encoding='utf-8') as fh:
            bas = fh.readline().rstrip('\n').split('\t')
            ik, ise, ig = bas.index('bin'), bas.index('CHOSEN set'), bas.index('trusted')
            for sat in fh:
                p = sat.rstrip('\n').split('\t')
                if len(p) > max(ik, ise, ig):
                    secim[p[ik]], guven[p[ik]] = p[ise], p[ig]
        print(u'selection table : %d bins (%d untrusted)' % (len(secim), sum(1 for v in guven.values() if v.startswith('NO'))))
    taninmayan = []
    for etiket, kl in ONCELIK[a.oncelik]:
        yollar = sorted(glob.glob(os.path.join(a.kok, kl, '**', '*.fasta'), recursive=True))
        for y in yollar:
            k = kutu_adi(y)
            if not k:
                taninmayan.append(os.path.relpath(y, a.kok))
                continue
            if secim and k in secim and secim[k] != etiket:
                continue          # the selection table chose another set for this bin
            sn = orientation.sinifi(os.path.basename(y))
            if sn == '?':
                sn = orientation.sinifi(y)
            if sn == '?':
                continue
            if k in gorulen:
                manifest.append(dict(kutu=k, sinif=sn, kaynak=etiket, dosya=os.path.relpath(y, a.kok),
                                     eski_yon='', cevrildi='', uzunluk='', durum='skipped (%s won)' % gorulen[k], secim_guven=guven.get(k, '')))
                continue
            kayitlar, _ = orientation.dosya_kanonik(y)
            kayitlar = [r for r in kayitlar if len(r[1]) >= 200]
            if not kayitlar:
                manifest.append(dict(kutu=k, sinif=sn, kaynak=etiket, dosya=os.path.relpath(y, a.kok),
                                     eski_yon='', cevrildi='', uzunluk=0, durum='empty or too short, skipped', secim_guven=guven.get(k, '')))
                continue
            ad, dizi, karar, cev = max(kayitlar, key=lambda r: len(r[1]))
            if karar == 'BELIRSIZ':
                belirsiz.append(dict(kutu=k, sinif=sn, kaynak=etiket,
                                     dosya=os.path.relpath(y, a.kok), uzunluk=len(dizi),
                                     N_yuzde=round(100.0 * dizi.count('N') / len(dizi), 1),
                                     not_='yon belirlenemedi - KANONIGE ALINMADI, maskeli'))
                manifest.append(dict(kutu=k, sinif=sn, kaynak=etiket, dosya=os.path.relpath(y, a.kok),
                                     eski_yon='BELIRSIZ', cevrildi='', uzunluk=len(dizi),
                                     durum='UNDECIDED, not written', secim_guven=guven.get(k, '')))
                continue
            gorulen[k] = etiket
            if k in guven and guven[k].startswith('NO'):
                print(u'  %-14s taken from %s although the selection marks it UNTRUSTED (%s)' % (k, etiket, guven[k]))
            cy = os.path.join(cik, '%s.canonical.fa' % k)
            with open(cy, 'w', encoding='utf-8') as fh:
                fh.write('>%s kanonik=%s kaynak=%s eski_yon=%s cevrildi=%s\n'
                         % (k, orientation.KANONIK_YON, etiket, karar, int(cev)))
                for i in range(0, len(dizi), 70):
                    fh.write(dizi[i:i + 70] + '\n')
            manifest.append(dict(kutu=k, sinif=sn, kaynak=etiket, dosya=os.path.relpath(y, a.kok),
                                 eski_yon=karar, cevrildi='EVET' if cev else 'hayir',
                                 uzunluk=len(dizi), durum='yazildi', secim_guven=guven.get(k, '')))

    def yaz(ad, rows):
        if not rows:
            open(os.path.join(cik, ad), 'w', encoding='utf-8').write('(bos)\n')
            return
        with open(os.path.join(cik, ad), 'w', newline='', encoding='utf-8') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t')
            w.writeheader()
            for r in rows:
                w.writerow(r)
    yaz('MANIFEST.tsv', manifest)
    yaz('UNDECIDED.tsv', belirsiz)
    yaz('INDEX.tsv', [dict(kutu=m['kutu'], sinif=m['sinif'],
                            dosya='%s.canonical.fa' % m['kutu'], kaynak=m['kaynak'],
                            eski_yon=m['eski_yon'], cevrildi=m['cevrildi'],
                            uzunluk=m['uzunluk'])
                       for m in manifest if m['durum'] == 'yazildi'])
    open(os.path.join(cik, 'README.txt'), 'w', encoding='utf-8').write(
        u'CANONICAL CONSENSUS DIRECTORY\nValid files: *.canonical.fa  (listed in INDEX.tsv)\nThis directory can hold LEFTOVER files ending in *_kanonik.fasta. Those are the\nmisnamed output of the first run and could not be deleted on a mounted drive. IGNORE THEM.\nEvery consumer must read INDEX.tsv and must NOT use glob.\nCanonical orientation: SENSE. Definition and criterion: screening/orientation.py\n')

    yazilan = [m for m in manifest if m['durum'] == 'yazildi']
    cevrilen = [m for m in yazilan if m['cevrildi'] == 'EVET']
    print(u'\ncanonical directory : %s' % cik)
    print(u'bins written        : %d' % len(yazilan))
    print(u'  converted           : %d (ANTISENSE -> SENSE)' % len(cevrilen))
    print(u'  already sense       : %d' % (len(yazilan) - len(cevrilen)))
    print(u'BELIRSIZ       : %d (not written, UNDECIDED.tsv)' % len(belirsiz))
    if taninmayan:
        print(u'unrecognised file names (no bin id found, skipped): %d' % len(taninmayan))
        for t in taninmayan[:10]:
            print(u'    %s' % t)
    kay = {}
    for m in yazilan:
        kay[m['kaynak']] = kay.get(m['kaynak'], 0) + 1
    print(u'source distribution:', ', '.join('%s=%d' % x for x in sorted(kay.items())))

    # VERIFICATION: is every file written really SENSE
    kotu = 0
    for y in sorted(glob.glob(os.path.join(cik, '*.canonical.fa'))):
        kayitlar, sn = orientation.dosya_kanonik(y)
        for ad, dizi, karar, cev in kayitlar:
            if cev:
                kotu += 1
    print(u'\nVERIFICATION: files in the canonical directory that still need flipping =', kotu,
          '(0 olmali)' if kotu == 0 else '*** SORUN ***')
    return 0 if kotu == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
