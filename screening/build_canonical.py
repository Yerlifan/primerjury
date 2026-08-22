# -*- coding: utf-8 -*-
"""
build_canonical.py - it produces THE ONE CANONICAL SOURCE: konsensus_kanonik/

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
    konsensus_kanonik/<bin>_kanonik.fasta
    konsensus_kanonik/MANIFEST.tsv           per file: source, old orientation, flipped
    konsensus_kanonik/BELIRSIZ.tsv           files whose orientation could not be settled

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
# OUTPUT : <bin>.kanonik.fa per bin under konsensus_kanonik/; besides that
#          INDEKS.tsv (the one list consumers must read), MANIFEST.tsv (the source,
#          the old orientation, whether it was flipped), BELIRSIZ.tsv and
#          OKUBENI.txt. Exit code 0 = no file needing a flip is left in the
#          canonical directory.
# CALLED BY: hepsi.kanonik_kos() runs it as a separate process, as the 2nd stage of
#          key 9 (--priority ozgun) and the 4th stage (--priority yeni). It is also
#          the command suggested for running by hand in the error message of every
#          stage whose orientation gate fails.
#
# WHY INDEKS.tsv EXISTS: leftover old files on the mounted directory CANNOT BE
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
}


_KUTU = re.compile(r'(?:^|[^A-Za-z0-9])(A1|A2|F1|F2|B)[-_](\d)(?!\d)')
_TAX = re.compile(r'(?<![0-9])(\d{3,7})(?![0-9])')


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
    adaylar = [x for x in _TAX.findall(kalan)]
    if not adaylar:
        return None
    return '%s_%s' % (kutu, adaylar[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', dest='kok', required=True)
    ap.add_argument('--output', dest='cikti', default='konsensus_kanonik')
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
    # permitted). That is why valid files are written to the *.kanonik.fa pattern and
    # recorded in INDEKS.tsv. Consumers must read THE INDEX, NOT a glob; the old
    # *_kanonik.fasta leftovers are inert.

    manifest, belirsiz, gorulen = [], [], {}
    for etiket, kl in ONCELIK[a.oncelik]:
        yollar = sorted(glob.glob(os.path.join(a.kok, kl, '**', '*.fasta'), recursive=True))
        for y in yollar:
            k = kutu_adi(y)
            if not k:
                continue
            sn = orientation.sinifi(os.path.basename(y))
            if sn == '?':
                sn = orientation.sinifi(y)
            if sn == '?':
                continue
            if k in gorulen:
                manifest.append(dict(kutu=k, sinif=sn, kaynak=etiket, dosya=os.path.relpath(y, a.kok),
                                     eski_yon='', cevrildi='', uzunluk='', durum='atlandi (%s kazandi)' % gorulen[k]))
                continue
            kayitlar, _ = orientation.dosya_kanonik(y)
            kayitlar = [r for r in kayitlar if len(r[1]) >= 200]
            if not kayitlar:
                manifest.append(dict(kutu=k, sinif=sn, kaynak=etiket, dosya=os.path.relpath(y, a.kok),
                                     eski_yon='', cevrildi='', uzunluk=0, durum='bos/kisa - atlandi'))
                continue
            ad, dizi, karar, cev = max(kayitlar, key=lambda r: len(r[1]))
            if karar == 'BELIRSIZ':
                belirsiz.append(dict(kutu=k, sinif=sn, kaynak=etiket,
                                     dosya=os.path.relpath(y, a.kok), uzunluk=len(dizi),
                                     N_yuzde=round(100.0 * dizi.count('N') / len(dizi), 1),
                                     not_='yon belirlenemedi - KANONIGE ALINMADI, maskeli'))
                manifest.append(dict(kutu=k, sinif=sn, kaynak=etiket, dosya=os.path.relpath(y, a.kok),
                                     eski_yon='BELIRSIZ', cevrildi='', uzunluk=len(dizi),
                                     durum='BELIRSIZ - yazilmadi'))
                continue
            gorulen[k] = etiket
            cy = os.path.join(cik, '%s.kanonik.fa' % k)
            with open(cy, 'w', encoding='utf-8') as fh:
                fh.write('>%s kanonik=%s kaynak=%s eski_yon=%s cevrildi=%s\n'
                         % (k, orientation.KANONIK_YON, etiket, karar, int(cev)))
                for i in range(0, len(dizi), 70):
                    fh.write(dizi[i:i + 70] + '\n')
            manifest.append(dict(kutu=k, sinif=sn, kaynak=etiket, dosya=os.path.relpath(y, a.kok),
                                 eski_yon=karar, cevrildi='EVET' if cev else 'hayir',
                                 uzunluk=len(dizi), durum='yazildi'))

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
    yaz('BELIRSIZ.tsv', belirsiz)
    yaz('INDEKS.tsv', [dict(kutu=m['kutu'], sinif=m['sinif'],
                            dosya='%s.kanonik.fa' % m['kutu'], kaynak=m['kaynak'],
                            eski_yon=m['eski_yon'], cevrildi=m['cevrildi'],
                            uzunluk=m['uzunluk'])
                       for m in manifest if m['durum'] == 'yazildi'])
    open(os.path.join(cik, 'OKUBENI.txt'), 'w', encoding='utf-8').write(
        u'CANONICAL CONSENSUS DIRECTORY\nValid files: *.kanonik.fa  (listed in INDEKS.tsv)\nThis directory can hold LEFTOVER files ending in *_kanonik.fasta. Those are the\nmisnamed output of the first run and could not be deleted on a mounted drive. IGNORE THEM.\nEvery consumer must read INDEKS.tsv and must NOT use glob.\nCanonical orientation: SENSE. Definition and criterion: screening/orientation.py\n')

    yazilan = [m for m in manifest if m['durum'] == 'yazildi']
    cevrilen = [m for m in yazilan if m['cevrildi'] == 'EVET']
    print(u'\ncanonical directory : %s' % cik)
    print(u'bins written        : %d' % len(yazilan))
    print('  cevrildi     : %d (ANTISENSE -> SENSE)' % len(cevrilen))
    print('  zaten sense  : %d' % (len(yazilan) - len(cevrilen)))
    print(u'BELIRSIZ       : %d (not written, BELIRSIZ.tsv)' % len(belirsiz))
    kay = {}
    for m in yazilan:
        kay[m['kaynak']] = kay.get(m['kaynak'], 0) + 1
    print(u'source distribution:', ', '.join('%s=%d' % x for x in sorted(kay.items())))

    # VERIFICATION: is every file written really SENSE
    kotu = 0
    for y in sorted(glob.glob(os.path.join(cik, '*.kanonik.fa'))):
        kayitlar, sn = orientation.dosya_kanonik(y)
        for ad, dizi, karar, cev in kayitlar:
            if cev:
                kotu += 1
    print(u'\nVERIFICATION: files in the canonical directory that still need flipping =', kotu,
          '(0 olmali)' if kotu == 0 else '*** SORUN ***')
    return 0 if kotu == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
