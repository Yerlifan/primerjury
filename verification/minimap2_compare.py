# -*- coding: utf-8 -*-
"""
minimap2_compare.py - iki hizalayiciyi AYNI girdiyle yan yana olcer.

One question alone: does minimap2 find the same identity the pure Python engine finds
mu ve ne kadar hizli. Cevap "evet" ise minimap2 varsayilan yapilabilir; "hayir"
ise mevcut motorda kalinir. Hizli olan DOGRU SAYILMAZ.
"""
# -------------------------------------------------------------------------
# minimap2_compare.py
#
# INPUT  : the real databases under REFERANS_DB/ and the real bin consensuses
#          (consensus sequences/ or konsensus_kanonik/)
# OUTPUT : MINIMAP2_KARSILASTIRMA.md (in the root directory) and the .tsv of the
#          same name
# CALLED BY: by hand. It IS NOT WIRED into the menu, because it produces a decision
#          document rather than being a stage of the chain.
#
# WHAT IT MEASURES
#   1) IS THE BEST HIT THE SAME. Do the two engines bring the same record first.
#      That is the most important criterion: the identity decision comes out of the
#      best hit.
#   2) THE IDENTITY PERCENT DEVIATION. The difference between the percentages the
#      two engines give for the same record.
#   3) IS THE RANKING KEPT. How far the sets of the first 5 hits overlap.
#   4) THE SPEED-UP. The ratio of the time taken for the same work.
#
# WHY THESE FOUR
# Looking at speed alone would be wrong: an engine that is fast but answers
# differently does not save us time, it makes us order the wrong primer. Looking at
# the best hit alone is not enough either: if the ranking breaks, the "at least two
# independent databases must agree" rule runs over different records and the
# verdict can change.
#
# IF THEY DISAGREE
# The report lists the rows that DISAGREE separately. On those rows the decision is
# not left to minimap2; which one is right is settled by hand, with an alignment.
# This script DOES NOT SAY on its own which one is right, it only shows where they
# came apart.
# -------------------------------------------------------------------------

from __future__ import print_function
import argparse
import glob
import io
import os
import sys
import time


def _kd_yukle(kok):
    import importlib.util as u
    y = os.path.join(kok, 'verification', 'identity_verification.py')
    if not os.path.exists(y):
        sys.stderr.write(u'ERROR: %s does not exist. --root must point at the project directory.\n' % y)
        return None
    sp = u.spec_from_file_location('kd', y)
    m = u.module_from_spec(sp)
    try:
        sp.loader.exec_module(m)
    except SystemExit:
        pass
    return m


def kutu_konsensuslari(kok, en_fazla):
    """Gercek kutu konsensuslerini okur. Kanonik olanlar tercih edilir."""
    adaylar = []
    kn = os.path.join(kok, 'konsensus_kanonik')
    if os.path.isdir(kn):
        adaylar = sorted(glob.glob(os.path.join(kn, '*.kanonik.fa')))
    if not adaylar:
        adaylar = sorted(glob.glob(os.path.join(kok, 'consensus sequences', '*', '*.fasta')))
    cikti = []
    for y in adaylar[:en_fazla]:
        try:
            satir = io.open(y, encoding='utf-8', errors='ignore').read().split('\n')
        except Exception:
            continue
        diz = ''.join(s.strip() for s in satir if s and not s.startswith('>'))
        if len(diz) >= 300:
            cikti.append((os.path.basename(y), diz.upper()))
    return cikti


def main():
    p = argparse.ArgumentParser(description=u'Iki hizalayiciyi yan yana olcer')
    p.add_argument('--root', dest='kok', default='.')
    p.add_argument('--bin', dest='kutu', type=int, default=4, help='how many bins to try')
    p.add_argument('--records', dest='kayit', type=int, default=120,
                   help='how many records per database to align')
    p.add_argument('--db', dest='vtb', default='', help=u'only databases whose name contains this')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)

    sys.path.insert(0, os.path.join(kok, 'verification'))
    try:
        import minimap2_aligner as mm
    except ImportError:
        sys.stderr.write(u'ERROR: verification/minimap2_aligner.py was not found.\n')
        return 1

    if not mm.var_mi():
        sys.stderr.write(
            u'\nCOMPARISON NOT POSSIBLE: mappy does not work.\n  Reason: %s\n  Install: pip install mappy\n           (or: micromamba install -n mikro -c bioconda minimap2\n                micromamba run -n mikro pip install mappy)\n\n  The current engine keeps working and the chain is NOT AFFECTED.\n  The default aligner is unchanged: pure Python.\n\n' % mm.sebep())
        return 2

    kd = _kd_yukle(kok)
    if kd is None:
        return 1

    kutular = kutu_konsensuslari(kok, a.kutu)
    if not kutular:
        sys.stderr.write(u'ERROR: no bin consensus found.\n')
        return 1

    vtb = [(e, d) for e, d, _t, kullan, _n in kd.VTB
           if kullan and os.path.exists(os.path.join(kok, 'REFERANS_DB', d))]
    if a.vtb:
        vtb = [v for v in vtb if a.vtb.lower() in v[1].lower()]
    if not vtb:
        sys.stderr.write(u'ERROR: there is no usable database under REFERANS_DB.\n')
        return 1

    satirlar = []
    t_py_top = t_mm_top = 0.0
    print(u'%d bins x %d databases, %d records per database'
          % (len(kutular), len(vtb), a.kayit))

    for kutu, q in kutular:
        for etiket, dosya in vtb:
            yol = os.path.join(kok, 'REFERANS_DB', dosya)
            kayitlar = []
            for bas, diz in kd.fasta_akisi(yol):
                kayitlar.append((bas, diz))
                if len(kayitlar) >= a.kayit:
                    break
            if not kayitlar:
                continue

            # --- motor 1: saf Python ---
            t0 = time.time()
            py = {}
            for bas, diz in kayitlar:
                py[bas] = kd.hizala(q, diz)
            t_py = time.time() - t0

            # --- engine 2: minimap2, a bulk index ---
            t0 = time.time()
            mmr = mm.toplu_hizala(q, [(bas, diz) for bas, diz in kayitlar])
            t_mm = time.time() - t0

            t_py_top += t_py
            t_mm_top += t_mm

            py_s = sorted(py.items(), key=lambda x: -x[1][0])
            mm_s = sorted(mmr.items(), key=lambda x: -x[1][0])
            py1, mm1 = py_s[0][0], mm_s[0][0]
            ust5_py = set(x[0] for x in py_s[:5])
            ust5_mm = set(x[0] for x in mm_s[:5])
            ortusme = len(ust5_py & ust5_mm)
            sapma = abs(py[py1][0] - mmr.get(py1, (0.0, 0))[0])

            satirlar.append(dict(
                kutu=kutu, vtb=etiket, kayit=len(kayitlar),
                ayni_isabet='EVET' if py1 == mm1 else 'HAYIR',
                py_en_iyi=py1[:52], py_yuzde=py[py1][0],
                mm_en_iyi=mm1[:52], mm_yuzde=mmr[mm1][0],
                ayni_kayitta_sapma=round(sapma, 2),
                ust5_ortusme='%d/5' % ortusme,
                py_sn=round(t_py, 2), mm_sn=round(t_mm, 2),
                hizlanma=round(t_py / t_mm, 1) if t_mm > 0.001 else 'olculemedi'))
            print(u'  %-28s %-22s isabet=%s  sapma=%.2f  ust5=%d/5  %.1fx'
                  % (kutu[:28], etiket[:22], satirlar[-1]['ayni_isabet'],
                     sapma, ortusme,
                     (t_py / t_mm) if t_mm > 0.001 else 0))

    if not satirlar:
        sys.stderr.write(u'ERROR: not one comparison could be made.\n')
        return 1

    ayni = sum(1 for s in satirlar if s['ayni_isabet'] == 'EVET')
    en_buyuk_sapma = max(s['ayni_kayitta_sapma'] for s in satirlar)
    tam_ortusen = sum(1 for s in satirlar if s['ust5_ortusme'] == '5/5')
    hiz = (t_py_top / t_mm_top) if t_mm_top > 0.001 else 0

    # --- THE DECISION ---
    # minimap2 IS NOT made the default unless all three conditions are met.
    sartlar = [
        (u'butun satirlarda ayni en iyi isabet', ayni == len(satirlar)),
        (u'kimlik sapmasi her satirda 0,5 puanin altinda', en_buyuk_sapma < 0.5),
        (u'ilk bes isabet her satirda birebir ortusuyor', tam_ortusen == len(satirlar)),
    ]
    gecti = all(x[1] for x in sartlar)

    tsv = os.path.join(kok, 'MINIMAP2_KARSILASTIRMA.tsv')
    bas = ['kutu', 'vtb', 'kayit', 'ayni_isabet', 'py_en_iyi', 'py_yuzde',
           'mm_en_iyi', 'mm_yuzde', 'ayni_kayitta_sapma', 'ust5_ortusme',
           'py_sn', 'mm_sn', 'hizlanma']
    with io.open(tsv, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'\t'.join(bas) + u'\n')
        for s in satirlar:
            fh.write(u'\t'.join(unicode(s[k]) if sys.version_info[0] < 3
                                else str(s[k]) for k in bas) + u'\n')

    md = os.path.join(kok, 'MINIMAP2_KARSILASTIRMA.md')
    L = []
    A = L.append
    A(u'# minimap2 against the pure Python aligner\n')
    A(u'Generated: %s · mappy %s\n' % (time.strftime('%Y-%m-%d %H:%M'), mm.surum()))
    A(u'The measurement: %d bins x %d databases = %d comparisons, %d records per database.\n' % (len(kutular), len(vtb), len(satirlar), a.kayit))
    A(u'## Karar\n')
    A(u'**%s**\n' % (u'minimap2 CAN BE MADE THE DEFAULT.' if gecti
                     else u'WE STAY ON THE CURRENT ENGINE. minimap2 was not made the default.'))
    A(u'| Condition | Result |')
    A(u'|---|---|')
    for ad, ok in sartlar:
        A(u'| %s | %s |' % (ad, u'PASSED' if ok else u'**FAILED**'))
    A(u'')
    A(u'Speed-up: **%.1f fold** (pure Python %.1f s, minimap2 %.1f s).\n'
      % (hiz, t_py_top, t_mm_top))
    if not gecti:
        A(u'> Being fast does not make it right. If even one of the conditions above is not met, minimap2 is not made the default until it is settled by hand, with an alignment, which side is right on the rows that disagree.')
    ayrilan = [s for s in satirlar if s['ayni_isabet'] == 'HAYIR'
               or s['ust5_ortusme'] != '5/5' or s['ayni_kayitta_sapma'] >= 0.5]
    if ayrilan:
        A(u'## The rows that disagree, they need confirmation by hand\n')
        A(u'| Bin | Database | Python best | % | minimap2 best | % | First 5 |')
        A(u'|---|---|---|---|---|---|---|')
        for s in ayrilan:
            A(u'| %s | %s | %s | %s | %s | %s | %s |'
              % (s['kutu'], s['vtb'], s['py_en_iyi'], s['py_yuzde'],
                 s['mm_en_iyi'], s['mm_yuzde'], s['ust5_ortusme']))
        A(u'')
    A(u'## Every measurement\n')
    A(u'The raw table: `MINIMAP2_KARSILASTIRMA.tsv`\n')
    A(u'| Bin | Database | Same hit | Deviation | First 5 | Python s | minimap2 s | Speed-up |')
    A(u'|---|---|---|---|---|---|---|---|')
    for s in satirlar:
        A(u'| %s | %s | %s | %s | %s | %s | %s | %s |'
          % (s['kutu'], s['vtb'], s['ayni_isabet'], s['ayni_kayitta_sapma'],
             s['ust5_ortusme'], s['py_sn'], s['mm_sn'], s['hizlanma']))
    A(u'\n## Where it is not used\n')
    A(u'This comparison covers **the database scan of the identity stage** only. minimap2 **is not used** in the primer binding search or in the in-silico PCR product calculation: primers are 18 to 25 bases, minimap2 was not designed for queries of that length, and when it finds no seed it misses the binding site silently. The pigeonhole engine stays there, because its losslessness is a guarantee.')
    io.open(md, 'w', encoding='utf-8').write(u'\n'.join(L) + u'\n')

    print(u'\nwritten: %s' % md)
    print(u'         %s' % tsv)
    print(u'VERDICT: %s' % (u'minimap2 as the default IS POSSIBLE'
                            if gecti else u'MEVCUT MOTORDA KALINIR'))
    print(u'HIZ    : %.1f kat' % hiz)
    return 0


if __name__ == '__main__':
    sys.exit(main() or 0)
