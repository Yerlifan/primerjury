# -*- coding: utf-8 -*-
"""ACCESS VERIFICATION - is the search REALLY using every database?

WHY IT EXISTS
-------------
Having 12 databases in the configuration does not mean they are being LOOKED AT
properly. A measured example: for the F2-1_101201 bin, the Petriella records WERE
THERE in the fungi.28S database (P. guttulata 99.65%, P. musispora 99.00%) but
they WERE NOT ENTERING the short list; the tool reported Parascedosporium at
98.16% as the best hit. The record was there and the search could not bring it
back.

THE METHOD - a recall measurement
  1) A few records SPREAD across the file are chosen from each database.
  2) That record's OWN sequence is given as the query -> does the short list bring
     it back?
  3) The same sequence is corrupted with 8% random mutation and asked again (a
     nanopore consensus carries errors, and that is real use). Does it still bring
     it back?
  If (2) fails, the search IS NOT ACTUALLY USING that database.
  If (3) fails, the search works only on flawless sequences - while in our use the
     consensuses carry errors.

The output: GECTI / DUSTU per database, plus the reason. Without this table no
identity result should be trusted.

"""

# -------------------------------------------------------------------------
# access_check.py proves whether each reference database is REALLY BEING USED BY THE
# SEARCH, by measuring recall.
#
# INPUT  : the FASTA sets under REFERENCE_DB/ (which of them to ask is read from the
#          VTB list inside identity_verification.py) and those sets' OWN records; the
#          query does not come from outside, it is chosen from inside the database.
# OUTPUT : ACCESS_RESULT/erisim_dogrulama.tsv (GECTI / KISMEN / DUSTU per database,
#          plus the reason; it is APPENDED to the file, so earlier runs are kept).
# CALLED BY: verification/full_chain.py -> key E
#          (python3 verification/access_check.py --root .)
#
# WHY A SEPARATE MEASUREMENT - THE ROOT OF THE SHORT LIST STORY
# This script produced the evidence for the bug of the period when the short list cut
# off WAS BINDING: the Petriella records WERE THERE in the fungi.28S file
# (P. guttulata 99.65%, P. musispora 99.00%) but they NEVER entered the short list,
# because the list was cut at 60 records and the ranking criterion (the seed count)
# and the deciding criterion (alignment identity) were not the same thing. So the
# tool reported the wrong genus (Parascedosporium at 98.16%) as the best hit. The fix
# was NOT to improve the ranking criterion but to grow the list to 500 and make the
# cut off non-binding. This script re-measures on every run that the fix still holds.
# -------------------------------------------------------------------------
import os, sys, csv, time, random, argparse


# -------------------------------------------------------------------------
# THE RECALL MEASUREMENT: records SPREAD across the file are chosen from each
# database, that record's OWN sequence is given as the query, and we look at whether
# the short list brings it back. If a database cannot bring back its own record, that
# database is effectively not being scanned.
#
# The second test uses a copy with 8% mutation, and THAT is the real test: our
# consensuses come from nanopore, that is, they carry errors. A search that works
# only on a flawless sequence MISSES records in real use, and the "KISMEN" verdict
# marks exactly that case.
# -------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description='Confirm that the databases can '
                                            'be read')
    p.add_argument('--root', dest='kok', default='.')
    p.add_argument('--db', dest='vtb', default=None, help='this database only (a fragment of the file name)')
    p.add_argument('--records', dest='kayit', type=int, default=3, help='how many test records per '
                                                                        'database')
    p.add_argument('--cap', dest='tavan', type=int, default=0,
                   help='at most this many records scanned (0 means all of '
                        'them, which is the default)')
    p.add_argument('--mutation', dest='mutasyon', type=float, default=0.08)
    a = p.parse_args()
    kok = os.path.abspath(a.kok)

    # ---- THE ROOT CHECK (added 2026-08-04) --------------------------------
    # This script used to print a raw Python traceback (FileNotFoundError plus the
    # traceback) when --root was given wrongly. The traceback was technically correct but
    # did not tell the user what to do. Now it looks first and loads second; what is
    # missing is named and the fix is written out.
    # NOTE: this is a change of behaviour - it stopped with a non-zero code before too,
    # and what changed is only THE MESSAGE. The measurement logic was not touched.
    _kd_yolu = os.path.join(kok, 'verification', 'identity_verification.py')
    _eksik = []
    if not os.path.isdir(os.path.join(kok, 'verification')):
        _eksik.append('verification klasoru')
    elif not os.path.exists(_kd_yolu):
        _eksik.append('verification/identity_verification.py')
    if not os.path.isdir(os.path.join(kok, 'REFERENCE_DB')):
        _eksik.append('REFERENCE_DB klasoru')
    if _eksik:
        sys.stderr.write(
            u'ERROR: %s does not contain %s.\n  This script runs from the project root. The root is the same directory as\n  verification/full_chain.py, and it holds verification/ and REFERENCE_DB/.\n  Correct use:  python3 verification/access_check.py --root <project directory>\n  If you come from the menu the root is supplied correctly on its own (key E).\n'
            % (kok, u' and '.join(_eksik)))
        return 1
    if not [f for f in os.listdir(os.path.join(kok, 'REFERENCE_DB'))
            if f.endswith(('.fasta', '.fna'))]:
        sys.stderr.write(
            u'ERROR: there is no FASTA file in the REFERENCE_DB directory (%s).\n  The access check reads the databases\' OWN rec'
            % os.path.join(kok, 'REFERENCE_DB'))
        return 1

    sys.path.insert(0, os.path.join(kok, 'verification'))
    import importlib.util as u
    sp = u.spec_from_file_location('kd', _kd_yolu)
    kd = u.module_from_spec(sp)
    try:
        sp.loader.exec_module(kd)
    except SystemExit:
        pass
    except Exception as e:
        sys.stderr.write(
            u'ERROR: identity_verification.py could not be loaded: %s: %s\n  This file may be corrupted. First check its syntax with  python3 -c "import ast; ast.parse(open(r\'%s\').read())"\n'
            % (type(e).__name__, e, _kd_yolu))
        return 1

    CIKTI = os.path.join(kok, 'ACCESS_RESULT')
    os.makedirs(CIKTI, exist_ok=True)
    vtb = [(e, d, t) for e, d, t, kullan, _n in kd.VTB
           if kullan and os.path.exists(os.path.join(kok, 'REFERENCE_DB', d))]
    if a.vtb:
        vtb = [v for v in vtb if a.vtb.lower() in v[1].lower()]

    satirlar = []
    for etiket, dosya, lokus in vtb:
        yol = os.path.join(kok, 'REFERENCE_DB', dosya)
        boyut = os.path.getsize(yol)
        print('\n=== %s (%s, %.0f MB) ===' % (etiket, dosya, boyut / 1e6), flush=True)
        # 1) yayilmis kayitlari topla
        # REZERVUAR ORNEKLEME: tek gecisde, dosyanin TAMAMINA yayilmis ornek.
        # Eski hali ilk N kaydi aliyordu; tavanin otesi HIC sinanmiyordu -
        # tam da avladigimiz kor nokta deseniydi (kullanici yakaladi).
        havuz, n, uygun = [], 0, 0
        K_ORNEK = max(a.kayit * 8, 24)
        rnd0 = random.Random(20260804)
        t0 = time.time()
        for bas, diz in kd.fasta_akisi(yol):
            n += 1
            if a.tavan and n > a.tavan:
                break
            if len(diz) < 800:
                continue
            uygun += 1
            if len(havuz) < K_ORNEK:
                havuz.append((n, bas, diz))
            else:
                j = rnd0.randrange(uygun)
                if j < K_ORNEK:
                    havuz[j] = (n, bas, diz)
        havuz.sort(key=lambda x: x[0])
        if not havuz:
            satirlar.append(dict(vtb=etiket, dosya=dosya, sonuc='OLCULEMEDI',
                                 sebep='there is no record above 800 bp (%d were scanned)' % n))
            print(u'  NOT MEASURED: no suitable record'); continue
        # bas / orta / son bolgelerden birer tane garanti
        secilen = []
        if havuz:
            for oran in [i / float(max(a.kayit - 1, 1)) for i in range(a.kayit)]:
                idx = min(len(havuz) - 1, int(oran * (len(havuz) - 1)))
                if havuz[idx] not in secilen:
                    secilen.append(havuz[idx])
        kapsam = ('TAMAMI' if not a.tavan or n <= a.tavan else
                  'THE FIRST %d RECORDS (the cap)' % a.tavan)
        print(u'  %d records scanned (%.0f s) - coverage: %s | test record: %s'
              % (n, time.time() - t0, kapsam,
                 ', '.join('#%d' % x[0] for x in secilen)), flush=True)

        tam_ok = mut_ok = 0
        ayrinti = []
        for idx, bas, diz in secilen:
            q = diz[:4000]
            kl = kd.kisa_liste(yol, q, ilerle=None)
            # THE 2026-08-09 FIX: kisa_liste() returns a list of DICTS rather than 3-tuples
            # (the keys: tohum, skor, baslik, dizi, sira, kaynak). The old line gave
            # "too many values to unpack (expected 3, got 6)" and stage E was failing with
            # exit code 1 because of it.
            bulundu = any(r['baslik'] == bas for r in kl)
            tam_ok += 1 if bulundu else 0
            rnd = random.Random(idx)
            lq = list(q)
            for _ in range(int(len(lq) * a.mutasyon)):
                i = rnd.randrange(len(lq)); lq[i] = rnd.choice('ACGT')
            kl2 = kd.kisa_liste(yol, ''.join(lq), ilerle=None)
            bulundu2 = any(r['baslik'] == bas for r in kl2)
            mut_ok += 1 if bulundu2 else 0
            ayrinti.append('#%d %s: tam=%s mutasyonlu=%s'
                           % (idx, bas[:40], 'EVET' if bulundu else 'HAYIR',
                              'EVET' if bulundu2 else 'HAYIR'))
            print('    %s' % ayrinti[-1], flush=True)
        k = len(secilen)
        if tam_ok == k and mut_ok == k:
            sonuc, sebep = 'GECTI', u'tam ve %%%d mutasyonlu sorguda %d/%d geri getirildi' % (
                int(a.mutasyon * 100), k, k)
        elif tam_ok == k:
            sonuc, sebep = 'KISMEN', ('it returned %d of %d on the exact '
                                      'sequence but only %d of %d on a query '
                                      'with %d per cent mutations, so it '
                                      'MISSES with an imperfect consensus'
                                      % (k, k, int(a.mutasyon * 100), mut_ok, k))
        else:
            sonuc, sebep = 'DUSTU', ('it could not even return its own record '
                                     '(%d of %d exact, %d of %d mutated), so '
                                     'the search IS NOT ACTUALLY USING this '
                                     'database'
                                     % (tam_ok, k, mut_ok, k))
        print('  -> %s: %s' % (sonuc, sebep), flush=True)
        satirlar.append(dict(vtb=etiket, dosya=dosya, mb=round(boyut / 1e6),
                             taranan=n, kapsam=kapsam, sinama=k, tam=tam_ok, mutasyonlu=mut_ok,
                             sonuc=sonuc, sebep=sebep, ayrinti='; '.join(ayrinti)))

    yol = os.path.join(CIKTI, 'erisim_dogrulama.tsv')
    yeni = not os.path.exists(yol)
    with open(yol, 'a', encoding='utf-8', newline='') as fh:
        if yeni:
            fh.write(u'# ACCESS VERIFICATION - does the search really use every database?\n')
            fh.write(u'# PASSED: it retrieves its own record from both the exact and the flawed sequence.\n')
            fh.write(u'# PARTIAL: it only retrieves on a flawless sequence, and our consensuses are imperfect!\n')
            fh.write(u'# FAILED: it cannot even retrieve its own record, so the database is effectively unused.\n')
        w = csv.DictWriter(fh, delimiter='\t',
                           fieldnames=['vtb', 'dosya', 'mb', 'taranan', 'kapsam', 'sinama', 'tam',
                                       'mutasyonlu', 'sonuc', 'sebep', 'ayrinti'])
        if yeni:
            w.writeheader()
        for s in satirlar:
            w.writerow({k: s.get(k, '') for k in w.fieldnames})
    print('\nwritten: %s' % yol)
    return 0


if __name__ == '__main__':
    sys.exit(main() or 0)
