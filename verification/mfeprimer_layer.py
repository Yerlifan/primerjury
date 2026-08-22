# -*- coding: utf-8 -*-
"""MFEprimer layer, the independent third tool.

Wraps the MFEprimer 4.4 binary for off-target amplicons, hairpins and dimers.
Its value is that it is NOT our code: a bug in our engine cannot confirm itself
here.

Two traps this module handles, both measured the hard way:
  * MFEprimer does not overwrite an existing output file, it returns 1 and
    leaves the old result in place, so a re-run would silently reuse stale
    numbers. Old outputs are deleted first.
  * "Off-target" in MFEprimer means BY SIZE ONLY. For group and universal
    primers, members of the target clade legitimately amplify at other lengths:
    of 1,605 amplicons flagged off-target, 95.7% were inside the target clade.
    klad_siniflandir() therefore classifies each hit taxonomically, and only
    (same domain outside clade) + (different domain) reach the verdict.

"""

# -------------------------------------------------------------------------
# mfeprimer_layer.py drives the MFEprimer 4.4.0 binary (off-target amplicons,
# hairpins, dimers); it is the INDEPENDENT third layer of evidence in the
# verification chain.
#
# INPUT  : the .fna/.fasta sets under REFERANS_DB/ and their .primerqc.bin
#          indexes (the MFE_DB list), the tools/mfeprimer binary, and the primer
#          pairs to verify (which come from the calling script).
# OUTPUT : the input and raw output files under VERIFICATION_RESULT/mfe/
#          (girdi_*.tsv, spec_*.txt, hairpin.txt, dimer.txt) plus JSON checkpoints
#          per pair and database. It writes no report of its own; the result is
#          combined by specificity_round.py.
# CALLED BY: verification/full_chain.py -> key D (indirectly: specificity_round.py
#          loads this module as "LAYER 3"; it is skipped when --no-mfe is given).
#
# WHY AN INDEPENDENT TOOL IS REQUIRED: the sample measurement and the local database
# scan use OUR same engine. If that engine has a bug, both layers go wrong in the
# same direction and "confirm" one another. Unless a tool from outside is put to the
# same questions, the contradiction stays invisible.
# -------------------------------------------------------------------------
import os, re, csv, json, subprocess, time

MFE_ADAYLARI = [os.path.join('tools', 'mfeprimer'), 'mfeprimer']
# the indexes used in the spec scan (smallest to largest, to keep the time in hand)
MFE_DB = ['archaea.16S.fna', 'bacteria.16S.fna', 'fungi.ITS.fna',
          'fungi.28SrRNA.fna', 'fungi.18SrRNA.fna',
          'SILVA_138.2_SSURef_NR99.fasta']
URUN_ALT, URUN_UST = 70, 400
HAIRPIN_TM_UST = 45.0          # the panel's own geometry rule
DIMER_TM_UST = 45.0
# D-5: MFEprimer does not report a Tm on dimer and hairpin records; the measure it
# reports is Delta G. This threshold is taken from the -9 kcal/mol limit commonly
# used for 3' dimers (IDT and Primer3 practice). The Tm<45 rule CANNOT BE APPLIED
# to this tool, and pretending it had been applied would produce a silent 'no
# violation' vote.
DG_ESIGI = -9.0


# tools/mfeprimer is looked for inside the project first and then on PATH. The copy
# inside the project is given execute permission: a file coming out of a repository
# or an archive often arrives without it, and the layer would then be skipped
# silently as "the tool is missing".
def mfe_bul(kok):
    for a in MFE_ADAYLARI:
        y = a if os.path.isabs(a) else os.path.join(kok, a)
        if os.path.exists(y):
            try:
                os.chmod(y, 0o755)
            except OSError:
                pass
            return y
        if os.path.dirname(a) == '':
            from shutil import which
            w = which(a)
            if w:
                return w
    return None


def indeks_sina(kok, mfe, dosya, yaz, sure_tavani=120):
    """Is the index REALLY being read? It tests with a single synthetic pair."""
    db = os.path.join(kok, 'REFERANS_DB', dosya)
    if not os.path.exists(db):
        return (False, u'fasta yok')
    if not os.path.exists(db + '.primerqc.bin'):
        return (False, u'.primerqc.bin indeksi yok - "mfeprimer index -i %s" ile kurulmali' % dosya)
    import tempfile
    with tempfile.TemporaryDirectory() as t:
        gi = os.path.join(t, 'sina.tsv'); co = os.path.join(t, 'sina.txt')
        open(gi, 'w').write('SINAMA\tCTGCGGTTTAATTGGATTCAACGC\tGAACTGACGACGGCCATGC\n')
        try:
            pr = subprocess.run([mfe, 'spec', '-i', gi, '-o', co, '-d', db,
                                 '-c', '2', '-S', str(URUN_UST)],
                                capture_output=True, text=True, timeout=sure_tavani)
        except subprocess.TimeoutExpired:
            return (False, u'indeks sinamasi zaman asimina ugradi (%d sn)' % sure_tavani)
        if pr.returncode != 0:
            return (False, u'mfeprimer hata verdi: %s' % ((pr.stderr or '')[:160]))
        if not os.path.exists(co):
            return (False, u'cikti uretilmedi')
        s = open(co, encoding='utf-8', errors='replace').read()
        if 'RecordCount=0' in s or 'record count' in s.lower():
            return (False, u'indeks RecordCount=0 bildiriyor - yeniden kurulmali')
        if 'Primer ID' not in s:
            return (False, u'cikti ayristirilamadi')
        # THE D-4 BUG FIX (2026-08-06): this is where the real silent failure was.
        # Although SILVA_138.2_SSURef_NR99.fasta.primerqc.bin is 515 MB, it returned
        # 'Binding Number  Plus 0  Minus 0' for EVERY primer and the amplicon count came
        # out 0. On the same FASTA our local scanner finds 419,090 hits for
        # Bakteri_universal, which means the index is not being read (it was probably
        # built with a different MFEprimer version). The old gate MISSED this: it looked
        # only for the string 'RecordCount=0' and the presence of 'Primer ID', and both
        # look sound in this situation. The result: the largest database was silently
        # voting "0 off-target amplicons". Now, if the TOTAL binding count is zero the
        # index counts as BROKEN and that database goes onto the skipped list.
        _bag = re.findall(r'^\S+\s+[ACGTUNRYKMSWBDHVacgtu]{10,}\s+\d+\s+[\d.]+\s+'
                          r'[\d.]+\s+[-\d.]+\s+(\d+)\s+(\d+)\s*$', s, re.M)
        if _bag and all(int(a) == 0 and int(b) == 0 for a, b in _bag):
            return (False, u'indeks BOZUK: butun primerler icin Binding Number 0/0 '
                           u'donuyor (dosya var ama okunmuyor). "mfeprimer index -i %s" '
                           u'ile YENIDEN kurulmali.' % dosya)
    return (True, u'okunuyor')


def _spec_ayristir(yol, beklenen_bp=None, tolerans=10):
    """Parse the amplicons.

        IMPORTANT: MFEprimer counts the target's OWN amplicon too. When the expected
        product length is given, the ones at that length are counted SEPARATELY
        ('ayni_boyda') and 'hedef_disi' covers only the ones of a DIFFERENT length. An
        amplicon coming out at the same length may also be in another organism, which
        is why it is reported as a separate column rather than ignored.

    """
    if not os.path.exists(yol):
        return {}
    s = open(yol, encoding='utf-8', errors='replace').read()
    m = re.search(r'Descriptions of \[\s*(\d+)\s*\] potential amplicons', s)
    toplam = int(m.group(1)) if m else 0
    boylar = [int(x) for x in re.findall(r'^\d+\s+\S+\s+(\d+)\s+', s, re.M)]
    ayni = hedef_disi = 0
    if beklenen_bp:
        for b in boylar:
            if abs(b - int(beklenen_bp)) <= tolerans:
                ayni += 1
            else:
                hedef_disi += 1
    else:
        hedef_disi = toplam
    return dict(_toplam=toplam, _boylar=boylar[:50],
                ayni_boyda=ayni, hedef_disi=hedef_disi)


def spec_kos(kok, mfe, ciftler, CIKTI, yaz, kontrol, sure_tavani=1800):
    """The off-target amplicon scan for each pair. A checkpoint per pair."""
    d = os.path.join(CIKTI, 'mfe')
    os.makedirs(d, exist_ok=True)
    # EVERY PAIR IS RUN SEPARATELY. Given more than one pair in a single file,
    # MFEprimer writes one combined number and that number CANNOT BE DISTRIBUTED across
    # the pairs; every pair would get the same value.
    girdiler = []
    for c in ciftler:
        ad = re.sub(r'\W+', '_', c['hedef'])
        gi = os.path.join(d, 'girdi_%s.tsv' % ad)
        with open(gi, 'w', encoding='utf-8') as fh:
            fh.write('%s\t%s\t%s\n' % (ad, c['F'], c['R']))
        girdiler.append((c, ad, gi))

    kullanilan, atlanan = [], []
    for dosya in MFE_DB:
        ok, not_ = indeks_sina(kok, mfe, dosya, yaz)
        if ok:
            kullanilan.append(dosya)
            yaz(u'    index OK : %s' % dosya)
        else:
            atlanan.append((dosya, not_))
            yaz(u'    index SKIPPED: %s - %s' % (dosya, not_))
    if not kullanilan:
        return dict(durum='ATLANDI', sebep=u'hicbir MFEprimer indeksi okunamadi',
                    atlanan=atlanan), {}

    sonuc = {}
    for c, ad, gi in girdiler:
        sonuc[c['hedef']] = {}
        for dosya in kullanilan:
            # THE 2026-08-10 SEQUENCE SEAL: the file name carried only the target name, so when
            # the sequence changed the same file was read and the old MFEprimer result came
            # back. The sequence digest is now part of the name.
            import hashlib as _hl
            _dz = _hl.md5(((c.get('F') or '') + '|' + (c.get('R') or ''))
                          .encode('utf-8')).hexdigest()[:10]
            kp = os.path.join(kontrol, 'mfe_%s_%s_%s.json'
                              % (ad, _dz, re.sub(r'\W+', '_', dosya)))
            # THE D-9 BUG FIX (2026-08-07): A POISONED CHECKPOINT.
            # The checkpoints were only being checked for "does the file exist". The SILVA
            # index had run WHILE BROKEN at 13:0x on 2026-08-06 and
            # {"_toplam":0,"hedef_disi":0,"sure":2.5} had been written for all 16 pairs. The
            # index was REBUILT at 23:30, but because those files were still in place every
            # later run would read those zeros back WITHOUT RUNNING SILVA AT ALL, that is, it
            # would produce a SILENT and WRONG assurance that "there is nothing off-target in
            # the largest database". Exactly the failure pattern we have to avoid.
            # The rule: if the index is NEWER than the checkpoint, the checkpoint is INVALID.
            if os.path.exists(kp):
                _ix = os.path.join(kok, 'REFERANS_DB', dosya + '.primerqc.bin')
                _bayat = False
                try:
                    if os.path.exists(_ix) and os.path.getmtime(_ix) > os.path.getmtime(kp):
                        _bayat = True
                except OSError:
                    pass
                if _bayat:
                    yaz(u'    checkpoint STALE (the index is newer), re-running: %s / %s' % (ad[:28], dosya))
                    try:
                        os.remove(kp)
                    except OSError:
                        pass
                else:
                    try:
                        sonuc[c['hedef']][dosya] = json.load(open(kp, encoding='utf-8'))
                        continue
                    except Exception:
                        pass
            co = os.path.join(d, 'spec_%s_%s.txt' % (ad, re.sub(r'\W+', '_', dosya)))
            # THE D-10 BUG FIX (2026-08-07): MFEprimer 4.4.0 REFUSES TO RUN when the output
            # file EXISTS:
            #   "... is already exists, please remove it first"  (returncode 1)
            # This was a trap that never saw daylight: the spec_*.txt files of earlier runs sat
            # under VERIFICATION_RESULT/mfe/ (16 SILVA outputs, all of them written while the index
            # was broken, so all of them "0 potential amplicons"). While the checkpoint
            # matched, mfeprimer was never called and the error was invisible; and when the
            # stale checkpoint was cleared and a re-run attempted, mfeprimer returned 1, we
            # wrote dict(hata=...) and moved on, and SILVA was UNMEASURED again. So even after
            # we said "we fixed the index", the result came out silently incomplete once more.
            # The fix: DELETE the old output (and the .spec.tsv/.mfe.log beside it) before
            # running. If it cannot be deleted, DO NOT PASS SILENTLY; write a clear error.
            _silinemedi = None
            for _es in (co, co + '.spec.tsv', co + '.mfe.log'):
                if os.path.exists(_es):
                    try:
                        os.remove(_es)
                    except OSError as _e:
                        if _es == co:
                            _silinemedi = str(_e)
            if _silinemedi:
                sonuc[c['hedef']][dosya] = dict(
                    hata=u'eski cikti dosyasi silinemedi, MFEprimer uzerine YAZMAZ: %s'
                         % _silinemedi)
                yaz(u'    ERROR: the old output could not be deleted, so this database was not measured: %s'
                    % co)
                continue
            t0 = time.time()
            # THE D-11 SPEED FIX (2026-08-07, MEASURED): MFEprimer produces its amplicon report
            # with A LARGE NUMBER OF SMALL writes. Writing to a mounted Windows directory
            # (/mnt/c/...) accumulates the latency of every write and the job slows to the
            # point of being unrecognisable. THE MEASUREMENT:
            #   Proteolitik_Synergistaceae, SILVA, the same command:
            #     to a local disk (/tmp)      : 5.2 s   (11.5 MB of output)
            #     straight to the mounted dir : DID NOT FINISH in 3 minutes (2.35 MB written)
            #   Bakteri_universal, SILVA:
            #     to a local disk (/tmp)      : 49.8 s  (282 MB of output)
            #     straight to the mounted dir : ~38 minutes in the 2026-08-06 run
            # Against that, a BULK copy is fast (233 MB / 26 s measured), because it is one
            # large sequential write. So: mfeprimer runs on the LOCAL disk and the evidence
            # files are copied over IN BULK when the job finishes. No evidence file is lost;
            # only the write pattern changes.
            import tempfile, shutil
            _yerel = tempfile.mkdtemp(prefix='mfe_spec_')
            _co = os.path.join(_yerel, os.path.basename(co))
            try:
                try:
                    pr = subprocess.run([mfe, 'spec', '-i', gi, '-o', _co, '-d',
                                         os.path.join(kok, 'REFERANS_DB', dosya),
                                         '-c', '4', '-s', str(URUN_ALT), '-S', str(URUN_UST)],
                                        capture_output=True, text=True, timeout=sure_tavani)
                    if pr.returncode != 0:
                        sonuc[c['hedef']][dosya] = dict(hata=(pr.stderr or '')[:200]); continue
                except subprocess.TimeoutExpired:
                    sonuc[c['hedef']][dosya] = dict(hata='zaman asimi'); continue
                v = _spec_ayristir(_co, c.get('urun'))
                v['sure_hesap'] = round(time.time() - t0, 1)
                # kanit dosyalarini bagli klasore TOPLU yaz
                for _ek in ('', '.spec.tsv', '.mfe.log'):
                    if os.path.exists(_co + _ek):
                        try:
                            shutil.copyfile(_co + _ek, co + _ek)
                        except OSError as _e:
                            yaz(u'    WARNING: the evidence file could not be copied (%s): %s'
                                % (os.path.basename(co + _ek), _e))
            finally:
                shutil.rmtree(_yerel, ignore_errors=True)
            v['sure'] = round(time.time() - t0, 1)
            sonuc[c['hedef']][dosya] = v
            json.dump(v, open(kp, 'w', encoding='utf-8'), ensure_ascii=False, default=str)
        t = sum(x.get('hedef_disi', 0) for x in sonuc[c['hedef']].values()
                if isinstance(x, dict) and not x.get('hata'))
        a = sum(x.get('ayni_boyda', 0) for x in sonuc[c['hedef']].values()
                if isinstance(x, dict) and not x.get('hata'))
        yaz(u'    %-34s off-target %d, same length as the target %d' % (c['hedef'][:34], t, a))
    return dict(durum='TAMAM', kullanilan=kullanilan, atlanan=atlanan), sonuc


def yapi_kos(kok, mfe, ciftler, CIKTI, yaz):
    """hairpin plus dimer. Evaluated by the panel's own rules (Tm < 45)."""
    d = os.path.join(CIKTI, 'mfe')
    os.makedirs(d, exist_ok=True)
    fa = os.path.join(d, 'oligolar.fa')
    with open(fa, 'w', encoding='utf-8') as fh:
        for c in ciftler:
            ad = re.sub(r'\s+', '_', c['hedef'])
            fh.write('>%s_fp\n%s\n>%s_rp\n%s\n' % (ad, c['F'], ad, c['R']))
    out = {}
    for komut, ust in (('hairpin', HAIRPIN_TM_UST), ('dimer', DIMER_TM_UST)):
        co = os.path.join(d, '%s.txt' % komut)
        # D-10 (2026-08-07): the same trap as spec_kos - MFEprimer DOES NOT OVERWRITE an
        # existing output file, it returns 1 ("is already exists"). hairpin.txt and
        # dimer.txt were left over from earlier runs; without the fix these two measurements
        # would silently become dict(hata=...) on every re-run.
        for _es in (co, co + '.mfe.log'):
            if os.path.exists(_es):
                try:
                    os.remove(_es)
                except OSError:
                    pass
        if os.path.exists(co):
            out[komut] = dict(hata=u'eski cikti silinemedi, MFEprimer uzerine yazmaz: %s' % co)
            yaz(u'    ERROR: the old %s output could not be deleted, so it was not measured' % komut); continue
        try:
            pr = subprocess.run([mfe, komut, '-i', fa, '-o', co],
                                capture_output=True, text=True, timeout=600)
            if pr.returncode != 0:
                out[komut] = dict(hata=(pr.stderr or '')[:160]); continue
        except subprocess.TimeoutExpired:
            out[komut] = dict(hata='zaman asimi'); continue
        s = open(co, encoding='utf-8', errors='replace').read()
        m = re.search(r'%s List \((\d+)\)' % komut.capitalize(), s)
        n = int(m.group(1)) if m else 0
        # THE D-5 BUG FIX (2026-08-06): the old code looked for 'Tm:\s*([0-9.]+)'.
        # MFEprimer 3.0/4.4 report NO Tm ON dimer and hairpin RECORDS; every record says
        # 'Score: N, Delta G = -X.XX kcal/mol' ('Tm' appears only in the HEADING of the
        # oligo summary table above, without a colon). The result: tms stayed empty,
        # en_yuksek_tm was None, kural_ihlali was always False, and so this layer was
        # SILENTLY voting "no violation" despite 66 dimer records.
        # Now the measure that is actually reported (Delta G) is parsed, and the fact that
        # the Tm IS NOT COMPUTED is marked openly.
        dgs = [float(x) for x in re.findall(r'Delta G\s*=\s*(-?[0-9.]+)', s)]
        skor = [int(x) for x in re.findall(r'Score:\s*(\d+)', s)]
        tms = [float(x) for x in re.findall(r'Tm:\s*([0-9.]+)', s)]
        out[komut] = dict(sayi=n, en_yuksek_tm=(max(tms) if tms else None),
                          tm_hesaplanmadi=(not tms),
                          en_dusuk_dg=(min(dgs) if dgs else None),
                          en_yuksek_skor=(max(skor) if skor else None),
                          dg_esigi=DG_ESIGI,
                          kural_ihlali=bool(dgs and min(dgs) <= DG_ESIGI),
                          hukum_verilebilir=bool(dgs),
                          tm_kurali_uygulanamadi_ust=ust, dosya=co)
        yaz(u'    %s: %d records, lowest Delta G %s kcal/mol (a Tm per record is NOT '
            u'COMPUTED by this tool, so the Tm<%.0f rule could not be applied)'
            % (komut, n, ('%.2f' % min(dgs)) if dgs else '-', ust))
    return out


def hedef_disi_kimlikleri(CIKTI, ciftler, yaz, tolerans=10):
    """D-7 (2026-08-06): THE COUNT of 'N off-target amplicons' does not decide on its own.

        MFEprimer's off-target measure rests ON LENGTH ALONE: every amplicon outside
        the expected product length plus or minus a tolerance counts as 'off-target'.
        But for group and UNIVERSAL primers, an amplicon from a member inside the
        target clade can also come out at a different length (natural indel length
        polymorphism in 16S/18S). An example: ALL 31 of Bakteri_universal's
        'off-target' amplicons are bacterial 16S (Thermodesulfovibrio, Desulfurella,
        Gemmata, Planctomycetes and so on), and for a universal bacterial primer those
        are INSIDE the target BY DESIGN.

        So the IDENTITY is written beside the count: which accession, which organism,
        how many bp. The decision is made by looking at that file.

    """
    d = os.path.join(CIKTI, 'mfe')
    yol = os.path.join(CIKTI, 'mfe_hedef_disi_kimlikleri.tsv')
    # THE A4 FIX (2026-08-21): this used to be a silent 'pass'.
    # If 'urun' could not be converted to a number, that target NEVER entered the 'bek'
    # dictionary; and because beklenen_bp was then unknown, ALL of that target's
    # amplicons counted as off-target. So ONE malformed cell SILENTLY hardened that
    # row's verdict, which is this project's chief kind of bug: a wrong answer with no
    # error. Every target that drops is now written to the screen and to the file.
    bek = {}
    bek_dusen = []
    for c in ciftler:
        try:
            bek[re.sub(r'\W+', '_', c['hedef'])] = (c['hedef'], int(c['urun']))
        except (TypeError, ValueError, KeyError) as e:
            bek_dusen.append((c.get('hedef', '?'), repr(c.get('urun')), type(e).__name__))
    if bek_dusen:
        yaz(u'  WARNING: the expected product length could not be read for %d targets. For those, the "same length" separation CANNOT be made and every amplicon'
            % len(bek_dusen))
        for h, v, e in bek_dusen:
            yaz(u'    %-40s product=%s (%s)' % (h[:40], v, e))
    satir = 0
    with open(yol, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# The IDENTITY of MFEprimer "off-target" amplicons (D-7).\n# The off-target measure rests on LENGTH ALONE; below')
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'veritabani', 'amplikon_bp', 'beklenen_bp', 'erisim',
                    'organizma', 'FpTm', 'RpTm', 'Ta', 'termodinamik_notu'])
        for dosya in sorted(os.listdir(d)) if os.path.isdir(d) else []:
            if not dosya.startswith('spec_') or not dosya.endswith('.txt'):
                continue
            ad = None
            for k in bek:
                if dosya.startswith('spec_%s_' % k) and (ad is None or len(k) > len(ad)):
                    ad = k
            if ad is None:
                continue
            hedef, b = bek[ad]
            db = dosya[len('spec_%s_' % ad):-4]
            s2 = open(os.path.join(d, dosya), encoding='utf-8', errors='replace').read()
            tab = {}
            for m in re.finditer(r'^\s*(\d+)\s+(\S+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+'
                                 r'[-\d.]+\s+[-\d.]+\s+[\d.]+\s+([\d.]+)\s*$', s2, re.M):
                if abs(int(m.group(3)) - b) > tolerans:
                    tab[int(m.group(1))] = (int(m.group(3)), m.group(2), m.group(4),
                                            m.group(5), m.group(6))
            if not tab:
                continue
            org = {}
            for m in re.finditer(r'^Amp (\d+): .*?==> (\S+?):(\d+)-(\d+) (.*)$', s2, re.M):
                org[int(m.group(1))] = m.group(5).strip()
            for i in sorted(tab):
                bp, er, ftm, rtm, ta = tab[i]
                try:
                    _not = (u'primer baglanma Tm (%s/%s C) tavlama sicakligindan (%s C) '
                            u'DUSUK - bu urun standart kosulda OLUSMAZ'
                            % (ftm, rtm, ta)) if min(float(ftm), float(rtm)) < float(ta) - 5 \
                        else u'primer Tm tavlama sicakligina yakin - GERCEK risk'
                except ValueError:
                    _not = u''
                w.writerow([hedef, db, bp, b, er, org.get(i, ''), ftm, rtm, ta, _not])
                satir += 1
    yaz(u'  written: %s (the identity of %d off-target amplicons)' % (yol, satir))
    return yol, satir


def hedef_disi_say(spec_sonuc, hedef_ad, ciftler=None):
    """The TOTAL amplicon count MFEprimer reports for this pair.

        NOTE: MFEprimer counts the target's OWN amplicon too. For the off-target
        measure the TOTAL is used rather than 'the ones not at the expected product
        length', and that is reported openly, so that the contradiction flag catches it
        when two tools are not counting the same thing. The raw files sit under mfe/.

    """
    v = (spec_sonuc or {}).get(hedef_ad) or {}
    return sum(x.get('hedef_disi', 0) for x in v.values()
               if isinstance(x, dict) and not x.get('hata'))


# -------------------------------------------------------------------------
# D-12 (2026-08-07): THE COUNT of "N off-target amplicons" DOES NOT DECIDE.
#
# THE MEASURED REASONING: in the 2026-08-07 12:11 run, MFEprimer found 1605
# amplicons of a different length. Classifying the taxonomy strings in
# mfe_hedef_disi_kimlikleri.tsv gives this distribution:
#     (a) INSIDE the target clade, different length   1536   95.7%
#    (ao) inside the target domain but an ORGANELLE     31   1.9%   (plant mitochondrion/chloroplast)
#     (b) the same domain, outside the clade            24   1.5%
#     (c) a different domain (a real cross reaction)    14   0.9%
# So 95.7% of the raw count is a harmless length variant. Of Bakteri_universal's
# 1550 "off-target" hits, 1519 begin with "Bacteria;..." in SILVA, which for a
# universal bacterial primer is INSIDE the target BY DESIGN.
#
# This function produces the measure that enters the verdict: klad_disi = (b) + (c).
# The primer Tm of every (b)/(c) record is also compared against the panel's
# annealing temperature; if the Tm is clearly below it (by >=TM_MARJI C) that
# product does not form under standard conditions and does not enter the "can form"
# count.
# -------------------------------------------------------------------------
TM_MARJI = 5.0          # min(FpTm,RpTm) < Ta - TM_MARJI ise urun olusmaz sayilir
ORGANEL_JETONLARI = ('Chloroplast', 'Mitochondria')


def klad_tablosu(kok):
    """screening/hedef_klad.tsv -> {target: (domain, [clade tokens], source)}"""
    y = os.path.join(kok, 'screening', 'hedef_klad.tsv')
    if not os.path.exists(y):
        return {}
    out = {}
    for r in csv.DictReader((l for l in open(y, encoding='utf-8')
                             if not l.startswith('#')), delimiter='\t'):
        h = (r.get('hedef') or '').strip()
        if not h:
            continue
        out[h] = ((r.get('alan') or '').strip(),
                  [x.strip() for x in (r.get('klad') or '').split(',') if x.strip()],
                  (r.get('kaynak') or '').strip())
    return out


# The database in the spec file name -> its domain. RefSeq headers carry NO taxonomy
# string; the domain of those records comes from the DEFINITION of the database
# (every record inside bacteria.16S.fna is bacterial 16S by definition). Not a guess.
VTB_ALAN = {'archaea_16S_fna': 'Archaea', 'bacteria_16S_fna': 'Bacteria',
            'fungi_ITS_fna': 'Eukaryota', 'fungi_28SrRNA_fna': 'Eukaryota',
            'fungi_18SrRNA_fna': 'Eukaryota'}
VTB_KLAD = {'fungi_ITS_fna': 'Fungi', 'fungi_28SrRNA_fna': 'Fungi',
            'fungi_18SrRNA_fna': 'Fungi'}


def _kayit_coz(org, vtb):
    """(alan, jetonlar, organel_mi)"""
    if ';' in org and org.split(';')[0].strip() in ('Bacteria', 'Archaea', 'Eukaryota'):
        t = [x.strip() for x in org.split(';')]
        return (t[0], t, any(o in t for o in ORGANEL_JETONLARI))
    a = VTB_ALAN.get(vtb, '?')
    t = [org, a] + ([VTB_KLAD[vtb]] if vtb in VTB_KLAD else [])
    return (a, t, False)


def klad_siniflandir(kok, CIKTI, ciftler, ta_panel, yaz=None):
    """Reads mfe_hedef_disi_kimlikleri.tsv and counts a/ao/b/c per pair.

        Returns: {target: dict(a, ao, b, c, klad_disi, olusabilir, olusmaz,
                            klad, alan, kaynak, kayitlar=[...])}
        'klad_disi' is the measure that enters the verdict. 'olusabilir' is the number
        of those records whose primer Tm is close to the annealing temperature.
        If the table or the file is missing it returns EMPTY, and the caller then
        CARRIES ON using the raw count but MUST mark that openly as "the clade
        separation COULD NOT BE MADE" (counting it clean silently is forbidden).

    """
    tab = klad_tablosu(kok)
    yol = os.path.join(CIKTI, 'mfe_hedef_disi_kimlikleri.tsv')
    if not tab or not os.path.exists(yol):
        if yaz:
            yaz(u'    clade separation NOT POSSIBLE (%s missing); the raw count will be used'
                % ('hedef_klad.tsv' if not tab else os.path.basename(yol)))
        return {}
    out = {}
    for c in ciftler:
        h = c['hedef']
        if h in tab:
            out[h] = dict(a=0, ao=0, b=0, c=0, klad_disi=0, olusabilir=0,
                          olusmaz=0, alan=tab[h][0], klad=tab[h][1],
                          kaynak=tab[h][2], kayitlar=[])
    eksik = set()
    for r in csv.DictReader((l for l in open(yol, encoding='utf-8')
                             if not l.startswith('#')), delimiter='\t'):
        h = r.get('hedef')
        if h not in out:
            if h:
                eksik.add(h)
            continue
        d = out[h]
        k_alan, jet, organel = _kayit_coz(r.get('organizma') or '', r.get('veritabani') or '')
        ic = any(j in jet for j in d['klad'])
        if not ic and ';' not in (r.get('organizma') or ''):
            ic = any(re.search(r'\b%s' % re.escape(j), r.get('organizma') or '', re.I)
                     for j in d['klad'])
        if ic and organel:
            s = 'ao'
        elif ic:
            s = 'a'
        elif k_alan == d['alan']:
            s = 'b'
        else:
            s = 'c'
        d[s] += 1
        if s in ('b', 'c'):
            d['klad_disi'] += 1
        if s in ('b', 'c', 'ao'):
            try:
                dtm = min(float(r['FpTm']), float(r['RpTm'])) - float(ta_panel)
            except (ValueError, KeyError, TypeError):
                dtm = None
            if dtm is None:
                pass
            elif dtm < -TM_MARJI:
                d['olusmaz'] += 1
            else:
                d['olusabilir'] += 1
            d['kayitlar'].append(dict(sinif=s, bp=r.get('amplikon_bp'),
                                      erisim=r.get('erisim'),
                                      organizma=r.get('organizma'),
                                      FpTm=r.get('FpTm'), RpTm=r.get('RpTm'),
                                      dTm=dtm))
    if eksik and yaz:
        yaz(u'    WARNING: target(s) with no entry in hedef_klad.tsv: %s'
            % ', '.join(sorted(eksik))[:200])
    return out
