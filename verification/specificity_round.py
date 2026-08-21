# -*- coding: utf-8 -*-
"""Specificity verification, four independent evidence layers before ordering.

WHY A SEPARATE STAGE
    A pair that looks clean in your sample can still fail in the lab: a sample
    is 99 bins, not the world. A pair may face no competitor in the sample and
    still amplify something abundant that the sample never captured.

THE FOUR LAYERS, AND WHY THEY ARE BUILT DIFFERENTLY
    1  in-sample    in-silico PCR against the raw reads.
                    Independent of every reference database.
    2  local DB     scan of the local reference sets. Uses OUR code and OUR
                    engine, so it shares our bugs, and cannot corroborate
                    itself. That is exactly why layers 3 and 4 exist.
    3  MFEprimer    an external binary. Different implementation, different
                    thermodynamic model, written by other people.
    4  NCBI         Primer-BLAST against nt. Independent of our choice of
                    databases, which is layer 2's main blind spot.

    If the layers disagree the row is marked CELISKILI and is NOT orderable.
    Disagreement is treated as information, not noise: a contradiction means at
    least one measurement is wrong, and shipping either would be a gamble.

THREE STATES, NOT TWO
    BILINMIYOR (unknown) is distinct from TEMIZ (clean). A layer that did not
    run never votes in favour. This is the single most common way pipelines of
    this kind mislead: an unmeasured check silently reads as a passed check.

    Related: the in-sample column is DISPLAYED but does not VOTE. It is the
    admission criterion that put the pair on this list, so by construction it
    can never say RISKLI. Letting a constant vote made 16 of 16 rows come out
    CELISKILI mechanically (fixed 2026-08-06).

"""

# -------------------------------------------------------------------------
# specificity_round.py takes the NEW or CHANGED pairs that passed the threshold in
# the recovery round and tests them side by side against four independent layers of
# evidence, before any of them is ordered.
#
# INPUT  : KURTARMA_SONUC/kurtarma_satirlari.tsv (which pairs to verify),
#          TEK_PROTOKOL_SONUC/panel_tek_protokol.tsv (the primer sequences),
#          the KUMELER list under REFERANS_DB/, the tools/mfeprimer indexes,
#          NCBI Primer-BLAST (over the network) or a hand filled
#          NCBI_SONUC_SABLONU.tsv.
# OUTPUT : DOGRULAMA_SONUC/dogrulama_uc_sutun.tsv (the main table),
#          DOGRULAMA_SONUC/CELISKILER.md (read this one first),
#          DOGRULAMA_SONUC/yerel_vuruslar.tsv, DOGRULAMA_RAPORU.md,
#          NCBI_PRIMER_BLAST_GIRDI.tsv + NCBI_SONUC_SABLONU.tsv (the manual route),
#          DOGRULAMA_SONUC/kontrol/ (a checkpoint per set).
# CALLED BY: verification/full_chain.py -> key D
#          (python3 verification/specificity_round.py --root . ...)
#
# WHY FOUR LAYERS: a pair that looks good in the sample may not hold in the lab,
# because the sample is 99 bins, not the world. Layers 1 and 2 are OUR code and
# use the same engine, so if that engine has a bug both go wrong in the same
# direction. Layer 3 (MFEprimer) is an independent tool from outside, and layer 4
# (NCBI) is a source independent of our own database choices. When the layers
# disagree the row is marked CONTRADICTORY and CANNOT BE ORDERED.
# -------------------------------------------------------------------------
import os, sys, csv, json, time, argparse, re

VERSIYON = '1.0 (2026-08-03)'

URUN_ALT, URUN_UST = 70, 400          # yalanci urun aranan mesafe araligi

# VURUS_ESIGI = 0  -> even ONE off-target product makes the layer RISKLI.
#
# A3 (2026-08-21): this value stood without a reason, which made it the exception
# in a project where every constant carries one. The reason is this, and the
# strictness is DELIBERATE:
#   * The panel's acceptance rule is "not one wrong primer goes to order". A wrong
#     order is both expensive and slow; the cost of missing a false product is
#     larger than the cost of unfairly rejecting a candidate.
#   * The threshold DOES NOT DECIDE ON ITS OWN. A RISKLI layer is only a VOTE; the
#     decision is tied to the AGREEMENT of the four layers (see birlestir()). So a
#     zero threshold does not block an order by itself, it requires a second source
#     to confirm.
#   * What is being counted has already passed A FILTER: products inside the
#     URUN_ALT..URUN_UST range that satisfy the binding rule. Random noise does not
#     get this far.
#
# WHEN TO CHANGE IT: in exploration mode (the question "which candidate has any
# chance at all") it can be raised to 1 or 2. It MUST NOT BE RAISED in a pre-order
# verdict. If you change it, the report says which threshold the run used, so it
# never changes silently.
VURUS_ESIGI = int(os.environ.get('PT_VURUS_ESIGI', '0'))
NCBI_SONUC_TAVANI = 1000               # the product cap of the Primer-BLAST result page.
                                       # If n == the cap, the value IS NOT A COUNT. D-3.
# D-12 (2026-08-07): the panel's real annealing temperature. Whether an off-target
# amplicon counts as "able to form" is judged against this.
# The source: the panel's own Ta (a panel decision, 2026-08-07). MFEprimer writes
# ITS OWN Ta for every amplicon; that value is derived from the amplicon's GC and
# IS NOT the temperature we will use in our thermocyclers.
TA_PANEL = 57.9
BOY_TOL = 10                           # beklenen urun boyuna bu kadar yakin vurus
                                       # It is THE TARGET'S OWN product (the same tolerance as the
                                       # MFEprimer layer). The D-1 correction.

# The local sets to scan: (label, file name, description)
KUMELER = [
    ('SILVA SSU NR99', 'SILVA_138.2_SSURef_NR99.fasta', u'510 495 kayit; SSU (16S/18S)'),
    ('SILVA LSU NR99', 'SILVA_138.2_LSURef_NR99.fasta', u'95 279 kayit; LSU (23S/28S)'),
    ('UNITE ITS', 'UNITE_ITS.fasta', u'2 069 189 kayit; mantar ITS'),
    ('PR2 SSU', 'PR2_SSU_taxo_long.fasta', u'240 201 kayit; okaryot 18S'),
    ('ROD operon', 'ROD_v1.2_operon_variants.fasta', u'60 320 kayit; rRNA operon varyantlari'),
    ('RefSeq bakteri 16S', 'bacteria.16S.fna', u'26 877 kayit'),
    ('RefSeq arke 16S', 'archaea.16S.fna', u'1 160 kayit'),
    ('RefSeq mantar ITS', 'fungi.ITS.fna', u'20 394 kayit'),
    ('RefSeq mantar 28S', 'fungi.28SrRNA.fna', u'12 890 kayit'),
    ('RefSeq mantar 18S', 'fungi.18SrRNA.fna', u'4 037 kayit'),
    ('RefSeq ref_all2', 'ref_all2.fna', u'65 358 kayit; RefSeq birlesik'),
]
# SILVA Parc is DELIBERATELY not used in the SPECIFICITY scan: a non-dereplicated
# set of 1.3 million records adds no new information to the false product question
# over NR99, but it lengthens the run by hours. THE IDENTITY question is different;
# there Parc is REQUIRED (see identity_verification.py), because NR99 deletes rare
# genera. It can be turned on if wanted:
PARC_ISTEGE_BAGLI = ('SILVA LSU Parc', 'SILVA_138.2_LSUParc.fasta',
                     u'1 312 521 kayit; tekrarsizlastirilmamis')

OLCUT_NOTU = u"""
YEREL TARAMANIN OLCUTU - ACIKCA YAZILIYOR
=========================================
Kullanilan kod: screening/global_scan.py  (AYNEN, degistirilmedi)

  * iki primer KARSILIKLI YONELIMDE baglanacak
  * aralarindaki mesafe %d-%d bp
  * F ve R uyumsuzluklari TOPLAM en cok %s

3' SON IKI BAZ SARTI BU KATMANDA UYGULANMADI. Sebep: mevcut kuresel tarama
kodu bu sarti tasimiyor ve o kodu yeniden yazmamak icin oldugu gibi kullanildi.
Bu olcut, 3' son iki baz sarti olan olcutten DAHA GEVSEKTIR - yani bulunan
vuruslarin bir kismi gercekte urun VERMEYEBILIR. Risk taramasinda bu GUVENLI
taraftir: gercek riski gozden kacirmaktansa fazladan uyari uretir.
Bir vurus ciddiye alinacaksa 3' ucunun tuttugu ayrica bakilmalidir.
"""


def sure_metni(sn):
    sn = int(sn)
    if sn < 90:
        return '%d saniye' % sn
    if sn < 5400:
        return '%d dakika' % round(sn / 60.0)
    return '%.1f saat' % (sn / 3600.0)


def vir(x, b=2):
    if x is None or x == '':
        return '-'
    try:
        return ('%.*f' % (b, float(str(x).replace(',', '.')))).replace('.', ',')
    except (TypeError, ValueError):
        return str(x)


# ---------------------------------------------------------------- input
_ATLANAN = []


# -------------------------------------------------------------------------
# Builds the set of pairs to verify. There are two kinds of row:
#   A NEW PAIR  - the recovery round found a new F/R; the sequences are taken
#                 straight out of the row's text.
#   CHANGED     - the primers are the same and what changed is the MEASUREMENT or
#                 the MEMBERSHIP; the sequences come from stage P's panel table.
# A row whose primer sequence cannot be found in any source is NOT skipped
# silently; it goes into the _ATLANAN list and is shown in the run header.
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# THE ORDER LIST MODE (--siparis) - 2026-08-06
#
# WHY IT WAS ADDED: by default this script only verified the RECOVERED pairs
# (kurtarma_satirlari.tsv). In the last run that list was a single row, which means
# 14 of the 15 pairs going to order NEVER SAW the Primer-BLAST layer. Since the
# question asked before an order is "does every pair on the list carry a false
# product risk", the input set has to be the order list.
#
# The same four layers run (sample / local DB / MFEprimer / NCBI Primer-BLAST); the
# only thing that changes is WHICH pairs are tested.
# -------------------------------------------------------------------------
def siparistekiler(kok, hepsi=False):
    """SIPARIS_LISTESI.tsv -> siparise giden ciftler (KESIN + EVRENSEL).

    hepsi=True ise KOSULLU ve ONERILMEZ satirlar da alinir - kullanici onlari
    da sinatmak isterse. Varsayilan: yalniz siparise gidecekler.
    """
    yol = os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'SIPARIS_LISTESI.tsv')
    if not os.path.exists(yol):
        yol = os.path.join(kok, 'SIPARIS_LISTESI.tsv')
    if not os.path.exists(yol):
        sys.exit('HATA: SIPARIS_LISTESI.tsv yok.\n'
                 '      Once verification/full_chain.py -> secenek (T) kosulmalidir.')
    with open(yol, encoding='utf-8') as fh:
        satirlar = list(csv.DictReader(
            (s for s in fh if not s.startswith('#')), delimiter='\t'))
    kabul = ('KESIN', 'EVRENSEL') if not hepsi else ('KESIN', 'EVRENSEL',
                                                     'KOSULLU', 'ONERILMEZ')
    out = []
    for s in satirlar:
        sn = (s.get('SINIF') or '').strip().upper()
        if sn not in kabul:
            continue
        F, R = (s.get('F') or '').strip(), (s.get('R') or '').strip()
        if not F or not R:
            _ATLANAN.append(s.get('hedef', '?'))
            continue
        out.append(dict(hedef=s['hedef'], F=F, R=R,
                        urun=s.get('urun_bp', ''),
                        sinif=sn,
                        # SIPARISTE_MI: KESIN/EVRENSEL satirlar siparise gider.
                        # --ncbi-yalniz-siparis bu bayragi kullanir (2026-08-07).
                        sipariste=(sn in ('KESIN', 'EVRENSEL')),
                        tur=u'SIPARIS LISTESI (%s)' % sn,
                        numune_deger=s.get('ayrim_mm1', ''),
                        numune_olcu=u'dCq %s' % (s.get('dCq_karsiligi') or '-'),
                        yol=s.get('siparis_sarti', '')))
    return out, yol


def kurtarilanlar(kok):
    """KURTARMA_SONUC/kurtarma_satirlari.tsv -> esigi gecen YENI/DEGISMIS ciftler."""
    yol = os.path.join(kok, 'KURTARMA_SONUC', 'kurtarma_satirlari.tsv')
    if not os.path.exists(yol):
        sys.exit('HATA: %s yok.\n      Once verification/full_chain.py -> secenek (K) kosulmalidir.' % yol)
    with open(yol, encoding='utf-8') as fh:
        satirlar = list(csv.DictReader((s for s in fh if not s.startswith('#')), delimiter='\t'))

    tp = os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'panel_tek_protokol.tsv')
    ciftler = {}
    if os.path.exists(tp):
        with open(tp, encoding='utf-8') as fh:
            for r in csv.DictReader((s for s in fh if not s.startswith('#')), delimiter='\t'):
                ciftler[r['hedef']] = (r.get('F', ''), r.get('R', ''), r.get('urun_bp', ''))

    out = []
    for s in satirlar:
        if not (s.get('esigi_gecti_mi') or '').startswith('EVET'):
            continue
        yeni = s.get('yeni_deger') or ''
        m = re.search(r'YENI CIFT\s+([ACGT]+)\s*/\s*([ACGT]+)\s*\((\d+)\s*bp\)', yeni)
        if m:
            F, R, bp = m.group(1), m.group(2), m.group(3)
            tur = 'YENI CIFT'
        else:
            F, R, bp = ciftler.get(s['hedef'], ('', '', ''))
            tur = 'DEGISMIS (ayni primerler, olcu/uyelik degisti)'
        if not F or not R:
            _ATLANAN.append(s['hedef'])
            continue
        out.append(dict(hedef=s['hedef'], F=F, R=R, urun=bp, tur=tur,
                        numune_deger=s.get('yeni_deger') or s.get('eski_deger'),
                        numune_olcu=s.get('olcu', ''), yol=s.get('denenen_yol', '')))
    return out, yol


# ---------------------------------------------------------------- layer 1
# -------------------------------------------------------------------------
# LAYER 2 - the local database scan. What it looks for: is there a place in an
# off-target organism where the two primers bind facing one another and 70-400 bp
# apart. The existing global scan code (screening/global_scan.py) is used AS IT IS;
# writing a second scanner would mean two different criteria.
#
# This layer's criterion DOES NOT CARRY the last-two-bases-at-the-3'-end condition,
# so it is LOOSER than the real PCR criterion: some of the hits it finds may give
# no product in reality. In a risk scan that is the safe side. An extra warning is
# better than a missed real risk.
#
# SILVA Parc is DELIBERATELY out: a non-dereplicated set of 1.3 million records adds
# nothing to the false product question over NR99 and lengthens the run by hours.
# For THE IDENTITY question it is the other way round (there Parc is REQUIRED); the
# two questions are not the same.
# -------------------------------------------------------------------------
def katman1_yerel(kok, ciftler, yaz, kontrol_dizin, parc=False, kume_ust=0):
    """MEVCUT kuresel tarama kodunu kullanir. Her kume icin ayri kontrol noktasi."""
    sys.path.insert(0, kok)
    from screening import global_scan as KT, yapilandirma as C

    adaylar = [dict(ad=c['hedef'], F=c['F'], R=c['R'], lo=URUN_ALT, hi=URUN_UST)
               for c in ciftler]
    # THE 2026-08-06 BUG FIX (D-1): only 'urun' (the TOTAL hits) used to be kept, and
    # that number was written into the '2_hedef_disi_urun' column. IT WAS WRONG:
    # kuresel_tarama also counts THE TARGET'S OWN MEMBERS (the Methanosarcina pair
    # finds 485 Methanosarcina sequences in SILVA; those are not off-target, they are
    # THE TARGET ITSELF). The MFEprimer layer already made that distinction as
    # 'ayni_boyda' / 'hedef_disi'; because the two layers were not asking THE SAME
    # question, birlestir() kept producing CONTRADICTORY. The local layer now counts
    # the expected product length (+-BOY_TOL) separately, and only 'hedef_disi' enters
    # the verdict.
    bek = {c['hedef']: int(c['urun']) for c in ciftler
           if str(c.get('urun', '')).strip().isdigit()}

    # ----------------------------------------------------------------- A2
    # THE TAXONOMIC SEPARATION (2026-08-21). This layer used to make the in-target /
    # off-target distinction by LENGTH ALONE, while layer 3 (MFEprimer, D-12) answered
    # the same question by TAXON. Measured under D-12: 95.7% of the 1,605 amplicons
    # counted as "off-target" were from INSIDE the target clade itself and differed
    # only in length. Unless two layers ask THE SAME question, birlestir() produces
    # contradictions structurally, and that is exactly what the D-1 comment warns about.
    #
    # The classification happens DURING THE SCAN and is COUNTED; the identities are not
    # stored. Measured: Bakteri_universal gives 483,098 hits, and storing their
    # identities comes to ~100 MB. A counter holds constant memory per candidate and
    # the count stays COMPLETE.
    try:
        sys.path.insert(0, kok)
        from screening import taxonomy as TX
        import verification.mfeprimer_layer as _MK
        _klad = _MK.klad_tablosu(kok)
    except Exception as e:
        TX, _klad = None, {}
        yaz(u'  WARNING: taxonomic separation is NOT POSSIBLE (%s). Only the size criterion will be used, so this layer will not ask the SAME question as layer 3.' % e)

    _klad_yok = [c['hedef'] for c in ciftler if c['hedef'] not in _klad]
    if _klad and _klad_yok:
        yaz(u'  WARNING: %d targets have no entry in hedef_klad.tsv, so NO taxonomic separation is done for them: %s'
            % (len(_klad_yok), ', '.join(_klad_yok)[:160]))

    def _siniflandirici(aday_ad, baslik, db_ad):
        if TX is None or aday_ad not in _klad:
            return 'bilinmiyor'
        alan, jetonlar, _kaynak = _klad[aday_ad]
        return TX.sinifla(baslik, db_ad, jetonlar, alan)

    toplam = {c['hedef']: dict(urun=0, ayni_boyda=0, hedef_disi=0, kume={},
                               boy={}, vurus=[], tarandi=0, atlanan=0,
                               sinif={k: 0 for k in ('a', 'ao', 'b', 'c', 'bilinmiyor')},
                               siniflandirildi=False)
              for c in ciftler}
    kumeler = list(KUMELER) + ([PARC_ISTEGE_BAGLI] if parc else [])
    if kume_ust:
        # only the N smallest sets (for a quick test; this is evidence that it RUNS, not coverage)
        var = [(e, d, a) for e, d, a in kumeler
               if os.path.exists(os.path.join(kok, 'REFERANS_DB', d))]
        var.sort(key=lambda t: os.path.getsize(os.path.join(kok, 'REFERANS_DB', t[1])))
        kumeler = var[:kume_ust]
        yaz(u'  (cluster-max=%d: scanning only %s. This is evidence that it RUNS, not coverage)'
            % (kume_ust, ', '.join(e for e, _, _ in kumeler)))
    for etiket, dosya, aciklama in kumeler:
        db = os.path.join(kok, 'REFERANS_DB', dosya)
        if not os.path.exists(db):
            yaz(u'  [%s] SKIPPED, file missing: %s' % (etiket, dosya))
            for h in toplam:
                toplam[h]['kume'][etiket] = 'dosya yok'
                toplam[h]['atlanan'] += 1
            continue
        # K-7: the pickle holds per-candidate results, and the candidate set changes on
        # every run. Without a signature the second night died with KeyError.
        import hashlib
        # THE 2026-08-10 SEQUENCE SEAL: the signature was built from the candidate NAMES
        # alone, so when a sequence changed the same signature came out and the old scan
        # was read back.
        imza = hashlib.md5(
            '|'.join(sorted('%s>%s<%s' % (a['ad'], a.get('F', ''), a.get('R', ''))
                            for a in adaylar)).encode('utf-8')).hexdigest()[:10]
        dy = os.path.join(kontrol_dizin, 'yerel_%s_%s.pkl'
                          % (''.join(ch if ch.isalnum() else '_' for ch in etiket), imza))
        t0 = time.time()
        yaz(u'  [%s] scanning (%s)...' % (etiket, aciklama))

        def ilerle(pi, kayit, gecen):
            print(u'     ... chunk %d, %d records (%s)'
                  % (pi, kayit, sure_metni(gecen)), end='\r', flush=True)
        try:
            res = KT.tara(adaylar, db=db, durum_yolu=dy, ilerle=ilerle,
                          siniflandirici=_siniflandirici)
        except TypeError:
            # An old-signature kuresel_tarama (no classifier parameter).
            # This must not turn into a SILENT fallback: that no taxonomic separation was made
            # is visible in the report (the '2_klad_ayrimi' column reads HAYIR).
            yaz(u'  WARNING: kuresel_tarama has the old signature, so there is NO taxonomic separation')
            res = KT.tara(adaylar, db=db, durum_yolu=dy)
        for h, r in res.items():
            if r.get('hata'):
                toplam[h]['kume'][etiket] = r['hata']
                continue
            toplam[h]['kume'][etiket] = r.get('urun', 0)
            toplam[h]['urun'] += r.get('urun', 0)
            toplam[h]['tarandi'] += 1
            # D-1: boy histogramini biriktir ve beklenen boy +-BOY_TOL icinde
            # kalanlari 'ayni_boyda' say. Bu, MFEprimer katmaninin kullandigi
            # olcutun AYNISI (mfe_katmani._spec_ayristir, tolerans=10).
            _b = bek.get(h)
            for _sz, _n2 in (r.get('boy') or {}).items():
                try:
                    _sz = int(_sz)
                except (TypeError, ValueError):
                    continue
                toplam[h]['boy'][_sz] = toplam[h]['boy'].get(_sz, 0) + _n2
                if _b is not None and abs(_sz - _b) <= BOY_TOL:
                    toplam[h]['ayni_boyda'] += _n2
                else:
                    toplam[h]['hedef_disi'] += _n2
            if _b is None:
                # beklenen boy bilinmiyorsa ayrim YAPILAMAZ - hukum bilinmiyor olsun
                toplam[h]['boy_ayrimi_yok'] = True
            # A2: taksonomik sayaclari veritabanlari boyunca biriktir
            _s = r.get('sinif') or {}
            for _k, _v in _s.items():
                if _k in toplam[h]['sinif']:
                    toplam[h]['sinif'][_k] += _v
            if r.get('siniflandirildi'):
                toplam[h]['siniflandirildi'] = True
            for v in (r.get('vurus') or [])[:20]:
                toplam[h]['vurus'].append((etiket,) + tuple(v))
        yaz(u'     done (%s): %s' % (sure_metni(time.time() - t0),
                                      ', '.join('%s=%s' % (h, toplam[h]['kume'][etiket])
                                                for h in list(toplam)[:4])))
    return toplam


# ---------------------------------------------------------------- katman 2
PB_URL = 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi'

# The minimum wait between consecutive submissions to NCBI (seconds). Submitting
# jobs to Primer-BLAST in quick succession leads to an IP block. The old code moved
# straight to the next pair with no wait at all WHEN THE JOB KEY COULD NOT BE
# OBTAINED, so 16 pairs were submitted back to back within seconds. Now every
# submission is spaced.
PB_GONDERIM_ARASI = 10

# -------------------------------------------------------------------------
# THE D-8 BUG FIX (2026-08-07) - TWO SEPARATE BUGS, AND EITHER ONE LEFT THE
# PRE-ORDER NCBI LAYER COMPLETELY USELESS.
#
# BUG 1 - "the job key could not be obtained" (16/16 pairs).
#   The value sent: ORGANISM='Bacteria (taxid:2) OR Archaea (taxid:2157) OR
#   Fungi (taxid:4751)'. NCBI's RAW reply (read, not guessed):
#     "Exception error: Invalid organism or taxonomy id input: 2 OR Archaea .
#      Please check the spelling and make sure it is on the suggested organism
#      list in organism input field"
#   So the ORGANISM field takes ONE organism; there is NO 'OR' syntax. The
#   Primer-BLAST page's OWN javascript (js/primerInit.js) does multiple organisms
#   like this:
#     function AddOneOrgField(e, orgName, orgVal) { ... name=\"ORGANISM\" ... }
#     function AddOrgField(e) { AddOneOrgField(e,"ORGANISM"); ... }
#     function GetOrganismURL(){ jQuery(".multiOrg").each(function(){
#         url += "&ORGANISM=" + $(this).value; }); }
#   The new fields are named "ORGANISM" too, so multiple organisms means a
#   REPEATED ORGANISM field. That CANNOT BE DONE with a Python dict (one key per
#   name); it needs a list of pairs plus urlencode(list). That is the fix.
#
# BUG 2 - "EMPTY RESULT" (9/15 pairs in the earlier run; INDEPENDENT of the
#   organism restriction).
#   The result page said "No target templates were found in selected database" and
#   listed no product at all. The cause: we were NOT SENDING most of the Primer3
#   fields. When they are not sent, the CGI reads those fields out of UNINITIALISED
#   memory. The "Search Summary" table on the raw result page showed it plainly:
#         Opt Primer size      1086305756
#         Min Tm               4.94733e-316
#         Opt Tm               2.18186e+243
#         Max Tm               0            <-- THIS IS THE FATAL ONE
#         Max Tm difference    6.95299e-310
#   With "Max Tm = 0" no product can meet the criterion, the page comes back empty,
#   and on our side that was marked not "clean" but "no data", which means the
#   layer measured nothing at all. The fix: send ALL the default values from NCBI's
#   OWN form (the primer-blast/ page, the defVal attribute). After the fix the same
#   request returned 1000 product rows for the same pair, and the Search Summary
#   read Max Tm=75, Min Tm=45, Opt Primer size=20.
#
# NOTE: the primer size and Tm limits are WIDENED DELIBERATELY. In this round we
# are not DESIGNING primers, we are testing the oligos WE ALREADY HAVE, and the
# limits are kept wide so that Primer3's design filter does not reject our fixed
# oligos.
# -------------------------------------------------------------------------
PB_VARSAYILAN = {
    # --- NCBI's own form defaults (taken from the defVal attribute) ---
    'PRIMER_NUM_RETURN': '10', 'PRIMER_MAX_DIFF_TM': '20',
    'PRIMER_ON_SPLICE_SITE': '0',
    'SPLICE_SITE_OVERLAP_5END': '7', 'SPLICE_SITE_OVERLAP_3END': '4',
    'SPLICE_SITE_OVERLAP_3END_MAX': '8',
    'MIN_INTRON_SIZE': '1000', 'MAX_INTRON_SIZE': '1000000',
    'SEARCHMODE': '0', 'MAX_TARGET_SIZE': '4000',
    # D-13 (2026-08-07, MEASURED): 1000 is NOT A LIMIT, it is the form's DEFAULT.
    # The NCBI form source: <input name="NUM_TARGETS_WITH_PRIMERS" defVal="1000">
    # ("Max targets to show (for pre-designed primers)"). There is NO upper bound
    # validation on the client side. Running the same pair against the same database:
    #     sent 3000  -> 3000 rows returned
    #     sent 8000  -> 8000 rows returned
    #     sent 20000 -> 11999 rows returned (the limit did not bind, the targets ran out)
    # Three settings trim together: NUM_TARGETS_WITH_PRIMERS, MAX_TARGET_PER_TEMPLATE
    # and HITSIZE. All three were raised. Because the page prints NO warning when it
    # truncates, the old 1000 looked like a "cap".
    'NUM_TARGETS': '20', 'NUM_TARGETS_WITH_PRIMERS': '20000',
    'MAX_TARGET_PER_TEMPLATE': '1000',
    'TOTAL_MISMATCH_IGNORE': '6', 'HITSIZE': '100000', 'EVALUE': '30000',
    'WORD_SIZE': '7', 'MAX_CANDIDATE_PRIMER': '500',
    'PRIMER_MIN_GC': '10.0', 'PRIMER_MAX_GC': '90.0',
    'GC_CLAMP': '0', 'POLYX': '5',
    'PRIMER_MAX_END_STABILITY': '9', 'PRIMER_MAX_END_GC': '5',
    'PRIMER_MAX_TEMPLATE_MISPRIMING_TH': '40.00',
    'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH': '70.00',
    'PRIMER_MAX_SELF_ANY_TH': '45.0', 'PRIMER_MAX_SELF_END_TH': '35.0',
    'PRIMER_PAIR_MAX_COMPL_ANY_TH': '45.0', 'PRIMER_PAIR_MAX_COMPL_END_TH': '35.0',
    'PRIMER_MAX_HAIRPIN_TH': '24.0',
    'PRIMER_MAX_TEMPLATE_MISPRIMING': '12.00',
    'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING': '24.00',
    'SELF_ANY': '8.00', 'SELF_END': '3.00',
    'PRIMER_PAIR_MAX_COMPL_ANY': '8.00', 'PRIMER_PAIR_MAX_COMPL_END': '3.00',
    'OVERLAP_5END': '7', 'OVERLAP_3END': '4',
    'MONO_CATIONS': '50.0', 'DIVA_CATIONS': '1.5',
    'CON_DNTPS': '0.6', 'CON_ANEAL_OLIGO': '50.0',
    'SALT_FORMULAR': '1', 'TM_METHOD': '1',
    'PRIMER_MISPRIMING_LIBRARY': 'AUTO',
    'ALLOW_NO_ORGANISM': 'NO', 'UNGAPPED_BLAST': 'on',
    'LOW_COMPLEXITY_FILTER': 'on', 'SHOW_SVIEWER': 'on',
    'SEARCH_SPECIFIC_PRIMER': 'on',
    # --- widened deliberately for FIXED OLIGO TESTING mode ---
    'PRIMER_MIN_SIZE': '15', 'PRIMER_OPT_SIZE': '20', 'PRIMER_MAX_SIZE': '30',
    'PRIMER_MIN_TM': '45.0', 'PRIMER_OPT_TM': '60.0', 'PRIMER_MAX_TM': '75.0',
    # --- bos birakilanlar (formda da bos) ---
    'PRIMER5_START': '', 'PRIMER5_END': '', 'PRIMER3_START': '', 'PRIMER3_END': '',
    'ENTREZ_QUERY': '', 'PRODUCT_MIN_TM': '', 'PRODUCT_OPT_TM': '', 'PRODUCT_MAX_TM': '',
    'CMD': 'request',
}


def pb_ac(url, veri=None, deneme=4, timeout=90, yaz=None):
    """A request to Primer-BLAST. It retries on transient network and rate limit errors.

        Under load NCBI can close the connection without answering
        (RemoteDisconnected) or return 429/502/503. The old code wrote that pair off
        as FAILED on the first error and moved on; that is exactly what happened in
        testing. The wait doubles on every attempt (10, 20, 40 s).

    """
    import urllib.request, urllib.error
    son = None
    for i in range(deneme):
        try:
            req = urllib.request.Request(
                url, veri, headers={'User-Agent': 'PrimerJury-primer-QC/1.0'})
            with urllib.request.urlopen(req, timeout=timeout) as f:
                return f.read().decode('utf-8', 'replace')
        except Exception as e:
            son = e
            kod = getattr(e, 'code', None)
            # On PERMANENT errors such as 400/404 there is no point retrying
            if kod is not None and kod not in (429, 500, 502, 503, 504):
                raise
            if i == deneme - 1:
                break
            b = PB_GONDERIM_ARASI * (2 ** i)
            if yaz:
                yaz(u'    network error (%s), retrying in %d s (%d/%d)'
                    % (type(e).__name__, b, i + 2, deneme))
            time.sleep(b)
    raise son


def organizma_listesi(s):
    """Splits the '--organizma' string into SEPARATE organisms.

        Primer-BLAST's ORGANISM field takes ONE organism (see the D-8 note above).
        The 'A OR B OR C' form the user writes is accepted, but it is NOT SENT to the
        request in that form: it is broken up here and turned into REPEATED ORGANISM
        fields. ';' and a newline count as separators too.

    """
    if not s or not s.strip():
        return []
    parcalar = re.split(r'\s+OR\s+|;|\n', s, flags=re.I)
    return [x.strip() for x in parcalar if x.strip()]


# LAYER 4 - NCBI. blastn -remote IS NOT USED: it exceeds the 45 second process cap
# and the job was being cut off half way. The Primer-BLAST URL API is used instead:
# the job is submitted, a job key is obtained, and it is polled until it finishes.
# Every reply is written to disk raw, so that when the numbers are argued about the
# source text can be looked at.

# D-13b (2026-08-07): on the Primer-BLAST result page the products are split into
# <div class="prPairTl">HEADER</div> ... blocks. Counting the whole page puts "a
# product on the target" and "an off-target product" in the same bucket.
# MEASURED 2026-08-10 (a live single pair trial, Proteolitik_Cloacimonas):
# even when the target's own taxon is excluded with ENTREZ_QUERY, Primer-BLAST does
# NOT OPEN the "potentially unintended templates" section. Unless a template
# sequence is declared, it puts every product it finds under "target templates".
# The "Products on potentially unintended templates" text on the page is not a
# SECTION HEADING but the LINK LIST at the top of the page, and it is printed for
# an empty section too. So a verdict based on counting sections CANNOT WORK
# STRUCTURALLY.
#
# The measured numbers (Proteolitik_Cloacimonas, txid112 excluded):
#   42 products -> 40 of them "uncultured bacterium clone ...", 2 "Methanogenic
#   prokaryote enrichment culture ..." -> ZERO named taxa.
# Unnamed environmental clones belong to no taxon, so ENTREZ_QUERY cannot exclude
# them; and since our targets are themselves unnamed lineages, most of those clones
# could BE THE TARGET. The label cannot decide it.
#
# So the verdict now looks at THE HEADER rather than the section: a product with a
# NAMED identity is evidence of an off-target, while an unnamed environmental clone
# goes into the "cannot decide" column and its identity is settled by sequence
# comparison (layers 2 and 3).
_ADSIZ_IZLERI = (u'uncultured', u'unidentified', u'unclassified', u'metagenome',
                 u'environmental sample', u'enrichment culture', u'clone',
                 u'synthetic construct')


def _ncbi_urunleri(html):
    """Returns every product on the result page as (accession, header, product_length)."""
    d = re.sub(r'<[^>]+>', ' ', html)
    d = re.sub(r'&nbsp;?', ' ', d)
    d = re.sub(r'\s+', ' ', d)
    return [(m.group(1), m.group(2).strip(), int(m.group(3)))
            for m in re.finditer(
                r'>\s*([A-Z]{1,2}[_A-Z]*\d{5,}\.\d)\s+([^>]{5,200}?)\s+'
                r'product length\s*=\s*(\d+)', d)]


def _adsiz_mi(baslik):
    """Is the record UNNAMED? A SINGLE DEFINITION: ncbi_yeniden_siniflandir.adli_mi().

        The overnight fix of 2026-08-10. The first rule looked only for keywords
        (uncultured, clone, metagenome and so on) and counted these headers as NAMED:
            "Bacterium LC2012 16S ribosomal RNA gene"
            "Archaeon 2022-TM-MRBT1 gene for 16S rRNA"
            "anaerobic methanogenic archaeon E15-5 16S rRNA gene"
            "Environmental 16s rDNA sequence from Evry wastewater treatment plant"
        None of them has a genus name. That is what inflated the off-target counts:
        Bacteroidales 650 -> 82, Nitrosocosmicus 170 -> 9, the Methanothrix genus
        22 -> 1. The new rule does not look for keywords, it looks for A NAME.

    """
    try:
        from ncbi_reclassify import adli_mi as _adli
    except ImportError:
        b = baslik.lower()
        return any(iz in b for iz in _ADSIZ_IZLERI)
    return not _adli(baslik)


def _ncbi_bolum_say(html, baslik):
    """Counts the 'product length' lines following the given section heading."""
    par = re.split(r'<div class="prPairTl">(.*?)</div>', html)
    t = 0
    for i in range(1, len(par), 2):
        b = re.sub(r'<[^>]+>', '', par[i]).strip().lower()
        if baslik.lower() in b:
            t += len(re.findall(r'product length\s*=\s*\d+', par[i + 1], re.I))
    return t


def katman2_oto(ciftler, cikti, yaz, organizma='', bekleme=20, tur_ust=60,
                haric_taxid=''):
    """The NCBI Primer-BLAST URL API. blastn -remote IS NOT USED (it exceeds the 45 s cap).

        Submit -> get a job key -> poll until it finishes. Every reply is written to
        disk raw; on a machine with no network, or when NCBI is queued, it gives up
        CLEANLY and says to fall back to the manual route.

    """
    import urllib.request, urllib.parse
    ham = os.path.join(cikti, 'ncbi_ham')
    os.makedirs(ham, exist_ok=True)
    orgs = organizma_listesi(organizma)
    if orgs:
        yaz(u'  organism restriction (sent as %d separate ORGANISM fields): %s'
            % (len(orgs), ' | '.join(orgs)))
    else:
        yaz(u'  NO organism restriction, the whole of nt will be scanned (broad targets are likely to hit the result cap)')
    # D-13c (2026-08-07, MEASURED): a bare 'NOT txidN[Organism]' INVERTS the filter. It
    # returns ONLY that taxon instead of EXCLUDING it. The measured evidence (the same
    # pair, the same database, the genus distribution):
    #   ENTREZ_QUERY empty                    -> Escherichia 45, Bacillus 34,
    #                                            Pseudomonas 32, Staphylococcus 25
    #   'NOT txid1279[Organism]'              -> ONLY Staphylococcus 196  (INVERTED!)
    #   'all[filter] NOT txid1279[Organism]'  -> Escherichia 14, Bacillus 13,
    #                                            ..., Staphylococcus 1  (CORRECT)
    # So the prefix is enforced IN THE CODE; even if the user writes a bare NOT, it is
    # corrected.
    # 2026-08-10 PER-TARGET EXCLUSION. A single general taxid was being applied to the
    # whole run, when each target's OWN taxon has to be excluded. Under a general
    # organism restriction (Bacteria OR Archaea OR Fungi) Primer-BLAST puts every
    # product it finds under "target templates" and leaves the "unintended" section
    # empty; measured on 22 of 22 pages. The fix: exclude the target's own taxon with
    # ENTREZ_QUERY, and then every remaining product is off-target by definition. The
    # map file: screening/hedef_taxid.tsv
    def _ent_of(_tx):
        # Several taxids can be given comma separated; each becomes a SEPARATE NOT term.
        # For universal primers (Metanojen_universal, for example) the target is not one
        # taxon but the union of several orders, and one term is not enough.
        _ler = [x.strip().lstrip('txid') for x in str(_tx).split(',') if x.strip()]
        if not _ler:
            return ''
        return 'all[filter]' + ''.join(' NOT txid%s[Organism]' % x for x in _ler)

    HARITA = {}
    _hy = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'screening', 'hedef_taxid.tsv')
    if os.path.exists(_hy):
        for _l in open(_hy, encoding='utf-8'):
            _l = _l.rstrip('\n')
            if not _l.strip() or _l.startswith('#'):
                continue
            _p = _l.split('\t')
            if len(_p) >= 2 and _p[1].strip():
                HARITA[_p[0].strip()] = _p[1].strip()
        yaz(u'  per-target exclusion map loaded: %d targets (%s)'
            % (len(HARITA), os.path.basename(_hy)))
    ent = ''
    if haric_taxid:
        ent = _ent_of(haric_taxid)
        yaz(u'  ENTREZ_QUERY (GLOBAL exclusion): %s' % ent)
    if HARITA:
        yaz(u'    NOTE: the "all[filter]" prefix is MANDATORY. Without it NCBI inverts the filter and returns only that taxon (measured).')
    out = {}
    ilk = True
    for c in ciftler:
        ad = c['hedef']
        if not ilk:
            time.sleep(PB_GONDERIM_ARASI)      # NCBI hiz siniri - bkz. PB_GONDERIM_ARASI
        ilk = False
        p = dict(PB_VARSAYILAN)
        p.update(dict(PRIMER_LEFT_INPUT=c['F'], PRIMER_RIGHT_INPUT=c['R'],
                      PRIMER_PRODUCT_MIN=str(URUN_ALT), PRIMER_PRODUCT_MAX=str(URUN_UST),
                      PRIMER_SPECIFICITY_DATABASE='nt',
                      TOTAL_PRIMER_SPECIFICITY_MISMATCH='5',
                      PRIMER_3END_SPECIFICITY_MISMATCH='2',
                      MISMATCH_REGION_LENGTH='5'))
        _ent_c = _ent_of(HARITA[ad]) if ad in HARITA else ent
        _kendi_dislandi = ad in HARITA
        if _ent_c:
            p['ENTREZ_QUERY'] = _ent_c
            if ad in HARITA:
                yaz(u'  [%s] its own taxon was excluded: txid%s' % (ad[:40], HARITA[ad]))
        # ORGANISM IS NOT PUT IN THE DICT: multiple organisms means a REPEATED field, and a
        # dict holds one key. It is added as a list of pairs instead.
        alanlar = list(p.items()) + [('ORGANISM', o) for o in orgs]
        try:
            veri = urllib.parse.urlencode(alanlar).encode()
            s = pb_ac(PB_URL, veri, yaz=yaz)
            m = re.search(r'job_key=([A-Za-z0-9_\-]+)', s)
            if not m:
                # DO NOT GUESS THE REASON - read NCBI's own error text and write it.
                _h = re.search(r'(?:Exception error|Error)\s*:\s*([^<\n]{5,300})', s, re.I)
                _sebep = _h.group(1).strip() if _h else u'NCBI hata metni bulunamadi'
                open(os.path.join(ham, '%s_ANAHTARSIZ.html' % re.sub(r'\W+', '_', ad)),
                     'w', encoding='utf-8').write(s)
                out[ad] = dict(durum='BASARISIZ',
                               not_=u'is anahtari alinamadi - NCBI yaniti: %s' % _sebep)
                yaz(u'  [%s] NCBI: is anahtari alinamadi - NCBI diyor ki: %s'
                    % (ad, _sebep)); continue
            anahtar = m.group(1)
            yaz(u'  [%s] NCBI isi gonderildi (%s), bekleniyor...' % (ad, anahtar))
            son = ''
            for i in range(tur_ust):
                time.sleep(bekleme)
                u2 = PB_URL + '?job_key=' + anahtar
                son = pb_ac(u2, yaz=yaz)
                # THE 2026-08-06 BUG FIX - caught in the pre-order measurement.
                # The old condition looked only for 'still running' and 'please wait'.
                # Primer-BLAST's waiting page CONTAINS NEITHER string; the page says
                # "Status Running" and "Time since submission". The result: the loop broke
                # on the FIRST poll, 'product length' could not be found on the page of a job
                # that was still RUNNING, and the row was reported as "0 off-target", that is,
                # CLEAN. A false assurance before an order is exactly what we have to avoid.
                # Completion is now confirmed POSITIVELY.
                _d = son.lower()
                _kosuyor = ('still running' in _d or 'please wait' in _d
                            or 'status</th><td>running' in _d.replace(' ', '')
                            or re.search(r'status[^<]*<[^>]*>\s*running', _d)
                            or 'time since submission' in _d)
                _bitti = bool(re.search(r'primer pair \d|no significant|not find any target'
                                        r'|unintended template', _d))
                if _bitti and not _kosuyor:
                    break
                print('     ... %d. yoklama (%s)          '
                      % (i + 1, 'kosuyor' if _kosuyor else 'bekleniyor'),
                      end='\r', flush=True)
            open(os.path.join(ham, '%s.html' % re.sub(r'\W+', '_', ad)), 'w',
                 encoding='utf-8').write(son)
            _d = son.lower()
            _bitti = bool(re.search(r'primer pair \d|no significant|not find any target'
                                    r'|unintended template', _d))
            if not _bitti:
                # Is bitmeden tavana carptik. TEMIZ demek YASAK - bilinmiyor denir.
                out[ad] = dict(durum='BASARISIZ',
                               not_=u'NCBI isi %d yoklamada bitmedi (hala kuyrukta). '
                                    u'Sonuc BILINMIYOR - temiz sayilmadi.' % tur_ust)
                yaz(u'  [%s] NCBI: the job did not finish (queued). Fall back to the manual route.' % ad)
                continue
            n = len(re.findall(r'product length\s*=\s*\d+', son, re.I))
            hedefsiz = bool(re.search(r'no significant|not find any target', son, re.I))
            # D-13b (2026-08-07, MEASURED): the page puts products in SEPARATE
            # SECTIONS ('Products on intended targets', '... on potentially
            # unintended templates', '... on target templates'). Counting the whole
            # page destroys that distinction. In 22 of the 22 raw HTML pages from the
            # 2026-08-07 run the "potentially unintended templates" section was EMPTY
            # and every product sat under "target templates", which means max(0,n-1)
            # IS NOT AN OFF-TARGET COUNT but the TOTAL number of products found in nt.
            # On top of that, because the organism restriction was Bacteria/Archaea/
            # Fungi, the target's own members are in that list too. A measured example:
            # Proteiniphilum_cinsi 876 products, 110 of the headers "Proteiniphilum"
            # and 760 unnameable environmental clones ("uncultured bacterium clone
            # ..."), so the NCBI headers DO NOT CARRY the answer to this question.
            n_unint = _ncbi_bolum_say(son, 'potentially unintended templates')
            n_target = _ncbi_bolum_say(son, 'target templates')
            # THE D-3 BUG FIX (2026-08-06): 'n' was producing a false value in two
            # separate ways, and both of them counted as FINE.
            #   (a) THE CAP: the Primer-BLAST result page lists at most 1000 products.
            #       If n==1000 the real number is 1000 OR MORE, and writing 999 with
            #       max(0,n-1) is not a COUNT, it is the sign of hitting the cap.
            #       In that run exactly five targets came out at 999 (Metanomikrobiyales,
            #       Nitrosocosmicus, Microascaceae, Metanojen_universal, Mantar F2), all
            #       of them broad targets, all of them at the cap.
            #   (b) EMPTY: on a page where the 'Products on ...' sections come back EMPTY,
            #       n==0 and max(0, 0-1)==0, so it was reported CLEAN. A page that lists
            #       no product at all, not even the intended one, is not 'clean', it is
            #       'no data'. That happened on nine targets in that run.
            _ur = _ncbi_urunleri(son)
            _adli = [(a, b, L) for a, b, L in _ur if not _adsiz_mi(b)]
            _adsiz = [(a, b, L) for a, b, L in _ur if _adsiz_mi(b)]
            if not hedefsiz and n >= NCBI_SONUC_TAVANI:
                # The page is truncated. But if the truncated list holds a NAMED
                # off-target taxon, that is a LOWER BOUND and it is valid: "at least this
                # many" can be said. If it comes out zero, nothing can be said; the list is
                # incomplete.
                if _kendi_dislandi and _adli:
                    out[ad] = dict(durum='TAMAM (alt sinir)', hedef_disi=len(_adli),
                                   ncbi_toplam_urun=n, ncbi_adsiz_klon=len(_adsiz),
                                   ncbi_ornek=u'; '.join(b[:70] for _a, b, _L in _adli[:3]),
                                   not_=u'Sayfa tavana carpti (%d urun, gercek sayi daha '
                                        u'fazla). Hedefin kendi taksonu (txid%s) dislandi. '
                                        u'Kirpilmis listede ADLI hedef disi takson: %d - bu '
                                        u'bir ALT SINIRDIR, kesin sayi degildir. Adsiz cevre '
                                        u'klonu: %d (etiketten karar verilemez, katman 2-3 '
                                        u'karar verir).'
                                        % (n, HARITA.get(ad, '?'), len(_adli), len(_adsiz)))
                    yaz(u'  [%s] NCBI: hit the cap (%d) but NAMED off-target >= %d (a lower bound)'
                        % (ad, n, len(_adli)))
                    continue
                out[ad] = dict(durum='BASARISIZ - SONUC TAVANI',
                               ncbi_toplam_urun=n, ncbi_adsiz_klon=len(_adsiz),
                               not_=u'Primer-BLAST %d urun listeledi (sayfa tavani). '
                                    u'Gercek sayi >= %d. Bu bir SAYIM DEGIL; hukum '
                                    u'icin kullanilamaz. Organizma kisiti (--organizma) '
                                    u'ile daraltip yeniden kosun.'
                                    % (n, NCBI_SONUC_TAVANI))
                yaz(u'  [%s] NCBI: RESULT CAP (%d) - not a count, not tested' % (ad, n))
                continue
            if not hedefsiz and n == 0:
                out[ad] = dict(durum='BASARISIZ - BOS SONUC',
                               not_=u'Sayfa bitti ama hicbir "product length" satiri yok '
                                    u'(hedefteki urun bile listelenmemis). Bu TEMIZ degil, '
                                    u'VERI YOK. Sinanmadi sayilir.')
                yaz(u'  [%s] NCBI: EMPTY result page - not tested' % ad)
                continue
            # IF THERE IS NO ORGANISM RESTRICTION (--organizma empty), Primer-BLAST also
            # lists the target's OWN members under "unintended template"; max(0,n-1) is
            # only correct under the assumption that "there is exactly one intended
            # product", and for group and universal primers that assumption IS INVALID.
            # Mark this openly.
            # D-13b: the value that enters the verdict is now the count of the
            # "unintended templates" SECTION. If that section is empty and every product
            # sits under "target templates", the page IS NOT ANSWERING this question
            # (because no target template was declared, Primer-BLAST counts no product as
            # 'unintended'). In that case the layer DOES NOT VOTE; it does not count as
            # 'clean'.
            if hedefsiz:
                out[ad] = dict(durum='TAMAM', hedef_disi=0, ncbi_toplam_urun=0,
                               not_=u'Primer-BLAST hic urun bulamadi.')
                yaz(u'  [%s] NCBI: no products at all -> off-target 0' % ad)
                continue
            if n_unint == 0 and n_target > 0:
                if not _kendi_dislandi:
                    out[ad] = dict(
                        durum='BASARISIZ - DISLAMA HARITASINDA YOK',
                        ncbi_toplam_urun=n_target, ncbi_adsiz_klon=len(_adsiz),
                        not_=u'Sayfa %d urun listeledi ama bu hedef icin '
                             u'screening/hedef_taxid.tsv icinde dislanacak takson '
                             u'yazili degil. Hedefin kendi uyeleri de listede olabilir, '
                             u'ayirt edilemez. SINANMADI.' % n_target)
                    yaz(u'  [%s] NCBI: not in the exclusion map, not tested' % ad)
                    continue
                # Its own taxon was excluded. Even if the section heading never opens, every
                # remaining NAMED product is off-target by definition.
                out[ad] = dict(
                    durum='TAMAM', hedef_disi=len(_adli),
                    ncbi_toplam_urun=n_target, ncbi_adsiz_klon=len(_adsiz),
                    ncbi_ornek=u'; '.join(b[:70] for _a, b, _L in _adli[:3]),
                    not_=u'Hedefin kendi taksonu (txid%s) ENTREZ_QUERY ile dislandi. '
                         u'%d urunun %d tanesi ADLI takson (hedef disi kaniti), %d '
                         u'tanesi adsiz cevre klonu ("uncultured ...") - adsizlar '
                         u'hicbir taksona bagli olmadigi icin dislama suzgeci onlara '
                         u'islemez ve hedefin KENDISI olabilirler; kimliklerine dizi '
                         u'karsilastirmasi (katman 2-3) karar verir, hukme girmezler. '
                         u'Bolum basligina bakilmadi - sablon dizi bildirilmedigi surece '
                         u'Primer-BLAST "unintended" bolumunu hic acmiyor (olculdu).'
                         % (HARITA.get(ad, '?'), n_target, len(_adli), len(_adsiz)))
                yaz(u'  [%s] NCBI: named off-target %d / unnamed clones %d / total %d'
                    % (ad, len(_adli), len(_adsiz), n_target))
                continue
            _kusur = (u'ORGANIZMA KISITI YOK: hedefin kendi uyeleri de "unintended" '
                      u'altinda sayilmis olabilir. ' if not organizma else u'')
            out[ad] = dict(durum='TAMAM', hedef_disi=n_unint,
                           ncbi_toplam_urun=n_target,
                           not_=_kusur + u'"unintended templates" bolumunun sayimi; '
                                u'ham yanit ncbi_ham/ altinda')
            yaz(u'  [%s] NCBI: off-target (unintended section) %s / total products %s%s'
                % (ad, n_unint, n_target, u' (no organism restriction)' if _kusur else u''))
        except Exception as e:
            out[ad] = dict(durum='BASARISIZ', not_=u'%s: %s' % (type(e).__name__, e))
            yaz(u'  [%s] NCBI FAILED (%s), fall back to the manual route' % (ad, type(e).__name__))
    return out


def katman2_elle_girdi(ciftler, cikti, yaz, organizma=''):
    """A ready input plus a result template, so the user can run Primer-BLAST in a browser."""
    g = os.path.join(cikti, 'NCBI_PRIMER_BLAST_GIRDI.tsv')
    with open(g, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# READY-MADE INPUT for NCBI Primer-BLAST, paste it straight in.\n')
        fh.write(u'# Address: https://www.ncbi.nlm.nih.gov/tools/primer-blast/\n')
        fh.write(u'# One pair per row. Paste into these fields on the page:\n')
        fh.write(u'#   "Primer Parameters > Forward primer"  <- the F column\n')
        fh.write(u'#   "Primer Parameters > Reverse primer"  <- the R column\n')
        fh.write(u'#   "Exon/intron selection > PCR product size"  Min/Max <- urun_min / urun_max\n')
        fh.write(u'#   "Primer Pair Specificity Checking Parameters":\n')
        fh.write(u'#       Database = nt ; Organism = the organizma_kisiti column\n')
        fh.write(u'#       Total mismatches = 5 ; 3\' end mismatches = 2\n')
        fh.write(u'# Write the results into NCBI_SONUC_SABLONU.tsv and load them back with --ncbi-load.\n')
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'F', 'R', 'urun_min', 'urun_max', 'organizma_kisiti', 'not'])
        for c in ciftler:
            w.writerow([c['hedef'], c['F'], c['R'], URUN_ALT, URUN_UST,
                        organizma or '(bos = tum nt)', c['tur']])
    yaz(u'  written: %s' % g)

    s = os.path.join(cikti, 'NCBI_SONUC_SABLONU.tsv')
    with open(s, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# Write the NCBI results HERE, then:\n')
        fh.write(u'#   verification/full_chain.py -> (D) -> "load manual results", or\n')
        fh.write(u'#   python3 verification/specificity_round.py --root . --ncbi-load DOGRULAMA_SONUC/NCBI_SONUC_SABLONU.tsv\n')
        fh.write(u'# hedef_disi_urun_sayisi: how many products Primer-BLAST counts under "Products on potentially unintended templates".\n')
        fh.write(u'# Write 0 if there are none. Leave it empty if you did not look; that row counts as "NCBI not done".\n')
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'hedef_disi_urun_sayisi', 'en_yakin_hedef_disi_organizma', 'notunuz'])
        for c in ciftler:
            w.writerow([c['hedef'], '', '', ''])
    yaz(u'  written: %s' % s)
    return g, s


def ncbi_yukle(yol, yaz=None):
    """Reads the hand filled NCBI template.

        THE A5 FIX (2026-08-21): a malformed number field used to be silently
        'passed'. When a person filling in the template wrote '~3' or '3 (maybe)',
        the NCBI layer for that target was NEVER FORMED and the verdict table showed
        'BILINMIYOR'. The user thought they had entered a value while the layer had
        dropped. Leaving a field EMPTY and writing something MALFORMED are not the
        same thing: the first means "I did not look" and is legitimate, the second is
        a typing mistake and has to be visible.

    """
    out = {}
    bozuk = []
    with open(yol, encoding='utf-8') as fh:
        for r in csv.DictReader((s for s in fh if not s.startswith('#')), delimiter='\t'):
            n = (r.get('hedef_disi_urun_sayisi') or '').strip()
            if n == '':
                continue                      # bilerek bos: "bakmadim", mesru
            try:
                out[r['hedef'].strip()] = dict(durum='TAMAM (elle)', hedef_disi=int(n),
                                               not_=(r.get('en_yakin_hedef_disi_organizma') or ''))
            except ValueError:
                bozuk.append((r.get('hedef', '?').strip(), n))
    if bozuk and yaz:
        yaz(u'  WARNING: "hedef_disi_urun_sayisi" is NOT a number on %d rows. The NCBI layer was NOT built for those targets and the verdict will show BILINMIYOR (unknown):' % len(bozuk))
        for h, n in bozuk:
            yaz(u'    %-40s value: %r  (an integer is expected; write 0 if there are none, leave EMPTY if you did not look)' % (h[:40], n))
    return out


# ---------------------------------------------------------------- combine
# Reduces a layer's raw count to three states. "BILINMIYOR" IS A SEPARATE STATE: an
# unmeasured layer does not count as "clean", or layers that never ran would quietly
# vote in favour.
# THE C-1 BUG FIX (2026-08-07): the summary line was unreadable.
#   "INCELEME - gevsek olcut vurusu (1 adet); 3' son iki baz sinanmali: 1
#    KOSULLU ...: 6   RISKLI ...: 3
#    INCELEME - gevsek olcut vurusu (11 adet)...: 1
#    INCELEME - gevsek olcut vurusu (22 adet)...: 1 ..."
# The cause: the counter was using THE WHOLE VERDICT STRING as its key. Because the
# hit COUNT appears inside the string ("(11 adet)", "(22 adet)"), every row in the
# same category became A SEPARATE key; 16 pairs spread across 10+ pseudo-categories
# and the summary summarised nothing. The category is now SEPARATED FROM THE NUMBER:
# the first word of the verdict string is the category, and the rest stands on ITS
# OWN LINE as detail.
KATEGORILER = ('KESIN', 'KOSULLU', 'INCELEME', 'RISKLI', 'CELISKILI', 'EKSIK')
# These four are written even when they are ZERO: the reader has to see the
# difference between "no INCELEME appears" and "INCELEME is zero".
ANA_KATEGORILER = ('KESIN', 'KOSULLU', 'INCELEME', 'RISKLI')


def karar_kategorisi(karar):
    """Hukum dizgesinden SAYISIZ kategori anahtarini cikarir."""
    k = (karar or '').strip()
    # D-16 (2026-08-07): maxsplit was passed as a POSITIONAL argument; from
    # Python 3.13 that is a DeprecationWarning, and later an error.
    ilk = re.split(r'[\s\-:(]', k, maxsplit=1)[0].upper()
    return ilk if ilk in KATEGORILER else (k.split(' ')[0].upper() or 'BILINMIYOR')


# --------------------------------------------------------------- the D-18 measure
# The measure that REPLACES the yes/no question "do the last two bases at the 3'
# end match exactly".
#
# WHY IT CHANGED. The old question was not a measurement: MFEprimer 3.0 and later
# do not allow a mismatch at the 3' terminal base by definition [K16], so the
# terminal base matching in 69/69 records was a forced output of the algorithm.
# On top of that, the yes/no question systematically discarded terminal mismatches
# that the literature has MEASURED as amplifying EFFICIENTLY (terminal T in
# particular) [K13].
#
# THE NEW MEASURE. Not mismatch present or absent, but THE EXPECTED CYCLE PENALTY.
# A binding region counts as a real competitor only when its expected penalty is
# smaller than the dCq REQUIRED for that target. The measure then speaks in the
# same currency as the abundance sensitive threshold (cycles), and the two can be
# set against one another.
#
# THE NUMBERS COME FROM THE LITERATURE, they are not invented
# (LITERATUR_2026-08-07.md, section 3):
#   [K13] Kwok 1990 : 3' terminal A:G, G:A, C:C  -> ~100 fold  = 6.6 cycles
#                     3' terminal A:A            -> ~20 fold   = 4.3 cycles
#                     other 3' terminal mismatches amplified EFFICIENTLY = 0
#                     terminal T (T:G, T:C, T:T) is the least consequential; even
#                     with an adjacent mismatch, meaningful amplification = 0
#   [K14] Bru 2008  : the position ONE BEFORE the 3' end -> ~3 log = 10 cycles
#                     positions 5, 6 and 8 from the 3' end -> ~1 log = 3.3 cycles
#                     further than 4 bases from the 3' end on the reverse primer
#                     -> no effect = 0
#   [K17] Sozhamannan 2025 : THREE mismatches within the last FOUR bases -> >15 cycles
#
# THE UNMEASURED RANGE IS MARKED OPENLY: there is NO direct measurement for
# positions 3 and 4; Bru's 5-8 value is used as a LOWER BOUND and 'tahmin' is
# returned.
#
# THIS FUNCTION DOES NOT DECIDE ON ITS OWN. Its output is a cycle estimate; the
# verdict comes from comparing it against the required dCq = log2(R) + 4.3.

UC3_TERMINAL_CEZA = {              # (primer_bazi, sablon_bazi) -> dongu
    ('A', 'G'): 6.6, ('G', 'A'): 6.6, ('C', 'C'): 6.6,
    ('A', 'A'): 4.3,
    ('T', 'G'): 0.0, ('T', 'C'): 0.0, ('T', 'T'): 0.0,
}
UC3_TERMINAL_VARSAYILAN = 0.0      # Kwok: geri kalani verimli cogaldi


def uc3_ceza_dongu(uyumsuz_konumlar, terminal_ciftler=None):
    """Returns the EXPECTED cycle penalty of mismatches near the 3' end.

        uyumsuz_konumlar : the list of 1-based positions from the 3' end (1 = the
                           terminal base).
        terminal_ciftler : {position: (primer_base, template_base)} - used only for
                           position 1; if it is not given, Kwok's worst value
                           (6.6 cycles) is taken, which is the CAUTIOUS side.

        Returns: (penalty_cycles, the_basis_text, was_it_measured)

    """
    if not uyumsuz_konumlar:
        return (0.0, 'uyumsuzluk yok', True)
    kon = sorted(set(int(x) for x in uyumsuz_konumlar))
    tahmin = False
    parcalar = []
    toplam = 0.0

    # Sozhamannan: son dort baz icinde uc ve uzeri uyumsuzluk
    son4 = [k for k in kon if k <= 4]
    if len(son4) >= 3:
        return (15.0, 'son 4 baz icinde %d uyumsuzluk -> >15 dongu [K17]'
                % len(son4), True)

    for k in kon:
        if k == 1:
            cift = (terminal_ciftler or {}).get(1)
            if cift:
                c = UC3_TERMINAL_CEZA.get(tuple(cift), UC3_TERMINAL_VARSAYILAN)
                parcalar.append('terminal %s:%s = %.1f [K13]' % (cift[0], cift[1], c))
            else:
                c = 6.6
                parcalar.append('terminal (baz cifti bilinmiyor, en kotu) = 6,6 [K13]')
            toplam += c
        elif k == 2:
            toplam += 10.0
            parcalar.append('sondan ikinci = 10,0 [K14]')
        elif k in (3, 4):
            toplam += 3.3
            tahmin = True
            parcalar.append('konum %d = 3,3 (OLCULMEDI, 5-8 degerinden alt sinir) [K14]' % k)
        elif k <= 8:
            toplam += 3.3
            parcalar.append('konum %d = 3,3 [K14]' % k)
        else:
            parcalar.append('konum %d = 0,0 (3\' uctan uzak) [K13, K14]' % k)
    return (toplam, ' + '.join(parcalar), not tahmin)


def hukum(v):
    """Bir katmanin sonucunu TEMIZ / RISKLI / BILINMIYOR'a indirger."""
    if v is None:
        return 'BILINMIYOR'
    if isinstance(v, str):
        return 'BILINMIYOR'
    return 'TEMIZ' if v <= VURUS_ESIGI else 'RISKLI'


# -------------------------------------------------------------------------
# FOUR SOURCES SIDE BY SIDE. The decision is tied to the sources AGREEING, not to
# any one source's own number:
#   all four sources CLEAN            -> KESIN
#   three sources CLEAN, one missing  -> KOSULLU
#   the sources DISAGREE              -> CELISKILI (cannot be ordered)
#   none of them gave a result        -> EKSIK
# "Zero databases were scanned" and "all of them were scanned and came out clean"
# are not the same thing, which is why the local layer's value is accepted as a
# number only if a scan really happened (tarandi > 0); otherwise it is BILINMIYOR.
# -------------------------------------------------------------------------
def birlestir(ciftler, yerel, ncbi, mfe=None, klad=None):
    """THREE MEASUREMENT layers side by side (local DB / MFEprimer / NCBI).

        The sample layer DELIBERATELY does not vote: it produces a constant 'TEMIZ',
        so it is a tautology rather than a measurement (the D-2 fix, 2026-08-06).
        If all three layers are TEMIZ and all three measured, KESIN; if two measured,
        KOSULLU; if the layers that measured DISAGREE, CELISKILI.

    """
    out = []
    for c in ciftler:
        h = c['hedef']
        # THE D-2 BUG FIX (2026-08-06): n_ok was a CONSTANT 'TEMIZ' and it still joined
        # the set of 'known' votes. That is a tautology: the sample layer is not a
        # MEASUREMENT, it is the acceptance criterion that PUT the pair on this list, and
        # by definition it can never say RISKLI. A constant TEMIZ vote, once put into the
        # agreement test, made a SINGLE RISKLI reading from any layer turn the set into
        # {TEMIZ,RISKLI}, and the row came out NECESSARILY CONTRADICTORY. All 16 of the 16
        # pairs were marked contradictory for that reason. The sample is now shown AS A
        # COLUMN but IT DOES NOT VOTE.
        n_ok = 'TEMIZ'
        yv = yerel.get(h, {}) if yerel else {}
        # O-1: 'zero databases were scanned' is not 'all scanned, clean'
        # D-1: what enters the verdict is not the TOTAL hits but the ones DIFFERING from
        # the expected length.
        y_tum = yv.get('urun') if yv.get('tarandi') else None
        if yv.get('tarandi') and not yv.get('boy_ayrimi_yok'):
            y = yv.get('hedef_disi')
        else:
            y = None            # boy ayrimi yapilamadiysa bu katman hukum vermez
        y_ok = hukum(y)
        nb = ncbi.get(h) if ncbi else None
        nb_v = nb.get('hedef_disi') if nb and nb.get('durum', '').startswith('TAMAM') else None
        nb_ok = hukum(nb_v)
        mv = (mfe or {}).get(h)
        m_ham = mv.get('hedef_disi') if mv else None
        # D-12 (2026-08-07): THE MEASURE ENTERING THE VERDICT CHANGED. MFEprimer's
        # "off-target" count rests ON LENGTH ALONE; for universal and group primers the
        # target clade's own members also give amplicons of a different length. THE
        # MEASUREMENT (the 2026-08-07 run, 1605 amplicons, the taxonomy strings in
        # mfe_hedef_disi_kimlikleri.tsv):
        #   (a) inside the target clade, different length  1536  (b) same domain, outside the clade  24
        #  (ao) inside the target domain but an ORGANELLE     31  (c) a different domain               14
        # So 95.7% of the raw count is a harmless length variant. What enters the verdict
        # now is klad_disi = (b)+(c); the raw count STAYS as a column.
        kl = (klad or {}).get(h)
        if kl is not None:
            m_urun = kl['klad_disi']
            m_klad_ayrimi = True
        else:
            m_urun = m_ham
            m_klad_ayrimi = False
        m_ok = hukum(m_urun) if mv else 'BILINMIYOR'
        # D-2: the sample is DELIBERATELY out; a constant value cannot vote.
        kaynaklar = dict(yerel=y_ok, mfeprimer=m_ok, ncbi=nb_ok)
        bilinen = [x for x in kaynaklar.values() if x != 'BILINMIYOR']
        uyusan = 0
        if bilinen:
            en_cok = max(bilinen.count(x) for x in set(bilinen))
            uyusan = en_cok
        n_kaynak = len(bilinen)
        # D-6 (2026-08-06): 'local RISKLI plus MFEprimer TEMIZ' IS NOT A CONTRADICTION.
        # The two criteria are NESTED: the local scan allows up to 5 total mismatches and
        # DOES NOT APPLY the last-two-bases-at-the-3'-end condition (global_scan.py,
        # need_tail=False), while MFEprimer uses a thermodynamic criterion (a Tm cut at
        # 30 C). In other words the local criterion CONTAINS the MFEprimer criterion. The
        # looser criterion finding MORE hits than the stricter one is the EXPECTED result,
        # not a conflict between sources. A real contradiction is THE OTHER WAY ROUND: a hit
        # the strict one finds and the loose one misses.
        # These rows do not count as 'clean' either; they do not go to order, and the 3' end
        # is tested by hand.
        # D-15 (2026-08-07): the local layer's off-target count CANNOT be filtered
        # TAXONOMICALLY. The reason was measured: yerel_vuruslar.tsv keeps at most 20 hits
        # per set (raporla(), the '_vurus' list), so the identity of 4702 different-length
        # hits IS NOT on disk. For Bakteri_universal only 2 of the 100 sampled hits are of a
        # different length, so the sample does not answer this question.
        # For that reason the local layer's RISKLI vote cannot produce RISKLI on its own
        # while MFEprimer's CLADE FILTERED count is 0 -> INCELEME.
        _gevsek_fazla = (y_ok == 'RISKLI' and m_ok == 'TEMIZ'
                         and nb_ok != 'RISKLI')
        # D-17 (2026-08-07, MEASURED): do not let ORGANELLE products be hidden.
        # SILVA starts its chloroplast and mitochondrion records with "Bacteria;..."
        # (mitochondria under Rickettsiales, chloroplasts under Cyanobacteriota). The domain
        # test therefore counts them as INSIDE the target, which is taxonomically correct
        # for Bakteri_universal but wrong IN PRACTICE: these are plant organelles.
        # THE MEASUREMENT: in all 31 of Bakteri_universal's 31 organelle amplicons the F and
        # R mismatch count is ZERO and FpTm is 62.97 / RpTm 61.33, both ABOVE the Ta of
        # 57.9 C. So all 31 of the 31 products DO FORM under standard conditions (91-302 bp,
        # expected 130 bp). The hosts: Azolla, Isoetes, Equisetum, Ipomoea, Welwitschia,
        # Silene, and so on, that is, products that can occur in a digester fed with plant
        # material. This row cannot count as 'clean'.
        _organel = (kl or {}).get('ao') or 0
        _organel_notu = ''
        if _organel:
            _organel_notu = (u' | ORGANEL UYARISI: %d konak organel (kloroplast/'
                             u'mitokondri) amplikonu; olusabilir %s'
                             % (_organel, (kl or {}).get('olusabilir')))
        if _gevsek_fazla:
            # D-18 (2026-08-09): "the last two bases at the 3' end" DOES NOT DECIDE and is no
            # longer required. The reason was MEASURED and it has two layers:
            #  (a) MFEprimer 3.0 and later do not allow a mismatch at the 3' terminal base BY
            #      DEFINITION; the terminal base matching in 69/69 records is a forced output of
            #      the algorithm, not a finding in the data.
            #  (b) The bin scan was re-measured WITH the last-two condition REMOVED
            #      (2026-08-09, 17 targets): the dCq change was at most 0.41 cycles
            #      (Proteolitik_Synergistaceae -0.41; Petrimonas +0.36; on the remaining 15
            #      targets |difference| <= 0.09). No verdict changed.
            # What is asked for instead: the DISTANCE of the mismatch from the 3' end, and its
            # TYPE (see uc3_ceza_dongu). Detail: ESIK_VE_OLCUT_2026-08-08.md.
            karar = ('INCELEME - gevsek olcut vurusu (%s adet); 3\' uca yakin '
                     'uyumsuzlugun KONUMU ve TIPI degerlendirilmeli '
                     '(son iki baz sarti hukum vermez)' % y)
        elif len(set(bilinen)) > 1:
            karar = 'CELISKILI'
        elif n_kaynak == 0:
            karar = 'EKSIK - hicbir kaynak sonuc vermedi'
        elif set(bilinen) == {'TEMIZ'}:
            if _organel:
                # D-17: even when every layer looks clean, if there is an organelle
                # product that can form, the row is NOT clean; a human decision is needed.
                karar = ('INCELEME - katmanlar temiz ama %d organel amplikonu var'
                         % _organel)
            elif n_kaynak >= 3:
                karar = 'KESIN - uc olcum katmani da uyusuyor'
            elif n_kaynak == 2:
                karar = 'KOSULLU - iki katman uyusuyor, biri eksik'
            else:
                karar = 'EKSIK - yalnizca %d kaynak sonuc verdi' % n_kaynak
        else:
            karar = 'RISKLI - siparis edilmez'
        if _organel_notu and 'ORGANEL' not in karar:
            karar = karar + _organel_notu
        out.append(dict(c, kategori=karar_kategorisi(karar),
                        numune=n_ok, numune_deger=c['numune_deger'],
                        yerel=y_ok, yerel_urun=y, yerel_tum=y_tum,
                        yerel_ayni_boyda=yv.get('ayni_boyda'),
                        yerel_kume=(yerel.get(h, {}) or {}).get('kume', {}),
                        # A2: the taxonomic counters. THEY DO NOT ENTER THE VERDICT
                        # YET; the length based 'yerel_urun' still stands as the
                        # source of the verdict. The two are written side by side so
                        # that the difference is MEASURED first and decided after.
                        yerel_klad_ayrimi=(yv.get('siniflandirildi') or False),
                        yerel_a=(yv.get('sinif') or {}).get('a'),
                        yerel_ao=(yv.get('sinif') or {}).get('ao'),
                        yerel_b=(yv.get('sinif') or {}).get('b'),
                        yerel_c=(yv.get('sinif') or {}).get('c'),
                        yerel_bilinmiyor=(yv.get('sinif') or {}).get('bilinmiyor'),
                        yerel_klad_disi=((yv.get('sinif') or {}).get('b', 0)
                                         + (yv.get('sinif') or {}).get('c', 0))
                        if yv.get('siniflandirildi') else None,
                        mfeprimer=m_ok, mfe_urun=m_urun, mfe_ham=m_ham,
                        mfe_klad_ayrimi=m_klad_ayrimi,
                        mfe_a=(kl or {}).get('a'), mfe_ao=(kl or {}).get('ao'),
                        mfe_b=(kl or {}).get('b'), mfe_c=(kl or {}).get('c'),
                        mfe_olusabilir=(kl or {}).get('olusabilir'),
                        mfe_olusmaz=(kl or {}).get('olusmaz'),
                        mfe_durum=(mv or {}).get('durum', 'YAPILMADI'),
                        ncbi=nb_ok, ncbi_urun=nb_v,
                        ncbi_durum=(nb or {}).get('durum', 'YAPILMADI'),
                        kaynak_sayisi=n_kaynak, uyusan=uyusan, karar=karar))
    return out


# -------------------------------------------------------------------------
# It writes four files. CELISKILER.md is deliberately a separate file: the
# contradictions are the most valuable output of this round and must not get lost
# inside a long table.
#
# When there is no contradiction, the sentence written is NOT "everything is clean":
# layers that never ran count as MISSING and by definition produce no contradiction.
# That distinction is stated openly inside the report, because "no contradiction
# came up" is very easy to misread.
# -------------------------------------------------------------------------
def raporla(cikti, satirlar, yaz):
    t = os.path.join(cikti, 'dogrulama_uc_sutun.tsv')
    with open(t, 'w', encoding='utf-8', newline='') as fh:
        # 2026-08-21: English legend for human readers. Column NAMES are NOT
        # translated -- they are the machine-readable contract that other
        # stages and existing result files depend on.
        try:
            import sys as _s, os as _o
            _k = _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))
            if _k not in _s.path:
                _s.path.insert(0, _k)
            from screening import labels as _L
            fh.write(u'# FOUR INDEPENDENT MEASUREMENT LAYERS SIDE BY SIDE.\n')
            fh.write(u'# If the measuring layers disagree, the row is CELISKILI\n')
            fh.write(u'# (contradictory) and is NOT orderable.\n')
            fh.write(_L.verdict_legend(_L.OZGULLUK, 'VERDICTS'))
            fh.write(_L.verdict_legend(_L.KATMAN, 'PER-LAYER READINGS'))
            fh.write(_L.legend(_L.SUTUN.keys()))
        except Exception:
            pass
        fh.write(u'# THE MEASUREMENT LAYERS SIDE BY SIDE. If they disagree the row is CELISKILI (contradictory) and cannot be ordered.\n')
        # A3 (2026-08-21): the threshold can now be changed with PT_VURUS_ESIGI.
        # A criterion that can be changed makes its own measurement unreadable UNLESS the
        # value it ran with is written down, so it stands at the top of the output.
        fh.write(u'# VURUS_ESIGI = %d%s  (a layer that sees MORE off-target products than\n#   this votes RISKLI. It does not decide on its own; the verdict\n#   depends on the layers AGREEING)\n'
                 % (VURUS_ESIGI,
                    u'  [CHANGED via PT_VURUS_ESIGI, default is 0]'
                    if VURUS_ESIGI != 0 else u'  [default]'))
        fh.write(u'# 1_NUMUNE DOES NOT VOTE: it always reads TEMIZ (clean). It is the admission criterion, not a measurement (D-2).\n')
        fh.write(u'# 2_hedef_disi_urun: hits whose length DIFFERS from the expected product by more than +-%d bp.\n' % BOY_TOL)
        fh.write(u'# A2 (2026-08-21): the 2_klad_* columns are a TAXONOMIC separation and\n#   do NOT measure the same thing as 2_hedef_disi_urun. The size criterion\n#   is misleading: D-12 measured that 95.7% of the 1,605 amplicons called\n#   "off-target" were INSIDE the target clade (class a) and merely differed\n#   in length. 2_hedef_disi_urun STILL casts the vote; the taxonomic counts\n#   are printed beside it so the difference is MEASURED first. If the gap is\n#   large, which criterion should decide is a SEPARATE decision.\n#   2_klad_disi_b_c = b + c  (candidate real cross-reaction)\n#   2_bilinmiyor: the domain could not be resolved. NOT evidence, not a cross-reaction.\n#   If 2_klad_ayrimi_yapildi is HAYIR the class counts are ZERO, but that\n#   means "not measured", NOT "no cross-reaction".\n')
        fh.write(u'# 3_hedef_disi_amplikon: since D-12 this column is TAXONOMICALLY filtered:\n#   (b) same domain, outside the clade + (c) different domain. The RAW\n#   size-based count stays in 3_HAM_boya_dayali. The difference between the\n#   two is the target clade\'s own length variants (3_a), which are BY DESIGN\n#   for universal and group primers.\n# 3_olusabilir_Tm_yakin: of the (b)+(c)+(ao) records, those whose\n#   min(FpTm,RpTm) does NOT fall more than 5 C below the panel\'s Ta (%.1f C),\n#   that is, the ones that could form a product under standard conditions.\n' % TA_PANEL)
        fh.write(u'# 2_ayni_boyda_HEDEFIN_KENDISI: hits at the expected length. For group and\n#   universal primers these are BY DESIGN and are NOT counted as off-target (D-1).\n')
        fh.write(u'# 4_NCBI: if the status is "SONUC TAVANI" (result cap) or "BOS SONUC" (empty\n#   result) the value is NOT a count; that cell was NOT TESTED (D-3).\n')
        fh.write('# ' + OLCUT_NOTU.strip().replace('\n', '\n# ') + '\n')
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'cift_turu', 'F', 'R', 'urun_bp',
                    '1_NUMUNE_oy_vermez', '1_numune_deger',
                    '2_YEREL_DB', '2_hedef_disi_urun', '2_ayni_boyda_HEDEFIN_KENDISI',
                    '2_tum_vurus', '2_kume_dagilimi',
                    '2_klad_ayrimi_yapildi', '2_klad_disi_b_c',
                    '2_a_klad_ici', '2_ao_organel', '2_b_ayni_alan_klad_disi',
                    '2_c_farkli_alan', '2_bilinmiyor',
                    '3_MFEPRIMER', '3_hedef_disi_amplikon', '3_durum',
                    '3_HAM_boya_dayali', '3_klad_ayrimi_yapildi',
                    '3_a_klad_ici_uzunluk_varyanti', '3_ao_organel',
                    '3_b_ayni_alan_klad_disi', '3_c_farkli_alan',
                    '3_olusabilir_Tm_yakin', '3_olusmaz_Tm_dusuk',
                    '4_NCBI', '4_hedef_disi_urun', '4_durum',
                    'kaynak_sayisi', 'uyusan_kaynak', 'KARAR'])
        for s in satirlar:
            w.writerow([s['hedef'], s['tur'], s['F'], s['R'], s['urun'],
                        s['numune'], s['numune_deger'],
                        s['yerel'], s['yerel_urun'],
                        s.get('yerel_ayni_boyda', ''), s.get('yerel_tum', ''),
                        '; '.join('%s=%s' % kv for kv in (s['yerel_kume'] or {}).items()),
                        ('EVET' if s.get('yerel_klad_ayrimi') else 'HAYIR - yalniz boy'),
                        s.get('yerel_klad_disi', ''),
                        s.get('yerel_a', ''), s.get('yerel_ao', ''),
                        s.get('yerel_b', ''), s.get('yerel_c', ''),
                        s.get('yerel_bilinmiyor', ''),
                        s.get('mfeprimer', '-'), s.get('mfe_urun', ''), s.get('mfe_durum', ''),
                        s.get('mfe_ham', ''),
                        ('EVET' if s.get('mfe_klad_ayrimi') else 'HAYIR - ham sayi'),
                        s.get('mfe_a', ''), s.get('mfe_ao', ''),
                        s.get('mfe_b', ''), s.get('mfe_c', ''),
                        s.get('mfe_olusabilir', ''), s.get('mfe_olusmaz', ''),
                        s['ncbi'], s['ncbi_urun'], s['ncbi_durum'],
                        s.get('kaynak_sayisi', ''), s.get('uyusan', ''), s['karar']])
    yaz(u'  written: %s' % t)

    celiskili = [s for s in satirlar if s['karar'] == 'CELISKILI']
    c = os.path.join(cikti, 'CELISKILER.md')
    with open(c, 'w', encoding='utf-8') as fh:
        fh.write(u'# Contradictions, the most valuable output of this round\n\n')
        if not celiskili:
            fh.write(u'In this run the layers **did not disagree on a single row**.\n\nThat does NOT mean "everything is clean": the NCBI or the local layer')
        for s in celiskili:
            fh.write(u'## %s\n\n' % s['hedef'])
            fh.write(u'| source | result | value |\n|---|---|---|\n\n')
            fh.write(u'| 1 in-sample measurement | (DOES NOT VOTE - fixed value) | %s |\n\n' % s['numune_deger'])
            fh.write(u'| 2 local databases (ours) | %s | %s off-target (+ %s at the target\'s OWN length, %s in total) |\n\n'
                     % (s['yerel'], s['yerel_urun'], s.get('yerel_ayni_boyda', '-'),
                        s.get('yerel_tum', '-')))
            fh.write(u'| 3 MFEprimer (BAGIMSIZ) | %s | %s amplikon |\n' % (s.get('mfeprimer', '-'), s.get('mfe_urun', '-')))
            fh.write(u'| 4 NCBI (BAGIMSIZ) | %s | %s |\n\n' % (s['ncbi'], s['ncbi_urun']))
            fh.write(u'**What it means:** a sample is 99 bins, not the world. A pair that looks clean in the sample but hits in a database')
            fh.write(u'**What to do:** this row CANNOT BE ORDERED. Which organism the hit is in is shown in `dogrulama_uc_sutun.tsv` and')
    yaz(u'  written: %s' % c)

    v = os.path.join(cikti, 'yerel_vuruslar.tsv')
    with open(v, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'kume', 'kayit_basligi', 'urun_bp', 'F_uyumsuzluk', 'R_uyumsuzluk'])
        for s in satirlar:
            for vv in (s.get('_vurus') or []):
                w.writerow([s['hedef']] + list(vv))
    yaz(u'  written: %s' % v)

    r = os.path.join(cikti, 'DOGRULAMA_RAPORU.md')
    # KATEGORI bazinda say (sayidan arindirilmis anahtar - bkz. C-1 notu).
    say = {}
    for s in satirlar:
        kg = s.get('kategori') or karar_kategorisi(s.get('karar'))
        say[kg] = say.get(kg, 0) + 1
    # ayrinti ayri tutulur: kategori -> {tam hukum dizgesi: adet}
    ayrinti = {}
    for s in satirlar:
        kg = s.get('kategori') or karar_kategorisi(s.get('karar'))
        ayrinti.setdefault(kg, {})
        ayrinti[kg][s['karar']] = ayrinti[kg].get(s['karar'], 0) + 1
    with open(r, 'w', encoding='utf-8') as fh:
        fh.write(u'# Verification round\n\nGenerated: %s, script %s\n\n'
                 % (time.strftime('%Y-%m-%d %H:%M'), VERSIYON))
        if VURUS_ESIGI != 0:
            fh.write(u'> **WARNING: the hit threshold is not the default.** This run used `PT_VURUS_ESIGI=%d` (the default is 0). Before ordering' % VURUS_ESIGI)
        fh.write(u'## Result\n\n')
        for k in KATEGORILER:
            if k in say or k in ANA_KATEGORILER:
                fh.write(u'- **%s**: %d\n' % (k, say.get(k, 0)))
        for k in sorted(x for x in say if x not in KATEGORILER):
            fh.write(u'- **%s**: %d\n' % (k, say[k]))
        fh.write(u'- _toplam_: %d\n' % sum(say.values()))
        fh.write(u'\n### Detail (reasons within each category)\n\n')
        for k in list(KATEGORILER) + sorted(x for x in ayrinti if x not in KATEGORILER):
            if k not in ayrinti:
                continue
            for gerekce, n in sorted(ayrinti[k].items(), key=lambda t: -t[1]):
                fh.write(u'- %s (%d): %s\n' % (k, n, gerekce))
        fh.write(u'\n> Read the contradictions first: `CELISKILER.md`\n\n')
        fh.write(u'```' + OLCUT_NOTU + u'```\n\n## The layers side by side\n\n')
        fh.write(u'| target | 1 in-sample | 2 local DB | 3 MFEprimer | 4 NCBI | agreeing | verdict |\n|---|---|---|---|---|---|---|\n')
        for s in satirlar:
            fh.write(u'| %s | %s | %s (%s) | %s (%s) | %s (%s) | %s/%s | **%s** |\n'
                     % (s['hedef'], s['numune'], s['yerel'], s['yerel_urun'],
                        s.get('mfeprimer', '-'), s.get('mfe_urun', '-'),
                        s['ncbi'], s['ncbi_urun'],
                        s.get('uyusan', '-'), s.get('kaynak_sayisi', '-'), s['karar']))
        fh.write(u'\n## How to complete the NCBI layer\n\nIf the automatic route did not run, or failed:\n\n1. Paste the rows from `NCBI_PRIMER_BLAST_GIRDI.tsv` into https://www.ncbi.nlm.nih.gov/tools/primer-blast/\n2. Write the results into `NCBI_SONUC_SABLONU.tsv`.\n3. Run `verification/full_chain.py` -> (D) -> "load manual results".\n')
    yaz(u'  written: %s' % r)
    yaz('')
    # SUMMARY: the four main categories on one line, with keys stripped OF THE NUMBER.
    sirali = [k for k in KATEGORILER if k in say or k in ANA_KATEGORILER] + \
             sorted(x for x in say if x not in KATEGORILER)
    yaz(u'  SUMMARY   ' + '   '.join('%s: %d' % (k, say.get(k, 0)) for k in sirali)
        + u'   |   total: %d' % sum(say.values()))
    yaz(u'  DETAIL:')
    for k in sirali:
        for gerekce, n in sorted(ayrinti.get(k, {}).items(), key=lambda t: -t[1]):
            yaz(u'    %-10s %2d  %s' % (k, n, gerekce))



# --------------------------------------------------------------- guvenlik agi
def cikti_denetle(yaz, ad, dosyalar, asgari=1):
    """When the stage ends, it audits ITS OWN output.

        If the expected row count is zero, or the file is missing entirely, it DOES
        NOT CARRY ON SILENTLY: it prints a clear error and returns a non-zero code.
        This is so that it cannot produce an empty result overnight and then say
        "nothing was found" in the morning.

    """
    sorun = []
    for yol, etiket in dosyalar:
        if not os.path.exists(yol):
            sorun.append(u'%s URETILMEDI (%s)' % (etiket, yol)); continue
        try:
            with open(yol, encoding='utf-8') as fh:
                n = sum(1 for s in fh if s.strip() and not s.startswith('#'))
            n = max(0, n - 1)          # baslik satiri
        except OSError as e:
            sorun.append(u'%s OKUNAMADI (%s)' % (etiket, e)); continue
        if n < asgari:
            sorun.append(u'%s BOS - %d veri satiri (en az %d bekleniyordu)'
                         % (etiket, n, asgari))
    if not sorun:
        return 0
    yaz('')
    yaz('  ' + '!' * 70)
    yaz(u'  STAGE %s PRODUCED EMPTY OUTPUT - THE CHAIN WAS STOPPED HERE' % ad)
    for x in sorun:
        yaz(u'    - %s' % x)
    yaz('')
    yaz(u'  WHY IT STOPPED: the next stage would have read this file as input.')
    yaz(u'  Continuing with empty input does not crash; it produces a MEANINGLESS BUT')
    yaz(u'  CONVINCING summary, which is exactly the silent failure we hunt for.')
    yaz(u'  Read the run log above, fix the cause, then run the same command')
    yaz(u'  again; finished work is skipped from its checkpoints.')
    yaz('  ' + '!' * 70)
    return 4


def girdi_denetle(yaz, ad, dosyalar):
    """Before the stage STARTS: do the files it needs exist, and are they non-empty?"""
    eksik = []
    for yol, etiket, uretici in dosyalar:
        if not os.path.exists(yol):
            eksik.append(u'%s yok (%s) -> once %s asamasini kosun' % (etiket, yol, uretici))
            continue
        with open(yol, encoding='utf-8') as fh:
            n = sum(1 for s in fh if s.strip() and not s.startswith('#'))
        if n <= 1:
            eksik.append(u'%s BOS (%s) -> %s asamasi sonuc uretmemis'
                         % (etiket, yol, uretici))
    if not eksik:
        return 0
    yaz('')
    yaz('  ' + '!' * 70)
    yaz(u'  STAGE %s WAS NOT STARTED - INPUT MISSING' % ad)
    for x in eksik:
        yaz(u'    - %s' % x)
    yaz('  ' + '!' * 70)
    return 5

# -------------------------------------------------------------------------
# In this script the driver sits directly inside main(). The layers run in order:
#   2) the local database  ->  3) MFEprimer  ->  4) NCBI  ->  combine and report
# (Layer 1 is the sample measurement and arrives ready from round K; it is not run
# here.)
#
# EXIT CODES: 7 = A CASCADED FAILURE. Stage K produced no rows at all, so the input
# is empty; this IS NOT D's fault and it says so. 5 = its own input is missing,
# 4 = its own output is empty. The distinction matters: debugging the wrong stage
# costs hours.
# -------------------------------------------------------------------------

# --- CLI value normalisation ------------------------------------------------
# English option values are accepted alongside the original Turkish ones and
# mapped back here. The internal values are unchanged on purpose: they are
# compared in dozens of places and, in some cases, written to output files.
# Translating the interface must not translate the data.
_DEGER = {'auto': 'oto', 'manual': 'elle', 'none': 'yok', 'quick': 'hizli',
          'full': 'tam', 'member': 'uye', 'competitor': 'rakip',
          'exclude': 'disla'}


def _ing_deger(a):
    for _ad in ('nt', 'literatur', 'ncbi', 'karisik', 'moduller', 'mod'):
        _v = getattr(a, _ad, None)
        if isinstance(_v, str) and _v in _DEGER:
            setattr(a, _ad, _DEGER[_v])
    return a

def main():
    p = argparse.ArgumentParser(description='Kurtarilan ciftlerin uc katmanli dogrulanmasi')
    p.add_argument('--root', '--kok', dest='kok', default='.')
    p.add_argument('--ncbi', choices=['auto', 'manual', 'none', 'oto', 'elle', 'yok'], default='elle',
                   help='auto: NCBI URL API; manual: write pasteable input; none: skip')
    p.add_argument('--ncbi-load', '--ncbi-yukle', dest='ncbi_yukle', default=None, help='doldurulmus NCBI_SONUC_SABLONU.tsv')
    p.add_argument('--organism', '--organizma', dest='organizma', default='', help='NCBI organizma kisiti (bos = tum nt)')
    # D-13c (2026-08-07, MEASURED): the target's OWN taxon can be excluded with
    # ENTREZ_QUERY, and then every product left on the page is off-target by definition.
    # CAUTION, a measured trap: a bare 'NOT txidN[Organism]' INVERTS the filter (it
    # returns only that taxon). The correct form is 'all[filter] NOT txidN[Organism]'.
    # The code adds that prefix itself.
    p.add_argument('--ncbi-exclude-taxid', '--ncbi-haric-taxid', dest='ncbi_haric_taxid', default='',
                   help="taxid to EXCLUDE at NCBI (example: 2157). Added to ENTREZ_QUERY "
                        "'all[filter] NOT txid<N>[Organism]' olarak gonderilir.")
    p.add_argument('--local-only', '--yalniz-yerel', dest='yalniz_yerel', action='store_true', help='only katman 2 (yerel DB)')
    p.add_argument('--no-mfe', '--mfe-yok', dest='mfe_yok', action='store_true',
                   help='MFEprimer katmanini skip')
    p.add_argument('--cluster-max', '--kume-ust', dest='kume_ust', type=int, default=0,
                   help='only en kucuk N veritabani (hizli test)')
    p.add_argument('--parc-set', '--parc', dest='parc', action='store_true',
                   help='also scan the SILVA LSU Parc set (slow; not needed for specificity)')
    p.add_argument('--order', '--siparis', dest='siparis', action='store_true',
                   help='kurtarilanlar yerine SIPARIS LISTESINDEKI ciftleri dogrula '
                        '(siparis oncesi Primer-BLAST kontrolu icin)')
    p.add_argument('--order-all', '--siparis-hepsi', dest='siparis_hepsi', action='store_true',
                   help='with --order: include CONDITIONAL and NOT-RECOMMENDED rows as well')
    # -----------------------------------------------------------------------
    # --tumu  (2026-08-07)
    # A user request: "search all of them in that database as well". --siparis tests
    # only KESIN and EVRENSEL (16 pairs); the KOSULLU and ONERILMEZ rows (6 pairs) were
    # NEVER put through SILVA. --tumu adds those too -> 22 pairs, against every indexed
    # database, SILVA included.
    # The --siparis mode IS UNCHANGED; --tumu is a separate flag built on top of it.
    # -----------------------------------------------------------------------
    p.add_argument('--all', '--tumu', dest='tumu', action='store_true',
                   help='EVERY pair in the panel (CERTAIN + UNIVERSAL + CONDITIONAL +'
                        'ONERILMEZ) butun indeksli veritabanlarina, SILVA dahil')
    p.add_argument('--ncbi-order-only', '--ncbi-yalniz-siparis', dest='ncbi_yalniz_siparis', action='store_true',
                   help='--tumu with: KATMAN 4 (NCBI) only siparis listesindeki '
                        '(KESIN/EVRENSEL) ciftlere kosar. Listede olmayanlar '
                        'yalniz yerel + MFEprimer katmanlarini gorur. NCBI cift '
                        'basina ~75 sn + 10 sn bekleme oldugu icin sure kalemi.')
    p.add_argument('--reset', '--sifirla', dest='sifirla', action='store_true')
    a = p.parse_args()
    a = _ing_deger(a)

    # --tumu = --siparis plus --siparis-hepsi. NOT a separate route, the existing route
    # with a wider input set. That way the behaviour of --siparis mode does not change.
    if getattr(a, 'tumu', False):
        a.siparis = True
        a.siparis_hepsi = True
    if getattr(a, 'ncbi_yalniz_siparis', False) and not getattr(a, 'siparis_hepsi', False):
        # There is nothing to restrict: the input set is already the order list alone.
        a.ncbi_yalniz_siparis = False

    kok = os.path.abspath(a.kok)
    if not os.path.isdir(os.path.join(kok, 'screening')):
        sys.exit('HATA: %s icinde screening yok.' % kok)
    CIKTI = os.path.join(kok, 'DOGRULAMA_SONUC')
    KONTROL = os.path.join(CIKTI, 'kontrol')
    os.makedirs(KONTROL, exist_ok=True)
    if a.sifirla:
        for f in os.listdir(KONTROL):
            try:
                os.remove(os.path.join(KONTROL, f))
            except OSError as e:
                print('  silinemedi: %s (%s)' % (f, e))
    g = open(os.path.join(CIKTI, 'kosu_gunlugu.txt'), 'a', encoding='utf-8')

    def yaz(s=''):
        print(s, flush=True); g.write(s + '\n'); g.flush()

    yaz('=' * 78)
    yaz(u'  SPECIFICITY ROUND - recovered pairs tested against independent evidence layers')
    yaz(u'  version %s   %s' % (VERSIYON, time.strftime('%Y-%m-%d %H:%M')))
    yaz('=' * 78)

    kyol_ = os.path.join(kok, 'KURTARMA_SONUC', 'kurtarma_satirlari.tsv')
    if not os.path.exists(kyol_) or sum(
            1 for x in open(kyol_, encoding='utf-8') if x.strip() and not x.startswith('#')) <= 1:
        yaz('')
        yaz('  ' + '!' * 70)
        yaz(u'  STAGE D WAS NOT RUN - THE INPUT IS EMPTY')
        yaz(u'  Cause: stage K (recovery) produced no rows at all, so')
        yaz(u'  NO pairs to verify. This is NOT a failure of D; it is a knock-on from K,')
        yaz(u'  falling is a KNOCK-ON effect of that.')
        yaz(u'  What to do: first find out why stage K produced no rows.')
        yaz(u'  D is sound on its own; it was tested with a hand prepared input.')
        yaz('  ' + '!' * 70)
        g.close()
        return 7           # 7 = A CASCADED FAILURE (empty input), 5 = its own input is missing
    if getattr(a, 'siparis', False):
        ciftler, kyol = siparistekiler(kok, hepsi=getattr(a, 'siparis_hepsi', False))
    else:
        ciftler, kyol = kurtarilanlar(kok)
    if not ciftler:
        yaz(u'  The recovery round produced no NEW or CHANGED pair above the threshold, so there is nothing to verify.')
        return 0
    yaz(u'  kaynak            : %s' % os.path.basename(kyol))
    yaz(u'  source path       : %s' % kyol)
    yaz(u'  pairs to verify   : %d' % len(ciftler))
    if getattr(a, 'tumu', False):
        _sip = sum(1 for c in ciftler if c.get('sipariste'))
        yaz(u'  MODE: --all  (EVERY pair in the panel; %d on the order list, %d outside it)' % (_sip, len(ciftler) - _sip))
    if _ATLANAN:
        yaz(u'  SKIPPED (primer sequence not found; absent from the protocol output): %s'
            % ', '.join(_ATLANAN))
    for c in ciftler:
        yaz(u'     - %-42s %s' % (c['hedef'][:42], c['tur']))
    n_kume = sum(1 for _, d, _ in KUMELER if os.path.exists(os.path.join(kok, 'REFERANS_DB', d)))
    yaz('')
    # THE TIME STATEMENT (2026-08-07): the old line ESTIMATED 240 s per set and printed
    # "11.7 hours". MEASUREMENT did not bear that out: with the checkpoints ready, the
    # local layer took 26 s for 11 sets x 16 pairs (the 2026-08-07 measurement).
    # Instead of an estimate, THE NUMBER OF READY CHECKPOINTS is reported, so the reader
    # can see what will actually be re-run.
    import hashlib as _h
    _imza = _h.md5('|'.join(sorted('%s>%s<%s' % (a['hedef'], a.get('F', ''),
                   a.get('R', '')) for a in ciftler)).encode()).hexdigest()[:10]
    _hazir = sum(1 for e, d, _a in KUMELER
                 if os.path.exists(os.path.join(KONTROL, 'yerel_%s_%s.pkl'
                                                % (re.sub(r'\W+', '_', e), _imza))))
    yaz(u'  LAYER 2 (local): %d sets, %d from checkpoints, %d to be scanned from scratch.' % (n_kume, _hazir, max(0, n_kume - _hazir)))
    yaz(u'     (measured with everything cached: ~30 s; a fresh scan takes minutes per set)')
    yaz(u'  LAYER 3 (MFEprimer): 6 indexes including SILVA. Measured for SILVA: 16 pairs took ~85 s of spec time plus ~40 s copying evidence.')
    # KATMAN 4 kapsami: --ncbi-yalniz-siparis verilmisse NCBI yalniz siparis
    # listesindeki ciftlere kosar. Sure beyani da o sayidan hesaplanir.
    _ncbi_ciftler = ([c for c in ciftler if c.get('sipariste')]
                     if getattr(a, 'ncbi_yalniz_siparis', False) else list(ciftler))
    _n4 = len(_ncbi_ciftler)
    yaz(u'  LAYER 4 (NCBI): MEASURED ~75 s per pair (submit plus poll) plus %d s between submissions, so ~%s for %d pairs.'
        % (PB_GONDERIM_ARASI, _n4,
           sure_metni(_n4 * 75 + max(0, _n4 - 1) * PB_GONDERIM_ARASI)))
    if getattr(a, 'ncbi_yalniz_siparis', False):
        _dis = [c['hedef'] for c in ciftler if not c.get('sipariste')]
        yaz(u'  LAYER 4 RESTRICTED (--ncbi-order-only): %d pairs will NOT be seen by NCBI.'
            % len(_dis))
        yaz(u'     The NCBI column stays BILINMIYOR (unknown) on those rows. That is NOT "clean".')
        for _d in _dis:
            yaz(u'       - %s' % _d)
    yaz(u'  Resumable: state is saved after every set.')
    yaz('')

    yaz(u'--- LAYER 2: LOCAL DATABASE SCAN ---')
    yerel = katman1_yerel(kok, ciftler, yaz, KONTROL, a.parc, a.kume_ust)

    # --- KATMAN 3: MFEprimer (BAGIMSIZ ARAC) ---
    mfe_sonuc = {}
    if not a.mfe_yok:
        yaz(u'--- LAYER 3: MFEprimer (INDEPENDENT TOOL) ---')
        yaz(u'  The first two layers are also OUR code and share the same engine; if that')
        yaz(u'  engine has a bug, both will be wrong in the same direction. This layer puts')
        yaz(u'  an independent, external tool to the same questions.')
        import importlib.util as _u
        _sp = _u.spec_from_file_location(
            'mfe_katmani', os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        'mfeprimer_layer.py'))
        MK = _u.module_from_spec(_sp); _sp.loader.exec_module(MK)
        mfe = MK.mfe_bul(kok)
        if not mfe:
            yaz(u'  MFEprimer binary not found (tools/mfeprimer); the layer was SKIPPED.')
        else:
            yaz(u'  ikili: %s' % mfe)
            dur, spec = MK.spec_kos(kok, mfe, ciftler, CIKTI, yaz, KONTROL)
            yapi = MK.yapi_kos(kok, mfe, ciftler, CIKTI, yaz)
            # D-7: an IDENTITY file beside the 'off-target N' count. The number on its own
            # is read wrongly (see mfe_katmani.hedef_disi_kimlikleri).
            try:
                MK.hedef_disi_kimlikleri(CIKTI, ciftler, yaz)
            except Exception as _e:
                yaz(u'  WARNING: could not write the off-target identity file: %s' % _e)
            if dur.get('durum') == 'TAMAM':
                for c in ciftler:
                    mfe_sonuc[c['hedef']] = dict(
                        durum='TAMAM', hedef_disi=MK.hedef_disi_say(spec, c['hedef'], ciftler),
                        yapi=yapi, kullanilan_db=dur.get('kullanilan', []))
            else:
                yaz(u'  The MFEprimer layer returned no result: %s' % dur.get('sebep', ''))

    ncbi = {}
    if a.ncbi_yukle:
        yaz(u'--- LAYER 4: NCBI (loading manually entered results) ---')
        ncbi = ncbi_yukle(a.ncbi_yukle, yaz)
        yaz(u'  %d rows loaded' % len(ncbi))
    elif a.yalniz_yerel or a.ncbi == 'yok':
        yaz(u'--- LAYER 4: NCBI SKIPPED (by request) ---')
    elif a.ncbi == 'oto':
        yaz(u'--- KATMAN 4: NCBI OTOMATIK (URL API) ---')
        yaz(u'  Not: blastn -remote KULLANILMIYOR (45 sn tavanini asiyor).')
        if _n4 != len(ciftler):
            yaz(u'  COVERAGE: %d/%d pairs (those on the order list). The rest were only seen by layers 2 and 3.' % (_n4, len(ciftler)))
        ncbi = katman2_oto(_ncbi_ciftler, CIKTI, yaz, a.organizma,
                          haric_taxid=getattr(a, 'ncbi_haric_taxid', '') or '')
        if not any(v.get('durum', '').startswith('TAMAM') for v in ncbi.values()):
            yaz(u'  The automatic route returned nothing; writing the manual route files.')
            katman2_elle_girdi(_ncbi_ciftler, CIKTI, yaz, a.organizma)
    else:
        yaz(u'--- LAYER 4: NCBI MANUAL (writing the input and the template) ---')
        katman2_elle_girdi(_ncbi_ciftler, CIKTI, yaz, a.organizma)

    # D-12: what should enter the verdict is not MFEprimer's raw (length based) count
    # but the TAXONOMICALLY filtered klad_disi. The filter reads
    # mfe_hedef_disi_kimlikleri.tsv and screening/hedef_klad.tsv; if either is missing
    # it returns EMPTY, the raw count is used, and the report says so openly.
    klad_sonuc = {}
    try:
        import verification.mfeprimer_layer as _MK2
    except ImportError:
        _MK2 = None
        try:
            sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
            import mfeprimer_layer as _MK2
        except ImportError:
            _MK2 = None
    if _MK2 is not None and hasattr(_MK2, 'klad_siniflandir'):
        try:
            klad_sonuc = _MK2.klad_siniflandir(kok, CIKTI, ciftler, TA_PANEL, yaz)
        except Exception as _e:
            yaz(u'  klad suzgeci koselemedi (%s: %s) - ham sayi kullanilacak'
                % (type(_e).__name__, _e))
            klad_sonuc = {}
    satirlar = birlestir(ciftler, yerel, ncbi, mfe_sonuc, klad_sonuc)
    for s in satirlar:
        s['_vurus'] = (yerel.get(s['hedef'], {}) or {}).get('vurus', [])
    raporla(CIKTI, satirlar, yaz)
    rc = cikti_denetle(yaz, 'D (DOGRULAMA)', [
        (os.path.join(CIKTI, 'dogrulama_uc_sutun.tsv'), 'dogrulama_uc_sutun.tsv')])
    g.close()
    return rc


if __name__ == '__main__':
    sys.exit(main() or 0)
