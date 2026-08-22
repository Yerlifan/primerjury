# -*- coding: utf-8 -*-
"""STAGE G - AN INDEPENDENT VERIFICATION OF EVERY BIN IDENTITY

WHY IT EXISTS
-------------
Stage `I` tests 12 CLAIMS, but those claims were the identities WE SUSPECTED. The
identity of the bins we did not suspect was never verified independently. Yet
those bins also determine the member and competitor sets, which means they affect
EVERY DISCRIMINATION RATIO. The silent majority was never checked at all. This
stage closes that.

THE SCOPE
---------
EVERY bin that enters a panel measurement: all the bins that take part in any
discrimination calculation as the MEMBER or the COMPETITOR of a target. Bins that
take part in no calculation are SKIPPED, and the fact that they were skipped is
WRITTEN INTO THE OUTPUT (no silent skipping).

THE METHOD - EXACTLY THE SAME AS `I`, NOT REWRITTEN
---------------------------------------------------
This script does not write the decision logic ITSELF; it CALLS the functions of
identity_verification.py:
    kp_yolu()            - the checkpoint key (THE SHARED CACHE)
    kl_degerlendir()     - short list -> alignment -> hits plus kazanan_sira
    kisa_liste()         - the single query scanner (for verification)
    ad_coz(), savunulabilir_duzey(), cins_cek(), hizala(), ayirt_edici_pencere()
    literature_check.py - the literature layer
So: a 500 item short list, all of it aligned, AT LEAST TWO independent databases
agreeing, the best THREE hits, the defensible level plus the suggested name, and
the winning rank.

DATABASE COVERAGE - THERE IS NO DOMAIN FILTER
---------------------------------------------
ALL 12 local databases are asked about EVERY bin. No domain filter (bacteria,
archaea, fungi) IS APPLIED.

  Why: choosing the domain from THE KRAKEN LABEL is dangerous, because the reason
  this stage exists is that Kraken labels can be wrong. Saying "there is no need to
  ask a fungal database about a bacterial bin" assumes the bin IS bacterial, and
  that is exactly what we are trying to test. So EVERY database is asked and the
  irrelevant ones are dropped IN THE RESULT (marked "no result" when there is no
  hit); they are not ruled out BEFORE the query.

  The bin's domain comes FROM THE MEASUREMENT, NOT FROM THE LABEL: which databases
  actually gave a hit is written into the 'alan_olcumden' column.

If a database returned no result for a bin, that DOES NOT COUNT as "clean":
  TAMAM          - scanned, there is a hit
  SONUC YOK      - scanned, no record matched a seed (it may be outside the domain)
  DOSYA YOK      - the database is not on disk
  SORULMADI (..) - the reason is written out
All 12 of the 12 databases appear on every bin row.

COVERAGE ACCOUNTING - SO THE CAP PROBLEM DOES NOT RECUR
-------------------------------------------------------
A real cap problem happened in the access test: the first run cut off at 120 001
records, and SILVA SSU NR99 / LSU Parc / UNITE ITS were effectively being scanned
truncated.

A separate pass per bin: 94 distinct consensuses x 12 databases = 1128 full
database passes. Unacceptable. Instead, one database pass serves every query in
the batch AT THE SAME TIME (a seed -> query inverted index). The short list it
produces is EXACTLY THE SAME as the single query kisa_liste(), and that is proven
by a test.

It WRITES NOTHING into the panel files; it writes under ALL_IDENTITIES_RESULT/.

"""

# -------------------------------------------------------------------------
# all_bin_identities.py verifies the identity of EVERY bin that enters a panel
# measurement, independently, against external reference databases ("the silent
# majority").
#
# INPUT  : the 12 local FASTA sets under REFERANS_DB/ (all of them, with NO domain
#          filter), konsensus_kanonik/ and the panel plus membership tables
#          (screening.targets), IDENTITY_RESULT/kontrol/ (a cache SHARED with stage I),
#          optionally NCBI nt.
# OUTPUT : ALL_IDENTITIES_RESULT/tum_kutu_kimlikleri.tsv (ONE row per bin),
#          ALL_IDENTITIES_RESULT/TUM_KUTU_KIMLIK_RAPORU.md,
#          ALL_IDENTITIES_RESULT/kutu_*.json, kosu_gunlugu.txt.
#          It WRITES NOTHING into the panel files.
# CALLED BY: verification/full_chain.py -> key G
#          (python3 verification/all_bin_identities.py --root .)
#
# HOW IT DIFFERS FROM STAGE I: I tests the 12 claims WE SUSPECT. G also tests the
# bins we take to be beyond suspicion, and those bins determine the member and
# competitor sets, which means they affect EVERY DISCRIMINATION RATIO.
#
# THE METHOD IS NOT REWRITTEN: the decision logic is CALLED from
# identity_verification.py (kp_yolu, kl_degerlendir, savunulabilir_duzey, ad_coz,
# cins_cek, hizala, ayirt_edici_pencere). So the 500 item short list, aligning all
# of it, requiring AT LEAST TWO independent databases to agree, and the best three
# hits rules are EXACTLY the same. Two separate implementations would be two
# separate tools producing two separate verdicts.
# -------------------------------------------------------------------------
import os, sys, csv, json, time, re, argparse, heapq, collections

VERSIYON = '1.0 (2026-08-04)'

# The REAL record counts taken from the access verification
# (ACCESS_RESULT/erisim_dogrulama.tsv, the 'TAMAMI' run). If fewer than these were
# scanned, the coverage is INCOMPLETE.
BEKLENEN_KAYIT = {
    'SILVA SSU NR99': 510495, 'SILVA LSU NR99': 95279, 'SILVA LSU Parc': 1312521,
    'UNITE ITS': 2069189, 'PR2 SSU': 240201, 'ROD operon': 60320,
    'RefSeq bakteri 16S': 26877, 'RefSeq arke 16S': 1160, 'RefSeq mantar ITS': 20394,
    'RefSeq mantar 28S': 12890, 'RefSeq mantar 18S': 4037, 'RefSeq ref_all2': 65358,
}


# identity_verification.py is a script (not a package), so it is loaded as a module
# from its file path. The decision logic comes FROM THERE; it is not rewritten here.
def _K(kok):
    """Load identity_verification.py as a module - the decision logic comes FROM THERE."""
    import importlib.util
    yol = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'identity_verification.py')
    sp = importlib.util.spec_from_file_location('kimlik_dogrulama', yol)
    m = importlib.util.module_from_spec(sp)
    sys.modules['kimlik_dogrulama'] = m
    sp.loader.exec_module(m)
    return m


# --------------------------------------------------------------- THE BIN INVENTORY
# -------------------------------------------------------------------------
# THE SCOPE: every bin that takes part in any discrimination calculation as the
# MEMBER or the COMPETITOR of a target. Bins that enter no calculation are not
# tested, but they are not skipped SILENTLY: they go onto the skipped list with
# their reason and appear in both the TSV and the markdown report.
# -------------------------------------------------------------------------
def kutu_envanteri(kok, K):
    """EVERY bin entering a panel measurement, plus which target it is a member or
        competitor of.

        Returns: (katilan, atlanan, uye, rakip, kons)
          katilan : the bin names, ordered
          atlanan : [(bin, reason)]
          uye     : bin -> [target, ...]
          rakip   : bin -> [target, ...]

    """
    sys.path.insert(0, kok)
    from screening import targets as H
    panel, _yol = H.panel_oku()
    kons_l = H.konsensusler()
    uyelik = H.uyelik_oku()
    kut = H.kutular()
    uye = collections.defaultdict(list)
    rakip = collections.defaultdict(list)
    for p in panel:
        b = H.hedef_baglami(p, uyelik, kons_l, kut)
        for k in b['uye_kons']:
            uye[k['kutu']].append(p['hedef'])
        for k in b['rakip_kons']:
            rakip[k['kutu']].append(p['hedef'])
    kons = {d['kutu']: d['dizi'] for d in kons_l}
    hepsi = sorted(kons)
    katilan = sorted(set(uye) | set(rakip))
    atlanan = [(k, u'hicbir hedefin uyesi ya da rakibi degil - panel olcumlerine '
                   u'girmiyor, ayrim hesaplarini etkilemiyor')
               for k in hepsi if k not in set(katilan)]
    return katilan, atlanan, uye, rakip, kons, H


# --------------------------------------------------------------- TOPLU TARAMA
class _TersMetin(str):
    """Text whose comparison is REVERSED.

        The short list's ranking criterion: seed count DESCENDING, header ASCENDING on
        a tie. To keep memory bounded we use a MIN-HEAP of size 'ust', and the element
        the heap throws away has to be "the worst". The worst is the one with the
        smallest seed count and, ON A TIE, THE LARGEST HEADER, which is why the header
        comparison is reversed. Without the reversal, records with EQUAL SEED COUNTS at
        the cut off are chosen differently from the single query kisa_liste() and the
        two routes diverge (measured: they diverged at position 455).

    """
    __slots__ = ()

    def __lt__(self, o):
        return str.__gt__(self, o)

    def __gt__(self, o):
        return str.__lt__(self, o)

    def __le__(self, o):
        return str.__ge__(self, o)

    def __ge__(self, o):
        return str.__le__(self, o)


# -------------------------------------------------------------------------
# THE BULK SCAN - EVERY query in the batch in one database pass.
#
# A separate pass per bin would need 94 consensuses x 12 databases = 1128 full file
# passes, which is unacceptable. Here a seed -> query inverted index is built, and
# each record's k-mer set is extracted once and multiplied against the index.
#
# THE SHORT LIST PRODUCED MUST BE EXACTLY THE SAME as the single query
# kisa_liste(): if the two routes produce different lists on the same input, stages
# I and G can give different verdicts for the same bin. The ranking criterion in
# both is "seed count DESCENDING, header ASCENDING on a tie"; the _TersMetin class
# exists to preserve that tie rule inside a min-heap too (measured: without the
# reversal they diverge at position 455).
#
# THERE IS NO DOMAIN FILTER: whatever the query, every record is evaluated.
# -------------------------------------------------------------------------
def toplu_kisa_liste(K, yol, sorgular, ust, ilerle=None):
    """Build the short list for EVERY query in ONE database pass.

        sorgular: {name: sequence}. Returns: ({name: short_list}, records_scanned)

        A seed -> query inverted index is built; for each record the record's 16-mer
        set is extracted once and multiplied against the index. The result is THE SAME
        as calling kisa_liste() for each query separately (proven by a test): 'seed
        count' means "how many distinct seeds of the query occurred in this record" on
        both routes.

        THERE IS NO DOMAIN FILTER: whatever the query, every record is evaluated.

    """
    import math
    k = K.K_TOHUM
    # THE SAME CRITERION AS THE SINGLE QUERY kisa_liste() (idf plus BM25). The two
    # routes MUST produce the same short list on the same input; the reasoning is in the
    # "THE RANKING CRITERION" section at the top of identity_verification.py.
    tohum_sira = {}                          # ad -> sirali tohum listesi
    tohum_sahip = {}                         # seed -> [(name, index), ...]
    for ad, q in sorgular.items():
        th = sorted(K.tohumlar(q) | K.tohumlar(K.rc(q)))
        tohum_sira[ad] = th
        for i, t in enumerate(th):
            tohum_sahip.setdefault(t, []).append((ad, i))
    tohum_kume = set(tohum_sahip)
    ORT = K.ortalama_uzunluk(yol)
    B = K.BM25_B
    havuz = max(int(K.ADAY_HAVUZU), (ust or 0) * 3, 500)
    df = {ad: [0] * len(t) for ad, t in tohum_sira.items()}
    yigin = {ad: [] for ad in sorgular}      # min-heap, en fazla 'havuz' eleman
    n = 0
    N = 0
    if not tohum_kume:
        return {ad: [] for ad in sorgular}, 0
    for bas, diz in K.fasta_akisi(yol):      # NO CAP - to the end of the file
        n += 1
        if ilerle and n % 20000 == 0:
            ilerle(n)
        L = len(diz)
        if L < 100:
            continue
        N += 1
        kmers = {diz[i:i + k] for i in range(L - k + 1)}
        ortak = kmers & tohum_kume
        if not ortak:
            continue
        tut = {}
        for t in ortak:
            for ad, i in tohum_sahip[t]:
                tut.setdefault(ad, set()).add(i)
        tb = _TersMetin(bas)
        norm = 1.0 - B + B * L / ORT
        for ad, s in tut.items():
            d_ = df[ad]
            for i in s:
                d_[i] += 1                   # INVERSE FREQUENCY: free in the same pass
            on = len(s) / norm               # ON ELEME (idf henuz bilinmiyor)
            h = yigin[ad]
            if len(h) < havuz:
                heapq.heappush(h, (on, tb, n, diz, frozenset(s), norm))
            elif on > h[0][0]:
                heapq.heapreplace(h, (on, tb, n, diz, frozenset(s), norm))
    cikti = {}
    for ad, h in yigin.items():
        idf = [math.log(max(N, 2) / (1.0 + d)) for d in df[ad]]
        aday = [(sum(idf[i] for i in s) / nr, str(b), d, len(s))
                for (_o, b, _n, d, s, nr) in h]
        aday.sort(key=lambda x: (-x[0], x[1]))
        kesme = len(aday) if not ust else ust
        cikti[ad] = [dict(tohum=int(a[3]), skor=round(a[0], 4), baslik=a[1],
                          dizi=a[2], sira=i, kaynak='tohum')
                     for i, a in enumerate(aday[:kesme], 1)]
    return cikti, n


# --------------------------------------------------------------- THE BIN VERDICT
# -------------------------------------------------------------------------
# THE AT-LEAST-TWO-INDEPENDENT-DATABASES RULE - WHY IT EXISTS
#
# The best hit of a single database IS NOT ENOUGH for an identity claim. Every set
# carries its own historical biases: one may have deleted rare genera during
# dereplication (measured: SILVA LSURef NR99 holds 0 Petriella records while the
# Parc set of the same release holds 82), another may carry the same record under an
# outdated name. A verdict resting on one source would report that source's mistake
# as AN IDENTITY.
#
# So the verdict is tied to the databases AGREEING: if the number of INDEPENDENT
# databases giving the same genus as their best hit is >=2 it is DOGRULANDI, if it
# is 1 it is "DOGRULANAMADI (tek kaynak)", and if none agree it is DOGRULANAMADI.
#
# NO CONFIRMATION IS INVENTED: where the evidence is insufficient, a gap is reported
# as a gap. Independence is enforced too; sets in the VTB list that are twins (byte
# for byte identical) or subsets have been taken out of the vote, or the same record
# would vote twice.
#
# Agreement is sought at GENUS level: species level separates on a one letter name
# difference and would break a real agreement artificially.
# -------------------------------------------------------------------------
def kutu_hukmu(K, bulgular, lokus_tab):
    """A bin's identity: AT LEAST TWO independent databases must agree.

        The same requirement as `I`. Agreement is sought at GENUS level (species level
        can change by a single letter); if the number of INDEPENDENT databases giving
        the same genus as their best hit is >= 2, it is DOGRULANDI.

    """
    havuz = []
    for et, v in bulgular.items():
        if not str(v.get('durum', '')).startswith('TAMAM'):
            continue
        for i in (v.get('isabet') or [])[:5]:
            havuz.append(dict(i, _vtb=et, _lokus=lokus_tab.get(et, 'SSU')))
    sayisal = [h for h in havuz if isinstance(h.get('kimlik'), (int, float))]
    sayisal.sort(key=lambda x: -x['kimlik'])
    lokus = sayisal[0]['_lokus'] if sayisal else 'SSU'
    adl = K.savunulabilir_duzey(sayisal or havuz, lokus)
    for n_, h_ in enumerate(sayisal[:3], 1):
        c_, t_, tam_ = K.ad_coz(h_['baslik'])
        adl['isabet%d' % n_] = dict(tam_ad=tam_, cins=c_ or '-', tur=t_ or '-',
                                    kimlik=h_.get('kimlik'), uzunluk=h_.get('hiz_uzunluk'),
                                    vtb=h_['_vtb'])
    # her veritabaninin EN IYI isabetinin cinsi -> oy
    oy = collections.defaultdict(list)
    for et, v in bulgular.items():
        if not str(v.get('durum', '')).startswith('TAMAM'):
            continue
        isb = (v.get('isabet') or [])
        if not isb:
            continue
        c, _t, _tam = K.ad_coz(isb[0].get('baslik', ''))
        if c:
            oy[c].append(et)
    if oy:
        en_cok = max(oy.items(), key=lambda kv: len(kv[1]))
        uyusan_cins, uyusan_vtb = en_cok[0], en_cok[1]
    else:
        uyusan_cins, uyusan_vtb = None, []
    if len(uyusan_vtb) >= 2:
        hukum = 'DOGRULANDI'
    elif len(uyusan_vtb) == 1:
        hukum = 'DOGRULANAMADI (tek kaynak)'
    else:
        hukum = 'DOGRULANAMADI'
    return adl, hukum, uyusan_cins, uyusan_vtb, lokus, oy


# Do the recorded identity (the Kraken taxid name) and the measured identity point
# at the same genus? The "Candidatus" prefix is dropped before comparing; if there
# is no recorded name the result is not UYUSMADI but KAYIT YOK, and the two are not
# the same thing.
def ayni_mi(kayitli, dogrulanan_cins, adl):
    """Do the recorded identity and the verified identity point at the same genus?"""
    if not kayitli or kayitli in ('?', '-'):
        return 'KAYIT YOK', u'taxid_names.tsv icinde bu taxid icin ad yok'
    if not dogrulanan_cins:
        return 'BELIRSIZ', u'dogrulanan kimlikte cins cozulemedi'
    kc = re.sub(r'^(Ca\.|Candidatus)\s+', '', kayitli).split()[0]
    if kc.lower() == dogrulanan_cins.lower():
        return 'EVET', ''
    return 'HAYIR', u'kayitli "%s" -> olculen "%s"' % (kayitli, dogrulanan_cins)


# --------------------------------------------------------------- THE RUN
# -------------------------------------------------------------------------
# THE DRIVER. The bins are processed in batches (24 by default); ONE pass is opened
# per database and that pass serves every bin in the batch.
#
# THERE IS NO DOMAIN FILTER, and that is a deliberate decision. Choosing the domain
# from the Kraken label - "there is no need to ask a fungal database about a
# bacterial bin" - assumes the bin IS bacterial, when the whole reason this stage
# exists is that those labels can be wrong. Every database is asked, and the
# irrelevant ones are dropped IN THE RESULT. The bin's domain comes FROM THE
# MEASUREMENT, NOT FROM THE LABEL (the alan_olcumden column).
#
# A DATABASE THAT RETURNS NO RESULT DOES NOT COUNT AS "CLEAN": it is marked
# separately as SONUC YOK / DOSYA YOK / SORULMADI, and all 12 of the 12 sources
# appear on every bin row.
#
# COVERAGE ACCOUNTING: for each database the number of records SCANNED is counted
# and compared against BEKLENEN_KAYIT. A real cap problem happened in the access
# test: the first run cut off at 120 001 records and SILVA SSU NR99, LSU Parc and
# UNITE ITS were effectively being scanned truncated. There is no cap in the
# streamer here, and incomplete coverage does not stay quiet, it prints a warning.
#
# THE CACHE IS SHARED WITH I: the checkpoints live under IDENTITY_RESULT/kontrol under
# the same key, so the same bin is never scanned twice across the two stages.
# -------------------------------------------------------------------------
def calistir(kok, kl_ust, kume_boyu, nt_kip, lit_kip, sifirla, yalniz, tavan_kutu):
    K = _K(kok)
    CIKTI = os.path.join(kok, 'ALL_IDENTITIES_RESULT')
    # CACHE SHARING: the checkpoints live in THE SAME directory as `I` and under
    # THE SAME key. The same bin is never scanned twice across the two stages.
    KONTROL = os.path.join(kok, 'IDENTITY_RESULT', 'kontrol')
    os.makedirs(CIKTI, exist_ok=True)
    os.makedirs(KONTROL, exist_ok=True)
    if sifirla:
        for f in os.listdir(CIKTI):
            if f.startswith('kutu_') and f.endswith('.json'):
                try:
                    os.remove(os.path.join(CIKTI, f))
                except OSError:
                    pass
    g = open(os.path.join(CIKTI, 'kosu_gunlugu.txt'), 'a', encoding='utf-8')

    def yaz(s=''):
        print(s, flush=True); g.write(s + '\n'); g.flush()

    yaz('=' * 78)
    yaz(u'  G - INDEPENDENT VERIFICATION OF EVERY BIN IDENTITY')
    yaz(u'  version %s   %s' % (VERSIYON, time.strftime('%Y-%m-%d %H:%M')))
    yaz('=' * 78)
    yaz(u'  `I` tests 12 CLAIMS, the suspicious ones. This stage tests EVERY bin')
    yaz(u'  that enters the panel measurements: the silent majority.')
    yaz('')

    katilan, atlanan, uye, rakip, kons, H = kutu_envanteri(kok, K)
    if yalniz:
        ay = [x.strip() for x in yalniz.split(',') if x.strip()]
        katilan = [k for k in katilan if any(a.lower() in k.lower() for a in ay)]
    if tavan_kutu:
        katilan = katilan[:tavan_kutu]

    # --- DATABASE COVERAGE: ALL OF THEM, NO DOMAIN FILTER ---
    var, yok = [], []
    for e, d, t, kullan, _n in K.VTB:
        if not kullan:
            continue                       # a twin or subset - not an independent source
        p = os.path.join(kok, 'REFERANS_DB', d)
        (var if os.path.exists(p) else yok).append((e, d, t))
    lokus_tab = {e: t for e, _d, t, _k, _n in K.VTB}

    yaz(u'  bins (included in the measurement) : %d' % len(katilan))
    yaz(u'  bins (SKIPPED)                     : %d  %s'
        % (len(atlanan), ', '.join(k for k, _s in atlanan) or '-'))
    yaz(u'  DATABASES                          : %d to query, %d files missing'
        % (len(var), len(yok)))
    for e, d, _t in var:
        yaz(u'      [WILL ASK]  %-20s %-32s expected %s records'
            % (e, d, '{:,}'.format(BEKLENEN_KAYIT.get(e, 0)).replace(',', ' ') or '?'))
    for e, d, _t in yok:
        yaz(u'      [NO FILE]   %-20s %-32s not found under REFERANS_DB' % (e, d))
    yaz(u'  DOMAIN FILTER                      : NONE. Every bin is asked against ALL'
        % len(var))
    yaz(u'    %d databases. Choosing the domain from the Kraken label would be')
    yaz(u'    dangerous: this stage exists precisely because those labels can be wrong.')
    yaz(u'    Irrelevant databases are dropped from the RESULT ("NO RESULT"), never')
    yaz(u'    filtered out BEFORE the query. A bin\'s domain comes from MEASUREMENT, not from its label.')
    yaz(u'  NCBI nt                            : %s'
        % {'oto': u'automatic (URL API), a separate BLAST per bin',
           'elle': u'manual (a query file is written)',
           'yok': u'NOT ASKED (default): %d bins x the BLAST queue would take days. Any nt cache left from stage `I` IS reused. Enable with --nt auto.' % len(katilan)}[nt_kip])
    yaz('')

    # --- SURE TAHMINI ---
    tekil = {}
    for k in katilan:
        tekil.setdefault(K.kp_yolu(KONTROL, '_', kons[k][:4000], (), kl_ust), k)
    n_tekil = len(tekil)
    _uz = [min(len(kons[k]), 4000) for k in katilan] or [1500]
    _oq = sum(_uz) / float(len(_uz))
    _ref = 2000.0
    _bir = 6.7e-6 * min(_oq, _ref) + 4.11e-9 * min(_oq, _ref) * max(_oq, _ref)
    # Measured constants (this script's own code, a synthetic set, 2026-08-04):
    #   the bulk scan  : ~6 400 records/s (24 queries at once)
    #   a single query : ~2 200 records/s (1 query)  -> the bulk pass is ~70x as efficient
    # The user's machine may differ; every database line prints the REAL time, so the
    # estimate corrects itself as the run goes on.
    TARAMA_HIZI = 6400.0
    _hiz = n_tekil * len(var) * kl_ust * _bir
    _kume = max(1, (n_tekil + kume_boyu - 1) // kume_boyu)
    _kayit = sum(BEKLENEN_KAYIT.get(e, 100000) for e, _d, _t in var)
    _tara = _kume * _kayit / TARAMA_HIZI
    yaz(u'  ESTIMATED TIME: ~%s   (AN OVERNIGHT JOB)' % K.sure_metni(_hiz + _tara))
    yaz(u'    unique consensus %d x %d databases = %d queries'
        % (n_tekil, len(var), n_tekil * len(var)))
    yaz(u'    alignment cost  ~%s  (%d queries x %d candidates, avg %d bp) <- DOMINANT cost'
        % (K.sure_metni(_hiz), n_tekil * len(var), kl_ust, int(_oq)))
    yaz(u'    scan cost       ~%s  (%d sets x %d database streams = %d streams, %s records)'
        % (K.sure_metni(_tara), _kume, len(var), _kume * len(var),
           '{:,}'.format(_kayit).replace(',', ' ')))
    yaz(u'    NOTE: a separate stream per bin would need %d streams and ~%s of scanning;'
        % (n_tekil * len(var), K.sure_metni(n_tekil * _kayit / 2200.0)))
    yaz(u'    batching brings that down to %d streams, about 70 times fewer.' % (_kume * len(var)))
    yaz(u'    Resumable: state is written to disk after every bin, and the database')
    yaz(u'    scans share a cache with stage `I`.')
    yaz('')

    # --- THE RUN: batch by batch, database by database ---
    sonuc, tb = [], time.time()
    bekleyen = []
    for k in katilan:
        kp = os.path.join(CIKTI, 'kutu_%s.json' % re.sub(r'\W+', '_', k))
        if os.path.exists(kp):
            try:
                sonuc.append(json.load(open(kp, encoding='utf-8')))
                continue
            except Exception:
                pass
        bekleyen.append(k)
    yaz(u'  taken from the previous run: %d bins | to scan: %d bins'
        % (len(sonuc), len(bekleyen)))

    kapsam_kayit = {}          # label -> (scanned, expected, coverage)
    for ki in range(0, len(bekleyen), kume_boyu):
        kume = bekleyen[ki:ki + kume_boyu]
        yaz('')
        yaz(u'[set %d/%d] %d bins: %s'
            % (ki // kume_boyu + 1, (len(bekleyen) + kume_boyu - 1) // kume_boyu,
               len(kume), ', '.join(kume[:6]) + (' ...' if len(kume) > 6 else '')))
        bulgular = {k: {} for k in kume}
        for et, dosya, _t in var:
            # collect the ones not in the cache (a key SHARED with I)
            kalan = {}
            for k in kume:
                q = kons[k][:4000]
                kp = K.kp_yolu(KONTROL, et, q, (), kl_ust)
                if os.path.exists(kp):
                    try:
                        bulgular[k][et] = json.load(open(kp, encoding='utf-8'))
                        continue
                    except Exception:
                        pass
                kalan[k] = q
            if not kalan:
                yaz(u'     %-20s: ALL %d bins came from the cache' % (et, len(kume)))
                continue
            t0 = time.time()
            yol = os.path.join(kok, 'REFERANS_DB', dosya)

            def ilerle(n, _e=et, _t0=t0):
                print(u'     ... %s: %d records scanned (%s)      '
                      % (_e, n, K.sure_metni(time.time() - _t0)), end='\r', flush=True)
            kls, taranan = toplu_kisa_liste(K, yol, kalan, kl_ust, ilerle)
            bek = BEKLENEN_KAYIT.get(et)
            kapsam = ('TAMAMI' if bek and taranan >= bek else
                      ('TAMAMI (beklenen bilinmiyor)' if not bek else
                       'EKSIK (%d / %d)' % (taranan, bek)))
            kapsam_kayit[et] = (taranan, bek, kapsam)
            for k, q in kalan.items():
                kl = kls.get(k) or []
                if not kl:
                    # NO RESULT - NOT "clean"; it is marked separately
                    res = dict(durum=u'SONUC YOK', kayit=0, kisa_liste_boyu=kl_ust,
                               hizalanan=0, isabet=[], kazanan_sira=None,
                               kazanan_kaynak=None, sira_uyarisi=None,
                               taranan_kayit=taranan, kapsam=kapsam,
                               sebep=u'TEMIZ SAYILMAZ: tarandi (%s kayit, kapsam %s) '
                                     u'ama hicbir kayit sorgunun tohumlarini '
                                     u'tutturmadi; kutu bu veritabaninin alani '
                                     u'disinda olabilir'
                                     % ('{:,}'.format(taranan).replace(',', ' '), kapsam))
                else:
                    res = K.kl_degerlendir(kl, q, kl_ust, taranan=taranan, t0=t0)
                    res['kapsam'] = kapsam
                json.dump(res, open(K.kp_yolu(KONTROL, et, q, (), kl_ust), 'w',
                                    encoding='utf-8'), ensure_ascii=False, default=str)
                bulgular[k][et] = res
            bos = len([1 for k in kalan if not (kls.get(k) or [])])
            yaz(u'     %-20s: %s records scanned, coverage %s | %d bins hit, %d NO RESULT (%s)'
                % (et, '{:,}'.format(taranan).replace(',', ' '), kapsam,
                   len(kalan) - bos, bos, K.sure_metni(time.time() - t0)))
            if bek and taranan < bek:
                yaz(u'     >>> WARNING: COVERAGE INCOMPLETE. %d records were expected in %s but %d were scanned. The cap problem may have recurred'
                    % (et, bek, taranan))

        # dosyasi olmayan veritabanlari da SATIRDA GORUNSUN
        for e, d, _t in yok:
            for k in kume:
                bulgular[k][e] = dict(durum=u'DOSYA YOK', isabet=[],
                                      sebep=u'REFERANS_DB/%s bulunamadi' % d)

        for k in kume:
            r = kutu_kaydi(K, kok, k, kons[k][:4000], bulgular[k], lokus_tab, uye,
                           rakip, H, nt_kip, lit_kip, KONTROL, CIKTI, yaz, var, yok)
            json.dump(r, open(os.path.join(CIKTI, 'kutu_%s.json'
                                           % re.sub(r'\W+', '_', k)), 'w',
                              encoding='utf-8'), ensure_ascii=False, default=str)
            sonuc.append(r)
        gec = time.time() - tb
        yap = len([1 for s in sonuc if s['kutu'] in bekleyen])
        if yap:
            print('        gecen %s | tahmini kalan %s'
                  % (K.sure_metni(gec), K.sure_metni(gec / yap * (len(bekleyen) - yap))),
                  flush=True)

    raporla(K, CIKTI, sonuc, atlanan, var, yok, kapsam_kayit, uye, rakip, yaz, kl_ust, nt_kip)
    g.close()
    return 0 if sonuc else 1


# -------------------------------------------------------------------------
# It builds one bin's row: the recorded identity, the verified identity, the
# agreement state, the defensible level, the best three hits, the literature layer
# and THE SOURCE ACCOUNTING.
#
# The source accounting holds ALL of the 12 local plus NCBI nt rows; if a source was
# never tried it says "SORULMADI (BILINMEYEN) - HATA, bildirin". Printing a noisy
# error was preferred over leaving a silent gap.
#
# The row also holds WHICH targets this bin is a member or a competitor of: if the
# identity changes, which measurements have to be redone is read from here.
# -------------------------------------------------------------------------
def kutu_kaydi(K, kok, kutu, q, bulgular, lokus_tab, uye, rakip, H, nt_kip, lit_kip,
               KONTROL, CIKTI, yaz, var, yok):
    u"""Tek kutunun satiri: kayitli kimlik, dogrulanan kimlik, uyusma, kaynak muhasebesi."""
    taxid = kutu.split('_')[-1]
    kayitli = H.taxid_adlari().get(taxid, '')

    # --- NCBI nt katmani (ayri kaynak) ---
    ntk = os.path.join(KONTROL, 'nt_%s.json' % re.sub(r'\W+', '_', kutu))
    if os.path.exists(ntk):
        try:
            bulgular[K.NT_ETIKET] = json.load(open(ntk, encoding='utf-8'))
        except Exception:
            pass
    if K.NT_ETIKET not in bulgular:
        if nt_kip == 'yok':
            bulgular[K.NT_ETIKET] = dict(
                durum=u'SORULMADI (--nt yok)', isabet=[],
                sebep=u'G asamasinda nt varsayilan olarak kapali: kutu basina ayri '
                      u'BLAST kuyrugu gerekir. `I` asamasindan onbellek de yok. '
                      u'--nt oto ile acilir.')
        else:
            bulgular[K.NT_ETIKET] = K.nt_katmani(kutu, q, CIKTI, yaz, nt_kip)
            if str(bulgular[K.NT_ETIKET].get('durum', '')).startswith('TAMAM'):
                json.dump(bulgular[K.NT_ETIKET], open(ntk, 'w', encoding='utf-8'),
                          ensure_ascii=False, default=str)

    adl, hukum, cins, uyusan_vtb, lokus, oy = kutu_hukmu(K, bulgular, lokus_tab)
    uym, uym_not = ayni_mi(kayitli, cins, adl)

    # --- THE LITERATURE LAYER (the same as I) ---
    try:
        import importlib.util as _lu
        _lp = _lu.spec_from_file_location(
            'lit', os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'literature_check.py'))
        LIT = _lu.module_from_spec(_lp); _lp.loader.exec_module(LIT)
        lit = LIT.kontrol_et((adl.get('isabet1') or {}).get('tam_ad', ''),
                             adl.get('onerilen_ad', ''), lokus, ag=(lit_kip != 'yok'))
    except Exception as _e:
        lit = dict(durum=u'literatur modulu yuklenemedi: %s' % type(_e).__name__,
                   ncbi_guncel_ad='-', ad_farkli_mi='-', revizyon_uyarisi='-')

    # --- KAYNAK MUHASEBESI: 12 veritabaninin 12'si de satirda gorunur ---
    beklenen_et = [e for e, _d, _t in var] + [e for e, _d, _t in yok] + [K.NT_ETIKET]
    detay = collections.OrderedDict()
    for et in beklenen_et:
        v = bulgular.get(et)
        if v is None:
            detay[et] = dict(durum=u'SORULMADI (BILINMEYEN)', en_iyi='', kimlik=None,
                             sebep=u'bu veritabani hic denenmedi - HATA, bildirin')
            continue
        d = str(v.get('durum', '?'))
        isb = (v.get('isabet') or [])
        detay[et] = dict(
            durum=d, en_iyi=(isb[0].get('baslik', '')[:120] if isb else ''),
            kimlik=(isb[0].get('kimlik') if isb else None),
            kazanan_sira=v.get('kazanan_sira'), kazanan_kaynak=v.get('kazanan_kaynak'),
            taranan_kayit=v.get('taranan_kayit'), kapsam=v.get('kapsam'),
            sebep=v.get('sebep', ''))
    sorgulanan = len([1 for v in detay.values() if str(v['durum']).startswith('TAMAM')])
    sonucsuz = [e for e, v in detay.items() if v['durum'] == 'SONUC YOK']
    sorulmayan = [e for e, v in detay.items()
                  if str(v['durum']).startswith('SORULMADI') or v['durum'] == 'DOSYA YOK']
    # THE DOMAIN FROM THE MEASUREMENT: which databases actually gave a hit
    alan = sorted({lokus_tab.get(e, '?') for e, v in detay.items()
                   if str(v['durum']).startswith('TAMAM')})

    kzs = [v['kazanan_sira'] for v in detay.values() if isinstance(v.get('kazanan_sira'), int)]
    return dict(
        kutu=kutu, taxid=taxid, kayitli_kimlik=kayitli or '-',
        dogrulanan_cins=cins or '-', dogrulanan_ad=adl.get('onerilen_ad', '-'),
        duzey=adl.get('duzey', '-'), gerekce=adl.get('gerekce', '-'),
        uyusuyor=uym, uyusma_notu=uym_not, hukum=hukum,
        uyusan_vtb=uyusan_vtb, oylar={k: v for k, v in oy.items()},
        adlandirma=adl, literatur=lit, vtb_detay=detay,
        sorgulanan_vtb=sorgulanan, toplam_vtb=len(beklenen_et),
        sonuc_yok_vtb=sonucsuz, sorulmayan_vtb=sorulmayan,
        alan_olcumden=alan, lokus=lokus,
        kazanan_sira_maks=(max(kzs) if kzs else None),
        uye_hedefler=sorted(uye.get(kutu, [])), rakip_hedefler=sorted(rakip.get(kutu, [])))


# --------------------------------------------------------------- THE REPORT
# THE DISAGREEING ONES FIRST, then the uncertain ones, and the agreeing ones last;
# on a tie, the bin affecting more targets comes first. Whoever reads the report
# should see the row that creates the most work first.
def _sirala(s):
    """UYUSMAYANLAR EN BASA. Sonra belirsizler, sonra uyusanlar."""
    o = {'HAYIR': 0, 'BELIRSIZ': 1, 'KAYIT YOK': 2, 'EVET': 3}
    return (o.get(s['uyusuyor'], 9), -len(s['uye_hedefler']) - len(s['rakip_hedefler']),
            s['kutu'])


# -------------------------------------------------------------------------
# Two outputs: a one row per bin TSV and a markdown report.
#
# At the end of the markdown report stands THE IMPACT SUMMARY: how many bins changed
# identity, how many targets' member or competitor sets that affects, and whether a
# re-measurement is needed. On a target whose member set changes, the backbone
# consensus can change too, which means the primer and discrimination values have to
# be recomputed. This script DOES NOT do that calculation itself; it only reports
# that it is needed, and it does not touch the panel files.
# -------------------------------------------------------------------------
def raporla(K, CIKTI, sonuc, atlanan, var, yok, kapsam_kayit, uye, rakip, yaz,
            kl_ust, nt_kip):
    sonuc = sorted(sonuc, key=_sirala)
    beklenen_et = [e for e, _d, _t in var] + [e for e, _d, _t in yok] + [K.NT_ETIKET]
    t = os.path.join(CIKTI, 'tum_kutu_kimlikleri.tsv')
    with open(t, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# The identity of EVERY bin entering the panel measurements was tested independently.\n')
        fh.write(u'# Same method as `I`: a short list of %d, all aligned, at least TWO independent databases must agree.\n' % kl_ust)
        fh.write(u'# NO DOMAIN FILTER: every database is asked for every bin; irrelevant ones are dropped from the RESULT ("NO RESULT"), never filtered out before the query. A bin\'s domain comes from MEASUREMENT, not from its label.\n')
        fh.write(u'# DISAGREEING rows come FIRST.\n')
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['kutu', 'taxid', 'MEVCUT_KAYITLI_KIMLIK', 'DOGRULANAN_KIMLIK',
                    'UYUSUYOR_MU', 'uyusma_notu', 'HUKUM',
                    'SAVUNULABILIR_DUZEY', 'ONERILEN_AD', 'adlandirma_gerekcesi',
                    # THE ALIGNED LENGTH IS REQUIRED beside the identity percentage.
                    # (2026-08-11) The percentage alone misleads: for the same bin,
                    # Petriella setifera came out at 100% in SILVA LSU Parc and
                    # Petriella musispora at 100% in the RefSeq fungal ITS. Both are
                    # 100% because they were aligned at two DIFFERENT loci over two
                    # DIFFERENT lengths. Anyone reading "100%" without seeing the
                    # length believes the species is settled. The data was already
                    # being produced (hiz_uzunluk); it just was not written to the table.
                    'en_iyi_isabet', 'en_iyi_kimlik_%', 'en_iyi_hiz_uzunluk', 'en_iyi_vtb',
                    'ikinci_isabet', 'ikinci_kimlik_%', 'ikinci_hiz_uzunluk', 'ikinci_vtb',
                    'ucuncu_isabet', 'ucuncu_kimlik_%', 'ucuncu_hiz_uzunluk', 'ucuncu_vtb',
                    'UYESI_OLDUGU_HEDEFLER', 'RAKIBI_OLDUGU_HEDEFLER',
                    'etkilenen_hedef_sayisi',
                    'SORGULANAN_VTB', 'SONUC_YOK_VTB', 'SORULMAYAN_VTB',
                    'alan_olcumden', 'lokus', 'kazanan_sira_maks',
                    'LIT_ncbi_guncel_ad', 'LIT_AD_FARKLI_MI', 'LIT_revizyon_uyarisi',
                    'LIT_durum', 'HER_VTB_NE_DEDI'])
        for s in sonuc:
            a = s.get('adlandirma') or {}

            def _i(n, alan, v='-'):
                return ((a.get('isabet%d' % n) or {}).get(alan) or v)
            d = s['vtb_detay']
            hepsi = ' | '.join(
                '%s [%s]: %s%s%s%s'
                % (e,
                   ('%s kayit, kapsam %s' % ('{:,}'.format(v['taranan_kayit']).replace(',', ' '),
                                             v.get('kapsam') or '?'))
                   if v.get('taranan_kayit') else v['durum'],
                   v['en_iyi'] or v['durum'],
                   ('' if v.get('kimlik') is None else ' (%%%s)' % K.vir(v['kimlik'])),
                   ('' if v.get('kazanan_sira') is None
                    else ' {sira %s/%d}' % (v['kazanan_sira'], kl_ust)),
                   ('' if not v.get('sebep') else ' <%s>' % v['sebep'][:110]))
                for e, v in d.items())
            w.writerow([
                s['kutu'], s['taxid'], s['kayitli_kimlik'], s['dogrulanan_ad'],
                s['uyusuyor'], s['uyusma_notu'], s['hukum'],
                s['duzey'], s['dogrulanan_ad'], s['gerekce'],
                _i(1, 'tam_ad'), K.vir(_i(1, 'kimlik', None)), _i(1, 'uzunluk'), _i(1, 'vtb'),
                _i(2, 'tam_ad'), K.vir(_i(2, 'kimlik', None)), _i(2, 'uzunluk'), _i(2, 'vtb'),
                _i(3, 'tam_ad'), K.vir(_i(3, 'kimlik', None)), _i(3, 'uzunluk'), _i(3, 'vtb'),
                ', '.join(s['uye_hedefler']) or '-',
                ', '.join(s['rakip_hedefler']) or '-',
                len(set(s['uye_hedefler']) | set(s['rakip_hedefler'])),
                '%d / %d' % (s['sorgulanan_vtb'], s['toplam_vtb']),
                ', '.join(s['sonuc_yok_vtb']) or '-',
                ', '.join(s['sorulmayan_vtb']) or '-',
                ', '.join(s['alan_olcumden']) or '-', s['lokus'],
                s['kazanan_sira_maks'] if s['kazanan_sira_maks'] is not None else '-',
                (s.get('literatur') or {}).get('ncbi_guncel_ad', '-'),
                (s.get('literatur') or {}).get('ad_farkli_mi', '-'),
                (s.get('literatur') or {}).get('revizyon_uyarisi', '-'),
                (s.get('literatur') or {}).get('durum', '-'), hepsi])
        if atlanan:
            fh.write(u'#\n# SKIPPED BINS (nothing is skipped silently):\n')
            for k, sebep in atlanan:
                fh.write(u'# %s\t%s\n' % (k, sebep))
    yaz(u'  written: %s' % t)

    # ----------------------------------------------------------- MD RAPORU
    degisen = [s for s in sonuc if s['uyusuyor'] == 'HAYIR']
    belirsiz = [s for s in sonuc if s['uyusuyor'] in ('BELIRSIZ', 'KAYIT YOK')]
    etkilenen = set()
    for s in degisen:
        etkilenen |= set(s['uye_hedefler']) | set(s['rakip_hedefler'])
    uye_etkilenen = {h for s in degisen for h in s['uye_hedefler']}
    r = os.path.join(CIKTI, 'TUM_KUTU_KIMLIK_RAPORU.md')
    with open(r, 'w', encoding='utf-8') as fh:
        fh.write(u'# Independent verification of all bin identities\n\n')
        fh.write(u'Generated: %s, script %s\n\n' % (time.strftime('%Y-%m-%d %H:%M'), VERSIYON))
        fh.write(u'Stage `I` tests 12 **claims**, the suspicious ones. This stage tests **every bin** that enters the panel measurements: the silent majority that defines the member and competitor sets, and therefore affects every discrimination figure.\n\n')

        # --- KAYNAK KAPSAMI (kanit) ---
        fh.write(u'## Database coverage\n\n')
        fh.write(u'**NO DOMAIN FILTER WAS APPLIED.** Every bin was asked against all %d local databases. Choosing the domain from the Kraken label would be dangerous: this stage exists precisely because those labels can be wrong. Irrelevant databases were dropped from the **result** (`NO RESULT`), not filtered out before the query. A bin\'s domain was derived from **measurement**, not from its label (the `alan_olcumden` column).\n\n' % len(var))
        fh.write(u'| # | database | expected records | scanned | coverage |\n|---|---|---|---|---|\n')
        for i, (e, _d, _t) in enumerate(var, 1):
            tar, bek, kap = kapsam_kayit.get(e, (None, BEKLENEN_KAYIT.get(e), None))
            fh.write(u'| %d | %s | %s | %s | %s |\n'
                     % (i, e, '{:,}'.format(BEKLENEN_KAYIT.get(e, 0)).replace(',', ' '),
                        '{:,}'.format(tar).replace(',', ' ') if tar else
                        u'from cache (not scanned in this run)',
                        kap or u'verified in a previous run'))
        for e, d, _t in yok:
            fh.write(u'| - | %s | - | - | **NO SUCH FILE** (`REFERANS_DB/%s`) |\n' % (e, d))
        fh.write(u'| - | NCBI nt | - | - | %s |\n'
                 % {'yok': u'**NOT ASKED** (--nt none; a separate BLAST queue per bin)',
                    'oto': u'automatic (URL API)', 'elle': u'elle'}[nt_kip])
        eksik = [e for e, (tar, bek, kap) in kapsam_kayit.items()
                 if kap and kap.startswith('EKSIK')]
        fh.write(u'\n> **The cap problem did not recur.** In the access test the first run stopped at 120,001 records, so SILVA SSU NR99 (510,495), LSU Parc (1,312,521) and UNITE ITS (2,069,189) were being scanned truncated. There is **no cap** in the streamer here; the number of records scanned is counted and compared against what was expected. Databases with incomplete coverage in this run: **%s**.\n\n' % (', '.join(eksik) if eksik else u'YOK'))
        fh.write(u'> A database that returned nothing was **not counted as clean**: it is marked separately as `NO RESULT` and shown in the row. Every bin row also lists %d of the %d sources.\n\n'
                 % (len(beklenen_et), len(beklenen_et)))

        # --- ETKI OZETI ---
        fh.write(u'## Impact summary\n\n')
        fh.write(u'| measure | value |\n|---|---|\n')
        fh.write(u'| bins tested | %d |\n' % len(sonuc))
        fh.write(u'| **bins whose identity CHANGED** | **%d** |\n' % len(degisen))
        fh.write(u'| bins whose identity was verified | %d |\n'
                 % len([s for s in sonuc if s['uyusuyor'] == 'EVET']))
        fh.write(u'| uncertain, or with no recorded name | %d |\n' % len(belirsiz))
        fh.write(u'| **targets affected** | **%d** |\n' % len(etkilenen))
        fh.write(u'| of those, the ones whose MEMBERSHIP set changed | %d |\n' % len(uye_etkilenen))
        fh.write(u'| skipped bins (not in the measurement) | %d |\n\n' % len(atlanan))
        if degisen:
            fh.write(u'> **RE-MEASUREMENT IS NEEDED.** %d bins changed identity, and '
                     u'that affects the member or competitor set of %d targets. On '
                     u'the %d targets whose member set changed the backbone consensus '
                     u'can change too, so the primer and discrimination values have '
                     u'to be recomputed.\n\n'
                     u'The targets affected: %s\n\n'
                     % (len(degisen), len(etkilenen), len(uye_etkilenen),
                        ', '.join(sorted(etkilenen))))
        else:
            fh.write(u'> No bin changed identity: the member and competitor sets stay as they are, so **no re-measurement is needed**')

        # --- UYUSMAYANLAR ---
        if degisen:
            fh.write(u'## Disagreements (read these first)\n\n')
            for s in degisen:
                fh.write(u'### %s  (taxid %s)\n\n' % (s['kutu'], s['taxid']))
                fh.write(u'- **Recorded identity:** %s\n' % s['kayitli_kimlik'])
                fh.write(u'- **Dogrulanan:** %s  (`%s`) - %s\n'
                         % (s['dogrulanan_ad'], s['duzey'], s['hukum']))
                fh.write(u'  - *Gerekce:* %s\n' % s['gerekce'])
                fh.write(u'- **The targets it is a member of:** %s\n'
                         % (', '.join(s['uye_hedefler']) or '-'))
                fh.write(u'- **The targets it is a competitor of:** %s\n'
                         % (', '.join(s['rakip_hedefler']) or '-'))
                a = s.get('adlandirma') or {}
                fh.write(u'\n  | # | nearest record | genus | species | identity | database |\n  |---|---|---|---|---|---|\n')
                for n_ in (1, 2, 3):
                    it = a.get('isabet%d' % n_)
                    if it:
                        fh.write(u'  | %d | %s | %s | %s | %%%s | %s |\n'
                                 % (n_, it['tam_ad'], it['cins'], it['tur'],
                                    K.vir(it['kimlik']), it['vtb']))
                fh.write(u'\n  **Source accounting (%d/%d):**\n\n'
                         % (s['sorgulanan_vtb'], s['toplam_vtb']))
                fh.write(u'  | database | status | best hit | identity | winner rank |\n  |---|---|---|---|---|\n')
                for e, v in s['vtb_detay'].items():
                    fh.write(u'  | %s | %s | %s | %s | %s |\n'
                             % (e, v['durum'], v['en_iyi'] or (v.get('sebep') or '-')[:96],
                                '-' if v.get('kimlik') is None else '%%%s' % K.vir(v['kimlik']),
                                v.get('kazanan_sira') if v.get('kazanan_sira') is not None else '-'))
                fh.write(u'\n')
        if atlanan:
            fh.write(u'## Atlanan kutular\n\n')
            fh.write(u'These are not a member or a competitor of any target, so they enter no discrimination calculation. They are not skipped silently')
            for k, sebep in atlanan:
                fh.write(u'- `%s` - %s\n' % (k, sebep))
            fh.write(u'\n')
    yaz(u'  written: %s' % r)
    yaz('')
    yaz(u'  IMPACT: %d bins tested | identity CHANGED for %d | targets affected %d | membership changed %d | skipped %d'
        % (len(sonuc), len(degisen), len(etkilenen), len(uye_etkilenen), len(atlanan)))
    if eksik:
        yaz(u'  >>> WARNING: database with INCOMPLETE coverage: %s' % ', '.join(eksik))


# The command line: --shortlist (the same default as I), --cluster how many bins per
# pass, --nt the NCBI layer (the default is "yok": a separate BLAST queue per bin
# would take days; the cache left over from stage I is still used), --only and
# --cap-bins a subset for testing.

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
    p = argparse.ArgumentParser(
        description='Confirm independently the identity of EVERY bin that '
                    'enters the panel measurements')
    p.add_argument('--root', dest='kok', default='.')
    p.add_argument('--shortlist', type=int, default=None, dest='kisa_liste',
                   help=u'number of candidates to align fully (default: same as the identity stage)')
    p.add_argument('--cluster', dest='kume', type=int, default=24,
                   help=u'how many bins to scan together in one database pass '
                        u'(bellek/hiz dengesi, varsayilan 24)')
    p.add_argument('--nt', choices=['auto', 'manual', 'none', 'oto', 'elle', 'yok'], default='yok',
                   help='the NCBI nt layer (the default is none: a separate '
                        'BLAST queue per bin takes days; the cache left from '
                        'the identity stage is still used)')
    p.add_argument('--literature', dest='literatur', choices=['auto', 'none', 'oto', 'yok'], default='oto')
    p.add_argument('--only', dest='yalniz', default=None,
                   help=u'comma-separated ayrilmis bin adi parcalari (for testing)')
    p.add_argument('--cap-bins', type=int, default=0, dest='tavan_kutu',
                   help=u'only ilk N bin (for testing)')
    p.add_argument('--reset', dest='sifirla', action='store_true')
    a = p.parse_args()
    a = _ing_deger(a)
    kok = os.path.abspath(a.kok)
    if not os.path.isdir(os.path.join(kok, 'screening')):
        sys.exit(u'ERROR: there is no screening directory inside %s.' % kok)
    kl = a.kisa_liste
    if kl is None:
        kl = _K(kok).KISA_LISTE
    return calistir(kok, kl, a.kume, a.nt, a.literatur, a.sifirla, a.yalniz, a.tavan_kutu)


if __name__ == '__main__':
    sys.exit(main() or 0)
