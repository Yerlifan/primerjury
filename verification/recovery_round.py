# -*- coding: utf-8 -*-
"""THE RECOVERY ROUND - it tries to recover the rows the SINGLE PROTOCOL run left
below the threshold, by five separate routes.

THE BASIC FRAME
---------------
The target chosen in the meeting is really A BIN; the name is that day's Kraken
label and IT MAY BE WRONG. The aim is to amplify those bins. So a target's member
set = the bin that was pointed at, plus every bin the measurement shows to be THE
SAME ORGANISM. That is what rescued Petriella: the organism was spread over nine
bins, and because the others sat in the "competitor" column the metric was
comparing the target against itself.

FIVE ROUTES (all in this one option, in order)
  1) fix THE MEASURE on universal targets   (the discrimination ratio is undefined
     on those rows)
  2) NARROW THE MEMBERSHIP                  (every target whose coverage is partial)
  3) REDESIGN plus ARMS                     (the ones that fell just short)
  4) clean up rows that still have A COUNTERPART (a better pair took their place)
  5) A MULTI LOCUS search                   (split the whole consensus into regions)

THE THRESHOLD IS NOT LOWERED. Getting a row through by lowering the threshold is
FORBIDDEN. Fixing the measure (route 1) and loosening the threshold ARE DIFFERENT
THINGS; this script does only the first, and it writes on every row which measure
was used and why.

It WRITES NOTHING into the panel files. It only reads, and writes under
KURTARMA_SONUC/.

"""

# -------------------------------------------------------------------------
# recovery_round.py takes the rows the SINGLE PROTOCOL (P) run left below the
# threshold and tries to recover them by five separate routes. THE THRESHOLD IS
# NEVER LOWERED.
#
# INPUT  : TEK_PROTOKOL_SONUC/panel_tek_protokol.tsv (the below-threshold rows),
#          uyelik_yeniden_turetme_uyelik_*.tsv (stage U's measured membership),
#          protocol/ek_ciftler.tsv (target name aliases),
#          konsensus_kanonik/ and "fastq files"/ (the measurement sources).
# OUTPUT : KURTARMA_SONUC/kurtarma_satirlari.tsv (ONE row per target),
#          KURTARMA_SONUC/yeni_adaylar.tsv (the candidates from the route 3 and
#          route 5 scans),
#          KURTARMA_SONUC/KURTARMA_RAPORU.md, kontrol/ (one JSON per target).
#          It WRITES NOTHING into the panel files.
# CALLED BY: verification/full_chain.py -> key K
#          (python3 verification/recovery_round.py --root .)
#
# FIVE ROUTES, tried in order:
#   route 1  fix THE MEASURE on universal targets (the discrimination ratio is
#            undefined there)
#   route 2  NARROW THE MEMBERSHIP (by measured sequence identity, UNCONDITIONALLY)
#   route 3  REDESIGN plus ARMS (for rows that fell just short)
#   route 4  move rows that still have A COUNTERPART onto the ones that dropped
#   route 5  A MULTI LOCUS search (split the whole consensus into regions and
#            design separately in each; this is what rescued Petriella, where the
#            answer was not in ITS but in LSU)
#
# FIXING THE MEASURE and LOOSENING THE THRESHOLD ARE DIFFERENT THINGS. This script
# does only the first, and it writes which measure was used on which row openly in
# the 'olcu' column.
# -------------------------------------------------------------------------
import os, sys, csv, json, time, argparse, math

VERSIYON = '1.0 (2026-08-03)'

# THE THRESHOLD COMES FROM ONE SOURCE: screening/config.py -> ESIK_DCQ = 3.0
# Its fold equivalent is 2 ** ESIK_DCQ = 8.00. NO constant is EMBEDDED; if dCq
# changes it changes in one place. The reasoning and the efficiency warning are
# written in that file.
# LOWERING THE THRESHOLD FROM INSIDE THE CODE IS STILL FORBIDDEN, see the module
# header. A change is made only through ESIK_DCQ and DELIBERATELY; on 2026-08-06
# dCq was fixed at 3.
def _esik_yukle():
    """Reads the threshold from ONE SOURCE: screening/config.py.
        Since verification/ and screening/ are sibling directories the root is derived
        from here, so the script finds it whatever working directory it is called from.

    """
    import os as _o, sys as _s
    _kok = _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))
    if _kok not in _s.path:
        _s.path.insert(0, _kok)
    from screening import config as _y
    return _y

_C = _esik_yukle()
ESIK = _C.AYRIM_ESIK
KAPSAM_ESIGI = 0.20         # a member bin counts as 'covered' at >=20% product
OKUMA_TAVANI = 3000         # the same as in the SINGLE PROTOCOL run
ENKOTU_ASGARI = 150
URUN_ALT, URUN_UST = 60, 400
KIMLIK_ESIGI = 99.0         # the consensus identity for counting as the same organism
KIL_PAYI_ALT = 5.0          # yol 3'e girecek satirin asgari mevcut kati

# --- a SEPARATE measure for universal targets ----------------------------
EVRENSEL_KAPSAMA_ESIGI = 0.80   # uye kutularin en az %80'i urun vermeli
EVRENSEL_ALANDISI_UST = 5.0     # alan disi kutularda Wilson ust sinir en cok %5

GEREKCE_EVRENSEL = u"""
YOL 1 - EVRENSEL HEDEFLERDE OLCU NEDEN DEGISTI (esik DUSURULMEDI)

Ayrim kati = (uye alt siniri) / (rakip ust siniri). Evrensel bir hedefte
"rakip" diye bir kume YOKTUR: Bakteri_universal butun bakterileri, Arke_universal
butun arkeleri cogaltmak icin tasarlandi. Payda sifira giderse oran ya 0/0
(tanimsiz, betik 0,00 yazar) ya da devasa bir sayi olur - nitekim ayni sutunda
0,00 ile 117 056 685 yan yana duruyor. Bu sayilar bir SEYI OLCMUYOR.

Bu satirlarda iddia da farklidir: "her seyi ayirt ederim" degil, "grubun
tamamini gorurum, grup disina tasmam". Dogru olcu bu iddiayi olcendir:

  KAPSAMA      = uye kutularin kaci >=%%%d urun veriyor
  ALAN DISI    = hedef grubun DISINDAKI kutularda urun veren okuma orani
                 (Wilson UST siniri - muhafazakar taraf)

GECME OLCUTU (ikisi birden):
  KAPSAMA   >= %%%d  ve  ALAN DISI <= %%%.0f

Bu, 10x esiginin gevsetilmesi DEGILDIR: 10x oranini bu satirlarda uygulamak
zaten mumkun degil, cunku oranin paydasi tanimsiz. Diger butun satirlarda 10x
esigi AYNEN durur.
""" % (int(KAPSAM_ESIGI * 100), int(EVRENSEL_KAPSAMA_ESIGI * 100), EVRENSEL_ALANDISI_UST)

# --- rows ALREADY MEASURED, not to be tried again ------------------------
BILINEN = {
    'Proteiniphilum_cinsi': dict(
        sonuc='KURTARILAMAZ',
        sebep=u'Hedef organizma numunede beyan edildigi gibi MEVCUT DEGIL. '
              u'Uye kutularin 2/3\'u olculen kimlikte Fermentimonas caenicola '
              u'(%95,33 ve %97,13, cins FARKLI) ve cift Fermentimonas\'i bilerek '
              u'disliyor (0/137). Uye kutularda urun HIC yok (0/3 kapsam). '
              u'Uyeligi daraltmak burada bir sey kurtarmaz - sorun uyelik degil '
              u'HEDEF TANIMI. Yeniden olculmedi, zaman harcanmadi.',
        yol=u'atlandi (onceden olculdu)'),
    # 2026-08-06: THIS RECORD WAS REMOVED (left here as a comment).
    # Because of a hand written 'atlandi' stamp, the row was NOT ENTERING the recovery
    # ladder at all, so route 5 (the multi locus search) never ran either.
    # The reasoning rested on 16S and it was correct, but it did not contain the
    # CONCLUSION 'it will not work at another locus either' - that had not been
    # measured. It has now: the A2 bins (a 4309 bp full operon) are the same organism
    # (98.62-99.38% identity to A1), all THREE regions of the operon were scanned, and
    # no candidate came out of any of them. The row now enters the ladder and its
    # result comes from measurement.
    #     'Methanosarcina mazei / M. soligelidi grubu': dict(
    #         sonuc='KURTARILAMAZ',
    #         sebep=u'The limiting competitor A1-4_3078083 (M. hadiensis) is A SEPARATE '
    #               u'ORGANISM: its consensus identity to the mazei bins is 98.61-98.75%, '
    #               u'while the mazei bins are 99.79-99.93% among themselves; the 16S '
    #               u'species threshold is ~98.7%. The read level probe test gave the same '
    #               u'answer (the hadiensis reads stay on their own probe, 171 to 5). So '
    #               u'narrowing the membership does not rescue this row. A recovery route '
    #               u'DOES EXIST but it needs a primer change: mazei and hadiensis differ '
    #               u'at ~19 positions in 16S and the current pair holds none of them '
    #               u'(that goes to route 3, a separate design job).',
    # 
    #     yol=u'atlandi (onceden olculdu)'),

}


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
        return ('%.*f' % (b, float(x))).replace('.', ',')
    except (TypeError, ValueError):
        return str(x)


def kutu_adi_normalize(k):
    if '_' not in k:
        return k
    bas, _, son = k.rpartition('_')
    return bas.replace('_', '-') + '_' + son


# -------------------------------------------------------------------------
# THE WILSON SCORE INTERVAL - WHY NOT THE RAW PROPORTION
# The raw proportion k/n misleads on a small sample: if 3 of 3 reads gave a product
# the raw proportion is 100% with no evidence behind it; if 0 of 200 did, the raw
# proportion is 0% while the real one could be 1.5%. The Wilson interval turns that
# uncertainty into a number, and the CONSERVATIVE side is always taken: the LOWER
# bound on the member side, the UPPER bound on the competitor side. The "outside
# the domain" proportion in route 1 is an UPPER bound for the same reason; it
# measures how much spill outside the target group there COULD BE.
# -------------------------------------------------------------------------
def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 1.0)
    p = k / float(n); d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    s = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - s), min(1.0, c + s))


# ---------------------------------------------------------------- inputs
# Stage P's panel table is the ONE entrance to this round. If the file is missing
# the run does not start: finishing silently with "there is no row to recover"
# would give the impression that the job was done when in fact nothing was tried.
def tek_protokol_oku(kok):
    """TEK_PROTOKOL_SONUC/panel_tek_protokol.tsv -> [{hedef, kaynak, karar, ...}]"""
    yol = os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'panel_tek_protokol.tsv')
    if not os.path.exists(yol):
        sys.exit(u'ERROR: %s is missing.\n      verification/full_chain.py -> option (P) has to be run first.' % yol)
    with open(yol, encoding='utf-8') as fh:
        satirlar = [s for s in fh if not s.startswith('#')]
    return list(csv.DictReader(satirlar, delimiter='\t')), yol


def _f(s):
    """'8,45' -> 8.45 ; '-' -> None"""
    s = (s or '').strip().replace(',', '.')
    try:
        return float(s)
    except ValueError:
        return None


def uyelik_dosyasi(kok):
    import glob
    # 2026-08-10: "a[-1]" DID NOT MEAN THE NEWEST. Two globs were sorted
    # alphabetically and concatenated, so engine_SONUC entries beat the ones in the
    # root regardless of date. single_protocol_measure.py carried the same trap; had
    # the two picked different files, K and P would have measured with different
    # memberships and their dCq values would not have been comparable. Both now select
    # BY TIME and find the same file.
    a = glob.glob(os.path.join(kok, 'uyelik_yeniden_turetme_uyelik_*.tsv'))
    a += glob.glob(os.path.join(kok, 'engine_SONUC', '*uyelik*.tsv'))
    if not a:
        return None
    a.sort(key=lambda p: (os.path.getmtime(p), os.path.basename(p)))
    return a[-1]


# Reads the membership table. If yeni_uye_kutular is missing it falls back to
# eski_uye_kutular, so the row is never left empty whichever version the table is.
def uyelik_oku(yol):
    out = {}
    with open(yol, encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            bol = lambda s: [kutu_adi_normalize(x.strip())
                             for x in (r.get(s) or '').split(';') if x.strip()]
            out[r['hedef'].strip()] = dict(
                sinif=(r.get('sinif') or '').strip(),
                uye=bol('yeni_uye_kutular') or bol('eski_uye_kutular'),
                karisik=bol('karisik_kutular'), rakip=bol('rakip_kutular'))
    return out


# The pair name in ek_ciftler.tsv and the target name in the membership table may
# not be the same; this table maps the two. Without the mapping the extra pairs
# would be skipped as "no membership".
def takma_adlar(kok):
    """protocol/ek_ciftler.tsv: target -> membership target (the name in the membership table)"""
    yol = os.path.join(kok, 'protocol', 'ek_ciftler.tsv')
    out = {}
    if not os.path.exists(yol):
        return out
    with open(yol, encoding='utf-8') as fh:
        for r in csv.DictReader((s for s in fh if not s.startswith('#')), delimiter='\t'):
            h = (r.get('hedef') or '').strip(); u = (r.get('uyelik_hedefi') or '').strip()
            if h and u:
                out[h] = u
    return out


# Bir hedefin evrensel (alan duzeyi) olup olmadigi. Bu satirlarda ayrim kati
# tanimsizdir; yol 1 devreye girer.
def evrensel_mi(hedef, duzey=''):
    ad = hedef.lower()
    return ('universal' in ad or 'evrensel' in ad
            or (duzey or '').strip().lower() == 'alan')


# --- REQUESTS ASKED FOR IN THE MEETING THAT NEVER ENTERED THE PANEL ----
# The chain starts from the panel, so a request with no panel row is NEVER seen by
# the chain and is left as "could not be done, never attempted". For those requests
# at least ONE design attempt is made, using the bin's OWN consensus as the
# backbone (the same method applied to the M. barkeri target).
PANELSIZ_TALEPLER = [
    dict(hedef='Podospora_pseudopauciseta (PANELSIZ TALEP)', karar='Karar 1',
         sinif='F1', uye=['F1-1_2093780', 'F1-4_2093780'],
         not_=u'Toplantida TUR ozgul istendi. Organizmanin KENDISI numunede yok '
               u'(bes referans ciftinden ucu F1 sinifinin 85 804 okumasinin '
               u'tamamina karsi tarandi, sifir urun; bolluk ust siniri %0,011). '
               u'AMA KUTU VAR ve olculen kimligi Petriella (F1-4_2093780, %99,58). '
               u'Bu deneme adlandirilmis turu degil KUTUYU hedefler.'),
    dict(hedef='Dictyostelium_discoideum_44689 (PANELSIZ TALEP)', karar='Karar 1',
         sinif='F1', uye=['F1-1_44689', 'F1-2_44689', 'F1-3_44689', 'F1-4_44689'],
         not_=u'Toplantida TUR ozgul istendi. Kraken etiketi olcumle curutuldu ve '
               u'kutuya YENI bir ad konulamadi; hedef dizisi tanimlanamadigi icin '
               u'panele hic girmedi. Bu deneme kutunun KENDI konsensusunu omurga '
               u'alir - ad bilinmese de kutuyu cogaltan bir cift bulunabilir.'),
]


# ---------------------------------------------------------------- routes
# -------------------------------------------------------------------------
# ROUTE 1 - CHANGING THE MEASURE ON UNIVERSAL TARGETS (the threshold IS NOT LOWERED)
#
# The discrimination ratio = (the member lower bound) / (the competitor upper
# bound). On a universal target there IS NO set called "competitors":
# Bakteri_universal was designed to amplify all bacteria and Arke_universal all
# archaea. As the competitor set approaches empty the denominator goes to zero and
# the ratio is either 0/0 (undefined) or a huge number. Indeed, the old panel had
# 0.00 and 117 million standing side by side in the same column. Those numbers ARE
# NOT MEASURING ANYTHING.
#
# The claim on those rows is different too: not "I can tell everything apart" but
# "I see the whole group and I do not spill outside it". The right measure is the
# one that measures that claim:
#   COVERAGE     = how many of the member bins give >=20% product
#   OUTSIDE THE DOMAIN = the proportion of reads giving a product in bins OUTSIDE
#               the target group, the Wilson UPPER bound (the conservative side)
# The passing criterion is BOTH AT ONCE: coverage >= 80% and outside <= 5%.
#
# THIS IS NOT LOWERING THE THRESHOLD. Applying the 10x ratio on these rows is not
# even possible, because its denominator is undefined. On every other row the 10x
# stands EXACTLY as it was.
#
# If there is no competitor read at all, 'outside the domain' CANNOT BE MEASURED
# and None is returned. It used to be written as 0.0, that is, the MOST FAVOURABLE
# possible value was produced and the threshold passed itself; an absence of
# measurement counted as a success.
# -------------------------------------------------------------------------
def yol1_evrensel(nm, uye, rakip, F, R):
    """Kapsama + alan disi orani. Ayrim kati KULLANILMAZ (paydasi tanimsiz)."""
    ka = na = 0
    for k in uye:
        h = nm.havuz.get(k['kutu'])
        if h is None:
            continue
        p, n, _ = h.urun_veren(F, R, URUN_ALT, URUN_UST, 1)
        na += 1
        if n and p / float(n) >= KAPSAM_ESIGI:
            ka += 1
    rp = rn = 0
    for k in rakip:
        h = nm.havuz.get(k['kutu'])
        if h is None:
            continue
        p, n, _ = h.urun_veren(F, R, URUN_ALT, URUN_UST, 1)
        rp += p; rn += n
    kapsama = (ka / float(na)) if na else 0.0
    # O-7: with no competitor read, 'outside the domain' CANNOT BE MEASURED. It used to
    # be written as 0.0, so the MOST FAVOURABLE possible value was produced and the
    # threshold passed itself.
    alandisi = (100.0 * wilson(rp, rn)[1]) if rn else None
    gecti = (kapsama >= EVRENSEL_KAPSAMA_ESIGI
             and alandisi is not None and alandisi <= EVRENSEL_ALANDISI_UST)
    return dict(kapsama=kapsama, kapsam_pay='%d/%d' % (ka, na), alandisi=alandisi,
                alandisi_pay='%d/%d' % (rp, rn), gecti=gecti)


_KOD = None


def _enc(s):
    """THE SAME encoding as engine/rederive_membership.py."""
    global _KOD
    import numpy as np
    if _KOD is None:
        m = np.full(256, 4, dtype=np.uint8)
        for i, c in enumerate('ACGT'):
            m[ord(c)] = i
        _KOD = m
    return _KOD[np.frombuffer(s.encode(), dtype=np.uint8)]


# The percent identity of two consensuses. Infix (HW) alignment: the short sequence
# is aligned INSIDE the long one and the overhang at the ends is not penalised.
# Because consensus lengths differ so much (1.5 kb beside 4.5 kb), a global
# alignment would make the same organism look different.
def hw_kimlik(a, b):
    """Align the short one as the query INSIDE the long one (infix/HW); return the percent
        identity.
        The same definition as hw_kimlik in engine/rederive_membership.py; because that
        file is a script rather than a package it was rewritten here instead of imported.

    """
    import numpy as np
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    if not q or not t:
        return 0.0
    Q, T = _enc(q), _enc(t)
    onceki = np.zeros(len(T) + 1, dtype=np.int32)   # HW: bastan bosluk bedava
    for i in range(len(Q)):
        simdi = np.empty_like(onceki)
        simdi[0] = i + 1
        esit = (T != Q[i]) | (T == 4) | (Q[i] == 4)
        # The left neighbour dependency (insertion) is VECTORISED:
        #   now[j] = min(cand[j], now[j-1]+1)
        # Setting a[j] = now[j]-j gives a[j] = min(cand[j]-j, a[j-1]), which is a running
        # minimum, so np.minimum.accumulate does it in one pass. The Python inner loop is
        # gone: a 1.5 kb x 1.5 kb alignment drops from minutes to seconds.
        aday = np.minimum(onceki[:-1] + esit, onceki[1:] + 1)
        aday = np.concatenate(([i + 1], aday))
        idx = np.arange(len(aday))
        simdi = np.minimum.accumulate(aday - idx) + idx
        onceki = simdi
    d = int(onceki.min())
    return round(100.0 * (1 - d / float(max(len(q), 1))), 2)


# -------------------------------------------------------------------------
# ROUTE 2 - NARROWING THE MEMBERSHIP
#
# The target chosen in the meeting is really A BIN; the name is that day's Kraken
# label and it may be wrong. A target's member set = the bin that was pointed at,
# plus every bin MEASUREMENT shows to be the same organism. The member bins are
# clustered by single linkage on consensus sequence identity (threshold 99%); those
# that do not belong to the target move to the competitor column and the row is
# re-measured.
#
# THE DECISION IS MADE WITHOUT LOOKING AT THE PRIMER'S RESULT. The narrowing rests
# only on measured sequence identity and is adopted UNCONDITIONALLY, even when it
# makes the row drop. The reason "it was adopted because it got past the threshold"
# appears NOWHERE in this script.
#
# Treating the largest cluster as the member set is AN ASSUMPTION and it is flagged
# as one: which cluster is the target was chosen by cluster size rather than by
# in-sample sequence evidence, and confirming it against an external reference is
# the job of stages I and G.
#
# If no two bins come out >=99% together, the narrowing IS NOT APPLIED: picking one
# would be arbitrary. That target is marked HETEROJEN and the decision is left to
# external reference confirmation.
# -------------------------------------------------------------------------
def yol2_uyelik_daralt(kons, uye_adlari, capa=None):
    """Cluster the member bins by consensus identity; the cluster holding the anchor is
        the real member set. Returns: (new_members, removed, evidence_text)

    """
    d = {k: kons[k] for k in uye_adlari if k in kons and len(kons[k]) > 200}
    if len(d) < 2:
        return (list(uye_adlari), [], u'uye kutu sayisi 2\'den az - kumeleme yapilamaz')
    adlar = sorted(d)
    kimlik = {}
    for i, a in enumerate(adlar):
        for b in adlar[i + 1:]:
            kimlik[(a, b)] = kimlik[(b, a)] = hw_kimlik(d[a], d[b])
    # single linkage clustering, threshold KIMLIK_ESIGI
    kume = {a: {a} for a in adlar}
    for (a, b), v in kimlik.items():
        if v >= KIMLIK_ESIGI and kume[a] is not kume[b]:
            yeni = kume[a] | kume[b]
            for x in yeni:
                kume[x] = yeni
    kumeler = []
    for a in adlar:
        if not any(kume[a] is k for k in kumeler):
            kumeler.append(kume[a])
    en_buyuk = max(kumeler, key=len)
    if capa and capa in kume:
        secilen = kume[capa]
    elif len(en_buyuk) >= 2:
        secilen = en_buyuk
    else:
        # no two bins came out as the same organism -> narrowing is MEANINGLESS
        ozet = '; '.join('%d kutu' % len(k) for k in kumeler)
        return (list(uye_adlari), [],
                u'DARALTMA UYGULANMADI: %d uye kutunun hicbiri birbiriyle >=%%%s '
                u'kimlikte degil (kumeler = %s). Kutular ayni organizma degil; '
                u'birini secmek keyfi olurdu. Bu hedef HETEROJEN - once kutu '
                u'kimlikleri referansla dogrulanmali.'
                % (len(adlar), vir(KIMLIK_ESIGI, 1), ozet))
    cikan = [a for a in adlar if a not in secilen]
    ic = [kimlik[(a, b)] for i, a in enumerate(sorted(secilen))
          for b in sorted(secilen)[i + 1:]]
    ozet = '; '.join('%d kutu' % len(k) for k in sorted(kumeler, key=len, reverse=True))
    ick = (u'kume ici kimlik %%%s-%%%s' % (vir(min(ic)), vir(max(ic)))) if ic else \
          u'kumede tek kutu kaldi, kume ici karsilastirma yok'
    kanit = (u'%d uye kutu KONSENSUS DIZI KIMLIGINE gore kumelendi (esik %%%s): '
             u'kumeler = %s. En buyuk kume UYE sayildi (%d kutu, %s) - bu bir '
             u'VARSAYIMDIR: hangi kumenin hedef oldugu numune ici dizi kanitiyla '
             u'degil kume buyuklugu ile secildi; dis referansla teyidi I asamasinin '
             u'isidir. Kume disinda kalan ve RAKIP hanesine tasinan: %s. '
             u'(Bu karar primerin sonucuna HIC bakmadan verildi.)'
             % (len(adlar), vir(KIMLIK_ESIGI, 1), ozet, len(secilen), ick,
                ', '.join(cikan) if cikan else 'yok'))
    return (sorted(secilen) + [a for a in uye_adlari if a not in d], cikan, kanit)


# ROUTE 4 - ROWS THAT STILL HAVE A COUNTERPART. If another pair targets the same
# member set (>=80% overlap) and DOES pass the threshold, the below-threshold row is
# unnecessary: keeping it in the panel wastes a plate position. This is not a
# recovery but a CLEAN-UP; the row is marked "moved onto the ones that dropped".
def yol4_eslenik_bul(satirlar, uyelik, alias=None):
    """Is there another pair targeting the same member set that DOES pass the threshold?"""
    alias = alias or {}
    U = lambda h: uyelik.get(alias.get(h, h), {})
    gecen = [r for r in satirlar if (r.get('esik_gecti_mi') or '').startswith('ESIK USTU')]
    out = {}
    for r in satirlar:
        if (r.get('esik_gecti_mi') or '').startswith('ESIK USTU'):
            continue
        u1 = set(U(r['hedef']).get('uye', []))
        if not u1:
            continue
        for g in gecen:
            u2 = set(U(g['hedef']).get('uye', []))
            if not u2:
                continue
            ortak = len(u1 & u2) / float(max(len(u1 | u2), 1))
            if ortak >= 0.80:
                out[r['hedef']] = (g['hedef'], _f(g.get('ASIL_ayrim_mm1')), ortak)
                break
    return out


def _ayirt_onbellekli(U, uye_diz, rak_diz):
    """ayirt_edici_mi is cached PER PRIMER.

        A BOTTLENECK FIX: cift_akisi produces N forward x M reverse pairs, but
        ayirt_edici_mi depends on THE PRIMER, so the same primer is asked about in
        hundreds of pairs. Measured: 0.030 s per call; a 60x60 grid = 3600 pairs x 2
        calls = 216 s. With the cache only 120 distinct primers are asked = 3.6 s. On a
        full run the effect is larger still.

    """
    bellek = {}

    def sor(primer, geri=False):
        anahtar = (primer, geri)
        if anahtar not in bellek:
            bellek[anahtar] = bool(U.ayirt_edici_mi(primer, uye_diz, rak_diz, geri=geri)[0])
        return bellek[anahtar]
    return sor, bellek


# -------------------------------------------------------------------------
# ROUTE 3 - REDESIGN plus ARMS
#
# For rows that fell just short (the current fold >= KIL_PAYI_ALT), a new pair is
# searched for on the backbone consensus. The candidate primers go through the
# geometry filter and are then asked, against the member and competitor consensuses,
# whether they discriminate.
#
# ARMS = a DELIBERATE mismatch added at the 2nd and 3rd base from the 3' end of the
# forward primer. This IS NOT A DEGENERATE BASE: it is one fixed sequence, it does
# not increase the oligo count, and it does not violate the decision that no
# degenerate base be used. Its purpose is to make a single base difference at the 3'
# end decisive for extension.
#
# ELIMINATION IN TWO STAGES: every candidate is first measured under THE PRIMARY
# criterion (mm<=1), and THE SECONDARY criterion (mm<=3) runs only on the 25
# candidates in front. The work halves and the deciding column is still the primary
# criterion.
#
# The "yalniz_ileri" mode KEEPS the existing reverse primer and changes only the
# forward one (the NL1 problem): in that case the reverse primer's place on the
# backbone is looked for exactly first, and with <=1 mismatch if that fails.
# -------------------------------------------------------------------------
def yol3_yeniden_tasarim(kok, nm, hedef, uye, rakip, kons, mevcut_F, mevcut_R,
                         yalniz_ileri=False, aday_ust=400, tarama_ust=3000,
                         arms_ust=5, yaz=print):
    """Yeni cift ara + ARMS varyantlari. primer3 yoksa duzgunce atlar."""
    import importlib.util
    if importlib.util.find_spec('primer3') is None:
        return dict(durum='ATLANDI', adaylar=[],
                    sebep=u'primer3-py bu makinede kurulu degil - yeniden tasarim '
                          u'taramasi yapilamadi. Kurulum: '
                          u'pip3 install primer3-py --break-system-packages')
    try:
        from screening import config as C, motor, uretec as U, geometri as G
        G.tm('ACGTACGTACGTACGTAC')
    except SystemExit:
        return dict(durum='ATLANDI', adaylar=[],
                    sebep=u'geometri modulu primer3 bulamadigi icin durdu', adaylar2=[])
    except Exception as e:
        return dict(durum='ATLANDI', adaylar=[],
                    sebep=u'yeniden tasarim baslatilamadi (%s)' % type(e).__name__)

    capa = None
    for k in uye:
        if k['kutu'] in kons and len(kons[k['kutu']]) > 500:
            capa = k['kutu']; break
    if not capa:
        return dict(durum='ATLANDI', sebep=u'kullanilabilir omurga konsensusu yok', adaylar=[])
    omurga = kons[capa]
    uye_diz = [kons[k['kutu']] for k in uye if k['kutu'] in kons]
    rak_diz = [kons[k['kutu']] for k in rakip if k['kutu'] in kons]

    yaz(u'      backbone: %s (%d bp), member consensus %d, competitor consensus %d'
        % (capa, len(omurga), len(uye_diz), len(rak_diz)))
    ad = U.aday_primerler(omurga)
    yaz(u'      candidates passing geometry: %d forward / %d reverse' % (len(ad['F']), len(ad['R'])))

    # A SPEED LIMIT: ayirt_edici_mi scans every consensus for each candidate. When the
    # candidate count runs into the tens of thousands that takes hours. We bring it down
    # to a cap by sampling at even intervals ALONG the backbone; taking the first N
    # would mean choosing every candidate from the 5' end of the backbone.
    def sey(liste, tavan):
        if len(liste) <= tavan:
            return liste
        adim = len(liste) / float(tavan)
        return [liste[int(i * adim)] for i in range(tavan)]
    once = (len(ad['F']), len(ad['R']))
    ad['F'] = sey(ad['F'], tarama_ust)
    if not yalniz_ileri:
        ad['R'] = sey(ad['R'], tarama_ust)
    if (len(ad['F']), len(ad['R'])) != once:
        yaz(u'      scan cap: %d forward / %d reverse candidates will be measured (evenly spaced along the backbone)' % (len(ad['F']), len(ad['R'])))

    if yalniz_ileri:
        # mevcut geri primeri KORU, yalniz ileriyi degistir (NL1 sorunu).
        # cift_akisi geri primerin OMURGADAKI yerini ister; rc(R)'yi omurgada ara.
        hedef_diz = motor.rc(mevcut_R)
        iR = omurga.find(hedef_diz)
        if iR < 0:                       # <=1 uyumsuzlukla ara
            L = len(hedef_diz)
            for j in range(len(omurga) - L + 1):
                if sum(1 for a, b in zip(omurga[j:j + L], hedef_diz) if a != b) <= 1:
                    iR = j; break
        if iR < 0:
            return dict(durum='ATLANDI', adaylar=[],
                        sebep=u'mevcut geri primerin omurgadaki baglanma yeri '
                              u'bulunamadi - "yalniz ileri primeri degistir" '
                              u'kipi uygulanamadi')
        ad['R'] = [(iR, len(mevcut_R), mevcut_R, G.olc(mevcut_R))]
        yaz(u'      keeping the existing reverse primer, its position on the backbone: %d' % iR)

    secilen = []
    bakilan = 0
    t0 = time.time()
    sor, bellek = _ayirt_onbellekli(U, uye_diz, rak_diz)
    BAKILAN_UST = max(20000, aday_ust * 200)   # the cap on pairs LOOKED AT, not on pairs accepted
    for t in U.cift_akisi(ad):
        bakilan += 1
        if bakilan > BAKILAN_UST:
            yaz(u'      (the cap of %d examined pairs was reached; the scan stopped)' % BAKILAN_UST)
            break
        if bakilan % 5000 == 0:
            print(u'      ... %d pairs scanned, %d discriminating (%s)          '
                  % (bakilan, len(secilen), sure_metni(time.time() - t0)), end='\r', flush=True)
        c = U.cift_yap(t)
        if not sor(c['F']):
            continue
        if not sor(c['R'], True):
            continue
        secilen.append(c)
        if len(secilen) >= aday_ust:
            break
    yaz(u'      %d pairs scanned -> %d discriminating candidates (%s, %d distinct primers queried)'
        % (bakilan, len(secilen), sure_metni(time.time() - t0), len(bellek)))

    # The ARMS variants: on the forward primer of the best few candidates, plus on the
    # CURRENT pair. (A deliberate mismatch IS NOT a degenerate base; it does not
    # increase the oligo count.)
    arms = []
    for c in secilen[:arms_ust]:
        for v, etiket in U.arms_varyantlari(c['F']):
            arms.append(dict(F=v, R=c['R'], urun=c['urun'], arms='F ' + etiket))
    for v, etiket in U.arms_varyantlari(mevcut_F):
        arms.append(dict(F=v, R=mevcut_R, urun=0, arms='F ' + etiket + ' (mevcut cift)'))

    # TWO STAGES: everything is first eliminated under THE PRIMARY criterion (mm<=1);
    # THE SECONDARY criterion (mm<=3) is measured only for the ones in front, so the
    # work halves.
    hepsi = secilen + arms
    yaz(u'      candidates to measure: %d (%d pairs + %d ARMS variants)'
        % (len(hepsi), len(secilen), len(arms)))
    ilk = []
    for j, c in enumerate(hepsi, 1):
        if j % 10 == 0:
            print(u'      ... measuring candidate %d/%d          ' % (j, len(hepsi)), end='\r', flush=True)
        o = nm.olc(c['F'], c['R'], uye, rakip, lo=URUN_ALT, hi=URUN_UST, mm=1)
        if not o or o.get('kat_enkotu') is None:
            continue
        ilk.append((o['kat_enkotu'], c, o))
    ilk.sort(key=lambda x: -x[0])
    sonuc = []
    for kat1, c, o in ilk[:25]:
        o3 = nm.olc(c['F'], c['R'], uye, rakip, lo=URUN_ALT, hi=URUN_UST, mm=3)
        sonuc.append(dict(F=c['F'], R=c['R'], urun=c.get('urun', 0),
                          arms=c.get('arms', ''), kat1=kat1,
                          kat3=(o3 or {}).get('kat_enkotu'),
                          kapsam=o.get('uye_kapsam_pay', '')))
    return dict(durum='TARANDI', sebep='', adaylar=sonuc,
                taranan=len(hepsi), omurga=capa)


# --------------------------------------------------------------- ROUTE 5
# The multi locus search: split THE WHOLE consensus into regions and try A SEPARATE
# design in each. This is what rescued Petriella; there was no answer in ITS and
# there was one in LSU. It stops the search from being confined to one backbone
# window.
#
# HOW THE REGION BOUNDARIES ARE FOUND
# -----------------------------------
# With conserved ANCHOR sequences. The anchors are regions that have been used as
# universal primers for decades and were chosen precisely because they are
# conserved; since a nanopore consensus carries some error, they are searched for
# with <=3 mismatches and with IUPAC awareness (screening/motor.find_sites). If an
# anchor is not found, that boundary is marked "not found" and the region is built
# BY THE FALLBACK ROUTE (a proportional window); the report says OPENLY which region
# was built from an anchor and which from the fallback.
CAPALAR = [
    # (ad, dizi, aciklama)  - hepsi 5'->3', sense yonunde aranir
    ('SSU_baslangic', 'AGAGTTTGATCMTGGCTCAG',  u'27F - bakteri/arke 16S basi'),
    ('SSU_orta',      'GTGYCAGCMGCCGCGGTAA',   u'515F - SSU V4 basi (16S ve 18S)'),
    ('SSU_son',       'TACGGYTACCTTGTTACGACTT', u'1492R (sense) - SSU sonu'),
    ('ITS1_baslangic','CTTGGTCATTTAGAGGAAGTAA', u'ITS1F - ITS1 basi'),
    ('58S_baslangic', 'GCATCGATGAAGAACGCAGC',  u'ITS3 - 5.8S sonu / ITS2 basi'),
    ('58S_son',       'GCTGCGTTCTTCATCGATGC',  u'ITS2 (ters tumleyen) - 5.8S basi'),
    ('LSU_baslangic', 'GCATATCAATAAGCGGAGGAAAAG', u'NL1 - LSU D1 basi'),
    ('LSU_D2_son',    'GGTCCGTGTTTCAAGACGG',   u'NL4 - D2 sonu'),
    ('LSU_ic',        'TCCTCCGCTTATTGATATGC',  u'ITS4 - LSU basi (5.8S sonrasi)'),
]


def capa_bul(motor, dizi, max_mm=3):
    """Capalari konsensuste ara. Donen: {ad: konum} (bulunanlar)."""
    out = {}
    try:
        enc = motor.encode(dizi)
    except Exception:
        return out
    for ad, d, _acik in CAPALAR:
        try:
            y = motor.find_sites(enc, d, max_mm, need_tail=False)
        except Exception:
            y = None
        if y:
            out[ad] = int(sorted(y, key=lambda x: x[1])[0][0])
    return out


def bolgeler_kur(motor, dizi, yaz):
    """Konsensusu bolgelere ayir. Donen: [(ad, bas, son, kaynak)]"""
    L = len(dizi)
    c = capa_bul(motor, dizi)
    b = []

    def ekle(ad, bas, son, kaynak):
        bas, son = max(0, int(bas)), min(L, int(son))
        if son - bas >= 200:
            b.append((ad, bas, son, kaynak))

    ssu_bas = c.get('SSU_baslangic', c.get('SSU_orta'))
    ssu_son = c.get('SSU_son')
    its1 = c.get('ITS1_baslangic')
    s58_bas = c.get('58S_son')
    s58_son = c.get('58S_baslangic')
    lsu = c.get('LSU_baslangic', c.get('LSU_ic'))
    d2son = c.get('LSU_D2_son')

    if ssu_bas is not None:
        ekle('SSU (16S/18S)', ssu_bas, ssu_son if ssu_son else (its1 or ssu_bas + 1600), 'capa')
    if its1 is not None:
        ekle('ITS1', its1, s58_bas if s58_bas else its1 + 300, 'capa')
    if s58_bas is not None and s58_son is not None:
        ekle('5.8S', s58_bas, s58_son + 20, 'capa')
    if s58_son is not None:
        ekle('ITS2', s58_son, lsu if lsu else s58_son + 350, 'capa')
    if lsu is not None:
        ekle('LSU D1-D2', lsu, d2son + 20 if d2son else lsu + 650, 'capa')
        if d2son is not None and L - d2son > 300:
            ekle('LSU kalani', d2son, L, 'capa')
    if not b:
        # THE FALLBACK ROUTE: no anchor matched - split the consensus into equal parts
        n = max(2, min(6, L // 600))
        adim = L // n
        for i in range(n):
            ekle('bolge %d/%d (YEDEK - capa bulunamadi)' % (i + 1, n),
                 i * adim, min(L, (i + 1) * adim + 100), 'yedek')
    else:
        # THE PLACES THE ANCHORS DO NOT COVER ARE SCANNED TOO. The whole point of
        # route 5 is to try ALL of the consensus; if the anchor based regions leave
        # part of the operon outside (in the archaeal A2 consensus, for instance,
        # only the SSU anchor matches and 1.6 kb of the 4.3 kb is covered), the
        # remaining pieces are added under the name "kapsanmayan".
        kapsanan = sorted((x, y) for _a, x, y, _k in b)
        bosluk, imlec = [], 0
        for x, y in kapsanan:
            if x - imlec >= 400:
                bosluk.append((imlec, x))
            imlec = max(imlec, y)
        if L - imlec >= 400:
            bosluk.append((imlec, L))
        for i, (x, y) in enumerate(bosluk, 1):
            ekle('kapsanmayan %d (capa disi)' % i, x, y, 'yedek')
        b.sort(key=lambda t: t[1])
    yaz(u'      bolgeler: %s' % '; '.join('%s %d-%d(%s)' % (a, x, y, k) for a, x, y, k in b))
    yaz(u'      bulunan capa: %s' % (', '.join(sorted(c)) or 'YOK'))
    return b, c


# ROUTE 5 - THE MULTI LOCUS SEARCH. Route 3 is confined to one backbone window;
# this route splits THE WHOLE consensus into regions and tries A SEPARATE design in
# each. Every locus has been tried before anything says "cannot be done".
def yol5_cok_lokuslu(kok, nm, hedef, uye, rakip, kons, aday_ust=150,
                     tarama_ust=800, yaz=print):
    """A SEPARATE design attempt in each region. Returns: a report, region by region."""
    import importlib.util
    if importlib.util.find_spec('primer3') is None:
        return dict(durum='ATLANDI', bolge=[],
                    sebep=u'primer3-py kurulu degil - cok lokuslu arama yapilamadi')
    from screening import engine_gateway, uretec as U

    # THE BACKBONE = the LONGEST member consensus (the 2026-08-06 fix).
    # It used to take the FIRST bin over 800 bp. When the member set held both A1 bins
    # (a 16S amplicon, ~1.4 kb) and A2 bins (the full operon, ~4.3 kb), the list order
    # put A1 first and route 5 NEVER SAW 2.9 kb of the operon: while saying "every locus
    # was tried" it had in fact scanned only 16S. Route 5 exists to split the widest
    # sequence into regions, so it takes THE LONGEST.
    capa_kutu = None
    _en = 0
    for k in uye:
        _k = k['kutu']
        if _k in kons and len(kons[_k]) > 800 and len(kons[_k]) > _en:
            capa_kutu, _en = _k, len(kons[_k])
    if not capa_kutu:
        return dict(durum='ATLANDI', bolge=[],
                    sebep=u'800 bp ustu konsensus yok - bolgelere ayrilamaz')
    omurga = kons[capa_kutu]
    yaz(u'      omurga: %s (%d bp)' % (capa_kutu, len(omurga)))
    bolge, capalar = bolgeler_kur(motor, omurga, yaz)

    uye_diz = [kons[k['kutu']] for k in uye if k['kutu'] in kons]
    rak_diz = [kons[k['kutu']] for k in rakip if k['kutu'] in kons]
    rapor = []
    for ad, bas, son, kaynak in bolge:
        alt = omurga[bas:son]
        ad_p = U.aday_primerler(alt)
        def sey(l, t):
            if len(l) <= t:
                return l
            a = len(l) / float(t)
            return [l[int(i * a)] for i in range(t)]
        ad_p['F'] = sey(ad_p['F'], tarama_ust)
        ad_p['R'] = sey(ad_p['R'], tarama_ust)
        secilen = []
        sor, _b = _ayirt_onbellekli(U, uye_diz, rak_diz)
        bakilan = 0
        BAKILAN_UST = max(8000, aday_ust * 100)
        for t in U.cift_akisi(ad_p):
            bakilan += 1
            if bakilan > BAKILAN_UST:
                break
            c = U.cift_yap(t)
            if not sor(c['F']):
                continue
            if not sor(c['R'], True):
                continue
            secilen.append(c)
            if len(secilen) >= aday_ust:
                break
        en_iyi = None
        for c in secilen:
            o = nm.olc(c['F'], c['R'], uye, rakip, lo=URUN_ALT, hi=URUN_UST, mm=1)
            if not o or o.get('kat_enkotu') is None:
                continue
            if en_iyi is None or o['kat_enkotu'] > en_iyi['kat1']:
                o3 = nm.olc(c['F'], c['R'], uye, rakip, lo=URUN_ALT, hi=URUN_UST, mm=3)
                en_iyi = dict(F=c['F'], R=c['R'], urun=c['urun'], kat1=o['kat_enkotu'],
                              kat3=(o3 or {}).get('kat_enkotu'),
                              kapsam=o.get('uye_kapsam_pay', ''))
        rapor.append(dict(bolge=ad, bas=bas, son=son, kaynak=kaynak,
                          uzunluk=son - bas, aday=len(secilen), en_iyi=en_iyi))
        yaz(u'      %-32s %4d-%4d  candidates %3d  best %s'
            % (ad, bas, son, len(secilen),
               ('%s x' % vir(en_iyi['kat1'])) if en_iyi else '-'))
    return dict(durum='TARANDI', bolge=rapor, omurga=capa_kutu,
                capalar=sorted(capalar), sebep='')


# ---------------------------------------------------------------- driver
# -------------------------------------------------------------------------
# PREPARATION AND THE DRIVER. The order: check the input -> read the panel, the
# membership and the consensuses -> select the below-threshold rows -> add the
# requests with no panel row -> build the read pools -> _tur().
#
# The read cap (--okuma) must be THE SAME as in the single protocol run: fold values
# measured at different depths cannot be compared, and a "recovered" decision could
# come out of the depth difference alone.
#
# REQUESTS WITH NO PANEL ROW: targets asked for in the meeting that never entered
# the panel. Because the chain starts from the panel these are never seen and would
# be left as "could not be done, never attempted". Here the bin's OWN consensus is
# taken as the backbone and at least one design attempt is made.
# -------------------------------------------------------------------------
def calistir(kok, aday_ust, yalniz, sifirla, tarama_ust=3000, okuma=OKUMA_TAVANI,
             arms_ust=5, panelsiz_atla=False):
    os.environ['_KURTARMA_KOK'] = kok
    sys.path.insert(0, kok)
    from screening import sample as N, hedefler as H

    CIKTI = os.path.join(kok, 'KURTARMA_SONUC')
    KONTROL = os.path.join(CIKTI, 'kontrol')
    os.makedirs(KONTROL, exist_ok=True)
    if sifirla:
        for f in os.listdir(KONTROL):
            try:
                os.remove(os.path.join(KONTROL, f))
            except OSError as e:
                print('  silinemedi: %s (%s)' % (f, e))
    gunluk = open(os.path.join(CIKTI, 'kosu_gunlugu.txt'), 'a', encoding='utf-8')

    def yaz(s=''):
        print(s, flush=True); gunluk.write(s + '\n'); gunluk.flush()

    yaz('=' * 78)
    yaz(u'  RECOVERY ROUND - rows that fell below the threshold')
    yaz(u'  version %s   %s' % (VERSIYON, time.strftime('%Y-%m-%d %H:%M')))
    yaz('=' * 78)
    yaz(u'  THE %0.0fx THRESHOLD IS NOT CHANGED. Route 1 does not lower it, it changes the' % ESIK)
    yaz(u'  MEASURE; the reasoning is at the top of the report.')
    yaz('')

    rc = girdi_denetle(yaz, 'K (verification)', [
        (os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'panel_tek_protokol.tsv'),
         'P asamasinin panel tablosu', 'P')])
    if rc:
        return rc
    satirlar, tp_yolu = tek_protokol_oku(kok)
    uy_yol = uyelik_dosyasi(kok)
    if not uy_yol:
        sys.exit(u'ERROR: uyelik_yeniden_turetme_uyelik_*.tsv is missing. Run option (U) first.')
    uyelik = uyelik_oku(uy_yol)
    kons = {d['kutu']: d['dizi'] for d in H.konsensusler()}
    kut = {k['kutu']: k for k in H.kutular()}

    hedefler = [r for r in satirlar
                if not (r.get('esik_gecti_mi') or '').startswith('ESIK USTU')]
    # panelde satiri OLMAYAN toplanti talepleri de kurtarma kapsamina alinir
    panelsiz = []
    for t in ([] if panelsiz_atla else PANELSIZ_TALEPLER):
        panelsiz.append(dict(hedef=t['hedef'], sinif=t['sinif'], F='', R='',
                             urun_bp='', ASIL_ayrim_mm1='', ASIL_kapsam_mm1='',
                             esik_gecti_mi='PANELDE SATIR YOK', _panelsiz=t))
    hedefler = hedefler + panelsiz
    if yalniz:
        hedefler = [r for r in hedefler
                    if any(y.strip().lower() in r['hedef'].lower()
                           for y in yalniz.split(','))]
    yaz(u'  single-protocol source : %s' % os.path.basename(tp_yolu))
    yaz(u'  membership source      : %s' % os.path.basename(uy_yol))
    yaz(u'  rows below threshold   : %d' % len(hedefler))
    yaz('')

    alias = takma_adlar(kok)
    eslenik = yol4_eslenik_bul(satirlar, uyelik, alias)

    # okunacak kutular
    gerekli = {}
    for r in hedefler:
        u = uyelik.get(alias.get(r['hedef'], r['hedef']))
        if not u:
            continue
        for ad in u['uye'] + u['rakip'] + u['karisik']:
            if ad in kut:
                gerekli[ad] = kut[ad]
        if not u['rakip']:
            for k in kut.values():
                if k['sinif'] == (u['sinif'] or r.get('sinif')):
                    gerekli[k['kutu']] = k
    # panelsiz taleplerin kutulari da havuza girmeli - yoksa nm.olc onlari
    # goremez ve tasarim denemesi sessizce sifir olcumle biter
    for t in PANELSIZ_TALEPLER:
        for ad in t['uye']:
            if ad in kut:
                gerekli[ad] = kut[ad]
        for k in kut.values():
            if k['sinif'] == t['sinif']:
                gerekli[k['kutu']] = k
    yaz(u'  bins to read           : %d' % len(gerekli))

    def ilerK(i, n, ad):
        print(u'   ... read pool %d/%d  %s          ' % (i, n, ad), end='\r', flush=True)

    t0 = time.time()
    yaz('')
    yaz(u'Building read pools. Only bin names scroll past on screen during this step;')
    yaz(u'it is NOT stuck. The real work starts after this and each target is saved separately.')
    nm = N.Numune(list(gerekli.values()), n=okuma, ilerle=ilerK, otorite=True)
    top = sum(h.n_okuma for h in nm.havuz.values())
    yaz(u'\nPools ready: %d bins, %d reads (%s)' % (len(gerekli), top, sure_metni(time.time() - t0)))
    tasarim_n = sum(1 for r in hedefler
                    if (_f(r.get('ASIL_ayrim_mm1')) or 0) >= KIL_PAYI_ALT
                    and r['hedef'] not in BILINEN)
    tahmin = len(hedefler) * 20 + tasarim_n * max(300, top / 60.0)
    yaz(u'ESTIMATED TIME: ~%s  (the route 3 design scan will run on %d rows)'
        % (sure_metni(tahmin), tasarim_n))
    yaz('')
    return _tur(kok, CIKTI, KONTROL, yaz, nm, hedefler, uyelik, kons, kut,
                eslenik, aday_ust, gunluk, alias, tarama_ust, arms_ust, okuma)


# -------------------------------------------------------------------------
# THE verification LOOP, TARGET BY TARGET. For each target the order is fixed:
#   a request with no panel row -> KNOWN (an outcome already measured) -> route 4 ->
#   route 1 -> route 2 -> route 3 -> route 5
#
# The rows in the KNOWN table (Proteiniphilum, the M. mazei group) ARE NOT TRIED
# AGAIN: their reasons have already been measured, they are shown on the row, and no
# time is spent.
#
# THE CHECKPOINT SEAL (_ayar): the read depth, the candidate cap, the scan cap, the
# ARMS cap and the script version are all part of the seal. Silently reusing a result
# measured at a different depth is the very mistake this chain is trying to correct.
# -------------------------------------------------------------------------
def _tur(kok, CIKTI, KONTROL, yaz, nm, hedefler, uyelik, kons, kut, eslenik,
         aday_ust, gunluk, alias=None, tarama_ust=3000, arms_ust=5,
         okuma=OKUMA_TAVANI):
    alias = alias or {}
    def kp(ad):
        t = ''.join(ch if ch.isalnum() else '_' for ch in ad)
        return os.path.join(KONTROL, 'hedef_%s.json' % t)

    AYAR = dict(okuma=okuma, aday_ust=aday_ust, tarama_ust=tarama_ust,
                arms_ust=arms_ust, surum=VERSIYON)
    # THE 2026-08-10 SEQUENCE SEAL: the primer SEQUENCE was not included in the
    # checkpoint key. When a pair's forward or reverse sequence changed, the old result
    # was silently reused. The same bug wasted two full runs at stage P (5 h 29 min plus
    # 2 h 0 min). The sequence now goes into the key; if it changes the checkpoint
    # becomes invalid.
    import hashlib as _hl

    def _ayar_of(r):
        d = dict(AYAR)
        d['dizi'] = _hl.md5(((r.get('F') or '') + '|' + (r.get('R') or ''))
                            .encode('utf-8')).hexdigest()[:12]
        return d
    sonuc = []
    tb = time.time()
    for i, r in enumerate(hedefler, 1):
        hedef = r['hedef']
        yol = kp(hedef)
        if os.path.exists(yol):
            try:
                _v = json.load(open(yol, encoding='utf-8'))
                # O-9: the settings seal - a result measured at a different depth must not be
                # reused (that is the reason the single protocol module exists).
                if _v.get('_ayar') == _ayar_of(r):
                    sonuc.append(_v)
                    yaz(u'[%2d/%2d] %-44s (taken from the previous run)'
                        % (i, len(hedefler), hedef[:44]))
                    continue
                yaz(u'[%2d/%2d] %-44s (settings changed, re-measuring)'
                    % (i, len(hedefler), hedef[:44]))
            except Exception:
                pass

        eski = _f(r.get('ASIL_ayrim_mm1'))
        s = dict(_ayar=_ayar_of(r), hedef=hedef, kaynak=r.get('kaynak', ''), eski=eski,
                 eski_kapsam=r.get('ASIL_kapsam_mm1', ''), yollar=[], yeni=None,
                 gecti='HAYIR', olcu=u'ayrim kati (%s)' % _C.esik_metni(),
                 sebep='', ayrinti={})
        yaz(u'[%2d/%2d] %s  (was %s x, coverage %s)'
            % (i, len(hedefler), hedef, vir(eski), s['eski_kapsam']))

        pz = r.get('_panelsiz')
        if pz:
            u = dict(uye=pz['uye'], rakip=[], karisik=[], sinif=pz['sinif'])
        else:
            u = uyelik.get(alias.get(hedef, hedef)) or dict(uye=[], rakip=[], karisik=[],
                                                            sinif=r.get('sinif', ''))
        coz = lambda ad: [kut[a] for a in ad if a in kut]
        uye = coz(u['uye']); rakip = coz(u['rakip']) + coz(u['karisik'])
        if not rakip:
            rakip = [k for k in kut.values()
                     if k['sinif'] == (u['sinif'] or r.get('sinif'))
                     and k['kutu'] not in set(u['uye'])]

        # --- PANELSIZ TALEP: dogrudan tasarim denemesi ---
        if pz:
            s.update(olcu=u'panelde satiri yok - kutudan tasarim denemesi',
                     yollar=["route 3 - a request with no panel row, designed from the bin's own consensus"])
            yaz(u'      -> PANELSIZ TALEP (%s): %s' % (pz['karar'], pz['not_'][:90]))
            t = yol3_yeniden_tasarim(kok, nm, hedef, uye, rakip, kons, '', '',
                                     False, aday_ust, tarama_ust, arms_ust, yaz)
            s['ayrinti']['yol3'] = t
            en_iyi = t['adaylar'][0] if t.get('adaylar') else None
            if en_iyi and en_iyi['kat1'] >= ESIK:
                s.update(yeni=u'YENI CIFT %s / %s (%d bp) %s x'
                              % (en_iyi['F'], en_iyi['R'], en_iyi['urun'], vir(en_iyi['kat1'])),
                         gecti='EVET (yeni cift)',
                         sebep=u'Panelde satiri yoktu; kutu konsensusundan tasarim '
                               u'denendi ve esigi gecen aday bulundu. %s' % pz['not_'])
                yaz(u'         BULUNDU: %s x' % vir(en_iyi['kat1']))
            else:
                yaz(u'         not available in a single window - ROUTE 5: multi-locus search')
                t5 = yol5_cok_lokuslu(kok, nm, hedef, uye, rakip, kons,
                                      aday_ust=min(aday_ust, 40),
                                      tarama_ust=min(tarama_ust, 200), yaz=yaz)
                s['yollar'].append(u'yol 5 - cok lokuslu arama (%s)' % t5['durum'])
                s['ayrinti']['yol5'] = t5
                iyi = [b for b in t5.get('bolge', [])
                       if b.get('en_iyi') and b['en_iyi']['kat1'] >= ESIK]
                if iyi:
                    en = max(iyi, key=lambda b: b['en_iyi']['kat1']); e = en['en_iyi']
                    s.update(yeni=u'YENI CIFT (%s bolgesi) %s / %s (%d bp) %s x'
                                  % (en['bolge'], e['F'], e['R'], e['urun'], vir(e['kat1'])),
                             gecti='EVET (yeni cift)',
                             sebep=u'%s Panelde satiri yoktu; %s bolgesinde cozum bulundu.'
                                   % (pz['not_'], en['bolge']))
                else:
                    s['sebep'] = (u'%s DENENDI (tek pencere + %d bolgede cok lokuslu '
                                  u'arama): esigi gecen aday yok.'
                                  % (pz['not_'], len(t5.get('bolge', []))))

        # --- bilinen, tekrar denenmeyecek ---
        elif hedef in BILINEN:
            b = BILINEN[hedef]
            s.update(gecti='HAYIR', sebep=b['sebep'], yollar=[b['yol']])
            yaz(u'      -> %s  (%s)' % (b['sonuc'], b['yol']))

        # --- YOL 4: eslenik ---
        elif hedef in eslenik:
            g, kat, ort = eslenik[hedef]
            s.update(gecti='DUSENLERE TASINDI', olcu='eslenik',
                     sebep=u'Ayni uye kumesini hedefleyen ve esigi GECEN baska bir '
                           u'cift var: "%s" (%s x, uye kumesi ortusme %%%d). Bu satir '
                           u'artik gereksiz; panelde tutulmasi plaka yeri israfidir.'
                           % (g, vir(kat), int(100 * ort)),
                     yollar=['route 4 - a row whose counterpart is still there'])
            yaz(u'      -> ROUTE 4: an equivalent exists (%s, %s x), moved to the failed list' % (g, vir(kat)))

        # --- YOL 1: evrensel ---
        elif evrensel_mi(hedef, r.get('duzey', '')):
            o = yol1_evrensel(nm, uye, rakip, r['F'], r['R'])
            s.update(olcu=u'KAPSAMA + ALAN DISI (ayrim kati bu satirda tanimsiz)',
                     yollar=[u'yol 1 - olcu duzeltildi'],
                     yeni=u'kapsama %s (%%%d), alan disi %%%s'
                          % (o['kapsam_pay'], int(100 * o['kapsama']), vir(o['alandisi'])),
                     gecti='EVET' if o['gecti'] else 'HAYIR',
                     ayrinti=o)
            if not o['gecti']:
                s['sebep'] = (u'Kapsama %%%d (olcut %%%d) / alan disi %%%s (olcut en cok %%%.0f).'
                              % (int(100 * o['kapsama']), int(100 * EVRENSEL_KAPSAMA_ESIGI),
                                 vir(o['alandisi']), EVRENSEL_ALANDISI_UST))
            yaz(u'      -> ROUTE 1: coverage %s, outside the domain %%%s  => %s'
                % (o['kapsam_pay'], vir(o['alandisi']), s['gecti']))

        else:
            # --- YOL 2: uyelik daraltma ---
            kapsam_tam = (r.get('ASIL_kapsam_mm1') or '').split('/')
            tam = (len(kapsam_tam) == 2 and kapsam_tam[0] == kapsam_tam[1])
            # ---------------------------------------------------------------
            # THE MEMBERSHIP IS ADOPTED UNCONDITIONALLY - AND A DROP IS NOT A LOSS.
            # After narrowing, the measured value can come out LOWER than before. That is
            # not the primer getting worse but THE MEASURE being corrected: the old value
            # came from a wrong membership (bins not belonging to the target had been
            # counted as members, or bins that are the same organism had been left in the
            # competitor column) and was never valid. Letting the row drop is the proof
            # that the rule does not work in one direction only.
            # ---------------------------------------------------------------
            # --- ROUTE 2: NARROWING THE MEMBERSHIP ---
            # THE CRITICAL FIX (design review item 1): the membership is NO LONGER adopted
            # according to the primer's result. The narrowing is done only by MEASURED
            # SEQUENCE IDENTITY and is adopted UNCONDITIONALLY, even when the result makes
            # the target drop. The reason "it was adopted because it got past the
            # threshold" appears NOWHERE.
            if len(u['uye']) > 1:
                yeni_uye, cikan, kanit = yol2_uyelik_daralt(kons, u['uye'])
                s['ayrinti']['yol2'] = dict(kanit=kanit, cikan=cikan,
                                            uygulandi=bool(cikan))
                if cikan:
                    uye2 = coz(yeni_uye)
                    rakip2 = rakip + coz(cikan)
                    o = nm.olc(r['F'], r['R'], uye2, rakip2, lo=URUN_ALT, hi=URUN_UST, mm=1)
                    o3 = nm.olc(r['F'], r['R'], uye2, rakip2, lo=URUN_ALT, hi=URUN_UST, mm=3)
                    kat = (o or {}).get('kat_enkotu')
                    kat3 = (o3 or {}).get('kat_enkotu')
                    # KOSULSUZ BENIMSEME - yon ne olursa olsun
                    uye, rakip = uye2, rakip2
                    s['yollar'].append(u'yol 2 - uyelik daraltildi (KOSULSUZ benimsendi)')
                    s['ayrinti']['yol2'].update(kat1=kat, kat3=kat3,
                                                kapsam=(o or {}).get('uye_kapsam_pay'))
                    dus = (kat is not None and eski is not None and kat < eski)
                    yon = (u'DEGISMEDI' if (kat is None or eski is None) else
                           (u'YUKSELDI' if kat > eski else
                            u'DUSTU' if kat < eski else u'DEGISMEDI'))
                    yaz(u'      -> YOL 2: %s' % kanit)
                    yaz(u'         with the narrowed membership: %s x (previously %s x, %s), adopted UNCONDITIONALLY' % (vir(kat), vir(eski), yon))
                    if dus:
                        s['dusus_notu'] = (
                            u'DIKKAT - BU BIR KAYIP DEGIL, DUZELTMEDIR. Bu satirin '
                            u'eski %s x degeri YANLIS UYELIKTEN geliyordu: hedefe ait '
                            u'olmayan kutular uye sayilmis, ya da ayni organizma olan '
                            u'kutular rakip hanesinde birakilmisti. Uyelik olculen dizi '
                            u'kimligine gore duzeltilince gercek deger %s x cikti. '
                            u'Dusus, primerin kotulesmesi degil OLCUNUN duzelmesidir; '
                            u'eski deger hicbir zaman gecerli degildi. Uyelik karari '
                            u'primerin sonucuna BAKILMADAN verildi - bu satirin dusmesi '
                            u'kuralın tek yonlu calismadiginin kanitidir.'
                            % (vir(eski), vir(kat)))
                        yaz(u'         NOTE: the drop is NOT A LOSS, it is the measurement being corrected. The old %s x came from a wrong membership.' % vir(eski))
                    s['uyelik_gerekcesi'] = (
                        u'Uye kumesi YALNIZ olculen konsensus kimligine gore '
                        u'belirlendi (esik %%%s), primerin sonucundan BAGIMSIZ olarak '
                        u've kosulsuz benimsendi. Kanit: %s. Yeni deger %s x '
                        u'(eski %s x, %s) - bu deger benimseme kararini ETKILEMEDI.'
                        % (vir(KIMLIK_ESIGI, 1), kanit, vir(kat), vir(eski), yon))
                    s['eski'] = eski = kat      # bundan sonrasi duzeltilmis uyelikle
                    if kat is not None and kat >= ESIK:
                        s.update(yeni='%s x' % vir(kat), gecti='EVET',
                                 sebep=s['uyelik_gerekcesi'])
                    else:
                        s['sebep'] = s['uyelik_gerekcesi']
                else:
                    s['yollar'].append(u'yol 2 - daraltma uygulanmadi')
                    s['sebep'] = kanit
                    s['uyelik_gerekcesi'] = kanit
                    yaz(u'      -> YOL 2: %s' % kanit)

            # --- ROUTE 3: redesign plus ARMS ---
            # The KIL_PAYI_ALT gate is right FOR ROUTE 3: moving a primer inside the same
            # backbone window does not take 0.8x to 8x. It was WRONG for ROUTE 5, see the
            # note below.
            if s['gecti'] != 'EVET' and (eski or 0) >= KIL_PAYI_ALT:
                yalniz_ileri = ('microasca' in hedef.lower())
                yaz(u'      -> ROUTE 3: redesign%s' %
                    (u' (only the FORWARD primer will be changed - NL1)' if yalniz_ileri else ''))
                t = yol3_yeniden_tasarim(kok, nm, hedef, uye, rakip, kons,
                                         r['F'], r['R'], yalniz_ileri, aday_ust,
                                         tarama_ust, arms_ust, yaz)
                s['yollar'].append(u'yol 3 - yeniden tasarim + ARMS (%s)' % t['durum'])
                s['ayrinti']['yol3'] = t
                en_iyi = t['adaylar'][0] if t.get('adaylar') else None
                if en_iyi and en_iyi['kat1'] >= ESIK:
                    s.update(yeni=u'YENI CIFT %s / %s (%d bp) %s x%s'
                                  % (en_iyi['F'], en_iyi['R'], en_iyi['urun'],
                                     vir(en_iyi['kat1']),
                                     (u' [ARMS: %s]' % en_iyi['arms']) if en_iyi['arms'] else ''),
                             gecti='EVET (yeni cift)',
                             sebep=u'Mevcut cift esigi gecmiyor; taramada esigi gecen aday bulundu.')
                    yaz(u'         BULUNDU: %s x  %s / %s'
                        % (vir(en_iyi['kat1']), en_iyi['F'], en_iyi['R']))
                elif t['durum'] == 'ATLANDI':
                    s['sebep'] = s['sebep'] or t['sebep']
                    yaz(u'         %s' % t['sebep'])
                else:
                    s['sebep'] = s['sebep'] or (
                        u'Tarandi, esigi gecen aday yok (en iyi %s x).'
                        % (vir(en_iyi['kat1']) if en_iyi else '-'))
                    yaz(u'         no candidate passes the threshold')

                # --- ROUTE 5: THE MULTI LOCUS SEARCH ---
                # Route 3 rests on a single backbone window. For Petriella the answer was not
                # in ITS but in LSU; EVERY locus is tried before anything says "cannot be done".
                yaz(u'      -> YOL 5: cok lokuslu arama (bolge bolge)')
                t5 = yol5_cok_lokuslu(kok, nm, hedef, uye, rakip, kons,
                                      aday_ust=min(aday_ust, 40),
                                      tarama_ust=min(tarama_ust, 200), yaz=yaz)
                s['yollar'].append(u'yol 5 - cok lokuslu arama (%s)' % t5['durum'])
                s['ayrinti']['yol5'] = t5
                iyi = [b for b in t5.get('bolge', [])
                       if b.get('en_iyi') and b['en_iyi']['kat1'] >= ESIK]
                if iyi:
                    en = max(iyi, key=lambda b: b['en_iyi']['kat1'])
                    e = en['en_iyi']
                    s.update(yeni=u'YENI CIFT (%s bolgesi) %s / %s (%d bp) %s x'
                                  % (en['bolge'], e['F'], e['R'], e['urun'], vir(e['kat1'])),
                             gecti='EVET (yeni cift)',
                             sebep=u'Tek omurga penceresinde cozum yoktu; %s bolgesinde '
                                   u'bulundu. Taranan bolge: %d.'
                                   % (en['bolge'], len(t5.get('bolge', []))))
                    yaz(u'         FOUND: the %s region, %s x' % (en['bolge'], vir(e['kat1'])))
                elif t5.get('bolge'):
                    s['sebep'] += (u' COK LOKUSLU ARAMA: %d bolgenin hepsi tarandi '
                                   u'(%s), hicbirinde esigi gecen aday yok.'
                                   % (len(t5['bolge']),
                                      ', '.join(b['bolge'] for b in t5['bolge'])))
            elif s['gecti'] != 'EVET':
                # ---------------------------------------------------------------
                # THE 2026-08-06 LOGIC BUG FIX.
                # NOTHING AT ALL used to be tried in this branch: if the floor was below
                # KIL_PAYI_ALT the row was left with the note "redesign did not run".
                # Right FOR ROUTE 3, WRONG for ROUTE 5:
                #   - route 3 moves a primer inside THE SAME backbone window; if the current
                #     fold is 0.8x there is no answer in that window and the gate belongs.
                #   - route 5 moves to ANOTHER LOCUS. The current fold is the measure of the
                #     locus BEING ABANDONED; it says nothing about the new one.
                # The result: on five targets with a low floor the multi locus search NEVER RAN
                # and the report said "the discrimination was not enough", when what we needed
                # to be able to say was "it was tried at every locus". The gate now applies to
                # route 3 only.
                # ---------------------------------------------------------------
                _uzun = max([len(kons[k['kutu']]) for k in uye if k['kutu'] in kons] or [0])
                if _uzun >= 2000:
                    yaz(u'      -> ROUTE 5: the floor is low but the consensus is %d bp, so OTHER LOCI are tried' % _uzun)
                    t5 = yol5_cok_lokuslu(kok, nm, hedef, uye, rakip, kons,
                                          aday_ust=min(aday_ust, 40),
                                          tarama_ust=min(tarama_ust, 200), yaz=yaz)
                    s['yollar'].append(u'yol 5 - cok lokuslu arama (%s)' % t5['durum'])
                    s['ayrinti']['yol5'] = t5
                    iyi = [b for b in t5.get('bolge', [])
                           if b.get('en_iyi') and b['en_iyi']['kat1'] >= ESIK]
                    if iyi:
                        en = max(iyi, key=lambda b: b['en_iyi']['kat1']); e = en['en_iyi']
                        s.update(yeni=u'YENI CIFT (%s bolgesi) %s / %s (%d bp) %s x'
                                      % (en['bolge'], e['F'], e['R'], e['urun'], vir(e['kat1'])),
                                 gecti='EVET (yeni cift)',
                                 sebep=u'Mevcut lokusta taban cok dusuktu; %s bolgesinde '
                                       u'esigi gecen aday bulundu. Taranan bolge: %d.'
                                       % (en['bolge'], len(t5.get('bolge', []))))
                        yaz(u'         FOUND: the %s region, %s x' % (en['bolge'], vir(e['kat1'])))
                    else:
                        s['sebep'] = ((s['sebep'] + u'  ') if s['sebep'] else u'') + (
                            u'BUTUN LOKUSLAR DENENDI: %d bolge tarandi (%s), hicbirinde '
                            u'esigi gecen aday yok.'
                            % (len(t5.get('bolge', [])),
                               ', '.join(b['bolge'] for b in t5.get('bolge', [])) or '-'))
                elif not s['sebep']:
                    s['sebep'] = (u'Taban cok dusuk (%s x < %s x) ve konsensus %d bp - '
                                  u'tek lokus (16S) var, gidilecek baska bolge YOK.'
                                  % (vir(eski), vir(KIL_PAYI_ALT), _uzun))

        json.dump(s, open(yol, 'w', encoding='utf-8'), ensure_ascii=False, indent=1, default=str)
        sonuc.append(s)
        g = time.time() - tb
        print(u'        elapsed %s | estimated remaining %s'
              % (sure_metni(g), sure_metni(g / i * (len(hedefler) - i))), flush=True)

    raporla(CIKTI, sonuc, yaz)
    rc = cikti_denetle(yaz, 'K (verification)', [
        (os.path.join(CIKTI, 'kurtarma_satirlari.tsv'), 'kurtarma_satirlari.tsv')])
    gunluk.close()
    return rc


# -------------------------------------------------------------------------
# Three outputs: a one row per target table, a candidate primer table and a markdown
# report. The 'olcu' column says on every row WHICH measure was applied: coverage
# plus outside-the-domain on the universal rows, the discrimination ratio (10x) on
# the rest. That the threshold was not lowered, and why, is written in the report
# header.
# -------------------------------------------------------------------------
def raporla(CIKTI, sonuc, yaz):
    yol = os.path.join(CIKTI, 'kurtarma_satirlari.tsv')
    with open(yol, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# RECOVERY ROUND - ONE row per target.\n')
        fh.write(u'# MEMBERSHIP RULE: the member set is determined ONLY by measured sequence\n')
        fh.write(u'# identity, and is adopted UNCONDITIONALLY. Whether the primer passed the\n')
        fh.write(u'# threshold does NOT affect the membership decision. A reason such as\n')
        fh.write(u'# "adopted because it passed the threshold" does NOT exist in this file.\n')
        fh.write(u'# THE %0.0fx THRESHOLD WAS NOT CHANGED. The "olcu" column says which measure was applied;\n' % ESIK)
        fh.write(u'# on universal rows the discrimination ratio is UNDEFINED, so COVERAGE + OUT-OF-DOMAIN is used instead.\n')
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'eski_deger', 'eski_kapsam', 'denenen_yol', 'olcu',
                    'yeni_deger', 'esigi_gecti_mi', 'UYELIK_GEREKCESI',
                    'DUSUS_KAYIP_MI_DUZELTME_MI', 'sebep'])
        for s in sonuc:
            w.writerow([s['hedef'], vir(s['eski']), s['eski_kapsam'],
                        ' + '.join(s['yollar']) or '-', s['olcu'],
                        s['yeni'] or '-', s['gecti'],
                        s.get('uyelik_gerekcesi', '-'), s['sebep']])
    yaz(u'  written: %s' % yol)

    ay = os.path.join(CIKTI, 'yeni_adaylar.tsv')
    with open(ay, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# Candidates produced by the route 3 scan (the best 25 for each target).\n')
        fh.write(u'# ARMS = a DELIBERATE mismatch at the 2nd and 3rd base from the 3\' end. It is NOT a\n')
        fh.write(u'# degenerate base, does not increase the oligo count and does not break the agreed rules, but it is a separate item.\n')
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'bolge', 'F', 'R', 'urun_bp', 'arms',
                    'ayrim_mm1', 'ayrim_mm3', 'kapsam'])
        for s in sonuc:
            for a in (s.get('ayrinti', {}).get('yol3', {}) or {}).get('adaylar', []):
                w.writerow([s['hedef'], 'tek pencere (yol 3)', a['F'], a['R'],
                            a['urun'], a['arms'], vir(a['kat1']), vir(a['kat3']), a['kapsam']])
            for b in (s.get('ayrinti', {}).get('yol5', {}) or {}).get('bolge', []):
                e = b.get('en_iyi')
                if e:
                    w.writerow([s['hedef'], '%s (%d-%d, %s)' % (b['bolge'], b['bas'],
                                                               b['son'], b['kaynak']),
                                e['F'], e['R'], e['urun'], '', vir(e['kat1']),
                                vir(e['kat3']), e['kapsam']])
                else:
                    w.writerow([s['hedef'], '%s (%d-%d, %s)' % (b['bolge'], b['bas'],
                                                               b['son'], b['kaynak']),
                                '', '', '', '', 'aday yok (%d taranan)' % b['aday'], '', ''])
    yaz(u'  written: %s' % ay)

    gecen = [s for s in sonuc if s['gecti'].startswith('EVET')]
    tasi = [s for s in sonuc if s['gecti'] == 'DUSENLERE TASINDI']
    kalan = [s for s in sonuc if s not in gecen and s not in tasi]
    rp = os.path.join(CIKTI, 'KURTARMA_RAPORU.md')
    with open(rp, 'w', encoding='utf-8') as fh:
        fh.write(u'# Recovery round\n\nGenerated: %s, script %s\n\n'
                 % (time.strftime('%Y-%m-%d %H:%M'), VERSIYON))
        fh.write(u'## Result\n\n- Recovered: **%d**\n- Moved to failed (an equivalent exists): **%d**\n- Not recoverable: **%d**\n\n' % (len(gecen), len(tasi), len(kalan)))
        fh.write(u'> **The threshold was not lowered.** On universal rows route 1 replaces the discrimination ratio with a coverage plus out-of-domain measure. That is not a relaxation of the threshold: on those rows the denominator of the ratio is undefined. On every other row 10x was applied exactly as before.\n\n')
        fh.write(u'```' + GEREKCE_EVRENSEL + u'```\n\n')
        fh.write(u'## Row by row\n\n')
        fh.write(u'| target | before | route | after | passed |\n|---|---|---|---|---|\n')
        for s in sonuc:
            fh.write(u'| %s | %s | %s | %s | %s |\n'
                     % (s['hedef'], vir(s['eski']), ' + '.join(s['yollar']) or '-',
                        (s['yeni'] or '-')[:70], s['gecti']))
        fh.write(u'\n## Kurtarilamayanlarin sebebi\n\n')
        for s in kalan:
            fh.write(u'**%s** — %s\n\n' % (s['hedef'], s['sebep'] or '-'))
        fh.write(u'\n## Reading order\n\n1. This file. 2. `kurtarma_satirlari.tsv` (one row per target). 3. `yeni_adaylar.tsv` ')
    yaz(u'  written: %s' % rp)
    yaz('')
    yaz(u'  RECOVERED: %d   MOVED TO FAILED: %d   NOT RECOVERABLE: %d'
        % (len(gecen), len(tasi), len(kalan)))



# --------------------------------------------------------------- guvenlik agi
def cikti_denetle(yaz, ad, dosyalar, asgari=1):
    """When the stage ends, it audits ITS OWN output.

        If the expected row count is zero, or the file is missing entirely, it DOES NOT
        CARRY ON SILENTLY: it prints a clear error and returns a non-zero code. This is
        so that it cannot produce an empty result overnight and then say "nothing was
        found" in the morning.

    """
    sorun = []
    for yol, etiket in dosyalar:
        if not os.path.exists(yol):
            sorun.append(u'%s WAS NOT PRODUCED (%s)' % (etiket, yol)); continue
        try:
            with open(yol, encoding='utf-8') as fh:
                n = sum(1 for s in fh if s.strip() and not s.startswith('#'))
            n = max(0, n - 1)          # baslik satiri
        except OSError as e:
            sorun.append(u'%s COULD NOT BE READ (%s)' % (etiket, e)); continue
        if n < asgari:
            sorun.append(u'%s IS EMPTY, %d data rows (at least %d were expected)'
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
            eksik.append(u'there is no %s (%s) -> run stage %s first' % (etiket, yol, uretici))
            continue
        with open(yol, encoding='utf-8') as fh:
            n = sum(1 for s in fh if s.strip() and not s.startswith('#'))
        if n <= 1:
            eksik.append(u'%s IS EMPTY (%s) -> stage %s produced no result'
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

# The command line: --aday-ust and --tarama-ust set the size of the route 3 scan,
# --arms-ust how many candidates get an ARMS variant, --okuma the depth cap (it must
# match P), and --panelsiz-atla skips the requests with no panel row in a quick test.
def main():
    p = argparse.ArgumentParser(description='Esik alti satirlar icin kurtarma turu')
    p.add_argument('--root', '--kok', dest='kok', default='.')
    p.add_argument('--candidate-max', '--aday-ust', dest='aday_ust', type=int, default=400,
                   help='yol 3 taramasinda olculecek en fazla aday cift')
    p.add_argument('--scan-max', '--tarama-ust', dest='tarama_ust', type=int, default=3000,
                   help='yol 3 taramasinda omurgadan examplele nen en fazla primer adayi')
    p.add_argument('--only', '--yalniz', dest='yalniz', default=None, help='only targets whose name contains this (test)')
    p.add_argument('--arms-max', '--arms-ust', dest='arms_ust', type=int, default=5,
                   help='kac adayin ARMS varyantlari uretilsin')
    p.add_argument('--reads', '--okuma', dest='okuma', type=int, default=OKUMA_TAVANI,
                   help='cap on reads per bin (must match the single-protocol measurement)')
    p.add_argument('--skip-if-no-panel', '--panelsiz-atla', dest='panelsiz_atla', action='store_true',
                   help='skip requests with no row in the panel (quick testing only)')
    p.add_argument('--reset', '--sifirla', dest='sifirla', action='store_true')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)
    if not os.path.isdir(os.path.join(kok, 'screening')):
        sys.exit(u'ERROR: there is no screening directory inside %s. Give the project directory with --kok.' % kok)
    return calistir(kok, a.aday_ust, a.yalniz, a.sifirla, a.tarama_ust, a.okuma,
                    a.arms_ust, a.panelsiz_atla)


if __name__ == '__main__':
    sys.exit(main() or 0)
