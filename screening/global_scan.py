# -*- coding: utf-8 -*-
"""Global specificity scan, the most expensive step, run last.

WHAT IT DOES
    Scans a full reference database (SILVA SSURef NR99, ~500k records; UNITE or
    LSU when needed) for places where both primers bind in opposing orientation
    within the product window. The criterion matches the panel's global rule:
    at most 5 mismatches in total, F and R reported separately.

    This is a RAW SCANNER. It counts every product it finds; deciding which of
    those are genuinely off-target is the caller's job (see
    verification/specificity_round.py). Keeping the two apart is deliberate, a
    scanner that also judged would make its own bugs invisible.

EFFICIENCY
    The database is read ONCE and all candidates are measured against each
    chunk together. A per-candidate pass would re-read ~500k records once per
    candidate. Chunk size is bounded because peak memory is roughly six times
    the chunk.

RESUMABILITY
    Every chunk writes its state to disk; restarting continues from the last
    completed chunk.

OPTIONAL TAXONOMIC CLASSIFICATION
    Pass `siniflandirici` to have every hit classified into the four D-12
    classes (inside clade / organelle / same domain outside clade / different
    domain) and COUNTED. Identities are not stored: one universal primer alone
    produced 483,098 hits, and keeping headers for those would cost ~100 MB.
    Counters are constant memory per candidate and the count is complete; the
    `vurus` list stays capped at 300 as an evidence sample.

    'siniflandirildi' False means the classification was NOT RUN. It does not
    mean "no cross-reaction found". Callers must not conflate the two.

"""
# -------------------------------------------------------------------------
# global_scan.py scans the last surviving candidates against the full database in
#                     REFERENCE_DB; it is the most expensive step of the search.
#
# INPUT  : the candidate list ({'ad','F','R','lo','hi'}); as the database,
#          yapilandirma.SILVA_SSU (or another fasta the caller gives); the
#          kontrol/kuresel_<target>.pkl intermediate state if there is one. The scan
#          is done with engine_gateway.encode and engine_gateway.find_sites (a numpy
#          vector search).
# OUTPUT : if durum_yolu is given it writes a pickled intermediate state at the end
#          of every chunk. tara() returns a
#          {candidate_name: {'urun': n, 'boy': {...}, 'vurus': [...]}} dictionary.
# CALLED BY: stage E inside __main__.hedefi_isle, that is verification/full_chain.py
#          keys 1, 2, 3, 7 and the 7th stage of 9 (skipped if --light is given).
#          From outside, verification/specificity_round.py (key D) uses this module
#          AS IT IS; it is the local database layer of the verification round.
#
# THE EFFICIENCY REASON: the database is read ONCE and ALL the candidates are
# measured together on each chunk. Had a separate pass been made per candidate, the
# roughly 500 thousand records would have been scanned once per candidate. The chunk
# size is bounded by KURESEL_PARCA because the memory ceiling is about six times
# that size.
# -------------------------------------------------------------------------
import os, json, pickle, time
import numpy as np
from . import config as C
from . import engine_gateway

AYIRAC = 'N' * 60


def _parcalar(db, parca_baz=C.KURESEL_PARCA):
    """Split the FASTA into blocks of about the given size; each block is
        (a header list, a length list, the concatenated sequence)."""
    ad_l, uz_l, buf, tot = [], [], [], 0
    for ad, seq in engine_gateway.read_fasta(db):
        s = engine_gateway.clean(seq.upper())
        # A2 (2026-08-21): the header IS NO LONGER CUT (it used to be ad[:150]).
        # Measured: 16.6 percent of SILVA SSU headers go over 150 characters and
        # the tail that gets cut is exactly the GENUS and SPECIES tokens, that is,
        # the ones most likely to match the target clade. Classifying with a cut
        # header makes a record INSIDE the clade count as OUTSIDE it; it produces
        # the very fault we are trying to correct.
        # The cost was measured and is negligible: a 40 MB sequence chunk is about
        # 28,500 records at an average header of 134 characters, so roughly 3.8 MB
        # per chunk. ad_l is local to the chunk and IS NOT WRITTEN to the pickle.
        ad_l.append(ad); uz_l.append(len(s)); buf.append(s)
        tot += len(s) + len(AYIRAC)
        if tot >= parca_baz:
            yield ad_l, uz_l, AYIRAC.join(buf)
            ad_l, uz_l, buf, tot = [], [], [], 0
    if buf:
        yield ad_l, uz_l, AYIRAC.join(buf)


def _kayit_indeksi(uz_l):
    'The start of every record inside the concatenated sequence, separator included.'
    off = np.zeros(len(uz_l), dtype=np.int64)
    c = 0
    for i, u in enumerate(uz_l):
        off[i] = c
        c += u + len(AYIRAC)
    return off


# The checkpoint format version. A2 added the 'sinif' counters; old pickles that DO
# NOT CARRY that field ARE INVALID. Without the version check an old checkpoint
# would be read back silently and the taxonomic counters would stay zero, and
# "measured, no cross reaction found" could not be told apart from "never
# measured".
DURUM_SURUMU = 2

SINIFLAR = ('a', 'ao', 'b', 'c', 'bilinmiyor')


def _bos_sonuc():
    return dict(urun=0, boy={}, vurus=[],
                sinif={k: 0 for k in SINIFLAR}, siniflandirildi=False)


def tara(adaylar, db=None, durum_yolu=None, ilerle=None, max_mm=C.KURESEL_MAX_MM,
         siniflandirici=None):
    """adaylar: [{'ad':..,'F':..,'R':..,'lo':..,'hi':..}]  (there should be few of them)

    siniflandirici : None, or  f(candidate_name, header, db_file_name) -> a class string
        If it is given, every hit is put into one of D-12's classes and COUNTED.
        A2 (2026-08-21): COUNTING rather than STORING the identities is a deliberate
        decision. Measured: Bakteri_universal alone gives 483,098 hits; keeping a
        header for each would come to about 100 MB, and the 'vurus' list is capped at
        300 anyway. A counter, on the other hand, holds CONSTANT memory per candidate
        and the count is COMPLETE. The 'vurus' list stays at 300 as an evidence
        sample.

    Returns: {candidate_name: {'urun':n, 'boy':{}, 'vurus':[(header,length,mmF,mmR)],
                      'sinif':{'a':n,'ao':n,'b':n,'c':n,'bilinmiyor':n},
                      'siniflandirildi':bool}}
    If 'siniflandirildi' is False the class counters are ZERO, but that DOES NOT MEAN
    "no cross reaction"; it means no measurement was made. The caller must not confuse
    the two.

    """
    db = db or C.SILVA_SSU
    if not os.path.exists(db):
        return {a['ad']: dict(hata='there is no such database: %s' % db) for a in adaylar}
    db_ad = os.path.basename(db)

    durum = dict(surum=DURUM_SURUMU, parca=0, toplam_kayit=0,
                 res={a['ad']: _bos_sonuc() for a in adaylar})
    if durum_yolu and os.path.exists(durum_yolu):
        try:
            eski = pickle.load(open(durum_yolu, 'rb'))
            if eski.get('surum') == DURUM_SURUMU:
                durum = eski
            # if the version does not match it is scanned FROM SCRATCH; the old result is NOT returned silently.
        except Exception:
            pass

    baslangic = durum['parca']
    pi = 0
    t0 = time.time()
    for ad_l, uz_l, big in _parcalar(db):
        pi += 1
        if pi <= baslangic:
            continue
        enc = engine_gateway.encode(big)
        off = _kayit_indeksi(uz_l)
        durum['toplam_kayit'] += len(ad_l)
        for a in adaylar:
            F, R = a['F'], a['R']
            revrc = engine_gateway.rc(R)
            fs = engine_gateway.find_sites(enc, F, max_mm, need_tail=False)
            if not fs:
                continue
            rs = engine_gateway.find_sites(enc, revrc, max_mm, need_tail=False)
            if not rs:
                continue
            fpos = np.array([x[0] for x in fs]); fmm = np.array([x[1] for x in fs])
            rpos = np.array([x[0] for x in rs]); rmm = np.array([x[1] for x in rs])
            frec = np.searchsorted(off, fpos, 'right') - 1
            rrec = np.searchsorted(off, rpos, 'right') - 1
            # kayit ici konum
            fin = fpos - off[frec]; rin = rpos - off[rrec]
            gecerli_f = fin + len(F) <= np.array(uz_l)[frec]
            gecerli_r = rin + len(revrc) <= np.array(uz_l)[rrec]
            fpos, fmm, frec, fin = fpos[gecerli_f], fmm[gecerli_f], frec[gecerli_f], fin[gecerli_f]
            rpos, rmm, rrec, rin = rpos[gecerli_r], rmm[gecerli_r], rrec[gecerli_r], rin[gecerli_r]
            ort = set(frec.tolist()) & set(rrec.tolist())
            if not ort:
                continue
            r = durum['res'][a['ad']]
            for kid in sorted(ort):
                fi = fin[frec == kid]; fm = fmm[frec == kid]
                ri = rin[rrec == kid]; rm = rmm[rrec == kid]
                en = None
                for x, xm in zip(fi, fm):
                    for yy, ym in zip(ri, rm):
                        bp = int(yy + len(revrc) - x)
                        if a['lo'] <= bp <= a['hi'] and yy >= x + len(F) and xm + ym <= max_mm:
                            if en is None or xm + ym < en[1] + en[2]:
                                en = (bp, int(xm), int(ym))
                if en:
                    r['urun'] += 1
                    r['boy'][en[0]] = r['boy'].get(en[0], 0) + 1
                    # A2: SAYAC tavansiz, KIMLIK listesi tavanli.
                    # Taksonomik hukum sayaclardan uretilir; 'vurus' yalnizca
                    # insana gosterilecek kanit ornegidir.
                    if siniflandirici is not None:
                        try:
                            s = siniflandirici(a['ad'], ad_l[kid], db_ad)
                        except Exception:
                            s = 'bilinmiyor'
                        if s not in r['sinif']:
                            s = 'bilinmiyor'
                        r['sinif'][s] += 1
                        r['siniflandirildi'] = True
                    if len(r['vurus']) < 300:
                        r['vurus'].append((ad_l[kid], en[0], en[1], en[2]))
        durum['parca'] = pi
        if durum_yolu:
            pickle.dump(durum, open(durum_yolu, 'wb'))
        if ilerle:
            ilerle(pi, durum['toplam_kayit'], time.time() - t0)
        del enc, big
    return durum['res']
