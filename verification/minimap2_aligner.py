# -*- coding: utf-8 -*-
"""
minimap2_aligner.py - kimlik asamasi icin IKINCI, SECILEBILIR hizalayici.

Bu dosya mevcut motoru DEGISTIRMEZ, yanina durur. identity_verification.py icindeki
saf Python + numpy hizalayicisi (hizala) yerinde kalir ve varsayilan olmaya
devam eder. Buradaki minimap2 yolu ancak iki sart birden saglanirsa devreye
girer: mappy kurulu olacak VE kullanici acikca secmis olacak.
"""
# -------------------------------------------------------------------------
# minimap2_aligner.py
#
# INPUT  : a query sequence (a bin consensus) and the target records (the database
#          records)
# OUTPUT : (percent_identity, distance) pairs - in THE SAME form as
#          kimlik_dogrulama.hizala, so that the two can be swapped
# CALLED BY: verification/identity_verification.py (key I) and
#          verification/all_bin_identities.py (key G), only when HIZALAYICI=minimap2
#          is given
#
# WHY IT EXISTS
# Our slowest step is the identity stage: one bin consensus is scanned against twelve
# databases and 500 candidates are FULLY aligned in each. The pure Python route does
# O(len(q) x len(t)) dynamic programming for every pair, and even vectorised with
# numpy the work is still quadratic. minimap2 was designed for exactly this job: it
# aligns long, error prone reads against large reference sets with seed, chain and
# extend, and does the quadratic DP only in a narrow band around the chain.
#
# WHERE IT IS NOT USED - this boundary matters
#   * SEARCHING FOR PRIMER BINDINGS. Primers are 18-25 bases. minimap2's default
#     seed (k=15, w=10) and its chaining logic were not designed for queries of that
#     length; on short queries it fails to find a seed and SILENTLY misses the
#     binding site. There the pigeonhole engine (screening/read_engine.py) stays,
#     and must stay, because that engine's losslessness is A GUARANTEE.
#   * THE IN-SILICO PCR PRODUCT CALCULATION. The same reason: finding the places
#     where two primers bind facing one another in the right direction is a job of
#     short, near exact search. minimap2's approximate answer is not acceptable here.
#
# INSTALLATION
#     pip install mappy
#   or into the project environment:
#     micromamba install -n mikro -c bioconda minimap2
#     micromamba run -n mikro pip install mappy
#
# SAFETY: if mappy is missing, this file BREAKS NOTHING. var_mi() returns False, the
# caller carries on with the existing engine, and the chain is not broken.
# -------------------------------------------------------------------------

from __future__ import print_function

_MAPPY = None
_DENENDI = False
_SEBEP = u''


def var_mi():
    """Is mappy installed and working? It measures once, then serves from a cache.

        'it can be imported' is not enough: whether an index can really be built is
        tested too, because a broken installation can pass the import stage and blow up
        on first use, and that moment is in the middle of a run lasting hours.

    """
    global _MAPPY, _DENENDI, _SEBEP
    if _DENENDI:
        return _MAPPY is not None
    _DENENDI = True
    try:
        import mappy
    except ImportError:
        _SEBEP = u'mappy kurulu degil (pip install mappy)'
        return False
    try:
        deneme = mappy.Aligner(seq='ACGT' * 40, preset='map-ont', n_threads=1)
        if not deneme:
            _SEBEP = u'mappy indeksi bos dondu'
            return False
        list(deneme.map('ACGT' * 40))
    except Exception as e:
        _SEBEP = u'mappy kurulu ama calismadi: %s: %s' % (type(e).__name__, e)
        return False
    _MAPPY = mappy
    return True


def sebep():
    return _SEBEP


def surum():
    if not var_mi():
        return u'yok'
    return getattr(_MAPPY, '__version__', u'bilinmiyor')


# -------------------------------------------------------------------------
# HOW THE IDENTITY PERCENTAGE IS COMPUTED
#
# The existing engine (kimlik_dogrulama.hizala) returns this:
#     percent = 100 * (1 - edit_distance / len(query))
# that is, THE WHOLE query is seated inside the target (an infix/HW alignment) and
# the cost is divided by the query length.
#
# minimap2, by contrast, does local alignment: it reports only the aligned part of
# the query. Those two numbers ARE NOT THE SAME THING and cannot be compared
# directly. To make them comparable, the local result is rescaled against the query
# length: every unaligned base counts as a mismatch.
#
#     aligned_correct = blen - NM
#     unaligned       = len(q) - (q_en - q_st)
#     distance        = NM + unaligned
#     percent         = 100 * (1 - distance / len(q))
#
# That way the two engines measure the same definition and the comparison means
# something. This conversion is REQUIRED FOR THE COMPARISON TO BE VALID; skipped,
# minimap2 systematically gives a higher percentage and the two are taken to agree.
# -------------------------------------------------------------------------
def hizala_mm(q, t):
    """THE SAME signature and THE SAME return form as kimlik_dogrulama.hizala."""
    if not var_mi():
        raise RuntimeError(u'mappy yok: %s' % _SEBEP)
    if not q or not t:
        return (0.0, len(q or t or ' '))
    try:
        ind = _MAPPY.Aligner(seq=t, preset='map-ont', n_threads=1)
        if not ind:
            return (0.0, len(q))
        en_iyi = None
        for h in ind.map(q):
            if en_iyi is None or h.mlen > en_iyi.mlen:
                en_iyi = h
        if en_iyi is None:
            return (0.0, len(q))
        hizalanmayan = len(q) - (en_iyi.q_en - en_iyi.q_st)
        uzaklik = int(en_iyi.NM) + max(hizalanmayan, 0)
        yuzde = round(100.0 * (1 - uzaklik / float(len(q))), 2)
        return (max(yuzde, 0.0), uzaklik)
    except Exception:
        # One record failing to align must not fail the whole run.
        return (0.0, len(q))


def toplu_hizala(q, hedefler, iplik=3):
    """THE REAL SPEED GAIN IS HERE: ONE index, ONE mapping.

        THE 2026-08-05 FIX - A SILENT CHANGE OF ROUTE, CLOSED
        The first version called mappy.Aligner(seq=<list>). mappy DOES NOT ACCEPT that
        and raises TypeError; the code caught it and fell silently back to the SLOW route
        that builds an index per target one at a time. The measurements were still
        correct, but the claim of a "bulk index" WAS NOT TRUE and the comment said the
        opposite. Exactly the kind of bug this project chases: the code does SOMETHING
        ELSE without raising an error.

        The right route: the sequences are written to a temporary FASTA and the index is
        built from the file. That way h.ctg carries THE REAL record name and the name
        matching is not left to guesswork.

        hedefler: [(key, sequence), ...]
        Returns : {key: (percent, distance)}   - the unaligned ones as (0.0, len(q))

    """
    if not var_mi():
        raise RuntimeError(u'mappy yok: %s' % _SEBEP)
    sonuc = dict((a, (0.0, len(q))) for a, _d in hedefler)
    if not q or not hedefler:
        return sonuc

    # Since the keys become FASTA headers, they cannot carry spaces or newlines.
    # The temporary name -> real key mapping is kept separately, so that two records
    # carrying the same header are not confused with one another.
    import os as _os
    import tempfile as _tf
    esleme = {}
    fd, yol = _tf.mkstemp(suffix='.fa', prefix='mm2_')
    try:
        with _os.fdopen(fd, 'w') as fh:
            for i, (a, d) in enumerate(hedefler):
                ad = 's%d' % i
                esleme[ad] = a
                fh.write('>%s\n%s\n' % (ad, d))
        # best_n IS KEPT HIGH. By default minimap2 reports only the PRIMARY hit; in a bulk
        # index that would mean "everything except the best record is invisible", and would
        # bring back exactly THE CUT OFF we are trying to remove. If we are to use it as a
        # candidate finder, every reasonable hit has to come through; the elimination
        # decision is made by the Python scorer, not by minimap2.
        ind = _MAPPY.Aligner(fn_idx_in=yol, preset='map-ont', n_threads=iplik,
                             best_n=max(len(hedefler), 50))
        if not ind:
            return sonuc
        en_iyi = {}
        for h in ind.map(q):
            ad = h.ctg
            if ad not in en_iyi or h.mlen > en_iyi[ad].mlen:
                en_iyi[ad] = h
        for ad, h in en_iyi.items():
            anahtar = esleme.get(ad)
            if anahtar is None:
                continue
            hizalanmayan = len(q) - (h.q_en - h.q_st)
            uzaklik = int(h.NM) + max(hizalanmayan, 0)
            yuzde = round(100.0 * (1 - uzaklik / float(len(q))), 2)
            sonuc[anahtar] = (max(yuzde, 0.0), uzaklik)
        return sonuc
    finally:
        try:
            _os.unlink(yol)
        except OSError:
            pass



def secili_mi():
    """True when the user has chosen minimap2 EXPLICITLY.

        Two conditions at once: the environment variable HIZALAYICI=minimap2 AND mappy
        working. It is OFF by default. The reason: until the comparison report (see
        MINIMAP2_KARSILASTIRMA.md) shows the two engines giving the same answer, we do
        not take the fast one for the correct one. Project rule 1: no decision is left to
        a single code path, and if two measurements diverge it is not the fast one that
        wins but the one VERIFIED BY HAND.

    """
    import os
    if os.environ.get('HIZALAYICI', '').strip().lower() != 'minimap2':
        return False
    return var_mi()


if __name__ == '__main__':
    import sys
    print(u'is mappy working : %s' % (u'EVET' if var_mi() else u'HAYIR'))
    if not var_mi():
        print(u'reason           : %s' % _SEBEP)
        print(u'install            : pip install mappy')
        sys.exit(1)
    print(u'mappy surumu       : %s' % surum())
    print(u'ALIGNER selected : %s' % (u'minimap2' if secili_mi() else u'python (varsayilan)'))


# -------------------------------------------------------------------------
# THE HYBRID ROUTE - the recommended use (the 2026-08-05 measurement)
#
# WHY NOT A FULL REPLACEMENT
# The measurement showed this: minimap2 and the pure Python engine ARE NOT ASKING
# THE SAME QUESTION.
#   * Python's hizala() is an infix (HW) alignment: it FORCES the query INSIDE the
#     target and produces a number for EVERY target. Even in an irrelevant database
#     it returns values between 50 and 65 percent. Those values are NOISE, coming
#     from conserved regions and from the forced alignment.
#   * minimap2 is local alignment: with no real homology it reports NO alignment at
#     all. So in an irrelevant database it returns 0 candidates.
#
# The measured example (A1-1_1826872, 150 records per database):
#     SILVA SSU NR99 (the right locus) : minimap2 found 2 candidates, and the best
#                                        hit the hybrid chose came out THE SAME as
#                                        Python's full scan (73.2%)
#     SILVA LSU NR99 / LSU Parc / UNITE ITS (the wrong locus):
#                                        minimap2 0 candidates, Python produced
#                                        60.4% / 52.1% / 51.3% - none of them a real
#                                        match
#
# So the "divergence" is not minimap2 being wrong, it is minimap2 REFUSING TO REPORT
# RUBBISH. Those values fall far below the verdict thresholds anyway (species 98.7%,
# genus 90%) and enter no decision.
#
# THE RIGHT INTEGRATION IS THEREFORE HYBRID:
#     minimap2 FINDS THE CANDIDATES  ->  Python SCORES those candidates
# The definition of the identity percentage does not change (the existing engine's
# definition is kept), but instead of fully aligning thousands of records in every
# database, only the handful that really align are scored.
#
# AN EXTRA GAIN: THE SHORT LIST PROBLEM
# On the existing route the candidates were ranked BY SEED COUNT and the first 500
# taken, and measurement showed the cut off WAS BINDING (in one query the winner came
# from position 4171, and in 13 of 118 queries from beyond 400). minimap2 searches by
# minimizer chaining; there is no "first N" cut, and every record that aligns comes
# through. So the hybrid route REMOVES the short list cut ENTIRELY. That may be worth
# more than the speed gain.
# -------------------------------------------------------------------------
def hibrit_adaylar(q, hedefler):
    """Selects the targets minimap2 REALLY aligns.

        Returns: [(key, sequence), ...] - only the ones with an alignment.
        With no mappy it returns None rather than an empty list; the caller must then use
        the existing short list route. None and an empty list ARE DIFFERENT THINGS:
          None     -> there is no minimap2, no decision could be made, fall back
          []       -> minimap2 is there and nothing aligned in this database
        That distinction matters; confused, a database is silently skipped.

    """
    if not var_mi():
        return None
    # MEASURED 2026-08-05: toplu_hizala builds ONE index and in that mode minimap2
    # DISCARDS SECONDARY alignments. In a four target test two real homologues at 75.76%
    # appeared as "0". It cannot be used as a candidate finder: it puts a more OPAQUE cut
    # in place of the short list cut we are trying to remove. So a SEPARATE index is built
    # per target here; slower, but lossless.
    return [(a, d) for a, d in hedefler if hizala_mm(q, d)[0] > 0]
