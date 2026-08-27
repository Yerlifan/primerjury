# -*- coding: utf-8 -*-
"""MULTI LOCUS IDENTIFICATION: a shared decision module for fungal bins.

THE PROBLEM
  The F1 and F2 libraries carry A WHOLE OPERON: 18S + ITS + 28S. In F2 the median
  length is 3715 bp. Against that, identification was being asked only of
  `fungi.ITS.fna`, so barely 500 bases of the sequence were used and the remaining
  3200 were never weighed at all.

  MEASURED (2026-08-27, 44 fungal bins): with the three loci asked separately, 16
  bins produced information ITS alone could not give. For example:
    F2-1_40559  ITS: Geomyces sp.                 95.86 per cent /  483 bp
                18S: Pseudogymnoascus destructans 99.77 per cent / 1729 bp
                28S: Pseudogymnoascus sp.         99.09 per cent / 1315 bp
  18S gave a SPECIES name over 1729 bases while ITS stopped at genus over 483.

THE PRINCIPLE: RAW IDENTITIES ARE NOT MIXED TOGETHER
  99 per cent in 28S and 99 per cent in ITS are not the same thing. So each locus
  is FIRST decided on its own, with its own threshold; the combination is made over
  DECISIONS, not over percentages. For the same reason the hits of the three loci
  are not thrown into one pool.

THE COMBINATION RULE
  1. The loci that gave a usable answer are collected.
  2. A GENUS vote: if a genus comes from at least two loci, that genus is taken. If
     all of them differ, the genus is uncertain and all three are written down.
  3. THE SPECIES is taken only from ITS or 28S (18S cannot separate fungal
     species). If both give a species, ITS has priority; if they say DIFFERENT
     species, "cf." goes in.
  4. WHICH locus produced the name and over how many bases IS RECORDED. Without
     that the name cannot be defended.

SYNONYMS
  This module does not resolve synonyms. When two loci give a different genus that
  MAY be a real divergence, but it may equally be a difference of nomenclature. It
  was measured: *Geomyces destructans* and *Pseudogymnoascus destructans* are the
  same organism. So a divergence is not reported as "AN ERROR" but as "the loci
  diverge, they may be synonyms", and the decision is left to the reader.
"""
from __future__ import print_function
import os
import sys

_BURA = os.path.dirname(os.path.abspath(__file__))
if _BURA not in sys.path:
    sys.path.insert(0, _BURA)
from identity_verification import (TUR_ESIGI, CINS_ESIGI, AYRIM_PAYI,   # noqa: E402
                                   EN_AZ_HIZALAMA)

# The loci of a fungal bin: (label, database files, threshold key). The order does
# not matter; the decision goes by length and by vote, not by position.
MANTAR_LOKUSLARI = [
    (u'18S', ['fungi.18SrRNA.fna'], 'SSU'),
    (u'ITS', ['fungi.ITS.fna', 'UNITE_ITS.fasta'], 'ITS'),
    (u'28S', ['fungi.28SrRNA.fna'], 'LSU_MANTAR'),
]
EN_AZ_KANIT = 250       # below this many bases no name is given at all

# THE DISCRIMINATING POWER OF A LOCUS, WHICH IS NOT THE ALIGNMENT LENGTH.
#
# The first version said "take the longest alignment" and it was wrong. MEASURED:
#   F2-1_2034170  18S: Lomentospora cf. prolificans 99.83 per cent / 1719 bp
#                 ITS: Petriella sp.                97.25 per cent /  509 bp
#                 28S: Parascedosporium sp.         98.24 per cent / 1193 bp
# The rule picked 18S and made the name *Lomentospora*, because its alignment was
# the longest. But 18S is a highly CONSERVED region in fungi; being long does not
# make it discriminating, quite the opposite, similar sequences stay similar all
# the way along. The formal barcode of fungi is ITS (Schoch et al. 2012, PNAS);
# 28S (D1/D2) comes second; 18S is weak for separating genus and species.
#
# Therefore:
#   * a tie in the GENUS vote is broken by DISCRIMINATING POWER.
#   * a SPECIES name can only be taken from ITS or 28S. 18S contributes to the
#     genus agreement but CANNOT PRODUCE a species name on its own.
LOKUS_GUCU = {u'ITS': 3, u'28S': 2, u'18S': 1}
TUR_VEREBILEN = (u'ITS', u'28S')


def lokus_karari(isabetler, anahtar, ad_ayikla, cins_epitet, adsiz_izler):
    """The decision one locus reaches with its own threshold.

    Returns a dict: ad, cins, tur (when species level was reached), kimlik,
    hizalama, duzey ('tur' | 'cins' | 'yok'), notu.
    """
    bos = dict(ad=u'no match', cins=None, tur=None, kimlik=0.0,
               hizalama=0, duzey=u'yok', notu=u'no hit at this locus')
    if not isabetler:
        return bos
    h = sorted(isabetler, key=lambda x: (-x[0], -x[2]))
    pid, aln, _b, _q, tit = h[0]
    if aln < EN_AZ_KANIT:
        return dict(bos, ad=u'cannot be named', kimlik=pid, hizalama=aln,
                    notu=u'the evidence is %d bp, the floor is %d bp'
                         % (aln, EN_AZ_KANIT))
    ad = ad_ayikla(tit)
    if not ad or any(j in ad.lower() for j in adsiz_izler):
        return dict(bos, ad=u'cannot be named', kimlik=pid, hizalama=aln,
                    notu=u'the record matched is unnamed')
    te, ce = TUR_ESIGI.get(anahtar, 98.7), CINS_ESIGI.get(anahtar, 94.5)
    cins, epitet = cins_epitet(ad)
    if pid < ce:
        return dict(bos, ad=u'cannot be named', kimlik=pid, hizalama=aln,
                    notu=u'the identity is %.2f per cent, below the genus threshold '
                         u'of %.2f per cent' % (pid, ce))
    if pid < te or not epitet:
        return dict(ad=u'%s sp.' % cins, cins=cins, tur=None, kimlik=pid,
                    hizalama=aln, duzey=u'cins',
                    notu=u'below the species threshold of %.2f per cent' % te)
    if aln < EN_AZ_HIZALAMA.get(anahtar, 600):
        return dict(ad=u'%s sp.' % cins, cins=cins, tur=None, kimlik=pid,
                    hizalama=aln, duzey=u'cins',
                    notu=u'the alignment is %d bp and the floor for %s is %d bp'
                         % (aln, anahtar, EN_AZ_HIZALAMA.get(anahtar, 600)))
    # species level; a close rival WITHIN THE SAME LOCUS means "cf."
    for x in h[1:]:
        a2 = ad_ayikla(x[4])
        if a2 and a2 != ad and cins_epitet(a2)[1] and (pid - x[0]) < AYRIM_PAYI:
            return dict(ad=u'%s cf. %s' % (cins, epitet), cins=cins, tur=None,
                        kimlik=pid, hizalama=aln, duzey=u'cins',
                        notu=u'%s is only %.2f per cent behind, at %.2f per cent'
                             % (a2, pid - x[0], x[0]))
    return dict(ad=ad, cins=cins, tur=ad, kimlik=pid, hizalama=aln,
                duzey=u'tur', notu=u'')


def birlestir(kararlar):
    """Tie the per locus decisions into one result.

    kararlar: {locus_label: the output of lokus_karari()}
    Returns: (ad, kimlik_metni, notu)
    """
    kullanilir = {k: v for k, v in kararlar.items() if v['cins']}
    if not kullanilir:
        # none of them gave a name, so carry the most informative reason
        for k in (u'ITS', u'28S', u'18S'):
            if k in kararlar and kararlar[k]['notu']:
                return (u'cannot be named', u'%.2f' % kararlar[k]['kimlik'],
                        u'none of the three loci could give a name (%s: %s)'
                        % (k, kararlar[k]['notu']))
        return u'cannot be named', u'-', u'no hit at any of the three loci'

    # --- 2. the genus vote ---
    oy = {}
    for k, v in kullanilir.items():
        oy.setdefault(v['cins'], []).append(k)
    # The VOTE count first, then DISCRIMINATING POWER on a tie, and the alignment
    # length last.
    kazanan = sorted(oy.items(),
                     key=lambda x: (-len(x[1]),
                                    -max(LOKUS_GUCU.get(l, 0) for l in x[1]),
                                    -max(kullanilir[l]['hizalama'] for l in x[1])))[0]
    cins, veren = kazanan[0], kazanan[1]
    ayrisma = u''
    if len(oy) > 1:
        ayrisma = (u'the loci diverge (%s), they may be synonyms and this has to be '
                   u'checked'
                   % u' / '.join(u'%s: %s' % (u'+'.join(sorted(l)), c)
                                 for c, l in sorted(oy.items())))

    # --- 3. choosing the species: the species level with the longest alignment ---
    # A species name is taken ONLY from ITS or 28S; 18S cannot separate fungal
    # species.
    turler = [(LOKUS_GUCU.get(k, 0), v['hizalama'], k, v)
              for k, v in kullanilir.items()
              if v['duzey'] == u'tur' and v['cins'] == cins
              and k in TUR_VEREBILEN]
    _atlanan18s = [k for k, v in kullanilir.items()
                   if v['duzey'] == u'tur' and k not in TUR_VEREBILEN]
    if not turler:
        notlar = [u'the genus is %s (from %s)' % (cins, u'+'.join(sorted(veren)))]
        if _atlanan18s:
            notlar.append(u'%s did give a species name, but 18S cannot separate fungal '
                          u'species, so it was not accepted'
                          % u'+'.join(sorted(_atlanan18s)))
        if ayrisma:
            notlar.append(ayrisma)
        for k in sorted(kullanilir):
            if kullanilir[k]['notu']:
                notlar.append(u'%s: %s' % (k, kullanilir[k]['notu']))
        en = max(kullanilir.values(), key=lambda v: v['hizalama'])
        return (u'%s sp.' % cins, u'%.2f' % en['kimlik'], u'; '.join(notlar))

    turler.sort(reverse=True)
    _guc, aln, lok, v = turler[0]
    farkli = [(k2, v2['tur']) for _g, _a, k2, v2 in turler[1:]
              if v2['tur'] != v['tur']]
    notlar = [u'from the %s locus, over an alignment of %d bp' % (lok, aln)]
    if farkli:
        notlar.append(u'%s says a different species (%s)'
                      % (u'+'.join(k for k, _ in farkli),
                         u', '.join(t for _, t in farkli)))
        if ayrisma:
            notlar.append(ayrisma)
        cins_, epitet_ = v['cins'], v['tur'].split()[-1]
        return (u'%s cf. %s' % (cins_, epitet_), u'%.2f' % v['kimlik'],
                u'; '.join(notlar))
    if len(veren) > 1:
        notlar.append(u'%s confirms the same genus' % u'+'.join(sorted(veren)))
    if ayrisma:
        notlar.append(ayrisma)
    return v['tur'], u'%.2f' % v['kimlik'], u'; '.join(notlar)


# =====================================================================
# THE METHOD AS REPORTED
# =====================================================================
# This is the method the study PRESENTED. The threshold based decision chain above
# is DIFFERENT from it and does not replace it.
#
# From the document, in its own words:
#   "The matches were filtered so that the alignment was at least 250 bases long
#    and the identity at least 90 per cent. Among the matches that met that
#    condition THERE IS NO THRESHOLD FOR NAMING; whichever named species gave the
#    highest percentage was written down. The second nearest species and the gap
#    between them were kept on every row as well."
#   "The discriminating region of the class is looked at first: 16S in archaea and
#    bacteria, ITS in fungi. If nothing comes back from ITS then 28S is used, and
#    if that gives nothing either then 18S."
#
# So the locus rule is A LADDER, not a vote; and there is no threshold for a
# species name. The sharpness of a decision is shown not by a threshold but by THE
# GAP TO THE SECOND NEAREST SPECIES, so that a reader can make up their own mind.
#
# A note from 2026-08-27: the first form of this module applied a method with
# thresholds and voting. That method WAS NOT the one that had been reported; it had
# been written without looking at the document. The function below was added so
# that the result produced and the result reported do not part company, and it
# became THE PRIMARY column in the table.
BASAMAK = [u'ITS', u'28S', u'18S']
EN_AZ_OZDESLIK = 90.0


def raporlanan_yontem(lokus_isabet, ad_ayikla, cins_epitet, adsiz_izler,
                      basamak=None):
    """The method as reported. Returns: (ad, kimlik, hizalama, lokus, ikinci, fark).

    lokus_isabet: {locus: hits}. For a single locus library one key is given and
    the ladder consists of that key alone.
    """
    for lok in (basamak or BASAMAK):
        isb = lokus_isabet.get(lok) or []
        adli = []
        for pid, aln, _bit, _q, tit in sorted(isb, key=lambda x: (-x[0], -x[1])):
            if aln < EN_AZ_KANIT or pid < EN_AZ_OZDESLIK:
                continue
            ad = ad_ayikla(tit)
            if not ad or not cins_epitet(ad)[1]:
                continue
            if any(j in ad.lower() for j in adsiz_izler):
                continue
            adli.append((pid, aln, ad))
        if not adli:
            continue
        pid, aln, ad = adli[0]
        ikinci = next((x for x in adli[1:] if x[2] != ad), None)
        return (ad, pid, aln, lok,
                ikinci[2] if ikinci else None,
                (pid - ikinci[0]) if ikinci else None)
    return None, 0.0, 0, None, None, None
