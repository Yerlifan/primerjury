# -*- coding: utf-8 -*-
"""
orientation.py - THE DEFINITION AND NORMALISATION OF THE CANONICAL ORIENTATION.
A single source.

THE PROJECT'S CANONICAL ORIENTATION: SENSE (the reference, or plus, strand).
  SSU rRNA and ITS consensuses are stored in the direction the reference databases
  (SILVA/RefSeq) use. The reason: every in-silico PCR engine in this project
  (ispcr.amplify, okuma_motoru.Sonda) scans the given sequence ON THE PLUS STRAND
  ONLY. On a consensus stored in reverse, neither the forward primer nor the
  complement of the reverse primer is found, and the engine reports no product
  without raising anything.

"""
# -------------------------------------------------------------------------
# orientation.py - the ONE place the canonical orientation (SENSE) is defined: it
#          measures which direction a sequence is stored in and, where needed,
#          converts it to canonical by taking its reverse complement.
#
# INPUT  : it takes a sequence and an amplicon class directly; dosya_kanonik()
#          reads a fasta file. okuma_motoru.Sonda is used for the binding site
#          search (lossless, tolerating <=2 mismatches).
# OUTPUT : it writes no file. tespit() returns the pair (verdict, detail);
#          kanonik() the triple (canonical_sequence, verdict, was_flipped);
#          dosya_kanonik() a record list; kendini_sina() an empty list (passed) or
#          error texts. Run directly, it prints the test result to the screen.
# CALLED BY: build_canonical.py (producing canonical_consensus), hepsi.yon_kapisi (the
#          gate at the head of every stage), kendini_sina.yon_sinamasi and
#          build_consensus.py. So it runs indirectly on all of full_chain.py's
#          stages 1 to 9.
#
# WHY A SEPARATE AND SINGLE MODULE: the orientation bug was found and patched in
# three separate places on the same night. Three patches means there is no single
# canonical fix, and the bug escapes again on the next change. The decision is tied
# to two INDEPENDENT criteria (the panel's own universal pairs and the literature
# motifs); if the two disagree the sequence counts as UNCERTAIN and IS NOT
# normalised, it is flagged - because a sequence flipped the wrong way does as much
# silent damage as one never flipped at all.
# -------------------------------------------------------------------------
import os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import read_engine as om

KANONIK_YON = 'SENSE'
TOLERANS_MM = 2          # a nanopore consensus carries the occasional error

PANEL_CIFT = {
    'A':  ('Arke_universal',        'CTGCGGTTTAATTGGATTCAACGC', 'GAACTGACGACGGCCATGC'),
    'B':  ('Bakteri_universal',     'ACAAGCGGTGGAGCATGTG',      'ACGACAGCCATGCAGCAC'),
    'F1': ('Mantar_universal (F1)', 'GGTTACCCGCTGAACTTAAGC',    'CGCTTCACTCGCCGTTAC'),
    'F2': ('Mantar_universal (F2)', 'GTGCATGGCCGTTCTTAGTTG',    'CAAACTTCCATCGGCTTGAGC'),
}
MOTIF = {
    'SSU': ['GTGYCAGCMGCCGCGGTAA', 'GGATTAGATACCC', 'AGTCCCGCAACGAGCGCAACCC'],
    'ITS': ['TCCGTAGGTGAACCTGCGG', 'GCATATCAATAAGCGGAGGA'],
}
rc = om.rc
temizle = om.temizle


def sinifi(ad):
    for s, onek in (('F1', ('F1-', 'F1_')), ('F2', ('F2-', 'F2_')),
                    ('A', ('A1-', 'A2-', 'A1_', 'A2_')), ('B', ('B-', 'B_'))):
        if any(x in ad for x in onek):
            return s
    return '?'


def _var(dizi, desen, mm=TOLERANS_MM):
    return len(om.Sonda(desen, uc5=False, max_mm=mm, son2=False).bul(dizi)) > 0


def _olcut1(dizi, sinif):
    if sinif not in PANEL_CIFT:
        return 'BELIRSIZ', 0, 0
    _, F, R = PANEL_CIFT[sinif]
    d = int(_var(dizi, F)) + int(_var(dizi, rc(R)))
    t = int(_var(dizi, rc(F))) + int(_var(dizi, R))
    return ('SENSE' if d > t else 'ANTISENSE' if t > d else 'BELIRSIZ'), d, t


def _olcut2(dizi, sinif):
    tip = 'ITS' if sinif in ('F1', 'F2') else 'SSU'
    d = sum(1 for m in MOTIF[tip] if _var(dizi, m))
    t = sum(1 for m in MOTIF[tip] if _var(rc(dizi), m))
    return ('SENSE' if d > t else 'ANTISENSE' if t > d else 'BELIRSIZ'), d, t


def tespit(dizi, sinif):
    """Returns: (verdict, detail). verdict in {'SENSE','ANTISENSE','BELIRSIZ'}"""
    # The two criteria are chosen INDEPENDENTLY OF ONE ANOTHER: criterion 1 uses
    # the panel's own universal primers, and criterion 2 uses literature motifs that
    # know nothing of the panel (515F, 806R-sense, ITS1, rc(ITS4)). Two criteria
    # derived from the same source could go wrong in the same direction.
    # The decision rule is conservative: if the two agree the verdict is definite,
    # if one stays silent the speaking one holds, and IF THE TWO CONTRADICT the
    # verdict is UNCERTAIN and the sequence is not flipped. Flipping the wrong way
    # is more harmful than not flipping, because the later stages assume the file
    # is canonical.
    # Sequences shorter than 200 bp count as UNCERTAIN outright: in a fragment that
    # cannot hold all the motifs, the counts become accidental.
    dizi = temizle(dizi)
    if len(dizi) < 200:
        return 'BELIRSIZ', dict(sebep='the sequence is shorter than 200 bp', uzunluk=len(dizi))
    o1, d1, t1 = _olcut1(dizi, sinif)
    o2, d2, t2 = _olcut2(dizi, sinif)
    ay = dict(olcut1=o1, olcut1_duz=d1, olcut1_ters=t1,
              olcut2=o2, olcut2_duz=d2, olcut2_ters=t2)
    if o1 == o2 and o1 != 'BELIRSIZ':
        ay['sebep'] = 'both criteria give the same orientation'
        return o1, ay
    if o1 != 'BELIRSIZ' and o2 == 'BELIRSIZ':
        ay['sebep'] = 'criterion 1 decided on its own; criterion 2 is silent'
        return o1, ay
    if o2 != 'BELIRSIZ' and o1 == 'BELIRSIZ':
        ay['sebep'] = 'criterion 2 decided on its own; criterion 1 is silent'
        return o2, ay
    if o1 != o2:
        ay['sebep'] = 'THE TWO CRITERIA DIVERGED, so it is not normalised but masked'
        return 'BELIRSIZ', ay
    ay['sebep'] = 'both criteria are silent; no motif was found'
    return 'BELIRSIZ', ay


def kanonik(dizi, sinif):
    """Returns: (canonical_sequence, verdict, was_flipped).
        On BELIRSIZ the sequence IS NOT CHANGED - the caller must flag it.

    """
    dizi = temizle(dizi)
    karar, ay = tespit(dizi, sinif)
    if karar == 'ANTISENSE':
        return rc(dizi), karar, True
    return dizi, karar, False


def dosya_kanonik(yol):
    """Converts every record in a fasta file to canonical.
        Returns: [(header, canonical_sequence, verdict, flipped)]

    """
    sn = sinifi(os.path.basename(yol)) or sinifi(yol)
    if sn == '?':
        sn = sinifi(yol)
    out, ad, buf = [], None, []
    with open(yol, errors='ignore') as fh:
        for l in fh:
            if l.startswith('>'):
                if ad is not None:
                    out.append((ad, ''.join(buf)))
                ad, buf = l[1:].rstrip('\n'), []
            else:
                buf.append(l.strip())
    if ad is not None:
        out.append((ad, ''.join(buf)))
    son = []
    for a, s in out:
        k, karar, cev = kanonik(s, sn)
        son.append((a, k, karar, cev))
    return son, sn


def kendini_sina():
    """The module's own test. It runs before the main work starts (project rule 2)."""
    hata = []
    # 1) bilinen sense bir SSU parcasi: 515F + 806R-sense icerir
    s = ('GG' * 60 + 'GTGCCAGCAGCCGCGGTAA' + 'AC' * 120 + 'GGATTAGATACCC' + 'TT' * 60)
    k, karar, cev = kanonik(s, 'A')
    if karar != 'SENSE' or cev:
        hata.append('sentetik sense SSU yanlis: %s' % karar)
    # 2) the reverse of the same sequence must come out ANTISENSE and be flipped back
    k2, karar2, cev2 = kanonik(rc(s), 'A')
    if karar2 != 'ANTISENSE' or not cev2:
        hata.append('sentetik antisense SSU yanlis: %s' % karar2)
    elif k2 != temizle(s):
        hata.append('flipping it did not bring the sequence back')
    # 3) motifsiz dizi BELIRSIZ olmali
    _, karar3, _ = kanonik('A' * 400, 'A')
    if karar3 != 'BELIRSIZ':
        hata.append('a sequence with no motif is not marked UNDECIDED: %s' % karar3)
    # 4) idempotans: kanonigi tekrar kanoniklestirmek degistirmemeli
    k4, _, cev4 = kanonik(k2, 'A')
    if cev4 or k4 != k2:
        hata.append('idempotence is broken')
    # 5) sinif tespiti
    if sinifi('A1-4_2209_konsensus.fasta') != 'A' or sinifi('F2-1_101201.fasta') != 'F2':
        hata.append('the class detection is broken')
    return hata


if __name__ == '__main__':
    h = kendini_sina()
    if h:
        print(u'THE SELF TEST FAILED:')
        for x in h:
            print('  -', x)
        sys.exit(1)
    print(u'orientation.py self test: PASSED (5/5)  | the canonical orientation =', KANONIK_YON)
