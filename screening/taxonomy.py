# -*- coding: utf-8 -*-
"""Extracts the TAXONOMY from a FASTA header, the ONE place for every reference
database.

WHY A SEPARATE MODULE (2026-08-21, A2)
--------------------------------------
Layer 2 (the local database) was telling whether a product came from INSIDE or
OUTSIDE the target by LENGTH alone. Layer 3 (MFEprimer, D-12) answers the same
question by TAXON. It was measured in D-12: of the 1,605 amplicons MFEprimer counted
as "off target", 95.7 percent were from INSIDE the target clade and only their
length differed. So the length criterion is misleading; layer 2 has to look at the
taxon too.

The real obstacle to that is this: the files in REFERANS_DB use FIVE DIFFERENT
header formats. Assuming a single parser means counting most records as "a different
domain" (class c), that is, producing a larger version of the very fault being
fixed. The formats WERE MEASURED (2026-08-21, the first record of each file):

  SILVA  >AY846379.1.1791 Eukaryota;Archaeplastida;...;Monoraphidium sp.
         the accession is separated BY A SPACE, the taxonomy by ';'
  UNITE  >UDB016649|k__Fungi;p__Basidiomycota;...;s__Thelephora_albomarginata|SH1281904.10FU
         split by '|', the taxonomy in field 2, the tokens prefixed 'x__'
  PR2    >AB353770.1.1740_U|18S_rRNA|nucleus||Eukaryota|TSAR|Alveolata|...
         split by '|', the taxonomy from field 5 on
  RefSeq >NR_201932.1 Sphingosinicella wutangchuni strain LY54 16S ribosomal RNA...
         NO TAXONOMY. The domain comes from the database's DEFINITION (every record
         in bacteria.16S.fna is bacterial by definition). Not a guess.
  ROD    >GCA_000001215|AE014298.5/23211192-23217141|Eukaryota;Opisthokonta;...
         split by '|', the taxonomy in field 3, separated by ';'

CUTTING THE HEADER
------------------
The caller MUST NOT CUT the header. Measured: 16.6 percent of SILVA SSU headers go
over 150 characters and the tail that gets cut is exactly the GENUS and SPECIES
tokens, that is, the ones most likely to match the target clade. Classifying with a
cut header makes a record INSIDE the clade count as OUTSIDE it.

UNKNOWN IS A STATE OF ITS OWN
-----------------------------
If the domain cannot be resolved, '?' is returned and the caller COUNTS that as
'unknown', not as 'a different domain'. A record that could not be resolved IS NOT
EVIDENCE of a cross reaction.

"""
from __future__ import unicode_literals

import re

ALANLAR = ('Bacteria', 'Archaea', 'Eukaryota')
ORGANEL_JETONLARI = ('Chloroplast', 'Mitochondria')

# In databases that carry no taxonomy the domain comes from the file's DEFINITION.
# The keys are derived from the file name (non letter characters become '_').
VTB_ALAN = {
    'archaea_16S_fna': 'Archaea',
    'bacteria_16S_fna': 'Bacteria',
    'fungi_ITS_fna': 'Eukaryota',
    'fungi_28SrRNA_fna': 'Eukaryota',
    'fungi_18SrRNA_fna': 'Eukaryota',
    'ref_all2_fna': '?',        # karisik kume - alan kayittan cozulemez
    'ref_all_fna': '?',
}
VTB_KLAD = {
    'fungi_ITS_fna': 'Fungi',
    'fungi_28SrRNA_fna': 'Fungi',
    'fungi_18SrRNA_fna': 'Fungi',
}

# UNITE jeton oneki: k__Fungi -> Fungi
_UNITE_ONEK = re.compile(r'^[kpcofgs]__')
# UNITE and PR2 use '_' instead of a space in species names
_ALT_CIZGI = re.compile(r'_+')


def vtb_anahtari(dosya_adi):
    """Turns a file name into a VTB_ALAN or VTB_KLAD key."""
    import os
    t = os.path.basename(dosya_adi or '')
    return re.sub(r'\W+', '_', t).strip('_')


def _temizle(jetonlar):
    out = []
    for j in jetonlar:
        j = (j or '').strip()
        if not j:
            continue
        j = _UNITE_ONEK.sub('', j)
        out.append(j)
        # 'Thelephora_albomarginata' -> ayrica 'Thelephora albomarginata'
        if '_' in j:
            out.append(_ALT_CIZGI.sub(' ', j))
    return out


def _alan_bul(jetonlar):
    for j in jetonlar:
        if j in ALANLAR:
            return j
    # the UNITE kingdom is not a domain; Fungi -> Eukaryota
    for j in jetonlar:
        if j in ('Fungi', 'Metazoa', 'Viridiplantae', 'Archaeplastida'):
            return 'Eukaryota'
    return '?'


def coz(baslik, vtb=''):
    """A FASTA header -> (domain, [tokens], is_organelle, has_taxonomy)

    baslik : the header line WITHOUT '>' and NOT CUT
    vtb    : the database file name (needed for sources that carry no taxonomy)

    The domain can come back as '?'; that means 'unknown', it DOES NOT MEAN 'a different
    domain'.

    has_taxonomy : was a REAL taxonomy string found in the record.
      If it is False, all we have is the ORGANISM NAME (the RefSeq format:
      'NR_201932.1 Sphingosinicella wutangchuni strain LY54 16S ...').
      That distinction IS REQUIRED and was added after being measured (2026-08-21): an
      organism name carries only the GENUS and the SPECIES. A GENUS target like
      'Petrimonas' matches in the name, but a target ABOVE GENUS such as 'Bacteroidales'
      (an order) or 'Microascaceae' (a family) APPEARS in no species name. Classifying
      without knowing that, every Bacteroidales member in RefSeq was counted as "outside
      the clade": measured, 3,646 false cross reactions in bacteria.16S.fna. Where we
      cannot decide we have to say 'unknown'.

    """
    b = (baslik or '').strip()
    if b.startswith('>'):
        b = b[1:].strip()
    if not b:
        return ('?', [], False, False)

    jet = []

    if '|' in b:
        alanlar = [x.strip() for x in b.split('|')]
        # UNITE: one field holds a ';' list with 'x__' prefixes
        for a in alanlar:
            if ';' in a and _UNITE_ONEK.search(a):
                jet = _temizle(a.split(';'))
                break
        if not jet:
            # ROD: one field holds a ';' separated list starting with a DOMAIN name
            for a in alanlar:
                if ';' in a and a.split(';')[0].strip() in ALANLAR:
                    jet = _temizle(a.split(';'))
                    break
        if not jet:
            # PR2: the DOMAIN name appears among the '|' separated fields
            for i, a in enumerate(alanlar):
                if a in ALANLAR:
                    jet = _temizle(alanlar[i:])
                    break
    elif ';' in b:
        # SILVA: the accession is separated BY A SPACE, the taxonomy follows with ';'
        govde = b.split(' ', 1)[1] if ' ' in b.split(';', 1)[0] else b
        jet = _temizle(govde.split(';'))

    if jet:
        alan = _alan_bul(jet)
        if alan != '?':
            return (alan, jet, any(o in jet for o in ORGANEL_JETONLARI), True)

    # The RefSeq format: there is NO taxonomy. The domain comes from the database's definition.
    k = vtb_anahtari(vtb)
    alan = VTB_ALAN.get(k, '?')
    ad = b.split(' ', 1)[1] if ' ' in b else b
    # 'Sphingosinicella wutangchuni strain LY54 16S ...' -> ilk iki kelime cins+tur
    kelime = ad.split()
    jet = _temizle(jet + [ad] + kelime[:2] + ([' '.join(kelime[:2])] if len(kelime) > 1 else []))
    if k in VTB_KLAD:
        jet.append(VTB_KLAD[k])
    if alan == '?':
        alan = _alan_bul(jet)
    return (alan, jet, any(o in jet for o in ORGANEL_JETONLARI), False)


def sinifla(baslik, vtb, hedef_klad, hedef_alan):
    """Puts a hit into one of D-12's four classes.

      'a'          from INSIDE the target clade                 -> NOT off target
      'ao'         inside the target domain but an ORGANELLE (chloroplast or
                   mitochondrion)
      'b'          THE SAME domain, OUTSIDE the clade           -> a real cross reaction
      'c'          A DIFFERENT domain                           -> a real cross reaction
      'bilinmiyor' IT COULD NOT BE DECIDED, IT DOES NOT COUNT AS EVIDENCE

    hedef_klad : the target's clade tokens (screening/target_clades.tsv)
    hedef_alan : the target's domain ('Bacteria' / 'Archaea' / 'Eukaryota')

    'bilinmiyor' is returned in THREE cases and none of them MEANS 'no cross reaction':
      1) the domain could not be resolved at all
      2) the target's domain is unknown
      3) the record has NO taxonomy (RefSeq) and the target clade DOES NOT APPEAR in the
         organism name. In that case the record may be INSIDE the clade or OUTSIDE it;
         because an organism name carries only genus and species, membership at order or
         family level CANNOT BE READ FROM THE NAME. Measured: without this distinction,
         3,646 FALSE cross reactions were produced for Bacteroidales_kumesi in
         bacteria.16S.fna (a came out 0, which is impossible).

    """
    alan, jet, organel, taksonomi_var = coz(baslik, vtb)

    ic = any(j in jet for j in (hedef_klad or []))
    if not ic and hedef_klad:
        # Taksonomi dizgesi olmayan kayitlarda (RefSeq) jeton tam eslesmeyebilir;
        # kelime sinirinda ara. 'Methanothrix' -> 'Methanothrix soehngenii'
        duz = ' '.join(jet)
        ic = any(re.search(r'\b%s' % re.escape(j), duz, re.I) for j in hedef_klad)

    if ic:
        return 'ao' if organel else 'a'
    if alan == '?' or not hedef_alan:
        return 'bilinmiyor'
    if alan != hedef_alan:
        # A different domain: this decision NEEDS NO TAXONOMY. Even in RefSeq
        # the database's definition settles the domain (archaea.16S.fna ->
        # Archaea), so it IS CERTAIN that an archaeal record is in a different
        # domain from a bacterial target. 'c' can be given safely.
        return 'c'
    # From here on: the same domain, and NO clade match.
    # If there is a taxonomy, 'outside the clade' IS CERTAIN -> 'b'.
    # With no taxonomy WE CANNOT KNOW: an organism name carries no membership
    # above genus.
    return 'b' if taksonomi_var else 'bilinmiyor'
