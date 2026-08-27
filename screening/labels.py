# -*- coding: utf-8 -*-
"""English labels for verdicts and output columns.

WHY THIS IS A LOOKUP AND NOT A RENAME
-------------------------------------
The verdict values (DOGRULANDI, CELISKILI, TEMIZ, ...) and the TSV column names
are not merely identifiers. They are DATA:

  * compared as strings in dozens of places across the pipeline
  * written into TSV files that other stages read back
  * stored inside checkpoint files from previous runs

Renaming them would therefore be a schema migration, not a translation: every
existing result file and checkpoint would stop being readable, silently, and
the failure would surface as an empty table rather than as an error.

So the internal values stay exactly as they are, and this module supplies the
English wording where a HUMAN reads the output: Markdown reports, screen
messages, and the legend printed at the top of each TSV. Machine-readable
column names are left untouched on purpose; a program that reads these files
keeps working.

Translating an interface must not translate the data underneath it.
"""
from __future__ import unicode_literals


# --- verdicts ---------------------------------------------------------------
# Identity verification
KIMLIK = {
    'DOGRULANDI':    'VERIFIED',
    'DUZELTILMELI':  'NEEDS CORRECTION',
    'DOGRULANAMADI': 'UNVERIFIED',
}

# Specificity verification
OZGULLUK = {
    'KESIN':     'CERTAIN',
    'KOSULLU':   'CONDITIONAL',
    'INCELEME':  'NEEDS REVIEW',
    'RISKLI':    'RISKY',
    'CELISKILI': 'CONTRADICTORY',
    'EKSIK':     'INCOMPLETE',
}

# Per-layer readings
KATMAN = {
    'TEMIZ':      'CLEAN',
    'RISKLI':     'RISKY',
    'BILINMIYOR': 'UNKNOWN',
}

# Taxonomic classes (D-12)
SINIF = {
    'a':          'inside target clade',
    'ao':         'inside clade, organelle',
    'b':          'same domain, outside clade',
    'c':          'different domain',
    'bilinmiyor': 'undecidable',
}

# Taxonomic levels a claim can defend
DUZEY = {
    'TUR':                  'SPECIES',
    'CINS':                 'GENUS',
    'CINS (tur belirsiz)':  'GENUS (species uncertain)',
    'CINS (hizalama kisa)': 'GENUS (the alignment is short)',
    'CINS (ikinci isabet adsiz, ayni yakinlikta)':
        'GENUS (the second hit is unnamed and just as close)',
    'AILE ve USTU (ad VERILEMEZ)': 'FAMILY OR ABOVE (cannot be named)',
    'ADLANDIRILAMIYOR (referans adsiz)': 'CANNOT BE NAMED (reference is unnamed)',
    'BELIRLENEMEDI':        'UNDETERMINED',
}

_HEPSI = {}
for _d in (KIMLIK, OZGULLUK, KATMAN, DUZEY):
    _HEPSI.update(_d)


def en(deger, varsayilan=None):
    """Internal verdict -> English wording. Unknown values pass through.

    Passing an unrecognised value through unchanged is deliberate: a verdict
    this module has not heard of must stay visible, not be replaced by a
    plausible-looking guess.
    """
    if deger is None:
        return varsayilan
    d = str(deger).strip()
    # 'RISKLI - siparis edilmez' gibi ekli hukumler: bas kismi cevrilir
    for k in sorted(_HEPSI, key=len, reverse=True):
        if d == k:
            return _HEPSI[k]
        if d.startswith(k + ' '):
            return _HEPSI[k] + d[len(k):]
    return d


# --- output columns ---------------------------------------------------------
# Column NAMES are not translated in the files (they are the machine contract).
# These glosses are printed as a legend so a reader can follow the table.
SUTUN = {
    'hedef':            'target',
    'cift_turu':        'pair type',
    'urun_bp':          'product length (bp)',
    'KARAR':            'VERDICT',
    'kaynak_sayisi':    'number of sources that measured',
    'uyusan_kaynak':    'sources that agree',

    '1_NUMUNE_oy_vermez': 'layer 1, in-sample (shown, does NOT vote)',
    '1_numune_deger':     'layer 1 value',

    '2_YEREL_DB':                     'layer 2, local databases',
    '2_hedef_disi_urun':              'off-target products BY SIZE',
    '2_ayni_boyda_HEDEFIN_KENDISI':   'products at the expected size (the target itself)',
    '2_tum_vurus':                    'all hits',
    '2_kume_dagilimi':                'distribution across database sets',
    '2_klad_ayrimi_yapildi':          'was taxonomic separation possible?',
    '2_klad_disi_b_c':                'outside clade (b + c), TAXONOMIC measure',
    '2_a_klad_ici':                   'a: inside the target clade',
    '2_ao_organel':                   'ao: inside clade, organelle',
    '2_b_ayni_alan_klad_disi':        'b: same domain, outside clade',
    '2_c_farkli_alan':                'c: different domain',
    '2_bilinmiyor':                   'undecidable, NOT evidence of cross-reaction',

    '3_MFEPRIMER':                    'layer 3, MFEprimer (external tool)',
    '3_hedef_disi_amplikon':          'off-target amplicons, taxonomically filtered',
    '3_HAM_boya_dayali':              'raw count, size-based',
    '3_klad_ayrimi_yapildi':          'was clade separation applied?',
    '3_a_klad_ici_uzunluk_varyanti':  'a: length variant inside the clade',
    '3_ao_organel':                   'ao: organelle',
    '3_b_ayni_alan_klad_disi':        'b: same domain, outside clade',
    '3_c_farkli_alan':                'c: different domain',
    '3_olusabilir_Tm_yakin':          'can form (Tm near the annealing temperature)',
    '3_olusmaz_Tm_dusuk':             'cannot form (Tm too low)',

    '4_NCBI':            'layer 4, NCBI Primer-BLAST',
    '4_hedef_disi_urun': 'off-target products',
    '4_durum':           'status (a capped or empty page is NOT a count)',

    'HUKUM':                'VERDICT',
    'SAVUNULABILIR_DUZEY':  'defensible taxonomic level',
    'ONERILEN_AD':          'proposed name',
    'adlandirma_gerekcesi': 'reasoning for the name',
    'EN_YAKIN_5_ORGANIZMA': 'five nearest organisms (deduplicated by organism)',
    'sorgulanan_vtb_sayisi': 'databases queried',
    'sonuc_veren_vtb':      'databases that returned a result',
    'HER_VTB_NE_DEDI':      'what each database said',
    'kanit':                'evidence',
}


def legend(sutunlar):
    """English legend for a TSV, as '# ' comment lines.

    Only columns this module knows are listed; an unexplained column is better
    than an invented explanation.
    """
    satir = ['# COLUMN LEGEND (English). Column names themselves are left in the',
             '# original language on purpose: they are the machine-readable',
             '# contract that other stages and existing result files depend on.']
    for s in sutunlar:
        if s in SUTUN:
            satir.append('#   %-32s %s' % (s, SUTUN[s]))
    return '\n'.join(satir) + '\n'


def verdict_legend(sozluk, baslik):
    """English legend for a set of verdict values."""
    satir = ['# %s' % baslik]
    for k in sozluk:
        satir.append('#   %-22s %s' % (k, sozluk[k]))
    return '\n'.join(satir) + '\n'
