# -*- coding: utf-8 -*-
"""THE ORDER CLASSES plus THE OLD AND NEW IDENTITY COLUMNS, from one source.

WHY A SEPARATE MODULE: both pieces of information have to appear in more than one
output (the order list, the panel table, the decision table, the combined summary).
Recomputed in each output, one would drift from the other, and in this project that
is exactly what happened once.

1) THE CLASSIFICATION (2026-08-06, the user's decision)
   Rows that cannot pass the threshold ARE NOT REMOVED from the list. Instead of
   deleting them silently they are stamped and what has to be done in the laboratory
   is written down; the decision belongs to the user.

     KESIN       dCq >= 3,0  (>= 8,00x)   -> it can be ordered
     KOSULLU     dCq 2,0-3,0 (4,00-8,00x) -> it can be ordered BUT confirmation is
                                             required
     ONERILMEZ   dCq <  2,0  (<  4,00x)   -> it appears in the list with its reason
     EVRENSEL    the ratio is undefined   -> it is judged by coverage

2) THE OLD AND NEW IDENTITY
   For each target's member bins three columns are put side by side:
     kraken_etiketi   the name discussed at the meeting (MEVCUT_KAYITLI_KIMLIK)
     olculen_kimlik   the name we confirmed (DOGRULANAN_KIMLIK)
     savunulabilir_duzey  species / genus / family / unnameable
   Where there is a difference the row is marked with ">>", so that a match like
   "what was taken for Trichoderma is really Petriella" shows at a glance.
   THE TARGET NAMES ARE NOT CHANGED, they are used as codes; a column IS ADDED.

"""
import os, csv, re, collections
from . import config as C

KOSULLU_ALT_DCQ = 2.0                       # bunun altinda ONERILMEZ
KOSULLU_ALT_KAT = 2.0 ** KOSULLU_ALT_DCQ    # 4,00x

FARK_ISARETI = u'>> FARKLI'


def _f(v):
    try:
        return float(str(v).replace(',', '.'))
    except (TypeError, ValueError):
        return None


def sinifla(kat, olculemedi=False):
    'Returns (class, dcq, how far it is from the threshold in dCq, the laboratory note).'
    if olculemedi or kat is None:
        return (u'EVRENSEL', None, None,
                'The separation ratio is undefined, since the competitor set '
                'is empty. The verdict goes by COVERAGE. In the laboratory a '
                'melting curve plus a no template control is enough, and a '
                'single product peak is expected.')
    d = C.kat_dcq(kat)
    if d is None:
        # fold = 0 -> the log is undefined. NO sentinel number IS PRODUCED (it used to
        # write -99 and that appeared in the output as "-99,00 dCq"); dCq is left empty.
        return (u'ONERILMEZ', None, None,
                'The separation ratio is ZERO: the competitor bins are '
                'amplified as well as the target or better, and no dCq can be '
                'computed. A signal measured with this pair CANNOT BE '
                'ATTRIBUTED to the target. Laboratory confirmation WILL NOT '
                'rescue this row; the primer has to be redesigned.')
    fark = round(d - C.ESIK_DCQ, 2)
    if kat >= C.AYRIM_ESIK:
        return (u'KESIN', d, fark,
                'A melting curve plus a no template control is enough. With a '
                'single product peak and a clean control, the order counts as '
                'confirmed.')
    if kat >= KOSULLU_ALT_KAT:
        return (u'KOSULLU', d, fark,
                'CLOSE TO THE THRESHOLD (dCq %.2f, which is %.2f dCq short). '
                'A melting curve IS NOT ENOUGH ON ITS OWN, because a cross '
                'product melts at the same temperature. A GEL IS NEEDED to '
                'confirm the product length, and a no template control is '
                'required. If the gel does not give a single band, or the '
                'length does not match, amplicon SEQUENCING is needed.'
                % (d, fark))
    return (u'ONERILMEZ', d, fark,
            'The separation is far too low (dCq %.2f, which is %.2f dCq short '
            'of the threshold): the competitor bins are amplified as well as '
            'the target. A signal measured with this pair CANNOT BE '
            'ATTRIBUTED to the target. If it is to be ordered, AMPLICON '
            'SEQUENCING is mandatory; a gel and a melting curve cannot decide '
            'this row. The advice: redesign the primer first.' % (d, fark))


def kimlik_tablosu(kok):
    """target -> dict(kraken, olculen, duzey, farkli_mi, uye_sayisi, ayrinti)

    The source: ALL_IDENTITIES_RESULT/tum_kutu_kimlikleri.tsv (stage G).
    If the file is missing it returns EMPTY and the caller writes '-' in the columns; a
    wrong name IS NOT INVENTED silently.

    """
    yol = os.path.join(kok, 'ALL_IDENTITIES_RESULT', 'tum_kutu_kimlikleri.tsv')
    if not os.path.exists(yol):
        return {}
    with open(yol, encoding='utf-8', errors='ignore') as fh:
        sat = [s for s in fh if not s.startswith('#')]
    r = list(csv.reader(sat, delimiter='\t'))
    if len(r) < 2:
        return {}
    h = r[0]
    kutular = [dict(zip(h, x)) for x in r[1:] if len(x) > 3]
    hed = collections.defaultdict(list)
    for d in kutular:
        for t in (d.get('UYESI_OLDUGU_HEDEFLER') or '').split(','):
            t = t.strip()
            if t:
                hed[t].append(d)

    def _kisalt(v):
        """Makes the long 'an unnameable lineage - THE NEAREST RECORD: <a very
        long taxonomy> (N per cent)' string readable. The number and the meaning
        are kept and the taxonomy is dropped; the full form is already in the
        bin identity report.

        """
        v = (v or '').strip()
        if v.lower().startswith('an unnameable'):
            m = re.search(r'\(\s*([\d,\.]+) per cent\)\s*$', v)
            return (u'an unnameable lineage (%s per cent)' % m.group(1)) if m \
                else u'an unnameable lineage'
        return v

    def _teklestir(vals):
        """Reduces repeated values to one, keeps the order and drops the empty ones."""
        out = []
        for v in vals:
            v = _kisalt(v)
            if v and v not in out:
                out.append(v)
        return out

    tablo = {}
    for t, ds in hed.items():
        kr = _teklestir(d.get('MEVCUT_KAYITLI_KIMLIK') for d in ds)
        ol = _teklestir(d.get('DOGRULANAN_KIMLIK') for d in ds)
        dz = _teklestir(d.get('SAVUNULABILIR_DUZEY') for d in ds)
        farkli = any((d.get('UYUSUYOR_MU') or '').upper().startswith('HAYIR')
                     for d in ds)
        tablo[t] = dict(
            kraken=' | '.join(kr) if kr else '-',
            olculen=' | '.join(ol) if ol else '-',
            duzey=' | '.join(dz) if dz else '-',
            farkli_mi=farkli,
            uye_sayisi=len(ds),
            farkli_kutu=sum(1 for d in ds
                            if (d.get('UYUSUYOR_MU') or '').upper().startswith('HAYIR')),
        )
    return tablo


def kimlik_sutun_basliklari():
    return ['kraken_etiketi', 'olculen_kimlik', 'savunulabilir_duzey', 'ad_farkli_mi']


def kimlik_sutunlari(tablo, hedef):
    k = tablo.get(hedef)
    if not k:
        # THE REASON is written instead of a silent '-': it means the target has no
        # membership definition at stage G (an example: the 'ek' pairs outside the panel).
        # Missing data and filled data must not be mixed up.
        y = 'there is NO membership definition at the membership stage'
        return [y, y, y, u'karsilastirilamadi']
    return [k['kraken'], k['olculen'], k['duzey'],
            ('%s (%d of %d bins)' % (FARK_ISARETI, k['farkli_kutu'], k['uye_sayisi']))
            if k['farkli_mi'] else u'ayni']


# -------------------------------------------------------------------------
# THE SOURCE OF THE IDENTITY, the answer to "how do you know that?"
#
# The HER_VTB_NE_DEDI column already carried this information, but it was buried in
# one enormous cell and nobody read it. It is split into separate columns:
#
#   karar_veren_vtb    the database or databases that gave the identity
#   kimlik_yuzdesi     the alignment identity in that database
#   hizalanan_uzunluk  over how many bases it was measured (the reference record's
#                      length)
#   uyusan_vtb_sayisi  how many databases said THE SAME genus / how many answered
#   uyusmayan_vtb      WHAT the ones that disagreed said; it is not hidden, it is
#                      written beside them
#
# The aim: when a reviewer says "how do you know that", the row can be shown.
# -------------------------------------------------------------------------
_VTB_PAT = re.compile(r'([A-Za-z0-9_ .]+?) \[([^\]]*)\]: (.*?)(?=(?: \| [A-Za-z0-9_ .]+? \[)|$)')


def _cins(ad):
    """Extracts the GENUS name from a hit header."""
    ad = (ad or '').strip()
    if not ad or ad.startswith('-'):
        return ''
    if ad.lower().startswith('an unnameable'):
        return u'(unnameable)'
    m = re.search(r'[;|]g__([A-Za-z0-9_]+)', ad)          # UNITE
    if m:
        return m.group(1).replace('_', ' ').split()[0]
    if ';' in ad:                                          # SILVA / PR2 / ROD
        parca = [p.strip() for p in ad.split(';') if p.strip()]
        for p in reversed(parca):
            p = p.split('(')[0].strip()
            k = re.match(r'([A-Z][A-Za-z0-9_-]{2,})', p)
            if k and k.group(1).lower() not in (
                    'uncultured', 'unidentified', 'metagenome'):
                return k.group(1)
        return ''
    m = re.match(r'\s*[A-Z]{1,2}[A-Z_]*\d+\.\d+\s+([A-Z][a-z]{3,})', ad)   # RefSeq/nt
    if m:
        return m.group(1)
    m = re.search(r'\b([A-Z][a-z]{3,})\b', ad)
    return m.group(1) if m else ''


def vtb_ayristir(hucre):
    """The per database cell -> [(database, record_count, result_text, per cent), ...]"""
    out = []
    for m in _VTB_PAT.finditer(hucre or ''):
        vtb, kapsam, sonuc = m.group(1).strip(), m.group(2), m.group(3).strip()
        y = re.search(r'\(%\s*([\d,\.]+)\)', sonuc)
        out.append((vtb, kapsam, sonuc, y.group(1) if y else ''))
    return out


def kaynak_sutun_basliklari():
    return ['karar_veren_vtb', 'kimlik_yuzdesi', 'hizalanan_uzunluk',
            'uyusan_vtb_sayisi', 'uyusmayan_vtb_ne_dedi']


def kaynak_sutunlari_kutu(d):
    """Produces the source columns from ONE BIN row of stage G."""
    if not d:
        return ['-'] * 5
    kazanan = (d.get('en_iyi_vtb') or '-').strip()
    yuzde = (d.get('en_iyi_kimlik_%') or '-').strip()
    basl = d.get('en_iyi_isabet') or ''
    # hizalanan uzunluk: SILVA basliklarinda kayit adi '.<bas>.<son>' tasir
    m = re.match(r'\S*?\.(\d+)\.(\d+)\s', basl)
    uzun = str(abs(int(m.group(2)) - int(m.group(1))) + 1) if m else '-'
    # If the bin's confirmed name is 'an unnameable lineage' the genus comparison CANNOT
    # BE MADE from it; the winning hit's own genus is used.
    _da = (d.get('DOGRULANAN_KIMLIK') or '').strip()
    hedef_cins = _cins(basl if (not _da or _da.lower().startswith('an unnameable'))
                       else _da)
    uy, uym = [], []
    for vtb, _kap, sonuc, _y in vtb_ayristir(d.get('HER_VTB_NE_DEDI', '')):
        if re.match(r'(SONUC YOK|SORULMADI|ELLE|TAMAM\s*$|-\s*$)', sonuc):
            continue
        c = _cins(sonuc)
        if not c:
            continue
        # THE PREFIX COMPARISON: the source TSV headers CUT at 110 characters, so the last
        # piece of the taxonomy line can be left half written ('Methanosarcina' ->
        # 'Methanosarcin'). Had exact equality been required, the same database would have
        # looked as though it 'disagreed' with the identity it gave itself, and that is
        # exactly what was seen. The comparison is made over the first 8 characters.
        _e = lambda x: (x or '')[:8].lower()
        (uy if (hedef_cins and _e(c) == _e(hedef_cins)) else uym).append(
            u'%s: %s' % (vtb, sonuc[:52]))
    cevap = len([1 for _v, _k, s, _y in vtb_ayristir(d.get('HER_VTB_NE_DEDI', ''))
                 if not re.match(r'(SONUC YOK|SORULMADI|ELLE)', s)])
    return [kazanan, yuzde, uzun,
            u'%d / %d cevap veren' % (len(uy), cevap) if cevap else '-',
            (u' || '.join(uym[:6]) if uym else 'none; no database disagreed')]


def kaynak_tablosu(kok):
    """target -> the source columns (from the row giving the BEST identity among the member bins)."""
    yol = os.path.join(kok, 'ALL_IDENTITIES_RESULT', 'tum_kutu_kimlikleri.tsv')
    if not os.path.exists(yol):
        return {}
    with open(yol, encoding='utf-8', errors='ignore') as fh:
        sat = [s for s in fh if not s.startswith('#')]
    r = list(csv.reader(sat, delimiter='\t'))
    if len(r) < 2:
        return {}
    h = r[0]
    kutular = [dict(zip(h, x)) for x in r[1:] if len(x) > 3]
    hed = collections.defaultdict(list)
    for d in kutular:
        for t in (d.get('UYESI_OLDUGU_HEDEFLER') or '').split(','):
            t = t.strip()
            if t:
                hed[t].append(d)

    def _yuz(d):
        try:
            return float((d.get('en_iyi_kimlik_%') or '0').replace(',', '.'))
        except ValueError:
            return 0.0

    out = {}
    for t, ds in hed.items():
        # the bin representing the target: the one with the HIGHEST identity (it decides)
        out[t] = kaynak_sutunlari_kutu(sorted(ds, key=_yuz, reverse=True)[0])
    return out


def kaynak_sutunlari(tablo, hedef):
    return tablo.get(hedef, ['there is NO membership definition at the '
                             'membership stage'] * 3 + ['-', '-'])
