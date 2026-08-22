# -*- coding: utf-8 -*-
"""The output: TSV for a machine, a readable Markdown report for a person.

The question the report has to answer:
  "For this target, under which parameter setting is there a solution, is there
  none, and if there is, what does it cost?"

"""
# -------------------------------------------------------------------------
# report.py writes the search results as TSV for a machine and as Markdown for a
#            person; for each target it produces a one paragraph decision sentence.
#
# INPUT  : the kontrol/hedef_*.json results kontrol.hepsi() reads from disk (each
#          one is the dictionary __main__.hedefi_isle returns: the baseline
#          measurement, the grid table, the candidate list, the global scan result);
#          the panel rows and the panel path hedefler.panel_oku() gives.
# OUTPUT : SCREENING_RESULT/adaylar.tsv, parametre_izgarasi.tsv and
#          KAPSAMLI_ARAMA_RAPORU.md. uret() returns the list of those three paths.
# CALLED BY: inside __main__.aramayi_kos, once after each target finishes (quietly)
#          and once more at the end of the run. So it runs through
#          verification/full_chain.py keys 1, 2, 3, 7 and 9 (the 7th stage inside
#          it).
# -------------------------------------------------------------------------
import os, time, csv
from . import config as C


def aday_ozet(c):
    d = dict(
        F=c['F'], R=c['R'], urun=c['urun'],
        F_uz=c['mF']['uz'], F_tm=c['mF']['tm'], F_gc=c['mF']['gc'],
        F_uc=c['mF']['uc'], F_son5=c['mF']['son5'],
        F_hp_tm=c['mF']['hp_tm'], F_hd_tm=c['mF']['hd_tm'],
        F_hp_dg=c['mF']['hp_dg'], F_hd_dg=c['mF']['hd_dg'],
        R_uz=c['mR']['uz'], R_tm=c['mR']['tm'], R_gc=c['mR']['gc'],
        R_uc=c['mR']['uc'], R_son5=c['mR']['son5'],
        R_hp_tm=c['mR']['hp_tm'], R_hd_tm=c['mR']['hd_tm'],
        R_hp_dg=c['mR']['hp_dg'], R_hd_dg=c['mR']['hd_dg'],
        arms=c.get('arms', ''),
        izgara_hucresi=c.get('izgara_hucresi', ''), sikilik=c.get('sikilik', ''),
    )
    cf = c.get('cift', {})
    d.update({'cift_' + k: v for k, v in cf.items()})
    nm = c.get('numune', {})
    d['numune_olcut'] = nm.get('olcut', '')
    n3 = c.get('numune_mm3') or {}
    d['numune_olcut_2'] = n3.get('olcut', '')
    d['numune_kat_enkotu_mm3'] = n3.get('kat_enkotu', '')
    d['numune_kat_havuz_mm3'] = n3.get('kat_havuz', '')
    d['numune_uye_kapsam_mm3'] = n3.get('uye_kapsam_pay', '')
    for k in ('uye_alt', 'uye_min', 'uye_max', 'uye_kutu_sayisi', 'havuz',
              'havuz_ust', 'kat_havuz', 'kat_enkotu', 'enkotu_kutu',
              'uye_kapsam', 'uye_kapsam_pay', 'uye_alt_kapsayan',
              'kat_havuz_kapsayan', 'kat_enkotu_kapsayan'):
        d['numune_' + k] = nm.get(k, '')
    d['numune_urun_boylari'] = ';'.join('%s:%s' % (x[0], x[1]) for x in nm.get('urun_boylari', []))
    ru = c.get('ref_uye'); rr = c.get('ref_rakip')
    d['ref_uye'] = '%d/%d' % (ru['veren'], ru['toplam']) if ru else ''
    d['ref_rakip'] = '%d/%d' % (rr['veren'], rr['toplam']) if rr else ''
    kg = c.get('kuresel')
    if kg and 'urun' in kg:
        d['kuresel_urun'] = kg['urun']
        d['kuresel_boy'] = ';'.join('%s:%s' % (x[0], x[1]) for x in
                                    sorted(kg.get('boy', {}).items(), key=lambda y: -y[1])[:5])
    else:
        d['kuresel_urun'] = ''
        d['kuresel_boy'] = ''
    return d


SUTUNLAR = ['hedef', 'sira', 'numune_olcut', 'F', 'R', 'urun', 'cift_urun_sinifi', 'arms',
            'izgara_hucresi', 'sikilik',
            'numune_uye_alt', 'numune_uye_min', 'numune_uye_max', 'numune_uye_kutu_sayisi',
            'numune_uye_kapsam_pay', 'numune_havuz', 'numune_kat_havuz',
            'numune_kat_enkotu', 'numune_kat_havuz_kapsayan',
            'numune_kat_enkotu_kapsayan', 'numune_enkotu_kutu',
            'numune_olcut_2', 'numune_kat_enkotu_mm3', 'numune_kat_havuz_mm3',
            'numune_uye_kapsam_mm3',
            'ref_uye', 'ref_rakip', 'kuresel_urun', 'kuresel_boy',
            'cift_dTm', 'cift_het_tm', 'cift_het_dg', 'cift_uc_dg',
            'cift_Ta_kural', 'cift_Ta60_marj', 'cift_Ta60_uygun',
            'F_uz', 'F_tm', 'F_gc', 'F_uc', 'F_son5', 'F_hp_tm', 'F_hp_dg', 'F_hd_tm', 'F_hd_dg',
            'R_uz', 'R_tm', 'R_gc', 'R_uc', 'R_son5', 'R_hp_tm', 'R_hp_dg', 'R_hd_tm', 'R_hd_dg',
            'numune_urun_boylari']


def uret(sonuclar, panel, panel_yolu, yaz):
    os.makedirs(C.CIKTI, exist_ok=True)
    yollar = []
    yollar.append(_adaylar_tsv(sonuclar))
    yollar.append(_izgara_tsv(sonuclar))
    yollar.append(_rapor_md(sonuclar, panel, panel_yolu))
    if yaz:
        yaz('')
    return yollar


def _adaylar_tsv(sonuclar):
    p = os.path.join(C.CIKTI, 'adaylar.tsv')
    with open(p, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(SUTUNLAR)
        for s in sonuclar:
            for i, a in enumerate(s.get('adaylar', []), 1):
                a = dict(a); a['hedef'] = s['hedef']; a['sira'] = i
                w.writerow([a.get(k, '') for k in SUTUNLAR])
    return p


def _izgara_tsv(sonuclar):
    p = os.path.join(C.CIKTI, 'parametre_izgarasi.tsv')
    with open(p, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'GC_alt', 'GC_ust', 'Tm_alt', 'Tm_ust',
                    'urun_alt', 'urun_ust', "3'_uc_GC_sart", 'son5_GC_<=3_sart',
                    'gevseklik_puani', 'hayatta_kalan_aday', 'ornek_cift'])
        for s in sonuclar:
            for x in s.get('izgara', []):
                h = x['hucre']
                w.writerow([s['hedef'], h['gc'][0], h['gc'][1], h['tm'][0], h['tm'][1],
                            h['urun'][0], h['urun'][1],
                            'EVET' if h['uc_gc'] else 'HAYIR',
                            'EVET' if h['son5'] else 'HAYIR',
                            x['sikilik'], x['hayatta'], x['ornek']])
    return p


def _f(v, k=1):
    try:
        return ('%.' + str(k) + 'f') % float(v)
    except Exception:
        return str(v)


def _taban(s):
    t = s.get('panel_olcum') or {}
    return t.get('kat_enkotu'), t.get('kat_havuz')


def _hedef_karari(s):
    """'Is there a solution, and if so what does it cost' - it produces a one paragraph answer."""
    # The decision has THREE steps and their order matters. THE BASELINE is written
    # first: the value of the panel's current pair measured with the same engine on the
    # same bins; a candidate can count as "better" only if it beats that. The second
    # step is whether the candidate beat the current pair. The third is whether, even
    # having beaten it, it passes the 10x threshold: "better" and "enough" are not the
    # same thing, and a candidate that rises to 1.04x from a baseline of 0.19x is
    # better and still unusable.
    #
    # The COST list is deliberately tied to the decision sentence: even if a candidate
    # passes the threshold, if its product length goes over 250 bp, if it cannot be run
    # at 60 C, if it needs a deliberate mismatch (ARMS), or if it survives only in a
    # relaxed grid cell, that has to stand beside the number. Otherwise the report says
    # "there is a solution" and the order is placed without the cost being read.
    ad = s.get('adaylar', [])
    if s.get('durum') != 'TAMAMLANDI':
        return s.get('durum', 'bilinmiyor'), ''
    if not ad:
        return 'COZUM YOK', 'No candidate passed the in-sample measurement.'
    izg = s.get('izgara', [])
    dolu = [x for x in izg if x['hayatta'] > 0]
    en_siki = dolu[0] if dolu else None

    olculu = [a for a in ad if a.get('numune_kat_enkotu') not in ('', None)]
    olculu.sort(key=lambda a: -(a['numune_kat_enkotu'] or 0))
    if not olculu:
        return 'COZUM YOK', 'Not one of the candidates has a measurable separation.'
    iyi = olculu[0]
    arms_li = [a for a in olculu if a.get('arms')]
    duz_li = [a for a in olculu if not a.get('arms')]

    p = []
    tb_kotu, tb_havuz = _taban(s)
    if tb_kotu is not None or tb_havuz is not None:
        p.append(u'THE BASELINE, the panel\'s current pair under the same engine: discrimination %sx (the worst bin) / %sx (the pool). The numbers below are to be compared against it.'
                 % (_f(tb_kotu), _f(tb_havuz)))
    if en_siki:
        p.append(u'Under the strictest setting (%s), %d candidates survive.' % (en_siki['ad'], en_siki['hayatta']))
    bos_siki = [x for x in izg if x['sikilik'] == 0 and x['hayatta'] == 0]
    if bos_siki:
        p.append(u'Under the FULLY STRICT setting (GC 40-60, Tm 58-62, product 60-150, a G or C required at the 3\' end, <=3 G or C in the last 5 bases) THERE IS NO SOLUTION.')
    if duz_li:
        b = duz_li[0]
        p.append(u'The best candidate without ARMS is %s / %s (%s bp), discrimination %sx (the worst bin), member %%%s-%%%s, grid cell: %s.'
                 % (b['F'], b['R'], b['urun'], _f(b['numune_kat_enkotu']),
                    _f(b['numune_uye_min']), _f(b['numune_uye_max']), b['izgara_hucresi']))
    if arms_li:
        b = arms_li[0]
        p.append('ARMS varyantiyla en iyi: %s / %s (%s bp), ayrim %sx  [%s].'
                 % (b['F'], b['R'], b['urun'], _f(b['numune_kat_enkotu']), b['arms']))
    if iyi.get('kuresel_urun') not in ('', None):
        p.append(u'The product count of the best candidate in the global scan: %s.' % iyi['kuresel_urun'])
    # bedel
    bedeller = []
    if iyi.get('cift_urun_sinifi', '').startswith('kabul'):
        bedeller.append(u'a product of 150-250 bp -> a 30 s annealing and extension in the protocol')
    if iyi.get('cift_urun_sinifi', '').startswith('ONERILMEZ'):
        bedeller.append(u'a product above 250 bp -> not recommended for QuantiNova')
    if iyi.get('cift_Ta60_uygun') in (False, 'False'):
        bedeller.append(u'under the Ta = min(Tm)-3 rule this comes out at %s C, below the 60 C aim'
                        % _f(iyi.get('cift_Ta_kural')))
    if iyi.get('arms'):
        bedeller.append(u'a deliberate mismatch (ARMS) is needed, a separate meeting item')
    if 'serbest' in str(iyi.get('izgara_hucresi', '')):
        bedeller.append(u'the geometry rule has to be relaxed: %s' % iyi['izgara_hucresi'])
    if bedeller:
        p.append('BEDELI: ' + '; '.join(bedeller) + '.')
    en = iyi.get('numune_kat_enkotu') or 0
    if tb_kotu and en <= tb_kotu:
        karar = 'NO SOLUTION (the existing pair is better)'
        p.append(u'RESULT: not one candidate scanned BEAT the current pair, so the current pair should be kept.')
    elif en >= 10:
        karar = 'COZUM VAR'
    else:
        karar = 'KISMI COZUM'
    return karar, ' '.join(p)


def _rapor_md(sonuclar, panel, panel_yolu):
    p = os.path.join(C.CIKTI, 'KAPSAMLI_ARAMA_RAPORU.md')
    L = []
    A = L.append
    A(u'# The comprehensive primer search, the report')
    A('')
    A(u'Generated: %s' % time.strftime('%Y-%m-%d %H:%M:%S'))
    A('')
    A(u'Source panel: `%s`' % os.path.basename(panel_yolu))
    A('')
    A(u'## How to read it')
    A('')
    A(u'The question asked for every target: **under which parameter setting is there a solution, is there none, and if there is, what does it cost?** The `Karar` column in the summary table answers it in one word, and the paragraph under it in detail.')
    A('')
    A('Sabit qPCR kisitlari (QIAGEN Rotor-Gene Q + QuantiNova SYBR Green):')
    A('')
    A(u'- An amplicon of **60-150 bp is preferred**; 150-250 is acceptable but needs a **30 s** annealing and extension in the protocol; **above 250 is not recommended**.')
    A(u'- The Rotor-Gene runs a single cycling program per run: the aim is for the whole panel to run at **the same Ta**, with **60 C preferred**. Every candidate was judged at 60 C as well (the `cift_Ta60_marj` and `cift_Ta60_uygun` columns).')
    A(u'- Because this is SYBR Green, **primer dimer and hairpin are an ELIMINATING criterion**, not a warning: a candidate whose hairpin, homodimer or heterodimer Tm is >= 45 C, or whose dG is < -9 kcal/mol, is eliminated.')
    A('')
    A(u'> **The 60 C warning (a structural consequence that comes out of the measurement).** The panel\'s rule is `Ta = min(Tm) - 3`. Under that rule, for Ta to be 60 C the primer with the lower Tm has to be at **63 C**. Since the widest Tm window in the grid is 56-64, a shared Ta of 60 C is possible only with candidates **inside the Tm 56-64 window and with Tm >= 63**. Whether a candidate suits 60 C is marked separately in the `cift_Ta60_uygun` column; for those that do not, the choice is either to move the Tm window up or to relax the `Ta = min(Tm) - 3` rule for this panel. That is a **meeting decision**; the tool does not choose by itself.')
    A('')
    A(u'> **About ARMS.** A deliberate mismatch **is not a degenerate base**: it is one defined base, one oligo stays in the tube, the synthesis cost does not rise and it does not spoil the "no degenerate base in the panel" record. But it **does not match the template exactly**: it lowers the yield and it is a separate meeting item. The report gives the best candidate with ARMS and without ARMS **separately**, so that the decision stays with the user.')
    A('')

    A(u'> **About the baseline values.** For every target the CURRENT pair in the panel is measured again with the same engine, on the same bins, under the same criterion, and the candidates are compared **against that baseline**. On some targets that baseline departs from the number the panel published (a different read depth, a different subset of member bins, or a different relaxation setting). Where the departure is large it is written under the target as a **WARNING**. The comparison still holds: candidate and baseline are measured under **the same** conditions. But the panel\'s number and the number in this report **must not be compared directly**.')
    A('')
    A(u'## Summary')
    A('')
    A(u'The letter codes in the `Neden sorunlu` column: **G** a geometry violation, **K** conditional or a preliminary decision, **A** discrimination or coverage below the threshold, **U** the product length is outside the qPCR ideal, **C** removed from the panel, **P** it cannot be separated on a gel within the plate.')
    A('')
    A(u'| Target | Why it is a problem | Decision | Current pair (x) | Best candidate (x) | Product (bp) | Was ARMS needed |')
    A('|---|---|---|---|---|---|---|')
    for s in sonuclar:
        kar, _ = _hedef_karari(s)
        ad = s.get('adaylar', [])
        olculu = sorted([a for a in ad if a.get('numune_kat_enkotu') not in ('', None)],
                        key=lambda a: -(a['numune_kat_enkotu'] or 0))
        b = olculu[0] if olculu else {}
        tk, _th = _taban(s)
        A('| %s | %s | **%s** | %s | %s | %s | %s |' % (
            s['hedef'], s.get('etiketler', ''), kar, _f(tk) if tk else '-',
            _f(b.get('numune_kat_enkotu', '-')), b.get('urun', '-'),
            'EVET' if b.get('arms') else 'hayir'))
    A('')

    for s in sonuclar:
        A('---')
        A('')
        A('## %s' % s['hedef'])
        A('')
        if s.get('durum') != 'TAMAMLANDI':
            A('Durum: **%s**' % s.get('durum'))
            if s.get('hata'):
                A('')
                A('```')
                A(s['hata'][-1200:])
                A('```')
            A('')
            continue
        pn = s.get('panel', {})
        A(u'**The pair in the panel:** `%s` / `%s`, %s bp, plate %s, Ta %s'
          % (pn.get('F'), pn.get('R'), pn.get('urun'), pn.get('plaka'), pn.get('ta')))
        A('')
        t = s.get('panel_olcum') or {}
        if t:
            A('')
            A(u'**The current pair\'s values measured with THE SAME engine (the comparison baseline):** member %%%s-%%%s (%s bins), competitor pool %s, **discrimination %sx (pool) / %sx (worst bin: %s)**'
              % (_f(t.get('uye_min')), _f(t.get('uye_max')), t.get('uye_kutu_sayisi'),
                 t.get('havuz'), _f(t.get('kat_havuz')), _f(t.get('kat_enkotu')),
                 t.get('enkotu_kutu')))
        if s.get('uyelik_uyarisi'):
            A('')
            A(u'> **WARNING, the membership definition has to be checked.**')
            for u in s['uyelik_uyarisi']:
                A('> ' + u.replace('!! ', ''))
        A('')
        A('**Neden arandi:** ' + '; '.join(s.get('gerekceler', [])))
        A('')
        A(u'**The backbone:** `%s` (%s bp), the membership source is `%s`, member taxids: %s'
          % (s['omurga']['kutu'], s['omurga']['uzunluk'], s.get('uyelik_kaynagi'),
             ', '.join(s.get('uye_tax', []))))
        A('')
        sy = s.get('sayilar', {})
        A(u'**The length of the search:** %s windows -> %s forward + %s reverse candidates -> %s pairs (+%s ARMS variants) -> %s measured in the sample -> %s passed the pair structure.'
          % (sy.get('pencere'), sy.get('ileri'), sy.get('geri'), sy.get('cift'),
             sy.get('arms'), sy.get('numune_olculen'), sy.get('cift_yapisi_gecen')))
        A('')
        kar, aciklama = _hedef_karari(s)
        A(u'### The decision: %s' % kar)
        A('')
        A(aciklama)
        A('')
        A(u'### The parameter grid, how many candidates each setting leaves')
        A('')
        A(u'| GC | Tm | Product | 3\' end G or C | G or C in the last 5 | relaxation | survivors |')
        A('|---|---|---|---|---|---|---|')
        for x in s.get('izgara', [])[:24]:
            h = x['hucre']
            A('| %d-%d | %d-%d | %d-%d | %s | %s | %d | **%d** |' % (
                h['gc'][0], h['gc'][1], h['tm'][0], h['tm'][1], h['urun'][0], h['urun'][1],
                'sart' if h['uc_gc'] else 'serbest', '<=3' if h['son5'] else 'serbest',
                x['sikilik'], x['hayatta']))
        A('')
        A('(Tam 144 hucre: `parametre_izgarasi.tsv`)')
        A('')
        A(u'### The best candidates')
        A('')
        A(u'| # | Forward | Reverse | bp | member coverage | discrimination x (worst bin) | pool x | member % | ref member | ref competitor | global product | ARMS | grid cell |')
        A('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
        ad = sorted([a for a in s.get('adaylar', [])],
                    key=lambda a: -(a.get('numune_kat_enkotu') or 0))
        for i, a in enumerate(ad[:15], 1):
            A('| %d | `%s` | `%s` | %s | %s | %s | %s | %s-%s | %s | %s | %s | %s | %s |' % (
                i, a['F'], a['R'], a['urun'], a.get('numune_uye_kapsam_pay', ''),
                _f(a.get('numune_kat_enkotu', '')), _f(a.get('numune_kat_havuz', '')),
                _f(a.get('numune_uye_min', '')), _f(a.get('numune_uye_max', '')),
                a.get('ref_uye', ''), a.get('ref_rakip', ''), a.get('kuresel_urun', ''),
                a.get('arms', '') or '-', a.get('izgara_hucresi', '')))
        A('')
        A(u'Every column (Tm, GC, hairpin and dimer dG, the 60 C margin, the product length distribution) is in `adaylar.tsv`')
        A('')

    A('---')
    A('')
    A('## Yontem ve sinirlar')
    A('')
    A(u'- The measurement engine **was not rewritten**: `engine/ispcr.py` (`find_sites` and `amplify`), `engine/scanner.py` (`Havuz`) and `engine/pair.py` (`urunler`) are imported directly. The geometry thresholds are tested on every run to be identical to `engine/geometry_core.py` (see the self test).')
    A(u'- **The criterion label is written on every row** (the `numune_olcut` column). The elimination is done with **<=1 mismatch** (the panel\'s sample criterion); the best candidates are also measured with **<=3** and written into the `numune_olcut_2` and `*_mm3` columns. The two criteria are separate and cannot stand in for one another.')
    A(u'- The sample criterion is the same as the panel\'s: **<=1 mismatch plus an EXACT last 2 bases at the 3\' end**. The global criterion: **<=5 mismatches in total**, F and R separately. The two criteria are separate.')
    A(u'- The discrimination ratios are conservative through **Wilson**: the LOWER bound for a member and the UPPER bound for a competitor. The ratios of different targets are measured at different read depths and **cannot be compared directly** (the same warning as the panel\'s A26).')
    A(u'- On the completeness of the search: the backbone is scanned on ONE strand; since a pair is fully defined on a double stranded template (F on the + strand, R on the -), the reverse strand gives the same set. Some of the consensuses being stored in reverse does not affect the coverage.')
    A(u'- The funnel: window -> geometry -> pair -> sample -> reference -> **global** (the most expensive step last, and only on the candidates that passed every other filter).')
    A(u'- This tool **does not decide**: it measures, it writes down the cost, and it leaves the choice to the user.')
    A('')
    with open(p, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(L))
    return p
