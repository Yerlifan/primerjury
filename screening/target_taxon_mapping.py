# -*- coding: utf-8 -*-
"""target_taxon_mapping.py produces two tables:
  Table A  the targets in the panel -> which decision each one comes from
  Table B  the taxa in the sample -> which target covers each one, measured

Output: a markdown mapping table, plus the sheet "18 The target to taxon
mapping" inside the panel workbook

"""
# -------------------------------------------------------------------------
# target_taxon_mapping.py ties two questions to a table: which meeting decision
#                          each target in the panel comes from, and which target
#                          covers each taxon in the sample.
#
# INPUT  : the measurement result pattern given with --r2, the json cross_coverage.py
#          produces given with --cross, the pair table with --pairs,
#          taxid_names.tsv with --taxid and the delivery panel with --xlsx. The
#          target to decision mapping (the KARAR dictionary) is fixed inside the
#          file; the rows that do not come from a meeting decision are marked there
#          plainly as "Karar 5 - derived from the measurement".
# OUTPUT : the Markdown file given with --md and the sheet
#          the sheet "18 The target to taxon mapping" added to the --xlsx file.
# CALLED BY: IT IS NOT IN THE MENU, it is run by hand; it changes the delivery
#          xlsx. cross_coverage.py, one of its inputs, is run by hand as well.
#
# The value of Table B is that it makes visible the taxa no specific target covers:
# if a taxon amplifies only with a universal or control primer, there is no specific
# measurement in the panel for that organism.
# -------------------------------------------------------------------------
import sys, os, json, glob, csv, argparse
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment

KIRMIZI = PatternFill('solid', fgColor='FFC7CE')
SARI    = PatternFill('solid', fgColor='FFEB9C')
GRI     = PatternFill('solid', fgColor='D9D9D9')
YESIL   = PatternFill('solid', fgColor='C6EFCE')
KALIN   = Font(bold=True)
SAR     = Alignment(wrap_text=True, vertical='top')

# the panel row -> (target name, level, decision, an explanation of the source)
KARAR = {
 2:  ('Metanojen_universal',              'grup (islevsel)', 'Decision 4, domain and '
                                                             'universal',
      'Asked for: a universal methanogen control primer.'),
 3:  ('Methanothrix_cinsi',               'cins',            'Decision 5, derived from a '
                                                             'measurement',
      'NOT in the decisions. It was designed in this session; genus level was '
      'put in when the M. soehngenii species claim was dropped.'),
 4:  ('Petrimonas_cinsi',                 'cins',            'Decision 2, four genus '
                                                             'specific',
      'Asked for: one of the four genera.'),
 5:  ('Metanomikrobiyales_hidrojenotrof', 'takim',           'Decision 3, a functional or '
                                                             'ecological group',
      'Asked for: hydrogenotrophic methanogens. The name was narrowed to the '
      'order Methanomicrobiales according to the measured coverage.'),
 6:  ('Mantar_universal (F1)',            'alan',            'Decision 4, domain and '
                                                             'universal',
      'Asked for: a universal fungal control primer, the F1 barcode class.'),
 7:  ('Methanosarcina_cinsi',             'cins',            'Decision 5, derived from a '
                                                             'measurement',
      'NOT directly in the decisions. It STANDS IN FOR the species specific '
      'M. barkeri target of decision 1: the species is not in the sample, '
      'since the 2208 bin is 97.4 to 97.9 per cent to M. vacuolata, so genus '
      'level was given.'),
 8:  ('Metilotrofik_metanojen',           'grup (islevsel)', 'Decision 3, a functional or '
                                                             'ecological group',
      'Asked for: methylotrophic methanogens.'),
 9:  ('Nitrosocosmicus_AOA',              'grup (ekolojik)', 'Decision 3, a functional or '
                                                             'ecological group',
      'Asked for: ammonia oxidising archaea. It was not added later.'),
 10: ('Proteiniphilum_cinsi',             'cins',            'Decision 2, four genus '
                                                             'specific',
      'Asked for: one of the four genera.'),
 11: ('Proteolitik_Synergistaceae',       'soy',             'Decision 5, derived from a '
                                                             'measurement, under decision 3',
      'What was asked for was proteolytic and syntrophic bacteria; the '
      'measured lineage came out Synergistaceae. Formerly named '
      'Proteolitik_Cloacibacillus.'),
 12: ('Bacteroidales_kumesi',             'soy',             'Decision 5, derived from a '
                                                             'measurement',
      'It STANDS IN FOR the Bacteroides and Alistipes targets of decision 2: '
      'neither genus is in the sample, the bins are 95.3 to 96.8 per cent '
      'similar to one another, and they are a single unnameable Bacteroidales '
      'lineage.'),
 13: ('Mantar_universal (F2)',            'alan',            'Decision 4, domain and '
                                                             'universal',
      'Asked for: a universal fungal control primer, the F2 barcode class.'),
 14: ('Microascaceae_askomikot',          'aile',            'Decision 5, derived from a '
                                                             'measurement',
      'A FAMILY LEVEL TARGET. It STANDS IN FOR the species specific '
      'Trichoderma asperellum target of decision 1: the organism in the '
      '101201 bin is not Trichoderma but Petriella and Microascaceae in ITS. '
      'TAKEN OUT OF THE PANEL, with a separation of 0.7x.'),
 15: ('Sakarolitik_Sphaerochaeta',        'cins',            'Decision 5, derived from a '
                                                             'measurement, under decision 3',
      'What was asked for was saccharolytic bacteria; the measured member '
      'came out Sphaerochaeta associata.'),
 16: ('Asetoklastik_metanojenler',        'grup (islevsel)', 'Decision 3, a functional or '
                                                             'ecological group',
      'Asked for: acetoclastic methanogens.'),
 17: ('Bakteri_universal',                'alan',            'Decision 4, domain and '
                                                             'universal',
      'Asked for: a universal bacterial control primer.'),
 18: ('Petriella_musispora',              'cins -> sinif',   'Decision 5, derived from a '
                                                             'measurement',
      'It came out of the measurement in place of the species specific '
      'Podospora pseudopauciseta target of decision 1. TAKEN OUT OF THE '
      'PANEL, with a separation of 0.7x.'),
 19: ('Proteolitik_Cloacimonas',          'cins',            'Decision 5, derived from a '
                                                             'measurement, under decision 3',
      'What was asked for was proteolytic and syntrophic bacteria; the '
      'measured member came out Ca. Cloacimonas acidaminovorans.'),
 20: ('Arke_universal',                   'alan',            'Decision 4, domain and '
                                                             'universal',
      'Asked for: a universal archaeal control primer.'),
 21: ('Methanothrix_soehngenii_turu',     'TUR (kosullu)',   'Decision 1, six species '
                                                             'specific',
      'Asked for: one of the six species. SPECIES LEVEL WAS GIVEN, '
      'conditionally.'),
 22: ('Methanosarcina_mazei_turu',        'TUR GRUBU',       'Decision 1, six species '
                                                             'specific',
      'Asked for: one of the six species. A SPECIES GROUP was given, because '
      'M. mazei and M. soligelidi do not separate. BELOW THE THRESHOLD after '
      'the read engine fix; see the read engine sheet.'),
}

# requests asked for under Decisions 1 and 2 that COULD NOT ENTER THE PANEL AS A TARGET
KARSILANMAYAN = [
 ('Decision 1, species specific', 'Methanosarcina barkeri',
  'The organism is not in the sample: the 2208 bin is 97.4 to 97.9 per cent '
  'to M. vacuolata. Row 7, Methanosarcina_cinsi, stands in for it.'),
 ('Decision 1, species specific', 'Podospora pseudopauciseta',
  'The organism is not in the sample. Row 18, Petriella_musispora, stands in '
  'for it, and that one was taken out of the panel too.'),
 ('Decision 1, species specific', 'Dictyostelium discoideum (44689)',
  'The Kraken2 label was refuted and the bin is heterogeneous. THERE IS NO '
  'PANEL TARGET for it; only the universal fungal F1 pair amplifies it.'),
 ('Decision 1, species specific', 'Trichoderma asperellum (101201)',
  'The organism in the bin is not Trichoderma; in ITS it is Petriella and '
  'Microascaceae. Row 14 stands in for it, and that one was taken out too.'),
 ('Decision 2, genus specific', 'Bacteroides',
  'The genus is not in the sample; the best match, Alistipes putredinis, is '
  '85.5 per cent. Row 12, Bacteroidales_kumesi, stands in for it.'),
 ('Decision 2, genus specific', 'Alistipes',
  'The same organism as Bacteroides; they were merged under row 12.'),
]

GENIS_AD = {2: 'Metanojen_universal', 6: 'Mantar_universal (F1)', 13: 'Mantar_universal (F2)',
            17: 'Bakteri_universal', 20: 'Arke_universal'}


def yukle(r2_desen, capraz, ciftler_tsv, taxid_tsv):
    R = {}
    for f in glob.glob(r2_desen):
        R.update(json.load(open(f)))
    CP = json.load(open(capraz))
    ad = {}
    for l in open(taxid_tsv, encoding='utf-8'):
        p = l.rstrip('\n').split('\t')
        if len(p) >= 2:
            ad[p[0]] = p[1]
    uye = {}
    for c in csv.DictReader(open(ciftler_tsv, encoding='utf-8'), delimiter='\t'):
        uye[int(c['satir'])] = set(x for x in c['uye_taksonlar'].split(',') if x)
    return R, CP, ad, uye


def takson_tablosu(R, CP, ad, uye):
    tax = {}
    for kaynak, veri in (('sinif', R), ('capraz', CP)):
        for k, v in veri.items():
            s, kutu = k.split('|')
            s = int(s)
            tx = kutu.split('_')[1]
            if kaynak == 'sinif':
                y1, n = v[1], v[4]
                y3 = v[3]
            else:
                y1, y3, n = v[0], v[1], v[2]
            if not n:
                continue
            d = tax.setdefault(tx, {})
            a = d.get(s, (0.0, 0.0))
            d[s] = (max(a[0], 100.0 * y1 / n), max(a[1], 100.0 * y3 / n))
    satirlar = []
    for tx in sorted(tax, key=lambda x: ad.get(x, 'zz')):
        d = tax[tx]
        u = sorted(KARAR[s][0] for s in d if tx in uye.get(s, set()))
        c1 = sorted(KARAR[s][0] for s in d if d[s][0] >= 10)
        c3 = sorted(KARAR[s][0] for s in d if d[s][1] >= 10 and d[s][0] < 10)
        ozgul1 = [x for x in c1 if x not in GENIS_AD.values()]
        if c1:
            durum = 'covered by a specific target' if ozgul1 else 'a universal or control primer only'
        elif c3:
            durum = 'amplified at mm<=3 only'
        else:
            durum = 'A GAP: no target amplifies it'
        satirlar.append(dict(taxid=tx, ad=ad.get(tx, '?'),
                             uye='; '.join(u) or '-',
                             cog1='; '.join(c1) or '-',
                             cog3='; '.join(c3) or '-', durum=durum))
    return satirlar


def md_yaz(yol, satirlar):
    L = []
    L.append(u'# The target to taxon mapping, 2 August 2026\n')
    L.append(u'The 21 targets in the panel are not the 44 taxa in the sample themselves. The targets are derived from the meeting decisions and from the measurement.\n')
    L.append(u'\n## Table A, the 21 targets in the panel by their source\n')
    L.append(u'| # | Target | Level | Source | Explanation |')
    L.append('|---|---|---|---|---|')
    for s in sorted(KARAR):
        h, dz, kr, ac = KARAR[s]
        L.append('| %d | %s | %s | %s | %s |' % (s, h, dz, kr, ac))
    L.append(u'\n### The count per decision\n')
    say = {}
    for s in KARAR:
        k = KARAR[s][2].split(' (')[0]
        say[k] = say.get(k, 0) + 1
    L.append(u'| Source | Targets | Which ones |')
    L.append('|---|---|---|')
    for k in sorted(say):
        h = [KARAR[s][0] for s in sorted(KARAR) if KARAR[s][2].split(' (')[0] == k]
        L.append('| %s | %d | %s |' % (k, say[k], ', '.join(h)))
    L.append(u'\n### Asked for in the meeting decisions but not admitted to the panel as a target\n')
    L.append(u'| Decision | Asked for | What happened |')
    L.append('|---|---|---|')
    for a, b, c in KARSILANMAYAN:
        L.append('| %s | %s | %s |' % (a, b, c))
    L.append(u'\n### The family level target\n')
    L.append(u'The one family level target: **Microascaceae_askomikot** (row 14, Decision 5). It was removed from the panel (discrimination 0,7x).\n')
    L.append(u'\n## Table B, the 44 taxa in the sample and which target covers them\n')
    L.append(u'The measurement: the corrected read engine, at most 3000 reads per bin, a threshold of >=%10 product. The five universal and broad targets were measured in all 99 bins with no class boundary.\n')
    L.append(u'| taxid | Taxon | The target(s) it is a member of | The targets that amplify it (mm<=1) | Extra (mm<=3 only) | State |')
    L.append('|---|---|---|---|---|---|')
    for r in satirlar:
        L.append('| %s | %s | %s | %s | %s | %s |' % (r['taxid'], r['ad'], r['uye'],
                                                      r['cog1'], r['cog3'], r['durum']))
    bos = [r for r in satirlar if r['durum'].startswith('BOSLUK')]
    yev = [r for r in satirlar if r['durum'] == 'a universal or control primer only']
    m3 = [r for r in satirlar if r['durum'] == 'amplified at mm<=3 only']
    L.append(u'\n### The gaps\n')
    L.append(u'| State | Count | Taxa |')
    L.append('|---|---|---|')
    L.append(u'| No target amplifies it | %d | %s |' % (len(bos), ', '.join(r['ad'] for r in bos) or '-'))
    L.append(u'| It amplifies only under the mm<=3 criterion | %d | %s |' % (len(m3), ', '.join(r['ad'] for r in m3) or '-'))
    L.append(u'| Only a universal or control primer amplifies it, it has no specific target | %d | %s |'
             % (len(yev), ', '.join(r['ad'] for r in yev) or '-'))
    L.append(u'| It is covered by a specific target | %d | - |'
             % len([r for r in satirlar if r['durum'].startswith('kapsanan')]))
    L.append(u'\nTaxa in total: **%d**\n' % len(satirlar))
    open(yol, 'w', encoding='utf-8').write('\n'.join(L) + '\n')


def xlsx_yaz(xlsx, satirlar):
    wb = openpyxl.load_workbook(xlsx)
    adx = '18 The target to taxon mapping'
    if adx in wb.sheetnames:
        del wb[adx]
    ws = wb.create_sheet(adx)
    for w, c in zip((8, 34, 20, 30, 44, 44, 26), 'ABCDEFG'):
        ws.column_dimensions[c].width = w
    n = 1

    def yaz(r, c, v, fill=None, bold=False):
        h = ws.cell(r, c, v)
        if fill: h.fill = fill
        if bold: h.font = KALIN
        h.alignment = SAR
        return h

    yaz(n, 1, u'TARGET TO TAXON MAPPING', bold=True); n += 1
    yaz(n, 1, u'The 21 targets in the panel are NOT the same thing as the 44 taxa in the sample. The targets are derived from the meeting decisions and from measurement.'); n += 2

    yaz(n, 1, u'TABLE A - the 21 targets in the panel, by where each came from', bold=True, fill=GRI); n += 1
    for j, h in enumerate(['#', 'Hedef', 'Duzey', 'Kaynak', 'Aciklama'], 1):
        yaz(n, j, h, bold=True, fill=GRI)
    n += 1
    for s in sorted(KARAR):
        h, dz, kr, ac = KARAR[s]
        f = SARI if kr.startswith('Karar 5') else None
        yaz(n, 1, s, fill=f); yaz(n, 2, h, fill=f); yaz(n, 3, dz, fill=f)
        yaz(n, 4, kr, fill=f); yaz(n, 5, ac, fill=f)
        ws.row_dimensions[n].height = 30
        n += 1
    n += 1

    yaz(n, 1, u'THE COUNT PER DECISION', bold=True, fill=GRI); n += 1
    say = {}
    for s in KARAR:
        k = KARAR[s][2].split(' (')[0]
        say.setdefault(k, []).append(KARAR[s][0])
    for k in sorted(say):
        yaz(n, 1, len(say[k]), bold=True); yaz(n, 2, k, bold=True)
        yaz(n, 4, ', '.join(say[k])); ws.merge_cells(start_row=n, start_column=4, end_row=n, end_column=5)
        ws.row_dimensions[n].height = 30
        n += 1
    n += 1

    yaz(n, 1, u'REQUESTED IN THE DECISIONS BUT COULD NOT ENTER THE PANEL AS A TARGET', bold=True, fill=KIRMIZI); n += 1
    for j, h in enumerate(['Karar', 'Istenen', 'Ne oldu'], 1):
        yaz(n, j + 1, h, bold=True, fill=GRI)
    n += 1
    for a, b, c in KARSILANMAYAN:
        yaz(n, 2, a); yaz(n, 3, b); yaz(n, 4, c)
        ws.merge_cells(start_row=n, start_column=4, end_row=n, end_column=5)
        ws.row_dimensions[n].height = 30
        n += 1
    n += 1

    yaz(n, 1, u'TABLE B - the 44 taxa in the sample, and which target covers each', bold=True, fill=GRI); n += 1
    yaz(n, 1, u'Measurement: the corrected read engine, <=3000 reads per bin, threshold >=%10 product. The five universal or broad targets were measured across all 99 bins WITH NO CLASS BOUNDARY (the panel\'s own measurements were class based).')
    ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=7)
    ws.row_dimensions[n].height = 30
    n += 1
    for j, h in enumerate(['taxid', 'Takson', 'The target or targets it is a '
                                              'member of', 'The targets that amplify it '
                                                                         '(mm<=1)',
                           'Extra, at mm<=3 only', 'Durum'], 1):
        yaz(n, j, h, bold=True, fill=GRI)
    n += 1
    for r in satirlar:
        f = (KIRMIZI if r['durum'].startswith('BOSLUK') else
             SARI if r['durum'] in ('a universal or control primer only', 'amplified at mm<=3 only') else None)
        yaz(n, 1, r['taxid'], fill=f); yaz(n, 2, r['ad'], fill=f); yaz(n, 3, r['uye'], fill=f)
        yaz(n, 4, r['cog1'], fill=f); yaz(n, 5, r['cog3'], fill=f); yaz(n, 6, r['durum'], fill=f)
        ws.row_dimensions[n].height = 26
        n += 1
    n += 1

    yaz(n, 1, 'BOSLUKLAR', bold=True, fill=KIRMIZI); n += 1
    gr = [('No target amplifies it', [r for r in satirlar if r['durum'].startswith('BOSLUK')], KIRMIZI),
          ('It is amplified at the mm<=3 criterion only', [r for r in satirlar if r['durum'] == 'amplified at mm<=3 only'], SARI),
          ('Only a universal or control primer amplifies it; it has no '
           'specific target',
           [r for r in satirlar if r['durum'] == 'a universal or control primer only'], SARI),
          ('It is covered by a specific target',
           [r for r in satirlar if r['durum'].startswith('kapsanan')], YESIL)]
    for ad2, lst, f in gr:
        yaz(n, 1, len(lst), bold=True, fill=f); yaz(n, 2, ad2, fill=f)
        yaz(n, 4, ', '.join(r['ad'] for r in lst) or '-', fill=f)
        ws.merge_cells(start_row=n, start_column=4, end_row=n, end_column=6)
        ws.row_dimensions[n].height = 30
        n += 1
    yaz(n, 1, u'TOTAL', bold=True); yaz(n, 2, len(satirlar), bold=True)
    wb.save(xlsx)
    return len(satirlar)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--r2', required=True)
    ap.add_argument('--cross', dest='capraz', required=True)
    ap.add_argument('--pairs', dest='ciftler', required=True)
    ap.add_argument('--taxid', required=True)
    ap.add_argument('--md', required=True)
    ap.add_argument('--xlsx', required=True)
    a = ap.parse_args()
    R, CP, ad, uye = yukle(a.r2, a.capraz, a.ciftler, a.taxid)
    satirlar = takson_tablosu(R, CP, ad, uye)
    md_yaz(a.md, satirlar)
    n = xlsx_yaz(a.xlsx, satirlar)
    print('md:', a.md)
    print(u'xlsx sheet added, taxon:', n)
    for r in satirlar:
        if not r['durum'].startswith('kapsanan'):
            print(' ', r['durum'], '|', r['ad'])


if __name__ == '__main__':
    main()
