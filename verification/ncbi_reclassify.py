# -*- coding: utf-8 -*-
"""RECLASSIFY THE NCBI RESULTS, from the saved pages, without going to the network.

WHY
---
Layer 4's first classifier looked for keywords to count a record as "unnamed"
(uncultured, clone, metagenome and so on). Once the run finished, a look at the
example headers showed that the records below are UNNAMED as well, and yet they had
been counted as named:

    Bacterium LC2012 16S ribosomal RNA gene
    Archaeon 2022-TM-MRBT1 gene for 16S rRNA
    anaerobic methanogenic archaeon E15-5 16S rRNA gene

"""
from __future__ import print_function

import argparse
import io
import os
import re
import sys
import time

# The words that are NOT a taxon name and mean "something". If a header starts with
# one of them the record is unnamed.
AD_DEGIL = {
    'bacterium', 'archaeon', 'organism', 'prokaryote', 'eukaryote',
    'uncultured', 'unidentified', 'unclassified', 'environmental',
    'anaerobic', 'aerobic', 'methanogenic', 'thermophilic', 'halophilic',
    'marine', 'soil', 'rumen', 'gut', 'sludge', 'compost', 'fungal', 'fungus',
    'yeast', 'synthetic', 'cloning', 'expression', 'mixed', 'candidate',
    'metagenome', 'mag:', 'mag', 'toluene-degrading', 'sulfate-reducing',
    'iron-reducing', 'nitrate-reducing', 'ammonia-oxidizing', 'symbiont',
    'endophyte', 'contaminant', 'unverified', 'predicted', 'putative',
}

# The keyword filter (the old rule); it is still applied, only it is no longer the
# ONLY criterion.
ADSIZ_IZ = ('uncultured', 'unidentified', 'unclassified', 'metagenome',
            'environmental', 'enrichment culture', 'clone', 'synthetic construct',
            'consortium', 'isolate ')

CINS = re.compile(r'^(?:Candidatus\s+)?(\[?[A-Z][a-z]{2,}\]?)\s+'
                  r'(?:cf\.\s+|aff\.\s+|sp\.|subsp\.|var\.|[a-z][a-z\-]{2,})')


def adli_mi(baslik):
    """Baslik gercek bir 'Cins tur' adiyla basliyor mu."""
    b = (baslik or '').strip()
    if not b:
        return False
    bl = b.lower()
    for iz in ADSIZ_IZ:
        if iz in bl:
            return False
    ilk = re.split(r'[\s,;:]+', b)[0].strip('[]').lower()
    if ilk in AD_DEGIL:
        return False
    m = CINS.match(b)
    if not m:
        return False
    if m.group(1).strip('[]').lower() in AD_DEGIL:
        return False
    return True


def urunler(html):
    d = re.sub(r'<[^>]+>', ' ', html)
    d = re.sub(r'&nbsp;?', ' ', d)
    d = re.sub(r'\s+', ' ', d)
    return [(m.group(1), m.group(2).strip())
            for m in re.finditer(
                r'>\s*([A-Z]{1,2}[_A-Z]*\d{5,}\.\d)\s+([^>]{5,200}?)\s+'
                r'product length\s*=\s*\d+', d)]



# --- BILINEN CEVAPLI SINAV -------------------------------------------------
# Kural degistiginde sessizce bozulmasin diye. Hepsi gercek NCBI basliklari;
# 2026-08-10 kosusunun ham sayfalarindan alindi.
SINAV = [
    # (baslik, ADLI mi)
    (u'Bacterium LC2012 16S ribosomal RNA gene, partial sequence', False),
    (u'Bacterium strain ASV8595 16S ribosomal RNA gene', False),
    (u'Archaeon 2022-TM-MRBT1 gene for 16S rRNA, partial sequence', False),
    (u'anaerobic methanogenic archaeon E15-5 16S rRNA gene, partial', False),
    (u'Environmental 16s rDNA sequence from Evry wastewater treatment plant', False),
    (u'Uncultured bacterium clone OTU639 16S ribosomal RNA gene', False),
    (u'Methanogenic prokaryote enrichment culture B31_1_13 16S ribosomal RNA', False),
    (u'Toluene-degrading methanogenic consortium archaeon M2 16S ribosomal RNA', False),
    (u'MAG: Methanothrix sp. isolate d1628e35 genome assembly', False),
    (u'Archaeoglobus fulgidus DSM 8774, complete genome', True),
    (u'Euplotes octocarinatus genome assembly, organelle: macronuclear', True),
    (u'Candidatus Methanocrinis natronophilus strain Mx 16S ribosomal RNA gene', True),
    (u'[Petriella] asymmetrica var. cypria culture CBS:258.31', True),
    (u'Methanosarcina sp. 1H1 gene for 16S ribosomal RNA, partial sequence', True),
    (u'Porphyromonas gingivalis strain UTM FZZ12 16S ribosomal RNA gene', True),
    (u'Albugo laibachii Nc14, genomic contig CONTIG_1449', True),
]


def sinav():
    gecti = dusen = 0
    for baslik, bek in SINAV:
        v = adli_mi(baslik)
        if v == bek:
            gecti += 1
        else:
            dusen += 1
            print('  DUSTU  bekleniyor %-5s cikan %-5s  %s'
                  % ('ADLI' if bek else 'ADSIZ', 'ADLI' if v else 'ADSIZ', baslik[:70]))
    print('  ad kurali sinavi: %d gecti, %d dustu' % (gecti, dusen))
    return 0 if dusen == 0 else 1


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', dest='kok', default='.')
    p.add_argument('--selftest', dest='sina', action='store_true',
                   help='only ad kurali sinavini kos, file reads')
    a = p.parse_args()
    if a.sina:
        return sinav()
    kok = os.path.abspath(a.kok)
    # Sinav HER kosuda once gecer. Kural bozuksa sayimlar da bozuktur.
    if sinav() != 0:
        print(u'  THE NAME RULE TEST FAILED, no classification was done.')
        return 1
    ham = os.path.join(kok, 'DOGRULAMA_SONUC', 'ncbi_ham')
    if not os.path.isdir(ham):
        sys.exit(u'ERROR: %s is missing.' % ham)

    eski = {}
    ey = os.path.join(kok, 'DOGRULAMA_SONUC', 'ncbi_katman4.tsv')
    if os.path.exists(ey):
        bas = None
        for l in io.open(ey, encoding='utf-8'):
            if l.startswith('#'):
                continue
            p_ = l.rstrip('\n').split('\t')
            if bas is None:
                bas = p_
                continue
            r = dict(zip(bas, p_))
            eski[r.get('hedef', '')] = r

    print('=' * 78)
    print(u'  RECLASSIFYING THE NCBI RESULTS (no network access)  %s'
          % time.strftime('%Y-%m-%d %H:%M'))
    print('=' * 78)
    print('  %-44s %8s %8s %8s' % ('hedef', 'eski adli', 'yeni adli', 'toplam'))
    satir = []
    for f in sorted(os.listdir(ham)):
        if not f.endswith('.html'):
            continue
        yol = os.path.join(ham, f)
        h = io.open(yol, encoding='utf-8', errors='replace').read()
        ur = urunler(h)
        if not ur:
            continue
        adli = [x for x in ur if adli_mi(x[1])]
        ad = f[:-5]
        # the mapping between the file name and the target name: the best match among the
        # targets in the TSV is chosen (the punctuation in the file name is simplified).
        hedef = ad
        for h2 in eski:
            if re.sub(r'\W+', '_', h2) == ad:
                hedef = h2
                break
        e = eski.get(hedef, {})
        try:
            eski_adli = int(e.get('adli_hedef_disi') or -1)
        except ValueError:
            eski_adli = -1
        print('  %-44s %8s %8d %8d'
              % (hedef[:44], eski_adli if eski_adli >= 0 else '-', len(adli), len(ur)))
        satir.append((hedef, eski_adli, len(adli), len(ur) - len(adli), len(ur),
                      [x[1] for x in adli[:5]]))

    cy = os.path.join(kok, 'DOGRULAMA_SONUC', 'ncbi_katman4_siki.tsv')
    with io.open(cy, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# NCBI layer 4, the STRICT name rule. Generated %s\n'
                 % time.strftime('%Y-%m-%d %H:%M'))
        fh.write(u'# A record counts as NAMED only if its header starts with a real\n# "Genus species" name. Headers such as "Bacterium LC2012", "Archaeon\n# 2022-TM-MRBT1", "anaerobic methanogenic archaeon E15-5" and\n# "Environmental 16s rDNA sequence ..." carry NO genus name; the old\n# loose rule counted all of them as named.\n')
        fh.write(u'hedef\tgevsek_kural_adli\tsiki_kural_adli\tsiki_adsiz\ttoplam\t'
                 u'ornek_adli_basliklar\n')
        for h, e, y, ad, t, orn in satir:
            fh.write(u'%s\t%s\t%d\t%d\t%d\t%s\n'
                     % (h, e if e >= 0 else '', y, ad, t,
                        ' | '.join(o[:70] for o in orn)))
    print()
    print(u'  written: %s' % cy)
    dus = [(h, e, y) for h, e, y, _a, _t, _o in satir if e >= 0 and y < e]
    print(u'  targets whose count DROPS under the strict rule: %d' % len(dus))
    for h, e, y in sorted(dus, key=lambda x: -(x[1] - x[2]))[:8]:
        print('    %-44s %d -> %d' % (h[:44], e, y))
    return 0


if __name__ == '__main__':
    sys.exit(main())
