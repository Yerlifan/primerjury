# -*- coding: utf-8 -*-
"""screening/taxonomy.py'nin KANITI - bes baslik bicimi.

NE SINANIYOR
------------
Katman 2'ye taksonomik ayrim eklemenin onundeki gercek engel, REFERANS_DB'deki
dosyalarin BES AYRI baslik bicimi kullanmasidir. Tek bir ayristirici varsaymak
kayitlarin cogunu "farkli alan" saydirir - yani duzeltmeye calistigimiz hatanin
daha buyugunu uretir. Bu sinama her bicimin DOGRU cozuldugunu gosterir.

Iki bolum var:
  1) SENTETIK - her bicimden birer baslik, beklenen alan ve klad uyelikleri
  2) GERCEK   - REFERANS_DB varsa her dosyadan ilk N baslik cozulur ve
                'alan cozulemedi' orani olculur. Yuksek oran, ayristiricinin
                o bicimi tanimadigini gosterir.

KOSMA
-----
    python3 tests/test_taxonomy.py
    python3 tests/test_taxonomy.py --ornek 5000
"""
from __future__ import print_function
import argparse
import io
import os
import sys

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, KOK)
from screening import taxonomy as T           # noqa: E402

# --- 1) SENTETIK: (ad, baslik, vtb, beklenen_alan, iceride_olmali_klad) -----
SENTETIK = [
    ('SILVA SSU',
     'AY846379.1.1791 Eukaryota;Archaeplastida;Chloroplastida;Chlorophyta;'
     'Chlorophyceae;Sphaeropleales;Monoraphidium;Monoraphidium sp. Itas 9/21 14-6w',
     'SILVA_138.2_SSURef_NR99.fasta', 'Eukaryota', 'Chlorophyta'),

    ('SILVA LSU',
     'AC152122.8318.11813 Eukaryota;Amorphea;Obazoa;Opisthokonta;Nucletmycea;'
     'Fungi;Dikarya;Ascomycota;Saccharomycotina;Saccharomycetes',
     'SILVA_138.2_LSURef_NR99.fasta', 'Eukaryota', 'Ascomycota'),

    ('UNITE',
     'UDB016649|k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Thelephorales;'
     'f__Thelephoraceae;g__Thelephora;s__Thelephora_albomarginata|SH1281904.10FU',
     'UNITE_ITS.fasta', 'Eukaryota', 'Basidiomycota'),

    ('PR2',
     'AB353770.1.1740_U|18S_rRNA|nucleus||Eukaryota|TSAR|Alveolata|'
     'Dinoflagellata|Dinophyceae|Peridiniales|Kryptoperidiniaceae|Unruhdinium',
     'PR2_SSU_taxo_long.fasta', 'Eukaryota', 'Alveolata'),

    ('ROD',
     'GCA_000001215|AE014298.5/23211192-23217141|Eukaryota;Opisthokonta;Metazoa;'
     'Arthropoda;Insecta;Diptera;Drosophilidae;Drosophila;Drosophila_melanogaster',
     'ROD_v1.2_operon_variants.fasta', 'Eukaryota', 'Arthropoda'),

    ('RefSeq bakteri (taksonomi YOK - alan dosyadan)',
     'NR_201932.1 Sphingosinicella wutangchuni strain LY54 16S ribosomal RNA, partial sequence',
     'bacteria.16S.fna', 'Bacteria', 'Sphingosinicella'),

    ('RefSeq arke (taksonomi YOK - alan dosyadan)',
     'NR_201921.1 Methanothermococcus jasoni strain Ax23 16S ribosomal RNA, complete sequence',
     'archaea.16S.fna', 'Archaea', 'Methanothermococcus'),

    ('RefSeq mantar (taksonomi YOK - alan dosyadan)',
     'NR_202962.1 Tremella indecorata ITS region; from TYPE material',
     'fungi.ITS.fna', 'Eukaryota', 'Tremella'),
]

# --- 2) SINIFLANDIRMA: hedef klad/alan verilince a/ao/b/c dogru mu ---------
SINIF = [
    # (ad, baslik, vtb, hedef_klad, hedef_alan, beklenen_sinif)
    ('hedefin KENDI kladi -> a',
     'AY846379.1.1791 Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;'
     'Marinifilaceae;Petrimonas;Petrimonas sulfuriphila',
     'SILVA_138.2_SSURef_NR99.fasta', ['Petrimonas'], 'Bacteria', 'a'),

    ('ayni alan, klad DISI -> b',
     'AY000001.1.1500 Bacteria;Pseudomonadota;Gammaproteobacteria;'
     'Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli',
     'SILVA_138.2_SSURef_NR99.fasta', ['Petrimonas'], 'Bacteria', 'b'),

    ('FARKLI alan -> c',
     'AY000002.1.1700 Eukaryota;Amorphea;Obazoa;Opisthokonta;Nucletmycea;Fungi',
     'SILVA_138.2_SSURef_NR99.fasta', ['Petrimonas'], 'Bacteria', 'c'),

    ('hedef alan icinde ORGANEL -> ao',
     'AY000003.1.1600 Bacteria;Cyanobacteriota;Chloroplast;Streptophyta',
     'SILVA_138.2_SSURef_NR99.fasta', ['Chloroplast'], 'Bacteria', 'ao'),

    ('alan cozulemez -> bilinmiyor (KANIT DEGIL)',
     'bilinmeyen_bir_baslik_taksonomisiz',
     'gizemli_dosya.fasta', ['Petrimonas'], 'Bacteria', 'bilinmiyor'),

    ('RefSeq, hedef cinsi adda geciyor -> a',
     'NR_201921.1 Methanothrix soehngenii strain X 16S ribosomal RNA',
     'archaea.16S.fna', ['Methanothrix'], 'Archaea', 'a'),

    # A2 duzeltmesi (2026-08-21): taksonomi TASIMAYAN kayitta, CINS USTU bir
    # hedef klad organizma adinda GECMEZ. Eskiden 'b' (sahte capraz) donuyordu.
    ('RefSeq + CINS USTU hedef -> bilinmiyor (b DEGIL)',
     'NR_201932.1 Bacteroides fragilis strain X 16S ribosomal RNA',
     'bacteria.16S.fna', ['Bacteroidales'], 'Bacteria', 'bilinmiyor'),

    ('RefSeq + FARKLI alan -> c (taksonomi gerekmez, vtb tanimi yeter)',
     'NR_201921.1 Methanothermococcus jasoni strain Ax23 16S ribosomal RNA',
     'archaea.16S.fna', ['Bacteroidales'], 'Bacteria', 'c'),

    ('SILVA + CINS USTU hedef -> b (taksonomi VAR, karar kesin)',
     'AY000001.1.1500 Bacteria;Pseudomonadota;Gammaproteobacteria;'
     'Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli',
     'SILVA_138.2_SSURef_NR99.fasta', ['Bacteroidales'], 'Bacteria', 'b'),
]


def bolum1():
    print('=== 1) BASLIK COZUMU (bes bicim) ===')
    print('%-46s | %-10s | %-10s | %s' % ('bicim', 'alan', 'beklenen', 'klad uyeligi'))
    print('-' * 100)
    ok = True
    for ad, baslik, vtb, bek_alan, bek_jeton in SENTETIK:
        alan, jet, _org, _tv = T.coz(baslik, vtb)
        jeton_var = bek_jeton in jet or any(bek_jeton in j for j in jet)
        gecti = (alan == bek_alan) and jeton_var
        ok = ok and gecti
        print('%-46s | %-10s | %-10s | %-5s %s'
              % (ad[:46], alan, bek_alan,
                 'VAR' if jeton_var else 'YOK',
                 '<-- DOGRU' if gecti else '<-- YANLIS (%s)' % bek_jeton))
    return ok


def bolum2():
    print()
    print('=== 2) SINIFLANDIRMA (a / ao / b / c / bilinmiyor) ===')
    print('%-46s | %-12s | %-12s | %s' % ('durum', 'sonuc', 'beklenen', ''))
    print('-' * 100)
    ok = True
    for ad, baslik, vtb, klad, alan, bek in SINIF:
        s = T.sinifla(baslik, vtb, klad, alan)
        gecti = (s == bek)
        ok = ok and gecti
        print('%-46s | %-12s | %-12s | %s'
              % (ad[:46], s, bek, '<-- DOGRU' if gecti else '<-- YANLIS'))
    return ok


def bolum3(ornek):
    """GERCEK dosyalar: cozulemeyen baslik orani olculur."""
    refdb = os.path.join(KOK, 'REFERANS_DB')
    if not os.path.isdir(refdb):
        print('\n(REFERANS_DB yok - gercek dosya bolumu atlandi)')
        return True
    print()
    print('=== 3) GERCEK DOSYALAR (ilk %d baslik) ===' % ornek)
    print('%-40s | %8s | %8s | %s' % ('dosya', 'cozuldu', 'cozulmedi', 'baskin alan'))
    print('-' * 100)
    ok = True
    dosyalar = [
        'SILVA_138.2_SSURef_NR99.fasta', 'SILVA_138.2_LSURef_NR99.fasta',
        'UNITE_ITS.fasta', 'PR2_SSU_taxo_long.fasta',
        'bacteria.16S.fna', 'archaea.16S.fna', 'fungi.ITS.fna',
        'fungi.18SrRNA.fna', 'fungi.28SrRNA.fna',
        'ROD_v1.2_operon_variants.fasta',
    ]
    for d in dosyalar:
        y = os.path.join(refdb, d)
        if not os.path.exists(y):
            print('%-40s | (dosya yok)' % d)
            continue
        say = {}
        n = 0
        with io.open(y, encoding='utf-8', errors='ignore') as fh:
            for satir in fh:
                if not satir.startswith('>'):
                    continue
                alan, _jet, _o, _tv = T.coz(satir.rstrip(), d)
                say[alan] = say.get(alan, 0) + 1
                n += 1
                if n >= ornek:
                    break
        cozulmedi = say.get('?', 0)
        cozuldu = n - cozulmedi
        baskin = sorted(((v, k) for k, v in say.items() if k != '?'), reverse=True)
        oran = 100.0 * cozulmedi / max(n, 1)
        # %5'ten fazla cozulememe, ayristiricinin o bicimi tanimadigini gosterir
        gecti = oran <= 5.0
        ok = ok and gecti
        print('%-40s | %8d | %8d | %s %s'
              % (d, cozuldu, cozulmedi,
                 ('%s %%%.1f' % (baskin[0][1], 100.0 * baskin[0][0] / max(n, 1))) if baskin else '-',
                 '' if gecti else '<-- COZULEMEYEN ORANI YUKSEK (%%%.1f)' % oran))
    return ok


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--sample', '--ornek', dest='ornek', type=int, default=2000)
    a = p.parse_args()
    s = [bolum1(), bolum2(), bolum3(a.ornek)]
    print()
    print('SONUC: ' + ('UC BOLUMUN UCU DE GECTI' if all(s) else 'BASARISIZ'))
    return 0 if all(s) else 1


if __name__ == '__main__':
    sys.exit(main())
