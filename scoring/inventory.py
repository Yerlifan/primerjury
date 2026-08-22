# -*- coding: utf-8 -*-
"""THE PHYSICAL LOCUS INVENTORY: which loci are REALLY in the sample.

The sentence "there is no other locus" can only be established by measurement. Three
INDEPENDENT measures:

  1) THE RAW READ LENGTH DISTRIBUTION (unfiltered). If the longest read in a library
     is 1650 bp, then 23S and ITS are PHYSICALLY absent there, whatever anchor is
     found. This measure DOES NOT depend on the anchor choice, and it is the
     strongest evidence.
  2) THE CANONICAL CONSENSUS LENGTH.
  3) THE ANCHOR SCAN: conserved rRNA motifs (the SSU / 5.8S / LSU / ITS boundaries).
     The engine used is the project's OWN okuma_motoru.Sonda class (IUPAC aware plus
     pigeonhole seeding, lossless). The son2 rule is OFF here: an anchor search is a
     detection of PRESENCE, not a PCR binding.

Failing to find an anchor DOES NOT on its own mean "it is absent" (the anchor
sequence may have changed in that lineage). So the decision is made with measure 1
and the confirmation with measure 3.

"""
import os, sys, json, statistics, argparse, gzip
sys.path.insert(0, '/tmp/mrb')
from shared_scorer import Puanlayici

KOK = '/sessions/dreamy-elegant-wozniak/mnt/PrimerTasarlama'
KP = '/tmp/mrb/kontrol/envanter.json'

HEDEFLER = [
    'Petrimonas_cinsi',
    'Methanosarcina mazei / M. soligelidi grubu',
    'Bacteroidales_kumesi',
    'Petriella_musispora',
    'Proteiniphilum_cinsi',
    'Mantar_universal (F1)',
    'Bakteri_universal',
]

# The conserved rRNA anchors. Name -> (sequence, which part of the operon)
CAPALAR = [
    ('27F',     'AGAGTTTGATCMTGGCTCAG', 'SSU 5\' ucu (bakteri/arke)'),
    ('515F',    'GTGYCAGCMGCCGCGGTAA',  'SSU V4 (evrensel)'),
    ('806R_rc', 'ATTAGAWACCCBNGTAGTCC', 'SSU V4 sonu (evrensel)'),
    ('1492R',   'GGTTACCTTGTTACGACTT',  'SSU 3\' ucu (bakteri/arke)'),
    ('23S_129F','CCGAATGGGGRAACCC',     'LSU/23S 5\' bolgesi (bakteri)'),
    ('23S_2241R','ACCGCCCCAGTHAAACT',   'LSU/23S 3\' bolgesi (bakteri)'),
    ('NS1',     'GTAGTCATATGCTTGTCTC',  'okaryot SSU 5\' ucu'),
    ('ITS1',    'TCCGTAGGTGAACCTGCGG',  'SSU sonu / ITS1 basi (mantar)'),
    ('ITS3',    'GCATCGATGAAGAACGCAGC', '5.8S (mantar)'),
    ('ITS4',    'TCCTCCGCTTATTGATATGC', 'ITS2 sonu / LSU basi (mantar)'),
    ('NL1',     'GCATATCAATAAGCGGAGGAAAAG', 'LSU D1 basi (mantar)'),
    ('NL4',     'GGTCCGTGTTTCAAGACGG',  'LSU D2 sonu (mantar)'),
]


def ham_boylar(yol, tavan=4000):
    """Read lengths WITH NO FILTER. The sample filter (200-6000) IS NOT APPLIED here;
    the aim is precisely to see the long and short reads the filter could hide.

    """
    ac = gzip.open if yol.endswith('.gz') else open
    boy = []
    with ac(yol, 'rt', errors='ignore') as fh:
        for i, s in enumerate(fh):
            if i % 4 == 1:
                boy.append(len(s.strip()))
                if tavan and len(boy) >= tavan:
                    break
    return boy


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--start', '--bas', dest='bas', type=int, default=0)
    p.add_argument('--count', '--adet', dest='adet', type=int, default=1)
    a = p.parse_args()

    P = Puanlayici(KOK)
    om = P.motor.okuma_motoru
    kons = {}
    ix = os.path.join(KOK, 'konsensus_kanonik', 'INDEKS.tsv')
    import csv as _csv
    for r in _csv.DictReader(open(ix, encoding='utf-8'), delimiter='\t'):
        kons[r['kutu']] = dict(uzunluk=int(r['uzunluk']),
                               yol=os.path.join(os.path.dirname(ix), r['dosya']))

    os.makedirs('/tmp/mrb/kontrol', exist_ok=True)
    durum = json.load(open(KP, encoding='utf-8')) if os.path.exists(KP) else {}
    sondalar = [(ad, om.Sonda(dz, False, 2, False), parca) for ad, dz, parca in CAPALAR]

    for hedef in HEDEFLER[a.bas:a.bas + a.adet]:
        b = P.baglam(hedef)
        if not b:
            durum[hedef] = dict(hata='uyelik yok'); continue
        kutular = []
        for k in b['uye']:
            boy = ham_boylar(k['yol'])
            kseq = ''
            ku = kons.get(k['kutu'])
            if ku and os.path.exists(ku['yol']):
                kseq = ''.join(l.strip() for l in open(ku['yol'], encoding='utf-8',
                                                       errors='ignore')
                               if not l.startswith('>')).upper()
            # the anchor scan: the consensus plus the first 300 reads (in both directions)
            capa = {}
            ac = gzip.open if k['yol'].endswith('.gz') else open
            ornek = []
            with ac(k['yol'], 'rt', errors='ignore') as fh:
                for i, s in enumerate(fh):
                    if i % 4 == 1:
                        ornek.append(om.temizle(s.strip()))
                        if len(ornek) >= 300:
                            break
            for ad, sd, parca in sondalar:
                n_kons = 0
                if kseq:
                    n_kons = len(sd.bul(om.temizle(kseq))) + len(sd.bul(om.rc(om.temizle(kseq))))
                n_ok = sum(1 for s in ornek if sd.bul(s) or sd.bul(om.rc(s)))
                if n_kons or n_ok:
                    capa[ad] = dict(konsensus=n_kons, okuma=n_ok, okuma_taban=len(ornek),
                                    parca=parca)
            kutular.append(dict(
                kutu=k['kutu'], sinif=k['sinif'],
                konsensus_bp=(ku['uzunluk'] if ku else None),
                okuma_sayisi=len(boy),
                boy_medyan=int(statistics.median(boy)) if boy else 0,
                boy_min=min(boy) if boy else 0, boy_maks=max(boy) if boy else 0,
                boy_ust_2500=sum(1 for x in boy if x > 2500),
                boy_ust_1800=sum(1 for x in boy if x > 1800),
                capa=capa))
        durum[hedef] = dict(uye_kutu=len(b['uye']), rakip_kutu=len(b['rakip']),
                            karisik_kutu=len(b['karisik']), kutular=kutular)
        mx = max([x['boy_maks'] for x in kutular] or [0])
        md = [x['boy_medyan'] for x in kutular]
        bulunan = sorted({a2 for x in kutular for a2 in x['capa']})
        print(u'%-46s member=%2d  longest read=%5d  median %s..%s' % (
            hedef[:46], len(b['uye']), mx, min(md or [0]), max(md or [0])))
        print('    capalar: %s' % (', '.join(bulunan) or 'YOK'))
        json.dump(durum, open(KP, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    print('--- kayit: %d / %d' % (len(durum), len(HEDEFLER)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
