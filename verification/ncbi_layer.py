# -*- coding: utf-8 -*-
"""KATMAN 4 (NCBI Primer-BLAST) TEK BASINA - duzeltilmis hukum kuraliyla.

NEDEN AYRI BIR BETIK
--------------------
Dogrulama turunun tamami 1 sa 55 dk suruyor ve bunun 1 sa 40 dk'si katman 2
(yerel veritabani taramasi). Katman 4'un hukum kurali 2026-08-10'da degisti;
katman 2-3'un sonuclari degismedi. Butun turu bastan kosmak, degismeyen 1 sa
40 dk'yi bosa harcamak demek. Bu betik YALNIZ katman 4'u kosar (~30 dk) ve
sonucu ayri bir dosyaya yazar; birlestirme adimi onu okur.

HUKUM KURALI NEDEN DEGISTI (olculdu, tahmin degil)
--------------------------------------------------
Primer-BLAST'a SABLON DIZI bildirilmedigi surece bulunan her urunu
"Products on target templates" hanesine koyar; "potentially unintended
templates" bolumu hic acilmaz. 2026-08-07 kosusunun 22 ham sayfasinin
22'sinde boyleydi. Sayfanin ust kismindaki "Products on potentially
unintended templates" yazisi bir bolum basligi degil, BAG LISTESIDIR.
Bolum sayimina dayali eski hukum bu yuzden yapisal olarak calismiyordu.

2026-08-10, canli tek cift denemesi (Proteolitik_Cloacimonas, hedefin kendi
taksonu txid112 ENTREZ_QUERY ile dislanmis): 42 urun dondu, 42'si de
"uncultured bacterium clone ..." yani adsiz cevre klonu, ADLI takson sifir.
Adsiz klonlar hicbir taksona bagli olmadigi icin dislama suzgeci onlara
islemez; ustelik hedeflerimiz zaten adlandirilmamis soylar oldugundan bu
klonlarin buyuk kismi HEDEFIN KENDISI olabilir. Yani NCBI etiketi bizim
hedeflerimiz icin karar veremez.

Yeni kural: hukum BASLIGA bakar.
  ADLI takson  -> hedef disi kaniti (hedefin kendi taksonu zaten dislandi)
  adsiz klon   -> "karar veremez" hanesi; hukme GIRMEZ, sayisi raporda durur.
                  Kimligine dizi karsilastirmasi (katman 2-3) karar verir.

Kosum:
    python verification/ncbi_layer.py --kok .
Cikti:
    DOGRULAMA_SONUC/ncbi_katman4.tsv
    DOGRULAMA_SONUC/NCBI_KATMAN4_RAPORU.md
"""
from __future__ import print_function

import argparse
import io
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import specificity_round as D                                   # noqa: E402

SURUM = '1.0 (2026-08-10)'


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', '--kok', dest='kok', default='.')
    p.add_argument('--organism', '--organizma', dest='organizma', default='',
                   help='NCBI organizma kisiti (bos = tum nt). BOS BIRAKIN: '
                        'genel kisit butun urunleri "target templates" hanesine '
                        'itiyor (olculdu).')
    p.add_argument('--wait', '--bekleme', dest='bekleme', type=int, default=20)
    p.add_argument('--species-max', '--tur-ust', dest='tur_ust', type=int, default=60)
    a = p.parse_args()

    kok = os.path.abspath(a.kok)
    cikti = os.path.join(kok, 'DOGRULAMA_SONUC')
    os.makedirs(cikti, exist_ok=True)
    gun = io.open(os.path.join(cikti, 'ncbi_katman4_gunluk.txt'), 'a',
                  encoding='utf-8')

    def yaz(s=''):
        print(s, flush=True)
        gun.write(s + u'\n')
        gun.flush()

    yaz(u'=' * 78)
    yaz(u'  LAYER 4 - NCBI Primer-BLAST (on its own)   version %s' % SURUM)
    yaz(u'  %s' % time.strftime('%Y-%m-%d %H:%M'))
    yaz(u'=' * 78)

    # --- girdi: panelin TAMAMI (kosullu ve onerilmezler dahil) ---
    ciftler, yol = D.siparistekiler(kok, hepsi=True)
    yaz(u'  input: %s' % yol)
    yaz(u'  pairs: %d' % len(ciftler))
    if not ciftler:
        yaz(u'  INPUT EMPTY - nothing was run.')
        return 2

    # --- dislama haritasi ONCEDEN denetlenir: eksik hedef varsa BASTAN soyle ---
    hy = os.path.join(kok, 'screening', 'hedef_taxid.tsv')
    harita = {}
    if os.path.exists(hy):
        for l in io.open(hy, encoding='utf-8'):
            l = l.rstrip('\n')
            if not l.strip() or l.startswith('#'):
                continue
            pp = l.split('\t')
            if len(pp) >= 2 and pp[1].strip():
                harita[pp[0].strip()] = pp[1].strip()
    eksik = [c['hedef'] for c in ciftler if c['hedef'] not in harita]
    yaz(u'  exclusion map: %d targets listed (%s)' % (len(harita), hy))
    if eksik:
        yaz(u'  WARNING: %d targets are NOT in the map. For those, the target\'s own members cannot be excluded and the cell is written as NOT TESTED:' % len(eksik))
        for e in eksik:
            yaz(u'      - %s' % e)

    t0 = time.time()
    sonuc = D.katman2_oto(ciftler, cikti, yaz, a.organizma,
                          bekleme=a.bekleme, tur_ust=a.tur_ust, haric_taxid='')
    gecen = time.time() - t0

    # --- TSV ---
    ty = os.path.join(cikti, 'ncbi_katman4.tsv')
    bas = ['hedef', 'durum', 'adli_hedef_disi', 'adsiz_klon', 'toplam_urun',
           'dislanan_taxid', 'ornek_adli_baslik', 'not']
    with io.open(ty, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# LAYER 4 - NCBI Primer-BLAST. Generated %s, version %s\n'
                 % (time.strftime('%Y-%m-%d %H:%M'), SURUM))
        fh.write(u'# adli_hedef_disi = the number that ENTERS the verdict. These are the NAMED\n#   taxa left after the target\'s own taxon is excluded via ENTREZ_QUERY.\n#   They are evidence of off-target product.\n')
        fh.write(u'# adsiz_klon = records whose title carries "uncultured / unidentified /\n#   clone / enrichment culture". They are attached to no taxon, so the\n#   exclusion filter does not reach them, and because our targets are\n#   unnamed lineages they may BE THE TARGET ITSELF.\n#   THEY DO NOT ENTER THE VERDICT. Layers 2 and 3 (sequence comparison)\n#   decide what they are.\n')
        fh.write(u'\t'.join(bas) + u'\n')
        for c in ciftler:
            ad = c['hedef']
            v = sonuc.get(ad, {}) or {}
            fh.write(u'\t'.join([
                ad,
                str(v.get('durum', 'YOK')),
                str(v.get('hedef_disi', '')),
                str(v.get('ncbi_adsiz_klon', '')),
                str(v.get('ncbi_toplam_urun', '')),
                harita.get(ad, ''),
                str(v.get('ncbi_ornek', '')).replace('\t', ' '),
                str(v.get('not_', '')).replace('\t', ' ').replace('\n', ' '),
            ]) + u'\n')

    # --- rapor ---
    tamam = [k for k, v in sonuc.items() if str(v.get('durum', '')).startswith('TAMAM')]
    temiz = [k for k in tamam if (sonuc[k].get('hedef_disi') or 0) == 0]
    kirli = [k for k in tamam if (sonuc[k].get('hedef_disi') or 0) > 0]
    dusen = [k for k, v in sonuc.items() if not str(v.get('durum', '')).startswith('TAMAM')]

    ry = os.path.join(cikti, 'NCBI_KATMAN4_RAPORU.md')
    with io.open(ry, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# Layer 4 - NCBI Primer-BLAST\n\n')
        fh.write(u'Generated: %s, version %s, elapsed %s\n\n'
                 % (time.strftime('%Y-%m-%d %H:%M'), SURUM, D.sure_metni(gecen)))
        fh.write(u'## The verdict rule\n\n')
        fh.write(u'Unless a template sequence is declared to Primer-BLAST, it puts every product it finds under "target templates" and the "unintended" section never opens. The verdict therefore reads the product headings, not the section heading. After the target\'s own taxon is excluded via ENTREZ_QUERY, the **named** taxa that remain are evidence of off-target product. **Unnamed** environmental clones ("uncultured ...") are attached to no taxon, so the filter does not reach them, and because our targets are unnamed lineages they may be the target itself. They do not enter the verdict; layers 2 and 3 decide what they are.\n\n')
        fh.write(u'## The numbers\n\n')
        fh.write(u'| result | how many pairs |\n|---|---|\n')
        fh.write(u'| no named off-target | %d |\n' % len(temiz))
        fh.write(u'| named off-target present | %d |\n' % len(kirli))
        fh.write(u'| not testable | %d |\n\n' % len(dusen))
        if kirli:
            fh.write(u'## Pairs with a named off-target product\n\n')
            fh.write(u'| target | named | unnamed | total | example header |\n|---|---|---|---|---|\n')
            for k in sorted(kirli, key=lambda x: -(sonuc[x].get('hedef_disi') or 0)):
                v = sonuc[k]
                fh.write(u'| %s | %s | %s | %s | %s |\n'
                         % (k, v.get('hedef_disi'), v.get('ncbi_adsiz_klon'),
                            v.get('ncbi_toplam_urun'), (v.get('ncbi_ornek') or '')[:110]))
            fh.write(u'\n')
        if dusen:
            fh.write(u'## Could not be tested, and why\n\n')
            for k in sorted(dusen):
                fh.write(u'- **%s** — %s\n' % (k, (sonuc[k].get('not_') or '')[:300]))
            fh.write(u'\n')
        fh.write(u'## What these numbers are NOT: read this first\n\n')
        fh.write(u'This layer measures against **the whole of GenBank**, that is, against every record in the world. "650 named off-target" does NOT mean there are 650 organisms in the digester sample; it means the primer could also bind 650 named relatives worldwide.\n\n')
        fh.write(u'The layer that measures what is in the sample is a different one (layer 2, local databases and sample reads), and the dCq that decides the order comes from there. The two layers answer different questions:\n\n')
        fh.write(u'| layer | question | scope |\n|---|---|---|\n')
        fh.write(u'| 1-2 | does this primer discriminate its target IN THIS SAMPLE | the bins in the sample |\n')
        fh.write(u'| 4 (this one) | what else would this primer bind IN THE WORLD | GenBank nt |\n\n')
        fh.write(u'For broad targets (the Bacteroidales cluster, universal primers, fungal ITS) a high world count is EXPECTED behaviour: those primers were designed to bind a broad clade, and the only taxon excluded is the target\'s own clade. For narrow targets, a high count is a real warning.\n\n')
        fh.write(u'## What this layer cannot measure\n\n')
        fh.write(u'- The identity of unnamed environmental clones. They carry no label; this layer counts them but does not classify them.\n')
        fh.write(u'- For pairs where the page hit its cap, the total product count is a lower bound, not a count.\n')

    yaz(u'')
    yaz(u'  written: %s' % ty)
    yaz(u'  written: %s' % ry)
    yaz(u'  named off-target NONE: %d | PRESENT: %d | not testable: %d | time: %s'
        % (len(temiz), len(kirli), len(dusen), D.sure_metni(gecen)))
    gun.close()
    return 0 if tamam else 1


if __name__ == '__main__':
    sys.exit(main())
