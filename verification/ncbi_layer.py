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
    yaz(u'  KATMAN 4 - NCBI Primer-BLAST (tek basina)   surum %s' % SURUM)
    yaz(u'  %s' % time.strftime('%Y-%m-%d %H:%M'))
    yaz(u'=' * 78)

    # --- girdi: panelin TAMAMI (kosullu ve onerilmezler dahil) ---
    ciftler, yol = D.siparistekiler(kok, hepsi=True)
    yaz(u'  girdi: %s' % yol)
    yaz(u'  cift sayisi: %d' % len(ciftler))
    if not ciftler:
        yaz(u'  GIRDI BOS - kosu yapilmadi.')
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
    yaz(u'  dislama haritasi: %d hedef yazili (%s)' % (len(harita), hy))
    if eksik:
        yaz(u'  UYARI: haritada OLMAYAN %d hedef var. Bunlar icin hedefin kendi'
            u' uyeleri ayirt edilemez ve hucre SINANMADI yazilir:' % len(eksik))
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
        fh.write(u'# KATMAN 4 - NCBI Primer-BLAST. Uretim %s, surum %s\n'
                 % (time.strftime('%Y-%m-%d %H:%M'), SURUM))
        fh.write(u'# adli_hedef_disi = hukme GIREN sayi. Hedefin kendi taksonu\n'
                 u'#   ENTREZ_QUERY ile dislandiktan sonra geriye kalan ADLI\n'
                 u'#   taksonlar. Bunlar hedef disi urun kanitidir.\n')
        fh.write(u'# adsiz_klon = "uncultured/unidentified/clone/enrichment culture"\n'
                 u'#   basligi tasiyan kayitlar. Hicbir taksona bagli olmadiklari\n'
                 u'#   icin dislama suzgeci onlara islemez ve hedeflerimiz\n'
                 u'#   adlandirilmamis soylar oldugundan HEDEFIN KENDISI olabilirler.\n'
                 u'#   HUKME GIRMEZ. Kimliklerine katman 2-3 (dizi karsilastirmasi)\n'
                 u'#   karar verir.\n')
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
        fh.write(u'# Katman 4 - NCBI Primer-BLAST\n\n')
        fh.write(u'Üretim: %s · sürüm %s · geçen süre %s\n\n'
                 % (time.strftime('%Y-%m-%d %H:%M'), SURUM, D.sure_metni(gecen)))
        fh.write(u'## Hüküm kuralı\n\n')
        fh.write(u'Primer-BLAST\'a şablon dizi bildirilmediği sürece bulduğu her ürünü '
                 u'"target templates" hanesine koyuyor, "unintended" bölümü hiç açılmıyor. '
                 u'Bu yüzden hüküm bölüm başlığına değil ürün başlığına bakıyor. Hedefin '
                 u'kendi taksonu ENTREZ_QUERY ile dışlandıktan sonra geriye kalan **adlı** '
                 u'taksonlar hedef dışı ürün kanıtıdır. **Adsız** çevre klonları '
                 u'("uncultured ...") hiçbir taksona bağlı olmadığı için süzgeç onlara '
                 u'işlemez ve hedeflerimiz adlandırılmamış soylar olduğundan hedefin '
                 u'kendisi olabilirler; hükme girmezler, kimliklerine katman 2-3 karar verir.\n\n')
        fh.write(u'## Sayılar\n\n')
        fh.write(u'| sonuç | kaç çift |\n|---|---|\n')
        fh.write(u'| adlı hedef dışı YOK | %d |\n' % len(temiz))
        fh.write(u'| adlı hedef dışı VAR | %d |\n' % len(kirli))
        fh.write(u'| sınanamadı | %d |\n\n' % len(dusen))
        if kirli:
            fh.write(u'## Adlı hedef dışı ürün bulunan çiftler\n\n')
            fh.write(u'| hedef | adlı | adsız | toplam | örnek başlık |\n|---|---|---|---|---|\n')
            for k in sorted(kirli, key=lambda x: -(sonuc[x].get('hedef_disi') or 0)):
                v = sonuc[k]
                fh.write(u'| %s | %s | %s | %s | %s |\n'
                         % (k, v.get('hedef_disi'), v.get('ncbi_adsiz_klon'),
                            v.get('ncbi_toplam_urun'), (v.get('ncbi_ornek') or '')[:110]))
            fh.write(u'\n')
        if dusen:
            fh.write(u'## Sınanamayanlar ve sebebi\n\n')
            for k in sorted(dusen):
                fh.write(u'- **%s** — %s\n' % (k, (sonuc[k].get('not_') or '')[:300]))
            fh.write(u'\n')
        fh.write(u'## Bu sayılar ne DEĞİLDİR — önce bunu okuyun\n\n')
        fh.write(u'Bu katman **GenBank\'in tamamına** karşı ölçer, yani dünyadaki '
                 u'her kayda karşı. "Adlı hedef dışı 650" demek, çürütücü '
                 u'numunesinde 650 organizma var demek DEĞİLDİR; o primerin '
                 u'dünyada 650 adlandırılmış akrabayı da tutabileceği demektir.\n\n')
        fh.write(u'Numunede ne olduğunu ölçen katman ayrıdır (katman 2, yerel '
                 u'veritabanları ve numune okumaları) ve sipariş kararını veren '
                 u'dCq oradan gelir. İki katman farklı soruları yanıtlar:\n\n')
        fh.write(u'| katman | soru | kapsam |\n|---|---|---|\n')
        fh.write(u'| 1-2 | bu primer BU NUMUNEDE hedefini ayırt ediyor mu | '
                 u'numunedeki kutular |\n')
        fh.write(u'| 4 (bu) | bu primer DÜNYADA başka neyi tutar | GenBank nt |\n\n')
        fh.write(u'Geniş kapsamlı hedeflerde (Bacteroidales kümesi, evrensel '
                 u'primerler, mantar ITS) dünya sayısının yüksek çıkması '
                 u'BEKLENEN davranıştır: o primerler zaten geniş bir kladı '
                 u'tutmak için tasarlandı ve dışlanan takson yalnızca hedefin '
                 u'kendi kladıdır. Dar hedeflerde ise yüksek sayı gerçek bir '
                 u'uyarıdır.\n\n')
        fh.write(u'## Bu katmanın ölçemedikleri\n\n')
        fh.write(u'- Adsız çevre klonlarının kimliği. Etiket taşımıyorlar; '
                 u'bu katman onları sayar ama sınıflandırmaz.\n')
        fh.write(u'- Sayfa tavana çarpan çiftlerde toplam ürün sayısı bir alt sınırdır, '
                 u'sayım değildir.\n')

    yaz(u'')
    yaz(u'  yazildi: %s' % ty)
    yaz(u'  yazildi: %s' % ry)
    yaz(u'  adli hedef disi YOK: %d | VAR: %d | sinanamadi: %d | sure: %s'
        % (len(temiz), len(kirli), len(dusen), D.sure_metni(gecen)))
    gun.close()
    return 0 if tamam else 1


if __name__ == '__main__':
    sys.exit(main())
