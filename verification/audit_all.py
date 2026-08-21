# -*- coding: utf-8 -*-
"""HER KOSUDA CALISAN DENETIM KAPISI  -  "sormadan bak" isini kodun isi yapar.

NEDEN VAR
---------
Bu projede hatalarin cogu olcumde degil, olcumun DAYANDIGI tabloda cikti:
bayat onbellek, bayat referans, hedefin adina gore secilmis dislama taksidi,
NCBI'nin degistirdigi taksonomi, birbirine uymayan anahtar adlari. Hicbiri
kosuyu dusurmuyordu; hepsi sonucu SESSIZCE yanlis yapiyordu ve ancak birisi
"baska hata var mi" diye sordugunda bulunuyordu.

Sormaya bagli denetim denetim degildir. Bu betik o soruyu her kosuda kendisi
sorar. Hicbir olcum yapmaz, hicbir dosya degistirmez; yalnizca bakar ve
tutmayani yazar.

DENETLENENLER
  1  Dislama haritasi anahtarlari SIPARIS_LISTESI ile birebir mi
  2  Dislama taksidleri NCBI'da var mi ve uyeleri KAPSIYOR mu (ag gerekir)
  3  Hizli test referanslari hangi ciftten olculdu, cift o gunden beri degisti mi
  4  Panel kaynagi ile SIPARIS_LISTESI ayni dizileri mi tasiyor
  5  Kontrol noktasi muhurleri diziyi iceriyor mu (bayat onbellek tuzagi)
  6  Cikti dosyalari girdilerinden TAZE mi
  7  Paneldeki her cift geometri kapisindan AYNI diziyle gecti mi
  8  Ayni plakadaki urunler jelde ayrilir mi, bant sinifi uygun mu
  9  Yol gosterici belgelerdeki cift sayisi panelin bugunku sayisi mi
 10  NCBI ad kurali bilinen cevapli sinavi geciyor mu
 11  Siparis sayfasindaki (xlsx) diziler panelin su anki dizileri mi
 12  KOSULSUZ siparis edilecek her ciftin kaniti tam mi
 13  Tablodaki urun boyu, ciftin konsensusta URETTIGI boy mu
 14  .bat dosyalari CRLF ve saf ASCII mi, goto hedefleri var mi
 15  Kanonik konsensus klasorunde indekste olmayan kalinti var mi
 16  P ile K ayni uyelik dosyasini mi okuyor
 17  Indeksli veritabanlari DNA alfabesinde mi, ikizler gercekten ikiz mi
 18  Evrensel/kontrol primerleri kendi hedeflerimizi goruyor mu

Cikis kodu: 0 hepsi temiz, 1 en az bir denetim dustu, 2 ag gerektiren
denetim atlandi ama yerel denetimler temiz.

Kosum:
    python verification/audit_all.py --kok .
    python verification/audit_all.py --kok . --agsiz     (NCBI adimlarini atla)
"""
from __future__ import print_function

import argparse
import csv
import hashlib
import io
import os
import subprocess
import sys
import time

BULGU = []
ATLANAN = []


# Onem dereceleri. Sozcukler bilerek acik: "BLOKE" bir yargi degil, "bu
# giderilmeden siparis verilirse para ve zaman kaybi kesin" demek.
BLOKE, DIKKAT, BILGI = u'SIPARISI DURDURUR', u'DIKKAT', u'BILGI'


def bulgu(baslik, ayrinti, onem=None):
    BULGU.append((baslik, ayrinti, onem or DIKKAT))


def _tsv(yol):
    if not os.path.exists(yol):
        return []
    with io.open(yol, encoding='utf-8') as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith('#')),
                                   delimiter='\t'))


def _harita(yol):
    h = {}
    if not os.path.exists(yol):
        return h
    for l in io.open(yol, encoding='utf-8'):
        l = l.rstrip('\n')
        if not l.strip() or l.startswith('#') or l.startswith('hedef\t'):
            continue
        p = l.split('\t')
        if len(p) >= 2:
            h[p[0].strip()] = p[1].strip()
    return h


# --- 1 ------------------------------------------------------------------
def d1_harita_anahtarlari(kok, yaz):
    sl = os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'SIPARIS_LISTESI.tsv')
    hy = os.path.join(kok, 'screening', 'hedef_taxid.tsv')
    if not os.path.exists(sl) or not os.path.exists(hy):
        ATLANAN.append(u'1 harita anahtarlari (dosya yok)')
        return
    adlar = set((r.get('hedef') or '').strip() for r in _tsv(sl) if r.get('hedef'))
    H = _harita(hy)
    eksik = sorted(adlar - set(H))
    fazla = sorted(set(H) - adlar)
    yaz(u'  [1] dislama haritasi anahtarlari: %d hedef, %d harita satiri'
        % (len(adlar), len(H)))
    if eksik:
        bulgu(u'Dislama haritasinda OLMAYAN hedef',
              u'%s\n      Bu hedefler NCBI katmaninda SINANMADI yazilir. Anahtar '
              u'adi SIPARIS_LISTESI ile BIREBIR ayni olmalidir.' % ', '.join(eksik))
    if fazla:
        bulgu(u'Haritada olup panelde OLMAYAN anahtar',
              u'%s\n      Muhtemelen ad degisti ve eski anahtar kaldi; olu satir '
              u'yanlis guven verir.' % ', '.join(fazla))
    bos = sorted(k for k, v in H.items() if not v)
    if bos:
        bulgu(u'Dislama taksidi BOS birakilmis hedef',
              u'%s\n      Bos taxid = o hedef NCBI katmaninda hic sinanmaz.'
              % ', '.join(bos))


# --- 2 ------------------------------------------------------------------
def d2_kapsama(kok, yaz, agsiz):
    bet = os.path.join(kok, 'screening', 'exclusion_coverage_check.py')
    if not os.path.exists(bet):
        ATLANAN.append(u'2 kapsama denetimi (betik yok)')
        return
    if agsiz:
        ATLANAN.append(u'2 kapsama denetimi (--agsiz verildi; NCBI Taxonomy gerekir)')
        return
    yaz(u'  [2] dislama kapsama denetimi kosuyor (NCBI Taxonomy)...')
    try:
        p = subprocess.run([sys.executable, bet, '--kok', kok],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           timeout=900)
        cik = p.stdout.decode('utf-8', 'replace')
    except Exception as e:
        ATLANAN.append(u'2 kapsama denetimi (kosturulamadi: %s)' % e)
        return
    for satir in cik.splitlines():
        if 'KAPSAMIYOR' in satir or 'denetlenemedi' in satir:
            yaz(u'      %s' % satir.strip()[:110])
    if p.returncode != 0:
        kotu = [s.strip() for s in cik.splitlines() if 'KAPSAMIYOR' in s]
        if kotu:
            bulgu(u'Dislama taksonu uyeleri KAPSAMIYOR',
                  u'\n      '.join(kotu[:8]) +
                  u'\n      Kapsanmayan uye hedefin KENDISI oldugu halde hedef disi '
                  u'sayilir. Ayrinti: screening/exclusion_coverage_check.py', BLOKE)
        else:
            yaz(u'      (kapsama tamam; yalniz denetlenemeyen hedef var)')


# --- 3 ------------------------------------------------------------------
def d3_referans_bayat(kok, yaz):
    ry = os.path.join(kok, 'HIZLI_TEST', 'referans_degerler.tsv')
    py_ = os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'panel_tek_protokol.tsv')
    if not os.path.exists(ry):
        bulgu(u'Hizli test referans dosyasi YOK',
              u'HIZLI_TEST/referans_degerler.tsv uretilmemis. Referans koda gomulu '
              u'sabitlerden okunur ve hangi ciftten olculdugu bilinmez. Uretmek '
              u'icin: python verification/refresh_reference.py --kok .')
        return
    if not os.path.exists(py_):
        ATLANAN.append(u'3 referans tazeligi (panel_tek_protokol.tsv yok)')
        return
    ref = {}
    for r in _tsv(ry):
        ref[(r.get('hedef') or '').strip()] = ((r.get('F') or '').strip().upper(),
                                               (r.get('R') or '').strip().upper())
    bayat = []
    for r in _tsv(py_):
        ad = (r.get('hedef') or '').strip()
        if ad not in ref:
            continue
        simdi = ((r.get('F') or '').strip().upper(), (r.get('R') or '').strip().upper())
        if simdi[0] and simdi != ref[ad]:
            bayat.append(ad)
    yaz(u'  [3] referans tazeligi: %d satir referansta, %d cift degismis'
        % (len(ref), len(bayat)))
    if bayat:
        bulgu(u'Referans BAYAT - primer cifti degismis',
              u'%s\n      Bu satirlarda eski cifte ait sayi yeni ciftle '
              u'karsilastirilir ve sahte "gerileme" uretir. Yenile: '
              u'python verification/refresh_reference.py --kok .' % ', '.join(bayat), DIKKAT)


# --- 4 ------------------------------------------------------------------
def d4_kaynak_tutarliligi(kok, yaz):
    """Panelin dizi kaynagi ile siparis listesi ayni diziyi mi soyluyor."""
    sl = os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'SIPARIS_LISTESI.tsv')
    pk = os.path.join(kok, 'primer_final',
                      'devir_ciftleri_20260802_sonrotus_TESLIM.tsv')
    if not os.path.exists(sl) or not os.path.exists(pk):
        ATLANAN.append(u'4 kaynak tutarliligi (dosya yok)')
        return
    panel = {}
    with io.open(pk, encoding='utf-8') as fh:
        sat = [l.rstrip('\n').split('\t') for l in fh]
    if not sat:
        ATLANAN.append(u'4 kaynak tutarliligi (panel kaynagi bos)')
        return
    bas = sat[0]
    try:
        iH = bas.index('Hedef')
        iF = [i for i, b in enumerate(bas) if b.startswith('Ileri primer')][0]
        iR = [i for i, b in enumerate(bas) if b.startswith('Geri primer')][0]
    except (ValueError, IndexError):
        ATLANAN.append(u'4 kaynak tutarliligi (panel sutunlari taninmadi)')
        return
    for r in sat[1:]:
        if len(r) > max(iF, iR) and r[iH].strip():
            panel[r[iH].strip()] = (r[iF].strip().upper(), r[iR].strip().upper())
    fark = []
    liste_adlari = set()
    for r in _tsv(sl):
        ad = (r.get('hedef') or '').strip()
        if not ad:
            continue
        liste_adlari.add(ad)
        if ad not in panel:
            continue
        simdi = ((r.get('F') or '').strip().upper(), (r.get('R') or '').strip().upper())
        if simdi[0] and panel[ad][0] and simdi != panel[ad]:
            fark.append(ad)
    # EKSIK SATIR da bir farktir. Yalniz ORTAK hedefleri karsilastirmak,
    # bir dosyada olup otekinde hic olmayan cifti gormezden gelir; 2026-08-10'da
    # Petriella_cinsi siparis listesinde KESIN yaziliyken panel kaynaginda
    # satiri yoktu ve bu denetim "0 fark" diyordu.
    yalniz_listede = sorted(liste_adlari - set(panel))
    yalniz_panelde = sorted(set(panel) - liste_adlari)
    yaz(u'  [4] kaynak tutarliligi: %d ortak hedef, %d dizi farki, '
        u'%d yalniz listede, %d yalniz panelde'
        % (len(liste_adlari & set(panel)), len(fark),
           len(yalniz_listede), len(yalniz_panelde)))
    if yalniz_listede:
        bulgu(u'Siparis listesinde VAR, panel kaynaginda YOK',
              u'%s\n      Bu ciftlerin plaka ve Ta bilgisi panel tablosunda yok; '
              u'siparise giriyorlar ama deney duzeninde yerleri belirsiz.'
              % ', '.join(yalniz_listede), BLOKE)
    if fark:
        bulgu(u'Panel kaynagi ile SIPARIS_LISTESI AYRISIYOR',
              u'%s\n      Iki dosya farkli dizi soyluyor. Hangisinin siparis '
              u'edilecegi belirsiz - once bu giderilmeli.' % ', '.join(fark), BLOKE)


# --- 5 ------------------------------------------------------------------
def d5_muhur_diziyi_iceriyor_mu(kok, yaz):
    """Kontrol noktasi muhurleri primer DIZISINI iceriyor mu.

    Icermeyen muhur, dizi degistiginde ESKI olcumu taze sanar. Bu proje bu
    hatayi bir gunde uc ayri betikte yasadi; artik mekanik olarak aranir.
    """
    bakilacak = [
        ('protocol/single_protocol_measure.py', True),
        ('verification/recovery_round.py', True),
        ('verification/specificity_round.py', True),
        ('verification/mfeprimer_layer.py', False),
        ('engine/rederive_membership.py', True),
    ]
    eksik = []
    for yol, zorunlu in bakilacak:
        t = os.path.join(kok, yol)
        if not os.path.exists(t):
            continue
        s = io.open(t, encoding='utf-8', errors='replace').read()
        # muhur/imza hesabinda F ve R gecmeli
        imzali = ("'F'" in s or '"F"' in s or "get('F'" in s)
        anahtar = ('md5' in s or 'sha1' in s or 'hashlib' in s)
        if anahtar and not imzali and zorunlu:
            eksik.append(yol)
    yaz(u'  [5] kontrol noktasi muhurleri: %d betik bakildi, %d suphe'
        % (len(bakilacak), len(eksik)))
    if eksik:
        bulgu(u'Muhur diziyi ICERMIYOR olabilir',
              u'%s\n      Muhurde primer dizisi yoksa dizi degistiginde eski '
              u'olcum TAZE sanilir. Elle bakin.' % ', '.join(eksik))


# --- 6 ------------------------------------------------------------------
def d6_cikti_tazeligi(kok, yaz, uretilecek=()):
    """uretilecek: bu kosuda YENIDEN URETILECEK ciktilar - onlara bayat denmez.

    2026-08-10: NCBI kapisi, koşunun kendi uretecegi ncbi_katman4.tsv'yi
    "bayat" diye gosterip onay sordu. Kapinin kurt masali anlatmasi, kapiyi
    ciddiye alinmaz hale getirir. Bu yuzden kosunun uretecegi dosya hariç
    tutulur ve hariç tutuldugu EKRANA YAZILIR - sessizce gizlenmez.
    """
    ciftler = [
        ('primer_final/devir_ciftleri_20260802_sonrotus_TESLIM.tsv',
         'TEK_PROTOKOL_SONUC/panel_tek_protokol.tsv'),
        ('TEK_PROTOKOL_SONUC/panel_tek_protokol.tsv',
         'HIZLI_TEST/referans_degerler.tsv'),
        ('TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv',
         'DOGRULAMA_SONUC/dogrulama_uc_sutun.tsv'),
        ('screening/hedef_taxid.tsv',
         'DOGRULAMA_SONUC/ncbi_katman4.tsv'),
        # 2026-08-11: uyelik tablosu degisince G asamasi tablosu bayatlar.
        # Siparis listesindeki kimlik sutunlari (olculen_kimlik, ad_farkli_mi,
        # uye sayisi) dogrudan o tablonun UYESI_OLDUGU_HEDEFLER sutunundan
        # geliyor. Uyelik duzeltildi ama G yeniden kosulmadiysa, siparis satiri
        # ARTIK UYE OLMAYAN kutularin adlarini sayar: Petriella_musispora
        # satiri "10/10 kutu" deyip Microascus, Lomentospora ve Graphium'u
        # sayiyordu, oysa olcum 9 kutuyla yapilmisti.
        ('screening/hedef_uyelik.tsv',
         'TUM_KIMLIK_SONUC/tum_kutu_kimlikleri.tsv'),
    ]
    def _diziler(y2, ad_h='hedef', ad_f='F', ad_r='R'):
        out = {}
        for r in _tsv(y2):
            a = (r.get(ad_h) or '').strip()
            if a:
                out[a] = ((r.get(ad_f) or '').strip().upper(),
                          (r.get(ad_r) or '').strip().upper())
        return out

    bayat = []
    for g, c in ciftler:
        # ICERIK ONCE, ZAMAN SONRA. panel kaynagi ile P ciktisi arasinda
        # yalniz Tm/urun boyu sutunlari degistiginde zaman damgasi "bayat"
        # der ama olcum GECERLIDIR: dCq diziye baglidir, Tm'e degil.
        # (2026-08-11 07:02'de tam bu oldu.) Diziler ayniysa bayat sayilmaz.
        if (g.endswith('devir_ciftleri_20260802_sonrotus_TESLIM.tsv')
                and c.endswith('panel_tek_protokol.tsv')):
            gp2, cp2 = os.path.join(kok, g), os.path.join(kok, c)
            if os.path.exists(gp2) and os.path.exists(cp2):
                sat2 = [l.rstrip('\n').split('\t') for l in io.open(gp2, encoding='utf-8')]
                b2 = sat2[0]
                try:
                    iH2 = b2.index('Hedef')
                    iF2 = next(k for k, x in enumerate(b2) if x.startswith('Ileri primer'))
                    iR2 = next(k for k, x in enumerate(b2) if x.startswith('Geri primer'))
                except (ValueError, StopIteration):
                    iH2 = None
                if iH2 is not None:
                    kay = {}
                    for r2 in sat2[1:]:
                        if len(r2) > max(iF2, iR2) and r2[iH2].strip():
                            kay[r2[iH2].strip()] = (r2[iF2].strip().upper(),
                                                    r2[iR2].strip().upper())
                    cik = _diziler(cp2)
                    ortak = set(kay) & set(cik)
                    if ortak and all(kay[k2] == cik[k2] for k2 in ortak):
                        yaz(u'      (%s zaman damgasi eski ama DIZILER ayni - '
                            u'olcum gecerli, bayat sayilmadi)' % c)
                        continue
        if c in uretilecek:
            yaz(u'      (%s bu kosuda yeniden uretilecek - tazelik aranmadi)' % c)
            continue
        gp, cp = os.path.join(kok, g), os.path.join(kok, c)
        if not os.path.exists(gp) or not os.path.exists(cp):
            continue
        if os.path.getmtime(cp) < os.path.getmtime(gp):
            bayat.append(u'%s (%s) < %s (%s)'
                         % (c, time.strftime('%d.%m %H:%M', time.localtime(os.path.getmtime(cp))),
                            g, time.strftime('%d.%m %H:%M', time.localtime(os.path.getmtime(gp)))))
    yaz(u'  [6] cikti tazeligi: %d bagimlilik, %d bayat' % (len(ciftler), len(bayat)))
    if bayat:
        bulgu(u'Cikti girdisinden ESKI',
              u'\n      '.join(bayat) +
              u'\n      Bu cikti guncel girdiyi gormemis; yeniden uretilmeli.')



# --- 7 ------------------------------------------------------------------
def d7_geometri_kapisi(kok, yaz):
    """Paneldeki her cift, geometri denetiminden AYNI diziyle gecmis mi.

    2026-08-10: alti ciftin dizisi 2 Agustos'taki geometri denetiminden sonra
    degistirilmis ama geometri yeniden kosulmamisti. Yani o alti cift panelin
    kendi kurallarindan (uzunluk, GC, Tm penceresi, sac tokasi, dimer) HIC
    gecmemisti ve bunu hicbir sey soylemiyordu.
    """
    import glob
    pk = os.path.join(kok, 'primer_final',
                      'devir_ciftleri_20260802_sonrotus_TESLIM.tsv')
    adaylar = sorted(glob.glob(os.path.join(kok, 'primer_final',
                                            'geometri_denetimi_*.tsv')))
    adaylar = [x for x in adaylar if 'yedek' not in x]
    if not os.path.exists(pk) or not adaylar:
        ATLANAN.append(u'7 geometri kapisi (dosya yok)')
        return
    gy = max(adaylar, key=os.path.getmtime)
    g = {}
    for r in _tsv(gy):
        h = (r.get('Hedef') or '').strip()
        pr = (r.get('Primer') or '').strip()
        if pr in ('Ileri', 'Geri'):
            g.setdefault(h, {})[pr] = (r.get('Dizi') or '').strip().upper()
    sat = [l.rstrip('\n').split('\t') for l in io.open(pk, encoding='utf-8')]
    bas = sat[0]
    try:
        iH = bas.index('Hedef')
        iF = next(i for i, b in enumerate(bas) if b.startswith('Ileri primer'))
        iR = next(i for i, b in enumerate(bas) if b.startswith('Geri primer'))
    except (ValueError, StopIteration):
        ATLANAN.append(u'7 geometri kapisi (panel sutunlari taninmadi)')
        return
    import re as _re
    gecmemis = []
    for r in sat[1:]:
        if len(r) <= max(iF, iR) or not r[iH].strip():
            continue
        if not _re.match(r'^[A-Za-z]', r[iH].strip()):
            continue
        F, R = r[iF].strip().upper(), r[iR].strip().upper()
        if not _re.fullmatch(r'[ACGT]+', F or ''):
            continue
        ad = r[iH].strip()
        if ad not in g or (F, R) != (g[ad].get('Ileri'), g[ad].get('Geri')):
            gecmemis.append(ad)
    yaz(u'  [7] geometri kapisi: kaynak %s, gecmemis %d'
        % (os.path.basename(gy), len(gecmemis)))
    if gecmemis:
        bulgu(u'Geometri kapisindan GECMEMIS cift',
              u'%s\n      Bu ciftlerin dizisi son geometri denetiminden sonra '
              u'degismis. Uzunluk, GC, Tm penceresi, sac tokasi ve dimer '
              u'kurallari bu diziler icin HIC olculmedi. Kosun: '
              u'python verification/refresh_geometry.py --kok .' % ', '.join(gecmemis), BLOKE)



# --- 8 ------------------------------------------------------------------
def d8_plaka_jel_ve_bant(kok, yaz):
    """Ayni plakada urunler jelde ayrilir mi + QuantiNova bant sinifi.

    2026-08-10: bugun degistirilen ciftler plaka ici jel ayrimini yeniden
    bozdu. Bacteroidales urunu 241 bp'den 150 bp'ye dondu ve Mantar F2'nin
    145 bp'siyle 5 bp'ye yaklasti - %2 agarozda ayirt edilemez. Cift
    degistirmek yalniz o cifti degil, PLAKAYI da etkiliyor; bu yuzden her
    kosuda bakilir.
    """
    import itertools
    import re as _re
    pk = os.path.join(kok, 'primer_final',
                      'devir_ciftleri_20260802_sonrotus_TESLIM.tsv')
    if not os.path.exists(pk):
        ATLANAN.append(u'8 plaka/jel (panel dosyasi yok)')
        return
    sat = [l.rstrip('\n').split('\t') for l in io.open(pk, encoding='utf-8')]
    bas = sat[0]
    try:
        iP, iT, iH = bas.index('Plaka'), bas.index('Ta (C)'), bas.index('Hedef')
        iU = bas.index('Urun (bp)')
        iF = next(i for i, b in enumerate(bas) if b.startswith('Ileri primer'))
    except (ValueError, StopIteration):
        ATLANAN.append(u'8 plaka/jel (sutunlar taninmadi)')
        return
    pl = {}
    bant = []
    for r in sat[1:]:
        if len(r) <= max(iU, iF) or not r[iH].strip():
            continue
        if not _re.fullmatch(r'[ACGT]+', (r[iF] or '').strip().upper()):
            continue          # not satirlarini eler - primer dizisi yoksa cift degildir
        try:
            u = int(_re.sub(r'\D', '', r[iU]))
        except ValueError:
            continue
        pl.setdefault((r[iP].strip(), r[iT].strip()), []).append((r[iH].strip(), u))
        if u > 250:
            bant.append(u'%s %d bp (QuantiNova >250 ONERMIYOR)' % (r[iH].strip(), u))
    cak = []
    for k, v in sorted(pl.items()):
        for (a1, u1), (a2, u2) in itertools.combinations(sorted(v, key=lambda x: x[1]), 2):
            if abs(u1 - u2) < 10:
                cak.append(u'plaka %s Ta %s: %s (%d bp) / %s (%d bp) - fark %d bp'
                           % (k[0], k[1], a1, u1, a2, u2, abs(u1 - u2)))
    yaz(u'  [8] plaka/jel: %d plaka grubu, %d cakisma, %d bant disi urun'
        % (len(pl), len(cak), len(bant)))
    if cak:
        bulgu(u'Plaka ici JEL AYRIMI cakismasi',
              u'\n      '.join(cak) +
              u'\n      Ayni plakada 10 bp\'den yakin iki urun %2 agarozda '
              u'ayirt edilemez. Ya plaka yeniden atanmali ya da fark kabul '
              u'edildigi RAPORDA yazili olmali.', DIKKAT)
    if bant:
        bulgu(u'Amplikon bant sinifi disinda',
              u'\n      '.join(bant) +
              u'\n      QuantiNova SYBR Green icin ideal 60-150 bp, 150-250 bp '
              u'30 sn uzatma ister, 250 uzeri onerilmez.', DIKKAT)



# --- 9 ------------------------------------------------------------------
# Sayinin ELLE yazildigi her belge bir gun bayatlar. Bu denetim, yol gosterici
# belgelerdeki "N cift" iddiasini panelin BUGUNKU sayisiyla karsilastirir.
# 2026-08-10: uc belge uc farkli sayi soyluyordu (16, 16, 11); dogrusu 20 idi.
# TARIHLI denetim kayitlari (SON_KONTROL.md gibi) bu listede YOKTUR: onlarin
# sayisi yazildigi gunun sayisidir ve degistirilmemelidir. Listede yalniz
# "su an ne yapmaliyim" sorusuna cevap veren belgeler bulunur.
YOL_GOSTERICI = ('OKU_ONCE.md', 'NASIL_DEVAM_EDILIR.md', 'CALISTIRMA_KILAVUZU.md',
                 'GUNCEL_DURUM.md')


def d9_belgelerde_bayat_sayi(kok, yaz):
    import re as _re
    sl = _tsv(os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'SIPARIS_LISTESI.tsv'))
    if not sl:
        ATLANAN.append(u'9 belge sayilari (SIPARIS_LISTESI yok)')
        return
    kesin = sum(1 for r in sl if (r.get('SINIF') or '').strip().upper() == 'KESIN')
    evr = sum(1 for r in sl if (r.get('SINIF') or '').strip().upper() == 'EVRENSEL')
    dogru = {kesin, evr, kesin + evr, len(sl), len(sl) - kesin - evr}
    kalip = _re.compile(r'(?:KES[İI]N|[Ss]ipari[şs] edilebilir|[Ss]ipari[şs] edilecek)'
                        r'[^0-9\n]{0,40}?(\d{1,2})\s*[çc]ift')
    kotu = []
    bakilan = 0
    for ad in YOL_GOSTERICI:
        t = os.path.join(kok, ad)
        if not os.path.exists(t):
            continue
        bakilan += 1
        metin = io.open(t, encoding='utf-8', errors='replace').read()
        for m in kalip.finditer(metin):
            n = int(m.group(1))
            if n not in dogru:
                kotu.append(u'%s: "%s" (bugunku dogru sayi: %d siparis, %d hedef '
                            u'ozgul, %d evrensel)'
                            % (ad, m.group(0).strip()[:60], kesin + evr, kesin, evr))
    yaz(u'  [9] belgelerdeki cift sayilari: %d belge bakildi, %d bayat iddia'
        % (bakilan, len(kotu)))
    if kotu:
        bulgu(u'Belgede BAYAT cift sayisi',
              u'\n      '.join(kotu) +
              u'\n      Sayiyi elle yazan her belge bayatlar. Cumleyi '
              u'GUNCEL_DURUM.md\'e isaret edecek sekilde degistirin; o dosya '
              u'her kosuda uretilir.', DIKKAT)



# --- 10 -----------------------------------------------------------------
def d10_ad_kurali_sinavi(kok, yaz):
    """NCBI ad kurali (adli/adsiz) bilinen cevapli sinavi geciyor mu.

    Bu kural hedef disi SAYISINI belirliyor. Bozulursa rapora giren tabloda
    "650 hedef disi" gibi sisirilmis sayilar cikar (2026-08-10'da tam olarak
    boyle oldu: gevsek kural 650 diyordu, siki kural 82).
    """
    bet = os.path.join(kok, 'verification', 'ncbi_reclassify.py')
    if not os.path.exists(bet):
        ATLANAN.append(u'10 ad kurali sinavi (betik yok)')
        return
    try:
        p = subprocess.run([sys.executable, bet, '--sina'],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           timeout=120)
        cik = p.stdout.decode('utf-8', 'replace').strip()
    except Exception as e:
        ATLANAN.append(u'10 ad kurali sinavi (kosturulamadi: %s)' % e)
        return
    son = [l for l in cik.splitlines() if 'sinavi:' in l]
    yaz(u'  [10] NCBI ad kurali: %s' % (son[-1].strip() if son else cik[:60]))
    if p.returncode != 0:
        bulgu(u'NCBI ad kurali sinavi DUSTU',
              cik[-600:] +
              u'\n      Bu kural hedef disi SAYISINI belirler; bozuksa rapordaki '
              u'sayilar yanlistir.', BLOKE)



# --- 11 -----------------------------------------------------------------
def d11_siparis_dizileri(kok, yaz):
    """Rapora/tedarikciye giden xlsx ile panelin SU ANKI dizileri ayni mi.

    2026-08-10 gece: sabah ozeti "diziler buradan kopyalanacak" diye
    PrimerJury_..._TESLIM.xlsx "2 Panel" sayfasini gosteriyordu; o sayfada
    ALTI ciftin dizisi eskiydi. Oradan siparis verilseydi 20 ciftin 6'si
    yanlis oligo olarak gelirdi. Bu, projenin verebilecegi en pahali hata.
    """
    import glob as _glob
    import re as _re
    pk = os.path.join(kok, 'primer_final',
                      'devir_ciftleri_20260802_sonrotus_TESLIM.tsv')
    # 2026-08-11: eski teslim xlsx'i arsive tasindi. Artik her kosuda
    # URETILEN PrimerJury_PANEL_*.xlsx denetlenir; o da yoksa denetim
    # atlanir. Arsivdeki dosyaya BAKILMAZ.
    _p = sorted(_glob.glob(os.path.join(kok, 'PrimerJury_PANEL_*.xlsx')))
    xl = _p[-1] if _p else ''
    if not os.path.exists(pk):
        ATLANAN.append(u'11 siparis dizileri (panel kaynagi yok)')
        return
    sat = [l.rstrip('\n').split('\t') for l in io.open(pk, encoding='utf-8')]
    b = sat[0]
    try:
        iH = b.index('Hedef')
        iF = next(i for i, x in enumerate(b) if x.startswith('Ileri primer'))
        iR = next(i for i, x in enumerate(b) if x.startswith('Geri primer'))
    except (ValueError, StopIteration):
        ATLANAN.append(u'11 siparis dizileri (panel sutunlari taninmadi)')
        return
    tsv = {}
    for r in sat[1:]:
        if len(r) > max(iF, iR) and r[iH].strip():
            F = r[iF].strip().upper()
            if _re.fullmatch(r'[ACGT]+', F or ''):
                tsv[r[iH].strip()] = (F, r[iR].strip().upper())
    if not xl or not os.path.exists(xl):
        yaz(u'  [11] siparis dizileri: uretilmis Excel yok - '
            u'python verification/build_excel.py --kok . ile uretin')
        ATLANAN.append(u'11 siparis dizileri (uretilmis Excel yok)')
        return
    try:
        import openpyxl
    except ImportError:
        ATLANAN.append(u'11 siparis dizileri (openpyxl yok: '
                       u'pip3 install openpyxl --break-system-packages)')
        return
    try:
        wb = openpyxl.load_workbook(xl, data_only=True, read_only=True)
        ws = wb['1 Siparis']
        rows = [[('' if c is None else str(c)).strip() for c in r]
                for r in ws.iter_rows(values_only=True)]
    except Exception as e:
        ATLANAN.append(u'11 siparis dizileri (xlsx okunamadi: %s)' % e)
        return
    # YENI Excel bicimi: "1 Siparis" sayfasi, sutunlar "oligo adi",
    # "dizi (5→3)", "hedef", "yon". Eski teslim dosyasinin bicimi degildi;
    # basligi eski adlarla aramak "baslik bulunamadi" verip denetimi
    # sessizce atlatiyordu (2026-08-11, ilk denemede oldu).
    bas = None
    for i, r in enumerate(rows[:12]):
        if 'oligo adi' in r and 'hedef' in r:
            bas = i
            break
    if bas is None:
        ATLANAN.append(u'11 siparis dizileri (Excel "1 Siparis" basligi bulunamadi)')
        return
    h = rows[bas]
    xH, xD, xY = h.index('hedef'), h.index('dizi (5→3)'), h.index('yon')
    xls = {}
    for r in rows[bas + 1:]:
        if len(r) <= max(xH, xD, xY) or not r[xH]:
            continue
        d = (r[xD] or '').upper()
        if not _re.fullmatch(r'[ACGT]+', d or ''):
            continue
        g = xls.setdefault(r[xH].strip(), {})
        g['F' if r[xY].strip().startswith('ileri') else 'R'] = d
    xls = {k: (v.get('F', ''), v.get('R', '')) for k, v in xls.items() if len(v) == 2}

    fark = sorted(k for k in set(xls) & set(tsv) if xls[k] != tsv[k])
    # EKSIK karsilastirmasi yalniz SIPARISE GIRENLER uzerinden yapilir.
    # Siparis sayfasinda ONERILMEZ ciftler zaten YOKTUR; hepsiyle
    # karsilastirmak "Proteiniphilum eksik" gibi yanlis alarm uretir
    # (2026-08-11, ilk denemede oldu).
    _sl = _tsv(os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'SIPARIS_LISTESI.tsv'))
    _sip = set((r.get('hedef') or '').strip() for r in _sl
               if (r.get('SINIF') or '').strip().upper() in ('KESIN', 'EVRENSEL'))
    eksik = sorted((set(tsv) & _sip) - set(xls)) if _sip else sorted(set(tsv) - set(xls))
    yaz(u'  [11] siparis dizileri: xlsx %d cift, panel %d cift, ayrisan %d'
        % (len(xls), len(tsv), len(fark)))
    if fark:
        bulgu(u'SIPARIS DIZILERI AYRISIYOR - xlsx ESKI',
              u'%s\n      "%s" dosyasindaki "2 Panel" sayfasinda bu ciftlerin '
              u'dizisi panelin su anki dizisinden FARKLI. O sayfadan siparis '
              u'verilirse YANLIS OLIGO gelir. Dogru liste: SIPARIS_FORMU.tsv '
              u'(uretilmis, elle yazilmamis).'
              % (', '.join(fark), os.path.basename(xl)), BLOKE)
    if eksik:
        bulgu(u'xlsx\'te OLMAYAN cift',
              u'%s\n      Panelde var, siparis sayfasinda yok.' % ', '.join(eksik))

    # BUTUN xlsx dosyalari: icinde primer dizisi tasiyan her dosya bir gun
    # siparis kaynagi sanilabilir. Hangilerinin bayat oldugunu ADIYLA yaz.
    import glob as _glob
    guncel = set(v[0] for v in tsv.values()) | set(v[1] for v in tsv.values())
    bayat_dosya = []
    for xy in sorted(_glob.glob(os.path.join(kok, '*.xlsx'))):
        try:
            w2 = openpyxl.load_workbook(xy, data_only=True, read_only=True)
        except Exception:
            continue
        bulundu = set()
        var = 0
        for sn in w2.sheetnames:
            for row in w2[sn].iter_rows(values_only=True):
                for c in row:
                    if isinstance(c, str):
                        v = c.strip().upper()
                        if 15 <= len(v) <= 30 and _re.fullmatch(r'[ACGT]+', v):
                            var += 1
                            if v in guncel:
                                bulundu.add(v)
        if var and len(bulundu) < len(guncel):
            bayat_dosya.append(u'%s (guncel dizi %d/%d)'
                               % (os.path.basename(xy), len(bulundu), len(guncel)))
    if bayat_dosya:
        bulgu(u'Icinde ESKI primer dizisi tasiyan xlsx dosyalari',
              u'\n      '.join(bayat_dosya) +
              u'\n      Bunlarin hicbiri siparis kaynagi DEGILDIR. Tek yetkili '
              u'liste: SIPARIS_FORMU.tsv (her kosuda uretilir).', BILGI)



# --- 12 -----------------------------------------------------------------
def d12_kanitsiz_kosulsuz(kok, yaz):
    """KOSULSUZ siparis edilecek bir ciftin kaniti eksik olabilir mi.

    2026-08-10: Petriella_cinsi siparis listesinde SINIF=KESIN,
    siparis_sarti=KOSULSUZ yaziyor; ama olculen_kimlik alani "G asamasinda
    uyelik tanimi YOK" diyor, uyelik tablosunda satiri yok ve panel
    kaynaginda plaka/Ta bilgisi yok. Yani kimligi dogrulanmamis bir hedefe
    kosulsuz siparis hukmu verilmis. "Kosulsuz" sozcugu raporda "bunda
    tartisilacak bir sey yok" diye okunur; oysa var.
    """
    sl = _tsv(os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'SIPARIS_LISTESI.tsv'))
    if not sl:
        ATLANAN.append(u'12 kanitsiz kosulsuz (SIPARIS_LISTESI yok)')
        return
    supheli = []
    for r in sl:
        ad = (r.get('hedef') or '').strip()
        sinif = (r.get('SINIF') or '').strip().upper()
        sart = (r.get('siparis_sarti') or '').strip().upper()
        if sinif != 'KESIN' or not sart.startswith('KOSULSUZ'):
            continue
        kim = (r.get('olculen_kimlik') or '').strip()
        uy = (r.get('uyusan_vtb_sayisi') or '').strip()
        eksik = []
        if not kim or 'YOK' in kim.upper() or 'TANIMI' in kim.upper():
            eksik.append(u'olculen kimlik: "%s"' % (kim or 'bos'))
        if uy in ('', '-', '0'):
            eksik.append(u'uyusan veritabani sayisi: "%s"' % (uy or 'bos'))
        if eksik:
            supheli.append(u'%s -> %s' % (ad, '; '.join(eksik)))
    yaz(u'  [12] kosulsuz siparis kaniti: %d satir bakildi, %d suphe'
        % (sum(1 for r in sl if (r.get('siparis_sarti') or '').strip().upper()
               .startswith('KOSULSUZ')), len(supheli)))
    if supheli:
        bulgu(u'KOSULSUZ yazili ama kaniti eksik',
              u'\n      '.join(supheli) +
              u'\n      "Kosulsuz" raporda "tartisilacak bir sey yok" diye okunur. '
              u'Kaniti eksik bir satir ya KOSULLU yapilmali ya da eksigi '
              u'giderilmeli.', BLOKE)



# --- 13 -----------------------------------------------------------------
def d13_urun_boyu(kok, yaz):
    """Tablodaki urun boyu, ciftin konsensusta URETTIGI boy mu.

    Urun boyu elle yazilabilen bir sayidir ve elle yazilan her sayi bir gun
    kayar. 2026-08-10: Proteolitik_Synergistaceae satirinda 173 bp yaziyordu,
    olcum 172 bp diyor (B-2-1197717 konsensusu, F@228 R@381, 381+19-228=172).
    Bir bazlik fark jel ayrimi hesabina ve bant sinifina girer.

    Olcut panelin kendi olcutu: mm<=1 ve 3' son iki baz tam.
    """
    import glob as _glob
    import re as _re
    pk = os.path.join(kok, 'primer_final',
                      'devir_ciftleri_20260802_sonrotus_TESLIM.tsv')
    kd = os.path.join(kok, 'konsensus_kanonik')
    if not os.path.exists(pk) or not os.path.isdir(kd):
        ATLANAN.append(u'13 urun boyu (panel ya da konsensus klasoru yok)')
        return
    # GLOB KULLANILMAZ. konsensus_kanonik klasorunde 250 dosya var ama yalniz
    # 100'u gecerli; kalan 150'si bagli klasorde silinemeyen kalinti ve bir
    # kismi ayni kutunun FARKLI icerikli eski surumu (33 kutuda olculdu).
    # Panelin kendi yukleyicisi (hedefler.konsensusler) bu yuzden INDEKS.tsv
    # okuyor; denetim de ayni kaynagi okumazsa panelin gormedigi bir diziyle
    # olcum yapar ve uydurma bir "sapma" uretir. 2026-08-10 gece: ilk surumum
    # tam bu hatayi yapiyordu.
    ixy = os.path.join(kd, 'INDEKS.tsv')
    if not os.path.exists(ixy):
        ATLANAN.append(u'13 urun boyu (konsensus_kanonik/INDEKS.tsv yok - '
                       u'kalinti dosyalarla olcum YAPILMAZ)')
        return
    kons = {}
    for r in _tsv(ixy):
        f = os.path.join(kd, (r.get('dosya') or '').strip())
        if not os.path.exists(f):
            continue
        kons[os.path.basename(f)] = ''.join(
            l.strip() for l in io.open(f, encoding='utf-8', errors='replace')
            if not l.startswith('>')).upper()
    if not kons:
        ATLANAN.append(u'13 urun boyu (indekste okunabilir dosya yok)')
        return

    def _rc(x):
        return x.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

    def _yer(p, d):
        n = len(p)
        out = []
        for i in range(len(d) - n + 1):
            f = 0
            for a, c in zip(p, d[i:i + n]):
                if a != c:
                    f += 1
                    if f > 1:
                        break
            else:
                if p[-2:] == d[i + n - 2:i + n]:
                    out.append(i)
        return out

    sat = [l.rstrip('\n').split('\t') for l in io.open(pk, encoding='utf-8')]
    b = sat[0]
    try:
        iH = b.index('Hedef')
        iU = b.index('Urun (bp)')
        iF = next(i for i, x in enumerate(b) if x.startswith('Ileri primer'))
        iR = next(i for i, x in enumerate(b) if x.startswith('Geri primer'))
    except (ValueError, StopIteration):
        ATLANAN.append(u'13 urun boyu (sutunlar taninmadi)')
        return
    sapan = []
    bakilan = 0
    for r in sat[1:]:
        if len(r) <= max(iU, iF, iR) or not r[iH].strip():
            continue
        F = r[iF].strip().upper()
        R = r[iR].strip().upper()
        if not _re.fullmatch(r'[ACGT]+', F or ''):
            continue
        try:
            u = int(_re.sub(r'\D', '', r[iU]))
        except ValueError:
            continue
        bakilan += 1
        Rrc = _rc(R)
        boylar = set()
        for d in kons.values():
            fs = _yer(F, d)
            if not fs:
                continue
            rs = _yer(Rrc, d)
            for i in fs:
                for j in rs:
                    if j >= i:
                        L = j + len(R) - i
                        if 40 <= L <= 600:
                            boylar.add(L)
        if boylar and u not in boylar:
            sapan.append(u'%s: tabloda %d bp, olculen %s'
                         % (r[iH].strip(), u, sorted(boylar)[:5]))
    yaz(u'  [13] urun boyu: %d cift olculdu, %d sapma' % (bakilan, len(sapan)))
    if sapan:
        bulgu(u'Tablodaki urun boyu olculenle TUTMUYOR',
              u'\n      '.join(sapan) +
              u'\n      Urun boyu jel ayrimi hesabina ve QuantiNova bant sinifina '
              u'girer; yanlis sayi iki hesabi da bozar.', DIKKAT)



# --- 14 -----------------------------------------------------------------
def d14_bat_dosyalari(kok, yaz):
    """.bat dosyalari CRLF ve saf ASCII mi, goto hedefleri var mi.

    2026-08-09: SINA_BAT.bat LF satir sonluydu ve on bes sinamanin on ucu
    SESSIZCE atlanmisti - goto hedefini bulamiyor, hata da vermiyor.
    2026-08-10: verification/one_key.py ve verification/full_chain.py da LF cikti.
    Turkce karakter de yorumlayiciyi bozar; bu yuzden saf ASCII sarti var.
    """
    import glob as _glob
    import re as _re
    sorun = []
    n = 0
    for y2 in sorted(_glob.glob(os.path.join(kok, '*.bat'))):
        n += 1
        b = io.open(y2, 'rb').read()
        ad = os.path.basename(y2)
        try:
            t = b.decode('ascii')
        except UnicodeDecodeError:
            sorun.append(u'%s: saf ASCII DEGIL (Turkce karakter yorumlayiciyi bozar)' % ad)
            continue
        if b.count(b'\n') and b.count(b'\r\n') != b.count(b'\n'):
            sorun.append(u'%s: LF satir sonu (%d satirin %d\'i CRLF). goto '
                         u'hedefini bulamayabilir ve HATA VERMEZ.'
                         % (ad, b.count(b'\n'), b.count(b'\r\n')))
        et = set(m.group(1).lower()
                 for m in _re.finditer(r'^:([a-z0-9_]+)', t, _re.M | _re.I))
        git = set(m.group(1).lower()
                  for m in _re.finditer(r'^[^\n]*?\bgoto\s+:([a-z0-9_]+)', t, _re.M | _re.I))
        eksik = sorted(x for x in git - et if x != 'eof')
        if eksik:
            sorun.append(u'%s: karsiligi olmayan goto hedefi: %s' % (ad, ', '.join(eksik)))
        # TEKRAR EDEN ETIKET. goto her zaman ILK etikete atlar; ayni ad iki kez
        # geciyorsa menunun bir tusu sessizce yanlis yere gider. 2026-08-11'de
        # PANEL.bat'a yeni tus eklerken tam bu oldu: dagitim satirinin altina
        # ikinci bir ":sk" dustu ve K tusu kutu planini degil, dagitimin
        # devamini calistirir hale geldi. Hicbir hata mesaji cikmadi.
        tekrar = {}
        for m in _re.finditer(r'^:([a-z0-9_]+)', t, _re.M | _re.I):
            k = m.group(1).lower()
            tekrar[k] = tekrar.get(k, 0) + 1
        cift = sorted(k for k, v in tekrar.items() if v > 1)
        if cift:
            sorun.append(u'%s: AYNI etiket birden cok kez tanimli: %s. goto ilk '
                         u'etikete atlar, o tus sessizce yanlis yere gider.'
                         % (ad, ', '.join(cift)))
    yaz(u'  [14] bat dosyalari: %d dosya, %d sorun' % (n, len(sorun)))
    if sorun:
        bulgu(u'.bat dosyasinda bicim sorunu',
              u'\n      '.join(sorun) +
              u'\n      Duzeltme: satir sonlarini CRLF yapin, dosyayi saf ASCII '
              u'kaydedin. Sessiz atlama bu yuzden olur.', BLOKE)



# --- 15 -----------------------------------------------------------------
def d15_konsensus_kalintilari(kok, yaz):
    """Kanonik konsensus klasorunde indekste OLMAYAN kalinti dosya var mi.

    2026-08-10 gece: klasorde 250 dosya var, INDEKS.tsv 100 tanesini
    tanimliyor. Kalan 150'si silinemeyen kalinti ve 33 kutuda ayni kutunun
    FARKLI icerikli iki-uc surumu duruyor (A1-1_2223.kanonik.fa ile
    A1-1_2223_kanonik.fasta ayni kutu, ayri dizi).

    Panelin kendi yukleyicisi indeks okudugu icin bundan etkilenmiyor. Ama
    glob yazan HER yeni betik sessizce yanlis diziyi secer - benim ilk urun
    boyu denetimim tam bunu yapti. Bu madde, tuzagin durdugunu hatirlatir.
    """
    import glob as _glob
    d = os.path.join(kok, 'konsensus_kanonik')
    ix = os.path.join(d, 'INDEKS.tsv')
    if not os.path.isdir(d):
        ATLANAN.append(u'15 konsensus kalintilari (klasor yok)')
        return
    if not os.path.exists(ix):
        bulgu(u'Kanonik konsensus INDEKSI YOK',
              u'%s bulunamadi. Indeks olmadan hangi dosyanin gecerli oldugu '
              u'bilinmez ve olcumler dosya adi sirasina kalir. Uretin: '
              u'python screening/build_canonical.py --kok .' % ix)
        return
    gecerli = set()
    for r in _tsv(ix):
        f = (r.get('dosya') or '').strip()
        if f:
            gecerli.add(f)
    hepsi = set(os.path.basename(f) for f in _glob.glob(os.path.join(d, '*.fa'))
                + _glob.glob(os.path.join(d, '*.fasta')))
    kalinti = hepsi - gecerli
    yaz(u'  [15] konsensus kalintilari: indekste %d, klasorde %d, kalinti %d'
        % (len(gecerli), len(hepsi), len(kalinti)))
    if kalinti:
        bulgu(u'Kanonik konsensus klasorunde KALINTI dosya',
              u'%d dosya indekste yok (klasorde %d, gecerli %d). Ornek: %s\n'
              u'      Panelin yukleyicisi indeks okudugu icin etkilenmiyor, ama '
              u'konsensus_kanonik/*.fa* diye GLOB yazan her betik ayni kutunun '
              u'eski surumunu secebilir. Yeni betik yazarken INDEKS.tsv okuyun.'
              % (len(kalinti), len(hepsi), len(gecerli),
                 ', '.join(sorted(kalinti)[:3])), BILGI)



# --- 16 -----------------------------------------------------------------
def d16_uyelik_kaynagi(kok, yaz):
    """P (tek protokol) ile K (kurtarma) AYNI uyelik dosyasini mi okuyor.

    Ikisi ayri dosya secerse dCq'lari ayri zeminde olcer ve karsilastirilamaz
    hale gelir. 2026-08-10: her ikisi de iki globu birlestirip a[-1] aliyordu;
    bu "en yeni" degil "alt klasor her zaman kazanir" demekti. Duzeltildi ama
    kural kodun disinda da sinanmali - iki betikten biri yarin degisebilir.
    """
    sec = {}
    for ad, dizin, mod in (('P', 'protocol', 'tek_protokol_olc'),
                           ('K', 'verification', 'kurtarma_turu')):
        d = os.path.join(kok, dizin)
        if not os.path.isdir(d):
            continue
        kod = ('import sys,os;sys.path.insert(0,%r);import %s as M;'
               'print(os.path.abspath(M.uyelik_dosyasi(%r) or ""))'
               % (d, mod, kok))
        try:
            p = subprocess.run([sys.executable, '-c', kod],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               timeout=120)
            sec[ad] = p.stdout.decode('utf-8', 'replace').strip()
        except Exception as e:
            ATLANAN.append(u'16 uyelik kaynagi (%s cagrilamadi: %s)' % (ad, e))
    if len(sec) < 2:
        ATLANAN.append(u'16 uyelik kaynagi (iki betik de okunamadi)')
        return
    ayni = len(set(v for v in sec.values() if v)) == 1
    yaz(u'  [16] uyelik kaynagi: P ve K ayni dosyayi mi okuyor -> %s'
        % ('evet' if ayni else 'HAYIR'))
    if not ayni:
        bulgu(u'P ile K FARKLI uyelik dosyasi okuyor',
              u'\n      '.join(u'%s: %s' % (k, v or 'bulunamadi')
                               for k, v in sorted(sec.items())) +
              u'\n      Ayri uyelikle olculen dCq degerleri ayni zeminde '
              u'degildir ve karsilastirilamaz.', BLOKE)



# --- 17 -----------------------------------------------------------------
def d17_veritabani_alfabesi(kok, yaz):
    """Taranan her FASTA DNA alfabesinde mi, "ikiz" denen dosyalar gercekten ikiz mi.

    Iki ayri tuzak, ikisi de yasanmis:
      1) SILVA RNA saklar (U). RNA alfabeli bir indeks DNA sorgularini HIC
         tutturmaz ve sonuc "0 hedef disi" yani TEMIZ gorunur. 2026-08-09'da
         SILVA indeksi tam bunu yapiyordu (3^9 k-mer alfabesinden anlasildi).
      2) identity_verification.py iki dosyayi "BAYT BAYT AYNI (cmp ile dogrulandi)"
         diye isaretlemisti. 2026-08-10 olcumu: SSU ciftinin boyutlari ayni
         ama ICERIKLERI farkli - 138.2 surumu U->T cevrilmis, ikizi hala RNA.
         Not, donusumden onceye aitti ve kimse geri donup bakmamisti.
    """
    import hashlib as _h
    d = os.path.join(kok, 'REFERANS_DB')
    if not os.path.isdir(d):
        ATLANAN.append(u'17 veritabani alfabesi (REFERANS_DB yok)')
        return
    # indeksi olan FASTA'lar taranan veritabanlaridir; alfabeleri DNA olmali
    rna = []
    bakilan = 0
    for f in sorted(os.listdir(d)):
        if not f.endswith(('.fna', '.fasta', '.fa')):
            continue
        y2 = os.path.join(d, f)
        if not any(os.path.exists(y2 + e) for e in ('.primerqc.bin', '.nsq')):
            continue          # indekssiz dosya taranmiyor, alfabesi onemsiz
        bakilan += 1
        u = t = 0
        with io.open(y2, encoding='utf-8', errors='replace') as fh:
            for i, l in enumerate(fh):
                if l.startswith('>'):
                    continue
                u += l.count('U')
                t += l.count('T')
                if i > 200000:
                    break
        if u:
            rna.append(u'%s: U=%d T=%d' % (f, u, t))
    # indeks FASTA'dan TAZE mi: FASTA sonradan degistiyse indeks eski veriyi
    # tarar ve bunu hicbir sey soylemez.
    bayat_ix = []
    for f in sorted(os.listdir(d)):
        if not f.endswith(('.fna', '.fasta', '.fa')):
            continue
        y2 = os.path.join(d, f)
        for e in ('.primerqc.bin', '.nsq'):
            ix = y2 + e
            if os.path.exists(ix) and os.path.getmtime(ix) < os.path.getmtime(y2) - 60:
                bayat_ix.append(u'%s%s indeksi FASTA\'dan ESKI' % (f, e))
    yaz(u'  [17] veritabani alfabesi: %d indeksli FASTA, %d RNA, %d bayat indeks'
        % (bakilan, len(rna), len(bayat_ix)))
    if bayat_ix:
        bulgu(u'Indeks FASTA\'dan ESKI',
              u'\n      '.join(bayat_ix) +
              u'\n      Indeks eski veriyi tarar; sonuc sessizce yanlis olur. '
              u'Yeniden indeksleyin.', BILGI)
    if rna:
        bulgu(u'Indeksli veritabani RNA alfabesinde',
              u'\n      '.join(rna) +
              u'\n      RNA alfabeli bir indeks DNA sorgularini tutturmaz ve '
              u'sonuc "0 hedef disi" yani TEMIZ gorunur. U->T cevirip yeniden '
              u'indeksleyin: bash build_index.sh <dosya>', BLOKE)

    # "ikiz" iddialari gercekten dogru mu
    ikizler = [('SILVA_138.2_SSURef_NR99.fasta', 'SILVA_SSURef_NR99.fasta'),
               ('SILVA_138.2_LSURef_NR99.fasta', 'SILVA_LSURef_NR99.fasta')]
    bozuk = []
    for a, b in ikizler:
        pa, pb = os.path.join(d, a), os.path.join(d, b)
        if not (os.path.exists(pa) and os.path.exists(pb)):
            continue
        if os.path.getsize(pa) != os.path.getsize(pb):
            bozuk.append(u'%s / %s: boyutlar farkli' % (a, b))
            continue
        ha, hb = _h.md5(), _h.md5()
        with io.open(pa, 'rb') as fa, io.open(pb, 'rb') as fb:
            ha.update(fa.read(20000000))
            hb.update(fb.read(20000000))
        if ha.hexdigest() != hb.hexdigest():
            bozuk.append(u'%s ile %s: boyut ayni ama ICERIK farkli '
                         u'(ilk 20 MB md5 tutmuyor)' % (a, b))
    if bozuk:
        bulgu(u'"Ikiz" denen dosyalar artik ikiz DEGIL',
              u'\n      '.join(bozuk) +
              u'\n      identity_verification.py bu dosyalari "bayt bayt ayni" diye '
              u'isaretliyor. Oylamaya girmiyorlar ama not YANLIS; birisi bayragi '
              u'cevirirse sessiz sifir uretir.', BILGI)



# --- 18 -----------------------------------------------------------------
def d18_evrensel_kapsam(kok, yaz):
    """Evrensel/kontrol primerleri KENDI hedeflerimizi goruyor mu.

    Evrensel primerlerde dCq TANIMSIZDIR (rakip kumesinin paydasi sifira
    gider); olcu KAPSAMDIR. Ama projede yazili bir kapsam esigi YOK ve
    verilen etiketler olculen degerlerle tutmuyor: 8 Agustos tablosunda
    Metanojen_universal'e %88 kapsamla "kapsam dusuk" denmis,
    Arke_universal'e %74 kapsamla hicbir sey denmemis.

    Daha onemlisi: 2026-08-11 olcumu, Arke_universal'in kendi uye
    kutularinin altisinda urun VERMEDIGINI gosteriyor ve bunlarin ucu
    Nitrosocosmicus (panelin KENDI hedefi), biri Methanomassiliicoccus
    (metilotrofik metanojen hedefi). Normallestirme kontrolu, normalize
    edecegi hedefi gormuyor.

    Bu madde hukum vermez; sayilari yan yana koyar ve yazili bir olcut
    olmadigini soyler. Esigi insan koyar.
    """
    import csv as _csv
    import re as _re
    p = _tsv(os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'panel_tek_protokol.tsv'))
    sl = {(r.get('hedef') or '').strip(): r
          for r in _tsv(os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'SIPARIS_LISTESI.tsv'))}
    if not p or not sl:
        ATLANAN.append(u'18 evrensel kapsam (tablo yok)')
        return
    satir = []
    for r in p:
        ad = (r.get('hedef') or '').strip()
        if (sl.get(ad, {}).get('SINIF') or '').strip().upper() != 'EVRENSEL':
            continue
        k = (r.get('ASIL_kapsam_mm1') or '').strip()
        m = _re.match(r'(\d+)\s*/\s*(\d+)', k)
        oran = (100.0 * int(m.group(1)) / int(m.group(2))) if m and int(m.group(2)) else None
        satir.append((ad, k, oran))
    if not satir:
        ATLANAN.append(u'18 evrensel kapsam (evrensel cift bulunamadi)')
        return
    dusuk = [(a, k, o) for a, k, o in satir if o is not None and o < 90.0]
    yaz(u'  [18] evrensel kapsam: %d cift | %s'
        % (len(satir), ', '.join(u'%s %s' % (a.split('_')[0][:12], k) for a, k, _o in satir)))
    if dusuk:
        bulgu(u'Evrensel/kontrol primerinde kapsam %90 altinda',
              u'\n      '.join(u'%s: %s (%%%.0f)' % (a, k, o) for a, k, o in dusuk) +
              u'\n      Evrensel primerlerde olcu KAPSAMDIR ve projede yazili bir '
              u'kapsam esigi YOKTUR; buradaki %90 bu denetimin koydugu gecici '
              u'siniridir, panelin kurali degildir. Kontrol primeri kendi '
              u'hedeflerini gormuyorsa normallestirme yanlidir. Olcutu yazin ya '
              u'da bu satirlar icin raporda gerekce verin.', DIKKAT)


def d19_uyelik_icerigi(kok, yaz):
    u"""IKI UYELIK KAYNAGI AYNI KUTULARI MI SOYLUYOR.

    16. madde iki betigin AYNI DOSYAYI okudugunu sinar. Ama proje iki AYRI
    uyelik dosyasi kullaniyor ve bu bilerek boyle:
        screening/hedef_uyelik.tsv          - arama/tarama tarafi
        uyelik_yeniden_turetme_uyelik_*.tsv      - tek protokol olcumu
    Ikisi ayni dosya olmadigi icin 16. madde bu ikisini hic karsilastirmaz.
    Icerikleri sessizce ayrisirsa arama bir kumeyi hedef sanip optimize eder,
    olcum baska bir kumeye gore not verir.

    2026-08-11'de tam bu oldu: Petriella_cinsi icin arama tarafi F2-4_500148'i
    UYE sayiyordu (taxid 500148 uzerinden), olcum tarafi ayni kutuyu RAKIP
    sayiyordu - ve o kutu cifti dusuren karar kutusuydu (1655/3000 okuma).
    Yani lokus taramasi, olcumun ISTEDIGININ TERSINI optimize etti. On kadar
    hedefte buna benzer ayrisma vardi.

    Kural: fark varsa ve farktaki kutunun OLCULEN kimligi varsa -> SIPARISI
    DURDURUR. Kimligi hic olculmemis kutulardan kaynaklanan fark -> DIKKAT
    (once o kutu olculmeli, karar ondan sonra verilir).
    """
    import glob as _glob
    uy = [x for x in _glob.glob(os.path.join(kok, 'uyelik_yeniden_turetme_uyelik_*.tsv'))
          if '.yedek' not in x]
    if not uy:
        ATLANAN.append(u'19 uyelik icerigi (tek protokol uyelik dosyasi yok)')
        return
    uy.sort(key=lambda p: (os.path.getmtime(p), os.path.basename(p)))
    kod = (
        'import sys,os,io,csv\n'
        'sys.path.insert(0,%r)\n'
        'from screening import targets as H\n'
        'sat=[l.rstrip("\\n").split("\\t") for l in io.open(%r,encoding="utf-8")]\n'
        'bas=sat[0]; iu=bas.index("yeni_uye_kutular")\n'
        'TP={r[0].strip():set(x for x in r[iu].split(";") if x) '
        'for r in sat[1:] if r and r[0].strip()}\n'
        'panel,_=H.panel_oku(); kons=H.konsensusler(); kut=H.kutular()\n'
        'var={k["kutu"] for k in kut}\n'
        # KIMLIGI OLCULMUS kutular IKI kaynaktan toplanir. Yalniz G asamasi
        # tablosuna bakmak yaniltiyordu: F1-*_44689 kutulari o tabloda "hicbir
        # hedefin uyesi degil" diye ATLANMIS, ama bolluk calismasi (11 Agustos)
        # dordunu de olcmus (Zoopagomycota / Nucletmycea). Tek tabloya bakan
        # denetim "olculmemis" deyip yanlis onem veriyordu.
        'kim=set()\n'
        'y=os.path.join(%r,"TUM_KIMLIK_SONUC","tum_kutu_kimlikleri.tsv")\n'
        'if os.path.exists(y):\n'
        '    s2=[l.rstrip("\\n").split("\\t") for l in io.open(y,encoding="utf-8")]\n'
        '    bi=[i for i,r in enumerate(s2) if r and r[0].strip()=="kutu"][0]\n'
        '    kim={r[0].strip() for r in s2[bi+1:] if r and r[0].strip() '
        'and not r[0].startswith("#")}\n'
        'import glob as _g\n'
        'for y2 in _g.glob(os.path.join(os.path.dirname(y),"..","BOLLUK_OLCULEN_*",'
        '"karsilastirma_kutu.tsv")):\n'
        '    with io.open(y2,encoding="utf-8") as fh2:\n'
        '        for r2 in csv.DictReader((l for l in fh2 if not l.startswith("#")),'
        'delimiter="\\t"):\n'
        '            if (r2.get("YENI_olculen_kimlik") or "").strip():\n'
        '                kim.add((r2.get("kutu") or "").strip())\n'
        'for p in panel:\n'
        '    ad=p["hedef"]\n'
        '    b=H.hedef_baglami(p,kons=kons,kut=kut)\n'
        '    ka=set(k["kutu"] for k in b["uye_kutu"])\n'
        '    tp=TP.get(ad) or TP.get(H.AD_ESLEME.get(ad,""))\n'
        '    if tp is None: continue\n'
        '    tp={x for x in tp if x in var}\n'
        '    if ka!=tp:\n'
        '        f=sorted((ka-tp)|(tp-ka))\n'
        '        olculen=[x for x in f if x in kim]\n'
        '        print("%%s\\t%%s\\t%%s" %% (ad,",".join(f),",".join(olculen)))\n'
        % (kok, uy[-1], kok))
    try:
        p = subprocess.run([sys.executable, '-c', kod], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, timeout=600)
    except Exception as e:
        ATLANAN.append(u'19 uyelik icerigi (%s)' % e)
        return
    if p.returncode != 0:
        ATLANAN.append(u'19 uyelik icerigi (karsilastirma kosmadi: %s)'
                       % p.stderr.decode('utf-8', 'replace').strip()[-200:])
        return
    satirlar = [l for l in p.stdout.decode('utf-8', 'replace').splitlines() if l.strip()]
    yaz(u'  [19] uyelik icerigi: iki kaynagin ayristigi hedef -> %d' % len(satirlar))
    for l in satirlar:
        pr = l.split('\t')
        ad, fark = pr[0], pr[1]
        olculen = pr[2] if len(pr) > 2 else ''
        if olculen:
            bulgu(u'%s: iki uyelik kaynagi ayni kutularda anlasmiyor' % ad,
                  u'Ayrisan kutu: %s\n      Bunlarin OLCULEN kimligi var (%s), '
                  u'yani karar verilebilir bir veri duruyor ve iki taraf yine de '
                  u'farkli sayiyor. Arama bir kumeyi optimize ederken olcum baska '
                  u'kumeye not veriyor; dCq degerleri ayni zeminde degil.'
                  % (fark, olculen), BLOKE)
        else:
            bulgu(u'%s: uyelik ayrisiyor, ayrisan kutularin kimligi OLCULMEMIS' % ad,
                  u'Ayrisan kutu: %s\n      Bu kutularin olculen kimligi yok; '
                  u'once kimlikleri olculmeli, uye mi rakip mi ondan sonra '
                  u'yazilmali. Kraken etiketine gore karar VERILMEZ.' % fark,
                  DIKKAT)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', '--kok', dest='kok', default='.')
    p.add_argument('--offline', '--agsiz', dest='agsiz', action='store_true',
                   help='NCBI gerektiren denetimleri atla')
    p.add_argument('--generate', '--uretilecek', dest='uretilecek', default='',
                   help='bu kosuda yeniden uretilecek cikti yollari (virgulle); '
                        'bunlara tazelik denetimi uygulanmaz')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)

    def yaz(s=''):
        print(s, flush=True)

    yaz(u'=' * 78)
    yaz(u'  HER KOSUDA DENETIM   %s' % time.strftime('%Y-%m-%d %H:%M'))
    yaz(u'  Olcum yapilmaz, dosya degistirilmez. Yalnizca bakilir.')
    yaz(u'=' * 78)

    d1_harita_anahtarlari(kok, yaz)
    d3_referans_bayat(kok, yaz)
    d4_kaynak_tutarliligi(kok, yaz)
    d5_muhur_diziyi_iceriyor_mu(kok, yaz)
    d7_geometri_kapisi(kok, yaz)
    d8_plaka_jel_ve_bant(kok, yaz)
    d9_belgelerde_bayat_sayi(kok, yaz)
    d10_ad_kurali_sinavi(kok, yaz)
    d11_siparis_dizileri(kok, yaz)
    d12_kanitsiz_kosulsuz(kok, yaz)
    d13_urun_boyu(kok, yaz)
    d14_bat_dosyalari(kok, yaz)
    d15_konsensus_kalintilari(kok, yaz)
    d16_uyelik_kaynagi(kok, yaz)
    d17_veritabani_alfabesi(kok, yaz)
    d18_evrensel_kapsam(kok, yaz)
    d19_uyelik_icerigi(kok, yaz)
    d6_cikti_tazeligi(kok, yaz, tuple(x.strip() for x in a.uretilecek.split(',') if x.strip()))
    d2_kapsama(kok, yaz, a.agsiz)

    yaz('')
    if BULGU:
        sirali = sorted(BULGU, key=lambda x: (BLOKE, DIKKAT, BILGI).index(x[2]))
        say = {}
        for _b, _a, o in BULGU:
            say[o] = say.get(o, 0) + 1
        yaz(u'  %d BULGU  (%s)' % (len(BULGU), ', '.join(
            u'%s %d' % (o, say[o]) for o in (BLOKE, DIKKAT, BILGI) if o in say)))
        for b, ay, o in sirali:
            yaz(u'   [%s] %s' % (o, b))
            yaz(u'      %s' % ay)
    else:
        yaz(u'  Butun denetimler temiz.')
    if ATLANAN:
        yaz(u'  Atlanan denetim: %s' % '; '.join(ATLANAN))
    yaz(u'=' * 78)

    rapor = os.path.join(kok, 'TEK_TUS_SONUC')
    if os.path.isdir(rapor):
        with io.open(os.path.join(rapor, 'DENETIM_RAPORU.md'), 'w',
                     encoding='utf-8', newline='') as fh:
            fh.write(u'# Her koşuda denetim\n\nÜretim: %s\n\n'
                     % time.strftime('%Y-%m-%d %H:%M'))
            if BULGU:
                say = {}
                for _b, _a, o in BULGU:
                    say[o] = say.get(o, 0) + 1
                fh.write(u'## %d bulgu\n\n' % len(BULGU))
                for o in (BLOKE, DIKKAT, BILGI):
                    grup = [x for x in BULGU if x[2] == o]
                    if not grup:
                        continue
                    fh.write(u'### %s (%d)\n\n' % (o, len(grup)))
                    for b, ay, _o in grup:
                        fh.write(u'- **%s** — %s\n' % (b, ay.replace('\n', ' ')))
                    fh.write(u'\n')
            else:
                fh.write(u'## Bütün denetimler temiz\n\n')
            if ATLANAN:
                fh.write(u'\n## Atlananlar\n\n')
                for x in ATLANAN:
                    fh.write(u'- %s\n' % x)

    if any(o == BLOKE for _b, _a, o in BULGU):
        return 1
    if BULGU:
        return 3          # yalniz DIKKAT/BILGI - kosu durdurulmaz
    return 2 if ATLANAN else 0


if __name__ == '__main__':
    sys.exit(main())
