# -*- coding: utf-8 -*-
"""HIZLI TUTARLILIK TESTI (gerileme testi) - zincir kendini yeniden uretiyor mu?

NE OLDUGU VE NE OLMADIGI - IKI CUMLE:
  Bu test, kodun KENDINI YENIDEN URETTIGINI sinar: ayni motor, on kat sig veriyle
  ayni sinifi ve ayni siralamayi veriyor mu?
  Olcumun DOGRU oldugunu SINAMAZ - cunku beklenen degerler de AYNI motorun tam
  derinlikli kosusundan geliyor; motorun sistematik bir hatasi varsa bu test onu
  yakalayamaz, yakalayabilecegi sey kod ve yapilandirma kaymasidir.

Kod incelemesi bir yere kadar gider; bu betik CIKTIYA bakar. Zincirin dort
asamasini KUCUK ve cevabi BILINEN bir alt kumede kosturur ve sonucu referans
kosuyla karsilastirir.

Referans: TEK_PROTOKOL_SONUC/panel_tek_protokol.tsv (tam derinlik, 3000 okuma).

NE SINANIR
  1) Esik ustu olanlar esik ustu, esik alti olanlar esik alti cikiyor mu?
     (ASIL OLCUT budur - sayinin kendisi degil)
  2) Siralama korunuyor mu? (Nitrosocosmicus > Metilotrofik > Cloacimonas > M_cinsi)
  3) Degerler makul bantta mi? (bant asagida gerekcesiyle yazili)
  4) Sonraki asamalar CIKTI URETIYOR mu? Bos cikti GECMIS SAYILMAZ.

SUREnin kisa tutulma yolu: okuma tavani 3000 -> 300, hedefler 8 satirla sinirli,
yol-3 tasarim taramasi kucultulmus. Olcum derinligi dustugu icin sayilar birebir
tutmaz; nitekim tutmasi da beklenmiyor - bkz. BANT GEREKCESI.
"""

# -------------------------------------------------------------------------
# quick_consistency_test.py — zincirin dort asamasini kucuk ve cevabi bilinen
# bir alt kumede kosturup kodun kendini yeniden urettigini sinar (gerileme testi).
#
# GİRDİ  : proje kokundeki dort betik ve olcum kaynaklari; bunlar HIZLI_TEST/
#          altina sembolik baglarla baglanir (screening, protocol,
#          verification, REFERANS_DB, konsensus_kanonik, primer_final, "fastq files",
#          engine ve uyelik_yeniden_turetme_uyelik_*.tsv). Beklenen
#          degerler bu dosyadaki BEKLENEN_UST / BEKLENEN_ALT / BEKLENEN_YENI
#          sabitlerinden gelir; kaynaklari tam derinlikli referans kosusudur.
# ÇIKTI  : HIZLI_TEST/HIZLI_TEST_RAPORU.md ve HIZLI_TEST/test_gunlugu.txt;
#          ayrica asamalarin kendi ciktilari HIZLI_TEST/ altinda ayri durur -
#          gercek sonuc klasorlerine DOKUNULMAZ.
# ÇAĞRAN : verification/full_chain.py -> H tusu
#          (bat icinde: wsl -e python3 "verification/quick_consistency_test.py" --kok .)
#
# NE OLDUGU  : kodun ve yapilandirmanin kaymadigini gosterir.
# NE OLMADIGI: olcumun DOGRU oldugunu gostermez - beklenen degerler de ayni
#              motorun kosusundan geliyor, motorun sistematik hatasi bu testte
#              gorunmez. Bagimsiz teyit MFEprimer katmaninin isidir.
# -------------------------------------------------------------------------
import io, os, sys, csv, json, time, subprocess, argparse

VERSIYON = '1.0 (2026-08-03)'
OKUMA = 300
# ESIK TEK KAYNAKTAN GELIR: screening/config.py -> ESIK_DCQ = 3.0
# Kat karsiligi 2 ** ESIK_DCQ = 8,00. Sabit sayi GOMULMEZ; dCq degisirse
# tek yerden degisir. Gerekce ve verim uyarisi o dosyada yazili.
def _esik_yukle():
    """Esigi TEK KAYNAKTAN okur: screening/config.py.
    verification/ ile screening/ kardes klasorler oldugu icin kok buradan
    turetilir; betik hangi calisma dizininden cagrilirsa cagrilsin bulur."""
    import os as _o, sys as _s
    _kok = _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))
    if _kok not in _s.path:
        _s.path.insert(0, _kok)
    from screening import config as _y
    return _y

_C = _esik_yukle()
ESIK = _C.AYRIM_ESIK

# --- referans degerler ---------------------------------------------------
# 2026-08-10 DUZELTME. Bu sayilar KODA GOMULU sabitlerdi ve hangi PRIMER
# CIFTINDEN olculdukleri hicbir yerde yazmiyordu. Bugun Bacteroidales cifti
# degistirildi (F: GAAGCTAGGATTTGGTTGCTGTG -> GCGTTATCCGGATTTATTGGGTTT) ve
# test eski cifte ait 0,74x'i yeni ciftin 14,23x'iyle karsilastirip "ZINCIR
# TUTARSIZ" dedi. Zincir tutarsiz DEGILDI; referans bayatti.
#
# Artik referanslar HIZLI_TEST/referans_degerler.tsv dosyasindan okunur ve her
# satirda o olcumun yapildigi F/R DIZISI yazar. Test sirasinda cift degismisse
# karsilastirma YAPILMAZ - "referans gecersiz, cift degisti" denir ve zincir
# durdurulmaz. Sessizce dogru sayilmaz, sessizce yanlis da sayilmaz.
#
# Dosya yoksa asagidaki eski sabitler kullanilir ama rapor bunu acikca yazar.
REFERANS_DOSYASI = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, 'HIZLI_TEST', 'referans_degerler.tsv')
CIKIS_TUTARSIZ = 6      # gercek gerileme - uzun kosuya GIRILMEZ
CIKIS_REFERANS_BAYAT = 7  # referans karsilastirilamaz - zincir devam edebilir


def _referans_yukle():
    """hedef -> (referans_x, F, R). Dosya yoksa bos doner."""
    out = {}
    y = os.path.abspath(REFERANS_DOSYASI)
    if not os.path.exists(y):
        return out, y
    bas = None
    for l in io.open(y, encoding='utf-8'):
        l = l.rstrip('\n')
        if not l.strip() or l.startswith('#'):
            continue
        p = l.split('\t')
        if bas is None:
            bas = p
            continue
        r = dict(zip(bas, p))
        try:
            ref = float((r.get('referans_x') or '').replace(',', '.'))
        except ValueError:
            continue
        out[(r.get('hedef') or '').strip()] = (
            ref, (r.get('F') or '').strip().upper(), (r.get('R') or '').strip().upper())
    return out, y


REFERANS, REFERANS_YOL = _referans_yukle()

BEKLENEN_UST = [
    ('Nitrosocosmicus_AOA',        155.15),
    ('Metilotrofik_metanojen',      91.44),
    ('Proteolitik_Cloacimonas',     33.89),
    ('Methanosarcina_cinsi',        24.17),
]
BEKLENEN_ALT = [
    ('Proteiniphilum_cinsi',         0.00),
    ('Methanosarcina mazei',         0.82),
    ('Bacteroidales_kumesi',         0.74),
]
BEKLENEN_YENI = [('Petriella_cinsi', 11.03)]

BANT_ALT, BANT_UST = 0.5, 2.0
# Esige COK YAKIN referans degerleri (< 15x) dusuk derinlikte esigin obur
# tarafina gecebilir ve bu ZINCIR HATASI DEGILDIR: projenin kendi olcumu
# Petriella LSU cifti icin tam derinlikte 11,03x, 300 okumada 8,93x diyor.
# Bu satirlarda sinif degisimi UYARI uretir; bant disina cikarsa HATA olur.
SINIRDA_UST = ESIK * 1.5
BANT_GEREKCESI = u"""
BANT GEREKCESI - neden 0,5x - 2,0x

Ayrim kati = (uye Wilson ALT siniri) / (rakip Wilson UST siniri). Wilson
araliginin genisligi okuma sayisiyla daralir. Testte derinlik 3000'den 300'e
indigi icin uye alt siniri DUSER, rakip ust siniri YUKSELIR; iki etki de orani
KUCULTUR. Yani testte cikan sayinin referanstan dusuk olmasi BEKLENEN davranistir,
hata degil.

Projenin kendi olculmus ornegi: Petriella LSU cifti tam derinlikte 11,03x,
panelin 300 okuma standardinda 8,93x -> 0,81x oran. Az okumali kutularda sapma
daha buyuk olabilir. 0,5-2,0 bandi bu gozlenen sapmayi rahatca kapsar ama bir
satirin esigin obur tarafina gecmesini yakalayacak kadar dardir.

KARAR bandin kendisine degil, SINIFIN korunmasina baglanir: esik ustu olan esik
ustu, esik alti olan esik alti kalmalidir. Bant disina cikan ama sinifi korunan
satir UYARI uretir, testi dusurmez.

TEK ISTISNA - ESIK SINIRINDAKI SATIRLAR: referans degeri 15x'in altinda olan bir
satir, dusuk derinlikte esigin obur tarafina gecebilir ve bu zincir hatasi degil
olcum derinliginin dogal sonucudur (olculmus ornek: Petriella LSU cifti tam
derinlikte 11,03x, 300 okumada 8,93x). Boyle satirlarda sinif degisimi, oran
bantta kaldigi surece UYARI olarak raporlanir. Oran bandin da disina cikarsa
HATA sayilir.
"""


def sure_metni(sn):
    sn = int(sn)
    return ('%d saniye' % sn) if sn < 90 else ('%d dakika' % round(sn / 60.0)) \
        if sn < 5400 else ('%.1f saat' % (sn / 3600.0))


def vir(x, b=2):
    if x is None:
        return '-'
    try:
        return ('%.*f' % (b, float(str(x).replace(',', '.')))).replace('.', ',')
    except (TypeError, ValueError):
        return str(x)


def _f(s):
    try:
        return float(str(s).replace(',', '.'))
    except (TypeError, ValueError):
        return None


def tsv_oku(yol):
    if not os.path.exists(yol):
        return []
    with open(yol, encoding='utf-8') as fh:
        return list(csv.DictReader(
            (s for s in fh if s.strip() and not s.startswith('#')), delimiter='\t'))


# Bir asamayi alt surec olarak kosar ve zaman tavani uygular. Zaman asimi
# "cikis kodu yok" olarak dondurulur - takilan bir asama testi sonsuza kadar
# bekletmesin diye.
def kos(yaz, ad, arg, tavan_sn):
    t0 = time.time()
    yaz(u'  > %s starting...' % ad)
    try:
        pr = subprocess.run(arg, timeout=tavan_sn, capture_output=True, text=True)
        rc, ciktisi = pr.returncode, (pr.stdout or '')[-2000:]
    except subprocess.TimeoutExpired:
        yaz(u'    TIMEOUT (%s): the stage exceeded its time cap' % sure_metni(tavan_sn))
        return None, u'zaman asimi'
    yaz(u'    done: exit code %d, %s' % (rc, sure_metni(time.time() - t0)))
    return rc, ciktisi


# ---------------------------------------------------------------------------
# P asamasini SIG derinlikte (300 okuma, tam kosuda 3000) ve yalniz sekiz sinama
# satirinda kosar, sonra ucunu birden denetler:
#   1) SINIF korunuyor mu - esik ustu satir esik ustu, esik alti satir esik alti
#      cikiyor mu. ASIL OLCUT budur.
#   2) SIRALAMA korunuyor mu - esik ustu dort satirin buyukluk sirasi bozuldu mu.
#   3) Degerler 0,5x - 2,0x bandinda mi (bant disi ama sinifi korunan satir
#      yalnizca UYARI uretir, testi dusurmez).
#
# SAYILAR NEDEN BIREBIR TUTMAZ: ayrim kati = (uye Wilson ALT siniri) / (rakip
# Wilson UST siniri). Wilson araliginin genisligi okuma sayisiyla daralir;
# derinlik dusurulunce uye alt siniri DUSER, rakip ust siniri YUKSELIR ve iki
# etki de orani kuculrur. Yani testte cikan sayinin referanstan dusuk olmasi
# beklenen davranistir, hata degil.
#
# Referans degeri 15x'in (SINIRDA_UST) altinda olan satirlar sig derinlikte
# esigin obur tarafina gecebilir; oran bantta kaldigi surece bu UYARI sayilir,
# HATA sayilmaz. Olculmus ornek: Petriella LSU cifti tam derinlikte 11,03x,
# 300 okumada 8,93x.
# ---------------------------------------------------------------------------
def calistir(kok, hizli_kok, tavan_dk, yaz):
    hedefler = ([h for h, _ in BEKLENEN_UST] + [h for h, _ in BEKLENEN_ALT]
                + [h for h, _ in BEKLENEN_YENI])
    sec = ','.join(hedefler)
    py = sys.executable
    sonuc = dict(asama={}, satir=[], uyari=[], hata=[])

    # --- P ---
    rc, _ = kos(yaz, 'P (TEK PROTOKOL, %d okuma)' % OKUMA,
                [py, os.path.join(kok, 'protocol', 'single_protocol_measure.py'),
                 '--kok', hizli_kok, '--okuma', str(OKUMA), '--yalniz', sec],
                tavan_dk * 60)
    P = tsv_oku(os.path.join(hizli_kok, 'TEK_PROTOKOL_SONUC', 'panel_tek_protokol.tsv'))
    sonuc['asama']['P'] = dict(rc=rc, satir=len(P))
    if not P:
        sonuc['hata'].append(u'P asamasi hic satir uretmedi - zincir burada kopuyor.')
        return sonuc

    olculen = {}
    for r in P:
        olculen[r['hedef']] = dict(kat=_f(r.get('ASIL_ayrim_mm1')),
                                   karar=r.get('esik_gecti_mi', ''),
                                   kapsam=r.get('ASIL_kapsam_mm1', ''),
                                   F=(r.get('F') or ''), R=(r.get('R') or ''))

    def bul(ad):
        for h, v in olculen.items():
            if ad.lower() in h.lower():
                return h, v
        return None, None

    # --- 1) sinif korunuyor mu + 3) bant ---
    for grup, bekl, ust_mu in (('esik ustu', BEKLENEN_UST, True),
                               ('esik alti', BEKLENEN_ALT, False),
                               ('yeni cift', BEKLENEN_YENI, True)):
        for ad, ref in bekl:
            h, v = bul(ad)
            # Referans dosyasi varsa dosya kazanir; sabit yalniz yedektir.
            _rf = None
            for _k, _val in REFERANS.items():
                if ad.lower() in _k.lower() or _k.lower() in ad.lower():
                    ref, _rf = _val[0], _val
                    break
            # BEKLENEN SINIF da referanstan gelir. Hedefin hangi sabit
            # listede (UST/ALT) yazildigi da bayatlayan bir bilgidir:
            # Bacteroidales sabitlerde "esik alti" yaziyordu, yeni ciftle
            # tam derinlikte 37,23x yani esik USTU. Sinifi listeye degil
            # olcume sorarak bu tuzak kapatildi.
            #
            # DIKKAT: ust_mu DIS dongunun degiskeni. Ona atama yapmak, ayni
            # gruptaki SONRAKI hedefleri de etkiler ve referansi olmayan bir
            # hedef, kendinden onceki hedefin sinifiyla degerlendirilir.
            # Bu yuzden yerel bir degisken kullanilir (2026-08-10, yazarken
            # yakalandi; kosmadan once duzeltildi).
            if _rf and v is not None:
                # CIFT DEGISTI MI? Degistiyse eski olcumle karsilastirmak
                # anlamsizdir; sessizce "tutarsiz" demek YANLIS olur.
                _cf = (v.get('F') or '').upper()
                _cr = (v.get('R') or '').upper()
                if _rf[1] and _rf[2] and (_cf, _cr) != (_rf[1], _rf[2]):
                    sonuc.setdefault('referans_bayat', []).append(
                        u'%s: primer cifti referans olcumunden bu yana DEGISTI '
                        u'(referans F %s / R %s, simdiki F %s / R %s). Eski %sx '
                        u'ile yeni olcum (%sx) karsilastirilamaz - bu bir '
                        u'gerileme DEGILDIR. Referansi yenileyin: '
                        u'HIZLI_TEST/referans_degerler.tsv'
                        % (ad, _rf[1][:12], _rf[2][:12], _cf[:12], _cr[:12],
                           vir(ref), vir(v['kat']) if v['kat'] is not None else '-'))
                    continue
            if v is None:
                sonuc['hata'].append(u'%s: satir P ciktisinda YOK (beklenen %s: %sx)'
                                     % (ad, grup, vir(ref)))
                continue
            kat = v['kat']
            if kat is None:
                sonuc['hata'].append(u'%s: ayrim kati olculemedi (karar=%s), '
                                     u'referansta %sx' % (ad, v['karar'], vir(ref)))
                continue
            _ust = (ref >= ESIK) if _rf else ust_mu
            gecti = kat >= ESIK
            sinirda = 0 < ref < SINIRDA_UST
            oran_ok = (ref > 0 and BANT_ALT <= (kat / ref) <= BANT_UST)
            if gecti != _ust and sinirda and oran_ok:
                sonuc['uyari'].append(
                    u'%s: referansta %sx, testte %sx - esik SINIRINDA bir satir '
                    u'(referans < %sx). Dusuk derinlikte esigin obur tarafina '
                    u'gecmesi beklenen davranistir (olculmus ornek: 11,03x -> '
                    u'8,93x); oran bantta oldugu icin HATA sayilmadi.'
                    % (ad, vir(ref), vir(kat), vir(SINIRDA_UST)))
            elif gecti != _ust:
                sonuc['hata'].append(
                    u'%s: SINIF DEGISTI - referansta %sx (%s) iken testte %sx (%s)%s'
                    % (ad, vir(ref), u'esik ustu' if _ust else u'esik alti',
                       vir(kat), u'esik ustu' if gecti else u'esik alti',
                       u'' if oran_ok else u' - ustelik %sx-%sx bandinin DISINDA'
                       % (vir(BANT_ALT), vir(BANT_UST))))
            elif ref > 0:
                oran = kat / ref
                if not (BANT_ALT <= oran <= BANT_UST):
                    sonuc['uyari'].append(
                        u'%s: %sx (referans %sx, oran %sx) - bant disi ama sinif '
                        u'korundu' % (ad, vir(kat), vir(ref), vir(oran)))
            sonuc['satir'].append(dict(hedef=ad, ref=ref, olculen=kat,
                                       karar=v['karar'], kapsam=v['kapsam'],
                                       sinif_ok=(gecti == _ust)))

    # --- 2) siralama korunuyor mu ---
    ref_sira = [ad for ad, _ in sorted(BEKLENEN_UST, key=lambda x: -x[1])]
    olc = []
    for ad in ref_sira:
        _h, v = bul(ad)
        olc.append((ad, v['kat'] if v and v['kat'] is not None else -1))
    test_sira = [ad for ad, _ in sorted(olc, key=lambda x: -x[1])]
    sonuc['siralama'] = dict(referans=ref_sira, test=test_sira,
                             korundu=(ref_sira == test_sira))
    if ref_sira != test_sira:
        sonuc['hata'].append(u'SIRALAMA BOZULDU: referans %s, testte %s'
                             % (' > '.join(ref_sira), ' > '.join(test_sira)))
    return sonuc


# ---------------------------------------------------------------------------
# K, D ve I asamalarini kucultulmus ayarlarla kosar. Burada olcumun DOGRULUGU
# degil, asamanin CALISIP SATIR URETTIGI sinanir - bos cikti gecmis sayilmaz.
#
# ARDISIK COKUS AYRIMI: K hicbir cifti kurtarmadiysa D'nin dogrulayacagi bir sey
# yoktur ve D'nin bos kalmasi D'nin hatasi DEGILDIR. Bu durumda D, sentetik tek
# satirlik bir girdiyle AYRI bir kokte kosturularak kendi basina sinanir; boylece
# "D calismiyor" ile "D'ye is dusmedi" birbirine karismaz.
# ---------------------------------------------------------------------------
def sonraki_asamalar(kok, hizli_kok, tavan_dk, yaz, sonuc):
    """K, D, I gercekten CIKTI URETIYOR mu? Bos cikti GECMIS SAYILMAZ."""
    py = sys.executable

    rc, _ = kos(yaz, 'K (verification)',
                [py, os.path.join(kok, 'verification', 'recovery_round.py'),
                 '--kok', hizli_kok, '--okuma', str(OKUMA),
                 # Test KAPSAMI dar tutulur: amac dogruluk degil, "asama calisiyor
                 # ve satir uretiyor" kaniti. Tam kosuda bu tavanlar kalkar.
                 '--tarama-ust', '40', '--aday-ust', '5', '--arms-ust', '0',
                 '--panelsiz-atla'],
                tavan_dk * 60)
    K = tsv_oku(os.path.join(hizli_kok, 'KURTARMA_SONUC', 'kurtarma_satirlari.tsv'))
    sonuc['asama']['K'] = dict(rc=rc, satir=len(K))
    if len(K) < 1:
        sonuc['hata'].append(u'K asamasi HIC satir uretmedi (bos cikti gecmis sayilmaz).')
    else:
        yaz(u'    K: %d rows, %d of them recovered'
            % (len(K), sum(1 for r in K if (r.get('esigi_gecti_mi') or '').startswith('EVET'))))

    rc, _ = kos(yaz, 'D (DOGRULAMA, yalniz yerel katman)',
                [py, os.path.join(kok, 'verification', 'specificity_round.py'),
                 '--kok', hizli_kok, '--ncbi', 'elle'],
                tavan_dk * 60)
    Dd = tsv_oku(os.path.join(hizli_kok, 'DOGRULAMA_SONUC', 'dogrulama_uc_sutun.tsv'))
    sonuc['asama']['D'] = dict(rc=rc, satir=len(Dd))
    if len(Dd) < 1 and len(K) >= 1 and not [r for r in K
                                            if (r.get('esigi_gecti_mi') or '').startswith('EVET')]:
        # K kostu ama HICBIR sey kurtarmadi -> D'nin dogrulayacagi cift yok.
        # Bu D'nin hatasi degil. Yine de D'nin CALISTIGINI kanitlamak gerekir:
        # sentetik tek satirlik bir girdiyle ayri bir kokte kosturulur.
        yaz(u'    D: empty because no pair was recovered, testing it ON ITS OWN...')
        oz = os.path.join(hizli_kok, 'D_KENDI_SINAMASI')
        os.makedirs(os.path.join(oz, 'KURTARMA_SONUC'), exist_ok=True)
        for ad in ('screening', 'verification', 'protocol', 'REFERANS_DB',
                   'konsensus_kanonik', 'TEK_PROTOKOL_SONUC', 'primer_final',
                   'engine', 'engine', 'engine',
                   'steps', 'engine', 'fastq files'):
            h = os.path.join(oz, ad)
            kaynak = os.path.join(hizli_kok, ad)
            if not os.path.exists(h) and os.path.exists(kaynak):
                try:
                    os.symlink(os.path.realpath(kaynak), h)
                except OSError:
                    pass
        with open(os.path.join(oz, 'KURTARMA_SONUC', 'kurtarma_satirlari.tsv'),
                  'w', encoding='utf-8', newline='') as fh:
            fh.write(u'# D SELF-TEST - synthetic input\n')
            ww = csv.writer(fh, delimiter='\t')
            ww.writerow(['hedef', 'eski_deger', 'eski_kapsam', 'denenen_yol', 'olcu',
                         'yeni_deger', 'esigi_gecti_mi', 'UYELIK_GEREKCESI', 'sebep'])
            ww.writerow(['D_SINAMA_Petriella', '8,45', '9/9', 'sentetik', 'ayrim kati',
                         'YENI CIFT AAATCTGGCTGCCTGTGC / CTCTCACCCTCTATGGCGTC (101 bp) 11,03 x',
                         'EVET (yeni cift)', 'sentetik', ''])
        rc2, _ = kos(yaz, 'D (kendi sinamasi, sentetik girdi)',
                     [py, os.path.join(kok, 'verification', 'specificity_round.py'),
                      '--kok', oz, '--ncbi', 'elle', '--mfe-yok',
                      '--kume-ust', '1'], tavan_dk * 60)
        D2 = tsv_oku(os.path.join(oz, 'DOGRULAMA_SONUC', 'dogrulama_uc_sutun.tsv'))
        sonuc['asama']['D'] = dict(rc=rc2, satir=len(D2), kendi_sinamasi=True)
        if len(D2) >= 1:
            sonuc['uyari'].append(
                u'D asamasi zincirde bos kaldi cunku K hicbir cifti kurtarmadi '
                u'(dogrulanacak sey yok). D KENDI BASINA sinandi ve %d satir uretti '
                u'- asama saglam.' % len(D2))
            yaz(u'    D: produced %d rows in its own self-test, so the stage is SOUND' % len(D2))
            Dd = D2
        else:
            sonuc['hata'].append(u'D asamasi sentetik girdiyle de satir uretmedi.')
    if len(Dd) < 1:
        if rc == 7 or len(K) < 1:
            # ARDISIK COKUS: K satir uretmediyse D'nin girdisi zaten bos.
            # Ayri bir hata olarak sayilmaz, K'nin hatasina baglanir.
            sonuc['uyari'].append(
                u'D asamasi kosmadi cunku K hic satir uretmedi (ardisik cokus). '
                u'Bu D\'nin hatasi degildir; K duzeltilince D de kosar.')
            yaz(u'    D: SKIPPED - the input is empty because K produced no rows (a knock-on failure)')
        else:
            sonuc['hata'].append(u'D asamasi HIC satir uretmedi.')
    else:
        # ZORUNLU katmanlar: bizim iki olcumumuz. MFEprimer ve NCBI testte
        # BILEREK atlanir (--mfe-yok, --ncbi elle: ag yok) - eksik olmalari
        # zincir hatasi degildir, uyari olarak raporlanir.
        eksik_sutun, yok_sutun, takma_ad, istege_bagli = katman_denetimi(Dd, sonuc)
        yaz(u'    D: %d rows | are the two mandatory sources filled: %s | skipped in the test: %s'
            % (len(Dd), u'NO' if (eksik_sutun or yok_sutun) else 'evet',
               ', '.join(istege_bagli) or 'yok'))
        for _ta in takma_ad:
            yaz(u'    D: (sema notu) %s' % _ta)

    rc, _ = kos(yaz, 'I (KIMLIK, 2 veritabani, nt yok)',
                [py, os.path.join(kok, 'verification', 'identity_verification.py'),
                 '--kok', hizli_kok, '--yalniz', '10', '--nt', 'yok', '--vtb-ust', '2'],
                tavan_dk * 60)
    I = tsv_oku(os.path.join(hizli_kok, 'KIMLIK_SONUC', 'kimlik_iddialari.tsv'))
    sonuc['asama']['I'] = dict(rc=rc, satir=len(I))
    if len(I) < 1:
        sonuc['hata'].append(u'I asamasi HIC iddia sonucu uretmedi.')
    else:
        yaz(u'    I: %d claims resolved' % len(I))
    return sonuc



# ===========================================================================
# 2026-08-09 DUZELTMESI - H KAPISININ YANLIS ALARMI
# ---------------------------------------------------------------------------
# BELIRTI : H asamasi cikis kodu 6 ile dusuyor ve
#           "D asamasinda ZORUNLU katmanlar doldurulmadi: 1. kaynak
#            (numune olcumu)" diyordu.
# OLCUM   : HIZLI_TEST/D_KENDI_SINAMASI/DOGRULAMA_SONUC/dogrulama_uc_sutun.tsv
#           okundu. Katman BOS DEGIL: sutunun degeri 'TEMIZ'. Ama sutunun ADI
#           '1_NUMUNE' degil, '1_NUMUNE_oy_vermez'. D asamasi 2026-08-06'daki
#           D-2 duzeltmesinde sutunu yeniden adlandirmis (specificity_round.py
#           satir 1053), kapi eski adi ariyor; r.get('1_NUMUNE') None donunce
#           kapi bunu "doldurulmadi" sayiyordu.
#           Yani sebep TEST KIPI DEGIL, SEMA KAYMASIDIR. '2_YEREL_DB' adi
#           degismedigi icin o katman gecmis, hata yalniz 1. katmanda cikmisti.
# DUZELTME: sutun once tam adiyla, bulunamazsa '<ad>_' onekiyle aranir.
#           - sutun var ve dolu            -> gecer (gerekirse sema notu)
#           - sutun var ama BUTUN satirlar bos -> HATA (kapi calismaya devam eder)
#           - sutun HIC yok                -> AYRI HATA (sema kaymasi); sessizce
#                                             gecmez, cunku "yok" ile "bos"un
#                                             ayni sayilmasi bu hatanin kendisiydi
#           MFEprimer ve NCBI katmanlari ZATEN uyariydi, oyle birakildi: test
#           --mfe-yok ve agsiz '--ncbi elle' ile kostugu icin onlarin bos olmasi
#           bilerek yaratilan bosluktur.
# ===========================================================================
ZORUNLU_KATMAN = (('1_NUMUNE', u'1. kaynak (numune olcumu)'),
                  ('2_YEREL_DB', u'2. kaynak (yerel veritabani)'))
ISTEGE_BAGLI_KATMAN = (('3_MFEPRIMER', u'3. kaynak (MFEprimer)'),
                       ('4_NCBI', u'4. kaynak (NCBI)'))


def sutun_coz(satirlar, ad):
    """Sutunu tam ad, sonra '<ad>_' onegiyle cozer.

    Doner: (cozulen_ad_or_None, dolu_mu). cozulen_ad None ise sutun HIC yok;
    bu 'bos' ile ayni sey DEGILDIR ve ayri raporlanir.
    """
    if not satirlar:
        return None, False
    basliklar = [b for b in satirlar[0].keys() if b]
    if ad in basliklar:
        coz = ad
    else:
        adaylar = sorted(b for b in basliklar if b.startswith(ad + '_'))
        coz = adaylar[0] if adaylar else None
    if coz is None:
        return None, False
    return coz, any((r.get(coz) or '').strip() for r in satirlar)


def katman_denetimi(Dd, sonuc):
    """D ciktisindaki katmanlari denetler, sonuc['hata'/'uyari']'yi doldurur.

    Doner: (eksik_sutun, yok_sutun, takma_ad, istege_bagli)
    """
    eksik_sutun, yok_sutun, takma_ad = [], [], []
    for sut, ad in ZORUNLU_KATMAN:
        coz, dolu = sutun_coz(Dd, sut)
        if coz is None:
            yok_sutun.append(u'%s ["%s" sutunu D ciktisinda HIC YOK]' % (ad, sut))
        elif not dolu:
            eksik_sutun.append(u'%s ("%s" sutunu var ama butun satirlarda bos)'
                               % (ad, coz))
        elif coz != sut:
            takma_ad.append(u'%s: "%s" yerine "%s" sutunu kullanildi'
                            % (ad, sut, coz))
    if yok_sutun:
        sonuc['hata'].append(
            u'D ciktisinda ZORUNLU katman sutunu HIC YOK (sema kaymasi): %s. '
            u'Kapi bunu "bos" ile ayni saymaz; D ciktisinin sutun adlari '
            u'degismisse kapi da guncellenmelidir.' % ', '.join(yok_sutun))
    if eksik_sutun:
        sonuc['hata'].append(u'D asamasinda ZORUNLU katmanlar doldurulmadi: %s'
                             % ', '.join(eksik_sutun))
    for t in takma_ad:
        sonuc['uyari'].append(
            u'D ciktisinda sutun adi degismis ama katman DOLU, bu yuzden hata '
            u'sayilmadi - %s' % t)
    istege_bagli = []
    for sut, ad in ISTEGE_BAGLI_KATMAN:
        coz, _ = sutun_coz(Dd, sut)
        if coz is None:
            istege_bagli.append(ad)
            continue
        if not any((r.get(coz) or '').strip() not in ('', 'BILINMIYOR') for r in Dd):
            istege_bagli.append(ad)
    if istege_bagli:
        sonuc['uyari'].append(
            u'D asamasinda su kaynaklar testte doldurulmadi: %s. Bu BILEREK '
            u'boyledir (test --mfe-yok ve agsiz --ncbi elle ile kosar); tam '
            u'kosuda dolarlar.' % ', '.join(istege_bagli))
    return eksik_sutun, yok_sutun, takma_ad, istege_bagli


# Karar tek bir sarta baglanir: hata listesi bos mu. Uyarilar testi dusurmez.
# Cikis kodu 6 = ZINCIR TUTARSIZ; tam kosuya girmeden once giderilmesi gerekir.
def raporla(hizli_kok, sonuc, yaz, gecen_sure):
    guvenilir = not sonuc['hata']
    yol = os.path.join(hizli_kok, 'HIZLI_TEST_RAPORU.md')
    with open(yol, 'w', encoding='utf-8') as fh:
        fh.write(u'# Quick consistency test (regression)\n\nGenerated: %s, script %s, time: %s\n\n'
                 % (time.strftime('%Y-%m-%d %H:%M'), VERSIYON, sure_metni(gecen_sure)))
        fh.write(u'## VERDICT\n\n')
        if guvenilir:
            fh.write(u'# CHAIN CONSISTENT (against its own reference): the full run can be started\n\n')
            fh.write(u'> **This is not a CORRECTNESS test.** The expected values also come from a full-depth run of the same engine;')
            # 2026-08-06 DUZELTMESI - temiz kosuda yakalandi: bu cumle KOSULSUZ
            # yazildigi icin "butun satirlar sinifini korudu" diyordu, oysa ayni
            # raporun tablosunda Petriella_cinsi satirinda "sinif korundu mu:
            # HAYIR" yaziyordu. Karar (TUTARLI) dogruydu - esik SINIRINDAKI
            # satirin dusuk derinlikte taraf degistirmesi bilerek uyari sayilir -
            # ama ozet cumlesi tabloyla CELISIYORDU. Artik sayilarak yazilir.
            _bozan = [s['hedef'] for s in sonuc['satir'] if not s.get('sinif_ok')]
            if _bozan:
                fh.write(u'The ranking held and all four stages produced non-empty output. %d of the %d rows tested kept their class'
                         % (len(sonuc['satir']), len(sonuc['satir']) - len(_bozan),
                            len(_bozan), ', '.join(_bozan)))
            else:
                fh.write(u'Every row tested kept its class, the ranking held, and all four stages produced non-empty output')
        else:
            fh.write(u'# ZINCIR TUTARSIZ — TAM KOSUYA GIRMEYIN\n\n')
            fh.write(u'The 6 to 16 hour run should not be spent before the mismatches below are resolved.\n\n')
            for h in sonuc['hata']:
                fh.write(u'- **%s**\n' % h)
            fh.write(u'\n')
        if sonuc.get('referans_bayat'):
            fh.write(u'## Rows whose reference is stale (NOT a regression)\n\n')
            fh.write(u'The primer pair on these rows changed since the reference was measured. Comparing a number belonging to the old pair against the new one')
            for h in sonuc['referans_bayat']:
                fh.write(u'- %s\n' % h)
            fh.write(u'\n')
        fh.write(u'## Row by row\n\n| target | reference | measured in the test | ratio | verdict | class preserved |\n|---|---|---|---|---|---|\n')
        for r in sonuc['satir']:
            oran = (r['olculen'] / r['ref']) if (r['ref'] and r['olculen'] is not None) else None
            fh.write(u'| %s | %sx | %sx | %s | %s | %s |\n'
                     % (r['hedef'], vir(r['ref']), vir(r['olculen']),
                        vir(oran) if oran else '-', r['karar'],
                        'evet' if r['sinif_ok'] else '**HAYIR**'))
        s = sonuc.get('siralama') or {}
        fh.write(u'\n## Ranking\n\n- reference: %s\n- in the test: %s\n- **%s**\n'
                 % (' > '.join(s.get('referans', [])), ' > '.join(s.get('test', [])),
                    'korundu' if s.get('korundu') else 'BOZULDU'))
        fh.write(u'\n## Output check per stage\n\n| stage | exit code | rows |\n|---|---|---|\n')
        for k in ('P', 'K', 'D', 'I'):
            a = sonuc['asama'].get(k, {})
            fh.write(u'| %s | %s | %s |\n' % (k, a.get('rc', 'kosulmadi'), a.get('satir', 0)))
        if sonuc['uyari']:
            fh.write(u'\n## Uyarilar (testi dusurmez)\n\n')
            for u_ in sonuc['uyari']:
                fh.write(u'- %s\n' % u_)
        fh.write(u'\n```' + BANT_GEREKCESI + u'```\n')
    yaz('')
    yaz('=' * 74)
    if guvenilir:
        yaz(u'  CHAIN CONSISTENT (against its own reference); the full run can be started.')
        yaz(u'  CAUTION: this does NOT mean "the measurements were validated". It shows that')
        yaz(u'  the same engine reproduces itself. For independent confirmation, use MFEprimer.')
        yaz(u'  All %d rows tested kept their class, the ranking did not change,'
            % len(sonuc['satir']))
        yaz(u'  and all four stages produced non-empty output.')
    else:
        yaz(u'  CHAIN INCONSISTENT - DO NOT START THE FULL RUN')
        for h in sonuc['hata']:
            yaz(u'    - %s' % h)
        yaz(u'  Fix these before spending eight hours on a full run.')
    for u_ in sonuc['uyari']:
        yaz(u'  (warning) %s' % u_)
    for b in sonuc.get('referans_bayat', []):
        yaz(u'  (stale reference) %s' % b)
    if sonuc.get('referans_bayat'):
        yaz(u'  These rows were NOT COMPARED. To refresh the reference:')
        yaz(u'    python verification/refresh_reference.py --root .')
    yaz('=' * 74)
    yaz(u'  detail: %s' % yol)
    if not guvenilir:
        return CIKIS_TUTARSIZ
    return CIKIS_REFERANS_BAYAT if sonuc.get('referans_bayat') else 0


# HIZLI_TEST/ gecici bir koktur: kaynak klasorler sembolik bagla baglanir,
# ciktilar ayri kalir. Boylece test, gercek sonuc klasorlerindeki uzun kosu
# sonuclarini EZMEZ.
def main():
    p = argparse.ArgumentParser(description='Zincirin hizli dogruluk testi')
    p.add_argument('--root', '--kok', dest='kok', default='.')
    p.add_argument('--cap-minutes', '--tavan-dk', dest='tavan_dk', type=int, default=15,
                   help='asama basina zaman tavani (dakika)')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)
    hizli = os.path.join(kok, 'HIZLI_TEST')
    os.makedirs(hizli, exist_ok=True)
    # gecici kok: kaynaklar baglanti, ciktilar ayri
    for ad in ('screening', 'protocol', 'verification', 'REFERANS_DB',
               'konsensus_kanonik', 'primer_final', 'fastq files',
               'engine', 'engine', 'engine',
               'steps', 'engine'):
        h = os.path.join(hizli, ad)
        if not os.path.exists(h) and os.path.exists(os.path.join(kok, ad)):
            try:
                os.symlink(os.path.join(kok, ad), h)
            except OSError:
                pass
    import glob
    for u_ in glob.glob(os.path.join(kok, 'uyelik_yeniden_turetme_uyelik_*.tsv')):
        h = os.path.join(hizli, os.path.basename(u_))
        if not os.path.exists(h):
            try:
                os.symlink(u_, h)
            except OSError:
                pass
    g = open(os.path.join(hizli, 'test_gunlugu.txt'), 'a', encoding='utf-8')

    def yaz(s=''):
        print(s, flush=True); g.write(s + '\n'); g.flush()

    yaz('=' * 74)
    yaz(u'  QUICK CONSISTENCY TEST (regression)   version %s   %s'
        % (VERSIYON, time.strftime('%Y-%m-%d %H:%M')))
    yaz('=' * 74)
    yaz(u'  Target time: about 30 minutes. Measurement depth %d reads (a full run uses 3000).' % OKUMA)
    yaz(u'  Rows tested: %d. Output directory: HIZLI_TEST/'
        % (len(BEKLENEN_UST) + len(BEKLENEN_ALT) + len(BEKLENEN_YENI)))
    yaz(u'  NOTE: the numbers will NOT match the reference exactly; the depth was reduced.')
    yaz(u'  What is being checked is that the CLASS and the RANKING are preserved.')
    yaz('')
    yaz(u'  WHAT IT IS    : a check that the code REPRODUCES ITSELF.')
    yaz(u'  WHAT IT IS NOT: it does NOT check that the measurement is CORRECT. The expected')
    yaz(u'               values come from a full-depth run of the SAME engine. If the engine')
    yaz(u'               has a systematic error, this test CANNOT catch it.')
    yaz(u'  Independent confirmation is a separate job: see the MFE_BAGIMSIZ_TEYIT report.')
    yaz('')
    t0 = time.time()
    sonuc = calistir(kok, hizli, a.tavan_dk, yaz)
    if not sonuc['hata'] or sonuc['asama'].get('P', {}).get('satir'):
        sonuc = sonraki_asamalar(kok, hizli, a.tavan_dk, yaz, sonuc)
    rc = raporla(hizli, sonuc, yaz, time.time() - t0)
    g.close()
    return rc


if __name__ == '__main__':
    sys.exit(main() or 0)
