# -*- coding: utf-8 -*-
"""
kanonik_uret.py - TEK KANONIK KAYNAK uretir: konsensus_kanonik/

Butun konsensusleri (karisik yonlu klasorlerden) okur, yon.py ile SENSE yonune
cevirir ve TEK bir klasore yazar. Bundan sonra her betik BURAYI okur; hicbir
betigin kendi yon yamasi olmaz.

Girdi klasorleri (oncelik sirasiyla; ayni kutu birden fazla klasorde varsa
oncelikli olan kazanir ve digeri manifest'e "atlandi" yazilir):
    1. KAPSAMLI_ARAMA_SONUC/konsensus_yeni   (varsa - en yeni uretim)
    2. referans_konsensus/konsensus          (gece normalize edilmis set)
    3. consensus sequences                   (Ali'nin ozgun ciktisi)

Cikti:
    konsensus_kanonik/<kutu>_kanonik.fasta
    konsensus_kanonik/MANIFEST.tsv           her dosya: kaynak, eski yon, cevrildi mi
    konsensus_kanonik/BELIRSIZ.tsv           yonu belirlenemeyen dosyalar (MASKELI)

Kullanim:
    python kanonik_uret.py --kok ..
    python kanonik_uret.py --kok .. --yeniden        (varsa uzerine yaz)
"""
# ---------------------------------------------------------------------------
# kanonik_uret.py — butun konsensus klasorlerini tarar, her kutunun dizisini
#                   yon.py ile SENSE yonune cevirir ve tek bir kanonik klasore
#                   yazar; boylece yon sorusu tek yerde cozulmus olur.
#
# GIRDI  : --kok altindaki uc kaynak klasor, --oncelik ile secilen sirayla:
#          "consensus sequences" (panelin uzerine kuruldugu ozgun set),
#          "KAPSAMLI_ARAMA_SONUC/konsensus_yeni" (yeni uretim) ve
#          "referans_konsensus/konsensus". Yon karari yon.dosya_kanonik() ile.
# CIKTI  : konsensus_kanonik/ altina kutu basina <kutu>.kanonik.fa; ayrica
#          INDEKS.tsv (tuketicilerin okumasi gereken tek liste), MANIFEST.tsv
#          (kaynak, eski yon, cevrildi mi), BELIRSIZ.tsv ve OKUBENI.txt.
#          Cikis kodu 0 = kanonik klasorde cevrilmesi gereken dosya kalmadi.
# CAGRAN : hepsi.kanonik_kos() ayri bir surec olarak calistirir - tus 9'un
#          2. asamasi (--oncelik ozgun) ve 4. asamasi (--oncelik yeni). Yon
#          kapisi dusen her asamanin hata mesajinda da elle calistirilmasi
#          onerilen komut budur.
#
# INDEKS.tsv NEDEN VAR: bagli klasorde eski kalinti dosyalar SILINEMIYOR. glob
# ile okuyan bir tuketici o kalintilari da toplar ve karisik yonlu eski
# dosyalari kanonik sanir. Bu yuzden gecerli dosyalarin listesi ayri tutulur.
# ---------------------------------------------------------------------------
import os, sys, re, csv, glob, argparse, shutil

BURA = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BURA)
import yon

# Oncelik: ayni kutu birden fazla klasorde varsa ilk sirada olan kazanir.
#
# VARSAYILAN 'ozgun' (2026-08-02 duzeltmesi). Onceki varsayilan 'referans' idi ve
# SESSIZ BIR KAYNAK DEGISIKLIGI yapiyordu: referans_konsensus/ klasoru, panelin
# uzerine kuruldugu "consensus sequences/*_consensus_strict.fasta" dosyalarindan
# FARKLI bir konsensus yeniden-kurulumudur (uzunluklar bile farkli: 1503 vs 1534).
# Olculdu: 'referans' ile Bakteri_universal (EVRENSEL bakteri cifti) 20 B
# kutusunun yalniz 2'sinde urun veriyordu ve boy 130 yerine 135 cikiyordu;
# 'ozgun' ile 7 kutuda urun veriyor ve panelin 130 bp degeri yeniden uretiliyor.
# Panelin butun sayilari 'ozgun' set uzerinde olculmustur - taban o olmalidir.
# YON normalizasyonu ayri bir istir ve her iki kaynakta da uygulanir.
ONCELIK = {
 'referans': [('referans_konsensus', 'referans_konsensus/konsensus'),
              ('konsensus_yeni', 'KAPSAMLI_ARAMA_SONUC/konsensus_yeni'),
              ('ozgun', 'consensus sequences')],
 # 'yeni': gece uretimi bitince kullanilir. YEDEK SIRASI ONEMLI - konsensus_yeni
 # bir kutuyu uretememisse ONCE 'ozgun' (panelin tabani) gelir; referans_konsensus
 # en sona alindi cunku farkli bir yeniden-kurulumdur (bkz. yukaridaki not).
 'yeni':     [('konsensus_yeni', 'KAPSAMLI_ARAMA_SONUC/konsensus_yeni'),
              ('ozgun', 'consensus sequences'),
              ('referans_konsensus', 'referans_konsensus/konsensus')],
 'ozgun':    [('ozgun', 'consensus sequences'),
              ('konsensus_yeni', 'KAPSAMLI_ARAMA_SONUC/konsensus_yeni'),
              ('referans_konsensus', 'referans_konsensus/konsensus')],
}


_KUTU = re.compile(r'(?:^|[^A-Za-z0-9])(A1|A2|F1|F2|B)[-_](\d)(?!\d)')
_TAX = re.compile(r'(?<![0-9])(\d{3,7})(?![0-9])')


def kutu_adi(yol):
    """dosya adindan <sinif>-<barkod>_<taxid> kutusunu cikar.
    Kaynak klasorlerin adlandirmasi TUTARSIZ (A1-1-reads_2209_consensus_strict,
    A1_1_reads_1826872_consensus_strict, A1-1_2209_yeniden_konsensus ...) - bu
    yuzden ad parcalama degil, desen eslemesi kullanilir."""
    b = os.path.basename(yol).replace('.fasta', '')
    m = _KUTU.search(b)
    if not m:
        return None
    kutu = '%s-%s' % (m.group(1), m.group(2))
    kalan = b[m.end():]
    adaylar = [x for x in _TAX.findall(kalan)]
    if not adaylar:
        return None
    return '%s_%s' % (kutu, adaylar[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--kok', required=True)
    ap.add_argument('--cikti', default='konsensus_kanonik')
    ap.add_argument('--yeniden', action='store_true')
    ap.add_argument('--oncelik', default='ozgun', choices=sorted(ONCELIK))
    a = ap.parse_args()

    h = yon.kendini_sina()
    if h:
        sys.exit('yon.py kendini sinamayi GECEMEDI: %s' % h)
    print('yon.py kendini sinama: GECTI. Kanonik yon =', yon.KANONIK_YON)
    print('kaynak onceligi :', ' > '.join(e for e, _ in ONCELIK[a.oncelik]))

    cik = os.path.join(a.kok, a.cikti)
    os.makedirs(cik, exist_ok=True)
    # NOT: bagli klasorde dosya SILINEMIYOR (Operation not permitted). Bu yuzden
    # gecerli dosyalar '*.kanonik.fa' desenine yazilir ve INDEKS.tsv'ye kaydedilir.
    # Tuketiciler glob DEGIL INDEKS okumalidir; eski '*_kanonik.fasta' kalintilari inert.

    manifest, belirsiz, gorulen = [], [], {}
    for etiket, kl in ONCELIK[a.oncelik]:
        yollar = sorted(glob.glob(os.path.join(a.kok, kl, '**', '*.fasta'), recursive=True))
        for y in yollar:
            k = kutu_adi(y)
            if not k:
                continue
            sn = yon.sinifi(os.path.basename(y))
            if sn == '?':
                sn = yon.sinifi(y)
            if sn == '?':
                continue
            if k in gorulen:
                manifest.append(dict(kutu=k, sinif=sn, kaynak=etiket, dosya=os.path.relpath(y, a.kok),
                                     eski_yon='', cevrildi='', uzunluk='', durum='atlandi (%s kazandi)' % gorulen[k]))
                continue
            kayitlar, _ = yon.dosya_kanonik(y)
            kayitlar = [r for r in kayitlar if len(r[1]) >= 200]
            if not kayitlar:
                manifest.append(dict(kutu=k, sinif=sn, kaynak=etiket, dosya=os.path.relpath(y, a.kok),
                                     eski_yon='', cevrildi='', uzunluk=0, durum='bos/kisa - atlandi'))
                continue
            ad, dizi, karar, cev = max(kayitlar, key=lambda r: len(r[1]))
            if karar == 'BELIRSIZ':
                belirsiz.append(dict(kutu=k, sinif=sn, kaynak=etiket,
                                     dosya=os.path.relpath(y, a.kok), uzunluk=len(dizi),
                                     N_yuzde=round(100.0 * dizi.count('N') / len(dizi), 1),
                                     not_='yon belirlenemedi - KANONIGE ALINMADI, maskeli'))
                manifest.append(dict(kutu=k, sinif=sn, kaynak=etiket, dosya=os.path.relpath(y, a.kok),
                                     eski_yon='BELIRSIZ', cevrildi='', uzunluk=len(dizi),
                                     durum='BELIRSIZ - yazilmadi'))
                continue
            gorulen[k] = etiket
            cy = os.path.join(cik, '%s.kanonik.fa' % k)
            with open(cy, 'w', encoding='utf-8') as fh:
                fh.write('>%s kanonik=%s kaynak=%s eski_yon=%s cevrildi=%s\n'
                         % (k, yon.KANONIK_YON, etiket, karar, int(cev)))
                for i in range(0, len(dizi), 70):
                    fh.write(dizi[i:i + 70] + '\n')
            manifest.append(dict(kutu=k, sinif=sn, kaynak=etiket, dosya=os.path.relpath(y, a.kok),
                                 eski_yon=karar, cevrildi='EVET' if cev else 'hayir',
                                 uzunluk=len(dizi), durum='yazildi'))

    def yaz(ad, rows):
        if not rows:
            open(os.path.join(cik, ad), 'w', encoding='utf-8').write('(bos)\n')
            return
        with open(os.path.join(cik, ad), 'w', newline='', encoding='utf-8') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t')
            w.writeheader()
            for r in rows:
                w.writerow(r)
    yaz('MANIFEST.tsv', manifest)
    yaz('BELIRSIZ.tsv', belirsiz)
    yaz('INDEKS.tsv', [dict(kutu=m['kutu'], sinif=m['sinif'],
                            dosya='%s.kanonik.fa' % m['kutu'], kaynak=m['kaynak'],
                            eski_yon=m['eski_yon'], cevrildi=m['cevrildi'],
                            uzunluk=m['uzunluk'])
                       for m in manifest if m['durum'] == 'yazildi'])
    open(os.path.join(cik, 'OKUBENI.txt'), 'w', encoding='utf-8').write(
        'KANONIK KONSENSUS KLASORU\n'
        'Gecerli dosyalar: *.kanonik.fa  (liste: INDEKS.tsv)\n'
        'Bu klasorde *_kanonik.fasta uzantili KALINTI dosyalar olabilir - bunlar ilk\n'
        'kosunun hatali adlandirilmis ciktisidir, bagli klasorde silinemedi. YOKSAYIN.\n'
        'Her tuketici INDEKS.tsv okumalidir, glob KULLANMAMALIDIR.\n'
        'Kanonik yon: SENSE. Tanim ve olcut: KAPSAMLI_ARAMA/yon.py\n')

    yazilan = [m for m in manifest if m['durum'] == 'yazildi']
    cevrilen = [m for m in yazilan if m['cevrildi'] == 'EVET']
    print('\nkanonik klasor : %s' % cik)
    print('yazilan kutu   : %d' % len(yazilan))
    print('  cevrildi     : %d (ANTISENSE -> SENSE)' % len(cevrilen))
    print('  zaten sense  : %d' % (len(yazilan) - len(cevrilen)))
    print('BELIRSIZ       : %d (yazilmadi, BELIRSIZ.tsv)' % len(belirsiz))
    kay = {}
    for m in yazilan:
        kay[m['kaynak']] = kay.get(m['kaynak'], 0) + 1
    print('kaynak dagilimi:', ', '.join('%s=%d' % x for x in sorted(kay.items())))

    # DOGRULAMA: yazilan her dosya gercekten SENSE mi
    kotu = 0
    for y in sorted(glob.glob(os.path.join(cik, '*.kanonik.fa'))):
        kayitlar, sn = yon.dosya_kanonik(y)
        for ad, dizi, karar, cev in kayitlar:
            if cev:
                kotu += 1
    print('\nDOGRULAMA: kanonik klasorde hala cevrilmesi gereken dosya =', kotu,
          '(0 olmali)' if kotu == 0 else '*** SORUN ***')
    return 0 if kotu == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
