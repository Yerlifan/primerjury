# -*- coding: utf-8 -*-
"""ERISIM DOGRULAMASI - arama her veritabanini GERCEKTEN kullaniyor mu?

NEDEN VAR
---------
Ayarda 12 veritabani yazmasi, oralara duzgun BAKILDIGI anlamina gelmiyor.
Olculmus ornek: F2-1_101201 kutusu icin Petriella kayitlari fungi.28S
veritabaninda VARDI (P. guttulata %99,65, P. musispora %99,00) ama kisa listeye
GIRMIYORDU; arac en iyi isabet olarak %98,16'lik Parascedosporium'u bildirdi.
Kayit oradaydi, arama onu geri getiremiyordu.

YONTEM - geri cagirma (recall) olcumu
  1) Her veritabanindan, dosya boyunca YAYILMIS birkac kayit secilir.
  2) O kaydin KENDI dizisi sorgu olarak verilir -> kisa liste onu geri getiriyor mu?
  3) Ayni dizi %8 rastgele mutasyonla bozulup tekrar sorulur (nanopore konsensusu
     hatali olur; gercek kullanim bu). Yine geri getiriyor mu?
  (2) duserse arama o veritabanini FIILEN KULLANMIYOR demektir.
  (3) duserse arama yalnizca kusursuz dizilerde calisiyor demektir - bizim
     kullanimimizda ise konsensusler hatalidir.

Cikti: veritabani basina GECTI / DUSTU + sebep. Bu tablo olmadan hicbir kimlik
sonucuna guvenilmemelidir.
"""

# -------------------------------------------------------------------------
# access_check.py — her referans veritabaninin ARAMA TARAFINDAN gercekten
# kullanilip kullanilmadigini geri cagirma (recall) olcerek kanitlar.
#
# GİRDİ  : REFERANS_DB/ altindaki FASTA kumeleri (hangileri sorulacagi
#          identity_verification.py icindeki VTB listesinden okunur) ve o kumelerin
#          KENDI kayitlari - sorgu disaridan gelmez, veritabaninin icinden secilir.
# ÇIKTI  : ERISIM_SONUC/erisim_dogrulama.tsv  (veritabani basina GECTI / KISMEN /
#          DUSTU + sebep; dosyaya EKLENEREK yazilir, eski kosular korunur).
# ÇAĞRAN : verification/full_chain.py -> E tusu
#          (bat icinde: wsl -e python3 "verification/access_check.py" --kok .)
#
# NEDEN AYRI BIR OLCUM - KISA LISTE HIKAYESININ KOKU
# Bu betik, kisa liste kesme noktasinin BAGLAYICI oldugu doneme ait hatanin
# kanitini uretti: Petriella kayitlari fungi.28S dosyasinda VARDI (P. guttulata
# %99,65, P. musispora %99,00) ama kisa listeye HIC girmiyordu, cunku liste
# 60 kayitla kesiliyordu ve siralama olcutu (tohum sayisi) ile karar olcutu
# (hizalama kimligi) ayni sey degildi. Arac bu yuzden yanlis cinsi (%98,16
# Parascedosporium) en iyi isabet diye bildirdi. Duzeltme siralama olcutunu
# iyilestirmek DEGIL, listeyi 500'e buyutup kesme noktasini baglayici olmaktan
# cikarmaktir. Bu betik o duzeltmenin fiilen tuttugunu her kosuda yeniden olcer.
# -------------------------------------------------------------------------
import os, sys, csv, time, random, argparse


# ---------------------------------------------------------------------------
# GERI CAGIRMA OLCUMU: her veritabanindan dosya boyunca YAYILMIS kayitlar
# secilir, o kaydin KENDI dizisi sorgu olarak verilir ve kisa listenin onu geri
# getirip getirmedigine bakilir. Bir veritabani kendi kaydini geri getiremiyorsa
# o veritabani fiilen taranmiyor demektir.
#
# Ikinci sinama %8 mutasyonlu kopya iledir ve ASIL sinama odur: bizim
# konsensuslerimiz nanopore kaynaklidir, yani hatalidir. Yalniz kusursuz dizide
# calisan bir arama gercek kullanimda kayit KACIRIR - "KISMEN" hukmu tam olarak
# bu durumu isaretler.
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description='Veritabani erisim dogrulamasi')
    p.add_argument('--root', '--kok', dest='kok', default='.')
    p.add_argument('--db', '--vtb', dest='vtb', default=None, help='this database only (a fragment of the file name)')
    p.add_argument('--records', '--kayit', dest='kayit', type=int, default=3, help='veritabani basina test kaydi')
    p.add_argument('--cap', '--tavan', dest='tavan', type=int, default=0,
                   help='taranacak en fazla kayit (0 = all of them, default)')
    p.add_argument('--mutation', '--mutasyon', dest='mutasyon', type=float, default=0.08)
    a = p.parse_args()
    kok = os.path.abspath(a.kok)

    # ---- KOK DENETIMI (eklendi 2026-08-04) --------------------------------
    # Eskiden bu betik, --kok yanlis verildiginde ham bir Python izlemesi
    # (FileNotFoundError + traceback) basiyordu. Izleme teknik olarak dogruydu
    # ama kullaniciya ne yapacagini soylemiyordu. Artik once bakilir, sonra
    # yuklenir; eksik olan sey adiyla soylenir ve cozum yazilir.
    # DIKKAT: bu bir davranis degisikligidir - eskiden de sifirdan farkli kodla
    # duruyordu, degisen yalniz MESAJDIR. Olcum mantigina dokunulmadi.
    _kd_yolu = os.path.join(kok, 'verification', 'identity_verification.py')
    _eksik = []
    if not os.path.isdir(os.path.join(kok, 'verification')):
        _eksik.append('verification klasoru')
    elif not os.path.exists(_kd_yolu):
        _eksik.append('verification/identity_verification.py')
    if not os.path.isdir(os.path.join(kok, 'REFERANS_DB')):
        _eksik.append('REFERANS_DB klasoru')
    if _eksik:
        sys.stderr.write(
            u'ERROR: %s does not contain %s.\n  This script runs from the project root. The root is the same directory as\n  verification/full_chain.py, and it holds verification/ and REFERANS_DB/.\n  Correct use:  python3 verification/access_check.py --kok <project directory>\n  If you come from the menu the root is supplied correctly on its own (key E).\n'
            % (kok, u' and '.join(_eksik)))
        return 1
    if not [f for f in os.listdir(os.path.join(kok, 'REFERANS_DB'))
            if f.endswith(('.fasta', '.fna'))]:
        sys.stderr.write(
            u'ERROR: there is no FASTA file in the REFERANS_DB directory (%s).\n  The access check reads the databases\' OWN rec'
            % os.path.join(kok, 'REFERANS_DB'))
        return 1

    sys.path.insert(0, os.path.join(kok, 'verification'))
    import importlib.util as u
    sp = u.spec_from_file_location('kd', _kd_yolu)
    kd = u.module_from_spec(sp)
    try:
        sp.loader.exec_module(kd)
    except SystemExit:
        pass
    except Exception as e:
        sys.stderr.write(
            u'ERROR: identity_verification.py could not be loaded: %s: %s\n  This file may be corrupted. First check its syntax with  python3 -c "import ast; ast.parse(open(r\'%s\').read())"\n'
            % (type(e).__name__, e, _kd_yolu))
        return 1

    CIKTI = os.path.join(kok, 'ERISIM_SONUC')
    os.makedirs(CIKTI, exist_ok=True)
    vtb = [(e, d, t) for e, d, t, kullan, _n in kd.VTB
           if kullan and os.path.exists(os.path.join(kok, 'REFERANS_DB', d))]
    if a.vtb:
        vtb = [v for v in vtb if a.vtb.lower() in v[1].lower()]

    satirlar = []
    for etiket, dosya, lokus in vtb:
        yol = os.path.join(kok, 'REFERANS_DB', dosya)
        boyut = os.path.getsize(yol)
        print('\n=== %s (%s, %.0f MB) ===' % (etiket, dosya, boyut / 1e6), flush=True)
        # 1) yayilmis kayitlari topla
        # REZERVUAR ORNEKLEME: tek gecisde, dosyanin TAMAMINA yayilmis ornek.
        # Eski hali ilk N kaydi aliyordu; tavanin otesi HIC sinanmiyordu -
        # tam da avladigimiz kor nokta deseniydi (kullanici yakaladi).
        havuz, n, uygun = [], 0, 0
        K_ORNEK = max(a.kayit * 8, 24)
        rnd0 = random.Random(20260804)
        t0 = time.time()
        for bas, diz in kd.fasta_akisi(yol):
            n += 1
            if a.tavan and n > a.tavan:
                break
            if len(diz) < 800:
                continue
            uygun += 1
            if len(havuz) < K_ORNEK:
                havuz.append((n, bas, diz))
            else:
                j = rnd0.randrange(uygun)
                if j < K_ORNEK:
                    havuz[j] = (n, bas, diz)
        havuz.sort(key=lambda x: x[0])
        if not havuz:
            satirlar.append(dict(vtb=etiket, dosya=dosya, sonuc='OLCULEMEDI',
                                 sebep=u'800 bp ustu kayit yok (taranan %d)' % n))
            print(u'  NOT MEASURED: no suitable record'); continue
        # bas / orta / son bolgelerden birer tane garanti
        secilen = []
        if havuz:
            for oran in [i / float(max(a.kayit - 1, 1)) for i in range(a.kayit)]:
                idx = min(len(havuz) - 1, int(oran * (len(havuz) - 1)))
                if havuz[idx] not in secilen:
                    secilen.append(havuz[idx])
        kapsam = ('TAMAMI' if not a.tavan or n <= a.tavan else
                  'ILK %d KAYIT (tavan)' % a.tavan)
        print(u'  %d records scanned (%.0f s) - coverage: %s | test record: %s'
              % (n, time.time() - t0, kapsam,
                 ', '.join('#%d' % x[0] for x in secilen)), flush=True)

        tam_ok = mut_ok = 0
        ayrinti = []
        for idx, bas, diz in secilen:
            q = diz[:4000]
            kl = kd.kisa_liste(yol, q, ilerle=None)
            # DUZELTME 2026-08-09: kisa_liste() 3'lu demet degil, DICT listesi
            # donduruyor (anahtarlar: tohum, skor, baslik, dizi, sira, kaynak).
            # Eski satir "too many values to unpack (expected 3, got 6)" veriyordu
            # ve E asamasi bu yuzden cikis kodu 1 ile dusuyordu.
            bulundu = any(r['baslik'] == bas for r in kl)
            tam_ok += 1 if bulundu else 0
            rnd = random.Random(idx)
            lq = list(q)
            for _ in range(int(len(lq) * a.mutasyon)):
                i = rnd.randrange(len(lq)); lq[i] = rnd.choice('ACGT')
            kl2 = kd.kisa_liste(yol, ''.join(lq), ilerle=None)
            bulundu2 = any(r['baslik'] == bas for r in kl2)
            mut_ok += 1 if bulundu2 else 0
            ayrinti.append('#%d %s: tam=%s mutasyonlu=%s'
                           % (idx, bas[:40], 'EVET' if bulundu else 'HAYIR',
                              'EVET' if bulundu2 else 'HAYIR'))
            print('    %s' % ayrinti[-1], flush=True)
        k = len(secilen)
        if tam_ok == k and mut_ok == k:
            sonuc, sebep = 'GECTI', u'tam ve %%%d mutasyonlu sorguda %d/%d geri getirildi' % (
                int(a.mutasyon * 100), k, k)
        elif tam_ok == k:
            sonuc, sebep = 'KISMEN', (u'tam dizide %d/%d geri getirdi ama %%%d mutasyonlu '
                                      u'sorguda yalniz %d/%d - hatali konsensusle KACIRIYOR'
                                      % (k, k, int(a.mutasyon * 100), mut_ok, k))
        else:
            sonuc, sebep = 'DUSTU', (u'kendi kaydini bile geri getiremedi (%d/%d tam, %d/%d '
                                     u'mutasyonlu) - arama bu veritabanini FIILEN KULLANMIYOR'
                                     % (tam_ok, k, mut_ok, k))
        print('  -> %s: %s' % (sonuc, sebep), flush=True)
        satirlar.append(dict(vtb=etiket, dosya=dosya, mb=round(boyut / 1e6),
                             taranan=n, kapsam=kapsam, sinama=k, tam=tam_ok, mutasyonlu=mut_ok,
                             sonuc=sonuc, sebep=sebep, ayrinti='; '.join(ayrinti)))

    yol = os.path.join(CIKTI, 'erisim_dogrulama.tsv')
    yeni = not os.path.exists(yol)
    with open(yol, 'a', encoding='utf-8', newline='') as fh:
        if yeni:
            fh.write(u'# ACCESS VERIFICATION - does the search really use every database?\n')
            fh.write(u'# PASSED: it retrieves its own record from both the exact and the flawed sequence.\n')
            fh.write(u'# PARTIAL: it only retrieves on a flawless sequence, and our consensuses are imperfect!\n')
            fh.write(u'# FAILED: it cannot even retrieve its own record, so the database is effectively unused.\n')
        w = csv.DictWriter(fh, delimiter='\t',
                           fieldnames=['vtb', 'dosya', 'mb', 'taranan', 'kapsam', 'sinama', 'tam',
                                       'mutasyonlu', 'sonuc', 'sebep', 'ayrinti'])
        if yeni:
            w.writeheader()
        for s in satirlar:
            w.writerow({k: s.get(k, '') for k in w.fieldnames})
    print('\nyazildi: %s' % yol)
    return 0


if __name__ == '__main__':
    sys.exit(main() or 0)
