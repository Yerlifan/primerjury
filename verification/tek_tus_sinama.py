# -*- coding: utf-8 -*-
u"""tek_tus.py'nin SENARYO SINAMALARI - hicbiri gercek klasore dokunmaz.

NEDEN GOLGE KOK
---------------
Sinama, asamalarin DUSMESINI ve YARIM KALMASINI da denemek zorunda. Bunu
gercek klasorde yapmak siparis karariniN dayandigi TSV'leri bozardi. Bu
yuzden her senaryo /tmp altinda kurulan bir GOLGE KOK'te kosar:
  * agir ve salt-okunur klasorler SEMBOLIK BAGLANTI ile baglanir (kopyalanmaz)
  * asama betikleri SAHTE betiklerle degistirilir (davranisi biz belirleriz)
  * butun yazmalar /tmp icinde kalir
Bagli klasore tek bayt yazilmaz. Sinama sonunda golge kok /tmp'de kalir.

KOSUM:  python3 verification/tek_tus_sinama.py --kok .
"""

import os, sys, io, json, time, shutil, argparse, subprocess, tempfile

BU = os.path.dirname(os.path.abspath(__file__))
TEK_TUS = os.path.join(BU, 'tek_tus.py')

# Sahte betik sablonu. Cikis kodunu ve urettigi dosyalari disaridan aliyoruz;
# boylece "asama dustu", "cikti bos", "cikti hic uretilmedi" gibi durumlarin
# hepsi gercek bir alt surecle sinanabiliyor.
SAHTE = u'''# -*- coding: utf-8 -*-
import os, sys, time
KOK = os.path.dirname(os.path.abspath(__file__))
while KOK and not os.path.exists(os.path.join(KOK, '_SINAMA_AYAR.json')):
    ust = os.path.dirname(KOK)
    if ust == KOK:
        break
    KOK = ust
import json
ayar = json.load(open(os.path.join(KOK, '_SINAMA_AYAR.json'), encoding='utf-8'))
ad = %(ad)r
# screening/__main__.py IKI asama tarafindan cagriliyor: 8 (--sina) ve
# S (--mod ozet). Hangisi oldugunu argv soyler; sahte betik de ayni ayrimi
# yapmali, yoksa S asamasi hicbir zaman ciktisini uretmez.
if ad == 'sina' and '--mod' in sys.argv:
    ad = 'S'
d = ayar.get(ad, {})
print('[sahte %%s] basladi' %% ad)
sys.stdout.flush()
if d.get('bekle'):
    for i in range(int(d['bekle'])):
        print('[sahte %%s] calisiyor %%d' %% (ad, i)); sys.stdout.flush()
        time.sleep(1)
for y, icerik in (d.get('yaz') or {}).items():
    t = os.path.join(KOK, y)
    os.makedirs(os.path.dirname(t), exist_ok=True)
    open(t, 'w', encoding='utf-8').write(icerik)
    print('[sahte %%s] yazdi: %%s' %% (ad, y))
if d.get('selftest_metni'):
    print(d['selftest_metni'])
print('[sahte %%s] bitti, cikis kodu %%s' %% (ad, d.get('rc', 0)))
sys.exit(int(d.get('rc', 0)))
'''

# (gercek yol, sahte betik adi) - tek_tus.py'nin ASAMALAR cizelgesindeki
# 'betik' alanlariyla BIREBIR ayni olmali.
BETIKLER = [
    ('screening/__main__.py', 'sina'),
    ('verification/hizli_tutarlilik_testi.py', 'H'),
    ('verification/erisim_dogrulama.py', 'E'),
    ('engine/uyelik_yeniden_turet.py', 'U'),
    ('protocol/tek_protokol_olc.py', 'P'),
    ('verification/kurtarma_turu.py', 'K'),
    ('verification/dogrulama_turu.py', 'D'),
    ('verification/kimlik_dogrulama.py', 'I'),
    ('verification/tum_kutu_kimlikleri.py', 'G'),
]

# On kontrolun ZORUNLU saydigi klasorler. Golge kokte hepsi bos da olsa
# VAR olmali; yoksa on kontrol dogru sekilde durur ve asil sinamaya
# giremeyiz. (On kontrolun DURDUGUNU ayrica S4 senaryosu sinar.)
KLASORLER = ['fastq files', 'consensus sequences', 'primer_final', 'REFERANS_DB',
             'screening', 'protocol', 'engine',
             'engine', 'engine', 'engine',
             'konsensus_kanonik', 'tools', 'verification']

MFE_IX = ['archaea.16S.fna', 'bacteria.16S.fna', 'fungi.ITS.fna',
          'fungi.28SrRNA.fna', 'fungi.18SrRNA.fna', 'SILVA_138.2_SSURef_NR99.fasta']
KUMELER = MFE_IX + ['SILVA_138.2_LSURef_NR99.fasta', 'UNITE_ITS.fasta',
                    'PR2_SSU_taxo_long.fasta', 'ROD_v1.2_operon_variants.fasta',
                    'ref_all2.fna']


def golge_kur(taban):
    u"""Bos ama on kontrolu gecebilen bir golge kok kurar."""
    if os.path.exists(taban):
        shutil.rmtree(taban)
    os.makedirs(taban)
    for k in KLASORLER:
        os.makedirs(os.path.join(taban, k), exist_ok=True)
        # klasorler BOS olmamali (on kontrol oge sayisi yaziyor, bos da gecer
        # ama gercege yakin olsun)
        with io.open(os.path.join(taban, k, '_yer_tutucu.txt'), 'w',
                     encoding='utf-8') as fh:
            fh.write(u'sinama\n')
    # screening bir PAKET olarak cagriliyor (python3 -m screening)
    with io.open(os.path.join(taban, 'screening', '__init__.py'), 'w',
                 encoding='utf-8') as fh:
        fh.write(u'')
    # mfeprimer ikilisi
    mp = os.path.join(taban, 'tools', 'mfeprimer')
    with io.open(mp, 'w', encoding='utf-8') as fh:
        fh.write(u'#!/bin/sh\nexit 0\n')
    os.chmod(mp, 0o755)
    # SILVA: DNA alfabesi kapisini gecen kucuk bir sahte fasta
    with io.open(os.path.join(taban, 'REFERANS_DB',
                              'SILVA_138.2_SSURef_NR99.fasta'), 'w',
                 encoding='utf-8') as fh:
        for i in range(50):
            fh.write(u'>sahte%d test\nACGTACGTACGTTTTTACGTACGTACGTTTTT\n' % i)
    for f in KUMELER:
        p = os.path.join(taban, 'REFERANS_DB', f)
        if not os.path.exists(p):
            with io.open(p, 'w', encoding='utf-8') as fh:
                fh.write(u'>sahte\nACGT\n')
    for f in MFE_IX:
        with io.open(os.path.join(taban, 'REFERANS_DB', f + '.primerqc.bin'),
                     'wb') as fh:
            fh.write(b'0' * 100)
    # Sahte asama betikleri
    for yol, ad in BETIKLER:
        t = os.path.join(taban, yol)
        os.makedirs(os.path.dirname(t), exist_ok=True)
        with io.open(t, 'w', encoding='utf-8') as fh:
            fh.write(SAHTE % dict(ad=ad))
    # Nihai hukum tablosu (ozetteki siparis tablosu bunu okur).
    # Gercekteki bicimin aynisi: yorum satirlari + YENI_HUKUM sutunu.
    # Uc satir siparise girer (SIPARIS EDILEBILIR x2 + KOSULLU x1), iki satir
    # girmez (ESIK ALTI + ONERILMEZ). Beklenen sayi: 3.
    with io.open(os.path.join(taban, 'ESIK_VE_OLCUT_2026-08-08.tsv'), 'w',
                 encoding='utf-8') as fh:
        fh.write(u'# sinama dosyasi\n')
        fh.write(u'hedef\tESKI_HUKUM\tYENI_HUKUM\tdCq_olculen\n')
        fh.write(u'A_hedefi\tKOSULLU\tSIPARIS EDILEBILIR (kosullu)\t7,28\n')
        fh.write(u'B_hedefi\tKOSULLU\tSIPARIS EDILEBILIR (kosullu)\t4,01\n')
        fh.write(u'C_evrensel\tRISKLI\tKOSULLU (kontrol primeri)\t-\n')
        fh.write(u'D_hedefi\tKOSULLU\tESIK ALTI (kosullu, silinmedi)\t-0,43\n')
        fh.write(u'E_hedefi\tRISKLI\tONERILMEZ (hedef numunede yok)\t-\n')
    return taban


def ayar_yaz(taban, ayar):
    with io.open(os.path.join(taban, '_SINAMA_AYAR.json'), 'w',
                 encoding='utf-8') as fh:
        fh.write(json.dumps(ayar, ensure_ascii=False, indent=1))


def basarili_ayar():
    u"""Butun asamalari basarili yapan varsayilan davranis."""
    return {
        'sina': dict(rc=0, selftest_metni=u'TUM SINAMALAR GECTI.'),
        'H': dict(rc=0, yaz={'HIZLI_TEST/HIZLI_TEST_RAPORU.md':
                             u'# rapor\n\nZINCIR TUTARLI (kendi referansina gore)\n'}),
        'E': dict(rc=0, yaz={'ERISIM_SONUC/erisim_dogrulama.tsv':
                             u'vtb\tsonuc\narchaea\tGECTI\n'}),
        'U': dict(rc=0, yaz={'uyelik_yeniden_turetme_uyelik_29991231.tsv':
                             u'kutu\thedef\nA1\tX\n'}),
        'P': dict(rc=0, yaz={'TEK_PROTOKOL_SONUC/panel_tek_protokol.tsv':
                             u'hedef\tayrim\nX\t9\n',
                             'TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv':
                             u'hedef\tSINIF\tF\tR\nX\tKESIN\tAAA\tTTT\n'}),
        'K': dict(rc=0, yaz={'KURTARMA_SONUC/kurtarma_satirlari.tsv':
                             u'hedef\tgecti\nX\tEVET\n'}),
        'D': dict(rc=0, yaz={'DOGRULAMA_SONUC/dogrulama_uc_sutun.tsv':
                             u'hedef\tKARAR\nX\tKOSULLU\n'}),
        'I': dict(rc=0, yaz={'KIMLIK_SONUC/kimlik_iddialari.tsv':
                             u'iddia\tsonuc\n1\tDOGRULANDI\n'}),
        'G': dict(rc=0, yaz={'TUM_KIMLIK_SONUC/tum_kutu_kimlikleri.tsv':
                             u'kutu\tkimlik\nA1\tX\n'}),
        'S': dict(rc=0, yaz={'KAPSAMLI_ARAMA_SONUC/00_OZET_HEPSI.md':
                             u'# ozet\n\nSahte ozet dosyasi, sinama icin.\n'}),
    }


def kos(taban, ek=(), zaman_asimi=180, sinyal_sn=None):
    argv = [sys.executable, TEK_TUS, '--kok', taban, '--onayla',
            '--ncbi', 'yok', '--canlilik', '3'] + list(ek)
    p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if sinyal_sn:
        # KESINTI SENARYOSU: belirlenen saniyede SIGINT gonderilir; bu,
        # Windows'ta Ctrl+C'nin surece ulasmasinin Linux karsiligidir.
        import signal as _s, threading as _t
        def kes():
            time.sleep(sinyal_sn)
            try:
                p.send_signal(_s.SIGINT)
            except Exception:
                pass
        th = _t.Thread(target=kes)
        th.daemon = True
        th.start()
    out = p.communicate(timeout=zaman_asimi)[0].decode('utf-8', 'replace')
    return p.returncode, out


SONUC = []


def sina(ad, kosul, ayrinti=u''):
    SONUC.append((ad, bool(kosul), ayrinti))
    print(u'  [%s] %s%s' % (u'GECTI' if kosul else u'DUSTU', ad,
                            (u'   -> ' + ayrinti) if ayrinti else u''))
    return bool(kosul)


def durum_oku(taban):
    y = os.path.join(taban, 'TEK_TUS_SONUC', 'durum.json')
    if not os.path.exists(y):
        return {}
    return json.load(io.open(y, encoding='utf-8'))


# ===========================================================================
def s1_sifirdan(ana):
    print(u'\n--- S1: SIFIRDAN KOSU (butun asamalar basarili) ---')
    t = golge_kur(os.path.join(ana, 's1'))
    ayar_yaz(t, basarili_ayar())
    rc, out = kos(t)
    d = durum_oku(t)
    sina(u'S1 cikis kodu 0', rc == 0, u'rc=%d' % rc)
    sina(u'S1 on kontrol gecti', u'ON KONTROL GECTI' in out)
    for k in ('8', 'H', 'E', 'U', 'P', 'K', 'D', 'I', 'G', 'S'):
        sina(u'S1 asama %s bitti' % k, d.get(k, {}).get('durum') == 'bitti',
             d.get(k, {}).get('durum', 'YOK'))
    sina(u'S1 Kraken asamalari arac yok diye atlandi',
         all(d.get(k, {}).get('durum', '').startswith('atlandi') for k in 'WXZ'),
         u'W=%s X=%s Z=%s' % tuple(d.get(k, {}).get('durum', '?') for k in 'WXZ'))
    sina(u'S1 sira dogru: U, P den once kostu',
         out.index(u'>> [') >= 0 and
         out.find(u'\n>> [') >= 0 and
         out.find(u'] U ') < out.find(u'] P '),
         u'U konumu %d, P konumu %d' % (out.find(u'] U '), out.find(u'] P ')))
    sina(u'S1 sira dogru: P, K, D bu sirada',
         out.find(u'] P ') < out.find(u'] K ') < out.find(u'] D '))
    sina(u'S1 sira dogru: I, G den once',
         out.find(u'] I ') < out.find(u'] G '))
    ozet = os.path.join(t, 'TEK_TUS_SONUC', '00_SABAH_OZETI.md')
    sina(u'S1 sabah ozeti uretildi', os.path.exists(ozet))
    m = io.open(ozet, encoding='utf-8').read() if os.path.exists(ozet) else u''
    sina(u'S1 ozette "DUSEN ASAMA YOK" yaziyor', u'DÜŞEN AŞAMA YOK' in m)
    sina(u'S1 ozette siparis tablosu var', u'KOSULLU' in m and u'A_hedefi' in m)
    sina(u'S1 siparis sayimi dogru (3 girer, 2 girmez)',
         u'SİPARİŞ EDİLEBİLİR: 3 çift = 6 oligo' in m and u'Sipariş dışı 2 çift' in m,
         u'ozette bulunan: %s' % (u' | '.join(
             x.strip() for x in m.splitlines() if u'SİPARİŞ EDİLEBİLİR' in x) or u'yok'))
    sina(u'S1 siparis tablosu ESIK_VE_OLCUT dosyasindan okundu',
         u'ESIK_VE_OLCUT_2026-08-08.tsv' in m and u'YENI_HUKUM' in m)
    sina(u'S1 gunluk dosyasi zaman damgali',
         any(x.startswith('gunluk_') for x in
             os.listdir(os.path.join(t, 'TEK_TUS_SONUC'))))
    return t


def s2_devam(ana, s1_taban):
    print(u'\n--- S2: IKINCI KOSU - bitmis asamalar atlanmali (0 sn) ---')
    t = s1_taban
    t0 = time.time()
    rc, out = kos(t)
    gecen = time.time() - t0
    d = durum_oku(t)
    sina(u'S2 cikis kodu 0', rc == 0, u'rc=%d' % rc)
    sina(u'S2 P atlandi (kontrol noktasi gecerli)',
         u'>> P' in out and u'ATLANDI' in out.split(u'>> P')[1][:400],
         out.split(u'>> P')[1][:120].replace(u'\n', u' ') if u'>> P' in out else u'')
    sina(u'S2 hepsi ikinci kosuda hizli bitti (< 40 sn)', gecen < 40,
         u'%.1f sn' % gecen)
    sina(u'S2 durum.json korundu', d.get('D', {}).get('durum') == 'bitti')
    return t


def s3_bayat(ana, taban):
    print(u'\n--- S3: BAYAT KONTROL NOKTASI - girdi ciktidan yeni olunca ---')
    # D'nin bir girdisine dokunuyoruz: hedef_klad.tsv. Bu, 2026-08-07'de
    # yasanan D-9 hatasinin ta kendisi: indeks yenilendi, kontrol noktasi
    # eski sifirlari geri okudu.
    hk = os.path.join(taban, 'screening', 'hedef_klad.tsv')
    with io.open(hk, 'w', encoding='utf-8') as fh:
        fh.write(u'hedef\tklad\nX\tBacteria\n')
    os.utime(hk, (time.time() + 5, time.time() + 5))   # ciktidan KESIN yeni
    rc, out = kos(taban, ek=['--plan'])
    sina(u'S3 D "BAYAT" diye isaretlendi',
         u'D  KOSACAK' in out and u'BAYAT' in out.split(u'D  KOSACAK')[1][:300],
         out.split(u'D  KOSACAK')[1][:150].replace(u'\n', u' ')
         if u'D  KOSACAK' in out else u'D KOSACAK satiri yok')
    sina(u'S3 bayat olmayan asamalar hala ATLANIR',
         u'I  ATLANIR' in out and u'G  ATLANIR' in out)
    rc, out = kos(taban)
    d = durum_oku(taban)
    sina(u'S3 yeniden kosuldu ve bitti', d.get('D', {}).get('durum') == 'bitti'
         and d.get('D', {}).get('sure', 0) > 0, u'sure=%s' % d.get('D', {}).get('sure'))
    return taban


def s4_eksik_dosya(ana):
    print(u'\n--- S4: GEREKLI DOSYA EKSIK - on kontrol DURDURMALI ---')
    t = golge_kur(os.path.join(ana, 's4'))
    ayar_yaz(t, basarili_ayar())
    os.remove(os.path.join(t, 'REFERANS_DB',
                           'SILVA_138.2_SSURef_NR99.fasta.primerqc.bin'))
    shutil.rmtree(os.path.join(t, 'fastq files'))
    rc, out = kos(t)
    sina(u'S4 cikis kodu 2 (on kontrol kapisi)', rc == 2, u'rc=%d' % rc)
    sina(u'S4 "ON KONTROL DUSTU" yazdi', u'ON KONTROL DUSTU' in out)
    sina(u'S4 eksik SILVA indeksini adiyla yazdi',
         u'SILVA_138.2_SSURef_NR99.fasta' in out and u'MFE indeksi' in out)
    sina(u'S4 eksik fastq klasorunu adiyla yazdi', u'fastq files' in out)
    sina(u'S4 kurulum komutunu yazdi', u'mfeprimer index -i' in out)
    sina(u'S4 HICBIR ASAMA KOSMADI',
         u'>> [1/' not in out and not os.path.exists(
             os.path.join(t, 'TEK_PROTOKOL_SONUC', 'panel_tek_protokol.tsv')))
    return t


def s5_asama_dustu(ana):
    print(u'\n--- S5: BIR ASAMA DUSTU (rc != 0) - bagimlilar KOSMAMALI ---')
    t = golge_kur(os.path.join(ana, 's5'))
    a = basarili_ayar()
    a['P'] = dict(rc=3, yaz={  # cikti URETILIYOR ama cikis kodu SIFIR DEGIL.
        'TEK_PROTOKOL_SONUC/panel_tek_protokol.tsv': u'hedef\tayrim\nX\t9\n',
        'TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv': u'hedef\tSINIF\tF\tR\nX\tKESIN\tA\tT\n'})
    ayar_yaz(t, a)
    rc, out = kos(t)
    d = durum_oku(t)
    sina(u'S5 zincir cikis kodu 3', rc == 3, u'rc=%d' % rc)
    sina(u'S5 P DUSTU (cikti dolu olsa da rc!=0 maskelenmedi)',
         d.get('P', {}).get('durum') == 'DUSTU', d.get('P', {}).get('durum', 'YOK'))
    sina(u'S5 dusme sebebi cikis kodunu yaziyor',
         'CIKIS KODU 3' in (d.get('P', {}).get('sebep') or ''),
         (d.get('P', {}).get('sebep') or '')[:90])
    sina(u'S5 K atlandi (bagimli)', d.get('K', {}).get('durum') == 'atlandi (bagimli)',
         d.get('K', {}).get('durum', 'YOK'))
    sina(u'S5 D atlandi (bagimli)', d.get('D', {}).get('durum') == 'atlandi (bagimli)',
         d.get('D', {}).get('durum', 'YOK'))
    sina(u'S5 BAGIMSIZ asamalar (I, G) yine de kostu',
         d.get('I', {}).get('durum') == 'bitti' and d.get('G', {}).get('durum') == 'bitti',
         u'I=%s G=%s' % (d.get('I', {}).get('durum'), d.get('G', {}).get('durum')))
    m = io.open(os.path.join(t, 'TEK_TUS_SONUC', '00_SABAH_OZETI.md'),
                encoding='utf-8').read()
    sina(u'S5 ozet "DUSEN ASAMALAR" bolumu acti', u'DÜŞEN AŞAMALAR' in m)
    sina(u'S5 ozet DUSTU kelimesini P satirinda tasiyor', u'DUSTU' in m)
    return t


def s6_bos_cikti(ana):
    print(u'\n--- S6: CIKIS KODU 0 AMA CIKTI BOS - "bitti" SAYILMAMALI ---')
    t = golge_kur(os.path.join(ana, 's6'))
    a = basarili_ayar()
    a['K'] = dict(rc=0, yaz={'KURTARMA_SONUC/kurtarma_satirlari.tsv':
                             u'hedef\tgecti\n'})     # yalniz baslik: 0 veri satiri
    ayar_yaz(t, a)
    rc, out = kos(t)
    d = durum_oku(t)
    sina(u'S6 K DUSTU (bos cikti geciyor sayilmadi)',
         d.get('K', {}).get('durum') == 'DUSTU', d.get('K', {}).get('durum', 'YOK'))
    sina(u'S6 sebep "BOS" diyor', u'BOS' in (d.get('K', {}).get('sebep') or u''),
         (d.get('K', {}).get('sebep') or u'')[:90])
    sina(u'S6 D atlandi (bagimli)', d.get('D', {}).get('durum') == 'atlandi (bagimli)')
    sina(u'S6 zincir cikis kodu 3', rc == 3, u'rc=%d' % rc)
    return t


def s7_kesinti(ana):
    print(u'\n--- S7: YARIDA KESILDI (Ctrl+C) ve AYNI TUSLA DEVAM ---')
    t = golge_kur(os.path.join(ana, 's7'))
    a = basarili_ayar()
    a['D'] = dict(rc=0, bekle=25, yaz={'DOGRULAMA_SONUC/dogrulama_uc_sutun.tsv':
                                       u'hedef\tKARAR\nX\tKOSULLU\n'})
    ayar_yaz(t, a)
    rc, out = kos(t, sinyal_sn=14)
    d = durum_oku(t)
    sina(u'S7 kesme cikis kodu 130', rc == 130, u'rc=%d' % rc)
    sina(u'S7 "KESME ISTEGI ALINDI" yazildi', u'KESME ISTEGI ALINDI' in out)
    sina(u'S7 D kesildi olarak damgalandi', d.get('D', {}).get('durum') == 'kesildi',
         d.get('D', {}).get('durum', 'YOK'))
    sina(u'S7 D den ONCEKI asamalar bitti damgali',
         all(d.get(k, {}).get('durum') == 'bitti' for k in ('H', 'E', 'U', 'P', 'K')),
         u', '.join(u'%s=%s' % (k, d.get(k, {}).get('durum')) for k in 'HEUPK'))
    sina(u'S7 kesintide de sabah ozeti yazildi',
         os.path.exists(os.path.join(t, 'TEK_TUS_SONUC', '00_SABAH_OZETI.md')))
    # DEVAM: ayni komut yeniden kosuluyor
    a['D'] = dict(rc=0, bekle=0, yaz={'DOGRULAMA_SONUC/dogrulama_uc_sutun.tsv':
                                      u'hedef\tKARAR\nX\tKOSULLU\n'})
    ayar_yaz(t, a)
    t0 = time.time()
    rc2, out2 = kos(t)
    gecen = time.time() - t0
    d2 = durum_oku(t)
    sina(u'S7-devam cikis kodu 0', rc2 == 0, u'rc=%d' % rc2)
    sina(u'S7-devam bitmis asamalar ATLANDI',
         u'ATLANDI' in out2 and u'>> [1/' in out2)
    sina(u'S7-devam D artik bitti', d2.get('D', {}).get('durum') == 'bitti',
         d2.get('D', {}).get('durum', 'YOK'))
    sina(u'S7-devam hizli bitti (< 40 sn)', gecen < 40, u'%.1f sn' % gecen)
    return t


def s8_belirlenimci_imza(ana, taban):
    print(u'\n--- S8: KONTROL NOKTASI ANAHTARI BELIRLENIMCI MI (md5) ---')
    # Ayni girdilerle iki ayri surecte hesaplanan imza AYNI olmali. Python'un
    # hash() fonksiyonu kullanilsaydi PYTHONHASHSEED yuzunden farkli cikardi.
    kod = (u'import sys,os,json;sys.path.insert(0,%r);'
           u'import tek_tus as T;'
           u'a=[x for x in T.ASAMALAR(dict(ncbi="yok",karac=None)) if x["kod"]=="D"][0];'
           u'print(T.girdi_imzasi(%r,a))' % (os.path.dirname(TEK_TUS), taban))
    im = []
    for tohum in ('0', '1', '12345'):
        cev = dict(os.environ, PYTHONHASHSEED=tohum)
        r = subprocess.run([sys.executable, '-c', kod], capture_output=True, env=cev)
        im.append(r.stdout.decode().strip() or r.stderr.decode()[-200:])
    sina(u'S8 imza uc farkli PYTHONHASHSEED ile AYNI',
         len(set(im)) == 1 and len(im[0]) == 16, u' / '.join(im))
    return taban


def s9_yalniz_atla(ana):
    print(u'\n--- S9: --yalniz ve --atla suzgecleri ---')
    t = golge_kur(os.path.join(ana, 's9'))
    ayar_yaz(t, basarili_ayar())
    rc, out = kos(t, ek=['--yalniz', '8S'])
    d = durum_oku(t)
    sina(u'S9 --yalniz 8S: yalniz 8 ve S kaydedildi',
         set(d.keys()) <= {'8', 'S'} and '8' in d, u'anahtarlar: %s' % sorted(d.keys()))
    t2 = golge_kur(os.path.join(ana, 's9b'))
    ayar_yaz(t2, basarili_ayar())
    rc2, out2 = kos(t2, ek=['--atla', 'IGD'])
    d2 = durum_oku(t2)
    sina(u'S9 --atla IGD: I, G, D hic kaydedilmedi',
         not ({'I', 'G', 'D'} & set(d2.keys())), u'anahtarlar: %s' % sorted(d2.keys()))
    sina(u'S9 --atla sonrasi cikis kodu 0', rc2 == 0, u'rc=%d' % rc2)
    return t


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--taban', default=os.path.join(tempfile.gettempdir(),
                                                   'tek_tus_sinama'))
    p.add_argument('--yalniz', default='')
    A = p.parse_args()
    ana = A.taban
    if not os.path.exists(TEK_TUS):
        print(u'HATA: tek_tus.py bulunamadi: %s' % TEK_TUS)
        return 1
    os.makedirs(ana, exist_ok=True)
    print(u'=' * 78)
    print(u'  tek_tus.py SENARYO SINAMALARI')
    print(u'  golge kok: %s   (bagli klasore TEK BAYT yazilmaz)' % ana)
    print(u'=' * 78)

    t0 = time.time()
    hepsi = A.yalniz.split(',') if A.yalniz else None

    def calis(ad):
        return (not hepsi) or ad in hepsi

    taban1 = None
    if calis('s1'):
        taban1 = s1_sifirdan(ana)
    if calis('s2') and taban1:
        s2_devam(ana, taban1)
    if calis('s3') and taban1:
        s3_bayat(ana, taban1)
    if calis('s8') and taban1:
        s8_belirlenimci_imza(ana, taban1)
    if calis('s4'):
        s4_eksik_dosya(ana)
    if calis('s5'):
        s5_asama_dustu(ana)
    if calis('s6'):
        s6_bos_cikti(ana)
    if calis('s7'):
        s7_kesinti(ana)
    if calis('s9'):
        s9_yalniz_atla(ana)

    gecti = [x for x in SONUC if x[1]]
    dustu = [x for x in SONUC if not x[1]]
    print(u'\n' + u'=' * 78)
    print(u'  SONUC: %d sinamanin %d tanesi GECTI, %d tanesi DUSTU  (%.0f sn)'
          % (len(SONUC), len(gecti), len(dustu), time.time() - t0))
    if dustu:
        print(u'  DUSENLER:')
        for ad, _, ayr in dustu:
            print(u'    * %s   %s' % (ad, ayr))
    print(u'=' * 78)
    return 0 if not dustu else 1


if __name__ == '__main__':
    sys.exit(main())
