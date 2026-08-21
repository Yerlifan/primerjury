# -*- coding: utf-8 -*-
"""
order_and_status.py - iki teslim dosyasi uretir, CANLI panelden:
  SIPARIS_LISTESI_20260802.tsv          yalniz siparis edilecek ciftler
  TOPLANTI_KARARLARI_DURUM_20260802.md  karar bazinda YAPILDI / KISMEN / YAPILAMIYOR

TAZELIK KURALI: xlsx her calistirmada bagli klasordeki CANLI dosyadan okunur,
/tmp kopyasi kullanilmaz. Yazdiktan sonra dosyalar geri okunup panelle karsilastirilir.
"""
# ---------------------------------------------------------------------------
# order_and_status.py — canli panelden iki teslim dosyasi uretir: siparis
#                       edilecek ciftlerin listesi ve toplanti kararlarinin
#                       yapildi / kismen / yapilamiyor durumu.
#
# GIRDI  : --xlsx ile bagli klasordeki CANLI panel dosyasi (openpyxl ile
#          okunur, md5'i alinip ekrana basilir) ve --kok proje klasoru. Panel
#          satiri -> toplanti karari eslemesi (KARAR_ESLEME) dosyanin icinde
#          sabittir.
# CIKTI  : SIPARIS_LISTESI_<tarih>.tsv ve TOPLANTI_KARARLARI_DURUM_<tarih>.md.
#          Yazdiktan sonra iki dosya da geri okunur ve primer dizileri panelle
#          karsilastirilir; bulunan her fark ekrana basilir.
# CAGRAN : MENUDE DEGILDIR - elle calistirilan bir teslim uretecidir. Menudeki
#          (P) tusu ayri bir siparis listesi uretir (protocol klasorunde) ve
#          bu betigi kullanmaz.
#
# TAZELIK KURALI KRITIKTIR: xlsx her calistirmada canli dosyadan okunur, gecici
# kopya kullanilmaz. Siparis dizisi bayat bir kopyadan alinirsa YANLIS PRIMER
# SIPARIS EDILIR ve hata ancak laboratuvarda gorunur.
# ---------------------------------------------------------------------------
import os, sys, csv, argparse, hashlib
import shutil
import openpyxl

# panel satiri -> (karar, istenen hedef)  ; Karar 1-4 toplantida ISTENEN hedeflerdir
KARAR_ESLEME = {
    22: ('Karar 1 - tur ozgul', 'Methanosarcina mazei'),
    21: ('Karar 1 - tur ozgul', 'Methanothrix soehngenii'),
    10: ('Karar 2 - cins ozgul', 'Proteiniphilum'),
    4:  ('Karar 2 - cins ozgul', 'Petrimonas'),
    5:  ('Karar 3 - islev/ekolojik grup', 'Hidrojenotrofik metanojenler'),
    16: ('Karar 3 - islev/ekolojik grup', 'Asetoklastik metanojenler'),
    8:  ('Karar 3 - islev/ekolojik grup', 'Metilotrofik metanojen'),
    9:  ('Karar 3 - islev/ekolojik grup', 'Nitrosocosmicus AOA'),
    15: ('Karar 3 - islev/ekolojik grup', 'Sakarolitik bakteriler'),
    11: ('Karar 3 - islev/ekolojik grup', 'Proteolitik / sintrofik bakteriler'),
    19: ('Karar 3 - islev/ekolojik grup', 'Proteolitik / sintrofik bakteriler'),
    20: ('Karar 4 - alan / evrensel', 'Arke universal'),
    17: ('Karar 4 - alan / evrensel', 'Bakteri universal'),
    6:  ('Karar 4 - alan / evrensel', 'Mantar universal (F1)'),
    13: ('Karar 4 - alan / evrensel', 'Mantar universal (F2)'),
    2:  ('Karar 4 - alan / evrensel', 'Universal metanojen'),
    12: ('Karar 5 - olcumden turetilen', 'Bacteroidales kumesi'),
    7:  ('Karar 5 - olcumden turetilen', 'Methanosarcina cinsi'),
    3:  ('Karar 5 - olcumden turetilen', 'Methanothrix cinsi'),
    14: ('Karar 5 - olcumden turetilen', 'Microascaceae askomikot'),
    18: ('Karar 5 - olcumden turetilen', 'Petriella musispora'),
}

# toplantida ISTENIP hic cift verilemeyen hedefler (panel satiri yok)
YAPILAMIYOR = [
 ('Karar 1 - tur ozgul', 'Methanosarcina barkeri',
  'Organizma numunede yok: 2208 kutusunun en yakin referansi M. vacuolata %97,4-97,9 - tur esiginin altinda. Yerine CINS duzeyi verildi (Methanosarcina_cinsi).'),
 ('Karar 1 - tur ozgul', 'Podospora pseudopauciseta',
  'Organizma numunede yok: bes referans ciftinden ucu F1 sinifinin 85 804 okumasinin tamaminda 0 urun verdi.'),
 ('Karar 1 - tur ozgul', 'Dictyostelium discoideum (44689)',
  'Kraken2 etiketi olcumle curutuldu: SILVA LSU 28S D1-D2 testinde D. discoideum skoru 195, rastgele bir mantar 480. Kutu heterojen (dort barkod birbirine %70-76).'),
 ('Karar 1 - tur ozgul', 'Trichoderma asperellum (101201)',
  'Kutudaki organizma Trichoderma degil: ITS Petriella/Microascaceae veriyor. Tasarlanan cift panelden cikarildi (ayrim 0,7x).'),
 ('Karar 2 - cins ozgul', 'Bacteroides',
  'Cins numunede yok: bes kutunun en iyi Bacteroidales eslesmesi %84,6-85,9, 16S cins esigi ~%94-95. Yerine adlandirilamayan Bacteroidales soyu icin kume cifti verildi.'),
 ('Karar 2 - cins ozgul', 'Alistipes',
  'Bacteroides ile AYNI organizma (kutular birbirine %95,3-96,8 benziyor); Bacteroidales_kumesi altinda birlestirildi.'),
 ('Karar 3 - islev/ekolojik grup', 'Trichoderma cinsi',
  'Hedef numunede var ama cins degil: olculen kimlik Petriella/Microascaceae. Cift panelden cikarildi (ayrim 0,7x).'),
]

# KISMEN gerekceleri (panel satiri VAR ama istenen duzey verilemedi)
KISMEN_NOT = {
 15: 'Toplanti "sakarolitik bakteriler" grubunu istedi; grup capinda cift bulunamadi. Verilen: UYE BAZLI tek cins cifti (Sphaerochaeta associata). Grubun diger uyeleri kapsanmiyor.',
 11: 'Toplanti "proteolitik/sintrofik bakteriler" grubunu istedi; grup capinda cift yok. Verilen: uye bazli iki ayri cift (Synergistaceae soyu + Cloacimonas cinsi). Hedef adi olculen kimlige cekildi - numunedeki organizma Cloacibacillus degil, adlandirilamayan Synergistaceae (%99,39; Cloacibacillus %90,02, cins esigi %94,5).',
 19: 'Ayni grubun ikinci uyesi. Bkz. Proteolitik_Synergistaceae satiri.',
 2:  'Kapsam TAM DEGIL: 34 metanojen kutusunun 33\'u cogaliyor. Ca. Methanomassiliicoccus kutusu bu ciftle cogalmiyor, onu Metilotrofik_metanojen cifti kapsiyor.',
 21: 'TUR duzeyi verildi ama KOSULLU: amplikon dizilemesi sarti var. Capraz vurusun 52\'si ayni ailedeki diger Methanothrix kayitlari - erime egrisi ayirmiyor.',
 22: 'TUR GRUBU verildi (M. mazei / M. soligelidi ayrilamiyor). OKUMA MOTORU DUZELTMESI SONRASI ESIK ALTI: en kotu tek kutu 0,82x. M. hadiensis kutusu hedef kadar iyi cogaliyor (%47,22). YENIDEN KARAR GEREKIYOR.',
 12: 'Toplantinin istedigi Bacteroides/Alistipes cinsleri numunede yok; yerine adlandirilamayan Bacteroidales SOYU icin kume cifti verildi. Ayrim 5,9x - 10x esiginin ALTINDA, teslim kosullu.',
 4:  'Cins OZGUL ama cins KAPSAMLI DEGIL: dogrulanmis P. sulfuriphila kutusu %57,0; Petrimonas kutularinin tamami kapsanmiyor.',
}


def kisa(v, n=110):
    s = '' if v is None else str(v).replace('\n', ' ').strip()
    return s if len(s) <= n else s[:n - 1] + '…'


# ---------------------------------------------------------------------------
# DUZELTME 2026-08-04: cikti adlari SABIT DEGILDIR.
#
# Eskiden bu betik ciktiyi her kosuda 'SIPARIS_LISTESI_20260802.tsv' adiyla
# yaziyordu. Sonucu suydu: 3 Agustos'ta uretilen daha yeni listeler yan yana
# durdugu halde, betigin urettigi dosya adiyla 2 Agustos'ta kalmis gorunuyor ve
# hangisinin guncel oldugu dosya adindan anlasilamiyordu.
#
# Yeni davranis uc parcalidir:
#   1) Tarihli dosya, PANELIN kendi degistirilme tarihinden turetilir. Boylece
#      ad, ciktinin gercekten hangi panel surumunden geldigini tasir.
#   2) Ayrica KANONIK bir 'SIPARIS_LISTESI.tsv' yazilir. Siparis her zaman bu
#      addan verilir; tarihli dosyalar kayit icin durur.
#   3) Klasorde bu betigin URETMEDIGI, daha yeni tarihli bir siparis listesi
#      varsa UYARI basilir ve kanonik dosya USTUNE YAZILMAZ. Sessizce eskiye
#      donmek, bu projede tam olarak bir kez pahaliya mal oldu.
# ---------------------------------------------------------------------------

KANONIK_SIPARIS = 'SIPARIS_LISTESI.tsv'


def panel_tarihi(xlsx_yolu):
    """Ciktinin tarih etiketini PANELDEN turetir, bugunun tarihinden degil.

    Neden: ayni panelden iki kez uretilen liste ayni adi almalidir. Bugunun
    tarihi kullanilsaydi, hicbir sey degismedigi halde her gun yeni bir dosya
    olusur ve 'hangisi guncel' sorusu yeniden dogardi.
    """
    import datetime
    return datetime.datetime.fromtimestamp(os.path.getmtime(xlsx_yolu)).strftime('%Y%m%d')


def daha_yeni_liste_var_mi(kok, bizim_dosya, kaynak_panel):
    """Bu betigin URETMEDIGI, daha yeni bir siparis listesi var mi.

    Ada gore degil, DEGISTIRILME TARIHINE gore bakar; ad tahmin ettirir, tarih
    olcer. Bulursa (dosya_adi, tarih) doner, yoksa None.

    DIKKAT - karsilastirma KAYNAK PANELIN tarihine gore yapilir, bu betigin az
    once yazdigi dosyanin tarihine gore DEGIL. Ikincisi her kosuda "simdi" olur
    ve hicbir aday ondan yeni cikamazdi; uyari hicbir zaman tetiklenmezdi.
    Sorulan soru sudur: turedigim panelden SONRA uretilmis bir liste var mi?
    """
    import glob, datetime
    bizim_t = os.path.getmtime(kaynak_panel)
    aday = []
    for d in ('', '1_TESLIM'):
        for y in glob.glob(os.path.join(kok, d, 'SIPARIS_LISTESI*.tsv')):
            if os.path.basename(y) == KANONIK_SIPARIS:
                continue
            if os.path.abspath(y) == os.path.abspath(bizim_dosya):
                continue
            if os.path.getmtime(y) > bizim_t:
                aday.append((y, os.path.getmtime(y)))
    if not aday:
        return None
    y, t = max(aday, key=lambda x: x[1])
    return (os.path.relpath(y, kok),
            datetime.datetime.fromtimestamp(t).strftime('%Y-%m-%d %H:%M'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--xlsx', required=True)
    ap.add_argument('--root', '--kok', dest='kok', required=True)
    a = ap.parse_args()

    # ---- CANLI dosyadan oku ------------------------------------------------
    xy = os.path.abspath(a.xlsx)
    ozet = hashlib.md5(open(xy, 'rb').read()).hexdigest()
    wb = openpyxl.load_workbook(xy, data_only=True)
    ws = wb['2 Panel']
    print('CANLI panel : %s' % xy)
    print('md5         : %s' % ozet)
    print('sayfa       : %d' % len(wb.sheetnames))

    satirlar = []
    for i in range(2, 23):
        satirlar.append(dict(
            satir=i, hedef=str(ws.cell(i, 3).value or ''), duzey=str(ws.cell(i, 4).value or ''),
            plaka=str(ws.cell(i, 1).value or ''), ta=ws.cell(i, 2).value,
            F=str(ws.cell(i, 6).value or ''), Fuz=ws.cell(i, 7).value, FTm=ws.cell(i, 8).value,
            R=str(ws.cell(i, 10).value or ''), Ruz=ws.cell(i, 11).value, RTm=ws.cell(i, 12).value,
            urun=ws.cell(i, 15).value, aralik=kisa(ws.cell(i, 27).value, 60),
            ayrim=kisa(ws.cell(i, 18).value, 70), siparis=str(ws.cell(i, 23).value or ''),
            durum=str(ws.cell(i, 20).value or ''),
            duz1=kisa(ws.cell(i, 29).value, 24), esik=kisa(ws.cell(i, 32).value, 24)))

    sip = [r for r in satirlar if not r['siparis'].upper().startswith('HAYIR')]
    hayir = [r for r in satirlar if r['siparis'].upper().startswith('HAYIR')]

    # ---- 1) SIPARIS LISTESI ------------------------------------------------
    ETIKET = panel_tarihi(xy)                      # panelden turetilir, sabit degil
    syol = os.path.join(a.kok, 'SIPARIS_LISTESI_%s.tsv' % ETIKET)
    BAS = ['#', 'Hedef', 'Duzey', 'Plaka', 'Ta (C)',
           "Ileri primer (5'->3')", 'Ileri uz', 'Ileri Tm',
           "Geri primer (5'->3')", 'Geri uz', 'Geri Tm',
           'Urun (bp)', 'Urun boyu araligi (numunede)', 'Siparis notu']
    with open(syol, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['# SIPARIS YALNIZ BU DOSYADAN VERILIR.'])
        w.writerow(['# Kaynak: %s ("2 Panel"), md5 %s' % (os.path.basename(xy), ozet)])
        w.writerow(['# Uretim: %s | Bu dosyada YALNIZ siparis edilecek %d cift vardir. '
                    'Siparis edilmeyecek satir yoktur.'
                    % (ETIKET[:4] + '-' + ETIKET[4:6] + '-' + ETIKET[6:], len(sip))])
        w.writerow(['# Panelden cikarilan ve siparis EDILMEYECEK ciftler: %s'
                    % ', '.join(r['hedef'] for r in hayir)])
        w.writerow(['# Toplam oligo: %d (her cift 2 oligo, dejenere baz YOK - toplanti karari).'
                    % (2 * len(sip))])
        w.writerow([])
        w.writerow(BAS)
        for n, r in enumerate(sip, 1):
            w.writerow([n, r['hedef'], r['duzey'], r['plaka'], r['ta'],
                        r['F'], r['Fuz'], r['FTm'], r['R'], r['Ruz'], r['RTm'],
                        r['urun'], r['aralik'], r['siparis']])

    # ---- GERI OKU VE PANELLE KARSILASTIR -----------------------------------
    okunan = list(csv.reader(open(syol, encoding='utf-8'), delimiter='\t'))
    veri = [x for x in okunan if x and x[0].isdigit()]
    kar = fark = 0
    for x in veri:
        h = x[1]
        p = next((r for r in sip if r['hedef'] == h), None)
        if p is None:
            fark += 1; print(u'  !! not in the panel:', h); continue
        for sut, deg in ((5, p['F']), (8, p['R'])):
            kar += 1
            if x[sut] != deg:
                fark += 1
                print(u'  !! SEQUENCE DIFFERENCE %s column %d: file=%s panel=%s' % (h, sut, x[sut], deg))
    print('\nSIPARIS DOSYASI: %s' % syol)
    print(u'  pairs to be ordered : %d   (NOT to be ordered: %d)' % (len(sip), len(hayir)))
    print('  karsilastirilan dizi  : %d' % kar)
    print('  bulunan fark          : %d' % fark)

    # ---- 2) TOPLANTI KARARLARI DURUMU --------------------------------------
    def kat(r):
        if r['siparis'].upper().startswith('HAYIR'):
            return 'YAPILAMIYOR'
        if r['satir'] in KISMEN_NOT:
            return 'KISMEN'
        return 'YAPILDI'

    yapildi, kismen = [], []
    for r in satirlar:
        if r['satir'] not in KARAR_ESLEME:
            continue
        k = kat(r)
        if k == 'YAPILDI':
            yapildi.append(r)
        elif k == 'KISMEN':
            kismen.append(r)

    myol = os.path.join(a.kok, 'TOPLANTI_KARARLARI_DURUM_%s.md' % ETIKET)
    L = []
    A = L.append
    A('# Toplantı kararları — ne yapabildik, ne yapamadık\n')
    A('2 Ağustos 2026. Sayılar **canlı panelden** alındı: '
      '`%s` (`2 Panel`), md5 `%s`.\n' % (os.path.basename(xy), ozet))
    A('## Hedef sayısı — belgeler arasındaki çelişki\n')
    A('Bazı belgeler "yirmi hedef", bazıları "altı tür + dört cins" diyor. **Esas alınan:** '
      'panelin karar defteri `6 Karar Durumu`. Ona göre toplantıda **21 hedef** istendi:\n')
    A('| Karar | İstenen hedef sayısı |')
    A('|---|---|')
    A('| Karar 1 — tür özgül | 6 |')
    A('| Karar 2 — cins özgül | 4 |')
    A('| Karar 3 — işlev/ekolojik grup | 7 |')
    A('| Karar 4 — alan/evrensel | 4 |')
    A('| **Toplam istenen** | **21** |')
    A('| Karar 5 — ölçümden türetilen (toplantıda istenmedi) | 8 |\n')
    A('"Altı tür + dört cins" ifadesi yalnız **Karar 1 + Karar 2**\'yi (10 hedef) anlatır, '
      'toplantının tamamını değil. "Yirmi" ise eski bir sayımdır; karar defterindeki 21 esastır.\n')
    A('---\n')

    # uye kumesi guvenilirligi - CANLI ciftler.tsv'den
    uyedur = {}
    cp = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ciftler.tsv')
    if os.path.exists(cp):
        for c in csv.DictReader(open(cp, encoding='utf-8'), delimiter='\t'):
            uyedur[int(c['satir'])] = c['uye_kumesi_durumu']

    def ayrim_ve_not(r):
        ham = r['duz1'] if r['duz1'] and r['duz1'] != '- / -' else ''
        notlar = []
        deger = ham or r['ayrim']
        try:
            enk = float(ham.split('/')[0].strip()) if ham else None
        except ValueError:
            enk = None
        if uyedur.get(r['satir']) == 'YENIDEN_KURULDU' and ham:
            # duzeltilmis olcum yeniden kurulan uye kumesiyle yapildi - tek basina guvenilmez
            deger = 'panel: %s' % r['ayrim']
            notlar.append('düzeltilmiş motorla üye kümesi **doğrulanamadı** (panel değeri gösteriliyor)')
        elif enk is not None and enk < 10:
            notlar.append('**en kötü tek kutu 10x ALTINDA**')
        if str(r['esik']).upper().startswith('EVET'):
            if not notlar:
                notlar.append('**10x eşiğinin altında**')
        return deger, '; '.join(notlar) or '-'

    A('## YAPILDI — çift var, sipariş edilebilir (%d)\n' % len(yapildi))
    A('| Karar | Hedef | İleri / Geri primer | Ürün | Ayrım | Uyarı |')
    A('|---|---|---|---|---|---|')
    for r in sorted(yapildi, key=lambda x: KARAR_ESLEME[x['satir']][0]):
        k, ist = KARAR_ESLEME[r['satir']]
        ay, nt = ayrim_ve_not(r)
        A('| %s | %s | `%s` / `%s` | %s bp | %s | %s |' % (k.split(' - ')[0], r['hedef'],
                                                           r['F'], r['R'], r['urun'], ay, nt))
    A('')
    A('> Ayrım sütunu: düzeltilmiş okuma motoruyla **en kötü tek kutu / havuz** (mm≤1). '
      '"panel:" ön eki, düzeltilmiş ölçümün üye kümesinin doğrulanamadığını ve panelin '
      'kendi değerinin gösterildiğini bildirir. "x/y kutu" yazan hedefler kapsam ölçüsüdür, '
      'ayrım ölçülmez. **10x eşiğinin altındaki çiftler sipariş edilebilir ama koşulludur** — '
      'amplikon dizilemesi ya da jel doğrulaması gerekir.\n')

    A('## KISMEN — bir şey verildi ama istenen düzeyde değil (%d)\n' % len(kismen))
    A('| Karar | İstenen | Verilen | Neden |')
    A('|---|---|---|---|')
    for r in sorted(kismen, key=lambda x: KARAR_ESLEME[x['satir']][0]):
        k, ist = KARAR_ESLEME[r['satir']]
        A('| %s | %s | %s (%s bp) | %s |' % (k.split(' - ')[0], ist, r['hedef'],
                                             r['urun'], KISMEN_NOT[r['satir']]))
    A('')

    A('## YAPILAMIYOR — hiç verilemedi (%d)\n' % len(YAPILAMIYOR))
    A('| Karar | İstenen | Ölçüm gerekçesi |')
    A('|---|---|---|')
    for k, ist, ger in YAPILAMIYOR:
        A('| %s | %s | %s |' % (k.split(' - ')[0], ist, ger))
    A('')
    A('---\n')
    A('## Sipariş\n')
    A('Sipariş edilecek **%d çift = %d oligo**. Diziler yalnız '
      '`%s` dosyasından kopyalanmalıdır.\n' % (len(sip), 2 * len(sip), KANONIK_SIPARIS))
    A('Panelden çıkarılan ve sipariş edilmeyecek çiftler: %s.\n'
      % ', '.join(r['hedef'] for r in hayir))
    open(myol, 'w', encoding='utf-8').write('\n'.join(L) + '\n')

    # ---- KANONIK KOPYA -----------------------------------------------------
    # Siparis her zaman KANONIK_SIPARIS adindan verilir. Ama once, bu betigin
    # uretmedigi daha yeni bir liste var mi diye bakilir; varsa kanonik dosyaya
    # DOKUNULMAZ ve kullanici uyarilir. Sessizce eskiye donmek yasak.
    kyol = os.path.join(a.kok, KANONIK_SIPARIS)
    yeni_olan = daha_yeni_liste_var_mi(a.kok, syol, xy)
    if yeni_olan:
        print(u'\n  WARNING: there is a NEWER order list that this script did not produce:')
        print('           %s   (%s)' % yeni_olan)
        print(u'         list produced : %s' % os.path.basename(syol))
        print(u'         %s WAS NOT OVERWRITTEN. Decide for yourself which one holds;'
              % KANONIK_SIPARIS)
        print(u'         if you want to change the canonical one, copy that file by hand.')
    else:
        shutil.copyfile(syol, kyol)
        print(u'\nCANONICAL LIST : %s   (the order is placed from THIS file)' % kyol)

    print('\nDURUM DOSYASI  : %s' % myol)
    print('  YAPILDI %d | KISMEN %d | YAPILAMIYOR %d' % (len(yapildi), len(kismen), len(YAPILAMIYOR)))
    return 0 if fark == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
