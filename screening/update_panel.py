# -*- coding: utf-8 -*-
r"""
update_panel.py - okuma motoru duzeltmesini panele isler.

  * yeni sayfa: "16 Okuma Motoru Duzeltmesi"
  * "2 Panel"  : her satira olcut etiketi + duzeltilmis degerler sutunlari + NOT
  * "1 Rapora Ozet" ve "6 Karar Durumu": degisen satirlar ve uyari
  * degerlendiriciya_ozet TSV'si senkron

Kullanim: python update_panel.py --xlsx ..\PrimerJury_PCR_Paneli_2026-08-02_TESLIM.xlsx
                                   --tsv  ..\okuma_motoru_duzeltmesi_20260802.tsv
"""
# ---------------------------------------------------------------------------
# update_panel.py — okuma motoru duzeltmesinin sonuclarini teslim panelinin
#                     xlsx dosyasina isler (yeni sayfa acar, mevcut sayfalara
#                     olcut etiketi ve duzeltilmis deger sutunlari ekler).
#
# GIRDI  : --xlsx ile verilen teslim panel dosyasi (openpyxl ile acilir) ve
#          --tsv ile verilen okuma motoru duzeltmesi tablosu; istege bagli
#          --yedek yolu.
# CIKTI  : ayni xlsx dosyasinin uzerine yazar (wb.save) ve istenirse yedegini
#          alir; ayrica degerlendiriciya_ozet TSV'sini senkronlar. Ekrana guncellenen
#          dosya adini ve sayfa sayisini basar.
# CAGRAN : MENUDE DEGILDIR - bilerek elle calistirilir. Sebep 00_OZET_HEPSI.md
#          icinde yazili: teslim dosyasini DEGISTIRIR ve o dosyaya baska
#          oturumlar da yaziyor, bu yuzden otomatik zincire baglanmamistir.
# ---------------------------------------------------------------------------
import sys, os, csv, argparse, shutil, datetime
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment

KIRMIZI = PatternFill('solid', fgColor='FFC7CE')
SARI    = PatternFill('solid', fgColor='FFEB9C')
YESIL   = PatternFill('solid', fgColor='C6EFCE')
GRI     = PatternFill('solid', fgColor='D9D9D9')
KALIN   = Font(bold=True)
SAR     = Alignment(wrap_text=True, vertical='top')


def yaz(ws, r, c, v, fill=None, bold=False, wrap=True):
    h = ws.cell(r, c, v)
    if fill: h.fill = fill
    if bold: h.font = KALIN
    if wrap: h.alignment = SAR
    return h


def sayfa16(wb, rows):
    ad = '16 Okuma Motoru Duzeltmesi'
    if ad in wb.sheetnames:
        del wb[ad]
    ws = wb.create_sheet(ad)
    for w, c in zip((6, 30, 15, 15, 15, 15, 15, 15, 46),
                    'ABCDEFGHI'):
        ws.column_dimensions[c].width = w
    n = 1
    yaz(ws, n, 1, 'OKUMA MOTORU DUZELTMESI - 2026-08-02', bold=True); n += 1
    yaz(ws, n, 1, u'A SILENT bug was found and fixed in the panel\'s sample measurement engine. This page records the bug, the fix, which values changed and which pairs fell below the 10x threshold.'); n += 2

    yaz(ws, n, 1, '1. HATA NEYDI', bold=True, fill=GRI); n += 1
    for s in [
        'engine/reads.py -> sinif `Sonda` (ve engine/scb.py -> sinif `S`) '
        'primeri okumada ararken 13 BAZLIK TAM ESLESEN TEK BIR TOHUM kullaniyordu:  '
        's = primer[-13:] ;  i = seq.find(s)',
        'Olcut "toplam uyumsuzluk <= max_mm" oldugu halde, uyumsuzluk 13 bazlik tohumun ICINE '
        'dustugunde find() hicbir sey bulamiyor ve baglanma yeri SESSIZCE kayboluyor. '
        'Program hata atmiyor, "urun yok" diyor. Bes onceki olcum hatasi gibi bu da sessiz.',
        'EK BULGU: 3\' son 2 baz TAM ESLESME kurali kodda HICBIR YERDE acikca uygulanmiyordu - '
        '13 bazlik tohumun yan etkisiydi. Duzeltilmis motorda kural acikca uygulanir.',
    ]:
        yaz(ws, n, 1, s); ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=9)
        ws.row_dimensions[n].height = 42; n += 1
    n += 1

    yaz(ws, n, 1, u'2. EVIDENCE - A TEST WITH A KNOWN ANSWER (bin A1-4_3078083, first 400 reads, the M. mazei pair, the SAME criterion mm<=1)', bold=True, fill=GRI); n += 1
    bas = ['Motor', 'Urun veren okuma / 400', 'Yuzde', 'Not']
    for j, h in enumerate(bas, 1):
        yaz(ws, n, j, h, bold=True, fill=GRI)
    n += 1
    for r in [('okuma.Sonda (panelin kullandigi)', '2', '0,50%', 'tohumlu'),
              ('kaba kuvvet (saf python, tohumsuz, bagimsiz yazildi)', '174', '43,50%', 'dogru cevap'),
              ('ispcr.find_sites (numpy, tohumsuz, panelin kendi kodu)', '174', '43,50%', 'kaba kuvvetle BIREBIR'),
              ('read_engine.py (duzeltilmis)', '174', '43,50%', 'kaba kuvvetle BIREBIR')]:
        for j, v in enumerate(r, 1):
            yaz(ws, n, j, v, fill=(KIRMIZI if r[1] == '2' else YESIL))
        n += 1
    yaz(ws, n, 1, u'In 202 of the forward primer\'s 205 binding sites (98.5%) the single mismatch falls inside the seed, and 188 of those are at base 6 of the primer. That is not scattered noise, it is one recurring variant base. MISS RATE: 172/174 = 98.9%.')
    ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=9); ws.row_dimensions[n].height = 30; n += 1
    yaz(ws, n, 1, u'CORRECTION: the "188/400" figure in the earlier note file came from running ispcr.amplify with its DEFAULT max_mm=3, so that comparison changed the seed AND the criterion at the same time. The 2 vs 174 above changes only the seed, and that is the number that gives the size of the bug.', fill=SARI)
    ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=9); ws.row_dimensions[n].height = 30; n += 2

    yaz(ws, n, 1, '3. NASIL DUZELTILDI', bold=True, fill=GRI); n += 1
    for s in [
        'GUVERCIN YUVASI (pigeonhole) TOHUMLAMASI: primer max_mm+1 tane ORTUSMEYEN bloga bolunur. '
        'En fazla max_mm uyumsuzluk varsa bu bloklardan EN AZ BIRI tam tutmak zorundadir; dolayisiyla '
        'bloklardan herhangi birinin tam eslesmesini aramak KAYIPSIZDIR. Her aday yer sonra tam kuralla '
        '(toplam uyumsuzluk + 3\' son 2 baz) tek tek dogrulanir.',
        'KAYIPSIZLIK KANITLANDI: screening/engine_test.py duzeltilmis motoru, tohum kullanmayan ve her '
        'pozisyonu tek tek deneyen bagimsiz bir kaba kuvvet uygulamasiyla karsilastirir. '
        'T1 sentetik (2 439 baglanma yeri, fark 0) - T2 gercek okuma (fark 0) - T3 urun sayisi (fark 0). '
        'Ayni testte ESKI motor sentetikte 1 386 yeri (%56,8), gercek okumada %35,2\'sini kaciriyor.',
        'HIZ: 400 okuma/cift icin kaba kuvvet 0,30 sn, duzeltilmis motor 0,03 sn (~10x). Kaba kuvvetle ayni '
        'sonucu verdigi icin hiz kazanci dogruluktan odun vermez.',
        'DOSYALAR: screening/read_engine.py (motor), brute_force.py (referans), engine_test.py (kanit), '
        'measure_panel.py (panel olcumu), independent_check.py, '
        'python3 -m screening --mod panel-olc --full-depth',
    ]:
        yaz(ws, n, 1, s); ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=9)
        ws.row_dimensions[n].height = 46; n += 1
    n += 1

    yaz(ws, n, 1, '4. IKINCI BULGU - OLCUT TUTARSIZLIGI (en az tohum hatasi kadar onemli)',
        bold=True, fill=GRI); n += 1
    for s in [
        '"2 Panel" sayfasindaki olcut notu numune olcutunu "uyumsuzluk <=1, 3\' son 2 baz TAM eslesme" '
        'diye tanimliyor. Ama satirlar TEK BIR OLCUTLE OLCULMEMIS. Her satirin yayimlanmis uye % araligi, '
        'dort motor/olcut kombinasyonuyla yeniden uretilerek hangi olcutle olculdugu tespit edildi:',
        'mm<=1 (tohumlu motor): Metanomikrobiyales (uye %51,2-80,0 BIREBIR), Proteolitik_Synergistaceae '
        '(%81,5 -> %81,2), Methanosarcina_mazei_turu (%40,6-60,7 ve rakip %4,49 BIREBIR).',
        'mm<=3 (tohumsuz motor = ispcr.amplify varsayilani): Proteiniphilum (panelde %29,0 -> mm<=3 ile '
        '%28,6, mm<=1 ile %1,6), Metilotrofik (%71,0 -> %72,4), Cloacimonas (%78,5 -> %77,8), '
        'Sakarolitik_Sphaerochaeta, Methanosarcina_cinsi.',
        'SONUC: "Ayrim (x)" sutunundaki degerler birbirleriyle KARSILASTIRILAMAZ. Bu tabloda her satir '
        'TEK motorla ve IKI olcutle (mm<=1 ve mm<=3) yeniden olculdu; her deger olcut etiketi tasiyor. '
        'OLCUT ETIKETI OLMAYAN HICBIR SAYIYI KULLANMAYIN.',
    ]:
        yaz(ws, n, 1, s, fill=(SARI if s.startswith('SONUC') else None))
        ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=9)
        ws.row_dimensions[n].height = 42; n += 1
    n += 1

    yaz(ws, n, 1, u'5. WHICH VALUES CHANGED (21 pairs, corrected engine, at most 3000 reads per bin)',
        bold=True, fill=GRI); n += 1
    yaz(ws, n, 1, u'THIS IS A SUBSET MEASUREMENT: at most 3000 reads per bin (the panel\'s own protocol). Full depth verification is done with python3 -m screening --mod panel-olc --full-depth; if the trend changes, this page must be updated.', fill=SARI)
    ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=9); ws.row_dimensions[n].height = 30; n += 1
    bas = ['#', 'Hedef', 'Panelin olcutu', 'ESKI kat (en kotu/havuz)', 'YENI mm<=1 (en kotu/havuz)',
           'YENI mm<=3 (en kotu/havuz)', 'Uye kumesi', '10x alti?', 'Not']
    for j, h in enumerate(bas, 1):
        yaz(ws, n, j, h, bold=True, fill=GRI)
    n += 1
    for r in rows:
        deg = bool(r['DEGISTI'])
        f = KIRMIZI if (r['ESIK_ALTI_mm1'] == 'EVET' and deg) else (SARI if deg else None)
        yaz(ws, n, 1, int(r['satir']), fill=f)
        yaz(ws, n, 2, r['hedef'], fill=f)
        yaz(ws, n, 3, r['panel_olcut_tespiti'], fill=f)
        yaz(ws, n, 4, '%s / %s' % (r['ESKI_kat_enkotu'] or '-', r['ESKI_kat_havuz'] or '-'), fill=f)
        yaz(ws, n, 5, '%s / %s' % (r['YENI1_kat_enkotu'] or '-', r['YENI1_kat_havuz'] or '-'), fill=f)
        yaz(ws, n, 6, '%s / %s' % (r['YENI3_kat_enkotu'] or '-', r['YENI3_kat_havuz'] or '-'), fill=f)
        yaz(ws, n, 7, r['uye_kumesi'], fill=f)
        yaz(ws, n, 8, ('mm<=1 ' if r['ESIK_ALTI_mm1'] else '') + ('mm<=3' if r['ESIK_ALTI_mm3'] else ''), fill=f)
        yaz(ws, n, 9, r['DEGISTI'] or u'no change (below the 5 point / 2 fold threshold)', fill=f)
        ws.row_dimensions[n].height = 30
        n += 1
    n += 1

    yaz(ws, n, 1, u'6. THE ONE REAL REGRESSION - Methanosarcina_mazei_turu (row 22)', bold=True, fill=KIRMIZI); n += 1
    for s in [
        'Duzeltmenin bir paneli degeri ESIGIN ALTINA dusurdugu TEK satir budur. Digerlerinde duzeltme '
        'ya hicbir sey degistirmiyor ya da degeri YUKARI cekiyor (kapsam eksik olculmustu).',
        'Sebep tek bir rakip kutu: A1-4_3078083 (Methanosarcina hadiensis, n=2215). Eski motorla %0,72 '
        'olculuyordu; dogru deger %47,22. Ikinci kutu A2-4_3078083 %4,49 -> %33,33. Uye kutulari yalniz '
        '+1,4 ile +2,2 puan oynuyor. Yani duzeltme neredeyse tamamen RAKIP tarafta - yani hata '
        'ozgullugu OLDUGUNDAN IYI gosteriyordu. Tehlikeli yon.',
        'ESKI (mm<=1): uye %40,63-60,63 / rakip maks %4,49 / en kotu kat 4,23x / havuz kat 49,96x. '
        'Panelde yayimlanan deger 187,9x (cins ici havuz - daha dar havuz, ayni hatayi tasiyor).',
        'YENI (mm<=1): uye %42,23-62,20 / rakip maks %47,22 / EN KOTU KAT 0,82x / havuz kat 11,41x.',
        'ANLAMI: bu cift Methanosarcina hadiensis\'i M. mazei kadar iyi cogaltiyor. "M. mazei / '
        'M. soligelidi grubu" iddiasi bu haliyle TASINMIYOR: ya cogalttigi kume M. hadiensis\'i de '
        'kapsayacak sekilde genisletilmeli, ya cift dusurulmeli, ya da amplikon dizilemesi sarti '
        'konmalidir. KARAR SIZIN - panel bu satiri "ESIK ALTI" olarak isaretledi, sessizce birakilmadi.',
    ]:
        yaz(ws, n, 1, s, fill=(KIRMIZI if s.startswith('ANLAMI') else None))
        ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=9)
        ws.row_dimensions[n].height = 46; n += 1
    n += 1

    yaz(ws, n, 1, u'7. ROWS THAT IMPROVED (their coverage had been under measured)', bold=True, fill=YESIL); n += 1
    for s in [
        'Methanosarcina_cinsi (satir 7): yedi M. mazei UYE kutusu eski motorla %0,5-26,5 olculuyordu; '
        'dogru deger %79,4-82,9. En kotu kat 0,04x -> 4,37x (mm<=1) / 4,66x (mm<=3), havuz 2,51x -> 81,59x. '
        'Panelde yayimlanan 28,4x havuz olcusudur; EN KOTU TEK KUTU olcusu hala 10x altinda.',
        'Asetoklastik_metanojenler (satir 16): uye alt sinir %0,0 -> %58,6; en kotu kat ~0 -> 4,22x, '
        'havuz ~0 -> 50,84x.',
        'Kapsam olculeri: Arke_universal 11/39 -> 32/39 (mm<=1; panelin 39/39 iddiasi mm<=3 degeridir), '
        'Bakteri_universal 4/20 -> 13/20 (mm<=1), Mantar_universal F1 14/20 -> 16/20.',
        'Methanothrix_cinsi (satir 3): mm<=1\'de degismiyor (13,54x -> 13,74x) ama mm<=3\'te 0,86x\'e '
        'dusuyor (bir rakip kutu %76,92). Bu satirin kaderi tamamen OLCUT SECIMINE bagli - olcut '
        'kararlastirilmadan siparise gitmemeli.',
    ]:
        yaz(ws, n, 1, s); ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=9)
        ws.row_dimensions[n].height = 46; n += 1
    n += 1

    yaz(ws, n, 1, '8. BAGIMSIZ DOGRULAMA (duzeltilmis motoru KULLANMAYAN yol)', bold=True, fill=GRI); n += 1
    for s in [
        'Baslik sayilari UC ayri yolla olculdu: (A) ispcr.find_sites - numpy vektor tarama, tohumsuz, '
        'panelin KENDI kodu, bu oturumda degistirilmedi; (B) brute_force.py - saf python, her pozisyon '
        'tek tek, ortak kod paylasmaz; (C) duzeltilmis read_engine.py.',
        'Yedi baslik kutusunun YEDISINDE de ucu ayni sayiyi verdi: A1-4_3078083 1046/2215 (%47,22), '
        'A2-4_3078083 52/156 (%33,33), A2-2_2209 1866/3000 (%62,20), A1-3_2209 1267/3000 (%42,23), '
        'A1-2_2209 (Methanosarcina_cinsi) 2389/3000 (%79,63), A2-3_2223 15/199 (%7,54), '
        'A1-2_2209 (Asetoklastik) 1828/3000 (%60,93).',
        'Kosmak icin: python screening/independent_check.py --fastq "fastq files"',
    ]:
        yaz(ws, n, 1, s); ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=9)
        ws.row_dimensions[n].height = 46; n += 1
    n += 1

    yaz(ws, n, 1, '9. ACIK KALAN - SONRAKI KISININ BILMESI GEREKENLER', bold=True, fill=SARI); n += 1
    for s in [
        'TAM DERINLIK: bu sayfadaki sayilar kutu basina <=3000 okuma alt kumesiyle uretildi. '
        'python3 -m screening --mod panel-olc --full-depth runs with no sampling. '
        'Egilim degisirse bu sayfa guncellenmelidir.',
        'UYE KUMELERI: panelin ozgun uye takson kumeleri hicbir betikte saklanmamis (komut satiri '
        'argumaniydi). screening/ciftler.tsv\'de hedef adlarindan ve taxid_adlari.tsv\'den yeniden '
        'kuruldu. "PANELLE_TUTUYOR" isaretli satirlarda yeniden kurulan kume panelin yayimladigi uye % '
        'araligini yeniden uretiyor; "YENIDEN_KURULDU" isaretlilerde uretmiyor - o satirlarda eski<->yeni '
        'FARKI gecerlidir (iki motor ayni kutularda kosuldu) ama mutlak degerler panelin ozgun kume '
        'tanimiyla ortusmeyebilir. Ozellikle Proteiniphilum (satir 10) ve Bacteroidales (satir 12) '
        'kumeleri kontrol edilmelidir.',
        'OLCUT KARARI: panelin tek bir numune olcutu secmesi gerekiyor. Onerilen mm<=1 (yayimlanan olcut '
        'etiketiyle tutarli olan); mm<=3 tasarim boru hattinin olcutudur ve daha gevsektir. Karar '
        'verilmeden "Ayrim (x)" sutunu satirlar arasi karsilastirilamaz.',
        'KONSENSUS: bu duzeltme HAM OKUMA olcumleridir; konsensus dosyalarina dokunulmadi. Konsensus '
        'yeniden uretimi ayrica yapilacaktir (iki yontemle: kalite agirlikli esigi dusurulmus + cogunluk '
        'oyu; ayrilan pozisyonlar maskelenecek, dejenere baz KULLANILMAYACAK - toplanti karari).',
    ]:
        yaz(ws, n, 1, s); ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=9)
        ws.row_dimensions[n].height = 56; n += 1
    return ws


def panel_sutunlari(wb, rows):
    ws = wb['2 Panel']
    c0 = ws.max_column + 1
    basliklar = ['OKUMA MOTORU DUZELTMESI - olcut etiketi',
                 'DUZELTILMIS ayrim (mm<=1) en kotu / havuz',
                 'DUZELTILMIS ayrim (mm<=3) en kotu / havuz',
                 'DEGISTI MI',
                 '10x ESIGININ ALTINDA MI']
    for j, h in enumerate(basliklar):
        yaz(ws, 1, c0 + j, h, bold=True, fill=GRI)
    ix = {int(r['satir']): r for r in rows}
    for satir in range(2, 23):
        r = ix.get(satir)
        if not r:
            continue
        deg = bool(r['DEGISTI'])
        alt = r['ESIK_ALTI_mm1'] == 'EVET'
        f = KIRMIZI if (deg and alt) else (SARI if deg else None)
        yaz(ws, satir, c0 + 0, r['panel_olcut_tespiti'], fill=f)
        yaz(ws, satir, c0 + 1, '%s / %s' % (r['YENI1_kat_enkotu'] or '-', r['YENI1_kat_havuz'] or '-'), fill=f)
        yaz(ws, satir, c0 + 2, '%s / %s' % (r['YENI3_kat_enkotu'] or '-', r['YENI3_kat_havuz'] or '-'), fill=f)
        yaz(ws, satir, c0 + 3, r['DEGISTI'] or 'hayir', fill=f)
        yaz(ws, satir, c0 + 4, ('EVET (mm<=1)' if alt else '') +
            (' EVET (mm<=3)' if r['ESIK_ALTI_mm3'] == 'EVET' else '') or 'hayir', fill=f)
    for j, w in enumerate((30, 26, 26, 60, 22)):
        ws.column_dimensions[openpyxl.utils.get_column_letter(c0 + j)].width = w
    # NOT satiri
    n = ws.max_row + 2
    yaz(ws, n, 1, 'NOT', bold=True)
    yaz(ws, n, 3, u'READ ENGINE FIX (2026-08-02, the most recent item): a SILENT bug was found in the engine that produced the "Uye urun %", "Rakip maks %" and "Ayrim (x)" columns. The 13 base EXACT matching seed never saw a binding site when the mismatch fell inside the seed (174/400 instead of 2/400 under the same criterion). The engine was fixed and all 21 pairs were re-measured. ALSO: those three columns were not measured under a SINGLE CRITERION. Some rows used mm<=1 and some mm<=3, so values inside a column cannot be compared across rows. Detail, the changed values and the pairs that fell below 10x: the "16 Okuma Motoru Duzeltmesi" sheet. The new columns (see to the right) are a subset measurement; full depth will be verified with python3 -m screening --mod panel-olc --full-depth.')
    ws.cell(n, 3).fill = KIRMIZI
    ws.cell(n, 3).alignment = SAR
    ws.row_dimensions[n].height = 90


def ozet_ve_karar(wb):
    ws = wb['1 Rapora Ozet']
    n = ws.max_row + 2
    yaz(ws, n, 1, 'OKUMA MOTORU DUZELTMESI (2026-08-02 - EN ONEMLI ACIK MADDE)', bold=True, fill=KIRMIZI)
    yaz(ws, n, 7, u'A silent bug was found in the engine that produced the discrimination ratios in this table: the 13 base exact matching seed never saw a binding site when the mismatch fell inside the seed (174/400 instead of 2/400 under the same criterion, a 98.9% miss rate). The engine was fixed and all 21 pairs were re-measured. THE ONE REAL REGRESSION: the Methanosarcina mazei / M. soligelidi group. Pool discrimination went 49.96x -> 11.41x, and the WORST SINGLE BIN 4.23x -> 0.82x. The cause: bin A1-4_3078083 of M. hadiensis amplifies at 47.22% instead of 0.72%. In other words this pair amplifies M. hadiensis as well as it amplifies the target. THE 187.9x VALUE ON THIS ROW IS NO LONGER VALID. On the other rows the fix is either neutral or pulls the value UP (their coverage had been under measured). Detail: "16 Okuma Motoru Duzeltmesi".')
    ws.cell(n, 7).fill = KIRMIZI; ws.cell(n, 7).alignment = SAR
    ws.row_dimensions[n].height = 120
    n += 1
    yaz(ws, n, 1, 'IKINCI BULGU - OLCUT TUTARSIZLIGI', bold=True, fill=SARI)
    yaz(ws, n, 7, u'The panel\'s "Ayrim (x)" column was not produced under a SINGLE CRITERION: some rows were measured with mismatches <=1 and some with <=3 (for example Proteiniphilum reads 29.0% in the panel, which is the mm<=3 value; under mm<=1 it is 1.6%). Values inside the column therefore CANNOT be compared across rows. On the "16 Okuma Motoru Duzeltmesi" sheet every row was re-measured with one engine under both criteria, and every value carries its criterion label. The panel needs to settle on a single sample criterion, and that decision is yours.')
    ws.cell(n, 7).fill = SARI; ws.cell(n, 7).alignment = SAR
    ws.row_dimensions[n].height = 90

    # M. mazei satirlari
    for r in (8, 23):
        h = ws.cell(r, 6)
        h.value = 'GECERSIZ: 187,9x -> duzeltilmis 11,41x (havuz) / 0,82x (en kotu tek kutu). ' \
                  'Bkz. "16 Okuma Motoru Duzeltmesi".'
        h.fill = KIRMIZI; h.alignment = SAR
    h = ws.cell(9, 6)   # Methanosarcina cinsi
    h.value = '28,4x -> duzeltilmis havuz 81,59x, en kotu tek kutu 4,66x (10x ALTINDA). ' \
              'Kapsam eksik olculmustu, duzeltme degeri yukari cekti.'
    h.fill = SARI; h.alignment = SAR

    ws = wb['6 Karar Durumu']
    n = ws.max_row + 2
    yaz(ws, n, 1, 'OKUMA MOTORU DUZELTMESI', bold=True, fill=KIRMIZI)
    yaz(ws, n, 2, 'Methanosarcina mazei / M. soligelidi grubu')
    yaz(ws, n, 3, 'ESIK ALTI - YENIDEN KARAR GEREKIYOR', fill=KIRMIZI)
    yaz(ws, n, 4, u'The seed bug in the sample engine was fixed. Within genus pool discrimination went 49.96x -> 11.41x, and the worst single bin 4.23x -> 0.82x. The cause: bin A1-4_3078083 of M. hadiensis measured 0.72% under the old engine, and its correct value is 47.22%. The pair amplifies M. hadiensis as well as it amplifies the target. OPTIONS: (a) widen the amplified set so that it covers M. hadiensis as well, (b) drop the pair, (c) keep it conditional on amplicon sequencing. It was not left in the panel silently; it is marked "ESIK ALTI". Detail: "16 Okuma Motoru Duzeltmesi".')
    ws.cell(n, 4).fill = KIRMIZI; ws.cell(n, 4).alignment = SAR
    ws.row_dimensions[n].height = 90
    n += 1
    yaz(ws, n, 1, 'OKUMA MOTORU DUZELTMESI', bold=True, fill=YESIL)
    yaz(ws, n, 2, 'Methanosarcina cinsi / Asetoklastik metanojenler')
    yaz(ws, n, 3, u'IMPROVED (the coverage had been under measured)', fill=YESIL)
    yaz(ws, n, 4, u'Methanosarcina_cinsi: seven M. mazei member bins measured 0.5-26.5% under the old engine, and their correct value is 79.4-82.9%. Pool discrimination 2.51x -> 81.59x; the worst single bin 0.04x -> 4.66x (still below 10x). Asetoklastik_metanojenler: the member floor went 0.0% -> 58.6%, and the pool ~0 -> 50.84x. The coverage measures also rose: Arke_universal 11/39 -> 32/39 (mm<=1), Bakteri_universal 4/20 -> 13/20.')
    ws.cell(n, 4).fill = YESIL; ws.cell(n, 4).alignment = SAR
    ws.row_dimensions[n].height = 76
    n += 1
    yaz(ws, n, 1, 'OKUMA MOTORU DUZELTMESI', bold=True, fill=SARI)
    yaz(ws, n, 2, 'Methanothrix cinsi')
    yaz(ws, n, 3, 'OLCUTE ASIRI DUYARLI - OLCUT KARARI BEKLIYOR', fill=SARI)
    yaz(ws, n, 4, u'Under the mm<=1 criterion the discrimination is 13.74x (above the threshold); under mm<=3 it is 0.86x (one competitor bin amplifies at 76.92%). The 75.2x published in the panel could not be reproduced under either criterion. This should not go to order before the criterion is settled.')
    ws.cell(n, 4).fill = SARI; ws.cell(n, 4).alignment = SAR
    ws.row_dimensions[n].height = 60


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--xlsx', required=True)
    ap.add_argument('--tsv', required=True)
    ap.add_argument('--backup', '--yedek', dest='yedek', default='')
    a = ap.parse_args()
    yedek = a.yedek or a.xlsx.replace('.xlsx', '_YEDEK_motor_oncesi.xlsx')
    if not os.path.exists(yedek):
        shutil.copy2(a.xlsx, yedek)
        print('yedek:', yedek)
    rows = list(csv.DictReader(open(a.tsv, encoding='utf-8'), delimiter='\t'))
    wb = openpyxl.load_workbook(a.xlsx)
    sayfa16(wb, rows)
    panel_sutunlari(wb, rows)
    ozet_ve_karar(wb)
    wb.save(a.xlsx)
    print('guncellendi:', a.xlsx, '| sayfalar:', len(wb.sheetnames))


if __name__ == '__main__':
    main()
