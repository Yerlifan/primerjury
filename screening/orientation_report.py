# -*- coding: utf-8 -*-
"""orientation_report.py - panele "19 Yon Normalizasyonu" sayfasini yazar."""
# ---------------------------------------------------------------------------
# orientation_report.py — yon normalizasyonunun butun kanitlarini teslim panelinin
#                  xlsx dosyasina "19 Yon Normalizasyonu" sayfasi olarak yazar.
#
# GIRDI  : --xlsx teslim paneli; --yon ile orientation_audit.py'nin urettigi dosya
#          bazli yon tablosu; --kod ile orientation_code_scan.py'nin urettigi kod
#          yolu siniflandirmasi; --indeks ile konsensus_kanonik/INDEKS.tsv.
#          Hangi dosyanin nasil duzeltildigine dair notlar (KOD_NOT) bu dosyanin
#          icinde sabittir.
# CIKTI  : verilen xlsx dosyasina yeni bir sayfa ekleyip kaydeder (wb.save);
#          ekrana yazilan sayfa adini ve satir sayisini basar.
# CAGRAN : MENUDE DEGILDIR - elle calistirilir, cunku teslim dosyasini
#          degistirir. Girdilerini uretecek iki betik de (orientation_audit.py,
#          orientation_code_scan.py) elle calistirilir.
#
# Sayfanin varlik sebebi: yon hatasi ayni gece uc ayri yerde ayri ayri bulunup
# yamandi. Bu sayfa "hangi dosya duzeltildi, hangi konsensus cevrildi, hangi
# kod yolu hala risktedir" sorularinin cevabini panelin icinde, izlenebilir
# bicimde tutar.
# ---------------------------------------------------------------------------
import os, sys, csv, argparse, collections
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment

KIRMIZI = PatternFill('solid', fgColor='FFC7CE')
SARI = PatternFill('solid', fgColor='FFEB9C')
YESIL = PatternFill('solid', fgColor='C6EFCE')
GRI = PatternFill('solid', fgColor='D9D9D9')
KALIN = Font(bold=True)
SAR = Alignment(wrap_text=True, vertical='top')

KOD_NOT = {
 'screening/config.py':
   'DUZELTILDI. KONSENSUS artik konsensus_kanonik. "consensus sequences" yalniz '
   'KONSENSUS_HAM olarak duruyor ve SADECE kanonik uretimi okur.',
 'screening/targets.py':
   'DUZELTILDI. konsensusler() artik konsensus_kanonik/INDEKS.tsv okuyor. Indeks '
   'yoksa HATA verir - karisik klasore sessizce DUSMEZ.',
 'screening/build_consensus.py':
   'DUZELTILDI (uc yerde). (a) sablon kanonige cevriliyor, (b) oksuz kutu sablonu '
   '(ham okuma) da cevriliyor, (c) CIKTI yazilmadan once bir kez daha olculup '
   'cevriliyor. Ayrica kosunun basina KAPI kondu: yon sinamasi gecmeden uretim baslamaz.',
 'screening/self_test.py':
   'yon_sinamasi() eklendi: orientation.py kendi sinavi + kanonik sette ters dosya var mi + '
   'bilinen cevapli etki testi. primer3 gibi opsiyonel bagimliliga TAKILMAZ.',
 'screening/generator.py':
   'Yon bagimsiz: diziyi hem duz hem rc ile deniyor. Kanonik kaynakla da dogru calisir.',
 'screening/panel_measurement.py':
   'HAM OKUMA olcer, konsensus yonunden etkilenmez (okumalar zaten iki yonde taranir).',
 'screening/membership_check.py':
   'Uyelik/olcum; konsensus yolunu targets.py uzerinden alir - kanonige bagli.',
 'screening/orientation.py':
   'YENI. Kanonik yon TANIMI ve normalizasyonu. Iki bagimsiz olcut, kendini sinama.',
 'screening/build_canonical.py':
   'YENI. konsensus_kanonik/ uretir + MANIFEST/INDEKS/BELIRSIZ yazar + kendi dogrulamasi.',
 'screening/orientation_audit.py':
   'YENI (denetim araci). Karisik klasorleri BILEREK okur - amaci yonu olcmek.',
 'screening/orientation_impact_test.py':
   'YENI. Bilinen cevapli test: dogru yon vs ters yon urun sayisi.',
 'steps/split_clusters.py':
   'ESKI HAT. Karisik klasoru tek yon varsayarak okuyor. Panelin mevcut sayilarini '
   'URETMIYOR (devre disi hat) ama yeniden kosulursa yanlis sonuc verir - kanonige '
   'cevrilmeli ya da arsivlenmeli.',
}


def yaz(ws, r, c, v, fill=None, bold=False):
    h = ws.cell(r, c, v)
    if fill: h.fill = fill
    if bold: h.font = KALIN
    h.alignment = SAR
    return h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--xlsx', required=True)
    ap.add_argument('--orientation', '--yon', dest='yon', required=True)
    ap.add_argument('--code', '--kod', dest='kod', required=True)
    ap.add_argument('--index', '--indeks', dest='indeks', required=True)
    a = ap.parse_args()

    yr = list(csv.DictReader(open(a.yon, encoding='utf-8'), delimiter='\t'))
    kr = list(csv.DictReader(open(a.kod, encoding='utf-8'), delimiter='\t'))
    ix = list(csv.DictReader(open(a.indeks, encoding='utf-8'), delimiter='\t'))

    say = collections.defaultdict(collections.Counter)
    for x in yr:
        say[x['klasor']][x['karar'].split(' (')[0]] += 1

    wb = openpyxl.load_workbook(a.xlsx)
    ad = '19 Yon Normalizasyonu'
    if ad in wb.sheetnames:
        del wb[ad]
    ws = wb.create_sheet(ad)
    for w, c in zip((44, 14, 14, 14, 14, 16, 60), 'ABCDEFG'):
        ws.column_dimensions[c].width = w
    n = 1

    yaz(ws, n, 1, 'YON NORMALIZASYONU - 2026-08-02', bold=True); n += 1
    yaz(ws, n, 1, u'ANSWER: the orientation IS NOW CANONICAL. It was not before: there were four separate consensus sets and two of them were mixed orientation. A single', fill=YESIL)
    ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=7)
    ws.row_dimensions[n].height = 32; n += 2

    yaz(ws, n, 1, u'1. WHY IT WAS ASKED - three separate patches, no single canonical fix', bold=True, fill=GRI); n += 1
    for s in ['The orientation bug was found and patched separately in AT LEAST THREE SEPARATE PLACES over the course of one night: (a) the strand choice of 85 on the source study side, (b) "the consensuses are the reverse complement of SILVA" on the design side, (c) "39 of 58 consensuses are in the reverse orientation" during the B/F re-measurement. Three separate patches means there is no single canonical fix, which means it escapes again on the next change. This page ties "I think it is fixed" to a measurement.']:
        yaz(ws, n, 1, s); ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=7)
        ws.row_dimensions[n].height = 46; n += 1
    n += 1

    yaz(ws, n, 1, u'2. MEASUREMENT - the orientation of every consensus file BEFORE the fix', bold=True, fill=GRI); n += 1
    yaz(ws, n, 1, u'Method: two INDEPENDENT criteria. (1) the panel\'s own universal pairs: in the sense orientation F and rc(R) are found directly' % len(yr))
    ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=7)
    ws.row_dimensions[n].height = 46; n += 1
    for j, h in enumerate(['Klasor', 'SENSE', 'ANTISENSE', 'BELIRSIZ/bos', 'Toplam', 'Karisik mi?', 'Not'], 1):
        yaz(ws, n, j, h, bold=True, fill=GRI)
    n += 1
    NOT = {
     'consensus sequences': 'the original output of the source study. MIXED. Cause: the samtools/minimap2 line '
                            'sablonu veriden secilen bir OKUMAYA gore kuruyor; nanopore okumalari '
                            '~50/50 iki yonde geldigi icin cikti yonu rastgele.',
     'referans_konsensus/konsensus': 'Gece normalize edilen set. Temiz.',
     'referans_konsensus/baskin/konsensus': 'Baskin alel seti. Temiz (6 dosya olculemedi - N orani yuksek).',
     'referans_konsensus/self/konsensus': 'Self set. Temiz (1 dosya bos).',
     'KAPSAMLI_ARAMA_SONUC/konsensus_yeni': 'BU GECE KOSULACAK URETIMIN CIKTISI. KARISIK ve '
                                            'agirlikli ANTISENSE. Kok neden: sablonu karisik '
                                            'klasorden aliyordu (asagi bakiniz).',
    }
    for k in sorted(say):
        c = say[k]
        se, an = c.get('SENSE', 0), c.get('ANTISENSE', 0)
        bl = sum(v for kk, v in c.items() if kk not in ('SENSE', 'ANTISENSE'))
        kar = 'EVET - KARISIK' if se and an else 'hayir'
        f = KIRMIZI if (se and an) else None
        yaz(ws, n, 1, k, fill=f); yaz(ws, n, 2, se, fill=f); yaz(ws, n, 3, an, fill=f)
        yaz(ws, n, 4, bl, fill=f); yaz(ws, n, 5, se + an + bl, fill=f)
        yaz(ws, n, 6, kar, fill=f); yaz(ws, n, 7, NOT.get(k, ''), fill=f)
        ws.row_dimensions[n].height = 44
        n += 1
    n += 1

    yaz(ws, n, 1, u'3. THE ROOT CAUSE (the evidence in the code)', bold=True, fill=KIRMIZI); n += 1
    for s in ['screening/config.py:  KONSENSUS = y(\'consensus sequences\')  '
              '-> paket KARISIK klasoru tek kaynak sayiyordu.',
              'screening/build_consensus.py -> _sablon_sec(): the template was taken from the "current consensus". The reads were anchored to the template in both directions and NORMALISED (the code called this "the orientation is normalised"), but THE TEMPLATE\'S OWN ORIENTATION was normalised nowhere. The result: the output inherits the template\'s orientation exactly. Measured evidence: konsensus_yeni was 28 antisense / 7 sense.',
              "For orphan bins the template was chosen straight from a RAW READ, and a read's orientation is random, so the output orientation was random too."]:
        yaz(ws, n, 1, s, fill=KIRMIZI if 'konsensus_uret' in s else None)
        ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=7)
        ws.row_dimensions[n].height = 46; n += 1
    n += 1

    yaz(ws, n, 1, u'4. A TEST WITH A KNOWN ANSWER - how far each number moves when the orientation is wrong',
        bold=True, fill=GRI); n += 1
    yaz(ws, n, 1, u'Setup: the same primer pair, the same consensus, THE ONLY DIFFERENCE BEING ORIENTATION. Source: referans_konsensus/konsensus (99 files')
    ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=7)
    ws.row_dimensions[n].height = 46; n += 1
    for j, h in enumerate(['Cift', 'Sinif', 'Dogru yon', 'TERS yon', 'Kayip', 'Sonuc', ''], 1):
        yaz(ws, n, j, h, bold=True, fill=GRI)
    n += 1
    for ad2, sn, d, t, top in [('Arke_universal', 'A', 38, 0, 39),
                               ('Bakteri_universal', 'B', 20, 0, 20),
                               ('Mantar_universal (F1)', 'F1', 18, 0, 20),
                               ('Mantar_universal (F2)', 'F2', 15, 0, 20),
                               ('Methanosarcina_mazei_turu', 'A', 12, 0, 39),
                               ('Methanosarcina_cinsi', 'A', 13, 0, 39),
                               ('Proteolitik_Cloacimonas', 'B', 1, 0, 20)]:
        yaz(ws, n, 1, ad2, fill=KIRMIZI); yaz(ws, n, 2, sn, fill=KIRMIZI)
        yaz(ws, n, 3, '%d/%d' % (d, top), fill=KIRMIZI)
        yaz(ws, n, 4, '%d/%d' % (t, top), fill=KIRMIZI)
        yaz(ws, n, 5, d - t, fill=KIRMIZI)
        yaz(ws, n, 6, u'the reverse orientation ZEROES the product', fill=KIRMIZI)
        n += 1
    yaz(ws, n, 1, u'TOTAL: 117 products in the correct orientation, 0 in the reverse orientation. LOSS 100%.', bold=True, fill=KIRMIZI)
    ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=7); n += 1
    yaz(ws, n, 1, u'CROSS-CHECK: the same test was repeated with the panel\'s OWN engine (ispcr.amplify). Arke_universal gave 38 correct')
    ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=7)
    ws.row_dimensions[n].height = 30; n += 2

    yaz(ws, n, 1, u'5. THE CANONICAL FIX - one source, one definition', bold=True, fill=YESIL); n += 1
    for s in ['THE CANONICAL ORIENTATION = SENSE (the reference, or plus, strand). The definition lives in ONE PLACE: screening/orientation.py. Two independent criteria, and if they disagree the file counts as UNCERTAIN, is NOT normalised, and is flagged instead (the project rule: no decision is left to a single code path).',
              'TEK KAYNAK: konsensus_kanonik/ - %d kutu, hepsi SENSE. Uretici: '
              'screening/build_canonical.py. Yaninda MANIFEST.tsv (her dosyanin kaynagi, '
              'eski yonu, cevrilip cevrilmedigi), INDEKS.tsv (gecerli dosya listesi) ve '
              'BELIRSIZ.tsv. Bu koside %d dosya ANTISENSE -> SENSE cevrildi, BELIRSIZ 0.'
              % (len(ix), sum(1 for r in ix if r['cevrildi'] == 'EVET')),
              'EVERY SCRIPT READS THIS PLACE: targets.py -> konsensusler() now reads INDEKS.tsv and RAISES AN ERROR when the index is missing; it does NOT fall back SILENTLY to the mixed directory. No script carries its own orientation patch any more.',
              'CAUTION - files CANNOT BE DELETED on a mounted drive. konsensus_kanonik/ still holds the misnamed "*_kanonik.fasta" leftovers of the first run. The valid files match the pattern "*.kanonik.fa" and are listed in INDEKS.tsv. Consumers read the INDEX, NOT a glob. The OKUBENI.txt in the directory repeats this.']:
        yaz(ws, n, 1, s, fill=SARI if s.startswith('DIKKAT') else None)
        ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=7)
        ws.row_dimensions[n].height = 52; n += 1
    n += 1

    yaz(ws, n, 1, u'6. CODE PATHS - every script that reads a consensus, and how it handles orientation',
        bold=True, fill=GRI); n += 1
    yaz(ws, n, 1, u'Automatic scan (orientation_code_scan.py, with comments and docstrings stripped) plus a manual note. Full list: yon_kod_tara' % len(kr))
    ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=7)
    ws.row_dimensions[n].height = 26; n += 1
    for j, h in enumerate(['Betik', 'Sinif', 'Karisik kaynak', 'Normalize kaynak',
                           'Iki yon deniyor', 'Kendi yamasi', 'Not'], 1):
        yaz(ws, n, j, h, bold=True, fill=GRI)
    n += 1
    onem = {'HAM_KAYNAK_VARSAYIM': 0, 'KAYNAK_BELIRSIZ': 1, 'NORMALIZE_KAYNAK': 2,
            'IKI_YON_DENIYOR': 3, 'KENDI_YAMASI': 4, 'YON_ILGISIZ': 5}
    for r in sorted(kr, key=lambda x: (onem.get(x['sinif'], 9), x['betik'])):
        if r['betik'].startswith('screening/eski/'):
            continue
        f = (KIRMIZI if r['sinif'] == 'HAM_KAYNAK_VARSAYIM' and 'yapilandirma' not in r['betik']
             and 'yon_denetimi' not in r['betik'] else
             YESIL if r['betik'] in KOD_NOT and KOD_NOT[r['betik']].startswith(('DUZELTILDI', 'YENI')) else None)
        yaz(ws, n, 1, r['betik'], fill=f); yaz(ws, n, 2, r['sinif'], fill=f)
        yaz(ws, n, 3, r['karisik_kaynak'], fill=f); yaz(ws, n, 4, r['normalize_kaynak'], fill=f)
        yaz(ws, n, 5, r['iki_yon'], fill=f); yaz(ws, n, 6, r['kendi_yamasi'], fill=f)
        yaz(ws, n, 7, KOD_NOT.get(r['betik'], ''), fill=f)
        if KOD_NOT.get(r['betik']):
            ws.row_dimensions[n].height = 40
        n += 1
    n += 1

    yaz(ws, n, 1, u'7. THE PRODUCTION RUN TONIGHT - checked specifically', bold=True, fill=YESIL); n += 1
    for s in ['build_consensus.py was fixed in THREE places: (a) the template is converted to canonical, (b) the orphan bin template (a raw read) is converted too, (c) the OUTPUT is measured once more and converted before it is written (the last seat belt). The output fasta header carries kanonik=... cevrildi=..., and the columns cikti_yon and cikti_cevrildi were added to the TSV.',
              'A GATE WAS PUT IN: konsensus_uret.calistir() runs the orientation test first, and if it does not pass, GENERATION IS NOT STARTED and what to do is printed. The reason: if the output of this step comes out in the wrong orientation, the whole night is wasted.',
              "yon_sinamasi() was added inside self_test.py (3 items: orientation.py's own test, whether the canonical set still holds a reversed file, and a known-answer impact test). It does NOT TRIP over optional dependencies such as primer3, and it can be called on its own. It was run in this session: 3/3 PASSED (correct 39/40, reversed 0/40).",
              'THE ORDER: (1) python screening/build_canonical.py --root .   (2) consensus regeneration from the verification/full_chain.py menu   (3) once generation finishes, run kanonik_uret AGAIN with --oncelik yeni so that the new set becomes the canonical source.']:
        yaz(ws, n, 1, s, fill=SARI if s.startswith('SIRA') else None)
        ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=7)
        ws.row_dimensions[n].height = 56; n += 1
    n += 1

    yaz(ws, n, 1, '8. ACIK KALAN', bold=True, fill=SARI); n += 1
    for s in ['konsensus_kanonik currently comes from referans_konsensus/konsensus (99 bins) plus 1 orphan bin (A1-1_2209, from the original directory, which was ANTISENSE and was converted). konsensus_yeni was NOT taken into the canonical set because it is INCOMPLETE (35/99); mixing a half set with a full one creates a heterogeneous base. Once the overnight generation finishes it should be re-run with --oncelik yeni.',
              "steps/split_clusters.py reads the mixed directory assuming a single orientation. It does not produce the panel's CURRENT numbers (it is a disabled line), but if it is re-run it gives a wrong answer. It should either be converted to canonical or archived.",
              'The scripts flagged KAYNAK_BELIRSIZ in yon_kod_taramasi (most of them under steps and engine) take the consensus path from the command line and leave the orientation to the caller. Those are the old line. If they are to be re-run, the canonical directory must be given.',
              'The orientation measurement on this page applies to the consensus files. RAW READ measurements are unaffected by orientation (reads are scanned in both directions anyway), so the numbers on the "16 Okuma Motoru Duzeltmesi" sheet are NOT AFFECTED by this finding.']:
        yaz(ws, n, 1, s)
        ws.merge_cells(start_row=n, start_column=1, end_row=n, end_column=7)
        ws.row_dimensions[n].height = 56; n += 1

    wb.save(a.xlsx)
    print('sayfa yazildi:', ad, '| satir', n)


if __name__ == '__main__':
    main()
