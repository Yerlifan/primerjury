# -*- coding: utf-8 -*-
r"""
orientation_audit.py - her konsensus dosyasinin SAKLANDIGI YONU tespit eder.

Neden gerek: yon hatasi gece boyunca EN AZ UC AYRI YERDE ayri ayri bulunup ayri ayri
yamandi (kaynak calismada 85'in iplik secimi, tasarim tarafinda "konsensusler SILVA'ya gore ters",
B/F yeniden olcumunde "58 konsensusten 39'u ters"). Uc ayri yama = tek kanonik cozum yok.
Bu betik "sanirim duzeltildi"yi olcuye baglar.

YONTEM - iki BAGIMSIZ olcut, ikisi ayrilirsa dosya BELIRSIZ sayilir (proje kurali):

  Olcut 1 (panelin kendi evrensel primerleri):
      Sense (referans/artı) yonde saklanan bir SSU konsensusunde ileri primer F
      DOGRUDAN, ters primerin tumleyeni rc(R) de DOGRUDAN bulunur.
      Antisense saklanmissa bunun yerine rc(F) ve R bulunur.
  Olcut 2 (literatur evrensel motifleri - panelden bagimsiz):
      SSU  : 515F  GTGYCAGCMGCCGCGGTAA  ve  806R bolgesi sense hali GGATTAGATACCC
      ITS  : ITS1  TCCGTAGGTGAACCTGCGG  ve  rc(ITS4)  GCATATCAATAAGCGGAGGA

Her olcut icin: motiflerin DUZ dizide mi yoksa TERS TUMLEYENINDE mi bulundugu sayilir.
Uyumsuzluk toleransi <=2 (nanopore konsensusunde tek tuk hata olur), kayipsiz motor.

Kullanim:
    python orientation_audit.py --kok ..  --out ..\yon_denetimi_20260802.tsv
    python orientation_audit.py --kok .. --klasor "consensus sequences"
"""
# ---------------------------------------------------------------------------
# orientation_audit.py — bir konsensus klasorundeki HER dosyanin saklandigi yonu iki
#                   bagimsiz olcutle tespit eder ve tabloya doker.
#
# GIRDI  : --kok altinda --klasor ile verilen bir ya da birden fazla konsensus
#          klasoru (varsayilanlar dosyanin icinde tanimli). Baglanma yeri
#          aramasi okuma_motoru.Sonda ile, <=2 uyumsuzluk toleransiyla yapilir.
# CIKTI  : --out ile verilen TSV dosyasi: dosya basina sinif, iki olcutun
#          sonucu ve nihai karar (SENSE / ANTISENSE / BELIRSIZ).
# CAGRAN : MENUDE DEGILDIR - elle calistirilan bir denetimdir. Urettigi tablo
#          orientation_report.py ile panele "19 Yon Normalizasyonu" sayfasi olarak
#          islenir. Kanonik uretimin kendisi build_canonical.py'de yapilir.
#
# IKI OLCUT NEDEN: karar tek bir kod yoluna birakilmaz. Olcut 1 panelin kendi
# evrensel primerlerini, olcut 2 panelden bagimsiz literatur motiflerini
# kullanir. Ikisi ayrilirsa dosya BELIRSIZ sayilir - yanlis yone cevrilmis bir
# konsensus, hic cevrilmemis kadar sessiz zarar verir.
# ---------------------------------------------------------------------------
import sys, os, glob, csv, argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import read_engine as om

# --- olcut 1: panelin evrensel ciftleri (sinifa gore) ------------------------
PANEL = {
    'A':  ('Arke_universal',        'CTGCGGTTTAATTGGATTCAACGC', 'GAACTGACGACGGCCATGC'),
    'B':  ('Bakteri_universal',     'ACAAGCGGTGGAGCATGTG',      'ACGACAGCCATGCAGCAC'),
    'F1': ('Mantar_universal (F1)', 'GGTTACCCGCTGAACTTAAGC',    'CGCTTCACTCGCCGTTAC'),
    'F2': ('Mantar_universal (F2)', 'GTGCATGGCCGTTCTTAGTTG',    'CAAACTTCCATCGGCTTGAGC'),
}

# --- olcut 2: literatur motifleri, sense (artı) iplikte beklenen -------------
MOTIF = {
    'SSU': ['GTGYCAGCMGCCGCGGTAA',      # 515F
            'GGATTAGATACCC',            # 806R bolgesi, sense hali
            'AGTCCCGCAACGAGCGCAACCC'],  # 1100 bolgesi, sense (SSU korunmus)
    'ITS': ['TCCGTAGGTGAACCTGCGG',      # ITS1
            'GCATATCAATAAGCGGAGGA'],    # rc(ITS4)
}


def var_mi(dizi, desen, mm=2):
    """desen dizide (kayipsiz motor, <=mm uyumsuzluk) var mi"""
    s = om.Sonda(desen, uc5=False, max_mm=mm, son2=False)
    return len(s.bul(dizi)) > 0


def olcut1(dizi, sinif):
    ad, F, R = PANEL[sinif]
    duz = int(var_mi(dizi, F)) + int(var_mi(dizi, om.rc(R)))
    ters = int(var_mi(dizi, om.rc(F))) + int(var_mi(dizi, R))
    if duz > ters:
        return 'SENSE', duz, ters
    if ters > duz:
        return 'ANTISENSE', duz, ters
    return 'BELIRSIZ', duz, ters


def olcut2(dizi, tip):
    duz = sum(1 for m in MOTIF[tip] if var_mi(dizi, m))
    ters = sum(1 for m in MOTIF[tip] if var_mi(om.rc(dizi), m))
    if duz > ters:
        return 'SENSE', duz, ters
    if ters > duz:
        return 'ANTISENSE', duz, ters
    return 'BELIRSIZ', duz, ters


def sinifi(yol):
    b = os.path.basename(yol)
    for s in ('A1-', 'A2-'):
        if s in b or s in yol:
            return 'A'
    if 'F1-' in b or 'F1-' in yol:
        return 'F1'
    if 'F2-' in b or 'F2-' in yol:
        return 'F2'
    if 'B-' in b or 'B-' in yol:
        return 'B'
    return '?'


def oku(yol):
    ad, buf = None, []
    with open(yol, errors='ignore') as fh:
        for l in fh:
            if l.startswith('>'):
                if ad is not None:
                    yield ad, ''.join(buf)
                ad, buf = l[1:].strip(), []
            else:
                buf.append(l.strip())
    if ad is not None:
        yield ad, ''.join(buf)


def tara(kok, klasorler):
    satirlar = []
    for kl in klasorler:
        for yol in sorted(glob.glob(os.path.join(kok, kl, '**', '*.fasta'), recursive=True)):
            sn = sinifi(yol)
            if sn == '?':
                continue
            tip = 'ITS' if sn in ('F1', 'F2') else 'SSU'
            seqs = [om.temizle(s) for _, s in oku(yol)]
            if not seqs:
                satirlar.append(dict(klasor=kl, dosya=os.path.basename(yol), sinif=sn,
                                     uzunluk=0, N_yuzde='', olcut1='DOSYA_BOS', olcut2='DOSYA_BOS',
                                     karar='DOSYA_BOS', ayrinti=''))
                continue
            s = max(seqs, key=len)
            if not s:
                satirlar.append(dict(klasor=kl, dosya=os.path.basename(yol), sinif=sn,
                                     uzunluk=0, N_yuzde='', olcut1='DIZI_BOS', olcut2='DIZI_BOS',
                                     karar='DIZI_BOS', ayrinti=os.path.relpath(yol, kok)))
                continue
            nN = s.count('N')
            o1, d1, t1 = olcut1(s, sn)
            o2, d2, t2 = olcut2(s, tip)
            if o1 == o2 and o1 != 'BELIRSIZ':
                karar = o1
            elif 'BELIRSIZ' in (o1, o2) and o1 != o2:
                karar = (o1 if o2 == 'BELIRSIZ' else o2) + ' (tek olcut)'
            elif o1 != o2:
                karar = 'AYRILIK - MASKELE'
            else:
                karar = 'BELIRLENEMEDI'
            satirlar.append(dict(klasor=kl, dosya=os.path.basename(yol), sinif=sn,
                                 uzunluk=len(s), N_yuzde=round(100.0 * nN / len(s), 1),
                                 olcut1='%s (%d/%d)' % (o1, d1, t1),
                                 olcut2='%s (%d/%d)' % (o2, d2, t2),
                                 karar=karar, ayrinti=os.path.relpath(yol, kok)))
    return satirlar


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', '--kok', dest='kok', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--dir', '--klasor', dest='klasor', action='append', default=[])
    a = ap.parse_args()
    kl = a.klasor or ['consensus sequences',
                      'referans_konsensus/konsensus',
                      'referans_konsensus/baskin/konsensus',
                      'referans_konsensus/self/konsensus',
                      'KAPSAMLI_ARAMA_SONUC/konsensus_yeni']
    satirlar = tara(a.kok, kl)
    if satirlar:
        with open(a.out, 'w', newline='', encoding='utf-8') as fh:
            w = csv.DictWriter(fh, fieldnames=list(satirlar[0].keys()), delimiter='\t')
            w.writeheader()
            for r in satirlar:
                w.writerow(r)
    ozet = {}
    for r in satirlar:
        ozet.setdefault(r['klasor'], {}).setdefault(r['karar'], 0)
        ozet[r['klasor']][r['karar']] += 1
    print('%-44s %s' % (u'DIRECTORY', u'THE DISTRIBUTION OF VERDICTS'))
    for k in sorted(ozet):
        print('%-44s %s' % (k, ', '.join('%s=%d' % x for x in sorted(ozet[k].items()))))
    print(u'\ntotal files:', len(satirlar), '| TSV:', a.out)


if __name__ == '__main__':
    main()
