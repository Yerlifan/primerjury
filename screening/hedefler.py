# -*- coding: utf-8 -*-
"""Targets — read FROM FILE, never hard-coded.

Loads target definitions and the "problem" list from the panel tables. There is
no built-in list: a hard-coded target set silently diverges from the tables the
rest of the pipeline reads, and the two then disagree without any error.

Consensus sequences are served ONLY from the canonical, orientation-normalised
directory. The raw consensus directory is mixed-orientation (measured: 71
antisense / 27 sense) and a reverse-oriented consensus yields ZERO products in
in-silico PCR, silently — measured loss 100%.

--- ozgun aciklama ---
Panelden hedefleri ve 'sorunlu' olanlari DOSYADAN okur. Sabit liste yoktur.

Sorunlu sayilma gerekceleri (hepsi panel dosyasindan turetilir):
  G  geometri ihlali          -> "GEOMETRI (toplanti kurali)" sutununda IHLAL
  K  kosullu / on karar       -> "Durum" ya da "SIPARIS EDILEBILIR MI" sutunu
  A  ayrim esik alti (<10x)   -> "Ayrim (x)" sutunu sayisal ve 10'un altinda
  U  urun boyu qPCR disinda   -> >250 bp (onerilmez) ya da 150-250 (kosullu)
  C  panelden cikarilmis      -> SIPARIS EDILEBILIR MI = HAYIR
  P  plaka ici jel cakismasi  -> "PLAKA ICI JEL AYRIMI" sutununda AYRILMAZ
"""
# ---------------------------------------------------------------------------
# hedefler.py — paketin veri okuma katmani: panel, uyelik, kutu ve konsensus
#               listelerini dosyalardan cikarir ve her hedef icin uye/rakip
#               kumesini cozer.
#
# GIRDI  : primer_final/devir_ciftleri_*.tsv (panel; yoksa yedegi),
#          steps/hedefler.tsv (karar tablosu, grup uyeligi),
#          steps/taxid_adlari.tsv, screening/hedef_uyelik.tsv
#          (elle duzenlenebilen acik uyelik), "fastq files/*/reads_*.fastq"
#          (kutular) ve konsensus_kanonik/INDEKS.tsv (kanonik konsensusler).
# CIKTI  : dosyaya yazmaz. panel_oku() (satirlar, yol); sorunlu_hedefler()
#          (sorunlular, panel, yol); kutular() ve konsensusler() sozluk
#          listeleri; hedef_baglami() omurga + uye/rakip kutu + uye/rakip
#          konsensus sozlugu dondurur.
# CAGRAN : __main__.py, panel_olcum.py, uyelik_denetimi.py, konsensus_uret.py,
#          kendini_sina.py ve hepsi.yon_kapisi icinden - yani 1, 2, 3, 4, 5, 6,
#          7, 8, 9 tuslarinin hepsinde. Disaridan da cagrilir:
#          protocol/tek_protokol_olc.py (tus P), verification/kurtarma_turu.py
#          (tus K), verification/kimlik_dogrulama.py (tus I),
#          verification/tum_kutu_kimlikleri.py (tus G).
# ---------------------------------------------------------------------------
import os, csv, re, glob
from . import yapilandirma as C
from . import motor

# ESIK TEK KAYNAKTAN GELIR: screening/yapilandirma.py -> ESIK_DCQ = 3.0
# Kat karsiligi 2 ** ESIK_DCQ = 8,00. Sabit sayi GOMULMEZ; dCq degisirse
# tek yerden degisir. Gerekce ve verim uyarisi o dosyada yazili.
AYRIM_ESIK = C.AYRIM_ESIK


# ---------------------------------------------------------------- panel
def _panel_yolu():
    for p in (C.PANEL_TSV, C.PANEL_TSV_YEDEK):
        if os.path.exists(p):
            return p
    raise SystemExit('HATA: panel TSV bulunamadi: %s' % C.PANEL_TSV)


def _sutun(basliklar, *onek):
    for i, b in enumerate(basliklar):
        for o in onek:
            if b.strip().lower().startswith(o.lower()):
                return i
    return None


def _sayi(s):
    """'13,4 (ozgulluk) / havuz 47,0' -> 13.4 ; sayi yoksa None."""
    if not s:
        return None
    m = re.search(r'(\d+(?:[.,]\d+)?)', s)
    return float(m.group(1).replace(',', '.')) if m else None


def panel_oku():
    yol = _panel_yolu()
    with open(yol, encoding='utf-8') as fh:
        r = list(csv.reader(fh, delimiter='\t'))
    h = r[0]
    ix = dict(
        plaka=_sutun(h, 'Plaka'), ta=_sutun(h, 'Ta '), hedef=_sutun(h, 'Hedef'),
        duzey=_sutun(h, 'Duzey'), sinif=_sutun(h, 'Amplikon sinifi'),
        F=_sutun(h, 'Ileri primer'), R=_sutun(h, 'Geri primer'),
        urun=_sutun(h, 'Urun (bp)'), uye=_sutun(h, 'Uye urun'),
        rakip=_sutun(h, 'Rakip maks'), ayrim=_sutun(h, 'Ayrim'),
        durum=_sutun(h, 'Durum'), siparis=_sutun(h, 'SIPARIS'),
        geo=_sutun(h, 'GEOMETRI'), jel=_sutun(h, 'PLAKA ICI JEL'),
        not_=_sutun(h, 'Not'),
    )
    out = []
    for row in r[1:]:
        if len(row) <= ix['hedef']:
            continue
        if (row[ix['plaka']] or '').strip().upper() == 'NOT':
            continue
        ad = (row[ix['hedef']] or '').strip()
        F = (row[ix['F']] or '').strip().upper()
        R = (row[ix['R']] or '').strip().upper()
        if not ad or not F or not R:
            continue
        d = {k: (row[i].strip() if i is not None and i < len(row) else '')
             for k, i in ix.items()}
        d['hedef'] = ad
        d['F'], d['R'] = F, R
        d['urun_bp'] = int(_sayi(d['urun']) or 0)
        d['ayrim_sayi'] = _sayi(d['ayrim'])
        out.append(d)
    return out, yol


def sorun_gerekceleri(d):
    g = []
    # Ayrim esigi kontrolu, panelin "Ayrim" sutununda KAPSAM metni (ornegin
    # "13/13 kutu") gordugunde sayisal esikle karsilastirma YAPMAZ; onun yerine
    # kapsanan kutu / toplam kutu oranina bakar. Sebep evrensel ve genis
    # hedeflerdedir: oralarda rakip kumesi bosa yaklasir, ayrim katinin paydasi
    # sifira gider ve oran tanimsizlasir. Tanimsiz bir orani 10x esigiyle
    # karsilastirmak bir olcum degil, gurultu uretir. Bu yuzden o satirlarda
    # olcu kapsamadir. Esik degismemistir, olculen buyukluk degismistir.
    if 'IHLAL' in (d['geo'] or '').upper():
        g.append(('G', 'geometri ihlali: ' + d['geo'][:110]))
    dur = (d['durum'] or '').upper()
    sip = (d['siparis'] or '').upper()
    if 'ON KARAR' in dur or 'KOSULLU' in dur or 'KOSULLU' in sip:
        g.append(('K', 'kosullu / on karar: ' + (d['durum'] or d['siparis'])[:90]))
    if 'KAYITLI UYARI' in sip:
        g.append(('K', 'sipariste kayitli uyari: ' + d['siparis'][:90]))
    if sip.startswith('HAYIR'):
        g.append(('C', 'panelden cikarildi: ' + d['siparis'][:90]))
    a = d['ayrim_sayi']
    if a is not None and 'kapsam' not in (d['ayrim'] or '').lower() and 'kutu' not in (d['ayrim'] or '').lower():
        if a < AYRIM_ESIK:
            g.append(('A', 'ayrim %.1fx < %.0fx esigi' % (a, AYRIM_ESIK)))
    elif 'kutu' in (d['ayrim'] or '').lower():
        m = re.match(r'\s*(\d+)\s*/\s*(\d+)\s*kutu', d['ayrim'])
        if m and int(m.group(1)) < int(m.group(2)):
            g.append(('A', 'kapsam eksik: ' + d['ayrim'][:60]))
    bp = d['urun_bp']
    if bp > C.URUN_ONERILMEZ:
        g.append(('U', 'urun %d bp > %d (QuantiNova icin onerilmez)' % (bp, C.URUN_ONERILMEZ)))
    elif bp > C.URUN_IDEAL[1]:
        g.append(('U', 'urun %d bp ideal 60-150 disinda (30 sn ann/ext gerekir)' % bp))
    if 'AYRILMAZ' in (d['jel'] or '').upper():
        g.append(('P', 'plaka ici jel cakismasi: ' + d['jel'][:80]))
    # panel kapsam metninde eksik kutu
    m = re.match(r'\s*(\d+)\s*/\s*(\d+)\s*kutu', d.get('uye', '') or '')
    if m and int(m.group(1)) < int(m.group(2)):
        g.append(('A', 'uye kapsami %s' % d['uye'][:40]))
    return g


def sorunlu_hedefler():
    panel, yol = panel_oku()
    out = []
    for d in panel:
        g = sorun_gerekceleri(d)
        if g:
            d = dict(d)
            d['gerekceler'] = g
            d['etiketler'] = ''.join(sorted({x[0] for x in g}))
            out.append(d)
    return out, panel, yol


# ---------------------------------------------------------------- uyelik
UYELIK_TSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hedef_uyelik.tsv')


def acik_uyelik():
    """screening/hedef_uyelik.tsv - hedef adi -> uyelik. Elle duzenlenebilir."""
    m = {}
    if not os.path.exists(UYELIK_TSV):
        return m
    with open(UYELIK_TSV, encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            p = line.rstrip('\n').split('\t')
            if p[0] == 'hedef':
                continue
            if len(p) < 2:
                continue
            m[p[0].strip()] = dict(
                uye=[x.strip() for x in p[1].split(',') if x.strip()],
                haric=[x.strip() for x in (p[2] if len(p) > 2 else '').split(',') if x.strip()],
                kaynak=(p[3].strip() if len(p) > 3 else ''),
                not_=(p[4].strip() if len(p) > 4 else ''),
                duzey='')
    return m


def uyelik_oku():
    """hedefler.tsv -> hedef adi -> uye taxid listesi. '*A'/'*B'/'*F' = sinif tamami."""
    m = {}
    if not os.path.exists(C.HEDEFLER_TSV):
        return m
    with open(C.HEDEFLER_TSV, encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            p = line.rstrip('\n').split('\t')
            if p[0] == 'karar':
                continue
            if len(p) < 5:
                continue
            m[p[1].strip()] = dict(
                duzey=p[2].strip(),
                uye=[x for x in p[3].strip().split(',') if x],
                haric=[x for x in p[4].strip().split(',') if x],
                not_=p[5].strip() if len(p) > 5 else '',
            )
    return m


# panel adi -> hedefler.tsv adi (ad degisiklikleri; kaynak: KAPANIS_2026-08-02.md)
AD_ESLEME = {
    'Proteolitik_Synergistaceae': 'Proteolitik_Cloacibacillus',
    'Mantar_universal (F1)': 'Mantar_universal',
    'Mantar_universal (F2)': 'Mantar_universal',
    'Metanomikrobiyales_hidrojenotrof': 'Hidrojenotrofik_metanojenler',
    'Methanosarcina mazei / M. soligelidi grubu': 'Methanosarcina_mazei_turu',
    'Methanothrix_cinsi': 'Methanothrix_soehngenii_turu',
    'Methanosarcina_cinsi': 'Asetoklastik_metanojenler',
    'Metanojen_universal': 'Arke_universal',
    'Asetoklastik_metanojenler': 'Asetoklastik_metanojenler',
    'Bacteroidales_kumesi': 'Sakarolitik_bakteriler',
    'Sakarolitik_Sphaerochaeta': 'Sakarolitik_Sphaerochaeta',
    'Proteolitik_Cloacimonas': 'Proteolitik_sintrofik_bakteriler',
    'Petriella_musispora': 'Petriella_musispora',
    'Microascaceae_askomikot': 'Microascaceae_askomikot',
}


def taxid_adlari():
    m = {}
    if os.path.exists(C.TAXID_ADLARI):
        with open(C.TAXID_ADLARI, encoding='utf-8') as fh:
            for line in fh:
                p = line.rstrip('\n').split('\t')
                if len(p) >= 2:
                    m[p[0].strip()] = p[1].strip()
    return m


# ---------------------------------------------------------------- kutular
def kutular():
    """(kutu_adi, sinif, taxid, fastq_yolu) listesi. kutu_adi ornek: 'B-1_1129264'."""
    out = []
    for p in sorted(glob.glob(os.path.join(C.FASTQ, '*', '*.fastq*'))):
        d = os.path.basename(os.path.dirname(p))
        b = os.path.basename(p)
        m = re.search(r'reads_([0-9]+)\.fastq', b)
        if not m:
            continue
        tax = m.group(1)
        sinif = d.split('-')[0]          # A1 / A2 / B / F1 / F2
        out.append(dict(kutu='%s_%s' % (d, tax), grup=d, sinif=sinif, taxid=tax, yol=p))
    return out


def konsensusler():
    """KANONIK konsensusler. Kaynak: konsensus_kanonik/INDEKS.tsv (hepsi SENSE).
    2026-08-02 yon normalizasyonu: eskiden 'consensus sequences' klasoru dogrudan
    okunuyordu; o klasor KARISIK yonludur (71 antisense / 27 sense) ve ters yonlu
    bir konsensuste in-silico PCR SESSIZCE 0 urun verir. Artik tek kanonik kaynak
    okunur; glob DEGIL INDEKS kullanilir (bagli klasorde silinemeyen kalinti
    dosyalar var). Indeks yoksa hata verilir - sessizce eski klasore DUSULMEZ."""
    import csv as _csv
    ix = getattr(C, 'KONSENSUS_INDEKS', None)
    if not ix or not os.path.exists(ix):
        raise RuntimeError(
            'Kanonik konsensus indeksi yok: %s\n'
            'Once uretin:  python screening/kanonik_uret.py --kok .\n'
            'Karisik yonlu "consensus sequences" klasorune DUSULMEZ - yon hatasi '
            'sessizce 0 urun verir.' % ix)
    kok = os.path.dirname(ix)
    out = []
    for r in _csv.DictReader(open(ix, encoding='utf-8'), delimiter='\t'):
        yol = os.path.join(kok, r['dosya'])
        if not os.path.exists(yol):
            continue
        seq = ''.join(l.strip() for l in open(yol, encoding='utf-8', errors='ignore')
                      if not l.startswith('>'))
        kutu = r['kutu']
        out.append(dict(kutu=kutu, grup=kutu.split('_')[0],
                        sinif=r['sinif'], taxid=kutu.split('_')[-1],
                        yol=yol, yon='SENSE', kaynak=r['kaynak'],
                        dizi=motor.clean(seq.upper())))
    return out


def hedef_baglami(panel_satiri, uyelik=None, kons=None, kut=None):
    """Bir panel hedefi icin: omurga konsensus, uye kutular, rakip kutular."""
    # UYELIK KAYNAK ONCELIGI: once screening/hedef_uyelik.tsv (elle
    # duzenlenen acik tanim), yoksa steps/hedefler.tsv (projenin karar
    # tablosu). Hangi kaynagin kullanildigi 'uyelik_kaynagi' alaninda tasinir ve
    # her rapora yazilir - sessiz kaynak degisikligi olmaz. Tanim hic yoksa
    # 'YOK' yazilir ve hedef atlanir; varsayilan bir uyelik UYDURULMAZ, cunku
    # yanlis bir uyelik tanimi ayrim katini oldugundan kucuk ya da buyuk
    # gosterir (olculen ornek: ayni cift 0,0x ile 23,5x arasinda oynuyor).
    #
    # Omurga, uye konsensuslerin EN UZUNU secilir: aday primer pencereleri
    # omurga uzerinde uretildigi icin kisa bir omurga arama uzayini kirpar.
    # Konsensusler kanonik klasorden geldigi icin hepsi SENSE yondedir; ters
    # yonlu bir omurgada uretilen adaylar okumalarda hicbir urun vermezdi.
    uyelik = uyelik if uyelik is not None else uyelik_oku()
    kons = kons if kons is not None else konsensusler()
    kut = kut if kut is not None else kutular()
    ad = panel_satiri['hedef']
    acik = acik_uyelik()
    if ad in acik:
        u = acik[ad]
        anahtar = 'screening/hedef_uyelik.tsv (%s)' % (u.get('kaynak') or '?')
    else:
        anahtar = AD_ESLEME.get(ad, ad)
        u = uyelik.get(anahtar)
        anahtar = 'steps/hedefler.tsv:%s' % anahtar
    siniflar = [s.strip() for s in (panel_satiri['sinif'] or '').split('/') if s.strip()]
    if not siniflar:
        siniflar = ['A1', 'A2', 'B', 'F1', 'F2']

    if u is None:
        uye_tax, haric = [], []
    else:
        uye_tax, haric = list(u['uye']), list(u['haric'])

    # KUTU: oneki - uyelik KUTU duzeyinde yazilabilir.
    #
    # NEDEN (2026-08-11): taxid, uyeligi ANLATMAYA YETMIYOR. Ornek: F2-1_500148,
    # F2-2_500148 ve F2-4_500148 ayni Kraken taxid'ini tasir; olculen kimlikleri
    # ise Petriella setifera, Petriella setifera ve "adlandirilamayan soy - en
    # yakin kayit Lomentospora". Taxid ile yazilan bir uyelik ucunu de ayni
    # torbaya atar; oysa ucuncusu hedefin uyesi DEGIL, rakibidir - ve olcumde
    # cifti dusuren tam da o kutudur (1655/3000 okuma). Bu yuzden uyelik artik
    # "KUTU:F2-1_500148" biciminde, OLCULEN kimlige gore, kutu kutu yazilabilir.
    # Ayni kural haric sutununda da gecerlidir.
    kutu_uye = set(t[5:] for t in uye_tax if t.startswith('KUTU:'))
    kutu_haric = set(t[5:] for t in haric if t.startswith('KUTU:'))
    uye_tax = [t for t in uye_tax if not t.startswith('KUTU:')]
    haric = [t for t in haric if not t.startswith('KUTU:')]

    yildiz = [t for t in uye_tax if t.startswith('*')]
    if yildiz:
        onek = [t[1:] for t in yildiz]
        uye_kut = [k for k in kut if any(k['sinif'] == o or k['sinif'].startswith(o)
                                         for o in onek) and k['taxid'] not in haric]
        uye_tax = sorted({k['taxid'] for k in uye_kut})

    def _uye_mi(k):
        if k['kutu'] in kutu_haric:
            return False
        if k['kutu'] in kutu_uye:
            return True
        return (not kutu_uye) and k['taxid'] in uye_tax

    def _rakip_mi(k):
        if k['kutu'] in kutu_haric or _uye_mi(k):
            return False
        return k['taxid'] not in haric

    uye_kut = [k for k in kut if _uye_mi(k) and k['sinif'] in siniflar]
    if not uye_kut:
        uye_kut = [k for k in kut if _uye_mi(k)]
    rakip_kut = [k for k in kut if k['sinif'] in siniflar and _rakip_mi(k)]

    uye_kons = [k for k in kons if _uye_mi(k) and k['sinif'] in siniflar]
    if not uye_kons:
        uye_kons = [k for k in kons if _uye_mi(k)]
    rakip_kons = [k for k in kons if k['sinif'] in siniflar and _rakip_mi(k)]
    if kutu_uye:
        uye_tax = sorted({k['taxid'] for k in uye_kut} |
                         {k['taxid'] for k in uye_kons})

    omurga = max(uye_kons, key=lambda k: len(k['dizi'])) if uye_kons else None
    return dict(hedef=ad, anahtar=anahtar, siniflar=siniflar,
                uye_tax=uye_tax, haric=haric,
                uye_kutu=uye_kut, rakip_kutu=rakip_kut,
                uye_kons=uye_kons, rakip_kons=rakip_kons,
                omurga=omurga,
                uyelik_kaynagi=anahtar if u else 'YOK',
                uyelik_notu=(u.get('not_', '') if u else ''))
