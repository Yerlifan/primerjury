# -*- coding: utf-8 -*-
"""LITERATUR KONTROLU - onerilen ad gecerli mi, es anlamlisi var mi, cinsi
yakin zamanda revize edilmis mi?

NEDEN VAR
---------
Parascedosporium olayi gosterdi ki veritabanindaki ad ile GUNCEL KABUL EDILEN ad
ayni olmayabilir ve cins sinirlari revizyon altinda olabilir. Bu, sonucu degistirir.
Bu yuzden literatur kontrolu artik ZORUNLU BIR ADIMDIR, soru geldikce yapilan bir
sey degil.

UC KATMAN - ve sinirlari acikca yazili
  1) NCBI Taxonomy (OTOMATIK, E-utilities): erisim numarasi -> taxid -> guncel ad,
     es anlamlilar, rutbe, soy zinciri. HIZLI ve kuyruksuz. AMA NCBI Taxonomy bir
     NOMENKLATUR OTORITESI DEGILDIR - pratik bir siniflandirmadir.
  2) NOMENKLATUR OTORITESI (mantar: MycoBank / Index Fungorum; bakteri-arke: LPSN).
     Aga cikilabiliyorsa sorgulanir; cikilamiyorsa ELLE KONTROL LISTESI uretilir:
     hangi ad, hangi otoritede, dogrudan sorgu baglantisiyla.
  3) REVIZYON UYARISI (PubMed, hafif arama): cins adi + "comb. nov." / "gen. nov."
     / "revision". GURULTULU olabilir - KARAR VERDIRMEZ, yalniz satiri isaretler.

AG YOKSA adim ATLANMAZ: "literatur kontrolu yapilamadi" diye isaretlenir ve elle
kontrol listesi yine uretilir.
"""

# -------------------------------------------------------------------------
# literatur_kontrol.py — bir isabetin adinin GUNCEL kabul edilen ad olup
# olmadigini, es anlamlilarini ve cinsin yakin zamanda revize edilip
# edilmedigini uc katmanda sorgular.
#
# GİRDİ  : cagiran betikten gelen referans kaydi basligi (icinden erisim
#          numarasi cikarilir) ve onerilen ad; ag uzerinden NCBI E-utilities
#          (esummary / efetch / esearch) ve PubMed. Yerel dosya okumaz.
# ÇIKTI  : cagirana dondurulen dict (guncel ad, es anlamlilar, rutbe, soy,
#          revizyon uyarisi, otorite baglantilari) ve istege bagli olarak
#          <cikti>/LITERATUR_ELLE_KONTROL.tsv dosyasi.
# ÇAĞRAN : screening.bat -> I ve G tuslari (dolayli: kimlik_dogrulama.py
#          ve tum_kutu_kimlikleri.py bu modulu ZORUNLU adim olarak yukler).
#
# SINIR: NCBI Taxonomy pratik bir siniflandirmadir, NOMENKLATUR OTORITESI
# DEGILDIR. Bu yuzden otomatik katman hicbir zaman tek basina karar vermez;
# mantar icin MycoBank/Index Fungorum, bakteri-arke icin LPSN baglantisi her
# satirda uretilir ve "otorite kontrolu GEREKLI" damgasi dusurulmez.
# -------------------------------------------------------------------------
import os, re, json, time, urllib.parse

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
ARAC = 'PrimerJury'
POSTA = ''          # NCBI nezaket alani; kullanici doldurabilir

OTORITE = {
    'mantar': [('MycoBank', 'https://www.mycobank.org/page/Name%%20details%%20page/?Rec=%s'),
               ('Index Fungorum', 'http://www.indexfungorum.org/names/Names.asp?strGenus=%s')],
    'bakteri_arke': [('LPSN', 'https://lpsn.dsmz.de/search?word=%s')],
}
MANTAR_IPUCU = ('fungi', 'ascomyc', 'basidiomyc', 'its', '28s', '18s',
                'petriella', 'parascedosporium', 'trichoderma', 'metarhizium')


def _al(url, zaman=25):
    import urllib.request
    with urllib.request.urlopen(url, timeout=zaman) as f:
        return f.read().decode('utf-8', 'replace')


def erisim_no(baslik):
    """Referans basligindan erisim numarasini cikar (NG_064074.1, NR_1234.1, AY882356)."""
    m = re.search(r'\b([A-Z]{1,2}_?[0-9]{5,9}(?:\.[0-9]+)?)\b', baslik or '')
    return m.group(1) if m else None


# Erisim numarasi -> taxid. esummary yanitindaki "uids" anahtari kayit degil
# indeks oldugu icin atlanir; ilk taxid tasiyan kayit dondurulur.
def taxid_al(acc):
    u = (EUTILS + 'esummary.fcgi?db=nuccore&id=%s&retmode=json&tool=%s'
         % (urllib.parse.quote(acc), ARAC)) + (('&email=%s' % POSTA) if POSTA else '')
    d = json.loads(_al(u))
    r = d.get('result', {})
    for k, v in r.items():
        if k == 'uids':
            continue
        t = v.get('taxid')
        if t:
            return str(t)
    return None


# taxid -> guncel bilimsel ad, rutbe, soy zinciri ve es anlamlilar. Es anlamlilar
# hem <Synonym> hem <EquivalentName> etiketlerinden toplanir: veritabaninda gecen
# ad bunlardan biriyse "ad farkli" uyarisi gercek bir celiski degil sadece eskimis
# bir kullanimdir; ikisini ayirt edebilmek icin ikisi de raporlanir.
def taxonomy_cek(taxid):
    u = (EUTILS + 'efetch.fcgi?db=taxonomy&id=%s&retmode=xml&tool=%s' % (taxid, ARAC))
    x = _al(u)
    def bul(etiket):
        m = re.search(r'<%s>(.*?)</%s>' % (etiket, etiket), x, re.S)
        return (m.group(1).strip() if m else '')
    esan = re.findall(r'<Synonym>(.*?)</Synonym>', x, re.S)
    esan += re.findall(r'<EquivalentName>(.*?)</EquivalentName>', x, re.S)
    return dict(guncel_ad=bul('ScientificName'), rutbe=bul('Rank'),
                soy=bul('Lineage'), es_anlamlilar='; '.join(dict.fromkeys(esan)))


def revizyon_ara(cins, yil=8):
    """PubMed'de cins + nomenklatur terimleri. KARAR VERDIRMEZ, uyari uretir."""
    if not cins:
        return dict(sayi=0, pmid=[], not_='cins adi yok')
    terim = ('%s[Title/Abstract] AND ("comb. nov."[All Fields] OR "gen. nov."[All Fields] '
             'OR revision[Title] OR taxonomy[Title])' % cins)
    u = (EUTILS + 'esearch.fcgi?db=pubmed&term=%s&retmax=5&retmode=json&datetype=pdat'
         '&reldate=%d&tool=%s' % (urllib.parse.quote(terim), yil * 365, ARAC))
    d = json.loads(_al(u))
    r = d.get('esearchresult', {})
    pm = r.get('idlist', []) or []
    return dict(sayi=int(r.get('count', 0) or 0), pmid=pm,
                not_=(u'son %d yilda %s kayit' % (yil, r.get('count', '0'))))


# Hangi nomenklatur otoritesine gidilecegini belirler (mantar mi, bakteri-arke mi).
# Bu bir KIMLIK karari degildir, yalnizca hangi baglantinin uretilecegini secer;
# yanlis tahmin sonucu bozmaz, kullaniciya fazladan bir baglanti gosterir.
def alan_tahmin(baslik, lokus=''):
    s = ((baslik or '') + ' ' + (lokus or '')).lower()
    return 'mantar' if any(k in s for k in MANTAR_IPUCU) else 'bakteri_arke'


# Index Fungorum ve LPSN cins bazli sorgu alir, MycoBank tam ad alir; bu yuzden
# hedef dizge kalibin kendisine gore secilir.
def otorite_baglantilari(ad, alan):
    out = []
    if not ad or ad == '-':
        return out
    cins = ad.split()[0]
    for otr, kalip in OTORITE.get(alan, []):
        hedef = cins if 'indexfungorum' in kalip or 'lpsn' in kalip else ad
        try:
            out.append((otr, kalip % urllib.parse.quote(hedef)))
        except TypeError:
            out.append((otr, kalip))
    return out


def kontrol_et(baslik, onerilen_ad, lokus='', ag=True):
    """Tek bir isabet icin uc katmanli literatur kontrolu."""
    acc = erisim_no(baslik)
    alan = alan_tahmin(baslik, lokus)
    sonuc = dict(erisim_no=acc or '-', alan=alan, vtb_adi='-', ncbi_guncel_ad='-',
                 es_anlamlilar='-', rutbe='-', soy='-', ad_farkli_mi='-',
                 revizyon_uyarisi='-', revizyon_pmid='-',
                 otorite_kontrolu='GEREKLI', otorite_baglantilari='',
                 durum='YAPILAMADI')
    # veritabanindaki ad (baslikta gecen ikili ad)
    m = re.search(r'\b([A-Z][a-z]{2,})\s+([a-z]{2,})\b', baslik or '')
    if m:
        sonuc['vtb_adi'] = '%s %s' % (m.group(1), m.group(2))
    bag = otorite_baglantilari(onerilen_ad if onerilen_ad and onerilen_ad != '-'
                               else sonuc['vtb_adi'], alan)
    sonuc['otorite_baglantilari'] = ' | '.join('%s: %s' % (o, u) for o, u in bag)
    if not ag or not acc:
        sonuc['durum'] = ('AG YOK - literatur kontrolu YAPILAMADI' if not ag
                          else 'erisim numarasi cikarilamadi')
        return sonuc
    try:
        t = taxid_al(acc)
        if not t:
            sonuc['durum'] = 'taxid bulunamadi'
            return sonuc
        time.sleep(0.34)                       # NCBI hiz siniri (3/sn)
        tx = taxonomy_cek(t)
        time.sleep(0.34)
        sonuc.update(ncbi_guncel_ad=tx['guncel_ad'] or '-', rutbe=tx['rutbe'] or '-',
                     es_anlamlilar=tx['es_anlamlilar'] or '-', soy=(tx['soy'] or '-')[:200])
        if sonuc['vtb_adi'] != '-' and tx['guncel_ad'] and \
                sonuc['vtb_adi'].lower() != tx['guncel_ad'].lower():
            sonuc['ad_farkli_mi'] = ('EVET - veritabani "%s" diyor, NCBI guncel adi "%s"'
                                     % (sonuc['vtb_adi'], tx['guncel_ad']))
        else:
            sonuc['ad_farkli_mi'] = 'hayir'
        cins = (tx['guncel_ad'] or sonuc['vtb_adi']).split()[0] if (
            tx['guncel_ad'] or sonuc['vtb_adi'] != '-') else None
        rv = revizyon_ara(cins)
        time.sleep(0.34)
        if rv['sayi']:
            sonuc['revizyon_uyarisi'] = (u'UYARI - %s cinsi icin %s. Karar verdirmez, '
                                         u'ELLE bakilmali.' % (cins, rv['not_']))
            sonuc['revizyon_pmid'] = ','.join(rv['pmid'])
        else:
            sonuc['revizyon_uyarisi'] = 'yok'
        sonuc['durum'] = 'TAMAM (NCBI Taxonomy)'
        sonuc['otorite_kontrolu'] = (
            u'GEREKLI - NCBI Taxonomy nomenklatur otoritesi DEGILDIR; '
            u'%s icin %s' % (alan, ', '.join(o for o, _ in bag) or 'otorite yok'))
    except Exception as e:
        sonuc['durum'] = u'AG/SORGU HATASI: %s' % type(e).__name__
    return sonuc


def elle_liste_yaz(cikti, satirlar):
    """Aga cikilamayan ya da otorite gerektiren adlar icin hazir kontrol listesi."""
    yol = os.path.join(cikti, 'LITERATUR_ELLE_KONTROL.tsv')
    with open(yol, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# ELLE LITERATUR KONTROLU - her ad icin otoritede dogrulanmali.\n')
        fh.write(u'# NCBI Taxonomy pratik bir siniflandirmadir, NOMENKLATUR OTORITESI DEGILDIR.\n')
        fh.write(u'# Mantar: MycoBank / Index Fungorum  |  Bakteri-Arke: LPSN\n')
        import csv as _c
        w = _c.writer(fh, delimiter='\t')
        w.writerow(['iddia_no', 'onerilen_ad', 'veritabani_adi', 'ncbi_guncel_ad',
                    'alan', 'otoritede_kontrol_et', 'sorgu_baglantilari', 'durum'])
        for s in satirlar:
            w.writerow([s.get('no', ''), s.get('onerilen_ad', '-'), s.get('vtb_adi', '-'),
                        s.get('ncbi_guncel_ad', '-'), s.get('alan', '-'),
                        s.get('otorite_kontrolu', 'GEREKLI'),
                        s.get('otorite_baglantilari', ''), s.get('durum', '-')])
    return yol
