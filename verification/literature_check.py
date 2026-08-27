# -*- coding: utf-8 -*-
"""THE LITERATURE CHECK: is the suggested name valid, does it have a synonym, has
its genus been revised recently?

WHY IT EXISTS
-------------
The Parascedosporium case showed that the name in the database and the CURRENT
ACCEPTED name may not be the same, and that genus boundaries can be under revision.
That changes the result. So the literature check is now A REQUIRED STEP rather than
something done when the question comes up.

THREE LAYERS, and their limits written out plainly
  1) NCBI Taxonomy (AUTOMATIC, E-utilities): accession -> taxid -> the current name,
     synonyms, rank, lineage. FAST and with no queue. BUT NCBI Taxonomy IS NOT A
     NOMENCLATURAL AUTHORITY, it is a practical classification.
  2) THE NOMENCLATURAL AUTHORITY (fungi: MycoBank / Index Fungorum; bacteria and
     archaea: LPSN). It is queried if the network is reachable; if it is not, A
     MANUAL CHECK LIST is produced: which name, at which authority, with a direct
     query link.
  3) THE REVISION WARNING (PubMed, a light search): the genus name plus "comb. nov."
     / "gen. nov." / "revision". It can be NOISY, so IT DOES NOT DECIDE, it only
     marks the row.

WITH NO NETWORK the step IS NOT SKIPPED: it is marked "the literature check could
not be made" and the manual check list is produced all the same.

"""

# -----------------------------------------------------------------------
# literature_check.py asks, in three layers, whether a hit's name is the CURRENT
# accepted name, what its synonyms are, and whether the genus has been revised
# recently.
#
# INPUT  : the reference record header coming from the calling script (the accession
#          number is taken out of it) and the suggested name; over the network, the
#          NCBI E-utilities (esummary / efetch / esearch) and PubMed. It reads no
#          local file.
# OUTPUT : the dict returned to the caller (the current name, synonyms, rank,
#          lineage, a revision warning, authority links) and, optionally, the file
#          <output>/LITERATUR_ELLE_KONTROL.tsv.
# CALLED BY: verification/full_chain.py -> keys I and G (indirectly:
#          identity_verification.py and all_bin_identities.py load this module as a
#          REQUIRED step).
#
# THE LIMIT: NCBI Taxonomy is a practical classification, NOT A NOMENCLATURAL
# AUTHORITY. That is why the automatic layer never decides on its own; a MycoBank
# or Index Fungorum link for fungi and an LPSN link for bacteria and archaea is
# produced on every row, and the "an authority check IS NEEDED" stamp is not
# dropped.
# -----------------------------------------------------------------------
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
    'Pull the accession number out of a reference header (NG_064074.1, NR_1234.1, AY882356).'
    m = re.search(r'\b([A-Z]{1,2}_?[0-9]{5,9}(?:\.[0-9]+)?)\b', baslik or '')
    return m.group(1) if m else None


# An accession number -> taxid. The "uids" key in the esummary reply is an index
# rather than a record, so it is skipped; the first record carrying a taxid is
# returned.
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


# taxid -> the current scientific name, rank, lineage and synonyms. The synonyms
# are collected from both the <Synonym> and the <EquivalentName> tags: if the name
# in the database is one of those, the "the name differs" warning is not a real
# contradiction but merely an outdated usage; both are reported so the two can be
# told apart.
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
    """The genus plus nomenclature terms in PubMed. IT DOES NOT DECIDE, it produces a warning."""
    if not cins:
        return dict(sayi=0, pmid=[], not_='there is no genus name')
    terim = ('%s[Title/Abstract] AND ("comb. nov."[All Fields] OR "gen. nov."[All Fields] '
             'OR revision[Title] OR taxonomy[Title])' % cins)
    u = (EUTILS + 'esearch.fcgi?db=pubmed&term=%s&retmax=5&retmode=json&datetype=pdat'
         '&reldate=%d&tool=%s' % (urllib.parse.quote(terim), yil * 365, ARAC))
    d = json.loads(_al(u))
    r = d.get('esearchresult', {})
    pm = r.get('idlist', []) or []
    return dict(sayi=int(r.get('count', 0) or 0), pmid=pm,
                not_=('%s records in the last %d years' % (yil, r.get('count', '0'))))


# It decides which nomenclatural authority to go to (fungi, or bacteria and
# archaea). This IS NOT an identity decision, it only chooses which link is
# produced; a wrong guess does not spoil the result, it shows the user one link too
# many.
def alan_tahmin(baslik, lokus=''):
    s = ((baslik or '') + ' ' + (lokus or '')).lower()
    return 'mantar' if any(k in s for k in MANTAR_IPUCU) else 'bakteri_arke'


# Index Fungorum and LPSN take a genus level query while MycoBank takes the full
# name; that is why the target string is chosen according to the template itself.
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
    """The three layer literature check for a single hit."""
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
        sonuc['durum'] = ('THERE IS NO NETWORK, so the literature check COULD NOT BE MADE' if not ag
                          else 'erisim numarasi cikarilamadi')
        return sonuc
    try:
        t = taxid_al(acc)
        if not t:
            sonuc['durum'] = 'no taxid found'
            return sonuc
        time.sleep(0.34)                       # NCBI hiz siniri (3/sn)
        tx = taxonomy_cek(t)
        time.sleep(0.34)
        sonuc.update(ncbi_guncel_ad=tx['guncel_ad'] or '-', rutbe=tx['rutbe'] or '-',
                     es_anlamlilar=tx['es_anlamlilar'] or '-', soy=(tx['soy'] or '-')[:200])
        if sonuc['vtb_adi'] != '-' and tx['guncel_ad'] and \
                sonuc['vtb_adi'].lower() != tx['guncel_ad'].lower():
            sonuc['ad_farkli_mi'] = ('YES: the database says "%s" while the current NCBI name is "%s"'
                                     % (sonuc['vtb_adi'], tx['guncel_ad']))
        else:
            sonuc['ad_farkli_mi'] = 'hayir'
        cins = (tx['guncel_ad'] or sonuc['vtb_adi']).split()[0] if (
            tx['guncel_ad'] or sonuc['vtb_adi'] != '-') else None
        rv = revizyon_ara(cins)
        time.sleep(0.34)
        if rv['sayi']:
            sonuc['revizyon_uyarisi'] = ('A WARNING: for the genus %s, %s. It '
                                         'decides nothing and has to be '
                                         'looked at BY HAND.' % (cins, rv['not_']))
            sonuc['revizyon_pmid'] = ','.join(rv['pmid'])
        else:
            sonuc['revizyon_uyarisi'] = 'yok'
        sonuc['durum'] = 'TAMAM (NCBI Taxonomy)'
        sonuc['otorite_kontrolu'] = (
            'REQUIRED: NCBI Taxonomy IS NOT a nomenclatural authority; for '
            '%s, %s' % (alan, ', '.join(o for o, _ in bag) or 'otorite yok'))
    except Exception as e:
        sonuc['durum'] = 'A NETWORK OR QUERY ERROR: %s' % type(e).__name__
    return sonuc


def elle_liste_yaz(cikti, satirlar):
    """A ready check list for the names the network could not reach, or that need an authority."""
    yol = os.path.join(cikti, 'LITERATUR_ELLE_KONTROL.tsv')
    with open(yol, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# MANUAL LITERATURE CHECK - every name must be verified against an authority.\n')
        fh.write(u'# NCBI Taxonomy is a practical classification, NOT A NOMENCLATURAL AUTHORITY.\n')
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
