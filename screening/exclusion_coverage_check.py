# -*- coding: utf-8 -*-
"""DISLAMA HARITASI KAPSAMA DENETIMI  -  hedef_taxid.tsv dogru mu?

SORU
----
hedef_taxid.tsv her hedef icin "NCBI'da dislanacak takson"u soyluyor. Bu
takson, hedefin BUTUN uyelerini icermek ZORUNDADIR. Icermezse disarida kalan
uye hedefin KENDISI oldugu halde "hedef disi urun" diye sayilir ve cift haksiz
yere kirli gorunur.

Bu denetim tahmine yer birakmaz: her uye taksid'in NCBI soy zincirini ceker ve
dislanan taksid o zincirde GERCEKTEN var mi diye bakar.

NEDEN YAZILDI (2026-08-10, canli olcumde yakalandi)
---------------------------------------------------
Asetoklastik_metanojenler icin dislama 94695 (Methanosarcinales) yazilmisti;
gerekce "Methanothrix de bu takimin altinda" idi. NCBI Methanothrix'i
Methanosarcinales'ten CIKARIP kendi takimina (Methanotrichales, 2905377)
tasimis. Yani gerekce artik dogru degildi ve Methanothrix uyeleri hedef disi
sayilacakti. Ezberden yazilan taksonomi eskiyor; bu yuzden zincir NCBI'dan
CEKILIR, hatirlanmaz.

Kosum:
    python screening/exclusion_coverage_check.py --kok .
Cikis kodu: 0 = butun hedefler kapsandi, 1 = kapsanmayan uye var.
"""
from __future__ import print_function

import argparse
import io
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'


def oku_harita(yol):
    h = {}
    for l in io.open(yol, encoding='utf-8'):
        l = l.rstrip('\n')
        if not l.strip() or l.startswith('#') or l.startswith('hedef\t'):
            continue
        p = l.split('\t')
        if len(p) >= 2 and p[1].strip():
            h[p[0].strip()] = [x.strip() for x in p[1].split(',') if x.strip()]
    return h


def oku_uyelik(yol):
    u = {}
    for l in io.open(yol, encoding='utf-8'):
        l = l.rstrip('\n')
        if not l.strip() or l.startswith('#') or l.startswith('hedef\t'):
            continue
        p = l.split('\t')
        if len(p) >= 2 and p[1].strip():
            u[p[0].strip()] = [x.strip() for x in p[1].split(',') if x.strip()]
    return u


def oku_uye_kutular(yol):
    """hedef -> uye KUTU adlari (A1-1_2209 gibi). Kraken taxid'i degil KUTU."""
    u = {}
    if not os.path.exists(yol):
        return u
    bas = None
    for l in io.open(yol, encoding='utf-8'):
        l = l.rstrip('\n')
        if not l.strip() or l.startswith('#'):
            continue
        p = l.split('\t')
        if bas is None:
            bas = p
            continue
        r = dict(zip(bas, p))
        k = (r.get('yeni_uye_kutular') or r.get('eski_uye_kutular') or '').strip()
        if k:
            u[r[bas[0]].strip()] = [x.strip() for x in k.split(';') if x.strip()]
    return u


def oku_olculen_kimlik(yol):
    """kutu -> OLCULEN kimlik metni. Kraken etiketi DEGIL, bizim olcumumuz."""
    o = {}
    if not os.path.exists(yol):
        return o
    for l in io.open(yol, encoding='utf-8'):
        l = l.rstrip('\n')
        if not l.strip() or l.startswith('#'):
            continue
        p = l.split('\t')
        if len(p) >= 4 and p[0].strip():
            o[p[0].strip()] = p[3].strip()
    return o


_ATLA_AD = ('sp.', 'cf.', 'strain', 'isolate', 'group', 'incertae', 'sedis',
            'uncultured', 'unclassified', 'environmental', 'bacterium',
            'archaeon', 'fungal', 'sp')


def ad_adaylari(kimlik):
    """Olculen kimlik metninden NCBI'da aranabilecek ad adaylari - genisten dara.

    Iki bicim var:
      "Petriella setifera"                       -> dogrudan ad
      "adlandirilamayan soy - EN YAKIN KAYIT: ...Fungi;Dikarya;Ascomycota;..."
                                                 -> soy zincirinden en DAR ad
    """
    k = (kimlik or '').strip()
    if not k:
        return []
    if 'EN YAKIN KAYIT' in k:
        gov = k.split('EN YAKIN KAYIT:', 1)[1]
        gov = re.sub(r'\(%[^)]*\)', ' ', gov)
        parca = [x.strip() for x in re.split(r'[;|]', gov) if x.strip()]
        parca = [x for x in parca
                 if re.match(r'^[A-Z][A-Za-z\- ]{2,}$', x)
                 and x.lower() not in ('eukaryota', 'bacteria', 'archaea')]
        return list(reversed(parca))[:6]          # en dardan basla
    ad = re.sub(r'\s*\(%[^)]*\)', '', k).strip()
    out = [ad]
    kel = ad.split()
    if len(kel) >= 2 and kel[1].lower() not in _ATLA_AD:
        out.append(' '.join(kel[:2]))
    if kel:
        out.append(kel[0])
    return [x for x in out if x and x.lower() not in _ATLA_AD]


_AD_ONBELLEK = {}


def ad_taxid(ad):
    if ad in _AD_ONBELLEK:
        return _AD_ONBELLEK[ad]
    q = urllib.parse.urlencode(dict(db='taxonomy',
                                    term='"%s"[Scientific Name]' % ad, retmode='json'))
    try:
        r = json.loads(urllib.request.urlopen(EUTILS + 'esearch.fcgi?' + q,
                                              timeout=40).read().decode())
        t = (r.get('esearchresult', {}).get('idlist') or [None])[0]
    except Exception:
        t = None
    time.sleep(0.35)
    _AD_ONBELLEK[ad] = t
    return t


def soy_zincirleri(taxidler, yaz):
    """taxid -> (ad, zincirdeki taxid kumesi). NCBI'dan cekilir, ezberlenmez.

    XML duzgun ayristirilir. Duz metin regex'i BU ISE YARAMAZ: <LineageEx>
    icinde ic ice <Taxon> ogeleri var ve regex ile bolunce her soy halkasi
    ayri bir kayit sanilir; kapsama sinavi o zaman dogru cifleri bile
    "kapsamiyor" der (2026-08-10'da tam olarak bu oldu).
    """
    import xml.etree.ElementTree as ET
    out = {}
    tl = sorted(set(taxidler))
    for i in range(0, len(tl), 80):
        oba = tl[i:i + 80]
        q = urllib.parse.urlencode(dict(db='taxonomy', id=','.join(oba), retmode='xml'))
        x = urllib.request.urlopen(EUTILS + 'efetch.fcgi?' + q,
                                   timeout=60).read().decode('utf-8', 'replace')
        kok = ET.fromstring(x)
        for tx in kok.findall('Taxon'):            # YALNIZ ust duzey kayitlar
            tid = tx.findtext('TaxId') or ''
            ad = tx.findtext('ScientificName') or ''
            zincir = {tid}
            le = tx.find('LineageEx')
            if le is not None:
                for h in le.findall('Taxon'):
                    if h.findtext('TaxId'):
                        zincir.add(h.findtext('TaxId'))
            # AKA/merged taxid: sorulan numara ile donen numara farkli olabilir
            for aka in tx.findall('AkaTaxIds/TaxId'):
                if aka.text:
                    out[aka.text] = (ad, zincir)
            out[tid] = (ad, zincir)
        time.sleep(0.4)
        yaz(u'  lineages pulled: %d/%d' % (min(i + 80, len(tl)), len(tl)))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', '--kok', dest='kok', default='.')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)

    def yaz(s=''):
        print(s, flush=True)

    hy = os.path.join(kok, 'screening', 'hedef_taxid.tsv')
    uy = os.path.join(kok, 'screening', 'hedef_uyelik.tsv')
    for f in (hy, uy):
        if not os.path.exists(f):
            sys.exit(u'ERROR: %s is missing.' % f)
    H = oku_harita(hy)
    U = oku_uyelik(uy)

    # KABUL EDILMIS ISTISNALAR: bilerek disarida birakilan uyeler. Alarmi
    # susturur ama sebebi dosyada durur - susturma sessizce olmaz.
    IST = {}
    iy = os.path.join(kok, 'screening', 'kapsama_istisna.tsv')
    if os.path.exists(iy):
        for l in io.open(iy, encoding='utf-8'):
            l = l.rstrip('\n')
            if not l.strip() or l.startswith('#'):
                continue
            p = l.split('\t')
            if len(p) >= 2:
                IST.setdefault(p[0].strip(), {})[p[1].strip()] = (
                    p[2].strip() if len(p) > 2 else '')

    yaz(u'=' * 78)
    yaz(u'  THE COVERAGE AUDIT OF THE EXCLUSION MAP   %s' % time.strftime('%Y-%m-%d %H:%M'))
    yaz(u'=' * 78)
    yaz(u'  exclusion map    : %d targets' % len(H))
    yaz(u'  membership table : %d targets' % len(U))

    # ad eslesmesi: uyelik tablosunun anahtari kisa olabilir (ornek:
    # "Proteolitik_Synergistaceae" <-> "Proteolitik_Synergistaceae (eski ad: ...)").
    def esle(ad):
        if ad in U:
            return ad
        for k in U:
            if ad.startswith(k) or k.startswith(ad):
                return k
        return None

    # OLCULEN KIMLIK YOLU: uye KUTULARIN olculen kimligi Kraken etiketinden
    # ayrisiyor (100 kutunun 78'inde). Kapsama sinavi Kraken taxid'i uzerinden
    # yapilirsa dogru dislamalar "kapsamiyor" gorunur - 2026-08-10'da
    # Petriella, Microascaceae ve Cloacimonas'ta tam olarak bu oldu:
    # kutular Kraken'de "Trichoderma"/"Cloacimonas", olculende "Petriella"/
    # "Planctomycetales". Sinav OLCULEN kimlige gore yapilir.
    ESD = {}
    ey = os.path.join(kok, 'screening', 'kimlik_taxid_esdegeri.tsv')
    if os.path.exists(ey):
        for l in io.open(ey, encoding='utf-8'):
            l = l.rstrip('\n')
            if not l.strip() or l.startswith('#'):
                continue
            p = l.split('\t')
            if len(p) >= 2 and p[1].strip().isdigit():
                ESD[p[0].strip()] = p[1].strip()
    KUT = oku_uye_kutular(os.path.join(kok, 'uyelik_yeniden_turetme_uyelik_20260803.tsv'))
    OLC = oku_olculen_kimlik(os.path.join(kok, 'TUM_KIMLIK_SONUC', 'tum_kutu_kimlikleri.tsv'))
    yaz(u'  member bin table : %d targets' % len(KUT))
    yaz(u'  measured identity: %d bins' % len(OLC))
    olculen_yol = bool(KUT and OLC)
    if not olculen_yol:
        yaz(u'  WARNING: the measured-identity route is CLOSED (file missing). The check will')
        yaz(u'  use the Kraken taxid instead, which can produce FALSE ALARMS.')

    gerek = set()
    icin = {}
    kimlik_notu = {}
    for hedef, exc in sorted(H.items()):
        gerek.update(exc)
        k = esle(hedef)
        if olculen_yol:
            kk = k if (k and k in KUT) else (hedef if hedef in KUT else None)
            if kk:
                tl, notlar = [], []
                for kutu in KUT[kk]:
                    kim = OLC.get(kutu, '')
                    if not kim:
                        notlar.append((kutu, 'olculen kimlik yok', None))
                        continue
                    tx = ESD.get(kutu) or ESD.get(kim.strip())
                    if tx:
                        notlar.append((kutu, kim[:50] + ' [esdeger tablosu]', tx))
                        tl.append(tx)
                        gerek.add(tx)
                        continue
                    for aday in ad_adaylari(kim):
                        tx = ad_taxid(aday)
                        if tx:
                            notlar.append((kutu, aday, tx))
                            break
                    if tx:
                        tl.append(tx)
                    else:
                        notlar.append((kutu, kim[:60] + ' (taxid bulunamadi)', None))
                icin[hedef] = (kk, tl)
                kimlik_notu[hedef] = notlar
                gerek.update(tl)
                continue
        if k is None:
            continue
        icin[hedef] = (k, U[k])
        gerek.update(U[k])
    gerek = {t for t in gerek if t and t.isdigit()}
    yaz(u'  cekilecek taxid  : %d' % len(gerek))
    Z = soy_zincirleri(sorted(gerek), yaz)

    yaz('')
    yaz(u'%-46s %-9s %-7s %s' % ('hedef', 'dislanan', 'uye', 'sonuc'))
    yaz(u'-' * 78)
    kotu = []
    eksik_uyelik = []
    for hedef in sorted(H):
        exc = H[hedef]
        if hedef not in icin:
            eksik_uyelik.append(hedef)
            yaz(u'%-46s %-9s %-7s NOT IN THE MEMBERSHIP TABLE - the member set for this target was never derived (verdict: MEMBERSHIP NOT VERIFIED' % (hedef[:46], ','.join(exc)[:9], '-'))
            continue
        k, uyeler = icin[hedef]
        uyeler = [m for m in uyeler if m and str(m).isdigit()]
        if not uyeler:
            yaz(u'%-46s %-9s %-7s uye taxid cozulemedi - denetlenemedi'
                % (hedef[:46], ','.join(exc)[:9], '-'))
            eksik_uyelik.append(hedef)
            continue
        disarida = []
        bagisik = []
        for m in uyeler:
            ad, zincir = Z.get(m, ('', set()))
            if not zincir:
                disarida.append((m, 'NCBI zinciri cekilemedi'))
                continue
            if not any(e in zincir or e == m for e in exc):
                if m in IST.get(hedef, {}):
                    bagisik.append((m, ad or '?', IST[hedef][m]))
                else:
                    disarida.append((m, ad or '?'))
        if bagisik and not disarida:
            yaz(u'%-46s %-9s %-7d ok (%d accepted exceptions)'
                % (hedef[:46], ','.join(exc)[:9], len(uyeler), len(bagisik)))
            for m, ad, sb in bagisik:
                yaz(u'      istisna: taxid %-9s %-24s %s' % (m, ad[:24], sb[:60]))
            continue
        if disarida:
            kotu.append((hedef, exc, disarida))
            yaz(u'%-46s %-9s %-7d DOES NOT COVER (%d members left out)'
                % (hedef[:46], ','.join(exc)[:9], len(uyeler), len(disarida)))
        else:
            yaz(u'%-46s %-9s %-7d ok' % (hedef[:46], ','.join(exc)[:9], len(uyeler)))

    yaz('')
    if kotu:
        yaz(u'  TARGETS WITH INCOMPLETE COVERAGE - these members ARE the target itself')
        yaz(u'  yet are counted as "off-target product":')
        for hedef, exc, dis in kotu:
            yaz(u'   * %s  (dislanan: %s)' % (hedef, ','.join(exc)))
            for m, ad in dis:
                zincir_ad = Z.get(m, ('', set()))[0]
                yaz(u'       disarida: taxid %-9s %s' % (m, zincir_ad or ad))
    else:
        yaz(u'  For every target, the excluded taxon covers ALL of its members.')
    if eksik_uyelik:
        yaz(u'  NOT IN THE MEMBERSHIP TABLE (could not be audited): %s' % ', '.join(eksik_uyelik))
    yaz(u'=' * 78)
    return 1 if (kotu or eksik_uyelik) else 0


if __name__ == '__main__':
    sys.exit(main())
