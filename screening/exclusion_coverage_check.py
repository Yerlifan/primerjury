# -*- coding: utf-8 -*-
"""THE COVERAGE AUDIT OF THE EXCLUSION MAP: is hedef_taxid.tsv right?

THE QUESTION
------------
hedef_taxid.tsv names, for each target, the taxon to be excluded at NCBI. That taxon
MUST hold ALL of the target's members. If it does not, a member left outside counts
as an "off target product" although it IS the target itself, and the pair looks
dirty for no reason.

This audit leaves no room for guessing: it pulls the NCBI lineage of every member
taxid and asks whether the excluded taxid is REALLY in that lineage.

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
    """target -> the member BIN names (A1-1_2209 and the like). The BIN, not the Kraken taxid."""
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
    """bin -> the MEASURED identity text. NOT the Kraken label, our own measurement."""
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
    """Name candidates that can be looked up at NCBI, taken from the measured identity
    text, from broad to narrow.

    There are two forms:
      "Petriella setifera"                       -> a direct name
      "adlandirilamayan soy - EN YAKIN KAYIT: ...Fungi;Dikarya;Ascomycota;..."
                                                 -> the NARROWEST name in the lineage

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
    """taxid -> (name, the set of taxids in the lineage). It is pulled from NCBI and not
    memorised.

    The XML is parsed properly. A plain text regex DOES NOT WORK HERE: <LineageEx> holds
    nested <Taxon> elements and splitting with a regex makes each lineage ring look like
    a separate record; the coverage test then says even the correct pairs "do not cover"
    (which is exactly what happened on 2026-08-10).

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
            # an AKA or merged taxid: the number asked for and the number returned can differ
            for aka in tx.findall('AkaTaxIds/TaxId'):
                if aka.text:
                    out[aka.text] = (ad, zincir)
            out[tid] = (ad, zincir)
        time.sleep(0.4)
        yaz(u'  lineages pulled: %d/%d' % (min(i + 80, len(tl)), len(tl)))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', dest='kok', default='.')
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

    # a name match: the membership table's key can be short (an example:
    # "Proteolitik_Synergistaceae" against "Proteolitik_Synergistaceae (eski ad: ...)").
    def esle(ad):
        if ad in U:
            return ad
        for k in U:
            if ad.startswith(k) or k.startswith(ad):
                return k
        return None

    # THE MEASURED IDENTITY ROUTE: the measured identity of the member BINS departs from
    # the Kraken label (in 78 of 100 bins). If the coverage test is made over the Kraken
    # taxid, correct exclusions look as though they "do not cover"; on 2026-08-10 that is
    # exactly what happened with Petriella, Microascaceae and Cloacimonas: the bins are
    # "Trichoderma" and "Cloacimonas" in Kraken and "Petriella" and "Planctomycetales"
    # under measurement. The test is made against the MEASURED identity.
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
