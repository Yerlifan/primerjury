# -*- coding: utf-8 -*-
u"""İŞLEV GENİ REFERANS DİZİLERİ — NCBI nuccore (18.09.2026).

# ==== TOPLANTI NOTU ====
# NE ISE YARAR : 18.09 kullanıcı: "Evet, dört gen için tasarla ve başlat" (mcrA, amoA, pmoA, dsrB + soxB). Bizim okumalarımız 16S/ITS olduğundan
#                bu genler okumalarda yok; tasarım için numunede bulunan (PAK kutuları) işlevli cinslerin gen kayıtları NCBI nuccore'dan alınır
#                (E-utilities esearch + efetch; TOOL/EMAIL verilir, istekler arası ≥ 0,4 s). Gen + cins sorgusu, 150–5 000 bp (genom kaydı değil),
#                cins başına en çok --en-cok kayıt. Başlık: >erişim|cins|NCBI başlığı.
# NASIL KOSAR  : python3 function_gene_references.py --gene mcrA --genera "Methanosarcina;Methanothrix" --out <gene_ref.fna> [--max 80] [--email ...]
# CIKIS        : <gen_ref.fna> ve <gen_ref.fna>.tsv (erişim, cins, uzunluk, başlık, sorgu); dosya varsa üstüne yazılmaz
# =======================
"""
from __future__ import print_function
import argparse
import io
import json
import os
import sys
import time
import urllib.parse
import urllib.request

EU = u'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
ARAC = u'primerjury-function-gene'
ARA = 0.4
GEN_TERIM = {
    u'mcrA': u'(mcrA[Gene Name] OR "methyl coenzyme M reductase alpha"[Title] OR "methyl-coenzyme M reductase subunit alpha"[Title] OR '
             u'"methyl coenzyme M reductase subunit alpha"[Title] OR "methyl-coenzyme M reductase alpha"[Title])',
    u'amoA': u'(amoA[Gene Name] OR "ammonia monooxygenase subunit A"[Title] OR "ammonia monooxygenase alpha"[Title])',
    u'pmoA': u'(pmoA[Gene Name] OR "particulate methane monooxygenase subunit A"[Title] OR "particulate methane monooxygenase alpha"[Title] OR '
             u'"methane monooxygenase subunit A"[Title])',
    u'dsrB': u'(dsrB[Gene Name] OR "dissimilatory sulfite reductase beta"[Title] OR "dissimilatory sulfite reductase subunit beta"[Title] OR '
             u'"sulfite reductase, dissimilatory-type subunit beta"[Title])',
    u'soxB': u'(soxB[Gene Name] OR "sulfur oxidation protein SoxB"[Title] OR "thiosulfohydrolase SoxB"[Title] OR "SoxB"[Title])',
}


def istek(yol, param):
    param = dict(param, tool=ARAC)
    url = EU + yol + u'?' + urllib.parse.urlencode(param)
    for deneme in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as y:
                veri = y.read().decode('utf-8', 'replace')
            time.sleep(ARA)
            return veri
        except Exception as e:   # D3 notu: ağ hatası 4 kez denenir, sonra yükseltilir
            print(u'   ağ hatası (%s), deneme %d' % (type(e).__name__, deneme + 1))
            time.sleep(5 * (deneme + 1))
    raise SystemExit(u'DURDU: NCBI yanıt vermedi: %s' % yol)


def cds_al(gen, cins, ep, en_cok_kayit=2):
    u"""Gen yalnız genom kaydında geçiyorsa (soxB, dsrB'de bazı cinsler): genomun CDS'lerini al, [gene=<gen>] olanları süz."""
    terim = u'%s AND "%s"[Organism]' % (GEN_TERIM[gen], cins)
    js = json.loads(istek(u'esearch.fcgi', dict(ep, db=u'nuccore', term=terim, retmax=en_cok_kayit, retmode=u'json', sort=u'relevance')))
    ids = js.get(u'esearchresult', {}).get(u'idlist', [])
    if not ids:
        return []
    metin = istek(u'efetch.fcgi', dict(ep, db=u'nuccore', id=u','.join(ids), rettype=u'fasta_cds_na', retmode=u'text'))
    out = []
    for k in metin.split(u'>'):
        if not k.strip():
            continue
        bas, _, dizi = k.partition(u'\n')
        if (u'[gene=%s]' % gen).lower() not in bas.lower() and (u'[protein=%s]' % gen).lower() not in bas.lower():
            continue
        dizi = u''.join(dizi.split()).upper()
        acc = bas.split()[0].replace(u'lcl|', u'')
        if len(dizi) >= 150:
            out.append((acc, bas[len(bas.split()[0]):].strip().replace(u'|', u'/').replace(u'\t', u' '), dizi))
    return out


def main(argv=None):
    ay = argparse.ArgumentParser()
    ay.add_argument('--gene', '--gen', dest='gen', required=True, choices=sorted(GEN_TERIM), help=u'gene name (mcrA, amoA, pmoA, dsrB, soxB)')
    ay.add_argument('--genera', '--cinsler', dest='cinsler', required=True, help=u'noktalı virgüllü cins adları (NCBI [Organism] sorgusu)')
    ay.add_argument('--out', '--cikti', dest='cikti', required=True, help=u'output FASTA')
    ay.add_argument('--max', '--en-cok', dest='en_cok', type=int, default=80, help=u'records per genus')
    ay.add_argument('--email', '--eposta', dest='eposta', default=u'', help=u'NCBI E-utilities e-mail')
    ay.add_argument('--cds-min', '--cds-en-az', dest='cds_en_az', type=int, default=5, help=u'bu sayının altında kayıt bulunan cinste genom CDS\'lerinden alınır')
    a = ay.parse_args(argv)
    if os.path.exists(a.cikti):
        raise SystemExit(u'DURDU: %s zaten var — üstüne yazılmaz' % a.cikti)
    ep = {u'email': a.eposta} if a.eposta else {}
    fa, tb = [], [u'erişim\tcins\tuzunluk\tbaşlık\tsorgu']
    gorulen = set()
    for cins in [x.strip() for x in a.cinsler.split(u';') if x.strip()]:
        terim = u'%s AND "%s"[Organism] AND 150:5000[SLEN]' % (GEN_TERIM[a.gen], cins)
        js = json.loads(istek(u'esearch.fcgi', dict(ep, db=u'nuccore', term=terim, retmax=a.en_cok, retmode=u'json', sort=u'relevance')))
        ids = js.get(u'esearchresult', {}).get(u'idlist', [])
        toplam = js.get(u'esearchresult', {}).get(u'count', u'?')
        yeni = [i for i in ids if i not in gorulen]
        gorulen.update(yeni)
        n = 0
        if yeni:
            metin = istek(u'efetch.fcgi', dict(ep, db=u'nuccore', id=u','.join(yeni), rettype=u'fasta', retmode=u'text'))
            kayit = [k for k in metin.split(u'>') if k.strip()]
            for k in kayit:
                bas, _, dizi = k.partition(u'\n')
                dizi = u''.join(dizi.split()).upper()
                acc = bas.split()[0]
                baslik = bas[len(acc):].strip().replace(u'|', u'/').replace(u'\t', u' ')
                if len(dizi) < 150:
                    continue
                fa.append(u'>%s|%s|%s\n%s' % (acc, cins, baslik, dizi))
                tb.append(u'\t'.join([acc, cins, u'%d' % len(dizi), baslik, terim]))
                n += 1
        c = 0
        if n < a.cds_en_az:   # gen amplikon kaydı olarak yok (ör. soxB); genom CDS'lerinden alınır
            for acc, baslik, dizi in cds_al(a.gen, cins, ep):
                if acc in gorulen:
                    continue
                gorulen.add(acc)
                fa.append(u'>%s|%s|%s' % (acc, cins, baslik) + u'\n' + dizi)
                tb.append(u'\t'.join([acc, cins, u'%d' % len(dizi), baslik, u'genom CDS']))
                c += 1
        print(u'  %-6s %-34s NCBI kayıt %6s | alınan %3d%s' % (a.gen, cins, toplam, n, (u' + genom CDS %d' % c) if c else u''))
        sys.stdout.flush()
    io.open(a.cikti, 'w', encoding='utf-8', newline='\n').write(u'\n'.join(fa) + u'\n')
    io.open(a.cikti + u'.tsv', 'w', encoding='utf-8', newline='\n').write(u'\n'.join(tb) + u'\n')
    print(u'  yazıldı %s: %d dizi' % (a.cikti, len(fa)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
