# -*- coding: utf-8 -*-
"""OLCULMEYEN TAKSONLAR ICIN KUTU URET  -  "hicbir veri bosta kalmasin".

DURUM
-----
Simdiye kadar barkod basina yalnizca 5 kutu (Bracken'in ilk 5'i) ayristirilip
olculdu. Bu, tur duzeyi okumalarin ortanca %83'unu kapsiyor ama EN KOTU
barkodda %27'de kaliyor. Yani sorun ortalamada degil dagilimda.

Bu betik kalan taksonlari da kutuya cevirir: Kraken'in okuma duzeyi
ciktisindan (edited_barcodeNN_output) o taksonun KLADINA atanmis okumalari
ham fastq'tan cikarir ve "fastq files/<kutu>/" altina mevcut adlandirmayla
yazar. Boylece var olan konsensus ureteci ve kimlik dogrulama zinciri onlari
kendiliginden gorur - yeni bir olcum yolu ACILMAZ, var olan yol beslenir.

KLAD KURALI (olculdu, varsayilmadi)
-----------------------------------
Bir kutunun okumalari, o taksonun TAM KLADINA atanmis okumalardir; yalniz
taksonun kendisine atananlar DEGIL. Ornek: A1-1'deki 2223 kutusunda 19 757
okuma var ama Kraken'de tam olarak 2223'e atanmis okuma SIFIR - hepsi sus
duzeyindeki 990316'ya (M. soehngenii GP6) atanmis. Klad kurali uygulanmadan
o kutu bos cikardi.

KALIBRASYON KAPISI - ONCE KANITLA, SONRA URET
---------------------------------------------
Betik once ELDE HAZIR olan kutulari yeniden hesaplar ve okuma kimligi
kumelerinin BIREBIR tutup tutmadigina bakar. Tutmuyorsa hicbir sey uretmez ve
duser. 2026-08-10'da 14/14 birebir tutmustur; yeni bir barkod ya da yeni bir
Kraken kosusu bu varsayimi bozarsa kapi onu yakalar.

Kosum:
    python verification/build_bins.py --root . --calibration-only
    python verification/build_bins.py --root . --coverage 0.95
"""
from __future__ import print_function

import argparse
import collections
import io
import os
import sys
import time

FASTQ_KOK = 'fastq files'
KRAKEN_KOK = 'kraken results'


# --------------------------------------------------------------- agac ve klad
def agac_kur(report):
    """Kraken raporunun GIRINTISINDEN takson agacini kurar (ust haritasi)."""
    ust = {}
    yol = []
    for l in io.open(report, encoding='utf-8', errors='replace'):
        p = l.rstrip('\n').split('\t')
        if len(p) < 6:
            continue
        ad = p[5]
        girinti = (len(ad) - len(ad.lstrip(' '))) // 2
        tx = p[4].strip()
        while len(yol) > girinti:
            yol.pop()
        ust[tx] = yol[-1] if yol else None
        yol.append(tx)
    return ust


def klad_kumesi(ust, kok):
    cocuk = collections.defaultdict(list)
    for c, u in ust.items():
        if u:
            cocuk[u].append(c)
    out = {str(kok)}
    yigin = [str(kok)]
    while yigin:
        x = yigin.pop()
        for c in cocuk.get(x, []):
            if c not in out:
                out.add(c)
                yigin.append(c)
    return out


def tur_satirlari(report):
    """(taxid, ad, klad_okuma) - yalniz rutbe S."""
    out = []
    for l in io.open(report, encoding='utf-8', errors='replace'):
        p = l.rstrip('\n').split('\t')
        if len(p) < 6:
            continue
        if p[3].strip() != 'S':
            continue
        try:
            n = int(p[1])
        except ValueError:
            continue
        out.append((p[4].strip(), p[5].strip(), n))
    out.sort(key=lambda x: -x[2])
    return out


def kraken_okuma_taxid(output_yolu):
    """okuma_id -> atanan taxid (yalniz siniflandirilmis okumalar)."""
    d = {}
    for l in io.open(output_yolu, encoding='utf-8', errors='replace'):
        p = l.rstrip('\n').split('\t')
        if len(p) >= 3 and p[0] == 'C':
            d[p[1].split()[0]] = p[2].strip()
    return d


def fastq_idler(yol):
    ids = set()
    with io.open(yol, encoding='utf-8', errors='replace') as fh:
        for i, l in enumerate(fh):
            if i % 4 == 0:
                ids.add(l[1:].split()[0] if len(l) > 1 else '')
    ids.discard('')
    return ids


# ------------------------------------------------------- barkod <-> kutu esleme
def esleme_kur(kok, yaz):
    """sinif klasorlerinden kutu <-> barkod eslemesi. VARSAYILMAZ, DOGRULANIR."""
    esleme = {}
    kr = os.path.join(kok, KRAKEN_KOK)
    if not os.path.isdir(kr):
        return esleme
    for sinif in sorted(os.listdir(kr)):
        d = os.path.join(kr, sinif)
        if not os.path.isdir(d):
            continue
        barkodlar = sorted(set(
            f.replace('edited_', '').replace('_kraken2.report', '')
            for f in os.listdir(d) if f.endswith('_kraken2.report')))
        kutular = sorted(k for k in os.listdir(os.path.join(kok, FASTQ_KOK))
                         if k.startswith(sinif + '-'))
        if len(barkodlar) != len(kutular):
            yaz(u'  WARNING: class %s has %d barcodes but %d bin directories, so no mapping could be built and this class is SKIPPED.'
                % (sinif, len(barkodlar), len(kutular)))
            continue
        for kutu, bc in zip(kutular, barkodlar):
            esleme[kutu] = (sinif, bc)
    return esleme


# ------------------------------------------------------------- kalibrasyon
def kalibrasyon(kok, esleme, yaz):
    """Elde HAZIR kutular klad kuraliyla birebir yeniden uretilebiliyor mu."""
    tam = sapan = 0
    ayrinti = []
    for kutu, (sinif, bc) in sorted(esleme.items()):
        d = os.path.join(kok, FASTQ_KOK, kutu)
        rep = os.path.join(kok, KRAKEN_KOK, sinif, 'edited_%s_kraken2.report' % bc)
        out = os.path.join(kok, KRAKEN_KOK, sinif, 'edited_%s_output' % bc)
        if not (os.path.isdir(d) and os.path.exists(rep) and os.path.exists(out)):
            continue
        ust = agac_kur(rep)
        atama = kraken_okuma_taxid(out)
        ters = collections.defaultdict(set)
        for rid, tx in atama.items():
            ters[tx].add(rid)
        for f in sorted(os.listdir(d)):
            if not f.endswith('.fastq'):
                continue
            tx = f.split('_')[-1].replace('.fastq', '')
            if not tx.isdigit():
                continue
            a = fastq_idler(os.path.join(d, f))
            kume = klad_kumesi(ust, tx)
            b = set()
            for t in kume:
                b |= ters.get(t, set())
            if a == b:
                tam += 1
            else:
                sapan += 1
                ayrinti.append(u'%s taxid %s: fastq %d, klad %d, yalniz fastq %d, '
                               u'yalniz klad %d'
                               % (kutu, tx, len(a), len(b), len(a - b), len(b - a)))
    yaz(u'  kalibrasyon: %d birebir, %d sapan' % (tam, sapan))
    for x in ayrinti[:10]:
        yaz(u'      %s' % x)
    return tam, sapan


# ------------------------------------------------------------- ham okuma bulma
def ham_fastq_bul(kok, bc, ek_kokler):
    adaylar = []
    for k in ek_kokler:
        if not k or not os.path.isdir(k):
            continue
        for dirpath, _dn, fn in os.walk(k):
            for f in fn:
                if bc in f and (f.endswith('.fastq') or f.endswith('.fq')
                                or f.endswith('.fastq.gz')):
                    adaylar.append(os.path.join(dirpath, f))
    return sorted(set(adaylar))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', dest='kok', default='.')
    p.add_argument('--coverage', dest='kapsam', type=float, default=0.95,
                   help='species-read coverage to reach in each barcode (0-1)')
    p.add_argument('--min-reads', dest='asgari_okuma', type=int, default=50,
                   help='taxa with fewer reads than this are not turned into '
                        'a bin, because the consensus would not be '
                        'trustworthy')
    p.add_argument('--raw-root', dest='ham_kok', default='',
                   help='ham barkod fastq directory; bos ise bilinen yerlere bakilir')
    p.add_argument('--calibration-only', dest='yalniz_kalibrasyon', action='store_true')
    p.add_argument('--plan-only', dest='yalniz_plan', action='store_true',
                   help='write down what would be produced, and produce '
                        'nothing')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)

    def yaz(s=''):
        print(s, flush=True)

    yaz(u'=' * 78)
    yaz(u'  BUILDING BINS FOR THE UNMEASURED TAXA   %s'
        % time.strftime('%Y-%m-%d %H:%M'))
    yaz(u'=' * 78)

    esleme = esleme_kur(kok, yaz)
    yaz(u'  bin to barcode mapping: %d bin directories' % len(esleme))
    if not esleme:
        yaz(u'  ESLEME KURULAMADI - cikiliyor.')
        return 1

    tam, sapan = kalibrasyon(kok, esleme, yaz)
    if sapan:
        yaz(u'')
        yaz(u'  CALIBRATION FAILED. The clade rule will re-derive the bins already on disk')
        yaz(u'  cannot produce it; applying it to the remaining taxa WOULD BE WRONG.')
        yaz(u'  No file was produced.')
        return 1
    if tam == 0:
        yaz(u'  No ready bin to calibrate was found, so the gate never opened.')
        return 1
    yaz(u'  The calibration PASSED (%d/%d). Production is allowed.' % (tam, tam))
    if a.yalniz_kalibrasyon:
        return 0

    # --- plan -------------------------------------------------------------
    plan = []
    for kutu, (sinif, bc) in sorted(esleme.items()):
        rep = os.path.join(kok, KRAKEN_KOK, sinif, 'edited_%s_kraken2.report' % bc)
        if not os.path.exists(rep):
            continue
        d = os.path.join(kok, FASTQ_KOK, kutu)
        var = set()
        if os.path.isdir(d):
            for f in os.listdir(d):
                t = f.split('_')[-1].replace('.fastq', '')
                if t.isdigit():
                    var.add(t)
        turler = tur_satirlari(rep)
        toplam = sum(n for _t, _ad, n in turler) or 1
        birikim = 0
        for tx, ad, n in turler:
            birikim += n
            if tx in var:
                continue
            if n < a.asgari_okuma:
                continue
            plan.append(dict(kutu=kutu, sinif=sinif, barkod=bc, taxid=tx,
                             ad=ad, okuma=n, birikim=birikim / float(toplam)))
            if birikim / float(toplam) >= a.kapsam:
                break

    yaz(u'')
    yaz(u'  PLAN: coverage target %%%d, minimum reads %d'
        % (int(a.kapsam * 100), a.asgari_okuma))
    yaz(u'  new bins to create: %d' % len(plan))
    sinif_say = collections.Counter(x['sinif'] for x in plan)
    for s, n in sorted(sinif_say.items()):
        yaz(u'      %-4s %d' % (s, n))
    py = os.path.join(kok, 'KUTU_URETIM_PLANI.tsv')
    with io.open(py, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# Bins to be produced. Generated %s, coverage target %.2f\n'
                 % (time.strftime('%Y-%m-%d %H:%M'), a.kapsam))
        fh.write(u'kutu\tsinif\tbarkod\ttaxid\tkraken_adi\tokuma\tbirikimli_kapsam\n')
        for x in plan:
            fh.write(u'%s\t%s\t%s\t%s\t%s\t%d\t%.4f\n'
                     % (x['kutu'], x['sinif'], x['barkod'], x['taxid'],
                        x['ad'], x['okuma'], x['birikim']))
    yaz(u'  written: %s' % py)
    if a.yalniz_plan:
        return 0

    # --- uretim -----------------------------------------------------------
    ek = [a.ham_kok] if a.ham_kok else []
    ek += [os.path.join(kok, os.pardir, 'tools', 'SONUCLAR', 'fastq files'),
           os.path.join(kok, os.pardir, 'tools', 'SONUCLAR'),
           os.path.join(kok, 'HAM_OKUMA')]
    uretilen = 0
    atlanan = []
    bc_grup = collections.defaultdict(list)
    for x in plan:
        bc_grup[(x['sinif'], x['barkod'], x['kutu'])].append(x)
    for (sinif, bc, kutu), isler in sorted(bc_grup.items()):
        ham = ham_fastq_bul(kok, bc, ek)
        if not ham:
            atlanan.append(u'%s: ham fastq bulunamadi (%s)' % (bc, kutu))
            continue
        rep = os.path.join(kok, KRAKEN_KOK, sinif, 'edited_%s_kraken2.report' % bc)
        out = os.path.join(kok, KRAKEN_KOK, sinif, 'edited_%s_output' % bc)
        ust = agac_kur(rep)
        atama = kraken_okuma_taxid(out)
        ters = collections.defaultdict(set)
        for rid, tx in atama.items():
            ters[tx].add(rid)
        istenen = {}
        for x in isler:
            kume = klad_kumesi(ust, x['taxid'])
            s = set()
            for t in kume:
                s |= ters.get(t, set())
            if s:
                istenen[x['taxid']] = s
        if not istenen:
            continue
        hedef_dizin = os.path.join(kok, FASTQ_KOK, kutu)
        if not os.path.isdir(hedef_dizin):
            os.makedirs(hedef_dizin)
        acik = {}
        try:
            for tx in istenen:
                acik[tx] = io.open(
                    os.path.join(hedef_dizin, '%s_reads_%s.fastq'
                                 % (kutu.replace('-', '_'), tx)),
                    'w', encoding='utf-8', newline='')
            hangi = {}
            for tx, s in istenen.items():
                for rid in s:
                    hangi[rid] = tx
            for hy in ham:
                with io.open(hy, encoding='utf-8', errors='replace') as fh:
                    while True:
                        b1 = fh.readline()
                        if not b1:
                            break
                        b2 = fh.readline()
                        b3 = fh.readline()
                        b4 = fh.readline()
                        rid = b1[1:].split()[0] if len(b1) > 1 else ''
                        tx = hangi.get(rid)
                        if tx:
                            acik[tx].write(b1 + b2 + b3 + b4)
        finally:
            for f in acik.values():
                f.close()
        uretilen += len(istenen)
        yaz(u'  %-8s %-12s %d bins written' % (kutu, bc, len(istenen)))

    yaz(u'')
    yaz(u'  bins created: %d' % uretilen)
    for x in atlanan[:20]:
        yaz(u'  SKIPPED: %s' % x)
    yaz(u'  Next step: consensus generation and the identity verification chain')
    yaz(u'    verification/full_chain.py -> option 6 (consensus), then stage G')
    yaz(u'=' * 78)
    return 0


if __name__ == '__main__':
    sys.exit(main())
