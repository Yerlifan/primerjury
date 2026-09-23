# -*- coding: utf-8 -*-
u"""İŞLEV GENİ PRİMER TASARIMI (18.09.2026, v2: çoklu hizalama).

# ==== TOPLANTI NOTU ====
# NE ISE YARAR : İşlev grubunu tek çiftle ölçmenin yolu işin kendi genidir (mcrA metan üretimi, amoA amonyak oksidasyonu, pmoA metan
#                oksidasyonu, dsrB + soxB kükürt döngüsü); 16S akrabalığı gruplar, işlevi değil. rRNA amplikon okumalarında bu genler YOK:
#                tasarım ve kapsama ölçümü NCBI referans dizileri üzerinde yapılır (tools/function_gene_references.py); numunenin kendi
#                okumalarında sınanamaz — bu sınır çıktıya yazılır.
#                v1 (ikili hizalama) uzak akrabalarda primer yerini eşleştiremiyordu (payda 0/0, ölçüldü 18.09) → v2: bütün diziler
#                MAFFT ile birlikte hizalanır, primerler hizalamanın uzlaşı dizisinden (sütun doluluğu ≥ %50) üretilir, kapsama doğrudan
#                hizalama sütunlarından sayılır (boşluk = uyumsuzluk; 3' son 3 baz birebir; ≤ 1 uyumsuzluk).
#                Ölçüt: SYBR Green qPCR koşulu (P3 aşağıda; tek 60 °C, ürün 70–250 bp). Dejenere: ikili IUPAC, 3' son 3 baz sabit,
#                ≤ --degenerate varyant; bütün varyantlar Tm 57–64 °C, GC %35–65, saç tokası/dimer Tm ≤ 45 °C (QPCR_KOSUL).
#                Tek dosya: tasarım ölçütleri, IUPAC açılımı ve Tm ölçümü burada; dışarıdan yalnız primer3-py ve MAFFT gerekir.
# NASIL KOSAR  : python3 function_gene_design.py --gene mcrA --fasta <gene_ref.fna> --out <prefix> [--degenerate 8] [--max-product 250]
#                [--top 5] [--mafft <path>]   (MAFFT: --mafft, yoksa MAFFT ortam değişkeni, yoksa PATH'teki mafft)
# CIKIS        : <prefix>.aln (MAFFT; varsa yeniden kullanılır), <prefix>_adaylar.tsv (cins başına ürün veren/sınanabilen),
#                <prefix>_secim.tsv (en iyi --top)
# =======================
"""
from __future__ import print_function
import argparse
import collections
import io
import os
import subprocess
import sys

EN_COK_MM = 1        # ürün için primer başına en çok uyumsuzluk (3' son UC3 baz birebir)
SUTUN_DOLU = 0.5     # uzlaşı dizisine giren sütunun en az doluluğu
YER_DOLU = 0.8       # primer yerinin sınanabilir sayılması için o dizideki en az doluluk (boşluksuz baz oranı)
UC3 = 3              # 3' uçta birebir olması gereken baz sayısı

# SYBR Green qPCR (QuantiNova, Rotor-Gene Q; iki adımlı 60 °C) ölçütleri — primer3 tasarım ayarları
P3 = {
    'PRIMER_TASK': 'generic', 'PRIMER_PICK_LEFT_PRIMER': 1, 'PRIMER_PICK_RIGHT_PRIMER': 1, 'PRIMER_PICK_INTERNAL_OLIGO': 0,
    'PRIMER_OPT_SIZE': 20, 'PRIMER_MIN_SIZE': 18, 'PRIMER_MAX_SIZE': 24,
    'PRIMER_OPT_TM': 60.0, 'PRIMER_MIN_TM': 58.0, 'PRIMER_MAX_TM': 62.0, 'PRIMER_PAIR_MAX_DIFF_TM': 1.5,
    'PRIMER_MIN_GC': 40.0, 'PRIMER_MAX_GC': 60.0, 'PRIMER_GC_CLAMP': 1, 'PRIMER_MAX_END_GC': 3,
    'PRIMER_MAX_POLY_X': 4, 'PRIMER_MAX_NS_ACCEPTED': 0,
    'PRIMER_PRODUCT_SIZE_RANGE': [[70, 200]], 'PRIMER_PRODUCT_OPT_SIZE': 120,
    'PRIMER_SALT_MONOVALENT': 50.0, 'PRIMER_SALT_DIVALENT': 3.0, 'PRIMER_DNTP_CONC': 0.8, 'PRIMER_DNA_CONC': 250.0,
    'PRIMER_TM_FORMULA': 1, 'PRIMER_SALT_CORRECTIONS': 1,
    'PRIMER_MAX_SELF_ANY_TH': 45.0, 'PRIMER_MAX_SELF_END_TH': 35.0, 'PRIMER_MAX_HAIRPIN_TH': 40.0,
    'PRIMER_PAIR_MAX_COMPL_ANY_TH': 45.0, 'PRIMER_PAIR_MAX_COMPL_END_TH': 35.0,
    'PRIMER_MAX_END_STABILITY': 9.0, 'PRIMER_NUM_RETURN': 30,
}
PENCERE, ADIM = 600, 300
# Tepkime koşulu, P3 tuzlarıyla aynı. Koşulsuz primer3.calc_tm Mg 0 / 50 nM varsayar ve Tm'yi ~4 °C düşük verir (ölçüldü 14.09).
QPCR_KOSUL = dict(mv_conc=50.0, dv_conc=3.0, dntp_conc=0.8, dna_conc=250.0)

IUPAC = {u'A': u'A', u'C': u'C', u'G': u'G', u'T': u'T', u'K': u'GT', u'M': u'AC', u'R': u'AG', u'Y': u'CT', u'S': u'CG', u'W': u'AT', u'N': u'ACGT'}
IKILI = dict((frozenset(v), k) for k, v in IUPAC.items() if len(v) == 2)
_RC = {u'A': u'T', u'C': u'G', u'G': u'C', u'T': u'A', u'N': u'N', u'K': u'M', u'M': u'K', u'R': u'Y', u'Y': u'R', u'S': u'S', u'W': u'W'}


def rc(s):
    return u''.join(_RC.get(c, u'N') for c in reversed(s.upper()))


def gc(s):
    s = s.upper()
    return 100.0 * sum(1 for c in s if c in u'GC') / len(s)


def varyantlar(p, en_cok=64):
    u"""IUPAC'lı primer -> somut diziler (K -> G ve T). Dejenere olmayan primer [p] verir."""
    out = [u'']
    for c in p.upper():
        out = [x + y for x in out for y in IUPAC.get(c, c)]
        if len(out) > en_cok:
            raise ValueError(u'primer %d diziden fazla açılıyor: %s' % (en_cok, p))
    return out


def adaylar(sablon, p3):
    u"""primer3 ile 600 bp pencerelerde (300 bp adım) aday çiftler; aynı çift bir kez."""
    import primer3
    out = {}
    for b0 in range(0, max(1, len(sablon) - PENCERE + 1), ADIM):
        win = sablon[b0:b0 + PENCERE]
        if len(win) < 150:
            continue
        try:
            res = primer3.bindings.design_primers({'SEQUENCE_ID': 'w', 'SEQUENCE_TEMPLATE': win}, p3)
        except Exception as e:   # D3 notu: primer3 pencere hatası (N'li dizi) atlanır, sayısı bildirilir
            print(u'    primer3 pencere %d atlandı: %s' % (b0, e))
            continue
        for i in range(res.get('PRIMER_PAIR_NUM_RETURNED', 0)):
            F = res['PRIMER_LEFT_%d_SEQUENCE' % i]
            R = res['PRIMER_RIGHT_%d_SEQUENCE' % i]
            if (F, R) not in out:
                out[(F, R)] = dict(F=F, R=R, urun=res['PRIMER_PAIR_%d_PRODUCT_SIZE' % i])
    return list(out.values())


def olc(F, R):
    u"""Dejenere çiftin varyantlarında Tm / GC / saç tokası / dimer (QPCR_KOSUL). Ölçüt dışıysa None."""
    import primer3
    fv, rv = varyantlar(F), varyantlar(R)
    tmf = [primer3.calc_tm(v, **QPCR_KOSUL) for v in fv]
    tmr = [primer3.calc_tm(v, **QPCR_KOSUL) for v in rv]
    if min(tmf + tmr) < 57.0 or max(tmf + tmr) > 64.0:
        return None
    if abs(sum(tmf) / len(tmf) - sum(tmr) / len(tmr)) > 2.5:
        return None
    if any(not (35.0 <= gc(v) <= 65.0) for v in fv + rv):
        return None
    hp = max(primer3.calc_hairpin(v, **QPCR_KOSUL).tm for v in fv + rv)
    ho = max(primer3.calc_homodimer(v, **QPCR_KOSUL).tm for v in fv + rv)
    he = max(primer3.calc_heterodimer(f, r, **QPCR_KOSUL).tm for f in fv for r in rv)
    if hp > 45.0 or ho > 45.0 or he > 45.0:
        return None
    aralik = lambda t: (u'%.1f' % t[0]) if round(max(t), 1) == round(min(t), 1) else (u'%.1f–%.1f' % (min(t), max(t)))
    return dict(tmF=aralik(tmf), tmR=aralik(tmr), gcF=round(sum(gc(v) for v in fv) / len(fv), 1), gcR=round(sum(gc(v) for v in rv) / len(rv), 1),
                dimer=round(max(ho, he), 1), varyant=len(fv) * len(rv))


def fasta_oku(yol):
    u"""(erişim, cins, başlık, dizi) listesi."""
    out, bas, par = [], None, []
    for s in io.open(yol, encoding='utf-8'):
        s = s.rstrip(u'\n')
        if s.startswith(u'>'):
            if bas:
                out.append(bas + (u''.join(par),))
            p = s[1:].split(u'|')
            bas = (p[0], p[1] if len(p) > 1 else u'?', p[2] if len(p) > 2 else u'')
            par = []
        elif bas:
            par.append(s.strip().upper())
    if bas:
        out.append(bas + (u''.join(par),))
    return [(a, c, b, d) for a, c, b, d in out if len(d.replace(u'-', u'')) >= 150]


def hizala(fasta, cikti, mafft):
    u"""MAFFT çoklu hizalama (varsa hazır dosya kullanılır)."""
    if os.path.exists(cikti) and os.path.getsize(cikti) > 0:
        print(u'  hizalama hazır: %s' % cikti)
        return cikti
    with io.open(cikti, 'wb') as g:
        subprocess.check_call([mafft, u'--auto', u'--thread', u'2', u'--quiet', fasta], stdout=g)
    return cikti


def uzlasi(dizi):
    u"""Sütun doluluğu ≥ SUTUN_DOLU olan sütunlardan uzlaşı dizisi; (uzlaşı, uzlaşı konumu -> sütun)."""
    n = len(dizi[0])
    kons, sut = [], []
    for j in range(n):
        baz = [d[j] for d in dizi if d[j] in u'ACGT']
        if len(baz) < SUTUN_DOLU * len(dizi):
            continue
        kons.append(collections.Counter(baz).most_common(1)[0][0])
        sut.append(j)
    return u''.join(kons), sut


def yer_al(d, sutunlar):
    u"""Dizinin bu sütunlardaki bazları (boşluk '-' kalır)."""
    return u''.join(d[j] for j in sutunlar)


def uyumsuzluk(yer, p):
    u"""IUPAC'lı p ile yer arasındaki uyumsuzluk; 3' son UC3 baz birebir değilse ya da yer çok boşluklu ise None."""
    if sum(1 for c in yer if c in u'ACGT') < YER_DOLU * len(p):
        return None
    for k in range(len(p) - UC3, len(p)):
        if yer[k] not in IUPAC.get(p[k], p[k]):
            return None
    return sum(1 for k in range(len(p) - UC3) if yer[k] not in IUPAC.get(p[k], p[k]))


def dejenere_kur(p, yerler, en_cok_varyant):
    u"""Üye bağlanma yerlerine göre ikili IUPAC; en çok kazandıran konumlardan başlayarak, 3' son UC3 baz sabit."""
    if not yerler:
        return p
    sayim = [collections.Counter(y[k] for y in yerler if y[k] in u'ACGT') for k in range(len(p))]
    aday = []
    for k in range(len(p) - UC3):
        c = sayim[k]
        ikinci = [b for b, _n in c.most_common() if b != p[k]]
        if not ikinci:
            continue
        aday.append((c[ikinci[0]], k, ikinci[0]))
    out, varyant = list(p), 1
    for _n, k, b in sorted(aday, reverse=True):
        if varyant * 2 > en_cok_varyant or p[k] not in u'ACGT':
            break
        out[k] = IKILI[frozenset([p[k], b])]
        varyant *= 2
    return u''.join(out)


def degerlendir(F, R, sF, sR, kayit, hiz):
    u"""Cins -> [ürün veren, sınanabilen]; ürün boyları (hizalamadaki boşluksuz baz sayısı)."""
    say = collections.defaultdict(lambda: [0, 0])
    boy = []
    for (acc, cins, _b, _d) in kayit:
        d = hiz[acc]
        yF, yR = yer_al(d, sF), yer_al(d, sR)
        mF, mR = uyumsuzluk(yF, F), uyumsuzluk(rc(yR), R)
        if mF is None and mR is None:
            continue
        say[cins][1] += 1
        if mF is not None and mR is not None and mF <= EN_COK_MM and mR <= EN_COK_MM:
            say[cins][0] += 1
            boy.append(sum(1 for c in d[sF[0]:sR[-1] + 1] if c in u'ACGT'))
    return say, boy


def main(argv=None):
    ay = argparse.ArgumentParser(description=u'Function-gene qPCR primers from a MAFFT alignment of reference sequences')
    ay.add_argument('--gene', '--gen', dest='gen', required=True, help=u'gene label written to the output (mcrA, amoA, ...)')
    ay.add_argument('--fasta', required=True, help=u'reference FASTA, headers >accession|genus|title (function_gene_references.py)')
    ay.add_argument('--out', '--cikti', dest='cikti', required=True, help=u'output prefix')
    ay.add_argument('--degenerate', '--dejenere', dest='dejenere', type=int, default=8, help=u'most variants per primer pair')
    ay.add_argument('--max-product', '--urun-en-cok', dest='urun_en_cok', type=int, default=250, help=u'longest product, bp')
    ay.add_argument('--top', '--en-iyi', dest='en_iyi', type=int, default=5, help=u'pairs written to <prefix>_secim.tsv')
    ay.add_argument('--mafft', default=os.environ.get('MAFFT', u'mafft'), help=u'MAFFT executable')
    a = ay.parse_args(argv)
    kayit = fasta_oku(a.fasta)
    if len(kayit) < 3:
        raise SystemExit(u'DURDU: %s içinde 3\'ten az dizi' % a.fasta)
    cinsler = sorted(set(k[1] for k in kayit))
    print(u'%s: %d dizi, %d cins (%s)' % (a.gen, len(kayit), len(cinsler),
                                          u', '.join(u'%s %d' % (c, sum(1 for x in kayit if x[1] == c)) for c in cinsler)))
    aln = hizala(a.fasta, a.cikti + u'.aln', a.mafft)
    hk = fasta_oku(aln)
    hiz = dict((k[0], k[3].upper()) for k in hk)
    kayit = [k for k in kayit if k[0] in hiz]
    boylar = set(len(v) for v in hiz.values())
    if len(boylar) != 1:
        raise SystemExit(u'DURDU: hizalama satır boyları farklı: %s' % sorted(boylar)[:5])
    kons, sut = uzlasi([hiz[k[0]] for k in kayit])
    print(u'  hizalama %d sütun; uzlaşı %d baz (sütun doluluğu ≥ %%%d)' % (list(boylar)[0], len(kons), 100 * SUTUN_DOLU))
    if len(kons) < 150:
        raise SystemExit(u'DURDU: uzlaşı dizisi çok kısa (%d baz) — diziler çok farklı' % len(kons))
    ad_ = adaylar(kons, dict(P3, PRIMER_PRODUCT_SIZE_RANGE=[[70, a.urun_en_cok]]))
    print(u'  %d aday (uzlaşı dizisinden)' % len(ad_))
    sys.stdout.flush()
    sat = []
    for x in ad_:
        p0 = kons.find(x['F'])
        p1 = kons.rfind(rc(x['R']))
        if p0 < 0 or p1 < 0 or p1 <= p0:
            continue
        sF = sut[p0:p0 + len(x['F'])]
        sR = sut[p1:p1 + len(x['R'])]
        yerF = [yer_al(hiz[k[0]], sF) for k in kayit]
        yerR = [rc(yer_al(hiz[k[0]], sR)) for k in kayit]
        temiz = lambda ys: [y for y in ys if sum(1 for c in y if c in u'ACGT') >= YER_DOLU * len(y)]
        for kip, F, R in ((u'düz', x['F'], x['R']),
                          (u'dejenere', dejenere_kur(x['F'], temiz(yerF), a.dejenere), dejenere_kur(x['R'], temiz(yerR), a.dejenere))):
            if kip == u'dejenere' and (F, R) == (x['F'], x['R']):
                continue
            m = olc(F, R)
            if not m:
                continue
            say, boy = degerlendir(F, R, sF, sR, kayit, hiz)
            kapsanan = sum(1 for c in cinsler if say[c][1] and say[c][0] >= 0.5 * say[c][1])
            sinanan = sum(1 for c in cinsler if say[c][1])
            oran = [1.0 * say[c][0] / say[c][1] for c in cinsler if say[c][1]]
            sat.append(dict(gen=a.gen, kip=kip, F=F, R=R, duz_F=x['F'], duz_R=x['R'],
                            urun=u'%d' % (sum(boy) // len(boy)) if boy else u'%d' % x['urun'],
                            varyant=m['varyant'], tmF=m['tmF'], tmR=m['tmR'], gcF=m['gcF'], gcR=m['gcR'], dimer=m['dimer'],
                            kapsanan=kapsanan, sinanan=sinanan, ortalama=round(100.0 * sum(oran) / len(oran), 1) if oran else 0.0,
                            toplam=u'%d/%d' % (sum(say[c][0] for c in cinsler), sum(say[c][1] for c in cinsler)),
                            cins_metin=u'; '.join(u'%s %d/%d' % (c, say[c][0], say[c][1]) for c in cinsler),
                            eksik=u'; '.join(c for c in cinsler if not say[c][1] or say[c][0] < 0.5 * say[c][1]) or u'-'))
    if not sat:
        raise SystemExit(u'DURDU: hiçbir aday ölçütleri tutmadı')
    bas = [u'gen', u'kip', u'F', u'R', u'duz_F', u'duz_R', u'urun', u'varyant', u'tmF', u'tmR', u'gcF', u'gcR', u'dimer', u'kapsanan',
           u'sinanan', u'ortalama', u'toplam', u'cins_metin', u'eksik']
    sat.sort(key=lambda d: (-d['kapsanan'], -d['ortalama'], d['varyant']))
    yaz = lambda yol, satirlar: io.open(yol, 'w', encoding='utf-8', newline='\n').write(
        u'\n'.join([u'\t'.join(bas)] + [u'\t'.join(u'%s' % d[k] for k in bas) for d in satirlar]) + u'\n')
    yaz(a.cikti + u'_adaylar.tsv', sat)
    en_iyi, gorulen = [], set()
    for d in sat:
        if (d['F'], d['R']) in gorulen:
            continue
        gorulen.add((d['F'], d['R']))
        en_iyi.append(d)
        if len(en_iyi) >= a.en_iyi:
            break
    yaz(a.cikti + u'_secim.tsv', en_iyi)
    print(u'  %d aday; en iyi %d:' % (len(sat), len(en_iyi)))
    for d in en_iyi:
        print(u'    %s / %s  ürün %s bp | %s | kapsanan cins %d/%d, toplam %s | %s | eksik: %s'
              % (d['F'], d['R'], d['urun'], d['kip'], d['kapsanan'], d['sinanan'], d['toplam'], d['cins_metin'], d['eksik']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
