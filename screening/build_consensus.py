# -*- coding: utf-8 -*-
"""OPTION 5 - REGENERATE the consensuses from the raw reads.

THE METHOD
----------
For each bin the reads are seated on a TEMPLATE (the existing consensus, or the
best read if there is none) with shared k-mer anchors, and voted on column by
column. TWO SEPARATE methods are used:

  A) QUALITY WEIGHTED (with a lowered threshold): every vote is weighted by the
     read's Phred quality at that base; the base is called when the weighted
     majority is >= 50%.
  B) MAJORITY VOTE (unweighted): a flat plurality.

If the two methods give THE SAME base, the base is DEFINITE. Where they differ,
that column is MASKED (an N is written). NO DEGENERATE BASE IS PRODUCED; the
output holds only A/C/G/T/N.

KNOWN TRAPS (all of them handled deliberately)
----------------------------------------------
1. The read length filter: class A2 reads are 4.2-4.5 kb and class F2 ~3.7 kb.
   The old 200-3000 bp filter discarded both classes entirely. The filter here is
   200-%(max)d bp, and which range was used per class is reported.
2. The consensus ORIENTATION: reads come in both directions and some of the stored
   consensuses are in the reverse complement direction. Every read is anchored in
   both directions and whichever holds more anchors is used, so the orientation is
   NORMALISED.
3. THE 5' END IS NOT TRIMMED: the output is at template length; uncovered columns
   are not clipped but marked with N. That way the coordinates are not disturbed.

"""
# -------------------------------------------------------------------------
# build_consensus.py regenerates each bin's consensus sequence from the raw fastq
#                     reads by two independent methods; where the methods disagree,
#                     the column is masked with N.
#
# INPUT  : the reads under "fastq files" through hedefler.kutular() (at most
#          MAX_OKUMA reads per bin, with a fixed seed); as the template, the
#          canonical consensus hedefler.konsensusler() gives, and failing that the
#          best read; orientation normalisation through orientation.py.
# OUTPUT : one fasta per bin under SCREENING_RESULT/konsensus_yeni/;
#          SCREENING_RESULT/KONSENSUS_YENIDEN_URETIM.md and
#          konsensus_yeniden_uretim.tsv (calistir returns those two paths); plus a
#          checkpoint per bin. The files produced are NOT USED DIRECTLY; run_all.py
#          first takes them into the canonical directory through
#          build_canonical.py --priority yeni, passing them through orientation
#          normalisation.
# CALLED BY: verification/full_chain.py key 6 (--mode konsensus) and the 3rd stage
#          inside key 9 (hepsi.calistir -> konsensus_uret.calistir).
#
# ORIENTATION is handled at three separate points here (the gate at the start of the
# run, converting the template to canonical, and a final measurement before the
# output is written), because the new consensus produced becomes the next stage's
# backbone; a backbone produced in reverse silently zeroes every in-silico PCR result
# that follows.
# -------------------------------------------------------------------------
import os, gzip, glob, time, json, math
from collections import defaultdict, Counter
from . import config as C
from . import engine_gateway, targets as H, checks
from . import orientation

K = 15                      # capa k-mer boyu
MIN_DERINLIK = 5            # the minimum reads for a column to be called
AGIRLIK_ESIGI = 0.50        # kalite agirlikli yontemin esigi (DUSURULMUS)
COGUNLUK_ESIGI = 0.50
OKUMA_MIN, OKUMA_MAX = 200, 8000     # A2 (4,5 kb) ve F2 (3,7 kb) ELENMEZ
MAX_OKUMA = 600             # the most reads used per bin


def fastq_oku(yol, n=MAX_OKUMA, tohum=C.NUMUNE_TOHUM):
    '(sequence, quality) pairs; the length filter DOES NOT drop the long amplicon classes.'
    ac = gzip.open if yol.endswith('.gz') else open
    out = []
    with ac(yol, 'rt', errors='ignore') as fh:
        s = None
        for i, line in enumerate(fh):
            m = i % 4
            if m == 1:
                s = line.strip().upper()
            elif m == 3 and s is not None:
                q = line.rstrip('\n')
                if OKUMA_MIN <= len(s) <= OKUMA_MAX and len(q) >= len(s):
                    out.append((s, q[:len(s)]))
                s = None
    if n and len(out) > n:
        import random
        out = random.Random(tohum).sample(out, n)
    return out


def _kmer_indeksi(sablon):
    d = defaultdict(list)
    for i in range(len(sablon) - K + 1):
        d[sablon[i:i + K]].append(i)
    # repeated k-mers are UNRELIABLE as anchors and are discarded
    return {k: v[0] for k, v in d.items() if len(v) == 1}


def _rc_kalite(s, q):
    return engine_gateway.rc(s), q[::-1]


def _capala(okuma, kal, idx):
    'Fit the read onto the template. Returns segments of (template_pos, read_pos, length).'
    vur = []
    for i in range(0, len(okuma) - K + 1):
        t = idx.get(okuma[i:i + K])
        if t is not None:
            vur.append((t - i, t, i))
    if not vur:
        return [], 0
    say = Counter(d for d, _, _ in vur)
    segler = []
    for diag, n in say.most_common(8):
        if n < 3:
            continue
        noktalar = sorted((t, r) for d, t, r in vur if d == diag)
        segler.append((noktalar[0][0], noktalar[0][1],
                       noktalar[-1][0] - noktalar[0][0] + K))
    return segler, len(vur)


def kutu_konsensusu(reads, sablon):
    'A per column vote table, two methods, and the sequence they agree on.'
    # Every read is anchored to the template BOTH as itself AND as its reverse
    # complement, and the direction holding more anchors is chosen. That is required,
    # because nanopore reads come in both directions; anchored in one direction only,
    # roughly half the reads can vote on no column at all and the depth halves. The
    # counter of the chosen direction (ters) appears in the report.
    #
    # The two methods are kept SEPARATE and a base is called only when both give THE
    # SAME base; a column where they differ is masked with N. That way the uncertainty
    # is not hidden but written into the sequence. NO degenerate IUPAC base IS PRODUCED:
    # a degenerate base means more than one oligo at synthesis and is against the panel's
    # rules, so the uncertainty stays an N.
    idx = _kmer_indeksi(sablon)
    L = len(sablon)
    agir = [defaultdict(float) for _ in range(L)]
    duz = [Counter() for _ in range(L)]
    kullanilan = ters = 0
    for s, q in reads:
        s = engine_gateway.clean(s)
        ileri = _capala(s, q, idx)
        geri_s, geri_q = _rc_kalite(s, q)
        geri = _capala(geri_s, geri_q, idx)
        if ileri[1] >= geri[1]:
            segler, seq, qual = ileri[0], s, q
        else:
            segler, seq, qual = geri[0], geri_s, geri_q
            ters += 1
        if not segler:
            continue
        kullanilan += 1
        for tpos, rpos, uz in segler:
            for j in range(uz):
                t = tpos + j; r = rpos + j
                if not (0 <= t < L) or not (0 <= r < len(seq)):
                    continue
                b = seq[r]
                if b not in 'ACGT':
                    continue
                Q = ord(qual[r]) - 33 if r < len(qual) else 20
                w = 1.0 - 10 ** (-max(Q, 1) / 10.0)
                agir[t][b] += w
                duz[t][b] += 1

    A = []; B = []
    for i in range(L):
        ta = sum(agir[i].values()); td = sum(duz[i].values())
        if td < MIN_DERINLIK:
            A.append('N'); B.append('N'); continue
        ba, va = max(agir[i].items(), key=lambda x: x[1]) if agir[i] else ('N', 0)
        bb, vb = duz[i].most_common(1)[0] if duz[i] else ('N', 0)
        A.append(ba if ta and va / ta >= AGIRLIK_ESIGI else 'N')
        B.append(bb if td and vb / td >= COGUNLUK_ESIGI else 'N')
    uzlasi = ''.join(a if (a == b and a != 'N') else 'N' for a, b in zip(A, B))
    return dict(A=''.join(A), B=''.join(B), uzlasi=uzlasi,
                kullanilan=kullanilan, ters=ters, toplam=len(reads),
                derinlik_ort=round(sum(sum(d.values()) for d in duz) / max(L, 1), 1))


def _sablon_sec(kutu, kons_haritasi, reads):
    """ORIENTATION NORMALISATION (2026-08-02): the template is converted to CANONICAL
        (SENSE).
        How it was before: the template was taken from the existing consensus and that
        directory was MIXED orientation; the reads were normalised to the template, but
        THE TEMPLATE'S OWN ORIENTATION was never normalised. The result: the output
        inherited the template's orientation (measured: konsensus_yeni was 28 antisense /
        7 sense). On a reversed consensus, in-silico PCR silently gives 0 products.

    """
    sn = orientation.sinifi(kutu)
    k = kons_haritasi.get(kutu)
    if k:
        d, karar, cev = orientation.kanonik(k['dizi'], sn)
        return d, 'the existing consensus (orientation=%s%s)' % (karar, ', KANONIGE CEVRILDI' if cev else '')
    if not reads:
        return None, 'okuma yok'
    # a bin with no consensus (an orphan): the median of the longest reads becomes the template
    sirali = sorted(reads, key=lambda x: -len(x[0]))[:15]
    s = sirali[len(sirali) // 2][0]
    d, karar, cev = orientation.kanonik(engine_gateway.clean(s), orientation.sinifi(kutu))
    return d, ('there is no template; the median of the longest reads was used (orientation=%s%s)'
               % (karar, ', KANONIGE CEVRILDI' if cev else ''))


def calistir(yaz, sure, yalniz=None, yeniden=False):
    # THE GATE: no consensus is produced before the orientation normalisation test passes.
    # The reason: if the output of this step comes out in the wrong orientation the whole
    # night is wasted (on a reversed consensus, in-silico PCR SILENTLY gives 0 products).
    from . import self_test as _KS
    if not _KS.yon_sinamasi(yaz):
        yaz(u'   THE ORIENTATION TEST FAILED, so consensus generation WAS NOT STARTED.')
        yaz(u'   Run this first: python screening/build_canonical.py --root .')
        return None
    from .run_all import yon_kapisi
    _ok, _m = yon_kapisi(yaz, 'the consensus regeneration')
    for _x in _m:
        yaz('  ' + _x)
    if not _ok:
        yaz('')
        yaz(u'  *** INPUT VERIFICATION FAILED - THIS STAGE WAS NOT STARTED ***')
        yaz(u'  Cause: the consensus sequences to be read are not canonical. On a reverse-oriented')
        yaz(u'  consensus, in-silico PCR returns 0 products without any warning,')
        yaz(u'  so the whole run would silently produce a wrong result.')
        yaz(u'  Fix:    python3 screening/build_canonical.py --root . --rerun')
        raise SystemExit(2)

    checks.hazirla()
    cikti = os.path.join(C.CIKTI, 'konsensus_yeni')
    os.makedirs(cikti, exist_ok=True)
    kut = H.kutular()
    kons = {k['kutu']: k for k in H.konsensusler()}
    if yalniz:
        kut = [k for k in kut if yalniz.lower() in k['kutu'].lower()]
    yaz('=' * 78)
    yaz(u'  REGENERATING THE CONSENSUSES FROM THE RAW READS')
    yaz('=' * 78)
    yaz(u'  number of bins: %d' % len(kut))
    yaz(u'  read filter   : %d-%d bp  (A2 at ~4.5 kb and F2 at ~3.7 kb are NOT discarded)'
        % (OKUMA_MIN, OKUMA_MAX))
    yaz(u'  method        : (A) quality weighted, threshold %.2f   (B) majority vote, threshold %.2f'
        % (AGIRLIK_ESIGI, COGUNLUK_ESIGI))
    yaz(u'  a column with no agreement becomes N (masked). NO DEGENERATE BASE IS PRODUCED.')
    yaz(u'  the 5\' end is NOT TRIMMED. The output keeps the template length, and any uncovered column is N.')
    yaz('')
    tahmin = len(kut) * 12
    yaz(u'ESTIMATED TIME: ~%s   (each bin is written to disk as it finishes, so it resumes)\n' % sure(tahmin))

    satirlar = []
    t0 = time.time()
    for i, k in enumerate(kut, 1):
        kp = os.path.join(C.KONTROL, 'konsensus_%s.json'
                          % ''.join(c if c.isalnum() else '_' for c in k['kutu']))
        if os.path.exists(kp) and not yeniden:
            try:
                _v = json.load(open(kp, encoding='utf-8'))
                if not checks.ayar_uyuyor(_v):
                    raise ValueError('ayar degisti')
                satirlar.append(json.load(open(kp, encoding='utf-8')))
                yaz(u'[%d/%d] %-18s (from the previous run)' % (i, len(kut), k['kutu']))
                continue
            except Exception:
                pass   # bayat/bozuk: silmeye calisma, uzerine yazilacak
        reads = fastq_oku(k['yol'])
        sablon, kaynak = _sablon_sec(k['kutu'], kons, reads)
        if sablon is None or len(sablon) < 200:
            yaz(u'[%d/%d] %-18s SKIPPED (%s)' % (i, len(kut), k['kutu'], kaynak))
            continue
        r = kutu_konsensusu(reads, sablon)
        eski = kons.get(k['kutu'], {}).get('dizi', '')
        ayni = sum(1 for a, b in zip(r['uzlasi'], eski) if a == b) if eski else 0
        farkli = sum(1 for a, b in zip(r['uzlasi'], eski)
                     if a != b and a != 'N' and b != 'N') if eski else 0
        satir = dict(kutu=k['kutu'], sinif=k['sinif'], taxid=k['taxid'],
                     sablon_kaynagi=kaynak, sablon_uz=len(sablon),
                     okuma_toplam=r['toplam'], okuma_kullanilan=r['kullanilan'],
                     ters_yonde=r['ters'], derinlik_ort=r['derinlik_ort'],
                     N_sayisi=r['uzlasi'].count('N'),
                     N_yuzde=round(100.0 * r['uzlasi'].count('N') / max(len(sablon), 1), 2),
                     eski_ile_ayni=ayni, eski_ile_farkli=farkli,
                     eski_uz=len(eski),
                     yontem_A_N=r['A'].count('N'), yontem_B_N=r['B'].count('N'),
                     yontemler_ayrildi=sum(1 for a, b in zip(r['A'], r['B'])
                                           if a != b and 'N' not in (a, b)),
                     dizi=r['uzlasi'])
        # THE OUTPUT IS CONVERTED TO CANONICAL - the last seat belt. Even when the template
        # is canonical, the output is measured once more here; if it comes out UNCERTAIN that
        # is written into the header.
        _ck, _karar, _cev = orientation.kanonik(r['uzlasi'], orientation.sinifi(k['kutu']))
        r['uzlasi'] = _ck
        satir['cikti_yon'] = _karar
        satir['cikti_cevrildi'] = 'EVET' if _cev else 'hayir'
        with open(os.path.join(cikti, '%s_yeniden_konsensus.fasta' % k['kutu']),
                  'w', encoding='utf-8') as fh:
            fh.write('>%s yeniden_uretildi sablon=%s derinlik=%s N=%d kanonik=%s cevrildi=%d\n'
                     % (k['kutu'], kaynak, r['derinlik_ort'], satir['N_sayisi'],
                        _karar, int(_cev)))
            for j in range(0, len(r['uzlasi']), 70):
                fh.write(r['uzlasi'][j:j + 70] + '\n')
        satir['_ayar'] = dict(checks.AYAR)
        with open(kp, 'w', encoding='utf-8') as fh:
            json.dump(satir, fh, ensure_ascii=False, default=str)
        satirlar.append(satir)
        yaz(u'[%d/%d] %-18s %5d bp  depth %5s  N %%%-5s  bases differing from the old one: %d'
            % (i, len(kut), k['kutu'], len(sablon), r['derinlik_ort'],
               satir['N_yuzde'], farkli))
        gecen = time.time() - t0
        print('       gecen %s  tahmini kalan %s' % (
            sure(gecen), sure(gecen / i * (len(kut) - i))), flush=True)

    # the report is produced from ALL the checkpoint files (earlier runs included)
    hepsi = []
    for f in sorted(os.listdir(C.KONTROL)):
        if f.startswith('konsensus_') and f.endswith('.json'):
            try:
                hepsi.append(json.load(open(os.path.join(C.KONTROL, f), encoding='utf-8')))
            except Exception:
                pass
    yollar = rapor_yaz(hepsi or satirlar, cikti)
    yaz('')
    yaz('=' * 78)
    yaz(u'  CONSENSUS GENERATION FINISHED (%s)' % sure(time.time() - t0))
    for p in yollar:
        yaz('    %s' % p)
    yaz('=' * 78)
    return yollar


SUT = ['kutu', 'sinif', 'taxid', 'sablon_kaynagi', 'sablon_uz', 'eski_uz',
       'okuma_toplam', 'okuma_kullanilan', 'ters_yonde', 'derinlik_ort',
       'N_sayisi', 'N_yuzde', 'cikti_yon', 'cikti_cevrildi', 'yontem_A_N', 'yontem_B_N', 'yontemler_ayrildi',
       'eski_ile_ayni', 'eski_ile_farkli']


def rapor_yaz(satirlar, cikti):
    import csv
    tsv = os.path.join(C.CIKTI, 'konsensus_yeniden_uretim.tsv')
    with open(tsv, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(SUT)
        for s in satirlar:
            w.writerow([s.get(k, '') for k in SUT])
    md = os.path.join(C.CIKTI, 'KONSENSUS_YENIDEN_URETIM.md')
    L = []; A = L.append
    A(u'# The consensuses reproduced from the raw reads')
    A('')
    A(u'Generated: %s' % time.strftime('%Y-%m-%d %H:%M:%S'))
    A('')
    A(u'The new fasta files: `%s`' % cikti)
    A('')
    A('## Yontem')
    A('')
    A(u'For every bin the reads were seated on a template with shared %d-mer anchors, and each column was called **by two separate methods**:' % K)
    A('')
    A('- **A - kalite agirlikli (esigi dusurulmus):** her oy Phred kalitesiyle '
      'agirliklandi; agirlikli cogunluk >= %.0f%% ise baz cagrildi.' % (100 * AGIRLIK_ESIGI))
    A(u'- **B, the majority vote:** unweighted plurality, threshold %.0f%%.' % (100 * COGUNLUK_ESIGI))
    A('')
    A(u'If the two methods give the same base, the base is **definite**. If they diverge the column is **masked with N**. The output holds only `A/C/G/T/N`: **there is no degenerate base**.')
    A('')
    A(u'### The known traps, every one handled deliberately')
    A('')
    A(u'| trap | what was done |')
    A('|---|---|')
    A(u'| The read length filter was eliminating class A2 (4,2-4,5 kb) and F2 (~3,7 kb) | The filter is **%d-%d bp**; neither class is eliminated |' % (OKUMA_MIN, OKUMA_MAX))
    A(u'| The consensus orientation was mixed (some reverse complemented) | Every read was anchored **in both directions** and the direction holding more anchors was chosen, so the orientation is normalised (the `ters_yonde` column) |')
    A(u'| The 5\' end was being cut | The output is **the length of the template**; the uncovered columns were not trimmed but marked with `N`, so the coordinates are kept |')
    A('')
    A('## Sonuclar')
    A('')
    A(u'| bin | class | length | depth | N %% | the methods diverged | bases differing from the old one | template |')
    A('|---|---|---|---|---|---|---|---|')
    for s in satirlar:
        A('| %s | %s | %s | %s | %s | %s | %s | %s |' % (
            s['kutu'], s['sinif'], s['sablon_uz'], s['derinlik_ort'], s['N_yuzde'],
            s['yontemler_ayrildi'], s['eski_ile_farkli'], s['sablon_kaynagi']))
    A('')
    A('## Sinirlar - bu ONEMLI')
    A('')
    A(u'- This is a **template based recall**, NOT a de novo assembly. It corrects column level substitution errors; large structural differences that are not in the template (long insertions or deletions) it **cannot find**.')
    A(u'- The anchors are non-repeating %d-mers; in repetitive regions the coverage falls and those columns become `N`. The `N_yuzde` column shows that.' % K)
    A(u'- The `yontemler_ayrildi` column is the number of columns where the two methods contradict each other. Where that number is high the consensus of that bin is **not reliable**; look at it first.')
    A(u'- The new consensuses were **not written into the panel automatically**. Before they are used, the `eski_ile_farkli` and `N_yuzde` columns have to be reviewed.')
    A('')
    with open(md, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(L))
    return [md, tsv]
