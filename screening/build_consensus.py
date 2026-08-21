# -*- coding: utf-8 -*-
"""SECENEK 5 - Konsensusleri ham okumalardan YENIDEN URET.

YONTEM
------
Her kutu icin okumalar bir SABLONA (mevcut konsensus; yoksa en iyi okuma)
ortak k-mer capalariyla oturtulur ve sutun sutun oy verilir. Iki AYRI yontem
kullanilir:

  A) KALITE AGIRLIKLI (esigi dusurulmus): her oy, okumanin o bazdaki Phred
     kalitesiyle agirliklanir; agirlikli cogunluk >= %50 ise baz cagrilir.
  B) COGUNLUK OYU (agirliksiz): duz plurality.

Iki yontem AYNI bazi verirse baz KESINDIR. Ayrilirlarsa o sutun MASKELENIR
(N yazilir). DEJENERE BAZ URETILMEZ - cikti yalniz A/C/G/T/N icerir.

BILINEN TUZAKLAR (hepsi bilerek ele alindi)
-------------------------------------------
1. Okuma uzunluk filtresi: A2 sinifi okumalari 4,2-4,5 kb, F2 sinifi ~3,7 kb.
   Eski 200-3000 bp filtresi bu iki sinifi tumden eliyordu. Burada filtre
   200-%(max)d bp'dir ve sinif bazinda hangi araligin kullanildigi raporlanir.
2. Konsensus YONU: okumalar iki yonde gelir, saklanan konsensuslerin bir kismi
   ters tumleyen yonde. Her okuma iki yonde de capalanir, cok capa hangisinde
   tutuyorsa o yon kullanilir - yon NORMALIZE edilir.
3. 5' UC KESILMEZ: cikti sablon boyundadir; kapsanmayan sutunlar kirpilmaz,
   N ile isaretlenir. Boylece koordinatlar bozulmaz.
"""
# ---------------------------------------------------------------------------
# build_consensus.py — her kutunun konsensus dizisini ham fastq okumalarindan
#                     iki bagimsiz yontemle yeniden uretir; yontemler ayrilan
#                     sutunlari N ile maskeler.
#
# GIRDI  : hedefler.kutular() ile "fastq files" altindaki okumalar (kutu basina
#          en cok MAX_OKUMA okuma, sabit tohum); sablon olarak
#          hedefler.konsensusler()'in verdigi kanonik konsensus, o yoksa en iyi
#          okuma; orientation.py ile yon normalizasyonu.
# CIKTI  : KAPSAMLI_ARAMA_SONUC/konsensus_yeni/ altina kutu basina fasta;
#          KAPSAMLI_ARAMA_SONUC/KONSENSUS_YENIDEN_URETIM.md ve
#          konsensus_yeniden_uretim.tsv (calistir bu iki yolu dondurur);
#          kutu basina kontrol noktasi. Uretilen dosyalar DOGRUDAN kullanilmaz;
#          run_all.py once build_canonical.py --oncelik yeni ile bunlari kanonik
#          klasore yon normalizasyonundan gecirerek alir.
# CAGRAN : verification/full_chain.py tusu 6 (--mod konsensus) ve tusu 9 icindeki
#          3. asama (hepsi.calistir -> konsensus_uret.calistir).
#
# YON burada uc ayri noktada ele alinir (kosunun basindaki kapi, sablonun
# kanonige cevrilmesi ve cikti yazilmadan onceki son olcum) cunku uretilen yeni
# konsensus bir sonraki asamanin omurgasi olur; ters yonde uretilmis bir omurga
# sonraki butun in-silico PCR sonuclarini sessizce sifirlar.
# ---------------------------------------------------------------------------
import os, gzip, glob, time, json, math
from collections import defaultdict, Counter
from . import config as C
from . import engine_gateway, hedefler as H, kontrol
from . import orientation

K = 15                      # capa k-mer boyu
MIN_DERINLIK = 5            # bir sutunun cagrilabilmesi icin en az okuma
AGIRLIK_ESIGI = 0.50        # kalite agirlikli yontemin esigi (DUSURULMUS)
COGUNLUK_ESIGI = 0.50
OKUMA_MIN, OKUMA_MAX = 200, 8000     # A2 (4,5 kb) ve F2 (3,7 kb) ELENMEZ
MAX_OKUMA = 600             # kutu basina kullanilan en fazla okuma


def fastq_oku(yol, n=MAX_OKUMA, tohum=C.NUMUNE_TOHUM):
    """(dizi, kalite) ciftleri; uzunluk filtresi A2/F2'yi ELEMEZ."""
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
    # tekrar eden k-mer'ler capa olarak GUVENILMEZ - atilir
    return {k: v[0] for k, v in d.items() if len(v) == 1}


def _rc_kalite(s, q):
    return motor.rc(s), q[::-1]


def _capala(okuma, kal, idx):
    """Okumayi sablona oturt. Donen: [(sablon_pos, okuma_pos, uzunluk)] segmentleri."""
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
    """Sutun bazli oy tablosu -> iki yontem -> uzlasi dizisi."""
    # Her okuma HEM duz HEM ters tumleyen halinde sablona capalanir ve daha cok
    # capa tutan yon secilir. Bu sart, cunku nanopore okumalari cift yonlu gelir;
    # tek yonde capalansa okumalarin yaklasik yarisi hicbir sutuna oy veremez ve
    # derinlik yariya duser. Secilen yon sayaci (ters) raporda gorunur.
    #
    # Iki yontem AYRI tutulur ve ancak ikisi AYNI bazi verirse baz cagrilir;
    # ayrildiklari sutun N ile maskelenir. Boylece belirsizlik gizlenmez, diziye
    # yazilir. Dejenere IUPAC bazi URETILMEZ: dejenere baz sentezde birden fazla
    # oligo demektir ve panelin kurallarina aykiridir - belirsizlik N kalir.
    idx = _kmer_indeksi(sablon)
    L = len(sablon)
    agir = [defaultdict(float) for _ in range(L)]
    duz = [Counter() for _ in range(L)]
    kullanilan = ters = 0
    for s, q in reads:
        s = motor.clean(s)
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
    """YON NORMALIZASYONU (2026-08-02): sablon KANONIGE (SENSE) cevrilir.
    Onceki hal: sablon mevcut konsensusten aliniyordu ve o klasor KARISIK yonluydu;
    okumalar sablona normalize ediliyordu ama SABLONUN KENDI YONU hic normalize
    edilmiyordu. Sonuc: cikti sablonun yonunu miras aliyordu (olculdu:
    konsensus_yeni 28 antisense / 7 sense). Ters yonlu konsensuste in-silico PCR
    sessizce 0 urun verir."""
    sn = yon.sinifi(kutu)
    k = kons_haritasi.get(kutu)
    if k:
        d, karar, cev = yon.kanonik(k['dizi'], sn)
        return d, 'mevcut konsensus (yon=%s%s)' % (karar, ', KANONIGE CEVRILDI' if cev else '')
    if not reads:
        return None, 'okuma yok'
    # konsensusu olmayan kutu (oksuz): en uzun okumalarin ortancasi sablon olur
    sirali = sorted(reads, key=lambda x: -len(x[0]))[:15]
    s = sirali[len(sirali) // 2][0]
    d, karar, cev = yon.kanonik(motor.clean(s), yon.sinifi(kutu))
    return d, ('sablon yok - en uzun okumalarin ortancasi (yon=%s%s)'
               % (karar, ', KANONIGE CEVRILDI' if cev else ''))


def calistir(yaz, sure, yalniz=None, yeniden=False):
    # KAPI: yon normalizasyonu sinamasi gecmeden konsensus uretilmez.
    # Gerekce: bu adimin ciktisi yanlis yonde cikarsa butun gece bosa gider
    # (ters yonlu konsensuste in-silico PCR SESSIZCE 0 urun verir).
    from . import self_test as _KS
    if not _KS.yon_sinamasi(yaz):
        yaz(u'   THE ORIENTATION TEST FAILED, so consensus generation WAS NOT STARTED.')
        yaz('   Once calistirin: python screening/build_canonical.py --kok .')
        return None
    from .hepsi import yon_kapisi
    _ok, _m = yon_kapisi(yaz, 'konsensus yeniden uretim')
    for _x in _m:
        yaz('  ' + _x)
    if not _ok:
        yaz('')
        yaz('  *** GIRDI DOGRULAMASI BASARISIZ - BU ASAMA BASLATILMADI ***')
        yaz(u'  Cause: the consensus sequences to be read are not canonical. On a reverse-oriented')
        yaz(u'  consensus, in-silico PCR returns 0 products without any warning,')
        yaz(u'  so the whole run would silently produce a wrong result.')
        yaz(u'  Fix:    python3 screening/build_canonical.py --root . --rerun')
        raise SystemExit(2)

    kontrol.hazirla()
    cikti = os.path.join(C.CIKTI, 'konsensus_yeni')
    os.makedirs(cikti, exist_ok=True)
    kut = H.kutular()
    kons = {k['kutu']: k for k in H.konsensusler()}
    if yalniz:
        kut = [k for k in kut if yalniz.lower() in k['kutu'].lower()]
    yaz('=' * 78)
    yaz('  KONSENSUSLERIN HAM OKUMALARDAN YENIDEN URETIMI')
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
                if not kontrol.ayar_uyuyor(_v):
                    raise ValueError('ayar degisti')
                satirlar.append(json.load(open(kp, encoding='utf-8')))
                yaz('[%d/%d] %-18s (onceki kosudan)' % (i, len(kut), k['kutu']))
                continue
            except Exception:
                pass   # bayat/bozuk: silmeye calisma, uzerine yazilacak
        reads = fastq_oku(k['yol'])
        sablon, kaynak = _sablon_sec(k['kutu'], kons, reads)
        if sablon is None or len(sablon) < 200:
            yaz('[%d/%d] %-18s ATLANDI (%s)' % (i, len(kut), k['kutu'], kaynak))
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
        # CIKTI KANONIGE CEVRILIR - son emniyet kemeri. Sablon kanonik olsa bile
        # cikti burada bir kez daha olculur; BELIRSIZ cikarsa basliga yazilir.
        _ck, _karar, _cev = yon.kanonik(r['uzlasi'], yon.sinifi(k['kutu']))
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
        satir['_ayar'] = dict(kontrol.AYAR)
        with open(kp, 'w', encoding='utf-8') as fh:
            json.dump(satir, fh, ensure_ascii=False, default=str)
        satirlar.append(satir)
        yaz('[%d/%d] %-18s %5d bp  derinlik %5s  N %%%-5s  eskiyle farkli baz: %d'
            % (i, len(kut), k['kutu'], len(sablon), r['derinlik_ort'],
               satir['N_yuzde'], farkli))
        gecen = time.time() - t0
        print('       gecen %s  tahmini kalan %s' % (
            sure(gecen), sure(gecen / i * (len(kut) - i))), flush=True)

    # rapor BUTUN kontrol dosyalarindan uretilir (onceki kosular dahil)
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
    yaz('  KONSENSUS URETIMI BITTI (%s)' % sure(time.time() - t0))
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
    A('# Konsensuslerin ham okumalardan yeniden uretimi')
    A('')
    A('Uretim zamani: %s' % time.strftime('%Y-%m-%d %H:%M:%S'))
    A('')
    A('Yeni fasta dosyalari: `%s`' % cikti)
    A('')
    A('## Yontem')
    A('')
    A('Her kutu icin okumalar bir sablona ortak %d-mer capalariyla oturtuldu ve '
      'her sutun **iki ayri yontemle** cagrildi:' % K)
    A('')
    A('- **A - kalite agirlikli (esigi dusurulmus):** her oy Phred kalitesiyle '
      'agirliklandi; agirlikli cogunluk >= %.0f%% ise baz cagrildi.' % (100 * AGIRLIK_ESIGI))
    A('- **B - cogunluk oyu:** agirliksiz plurality, esik %.0f%%.' % (100 * COGUNLUK_ESIGI))
    A('')
    A('Iki yontem ayni bazi verirse baz **kesindir**. Ayrilirlarsa sutun **N ile '
      'maskelenmistir**. Cikti yalniz `A/C/G/T/N` icerir - **dejenere baz yoktur**.')
    A('')
    A('### Bilinen tuzaklar - hepsi bilerek ele alindi')
    A('')
    A('| tuzak | ne yapildi |')
    A('|---|---|')
    A('| Okuma uzunluk filtresi A2 (4,2-4,5 kb) ve F2 (~3,7 kb) sinifini eliyordu '
      '| Filtre **%d-%d bp**; iki sinif da elenmiyor |' % (OKUMA_MIN, OKUMA_MAX))
    A('| Konsensus yonu karisik (bazilari ters tumleyen) | Her okuma **iki yonde de** '
      'capalandi, cok capa tutan yon secildi - yon normalize edildi (`ters_yonde` sutunu) |')
    A("| 5' uc kesiliyordu | Cikti **sablon boyunda**; kapsanmayan sutunlar "
      'kirpilmadi, `N` ile isaretlendi - koordinatlar korunuyor |')
    A('')
    A('## Sonuclar')
    A('')
    A('| kutu | sinif | uzunluk | derinlik | N %% | yontemler ayrildi | eskiyle farkli baz | sablon |')
    A('|---|---|---|---|---|---|---|---|')
    for s in satirlar:
        A('| %s | %s | %s | %s | %s | %s | %s | %s |' % (
            s['kutu'], s['sinif'], s['sablon_uz'], s['derinlik_ort'], s['N_yuzde'],
            s['yontemler_ayrildi'], s['eski_ile_farkli'], s['sablon_kaynagi']))
    A('')
    A('## Sinirlar - bu ONEMLI')
    A('')
    A('- Bu bir **sablona dayali yeniden cagirma**dir, sifirdan (de novo) birlestirme '
      'DEGILDIR. Sutun bazli ikame hatalarini duzeltir; sablonda olmayan buyuk '
      'yapisal farklari (uzun ekleme/silme) **bulamaz**.')
    A('- Capalar tekrar etmeyen %d-mer\'lerdir; tekrarli bolgelerde kapsama duser ve '
      'o sutunlar `N` olur. `N_yuzde` sutunu bunu gosterir.' % K)
    A('- `yontemler_ayrildi` sutunu iki yontemin celistigi sutun sayisidir. Bu sayi '
      'yuksekse o kutunun konsensusu **guvenilmez**; once ona bakin.')
    A('- Yeni konsensuslar **panele otomatik islenmedi**. Kullanilmadan once '
      '`eski_ile_farkli` ve `N_yuzde` sutunlari gozden gecirilmelidir.')
    A('')
    with open(md, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(L))
    return [md, tsv]
