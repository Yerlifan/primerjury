# -*- coding: utf-8 -*-
"""HEPSINI SIRAYLA KOS - P -> K -> D -> I, bagimlilik sirasiyla, tek tusla.

Her asama bitince duruma yazilir; kesilirse BITMIS asamalar atlanir ve
kaldigi yerden devam eder. Sonunda tek bir birlesik ozet uretilir:
son siparis listesi + kurtarilan hedefler + celiskiler + kimlik tablosu.

Gece baslatilip sabah TEK DOSYAYA bakilmak icin tasarlandi:
    TUM_KOSU_SONUC/00_BIRLESIK_OZET.md
"""

# -------------------------------------------------------------------------
# run_all_stages.py — dort olcum asamasini bagimlilik sirasiyla (P -> K -> D -> I)
# kosar ve sabah tek dosyadan okunacak birlesik ozeti uretir.
#
# GİRDİ  : ASAMALAR tablosundaki dort betik (protocol/single_protocol_measure.py,
#          verification/recovery_round.py, verification/specificity_round.py,
#          verification/identity_verification.py) ve onlarin urettigi TSV dosyalari;
#          ayrica TUM_KOSU_SONUC/durum.json (onceki kosunun nerede kaldigi).
# ÇIKTI  : TUM_KOSU_SONUC/00_BIRLESIK_OZET.md (sabah okunacak tek dosya),
#          TUM_KOSU_SONUC/durum.json, TUM_KOSU_SONUC/kosu_gunlugu.txt.
# ÇAĞRAN : verification/full_chain.py -> T tusu
#          (bat icinde: wsl -e python3 "verification/run_all_stages.py" --kok .)
#
# SIRA ZORUNLUDUR, KEYFI DEGIL: K girdisini P'nin panel tablosundan, D girdisini
# K'nin kurtardigi ciftlerden alir. Bir asama satir uretmezse sonraki asamanin
# girdisi bostur; zincir orada BILEREK durdurulur (bkz. calistir icindeki cikti
# denetimi). Bos girdiyle devam etmek cokme uretmez, gece boyunca anlamsiz ama
# inandirici bir ozet uretir.
# -------------------------------------------------------------------------
import os, sys, csv, json, time, subprocess, argparse

VERSIYON = '1.0 (2026-08-03)'

ASAMALAR = [
    ('P', u'TEK PROTOKOL - panelin tamami tek kuralla olculur',
     ['protocol/single_protocol_measure.py'], 'TEK_PROTOKOL_SONUC/panel_tek_protokol.tsv',
     u'1-2 saat'),
    ('K', u'verification - esik alti satirlar dort yolla kurtarilir',
     ['verification/recovery_round.py'], 'KURTARMA_SONUC/kurtarma_satirlari.tsv',
     u'1-3 saat'),
    ('D', u'DOGRULAMA - kurtarilan ciftler uc katmanla sinanir',
     ['verification/specificity_round.py'], 'DOGRULAMA_SONUC/dogrulama_uc_sutun.tsv',
     u'1-4 saat'),
    ('I', u'KIMLIK DOGRULAMA - rapora giren iddialar bagimsiz sinanir',
     ['verification/identity_verification.py'], 'KIMLIK_SONUC/kimlik_iddialari.tsv',
     u'3-7 saat'),
]



# --- TOPLANTI TALEPLERININ TAM LISTESI ---------------------------------
# Kaynak: "6 Karar Durumu" sayfasi. O sayfada 29 satir var ama hepsi TALEP
# DEGIL: Karar 1-4 = 21 TOPLANTI TALEBI, Karar 5 = 8 OLCUMDEN TURETILEN hedef
# (toplantida istenmedi, olcum sirasinda ortaya cikti). Bu tablo 21 talebi
# tasir ve her birinin hangi panel satir(lar)iyla karsilandigini soyler.
TOPLANTI_TALEPLERI = [
 ('Karar 1', 'Methanosarcina mazei',              ['mazei']),
 ('Karar 1', 'Methanothrix soehngenii',           ['soehngenii', 'Methanothrix_cinsi']),
 ('Karar 1', 'Methanosarcina barkeri',            ['Methanosarcina_cinsi']),
 ('Karar 1', 'Podospora pseudopauciseta',         ['Podospora']),
 ('Karar 1', 'Dictyostelium discoideum (44689)',  ['Dictyostelium', '44689']),
 ('Karar 1', 'Trichoderma asperellum (101201)',   ['Petriella', 'Microasca']),
 ('Karar 2', 'Bacteroides',                       ['Bacteroidales_kumesi']),
 ('Karar 2', 'Alistipes',                         ['Bacteroidales_kumesi']),
 ('Karar 2', 'Proteiniphilum',                    ['Proteiniphilum']),
 ('Karar 2', 'Petrimonas',                        ['Petrimonas']),
 ('Karar 3', 'Hidrojenotrofik metanojenler',      ['Metanomikrobiyales']),
 ('Karar 3', 'Asetoklastik metanojenler',         ['Asetoklastik']),
 ('Karar 3', 'Metilotrofik metanojen',            ['Metilotrofik']),
 ('Karar 3', 'Nitrosocosmicus AOA',               ['Nitrosocosmicus']),
 ('Karar 3', 'Trichoderma cinsi',                 ['Petriella', 'Microasca']),
 ('Karar 3', 'Sakarolitik bakteriler',            ['Sphaerochaeta']),
 ('Karar 3', 'Proteolitik / sintrofik bakteriler',['Synergistaceae', 'Cloacimonas', 'Proteiniphilum']),
 ('Karar 4', 'Arke universal',                    ['Arke_universal']),
 ('Karar 4', 'Bakteri universal',                 ['Bakteri_universal']),
 ('Karar 4', 'Mantar universal',                  ['Mantar_universal']),
 ('Karar 4', 'Universal metanojen',               ['Metanojen_universal']),
]


def sure_metni(sn):
    sn = int(sn)
    return ('%d saniye' % sn) if sn < 90 else ('%d dakika' % round(sn / 60.0)) \
        if sn < 5400 else ('%.1f saat' % (sn / 3600.0))


def tsv_oku(yol):
    if not os.path.exists(yol):
        return []
    with open(yol, encoding='utf-8') as fh:
        # K-2 ikinci savunma hatti: BOS satir da atlanmali, yoksa
        # DictReader bos satiri baslik sanip fieldnames=[] yapiyor.
        return list(csv.DictReader(
            (s for s in fh if s.strip() and not s.startswith('#')), delimiter='\t'))


# ---------------------------------------------------------------------------
# Dort asamayi sirayla kosar. Kesintiye dayaniklidir: her asama bitince durumu
# durum.json icine yazar, tekrar kosuldugunda "bitti" damgali asamalari atlar.
#
# NEDEN "cikis kodu 0" YETMEZ: bir asama sifir donup hicbir cikti uretmemis
# olabilir. Yalniz cikis koduna bakilsaydi o asamaya "bitti" damgasi vurulur ve
# bir sonraki gece sessizce atlanirdi - hata kalicilasirdi. Bu yuzden her
# asamadan sonra beklenen TSV dosyasi ACILIP veri satiri sayilir; dosya yoksa ya
# da bossa durum "CIKTI YOK" / "CIKTI BOS" yazilir ve zincir durdurulur.
# ---------------------------------------------------------------------------
def calistir(kok, ncbi, yeniden, atla):
    CIKTI = os.path.join(kok, 'TUM_KOSU_SONUC')
    os.makedirs(CIKTI, exist_ok=True)
    dyol = os.path.join(CIKTI, 'durum.json')
    durum = {}
    if os.path.exists(dyol) and not yeniden:
        try:
            durum = json.load(open(dyol, encoding='utf-8'))
        except Exception:
            durum = {}
    g = open(os.path.join(CIKTI, 'kosu_gunlugu.txt'), 'a', encoding='utf-8')

    def yaz(s=''):
        print(s, flush=True); g.write(s + '\n'); g.flush()

    yaz('=' * 78)
    yaz(u'  RUN ALL STAGES IN ORDER   version %s   %s' % (VERSIYON, time.strftime('%Y-%m-%d %H:%M')))
    yaz('=' * 78)
    yaz(u'  TOTAL ESTIMATED TIME: 6-16 hours. Can be left running overnight.')
    yaz(u'  State is saved after every stage; if interrupted, FINISHED stages are skipped.')
    yaz('')
    for k, ad, _, _, sure in ASAMALAR:
        d = durum.get(k, {}).get('durum', 'bekliyor')
        yaz(u'   %s  %-58s %-10s [%s]' % (k, ad[:58], sure, d))
    yaz('')

    t0 = time.time()
    for k, ad, komut, ciktilar, sure in ASAMALAR:
        if k in atla:
            yaz(u'--- STAGE %s SKIPPED (on request) ---' % k); continue
        if durum.get(k, {}).get('durum') == 'bitti' and not yeniden:
            yaz(u'--- STAGE %s already finished, skipping ---' % k); continue
        betik = os.path.join(kok, komut[0])
        if not os.path.exists(betik):
            yaz(u'--- STAGE %s SKIPPED: script missing (%s) ---' % (k, komut[0]))
            durum[k] = dict(durum='betik yok'); continue
        yaz('')
        yaz('=' * 78)
        yaz(u'  STAGE %s STARTING: %s   (estimated %s)' % (k, ad, sure))
        yaz('=' * 78)
        arg = [sys.executable, betik, '--kok', kok]
        if k == 'D':
            arg += (['--ncbi', ncbi] if ncbi != 'yok' else ['--yalniz-yerel'])
        ta = time.time()
        rc = subprocess.call(arg)
        # O-14: rc==0 YETMEZ. Asama sifir donup hicbir cikti uretmemis olabilir;
        # o zaman 'bitti' damgasi vurulursa sonraki gece sessizce atlanir.
        bekl = os.path.join(kok, ciktilar)
        satir = len(tsv_oku(bekl)) if os.path.exists(bekl) else -1
        if rc != 0:
            dd = 'hata (%d)' % rc
        elif satir < 0:
            dd = 'CIKTI YOK (%s uretilmedi)' % ciktilar
        elif satir == 0:
            dd = 'CIKTI BOS (%s icinde 0 veri satiri)' % ciktilar
        else:
            dd = 'bitti'
        durum[k] = dict(durum=dd, sure=round(time.time() - ta, 1),
                        satir=satir, zaman=time.strftime('%Y-%m-%d %H:%M'))
        json.dump(durum, open(dyol, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
        yaz(u'  STAGE %s: %s  (%s, %d rows)'
            % (k, durum[k]['durum'], sure_metni(time.time() - ta), max(satir, 0)))
        if dd != 'bitti':
            yaz('')
            yaz(u'  ' + '!' * 70)
            yaz(u'  STOPPED. Stage %s did not produce the expected output: %s' % (k, dd))
            yaz(u'  The following stages would have used that file as input, and continuing')
            yaz(u'  with empty input produces a MEANINGLESS but CONVINCING summary overnight.')
            yaz(u'  The chain was stopped here deliberately.')
            yaz(u'  What to do: read the output of stage %s above, fix the cause,')
            yaz(u'  then run the SAME command again; finished stages will be skipped.')
            yaz(u'  ' + '!' * 70)
            ozet(kok, CIKTI, durum, yaz)
            g.close()
            return 3

    yaz('')
    yaz(u'All stages finished (%s). Writing the combined summary...' % sure_metni(time.time() - t0))
    ozet(kok, CIKTI, durum, yaz)
    g.close()
    return 0


# ---------------------------------------------------------------------------
# Dort asamanin ciktisini TEK markdown dosyasinda birlestirir. Olcum YAPMAZ,
# yalniz mevcut TSV dosyalarini okur; bu yuzden zincir yarida kesilse bile
# cagrilabilir ve elde ne varsa onu gosterir.
#
# Siparis listesi burada iki suzgecten gecer: P asamasinin esigi gecirdigi
# ciftlerden, D asamasinin CELISKILI ya da RISKLI isaretledikleri DUSULUR.
# "EKSIK" isaretliler sayidan dusulmez ama ayrica yazilir - dogrulamasi
# tamamlanmamis bir cift "temiz" ile ayni sey degildir.
# ---------------------------------------------------------------------------
def ozet(kok, CIKTI, durum, yaz):
    d_bitti = durum.get('D', {}).get('durum') == 'bitti'
    k_bitti = durum.get('K', {}).get('durum') == 'bitti'
    """Dort asamanin ciktisini TEK dosyada birlestirir."""
    P = tsv_oku(os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'panel_tek_protokol.tsv'))
    SIP = tsv_oku(os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'SIPARIS_LISTESI.tsv'))
    K = tsv_oku(os.path.join(kok, 'KURTARMA_SONUC', 'kurtarma_satirlari.tsv'))
    D = tsv_oku(os.path.join(kok, 'DOGRULAMA_SONUC', 'dogrulama_uc_sutun.tsv'))
    I = tsv_oku(os.path.join(kok, 'KIMLIK_SONUC', 'kimlik_iddialari.tsv'))

    yol = os.path.join(CIKTI, '00_BIRLESIK_OZET.md')
    with open(yol, 'w', encoding='utf-8') as fh:
        fh.write(u'# Combined summary, four stages in one file\n\n')
        fh.write(u'Generated: %s, version %s\n\n' % (time.strftime('%Y-%m-%d %H:%M'), VERSIYON))

        fh.write(u'## Stage status\n\n| stage | status | time |\n|---|---|---|\n')
        for k, ad, _, _, _ in ASAMALAR:
            d = durum.get(k, {})
            fh.write(u'| **%s** %s | %s | %s |\n'
                     % (k, ad, d.get('durum', u'not run'),
                        sure_metni(d.get('sure', 0)) if d.get('sure') else '-'))

        # --- 1) SIPARIS LISTESI ---
        fh.write(u'\n---\n\n## 1. Final order list\n\n')
        gecen = [r for r in SIP if (r.get('durum') or '').startswith('ESIK USTU')]
        kalan = [r for r in SIP if (r.get('durum') or '').startswith('ESIK ALTI')]
        olcx = [r for r in SIP if (r.get('durum') or '').startswith('OLCULEMEDI')]
        # dogrulamadan gecmeyenleri DUS
        celiskili = {r['hedef'] for r in D if r.get('KARAR') == 'CELISKILI'}
        riskli = {r['hedef'] for r in D if (r.get('KARAR') or '').startswith('RISKLI')}
        eksik = {r['hedef'] for r in D if (r.get('KARAR') or '').startswith('EKSIK')}
        # D-6: 'INCELEME' = gevsek yerel olcutte vurus var, MFEprimer temiz.
        # Celiski degil ama TEMIZ de degil - 3' son iki baz sinanana kadar
        # siparise GIRMEZ. Dusulenler ayri yazilir ki sayi kaybolmasin.
        inceleme = {r['hedef'] for r in D
                    if (r.get('KARAR') or '').startswith('INCELEME')}
        temiz = [r for r in gecen if r['hedef'] not in celiskili
                 and r['hedef'] not in riskli and r['hedef'] not in inceleme]
        # --- 2026-08-06: UC SAYI AYRI AYRI. SIP satirlari artik SINIF sutunu
        # tasiyor (KESIN / EVRENSEL / KOSULLU / ONERILMEZ). Esigi gecemeyen satir
        # listeden SILINMEZ; ozet de bunu boyle yazmali, yoksa siparis dosyasi ile
        # bu ozet farkli sey soyler - bu projede daha once tam olarak bu oldu.
        _sn = lambda r: (r.get('SINIF') or '').strip().upper()
        kesin = [r for r in SIP if _sn(r) in ('KESIN', 'EVRENSEL')]
        kosullu = [r for r in SIP if _sn(r) == 'KOSULLU']
        onerilmez = [r for r in SIP if _sn(r) == 'ONERILMEZ']
        kesin_temiz = [r for r in kesin
                       if r['hedef'] not in celiskili and r['hedef'] not in riskli
                       and r['hedef'] not in inceleme]
        if kesin or kosullu or onerilmez:
            _dus = len(kesin) - len(kesin_temiz)
            fh.write(u'- **CERTAIN: %d pairs**  (dCq >= 3.0, or universal/coverage)%s\n'
                     % (len(kesin),
                        (u' - **%d of these are CONTRADICTORY/RISKY/NEEDS REVIEW** in verification and do not go into the order' % _dus) if _dus else u''))
            fh.write(u'- **ORDERABLE: %d pairs = %d oligos**\n'
                     % (len(kesin_temiz), 2 * len(kesin_temiz)))
            fh.write(u'- **CONDITIONAL: %d pairs**  (dCq 2.0-3.0 - orderable BUT laboratory validation is REQUIRED)\n' % len(kosullu))
            fh.write(u'- **NOT ADVISED: %d pairs**  (dCq < 2.0 - stays on the list, with its reason written out)\n'
                     % len(onerilmez))
            fh.write(u'- A row that fails the threshold is NEVER DELETED SILENTLY; the decision is yours. Every row\'s laboratory note is in SIPARIS_LISTESI.tsv.\n')
            fh.write(u'\n')
        fh.write(u'- Pairs passing the discrimination ratio (tool count): **%d**\n' % len(gecen))
        fh.write(u'- Of those, marked CONTRADICTORY or RISKY during verification: **%d**\n'
                 % (len(gecen) - len(temiz)))
        if eksik:
            fh.write(u'- Verification NOT COMPLETE (NCBI missing): %d. These ARE included in the number above, but they cannot be ordered before the NCBI step finishes' % len(eksik))
        # ESKI/YENI KIMLIK yan yana: ">>" isareti toplantida konusulan ad ile
        # olculen kimligin AYNI OLMADIGINI gosterir.
        fh.write(u'\n| # | CLASS | target | Kraken label | Measured identity | Level | product bp | discrimination mm<=1 | dCq | verification |\n|---|---|---|---|---|---|---|---|---|---|\n')
        _sirali = (kesin + kosullu + onerilmez) if (kesin or kosullu or onerilmez) else gecen
        for i, r in enumerate(_sirali, 1):
            dg = ('CELISKILI' if r['hedef'] in celiskili else
                  'RISKLI' if r['hedef'] in riskli else
                  'eksik' if r['hedef'] in eksik else
                  ('temiz' if d_bitti else 'DOGRULANMADI'))
            _f = (r.get('ad_farkli_mi') or '')
            _im = u'**>>** ' if _f.startswith('>>') else u''
            # '|' MARKDOWN SUTUN AYIRICISIDIR. Kimlik dizeleri birden fazla adi
            # ' | ' ile birlestirdigi icin tabloyu kaydiriyordu; ' / ' yapilir.
            _t = lambda v, n: ((v or '-').replace('|', '/'))[:n]
            fh.write(u'| %d | %s%s | `%s` | %s | %s | %s | %s | %s | %s | %s |\n'
                     % (i, _im, _sn(r) or '-', r['hedef'],
                        _t(r.get('kraken_etiketi'), 44),
                        _t(r.get('olculen_kimlik'), 44),
                        _t(r.get('savunulabilir_duzey'), 28),
                        r.get('urun_bp', ''), r.get('ayrim_mm1', ''),
                        r.get('dCq_karsiligi', '-'), dg))
        _fs = sum(1 for r in _sirali if (r.get('ad_farkli_mi') or '').startswith('>>'))
        if _fs:
            fh.write(u'\nOn %d rows marked **>>** the name discussed in the meeting is NOT THE SAME as the measured identity. Example: `Microascaceae' % _fs)
        if kalan or olcx:
            fh.write(u'\n**Failing the discrimination ratio:** %d; ratio not measurable (universal): %d. Both appear in the CLASS column above'
                     % (len(kalan), len(olcx)))

        # --- 2) KURTARILANLAR ---
        fh.write(u'\n---\n\n## 2. Kurtarma turu\n\n')
        if not K:
            fh.write(u'Kurtarma turu kosulmadi ya da cikti uretmedi.\n')
        else:
            kg = [r for r in K if (r.get('esigi_gecti_mi') or '').startswith('EVET')]
            kt = [r for r in K if r.get('esigi_gecti_mi') == 'DUSENLERE TASINDI']
            fh.write(u'- Kurtarilan: **%d**  ·  Dusenlere tasinan: **%d**  ·  '
                     u'Kurtarilamayan: **%d**\n\n' % (len(kg), len(kt), len(K) - len(kg) - len(kt)))
            fh.write(u'| target | before | route | after | result |\n|---|---|---|---|---|\n')
            for r in K:
                fh.write(u'| %s | %s | %s | %s | **%s** |\n'
                         % (r['hedef'], r.get('eski_deger', ''), (r.get('denenen_yol') or '')[:46],
                            (r.get('yeni_deger') or '-')[:52], r.get('esigi_gecti_mi', '')))
            kk = [r for r in K if not (r.get('esigi_gecti_mi') or '').startswith('EVET')
                  and r.get('sebep')]
            if kk:
                fh.write(u'\n**Kurtarilamayanlarin sebebi**\n\n')
                for r in kk:
                    fh.write(u'- **%s** — %s\n' % (r['hedef'], r['sebep']))

        # --- 3) CELISKILER ---
        fh.write(u'\n---\n\n## 3. Contradictions (specificity round)\n\n')
        cel = [r for r in D if r.get('KARAR') == 'CELISKILI']
        if not D:
            fh.write(u'Dogrulama turu kosulmadi ya da cikti uretmedi.\n')
        elif not cel:
            fh.write(u'Uc katman hicbir satirda ayrilmadi. **Bu "her sey temiz" demek '
                     u'DEGILDIR**: NCBI ya da yerel katmani kosulmamis satirlar "EKSIK" '
                     u'sayilir ve celiski uretmez.\n')
        else:
            fh.write(u'**%d satir celiskili — hicbiri siparis edilemez.**\n\n' % len(cel))
            fh.write(u'| hedef | 1 numune | 2 yerel DB | 3 NCBI |\n|---|---|---|---|\n')
            for r in cel:
                fh.write(u'| %s | %s | %s (%s) | %s (%s) |\n'
                         % (r['hedef'], r.get('1_NUMUNE', ''), r.get('2_YEREL_DB', ''),
                            r.get('2_hedef_disi_urun', ''), r.get('3_NCBI', ''),
                            r.get('3_hedef_disi_urun', '')))
            fh.write(u'\nAyrinti ve ne yapilmasi gerektigi: `DOGRULAMA_SONUC/CELISKILER.md`\n')

        # --- 4) KIMLIK ---
        fh.write(u'\n---\n\n## 4. Identity claims\n\n')
        if not I:
            fh.write(u'Kimlik dogrulama turu kosulmadi ya da cikti uretmedi.\n')
        else:
            say = {}
            for r in I:
                say[r['HUKUM']] = say.get(r['HUKUM'], 0) + 1
            fh.write('  ·  '.join('**%s**: %d' % kv for kv in say.items()) + u'\n\n')
            fh.write(u'| # | claim | verdict | agreeing databases |\n|---|---|---|---|\n')
            for r in I:
                fh.write(u'| %s | %s | **%s** | %s |\n'
                         % (r['no'], r['iddia'][:78], r['HUKUM'], r.get('sonuc_veren_vtb', r.get('uyusan_vtb_sayisi', ''))))
            duz = [r for r in I if r['HUKUM'] == 'DUZELTILMELI']
            if duz:
                fh.write(u'\n### Rapora gonderilecek duzeltmeler\n\n')
                for r in duz:
                    fh.write(u'**%s. %s**\n\n' % (r['no'], r['iddia']))
                    fh.write(u'> %s\n\n' % r.get('DOGRU_IFADE (duzeltilmeliyse)', ''))
                    fh.write(u'Kanit: %s\n\n' % r.get('kanit', ''))


        # --- 5) TOPLANTI TALEPLERININ TAMAMI ---
        fh.write(u'\n---\n\n## 5. ALL requested targets\n\n')
        fh.write(u'The "6 Decision Status" sheet has 29 rows, but not all of them are requests: **decisions 1-4 are the 21 requested targets**, while decision 5 covers 8 targets derived from measurement (they were never requested). The table below shows ALL 21 requests, including those with no row in the panel.\n\n')
        p_ad = [r.get('hedef', '') for r in P]
        k_ad = {r.get('hedef', ''): r for r in K}
        sip_ad = {r.get('hedef', ''): r for r in SIP}
        fh.write(u'| decision | requested target | panel row(s) | measured? | result |\n|---|---|---|---|---|\n')
        panelsiz = 0
        for kar, ad, anah in TOPLANTI_TALEPLERI:
            bul = sorted({x for x in p_ad if any(a.lower() in x.lower() for a in anah)})
            kb = [k_ad[x] for x in k_ad if any(a.lower() in x.lower() for a in anah)]
            if not bul and not kb:
                panelsiz += 1
                fh.write(u'| %s | %s | **YOK** | hayir | **panelde satiri yok - '
                         u'kurtarma turunda tasarim denemesi kapsaminda** |\n'
                         % (kar, ad))
                continue
            if not bul and kb:
                fh.write(u'| %s | %s | *(panelsiz, K turunda denendi)* | evet | %s |\n'
                         % (kar, ad, kb[0].get('esigi_gecti_mi', '-')))
                continue
            durumlar = []
            for x in bul:
                sr = sip_ad.get(x)
                if sr:
                    durumlar.append(u'%s: %s' % (x, sr.get('durum', '-')))
                elif x in k_ad:
                    durumlar.append(u'%s: %s' % (x, k_ad[x].get('esigi_gecti_mi', '-')))
                else:
                    durumlar.append(u'%s: olculdu' % x)
            fh.write(u'| %s | %s | %s | evet | %s |\n'
                     % (kar, ad, ', '.join(bul), '; '.join(durumlar)[:120]))
        fh.write(u'\n**Requests with no row in the panel: %d.** ' % panelsiz)
        if panelsiz:
            fh.write(u'Bu talepler kurtarma turunun `PANELSIZ_TALEPLER` tablosunda '
                     u'tanimlidir ve kutunun KENDI konsensusundan tasarim denemesi '
                     u'yapilir; "hic denenmedi" diye kalmazlar.\n')
        else:
            fh.write(u'Zincir toplanti listesinin TAMAMINI kapsiyor.\n')

        fh.write(u'\n---\n\n## Where to look\n\n| Question | File |\n|---|---|\n| What should I order? | `TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv` |\n| Why did this row fail? | `KURTARMA_SONUC/KURTARMA_RAPORU.md` |\n| Which pair is suspect? | `DOGRULAMA_SONUC/CELISKILER.md` |\n| What do I write in the report? | `KIMLIK_SONUC/KIMLIK_DOGRULAMA_RAPORU.md` |\n')
    yaz(u'  written: %s' % yol)
    yaz('')
    yaz(u'  LOOK AT THIS FIRST: %s' % yol)


# Komut satiri: --kok proje klasoru, --ncbi D asamasinin NCBI kipi, --yeniden
# bitmis asamalari da tekrar kosar, --atla belirli asamalari atlar (orn. "DI").

# --- CLI value normalisation ------------------------------------------------
# English option values are accepted alongside the original Turkish ones and
# mapped back here. The internal values are unchanged on purpose: they are
# compared in dozens of places and, in some cases, written to output files.
# Translating the interface must not translate the data.
_DEGER = {'auto': 'oto', 'manual': 'elle', 'none': 'yok', 'quick': 'hizli',
          'full': 'tam', 'member': 'uye', 'competitor': 'rakip',
          'exclude': 'disla'}


def _ing_deger(a):
    for _ad in ('nt', 'literatur', 'ncbi', 'karisik', 'moduller', 'mod'):
        _v = getattr(a, _ad, None)
        if isinstance(_v, str) and _v in _DEGER:
            setattr(a, _ad, _DEGER[_v])
    return a

def main():
    p = argparse.ArgumentParser(description='P -> K -> D -> I sirayla')
    p.add_argument('--root', '--kok', dest='kok', default='.')
    p.add_argument('--ncbi', choices=['auto', 'manual', 'none', 'oto', 'elle', 'yok'], default='elle')
    p.add_argument('--rerun', '--yeniden', dest='yeniden', action='store_true', help='re-run finished stages as well')
    p.add_argument('--skip', '--atla', dest='atla', default='', help='stages to skip, e.g. I or DI')
    a = p.parse_args()
    a = _ing_deger(a)
    kok = os.path.abspath(a.kok)
    if not os.path.isdir(os.path.join(kok, 'screening')):
        sys.exit('HATA: %s icinde screening yok.' % kok)
    return calistir(kok, a.ncbi, a.yeniden, set(a.atla.upper()))


if __name__ == '__main__':
    sys.exit(main() or 0)
