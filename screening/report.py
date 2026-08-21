# -*- coding: utf-8 -*-
"""Cikti: makine icin TSV, insan icin okunabilir Markdown raporu.

Raporun cevaplamasi gereken soru:
  "Bu hedef icin hangi parametre ayarinda cozum var, yok mu, varsa bedeli ne?"
"""
# ---------------------------------------------------------------------------
# report.py — arama sonuclarini makine icin TSV, insan icin Markdown olarak
#            yazar; her hedef icin tek paragraflik bir karar cumlesi uretir.
#
# GIRDI  : kontrol.hepsi()'nin diskten okudugu kontrol/hedef_*.json sonuclari
#          (her biri __main__.hedefi_isle'nin dondurdugu sozluk: taban olcumu,
#          izgara tablosu, aday listesi, kuresel tarama sonucu);
#          hedefler.panel_oku()'nun verdigi panel satirlari ve panel yolu.
# CIKTI  : KAPSAMLI_ARAMA_SONUC/adaylar.tsv, parametre_izgarasi.tsv ve
#          KAPSAMLI_ARAMA_RAPORU.md. uret() bu uc yolun listesini dondurur.
# CAGRAN : __main__.aramayi_kos icinden - her hedef bittikten sonra bir kez
#          (sessiz) ve kosunun sonunda bir kez daha. Yani verification/full_chain.py
#          tuslari 1, 2, 3, 7 ve 9 (icindeki 7. asama) uzerinden calisir.
# ---------------------------------------------------------------------------
import os, time, csv
from . import config as C


def aday_ozet(c):
    d = dict(
        F=c['F'], R=c['R'], urun=c['urun'],
        F_uz=c['mF']['uz'], F_tm=c['mF']['tm'], F_gc=c['mF']['gc'],
        F_uc=c['mF']['uc'], F_son5=c['mF']['son5'],
        F_hp_tm=c['mF']['hp_tm'], F_hd_tm=c['mF']['hd_tm'],
        F_hp_dg=c['mF']['hp_dg'], F_hd_dg=c['mF']['hd_dg'],
        R_uz=c['mR']['uz'], R_tm=c['mR']['tm'], R_gc=c['mR']['gc'],
        R_uc=c['mR']['uc'], R_son5=c['mR']['son5'],
        R_hp_tm=c['mR']['hp_tm'], R_hd_tm=c['mR']['hd_tm'],
        R_hp_dg=c['mR']['hp_dg'], R_hd_dg=c['mR']['hd_dg'],
        arms=c.get('arms', ''),
        izgara_hucresi=c.get('izgara_hucresi', ''), sikilik=c.get('sikilik', ''),
    )
    cf = c.get('cift', {})
    d.update({'cift_' + k: v for k, v in cf.items()})
    nm = c.get('numune', {})
    d['numune_olcut'] = nm.get('olcut', '')
    n3 = c.get('numune_mm3') or {}
    d['numune_olcut_2'] = n3.get('olcut', '')
    d['numune_kat_enkotu_mm3'] = n3.get('kat_enkotu', '')
    d['numune_kat_havuz_mm3'] = n3.get('kat_havuz', '')
    d['numune_uye_kapsam_mm3'] = n3.get('uye_kapsam_pay', '')
    for k in ('uye_alt', 'uye_min', 'uye_max', 'uye_kutu_sayisi', 'havuz',
              'havuz_ust', 'kat_havuz', 'kat_enkotu', 'enkotu_kutu',
              'uye_kapsam', 'uye_kapsam_pay', 'uye_alt_kapsayan',
              'kat_havuz_kapsayan', 'kat_enkotu_kapsayan'):
        d['numune_' + k] = nm.get(k, '')
    d['numune_urun_boylari'] = ';'.join('%s:%s' % (x[0], x[1]) for x in nm.get('urun_boylari', []))
    ru = c.get('ref_uye'); rr = c.get('ref_rakip')
    d['ref_uye'] = '%d/%d' % (ru['veren'], ru['toplam']) if ru else ''
    d['ref_rakip'] = '%d/%d' % (rr['veren'], rr['toplam']) if rr else ''
    kg = c.get('kuresel')
    if kg and 'urun' in kg:
        d['kuresel_urun'] = kg['urun']
        d['kuresel_boy'] = ';'.join('%s:%s' % (x[0], x[1]) for x in
                                    sorted(kg.get('boy', {}).items(), key=lambda y: -y[1])[:5])
    else:
        d['kuresel_urun'] = ''
        d['kuresel_boy'] = ''
    return d


SUTUNLAR = ['hedef', 'sira', 'numune_olcut', 'F', 'R', 'urun', 'cift_urun_sinifi', 'arms',
            'izgara_hucresi', 'sikilik',
            'numune_uye_alt', 'numune_uye_min', 'numune_uye_max', 'numune_uye_kutu_sayisi',
            'numune_uye_kapsam_pay', 'numune_havuz', 'numune_kat_havuz',
            'numune_kat_enkotu', 'numune_kat_havuz_kapsayan',
            'numune_kat_enkotu_kapsayan', 'numune_enkotu_kutu',
            'numune_olcut_2', 'numune_kat_enkotu_mm3', 'numune_kat_havuz_mm3',
            'numune_uye_kapsam_mm3',
            'ref_uye', 'ref_rakip', 'kuresel_urun', 'kuresel_boy',
            'cift_dTm', 'cift_het_tm', 'cift_het_dg', 'cift_uc_dg',
            'cift_Ta_kural', 'cift_Ta60_marj', 'cift_Ta60_uygun',
            'F_uz', 'F_tm', 'F_gc', 'F_uc', 'F_son5', 'F_hp_tm', 'F_hp_dg', 'F_hd_tm', 'F_hd_dg',
            'R_uz', 'R_tm', 'R_gc', 'R_uc', 'R_son5', 'R_hp_tm', 'R_hp_dg', 'R_hd_tm', 'R_hd_dg',
            'numune_urun_boylari']


def uret(sonuclar, panel, panel_yolu, yaz):
    os.makedirs(C.CIKTI, exist_ok=True)
    yollar = []
    yollar.append(_adaylar_tsv(sonuclar))
    yollar.append(_izgara_tsv(sonuclar))
    yollar.append(_rapor_md(sonuclar, panel, panel_yolu))
    if yaz:
        yaz('')
    return yollar


def _adaylar_tsv(sonuclar):
    p = os.path.join(C.CIKTI, 'adaylar.tsv')
    with open(p, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(SUTUNLAR)
        for s in sonuclar:
            for i, a in enumerate(s.get('adaylar', []), 1):
                a = dict(a); a['hedef'] = s['hedef']; a['sira'] = i
                w.writerow([a.get(k, '') for k in SUTUNLAR])
    return p


def _izgara_tsv(sonuclar):
    p = os.path.join(C.CIKTI, 'parametre_izgarasi.tsv')
    with open(p, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'GC_alt', 'GC_ust', 'Tm_alt', 'Tm_ust',
                    'urun_alt', 'urun_ust', "3'_uc_GC_sart", 'son5_GC_<=3_sart',
                    'gevseklik_puani', 'hayatta_kalan_aday', 'ornek_cift'])
        for s in sonuclar:
            for x in s.get('izgara', []):
                h = x['hucre']
                w.writerow([s['hedef'], h['gc'][0], h['gc'][1], h['tm'][0], h['tm'][1],
                            h['urun'][0], h['urun'][1],
                            'EVET' if h['uc_gc'] else 'HAYIR',
                            'EVET' if h['son5'] else 'HAYIR',
                            x['sikilik'], x['hayatta'], x['ornek']])
    return p


def _f(v, k=1):
    try:
        return ('%.' + str(k) + 'f') % float(v)
    except Exception:
        return str(v)


def _taban(s):
    t = s.get('panel_olcum') or {}
    return t.get('kat_enkotu'), t.get('kat_havuz')


def _hedef_karari(s):
    """'Cozum var mi, varsa bedeli ne' - tek paragraflik cevap uretir."""
    # Karar UC basamaklidir ve sirasi onemlidir. Once TABAN yazilir: paneldeki
    # mevcut cift ayni motorla, ayni kutularda olculmus degeridir; aday ancak
    # bunu geciyorsa "daha iyi" sayilabilir. Ikinci basamak, adayin mevcut cifti
    # gecip gecmedigidir. Ucuncusu, gecse bile 10x esigini asip asmadigidir -
    # "daha iyi" ile "yeterli" ayni sey degildir, taban 0,19x iken 1,04x'e cikan
    # bir aday daha iyidir ama hala kullanilamaz.
    #
    # BEDEL listesi bilerek karar cumlesine baglanir: bir aday esigi gecse bile
    # urun boyu 250 bp'yi asiyorsa, 60 C'de kosulamiyorsa, kasitli uyumsuzluk
    # (ARMS) gerektiriyorsa ya da ancak gevsetilmis bir izgara hucresinde
    # hayattaysa bu bilgi sayinin yaninda durmalidir. Aksi halde rapor "cozum
    # var" der ve bedeli okunmadan siparise gidilir.
    ad = s.get('adaylar', [])
    if s.get('durum') != 'TAMAMLANDI':
        return s.get('durum', 'bilinmiyor'), ''
    if not ad:
        return 'COZUM YOK', 'Hicbir aday numune olcumunu gecmedi.'
    izg = s.get('izgara', [])
    dolu = [x for x in izg if x['hayatta'] > 0]
    en_siki = dolu[0] if dolu else None

    olculu = [a for a in ad if a.get('numune_kat_enkotu') not in ('', None)]
    olculu.sort(key=lambda a: -(a['numune_kat_enkotu'] or 0))
    if not olculu:
        return 'COZUM YOK', 'Adaylarin hicbirinde olculebilir ayrim yok.'
    iyi = olculu[0]
    arms_li = [a for a in olculu if a.get('arms')]
    duz_li = [a for a in olculu if not a.get('arms')]

    p = []
    tb_kotu, tb_havuz = _taban(s)
    if tb_kotu is not None or tb_havuz is not None:
        p.append('TABAN - paneldeki mevcut cift ayni motorla: ayrim %sx (en kotu kutu) / '
                 '%sx (havuz). Asagidaki sayilar bununla karsilastirilmalidir.'
                 % (_f(tb_kotu), _f(tb_havuz)))
    if en_siki:
        p.append('En siki ayar (%s) altinda %d aday hayatta.' % (en_siki['ad'], en_siki['hayatta']))
    bos_siki = [x for x in izg if x['sikilik'] == 0 and x['hayatta'] == 0]
    if bos_siki:
        p.append('TAM SIKI ayarda (GC 40-60, Tm 58-62, urun 60-150, 3\' uc G/C sart, '
                 'son 5 baz <=3 G/C) COZUM YOK.')
    if duz_li:
        b = duz_li[0]
        p.append('ARMS\'siz en iyi aday %s / %s (%s bp), ayrim %sx (en kotu kutu), '
                 'uye %%%s-%%%s, izgara hucresi: %s.'
                 % (b['F'], b['R'], b['urun'], _f(b['numune_kat_enkotu']),
                    _f(b['numune_uye_min']), _f(b['numune_uye_max']), b['izgara_hucresi']))
    if arms_li:
        b = arms_li[0]
        p.append('ARMS varyantiyla en iyi: %s / %s (%s bp), ayrim %sx  [%s].'
                 % (b['F'], b['R'], b['urun'], _f(b['numune_kat_enkotu']), b['arms']))
    if iyi.get('kuresel_urun') not in ('', None):
        p.append('Kuresel taramada en iyi adayin urun sayisi: %s.' % iyi['kuresel_urun'])
    # bedel
    bedeller = []
    if iyi.get('cift_urun_sinifi', '').startswith('kabul'):
        bedeller.append('urun 150-250 bp -> protokolde 30 sn annealing/extension')
    if iyi.get('cift_urun_sinifi', '').startswith('ONERILMEZ'):
        bedeller.append('urun >250 bp -> QuantiNova icin onerilmez')
    if iyi.get('cift_Ta60_uygun') in (False, 'False'):
        bedeller.append('Ta = min(Tm)-3 kurali ile %s C cikiyor, 60 C hedefinin altinda'
                        % _f(iyi.get('cift_Ta_kural')))
    if iyi.get('arms'):
        bedeller.append('kasitli uyumsuzluk (ARMS) gerekiyor - ayri toplanti maddesi')
    if 'serbest' in str(iyi.get('izgara_hucresi', '')):
        bedeller.append('geometri kuralinin gevsetilmesi gerekiyor: %s' % iyi['izgara_hucresi'])
    if bedeller:
        p.append('BEDELI: ' + '; '.join(bedeller) + '.')
    en = iyi.get('numune_kat_enkotu') or 0
    if tb_kotu and en <= tb_kotu:
        karar = 'COZUM YOK (mevcut cift daha iyi)'
        p.append('SONUC: taranan hicbir aday mevcut cifti GECEMEDI - mevcut cift korunmali.')
    elif en >= 10:
        karar = 'COZUM VAR'
    else:
        karar = 'KISMI COZUM'
    return karar, ' '.join(p)


def _rapor_md(sonuclar, panel, panel_yolu):
    p = os.path.join(C.CIKTI, 'KAPSAMLI_ARAMA_RAPORU.md')
    L = []
    A = L.append
    A('# Kapsamli primer aramasi - rapor')
    A('')
    A('Uretim zamani: %s' % time.strftime('%Y-%m-%d %H:%M:%S'))
    A('')
    A('Kaynak panel: `%s`' % os.path.basename(panel_yolu))
    A('')
    A('## Nasil okunur')
    A('')
    A('Her hedef icin sorulan soru: **hangi parametre ayarinda cozum var, yok mu, '
      'varsa bedeli ne?** Ozet tablodaki `Karar` sutunu bunu tek kelimeyle, '
      'altindaki paragraf ayrintisiyla soyler.')
    A('')
    A('Sabit qPCR kisitlari (QIAGEN Rotor-Gene Q + QuantiNova SYBR Green):')
    A('')
    A('- Amplikon **60-150 bp tercih edilir**; 150-250 kabul edilebilir ama protokolde '
      'annealing/extension **30 sn** gerektirir; **>250 onerilmez**.')
    A('- Rotor-Gene tek kosuda tek dongu programi calistirir: butun panelin **ayni Ta**\'da '
      'kosmasi hedeftir, **60 C oncelikli**. Her aday 60 C\'de de degerlendirildi '
      '(`cift_Ta60_marj`, `cift_Ta60_uygun` sutunlari).')
    A('- SYBR Green oldugu icin **primer-dimer ve hairpin ELEYICI olcuttur**, uyari degil: '
      'hairpin/homodimer/heterodimer Tm >= 45 C ya da dG < -9 kcal/mol olan aday elenir.')
    A('')
    A('> **60 C uyarisi (olcumden cikan yapisal sonuc).** Panelin kurali `Ta = min(Tm) - 3`. '
      'Bu kuralla Ta\'nin 60 C olmasi icin daha dusuk Tm\'li primerin **63 C** olmasi gerekir. '
      'Izgaradaki en genis Tm penceresi 56-64 oldugu icin 60 C\'lik ortak Ta yalnizca '
      '**Tm 56-64 penceresinde ve Tm >= 63** olan adaylarla mumkundur. Adaylarin 60 C '
      'uygunlugu `cift_Ta60_uygun` sutununda ayrica isaretlendi; uymayanlar icin secenek '
      'ya Tm penceresini yukari tasimak ya da `Ta = min(Tm) - 3` kuralini bu panel icin '
      'gevsetmektir. Bu bir **toplanti karari**dir, arac kendiliginden secmez.')
    A('')
    A('> **ARMS hakkinda.** Kasitli uyumsuzluk **dejenere baz degildir**: tek tanimli bir '
      'bazdir, tupte tek oligo kalir, sentez maliyeti artmaz ve "panelde dejenere baz yok" '
      'kaydini bozmaz. Ama sablonla **tam eslesmez**: verimi dusurur ve ayri bir toplanti '
      'maddesidir. Rapor ARMS\'li ve ARMS\'siz en iyi adayi **ayri ayri** verir ki karar '
      'kullanicida kalsin.')
    A('')

    A('> **Taban degerler hakkinda.** Her hedefte paneldeki MEVCUT cift de ayni motorla, '
      'ayni kutularda, ayni olcutle yeniden olculur ve adaylar **o tabanla** '
      'karsilastirilir. Bazi hedeflerde bu taban, panelin yayimladigi sayidan sapar '
      '(farkli okuma derinligi, farkli uye kutu alt kumesi ya da farkli gevseklik ayari '
      'yuzunden). Sapma buyukse hedefin altinda **UYARI** olarak yazilir. '
      'Karsilastirma yine de gecerlidir: aday ve taban **ayni** kosullarda olculur. '
      'Ama panel sayisiyla bu rapordaki sayi **dogrudan karsilastirilmamalidir**.')
    A('')
    A('## Ozet')
    A('')
    A('Sutun `Neden sorunlu` harf kodlari: '
      '**G** geometri ihlali, **K** kosullu/on karar, **A** ayrim ya da kapsam esik alti, '
      '**U** urun boyu qPCR ideali disinda, **C** panelden cikarilmis, '
      '**P** plaka ici jelde ayrilamiyor.')
    A('')
    A('| Hedef | Neden sorunlu | Karar | Mevcut cift (x) | En iyi aday (x) | Urun (bp) | ARMS gerekti mi |')
    A('|---|---|---|---|---|---|---|')
    for s in sonuclar:
        kar, _ = _hedef_karari(s)
        ad = s.get('adaylar', [])
        olculu = sorted([a for a in ad if a.get('numune_kat_enkotu') not in ('', None)],
                        key=lambda a: -(a['numune_kat_enkotu'] or 0))
        b = olculu[0] if olculu else {}
        tk, _th = _taban(s)
        A('| %s | %s | **%s** | %s | %s | %s | %s |' % (
            s['hedef'], s.get('etiketler', ''), kar, _f(tk) if tk else '-',
            _f(b.get('numune_kat_enkotu', '-')), b.get('urun', '-'),
            'EVET' if b.get('arms') else 'hayir'))
    A('')

    for s in sonuclar:
        A('---')
        A('')
        A('## %s' % s['hedef'])
        A('')
        if s.get('durum') != 'TAMAMLANDI':
            A('Durum: **%s**' % s.get('durum'))
            if s.get('hata'):
                A('')
                A('```')
                A(s['hata'][-1200:])
                A('```')
            A('')
            continue
        pn = s.get('panel', {})
        A('**Paneldeki cift:** `%s` / `%s` — %s bp, plaka %s, Ta %s'
          % (pn.get('F'), pn.get('R'), pn.get('urun'), pn.get('plaka'), pn.get('ta')))
        A('')
        t = s.get('panel_olcum') or {}
        if t:
            A('')
            A('**Mevcut ciftin AYNI motorla olculmus degerleri (karsilastirma tabani):** '
              'uye %%%s-%%%s (%s kutu), rakip havuz %s, **ayrim %sx (havuz) / %sx (en kotu kutu: %s)**'
              % (_f(t.get('uye_min')), _f(t.get('uye_max')), t.get('uye_kutu_sayisi'),
                 t.get('havuz'), _f(t.get('kat_havuz')), _f(t.get('kat_enkotu')),
                 t.get('enkotu_kutu')))
        if s.get('uyelik_uyarisi'):
            A('')
            A('> **UYARI - uyelik tanimi kontrol edilmeli.**')
            for u in s['uyelik_uyarisi']:
                A('> ' + u.replace('!! ', ''))
        A('')
        A('**Neden arandi:** ' + '; '.join(s.get('gerekceler', [])))
        A('')
        A('**Omurga:** `%s` (%s bp) — uyelik kaynagi `%s`, uye taxid: %s'
          % (s['omurga']['kutu'], s['omurga']['uzunluk'], s.get('uyelik_kaynagi'),
             ', '.join(s.get('uye_tax', []))))
        A('')
        sy = s.get('sayilar', {})
        A('**Arama boyu:** %s pencere -> %s ileri + %s geri aday -> %s cift '
          '(+%s ARMS varyanti) -> %s numunede olculdu -> %s cift yapisini gecti.'
          % (sy.get('pencere'), sy.get('ileri'), sy.get('geri'), sy.get('cift'),
             sy.get('arms'), sy.get('numune_olculen'), sy.get('cift_yapisi_gecen')))
        A('')
        kar, aciklama = _hedef_karari(s)
        A('### Karar: %s' % kar)
        A('')
        A(aciklama)
        A('')
        A('### Parametre izgarasi — hangi ayar kac aday birakiyor')
        A('')
        A('| GC | Tm | Urun | 3\' uc G/C | son 5 G/C | gevseklik | hayatta kalan |')
        A('|---|---|---|---|---|---|---|')
        for x in s.get('izgara', [])[:24]:
            h = x['hucre']
            A('| %d-%d | %d-%d | %d-%d | %s | %s | %d | **%d** |' % (
                h['gc'][0], h['gc'][1], h['tm'][0], h['tm'][1], h['urun'][0], h['urun'][1],
                'sart' if h['uc_gc'] else 'serbest', '<=3' if h['son5'] else 'serbest',
                x['sikilik'], x['hayatta']))
        A('')
        A('(Tam 144 hucre: `parametre_izgarasi.tsv`)')
        A('')
        A('### En iyi adaylar')
        A('')
        A('| # | Ileri | Geri | bp | uye kapsam | ayrim x (en kotu kutu) | havuz x | uye % | ref uye | ref rakip | kuresel urun | ARMS | izgara hucresi |')
        A('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
        ad = sorted([a for a in s.get('adaylar', [])],
                    key=lambda a: -(a.get('numune_kat_enkotu') or 0))
        for i, a in enumerate(ad[:15], 1):
            A('| %d | `%s` | `%s` | %s | %s | %s | %s | %s-%s | %s | %s | %s | %s | %s |' % (
                i, a['F'], a['R'], a['urun'], a.get('numune_uye_kapsam_pay', ''),
                _f(a.get('numune_kat_enkotu', '')), _f(a.get('numune_kat_havuz', '')),
                _f(a.get('numune_uye_min', '')), _f(a.get('numune_uye_max', '')),
                a.get('ref_uye', ''), a.get('ref_rakip', ''), a.get('kuresel_urun', ''),
                a.get('arms', '') or '-', a.get('izgara_hucresi', '')))
        A('')
        A('Butun sutunlar (Tm, GC, hairpin/dimer dG, 60 C marji, urun boyu dagilimi): `adaylar.tsv`')
        A('')

    A('---')
    A('')
    A('## Yontem ve sinirlar')
    A('')
    A('- Olcum motoru **yeniden yazilmadi**: `engine/ispcr.py` '
      '(`find_sites`/`amplify`), `engine/scanner.py` (`Havuz`) ve '
      '`engine/pair.py` (`urunler`) dogrudan ice aktarildi. '
      'Geometri esikleri `engine/geometry_core.py` ile birebir olacak sekilde '
      'her koşuda sinaniyor (bkz. kendini sinama).')
    A('- **Olcut etiketi her satirda yazilidir** (`numune_olcut` sutunu). Eleme '
      '**<=1 uyumsuzluk** ile yapilir (panelin numune olcutu); en iyi adaylar '
      'ayrica **<=3** ile de olculup `numune_olcut_2` / `*_mm3` sutunlarina '
      'yazilir. Iki olcut ayridir, birbirinin yerine kullanilamaz.')
    A('- Numune olcutu panelin olcutuyle ayni: **uyumsuzluk <=1 + 3\' son 2 baz TAM**. '
      'Kuresel olcut: **toplam <=5 uyumsuzluk**, F ve R ayri. Iki olcut ayridir.')
    A('- Ayrim oranlari **Wilson** ile muhafazakar yonde: uye icin ALT sinir, rakip icin '
      'UST sinir. Farkli hedeflerin oranlari farkli okuma derinliginde olculur, '
      '**dogrudan karsilastirilamaz** (panelin A26 uyarisinin aynisi).')
    A('- Aramanin tamligi: omurga TEK zincirde taranir; cift-zincirli sablonda bir cift '
      '(+ zincirde F, - zincirde R) ile tam tanimlandigi icin ters zincir ayni kumeyi '
      'verir. Konsensuslerin bir kisminin ters yonde saklanmis olmasi kapsami etkilemez.')
    A('- Huni yapisi: pencere -> geometri -> cift -> numune -> referans -> **kuresel** '
      '(en pahali adim en sonda, yalniz diger butun suzgeclerden gecen adaylara).')
    A('- Bu araç **karar vermez**: olcer, bedelini yazar, secimi kullaniciya birakir.')
    A('')
    with open(p, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(L))
    return p
