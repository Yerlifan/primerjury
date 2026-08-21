# -*- coding: utf-8 -*-
"""KAPSAMLI ARAMA - ana akis.

Kullanim (normalde verification/full_chain.py cagirir):
    python3 -m screening --mod tam
    python3 -m screening --mod sorunlu
    python3 -m screening --mod devam
    python3 -m screening --sina          (kendini sina, olcum yapmaz)
"""
# ---------------------------------------------------------------------------
# __main__.py — paketin tek giris noktasi; komut satiri secenegini okur ve
#               ilgili asamayi (arama, panel olcumu, uyelik denetimi, konsensus
#               uretimi, ozet, hepsi) sirasiyla calistirir.
#
# GIRDI  : config.py'deki butun yollar ve sabitler; hedefler.panel_oku()
#          ile panel TSV'si; hedefler.konsensusler() ile kanonik konsensusler;
#          hedefler.kutular() ile "fastq files" altindaki kutu dosyalari;
#          numune.Numune ile kutu basina ham okuma havuzlari; reference.py ile
#          SILVA/UNITE havuzlari; checks.py ile onceki kosunun kontrol
#          noktalari. Argumanlari argparse ile alir (--mod, --hedef, --okuma,
#          --hafif, --yeniden, --tam-derinlik, --sina, --sinama-atla).
# CIKTI  : KAPSAMLI_ARAMA_SONUC/ altina kontrol/hedef_*.json kontrol noktalari;
#          rapor.uret() araciligiyla adaylar.tsv, parametre_izgarasi.tsv ve
#          KAPSAMLI_ARAMA_RAPORU.md. main() surec cikis kodunu dondurur
#          (0 = basarili, 2 = kapi/sinama dusuklugu).
# CAGRAN : verification/full_chain.py icinde "wsl -e python3 -m screening
#          --mod %MOD%" satiri. Tuslar: 1 (--mod tam), 2 (--mod sorunlu),
#          3 (--mod devam), 4 (--mod panel-olc --tam-derinlik), 5 (--mod
#          uyelik), 6 (--mod konsensus), 7 (tek hedef; secime gore sorunlu /
#          panel-olc / uyelik + --hedef), 8 (--sina), 9 (--mod hepsi),
#          S (--mod ozet). P, K, D, T, I, G, E, R, U, H ve W/X/Y/Z tuslari
#          verification/, protocol/, engine/ ve tools/WSL klasorlerindeki
#          ayri betikleri calistirir; bu paketi dogrudan cagirmazlar.
# ---------------------------------------------------------------------------
import os, sys, time, argparse, traceback

from . import config as C
from . import engine_gateway, geometri as G, hedefler as H, uretec as U, numune as N
from . import reference as REF, kuresel_tarama as KT, kontrol, rapor

BASLANGIC = time.time()


def yaz(s=''):
    print(s, flush=True)


def sure(sn):
    sn = int(max(sn, 0))
    if sn < 60:
        return '%d sn' % sn
    if sn < 3600:
        return '%d dk %02d sn' % (sn // 60, sn % 60)
    return '%d sa %02d dk' % (sn // 3600, (sn % 3600) // 60)


def cizgi(ch='-'):
    yaz(ch * 78)


# ---------------------------------------------------------------- tek hedef
# Tek bir panel hedefi icin BUTUN arama hunisini bastan sona kosar. Asamalarin
# sirasi keyfi degil, maliyete gore dizilmistir: once ucuz ve eleyici olan
# (geometri), sonra pahali olan (ham okumalarda in-silico PCR), en sonda en
# pahali olan (yaklasik 500 bin kayitlik kuresel tarama) gelir. Boylece pahali
# adima yalniz onceki suzgeclerden gecen az sayida aday girer.
#
# Asama [0] mevcut panel cifti icin bir TABAN olcer. Bu taban olmadan "aday
# daha iyi" cumlesi olculemez: aday ile mevcut cift AYNI motorda, AYNI
# kutularda, AYNI olcutle olculmedikce karsilastirma anlamsizdir.
def hedefi_isle(satir, baglam, numune, sira, toplam, hafif=False):
    ad = satir['hedef']
    t0 = time.time()
    cizgi('=')
    yaz('[%d/%d] HEDEF: %s' % (sira, toplam, ad))
    yaz('       gerekce: ' + '; '.join(g[1] for g in satir['gerekceler'])[:200])
    yaz('       paneldeki cift: %s / %s  (%d bp)' % (satir['F'], satir['R'], satir['urun_bp']))
    cizgi('=')

    om = baglam['omurga']
    if om is None:
        yaz('  ATLANDI: bu hedefin uye konsensusu bulunamadi (uyelik kaynagi: %s).'
            % baglam['uyelik_kaynagi'])
        return dict(hedef=ad, durum='ATLANDI - uye konsensusu yok',
                    uyelik_kaynagi=baglam['uyelik_kaynagi'])

    yaz('  omurga konsensus : %s (%d bp)' % (om['kutu'], len(om['dizi'])))
    yaz('  uye kutu / rakip : %d / %d      uye konsensus / rakip: %d / %d'
        % (len(baglam['uye_kutu']), len(baglam['rakip_kutu']),
           len(baglam['uye_kons']), len(baglam['rakip_kons'])))

    # ---------------- TABAN: paneldeki mevcut cift AYNI motorla olculur
    yaz('\n  [0] Paneldeki mevcut cift ayni motorla olculuyor (karsilastirma tabani)')
    taban = numune.olc(satir['F'], satir['R'], baglam['uye_kutu'], baglam['rakip_kutu'],
                       lo=C.URUN_IDEAL[0], hi=C.URUN_MUTLAK_UST)
    if taban:
        yaz('      mevcut cift: uye %%%.1f-%%%.1f | kapsam %s kutu (>=%d%%) | havuz %s'
            % (taban['uye_min'], taban['uye_max'], taban['uye_kapsam_pay'],
               int(100 * C.KAPSAM_ESIGI), taban['havuz']))
        yaz('                   ayrim %s x (havuz) / %s x (en kotu kutu)   |   '
            'yalniz kapsanan kutular: %s x / %s x'
            % (taban['kat_havuz'], taban['kat_enkotu'],
               taban['kat_havuz_kapsayan'], taban['kat_enkotu_kapsayan']))
        yaz('      -> Bir adayin "daha iyi" sayilmasi icin BU sayilari gecmesi gerekir.')
        uyari = uyelik_uyarisi(satir, taban)
        if uyari:
            yaz('')
            for satir_u in uyari:
                yaz('      ' + satir_u)
    else:
        yaz('      mevcut cift olculemedi (uye kutusu yok).')

    # ---------------- ASAMA A: pencereler + geometri
    yaz('\n  [A] Pencere uretimi ve geometri olcumu (18-25 bp, her pozisyon, iki yon)')
    ta = time.time()

    def ilerA(n):
        el = time.time() - ta
        print('      ... %d pencere olculdu (%s)' % (n, sure(el)), end='\r', flush=True)

    ad_p = U.aday_primerler(om['dizi'], ilerle=ilerA)
    yaz('      taranan pencere : %d' % ad_p['taranan_pencere'])
    yaz('      degismez kurallari gecen: ileri %d, geri %d   (%s)'
        % (len(ad_p['F']), len(ad_p['R']), sure(time.time() - ta)))

    if not ad_p['F'] or not ad_p['R']:
        return dict(hedef=ad, durum='COZUM YOK - hicbir pencere degismez kurallari gecmedi',
                    pencere=ad_p['taranan_pencere'])

    # ---------------- ASAMA B: cift kurma + izgara
    yaz('\n  [B] Es primer kombinasyonlari (urun %d-%d bp) ve 144 hucreli parametre izgarasi'
        % (C.URUN_IDEAL[0], C.URUN_MUTLAK_UST))
    tb = time.time()

    def ilerB(n):
        print('      ... %d cift sayildi (%s)' % (n, sure(time.time() - tb)), end='\r', flush=True)

    top = U.tara_ve_topla(ad_p, hucre_basina=8, ilerle=ilerB)
    cl = top['temsilciler']
    yaz('      kurala uyan cift : %d  (tamami sayildi, ust sinir yok)   (%s)'
        % (top['toplam'], sure(time.time() - tb)))
    if top['toplam'] == 0:
        return dict(hedef=ad, durum='COZUM YOK - urun boyu penceresinde hic cift yok',
                    pencere=ad_p['taranan_pencere'],
                    sayilar=dict(pencere=ad_p['taranan_pencere'], ileri=len(ad_p['F']),
                                 geri=len(ad_p['R']), cift=0))

    izgara = U.izgara_tablosu_sayactan(top['sayac'])
    ornekler = {}
    for c in cl:
        s_, adh = U.hucre_etiketle(c)
        ornekler.setdefault(adh, '%s / %s (%d bp)' % (c['F'], c['R'], c['urun']))
    for x in izgara:
        x['ornek'] = ornekler.get(x['ad'], '')
    yaz('      izgara: %d/%d hucrede en az bir aday var'
        % (sum(1 for x in izgara if x['hayatta']), len(izgara)))
    for x in izgara[:4]:
        yaz('        %-62s %8d aday' % (x['ad'], x['hayatta']))
    yaz('      temsilci aday (her izgara hucresinden ornek): %d' % len(cl))

    # ---------------- ASAMA C: numune taramasi
    secili = sec_ornekle(cl, C.HUNI['numuneye_giden'])
    yaz('\n  [C] Numunenin ham okumalarinda in-silico PCR  (%d aday x %d kutu)'
        % (len(secili), len(baglam['uye_kutu']) + len(baglam['rakip_kutu'])))
    tc = time.time()
    olculen = numunede_olc(secili, numune, baglam, tc, yaz)
    yaz('\n      numune olcumu biten aday: %d   (%s)' % (len(olculen), sure(time.time() - tc)))
    olculen.sort(key=puan)

    # ---------------- ASAMA B2: ARMS - EN IYI adaylar uzerinde
    yaz("\n  [B2] ARMS varyantlari (ayirt edici 3' son baz + -2/-3 kasitli uyumsuzluk)")
    yaz('       Numunede en iyi olculen adaylar uzerinde denenir; varyantlar da')
    yaz('       AYNI olcume sokulur, boylece kazanc gercekten olculur.')
    arms = []
    if baglam['rakip_kons'] and olculen:
        uye_d = [k['dizi'] for k in baglam['uye_kons']][:4]
        rak_d = [k['dizi'] for k in baglam['rakip_kons']][:10]
        tab = time.time()
        for c in olculen[:C.HUNI['arms_taban']]:
            for yon, p in (('F', c['F']), ('R', c['R'])):
                ok, umm, rmm = U.ayirt_edici_mi(p, uye_d, rak_d)
                if not ok:
                    continue
                for v, et in U.arms_varyantlari(p):
                    mv = G.olc(v)
                    if not G.sabit_gecti(mv):
                        continue
                    yeni_c = dict(c)
                    yeni_c.pop('numune', None); yeni_c.pop('pm', None); yeni_c.pop('um', None)
                    if yon == 'F':
                        yeni_c['F'], yeni_c['mF'] = v, mv
                    else:
                        yeni_c['R'], yeni_c['mR'] = v, mv
                    yeni_c['arms'] = '%s %s (uye_mm=%s rakip_mm=%s)' % (yon, et, umm, rmm)
                    yeni_c['arms_taban'] = '%s / %s' % (c['F'], c['R'])
                    arms.append(yeni_c)
        yaz('      uretilen ARMS varyanti: %d  (%s)' % (len(arms), sure(time.time() - tab)))
        if arms:
            arms = arms[:C.HUNI['arms_ust']]
            ta2 = time.time()
            olculen_arms = numunede_olc(arms, numune, baglam, ta2, yaz)
            yaz('\n      numunede olculen ARMS varyanti: %d   (%s)'
                % (len(olculen_arms), sure(time.time() - ta2)))
            olculen += olculen_arms
            olculen.sort(key=puan)
    else:
        yaz('      atlandi (hafif mod, rakip konsensus yok ya da olculen aday yok)')
    yaz('      NOT: kasitli uyumsuzluk DEJENERE BAZ DEGILDIR (tek tanimli baz, 1 oligo);')
    yaz('           ama sablonla tam eslesmez - ayri bir toplanti maddesidir.')

    en_iyi = olculen[:C.HUNI['referansa_giden']]

    # ---------------- ASAMA D: geometri detayi + referans kapsam
    yaz('\n  [D] Cift geometrisi (hairpin/dimer dG, 60 C) + referans kapsam / rakip ayrimi')
    taxad = H.taxid_adlari()
    cinsler, rakip_cins = REF.uye_ve_rakip_anahtar(baglam, taxad)
    uye_havuz = REF.havuz_cikar(cinsler, baglam['siniflar'][0],
                                'uye_' + ad) if cinsler and not hafif else []
    rak_havuz = REF.havuz_cikar(rakip_cins, baglam['siniflar'][0],
                                'rakip_' + ad) if rakip_cins and not hafif else []
    yaz('      referans havuzu: uye %d dizi (%s) | rakip %d dizi (%s)'
        % (len(uye_havuz), ', '.join(cinsler[:4]), len(rak_havuz), ', '.join(rakip_cins[:4])))

    for i, c in enumerate(en_iyi, 1):
        # IKINCI OLCUT: panelin bazi satirlari <=3 ile olculmustu; en iyi
        # adaylar iki olcutle de verilir ki karsilastirilabilsinler.
        try:
            c['numune_mm3'] = numune.olc(c['F'], c['R'], baglam['uye_kutu'],
                                         baglam['rakip_kutu'], lo=C.URUN_IDEAL[0],
                                         hi=C.URUN_MUTLAK_UST, mm=C.KURESEL_MAX_MM
                                         if False else 3)
        except Exception:
            c['numune_mm3'] = {}
        c['cift'] = G.cift_olc(c['mF'], c['mR'], c['urun'])
        c['sikilik'], c['izgara_hucresi'] = U.hucre_etiketle(c)
        if uye_havuz:
            c['ref_uye'] = REF.kapsam(uye_havuz, c['F'], c['R'],
                                      C.URUN_IDEAL[0], C.URUN_MUTLAK_UST)
        if rak_havuz:
            c['ref_rakip'] = REF.kapsam(rak_havuz, c['F'], c['R'],
                                        C.URUN_IDEAL[0], C.URUN_MUTLAK_UST)
        if i % 10 == 0:
            print('      ... %d/%d' % (i, len(en_iyi)), end='\r', flush=True)

    gecen = [c for c in en_iyi if c['cift']['cift_gecti']]
    yaz('      cift yapisi (dTm/heterodimer/dG) gecen: %d/%d' % (len(gecen), len(en_iyi)))

    # ---------------- ASAMA E: kuresel ozgulluk
    son = gecen[:C.HUNI['kusele_giden']] if gecen else en_iyi[:3]
    if hafif or not son:
        yaz('\n  [E] Kuresel ozgulluk ATLANDI (hafif mod ya da aday yok)')
    else:
        yaz('\n  [E] KURESEL OZGULLUK - en pahali adim, %d aday, ~500 bin kayit' % len(son))
        te = time.time()
        durum_yolu = os.path.join(C.KONTROL, 'kuresel_%s.pkl'
                                  % ''.join(ch if ch.isalnum() else '_' for ch in ad))

        def ilerE(parca, kayit, gecen_sn):
            print('      ... parca %d, %d kayit, gecen %s          '
                  % (parca, kayit, sure(gecen_sn)), end='\r', flush=True)

        try:
            kr = KT.tara([dict(ad='a%d' % i, F=c['F'], R=c['R'],
                               lo=C.URUN_IDEAL[0], hi=C.URUN_MUTLAK_UST)
                          for i, c in enumerate(son)],
                         durum_yolu=durum_yolu, ilerle=ilerE)
            for i, c in enumerate(son):
                c['kuresel'] = kr.get('a%d' % i, {})
            yaz('\n      kuresel tarama bitti (%s)' % sure(time.time() - te))
        except Exception as e:
            yaz('\n      KURESEL TARAMA HATASI (diger sonuclar gecerli): %s' % e)

    sonuc = dict(
        hedef=ad, durum='TAMAMLANDI',
        gerekceler=[g[1] for g in satir['gerekceler']],
        etiketler=satir['etiketler'],
        panel=dict(F=satir['F'], R=satir['R'], urun=satir['urun_bp'],
                   plaka=satir['plaka'], ta=satir['ta'], ayrim=satir['ayrim'],
                   geo=satir['geo'], jel=satir['jel']),
        panel_olcum=taban,
        uyelik_uyarisi=uyelik_uyarisi(satir, taban),
        omurga=dict(kutu=om['kutu'], uzunluk=len(om['dizi']), yol=om['yol']),
        uyelik_kaynagi=baglam['uyelik_kaynagi'],
        uye_tax=baglam['uye_tax'],
        sayilar=dict(pencere=ad_p['taranan_pencere'], ileri=len(ad_p['F']),
                     geri=len(ad_p['R']), cift=top['toplam'],
                     temsilci=len(cl), arms=len(arms),
                     numune_olculen=len(olculen), cift_yapisi_gecen=len(gecen)),
        izgara=izgara,
        adaylar=[rapor.aday_ozet(c) for c in en_iyi],
        sure_sn=round(time.time() - t0, 1),
    )
    yaz('\n  BITTI: %s  (%s)' % (ad, sure(time.time() - t0)))
    return sonuc


def puan(c):
    """Siralama olcutu: ONCE uye kapsami (kac uye kutusu cogaliyor), SONRA ayrim.

    Tek bir uye kutusunun bos cikmasi uye_alt'i 0 yapar ve butun adaylari
    esitler; kapsam ekseni bunu onler.
    """
    # NEDEN KAPSAM ONCE GELIYOR: ayrim kati, uye kutularinin Wilson ALT sinirinin
    # en kucugu uzerinden hesaplanir. Tek bir uye kutusu hic urun vermezse o alt
    # sinir 0 olur ve oran da 0 cikar - iyi ile kotu aday ayni puani alir,
    # siralama coker. Kapsam (kac uye kutusu >=%20 urun veriyor) bu cokusten
    # etkilenmez, bu yuzden birinci anahtar odur. Rakip kumesi bos olan evrensel
    # hedeflerde ayrim kati zaten tanimsizdir; orada karsilastirmayi tasiyan tek
    # eksen kapsamdir. Sonraki iki anahtar once "yalniz kapsanan kutular"
    # uzerinden, sonra butun uye kutulari uzerinden hesaplanan ayrim katidir.
    n = c['numune']
    return (-n.get('uye_kapsam', 0),
            -(n.get('kat_enkotu_kapsayan') or n.get('kat_havuz_kapsayan') or 0),
            -(n.get('kat_enkotu') or n.get('kat_havuz') or 0))


def uyelik_uyarisi(satir, taban):
    """Olculen taban, panelin YAYIMLANMIS degerinden cok saparsa uyar.

    En sik sebep uyelik tanimidir (uye/rakip kutu listesi). Sessizce yanlis
    olcmektense yuksek sesle uyarmak dogrudur.
    """
    # UYELIK KOSULSUZ BENIMSENMEZ. Bu fonksiyon sapmayi yalniz BILDIRIR; uye
    # ya da rakip kutu listesini kendiliginden degistirmez, hedef_uyelik.tsv'ye
    # yazmaz. Ilke: kanit yoklugu kanit sayilmaz - bir kutunun yer degistirmesi
    # icin pozitif olcum kaniti gerekir, "sayi beklenenden farkli cikti" boyle
    # bir kanit degildir. Bu yuzden burada yapilan tek sey, kullaniciyi dosyanin
    # hangi satirina bakmasi gerektigine yonlendirmektir.
    # Bant 0,34x - 3,0x arasidir: olculen degerlerden HERHANGI BIRI panelin
    # yayimladigi degerin bu araliginda kaliyorsa uyari basilmaz. Bant genis
    # tutulmustur cunku panel satirlari farkli okuma derinliklerinde olculmustu
    # ve Wilson araliginin genisligi derinlige baglidir.
    if not taban:
        return []
    p = satir.get('ayrim_sayi')
    if p is None or 'kapsam' in (satir.get('ayrim') or '').lower() \
            or 'kutu' in (satir.get('ayrim') or '').lower():
        return []
    olculen = [x for x in (taban.get('kat_enkotu'), taban.get('kat_havuz'),
                           taban.get('kat_enkotu_kapsayan'),
                           taban.get('kat_havuz_kapsayan')) if x]
    if not olculen:
        return ['!! UYARI: mevcut cift icin ayrim orani HIC olculemedi - uye kutularinin',
                '   hicbiri urun vermiyor. Uyelik tanimi yanlis olabilir:',
                '   screening/hedef_uyelik.tsv -> satir "%s"' % satir['hedef']]
    if any(0.34 * p <= o <= 3.0 * p for o in olculen):
        return []
    return [
        '!! UYARI: panelin yayimladigi ayrim %.1fx, bu koşuda olculen %s.' % (
            p, ' / '.join('%.1fx' % o for o in olculen)),
        '   Sapma buyuk. En olasi sebep UYELIK TANIMI (hangi kutu uye, hangisi rakip).',
        '   Once su dosyaya bakin: screening/hedef_uyelik.tsv  ->  satir "%s"' % satir['hedef'],
        '   (Okuma sayisi dusukse Wilson araligi genisler, sapmanin bir kismi ondandir.)',
    ]


def numunede_olc(adaylar, numune, baglam, t0, yaz):
    """Aday listesini ham okumalarda olc, ilerlemeyi ve tahmini kalan sureyi bas."""
    out = []
    for i, c in enumerate(adaylar, 1):
        r = numune.olc(c['F'], c['R'], baglam['uye_kutu'], baglam['rakip_kutu'],
                       lo=C.URUN_IDEAL[0], hi=C.URUN_MUTLAK_UST)
        if r:
            c = dict(c); c['numune'] = r
            out.append(c)
        if i % 5 == 0 or i == len(adaylar):
            gecen = time.time() - t0
            kalan = gecen / i * (len(adaylar) - i)
            print('      ... aday %d/%d  gecen %s  tahmini kalan %s        '
                  % (i, len(adaylar), sure(gecen), sure(kalan)), end='\r', flush=True)
    return out


def sec_ornekle(adaylar, ust):
    """Her izgara hucresinden temsilci al - boylece 'hangi ayarda cozum var'
    sorusu butun ayarlar icin cevaplanabilir olur."""
    from collections import defaultdict
    kova = defaultdict(list)
    for c in adaylar:
        s, ad = U.hucre_etiketle(c)
        kova[(s, ad)].append(c)
    for k in kova:
        kova[k].sort(key=lambda c: (abs(c['urun'] - 105), abs(c['mF']['tm'] - c['mR']['tm'])))
    out, i = [], 0
    anahtarlar = sorted(kova, key=lambda k: k[0])
    while len(out) < ust:
        eklendi = False
        for k in anahtarlar:
            if i < len(kova[k]):
                out.append(kova[k][i]); eklendi = True
                if len(out) >= ust:
                    break
        if not eklendi:
            break
        i += 1
    return out


# ---------------------------------------------------------------- ana


def aramayi_kos(a, yaz, sure, cizgi, mod=None):
    """Kapsamli arama akisi. run_all.py de bunu cagirir."""
    # YON KAPISI - bu asamanin ilk isi. Nanopore okumalari cift yonludur ve
    # konsensus ters yonde uretilmis olabilir. Projenin butun in-silico PCR
    # motorlari (ispcr.amplify, okuma_motoru.Sonda) verilen diziyi YALNIZ arti
    # iplikte tarar; ters yonlu bir konsensuste ne ileri primer ne de ters
    # primerin tumleyeni bulunur. Motor hata atmaz, sessizce "urun yok" der -
    # yani butun gece kosan bir arama hicbir uyari vermeden "hicbir hedefte
    # cozum yok" uretir. Olculen kayip %100'dur (orientation_impact_test.py: dogru yonde
    # 117 urun, ters yonde 0). Bu yuzden konsensusler kanonik degilse asama
    # BASLAMAZ; SystemExit(2) ile durulur.
    if mod:
        a.mod = mod
    from .hepsi import yon_kapisi
    _ok, _m = yon_kapisi(yaz, 'kapsamli arama')
    for _x in _m:
        yaz('  ' + _x)
    if not _ok:
        yaz('')
        yaz('  *** GIRDI DOGRULAMASI BASARISIZ - BU ASAMA BASLATILMADI ***')
        yaz('  Sebep: okunacak konsensusler kanonik degil. Ters yonlu bir')
        yaz('  konsensuste in-silico PCR hicbir uyari vermeden 0 urun dondurur,')
        yaz('  yani butun kosu sessizce yanlis sonuc uretirdi.')
        yaz('  Cozum:  python3 screening/build_canonical.py --kok . --yeniden')
        raise SystemExit(2)

    sorunlu, panel, panel_yolu = H.sorunlu_hedefler()
    yaz('\nPanel dosyasi : %s' % panel_yolu)
    yaz('Paneldeki cift: %d     sorunlu bulunan: %d' % (len(panel), len(sorunlu)))
    for d in sorunlu:
        yaz('  [%-5s] %-36s %s' % (d['etiketler'], d['hedef'][:36],
                                   '; '.join(g[1] for g in d['gerekceler'])[:90]))

    if a.mod == 'tam':
        liste = panel
        for d in liste:
            d.setdefault('gerekceler', H.sorun_gerekceleri(d) or [('-', 'panelin tamami taraniyor')])
            d.setdefault('etiketler', '-')
    else:
        liste = sorunlu
    if a.hedef:
        liste = [d for d in liste if a.hedef.lower() in d['hedef'].lower()]

    if a.mod == 'devam':
        liste = [d for d in liste if not kontrol.bitti_mi(d['hedef'])]
        yaz('\nDEVAM MODU: bitmis hedefler atlaniyor. Kalan: %d' % len(liste))

    if not liste:
        yaz('\nIslenecek hedef yok. Rapor mevcut kontrol dosyalarindan uretiliyor.')
        rapor.uret(kontrol.hepsi(), panel, panel_yolu, yaz)
        return 0
    uyelik = H.uyelik_oku(); kons = H.konsensusler(); kut = H.kutular()
    baglamlar = {d['hedef']: H.hedef_baglami(d, uyelik, kons, kut) for d in liste}

    gerekli = {}
    for b in baglamlar.values():
        for k in b['uye_kutu'] + b['rakip_kutu']:
            gerekli[k['kutu']] = k
    yaz('\nHam okuma havuzlari kuruluyor: %d kutu x %d okuma' % (len(gerekli), a.okuma))

    def ilerK(i, n, ad):
        print('   ... %d/%d  %s          ' % (i, n, ad), end='\r', flush=True)

    tk = time.time()
    numune = N.Numune(list(gerekli.values()), n=a.okuma, ilerle=ilerK)
    yaz('\nHavuzlar hazir (%s)' % sure(time.time() - tk))

    if not a.hafif:
        taxad = H.taxid_adlari()
        istekler = []
        for d in liste:
            b = baglamlar[d['hedef']]
            cins, rakip = REF.uye_ve_rakip_anahtar(b, taxad)
            sf = b['siniflar'][0] if b['siniflar'] else 'B'
            istekler.append(('uye_' + d['hedef'], cins, sf))
            istekler.append(('rakip_' + d['hedef'], rakip, sf))
        yaz('\nReferans havuzlari cikariliyor (veritabani basina TEK gecis)...')

        def ilerR(db, n, bulunan):
            print('   ... %s: %d kayit tarandi, %d dizi toplandi        '
                  % (db, n, bulunan), end='\r', flush=True)

        tr = time.time()
        try:
            REF.toplu_cikar(istekler, ilerle=ilerR)
            yaz('\nReferans havuzlari hazir (%s)' % sure(time.time() - tr))
        except Exception as e:
            yaz('\nReferans havuzu cikarilamadi (%s) - referans adimi atlanacak.' % e)

    yaz('\nTAHMINI SURE: hedef basina ~%s; %d hedef -> ~%s'
        % (sure(600 if not a.hafif else 90), len(liste),
           sure((600 if not a.hafif else 90) * len(liste))))
    yaz('Kesintiye dayaniklidir: pencereyi kapatip yeniden acabilirsiniz, (3) ile devam eder.\n')

    for i, d in enumerate(liste, 1):
        try:
            s = hedefi_isle(d, baglamlar[d['hedef']], numune, i, len(liste), hafif=a.hafif)
            kontrol.yaz(d['hedef'], s)
            rapor.uret(kontrol.hepsi(), panel, panel_yolu, lambda *x: None)
        except KeyboardInterrupt:
            yaz('\n\nKULLANICI DURDURDU. Biten hedefler kayitli; (3) ile devam edin.')
            break
        except Exception:
            yaz('\n  HATA (%s) - bu hedef atlandi, digerleri suruyor:' % d['hedef'])
            traceback.print_exc()
            kontrol.yaz(d['hedef'], dict(hedef=d['hedef'], durum='HATA',
                                         hata=traceback.format_exc()[-2000:]))

    yolar = rapor.uret(kontrol.hepsi(), panel, panel_yolu, yaz)
    cizgi('=')
    yaz('  TOPLAM SURE: %s' % sure(time.time() - BASLANGIC))
    yaz('  SONUCLAR:')
    for p in yolar:
        yaz('    %s' % p)
    cizgi('=')
    return 0




def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--mod', default='sorunlu',
                    choices=['tam', 'sorunlu', 'devam', 'panel-olc', 'konsensus',
                             'uyelik', 'hepsi', 'ozet'])
    ap.add_argument('--sina', action='store_true')
    ap.add_argument('--hedefleri-listele', action='store_true')
    ap.add_argument('--sinama-atla', action='store_true',
                    help='kendini sinamayi atla (yalniz gelistirme/test icin)')
    ap.add_argument('--yeniden', action='store_true',
                    help='kontrol noktalarini yok say, bastan hesapla')
    ap.add_argument('--hafif', action='store_true', help='referans/kuresel adimlari atla')
    ap.add_argument('--tam-derinlik', action='store_true',
                    help='kutudaki BUTUN okumalari kullan (--okuma 0 ile ayni)')
    ap.add_argument('--okuma', type=int, default=C.NUMUNE_OKUMA_SAYISI,
                    help='kutu basina okuma; 0 = kutudaki BUTUN okumalar (yavas ama kesin)')
    ap.add_argument('--hedef', default=None, help='yalniz bu hedef (test icin)')
    ap.add_argument('--aday-ust', type=int, default=None)
    a = ap.parse_args(argv)

    if a.tam_derinlik:
        a.okuma = 0
    if a.aday_ust:
        C.HUNI['numuneye_giden'] = a.aday_ust

    if a.hedefleri_listele:
        try:
            sorunlu, panel, _ = H.sorunlu_hedefler()
            sor = {d['hedef'] for d in sorunlu}
            for d in panel:
                print('   %-46s %s' % (d['hedef'][:46],
                                       '<- sorunlu' if d['hedef'] in sor else ''))
        except Exception as e:
            print('Hedef listesi okunamadi:', e)
        return 0

    kontrol.ayar_kur(okuma=a.okuma, hafif=bool(a.hafif),
                     aday_ust=C.HUNI['numuneye_giden'])
    kontrol.hazirla()
    cizgi('=')
    yaz('  KAPSAMLI PRIMER ARAMASI - PrimerJury paneli')
    yaz('  baslangic: %s     mod: %s' % (time.strftime('%Y-%m-%d %H:%M:%S'), a.mod))
    yaz('  cikti klasoru: %s' % C.CIKTI)
    cizgi('=')

    from . import self_test
    if a.mod == 'ozet':
        from . import run_all as HP
        rc = HP.yalniz_ozet(yaz, sure, cizgi)
        yaz('  TOPLAM SURE: %s' % sure(time.time() - BASLANGIC))
        return rc

    if a.mod == 'hepsi':
        from . import run_all as HP
        rc = HP.calistir(yaz, sure, cizgi, a)
        yaz('  TOPLAM SURE: %s' % sure(time.time() - BASLANGIC))
        return rc

    if a.mod == 'uyelik':
        from . import membership_check
        if not a.sinama_atla and not kendini_sina.calistir(yaz):
            yaz('\nKENDINI SINAMA BASARISIZ - denetim baslatilmadi.')
            return 2
        uyelik_denetimi.calistir(yaz, sure,
                                 okuma_sayisi=(a.okuma or C.NUMUNE_OKUMA_SAYISI),
                                 yalniz=a.hedef, yeniden=a.yeniden)
        yaz('  TOPLAM SURE: %s' % sure(time.time() - BASLANGIC))
        return 0

    if a.mod == 'konsensus':
        from . import build_consensus
        konsensus_uret.calistir(yaz, sure, yalniz=a.hedef, yeniden=a.yeniden)
        yaz('  TOPLAM SURE: %s' % sure(time.time() - BASLANGIC))
        return 0

    if not a.sinama_atla and not kendini_sina.calistir(yaz):
        yaz('\nKENDINI SINAMA BASARISIZ - arama baslatilmadi.')
        return 2
    if a.sina:
        yaz('\nYalniz sinama istendi, arama yapilmadi.')
        return 0

    if a.mod == 'panel-olc':
        from . import panel_measurement
        panel_olcum.calistir(yaz, sure, okuma_sayisi=a.okuma if a.okuma else 0,
                             yalniz=a.hedef, yeniden=a.yeniden)
        yaz('  TOPLAM SURE: %s' % sure(time.time() - BASLANGIC))
        return 0

    return aramayi_kos(a, yaz, sure, cizgi)

if __name__ == '__main__':
    sys.exit(main())
