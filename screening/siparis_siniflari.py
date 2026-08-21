# -*- coding: utf-8 -*-
"""SIPARIS SINIFLARI + ESKI/YENI KIMLIK SUTUNLARI - tek kaynak.

NEDEN AYRI MODUL: iki bilgi de birden fazla ciktida gorunmek zorunda (siparis
listesi, panel tablosu, karar tablosu, birlesik ozet). Her ciktida yeniden
hesaplansaydi biri otekinden kayardi - bu projede daha once tam olarak bu oldu.

1) SINIFLANDIRMA (2026-08-06, kullanici karari)
   Esigi gecemeyen satirlar listeden CIKARILMAZ. Sessizce silmek yerine
   damgalanir ve laboratuvarda ne yapilmasi gerektigi yazilir; karar
   kullanicinindir.

     KESIN       dCq >= 3,0  (>= 8,00x)   -> siparis edilebilir
     KOSULLU     dCq 2,0-3,0 (4,00-8,00x) -> siparis edilebilir AMA dogrulama sart
     ONERILMEZ   dCq <  2,0  (<  4,00x)   -> listede gorunur, gerekcesi yazili
     EVRENSEL    oran tanimsiz            -> kapsama gore degerlendirilir

2) ESKI/YENI KIMLIK
   Her hedefin uye kutulari icin uc sutun yan yana konur:
     kraken_etiketi   toplantida konusulan ad (MEVCUT_KAYITLI_KIMLIK)
     olculen_kimlik   bizim dogruladigimiz ad (DOGRULANAN_KIMLIK)
     savunulabilir_duzey  tur / cins / aile / adlandirilamiyor
   Fark varsa satir ">>" ile isaretlenir; "Trichoderma sanilan sey aslinda
   Petriella" gibi eslesmeler tek bakista gorunsun diye.
   HEDEF ADLARI DEGISTIRILMEZ - kod olarak kullaniliyorlar; sutun EKLENIR.
"""
import os, csv, re, collections
from . import yapilandirma as C

KOSULLU_ALT_DCQ = 2.0                       # bunun altinda ONERILMEZ
KOSULLU_ALT_KAT = 2.0 ** KOSULLU_ALT_DCQ    # 4,00x

FARK_ISARETI = u'>> FARKLI'


def _f(v):
    try:
        return float(str(v).replace(',', '.'))
    except (TypeError, ValueError):
        return None


def sinifla(kat, olculemedi=False):
    """(sinif, dcq, esikten_uzaklik_dcq, laboratuvar_notu) dondurur."""
    if olculemedi or kat is None:
        return (u'EVRENSEL', None, None,
                u'Ayrim orani tanimsiz (rakip kumesi bos). Karar KAPSAMA gore '
                u'verilir. Laboratuvarda: erime egrisi + negatif kontrol (NTC) '
                u'yeterli; tek urun tepesi bekleniyor.')
    d = C.kat_dcq(kat)
    if d is None:
        # kat = 0 -> log tanimsiz. Sentinel sayi URETILMEZ (once -99 yaziliyordu
        # ve ciktida "-99,00 dCq" diye gorunuyordu); dCq bos birakilir.
        return (u'ONERILMEZ', None, None,
                u'Ayrim orani SIFIR: rakip kutular hedef kadar - ya da daha iyi - '
                u'cogaliyor, dCq hesaplanamiyor. Bu ciftle olculen sinyal hedefe '
                u'ATFEDILEMEZ. Laboratuvar dogrulamasi bu satiri KURTARMAZ; '
                u'primer yeniden tasarlanmali.')
    fark = round(d - C.ESIK_DCQ, 2)
    if kat >= C.AYRIM_ESIK:
        return (u'KESIN', d, fark,
                u'Erime egrisi + NTC yeterli. Tek urun tepesi ve NTC temizse '
                u'siparis dogrulanmis sayilir.')
    if kat >= KOSULLU_ALT_KAT:
        return (u'KOSULLU', d, fark,
                u'ESIGE YAKIN (dCq %.2f, esikten %.2f dCq geride). Erime egrisi '
                u'TEK BASINA YETMEZ: capraz urun ayni sicaklikta erir. JEL '
                u'GEREKIYOR (urun boyu dogrulanmali) + NTC SART. Jelde tek bant '
                u'cikmazsa ya da boy uymazsa amplikon DIZILEMESI gerekir.'
                % (d, fark))
    return (u'ONERILMEZ', d, fark,
            u'Ayrim cok dusuk (dCq %.2f, esikten %.2f dCq geride) - rakip kutular '
            u'hedef kadar iyi cogaliyor. Bu ciftle olculen sinyal hedefe '
            u'ATFEDILEMEZ. Siparis edilecekse AMPLIKON DIZILEMESI zorunludur; '
            u'jel ve erime egrisi bu satirda karar veremez. Onerimiz: once '
            u'primer yeniden tasarlansin.' % (d, fark))


def kimlik_tablosu(kok):
    """hedef -> dict(kraken, olculen, duzey, farkli_mi, uye_sayisi, ayrinti)

    Kaynak: TUM_KIMLIK_SONUC/tum_kutu_kimlikleri.tsv (G asamasi).
    Dosya yoksa BOS doner ve cagiran taraf sutunlari '-' yazar; sessizce
    yanlis ad UYDURULMAZ.
    """
    yol = os.path.join(kok, 'TUM_KIMLIK_SONUC', 'tum_kutu_kimlikleri.tsv')
    if not os.path.exists(yol):
        return {}
    with open(yol, encoding='utf-8', errors='ignore') as fh:
        sat = [s for s in fh if not s.startswith('#')]
    r = list(csv.reader(sat, delimiter='\t'))
    if len(r) < 2:
        return {}
    h = r[0]
    kutular = [dict(zip(h, x)) for x in r[1:] if len(x) > 3]
    hed = collections.defaultdict(list)
    for d in kutular:
        for t in (d.get('UYESI_OLDUGU_HEDEFLER') or '').split(','):
            t = t.strip()
            if t:
                hed[t].append(d)

    def _kisalt(v):
        """Uzun 'adlandirilamayan soy - EN YAKIN KAYIT: <upuzun taksonomi> (%X)'
        dizesini okunur hale getirir. Sayi ve anlam korunur, taksonomi atilir -
        tam hali zaten TUM_KUTU_KIMLIK_RAPORU.md icinde duruyor."""
        v = (v or '').strip()
        if v.lower().startswith('adlandirilamayan'):
            m = re.search(r'\(%\s*([\d,\.]+)\)\s*$', v)
            c = re.search(r'([A-Z][a-z]+);?[^;]*$', v.split('(%')[0]) if '(%' in v else None
            return (u'adlandirilamayan soy (%%%s)' % m.group(1)) if m else u'adlandirilamayan soy'
        return v

    def _teklestir(vals):
        """Ayni degerleri teke indirir, sirayi korur; bos olanlari atar."""
        out = []
        for v in vals:
            v = _kisalt(v)
            if v and v not in out:
                out.append(v)
        return out

    tablo = {}
    for t, ds in hed.items():
        kr = _teklestir(d.get('MEVCUT_KAYITLI_KIMLIK') for d in ds)
        ol = _teklestir(d.get('DOGRULANAN_KIMLIK') for d in ds)
        dz = _teklestir(d.get('SAVUNULABILIR_DUZEY') for d in ds)
        farkli = any((d.get('UYUSUYOR_MU') or '').upper().startswith('HAYIR')
                     for d in ds)
        tablo[t] = dict(
            kraken=' | '.join(kr) if kr else '-',
            olculen=' | '.join(ol) if ol else '-',
            duzey=' | '.join(dz) if dz else '-',
            farkli_mi=farkli,
            uye_sayisi=len(ds),
            farkli_kutu=sum(1 for d in ds
                            if (d.get('UYUSUYOR_MU') or '').upper().startswith('HAYIR')),
        )
    return tablo


def kimlik_sutun_basliklari():
    return ['kraken_etiketi', 'olculen_kimlik', 'savunulabilir_duzey', 'ad_farkli_mi']


def kimlik_sutunlari(tablo, hedef):
    k = tablo.get(hedef)
    if not k:
        # Sessiz '-' yerine SEBEBI yazilir: hedefin G asamasinda uyelik tanimi
        # yok demektir (ornek: panel disi 'ek' ciftler). Eksik veriyle
        # dolu veri birbirine karismasin.
        y = u'G asamasinda uyelik tanimi YOK'
        return [y, y, y, u'karsilastirilamadi']
    return [k['kraken'], k['olculen'], k['duzey'],
            (u'%s (%d/%d kutu)' % (FARK_ISARETI, k['farkli_kutu'], k['uye_sayisi']))
            if k['farkli_mi'] else u'ayni']


# ---------------------------------------------------------------------------
# KIMLIGIN KAYNAGI - "bunu nereden biliyorsunuz?" sorusunun cevabi
#
# HER_VTB_NE_DEDI sutunu bu bilgiyi zaten tasiyordu ama tek bir devasa hucrede
# gomuluydu; kimse okumuyordu. Ayrilip ayri sutunlara konur:
#
#   karar_veren_vtb    kimligi veren veritabani/veritabanlari
#   kimlik_yuzdesi     o veritabanindaki hizalama kimligi
#   hizalanan_uzunluk  kac baz uzerinden olculdu (referans kaydin uzunlugu)
#   uyusan_vtb_sayisi  kac veritabani AYNI cinsi soyledi / kac tanesi cevap verdi
#   uyusmayan_vtb      uyusmayanlar NE DEDI - gizlenmez, yan yana yazilir
#
# Amac: degerlendirici "bunu nereden biliyorsunuz" dediginde satir gosterilebilsin.
# ---------------------------------------------------------------------------
_VTB_PAT = re.compile(r'([A-Za-z0-9_ .]+?) \[([^\]]*)\]: (.*?)(?=(?: \| [A-Za-z0-9_ .]+? \[)|$)')


def _cins(ad):
    """Isabet basligindan CINS adini cikarir.

    Basliklar uc bicimde geliyor ve ilk buyuk harfli kelimeyi almak YANLISTI -
    SILVA basliklarinda o kelime alan adidir ("Archaea", "Bacteria", "Eukaryota")
    ve butun kayitlar ayni cikardi:
      SILVA/PR2 : '<erisim> Alan;Sube;...;Cins;Tur'        -> son anlamli parca
      UNITE     : '<erisim>|k__Fungi;...;g__Cins;s__Tur|SH' -> g__ etiketi
      RefSeq    : 'NR_xxxxx.1 Cins tur strain ...'          -> erisimden sonraki
    """
    ad = (ad or '').strip()
    if not ad or ad.startswith('-'):
        return ''
    if ad.lower().startswith('adlandirilamayan'):
        return u'(adlandirilamayan)'
    m = re.search(r'[;|]g__([A-Za-z0-9_]+)', ad)          # UNITE
    if m:
        return m.group(1).replace('_', ' ').split()[0]
    if ';' in ad:                                          # SILVA / PR2 / ROD
        parca = [p.strip() for p in ad.split(';') if p.strip()]
        for p in reversed(parca):
            p = p.split('(')[0].strip()
            k = re.match(r'([A-Z][A-Za-z0-9_-]{2,})', p)
            if k and k.group(1).lower() not in (
                    'uncultured', 'unidentified', 'metagenome'):
                return k.group(1)
        return ''
    m = re.match(r'\s*[A-Z]{1,2}[A-Z_]*\d+\.\d+\s+([A-Z][a-z]{3,})', ad)   # RefSeq/nt
    if m:
        return m.group(1)
    m = re.search(r'\b([A-Z][a-z]{3,})\b', ad)
    return m.group(1) if m else ''


def vtb_ayristir(hucre):
    """HER_VTB_NE_DEDI hucresi -> [(vtb, kayit_sayisi, sonuc_metni, yuzde), ...]"""
    out = []
    for m in _VTB_PAT.finditer(hucre or ''):
        vtb, kapsam, sonuc = m.group(1).strip(), m.group(2), m.group(3).strip()
        y = re.search(r'\(%\s*([\d,\.]+)\)', sonuc)
        out.append((vtb, kapsam, sonuc, y.group(1) if y else ''))
    return out


def kaynak_sutun_basliklari():
    return ['karar_veren_vtb', 'kimlik_yuzdesi', 'hizalanan_uzunluk',
            'uyusan_vtb_sayisi', 'uyusmayan_vtb_ne_dedi']


def kaynak_sutunlari_kutu(d):
    """G asamasinin BIR KUTU satirindan kaynak sutunlarini uretir."""
    if not d:
        return ['-'] * 5
    kazanan = (d.get('en_iyi_vtb') or '-').strip()
    yuzde = (d.get('en_iyi_kimlik_%') or '-').strip()
    basl = d.get('en_iyi_isabet') or ''
    # hizalanan uzunluk: SILVA basliklarinda kayit adi '.<bas>.<son>' tasir
    m = re.match(r'\S*?\.(\d+)\.(\d+)\s', basl)
    uzun = str(abs(int(m.group(2)) - int(m.group(1))) + 1) if m else '-'
    # Kutunun dogrulanan adi 'adlandirilamayan soy' ise cins karsilastirmasi
    # ondan YAPILAMAZ; kazanan isabetin kendi cinsi kullanilir.
    _da = (d.get('DOGRULANAN_KIMLIK') or '').strip()
    hedef_cins = _cins(basl if (not _da or _da.lower().startswith('adlandirilamayan'))
                       else _da)
    uy, uym = [], []
    for vtb, _kap, sonuc, _y in vtb_ayristir(d.get('HER_VTB_NE_DEDI', '')):
        if re.match(r'(SONUC YOK|SORULMADI|ELLE|TAMAM\s*$|-\s*$)', sonuc):
            continue
        c = _cins(sonuc)
        if not c:
            continue
        # ONEK KARSILASTIRMASI: kaynak TSV basliklari 110 karakterde KESIYOR,
        # yani taksonomi satirinin son parcasi yarim kalabiliyor
        # ('Methanosarcina' -> 'Methanosarcin'). Tam esitlik aransaydi ayni
        # veritabani kendi verdigi kimlikle 'uyusmuyor' gorunurdu - nitekim
        # goruldu. Ilk 8 karakter uzerinden karsilastirilir.
        _e = lambda x: (x or '')[:8].lower()
        (uy if (hedef_cins and _e(c) == _e(hedef_cins)) else uym).append(
            u'%s: %s' % (vtb, sonuc[:52]))
    cevap = len([1 for _v, _k, s, _y in vtb_ayristir(d.get('HER_VTB_NE_DEDI', ''))
                 if not re.match(r'(SONUC YOK|SORULMADI|ELLE)', s)])
    return [kazanan, yuzde, uzun,
            u'%d / %d cevap veren' % (len(uy), cevap) if cevap else '-',
            (u' || '.join(uym[:6]) if uym else u'yok - uyusmayan veritabani cikmadi')]


def kaynak_tablosu(kok):
    """hedef -> kaynak sutunlari (uye kutularin EN IYI kimligini veren satirdan)."""
    yol = os.path.join(kok, 'TUM_KIMLIK_SONUC', 'tum_kutu_kimlikleri.tsv')
    if not os.path.exists(yol):
        return {}
    with open(yol, encoding='utf-8', errors='ignore') as fh:
        sat = [s for s in fh if not s.startswith('#')]
    r = list(csv.reader(sat, delimiter='\t'))
    if len(r) < 2:
        return {}
    h = r[0]
    kutular = [dict(zip(h, x)) for x in r[1:] if len(x) > 3]
    hed = collections.defaultdict(list)
    for d in kutular:
        for t in (d.get('UYESI_OLDUGU_HEDEFLER') or '').split(','):
            t = t.strip()
            if t:
                hed[t].append(d)

    def _yuz(d):
        try:
            return float((d.get('en_iyi_kimlik_%') or '0').replace(',', '.'))
        except ValueError:
            return 0.0

    out = {}
    for t, ds in hed.items():
        # hedefi temsil eden kutu: kimligi EN YUKSEK olan (karari o veriyor)
        out[t] = kaynak_sutunlari_kutu(sorted(ds, key=_yuz, reverse=True)[0])
    return out


def kaynak_sutunlari(tablo, hedef):
    return tablo.get(hedef, [u'G asamasinda uyelik tanimi YOK'] * 3 + ['-', '-'])
