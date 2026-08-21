# -*- coding: utf-8 -*-
"""FASTA basligindan TAKSONOMI cikarir - butun referans veritabanlari icin TEK yer.

NEDEN AYRI BIR MODUL (2026-08-21, A2)
-------------------------------------
Katman 2 (yerel veritabani) bir urunun hedefin ICINDEN mi DISINDAN mi geldigini
yalniz BOYA bakarak ayiriyordu. Katman 3 (MFEprimer, D-12) ayni soruyu TAKSONA
bakarak cevapliyor. D-12'de olculdu: MFEprimer'in "hedef disi" saydigi 1.605
amplikonun %95,7'si hedef kladin KENDI ICINDENdi, yalnizca boyu farkliydi.
Yani boy olcutu yaniltici; katman 2 de taksona bakmali.

Bunun onundeki gercek engel su: REFERANS_DB'deki dosyalar BES AYRI baslik
bicimi kullaniyor. Tek bir ayristirici varsaymak, kayitlarin cogunu "farkli
alan" (sinif c) saymak demektir - yani duzeltmeye calistigimiz hatanin daha
buyugunu uretmek. Bicimler OLCULDU (2026-08-21, her dosyanin ilk kaydi):

  SILVA  >AY846379.1.1791 Eukaryota;Archaeplastida;...;Monoraphidium sp.
         erisim BOSLUKLA ayrilir, taksonomi ';' ile
  UNITE  >UDB016649|k__Fungi;p__Basidiomycota;...;s__Thelephora_albomarginata|SH1281904.10FU
         '|' ile alanlara bolunur, taksonomi 2. alanda, jetonlar 'x__' onekli
  PR2    >AB353770.1.1740_U|18S_rRNA|nucleus||Eukaryota|TSAR|Alveolata|...
         '|' ile bolunur, taksonomi 5. alandan itibaren
  RefSeq >NR_201932.1 Sphingosinicella wutangchuni strain LY54 16S ribosomal RNA...
         TAKSONOMI YOK. Alan, veritabaninin TANIMINDAN gelir (bacteria.16S.fna
         icindeki her kayit tanimi geregi bakteridir). Tahmin degil.
  ROD    >GCA_000001215|AE014298.5/23211192-23217141|Eukaryota;Opisthokonta;...
         '|' ile bolunur, taksonomi 3. alanda, ';' ile

BASLIK KESILMESI
----------------
Cagiran taraf baslıgi KESMEMELIDIR. Olculdu: SILVA SSU basliklarinin %16,6'si
150 karakteri asiyor ve kesilen kuyruk tam da CINS ve TUR jetonlaridir - yani
hedef kladla eslesme ihtimali en yuksek olanlar. Kesik baslikla siniflandirmak,
klad ICI bir kaydi klad DISI saydirir.

BILINMIYOR AYRI BIR DURUMDUR
----------------------------
Alan cozulemezse '?' doner ve cagiran taraf bunu 'bilinmiyor' olarak SAYAR,
'farkli alan' saymaz. Cozulemeyen bir kayit, capraz reaksiyon KANITI degildir.
"""
from __future__ import unicode_literals

import re

ALANLAR = ('Bacteria', 'Archaea', 'Eukaryota')
ORGANEL_JETONLARI = ('Chloroplast', 'Mitochondria')

# Taksonomi tasimayan veritabanlarinda alan, dosyanin TANIMINDAN gelir.
# Anahtarlar dosya adindan turetilir (harf disi karakterler '_' olur).
VTB_ALAN = {
    'archaea_16S_fna': 'Archaea',
    'bacteria_16S_fna': 'Bacteria',
    'fungi_ITS_fna': 'Eukaryota',
    'fungi_28SrRNA_fna': 'Eukaryota',
    'fungi_18SrRNA_fna': 'Eukaryota',
    'ref_all2_fna': '?',        # karisik kume - alan kayittan cozulemez
    'ref_all_fna': '?',
}
VTB_KLAD = {
    'fungi_ITS_fna': 'Fungi',
    'fungi_28SrRNA_fna': 'Fungi',
    'fungi_18SrRNA_fna': 'Fungi',
}

# UNITE jeton oneki: k__Fungi -> Fungi
_UNITE_ONEK = re.compile(r'^[kpcofgs]__')
# UNITE ve PR2 tur adlarinda bosluk yerine '_' kullanilir
_ALT_CIZGI = re.compile(r'_+')


def vtb_anahtari(dosya_adi):
    """Dosya adini VTB_ALAN/VTB_KLAD anahtarina cevirir."""
    import os
    t = os.path.basename(dosya_adi or '')
    return re.sub(r'\W+', '_', t).strip('_')


def _temizle(jetonlar):
    out = []
    for j in jetonlar:
        j = (j or '').strip()
        if not j:
            continue
        j = _UNITE_ONEK.sub('', j)
        out.append(j)
        # 'Thelephora_albomarginata' -> ayrica 'Thelephora albomarginata'
        if '_' in j:
            out.append(_ALT_CIZGI.sub(' ', j))
    return out


def _alan_bul(jetonlar):
    for j in jetonlar:
        if j in ALANLAR:
            return j
    # UNITE kingdom'u alan degil; Fungi -> Eukaryota
    for j in jetonlar:
        if j in ('Fungi', 'Metazoa', 'Viridiplantae', 'Archaeplastida'):
            return 'Eukaryota'
    return '?'


def coz(baslik, vtb=''):
    """FASTA basligi -> (alan, [jetonlar], organel_mi, taksonomi_var_mi)

    baslik : '>' OLMADAN, KESILMEMIS baslik satiri
    vtb    : veritabani dosya adi (taksonomi tasimayan kaynaklar icin gerekli)

    alan '?' donebilir; bu 'bilinmiyor' demektir, 'farkli alan' DEMEK DEGILDIR.

    taksonomi_var_mi : kayitta GERCEK bir taksonomi dizgesi bulundu mu.
      False ise elimizde yalniz ORGANIZMA ADI vardir (RefSeq bicimi:
      'NR_201932.1 Sphingosinicella wutangchuni strain LY54 16S ...').
      Bu ayrim ZORUNLUDUR ve sonradan olculerek eklendi (2026-08-21): organizma
      adi yalniz CINS ve TUR tasir. 'Petrimonas' gibi bir CINS hedefi adda
      eslesir, ama 'Bacteroidales' (takim) ya da 'Microascaceae' (aile) gibi
      CINS USTU bir hedef hicbir tur adinda GECMEZ. Bu bilinmeden sinifladigimizda
      RefSeq'teki butun Bacteroidales uyeleri "klad disi" sayiliyordu - olculdu,
      bacteria.16S.fna'da 3.646 sahte capraz. Karar veremedigimiz yerde
      'bilinmiyor' demek zorundayiz.
    """
    b = (baslik or '').strip()
    if b.startswith('>'):
        b = b[1:].strip()
    if not b:
        return ('?', [], False, False)

    jet = []

    if '|' in b:
        alanlar = [x.strip() for x in b.split('|')]
        # UNITE: bir alanda 'x__' onekli ';' listesi var
        for a in alanlar:
            if ';' in a and _UNITE_ONEK.search(a):
                jet = _temizle(a.split(';'))
                break
        if not jet:
            # ROD: bir alanda ';' ile ayrilmis ve bir ALAN adiyla baslayan liste
            for a in alanlar:
                if ';' in a and a.split(';')[0].strip() in ALANLAR:
                    jet = _temizle(a.split(';'))
                    break
        if not jet:
            # PR2: '|' ile ayrilmis alanlarin icinde ALAN adi geciyor
            for i, a in enumerate(alanlar):
                if a in ALANLAR:
                    jet = _temizle(alanlar[i:])
                    break
    elif ';' in b:
        # SILVA: erisim BOSLUKLA ayrilir, sonrasi ';' ile taksonomi
        govde = b.split(' ', 1)[1] if ' ' in b.split(';', 1)[0] else b
        jet = _temizle(govde.split(';'))

    if jet:
        alan = _alan_bul(jet)
        if alan != '?':
            return (alan, jet, any(o in jet for o in ORGANEL_JETONLARI), True)

    # RefSeq bicimi: taksonomi YOK. Alan veritabaninin tanimindan gelir.
    k = vtb_anahtari(vtb)
    alan = VTB_ALAN.get(k, '?')
    ad = b.split(' ', 1)[1] if ' ' in b else b
    # 'Sphingosinicella wutangchuni strain LY54 16S ...' -> ilk iki kelime cins+tur
    kelime = ad.split()
    jet = _temizle(jet + [ad] + kelime[:2] + ([' '.join(kelime[:2])] if len(kelime) > 1 else []))
    if k in VTB_KLAD:
        jet.append(VTB_KLAD[k])
    if alan == '?':
        alan = _alan_bul(jet)
    return (alan, jet, any(o in jet for o in ORGANEL_JETONLARI), False)


def sinifla(baslik, vtb, hedef_klad, hedef_alan):
    """Bir vurusu D-12'nin dort sinifindan birine sokar.

      'a'          hedef kladin ICINDEN                        -> hedef disi DEGIL
      'ao'         hedef alan icinde ama ORGANEL (kloroplast/mitokondri)
      'b'          AYNI alan, klad DISI                        -> gercek capraz
      'c'          FARKLI alan                                 -> gercek capraz
      'bilinmiyor' KARAR VERILEMEDI - KANIT SAYILMAZ

    hedef_klad : hedefin klad jetonlari (KAPSAMLI_ARAMA/hedef_klad.tsv)
    hedef_alan : hedefin alani ('Bacteria' / 'Archaea' / 'Eukaryota')

    'bilinmiyor' UC durumda doner ve ucu de 'capraz yok' DEMEK DEGILDIR:
      1) alan hic cozulemedi
      2) hedef alan bilinmiyor
      3) kayitta taksonomi YOK (RefSeq) ve hedef klad organizma adinda
         GECMIYOR. Bu durumda kayit klad ICINDE de olabilir, DISINDA da -
         organizma adi yalniz cins ve tur tasidigi icin takim/aile duzeyindeki
         bir uyelik ADDAN OKUNAMAZ. Olculdu: bu ayrim yapilmadiginda
         bacteria.16S.fna'da Bacteroidales_kumesi icin 3.646 SAHTE capraz
         uretiliyordu (a=0 cikiyordu, ki imkansizdir).
    """
    alan, jet, organel, taksonomi_var = coz(baslik, vtb)

    ic = any(j in jet for j in (hedef_klad or []))
    if not ic and hedef_klad:
        # Taksonomi dizgesi olmayan kayitlarda (RefSeq) jeton tam eslesmeyebilir;
        # kelime sinirinda ara. 'Methanothrix' -> 'Methanothrix soehngenii'
        duz = ' '.join(jet)
        ic = any(re.search(r'\b%s' % re.escape(j), duz, re.I) for j in hedef_klad)

    if ic:
        return 'ao' if organel else 'a'
    if alan == '?' or not hedef_alan:
        return 'bilinmiyor'
    if alan != hedef_alan:
        # Alan farkli: bu karar TAKSONOMI GEREKTIRMEZ. RefSeq'te bile
        # veritabaninin tanimi alani belirler (archaea.16S.fna -> Archaea),
        # yani bir arke kaydinin bakteri hedefinden farkli alanda oldugu
        # KESINDIR. 'c' guvenle verilebilir.
        return 'c'
    # Buradan sonrasi: ayni alan, klad eslesmesi YOK.
    # Taksonomi varsa 'klad disi' KESINDIR -> 'b'.
    # Taksonomi yoksa BILEMEYIZ: organizma adi cins ustu uyeligi tasimaz.
    return 'b' if taksonomi_var else 'bilinmiyor'
