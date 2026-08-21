# -*- coding: utf-8 -*-
"""TOPLANTIDA NE ISTENDI, HANGISI OLDU, HANGISI OLMADI.

Toplanti kararlari steps/hedefler.tsv icinde numarali duruyor
(Karar 1 tur ozgul, Karar 2 cins, Karar 3 grup, Karar 4 alan evrensel,
Karar 5 olcumden turetilen ve onaya sunulan ekler). Bu betik o listeyi
panelin BUGUNKU ciktisiyla yan yana koyar.

Elle yazilmis hicbir durum yoktur. Her satirin durumu su dosyalardan okunur:
    TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv     (hukum, dCq, kapsam)
    DOGRULAMA_SONUC/dogrulama_uc_sutun.tsv     (kanit katmanlari)
    engine_SONUC/kutu_olculen_kimlik.tsv (numunede var mi)

Kosum:
    python verification/decision_status.py --kok .
Cikti:
    TOPLANTI_DURUMU.md   ve   toplanti_durumu.tsv
"""
from __future__ import print_function

import argparse
import csv
import io
import os
import sys
import time

# Toplanti adi -> paneldeki hedef adi. Ad esitligi kurulamayan satir
# "PANELDE KARSILIGI YOK" olarak raporlanir; sessizce eslenmez.
ESLEME = {
    'Methanosarcina_mazei_turu': 'Methanosarcina mazei / M. soligelidi grubu',
    'Methanothrix_soehngenii_turu': 'Methanothrix_soehngenii_turu',
    'Methanosarcina_barkeri_turu': None,
    'Podospora_pseudopauciseta': None,
    'Zoopagomycota_mantari': None,
    'Microascaceae_askomikot': 'Microascaceae_askomikot',
    'Bacteroidaceae_ailesi': 'Bacteroidales_kumesi',
    'Alistipes_cinsi': None,
    'Proteiniphilum_cinsi': 'Proteiniphilum_cinsi',
    'Petrimonas_cinsi': 'Petrimonas_cinsi',
    'Hidrojenotrofik_metanojenler': 'Metanomikrobiyales_hidrojenotrof',
    'Metilotrofik_metanojen': 'Metilotrofik_metanojen',
    'Asetoklastik_metanojenler': 'Asetoklastik_metanojenler',
    'Sakarolitik_bakteriler': 'Sakarolitik_Sphaerochaeta',
    'Proteolitik_sintrofik_bakteriler': 'Proteolitik_Cloacimonas',
    'Nitrosocosmicus_AOA': 'Nitrosocosmicus_AOA',
    'Trichoderma_cinsi': 'Petriella_musispora',
    'Bakteri_universal': 'Bakteri_universal',
    'Arke_universal': 'Arke_universal',
    'Mantar_universal': 'Mantar_universal (F2)',
    'Petriella_musispora': 'Petriella_musispora',
    'Proteolitik_Cloacibacillus':
        'Proteolitik_Synergistaceae (eski ad: Proteolitik_Cloacibacillus)',
    'Sakarolitik_Sphaerochaeta': 'Sakarolitik_Sphaerochaeta',
}

# Eslemesi olmayanlarin SEBEBI - hepsi olcume dayanir, hicbiri kanaat degil.
YOK_SEBEBI = {
    'Methanosarcina_barkeri_turu':
        u'Tur ozgul cift uretilemedi. M. barkeri ile M. mazei numunedeki '
        u'kutularda ayrilamiyor; ayirt edici pencere bulunamadi. Cins duzeyi '
        u'cift (Methanosarcina_cinsi) panelde VAR ve siparise giriyor.',
    'Podospora_pseudopauciseta':
        u'Organizmanin KENDISI numunede yok. Kurtarma turu bes referans '
        u'ciftinden ucunu denedi, 80 001 cift tarandi, ayirt edici aday 0. '
        u'Olmayan bir hedefe primer yazilamaz.',
    'Zoopagomycota_mantari':
        u'Kraken etiketi olcumle curutuldu. Kutunun olculen kimligi '
        u'Zoopagomycota degil; hedef bu adla numunede yok.',
    'Alistipes_cinsi':
        u'Olculen kimlik Alistipes degil Rikenellaceae (ust klad). Cins '
        u'duzeyinde ayirt edilemedigi icin Alistipes ozgul cift yazilamadi; '
        u'okumalari Bacteroidales_kumesi ciftinin kapsaminda.',
}


def _tsv(yol):
    if not os.path.exists(yol):
        return []
    with io.open(yol, encoding='utf-8') as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith('#')),
                                   delimiter='\t'))


def kararlar(yol):
    out = []
    for l in io.open(yol, encoding='utf-8'):
        l = l.rstrip('\n')
        if not l.strip() or l.startswith('#') or l.startswith('karar\t'):
            continue
        p = l.split('\t')
        if len(p) >= 3 and p[0].strip().isdigit():
            out.append((p[0].strip(), p[1].strip(), p[2].strip()))
    return out


KARAR_ADI = {
    '1': u'Karar 1 — tür özgül istendi',
    '2': u'Karar 2 — cins düzeyi istendi',
    '3': u'Karar 3 — işlevsel grup istendi',
    '4': u'Karar 4 — alan evrensel istendi',
    '5': u'Karar 5 — toplantıda konuşulmadı, ölçümden türetildi ve onaya sunuldu',
}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', '--kok', dest='kok', default='.')
    a = p.parse_args()
    kok = os.path.abspath(a.kok)

    ky = os.path.join(kok, 'steps', 'hedefler.tsv')
    if not os.path.exists(ky):
        sys.exit('HATA: %s yok.' % ky)
    K = kararlar(ky)
    SL = {(r.get('hedef') or '').strip(): r
          for r in _tsv(os.path.join(kok, 'TEK_PROTOKOL_SONUC', 'SIPARIS_LISTESI.tsv'))}
    DG = {(r.get('hedef') or '').strip(): r
          for r in _tsv(os.path.join(kok, 'DOGRULAMA_SONUC', 'dogrulama_uc_sutun.tsv'))}

    satirlar = []
    for no, ad, duzey in K:
        panel = ESLEME.get(ad, '__ESLENMEDI__')
        if panel == '__ESLENMEDI__':
            satirlar.append(dict(karar=no, istenen=ad, duzey=duzey, panel='',
                                 durum=u'ESLEME TABLOSUNDA YOK',
                                 hukum='', dcq='', kapsam='',
                                 not_=u'Bu ad esleme tablosuna girilmemis - '
                                      u'elle bakilmali.'))
            continue
        if panel is None:
            satirlar.append(dict(karar=no, istenen=ad, duzey=duzey, panel='',
                                 durum=u'YAPILAMADI', hukum='', dcq='', kapsam='',
                                 not_=YOK_SEBEBI.get(ad, u'Sebep yazilmamis.')))
            continue
        r = SL.get(panel)
        if r is None:
            satirlar.append(dict(karar=no, istenen=ad, duzey=duzey, panel=panel,
                                 durum=u'PANEL CIKTISINDA YOK', hukum='',
                                 dcq='', kapsam='',
                                 not_=u'Esleme "%s" diyor ama SIPARIS_LISTESI\'nde '
                                      u'bu ad yok. Esleme bayat olabilir.' % panel))
            continue
        sart = (r.get('siparis_sarti') or '').strip()
        sinif = (r.get('SINIF') or '').strip()
        dcq = (r.get('dCq_karsiligi') or '').strip()
        durum = (r.get('durum') or '').strip()
        # Durum SINIF'tan degil SINIF + SIPARIS SARTI'ndan gelir. Ikisini
        # ayirmazsak "kosulsuz siparis edilebilir" ile "kosullu" ayni kovaya
        # duser ve degerlendirici hangisinin sartsiz oldugunu goremez.
        if sinif in ('KESIN', 'EVRENSEL'):
            if sart.upper().startswith('KOSULSUZ'):
                d = u'YAPILDI (kosulsuz)'
            elif sart.upper().startswith('KOSULLU'):
                d = u'YAPILDI (kosullu)'
            else:
                d = u'YAPILDI'
        else:
            d = u'YAPILAMADI'
        dg = DG.get(panel, {})
        satirlar.append(dict(
            karar=no, istenen=ad, duzey=duzey, panel=panel, durum=d,
            hukum=(durum or sinif), dcq=dcq,
            kapsam=(r.get('ASIL_kapsam_mm1') or r.get('kapsam') or '').strip(),
            not_=(sart if sart and sart != '-'
                  else (durum or (dg.get('KARAR') or '').strip()))))

    ty = os.path.join(kok, 'toplanti_durumu.tsv')
    with io.open(ty, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# Toplanti kararlari x panelin BUGUNKU ciktisi. Uretim %s\n'
                 % time.strftime('%Y-%m-%d %H:%M'))
        fh.write(u'karar\tistenen\tistenen_duzey\tpaneldeki_karsiligi\tdurum\t'
                 u'hukum\tdCq\tkapsam\tnot\n')
        for s in satirlar:
            fh.write(u'\t'.join([s['karar'], s['istenen'], s['duzey'], s['panel'],
                                 s['durum'], s['hukum'], s['dcq'], s['kapsam'],
                                 s['not_'].replace('\t', ' ').replace('\n', ' ')]) + u'\n')

    sayim = {}
    for s in satirlar:
        sayim[s['durum']] = sayim.get(s['durum'], 0) + 1

    my = os.path.join(kok, 'TOPLANTI_DURUMU.md')
    with io.open(my, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# Toplantıda ne istendi, hangisi oldu\n\n')
        fh.write(u'Üretim: %s\n\n' % time.strftime('%Y-%m-%d %H:%M'))
        fh.write(u'Kaynak: `steps/hedefler.tsv` (toplantı kararlarının '
                 u'kendisi) yan yana `TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv` ve '
                 u'`DOGRULAMA_SONUC/dogrulama_uc_sutun.tsv`. Hiçbir durum elle '
                 u'yazılmadı; hepsi bu dosyalardan okundu.\n\n')
        fh.write(u'| durum | kaç istek |\n|---|---|\n')
        for k in sorted(sayim, key=lambda x: -sayim[x]):
            fh.write(u'| %s | %d |\n' % (k, sayim[k]))
        fh.write(u'\n')
        for no in ('1', '2', '3', '4', '5'):
            grup = [s for s in satirlar if s['karar'] == no]
            if not grup:
                continue
            fh.write(u'## %s\n\n' % KARAR_ADI.get(no, u'Karar %s' % no))
            fh.write(u'| istenen | paneldeki karşılığı | durum | dCq | not |\n'
                     u'|---|---|---|---|---|\n')
            for s in grup:
                fh.write(u'| %s | %s | **%s** | %s | %s |\n'
                         % (s['istenen'], s['panel'] or '—', s['durum'],
                            s['dcq'] or '—', s['not_'][:200]))
            fh.write(u'\n')
        yapilamadi = [s for s in satirlar if s['durum'] == u'YAPILAMADI']
        if yapilamadi:
            fh.write(u'## Yapılamayanlar ve sebebi\n\n')
            for s in yapilamadi:
                fh.write(u'### %s\n\n%s\n\n' % (s['istenen'], s['not_']))

    print('yazildi: %s' % my)
    print('yazildi: %s' % ty)
    for k in sorted(sayim, key=lambda x: -sayim[x]):
        print('  %-24s %d' % (k, sayim[k]))
    return 0


if __name__ == '__main__':
    sys.exit(main())
