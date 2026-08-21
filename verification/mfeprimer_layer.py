# -*- coding: utf-8 -*-
"""MFEprimer layer, the independent third tool.

Wraps the MFEprimer 4.4 binary for off-target amplicons, hairpins and dimers.
Its value is that it is NOT our code: a bug in our engine cannot confirm itself
here.

Two traps this module handles, both measured the hard way:
  * MFEprimer does not overwrite an existing output file, it returns 1 and
    leaves the old result in place, so a re-run would silently reuse stale
    numbers. Old outputs are deleted first.
  * "Off-target" in MFEprimer means BY SIZE ONLY. For group and universal
    primers, members of the target clade legitimately amplify at other lengths:
    of 1,605 amplicons flagged off-target, 95.7% were inside the target clade.
    klad_siniflandir() therefore classifies each hit taxonomically, and only
    (same domain outside clade) + (different domain) reach the verdict.

--- ozgun aciklama ---
MFEprimer katmani - BAGIMSIZ ucuncu arac.

NEDEN VAR
---------
Ilk iki katman da BIZIM kodumuz: numune olcumu ve yerel veritabani taramasi ayni
motoru (screening/read_engine.py, ispcr.py) kullaniyor. O motorda bir hata
varsa iki katman da AYNI yonde yanilir ve birbirini "dogrular". Bu katman disaridan
gelen, bagimsiz yazilmis bir araci (MFEprimer 4.4.0) ayni sorulara sokar.

NE YAPAR
  * spec    : hedef disi amplikon taramasi (klasordeki .primerqc.bin indeksleri)
  * hairpin : her oligo icin sac tokasi
  * dimer   : oligo cifti arasi dimer

INDEKSLER: REFERANS_DB altinda .primerqc.bin dosyalari hazir. Gecmiste
"RecordCount=0, uyumsuz" notu dusulmustu; bu bir IKILI/BICIM uyusmazligiydi,
indeksler saglam. Bu modul her kosuda once indeksi SINAR, okunmuyorsa o
veritabanini atlar ve sebebini yazar - sessizce gecmez.
"""

# -------------------------------------------------------------------------
# mfeprimer_layer.py, MFEprimer 4.4.0 ikilisini (hedef disi amplikon, hairpin,
# dimer) surer; dogrulama zincirinin BAGIMSIZ ucuncu kanit katmanidir.
#
# GİRDİ  : REFERANS_DB/ altindaki .fna/.fasta kumeleri ve onlarin
#          .primerqc.bin indeksleri (MFE_DB listesi), tools/mfeprimer ikilisi,
#          dogrulanacak primer ciftleri (cagiran betikten gelir).
# ÇIKTI  : DOGRULAMA_SONUC/mfe/ altina girdi ve ham cikti dosyalari
#          (girdi_*.tsv, spec_*.txt, hairpin.txt, dimer.txt) + cift/veritabani
#          basina JSON kontrol noktalari. Kendi basina rapor yazmaz; sonucu
#          specificity_round.py birlestirir.
# ÇAĞRAN : verification/full_chain.py -> D tusu (dolayli: specificity_round.py bu modulu
#          "KATMAN 3" olarak yukler; --mfe-yok verilirse atlanir).
#
# NEDEN BAGIMSIZ ARAC SART: numune olcumu ile yerel veritabani taramasi BIZIM
# ayni motorumuzu kullanir. O motorda bir hata varsa iki katman da ayni yonde
# yanilir ve birbirini "dogrular". Disaridan gelen bir arac ayni sorulara
# sokulmadikca celiski gorunmez.
# -------------------------------------------------------------------------
import os, re, csv, json, subprocess, time

MFE_ADAYLARI = [os.path.join('tools', 'mfeprimer'), 'mfeprimer']
# spec taramasinda kullanilacak indeksler (kucukten buyuge - sure kontrolu icin)
MFE_DB = ['archaea.16S.fna', 'bacteria.16S.fna', 'fungi.ITS.fna',
          'fungi.28SrRNA.fna', 'fungi.18SrRNA.fna',
          'SILVA_138.2_SSURef_NR99.fasta']
URUN_ALT, URUN_UST = 70, 400
HAIRPIN_TM_UST = 45.0          # panelin kendi geometri kurali
DIMER_TM_UST = 45.0
# D-5: MFEprimer dimer/hairpin kayitlarinda Tm raporlamaz; raporladigi olcu
# Delta G'dir. Bu esik, 3' dimerleri icin yaygin kullanilan -9 kcal/mol
# sinirindan alinmistir (IDT/Primer3 pratigi). Tm<45 kurali bu araca
# UYGULANAMAZ - uygulanmis gibi yapmak sessiz bir 'ihlal yok' oyu uretirdi.
DG_ESIGI = -9.0


# tools/mfeprimer once proje icinde aranir, sonra PATH uzerinde. Proje icindeki
# kopyaya calistirma izni verilir: depodan/arsivden cikan dosya cogu zaman izinsiz
# gelir ve katman sessizce "arac yok" diye atlanirdi.
def mfe_bul(kok):
    for a in MFE_ADAYLARI:
        y = a if os.path.isabs(a) else os.path.join(kok, a)
        if os.path.exists(y):
            try:
                os.chmod(y, 0o755)
            except OSError:
                pass
            return y
        if os.path.dirname(a) == '':
            from shutil import which
            w = which(a)
            if w:
                return w
    return None


def indeks_sina(kok, mfe, dosya, yaz, sure_tavani=120):
    """Indeks GERCEKTEN okunuyor mu? Tek bir sahte ciftle sinar."""
    db = os.path.join(kok, 'REFERANS_DB', dosya)
    if not os.path.exists(db):
        return (False, u'fasta yok')
    if not os.path.exists(db + '.primerqc.bin'):
        return (False, u'.primerqc.bin indeksi yok - "mfeprimer index -i %s" ile kurulmali' % dosya)
    import tempfile
    with tempfile.TemporaryDirectory() as t:
        gi = os.path.join(t, 'sina.tsv'); co = os.path.join(t, 'sina.txt')
        open(gi, 'w').write('SINAMA\tCTGCGGTTTAATTGGATTCAACGC\tGAACTGACGACGGCCATGC\n')
        try:
            pr = subprocess.run([mfe, 'spec', '-i', gi, '-o', co, '-d', db,
                                 '-c', '2', '-S', str(URUN_UST)],
                                capture_output=True, text=True, timeout=sure_tavani)
        except subprocess.TimeoutExpired:
            return (False, u'indeks sinamasi zaman asimina ugradi (%d sn)' % sure_tavani)
        if pr.returncode != 0:
            return (False, u'mfeprimer hata verdi: %s' % ((pr.stderr or '')[:160]))
        if not os.path.exists(co):
            return (False, u'cikti uretilmedi')
        s = open(co, encoding='utf-8', errors='replace').read()
        if 'RecordCount=0' in s or 'record count' in s.lower():
            return (False, u'indeks RecordCount=0 bildiriyor - yeniden kurulmali')
        if 'Primer ID' not in s:
            return (False, u'cikti ayristirilamadi')
        # D-4 HATA DUZELTMESI (2026-08-06): asil sessiz basarisizlik burasi.
        # SILVA_138.2_SSURef_NR99.fasta.primerqc.bin 515 MB olmasina ragmen
        # HER primer icin 'Binding Number  Plus 0  Minus 0' donuyor ve amplikon
        # sayisi 0 cikiyor. Ayni FASTA'da bizim yerel tarayici Bakteri_universal
        # icin 419.090 vurus buluyor - yani indeks okunmuyor (muhtemelen baska
        # bir MFEprimer surumuyle kurulmus). Eski kapi bunu KACIRIYORDU: yalniz
        # 'RecordCount=0' dizgesini ve 'Primer ID' varligini ariyordu, ikisi de
        # bu durumda saglam gorunuyor. Sonuc: en buyuk veritabani sessizce
        # "0 hedef disi amplikon" oyu veriyordu. Artik baglanma sayisi TOPLAM
        # sifirsa indeks BOZUK sayilir ve o veritabani ATLANANLAR'a yazilir.
        _bag = re.findall(r'^\S+\s+[ACGTUNRYKMSWBDHVacgtu]{10,}\s+\d+\s+[\d.]+\s+'
                          r'[\d.]+\s+[-\d.]+\s+(\d+)\s+(\d+)\s*$', s, re.M)
        if _bag and all(int(a) == 0 and int(b) == 0 for a, b in _bag):
            return (False, u'indeks BOZUK: butun primerler icin Binding Number 0/0 '
                           u'donuyor (dosya var ama okunmuyor). "mfeprimer index -i %s" '
                           u'ile YENIDEN kurulmali.' % dosya)
    return (True, u'okunuyor')


def _spec_ayristir(yol, beklenen_bp=None, tolerans=10):
    """Amplikonlari ayristir.

    ONEMLI: MFEprimer hedefin KENDI amplikonunu da sayar. Beklenen urun boyu
    verilirse o boydakiler AYRI sayilir ('ayni_boyda') ve 'hedef_disi' yalnizca
    FARKLI boydakileri kapsar. Ayni boyda cikan bir amplikon baska bir
    organizmada da olabilir - o yuzden ayri sutun olarak raporlanir, yok
    sayilmaz."""
    if not os.path.exists(yol):
        return {}
    s = open(yol, encoding='utf-8', errors='replace').read()
    m = re.search(r'Descriptions of \[\s*(\d+)\s*\] potential amplicons', s)
    toplam = int(m.group(1)) if m else 0
    boylar = [int(x) for x in re.findall(r'^\d+\s+\S+\s+(\d+)\s+', s, re.M)]
    ayni = hedef_disi = 0
    if beklenen_bp:
        for b in boylar:
            if abs(b - int(beklenen_bp)) <= tolerans:
                ayni += 1
            else:
                hedef_disi += 1
    else:
        hedef_disi = toplam
    return dict(_toplam=toplam, _boylar=boylar[:50],
                ayni_boyda=ayni, hedef_disi=hedef_disi)


def spec_kos(kok, mfe, ciftler, CIKTI, yaz, kontrol, sure_tavani=1800):
    """Her cift icin hedef disi amplikon taramasi. Cift bazinda kontrol noktasi."""
    d = os.path.join(CIKTI, 'mfe')
    os.makedirs(d, exist_ok=True)
    # HER CIFT AYRI KOSULUR. Tek dosyada birden fazla cift verilirse MFEprimer
    # tek bir toplam sayi yaziyor ve o sayi ciftlere DAGITILAMIYOR - butun
    # ciftler ayni degeri alirdi.
    girdiler = []
    for c in ciftler:
        ad = re.sub(r'\W+', '_', c['hedef'])
        gi = os.path.join(d, 'girdi_%s.tsv' % ad)
        with open(gi, 'w', encoding='utf-8') as fh:
            fh.write('%s\t%s\t%s\n' % (ad, c['F'], c['R']))
        girdiler.append((c, ad, gi))

    kullanilan, atlanan = [], []
    for dosya in MFE_DB:
        ok, not_ = indeks_sina(kok, mfe, dosya, yaz)
        if ok:
            kullanilan.append(dosya)
            yaz(u'    indeks OK : %s' % dosya)
        else:
            atlanan.append((dosya, not_))
            yaz(u'    index SKIPPED: %s - %s' % (dosya, not_))
    if not kullanilan:
        return dict(durum='ATLANDI', sebep=u'hicbir MFEprimer indeksi okunamadi',
                    atlanan=atlanan), {}

    sonuc = {}
    for c, ad, gi in girdiler:
        sonuc[c['hedef']] = {}
        for dosya in kullanilan:
            # 2026-08-10 DIZI MUHRU: dosya adi yalnizca hedef adini tasiyordu,
            # dizi degisince ayni dosya okunup eski MFEprimer sonucu geri
            # veriliyordu. Dizi ozeti artik adin parcasi.
            import hashlib as _hl
            _dz = _hl.md5(((c.get('F') or '') + '|' + (c.get('R') or ''))
                          .encode('utf-8')).hexdigest()[:10]
            kp = os.path.join(kontrol, 'mfe_%s_%s_%s.json'
                              % (ad, _dz, re.sub(r'\W+', '_', dosya)))
            # D-9 HATA DUZELTMESI (2026-08-07): ZEHIRLI KONTROL NOKTASI.
            # Kontrol noktalari yalnizca "dosya var mi" diye bakiliyordu. SILVA
            # indeksi 2026-08-06 13:0x'te BOZUKKEN kosmus ve 16 ciftin hepsi icin
            # {"_toplam":0,"hedef_disi":0,"sure":2.5} yazilmisti. Indeks 23:30'da
            # YENIDEN KURULDU ama o dosyalar yerinde durdugu icin sonraki her
            # kosu SILVA'yi HIC KOSMADAN o sifirlari geri okuyacakti - yani
            # "en buyuk veritabaninda hedef disi yok" diye SESSIZ ve YANLIS bir
            # guvence uretecekti. Tam da kacinmamiz gereken hata deseni.
            # Kural: indeks kontrol noktasindan YENIYSE kontrol noktasi GECERSIZ.
            if os.path.exists(kp):
                _ix = os.path.join(kok, 'REFERANS_DB', dosya + '.primerqc.bin')
                _bayat = False
                try:
                    if os.path.exists(_ix) and os.path.getmtime(_ix) > os.path.getmtime(kp):
                        _bayat = True
                except OSError:
                    pass
                if _bayat:
                    yaz(u'    checkpoint STALE (the index is newer), re-running: %s / %s' % (ad[:28], dosya))
                    try:
                        os.remove(kp)
                    except OSError:
                        pass
                else:
                    try:
                        sonuc[c['hedef']][dosya] = json.load(open(kp, encoding='utf-8'))
                        continue
                    except Exception:
                        pass
            co = os.path.join(d, 'spec_%s_%s.txt' % (ad, re.sub(r'\W+', '_', dosya)))
            # D-10 HATA DUZELTMESI (2026-08-07): MFEprimer 4.4.0 cikti dosyasi
            # VARSA CALISMAYI REDDEDER:
            #   "... is already exists, please remove it first"  (returncode 1)
            # Bu gunes gormeyen bir tuzakti: eski kosularin spec_*.txt dosyalari
            # DOGRULAMA_SONUC/mfe/ altinda duruyor (16 adet SILVA ciktisi, hepsi
            # indeks bozukken yazildigi icin "0 potential amplicons"). Kontrol
            # noktasi tutarken mfeprimer HIC cagrilmadigi icin hata gorunmuyordu;
            # bayat kontrol noktasi temizlenip yeniden kosulmak istendiginde ise
            # mfeprimer 1 donduruyor, biz onu dict(hata=...) yazip geciyoruz ve
            # SILVA yine OLCULMEMIS oluyordu. Yani "indeksi duzelttik" dedigimiz
            # halde sonuc yine sessizce eksik cikardi.
            # Cozum: kosmadan once eski ciktiyi (ve yanindaki .spec.tsv/.mfe.log)
            # SIL. Silinemiyorsa SESSIZ GECME - acik hata yaz.
            _silinemedi = None
            for _es in (co, co + '.spec.tsv', co + '.mfe.log'):
                if os.path.exists(_es):
                    try:
                        os.remove(_es)
                    except OSError as _e:
                        if _es == co:
                            _silinemedi = str(_e)
            if _silinemedi:
                sonuc[c['hedef']][dosya] = dict(
                    hata=u'eski cikti dosyasi silinemedi, MFEprimer uzerine YAZMAZ: %s'
                         % _silinemedi)
                yaz(u'    ERROR: the old output could not be deleted, so this database was not measured: %s'
                    % co)
                continue
            t0 = time.time()
            # D-11 HIZ DUZELTMESI (2026-08-07, OLCULDU): MFEprimer amplikon
            # raporunu COK SAYIDA KUCUK yazma ile uretir. Bagli Windows klasorune
            # (/mnt/c/...) yazarken her yazmanin gecikmesi toplanir ve is
            # taninmayacak kadar yavaslar. OLCUM:
            #   Proteolitik_Synergistaceae, SILVA, ayni komut:
            #     yerel diske (/tmp)     : 5,2 sn   (11,5 MB cikti)
            #     bagli klasore dogrudan : 3 dakikada BITMEDI (2,35 MB yazilmisti)
            #   Bakteri_universal, SILVA:
            #     yerel diske (/tmp)     : 49,8 sn  (282 MB cikti)
            #     bagli klasore dogrudan : 2026-08-06 kosusunda ~38 dakika
            # Buna karsilik TOPLU kopyalama hizlidir (233 MB / 26 sn olculdu),
            # cunku tek buyuk ardisik yazmadir. Bu yuzden: mfeprimer YEREL
            # diskte kosar, kanit dosyalari is bitince TOPLU kopyalanir.
            # Kanit dosyalarindan hicbiri kaybolmaz, yalnizca yazma deseni degisir.
            import tempfile, shutil
            _yerel = tempfile.mkdtemp(prefix='mfe_spec_')
            _co = os.path.join(_yerel, os.path.basename(co))
            try:
                try:
                    pr = subprocess.run([mfe, 'spec', '-i', gi, '-o', _co, '-d',
                                         os.path.join(kok, 'REFERANS_DB', dosya),
                                         '-c', '4', '-s', str(URUN_ALT), '-S', str(URUN_UST)],
                                        capture_output=True, text=True, timeout=sure_tavani)
                    if pr.returncode != 0:
                        sonuc[c['hedef']][dosya] = dict(hata=(pr.stderr or '')[:200]); continue
                except subprocess.TimeoutExpired:
                    sonuc[c['hedef']][dosya] = dict(hata='zaman asimi'); continue
                v = _spec_ayristir(_co, c.get('urun'))
                v['sure_hesap'] = round(time.time() - t0, 1)
                # kanit dosyalarini bagli klasore TOPLU yaz
                for _ek in ('', '.spec.tsv', '.mfe.log'):
                    if os.path.exists(_co + _ek):
                        try:
                            shutil.copyfile(_co + _ek, co + _ek)
                        except OSError as _e:
                            yaz(u'    UYARI: kanit dosyasi kopyalanamadi (%s): %s'
                                % (os.path.basename(co + _ek), _e))
            finally:
                shutil.rmtree(_yerel, ignore_errors=True)
            v['sure'] = round(time.time() - t0, 1)
            sonuc[c['hedef']][dosya] = v
            json.dump(v, open(kp, 'w', encoding='utf-8'), ensure_ascii=False, default=str)
        t = sum(x.get('hedef_disi', 0) for x in sonuc[c['hedef']].values()
                if isinstance(x, dict) and not x.get('hata'))
        a = sum(x.get('ayni_boyda', 0) for x in sonuc[c['hedef']].values()
                if isinstance(x, dict) and not x.get('hata'))
        yaz(u'    %-34s off-target %d, same length as the target %d' % (c['hedef'][:34], t, a))
    return dict(durum='TAMAM', kullanilan=kullanilan, atlanan=atlanan), sonuc


def yapi_kos(kok, mfe, ciftler, CIKTI, yaz):
    """hairpin + dimer. Panelin kendi kurallariyla (Tm < 45) degerlendirilir."""
    d = os.path.join(CIKTI, 'mfe')
    os.makedirs(d, exist_ok=True)
    fa = os.path.join(d, 'oligolar.fa')
    with open(fa, 'w', encoding='utf-8') as fh:
        for c in ciftler:
            ad = re.sub(r'\s+', '_', c['hedef'])
            fh.write('>%s_fp\n%s\n>%s_rp\n%s\n' % (ad, c['F'], ad, c['R']))
    out = {}
    for komut, ust in (('hairpin', HAIRPIN_TM_UST), ('dimer', DIMER_TM_UST)):
        co = os.path.join(d, '%s.txt' % komut)
        # D-10 (2026-08-07): spec_kos ile ayni tuzak - MFEprimer var olan cikti
        # dosyasinin uzerine YAZMAZ, 1 dondurur ("is already exists"). hairpin.txt
        # ve dimer.txt onceki kosulardan kaliyordu; duzeltilmezse bu iki olcum
        # her yeniden kosuda sessizce dict(hata=...) olurdu.
        for _es in (co, co + '.mfe.log'):
            if os.path.exists(_es):
                try:
                    os.remove(_es)
                except OSError:
                    pass
        if os.path.exists(co):
            out[komut] = dict(hata=u'eski cikti silinemedi, MFEprimer uzerine yazmaz: %s' % co)
            yaz(u'    ERROR: the old %s output could not be deleted, so it was not measured' % komut); continue
        try:
            pr = subprocess.run([mfe, komut, '-i', fa, '-o', co],
                                capture_output=True, text=True, timeout=600)
            if pr.returncode != 0:
                out[komut] = dict(hata=(pr.stderr or '')[:160]); continue
        except subprocess.TimeoutExpired:
            out[komut] = dict(hata='zaman asimi'); continue
        s = open(co, encoding='utf-8', errors='replace').read()
        m = re.search(r'%s List \((\d+)\)' % komut.capitalize(), s)
        n = int(m.group(1)) if m else 0
        # D-5 HATA DUZELTMESI (2026-08-06): eski kod 'Tm:\s*([0-9.]+)' ariyordu.
        # MFEprimer 3.0/4.4 dimer ve hairpin KAYITLARINDA Tm YOK - her kayit
        # 'Score: N, Delta G = -X.XX kcal/mol' yaziyor ('Tm' yalniz ustteki oligo
        # ozet tablosunun BASLIGINDA gecer, iki nokta olmadan). Sonuc: tms hep
        # bos kaliyor, en_yuksek_tm None, kural_ihlali her zaman False - yani bu
        # katman 66 dimer kaydina ragmen SESSIZCE "ihlal yok" oyu veriyordu.
        # Artik gercekten raporlanan olcu (Delta G) ayristirilir ve Tm'in
        # HESAPLANMADIGI acikca isaretlenir.
        dgs = [float(x) for x in re.findall(r'Delta G\s*=\s*(-?[0-9.]+)', s)]
        skor = [int(x) for x in re.findall(r'Score:\s*(\d+)', s)]
        tms = [float(x) for x in re.findall(r'Tm:\s*([0-9.]+)', s)]
        out[komut] = dict(sayi=n, en_yuksek_tm=(max(tms) if tms else None),
                          tm_hesaplanmadi=(not tms),
                          en_dusuk_dg=(min(dgs) if dgs else None),
                          en_yuksek_skor=(max(skor) if skor else None),
                          dg_esigi=DG_ESIGI,
                          kural_ihlali=bool(dgs and min(dgs) <= DG_ESIGI),
                          hukum_verilebilir=bool(dgs),
                          tm_kurali_uygulanamadi_ust=ust, dosya=co)
        yaz(u'    %s: %d records, lowest Delta G %s kcal/mol (this tool does NOT COMPUTE a Tm per record; Tm<%'
            % (komut, n, ('%.2f' % min(dgs)) if dgs else '-', ust))
    return out


def hedef_disi_kimlikleri(CIKTI, ciftler, yaz, tolerans=10):
    """D-7 (2026-08-06): 'hedef disi N amplikon' SAYISI tek basina hukum verdirmez.

    MFEprimer'in hedef disi olcusu SADECE BOYA dayanir: beklenen urun boyundan
    +-tolerans disinda kalan her amplikon 'hedef disi' sayilir. Ama grup ve
    EVRENSEL primerlerde hedef klad icindeki bir uyenin amplikonu da farkli
    boyda cikabilir (16S/18S'te dogal indel uzunluk polimorfizmi). Ornek:
    Bakteri_universal'in 31 'hedef disi' amplikonunun HEPSI bakteri 16S'idir
    (Thermodesulfovibrio, Desulfurella, Gemmata, Planctomycetes...) - evrensel
    bir bakteri primeri icin bunlar TASARIM GEREGI hedef ICIDIR.

    Bu yuzden sayinin yanina KIMLIK yazilir: hangi erisim numarasi, hangi
    organizma, kac bp. Karar bu dosyaya bakilarak verilir.
    """
    d = os.path.join(CIKTI, 'mfe')
    yol = os.path.join(CIKTI, 'mfe_hedef_disi_kimlikleri.tsv')
    # A4 DUZELTMESI (2026-08-21): burasi eskiden sessizce 'pass' idi.
    # 'urun' sayiya cevrilemezse o hedef 'bek' sozlugune HIC girmiyordu; sonuc
    # olarak beklenen_bp bilinmedigi icin o hedefin BUTUN amplikonlari hedef
    # disi sayiliyordu. Yani bozuk TEK bir hucre, o satirin hukmunu SESSIZCE
    # sertlestiriyordu - projenin bas hata turu: hata vermeden yanlis cevap.
    # Artik dusen her hedef ekrana ve dosyaya yazilir.
    bek = {}
    bek_dusen = []
    for c in ciftler:
        try:
            bek[re.sub(r'\W+', '_', c['hedef'])] = (c['hedef'], int(c['urun']))
        except (TypeError, ValueError, KeyError) as e:
            bek_dusen.append((c.get('hedef', '?'), repr(c.get('urun')), type(e).__name__))
    if bek_dusen:
        yaz(u'  WARNING: the expected product length could not be read for %d targets. For those, the "same length" separation CANNOT be made and every amplicon'
            % len(bek_dusen))
        for h, v, e in bek_dusen:
            yaz(u'    %-40s product=%s (%s)' % (h[:40], v, e))
    satir = 0
    with open(yol, 'w', encoding='utf-8', newline='') as fh:
        fh.write(u'# The IDENTITY of MFEprimer "off-target" amplicons (D-7).\n# The off-target measure rests on LENGTH ALONE; below')
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['hedef', 'veritabani', 'amplikon_bp', 'beklenen_bp', 'erisim',
                    'organizma', 'FpTm', 'RpTm', 'Ta', 'termodinamik_notu'])
        for dosya in sorted(os.listdir(d)) if os.path.isdir(d) else []:
            if not dosya.startswith('spec_') or not dosya.endswith('.txt'):
                continue
            ad = None
            for k in bek:
                if dosya.startswith('spec_%s_' % k) and (ad is None or len(k) > len(ad)):
                    ad = k
            if ad is None:
                continue
            hedef, b = bek[ad]
            db = dosya[len('spec_%s_' % ad):-4]
            s2 = open(os.path.join(d, dosya), encoding='utf-8', errors='replace').read()
            tab = {}
            for m in re.finditer(r'^\s*(\d+)\s+(\S+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+'
                                 r'[-\d.]+\s+[-\d.]+\s+[\d.]+\s+([\d.]+)\s*$', s2, re.M):
                if abs(int(m.group(3)) - b) > tolerans:
                    tab[int(m.group(1))] = (int(m.group(3)), m.group(2), m.group(4),
                                            m.group(5), m.group(6))
            if not tab:
                continue
            org = {}
            for m in re.finditer(r'^Amp (\d+): .*?==> (\S+?):(\d+)-(\d+) (.*)$', s2, re.M):
                org[int(m.group(1))] = m.group(5).strip()
            for i in sorted(tab):
                bp, er, ftm, rtm, ta = tab[i]
                try:
                    _not = (u'primer baglanma Tm (%s/%s C) tavlama sicakligindan (%s C) '
                            u'DUSUK - bu urun standart kosulda OLUSMAZ'
                            % (ftm, rtm, ta)) if min(float(ftm), float(rtm)) < float(ta) - 5 \
                        else u'primer Tm tavlama sicakligina yakin - GERCEK risk'
                except ValueError:
                    _not = u''
                w.writerow([hedef, db, bp, b, er, org.get(i, ''), ftm, rtm, ta, _not])
                satir += 1
    yaz(u'  written: %s (the identity of %d off-target amplicons)' % (yol, satir))
    return yol, satir


def hedef_disi_say(spec_sonuc, hedef_ad, ciftler=None):
    """MFEprimer'in bu cift icin bildirdigi TOPLAM amplikon sayisi.

    NOT: MFEprimer hedefin KENDI amplikonunu da sayar. Hedef disi olcusu icin
    'beklenen urun boyunda olmayanlar' degil, TOPLAM kullaniliyor ve bu acikca
    boyle raporlaniyor - iki arac ayni seyi saymadiginda celiski isareti bunu
    yakalasin diye. Ham dosyalar mfe/ altinda duruyor.
    """
    v = (spec_sonuc or {}).get(hedef_ad) or {}
    return sum(x.get('hedef_disi', 0) for x in v.values()
               if isinstance(x, dict) and not x.get('hata'))


# ---------------------------------------------------------------------------
# D-12 (2026-08-07): "hedef disi N amplikon" SAYISI HUKUM VERDIRMEZ.
#
# OLCULMUS GEREKCE: 2026-08-07 12:11 kosusunda MFEprimer 1605 farkli-boy
# amplikon buldu. mfe_hedef_disi_kimlikleri.tsv'deki taksonomi dizgeleri
# siniflandirildiginda dagilim sudur:
#     (a) hedef kladin ICINDEN, boyu farkli   1536   %95,7
#    (ao) hedef alan icinde ama ORGANEL         31   %1,9   (bitki mitokondri/kloroplast)
#     (b) ayni alan, klad disi                  24   %1,5
#     (c) farkli alan (gercek capraz)           14   %0,9
# Yani ham sayinin %95,7'si zararsiz uzunluk varyanti. Bakteri_universal'in
# 1550 "hedef disi"sinin 1519'u SILVA'da "Bacteria;..." ile basliyor - evrensel
# bir bakteri primeri icin TASARIM GEREGI hedef ICI.
#
# Bu fonksiyon hukme girecek olcuyu uretir: klad_disi = (b) + (c).
# Ayrica her (b)/(c) kaydinin primer Tm'i panelin baglanma sicakligiyla
# karsilastirilir; Tm belirgin altindaysa (>=TM_MARJI C) o urun standart
# kosulda olusmaz ve "olusabilir" sayimina girmez.
# ---------------------------------------------------------------------------
TM_MARJI = 5.0          # min(FpTm,RpTm) < Ta - TM_MARJI ise urun olusmaz sayilir
ORGANEL_JETONLARI = ('Chloroplast', 'Mitochondria')


def klad_tablosu(kok):
    """screening/hedef_klad.tsv -> {hedef: (alan, [klad jetonlari], kaynak)}"""
    y = os.path.join(kok, 'screening', 'hedef_klad.tsv')
    if not os.path.exists(y):
        return {}
    out = {}
    for r in csv.DictReader((l for l in open(y, encoding='utf-8')
                             if not l.startswith('#')), delimiter='\t'):
        h = (r.get('hedef') or '').strip()
        if not h:
            continue
        out[h] = ((r.get('alan') or '').strip(),
                  [x.strip() for x in (r.get('klad') or '').split(',') if x.strip()],
                  (r.get('kaynak') or '').strip())
    return out


# spec dosya adindaki veritabani -> alan. RefSeq basliklarinda taksonomi
# dizgesi YOKTUR; o kayitlarin alani veritabaninin TANIMINDAN gelir (bacteria
# .16S.fna icindeki her kayit tanimi geregi bakteri 16S'idir). Tahmin degil.
VTB_ALAN = {'archaea_16S_fna': 'Archaea', 'bacteria_16S_fna': 'Bacteria',
            'fungi_ITS_fna': 'Eukaryota', 'fungi_28SrRNA_fna': 'Eukaryota',
            'fungi_18SrRNA_fna': 'Eukaryota'}
VTB_KLAD = {'fungi_ITS_fna': 'Fungi', 'fungi_28SrRNA_fna': 'Fungi',
            'fungi_18SrRNA_fna': 'Fungi'}


def _kayit_coz(org, vtb):
    """(alan, jetonlar, organel_mi)"""
    if ';' in org and org.split(';')[0].strip() in ('Bacteria', 'Archaea', 'Eukaryota'):
        t = [x.strip() for x in org.split(';')]
        return (t[0], t, any(o in t for o in ORGANEL_JETONLARI))
    a = VTB_ALAN.get(vtb, '?')
    t = [org, a] + ([VTB_KLAD[vtb]] if vtb in VTB_KLAD else [])
    return (a, t, False)


def klad_siniflandir(kok, CIKTI, ciftler, ta_panel, yaz=None):
    """mfe_hedef_disi_kimlikleri.tsv'yi okur, cift basina a/ao/b/c sayar.

    Doner: {hedef: dict(a, ao, b, c, klad_disi, olusabilir, olusmaz,
                        klad, alan, kaynak, kayitlar=[...])}
    'klad_disi' hukme girecek olcudur. 'olusabilir', o kayitlardan primer
    Tm'i baglanma sicakligina yakin olanlarin sayisidir.
    Tablo ya da dosya yoksa BOS doner - o zaman cagiran taraf ham sayiyi
    kullanmaya DEVAM eder ama bunu acikca "klad ayrimi YAPILAMADI" diye
    isaretlemek zorundadir (sessizce temiz saymak yasak).
    """
    tab = klad_tablosu(kok)
    yol = os.path.join(CIKTI, 'mfe_hedef_disi_kimlikleri.tsv')
    if not tab or not os.path.exists(yol):
        if yaz:
            yaz(u'    clade separation NOT POSSIBLE (%s missing); the raw count will be used'
                % ('hedef_klad.tsv' if not tab else os.path.basename(yol)))
        return {}
    out = {}
    for c in ciftler:
        h = c['hedef']
        if h in tab:
            out[h] = dict(a=0, ao=0, b=0, c=0, klad_disi=0, olusabilir=0,
                          olusmaz=0, alan=tab[h][0], klad=tab[h][1],
                          kaynak=tab[h][2], kayitlar=[])
    eksik = set()
    for r in csv.DictReader((l for l in open(yol, encoding='utf-8')
                             if not l.startswith('#')), delimiter='\t'):
        h = r.get('hedef')
        if h not in out:
            if h:
                eksik.add(h)
            continue
        d = out[h]
        k_alan, jet, organel = _kayit_coz(r.get('organizma') or '', r.get('veritabani') or '')
        ic = any(j in jet for j in d['klad'])
        if not ic and ';' not in (r.get('organizma') or ''):
            ic = any(re.search(r'\b%s' % re.escape(j), r.get('organizma') or '', re.I)
                     for j in d['klad'])
        if ic and organel:
            s = 'ao'
        elif ic:
            s = 'a'
        elif k_alan == d['alan']:
            s = 'b'
        else:
            s = 'c'
        d[s] += 1
        if s in ('b', 'c'):
            d['klad_disi'] += 1
        if s in ('b', 'c', 'ao'):
            try:
                dtm = min(float(r['FpTm']), float(r['RpTm'])) - float(ta_panel)
            except (ValueError, KeyError, TypeError):
                dtm = None
            if dtm is None:
                pass
            elif dtm < -TM_MARJI:
                d['olusmaz'] += 1
            else:
                d['olusabilir'] += 1
            d['kayitlar'].append(dict(sinif=s, bp=r.get('amplikon_bp'),
                                      erisim=r.get('erisim'),
                                      organizma=r.get('organizma'),
                                      FpTm=r.get('FpTm'), RpTm=r.get('RpTm'),
                                      dTm=dtm))
    if eksik and yaz:
        yaz(u'    WARNING: target(s) with no entry in hedef_klad.tsv: %s'
            % ', '.join(sorted(eksik))[:200])
    return out
