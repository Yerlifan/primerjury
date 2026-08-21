#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
batch_design.py
hedefler.tsv dosyasındaki bütün toplantı kararlarını design_group_primers.py ile
sırayla tasarlar ve tek bir özet tablo üretir.

Tasarım her AMPLİKON SINIFI için ayrı yapılır (A1, A2, B, F1, F2), çünkü
sınıflar rDNA'nın farklı pencerelerini kapsıyor. Bir hedef birden fazla
sınıfta bulunuyorsa her biri ayrı satır olur ve en iyi sonuç raporlanır.

Rakip kümesi: aynı amplikon sınıfındaki, hedefte olmayan ve 'haric'
sütununda belirtilmeyen bütün taksonlar. Hedefin kendi üyeleri asla
rakip listesine girmez.

Kullanım:
  python3 batch_design.py \
      --kons "/.../referans_konsensus/self/konsensus" \
      --hedefler hedefler.tsv \
      --out "/.../primer_adaylari" \
      [--only Karar1] [--jobs 1] [--extra "--degeneracy-budget 2"]
"""
import argparse, csv, glob, hashlib, json, os, re, subprocess, sys, collections, time, datetime
import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
MOTOR = os.path.join(HERE, "design_group_primers.py")


def _ayirt_modulu():
    """indistinguishable_targets.py'yi tek kaynak olarak yukler."""
    yol = os.path.join(HERE, "indistinguishable_targets.py")
    if not os.path.exists(yol):
        return None
    spec = importlib.util.spec_from_file_location("ayirt10", yol)
    m = importlib.util.module_from_spec(spec)
    yedek, sys.argv = sys.argv, [yol]
    try:
        spec.loader.exec_module(m)
    finally:
        sys.argv = yedek
    return m


AYIRT = _ayirt_modulu()

# Arama merdiveni. Sirayla denenir, ilk aday veren kademede durulur.
#   kati        : toplanti kararindaki kural, yetim primer rakiplerde HIC
#                 baglanmamali
#   yetim3      : yetim primerin rakiplerdeki EN IYI yerlesimi bile en az uc
#                 uyumsuzluk tasiyorsa yeterli sayilir
#   yetim3_genis: ayni gevseme, ustune daha genis oligo taramasi
MERDIVEN = [
    ("kati",         []),
    ("yetim3",       ["--yetim-min-uyumsuzluk", "3"]),
    ("yetim3_genis", ["--yetim-min-uyumsuzluk", "3", "--max-oligo", "1200"]),
]


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--kons", required=True,
                   help="consensus directory (output of the freeze-reference or anchored-consensus step)")
    p.add_argument("--hedefler", default=os.path.join(HERE, "hedefler.tsv"))
    p.add_argument("--out", required=True)
    p.add_argument("--only", default=None, help="this target name only")
    p.add_argument("--karar", default=None, help="this decision group only")
    p.add_argument("--extra", default="", help="04'e flags passed through unchanged")
    p.add_argument("--grup-extra", default="",
                   help="passed in addition for decision-group 3 and 4 targets "
                        "bayraklar. Dejenerelik bayraklari kaldirildi: "
                        "toplanti karari butun hedeflerde salt ACGT oligo "
                        "istiyor, kalip belirsizligi --iupac-max ile "
                        "yonetiliyor.")
    p.add_argument("--min-uye", type=int, default=1)
    p.add_argument("--ayirt-edilemez-cikar", type=int, default=1,
                   help="1: too close to one of the target members to be distinguished "
                        "ozdes olan taksonlar rakip listesinden cikarilir ve "
                        "gerekcesi loglanir. 0: cikarilmaz (toplanti kuralinin "
                        "mantiken saglanamadigi durumda hedef sifir aday verir)")
    p.add_argument("--yeniden", action="store_true",
                   help="ignore the checkpoint, redesign every target from scratch")
    return p.parse_args()


def _girdi_parmak_izi(a):
    """Kosunun bagli oldugu her girdinin ozeti. Degisirse checkpoint duser."""
    h = hashlib.sha256()
    h.update(("kons=%s\n" % os.path.abspath(a.kons)).encode())
    for f in sorted(glob.glob(os.path.join(a.kons, "*_konsensus.fasta"))):
        try:
            st = os.stat(f)
            h.update(("%s|%d|%d\n" % (os.path.basename(f), st.st_size,
                                      int(st.st_mtime))).encode())
        except OSError:
            pass
    for yol in (a.hedefler, MOTOR, os.path.join(HERE, "generate_primer_candidates.py"),
                os.path.join(HERE, "indistinguishable_targets.py")):
        try:
            st = os.stat(yol)
            h.update(("%s|%d|%d\n" % (os.path.basename(yol), st.st_size,
                                      int(st.st_mtime))).encode())
        except OSError:
            h.update(("%s|yok\n" % yol).encode())
    h.update(("extra=%s|grup_extra=%s|min_uye=%d|ayirt=%d\n"
              % (a.extra, a.grup_extra, a.min_uye,
                 a.ayirt_edilemez_cikar)).encode())
    return h.hexdigest()[:16]


TS = lambda: datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
_LOG = [None]


def log(msg):
    """Her satir tarih ve saatle hem ekrana hem log dosyasina yazilir."""
    line = "[%s] %s" % (TS(), msg)
    print(line, flush=True)
    if _LOG[0]:
        with open(_LOG[0], "a", encoding="utf-8") as fh:
            fh.write(line + "\n")


def sinif_of(tag):
    """A1-1_2209 -> A1 ; F2-3_101201 -> F2 ; B-4_818 -> B"""
    grp = tag.rsplit("_", 1)[0]
    base = re.split(r"[-_]", grp)[0]
    return base


def taxid_of(tag):
    m = re.search(r"_(\d+)$", tag)
    return m.group(1) if m else None


def main():
    a = get_args()
    os.makedirs(a.out, exist_ok=True)
    _LOG[0] = os.path.join(a.out, "toplu_tasarim.log")
    CKPT = os.path.join(a.out, "checkpoint.json")
    # Checkpoint'in gecerliligi GIRDIYE baglidir. Konsensus kumesi, hedef
    # tablosu ya da motor betikleri degistiyse eski sonuclar artik ayni
    # soruya ait degildir; sessizce yeniden kullanilirsa kosu bitmis gibi
    # gorunur ama sonuc eski girdiye aittir. Bu yuzden girdi parmak izi
    # checkpoint'e yazilir ve her aciliata karsilastirilir.
    parmak = _girdi_parmak_izi(a)
    ckpt = {}
    if os.path.exists(CKPT) and not a.yeniden:
        try:
            ham = json.load(open(CKPT, encoding="utf-8"))
            eski = ham.get("_girdi_parmak_izi") if isinstance(ham, dict) else None
            kayitlar = {k: v for k, v in ham.items()
                        if not k.startswith("_")} if isinstance(ham, dict) else {}
            if eski is None:
                log(u'the checkpoint CARRIES NO fingerprint (left over from an old version), it is ignored and every target will be designed from scratch')
            elif eski != parmak:
                log(u'the checkpoint fingerprint DOES NOT MATCH, it is ignored')
                log(u'   recorded: %s' % eski)
                log(u'   current : %s' % parmak)
                log(u'   The input has changed (the consensus set, hedefler.tsv or the engine scripts). The old results will not be reused.')
            else:
                ckpt = kayitlar
                log(u'the checkpoint was found and the input fingerprint matches: %d target and class combinations are already done, they will be skipped' % len(ckpt))
        except Exception as e:
            log(u'the checkpoint could not be read (%s), starting from scratch' % e)
            ckpt = {}
    t0 = time.time()
    log(u'start. consensus=%s  output=%s' % (a.kons, a.out))
    files = sorted(glob.glob(os.path.join(a.kons, "*_konsensus.fasta")))
    if not files:
        sys.exit(u'no consensus found: %s' % a.kons)
    envanter = {}
    bozuk = []
    for f in files:
        tag = re.sub(r"_[A-Za-z]+_konsensus\.fasta$", "", os.path.basename(f))
        tid = taxid_of(tag)
        if not tid:
            continue
        # Bos ya da tumu N olan konsensus tasarima giremez. 07 self kosusunda
        # taxid 2233851 icin uzunluk 0 cikti; boyle bir uye hedef kumesine
        # girerse hicbir oligo ona baglanamaz ve butun hedef coker, rakip
        # kumesine girerse ozgulluk denetimi sessizce zayiflar.
        seq = "".join(l.strip() for l in open(f, encoding="utf-8",
                                              errors="replace")
                      if not l.startswith(">")).upper()
        kapsanan = sum(1 for c in seq if c != "N")
        if kapsanan < 200:
            bozuk.append((tag, len(seq), kapsanan))
            continue
        envanter[tag] = dict(path=f, taxid=tid, sinif=sinif_of(tag),
                             uzunluk=len(seq), kapsanan=kapsanan)
    if not bozuk:
        with open(os.path.join(a.out, "dislanan_takson.tsv"), "w",
                  encoding="utf-8") as df:
            df.write("grup\ttaxid\tetiket\tuzunluk\tkapsanan\n")
    siniflar = collections.defaultdict(list)
    for tag, d in envanter.items():
        siniflar[d["sinif"]].append(tag)
    if bozuk:
        log(u'EXCLUDED consensus (fewer than 200 covered bases): %d' % len(bozuk))
        for t, L, k in bozuk:
            log("   %-26s uzunluk=%d kapsanan=%d" % (t, L, k))
        # 09 da ayni taksonlari dislamali; fastq envanterinden calistigi
        # icin bunlari kendiliginden bilemez.
        with open(os.path.join(a.out, "dislanan_takson.tsv"), "w",
                  encoding="utf-8") as df:
            df.write("grup\ttaxid\tetiket\tuzunluk\tkapsanan\n")
            for t, L, k in bozuk:
                df.write("%s\t%s\t%s\t%d\t%d\n"
                         % (t.rsplit("_", 1)[0], taxid_of(t) or "", t, L, k))
        log(u'   dislanan_takson.tsv was written (specificity.py will exclude the same set)')
        log(u'   These files will be used neither as a target nor as a competitor.')
        log(u'   The reason can be checked in referans_konsensus/self/log/<label>_mm2.log')
        log(u'   for each excluded label.')
    log(u'consensus: %d files, %d classes (%s)'
        % (len(envanter), len(siniflar), ", ".join(sorted(siniflar))))
    for s in sorted(siniflar):
        tids = sorted(set(envanter[t]["taxid"] for t in siniflar[s]), key=int)
        log(u'   %-3s %2d files, %2d taxa: %s'
            % (s, len(siniflar[s]), len(tids), ",".join(tids)))

    # --- ayirt edilemez takson ciftleri ------------------------------
    # Kraken2 tek bir populasyonun okumalarini kardes tur dugumlerine
    # dagitabiliyor. Boyle iki kutu ayni konsensusu uretir. Biri rakip
    # listesine girerse "rakipte urun olusmasin" kurali MANTIKEN
    # saglanamaz ve hedef sessizce sifir aday verir. Bu yuzden once
    # olculur, sonra rakip listesinden cikarilir ve her cikarma loglanir.
    ayirt = {}
    if a.ayirt_edilemez_cikar and AYIRT is not None:
        temsil = {}
        for tag, d in envanter.items():
            key = (d["sinif"], d["taxid"])
            if key not in temsil or d["kapsanan"] > temsil[key][1]:
                seq = "".join(l.strip() for l in open(d["path"], encoding="utf-8",
                                                      errors="replace")
                              if not l.startswith(">")).upper()
                temsil[key] = (tag, d["kapsanan"], seq)
        if not AYIRT.MAPPY:
            log(u'WARNING: mappy is not installed. The indistinguishable bin measurement will be done with k-mer coverage only, with no alignment measurement. To install it: pip install mappy')
        ciftler = AYIRT.ayirt_edilemezler(temsil)
        for sn, t1, t2, u, o, k, g, kp, kt in ciftler:
            ayirt.setdefault((sn, t1), set()).add(t2)
            ayirt.setdefault((sn, t2), set()).add(t1)
        log(u'indistinguishable taxon pairs: %d' % len(ciftler))
        for sn, t1, t2, u, o, k, g, kp, kt in sorted(ciftler, key=lambda r: (r[0], -r[4])):
            log(u'   %-3s %-9s ~ %-9s aligned=%5d intersect=%%%.2f strict=%%%.2f coverage=%.2f kmer=%.3f  %s' % (sn, t1, t2, u, o, kt, kp, k, g))
        with open(os.path.join(a.out, "ayirt_edilemez.tsv"), "w",
                  newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["sinif", "taxid1", "taxid2", "hizalanan_bp",
                        "kesisimli_ozdeslik", "kmer_kapsamasi", "gerekce",
                        "hizalama_kapsami", "kati_ozdeslik"])
            w.writerows(ciftler)
    elif a.ayirt_edilemez_cikar:
        log(u'WARNING: indistinguishable_targets.py was not found, the indistinguishable bin cleanup WILL NOT BE DONE')

    hedefler = []
    with open(a.hedefler, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if p[0] == "karar":
                continue
            hedefler.append(dict(karar=p[0], hedef=p[1], duzey=p[2],
                                 inn=p[3], haric=p[4] if len(p) > 4 else "",
                                 note=p[5] if len(p) > 5 else ""))
    if a.only:
        hedefler = [h for h in hedefler if h["hedef"] == a.only]
    if a.karar:
        hedefler = [h for h in hedefler if h["karar"] == a.karar]
    log(u'targets to process: %d' % len(hedefler))

    ozet = []
    for h in hedefler:
        spec = h["inn"]
        haric = set(x for x in h["haric"].split(",") if x)
        # alan hedefleri: *B, *A, *F  -> o harfle baslayan butun siniflar
        if spec.startswith("*"):
            harf = spec[1:]
            hedef_siniflar = [s for s in siniflar if s.startswith(harf)]
            in_tids = None            # sinif icindeki her sey
        else:
            in_tids = set(x for x in spec.split(",") if x)
            hedef_siniflar = sorted(set(
                envanter[t]["sinif"] for t in envanter if envanter[t]["taxid"] in in_tids))
            # ALAN TUTARLILIGI: hedefin taxidleri, ait olmadigi bir lokus
            # kitapliginda da gecebiliyor (ornegin bakteriyel bir taxid mantar
            # ITS/28S kutusunda). O siniftaki tasarim hedefi temsil etmez;
            # dahasi, dogru siniftaki tasarim cift bulamazsa geriye yalniz bu
            # sahte cift kalir ve hedef kapsanmis gorunur. Alan bilgisi elle
            # yazilmaz, kutu dagilimindan cikarilir; hedef tek alandaysa
            # hicbir sey degismez.
            _alan = {}
            for t in envanter:
                if envanter[t]["taxid"] in in_tids:
                    _a = envanter[t]["sinif"][0]
                    _alan[_a] = _alan.get(_a, 0) + 1
            if len(_alan) > 1:
                _bas = max(_alan, key=lambda x: _alan[x])
                _dis = [s for s in hedef_siniflar if s[0] != _bas]
                if _dis:
                    hedef_siniflar = [s for s in hedef_siniflar if s[0] == _bas]
                    log(u'   %-34s a mixture of domains (%s); the dominant domain is %s, the class skipped: %s'
                        % (h["hedef"],
                           ", ".join("%s=%d" % kv for kv in sorted(_alan.items())),
                           _bas, ", ".join(_dis)))
        if not hedef_siniflar:
            log(u'SKIPPED %-34s it has no counterpart in any class in the data' % h["hedef"])
            ozet.append(dict(h, sinif="-", uye=0, rakip=0, cift=0, durum="veride yok",
                             en_iyi="", tsv=""))
            continue

        for sinif in hedef_siniflar:
            tags = siniflar[sinif]
            if in_tids is None:
                ing = [t for t in tags if envanter[t]["taxid"] not in haric]
                outg = []
                # alan hedefi: rakip = diger alanlarin siniflari
                for s2 in siniflar:
                    if s2[0] != sinif[0]:
                        outg.extend(siniflar[s2])
            else:
                ing = [t for t in tags if envanter[t]["taxid"] in in_tids]
                outg = [t for t in tags if envanter[t]["taxid"] not in in_tids
                        and envanter[t]["taxid"] not in haric]
            in_tid_kume = set(envanter[t]["taxid"] for t in ing)
            atilan = []
            if ayirt:
                temiz = []
                for t in outg:
                    tx = envanter[t]["taxid"]
                    carpisan = ayirt.get((sinif, tx), set()) & in_tid_kume
                    if carpisan:
                        atilan.append((tx, sorted(carpisan)))
                    else:
                        temiz.append(t)
                outg = temiz
            if atilan:
                gor = sorted(set(t[0] for t in atilan))
                log(u'   %-34s %-3s taken out of the competitors (indistinguishable from the target): %s'
                    % (h["hedef"], sinif, ", ".join(gor)))
            if len(ing) < a.min_uye:
                log(u'SKIPPED %-34s %-3s members=%d (below --min-uye %d)'
                    % (h["hedef"], sinif, len(ing), a.min_uye))
                ozet.append(dict(h, sinif=sinif, uye=len(ing), rakip=len(outg),
                                 cift=0, durum="uye yetersiz", kademe="",
                                 en_iyi="", tsv="", engelleyen="",
                                 sure=0, zaman=TS()))
                continue
            etiket = "%s__%s" % (h["hedef"], sinif)
            tsv = os.path.join(a.out, "%s.tsv" % etiket)
            if etiket in ckpt and not a.yeniden:
                c = ckpt[etiket]
                log(u'SKIPPED (checkpoint) %-40s %-3s pairs=%s' % (h["hedef"], sinif, c.get("cift")))
                ozet.append(dict(h, sinif=sinif, uye=c.get("uye", 0),
                                 rakip=c.get("rakip", 0),
                                 cift=c.get("cift", 0), durum=c.get("durum", ""),
                                 kademe=c.get("kademe", ""),
                                 en_iyi=c.get("en_iyi", ""), tsv=c.get("tsv", ""),
                                 engelleyen=c.get("engelleyen", ""),
                                 sure=c.get("sure", ""), zaman=c.get("zaman", "")))
                continue
            cmd = [sys.executable, MOTOR,
                   "--in-group"] + [envanter[t]["path"] for t in sorted(ing)]
            if outg:
                cmd += ["--out-group"] + [envanter[t]["path"] for t in sorted(outg)]
            cmd += ["--label", etiket, "--out", tsv]
            if a.extra:
                cmd += a.extra.split()
            if h["karar"] in ("3", "4") and a.grup_extra:
                cmd += a.grup_extra.split()
            # Kademeli arama. Once toplanti kararindaki KATI kural denenir:
            # primerlerden biri rakiplerde hic baglanmamali. Bu saglanamazsa
            # kural bir kademe gevsetilir ve gevseme cikti tablosunda ACIKCA
            # isaretlenir; boylece "aday yok" ile "aday var ama daha zayif
            # guvenceyle" ayrimi kaybolmaz.
            th = time.time()
            n, txt, kademe = 0, "", "kati"
            for kad_ad, kad_bayrak in MERDIVEN:
                r = subprocess.run(cmd + kad_bayrak, capture_output=True, text=True)
                txt = r.stdout + r.stderr
                m = re.search(r"gecerli cift sayisi\s*:\s*(\d+)", txt)
                n = int(m.group(1)) if m else 0
                kademe = kad_ad
                if n > 0:
                    break
                if kad_ad != MERDIVEN[-1][0]:
                    log(u'      %-40s %-3s the \'%s\' step gave zero, relaxing'
                        % (h["hedef"], sinif, kad_ad))
            sure = time.time() - th
            if r.returncode != 0 and n == 0:
                log(u'      ERROR: design_group_primers.py finished with exit code %d. The last lines:'
                    % r.returncode)
                for satir in txt.strip().splitlines()[-4:]:
                    log("         %s" % satir[:150])
            # Sifir cift veren hedefte ONCEKI kosudan kalan aday dosyasi
            # silinmeli; yoksa 09 bayat primerleri dogrular ve teslimata
            # gecersiz cift girer.
            if n == 0 and os.path.exists(tsv):
                try:
                    os.remove(tsv)
                    log(u'      the old candidate file was deleted: %s'
                        % os.path.basename(tsv))
                except OSError as e:
                    log(u'      WARNING: the old candidate file could not be deleted (%s)' % e)
            best = ""
            mb = re.search(r"En iyi bes aday:\n\s*(.+)", txt)
            if mb:
                best = mb.group(1).strip()
            durum = "TAMAM" if n > 0 else "cift yok"
            if n > 0 and kademe != "kati":
                durum = "TAMAM(%s)" % kademe
            if "gecerli cift bulunamadi" in txt and n == 0:
                durum = "cift yok"
            eng = ""
            me = re.search(r"en cok engelleyen uyeler:\s*(.+)", txt)
            if me:
                eng = me.group(1).strip()[:70]
            log(u'%-40s %-3s members=%2d competitors=%2d pairs=%-7d %-14s %5.1f s'
                % (h["hedef"], sinif, len(ing), len(outg), n, durum, sure))
            if eng and n == 0:
                log("      engelleyen: %s" % eng)
            with open(os.path.join(a.out, "%s.log" % etiket), "w",
                      encoding="utf-8") as lf:
                lf.write(txt)
            kayit = dict(uye=len(ing), rakip=len(outg), cift=n, durum=durum,
                         kademe=kademe if n else "",
                         en_iyi=best, tsv=tsv if n else "", engelleyen=eng,
                         sure=round(sure, 1), zaman=TS())
            ozet.append(dict(h, sinif=sinif, **kayit))
            ckpt[etiket] = kayit
            with open(CKPT, "w", encoding="utf-8") as cf:
                json.dump(dict(ckpt, _girdi_parmak_izi=parmak), cf,
                          ensure_ascii=False, indent=1)
            yaz_ozet(a.out, ozet, kismi=bool(a.only or a.karar))

    yaz_ozet(a.out, ozet, kismi=bool(a.only or a.karar))
    tamam = sum(1 for o in ozet if o["durum"].startswith("TAMAM"))
    log(u'summary written: %s' % os.path.join(a.out, "ozet.tsv"))
    log(u'target and class combinations: %d, with a pair found: %d' % (len(ozet), tamam))
    kalan = sorted(set(o["hedef"] for o in ozet if o["durum"] != "TAMAM")
                   - set(o["hedef"] for o in ozet if o["durum"].startswith("TAMAM")))
    if kalan:
        log(u'targets with no pair found in any class (%d): %s'
            % (len(kalan), ", ".join(kalan)))
    log(u'total time: %.1f minutes' % ((time.time() - t0) / 60))


def yaz_ozet(out, ozet, kismi=False):
    """Ozet her hedeften sonra yeniden yazilir; kosu yarida kesilse de
    o ana kadarki sonuclar diskte kalir.

    kismi=True (--only ya da --karar ile calisildiginda) mevcut ozet.tsv
    okunup bu kosuda islenmeyen satirlar korunur; aksi halde kismi bir kosu
    butun raporu tek hedefe indiriyordu."""
    cols = ["karar", "hedef", "duzey", "sinif", "uye", "rakip", "cift",
            "durum", "kademe", "en_iyi", "engelleyen", "sure", "zaman",
            "tsv", "note"]
    satirlar = list(ozet)
    yol = os.path.join(out, "ozet.tsv")
    if kismi and os.path.exists(yol):
        gorulen = set((x.get("hedef"), x.get("sinif")) for x in ozet)
        try:
            for eski in csv.DictReader(open(yol, encoding="utf-8"),
                                       delimiter="\t"):
                if (eski.get("hedef"), eski.get("sinif")) not in gorulen:
                    satirlar.append(eski)
        except OSError:
            pass
    with open(yol, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(satirlar)


if __name__ == "__main__":
    main()
