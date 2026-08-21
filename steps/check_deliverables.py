#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_deliverables.py
TESLIM EDILEN TABLONUN BAGIMSIZ DENETIMI.

Bu betik tasarim kodunun hicbir fonksiyonunu ice aktarmaz. primer_final.tsv
dosyasini okur ve her satiri toplanti kararindaki kurallara gore SIFIRDAN
yeniden olcer. Amac, tasarim kodu ile teslim edilen tablo arasinda hicbir
sessiz kayma kalmadigini gostermektir: tabloda yazan Tm, GC, urun boyu ve
tm_farki degerleri de yeniden hesaplanip karsilastirilir.

Denetlenen kurallar
  Oligo   : yalniz A/C/G/T, uzunluk 18-25, GC %40-60 (sert 35-65),
            3' uc G ya da C, son bes bazda en fazla 3 G/C,
            en fazla 4 ardisik ayni baz
  Termo   : Tm 58-62 (sert 57-63) iki kutuphaneyle, hairpin >= -3000,
            self-dimer >= -6000, hetero-dimer >= -6000
  Cift    : |TmF - TmR| < 1,5 ; urun 70-250 (sert 300)
  Tablo   : yazili ileri_tm/geri_tm/ileri_gc/geri_gc/tm_farki degerleri
            yeniden hesaplananla ayni mi
  Kalip   : (--kons verilirse) geri primer gercekten kalibin ters
            tumleyeni mi, urun ileri primerle baslayip geri primerin ters
            tumleyeniyle bitiyor mu

Kullanim:
  python3 check_deliverables.py --final pr_final --kons pr_kons/konsensus
"""
import argparse, csv, glob, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import field_audit

try:
    import primer3
except ImportError:
    sys.exit("primer3-py gerekli: pip install primer3-py --break-system-packages")
try:
    from Bio.SeqUtils import MeltingTemp as mt
    from Bio.Seq import Seq
except ImportError:
    sys.exit("biopython gerekli: pip install biopython --break-system-packages")

TAM = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(TAM)[::-1]


def gc_yuzde(s):
    return 100.0 * (s.count("G") + s.count("C")) / len(s) if s else 0.0


def en_uzun_tekrar(s):
    en = k = 1
    for i in range(1, len(s)):
        k = k + 1 if s[i] == s[i - 1] else 1
        if k > en:
            en = k
    return en


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--final", required=True, help="09'un output directory")
    p.add_argument("--kons", default=None,
                   help="baskin alel consensus directory; if given the template is checked as well")
    p.add_argument("--hedefler", default=None,
                   help="targets.tsv; if given field consistency is checked as well")
    p.add_argument("--yalniz-gecti", type=int, default=1)
    p.add_argument("--out", default=None, help="findings TSV'si")
    # toplanti karari esikleri
    p.add_argument("--uz-min", type=int, default=18)
    p.add_argument("--uz-max", type=int, default=25)
    p.add_argument("--gc-min", type=float, default=40.0)
    p.add_argument("--gc-max", type=float, default=60.0)
    p.add_argument("--gc-sert-min", type=float, default=35.0)
    p.add_argument("--gc-sert-max", type=float, default=65.0)
    p.add_argument("--gc-son-pencere", type=int, default=5)
    p.add_argument("--gc-son-max", type=int, default=3)
    p.add_argument("--homopolimer-max", type=int, default=4)
    p.add_argument("--tm-min", type=float, default=58.0)
    p.add_argument("--tm-max", type=float, default=62.0)
    p.add_argument("--tm-sert-min", type=float, default=57.0)
    p.add_argument("--tm-sert-max", type=float, default=63.0)
    p.add_argument("--tm-fark-max", type=float, default=1.5)
    p.add_argument("--tm-capraz-tol", type=float, default=2.0)
    p.add_argument("--hairpin-min", type=float, default=-3000.0)
    p.add_argument("--selfdimer-min", type=float, default=-6000.0)
    p.add_argument("--heterodimer-min", type=float, default=-6000.0)
    p.add_argument("--urun-min", type=int, default=70)
    p.add_argument("--urun-max", type=int, default=250)
    p.add_argument("--urun-sert-max", type=int, default=300)
    # termodinamik kosullar, 03/04 ile ayni
    p.add_argument("--mv", type=float, default=50.0)
    p.add_argument("--dv", type=float, default=1.5)
    p.add_argument("--dntp", type=float, default=0.6)
    # 03 ve 04'un --dna-conc varsayilani ile ayni olmak zorunda; farkli olursa
    # yeniden olculen Tm sistematik olarak kayar ve sahte bulgu uretir.
    p.add_argument("--dna", type=float, default=50.0)
    return p.parse_args()


def oku_fasta(f):
    return "".join(l.strip() for l in open(f, encoding="utf-8", errors="replace")
                   if not l.startswith(">")).upper()


def kalip_yukle(kok):
    d = {}
    if not kok:
        return d
    for p in sorted(glob.glob(os.path.join(kok, "*.fasta"))):
        et = re.sub(r"_(baskin|ref|self)?_?konsensus\.fasta$", "",
                    os.path.basename(p))
        d[et] = oku_fasta(p)
    return d


def main():
    a = get_args()
    tsv = os.path.join(a.final, "primer_final.tsv")
    if not os.path.exists(tsv):
        sys.exit("bulunamadi: %s" % tsv)
    rows = list(csv.DictReader(open(tsv, encoding="utf-8"), delimiter="\t"))
    if a.yalniz_gecti:
        rows = [r for r in rows if r.get("ozgulluk_durum") == "GECTI"]
    if not rows:
        sys.exit("denetlenecek satir yok")

    kaliplar = kalip_yukle(a.kons)
    print("=" * 72)
    print("TESLIM DENETIMI  (bagimsiz yeniden olcum)")
    print("=" * 72)
    print("denetlenen satir : %d" % len(rows))
    print("kalip dosyasi    : %d" % len(kaliplar))
    print("termodinamik     : Na=%.0f Mg=%.1f dNTP=%.1f DNA=%.0f"
          % (a.mv, a.dv, a.dntp, a.dna))
    print()

    def tm3(s):
        return primer3.calc_tm(s, mv_conc=a.mv, dv_conc=a.dv,
                               dntp_conc=a.dntp, dna_conc=a.dna)

    def tmb(s):
        return float(mt.Tm_NN(Seq(s), nn_table=mt.DNA_NN3, Na=a.mv, Mg=a.dv,
                              dNTPs=a.dntp, dnac1=a.dna, dnac2=0, saltcorr=7))

    bulgu = []          # (agirlik, hedef, sinif, kural, ayrinti)

    def ekle(ag, r, kural, ayrinti):
        bulgu.append(dict(agirlik=ag, hedef=r.get("hedef", ""),
                          sinif=r.get("sinif", ""), kural=kural,
                          ayrinti=ayrinti,
                          ileri=r.get("ileri_dizi", ""),
                          geri=r.get("geri_dizi", "")))

    # sistematik kayma once olculur, sonra ondan sapan aranir
    tumu = []
    for r in rows:
        for s in (r["ileri_dizi"], r["geri_dizi"]):
            tumu.append(s)
    farklar = sorted(tm3(s) - tmb(s) for s in tumu)
    kayma = farklar[len(farklar) // 2]
    print("iki kutuphane arasi medyan kayma: %+.2f C  (tolerans %.2f C)"
          % (kayma, a.tm_capraz_tol))

    for r in rows:
        F, R = r["ileri_dizi"].upper(), r["geri_dizi"].upper()
        for ad, s in (("ileri", F), ("geri", R)):
            if set(s) - set("ACGT"):
                ekle("KRITIK", r, "alfabe",
                     "%s primerde ACGT disi baz: %s" % (ad, sorted(set(s) - set("ACGT"))))
                continue
            if not (a.uz_min <= len(s) <= a.uz_max):
                ekle("KRITIK", r, "uzunluk", "%s uzunluk %d" % (ad, len(s)))
            g = gc_yuzde(s)
            if not (a.gc_sert_min <= g <= a.gc_sert_max):
                ekle("KRITIK", r, "gc_sert", "%s GC %%%.1f" % (ad, g))
            elif not (a.gc_min <= g <= a.gc_max):
                ekle("UYARI", r, "gc_tercih", "%s GC %%%.1f" % (ad, g))
            if s[-1] not in "GC":
                ekle("KRITIK", r, "3p_gc_kilit", "%s 3' uc %s" % (ad, s[-1]))
            kuyruk = s[-a.gc_son_pencere:]
            if kuyruk.count("G") + kuyruk.count("C") > a.gc_son_max:
                ekle("KRITIK", r, "3p_asiri_sabit",
                     "%s son %d bazda %d G/C" % (ad, a.gc_son_pencere,
                                                 kuyruk.count("G") + kuyruk.count("C")))
            if en_uzun_tekrar(s) > a.homopolimer_max:
                ekle("KRITIK", r, "homopolimer",
                     "%s %d ardisik ayni baz" % (ad, en_uzun_tekrar(s)))
            t3, tb = tm3(s), tmb(s)
            if abs((t3 - tb) - kayma) > a.tm_capraz_tol:
                ekle("KRITIK", r, "tm_capraz",
                     "%s primer3 %.2f, biopython %.2f, kaymadan sapma %.2f"
                     % (ad, t3, tb, abs((t3 - tb) - kayma)))
            if not (a.tm_sert_min <= t3 <= a.tm_sert_max):
                ekle("KRITIK", r, "tm_sert", "%s Tm %.2f" % (ad, t3))
            elif not (a.tm_min <= t3 <= a.tm_max):
                ekle("UYARI", r, "tm_tercih", "%s Tm %.2f" % (ad, t3))
            hp = primer3.calc_hairpin(s, mv_conc=a.mv, dv_conc=a.dv,
                                      dntp_conc=a.dntp, dna_conc=a.dna).dg
            if hp < a.hairpin_min:
                ekle("KRITIK", r, "hairpin", "%s dG %.0f" % (ad, hp))
            hd = primer3.calc_homodimer(s, mv_conc=a.mv, dv_conc=a.dv,
                                        dntp_conc=a.dntp, dna_conc=a.dna).dg
            if hd < a.selfdimer_min:
                ekle("KRITIK", r, "self_dimer", "%s dG %.0f" % (ad, hd))
            # tabloda yazan degerle karsilastir
            try:
                yazili_tm = float(r["%s_tm" % ad])
                if abs(yazili_tm - t3) > 0.05:
                    ekle("KRITIK", r, "tablo_tm",
                         "%s tabloda %.2f, yeniden olculen %.2f" % (ad, yazili_tm, t3))
            except (KeyError, ValueError):
                ekle("UYARI", r, "tablo_tm", "%s Tm okunamadi" % ad)
            try:
                yazili_gc = float(r["%s_gc" % ad])
                if abs(yazili_gc - g) > 0.6:
                    ekle("KRITIK", r, "tablo_gc",
                         "%s tabloda %.1f, yeniden olculen %.1f" % (ad, yazili_gc, g))
            except (KeyError, ValueError):
                pass

        het = primer3.calc_heterodimer(F, R, mv_conc=a.mv, dv_conc=a.dv,
                                       dntp_conc=a.dntp, dna_conc=a.dna).dg
        if het < a.heterodimer_min:
            ekle("KRITIK", r, "hetero_dimer", "dG %.0f" % het)
        dfark = abs(tm3(F) - tm3(R))
        if dfark >= a.tm_fark_max:
            ekle("KRITIK", r, "tm_farki", "%.2f C" % dfark)
        try:
            if abs(float(r["tm_farki"]) - dfark) > 0.05:
                ekle("KRITIK", r, "tablo_tm_farki",
                     "tabloda %.2f, yeniden olculen %.2f"
                     % (float(r["tm_farki"]), dfark))
        except (KeyError, ValueError):
            pass
        try:
            umin, umax = int(r["urun_min"]), int(r["urun_maks"])
            if not (a.urun_min <= umin and umax <= a.urun_sert_max):
                ekle("KRITIK", r, "urun_boyu", "%d-%d bp" % (umin, umax))
            elif umax > a.urun_max:
                ekle("UYARI", r, "urun_tercih", "%d-%d bp" % (umin, umax))
        except (KeyError, ValueError):
            ekle("UYARI", r, "urun_boyu", "okunamadi")

        # kalip denetimi: F ve rc(R) gercekten kalipta yan yana mi
        if kaliplar:
            bulundu = False
            for et, kal in kaliplar.items():
                if r.get("sinif") and not et.startswith(r["sinif"].split("-")[0]):
                    pass
                i = kal.find(F)
                if i < 0:
                    continue
                j = kal.find(rc(R), i)
                if j < 0:
                    continue
                urun = kal[i:j + len(R)]
                if urun.startswith(F) and urun.endswith(rc(R)) \
                        and a.urun_min <= len(urun) <= a.urun_sert_max:
                    bulundu = True
                    break
            if not bulundu:
                ekle("BILGI", r, "kalipta_tam_eslesme_yok",
                     "hicbir konsensuste birebir F...rc(R) bulunamadi "
                     "(uyumsuzluga izin verilen baglanma ayri denetlenir)")

    # --- alan (domain) tutarliligi -------------------------------------
    # Kural ihlali degil, ama biyolojik olarak tutarsiz sonucu yakalar.
    # Bir hedefin uye kutulari birden fazla ALANDA (A arke, B bakteri,
    # F mantar) bulunuyorsa, azinlik alandaki tasarim o alanin lokus
    # kitapligina dusmus yabanci okumalardan yapilmis demektir. Kural
    # denetiminden gecer ama laboratuvarda hedefi temsil etmez.
    # Alan bilgisi elle yazilmaz, veriden cikarilir: hangi taxid hangi
    # sinif kutularinda geciyorsa o sayilir.
    # Olcu 13 ile AYNI modulden gelir; iki yerde iki ayri kural olursa
    # Excel'den cikarilan cift burada temiz gorunebilir ya da tersi olur.
    if a.hedefler and a.kons:
        _ta = alan_denetimi.taxid_alanlari(a.kons)
        _ht = alan_denetimi.hedef_taxidleri(a.hedefler)
        for r in rows:
            uyumsuz, dag, baskin = alan_denetimi.alan_dagilimi(
                r.get("hedef", ""), r.get("sinif", ""),
                taxid_alan=_ta, hedef_taxid=_ht)
            if uyumsuz:
                ekle("KRITIK", r, "alan_karisimi",
                     "hedefin kutulari %s; bu cift %s alaninda, oysa baskin alan %s"
                     % (", ".join("%s=%d" % kv for kv in sorted(dag.items())),
                        (r.get("sinif") or "")[:1], baskin))

    kritik = [b for b in bulgu if b["agirlik"] == "KRITIK"]
    uyari = [b for b in bulgu if b["agirlik"] == "UYARI"]
    bilgi = [b for b in bulgu if b["agirlik"] == "BILGI"]

    def dok(baslik, liste, sinir=25):
        print("\n%s: %d" % (baslik, len(liste)))
        say = {}
        for b in liste:
            say[b["kural"]] = say.get(b["kural"], 0) + 1
        for k, v in sorted(say.items(), key=lambda x: -x[1]):
            print("   %-24s %d" % (k, v))
        for b in liste[:sinir]:
            print("      %-30s %-4s %-20s %s"
                  % (b["hedef"][:29], b["sinif"], b["kural"], b["ayrinti"][:60]))
        if len(liste) > sinir:
            print("      ... %d kayit daha" % (len(liste) - sinir))

    dok("KRITIK bulgu", kritik)
    dok("UYARI (tercih araligi disi, kural ihlali degil)", uyari)
    if bilgi:
        dok("BILGI", bilgi, sinir=5)

    if a.out and bulgu:
        with open(a.out, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(bulgu[0].keys()),
                               delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(bulgu)
        print("\nyazildi: %s" % a.out)

    print("\n" + "=" * 72)
    if kritik:
        print("SONUC: %d KRITIK bulgu var. Teslim edilmemeli." % len(kritik))
        print("=" * 72)
        sys.exit(1)
    print("SONUC: %d satirin tamami kurallara uyuyor. "
          "Tercih araligi disi %d uyari var, bunlar kural ihlali degildir."
          % (len(rows), len(uyari)))
    print("=" * 72)


if __name__ == "__main__":
    main()
