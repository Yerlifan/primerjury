#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_deliverables.py
AN INDEPENDENT AUDIT OF THE DELIVERED TABLE.

This script imports not one function of the design code. It reads final_primers.tsv
and measures every row again FROM SCRATCH against the rules of the meeting decision.
The aim is to show that no silent drift is left between the design code and the
delivered table: the Tm, GC, product length and tm_farki values written in the table
are recomputed and compared as well.

The rules audited
  Oligo   : A/C/G/T only, length 18-25, GC 40-60 percent (hard 35-65),
            G or C at the 3' end, at most 3 G or C in the last five bases,
            at most 4 identical bases in a row
  Thermo  : Tm 58-62 (hard 57-63) with two libraries, hairpin >= -3000,
            self-dimer >= -6000, hetero-dimer >= -6000
  Pair    : |TmF - TmR| < 1,5 ; the product 70-250 (hard 300)
  Table   : are the ileri_tm/geri_tm/ileri_gc/geri_gc/tm_farki values written the
            same as the recomputed ones
  Template: (if --consensus is given) is the reverse primer really the reverse complement
            of the template, and does the product start with the forward primer and
            end with the reverse complement of the reverse primer

Usage:
  python3 check_deliverables.py --final pr_final --consensus pr_kons/konsensus

"""
import argparse, csv, glob, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import field_audit

try:
    import primer3
except ImportError:
    sys.exit(u'primer3-py is required: pip install primer3-py --break-system-packages')
try:
    from Bio.SeqUtils import MeltingTemp as mt
    from Bio.Seq import Seq
except ImportError:
    sys.exit(u'biopython is required: pip install biopython --break-system-packages')

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
    p.add_argument("--consensus", default=None,
                   help="baskin alel consensus directory; if given the template is checked as well")
    p.add_argument("--targets", default=None,
                   help="targets.tsv; if given field consistency is checked as well")
    p.add_argument("--passed-only", type=int, default=1)
    p.add_argument("--out", default=None, help="findings TSV'si")
    # toplanti karari esikleri
    p.add_argument("--len-min", type=int, default=18)
    p.add_argument("--len-max", type=int, default=25)
    p.add_argument("--gc-min", type=float, default=40.0)
    p.add_argument("--gc-max", type=float, default=60.0)
    p.add_argument("--gc-hard-min", type=float, default=35.0)
    p.add_argument("--gc-hard-max", type=float, default=65.0)
    p.add_argument("--gc-clamp-last", type=int, default=5)
    p.add_argument("--gc-clamp-max", type=int, default=3)
    p.add_argument("--homopolymer-max", type=int, default=4)
    p.add_argument("--tm-min", type=float, default=58.0)
    p.add_argument("--tm-max", type=float, default=62.0)
    p.add_argument("--tm-hard-min", type=float, default=57.0)
    p.add_argument("--tm-hard-max", type=float, default=63.0)
    p.add_argument("--pair-tm-diff-max", type=float, default=1.5)
    p.add_argument("--tm-cross-tol", type=float, default=2.0)
    p.add_argument("--hairpin-min", type=float, default=-3000.0)
    p.add_argument("--selfdimer-min", type=float, default=-6000.0)
    p.add_argument("--heterodimer-min", type=float, default=-6000.0)
    p.add_argument("--prod-min", type=int, default=70)
    p.add_argument("--prod-max", type=int, default=250)
    p.add_argument("--prod-hard-max", type=int, default=300)
    # the thermodynamic conditions, the same as in generate_primer_candidates.py and design_group_primers.py
    p.add_argument("--mv", type=float, default=50.0)
    p.add_argument("--dv", type=float, default=1.5)
    p.add_argument("--dntp", type=float, default=0.6)
    # It has to be the same as the --dna-conc default of generate_primer_candidates.py
    # and design_group_primers.py; if it differs, the remeasured Tm drifts
    # systematically and produces a false finding.
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
    tsv = os.path.join(a.final, "final_primers.tsv")
    if not os.path.exists(tsv):
        sys.exit(u'not found: %s' % tsv)
    rows = list(csv.DictReader(open(tsv, encoding="utf-8"), delimiter="\t"))
    if a.passed_only:
        rows = [r for r in rows if r.get("ozgulluk_durum") == "GECTI"]
    if not rows:
        sys.exit(u'there is no row to audit')

    kaliplar = kalip_yukle(a.consensus)
    print("=" * 72)
    print(u'DELIVERY AUDIT  (an independent re-measurement)')
    print("=" * 72)
    print(u'rows audited : %d' % len(rows))
    print(u'template files   : %d' % len(kaliplar))
    print("termodinamik     : Na=%.0f Mg=%.1f dNTP=%.1f DNA=%.0f"
          % (a.mv, a.dv, a.dntp, a.dna))
    print()

    def tm3(s):
        return primer3.calc_tm(s, mv_conc=a.mv, dv_conc=a.dv,
                               dntp_conc=a.dntp, dna_conc=a.dna)

    def tmb(s):
        return float(mt.Tm_NN(Seq(s), nn_table=mt.DNA_NN3, Na=a.mv, Mg=a.dv,
                              dNTPs=a.dntp, dnac1=a.dna, dnac2=0, saltcorr=7))

    bulgu = []          # (weight, target, class, rule, detail)

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
    print(u'the median offset between the two libraries: %+.2f C  (tolerance %.2f C)'
          % (kayma, a.tm_cross_tol))

    for r in rows:
        F, R = r["ileri_dizi"].upper(), r["geri_dizi"].upper()
        for ad, s in (("ileri", F), ("geri", R)):
            if set(s) - set("ACGT"):
                ekle("KRITIK", r, "alfabe",
                     "%s primerde ACGT disi baz: %s" % (ad, sorted(set(s) - set("ACGT"))))
                continue
            if not (a.len_min <= len(s) <= a.len_max):
                ekle("KRITIK", r, "uzunluk", "%s uzunluk %d" % (ad, len(s)))
            g = gc_yuzde(s)
            if not (a.gc_hard_min <= g <= a.gc_hard_max):
                ekle("KRITIK", r, "gc_sert", "%s GC %%%.1f" % (ad, g))
            elif not (a.gc_min <= g <= a.gc_max):
                ekle("UYARI", r, "gc_tercih", "%s GC %%%.1f" % (ad, g))
            if s[-1] not in "GC":
                ekle("KRITIK", r, "3p_gc_kilit", "%s 3' uc %s" % (ad, s[-1]))
            kuyruk = s[-a.gc_clamp_last:]
            if kuyruk.count("G") + kuyruk.count("C") > a.gc_clamp_max:
                ekle("KRITIK", r, "3p_asiri_sabit",
                     "%s son %d bazda %d G/C" % (ad, a.gc_clamp_last,
                                                 kuyruk.count("G") + kuyruk.count("C")))
            if en_uzun_tekrar(s) > a.homopolymer_max:
                ekle("KRITIK", r, "homopolimer",
                     u'%s %d identical bases in a row' % (ad, en_uzun_tekrar(s)))
            t3, tb = tm3(s), tmb(s)
            if abs((t3 - tb) - kayma) > a.tm_cross_tol:
                ekle("KRITIK", r, "tm_capraz",
                     "%s primer3 %.2f, biopython %.2f, kaymadan sapma %.2f"
                     % (ad, t3, tb, abs((t3 - tb) - kayma)))
            if not (a.tm_hard_min <= t3 <= a.tm_hard_max):
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
                         u'%s the table says %.2f, remeasured %.2f' % (ad, yazili_tm, t3))
            except (KeyError, ValueError):
                ekle("UYARI", r, "tablo_tm", u'%s the Tm could not be read' % ad)
            try:
                yazili_gc = float(r["%s_gc" % ad])
                if abs(yazili_gc - g) > 0.6:
                    ekle("KRITIK", r, "tablo_gc",
                         u'%s the table says %.1f, remeasured %.1f' % (ad, yazili_gc, g))
            except (KeyError, ValueError):
                pass

        het = primer3.calc_heterodimer(F, R, mv_conc=a.mv, dv_conc=a.dv,
                                       dntp_conc=a.dntp, dna_conc=a.dna).dg
        if het < a.heterodimer_min:
            ekle("KRITIK", r, "hetero_dimer", "dG %.0f" % het)
        dfark = abs(tm3(F) - tm3(R))
        if dfark >= a.pair_tm_diff_max:
            ekle("KRITIK", r, "tm_farki", "%.2f C" % dfark)
        try:
            if abs(float(r["tm_farki"]) - dfark) > 0.05:
                ekle("KRITIK", r, "tablo_tm_farki",
                     u'the table says %.2f, remeasured %.2f'
                     % (float(r["tm_farki"]), dfark))
        except (KeyError, ValueError):
            pass
        try:
            umin, umax = int(r["urun_min"]), int(r["urun_maks"])
            if not (a.prod_min <= umin and umax <= a.prod_hard_max):
                ekle("KRITIK", r, "urun_boyu", "%d-%d bp" % (umin, umax))
            elif umax > a.prod_max:
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
                        and a.prod_min <= len(urun) <= a.prod_hard_max:
                    bulundu = True
                    break
            if not bulundu:
                ekle("BILGI", r, "kalipta_tam_eslesme_yok",
                     u'no exact F...rc(R) was found in any consensus (binding that allows a mismatch is audited separately)')

    # --- domain consistency --------------------------------------------
    # Not a rule violation, but it catches a biologically inconsistent
    # result. If a target's member bins sit in more than one DOMAIN (A
    # archaea, B bacteria, F fungi), the design in the minority domain was
    # made from foreign reads that fell into that domain's locus library. It
    # passes the rule check but it does not represent the target in the
    # laboratory.
    # The domain is not written by hand, it is derived from the data:
    # whichever class bins a taxid appears in is what counts.
    # The measure comes from the SAME module as check 13; with two separate
    # rules in two places, a pair taken out of the Excel could look clean here
    # or the other way round.
    if a.targets and a.consensus:
        _ta = field_audit.taxid_alanlari(a.consensus)
        _ht = field_audit.hedef_taxidleri(a.targets)
        for r in rows:
            uyumsuz, dag, baskin = field_audit.alan_dagilimi(
                r.get("hedef", ""), r.get("sinif", ""),
                taxid_alan=_ta, hedef_taxid=_ht)
            if uyumsuz:
                ekle("KRITIK", r, "alan_karisimi",
                     u'the target\'s bins are %s; this pair is in the %s domain, while the dominant domain is %s'
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
            print(u'      ... %d more records' % (len(liste) - sinir))

    dok("KRITIK bulgu", kritik)
    dok('A WARNING: outside the preferred range, but no rule is broken', uyari)
    if bilgi:
        dok("BILGI", bilgi, sinir=5)

    if a.out and bulgu:
        with open(a.out, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(bulgu[0].keys()),
                               delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(bulgu)
        print(u'\nwritten: %s' % a.out)

    print("\n" + "=" * 72)
    if kritik:
        print(u'RESULT: there are %d CRITICAL findings. This should not be delivered.' % len(kritik))
        print("=" * 72)
        sys.exit(1)
    print(u'RESULT: all %d rows follow the rules. There are %d warnings outside the preferred range, and those are not rule violations.'
          % (len(rows), len(uyari)))
    print("=" * 72)


if __name__ == "__main__":
    main()
