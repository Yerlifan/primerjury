#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
split_clusters.py
Tek bir universal primer ile kapsanamayan bir hedef kumesini, veriden
turetilerek en az sayida alt kumeye boler.

Yontem, acgozlu birakma: kume 04 ile denenir. Gecerli cift cikmazsa 04'un
"en cok engelleyen uyeler" raporundan en cok engelleyen uye cikarilir ve
tekrar denenir. Ilk gecerli cift ciktiginda kalan uyeler SET-1 olur,
cikarilanlar yeni bir kume olarak ayni islemden gecer.

Boylece bolme, dizi benzerligine degil TASARLANABILIRLIGE gore yapilir.
ITS gibi cok ayrisan bolgelerde k-mer benzerligi anlamsiz derecede dusuk
kalir (F1 grubunda ikili Jaccard 0,028 ile 0,035 arasi), bu yuzden benzerlik
kumelemesi yanlis yol olurdu.

Kullanim:
  python3 split_clusters.py --in-group "consensus sequences/F1-*/*.fasta" \
      --out-group "consensus sequences/B-1/*.fasta" \
      --label Mantar_F1 --outdir primer_adaylari/F1_setleri
"""
import argparse, glob, os, re, subprocess, sys, json

HERE = os.path.dirname(os.path.abspath(__file__))
ENGINE = os.path.join(HERE, "design_group_primers.py")


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--in-group", nargs="+", required=True)
    p.add_argument("--out-group", nargs="*", default=[])
    p.add_argument("--label", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--max-sets", type=int, default=6)
    p.add_argument("--min-members", type=int, default=2,
                   help="leftover clusters with fewer members than this are designed one by one")
    p.add_argument("--extra", nargs=argparse.REMAINDER, default=[],
                   help="04'e aynen gecirilecek ek bayraklar")
    return p.parse_args()


def expand(pats):
    out, seen = [], set()
    for pat in pats:
        for f in sorted(glob.glob(pat)):
            if f not in seen:
                seen.add(f)
                out.append(f)
    return out


def tag_of(p):
    """04 ile AYNI etiketleme. Iki yerde farkli olursa engelleyen uye
    eslestirilemez ve bolme sessizce durur."""
    return os.path.basename(p).split("_consensus")[0]


def run_engine(members, outg, label, out_tsv, extra):
    cmd = [sys.executable, ENGINE, "--in-group"] + members
    if outg:
        cmd += ["--out-group"] + outg
    cmd += ["--label", label, "--out", out_tsv] + list(extra)
    r = subprocess.run(cmd, capture_output=True, text=True)
    txt = r.stdout + r.stderr
    n = 0
    m = re.search(r"gecerli cift sayisi\s*:\s*(\d+)", txt)
    if m:
        n = int(m.group(1))
    blockers = []
    m = re.search(r"en cok engelleyen uyeler:\s*(.+)", txt)
    if m:
        for part in m.group(1).split(","):
            part = part.strip()
            if "=" in part:
                t, c = part.rsplit("=", 1)
                try:
                    blockers.append((t.strip(), int(c)))
                except ValueError:
                    pass
    best = ""
    mb = re.search(r"En iyi bes aday:\n(.+)", txt, re.S)
    if mb:
        best = mb.group(1).strip().splitlines()[0].strip() if mb.group(1).strip() else ""
    return n, blockers, best, txt


def main():
    a = get_args()
    os.makedirs(a.outdir, exist_ok=True)
    pool = expand(a.in_group)
    outg = expand(a.out_group)
    if not pool:
        sys.exit("hedef kumesi bos")
    print("baslangic uye sayisi : %d" % len(pool))
    print(u'number of competitors : %d' % len(outg))
    print("motor                : %s" % ENGINE)

    sets, si, log = [], 0, []
    while pool and si < a.max_sets:
        si += 1
        cur = list(pool)
        removed = []
        while True:
            if len(cur) < a.min_members:
                break
            tsv = os.path.join(a.outdir, "%s_SET%d.tsv" % (a.label, si))
            n, blockers, best, txt = run_engine(cur, outg,
                                                "%s_SET%d" % (a.label, si), tsv, a.extra)
            print(u'  SET%d attempt: %d members -> %d valid pairs' % (si, len(cur), n))
            if n > 0:
                sets.append(dict(set=si, members=[tag_of(x) for x in cur],
                                 pairs=n, best=best, tsv=tsv))
                print(u'  SET%d done: %d members, %d pairs' % (si, len(cur), n))
                if best:
                    print("     en iyi: %s" % best)
                pool = removed
                break
            if not blockers:
                print(u'  SET%d: no pair, and no blocking member was reported, so it stopped' % si)
                cur = []
                break
            worst = blockers[0][0]
            hit = [x for x in cur if tag_of(x) == worst]
            if not hit:
                print("  engelleyen uye eslestirilemedi: %s" % worst)
                cur = []
                break
            cur.remove(hit[0])
            removed.append(hit[0])
            log.append((si, worst, blockers[0][1]))
            print("     cikarildi: %s (%d elemeden sorumlu)" % (worst, blockers[0][1]))
        else:
            pass
        if not sets or sets[-1]["set"] != si:
            # Bu turda hic cift cikmadi. `cur` icinde hala uye olabilir
            # (min_members'in altina dusuldugu, engelleyen raporlanmadigi ya
            # da eslestirilemedigi durumlar). Eski surumde bu uyeler ne bir
            # SET'e ne de 'kalan' listesine giriyordu, yani sessizce
            # kayboluyorlardi ve rapor grubun tamamen karsilandigini
            # gosteriyordu. Kalanlar havuza geri verilir.
            pool = removed + [x for x in cur if x not in removed]
            if not pool:
                break
            if len(pool) == len(expand(a.in_group)):
                # hic ilerleme yok, sonsuz donguye girmemek icin durulur
                break

    print("\n=== SONUC ===")
    for s in sets:
        print(u'SET%d : %d members, %d valid pairs' % (s["set"], len(s["members"]), s["pairs"]))
        print("        uyeler: %s" % ", ".join(s["members"]))
        if s["best"]:
            print("        %s" % s["best"])
    kalan = [tag_of(x) for x in pool]
    # Koruma denetimi: girdi kumesindeki her uye ya bir SET'te ya kalan
    # listesinde olmali. Tutmuyorsa sessiz kayip var demektir.
    girdi = set(tag_of(x) for x in expand(a.in_group))
    yazilan = set(kalan)
    for st in sets:
        yazilan |= set(st["members"])
    kayip = sorted(girdi - yazilan)
    if kayip:
        print(u'\nWARNING: an unreported member exists and is being added to the list: %s'
              % ", ".join(kayip))
        kalan = sorted(set(kalan) | set(kayip))
    if kalan:
        print("Kalan, tek basina tasarlanmasi gereken uye: %s" % ", ".join(kalan))
    with open(os.path.join(a.outdir, "%s_bolme.json" % a.label), "w",
              encoding="utf-8") as fh:
        json.dump(dict(sets=sets, kalan=kalan, cikarma_kaydi=log), fh,
                  ensure_ascii=False, indent=1)
    print("\nyazildi: %s" % os.path.join(a.outdir, "%s_bolme.json" % a.label))


if __name__ == "__main__":
    main()
