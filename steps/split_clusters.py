#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
split_clusters.py
Splits a target set that cannot be covered by a single universal primer into the
smallest number of subsets, derived from the data.

The method is greedy dropping: the set is tried with design_group_primers.py. If no
valid pair comes out, the member blocking most is taken from that script's "the
members blocking most" report and it is tried again. When the first valid pair comes
out, the members left become SET-1 and the ones dropped go through the same process
as a new set.

So the split rests not on sequence similarity but on THE DESIGN ITSELF: a set stays
together as long as one pair can cover it.

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
    """THE SAME labelling as design_group_primers.py. If the two differ, the blocking
    member cannot be matched and the split stops silently.

    """
    return os.path.basename(p).split("_consensus")[0]


def run_engine(members, outg, label, out_tsv, extra):
    cmd = [sys.executable, ENGINE, "--in-group"] + members
    if outg:
        cmd += ["--out-group"] + outg
    cmd += ["--label", label, "--out", out_tsv] + list(extra)
    r = subprocess.run(cmd, capture_output=True, text=True)
    txt = r.stdout + r.stderr
    n = 0
    m = re.search(r"valid pairs\s*:\s*(\d+)", txt)
    if m:
        n = int(m.group(1))
    blockers = []
    m = re.search(r"the members blocking most:\s*(.+)", txt)
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
    mb = re.search(r"Five best candidates:\n(.+)", txt, re.S)
    if mb:
        best = mb.group(1).strip().splitlines()[0].strip() if mb.group(1).strip() else ""
    return n, blockers, best, txt


def main():
    a = get_args()
    os.makedirs(a.outdir, exist_ok=True)
    pool = expand(a.in_group)
    outg = expand(a.out_group)
    if not pool:
        sys.exit(u'the target set is empty')
    print(u'starting member count : %d' % len(pool))
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
                print("  the blocking member could not be matched: %s" % worst)
                cur = []
                break
            cur.remove(hit[0])
            removed.append(hit[0])
            log.append((si, worst, blockers[0][1]))
            print("     cikarildi: %s (%d elemeden sorumlu)" % (worst, blockers[0][1]))
        else:
            pass
        if not sets or sets[-1]["set"] != si:
            # No pair came out in this round. `cur` may still hold members (the
            # cases where it fell below min_members, where no blocker was reported,
            # or where one could not be matched). In the old version those members
            # went into neither a SET nor the 'kalan' list, so they disappeared
            # silently and the report showed the group as fully met. The ones left
            # are given back to the pool.
            pool = removed + [x for x in cur if x not in removed]
            if not pool:
                break
            if len(pool) == len(expand(a.in_group)):
                # no progress at all; it stops so as not to enter an infinite loop
                break

    print(u'\n=== RESULT ===')
    for s in sets:
        print(u'SET%d : %d members, %d valid pairs' % (s["set"], len(s["members"]), s["pairs"]))
        print("        uyeler: %s" % ", ".join(s["members"]))
        if s["best"]:
            print("        %s" % s["best"])
    kalan = [tag_of(x) for x in pool]
    # The conservation check: every member of the input set has to be either in a SET
    # or in the 'kalan' list. If that does not hold there is a silent loss.
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
        print(u'The members left, which have to be designed on their own: %s' % ", ".join(kalan))
    with open(os.path.join(a.outdir, "%s_bolme.json" % a.label), "w",
              encoding="utf-8") as fh:
        json.dump(dict(sets=sets, kalan=kalan, cikarma_kaydi=log), fh,
                  ensure_ascii=False, indent=1)
    print("\nwritten: %s" % os.path.join(a.outdir, "%s_bolme.json" % a.label))


if __name__ == "__main__":
    main()
