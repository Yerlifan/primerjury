#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
design_group_primers.py
Produces a primer pair for a target set with more than one member. A genus
specific set, a function group and a domain universal all use the same engine;
the only difference is how large the target set is and how much degeneracy is
allowed.

The method: no alignment is used. Candidates are produced from an anchor
consensus, and every candidate is then scanned against every member and every
competitor with THE BINDING RULE:
    the last two bases have to match the target exactly (extension starts there)
    at most one mismatch in the last five bases
    at most three mismatches across the whole primer
    an overhang on the 5' side is free
    the two primers sit on opposite strands with their 3' ends facing each other

The acceptance criteria:
    a product has to form in EVERY targeted member, within the length range
    no product may form in ANY competitor
    the separation has to be solid: at least one of the primers must find no
    binding site at all in the competitors. A cleanliness that comes from both
    primers binding weakly and only failing together is not accepted.

WHICH DIRECTORY THE CONSENSUS COMES FROM
    Use the canonical directory alone. Older examples pointed at the raw
    consensus directory, which is MIXED IN ORIENTATION (measured: 71 antisense
    against 27 sense). On a reversed consensus an in-silico PCR SILENTLY gives 0
    products, a measured loss of 100 per cent, and the evidence is in
    screening/orientation_impact_test.py.
    The danger is larger in this script, because the input is a GLOB: one member
    stored the other way round is silently counted as "no product in this member"
    and the pair is dropped for nothing. The bin to file mapping is in the index
    inside the canonical directory.

Usage:
  python3 design_group_primers.py \
     --in-group  "canonical_consensus/*_2209.canonical.fa" \
                 "canonical_consensus/*_2223.canonical.fa" \
     --out-group "canonical_consensus/*_394967.canonical.fa" \
     --label acetoclastic_methanogens \
     --out primer_candidates/acetoclastic.tsv
"""
import argparse, csv, glob, importlib.util, os, re, statistics, sys, bisect

HERE = os.path.dirname(os.path.abspath(__file__))


def load_engine():
    'It uses the rule functions of the candidate generator as its one source.'
    p = os.path.join(HERE, "generate_primer_candidates.py")
    if not os.path.exists(p):
        sys.exit(u'generate_primer_candidates.py was not found in the same directory: %s' % HERE)
    spec = importlib.util.spec_from_file_location("engine03", p)
    m = importlib.util.module_from_spec(spec)
    sys.argv_backup, sys.argv = sys.argv, [p]
    try:
        spec.loader.exec_module(m)
    finally:
        sys.argv = sys.argv_backup
    return m


E = load_engine()
rc, gc_pct, composition_ok = E.rc, E.gc_pct, E.composition_ok
read_fasta, read_mask = E.read_fasta, E.read_mask
tm_primer3, tm_biopython = E.tm_primer3, E.tm_biopython
import primer3

IUPAC_SET = {"A": "A", "C": "C", "G": "G", "T": "T",
             "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
             "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def base_match(p, t):
    'Does the primer base p agree with the template base t. The intersection of their IUPAC sets.'
    if t == "N":
        return False          # an N in the template does not count as a match, there is no information
    return bool(set(IUPAC_SET.get(p, "")) & set(IUPAC_SET.get(t, "")))


def build_index(seq, k, azami_acilim=64):
    """k-mer -> baslangic pozisyonlari listesi.

    Kalipta IUPAC kodu bulunabilir (konsensus -A ile uretildiginde) ve
    seed_variants yalnizca A, C, G, T, N uretir. Indeks ham k-mer'i anahtar
    yapsaydi, IUPAC kodu iceren pencereler hicbir cekirdek varyantiyla
    eslesmez ve o baglanma yerleri TAMAMEN gozden kacardi. Olculdu: 12
    gercek konsensus uzerinde find_bindings kuralin dogrudan uygulanmasina
    gore 20 baglanma yerini kaciriyordu ve bunlardan bazilari 'rakipte urun
    yok' kararini tersine cevirecek yerlerdi. Bu yuzden IUPAC iceren her
    k-mer, temsil ettigi somut ACGT k-mer'lerinin hepsine kaydedilir;
    boylece indeks anahtarlari salt ACGT (ve N) olur.

    azami_acilim: bir k-mer'in acilim sayisi bunu asarsa acilim yapilmaz ve
    ham k-mer anahtar olarak birakilir; log'a dusmez ama boyle bir pencere
    zaten bes bazin ucunden fazlasi belirsiz demektir ve primer tasarimina
    girmez."""
    idx = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if all(c in "ACGTN" for c in kmer):
            idx.setdefault(kmer, []).append(i)
            continue
        acilim = [""]
        tasti = False
        for c in kmer:
            secenek = IUPAC_SET.get(c, "ACGT") if c != "N" else "N"
            acilim = [p + o for p in acilim for o in secenek]
            if len(acilim) > azami_acilim:
                tasti = True
                break
        if tasti:
            idx.setdefault(kmer, []).append(i)
            continue
        for v in acilim:
            idx.setdefault(v, []).append(i)
    return idx


def seed_variants(tail, max_mm_in_tail, exact_last):
    """The TEMPLATE k-mer variants allowed for the last 'tail' bases.
        The last 'exact_last' bases are fixed and the rest may carry at most
        max_mm_in_tail mismatches. The template can hold an N; base_match counts an N
        as a mismatch, so at the free positions N is a variant too. N is produced at
        the free positions alone, because the last two bases have to match exactly
        and an exact match with N is not possible."""
    n = len(tail)
    free = n - exact_last
    out = set()

    def rec(i, cur, mm):
        if i == n:
            out.add("".join(cur))
            return
        allowed = IUPAC_SET.get(tail[i], "ACGT")
        fixed = i >= free
        for b in "ACGTN":
            if b == "N":
                if fixed:
                    continue          # son iki bazda N kabul edilemez
                hit = False           # kalipta N her zaman uyumsuzluk
            else:
                hit = b in allowed
                if fixed and not hit:
                    continue
            nm = mm + (0 if hit else 1)
            if nm > max_mm_in_tail:
                continue
            rec(i + 1, cur + [b], nm)

    rec(0, [], 0)
    return out


def find_bindings(oligo, seq, idx, k, a):
    """The oligo's binding sites on seq.
        An overhang on the 5' side is free by rule: the oligo's 5' end may run past
        the template and the part that runs past does not count as a mismatch. The
        3' end has to be inside the template, because extension starts there.
        Returns: [(the 3' end position, the mismatches in the overlapping region)],
        zero based."""
    L, n = len(seq), len(oligo)
    tail = oligo[-k:]
    hits = []
    seen = set()
    for var in seed_variants(tail, a.tail_max_mm, a.exact_last):
        for pos in idx.get(var, ()):          # pos = tail baslangici
            start = pos - (n - k)
            end = start + n
            if end > L:
                continue                      # 3' uc kalip disinda kalamaz
            j0 = max(0, -start)               # 5' sarkma kadar atla
            if n - j0 < a.min_overlap:
                continue
            mm = 0
            ok = True
            for j in range(j0, n):
                if not base_match(oligo[j], seq[start + j]):
                    mm += 1
                    if mm > a.total_max_mm:
                        ok = False
                        break
            if ok:
                key = (end - 1, mm)
                if key not in seen:
                    seen.add(key)
                    hits.append(key)
    return hits


def load_set(patterns, mask_dir=None):
    'Loads a consensus set from the patterns. Returns [(label, sequence)].'
    out, seen = [], set()
    for pat in patterns:
        for p in sorted(glob.glob(pat)):
            if p in seen:
                continue
            seen.add(p)
            name, s = read_fasta(p)
            tag = os.path.basename(p).split("_consensus")[0]
            out.append((tag, s, p))
    return out


def get_args():
    p = argparse.ArgumentParser(description='A primer pair for a target with '
                                            'more than one member')
    p.add_argument("--in-group", nargs="+", required=True,
                   help='the consensus files of the targeted members; a glob '
                        'is accepted')
    p.add_argument("--out-group", nargs="*", default=[],
                   help='the competitor consensus files; a glob is accepted')
    p.add_argument("--anchor", default=None,
                   help="consensus the candidates are generated from; if omitted, the one with fewest Ns")
    p.add_argument("--mask-dir", default=None, help="02 betiginin maske directory")
    p.add_argument("--label", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--max-pairs", type=int, default=2000)
    p.add_argument("--max-oligo", type=int, default=400,
                   help='how many oligos to keep per strand after the '
                        'thermodynamic filter. A single taxon target yields '
                        'thousands of conserved oligos and the forward times '
                        'reverse product reaches millions. The choice is '
                        'stratified by position, so it spreads evenly along '
                        'the template; 0 means unlimited.')
    p.add_argument("--stop-after", type=int, default=20000,
                   help="stop pairing after this many valid pairs; 0 means unlimited")
    # baglanma kurali
    p.add_argument("--exact-last", type=int, default=2)
    p.add_argument("--tail-len", type=int, default=5)
    p.add_argument("--tail-max-mm", type=int, default=1)
    p.add_argument("--total-max-mm", type=int, default=3)
    p.add_argument("--min-overlap", type=int, default=15,
                   help="5' overhang is allowed, but the part overlapping the template must be at least "
                        "bu kadar baz olmali")
    p.add_argument("--competitor-prod-max", type=int, default=0,
                   help='the upper bound used when a product is looked for in '
                        'a competitor. 0 means unlimited, so every band along '
                        'the sequence counts. The rule says no product may '
                        'form in any competitor, so the default is unlimited.')
    # ozgulluk
    p.add_argument("--competitor-prod-min", type=int, default=1,
                   help='the smallest length that counts as a product in a '
                        'competitor; because no band at all is wanted there, '
                        'the default is 1')
    p.add_argument("--orphan-min-mismatch", type=int, default=0,
                   help='0: the strict rule, an orphan primer must not bind '
                        'in the competitors AT ALL. Above 0: a primer whose '
                        'best placement in the competitors carries this many '
                        'mismatches counts as an orphan too, which is the '
                        'relaxed step')
    p.add_argument("--require-orphan-primer", type=int, default=1,
                   help='1 means one of the primers has to find NO binding '
                        'site at all in the competitors')
    # the same oligo, thermodynamics and product rules as generate_primer_candidates.py
    p.add_argument("--len-min", type=int, default=18)
    p.add_argument("--len-max", type=int, default=25)
    p.add_argument("--gc-min", type=float, default=40.0)
    p.add_argument("--gc-max", type=float, default=60.0)
    p.add_argument("--gc-hard-min", type=float, default=35.0)
    p.add_argument("--gc-hard-max", type=float, default=65.0)
    p.add_argument("--gc-clamp-last", type=int, default=5)
    p.add_argument("--gc-clamp-max", type=int, default=3)
    p.add_argument("--homopolymer-max", type=int, default=4)
    p.add_argument("--require-3p-gc", type=int, default=1)
    p.add_argument("--degeneracy-budget", type=int, default=0)
    p.add_argument("--degeneracy-fold-max", type=int, default=4)
    p.add_argument("--keep-variants", action="store_true",
                   help="keep IUPAC sibling variants of the same locus as separate rows")
    p.add_argument("--iupac-max", type=int, default=2)
    p.add_argument("--iupac-clamp-forbidden", type=int, default=5)
    p.add_argument("--tm-min", type=float, default=58.0)
    p.add_argument("--tm-max", type=float, default=62.0)
    p.add_argument("--tm-hard-min", type=float, default=57.0)
    p.add_argument("--tm-hard-max", type=float, default=63.0)
    p.add_argument("--tm-cross-tol", type=float, default=2.0)
    p.add_argument("--hairpin-dg-min", type=float, default=-3000.0)
    p.add_argument("--homodimer-dg-min", type=float, default=-6000.0)
    p.add_argument("--heterodimer-dg-min", type=float, default=-6000.0)
    p.add_argument("--mv", type=float, default=50.0)
    p.add_argument("--dv", type=float, default=1.5)
    p.add_argument("--dntp", type=float, default=0.6)
    p.add_argument("--dna-conc", type=float, default=50.0)
    p.add_argument("--prod-min", type=int, default=70)
    p.add_argument("--prod-max", type=int, default=250)
    p.add_argument("--prod-hard-max", type=int, default=300)
    p.add_argument("--prod-best-min", type=int, default=90)
    p.add_argument("--prod-best-max", type=int, default=150)
    p.add_argument("--pair-tm-diff-max", type=float, default=1.5)
    return p.parse_args()


def main():
    a = get_args()
    ing = load_set(a.in_group)
    outg = load_set(a.out_group) if a.out_group else []
    if not ing:
        sys.exit(u'the target set is empty, check the --in-group patterns')
    # A label collision: because the seqs dictionary is keyed by tag, two members
    # falling on the same label used to collapse silently into one sequence while the
    # output still claimed to cover both.
    seen_tags = {}
    for tag, _, p_ in ing + outg:
        seen_tags.setdefault(tag, []).append(p_)
    dup = {t: v for t, v in seen_tags.items() if len(v) > 1}
    if dup:
        for t, v in dup.items():
            print(u'LABEL COLLISION: %s -> %s' % (t, v), file=sys.stderr)
        sys.exit(u'there are files falling on the same label. Stopped to prevent a silent loss of members; separate the glob patterns.')
    print(u'label             : %s' % a.label)
    print(u'target members    : %d' % len(ing))
    for t, s, p in ing:
        print(u'    target      %-28s %5d bp  N=%d' % (t, len(s), s.count("N")))
    print(u'competitors       : %d' % len(outg))
    for t, s, p in outg:
        print(u'    competitor  %-28s %5d bp  N=%d' % (t, len(s), s.count("N")))
    print(u'binding rule      : last %d bases exact, at most %d mismatches in the last %d, at most %d in total'
          % (a.exact_last, a.tail_len, a.tail_max_mm, a.total_max_mm))
    print(u'degeneracy budget : %d positions, at most %d fold'
          % (a.degeneracy_budget, a.degeneracy_fold_max))

    # --- capa secimi -------------------------------------------------
    if a.anchor:
        anchor = [x for x in ing if a.anchor in x[0] or a.anchor in x[2]]
        if not anchor:
            sys.exit(u'anchor not found: %s' % a.anchor)
        anchor = anchor[0]
    else:
        anchor = min(ing, key=lambda x: (x[1].count("N"), -len(x[1])))
    print(u'anchor consensus  : %s' % anchor[0])

    # --- orientation normalisation ------------------------------------
    # The consensuses were produced without orientation normalisation: consensus2.sh
    # picks a seed read for each taxon and builds the consensus in that read's
    # direction, so some of the members may be stored in reverse. A reversed member
    # makes the forward primer bind the minus strand instead of the plus strand there,
    # and no product forms in that member at all. Every sequence's direction is voted
    # on with conserved probes taken from the anchor, and the reversed ones are turned
    # into their reverse complement.
    K0 = a.tail_len
    probes = [anchor[1][i:i + 20] for i in range(0, len(anchor[1]) - 20, 40)
              if "N" not in anchor[1][i:i + 20]][:40]
    if not probes:
        sys.exit(u'no probe could be produced on the anchor, the consensus holds too many N')

    # A second, independent orientation criterion: SSU motifs conserved across all
    # living things. When the anchor probes bind a distant member not at all
    # (plus=minus=0) the vote cannot decide; that is when these motifs come in. If
    # neither can decide, the member is not passed over silently, it is reported BY
    # NAME.
    UNIV = ["GTGCCAGCMGCCGCGGTAA", "GGATTAGATACCC", "AAACTCAAAGGAATTGACGG",
            "GTGYCAGCMGCCGCGGTAA", "ATTAGATACCCBDGTAGTCC"]

    def _expand(m):
        alt = [""]
        for ch in m:
            opts = IUPAC_SET.get(ch, "ACGT")
            alt = [p + o for p in alt for o in opts]
            if len(alt) > 64:
                return alt[:64]
        return alt

    UNIV_EXP = []
    for m in UNIV:
        UNIV_EXP.extend(_expand(m))

    def univ_vote(seq_):
        r = rc(seq_)
        f = sum(1 for m in UNIV_EXP if m in seq_)
        b = sum(1 for m in UNIV_EXP if m in r)
        return f, b

    def orient(tag, seq_):
        ip = build_index(seq_, K0)
        r = rc(seq_)
        im = build_index(r, K0)
        np_ = sum(1 for o in probes if find_bindings(o, seq_, ip, K0, a))
        nm = sum(1 for o in probes if find_bindings(o, r, im, K0, a))
        if np_ == 0 and nm == 0:
            uf, ub = univ_vote(seq_)
            if uf == 0 and ub == 0:
                return (seq_, np_, nm, False, "KARARSIZ")
            return ((seq_, np_, nm, False, "motif") if uf >= ub
                    else (r, np_, nm, True, "motif"))
        return ((seq_, np_, nm, False, "prob") if np_ >= nm
                else (r, np_, nm, True, "prob"))

    flipped, undecided = [], []
    ing2, outg2 = [], []
    for grp, src, dst in (("hedef", ing, ing2), ("rakip", outg, outg2)):
        for tag, s_, p_ in src:
            s2, np_, nm, fl, how = orient(tag, s_)
            dst.append((tag, s2, p_))
            if fl:
                flipped.append((tag, np_, nm, how))
            if how == "KARARSIZ":
                undecided.append((grp, tag, np_, nm))
    ing, outg = ing2, outg2
    print(u'\norientation normalisation: %d anchor probes plus %d universal motif variants'
          % (len(probes), len(UNIV_EXP)))
    if flipped:
        for t, np_, nm, how in flipped:
            print(u'   found REVERSED, flipped: %-30s plus=%d minus=%d (%s)'
                  % (t, np_, nm, how))
    else:
        print(u'   no sequence disagrees with the anchor orientation')
    if undecided:
        print(u'   WARNING: %d sequences could not be oriented; both directions scored zero:'
              % len(undecided))
        for grp, t, np_, nm in undecided:
            print(u'      %-6s %-30s plus=%d minus=%d' % (grp, t, np_, nm))
        print(u'   That means those sequences are far from the anchor. If they are in the target')
        print(u'   set they will give no product; if they are in the competitor set the')
        print(u'   specificity check becomes unreliable. Review the sets.')

    masked = set()
    if a.mask_dir:
        # The label has the form "A1-1-reads_2209" or "A1_1_reads_2223"; the group and the
        # taxid are used together, so the coordinates of other taxa are not laid over the
        # anchor.
        m_ = re.search(r"reads[-_](\d+)", anchor[0])
        tid = m_.group(1) if m_ else None
        grp = re.split(r"[-_]reads", anchor[0])[0].replace("_", "-")
        pats = []
        if tid:
            pats = [os.path.join(a.mask_dir, "%s_%s_maske.bed" % (grp, tid)),
                    os.path.join(a.mask_dir, "*%s*%s*maske.bed" % (grp, tid))]
        hits = []
        for pat in pats:
            hits.extend(glob.glob(pat))
        hits = sorted(set(hits))
        if not hits:
            sys.exit(u'ERROR: --mask-dir was given but no mask file was found for the anchor \'%s\' (group=%s taxid=%s). Stopped to prevent a silent absence of masking.' % (anchor[0], grp, tid))
        if len(hits) > 1:
            print(u'WARNING: more than one mask file matched the anchor: %s' % hits)
        for cand in hits:
            m, _c, _s = read_mask(cand, len(anchor[1]))
            masked |= m
        print("capada yasak poz. : %d  (maske: %s)" % (len(masked), hits[0]))
    else:
        print(u'forbidden positions in anchor: 0  (no mask given)')

    # --- 1. capadan aday uretimi, kompozisyon suzgeci -----------------
    seq = anchor[1]
    L = len(seq)
    raw, reasons = [], {}
    for start in range(L):
        for ln in range(a.len_min, a.len_max + 1):
            end = start + ln
            if end > L:
                break
            if any(i in masked for i in range(start, end)):
                reasons["maskeli"] = reasons.get("maskeli", 0) + 1
                continue
            win = seq[start:end]
            # A separate check per strand: R's 3' end is the START of the window.
            for strand in ("F", "R"):
                varyant, why = E.iupac_varyantlar(win, a, uc=strand)
                if varyant is None:
                    reasons[why] = reasons.get(why, 0) + 1
                    continue
                for coz, kac in varyant:
                    oligo = coz if strand == "F" else rc(coz)
                    ok, why = composition_ok(oligo, a)
                    if ok:
                        raw.append((strand, start, ln, oligo, kac))
                    else:
                        reasons[why] = reasons.get(why, 0) + 1
    print(u'\noligos after composition filter: %d' % len(raw))
    for k, v in sorted(reasons.items(), key=lambda x: -x[1])[:6]:
        print("   elenen %-16s %d" % (k, v))
    if not raw:
        sys.exit(u'no oligo passed the composition filter')

    # --- 2. iki bagimsiz Tm ve termodinamik ---------------------------
    tm3 = [tm_primer3(o, a) for _, _, _, o, _ in raw]
    tmb = [tm_biopython(o, a) for _, _, _, o, _ in raw]
    offset = statistics.median(x - y for x, y in zip(tm3, tmb))
    print(u'\nprimer3 minus Biopython, median offset: %+.2f C (tolerance %.2f)'
          % (offset, a.tm_cross_tol))
    kept = []
    for (strand, start, ln, o, kac), t3, tb in zip(raw, tm3, tmb):
        if abs((t3 - tb) - offset) > a.tm_cross_tol:
            continue
        if not (a.tm_hard_min <= t3 <= a.tm_hard_max):
            continue
        hp = primer3.calc_hairpin(o, mv_conc=a.mv, dv_conc=a.dv,
                                  dntp_conc=a.dntp, dna_conc=a.dna_conc).dg
        if hp < a.hairpin_dg_min:
            continue
        hd = primer3.calc_homodimer(o, mv_conc=a.mv, dv_conc=a.dv,
                                    dntp_conc=a.dntp, dna_conc=a.dna_conc).dg
        if hd < a.homodimer_dg_min:
            continue
        kept.append(dict(strand=strand, start=start, ln=ln, oligo=o, tm=t3,
                         tmb=tb, hairpin=hp, homodimer=hd, gc=gc_pct(o),
                         iupac=kac))
    print(u'oligos after thermodynamics: %d' % len(kept))
    if not kept:
        sys.exit(u'no oligo passed the thermodynamics filter')

    # --- 3. the binding scan for every member and competitor ----------
    K = a.tail_len
    seqs = {}
    for tag, s, _ in ing + outg:
        seqs[tag] = dict(plus=s, minus=rc(s),
                         idx_plus=build_index(s, K), idx_minus=build_index(rc(s), K))
    ing_tags = [t for t, _, _ in ing]
    out_tags = [t for t, _, _ in outg]

    print(u'\nscanning %d oligos x %d sequences' % (len(kept), len(seqs)))
    bind = {}     # oligo -> {tag: {"plus":[(3'poz,mm)], "minus":[...]}}
    for k in kept:
        o = k["oligo"]
        if o in bind:
            continue
        rec = {}
        for tag, d in seqs.items():
            rec[tag] = dict(
                L=len(d["plus"]),
                plus=find_bindings(o, d["plus"], d["idx_plus"], K, a),
                minus=find_bindings(o, d["minus"], d["idx_minus"], K, a))
        bind[o] = rec

    # her uyede en az bir baglanma yeri olan oligolar
    universal = [k for k in kept
                 if all(bind[k["oligo"]][t]["plus"] or bind[k["oligo"]][t]["minus"]
                        for t in ing_tags)]
    print(u'oligos binding every target member: %d' % len(universal))
    if not universal:
        sys.exit(u'no oligo binds every member. Try raising the degeneracy budget (--degeneracy-budget) or narrowing the target set.')

    # rakiplerde hic baglanmayan oligolar (yetim primer adaylari)
    orphan = set()
    rakip_en_iyi = {}      # oligo -> rakiplerdeki EN IYI (en dusuk) uyumsuzluk
    if out_tags:
        n_tam = 0
        for k in universal:
            o = k["oligo"]
            mm = [m for t in out_tags
                  for _, m in (bind[o][t]["plus"] + bind[o][t]["minus"])]
            en_iyi = min(mm) if mm else None
            rakip_en_iyi[o] = en_iyi
            if en_iyi is None:
                orphan.add(o); n_tam += 1
            elif a.orphan_min_mismatch and en_iyi >= a.orphan_min_mismatch:
                # The step used when the strict rule cannot be met: if even
                # the BEST placement of the primer in the competitors carries
                # this many mismatches, the binding stays weak even when the
                # annealing temperature drops. This step is marked separately
                # in the output.
                orphan.add(o)
        print(u'oligos that bind nowhere in the competitors: %d' % n_tam)
        if a.orphan_min_mismatch:
            print(u'oligos whose best placement in the competitors carries >=%d mismatches: %d (the relaxed step)' % (a.orphan_min_mismatch,
                                            len(orphan) - n_tam))
        dag = {}
        for v in rakip_en_iyi.values():
            dag[v] = dag.get(v, 0) + 1
        print(u'   the best mismatch distribution in the competitors: %s'
              % ", ".join("%s=%d" % ("hic" if k is None else k, v)
                          for k, v in sorted(dag.items(),
                                             key=lambda x: (x[0] is not None, x[0]))))
    else:
        print(u'the competitor set is empty, so the orphan primer rule cannot be applied and NO specificity guarantee is given')
        if a.require_orphan_primer:
            print(u'   (--require-orphan-primer 1 but there is no competitor, so the rule is skipped)')

    # --- 4. ciftleme ve her uyede urun dogrulamasi --------------------
    Fs = [k for k in universal if k["strand"] == "F"]
    Rs = [k for k in universal if k["strand"] == "R"]
    print(u'forward candidates: %d   reverse candidates: %d' % (len(Fs), len(Rs)))

    def ozgulluk_skoru(k):
        """The mismatch count of the oligo's BEST placement in the competitors.
                When it binds nowhere in them it counts as 99. Larger is more specific."""
        v = rakip_en_iyi.get(k["oligo"], None)
        return 99 if v is None else v

    def tabakala(lst, n):
        """Stratified selection by position: the template is cut into n slices and
                one oligo is taken from each. Inside a slice the choice goes first by
                SPECIFICITY and then by how close the Tm is to the middle of the 58 to 62
                band. Going by Tm alone let thousands of non-discriminating oligos in the
                conserved regions push the discriminating ones out; because the
                conserved backbone outnumbers the variable regions so heavily in rDNA,
                this selection criterion decides the result."""
        if not n or len(lst) <= n:
            return lst
        lst = sorted(lst, key=lambda k: k["start"])
        L = lst[-1]["start"] - lst[0]["start"] + 1
        mid = (a.tm_min + a.tm_max) / 2.0
        anahtar = lambda k: (-ozgulluk_skoru(k), abs(k["tm"] - mid))
        kova = {}
        for k in lst:
            b = int((k["start"] - lst[0]["start"]) * n / L)
            cur = kova.get(b)
            if cur is None or anahtar(k) < anahtar(cur):
                kova[b] = k
        out = list(kova.values())
        secili = set(id(k) for k in out)
        # if the bucket count is below n, fill the rest by specificity first and Tm second
        if len(out) < n:
            kalan = [k for k in lst if id(k) not in secili]
            kalan.sort(key=anahtar)
            out.extend(kalan[:n - len(out)])
        return out

    if a.max_oligo:
        f0, r0 = len(Fs), len(Rs)
        Fs, Rs = tabakala(Fs, a.max_oligo), tabakala(Rs, a.max_oligo)
        if len(Fs) < f0 or len(Rs) < r0:
            print(u'layered selection by position (--max-oligo %d): forward %d -> %d, reverse %d -> %d' % (a.max_oligo, f0, len(Fs), r0, len(Rs)))
    pairs = []
    fail_member = {}
    n_noprod = n_comp = n_orph = n_tmd = n_het = 0
    durdu = False
    for f in Fs:
        if durdu:
            break
        for r in Rs:
            if a.stop_after and len(pairs) >= a.stop_after:
                durdu = True
                break
            if abs(f["tm"] - r["tm"]) >= a.pair_tm_diff_max:
                n_tmd += 1
                continue
            # is there a product in every target member
            prods = {}
            ok = True
            for t in ing_tags:
                p = product_len(bind[f["oligo"]][t], bind[r["oligo"]][t],
                                f["ln"], r["ln"], a)
                if p is None:
                    ok = False
                    fail_member[t] = fail_member.get(t, 0) + 1
                    break
                prods[t] = p
            if not ok:
                n_noprod += 1
                continue
            # there must be no product in any competitor
            bad = False
            for t in out_tags:
                # ANY band in a competitor is rejected, not only one inside the target length
                # window; a 370 bp band forms in a PCR too
                cpmax = a.competitor_prod_max if a.competitor_prod_max else 0
                # NO band at all is wanted in a competitor. Applying the lower
                # bound here meant counting cross bands under 70 bp as "no
                # product"; a short band visible on a gel is cross amplification
                # as well.
                if product_len(bind[f["oligo"]][t], bind[r["oligo"]][t],
                               f["ln"], r["ln"], a, pmax=cpmax,
                               pmin=a.competitor_prod_min) is not None:
                    bad = True
                    break
            if bad:
                n_comp += 1
                continue
            # ayrim saglamligi: en az biri rakiplerde hic baglanmamali
            if a.require_orphan_primer and out_tags:
                if f["oligo"] not in orphan and r["oligo"] not in orphan:
                    n_orph += 1
                    continue
            het = primer3.calc_heterodimer(f["oligo"], r["oligo"], mv_conc=a.mv,
                                           dv_conc=a.dv, dntp_conc=a.dntp,
                                           dna_conc=a.dna_conc).dg
            if het < a.heterodimer_dg_min:
                n_het += 1
                continue
            pl = list(prods.values())
            pen = abs(f["tm"] - r["tm"]) * 1.0
            pen += (max(pl) - min(pl)) * 0.02          # uyeler arasi urun tutarliligi
            for t in (f["tm"], r["tm"]):
                if not (a.tm_min <= t <= a.tm_max):
                    pen += min(abs(t - a.tm_min), abs(t - a.tm_max)) * 2.0
            avg = statistics.mean(pl)
            if not (a.prod_best_min <= avg <= a.prod_best_max):
                pen += min(abs(avg - a.prod_best_min), abs(avg - a.prod_best_max)) * 0.05
            # --prod-max yumusak ust sinir: asildiginda ceza, --prod-hard-max
            # ise mutlak sinir. Eski surumde --prod-max hic okunmuyordu.
            if max(pl) > a.prod_max:
                pen += (max(pl) - a.prod_max) * 0.05
            # In each member the mismatch of the BEST binding in that member
            # is taken (the primer will bind there), and then the WORST among
            # the members is reported: that is, "with how many mismatches does
            # it bind in the weakest member".
            mmF = max(min(m for _, m in (bind[f["oligo"]][t]["plus"] +
                                         bind[f["oligo"]][t]["minus"]))
                      for t in ing_tags)
            mmR = max(min(m for _, m in (bind[r["oligo"]][t]["plus"] +
                                         bind[r["oligo"]][t]["minus"]))
                      for t in ing_tags)
            pen += (mmF + mmR) * 0.5
            pairs.append(dict(
                hedef_grubu=a.label, capa=anchor[0],
                ileri_baslangic=f["start"] + 1, ileri_uzunluk=f["ln"],
                geri_baslangic=r["start"] + 1, geri_uzunluk=r["ln"],
                ileri_dizi=f["oligo"], ileri_tm=round(f["tm"], 2), ileri_gc=round(f["gc"], 1),
                ileri_iupac_cozulen=f.get("iupac", 0),
                geri_dizi=r["oligo"], geri_tm=round(r["tm"], 2), geri_gc=round(r["gc"], 1),
                geri_iupac_cozulen=r.get("iupac", 0),
                tm_farki=round(abs(f["tm"] - r["tm"]), 2),
                heterodimer_dg=round(het, 1),
                uye_sayisi=len(ing_tags), rakip_sayisi=len(out_tags),
                urun_min=min(pl), urun_maks=max(pl), urun_ortalama=round(avg, 1),
                en_kotu_uyumsuzluk_ileri=mmF, en_kotu_uyumsuzluk_geri=mmR,
                yetim_primer=("rakip_verilmedi" if not out_tags else
                              ("ileri" if f["oligo"] in orphan else
                               ("geri" if r["oligo"] in orphan else "yok"))),
                yetim_kademe=("rakip_verilmedi" if not out_tags else
                              ("kati" if (rakip_en_iyi.get(f["oligo"]) is None or
                                          rakip_en_iyi.get(r["oligo"]) is None)
                               else ("gevsetilmis" if (f["oligo"] in orphan or
                                                      r["oligo"] in orphan)
                                     else "yok"))),
                ileri_rakip_en_iyi_uyumsuzluk=("hic" if rakip_en_iyi.get(f["oligo"]) is None
                                               else rakip_en_iyi.get(f["oligo"], "")),
                geri_rakip_en_iyi_uyumsuzluk=("hic" if rakip_en_iyi.get(r["oligo"]) is None
                                              else rakip_en_iyi.get(r["oligo"], "")),
                uye_urunleri=";".join("%s=%d" % (t, prods[t]) for t in ing_tags),
                ceza=round(pen, 3)))
    print(u'\npair construction')
    print(u'   dropped, pair Tm difference        : %d' % n_tmd)
    print(u'   dropped, no product in one member  : %d' % n_noprod)
    if fail_member:
        worst = sorted(fail_member.items(), key=lambda x: -x[1])[:5]
        print("      the members blocking most: %s"
              % ", ".join("%s=%d" % kv for kv in worst))
    print(u'   dropped, product forms in a competitor : %d' % n_comp)
    print(u'   dropped, no orphan primer          : %d' % n_orph)
    print(u'   dropped, hetero-dimer dG           : %d' % n_het)
    print(u'valid pairs                          : %d%s'
          % (len(pairs), u'  (stopped early via --stop-after)' if durdu else ""))
    if not pairs:
        sys.exit(u'no valid pair was found')
    pairs.sort(key=lambda x: x["ceza"])
    # The IUPAC sibling variants of the same locus are reduced to a single row.
    # Otherwise the first ten candidates fill up with the allele variants of one
    # region and the table loses its diversity. What is kept is the variant with the
    # lowest penalty; how many siblings there were is reported in its own column.
    if not a.keep_variants:
        onceki = len(pairs)
        secili, gorulen = [], {}
        for pr in pairs:
            anahtar = (pr["ileri_baslangic"], pr["ileri_uzunluk"],
                       pr["geri_baslangic"], pr["geri_uzunluk"])
            if anahtar in gorulen:
                gorulen[anahtar]["lokus_varyant_sayisi"] += 1
                continue
            pr["lokus_varyant_sayisi"] = 1
            gorulen[anahtar] = pr
            secili.append(pr)
        pairs = secili
        print(u'   locus merging: %d pairs -> %d loci' % (onceki, len(pairs)))
    else:
        for pr in pairs:
            pr["lokus_varyant_sayisi"] = 1
    pairs = pairs[:a.max_pairs]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pairs[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(pairs)
    print(u'\nwritten: %s' % a.out)
    print(u'\nFive best candidates:')
    for p in pairs[:5]:
        print(u'  penalty=%.2f  product %d-%d bp  orphan=%s  F=%s (Tm %.1f)  R=%s (Tm %.1f)'
              % (p["ceza"], p["urun_min"], p["urun_maks"], p["yetim_primer"],
                 p["ileri_dizi"], p["ileri_tm"], p["geri_dizi"], p["geri_tm"]))


def _one_config(b_plus, b_minus, ln_plus, ln_minus, L, pmin, pmax):
    """One primer on the plus strand and the other on the minus strand. The
        coordinate conversion: position m on the minus strand corresponds to L-1-m on
        the plus strand."""
    best = None
    for fend, _ in b_plus:
        fstart = fend - ln_plus + 1
        for rend, _ in b_minus:
            r_left_plus = L - 1 - rend
            if fend >= r_left_plus:             # 3' uclari birbirine bakmiyor
                continue
            r_right_plus = r_left_plus + ln_minus - 1
            plen = r_right_plus - fstart + 1
            if plen < pmin:
                continue
            if pmax and plen > pmax:
                continue
            if best is None or plen < best:
                best = plen
    return best


def product_len(bf, br, lnf, lnr, a, pmax=None, pmin=None):
    """Because the template is double stranded, BOTH configurations give a real
        product:
        (1) the first primer on the plus strand and the second on the minus strand
        (2) the first primer on the minus strand and the second on the plus strand
        Looking at only one of them marks a sequence stored the other way round as
        "no product", and in a competitor that escapes the specificity check
        entirely.
        Returns: the shortest valid product length, or None."""
    L = bf.get("L") or br.get("L")
    if not L:
        return None
    if pmax is None:
        pmax = a.prod_hard_max
    if pmin is None:
        pmin = a.prod_min
    cands = [
        _one_config(bf["plus"], br["minus"], lnf, lnr, L, pmin, pmax),
        _one_config(br["plus"], bf["minus"], lnr, lnf, L, pmin, pmax),
    ]
    cands = [c for c in cands if c is not None]
    return min(cands) if cands else None


if __name__ == "__main__":
    main()
