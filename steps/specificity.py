#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
specificity.py
Puts the candidate pairs produced by batch_design.py through the specificity and
verification rules of the panel decision.

The rules applied:
  1. VERIFICATION ON RAW READS. A pair that passed on the consensus is tested
     again on the target's raw reads. If no product comes out, the candidate is
     flagged "numuneden_dogrulanamadi" and eliminated.
  2. THE COMPETITOR RATIO, THE WILSON LOWER BOUND. The proportion of competitor
     taxa reads giving a product is judged by its Wilson lower bound rather than
     the raw ratio, so a small sample cannot make a competitor look clean.

"""
import hashlib
import argparse, csv, datetime, glob, importlib.util, json, math, os, re
import subprocess, sys, time, collections

HERE = os.path.dirname(os.path.abspath(__file__))
TS = lambda: datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
_LOG = [None]


def log(msg):
    line = "[%s] %s" % (TS(), msg)
    print(line, flush=True)
    if _LOG[0]:
        with open(_LOG[0], "a", encoding="utf-8") as fh:
            fh.write(line + "\n")


def load(p, n):
    sp = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(sp)
    bak, sys.argv = sys.argv, [p]
    try:
        sp.loader.exec_module(m)
    finally:
        sys.argv = bak
    return m


E = load(os.path.join(HERE, "generate_primer_candidates.py"), "e03")
G = load(os.path.join(HERE, "design_group_primers.py"), "g04")
rc = E.rc


class Kural:
    """Exactly the same binding rule as steps/design_group_primers.py."""
    exact_last = 2
    tail_len = 5
    tail_max_mm = 1
    total_max_mm = 3
    min_overlap = 15


KURAL = Kural()


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import alignment

MAPPY = alignment.ARKA_UC is not None


def _girdi_parmak_izi(a):
    """A summary of every input the run depends on. When it changes the checkpoint
        drops. Reusing the old verification results silently when the candidate
        files had changed made the run look finished while the result belonged to
        the old input."""
    h = hashlib.sha256()
    h.update(("adaylar=%s\n" % os.path.abspath(a.candidates)).encode())
    for f in sorted(glob.glob(os.path.join(a.candidates, "*__*.tsv"))):
        try:
            st = os.stat(f)
            h.update(("%s|%d|%d\n" % (os.path.basename(f), st.st_size,
                                      int(st.st_mtime))).encode())
        except OSError:
            pass
    for ek in ("ayirt_edilemez.tsv",):
        yol = os.path.join(a.candidates, ek)
        try:
            st = os.stat(yol)
            h.update(("%s|%d|%d\n" % (ek, st.st_size, int(st.st_mtime))).encode())
        except OSError:
            h.update(("%s|yok\n" % ek).encode())
    try:
        st = os.stat(a.targets)
        h.update(("hedefler|%d|%d\n" % (st.st_size, int(st.st_mtime))).encode())
    except OSError:
        pass
    # The engine scripts go into the fingerprint too: when the binding rule or the scan
    # logic changes, the old verification results are no longer valid.
    _burada = os.path.dirname(os.path.abspath(__file__))
    for _b in ("design_group_primers.py", "generate_primer_candidates.py",
               "specificity.py", "alignment.py"):
        try:
            _st = os.stat(os.path.join(_burada, _b))
            h.update(("%s|%d|%d\n" % (_b, _st.st_size, int(_st.st_mtime))).encode())
        except OSError:
            h.update(("%s|yok\n" % _b).encode())
    h.update(("kons=%s|top=%d|maxokuma=%d|wilson=%.5f|susici=%.5f|"
              "minuye=%s|bulasma=%d|mfe=%d|blast=%d\n"
              % (os.path.abspath(a.consensus) if a.consensus else "-", a.top,
                 a.max_reads, a.competitor_wilson_max, a.within_strain_diff_max,
                 getattr(a, "min_uye_orani", ""), a.contamination_sample,
                 int(a.atla_mfe), int(a.atla_blast))).encode())
    h.update(("bulasma_min=%s|sizinti_tavan=%s\n"
              % (getattr(a, "bulasma_min_okuma", ""),
                 getattr(a, "sizinti_tavan", ""))).encode())
    return h.hexdigest()[:16]


def wilson_ust(k, n, z=1.96):
    """Oranin Wilson UST siniri."""
    if n == 0:
        return 1.0
    p = k / n
    d = 1 + z * z / n
    m = p + z * z / (2 * n)
    s_ = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return min(1.0, (m + s_) / d)


def wilson_alt(k, n, z=1.96):
    """Judges a ratio by its Wilson lower bound rather than by the raw count.
        k successes, n trials. Returns 0 when n=0.

    """
    if n <= 0:
        return 0.0
    p = k / n
    d = 1 + z * z / n
    merkez = p + z * z / (2 * n)
    yari = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return max(0.0, (merkez - yari) / d)


def fastq_oku(path, limit):
    """Yields the read sequences, up to the limit."""
    n = 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                yield line.strip().upper()
                n += 1
                if limit and n >= limit:
                    return


class _PArg:
    """The smallest set of arguments product_len needs."""
    prod_min = 50
    prod_hard_max = 400


_OKUMA_BELLEK = {}


def okumalari_al(path, limit):
    """The same fastq is scanned again and again for every candidate of a target;
        it is read from disk once and held in memory.

    """
    k = (path, limit)
    v = _OKUMA_BELLEK.get(k)
    if v is None:
        v = [(s, rc(s)) for s in fastq_oku(path, limit)]
        # The cache holds at most --onbellek-dosya fastq files. In the old version it was
        # emptied completely on every miss, so the hit rate was effectively zero and the
        # same fastq was read again for every candidate.
        while len(_OKUMA_BELLEK) >= 3:
            _OKUMA_BELLEK.pop(next(iter(_OKUMA_BELLEK)))          # tek hedefin dosyalari; bellek sismesin
        _OKUMA_BELLEK[k] = v
    return v


_BULASMA = {}


def capraz_bulasma(hedef_kons_yolu, rakip_kons_yolu, rakip_fq, ornek=400,
                   min_uzunluk=400):
    """How much of the reads in a competitor bin actually belong to THE TARGET.

        Kraken bins can leak into one another: some of the reads in a bin are the
        molecules of another taxon. Those reads naturally give a product with the target
        primer. So a 'there is a product in the competitor' decision cannot be made
        without comparing against this measured leakage ratio. Returns: (k, n) = the
        number of reads that fit the target better, and the number of reads examined.

    """
    if not MAPPY:
        return (0, 0)
    anahtar = (hedef_kons_yolu, rakip_kons_yolu, rakip_fq)
    if anahtar in _BULASMA:
        return _BULASMA[anahtar]
    def _oku(p):
        """Ns at the start and the end are trimmed; INNER Ns are left in place. Deleting
                the inner Ns joins the two sides of a coverage gap and produces a chimeric
                reference; the alignment results become meaningless at that junction.

        """
        d = "".join(l.strip() for l in open(p, encoding="utf-8",
                                            errors="replace")
                    if not l.startswith(">")).upper()
        return d.strip("N")
    try:
        A_h = alignment.Hizalayici(seq=_oku(hedef_kons_yolu), preset="map-ont")
        A_r = alignment.Hizalayici(seq=_oku(rakip_kons_yolu), preset="map-ont")
    except Exception:
        _BULASMA[anahtar] = (0, 0)
        return (0, 0)
    if not A_h or not A_r:
        _BULASMA[anahtar] = (0, 0)
        return (0, 0)
    okumalar = {}
    n = 0
    try:
        with open(rakip_fq, errors="replace") as fh:
            for i, line in enumerate(fh):
                if i % 4 != 1:
                    continue
                r = line.strip()
                if len(r) < min_uzunluk:
                    continue
                n += 1
                if n > ornek:
                    break
                okumalar["o%d" % n] = r
    except OSError:
        pass
    # Batch alignment: starting a process per read on the command line backend would be
    # unacceptably slow.
    skor_h = {ad: max((x.mlen for x in hl), default=0)
              for ad, hl in A_h.map_toplu(okumalar)}
    skor_r = {ad: max((x.mlen for x in hl), default=0)
              for ad, hl in A_r.map_toplu(okumalar)}
    k = sum(1 for ad in okumalar if skor_h.get(ad, 0) > skor_r.get(ad, 0))
    n = len(okumalar)
    _BULASMA[anahtar] = (k, n)
    return (k, n)


def okuma_taramasi(path, F, R, prod_min, prod_max, limit):
    """Pair verification on the raw reads.

        In two stages. First the last 12 bases of the 3' end are searched quickly with
        str.find, which runs at C speed and discards most of the reads. Then, only on the
        candidate reads, the tested find_bindings and product_len functions of
        steps/design_group_primers.py are run, so the binding rule and the product
        coordinate are evaluated with exactly the same code as in the design stage.

        Returns: (total_reads, F_bound, R_bound, with_product)

    """
    pa = _PArg()
    pa.prod_min, pa.prod_hard_max = prod_min, prod_max
    K = KURAL.tail_len
    tot = f_hit = r_hit = both = 0
    # In the old version the gate required an exact match of F[-12:]. That was STRICTER
    # than the rule's own text: the rule allows one mismatch in the last five bases and
    # three in total, so even a single mismatch in the -12..-3 range left a valid binding
    # outside the scan. Measured: at an ONT-like 5% error, 15% of the real bindings were
    # lost, and at 8% error, 30%; the competitor product ratio came out systematically
    # LOW, which means non-specific pairs were PASSING. The gate is now built from THE
    # SAME core set that find_noidx uses, so it can never discard a read the rule would
    # have let through.

    # Building a 5-mer index per read was the bottleneck: on a 3700 base read, a
    # dictionary of ~3700 entries was created on every call. Since a read is a single
    # string, scanning the core variants with str.find instead of an index is far
    # cheaper; the result must be exactly the same as find_bindings, and that is tested
    # below.
    varyant = {}

    def _var(oligo):
        v = varyant.get(oligo)
        if v is None:
            v = sorted(G.seed_variants(oligo[-K:], KURAL.tail_max_mm,
                                       KURAL.exact_last))
            varyant[oligo] = v
        return v

    def find_noidx(oligo, seq):
        L, n = len(seq), len(oligo)
        saf = all(c in "ACGT" for c in oligo)
        hits, gor = [], set()
        for v in _var(oligo):
            st = 0
            while True:
                pos = seq.find(v, st)
                if pos < 0:
                    break
                st = pos + 1
                start = pos - (n - K)
                end = start + n
                if end > L:
                    continue
                j0 = max(0, -start)
                if n - j0 < KURAL.min_overlap:
                    continue
                # The fast path: raw reads hold only A, C, G, T and N.
                # If the oligo is pure ACGT as well, a direct character comparison is
                # enough instead of base_match's set intersection, and it gives the same
                # result (N always counts as a mismatch). Because the set operation is
                # called hundreds of times per read, the difference is large.
                mm = 0
                ok = True
                if saf:
                    for j in range(j0, n):
                        if oligo[j] != seq[start + j]:
                            mm += 1
                            if mm > KURAL.total_max_mm:
                                ok = False
                                break
                else:
                    for j in range(j0, n):
                        if not G.base_match(oligo[j], seq[start + j]):
                            mm += 1
                            if mm > KURAL.total_max_mm:
                                ok = False
                                break
                if ok:
                    key = (end - 1, mm)
                    if key not in gor:
                        gor.add(key)
                        hits.append(key)
        return hits

    def baglanmalar(oligo, seq, seqrc):
        return dict(L=len(seq), plus=find_noidx(oligo, seq),
                    minus=find_noidx(oligo, seqrc))

    kapiF = _var(F)
    kapiR = _var(R)

    def _kapi_gecer(seq, seqrc):
        """THE SAME core set as find_noidx. The gate discards only the reads in which
                no core variant is found at all; on those reads find_noidx returns empty by
                definition too, so the gate does not change the result.

        """
        for v in kapiF:
            if v in seq or v in seqrc:
                return True
        for v in kapiR:
            if v in seq or v in seqrc:
                return True
        return False

    for seq, seqrc in okumalari_al(path, limit):
        tot += 1
        if not _kapi_gecer(seq, seqrc):
            continue
        bf = baglanmalar(F, seq, seqrc)
        br = baglanmalar(R, seq, seqrc)
        fb = bool(bf["plus"] or bf["minus"])
        rb = bool(br["plus"] or br["minus"])
        if fb:
            f_hit += 1
        if rb:
            r_hit += 1
        if fb and rb:
            p = G.product_len(bf, br, len(F), len(R), pa, pmax=prod_max)
            if p is not None:
                both += 1
    return tot, f_hit, r_hit, both


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--candidates", required=True, help="08'in output directory")
    p.add_argument("--pt", required=True, help="PrimerTasarlama kok directory")
    p.add_argument("--out", required=True)
    p.add_argument("--targets", default=os.path.join(HERE, "targets.tsv"))
    p.add_argument("--top", type=int, default=15,
                   help="number of best candidates tested per target")
    p.add_argument("--max-reads", type=int, default=20000,
                   help='at most this many reads scanned per taxon; 0 means '
                        'unlimited. When the reads are trimmed it is written '
                        'into the log.')
    p.add_argument("--min-member-fraction", type=float, default=0.5,
                   help="at least this fraction of target members must be confirmed in the raw reads")
    p.add_argument("--consensus", default=None,
                   help='the consensus directory. When it is given, the cross '
                        'contamination in the competitor bins is measured and '
                        'the "a product in a competitor" verdict is decided '
                        'against that measured leak.')
    p.add_argument("--prod-min", type=int, default=70,
                   help="must match the value used at the design stage")
    p.add_argument("--prod-hard-max", type=int, default=300,
                   help="must match the value used at the design stage")
    p.add_argument("--contamination-sample", type=int, default=400)
    p.add_argument("--contamination-min-reads", type=int, default=100,
                   help="when the carry-over measurement rests on fewer reads than this "
                        "sizintinin esigi acmasina izin verilmez")
    p.add_argument("--leak-cap", type=float, default=0.15,
                   help="the leakage threshold may open up to this value at most")
    p.add_argument("--competitor-wilson-max", type=float, default=0.02,
                   help='a candidate is dropped when the Wilson lower bound '
                        'of the fraction of competitor reads giving a product '
                        'exceeds this value')
    p.add_argument("--within-strain-diff-max", type=float, default=0.40,
                   help='the largest difference allowed between the binding '
                        'rates of the two primers; a candidate that exceeds '
                        'it is penalised')
    p.add_argument("--atla-mfe", action="store_true")
    p.add_argument("--atla-blast", action="store_true")
    p.add_argument("--mfe", default=None, help='the path of the mfeprimer '
                                               'binary')
    p.add_argument("--rerun", action="store_true")
    return p.parse_args()


def main():
    a = get_args()
    os.makedirs(a.out, exist_ok=True)
    _LOG[0] = os.path.join(a.out, "ozgulluk.log")
    CKPT = os.path.join(a.out, "checkpoint.json")
    parmak = _girdi_parmak_izi(a)
    ckpt = {}
    if os.path.exists(CKPT) and not a.rerun:
        try:
            ham = json.load(open(CKPT, encoding="utf-8"))
            eski = ham.get("_girdi_parmak_izi") if isinstance(ham, dict) else None
            kayitlar = {k: v for k, v in ham.items()
                        if not k.startswith("_")} if isinstance(ham, dict) else {}
            if eski is None:
                log(u'the checkpoint CARRIES NO fingerprint (left over from an old version), it is ignored')
            elif eski != parmak:
                log(u'the checkpoint fingerprint DOES NOT MATCH, it is ignored')
                log(u'   recorded: %s' % eski)
                log(u'   current : %s' % parmak)
                log(u'   The candidate files or the thresholds have changed; the old verification results will not be reused.')
            else:
                ckpt = kayitlar
                log(u'the checkpoint is valid: %d targets are already processed, they will be skipped'
                    % len(ckpt))
        except Exception as e:
            log(u'the checkpoint could not be read (%s)' % e)
    t0 = time.time()
    log(u'start. candidates=%s' % a.candidates)

    mfe = a.mfe or os.path.join(a.pt, "tools", "mfeprimer")
    kullan_mfe = (not a.atla_mfe) and os.path.exists(mfe) and os.access(mfe, os.X_OK)
    if not a.atla_mfe and not kullan_mfe:
        log(u'WARNING: mfeprimer was not found or is not executable (%s), the external database step will be skipped' % mfe)
    # External database specificity is NOW done in external_databases.py.
    # The flags here stand only for backward compatibility; which step runs where is
    # written openly into the log, so that no silent difference is left between the
    # documentation and the code.
    kullan_blast = (not a.atla_blast) and bool(
        subprocess.run(["bash", "-c", "command -v blastn"],
                       capture_output=True).stdout.strip())
    if not a.atla_blast and not kullan_blast:
        log(u'WARNING: blastn was not found, the second measurement will be skipped')

    # fastq envanteri: sinif_taxid -> yol
    fq = {}
    for p in glob.glob(os.path.join(a.pt, "fastq files", "*", "*.fastq")):
        grp = os.path.basename(os.path.dirname(p))
        m = re.search(r"reads[-_](\d+)", os.path.basename(p))
        if m:
            fq[(re.split(r"[-_]", grp)[0], grp, m.group(1))] = p
    log(u'fastq inventory: %d files' % len(fq))

    # the consensus inventory: (group, taxid) -> path. For the cross contamination measurement.
    kons = {}
    if a.consensus:
        for p2 in glob.glob(os.path.join(a.consensus, "*_konsensus.fasta")):
            m = re.match(r"((?:A1|A2|B|F1|F2)-\d+)_(\d+)_", os.path.basename(p2))
            if m:
                kons[(m.group(1), m.group(2))] = p2
        log(u'consensus inventory: %d files' % len(kons))
        log(alignment.durum())
        if not MAPPY:
            log(u'WARNING: there is no alignment backend, cross contamination cannot be measured. The fixed threshold (--competitor-wilson-max %.3f) will be used.'
                % a.competitor_wilson_max)
    else:
        log(u'--consensus was not given: cross contamination will not be measured, a fixed threshold will be used (--competitor-wilson-max %.3f)' % a.competitor_wilson_max)

    def _kons_of(fq_yolu):
        b = os.path.basename(fq_yolu)
        g = os.path.basename(os.path.dirname(fq_yolu))
        m = re.search(r"reads[-_](\d+)", b)
        return kons.get((g, m.group(1))) if m else None

    # the target definitions
    hedefler = {}
    for line in open(a.targets, encoding="utf-8"):
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        if p[0] == "karar":
            continue
        hedefler[p[1]] = dict(karar=p[0], duzey=p[2], inn=p[3],
                              haric=p[4] if len(p) > 4 else "")

    # The indistinguishable taxon pairs measured by batch_design.py. Unless the same
    # exclusion is applied here too, a taxon batch_design.py removed from the competitor
    # list comes back in this stage and the target's own sequence counts as a competitor.
    # The taxa batch_design.py excluded because of a broken consensus. Because this stage
    # works from the fastq inventory, it was taking the same taxon back as a competitor
    # (and as a member on domain targets); the two stages have to work on the same sets.
    # A (group, taxid) pair. Excluding by taxon alone would be wrong: a taxon's consensus
    # can be empty in one sample and sound in another.
    # The exclusion applies ONLY to the member set; in the competitor set the raw reads
    # of these bins are valuable and removing them weakens the specificity check. An
    # empty consensus does not mean the reads are gone, only that no consensus could be
    # built.
    dislanan = set()
    dl = os.path.join(a.candidates, "dislanan_takson.tsv")
    if os.path.exists(dl):
        # The columns are read BY HEADER, not BY POSITION. In earlier versions this file had
        # 4 columns (taxid, etiket, uzunluk, kapsanan); code reading by position took the
        # (taxid, etiket) pair for (group, taxid) in that file, built a set matching nothing,
        # and THE EXCLUSION WAS SILENTLY DISABLED. Stopping when we do not recognise the
        # header is better than carrying on with the wrong set.
        with open(dl, encoding="utf-8") as fh:
            satirlar = [l.rstrip("\n") for l in fh
                        if l.strip() and not l.startswith("#")]
        if satirlar:
            basliklar = satirlar[0].split("\t")
            if "grup" not in basliklar or "taxid" not in basliklar:
                sys.exit(u'ERROR: the header of %s was not recognised: %s\nExpected columns: grup, taxid, etiket, uzunluk, kapsanan.\nThis file may be left over from an old batch_design.py version; run batch_design.py again.' % (dl, basliklar))
            ig, it = basliklar.index("grup"), basliklar.index("taxid")
            for satir in satirlar[1:]:
                p2 = satir.split("\t")
                if len(p2) > max(ig, it) and p2[ig].strip() and p2[it].strip():
                    dislanan.add((p2[ig].strip(), p2[it].strip()))
        log(u'excluded by batch_design.py (group, taxon): %d -> taken out of the member set, KEPT in the competitor set with its raw reads' % len(dislanan))
        for g, t in sorted(dislanan):
            log("   %s %s" % (g, t))

    ayirt = {}
    ae = os.path.join(a.candidates, "ayirt_edilemez.tsv")
    if os.path.exists(ae):
        for r in csv.DictReader(open(ae, encoding="utf-8"), delimiter="\t"):
            ayirt.setdefault((r["sinif"], r["taxid1"]), set()).add(r["taxid2"])
            ayirt.setdefault((r["sinif"], r["taxid2"]), set()).add(r["taxid1"])
        log(u'the indistinguishable pair table was read: %s (%d records)'
            % (ae, sum(len(v) for v in ayirt.values()) // 2))
    else:
        log(u'WARNING: %s is missing, indistinguishable taxa will not be filtered out' % ae)

    sonuc = []
    tsvler = sorted(glob.glob(os.path.join(a.candidates, "*__*.tsv")))
    log(u'NOTE: external database specificity (mfeprimer and blastn) does NOT run in this script, it runs in the external_databases.py step.')
    log(u'candidate files to process: %d' % len(tsvler))

    for ti, tsv in enumerate(tsvler, 1):
        etiket = os.path.basename(tsv)[:-4]
        hedef, sinif = etiket.rsplit("__", 1)
        if etiket in ckpt and not a.rerun:
            log(u'[%d/%d] SKIPPED (checkpoint) %s' % (ti, len(tsvler), etiket))
            sonuc.extend(ckpt[etiket])
            continue
        th = time.time()
        rows = list(csv.DictReader(open(tsv, encoding="utf-8"), delimiter="\t"))
        if not rows:
            log(u'[%d/%d] %s is empty, skipped' % (ti, len(tsvler), etiket))
            ckpt[etiket] = []
            continue
        h = hedefler.get(hedef, {})
        in_t = set(x for x in h.get("inn", "").split(",") if x and x != "*")
        haric = set(x for x in h.get("haric", "").split(",") if x)
        alan = h.get("inn", "").startswith("*")
        # The member fastq files are grouped PER TAXON. In the old version the verification
        # counter incremented PER FILE; when the same taxon had files from several years, a
        # member taxon that was never verified could be lost inside the ratio. Measured: when
        # only two of five member taxa gave a product for Asetoklastik_metanojenler in A2,
        # the file based ratio PASSED at 7/10=0.70, while the taxon based ratio 2/5=0.40
        # should have failed.
        if alan:
            # A domain target (*A, *B, *F): the members are every class starting
            # with the same letter, the competitors are the classes of the OTHER
            # domains. batch_design.py builds it exactly that way; the old version
            # of this stage left the competitor list empty and the universal primers
            # PASSED without a single competitor being tested. The 'haric' column
            # also applies to the member set on domain targets.
            # EXACTLY the same setup as batch_design.py: the members are ONLY the
            # taxa in this amplicon class, the competitors are the classes of other
            # DOMAINS.
            harf = sinif[0]
            uye_ikili = [(t, v) for (s, g, t), v in fq.items()
                         if s == sinif and t not in haric
                         and (g, t) not in dislanan]
            rakip_fq = [v for (s, g, t), v in fq.items() if s[0] != harf]
        else:
            uye_ikili = [(t, v) for (s, g, t), v in fq.items()
                         if s == sinif and t in in_t and (g, t) not in dislanan]
            rakip_fq = [v for (s, g, t), v in fq.items()
                        if s == sinif and t not in in_t and t not in haric]
        uye_fq = [v for _, v in uye_ikili]
        uye_takson = {}
        for t, v in uye_ikili:
            uye_takson.setdefault(t, []).append(v)
        atilan_tx = set()
        if ayirt and in_t:
            temiz = []
            for v in rakip_fq:
                tx = None
                for (s2, g2, t2), v2 in fq.items():
                    if v2 == v:
                        tx = t2
                        break
                if tx and (ayirt.get((sinif, tx), set()) & in_t):
                    atilan_tx.add(tx)
                else:
                    temiz.append(v)
            rakip_fq = temiz
        if atilan_tx:
            log(u'      taken out of the competitors (indistinguishable from the target): %s'
                % ", ".join(sorted(atilan_tx)))
        log(u'[%d/%d] %-46s candidates=%d member_taxa=%d member_fastq=%d competitor_fastq=%d (at most %d reads per fastq)'
            % (ti, len(tsvler), etiket, len(rows), len(uye_takson),
               len(uye_fq), len(rakip_fq), a.max_reads))
        if not uye_fq:
            log(u'      no member fastq was found, the raw read verification will be skipped')

        uye_kons = [k for k in (_kons_of(x) for x in uye_fq) if k]

        def _bulasma_getir(rakip_yolu):
            """Measures how many of the reads in a competitor bin actually belong
                        to THE TARGET. The target member with the highest leakage is taken.

            """
            rk = _kons_of(rakip_yolu)
            if not (uye_kons and rk and MAPPY):
                return (0, 0)
            en = (0, 0)
            for hk in uye_kons:
                k, n = capraz_bulasma(hk, rk, rakip_yolu, a.contamination_sample)
                if n and (not en[1] or k / n > en[0] / max(1, en[1])):
                    en = (k, n)
            return en

        rows.sort(key=lambda r: float(r.get("ceza", 9e9)))
        gecen = []
        for ri, r in enumerate(rows[:a.top], 1):
            F, R = r["ileri_dizi"], r["geri_dizi"]
            # The product window must be THE SAME as in the design stage; if it
            # differs, borderline products are judged differently in the two stages.
            pmin, pmax = a.prod_min, a.prod_hard_max
            uye_dogru = 0
            f_or, r_or = [], []
            dogrulanmayan = []
            for tx, yollar in sorted(uye_takson.items()):
                tx_urun = 0
                for p in yollar:
                    tot, fh_, rh, both = okuma_taramasi(p, F, R, pmin, pmax,
                                                        a.max_reads)
                    if tot:
                        f_or.append(fh_ / tot)
                        r_or.append(rh / tot)
                        tx_urun += both
                if tx_urun > 0:
                    uye_dogru += 1
                else:
                    dogrulanmayan.append(tx)
            uye_orani = (uye_dogru / len(uye_takson)) if uye_takson else None
            # the proportion of competitor reads giving a product, the Wilson lower bound
            rak_w = 0.0
            rak_detay = []
            rakipte_gercek = False
            for p in rakip_fq:
                tot, fh_, rh, both = okuma_taramasi(p, F, R, pmin, pmax, a.max_reads)
                w = wilson_alt(both, tot)
                # The UPPER bound of the measured bin leakage. If the product ratio is
                # below this, the reads giving a product can be explained by TARGET
                # reads that fell into the wrong bin; it is not evidence that the
                # competitor itself amplified.
                bk, bn = _bulasma_getir(p)
                # The Wilson UPPER bound approaches 1 with few reads. Left unguarded,
                # on a competitor bin of 10 reads the threshold rises to 0.28 and on
                # one of 5 reads to 0.43, so a competitor giving a product in half its
                # reads counts as "clean". So leakage can only raise the threshold when
                # there are enough reads, and the amount it raises it by is capped.
                if bn >= a.contamination_min_reads:
                    sizinti_ust = min(wilson_ust(bk, bn), a.leak_cap)
                    sizinti_not = "%.4f" % sizinti_ust
                else:
                    sizinti_ust = 0.0
                    sizinti_not = "yetersiz_okuma(%d)" % bn
                esik = max(a.competitor_wilson_max, sizinti_ust)
                if w > esik:
                    rakipte_gercek = True
                rak_detay.append("%s=%d/%d(W%.4f,sizinti<%s)"
                                 % (os.path.basename(p)[:18], both, tot, w,
                                    sizinti_not))
                rak_w = max(rak_w, w)
            fmin = min(f_or) if f_or else 0.0
            rmin = min(r_or) if r_or else 0.0
            sus_fark = abs((sum(f_or) / len(f_or) if f_or else 0)
                           - (sum(r_or) / len(r_or) if r_or else 0))
            durum = []
            # A candidate on which not a single raw read was scanned cannot count as
            # "passed". In the old version, when uye_fq or rakip_fq was empty no flag
            # was added at all and the candidate was delivered as PASSED.
            if not uye_fq:
                durum.append("uye_okumasi_yok")
            elif uye_orani is None or uye_orani < a.min_member_fraction:
                durum.append("numuneden_dogrulanamadi")
            if not rakip_fq:
                durum.append("rakip_sinanmadi")
            elif rakipte_gercek:
                durum.append("rakipte_urun")
            if sus_fark > a.within_strain_diff_max:
                durum.append("sus_ici_degisken")
            r2 = dict(r)
            r2.update(hedef=hedef, sinif=sinif, karar=h.get("karar", ""),
                      uye_dogrulanan=uye_dogru, uye_toplam=len(uye_takson),
                      dogrulanmayan_uye=",".join(dogrulanmayan),
                      uye_orani=round(uye_orani, 3) if uye_orani is not None else "",
                      rakip_wilson=round(rak_w, 5),
                      rakip_detay=";".join(rak_detay)[:200],
                      ileri_baglanma_min=round(fmin, 3),
                      geri_baglanma_min=round(rmin, 3),
                      sus_ici_fark=round(sus_fark, 3),
                      ozgulluk_durum=",".join(durum) if durum else "GECTI")
            gecen.append(r2)
            if len(durum) == 0:
                log(u'      [%d] GECTI  F=%s R=%s  members=%d/%d taxa competitorW=%.4f'
                    % (ri, F, R, uye_dogru, len(uye_takson), rak_w))
            if sum(1 for g in gecen if g["ozgulluk_durum"] == "GECTI") >= 5:
                log(u'      five valid candidates were found, stopped for this target')
                break
        ckpt[etiket] = gecen
        sonuc.extend(gecen)
        ckpt["_girdi_parmak_izi"] = parmak
        with open(CKPT, "w", encoding="utf-8") as cf:
            json.dump(ckpt, cf, ensure_ascii=False)
        yaz(a.out, sonuc)
        log(u'      done, %.1f s, passed=%d/%d'
            % (time.time() - th,
               sum(1 for g in gecen if g["ozgulluk_durum"] == "GECTI"), len(gecen)))

    yaz(a.out, sonuc)
    ok = sum(1 for s in sonuc if s["ozgulluk_durum"] == "GECTI")
    log(u'TOTAL: %d candidates tested, %d of them passed every rule' % (len(sonuc), ok))
    log(u'total time: %.1f minutes' % ((time.time() - t0) / 60))
    log(u'output: %s' % os.path.join(a.out, "final_primers.tsv"))


def yaz(out, sonuc):
    """The file is rewritten even when the result is empty. In the old version it
        returned early, so the previous run's stale final_primers.tsv stayed on disk and
        the Excel export took it for a valid one and published it.

    """
    if not sonuc:
        with open(os.path.join(out, "final_primers.tsv"), "w",
                  encoding="utf-8") as fh:
            fh.write("karar\thedef\tsinif\tozgulluk_durum\n")
        return
    cols = ["karar", "hedef", "sinif", "ozgulluk_durum",
            "ileri_dizi", "ileri_tm", "ileri_gc", "ileri_baslangic", "ileri_uzunluk",
            "geri_dizi", "geri_tm", "geri_gc", "geri_baslangic", "geri_uzunluk",
            "tm_farki", "urun_min", "urun_maks",
            "ileri_iupac_cozulen", "geri_iupac_cozulen", "lokus_varyant_sayisi",
            "uye_dogrulanan", "uye_toplam", "uye_orani", "rakip_wilson",
            "ileri_baglanma_min", "geri_baglanma_min", "sus_ici_fark",
            "yetim_primer", "heterodimer_dg", "ceza", "rakip_detay"]
    cols = [c for c in cols if any(c in s for s in sonuc)]
    with open(os.path.join(out, "final_primers.tsv"), "w", newline="",
              encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(sorted(sonuc, key=lambda s: (s["ozgulluk_durum"] != "GECTI",
                                                 s.get("karar", ""), s["hedef"],
                                                 float(s.get("ceza", 9e9)))))


if __name__ == "__main__":
    main()
