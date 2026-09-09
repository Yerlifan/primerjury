#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==== MEETING NOTE ====
# WHAT IT IS FOR   : PAK polish for EVERY bin: split the bin's reads into genus populations (ITS
#                    window for fungi, 16S window for archaea and bacteria), polish a consensus from
#                    the DOMINANT population only (medoid read as template, minimap2 + samtools
#                    consensus, two rounds) and decide it with the group's own rule (three loci for
#                    fungi, the 16S ladder for the rest). Generalised from fungal_bin_identity.py on
#                    2026-09-06: "not only the fungi, every organism".
# INPUT            : --root; "fastq files"/<lib>/ bin reads; REFERENCE_DB (or a fast copy)
# OUTPUT           : referans_konsensus/pak_polish/consensus/<bin>.fasta (a new candidate set),
#                    referans_konsensus/pak_polish/population/<bin>.tsv, PAK_POLISH.tsv
# HOW TO RUN IT    : python3 verification/population_polish.py --root . [--groups A1,A2,B,F1,F2]
#                    python3 verification/population_polish.py --self-test
# WHY IT IS LIKE THIS : Measured in the study: 6 of 59 archaeal/bacterial bins had a dominant
#                    population under 80 per cent (Bacteroidales mixtures), and a ten-read archaeal
#                    bin went from "cannot be named" at 94.2 to Methanofollis ethanolicus at 99.86.
# =======================
u"""PAK polish for every library. The fungal module holds the machinery; this file adds
the group logic (which locus window, which databases, which decision rule) and the
batch driver over all groups.

  group      window located with        best record per read over          decision
  F1, F2     fungi.ITS.fna              RefSeq ITS + UNITE                 three loci (locus_decision)
  A1, A2     archaea.16S.fna            archaea.16S + SILVA SSU NR99       16S ladder (target_identity)
  B          bacteria.16S.fna           bacteria.16S + SILVA SSU NR99      16S ladder (target_identity)

The polished sequence lands as one more candidate set; select_consensus weighs it
with the same read-support rule as every other candidate.
"""
from __future__ import print_function
import argparse
import collections
import io
import os
import shutil
import sys
import tempfile
import time

KOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _d in (os.path.join(KOK, 'verification'), os.path.join(KOK, 'steps')):
    if _d not in sys.path:
        sys.path.insert(0, _d)
import fungal_bin_identity as F
from identity_verification import is_concatemer, CONCATEMER_FACTOR                                             # noqa: E402
from fungal_bin_identity import (sample_reads, write_fasta, read_fasta_seq, medoid,   # noqa: E402
                                 polish, best_hits, best_hits_mm2, its_window,
                                 window_file, bin_files, db_path, populations,
                                 MIN_POPULATION, ITS_DB, blast)
from locus_decision import (MANTAR_LOKUSLARI, lokus_karari, birlestir, raporlanan_yontem, kingdom_gate,   # noqa: E402
                            LOKUS_KUMELERI, BASAMAK_ARKE)
from identity_verification import ADSIZ_JETONLARI                           # noqa: E402
from target_identity import ad_ayikla, cins_epitet, basamaktan_sec          # noqa: E402

SET_DIR = os.path.join('referans_konsensus', 'pak_polish', 'consensus')
POP_DIR = os.path.join('referans_konsensus', 'pak_polish', 'population')
GROUPS = ('A1', 'A2', 'B', 'F1', 'F2')
SEP = u'@@'
SSU_DB = ['archaea.16S.fna', 'bacteria.16S.fna', 'SILVA_138.2_SSURef_NR99.fasta']
MULTI = {'F1': 'COK_LOKUS', 'F2': 'COK_LOKUS', 'A2': 'COK_LOKUS_ARKE'}   # groups decided over several loci
ALT_DIR = os.path.join('referans_konsensus', 'pak_polish', 'sub')
MIN_SUB = 10          # a sub-population needs this many reads ...
SUB_SHARE = 15.0      # ... and this share of its parent population


def level(name):
    name = name or u''
    if u'cf.' in name:
        return u'cf.'
    if name.endswith(u'sp.'):
        return u'genus'
    if u' ' in name and name not in (u'adlandırılamıyor', u'eşleşme yok', u'cannot be named', u'no match'):
        return u'species'
    return u'unnamed'


def final_polish(tpl, pop_fa, kd, a):
    u"""With --medaka (the default) medaka runs on top of the samtools consensus.
    Returns (sequence or None, polish label)."""
    if not getattr(a, 'medaka', False):
        return None, u'samtools'
    if not F.medaka_present():
        return None, u'samtools (no medaka environment)'
    seq = F.medaka_polish(tpl, pop_fa, os.path.join(kd, 'medaka.fa'), a.threads)
    if not seq or len(seq) < 200:
        return None, u'samtools (medaka failed)'
    return seq, u'samtools+medaka'


def group_locus(group):
    return u'ITS' if group in ('F1', 'F2') else u'SSU'


def window_db(group):
    if group in ('F1', 'F2'):
        return ITS_DB[0]
    return 'archaea.16S.fna' if group in ('A1', 'A2') else 'bacteria.16S.fna'


def population_dbs(group):
    if group in ('F1', 'F2'):
        return list(ITS_DB)
    return [('archaea.16S.fna' if group in ('A1', 'A2') else 'bacteria.16S.fna'),
            'SILVA_138.2_SSURef_NR99.fasta']


def batch_hits(bin_reads, root, a, tmp, blast_fn=None):
    u"""Per group: all reads in one file, the locus window, the best record per read."""
    per_bin = {label: {} for label in bin_reads}
    windowed = collections.Counter()
    groups = collections.defaultdict(list)
    for label in bin_reads:
        groups[label.split(u'-')[0]].append(label)
    for group in sorted(groups):
        allreads = [(u'%s%s%s' % (label, SEP, k), d) for label in groups[group] for k, d in bin_reads[label]]
        reads_fa = os.path.join(tmp, 'reads_%s.fa' % group)
        write_fasta(reads_fa, allreads)
        t0 = time.time()
        rdb = db_path(root, window_db(group), a.db_dir, quiet=True) if root is not None else window_db(group)
        window = its_window(reads_fa, rdb, a.threads) if rdb else {}
        win_fa = os.path.join(tmp, 'windows_%s.fa' % group)
        n_window = window_file(allreads, window, win_fa)
        dbs = population_dbs(group)
        t1 = time.time()
        best = (best_hits(win_fa, dbs, root, a.threads, a.db_dir, blast_fn) if blast_fn
                else best_hits_mm2(win_fa, dbs, root, a.threads, a.db_dir))
        print(u'  %s: %d bins, %d reads; %s window %d (%.0f s); hits %d (%.0f s; %s)'
              % (group, len(groups[group]), len(allreads), group_locus(group), n_window, t1 - t0,
                 len(best), time.time() - t1, u' + '.join(dbs)))
        sys.stdout.flush()
        for key, v in best.items():
            label, k = key.split(SEP, 1)
            per_bin[label][k] = v
        for key in window:
            windowed[key.split(SEP, 1)[0]] += 1
    return per_bin, windowed


def sub_populations(label, reads, count, members):
    u"""SUB-POPULATIONS INSIDE A BIN (2026-09-06). A bin can hold several species
    (Petriella musispora + guttulata; Proteiniphilum mixtures). Two sources: (a) species
    groups inside the dominant genus (by the best record's real epithet), (b) the second
    genus population. Each with >= MIN_SUB reads and >= SUB_SHARE per cent is polished on
    its own. The parent bin's result does not change; sub-populations go to a separate
    report (a primer target of their own). Returns [(tag, source, reads)]."""
    out = []
    seq_of = dict(reads)
    total = sum(count.values())
    if not total:
        return out
    g1, n1 = count.most_common(1)[0]
    groups = collections.defaultdict(list)
    for k, _p, _a, name in members[g1]:
        c, e = cins_epitet(name)
        if e and e.lower() not in F.PLACEHOLDER:
            groups[u'%s %s' % (c, e)].append(k)
    big = [(t, ks) for t, ks in groups.items() if len(ks) >= MIN_SUB and 100.0 * len(ks) / n1 >= SUB_SHARE]
    if len(big) >= 2:
        for t, ks in sorted(big, key=lambda x: -len(x[1])):
            out.append((u'%s~%s' % (label, t.replace(u' ', u'_')), u'species group inside the dominant genus (%d/%d)' % (len(ks), n1),
                        [(k, seq_of[k]) for k in ks]))
    if len(count) > 1:
        g2, n2 = count.most_common(2)[1]
        if n2 >= MIN_SUB and 100.0 * n2 / total >= SUB_SHARE:
            out.append((u'%s~%s' % (label, g2.replace(u' ', u'_')), u'second genus population (%d/%d)' % (n2, total),
                        [(k, seq_of[k]) for k, _p, _a, _n in members[g2]]))
    return out


def polish_sub(root, tag, group, pop, a, tmp):
    kd = os.path.join(tmp, tag.replace(u'~', u'__'))
    os.makedirs(kd, exist_ok=True)
    med = medoid(pop)
    pop_fa = os.path.join(kd, 'population.fa')
    write_fasta(pop_fa, pop)
    tpl = os.path.join(kd, 'template0.fa')
    write_fasta(tpl, [(med[0], med[1])])
    seq = med[1]
    for rnd in range(1, a.rounds + 1):
        out = os.path.join(kd, 'round%d.fa' % rnd)
        seq, _n = polish(tpl, pop_fa, out, a.threads)
        if len(seq) < 200:
            seq = read_fasta_seq(tpl)
            break
        write_fasta(out, [(tag, seq)])
        tpl = out
    seq2, how = final_polish(tpl, pop_fa, kd, a)
    if seq2:
        seq = seq2
    return dict(bin=tag, group=group, locus=group_locus(group), reads=len(pop), pop_n=len(pop),
                template=med[0], template_similarity=med[2], bp=len(seq), sequence=seq, polish=how)


def polish_bin(root, label, reads, best, n_window, a, tmp):
    group = label.split(u'-')[0]
    key = group_locus(group)
    kd = os.path.join(tmp, label)
    if not os.path.isdir(kd):
        os.makedirs(kd)
    count, members, scount, spid = populations(best, key)
    seq_of = dict(reads)
    pop_dir = os.path.join(root, POP_DIR)
    if not os.path.isdir(pop_dir):
        os.makedirs(pop_dir)
    rows = [u'read\tbp\tgenus\tname\tidentity %\talignment bp\tdatabase']
    for k, d in reads:
        if k in best:
            pid, aln, title, db = best[k]
            name = ad_ayikla(title) or u''
            rows.append(u'%s\t%d\t%s\t%s\t%.2f\t%d\t%s'
                        % (k, len(d), cins_epitet(name)[0] if name else u'', name, pid, aln, db))
        else:
            rows.append(u'%s\t%d\t\t\t\t\t' % (k, len(d)))
    io.open(os.path.join(pop_dir, label + '.tsv'), 'w', encoding='utf-8',
            newline='\n').write(u'\n'.join(rows) + u'\n')
    total = sum(count.values())
    r = dict(bin=label, group=group, locus=key, reads=len(reads), assigned=total, windowed=n_window)
    if scount:
        s1, n1 = scount.most_common(1)[0]
        pids = sorted(spid[s1])
        r.update(read_species=s1, species_share=100.0 * n1 / sum(scount.values()),
                 species_n=n1, species_median=pids[len(pids) // 2])
    if total:
        g1, n1 = count.most_common(1)[0]
        g2, n2 = count.most_common(2)[1] if len(count) > 1 else (u'-', 0)
        r.update(genus=g1, share=100.0 * n1 / total, genus2=g2, share2=100.0 * n2 / total)
        pop = [(k, seq_of[k]) for k, _p, _a, _n in members[g1]]
        pop_note = u'dominant population %s %d/%d' % (g1, n1, total)
    else:
        pop = list(reads)
        r.update(genus=u'-', share=0.0, genus2=u'-', share2=0.0)
        pop_note = u'no %s hit above the genus threshold; all reads form one population' % key
    if len(pop) < MIN_POPULATION and total < 0.1 * len(reads):
        pop = list(reads)
        pop_note = u'%s; %s hits %d/%d (<10%%), all reads form one population' % (pop_note, key, total, len(reads))
    if len(pop) < MIN_POPULATION:
        r['note'] = u'%s: %d reads < %d, no consensus made' % (pop_note, len(pop), MIN_POPULATION)
        return r
    if getattr(a, 'sub', True) and total:
        r['sub'] = sub_populations(label, reads, count, members)
    med = medoid(pop)
    r.update(template=med[0], template_bp=len(med[1]), template_similarity=med[2])
    pop_fa = os.path.join(kd, 'population.fa')
    write_fasta(pop_fa, pop)
    tpl = os.path.join(kd, 'template0.fa')
    write_fasta(tpl, [(med[0], med[1])])
    seq, n_in = med[1], 0
    for rnd in range(1, a.rounds + 1):
        out = os.path.join(kd, 'round%d.fa' % rnd)
        seq, n_in = polish(tpl, pop_fa, out, a.threads)
        if len(seq) < 200:
            r['note'] = u'%s; round %d left %d bp, previous round kept' % (pop_note, rnd, len(seq))
            seq = read_fasta_seq(tpl)
            break
        write_fasta(out, [(label, seq)])
        tpl = out
    seq2, how = final_polish(tpl, pop_fa, kd, a)
    if seq2:
        seq, n_in = seq2, seq2.count(u'N')
    r.update(bp=len(seq), n_internal=n_in, sequence=seq, pop_n=len(pop), polish=how)
    return r


def batch_decide(results, root, a, tmp, blast_fn=None):
    u"""Fungal bins: three loci; the others: the 16S ladder of target_identity."""
    multi = {k: r for k, r in results.items() if r.get('sequence') and r['group'] in MULTI}
    others = {k: r for k, r in results.items() if r.get('sequence') and r['group'] not in MULTI}
    for key in sorted(set(MULTI[r['group']] for r in multi.values())):
        fungal = {k: r for k, r in multi.items() if MULTI[r['group']] == key}
        loci = LOKUS_KUMELERI[key]
        ladder = None if key == 'COK_LOKUS' else BASAMAK_ARKE
        query = os.path.join(tmp, 'final_%s.fa' % key)
        write_fasta(query, [(k, r['sequence']) for k, r in fungal.items()])
        hits = {k: {} for k in fungal}
        for lok, dbs, _an in loci:
            for name in dbs:
                db = db_path(root, name, a.db_dir) if root is not None else name
                if db:
                    for k, v in (blast_fn or blast)(query, db, threads=a.threads, top=500).items():
                        hits.setdefault(k, {}).setdefault(lok, []).extend(v)
        for label, r in fungal.items():
            dec = {lok: lokus_karari(hits.get(label, {}).get(lok), an, ad_ayikla, cins_epitet, ADSIZ_JETONLARI)
                   for lok, _d, an in loci}
            name, ident, note = birlestir(dec)
            gate = kingdom_gate(hits.get(label, {})) if key == 'COK_LOKUS' else None   # 2026-09-09
            if gate:
                name, ident, note = gate
            ra, rp, rl, rlok, rsecond, rgap = raporlanan_yontem(hits.get(label, {}), ad_ayikla, cins_epitet,
                                                                ADSIZ_JETONLARI, ladder)
            r.update(name=name, identity=ident, decision_note=note,
                     detail=u' | '.join(u'%s: %s %.2f%%/%d' % (l, dec[l]['ad'], dec[l]['kimlik'], dec[l]['hizalama'])
                                        for l, _d, _a in loci if l in dec),
                     reported=(u'%s (%.2f%%, %d bp, %s)' % (ra, rp, rl, rlok)) if ra else u'-')
    groups = collections.defaultdict(dict)
    for k, r in others.items():
        groups[r['group']][k] = r
    for group, bins in sorted(groups.items()):
        query = os.path.join(tmp, 'final_%s.fa' % group)
        write_fasta(query, [(k, r['sequence']) for k, r in bins.items()])
        hits = collections.defaultdict(list)
        for name in population_dbs(group):
            db = db_path(root, name, a.db_dir) if root is not None else name
            if db:
                for k, v in (blast_fn or blast)(query, db, threads=a.threads, top=500).items():
                    for x in v:
                        # (pid, aln, title, db, slen): the tuple the ladder expects
                        hits[k].append((x[0], x[1], x[4], os.path.basename(name), x[5] if len(x) > 5 else None))
        for label, r in bins.items():
            isb = sorted(hits.get(label, []), key=lambda x: (-x[0], -x[1]))
            res = basamaktan_sec({'16S': isb}, ['16S'], True)
            if res:
                name, note, pid, aln = res[0], res[1], res[2], res[3]
                r.update(name=name, identity=u'%.2f' % pid if pid is not None else u'-',
                         decision_note=note, detail=u'16S: %s %.2f%%/%d' % (name, pid or 0.0, aln or 0), reported=u'-')
            else:
                r.update(name=u'cannot be named', identity=u'-', decision_note=u'no 16S hit', detail=u'16S: no hit', reported=u'-')


HEADER = [u'bin', u'group', u'locus', u'sampled reads', u'locus-assigned reads', u'windowed reads',
          u'dominant genus', u'share %', u'second genus', u'second share %', u'read species', u'species share %',
          u'reads at species', u'species median identity %', u'template read', u'template bp',
          u'template k-mer similarity', u'polishing reads', u'consensus bp', u'internal N',
          u'DECISION', u'identity %', u'decision note', u'locus detail', u'reported method', u'note']


def make_row(r):
    def f(x, fmt=u'%s'):
        return (fmt % x) if x not in (None, u'') else u'-'
    return u'\t'.join([r['bin'], f(r.get('group')), f(r.get('locus')), f(r.get('reads')), f(r.get('assigned')),
                       f(r.get('windowed')), f(r.get('genus')), f(r.get('share'), u'%.0f'), f(r.get('genus2')),
                       f(r.get('share2'), u'%.0f'), f(r.get('read_species')), f(r.get('species_share'), u'%.0f'),
                       f(r.get('species_n')), f(r.get('species_median'), u'%.2f'), f(r.get('template')),
                       f(r.get('template_bp')), f(r.get('template_similarity'), u'%.3f'), f(r.get('pop_n')),
                       f(r.get('bp')), f(r.get('n_internal')), f(r.get('name')), f(r.get('identity')),
                       f(r.get('decision_note')), f(r.get('detail')), f(r.get('reported')), f(r.get('note'))])


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', default='.')
    ap.add_argument('--reads', type=int, default=150)
    ap.add_argument('--rounds', type=int, default=2)
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--groups', default=u','.join(GROUPS))
    ap.add_argument('--bins', default='')
    ap.add_argument('--db-dir', default=None)
    ap.add_argument('--out', default='PAK_POLISH.tsv')
    ap.add_argument('--redo', action='store_true')
    ap.add_argument('--no-sub', dest='sub', action='store_false', default=True,
                    help='skip the sub-population stage')
    ap.add_argument('--medaka', dest='medaka', action='store_true', default=True,
                    help=u'medaka on top of the samtools consensus (default; environment: %s)' % F.MEDAKA_ENV)
    ap.add_argument('--no-medaka', dest='medaka', action='store_false')
    ap.add_argument('--self-test', action='store_true')
    a = ap.parse_args(argv)
    if a.self_test:
        return self_test()
    root = os.path.abspath(a.root)
    groups = set(a.groups.split(u','))
    files = bin_files(root)
    bins = [b for b in sorted(files) if b.split(u'-')[0] in groups]
    if a.bins:
        want = set(a.bins.split(u','))
        bins = [b for b in bins if b in want]
    set_dir = os.path.join(root, SET_DIR)
    if not os.path.isdir(set_dir):
        os.makedirs(set_dir)
    out_path = os.path.join(root, a.out)
    previous = {}
    # with --bins the other bins' rows are kept whatever --redo says
    if os.path.exists(out_path) and (not a.redo or a.bins):
        for s in io.open(out_path, encoding='utf-8', errors='replace').read().splitlines()[1:]:
            if s.strip():
                previous[s.split(u'\t')[0]] = s
    rows, todo = {}, []
    for label in bins:
        if os.path.exists(os.path.join(set_dir, label + '.fasta')) and not a.redo and label in previous:
            rows[label] = previous[label]
        else:
            todo.append(label)
    if a.bins:
        for k, line in previous.items():
            if k not in bins:
                rows[k] = line
    if a.medaka and not F.medaka_present():
        print(u'  WARNING: --medaka is on but the environment is missing (%s); the samtools consensus stands' % F.MEDAKA_ENV)
    print(u'  bins: %d (done %d, to do %d; groups %s; reads %d, rounds %d, threads %d)'
          % (len(bins), len(bins) - len(todo), len(todo), u','.join(sorted(groups)), a.reads, a.rounds, a.threads))
    tmp = tempfile.mkdtemp(prefix='pak_polish_')
    levels = collections.Counter()

    def flush():
        io.open(out_path, 'w', encoding='utf-8', newline='\n').write(
            u'\t'.join(HEADER) + u'\n' + u'\n'.join(rows[k] for k in sorted(rows)) + u'\n')
    try:
        bin_reads = {}
        for label in todo:
            reads = sample_reads(files[label], a.reads)
            median_bp = sorted(len(d) for _k, d in reads)[len(reads) // 2] if reads else 0
            conc, expected = is_concatemer(label.split(u'-')[0], median_bp)
            if len(reads) < MIN_POPULATION:
                rows[label] = make_row(dict(bin=label, group=label.split(u'-')[0], reads=len(reads),
                                            note=u'too few reads (%d)' % len(reads)))
            elif conc:
                # 2026-09-07: two amplicons ligated end to end; not an organism, not named (see identity_verification)
                rows[label] = make_row(dict(bin=label, group=label.split(u'-')[0], reads=len(reads),
                                            note=u'concatemer: median read %d bp > %.1f x %d bp, two amplicons ligated, not named' % (median_bp, CONCATEMER_FACTOR, expected)))
            else:
                bin_reads[label] = reads
        if bin_reads:
            per_bin, windowed = batch_hits(bin_reads, root, a, tmp)
            results = {}
            for n, label in enumerate(sorted(bin_reads), 1):
                r = polish_bin(root, label, bin_reads[label], per_bin[label], windowed[label], a, tmp)
                results[label] = r
                print(u'  [%3d/%d] %-16s %s %.0f%% (2nd %s %.0f%%) | polish %s reads -> %s bp'
                      % (n, len(bin_reads), label, r.get('genus'), r.get('share', 0.0), r.get('genus2'),
                         r.get('share2', 0.0), r.get('pop_n', 0), r.get('bp', 0)))
                sys.stdout.flush()
            subs = {}
            for label, r in sorted(results.items()):
                for tag, source, pop in r.get('sub', []):
                    rs = polish_sub(root, tag, r['group'], pop, a, tmp)
                    rs['source'] = source
                    subs[tag] = rs
            batch_decide(results, root, a, tmp)
            if subs:
                batch_decide(subs, root, a, tmp)
                sub_dir = os.path.join(root, ALT_DIR)
                if not os.path.isdir(sub_dir):
                    os.makedirs(sub_dir)
                sub_tsv = os.path.join(root, 'PAK_SUB_POPULATIONS.tsv')
                old = {}
                if os.path.exists(sub_tsv):
                    for line in io.open(sub_tsv, encoding='utf-8', errors='replace').read().splitlines()[1:]:
                        if line.strip() and line.split(u'\t')[0].split(u'~')[0] not in results:
                            old[line.split(u'\t')[0]] = line
                for tag, rs in sorted(subs.items()):
                    with io.open(os.path.join(sub_dir, tag.replace(u'~', u'__') + '.fasta'), 'w', encoding='utf-8') as g:
                        g.write(u'>%s pak_sub %s reads=%d template=%s\n%s\n' % (tag, rs['source'], rs['pop_n'], rs['template'], rs['sequence']))
                    old[tag] = u'\t'.join([tag, tag.split(u'~')[0], rs['group'], rs['source'], u'%d' % rs['pop_n'], u'%d' % rs['bp'],
                                           rs.get('name') or u'-', u'%s' % rs.get('identity', u'-'), rs.get('decision_note') or u'-', rs.get('detail') or u'-'])
                    print(u'  SUB %-34s %-44s -> %s' % (tag[:34], rs['source'][:44], rs.get('name')))
                io.open(sub_tsv, 'w', encoding='utf-8', newline='\n').write(
                    u'\t'.join([u'sub tag', u'bin', u'group', u'source', u'reads', u'consensus bp', u'DECISION', u'identity %', u'decision note', u'locus detail'])
                    + u'\n' + u'\n'.join(old[k] for k in sorted(old)) + u'\n')
            for label, r in sorted(results.items()):
                if r.get('sequence'):
                    with io.open(os.path.join(set_dir, label + '.fasta'), 'w', encoding='utf-8') as g:
                        g.write(u'>%s pak_polish group=%s locus=%s population=%s share=%.0f reads=%d template=%s rounds=%d polish=%s\n%s\n'
                                % (label, r['group'], r['locus'], r.get('genus'), r.get('share', 0.0),
                                   r.get('pop_n', 0), r.get('template'), a.rounds, (r.get('polish') or u'samtools').replace(u' ', u'_'), r['sequence']))
                    levels[level(r.get('name'))] += 1
                    print(u'  %-16s %s -> %s' % (label, (r.get('detail') or u'')[:60], r.get('name')))
                else:
                    print(u'  %-16s %s' % (label, r.get('note')))
                rows[label] = make_row(r)
            flush()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    flush()
    print(u'')
    print(u'  levels (this run): %s' % dict(levels))
    print(u'  set: %s' % set_dir)
    print(u'  report: %s' % out_path)
    return 0


def self_test():
    u"""The fungal machinery's self-test plus the group logic and the 16S decision on a
    synthetic pair."""
    errors = []
    if not F.tools_present():
        print(u'  self-test skipped: blastn / makeblastdb / minimap2 / samtools not all present')
        return 0
    if F.self_test() != 0:
        errors.append(u'fungal machinery self-test failed')
    for g, want in ((u'F1', u'ITS'), (u'A2', u'SSU'), (u'B', u'SSU')):
        if group_locus(g) != want:
            errors.append(u'group_locus(%s) = %s' % (g, group_locus(g)))
    if window_db('B') != 'bacteria.16S.fna' or population_dbs('A1')[0] != 'archaea.16S.fna':
        errors.append(u'database mapping')
    # sub-populations on a synthetic split: two species groups inside the dominant genus + a second genus
    count = collections.Counter({u'Aus': 70, u'Bus': 30})
    members = {u'Aus': [(u'a%d' % i, 99.0, 500, u'Aus alpha' if i < 40 else u'Aus beta') for i in range(70)],
               u'Bus': [(u'b%d' % i, 99.0, 500, u'Bus beta') for i in range(30)]}
    reads = [(u'a%d' % i, u'ACGT' * 100) for i in range(70)] + [(u'b%d' % i, u'ACGT' * 100) for i in range(30)]
    subs = sub_populations(u'X-1_1', reads, count, members)
    kinds = sorted(x[1].split(u' (')[0] for x in subs)
    if len(subs) != 3 or kinds != [u'second genus population', u'species group inside the dominant genus', u'species group inside the dominant genus']:
        errors.append(u'sub-populations: %s' % [(x[0], x[1], len(x[2])) for x in subs])
    r = dict(bin=u'B-9_1', group=u'B', reads=3, note=u'too few reads (3)')
    if len(make_row(r).split(u'\t')) != len(HEADER):
        errors.append(u'make_row field count does not match the header')
    # the 16S decision on a synthetic pair through the ladder
    hits = [(99.9, 1450, u'NR_1 Aus alpha strain X 16S ribosomal RNA', 'db', 1500),
            (95.0, 1450, u'NR_2 Bus beta strain Y 16S ribosomal RNA', 'db', 1500)]
    res = basamaktan_sec({'16S': hits}, ['16S'], True)
    if not res or res[0] != u'Aus alpha':
        errors.append(u'16S ladder: %s' % (res,))
    for e in errors:
        print(u'  ERROR: %s' % e)
    print(u'  SELF-TEST (population polish): %s' % (u'PASSED' if not errors else u'%d ERROR(S)' % len(errors)))
    return 0 if not errors else 1


if __name__ == '__main__':
    sys.exit(main())
