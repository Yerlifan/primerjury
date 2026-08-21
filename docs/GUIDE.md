# PrimerJury — User Guide

A step-by-step guide for running the pipeline on your own sequencing data.

**Contents**

1. [What this tool is for](#1-what-this-tool-is-for)
2. [Installing](#2-installing)
3. [Preparing your input](#3-preparing-your-input)
4. [Defining your targets](#4-defining-your-targets)
5. [Running the pipeline](#5-running-the-pipeline)
6. [Reading the output](#6-reading-the-output)
7. [Understanding the verdicts](#7-understanding-the-verdicts)
8. [Repository layout](#8-repository-layout)
9. [Troubleshooting](#9-troubleshooting)
10. [Extending the tool](#10-extending-the-tool)

---

## 1. What this tool is for

PrimerJury designs qPCR primers from amplicon sequencing data and — more
importantly — **refuses to trust any single source when doing it**.

The name is the design: every claim faces a jury of independent methods, and a
split jury blocks the decision rather than being averaged away.

Two places where this matters:

**Identity.** Before designing a species-specific primer, you must know what the
sequence actually is. A classifier label is one opinion. On the dataset this was
built for, Kraken2's lowest-common-ancestor calls and alignment-based identity
disagreed at genus and phylum level on a substantial share of read bins. Design a
"species-specific" primer on a wrong label and you get a primer that works
perfectly, for the wrong organism.

**Specificity.** A primer pair that looks clean against one database may amplify
something else entirely. PrimerJury checks four layers built on deliberately
different mechanisms, so a bug in one cannot confirm itself in another.

### When this tool is a good fit

- You have amplicon sequencing data (16S/18S/ITS/28S; nanopore or Illumina).
- You want qPCR primers for specific taxa, functional groups, or universal
  controls.
- Getting it wrong is expensive — you are ordering oligos, or publishing.

### When it is not

- You want a quick primer for a single well-characterised isolate. Use Primer3
  directly; this pipeline's machinery would be wasted.
- You need a GUI. This is a command-line tool.

---

## 2. Installing

### Requirements

You need **WSL2 (Windows) or Linux**, and Python 3.8+.

| Tool | Purpose | Required? |
|---|---|---|
| blastn, makeblastdb | specificity scans | yes |
| MFEprimer 4.4 | thermodynamics, off-target amplicons | yes |
| minimap2, samtools | read alignment, consensus | yes |
| seqkit | sequence handling | yes |
| Kraken2, Bracken | classification | optional |
| QIIME2, PICRUSt2 | community and function analysis | optional |

### One command

```bash
bash install.sh durum          # measure what you already have; changes nothing
bash install.sh araclar        # install the tools
bash install.sh veritabani     # download reference databases (~28 GB)
bash install.sh qiime          # QIIME2 + PICRUSt2 (optional)
bash install.sh hepsi          # everything
```

> Subcommands are still Turkish words; they are being translated. `durum` =
> status, `araclar` = tools, `veritabani` = databases, `hepsi` = all.

### Downloads are verified, not trusted

After every download the file is measured: is it really FASTA, how many records,
RNA or DNA alphabet. If the check fails the file is renamed `.SUPHELI` (suspect)
and is **not used**.

This matters because SILVA and UNITE rename their files every release. A
hard-coded URL silently goes stale and fetches nothing — or the wrong release.
UNITE's URL is not hard-coded at all; the installer points you to the current
release page and verifies whatever you give it:

```bash
bash install.sh veritabani --yalniz unite --url <URL-you-copied>
```

### Building indexes

MFEprimer needs an index per database:

```bash
bash build_index.sh --liste          # show candidate files
bash build_index.sh SILVA_138.2_SSURef_NR99.fasta
```

**Do not skip this for SILVA.** SILVA stores sequences as RNA (U instead of T).
MFEprimer builds its 9-mer index over {A,C,G,T}; with no T in the file the
effective alphabet collapses to {A,C,G} and only 3⁹ = 19,683 k-mers enter the
index instead of 4⁹ = 262,144. **MFEprimer reports success anyway.**
`build_index.sh` converts U→T on sequence lines only (headers contain capital U
in words like "Unknown", and a blind `sed 's/U/T/g'` would corrupt the taxonomy).

### Choosing a Kraken2 k-mer length

Prebuilt Kraken2 databases are fixed at `k=35, l=31`. To choose your own:

```bash
bash install.sh kraken-kur --kmer 31 --db ~/k2db_k31
```

Shorter *k* raises sensitivity on error-prone long reads (nanopore) but pushes
the LCA up the tree — you get genus where you wanted species. Longer *k* is more
specific but loses hits to single sequencing errors.

Build parameters are written to `$DB/KURULUM_BILGISI.txt`, because `opts.k2d` is
binary and six months later nobody remembers which *k* was used.

**k-mer choice alone does not fix identification.** The label disagreements come
from the LCA step and database coverage, not from *k*. Whatever *k* you build
with, verify identity independently.

---

## 3. Preparing your input

Put your sequences in `sequences/`:

```
sequences/
├── sample-A1/
│   └── barcode01.fastq.gz
└── sample-B2/
    └── barcode02.fastq.gz
```

Flat files work too. Both `.fasta` and `.fastq` are accepted, gzipped or not.

**Base-name collisions are refused, not resolved.** `barcode03.fastq` and
`barcode03.fastq.gz` would map to the same output name. In the original study
exactly this caused a silent double-count, so the pipeline now stops instead of
picking one.

### A trap worth knowing about

Consensus sequences must be **orientation-normalised** before in-silico PCR. In
the source dataset the raw consensus folder was mixed: 71 antisense, 27 sense.
A reverse-oriented consensus silently yields **zero products** — measured loss
100%. The candidate then looks like it fails, when in fact it was never tested.

The pipeline produces a canonical, single-orientation set and reads only from
there. If you add consensus sequences by hand, normalise them first.

---

## 4. Defining your targets

This is currently the least generalised part of the tool. Two tables drive it:

### `steps/hedefler.tsv` — what you want to amplify

| column | meaning |
|---|---|
| `karar` | decision group: 1 species, 2 genus, 3 functional group, 4 universal |
| `hedef` | target name (used as an identifier everywhere downstream) |
| `duzey` | level the primer must discriminate at: `tur` (species), `cins` (genus) |
| `in_taxid` | taxids of members — the product **must** form in all of them |
| `haric` | taxids deliberately excluded from the competitor set |
| `hedef_tur` | human-readable organism name |

### `screening/hedef_klad.tsv` — how to tell inside from outside

| column | meaning |
|---|---|
| `hedef` | must match `hedefler.tsv` |
| `alan` | domain: `Bacteria`, `Archaea`, `Eukaryota` |
| `klad` | comma-separated clade tokens that mark "inside the target" |

This second table is what lets the pipeline say *"this off-target hit is actually
a length variant of the target itself"* rather than counting it as cross-reaction.
Without it, group and universal primers look far worse than they are — measured:
of 1,605 amplicons flagged "off-target" by size alone, **95.7% were inside the
target clade**.

Working examples for both files are in `examples/`.

### The four target levels

| Level | `karar` | What it means | Hardest part |
|---|---|---|---|
| Species | 1 | discriminates one species from its siblings | sibling species differ by very few bases |
| Genus | 2 | amplifies a genus, excludes neighbours | genus boundaries move in the literature |
| Functional group | 3 | a metabolic guild, possibly polyphyletic | members may share no single conserved region |
| Universal | 4 | a whole domain, used as a control | coverage matters more than specificity |

Universal primers are scored differently: coverage is the criterion, not absence
of off-target products.

---

## 5. Running the pipeline

### The full chain

```bash
python3 verification/tam_zincir.py --kok . --onayla
```

Ten stages run in dependency order. Order is not arbitrary:

| Stage | What it does | Typical time |
|---|---|---|
| self-test | code verifies itself before any measurement | 1–2 min |
| quick consistency | runs the whole chain on a small subset with a known answer | 25–40 min |
| Kraken environment | reports what is installed and which database | 1–5 min |
| identity verification | tests every identity claim independently | 3–4 h |
| all bin identities | verifies every bin entering the panel | 4–6 h |
| full measurement | panel measurement, recovery, verification | 30–90 min |
| threshold scan | Kraken confidence threshold sweep | 1–2 h |
| second database | re-run with a different Kraken database | 2–8 h |
| comparison table | four methods side by side | 1–2 min |
| summary | single combined report | 1 min |

Kraken stages are **skipped, not failed**, if Kraken2 is absent. The reason is
recorded and the chain continues.

### Individual stages

```bash
python3 verification/kimlik_dogrulama.py --kok .      # identity verification
python3 verification/dogrulama_turu.py --kok .        # four-layer specificity
python3 protocol/tek_protokol_olc.py --kok .          # single-protocol panel measurement
python3 cross_check.py --kok .                        # read-only independent audit
```

### Two rules the runner enforces

**A zero exit code is not success.** After every stage the expected output must
exist, be non-empty, and in several stages its content is inspected. The
recurring failure mode in this domain is a program that produces a wrong or empty
answer *without erroring*.

**Everything is resumable.** Each stage records where it stopped. Re-running the
same command skips completed stages. Long scans checkpoint per chunk.

### Tests

```bash
python3 tests/test_repo_health.py        # the repository itself: imports, references, packaging
python3 tests/test_taxonomy.py           # 5 FASTA header formats, against real databases
python3 tests/test_unnamed_records.py    # unnamed records cannot become species names
python3 tests/test_orientation_trap.py   # the mixed-orientation detector
```

Each exits non-zero on failure and prints what it measured.

`test_repo_health.py` deserves a note. It checks what a syntax check cannot:
whether the packages actually *import* (side effects included), whether every
script referenced by another script exists, whether any pre-rename or personal
name survived, and whether what git ships is enough to run. It exists because
`screening/motor.py` once passed every parse check and still could not be
imported — it loads its engine from disk, and those files had not been copied.
"110 files, 0 errors" was true and useless.

---

## 6. Reading the output

| File | Contents |
|---|---|
| `KIMLIK_SONUC/kimlik_iddialari.tsv` | identity claims, verdicts, five nearest organisms |
| `DOGRULAMA_SONUC/dogrulama_uc_sutun.tsv` | the four layers side by side, per pair |
| `DOGRULAMA_SONUC/CELISKILER.md` | contradictions — the most valuable output |
| `TEK_PROTOKOL_SONUC/SIPARIS_LISTESI.tsv` | the order list |
| `GUNCEL_DURUM.md` | current panel state, regenerated every run |

### The five nearest organisms

Identity output lists the five closest organisms, **deduplicated by organism
rather than by record** — otherwise the same species from five databases fills
the list and tells you nothing.

Three kinds of entry appear, and the distinction is deliberate:

```
1) adsiz: Uncultured prokaryote clone A13 16S...  %99,00  [NCBI nt]
2) Methanosarcina vacuolata                       %97,57  [RefSeq archaea 16S]
3) uncultured bacterium (Petrimonas)              %96,79  [SILVA SSU NR99]
```

- A **name** — the record is identified to species or genus.
- **`adsiz:`** — the record carries no taxonomy at all (an environmental clone).
  A 99% match here is *not* an identification.
- **`deepest (parent)`** — the record has full taxonomy but no species binomial.
  Classified, just not to species level. This is not the same as unnamed.

---

## 7. Understanding the verdicts

### Identity verdicts

| Verdict | Meaning |
|---|---|
| `DOGRULANDI` | ≥2 independent databases support the claim |
| `DUZELTILMELI` | ≥2 databases agree on a *different* answer; the correct wording is given |
| `DOGRULANAMADI` | evidence insufficient, contradictory, or single-source |

A single database's best hit is never an identification. Deduplicated sets delete
rare genera — measured: 0 *Petriella* records in SILVA LSURef NR99, 82 in the
Parc set of the same release. A verdict resting on one source would inherit that
gap as fact.

### Specificity verdicts

| Verdict | Meaning |
|---|---|
| `KESIN` | three layers measured and all clean |
| `KOSULLU` | two layers measured and agreeing |
| `CELISKILI` | measuring layers disagree — **not orderable** |
| `RISKLI` | off-target products found |
| `EKSIK` | too few layers measured to judge |

### Three states, not two

`BILINMIYOR` (unknown) is a distinct state from `TEMIZ` (clean). A layer that did
not run never votes in favour. This sounds obvious and is the single most common
way these pipelines lie: an unmeasured check silently reads as a passed check.

---

## 8. Repository layout

```
primerjury/
├── install.sh              one entry point: tools, databases, Kraken2, QIIME2
├── build_index.sh          MFEprimer index builder (handles the SILVA U→T trap)
├── cross_check.py          independent read-only audit of a finished run
├── sequences/              YOUR INPUT GOES HERE
├── examples/               example target and clade tables
├── screening/              search engine, in-silico PCR, scoring, configuration
├── verification/           orchestration, identity verification, the four layers
├── steps/                  the numbered pipeline, in dependency order
├── protocol/               single-protocol panel measurement
├── scoring/                shared scoring
├── engine/                 low-level sequence engine (ispcr, reads, scanner)
├── tools/                  Kraken2 environment and database tooling
├── tests/                  tests
└── docs/                   this guide, plus the pre-release audit
```

### Pipeline order

File names no longer carry numbers — a number in a filename breaks the moment
you insert a step. The order lives here instead:

| # | Script (`steps/`) | Purpose |
|---|---|---|
| 1 | `check_environment.sh` | cores, RAM, tool versions, databases found |
| 2 | `reclassify_kraken2.sh` | Kraken2 reclassification (memory-mapping aware) |
| 3 | `analyze_ambiguous_bases.sh` | separates low coverage from real strain variation |
| 4 | `generate_primer_candidates.py` | candidate oligos from a consensus |
| 5 | `design_group_primers.py` | primers for a group of members |
| 6 | `split_clusters.py` | splits heterogeneous bins |
| 7 | `anchored_reference_consensus.sh` | reference-anchored consensus |
| 8 | `freeze_reference.sh` | pins the reference set |
| 9 | `batch_design.py` | batch design across targets |
| 10 | `specificity.py` | specificity against in-sample competitors |
| 11 | `indistinguishable_targets.py` | targets that cannot be separated |
| 12 | `check_bin_identity.py` | bin identity check |
| 13 | `dominant_allele_consensus.py` | dominant-allele consensus |
| 14 | `export_excel.py` | Excel deliverables |
| 15 | `external_databases.py` | external database scan |
| 16 | `design_from_reference.py` | design from reference sequences |
| 17 | `check_primer_geometry.py` | geometry gate |
| 18 | `regression_test.py` | regression suite |
| 19 | `check_deliverables.py` | deliverable audit |
| 20 | `mfeprimer_layer.py` | MFEprimer layer |
| 21 | `recover_bins.py` | recovers bins below threshold |
| 22 | `community_trends.py` | community trends |
| 23 | `target_identity.py` | target identity |
| 24 | `reassign_confidence.py` | confidence reassignment |
| 25 | `abundance_rank.py` | abundance-rank curves |
| 26 | `reference_identity.py` | reference identity |
| 27 | `check_taxonomic_level.py` | does it discriminate at the required level? |

---

## 9. Troubleshooting

**"MFEprimer index built successfully" but no off-target hits anywhere.**
The index is probably built over a 3-letter alphabet. Re-run
`bash build_index.sh <file>`; it measures RNA vs DNA and reports what it did.

**A candidate shows zero products in every member.**
Check consensus orientation first. A reverse-oriented consensus gives zero
products silently.

**The chain stops at the quick consistency test.**
That gate compares against reference values recorded with the primer sequences
they were measured from. If you changed a pair, the reference is stale, and the
run reports "reference invalid, pair changed" instead of failing the chain.

**Everything is `BILINMIYOR`.**
Usually missing databases. Run `bash install.sh durum` — it reports how many of
the independent sources are present and what that does to the verdicts.

**A rerun reports "taken from previous run" and nothing changed after I fixed
something.**
Expensive scans are cached; verdicts are re-derived every run. If you changed
scan *parameters*, the signature changes and the scan re-runs. If a checkpoint
predates this separation the run says so explicitly and rescans.

**The computer locks up during a long run.**
WSL memory. See `install.sh durum` output and lower thread counts; the Kraken
tooling picks `nproc - 2` by default.

---

## 10. Extending the tool

### Known gaps

Being honest about what is not yet general:

- **Target definitions are still study-shaped.** `hedefler.tsv` and
  `hedef_klad.tsv` work, but assume a structure close to the original study
  (amplicon groups, bin naming). Generalising this is the main open work.
- **`screening/yapilandirma.py` holds every path and constant.** It is meant to
  be edited and is genuinely the only place paths are defined — moving it to
  YAML/TOML is planned and should be straightforward.
- **Layer 2 taxonomic discrimination is new.** It is measured and reported
  alongside the size-based criterion, but the size-based one still casts the
  vote. Switching it changes verdicts and needs a before/after study.
- **Code comments and identifiers are Turkish.** Interface, docs and file names
  are English. Translating 57k lines of identifiers is planned but has to be
  staged carefully: several identifiers are also TSV column names and checkpoint
  keys, so renaming them changes output schemas.

### The design rules

If you contribute, these are the rules the codebase is held to. Each was paid for
with a real bug:

1. **No decision rests on a single code path.** Every measurement is taken two
   independent ways; divergence rejects the candidate and is recorded.
2. **A zero exit code is not success.** Check the output, not the return value.
3. **"Unknown" is distinct from "clean".** An unmeasured layer never votes yes.
4. **Cache expensive measurement, re-derive cheap judgment.** Caching a verdict
   means a fix to judgment logic silently does nothing — this happened, and was
   caught only by diffing before/after output.
5. **Long runs are resumable.**
6. **Nothing is skipped silently.** Anything not installed, read, or measured is
   listed at the end and changes the exit code.
7. **If you change something that produces a number, show the before/after.**

### Reporting problems

Issues and pull requests are welcome, in English or Turkish. For a suspected
wrong result, please include the relevant `*_SONUC/` tables — the verdicts carry
their own evidence columns, which usually make the cause visible.
