# PrimerJury

<https://github.com/Yerlifan/primerjury>

**qPCR primer design from amplicon sequencing data, with identity verification that does not trust a single database.**

This pipeline was built for one stubborn problem: a classifier label is not an
identification. On real nanopore rDNA amplicon data, Kraken2's
lowest-common-ancestor calls and alignment-based identity **disagreed on a large
fraction of read bins**. The disagreement was not marginal. It reached genus and
phylum level:

| Kraken2 label | alignment-based identity |
|---|---|
| *Bacteroides ovatus* | *Porphyromonas* |
| *Ca. Cloacimonas acidaminovorans* | *Spirochaeta* |
| *Ca. Nitrosocosmicus hydrocola* | *Nitrososphaera* |
| *Colletotrichum higginsianum* | *Ramicandelaber* |

Design a species-specific primer on top of a wrong label and you get a primer
that works perfectly on the wrong organism. So this project verifies identity
independently before designing anything, then checks the resulting primers
against four independent layers before anything is ordered.

> **Status: research code, being generalised.** It runs and produces the
> published panel, but parts of it still assume the original study's targets.
> See [Current limitations](#current-limitations) before using it on your own
> data. Contributions welcome.

---

## What changed in September 2026

The study was re-run from the raw reads, and the run found things the tests
had not. Everything below is measured and recorded in `docs/WORK_RECORD.md`,
sections 13 to 16.

- **Bins without a classifier.** `./primerjury bins` turns demultiplexed BAMs
  into bins by length peak and minimap2 clustering; no Kraken2 needed, no label
  carried into the identity step.
- **The consensus template fault.** The dominant-allele step could not seed
  minimap2 on its own IUPAC-coded template: 0 of 3,001 reads aligned in a mixed
  bin, and a bin that looked healthy carried twenty wrong bases. Fixed; four
  bins rose a level, one from unnamed to species.
- **Consensus choice by the bin's own reads**, with a floor that names an
  untrustworthy choice (`./primerjury consensus`).
- **The exact-match rate.** Deliberate 3'-end variants had scored high under
  "at most one mismatch"; the read-level check reports mm0 and flags them
  (`./primerjury readcheck`). ARMS variants are off by default.
- **Evidence, not a number chosen for 16S.** A short reference matched over
  at least 90 per cent of its length counts in full; hits under 250 bases no
  longer veto a genus.
- **Stray files cannot become bins**, and a route that raised `NameError` the
  moment it was reached is fixed.

---

## What it does

```
sequences/  (your .fasta / .fastq, or demultiplexed BAMs)
     │
     ├─ 0. bin                 length peaks + minimap2 clustering, no classifier
     │                         (steps/from_raw.sh; replaces step 1 when there is no Kraken2)
     ├─ 1. classify            Kraken2 + Bracken           (optional, see below)
     ├─ 2. consensus           per-bin consensus, N-analysis
     ├─ 2b. fungal polish      split a mixed bin by ITS population, polish the dominant one
     │                         distinguishes low coverage from real strain variation
     ├─ 3. identity            12 reference databases, seed + full alignment,
     │                         NO taxonomy tree, ≥2 databases must agree
     ├─ 4. primer design       species / genus / functional-group / universal
     ├─ 5. specificity         4 independent layers (below)
     ├─ 5b. qiime2 route       profile + reference taxonomy + PICRUSt2 (optional, independent)
     └─ 6. report              ranked order list + evidence tables
```

### The four verification layers

A pair is only reported as safe to order when independent layers agree. They are
deliberately built on **different mechanisms**, so that a bug in one does not
silently confirm itself in another:

| Layer | Method | Independent of |
|---|---|---|
| 1 · in-sample | in-silico PCR against the raw reads | reference databases |
| 2 · local DB | scan of 12 local reference sets | network, NCBI |
| 3 · MFEprimer | external binary, thermodynamic amplicon search | our own code |
| 4 · NCBI | Primer-BLAST against `nt` | our database choices |

If the layers disagree, the pair is marked `CELISKILI` (contradictory) and is
**not** orderable. Disagreement is treated as information, not as noise.

### The exact-match rate is reported, not only "at most one mismatch"

Two oligos of the delivered panel disagreed with the template at the third
base from the 3' end in **100 per cent** of the target reads. Not a typo: the
candidate generator had produced deliberate -2/-3 variants (ARMS-style), the
panel criterion "at most one mismatch" could not see them, and the variant
scored artificially high. Nothing in the delivery said ARMS.

`verification/read_level_primer_check.py` therefore reports **mm0** (both
primers exact) next to mm1 for every pair and every bin, and flags a pair as
carrying a SYSTEMATIC MISMATCH when its members amplify at mm1 but the pooled
mm0 rate is zero. The deliberate-mismatch variants are off by default in the
generator.

### Identity verification

`verification/identity_verification.py` deliberately shares no mechanism with Kraken2:

- **No taxonomy tree, no k-mer LCA, no primers.** Seeds are extracted from the
  query consensus, the database is streamed, a short list is built, and every
  short-listed record is fully aligned (Levenshtein DP, infix).
- **Identity is measured twice**, over the whole overlap, and over the
  *discriminating window*: the columns where the best reference records differ
  from each other. Conserved regions (18S, 5.8S, LSU core) fall outside it, so a
  claim resting on false-high conserved-region identity becomes visible.
- **At least two independent databases must agree.** A single database's best
  hit is never an identification: deduplicated sets delete rare genera
  (measured: 0 *Petriella* records in SILVA LSURef NR99, 82 in the Parc set of
  the same release).
- **Unnamed environmental records cannot become a name.** A 99% match to
  `Uncultured bacterium clone 4B-11` is evidence that your sequence overlaps
  environmental clones, it is not a species. Reported as
  `CANNOT BE NAMED (reference is unnamed)`, never as a taxon.
- Output includes the **five nearest organisms**, deduplicated by organism
  rather than by record, so the list shows what else is close instead of the
  same species from five databases.

### Measuring what a bin really is

`steps/target_identity.py` puts the name a decision carries next to the organism
the sequence shows. Two rules govern it.

**Every class sees every rDNA database.** Each class used to be asked only its
own shortlist: archaea and bacteria of RefSeq 16S alone. That set holds 1,160 and
26,877 records while SILVA SSU NR99, in the same directory, holds 510,495. Worse,
limiting by domain assumes the label is right, which is the very thing this step
exists to test: a bin labelled bacterial that is really fungal can only be caught
by asking the fungal databases. `--databases narrow` restores the old, faster
coverage.

**Identities from different loci are never raced against each other.** 99 per
cent in 28S and 99 per cent in ITS are not the same measurement. Each locus is
decided with its own threshold, and the loci are walked in the ladder of the
class, the discriminating region first: 16S for archaea and bacteria, then ITS,
28S, 18S for fungi. When a later locus reaches species level and the one taken
did not, the row says so rather than looking more settled than the data is.

The thresholds come from one place, `verification/identity_verification.py`, and
follow the locus: 98.7 per cent for SSU (Kim et al. 2014), 99.6 for ITS and 99.8
for the fungal LSU (Vu et al. 2018). A species name is also refused when the
alignment is too short to carry it: 100 per cent over 484 bases and 100 per cent
over 2,900 are not the same evidence. A record matched over at least 250 bases that
covers at least 90 per cent of the record counts in full, because more than
half of the ITS references are shorter than the ITS floor (measured). Hits
under 250 bases are removed before ranking, so a short junk hit cannot veto
the genus a long one would give.

### Fungal bins: a mixed bin has no single consensus

Fungal amplicon bins are mixed. Measured on the study's data (2026-09-05): in
one bin 78 per cent of the reads went to *Petriella* by ITS and the rest to
other Microascaceae; the chain's blended consensus sat at 97.25 per cent to the
*P. musispora* type record while the bin's own single reads reached a median of
98.20 and a best of 100.00. A per-position majority over two organisms is a
sequence between them, and in ITS that is close to nothing; the bin stayed at
genus while its reads said species.

`./primerjury fungi` (`verification/fungal_bin_identity.py`) splits each fungal
bin read by read (best ITS record over RefSeq ITS and UNITE), keeps the dominant
genus population, picks that population's medoid read as the template (chosen by
k-mer similarity among the reads, not by any reference) and polishes it with
minimap2 and `samtools consensus` in two rounds. The polished sequence is decided
over the three loci with the same thresholds as everything else and lands as one
more candidate set: `select_consensus` weighs it with the same read-support
criterion, so it is used only where it represents the reads better. The per-read
table doubles as the read-level witness, from the same databases as the
consensus route. Self-test: `./primerjury fungi --self-test`.

Two things were measured on the way and changed the code. BLAST with a small
`-max_target_seqs` missed the best UNITE record in the Microascaceae pile-up (a
read written as *Microascus* 96.6 per cent had a real best record, *Acaulium*, at
98.3), so the per-read step maps the ITS windows with minimap2 against an index
built once beside the database: 401 windows in 9 s instead of 26 min, and every
genus that differed was a higher-identity record. And the top record is not the
last word: when it cannot name a species (unnamed, genus only, or a placeholder
such as "sp.") and a species-named record sits within the separation margin,
the decision rests on that record, in the multi-locus and the single-locus
paths alike (`tests/test_locus_tie.py`, `tests/test_single_locus_tie.py`).

Measured on the three roots of the study (2026-09-06): the polished set won the
read-support selection in 35 of 40, 31 of 40 and 140 of 173 fungal bins, and
the fungal bins at species level went from 0 to 18, 0 to 19 and 9 to 44, with
the independent re-derivation of every name agreeing 99 of 99 and the
consistency audit at zero on all three.

### The QIIME2 route

`./primerjury qiime2 all` runs an independent community profile beside the
identity chain, so that the two can be compared rather than one trusted:

1. `steps/qiime2_profile.sh`: per library, a length window from the sample's
   own read-length distribution (10th to 90th percentile; a 4.3 kb operon and a
   1.5 kb 16S library do not share a threshold), a random subsample with a fixed
   seed (the first reads of a nanopore file are the first minutes of the run,
   not the sample), de novo clustering, chimera removal, two feature thresholds
   (min 2 and min 10 reads, exported separately), phylogeny and diversity. Every
   step's count goes to `STEP_COUNTS.tsv`.
2. `steps/qiime2_classify.sh`: the ready-made classifier is trained for 16S and
   gave meaningless labels on fungi and on the long archaeal operon, so each
   library is aligned to the reference that fits it (SILVA SSU, SILVA LSU,
   UNITE) and `qiime2_reference_taxonomy.py` takes the taxonomy the best hits
   agree on, with the identity chain's alignment-length rule and thresholds.
3. `steps/picrust2_run.sh`: functional prediction for the archaeal and bacterial
   libraries, an operon cut down to its SSU part first, every step under a
   memory guard (the trait-table check at start-up is what eats memory, not the
   calculation; measured on a 16 GB machine).

### Starting from raw reads

The order is `./primerjury bins` (barcodes and bins), then
`bash steps/anchored_reference_consensus.sh`, then `./primerjury fungi` for the
fungal bins, then `./primerjury consensus` and `./primerjury run`.

`steps/from_raw.sh` turns demultiplexed BAMs (or one fastq per barcode) into
bins the chain reads, **without a classifier**:

```bash
bash steps/from_raw.sh --map examples/barcodes_example.tsv --bam-dir basecalled --parallel 2
```

1. `split_barcodes.py` reads the `BC:Z` tag dorado writes and produces one
   fastq per barcode.
2. `bin_reads.py` finds the length peaks of each barcode (a barcode rarely
   carries a single amplicon: measured shares of 6, 12 and 10 per cent for
   three products in one archaeal library), clusters a seed sample of every
   peak with minimap2 `ava-ont` at 97 per cent identity, and assigns every
   read to the nearest cluster centre at >= 90 per cent identity over >= 80
   per cent of the read. Reads that fit nowhere are counted, not dropped.
3. `select_bins.py` chooses the bins that enter the chain (share >= 0.5 per
   cent, or the five largest of a window, at most 40 per barcode) and records
   why every other bin was left out.

The bins are named `BIN<n>`; the consensus and identity steps accept them
beside Kraken taxids. Validated against the study's Kraken bins: the three
dominant methanogen bins were recovered at 90 per cent recall and 100 per cent
purity. The bin definition carries no database label, so the identity step
is the first place a name appears.

Two faults in the consensus steps were found on this route and fixed
(`docs/WORK_RECORD.md`, section 13): the dominant-allele step could not seed
minimap2 on its IUPAC-coded template (0 of 3,001 reads aligned in a mixed
bin; a healthy-looking bin carried twenty wrong bases), and files of another
library left in a bin folder became phantom bins. Both had passed every test.

---

## Requirements

| Tool | Used for |
|---|---|
| WSL2 or Linux | everything (the pipeline is POSIX; Windows is supported through WSL) |
| Python ≥ 3.8 | pipeline and analysis |
| blastn / makeblastdb | specificity scans |
| MFEprimer 4.4 | thermodynamics, off-target amplicons (layer 3) |
| minimap2, samtools | read alignment, consensus |
| seqkit | sequence handling |
| Kraken2 + Bracken | classification (optional, see note) |
| QIIME2 + PICRUSt2 | community/function analysis (optional) |

Python packages are in `requirements.txt`.

> **Kraken2 is optional and deliberately so.** The pipeline runs without it; the
> classification stages are marked skipped and the chain continues. Given the
> label disagreements above, Kraken2 is treated as one opinion to be checked,
> not as ground truth.

---

## Quick start

Linux or WSL2. Everything runs through one file:

```bash
./primerjury                   # the built-in guide
./primerjury check             # what is installed? changes nothing
./primerjury install all       # tools + reference databases + QIIME2
./primerjury bins --map barcodes.tsv --bam-dir basecalled   # from raw reads, no classifier
./primerjury run               # the full chain
```

Reference databases total roughly **28 GB** and are never committed to git.

### Downloads are verified, not trusted

Every download is measured after it lands: is it really FASTA, how many records,
RNA or DNA alphabet. If the measurement fails the file is renamed `.SUPHELI`
(suspect) and **is not used**. Silently running against a truncated database is
worse than not running at all.

This also protects against the real failure mode of hard-coded URLs: SILVA and
UNITE rename files every release, so a stale URL would otherwise fetch nothing, or the wrong release, without any complaint. UNITE's URL is not hard-coded at all
because it changes per release DOI; the installer tells you where to get the
current one and verifies whatever you hand it.

### Choosing the Kraken2 k-mer length

Prebuilt Kraken2 databases are fixed at `k=35, l=31`. To choose your own:

```bash
./primerjury install kraken --kmer 31
```

Shorter *k* raises sensitivity on error-prone long reads (ONT) but pushes the
LCA up the tree; longer *k* is more specific but loses hits to single errors.
The build parameters are written to `$DB/BUILD_INFO.txt`, because
`opts.k2d` is binary and six months later nobody remembers which *k* was used.

**k-mer choice alone does not fix identification.** The disagreements at the top
of this README come from the LCA step and database coverage, not from *k*.
Whatever *k* you build with, verify identity independently.

---

## Run

Put your sequences in `sequences/`, then:

```bash
./primerjury run
```

The full chain runs ten stages in dependency order and **checks the output of
each one**. A zero exit code is not accepted as success: the expected file must
exist, be non-empty, and in several stages its contents are inspected. A failed
check stops the chain rather than quietly continuing, the recurring failure
mode in this domain is a program that produces a wrong or empty answer without
erroring.

Individual stages:

```bash
./primerjury identity          # identity verification
./primerjury specificity       # four-layer specificity
./primerjury panel             # single-protocol panel measurement
./primerjury audit             # independent read-only audit
./primerjury status            # current panel state
./primerjury check             # what is installed; changes nothing
```

### Tests

```bash
./primerjury test      # every suite: 157 regression tests plus the four
                       # test files. Exits non-zero if anything fails.
./primerjury doctor    # the repository checks itself
```

Or individually:

```bash
python3 steps/regression_test.py       # 157 tests, expectations derived from
                                       # the design decisions or known maths
python3 tests/test_repo_health.py      # imports, names, format strings, line
                                       # endings, packaging, option handoffs
python3 tests/test_taxonomy.py         # 5 header formats, real DBs
python3 tests/test_unnamed_records.py  # unnamed records cannot become names
python3 tests/test_orientation_trap.py # orientation-trap detector
```

Each exits non-zero on failure and prints what it measured.

Several groups exist to hold down a **contract between two files** that nothing
in the language enforces: one script prints a counter and four others read it
back with a regular expression, one writes a verdict token into a TSV and the
workbook builder looks it up. Rewording either side breaks the link in silence,
which has happened, so a test now fails instead. Each of those tests was
confirmed by breaking the contract on purpose and watching it fail.

---

## Design rules

These are not style preferences; each was paid for with a real bug.

1. **No decision rests on a single code path.** Every measurement is taken two
   independent ways. If they diverge, the candidate is rejected and the
   divergence is recorded.
2. **A zero exit code is not success.** Output is checked for existence,
   emptiness, and content.
3. **"Unknown" is a distinct state from "clean".** An unmeasured layer never
   votes in favour. `BILINMIYOR` is never folded into `TEMIZ`.
4. **Expensive measurement and cheap judgment are cached separately.** Scans are
   checkpointed; verdicts are re-derived every run. Caching a verdict means a
   fix to the judgment logic silently does nothing, this happened, and was
   caught only by diffing before/after outputs.
5. **Long runs are resumable.** Every stage records where it stopped.
6. **Nothing is skipped silently.** Anything that could not be installed, read,
   or measured is listed at the end and changes the exit code.

---

## Current limitations

Honest list; these are the gaps between "runs for the original study" and
"general-purpose tool":

- **Targets are still study-specific.** `steps/targets.tsv` and
  `screening/target_clades.tsv` describe the original 20 targets and 5
  amplicon groups. Samples in `examples/` show the format. Generalising the
  target definition is the main open work item.
- **`screening/config.py` holds every path and constant** and is
  meant to be edited. Moving it to YAML/TOML is planned, the good news is that
  it is genuinely the only place paths are defined.
- **Layer 2 taxonomic discrimination is new** and its effect on verdicts is
  still being measured; the size-based criterion remains the one that votes.
- **Some internal names are Turkish, on purpose.** Everything a reader meets is
  English: screen messages, help text, comments, file and directory names, and
  command line options. What stays Turkish is the layer underneath, and only
  because it is DATA rather than wording: TSV column names, verdict values such
  as `DOGRULANDI` and `CELISKILI`, and the taxonomic level tokens. Those are
  compared as strings across the pipeline, written into files later stages read
  back, and stored inside checkpoints from earlier runs. Renaming them would be
  a schema migration disguised as a translation, and it would fail as an empty
  table rather than as an error. `screening/labels.py` supplies the English
  wording wherever a person reads one.

---

## Repository layout

| Path | Contents |
|---|---|
| `screening/` | search engine, in-silico PCR, scoring, configuration |
| `verification/` | orchestration, identity verification, the four layers |
| `steps/` | the pipeline steps, run in dependency order |
| `protocol/` | single-protocol panel measurement |
| `scoring/` | shared scoring |
| `tests/` | tests |
| `tools/` | Kraken2 environment/database tooling |
| `primerjury` | **the single entry point, start here** |
| `docs/` | user guide, audit report, measurements |
| `sequences/` | **your input goes here** |

---

## Documentation

**[Full user guide → `docs/GUIDE.md`](docs/GUIDE.md)**, installation, input
preparation, defining your own targets, reading the output, and troubleshooting.

`docs/AUDIT_2026-08-21.md` is the pre-release code audit: what was measured,
what was broken and what was fixed. `docs/WORK_RECORD.md` is the record of the
design decisions and the faults found while the pipeline was built.

## Licence

**[PolyForm Noncommercial 1.0.0](LICENSE)**, free for any noncommercial purpose.

- **Anyone may use it**: researchers, students, universities, public research
  organisations, health and environmental organisations, government institutions,
  hobbyists, regardless of how they are funded.
- **You may modify and redistribute it**, keeping this licence.
- **You may not make money from it**: no commercial products, paid services, or
  commercial advantage built on this software.

Note that this is a source-available licence, not an OSI-approved open-source
one; GitHub will show it as "Other". That is the intended trade-off.

## Citation

Author: **Burak Aslancan Pak** ([ORCID 0000-0002-7793-2215](https://orcid.org/0000-0002-7793-2215))

If this is useful in published work, please cite the repository, see
[CITATION.cff](CITATION.cff).

## Contributing

Issues and pull requests are welcome, in English or Turkish. If you change
anything that produces a number, please include the before/after measurement: that is the standard the rest of the codebase is held to.
