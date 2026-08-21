# MicRhoBooster

**qPCR primer design for anaerobic digester microbiomes — with independent identity verification that does not trust a single database.**

This pipeline was built for one stubborn problem: a Kraken2 label is not an
identification. On real nanopore rDNA amplicon data from an anaerobic digester,
Kraken2's lowest-common-ancestor calls and alignment-based identity **disagreed
on a large fraction of read bins** — not marginally, but at genus and phylum level:

| Kraken2 label | alignment-based identity |
|---|---|
| *Bacteroides ovatus* | *Porphyromonas* |
| *Ca. Cloacimonas acidaminovorans* | *Spirochaeta* |
| *Ca. Nitrosocosmicus hydrocola* | *Nitrososphaera* |
| *Colletotrichum higginsianum* | *Ramicandelaber* |

Designing species-specific primers on top of a wrong label produces a primer
that works perfectly — for the wrong organism. So this project verifies identity
independently before designing anything, and verifies the resulting primers
against four independent layers before anything is ordered.

> **Status: research code, being generalised.** It runs and produces the
> published panel, but parts of it still assume the original study's targets.
> See [Current limitations](#current-limitations) before using it on your own
> data. Contributions welcome.

---

## What it does

```
sekanslar/  (your .fasta / .fastq)
     │
     ├─ 1. classify            Kraken2 + Bracken           (optional, see below)
     ├─ 2. consensus           per-bin consensus, N-analysis
     │                         distinguishes low coverage from real strain variation
     ├─ 3. identity            12 reference databases, seed + full alignment,
     │                         NO taxonomy tree, ≥2 databases must agree
     ├─ 4. primer design       species / genus / functional-group / universal
     ├─ 5. specificity         4 independent layers (below)
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

### Identity verification

`KURTARMA/kimlik_dogrulama.py` deliberately shares no mechanism with Kraken2:

- **No taxonomy tree, no k-mer LCA, no primers.** Seeds are extracted from the
  query consensus, the database is streamed, a short list is built, and every
  short-listed record is fully aligned (Levenshtein DP, infix).
- **Identity is measured twice** — over the whole overlap, and over the
  *discriminating window*: the columns where the best reference records differ
  from each other. Conserved regions (18S, 5.8S, LSU core) fall outside it, so a
  claim resting on false-high conserved-region identity becomes visible.
- **At least two independent databases must agree.** A single database's best
  hit is never an identification: deduplicated sets delete rare genera
  (measured: 0 *Petriella* records in SILVA LSURef NR99, 82 in the Parc set of
  the same release).
- **Unnamed environmental records cannot become a name.** A 99% match to
  `Uncultured bacterium clone 4B-11` is evidence that your sequence overlaps
  environmental clones — it is not a species. Reported as
  `ADLANDIRILAMIYOR (referans adsız)`, never as a taxon.
- Output includes the **five nearest organisms**, deduplicated by organism
  rather than by record, so the list shows what else is close instead of the
  same species from five databases.

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
| Kraken2 + Bracken | classification (optional — see note) |
| QIIME2 + PICRUSt2 | community/function analysis (optional) |

Python packages are in `requirements.txt`.

> **Kraken2 is optional and deliberately so.** The pipeline runs without it; the
> classification stages are marked skipped and the chain continues. Given the
> label disagreements above, Kraken2 is treated as one opinion to be checked,
> not as ground truth.

---

## Install

One entry point handles tools, databases, Kraken2 and QIIME2:

```bash
bash KUR.sh durum          # measure what is installed — changes nothing
bash KUR.sh araclar        # kraken2, bracken, minimap2, samtools, blast, seqkit
bash KUR.sh veritabani     # SILVA + PR2 + RefSeq (+ UNITE/ROD instructions)
bash KUR.sh qiime          # QIIME2 + PICRUSt2 + SILVA classifier
bash KUR.sh hepsi          # all of the above
```

Reference databases total roughly **28 GB** and are never committed to git.

### Downloads are verified, not trusted

Every download is measured after it lands: is it really FASTA, how many records,
RNA or DNA alphabet. If the measurement fails the file is renamed `.SUPHELI`
(suspect) and **is not used**. Silently running against a truncated database is
worse than not running at all.

This also protects against the real failure mode of hard-coded URLs: SILVA and
UNITE rename files every release, so a stale URL would otherwise fetch nothing —
or the wrong release — without complaint. UNITE's URL is not hard-coded at all
because it changes per release DOI; the installer tells you where to get the
current one and verifies whatever you hand it.

### Choosing the Kraken2 k-mer length

Prebuilt Kraken2 databases are fixed at `k=35, l=31`. To choose your own:

```bash
bash KUR.sh kraken-kur --kmer 31 --db ~/k2db_k31
```

Shorter *k* raises sensitivity on error-prone long reads (ONT) but pushes the
LCA up the tree; longer *k* is more specific but loses hits to single errors.
The build parameters are written to `$DB/KURULUM_BILGISI.txt`, because
`opts.k2d` is binary and six months later nobody remembers which *k* was used.

**k-mer choice alone does not fix identification.** The disagreements at the top
of this README come from the LCA step and database coverage, not from *k*.
Whatever *k* you build with, verify identity independently.

---

## Run

Put your sequences in `sekanslar/`, then:

```bash
bash KAPSAMLI_ARAMA.bat        # Windows menu (legacy)
python3 KURTARMA/tam_zincir.py --kok . --onayla    # full chain
```

The full chain runs ten stages in dependency order and **checks the output of
each one**. A zero exit code is not accepted as success: the expected file must
exist, be non-empty, and in several stages its contents are inspected. A failed
check stops the chain rather than quietly continuing — the recurring failure
mode in this domain is a program that produces a wrong or empty answer without
erroring.

Individual stages:

```bash
python3 KURTARMA/kimlik_dogrulama.py --kok .      # identity verification
python3 KURTARMA/dogrulama_turu.py --kok .        # 4-layer specificity
python3 TEK_PROTOKOL/tek_protokol_olc.py --kok .  # single-protocol panel measurement
python3 capraz_kontrol.py --kok .                 # read-only independent audit
```

### Tests

```bash
python3 DENETIM_SINAMALARI/sina_taksonomi.py         # 5 header formats, real DBs
python3 DENETIM_SINAMALARI/sina_adsiz_kayit.py       # unnamed records cannot become names
python3 DENETIM_SINAMALARI/sina_D9_karisik_klasor.py # orientation-trap detector
```

Each exits non-zero on failure and prints what it measured.

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
   fix to the judgment logic silently does nothing — this happened, and was
   caught only by diffing before/after outputs.
5. **Long runs are resumable.** Every stage records where it stopped.
6. **Nothing is skipped silently.** Anything that could not be installed, read,
   or measured is listed at the end and changes the exit code.

---

## Current limitations

Honest list; these are the gaps between "runs for the original study" and
"general-purpose tool":

- **Targets are still study-specific.** `WSL_betikleri/hedefler.tsv` and
  `KAPSAMLI_ARAMA/hedef_klad.tsv` describe the original 20 targets and 5
  amplicon groups. Samples in `ornek_veri/` show the format. Generalising the
  target definition is the main open work item.
- **`KAPSAMLI_ARAMA/yapilandirma.py` holds every path and constant** and is
  meant to be edited. Moving it to YAML/TOML is planned — the good news is that
  it is genuinely the only place paths are defined.
- **Layer 2 taxonomic discrimination is new** and its effect on verdicts is
  still being measured; the size-based criterion remains the one that votes.
- **The Windows `.bat` menu is legacy** (2000 lines) and will be replaced by a
  proper CLI.
- Code comments and internal reports are in Turkish. Interface, docs and issues
  are in English. Renaming 57k lines of identifiers was judged too risky to be
  worth it.

---

## Repository layout

| Path | Contents |
|---|---|
| `KAPSAMLI_ARAMA/` | search engine, in-silico PCR, scoring, configuration |
| `KURTARMA/` | orchestration, identity verification, the four layers |
| `WSL_betikleri/` | numbered pipeline, `00` → `27` |
| `TEK_PROTOKOL/` | single-protocol panel measurement |
| `ORTAK_PUANLAYICI/` | shared scoring |
| `DENETIM_SINAMALARI/` | tests |
| `ARACLAR/` | Kraken2 environment/database tooling |
| `docs/` | audit report and working log (Turkish) |
| `sekanslar/` | **your input goes here** |

---

## Citation and licence

Released under the MIT licence — see [LICENSE](LICENSE). If you need a different
licence, change that file before publishing; nothing else depends on it.

If this is useful in published work, please cite the repository. Add author and
DOI details to `CITATION.cff` before release.

## Contributing

Issues and pull requests are welcome, in English or Turkish. If you change
anything that produces a number, please include the before/after measurement —
that is the standard the rest of the codebase is held to.
