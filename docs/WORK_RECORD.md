# The detailed work record

Last updated: 2026-08-01
Scripts: `steps/`

This document records every piece of work done, every number measured, every
fault found and every decision taken. The point is that no result should have an
unclear origin, and that every claim can be traced back to the measurement it
came from.

---

## 1. The task and the acceptance criteria

ONT nanopore rDNA amplicon data from anaerobic digester samples was classified
with Kraken2. The task is to design PCR primer pairs for the taxa named in the
design decisions.

### 1.1 The oligo rules

- A, C, G and T only; no degenerate base
- A length of 18 to 25
- GC between 40 and 60 per cent, with a hard bound of 35 to 65
- A G or a C at the 3' end
- At most three G or C in the last five bases
- At most four of the same base in a row

### 1.2 The thermodynamics

- Tm between 58 and 62 degrees, with a hard bound of 57 to 63
- A hairpin dG of at least minus 3000
- A self dimer of at least minus 6000
- A hetero dimer of at least minus 6000
- Tm is measured **with two independent libraries** (primer3 and Biopython) and
  the systematic offset between them is measured from the data

### 1.3 The product

- 70 to 250 bp, at most 300
- 90 to 150 scores best
- Less than 1.5 degrees between the Tm of the two primers
- The start of the product has to equal the forward primer and its end the
  reverse complement of the reverse primer, **exactly**, and this is checked by
  machine

### 1.4 The binding rule, which is PCR physics

- The last two bases have to match exactly
- At most one mismatch in the last five bases
- At most three mismatches in total
- A 5' overhang is free
- The primers sit on opposite strands with their 3' ends facing each other

### 1.5 Specificity

- A product has to form in **every** targeted member
- No product may form in **any** competitor
- At least one primer has to find no binding site at all in the competitors,
  which is the orphan primer
- The pairs are tested again on the raw reads
- Competitor ratios are judged with a **Wilson lower bound**

### 1.6 The general rule

No decision is left to a single code path. Every step takes two independent
measurements, and when they diverge the result is rejected and recorded. A primer
BLAST is mandatory.

### 1.7 The targets and the level asked for

They are coded in `targets.tsv`: twenty targets under four decision headings.

**Decision 1, species specific (6 targets):** Methanosarcina mazei, Methanothrix
soehngenii, Methanosarcina barkeri, Podospora pseudopauciseta, Dictyostelium
discoideum, Trichoderma asperellum

**Decision 2, genus specific (4 targets):** Bacteroides, Alistipes,
Proteiniphilum, Petrimonas

**Decision 3, functional group (7 targets):** hydrogenotrophic methanogens, a
methylotrophic methanogen, acetoclastic methanogens, saccharolytic bacteria,
proteolytic syntrophic bacteria, Nitrosocosmicus AOA, the Trichoderma genus

**Decision 4, universal (3 targets):** bacteria, archaea, fungi

**The tolerance at species specificity:** 1 to 2 cross reacting **species** are
accepted. What is counted is the number of cross reacting species, not the number
of products formed in them.

### 1.8 The working constraints

- Heavy work runs in the user's own WSL, not in a cloud box
- Everything is logged with a date and a time
- Checkpoints are kept against an interruption
- Everything is cross-checked
- Nothing is read from a cache without its checksum being confirmed
- Values that depend on the data are not embedded in the code
- A self audit always runs before delivery

---

## 2. The pipeline

### 2.1 The main chain

| The script | Its job |
|---|---|
| `check_environment.sh` | the dependency check |
| `reclassify_kraken2.sh` | Kraken2 classification |
| `analyze_ambiguous_bases.sh` | the N fraction analysis |
| `generate_primer_candidates.py` | oligo candidates, composition and thermodynamics |
| `design_group_primers.py` | pair formation, the binding and product rule engine |
| `split_clusters.py` | splitting a bin set |
| `anchored_reference_consensus.sh` | the reference anchored consensus |
| `freeze_reference.sh` | pinning the reference |
| `batch_design.py` | batch design per target |
| `specificity.py` | specificity on raw reads, the Wilson lower bound |
| `indistinguishable_targets.py` | finding indistinguishable bins |
| `check_bin_identity.py` | the bin identity check |
| `dominant_allele_consensus.py` | the dominant allele consensus |
| `export_excel.py` | the Excel delivery |

### 2.2 Verification and outside measurement

| The script | Its job |
|---|---|
| `external_databases.py` | the outside database sweep with blastn, the coverage check, the taxon separation |
| `design_from_reference.py` | design from a reference database |
| `check_primer_geometry.py` | the geometry check |
| `regression_test.py` | the regression suite, **147 tests** |
| `check_deliverables.py` | the independent delivery check |
| `mfeprimer_layer.py` | a second independent specificity measurement |
| `recover_bins.py` | consensus recovery without a reference |
| `community_trends.py` | the community trend workbook |
| `target_identity.py` | comparing the target name against the data |
| `reassign_confidence.py` | the Kraken2 confidence threshold, offline |
| `abundance_rank.py` | the rank at which abundance can be read |
| `reference_identity.py` | what the reference primers amplify in the sample |
| `check_taxonomic_level.py` | testing species and genus specificity directly |
| `field_audit.py` | field consistency across archaea, bacteria and fungi, a shared module |

### 2.3 The drivers

| The script | Its job |
|---|---|
| `run.sh` | the main chain |
| `heavy_jobs.sh` | the heavy steps A to H, in WSL |
| `sync.sh` | the SHA256 manifest, verified from both sides |

---

## 3. The faults found and how they were fixed

This is the most important section of the work. Every item carries a fault that
was measured, its size and its fix. Most of them were mistakes made here, and
they are marked as such.

### 3.1 The product length forgot the length of the reverse primer

The product length was measured as the distance between the two 3' ends. The
right measure is the distance from the 5' end of the forward primer to the 5' end
of the reverse primer. The difference is typically 44 bases, so real off target
products of 70 to 94 bp fell below the lower bound and were never counted. It was
fixed and tied to the regression suite: the test now builds a real template and
cuts the primers out of it, so the expected value comes from the sequence itself
rather than from the code's own measure.

### 3.2 The delivery check used the wrong DNA concentration

`check_deliverables.py` used `--dna 250` nM when it measured Tm again, while the
default in the design scripts is 50 nM. The remeasured Tm shifted systematically
and produced **197 false CRITICAL findings**. The default was pulled back to
50 nM.

### 3.3 A field mixture: a bacterial target on a fungal locus

The bins of the saccharolytic bacteria target split 6 to 2 between the bacterial
and the fungal class. When no pair could be found in the correct class, only the
pair on the wrong locus was left, and the table showed the target as covered. A
rule check cannot catch that, because no rule is broken. The shared module
`field_audit.py` was written; the field is derived from the consensus file names
and the taxid list in `targets.tsv`, not from a table filled in by hand. The
batch design produces no cross field design at all, the Excel export leaves them
out of the workbook, and the delivery check marks them CRITICAL on every run.

### 3.4 The specificity step read the excluded taxon file by position

The columns of `dislanan_takson.tsv` were read by position, so a change in the
file format would have turned the exclusion off **silently**. It was changed to
read by header and it now stops loudly on an unrecognised header. Five tests were
added.

### 3.5 The fungal identity: not Trichoderma but Petriella

Kraken2 labelled the fungal bins Trichoderma, Metarhizium and Podospora, and at
one point the conclusion here was "it is really Trichoderma". **That was wrong.**
The dominant organism is from the Microascaceae family, closest to Petriella: in
the ITS region 98.2 per cent against RefSeq and 98.8 to 98.9 per cent against
UNITE.

An intermediate contradiction was resolved along the way. Asking UNITE with the
full consensus returned Hypocreales, with a bit score of 5280 against 872. That
was a length artefact: the full operon records of UNITE produce 3,700 bp
alignments and the conserved 18S and 28S dominate the alignment. Asking with the
ITS window alone, 15 of the first 15 hits came back Microascaceae.

### 3.6 An orientation artefact and a merged rate in bin recovery

`recover_bins.py` reported 72 per cent divergence between two halves, and the
cause was that the seeds came from opposite strands. It was changed to an
alignment based, orientation aware comparison. A second fault surfaced right
after: a coverage shortfall was being folded into the substitution difference.
The real substitution rate was 0.0021 while the merged number came out 0.0145.
Three numbers are given separately now: the substitution rate, the indel rate and
the coverage.

One Blastochloris bin was recovered: because its own reference was 19 per cent
IUPAC, minimap2 could find only 2 minimizers, so a fully covered 1,449 bp
consensus was built from the reads themselves.

### 3.7 A stale file was read without a checksum

A fetch brought back an **old** copy of the trend workbook, with a checksum of
e58d7414 against 62b027a7 on the device. From then on every fetched file is
confirmed against the device by checksum. On one occasion the copy in the mirror
was also edited by hand, which broke the parity; the corrected version was
written to both sides.

### 3.8 The rank percentages added up to 368 per cent

`abundance_rank.py` counted nested nodes of the same rank twice. In a Kraken2
report the real phylum is coded "P", a subphylum "P1" and the one below it "P2",
and all of them fold into the same parent rank; because the clade of the upper
node already contains that of the lower one, counting both is double counting.
Before the fix, the percentages of one barcode added up to **368.54 per cent**.
It was fixed with an ancestor based disjointness rule and all twenty barcodes
settled into the range 99.85 to 100.00 per cent. This was caught by an arithmetic
check of our own before delivery.

### 3.9 The "same organism" warning looked at plurality instead of a majority

The warning fired for the acetoclastic methanogens on 6 of 19 bins, which is a
plurality and not a majority. It was changed to a real majority condition
(`d * 2 > n`) and the support ratios are always shown now.

### 3.10 The level check took comment lines for a header

`reference_targets.tsv` opens with lines beginning with `#`, and `DictReader`
took them for the header and raised a KeyError. They are filtered out.

### 3.11 The level check said "matches" while looking at the genus

The M. barkeri products matched at genus level only and were still written as
matching. Species and genus agreement are reported separately now.

### 3.12 The identity ranking looked at length

`target_identity.py` sorted the hits by length first, so a 94.47 per cent match
over 524 bp came ahead of a 98.21 per cent match over 504 bp and the identity
came out wrong. It was changed to bit score.

### 3.13 The Kraken2 confidence threshold cannot be chosen from memory

The value of 0.1 often recommended for short read data leaves **69 per cent** of
the reads in this ONT data unclassified, because most of the k-mers of an ONT
read find no counterpart in the database and those k-mers go into the denominator
of the score. The threshold was chosen from the data with a scan: 0.02.

`reassign_confidence.py` removes the need to reclassify at all; Kraken2's
`--output` files already carry the k-mer LCA string of every read and the
confidence score is computed from exactly that. Neither the 106 GB database nor
the raw fastq is needed. Measured: **99.84 per cent** of the hit taxa are present
in the tree built from the report files. The self test: at a confidence of 0,
exactly 0 reads move.

### 3.14 A Bracken proposal was refused by our own measurement

Running Bracken was proposed here, and then withdrawn after measuring that only
**17.3 per cent** of the reads reach genus level once the confidence correction
is applied. Bracken assumes the database is complete, and this data does not meet
that assumption.

### 3.15 ROD is a eukaryotic database and had been assigned to archaea and bacteria

**One of the most serious findings of the whole session.**

The lineage in the headers of `ROD_v1.2_operon_variants.fasta` was counted:
**60,320 of the 60,320 records are Eukaryota**, 0 Bacteria, 0 Archaea, with 9,753
fungi among them. ROD is the Ribosomal Operon Database of Krabberod and
colleagues: full length **eukaryotic** rDNA operons extracted from genome
assemblies.

In spite of that, ROD had been assigned to the archaeal and bacterial classes in
the class to database table and not to the fungal ones. That cut both ways:

- **A false cleanliness:** for 71 archaeal and bacterial pairs the sweep wrote
  "no off target product in ROD". That is not evidence of specificity; there is
  not one sequence from that domain in that database. In the same sweep SILVA SSU
  found a product in all 71 records, 41,324 products in total, so the primers
  were not unmeasurable; what was measured was empty.
- **A missed measurement:** the 9,753 full fungal operons in ROD were never
  scanned for the fungal primers at all.

`fungi.18SrRNA.fna`, 4,037 sequences and indexed, was in no class list either.
That was an oversight here too.

**The fix:** the hand written database choice per class was removed. Every class
now sees all ten of the rDNA databases on hand and the list is derived from one
shared source. The reasoning: in a real PCR a primer meets not only the rDNA of
its own domain but all the DNA in the environment, and an archaeal primer
misbinding in a bacterial 23S is exactly the kind of fault being looked for.

### 3.16 An index claim was incomplete

`ref_all.fna` and `ref_all2.fna` had been described as having no index. Correctly
put: they have no BLAST index, but they **do have an mfeprimer index**, in the
older `.primerqc` format plus `.fai` and `.json`. The claim was wrong because
only `.nin` had been looked for. All eight BLAST components
(`.ndb .nhr .nin .njs .not .nsq .ntf .nto`) and the multi volume `.00.*` pattern
are each scanned separately now.

Building those indexes is not needed either. Comparing the identifier sets in the
`.fai` files:

```
archaea.16S     1160
bacteria.16S   26877
fungi.ITS      20394
fungi.18SrRNA   4037
fungi.28SrRNA  12890
in total       65358   <-- ref_all2.fna holds exactly 65358 as well
in ref_all2 but not in the five : 0
in the five but not in ref_all2 : 0
```

`ref_all2.fna` is exactly the union of those five files, and `ref_all.fna` is the
same without 28S and 18S (65,358 minus 12,890 minus 4,037 is 48,431). Both are
already a subset of the files being scanned, so the earlier `makeblastdb`
proposal was withdrawn.

`SILVA_138.2_LSUParc.fasta` is left out deliberately: the Parc set also holds
partial and low quality records, while `LSURef_NR99` is the curated 99 per cent
non redundant representative of the same coverage. On the SSU side NR99 is used
rather than Parc for the same reason.

A side finding: `SILVA_SSURef_NR99.fasta` at 824 MB and
`SILVA_LSURef_NR99.fasta` at 297 MB sat in two places at once and `cmp` confirmed
they were identical. That is 1.1 GB of needless copies.

### 3.17 An empty measurement looked like a clean one, so a coverage check was added

That was the real lesson of the ROD fault. A coverage check was added to
`external_databases.py`: for every pair of a class and a database, before the
sweep, that class's own consensus sequences are searched against the database
with megablast and the longest alignment is measured. The threshold is not
invented, it comes from the data: because the product being looked for is at most
`prod_max` bases long, when the longest similar region in the database is shorter
than `prod_max` that product forming there is already mechanically impossible.

The result: **474 of 970 records came back with no coverage**. So of the 674 rows
that could be read as "no off target product", only 201 were really measuring
anything.

### 3.18 The "off target product" count included the target itself

The hit headers were resolved from the databases. Some of the records giving the
highest counts were **the target itself**:

| The record | The raw count | What the hits really are |
|---|---|---|
| acetoclastic against archaea.16S | 306 | *Methanothrix soehngenii* GP6, Opfikon |
| acetoclastic against SILVA | 1143 | all Methanosaetaceae;Methanothrix |
| Methanothrix soehngenii against SILVA | 308 | all Methanothrix |
| **Nitrosocosmicus against SILVA** | **1119** | **Halobacteria, Methanoperedenaceae, Cenarchaeum** |
| **Petrimonas against SILVA** | **707** | **Clostridium, Bacteroides** |
| **Trichoderma against fungi.28S** | **1150** | **Calonectria, Acremonium, Trichothecium** |

The first three are not a fault but the primer doing its job. Sorting by the raw
count sent the reader off to **fix the wrong primers**.

The fix: every product is separated by the taxon of the reference it formed in,
into `urun_kendi_taksonda`, `urun_yabanci_taksonda` and
`urun_takson_bilinmiyor`. The taxon names are not written by hand; they come from
two sources, `targets.tsv` plus `taxid_names.tsv` for the declared name and
`hedef_kimlik.tsv` for the measured one. Universal targets are marked separately,
because amplifying many taxa at once is the point there and no product counts as
foreign.

### 3.19 "A foreign taxon" was not enough on its own either

Measured after the second round: some of the targets are not a single taxon but a
**functional group**.

| The target | The "foreign" hits | What they really are |
|---|---|---|
| hydrogenotrophic methanogens | *Methanobacterium alcaliphilum*, *Methanosphaera stadtmanae* | both hydrogenotrophic methanogens |
| Nitrosocosmicus AOA | *Nitrosotalea*, *Nitrosopumilus* | both ammonia oxidising archaea |
| Zoopagomycota | *Piptocephalis moniliformis* | itself a Zoopagomycota |
| acetoclastic | *Methanolobus*, *Methanosalsum*, *Methanohalophilus* | the same family |

Against that, Petrimonas hitting *Flavobacterium* and *Phocaeicola vulgatus*, and
Trichoderma hitting *Calonectria* and *Acremonium*, really are distant.

The fix: the distance is measured **from the data** and not from a hand written
table of functional groups. The SILVA, ROD, UNITE and PR2 headers carry the
lineage; the target's own lineage is derived from the references that give a
product in its own taxon; and the depth every foreign hit shares with that
lineage is measured. Sharing down to family level is `yakin`, and not sharing is
`uzak`.

### 3.20 The distance threshold shifted with the number of hits

The reference lineage was taken from the **common prefix** at first. With one hit
the depth came out 7 and with two varied hits 4, so the same foreign hit could
count as distant in one record and close in another. It was changed to the
**dominant** lineage, the threshold was fixed, and a test was tied to it.

### 3.21 A pair giving zero products looked clean, at pair level

Three pairs of the saccharolytic bacteria target gave no product in any of the
eight covered databases, **not even in its own taxon**, with `kendi=0`. The
summary counted those as "giving no product in any foreign taxon". That is the
ROD fault again, one level down. A pair with `kendi=0` is no longer counted
clean; it is reported separately as `OLCUM GECERSIZ`.

### 3.22 The level check silently dropped the reference pairs

In `primer_referans.tsv` the name is `Methanosarcina_barkeri_referans`, and
stripping the `_referans` suffix gives `Methanosarcina_barkeri`, while the name in
`targets.tsv` is `Methanosarcina_barkeri_turu`. When the match failed, the
target's **only** primer set, since it has no de novo pair at all, dropped
silently and the target showed as having no pair. For Proteiniphilum and Podospora
the names matched exactly, so the fault stood out on one target alone. An exact
match is tried first and then a prefix match, and when neither holds it says
loudly that the target was left out of the measurement.

### 3.23 "sp" was being counted as a species name

Records of the form `s__Trichoderma_sp` in UNITE counted as a species name and
stood in the panel like a separate species: 16,910 records for Trichoderma alone,
1,712 for Marasmius and 1,326 for Podospora. A product forming in a record whose
species is not named does not show that the species separation failed. `sp`,
`spp`, `cf` and `aff` no longer count as a species name. `Trichoderma sp.` in
RefSeq was never caught by this because of the full stop; only the underscored
forms were slipping through.

### 3.24 The genus level measurement hid what fell outside the genus

The old count only said in how many species of the panel there is a product. Two
of the five Proteiniphilum pairs also amplify *Fermentimonas caenicola*, which is
another genus, and the old label wrote both as covering three of three, so the
pair breaking genus specificity looked like the pair with the widest coverage.
The genus verdict was split into `CINS_OZGUL`, `CINS_AYRIMI_YOK` and
`CINS_ICINDE_URUN_YOK`. The measured identity does not count as the target here:
if Proteiniphilum was asked for, then Fermentimonas is a cross amplification even
when the organism in the bins is that one.

### 3.25 The RefSeq type strain records never entered the panel

The RefSeq ITS and 28S records look like `NR_172285.1 Petriella musispora CBS
745.69 ITS region; from TYPE material`. In the earlier version the RefSeq branch
did not run at all when the header held a semicolon, so the species name could
not be pulled out of a **type material** record and those records never entered
the panel. The panel is exactly where they were needed most: type strain
sequences are the gold standard of species separation. For Petriella alone, 35
records were falling out this way. SILVA's lineage headers cannot pass through
that branch anyway, so no separate guard is needed there.

### 3.26 In reference design the competitor set was narrower than the verification panel

**The root cause of the new designs.** `design_from_reference.py` gathered its
competitors from one database and six records per name. The Podospora design was
made with 29 competitor sequences, since `fungi.ITS.fna` holds 14 Podospora
records in total, while the verification panel of the level check holds 242
records and 50 species including UNITE. The design believed it was excluding
*P. anserina* and *P. comata* while amplifying both in the wider panel. So the
cause of the failure was not the design engine but the competitors that were
**not visible** at design time.

It was reproduced synthetically: with the narrow competitor set **3,548 pairs**
pass, and with the wide set **0**. The old behaviour gave 3,548 for both.

The fix: the database column takes more than one file, separated by commas; the
competitor list is not limited to what was written by hand, because the other
species in the target's genus are found **from the data** and added; and the
definition of a species name is taken from the level check so that there are not
two separate copies of it.

### 3.27 A genus name was being searched as a substring

Written as `ic='Bacteroides'`, the selector counted *Parabacteroides* and
*Acetobacteroides* records as **members** of the target. In `bacteria.16S.fna` a
word bounded Bacteroides gives 86 records, while the substring search added 17
Parabacteroides and 1 Acetobacteroides to them; the primer would then be forced to
give a product in other genera too, and genus specificity would be lost at design
time. It was changed to a word boundary match.

---

## 4. The scientific findings

### 4.1 The dominant fungal organism

Microascaceae, closest to *Petriella*; the Trichoderma, Metarhizium and Podospora
labels of Kraken2 are wrong. Two independent databases agree in the ITS region,
which is the region that separates them.

### 4.2 Abundance can only be read at the rank the data supports

- Archaea: genus level, 86 to 97 per cent of the reads
- Bacteria: phylum level, 8.8 per cent
- Fungi: class and phylum level, 0.1 to 1.5 per cent

That this is not an artefact of the threshold was confirmed with a control run at
confidence 0, and read length and read count were explicitly ruled out as
confounders.

### 4.3 What the reference designed primers amplify in the sample

- **Proteiniphilum is confirmed:** 5 of 5 pairs give
  *Proteiniphilum acetatigenes* at 98.3 to 100 per cent
- **M. barkeri does not hold at species level:** the products are
  *M. thermophila* at 97.89 per cent and *M. flavescens* at 98.53 per cent, so it
  has to be labelled at genus level
- **There is no Podospora:** 0 product reads, a third independent confirmation

### 4.4 One tension and how it resolves

The level check says the *M. barkeri* primers do separate at species level
against the reference sequences, with one pair having no cross reaction and four
pairs inside the tolerance in a panel of 16 species. The identity step measured
that the products those primers amplify **in this sample** come back as
*M. thermophila* at 97.9 per cent and *M. flavescens* at 98.5 per cent. That is
not a contradiction: the organism in the sample is identical to no named species
and sits about 2 per cent away from both. The primer design is sound and the
organism in the sample is not *M. barkeri*. The report has to state that
distinction plainly.

### 4.5 Where mfeprimer and blastn diverge

31 diverging records over 27 pairs, every one of them in the direction of
mfeprimer finding something blastn does not. In order of priority: saccharolytic
F2 (12,760, already excluded), Alistipes B (545, 485, 19), Zoopagomycota F1
(137, 135, 135, 104), methylotrophic A1 (60, 11), Nitrosocosmicus (44, 44, 2).

A measured mfeprimer behaviour: run with `--misEnd 3`, a primer whose last base
at the 3' end has been changed gives **the same** number of amplicons as the
untouched primer, 323 against 323. A silent mfeprimer failure is caught as well,
where it prints `no valid db` and still returns 0.

---

## 5. The measurement results

### 5.1 The wide outside database sweep, third round

970 records. 674 of them have no product at all, and of those the coverage is
confirmed for 201. Pairs measured with a covered database and found specific: 67.
Measurable because they give a product in their own taxon: 47. Giving no product
in any distant taxon: 2. Giving no product even in their own taxon, so the
measurement is invalid: 20.

**The two really clean pairs:**

```
Nitrosocosmicus_AOA  A1   own=26  close=10  distant=0
   F ACCGCGTGTCACTATCGC
   R TCGATAGTACCAATTAGGCACCAC

Microascaceae_askomikot  F2   own=1  close=0  distant=0
   F ATCAATAAGCGGAGGAAAAGAAACC
   R CCTCTTCAAATTACAACTCGGACTG
```

The second rests on `own=1`, that is on a single record, which is weak evidence.

**The best measurable pair per target, by product in a distant taxon:**

| The target | Best distant | Close | Own |
|---|---|---|---|
| Nitrosocosmicus_AOA A1 | 0 | 10 | 26 |
| Microascaceae_askomikot F2 | 0 | 0 | 1 |
| Metilotrofik_metanojen A1 | 1 | 2 | 5 |
| Nitrosocosmicus_AOA A2 | 2 | 97 | 43 |
| Methanothrix_soehngenii A1 | 8 | 0 | 72 |
| Hidrojenotrofik A2 | 8 | 0 | 44 |
| Methanothrix_soehngenii A2 | 10 | 0 | 116 |
| Hidrojenotrofik A1 | 11 | 0 | 34 |
| Methanosarcina_mazei A2 | 42 | 10 | 58 |
| Petrimonas_cinsi B | 43 | 26 | 19 |
| Zoopagomycota F1 | 50 | 0 | 6 |
| Asetoklastik A1 | 227 | 0 | 749 |
| Trichoderma_cinsi F2 | 2192 | 85 | 215 |

**The spread inside a target is critical:** of the five Nitrosocosmicus A1 pairs
one gives 0 and another gives 1,098. The choice of pair changes the result
completely, so the list of recommended pairs in the workbook has to be rebuilt
against this criterion.

**Never measurable, because there is no product in their own taxon:**
Alistipes_cinsi B with 5 pairs, Bacteroidaceae_ailesi B with 5,
Sakarolitik_bakteriler F2 with 5, Zoopagomycota F1 with 4 and Trichoderma_cinsi
F2 with 1.

### 5.2 The level check

The cross reacting species tolerance is 2.

**Species specificity:**

| The target | The result |
|---|---|
| *Methanothrix soehngenii* | 10 of 10 pairs species specific, no cross reaction at all. **It meets the level** |
| *Methanosarcina barkeri* | 1 pair species specific, 4 inside the tolerance (*M. baltica* and others). **It meets the level** |
| *Methanosarcina mazei* | 1 pair inside the tolerance, the single cross reaction *M. soligelidi*. **It meets the level** |
| *Dictyostelium discoideum* | 1 of 5 pairs inside the tolerance, the cross reaction *Marasmius rhyssophyllus*; what it amplifies in the target species really is *D. discoideum*. **It meets the level** |
| *Trichoderma asperellum* | 2 pairs with 3 to 4 cross reactions (*Petriella guttulata*, *setifera*), 1 pair with no product in the target species. What it amplifies in the target species is not *T. asperellum* but *Petriella musispora*. **It does not meet the level** |
| *Podospora pseudopauciseta* | 5 of 5 pairs give no product in the target species at all and amplify its sibling species. **It does not meet the level** |

**Genus specificity:**

| The target | The result |
|---|---|
| **Proteiniphilum** | 3 of 5 pairs genus specific; 2 also amplify *Fermentimonas caenicola*. **It meets the level** |
| **Petrimonas** | 5 of 5 pairs genus specific inside the panel (*P. mucosa*, *P. sulfuriphila*), but the wide sweep finds distant taxa outside the panel. **Partly** |
| Bacteroides | 0 of 65 species. **It does not meet the level** |
| Alistipes | 0 of 22 species. **It does not meet the level** |

**Six of the ten targets meet the level asked for.**

One correction: the first conclusion for Dictyostelium was "there is no panel".
That was wrong. It is true that it has zero records in the three fungal
databases, but PR2, the eukaryotic SSU database, holds 112 records and the pair
amplifies them. The observation "it is not in the fungal databases" had been
generalised into "it is nowhere".

### 5.3 The overall state of the targets

**Fifteen** of the twenty targets are covered by a usable pair. The ones that are
not:

- **Proteolytic syntrophic bacteria:** no pair at all, neither de novo nor from a
  reference. Of 17,664 drops, 17,623 are caused by *Cloacibacillus porcorum*
  alone.
- **Saccharolytic bacteria:** 5 pairs appear to pass but all of them are in the
  fungal class, that is on a fungal locus. The delivery check marks them CRITICAL
  on every run and they do not enter the workbook.
- **Podospora pseudopauciseta:** the organism is not in the sample, so the
  question does not arise.

*Methanosarcina barkeri* and *Proteiniphilum* are covered by the reference set
alone.

---

## 6. The state of the primer BLAST

**It was done**, and on the rDNA side more widely than NCBI's default:

- `external_databases.py`: `blastn -task blastn-short`, ten rDNA databases, five
  classes, 970 records, the coverage check, and the own, foreign, close and
  distant separation
- `mfeprimer_layer.py`: mfeprimer `spec`, a second independent method
- `check_taxonomic_level.py`: blastn against panels of sibling species
- `target_identity.py` and `reference_identity.py`: the identity blastn runs

What Primer-BLAST does algorithmically, which is to search the primers separately
and find hit pairs that meet in opposite orientations inside the product length
window, has been done.

**Two surfaces are still open:**

1. Every database on hand is rDNA. Whether a primer binds somewhere outside the
   rDNA of the genome has not been examined; the `nt` search of NCBI
   Primer-BLAST covers that and it has not been run.
2. `external_databases.py` and `mfeprimer_layer.py` read the final primer table
   only. The pairs designed from a reference never entered the outside database
   sweep, and those are exactly the deliverable pairs of Proteiniphilum.

---

## 7. The regression suite

`regression_test.py`, **147 tests, 147 of 147 passing.** The expected result of
every test is derived from the design decisions or from known mathematics; the
code's own helper functions are not trusted.

| The group | The subject |
|---|---|
| 1 | the reverse complement |
| 2 | the composition rules |
| 3 | the IUPAC chain rule |
| 4 | the binding rule, compared against brute force over 400 trials |
| 5 | product geometry, with a real template |
| 6 to 9 | Wilson, the consensus, indistinguishability, the alignment back end |
| 10 | the raw read sweep against the design rule, for consistency |
| 11 | the outside database product length |
| 12 | the file format contract between the group engine and the specificity step |
| 13 | the outside database set mapping and the coverage check |
| 14 | separating the own taxon from a foreign taxon |
| 15 | the distance of a foreign hit, by lineage depth |
| 16 | the decision level check, species and genus specificity |
| 17 | the competitor set in reference design |
| 18 | the log contract between the group engine and the four scripts that read its counters |

---

## 8. Synchronisation

`sync.sh` produces a SHA256 manifest with a script section and an output section,
and `--dogrula` compares the two sides. That it catches a difference was
confirmed by introducing one deliberately. All 42 script files are byte for byte
identical on both sides.

---

## 9. Preparing the new designs, 2026-08-01

The decision: for organisms that are not in the sample, design **both** from a
reference under the name in the decision **and** against the real organism, and
try all three of the groups that failed.

`reference_targets.tsv` was raised to nine rows:

| The name | Class | The database | The target |
|---|---|---|---|
| Methanosarcina_barkeri_referans | A1 | archaea.16S | *M. barkeri* |
| Proteiniphilum_cinsi_referans | B | bacteria.16S | *P. saccharofermentans*, *P. acetatigenes* |
| Podospora_pseudopauciseta_referans | F1 | fungi.ITS + UNITE | *P. pseudopauciseta* |
| Bacteroidaceae_ailesi_referans | B | bacteria.16S | the Bacteroides genus |
| Alistipes_cinsi_referans | B | bacteria.16S | the Alistipes genus |
| Microascaceae_askomikot_referans | F2 | fungi.ITS + UNITE | *T. asperellum* |
| Petriella_musispora_referans | F2 | fungi.ITS + 28S + UNITE | *P. musispora* |
| Proteolitik_Cloacibacillus_referans | B | bacteria.16S | *C. porcorum*, *C. evryensis* |
| Sakarolitik_Sphaerochaeta_referans | B | bacteria.16S | the Sphaerochaeta genus |

Hand written sibling species names were **removed** from the list: they are found
from the data now, and a hand written list can be incomplete and therefore gives
false confidence. What stays by hand is only the competitors from other genera.

The taxid of *Petriella musispora* could not be established with certainty and
was **not invented**. An optional seventh column, `hedef_tur`, was added to
`targets.tsv` and the species name is written there directly. That column
**replaces** the taxid name rather than adding to it: adding would have let
Trichoderma into the target species set, and a pair amplifying Trichoderma would
have counted as giving a product in the target species.

Three new rows were numbered as a fifth decision and their note column says they
were derived from a measurement and are put forward for approval; they were never
discussed in a meeting and must not look as though they were decided.

---

## 10. What is left to do

1. Rebuilding the report and the workbook with every new measurement: the
   coverage check, the own against foreign and close against distant separation,
   the level check, the rank coverage and the identity findings
2. Running the new reference design and confirming it with the level check
3. Bringing the reference set into the outside database and mfeprimer sweeps
4. An `nt` sweep with NCBI Primer-BLAST, for off target binding outside the rDNA
5. A full rerun for the saccharolytic target with the patched batch design, about
   2.5 hours, so that the field mixture rows leave the raw table
6. Wet lab confirmation

---

## 11. The commands

`$PT` below is the project directory and the commands are run from `steps/`.

```bash
python3 regression_test.py
```

```bash
bash sync.sh --dogrula
```

```bash
bash heavy_jobs.sh --only H
```

```bash
python3 check_taxonomic_level.py --targets targets.tsv --names taxid_names.tsv --final "$PT/primer_final" --reference "$PT/primer_referans/primer_referans.tsv" --db "$PT/REFERANS_DB" --identity "$PT/primer_final/hedef_kimlik.tsv" --threads 4 --out "$PT/primer_final/duzey_denetimi.tsv"
```

```bash
python3 design_from_reference.py --db "$PT/REFERANS_DB" --pt "$PT" --reference-targets reference_targets.tsv --out "$PT/primer_referans"
```

The wide outside database sweep takes about 12 minutes, the level check about 5
minutes, and the new reference design up to half an hour.

---

## 12. Archived intermediate outputs

Every one of these was archived rather than deleted, so that a number in an old
report can still be traced to the file it came from.

| The file | Why it was archived |
|---|---|
| the wide outside sweep, ROD variant | made while ROD was assigned to the wrong class |
| the wide outside sweep, no taxon variant | made before the own against foreign separation |
| the wide outside sweep, no distance variant | made before the close against distant separation |
| the level check, first round | made while the reference match fault and the "sp" fault were still there |
| the level check, no genus variant | made before the genus verdict was added |
| the reference primers, narrow competitor variant | designed with the narrow competitor set |
