# The step scripts: reclassification, ambiguous base analysis and design

This directory holds the shell and Python steps that run before and around the
panel design. Two branches run in parallel and neither waits for the other. The
Kraken2 branch reclassifies the raw reads. The ambiguous base branch measures the
N positions and produces a mask, and it works from the consensus and fastq files
that are already there.

## Kraken2, a 106 GB database and 16 GB of RAM

`reclassify_kraken2.sh` compares the database size against the available RAM and
turns the `--memory-mapping` flag on by itself, so a 106 GB database can be run
with 16 GB of RAM. The source study ran pluspf16 with memory mapping in the same
way. If you have a larger database, pass it with `--db` and the script makes the
same decision again against its size.

Speed is the one unknown, which is why the script has an `--only-benchmark` mode:
it processes the smallest file, measures the real time and prints a number for
the total. That is not a gate but planning data; you look at the measurement and
then start the full run with `--yes`. If you want the run to take less time, the
two most direct things to change are raising the thread count with `--threads`
and moving the database onto the fastest disk you have.

## The raw read files

The files under the project's `fastq files` directory are not raw reads. They are
the five most abundant taxa extracted from each sample, which is the shape the
downstream evaluation wanted. The reclassification is done over the raw files
instead: pass the raw fastq directory to `reclassify_kraken2.sh` with `--in`.
Reclassifying the subsets would be pointless, because that set was already chosen
by pluspf16's own decisions, and it narrows the search for anything those
decisions ruled out.

## The consensus script of the source study carries a silent fault

This was confirmed by installing samtools 1.19.2 in a sandbox and running a
synthetic BAM through it. It is a measurement, not a guess. `consensus2.sh` calls
`samtools consensus` twice:

```
samtools consensus -a --use-qual              # "degenerate consensus" in the comment
samtools consensus -a --use-qual -c 0.9       # "strict" in the comment
```

According to the help output of samtools 1.22.1 the default mode is `bayesian`,
and `-c/--call-fract` together with `-q/--use-qual` are options of the `simple`
mode alone. The test gave this:

| The call | The output at a position splitting 60 to 40 |
|---|---|
| `-a --use-qual` | **N** |
| `-a --use-qual -c 0.9` | **N**, byte for byte the same as the first |
| `-a -A` | **M**, that is A or C |
| `-a -m simple --use-qual -c 0.9` | N |
| `-a -m simple -A -c 0.9` | M |

So `-c 0.9` has no effect at all in bayesian mode and the output of the two calls
in the script is identical. More important still, because `-A/--ambig` is never
given, samtools cannot print an IUPAC code under any condition; it writes an N
even at a position that really is biallelic.

The consequence: the 1,897 N's in the `_consensus_strict.fasta` files merge two
completely different causes indistinguishably. Some of them are low read depth
and some of them are real within strain variability.
`analyze_ambiguous_bases.sh` separates the two. The "degenerate consensus" in the
directory was never degenerate, so it has to be produced again before any
comparison means anything.

## The order to run things in

### 1. The environment check, which changes nothing

```bash
bash check_environment.sh 2>&1 | tee environment_report.txt
```

It measures the cores and the RAM, lists the tool versions, writes the help
output of `kraken2`, `bracken`, `samtools consensus`, `minimap2` and `blastn`
into a directory, finds the Kraken2 databases on disk by looking for `hash.k2d`,
reports the size of each one and whether the Bracken `kmer_distrib` files are
present, looks for the raw barcode fastq files, and measures the median read
length out of the existing Kraken output files.

### 2. The Kraken2 speed measurement

```bash
bash reclassify_kraken2.sh --db ~/k2db --in /path/to/raw/fastq --out <project>/kraken_new --only-benchmark
```

Seeing 16 GB of RAM against a 106 GB database, the script turns `--memory-mapping`
on by itself and writes that into the log. It processes the smallest file,
measures the time, prints the estimate and stops. Once you have seen the
estimate, `--yes` in place of `--only-benchmark` starts the full run.

The safeguards in the script:

- Every flag it will use is looked for in the `kraken2 --help` output first, and
  it stops when a flag is not there.
- The memory mapping decision is derived from the RAM and the database size, not
  written by hand. It can be overridden with `--force-mmap` and `--no-mmap`.
- The thread count is `nproc` minus two, and a value given by hand that exceeds
  the core count raises a warning.
- A base name clash is detected up front and the script stops. Because
  `barcode03.fastq` and `barcode03.fastq.gz` come down to the same output name,
  one of them would be dropped without a word; duplicate copies in one directory
  caused exactly that kind of silent double counting once.
- An existing output is skipped, never overwritten.

### 3. The ambiguous base analysis, which does not wait for the classification

```bash
bash analyze_ambiguous_bases.sh --pt <project> --out <project>/N_analysis
```

For every target it aligns the reads onto their own consensus with
`minimap2 -ax map-ont`, then takes two independent measurements: a consensus that
allows IUPAC codes with `samtools consensus -a -A`, and the depth, the IUPAC call
and the base string per position with `samtools consensus -a -A -f PILEUP`. The
classification rests on the PILEUP position numbers, not on the FASTA length, and
when the two disagree the number of disagreements is logged.

Every N position falls into one of four classes:

| The class | The criterion | The consequence |
|---|---|---|
| `iki_allelli` | the depth is enough, the IUPAC code names two bases, and the second base is above the threshold | real within strain variability, permanently forbidden for a primer footprint |
| `dusuk_derinlik` | the depth is below the threshold | there is no read support, forbidden under the mask |
| `kurtarilabilir` | the depth is enough and a single base is clearly in the majority | the sequence is not changed, it is only noted |
| `belirsiz` | anything outside the three above | forbidden, and to be looked at separately |

The thresholds are derived from the data. The depth threshold is ten per cent of
the median of the per target median depth with a floor of 5, and the second base
ratio threshold is samtools' own default `--het-fract` value of 0.15. Both can be
overridden with `--min-depth` and `--het-fract`, and the value used together with
the way it was derived is written into the output.

If your nanopore data is R10.4, adding `--config r10.4_sup` or `--config r10.4_dup`
corrects the quality calibration. It is not added by default, because the
basecalling model cannot be known from the files alone; the script checks whether
the value asked for is defined in the samtools version it finds.

The outputs: `N_pozisyonlari.tsv` with one row per N position, `hedef_ozeti.tsv`
with one row per target, `maske/*.bed` with the regions forbidden to primer
placement, and `ambig/*.fa` plus `pileup/*.txt` with the raw measurements.

### 4. The Bracken re-estimation

Two pieces of information are needed before this step: which
`database<K>mers.kmer_distrib` files the new database holds, and the measured
median read length. The source study used `database300mers.kmer_distrib`, while
the read lengths in its Kraken output are around 4 kb, which means a distribution
file of a suitable length may have to be produced with `bracken-build`.
`check_environment.sh` measures both.

The `bracken_species.sh` of the source study carries two faults. The leading
slash of its paths is missing and patched over with `/"$VAR"`, and its
`basename ... _report.txt` pattern does not match files with a `.report`
extension.

## What was cleaned out of the working tree, and why

26 GB was moved into a `_to_delete` directory rather than deleted, because a
mounted drive gives no delete permission from inside WSL. What it held is worth
knowing, because each item is a trap somebody else can walk into:

- A half finished UNITE mfeprimer index workspace under the reference database
  directory, 24 GB, plus its 412 MB log. The last 1 MB of the log was left in
  place, because that is where the failure is visible.
- SILVA copies confirmed identical by md5: `SILVA_SSURef_NR99.fasta` and its
  `.gz`, `SILVA_LSURef_NR99.fasta` and its `.gz`. The `SILVA_138.2_*` versions
  were kept.
- Five fastq files of one barcode that had been copied into another barcode's
  directory by mistake. Those copies are the silent double counting the base name
  clash check now catches.
- A Bun runtime under `tools`, which is not a primer tool at all, a stray Bracken
  report copy, and a single loose fasta.

## mfeprimer

Do not use the Windows `mfeprimer.exe`. That build is from the 4.2 generation and
cannot read the `v5-single-strand-64` index format the reference database holds.
The Linux ELF `tools/mfeprimer` is from the 4.4 generation and is the binary that
produced those indexes. For a clean tagged release:

```bash
wget https://github.com/quwubin/MFEprimer-3.0/releases/download/v4.4.0/mfeprimer-4.4.0-linux-amd64.gz
```

```bash
gunzip mfeprimer-4.4.0-linux-amd64.gz && chmod +x mfeprimer-4.4.0-linux-amd64
```

The `--bind-amp-only` flag that v4.4.0 brought is useful in the specificity
filter: instead of doing a thermodynamic calculation on every genome wide k-mer
match in the competitor databases, it reports only the bindings inside the
predicted amplicon. With 16 GB of RAM and large references that difference
matters.

## The ambiguous nucleotide picture, in summary

All ninety nine consensus files were scanned from the source disk. The alphabet
is A (56,064), C (66,639), G (59,428), T (60,139) and N (1,897) and nothing else;
not one of R, Y, S, W, K, M, B, D, H or V occurs, and there is no gap character.
Across 244,167 bp that is an N rate of 7.8 per thousand. Nineteen of the ninety
eight files are completely clean.

The distribution shows this is mostly a coverage problem. The worst file is
`A2-3-reads_118126` with 137 N's in 4,381 bp, spread over 118 separate blocks
whose longest run is 3 bp, and the fastq of that taxon is only 376 KB. The one
structural exception is `F1-1-reads_2093779`, which carries a single N block of
43 bp.

After masking the N's plus a three base safety margin, all 98 files still hold at
least one clean 25 bp window, and in only two files (`B-4-reads_1642646` and
`B-4-reads_1642647`) does the longest clean block stop at 145 bp. Because only
the two primer footprints have to be clean, and an N inside the product is
acceptable, the real constraint is far looser than that.

---

# generate_primer_candidates.py, the design rules applied directly

This script applies the oligo, thermodynamic, product and region rules straight
from the design decisions. It waits for neither the classification nor the
ambiguous base analysis; without a mask file it works over the whole consensus,
and with one the masked positions never enter a primer footprint.

```bash
python3 generate_primer_candidates.py --consensus <consensus.fasta> --mask <mask.bed> --out <candidates.tsv>
```

## Two independent Tm measurements, tested against real data

The rule was this: the temperature is computed with two independent libraries and
an oligo that departs from the constant offset between them by more than the
tolerance is dropped. The script does not take the offset by hand, it measures it
from the median difference across every candidate on that target. Measured on a
real *Methanosarcina mazei* consensus (1,445 bp) over 5,782 candidates:

| The measurement | The value |
|---|---|
| primer3 minus Biopython, the median offset | minus 1.50 degrees |
| the standard deviation of the offset | 0.44 degrees |
| oligos outside the tolerance | 0 |

The offset is systematic and narrow, so the assumption the rule rests on holds in
this data. The tolerance defaults to 2.0 degrees and can be changed with
`--tm-cross-tol`.

## The rules that are applied

The oligo: a length of 18 to 25, a GC fraction preferred between 40 and 60 per
cent with a hard bound of 35 to 65, a G or a C at the 3' end, at most three G or
C in the last five bases, and at most four repeats of the same base. The
thermodynamics: a Tm preferred between 58 and 62 with a hard bound of 57 to 63, a
hairpin dG of at least minus 3000, a self dimer of at least minus 6000 and a
hetero dimer inside the pair of at least minus 6000. The product: 70 to 250 and
up to 300, with 90 to 150 scoring best, a product GC preferred between 40 and 60,
and less than 1.5 degrees between the Tm of the two primers. Every threshold can
be changed from the command line and the value used is written at the head of the
output.

The product is machine verified on every candidate: the piece cut out of the
template is compared against the forward primer at its start and against the
reverse complement of the reverse primer at its end. Because it is cut from the
same template, zero failures are expected at this stage; what the check is really
for shows up when it runs against competitor consensuses and raw reads.

## The degenerate base policy

`--degeneracy-budget 0` is the default, so there is no degenerate base in an
oligo. For a function group primer or an anaerobic universal primer, though, a
single ACGT sequence usually cannot catch every member. Limited degeneracy can be
opened for those two cases:

```bash
--degeneracy-budget 2 --degeneracy-fold-max 4
```

The script then allows at most two degenerate positions and a total fold of four,
and it forbids a degenerate base in the last five bases under every condition,
because that is where extension starts. For species and genus specific primers
the budget is better left at zero, since specificity matters more than anything
else there.

## What was measured on real targets

| The target | Length | N | Masked fraction | Valid pairs |
|---|---|---|---|---|
| taxid 2209, *M. mazei* | 1,445 bp | 0 | 0 per cent | 56,454 |
| taxid 118126, the worst N density | 4,381 bp | 137 | 3.13 per cent | 101,245 |
| taxid 1642646, the shortest clean block | 1,493 bp | 65 | 4.35 per cent | 20,197 |

Even on the worst targets tens of thousands of valid pairs remain, so the N's
block the design on no target at all. The bottleneck comes at the specificity
stage. The output is sorted by penalty score and trimmed to the first 5,000 rows
by default, which `--max-pairs` changes.

Processing one target takes about 14 seconds for 1.5 kb and about a minute for
4.4 kb, so under half an hour for 98 targets.

---

# design_group_primers.py, the multi member engine

A genus specific set, a function group and a universal primer all use the same
engine; the only difference is how large the target set is and how much
degeneracy is allowed. No alignment is used. Candidates are produced from an
anchor consensus, and then every candidate is scanned against every member and
every competitor with the binding rule.

The acceptance criteria: a product has to form in every targeted member, no
product may form in any competitor, and the separation has to be solid, meaning
at least one of the primers finds no binding site at all in the competitors. A
cleanliness that comes from both primers binding weakly and only failing together
is rejected.

## The consensus files carry no orientation normalisation, and the engine fixes it

`consensus2.sh` picks a seed read for each taxon and builds the consensus in that
read's orientation. As a result the consensuses are not normalised against one
another. Measured: one of the five members of the acetoclastic methanogen set,
`A1-4-reads_3078083`, is stored the other way round from the rest. Twenty six of
forty conserved probes bind its minus strand and not one binds its plus strand.

The effect is silent and completely destructive: because the forward primer binds
the minus strand instead of the plus strand in that member, no product forms
there at all, so every pair is dropped for "no product in one member". Before
normalisation was added the acetoclastic group returned zero valid pairs; after
it was added, 4,755 pairs.

On every run the engine produces 40 conserved probes from the anchor, votes each
sequence with them, turns the reversed ones into their reverse complement, and
writes into the log which sequence was turned. The source files are not touched.

## A confirmation: the engine found a known universal primer on its own

One of the conserved oligos the engine found is `GGTTACCTTGTTACGACTTA`, the core
region of the classic 1492R universal 16S primer. The engine found that sequence
without looking at any reference list, from the conservation across five member
consensuses alone. The binding scan was also tested against a brute force
comparison and gave the same mismatch count in every member on three sample
oligos.

## A worked example: the acetoclastic methanogens

The target set is taxid 2209, 2223, 2208, 3078083 and 1434102. The competitors
are taxid 394967, 83984, 224719 and 1406512, that is the hydrogenotrophic and
methylotrophic methanogens.

```
2,041 oligos x 9 sequences scanned
oligos binding every target member   : 897
oligos binding no competitor         : 188
dropped, no product in one member    : 77,746
dropped, a product forms in a competitor : 6,044
dropped, no orphan primer            : 1,064
dropped, hetero dimer dG             : 186
valid pairs                          : 4,755
```

The best candidate:

```
F  ATCTCCGGGCTCTTGCTCTC   Tm 61.4
R  TGGGTCTGCGGCCTATCAG    Tm 61.4
   the product is 88 bp in all five members, and the forward primer binds no competitor
```

The same product length coming out in all five members shows the region really is
conserved. One target group is processed in three seconds.

---

# The independent cross audit and the fixes it forced

Two independent auditors went through the candidate generator and the group
engine rule by rule and wrote tests with both synthetic and real data. Every
defect found was fixed and the run repeated afterwards.

## The heavy defects found in the group engine

| The defect | Its effect | The fix |
|---|---|---|
| `product_len` scanned only the "F on the plus strand, R on the minus strand" configuration | Because the template is double stranded the other configuration is a real product too. A competitor was wrongly counted clean and a member wrongly counted as having no product | Both configurations are scanned |
| When the orientation vote tied at plus = minus = 0 it silently said "not reversed"; in real data that happens on 21 sequences, 3 of which really are reversed | A reversed competitor escaped the specificity check completely | A second independent criterion was added: SSU motifs conserved across all life. When neither can decide, the sequence name is printed with a warning |
| `find_bindings` could not see an N in the template during the seed scan | Across 98 consensuses 3 per cent of the 3' end alignment positions were blind, and 21.8 per cent in the worst file. It produced false "orphan primer" claims | The seed variants now produce N as well |
| The rule allowing a free 5' overhang had never been applied | A blind band of about 24 bp at both ends of every sequence, exactly where the coverage falls off | The 5' overhang is free and the overlapping part has to be at least `--min-overlap` bases |
| A product in a competitor was rejected only when it fell inside the target length window | A 370 bp competitor band was being counted as clean | Every band in a competitor is now rejected, and `--competitor-prod-max 0` means unlimited |
| Two members falling on the same label collapsed silently into one sequence | The output claimed to cover two members | The script stops on a label clash |
| With no competitor given, every oligo was declared an orphan | A specificity guarantee with nothing behind it | `rakip_verilmedi` is written and a warning is printed |
| The mask search was both too broad and silently empty | The coordinates of other taxa were laid over the anchor, or the mask was ignored without a word | The group and the taxid are matched together, and the script stops when neither is found |
| `--degeneracy-budget 2` crashed the script | N counted as a legitimate degenerate base and primer3 blew up | N was taken out of the fold table and no oligo containing an N is produced at all |
| The "worst mismatch" was computed with `min` | The name and the computation contradicted each other | It was made explicit: the best inside a member, the worst across members |

## The defects found in the candidate generator

`--degeneracy-budget > 0` crashed it, which is the same N problem, now fixed.
`--prod-max` was a dead parameter and is now penalised as a soft upper bound.
There was no guard against the F and R footprints overlapping, and one was added.
A mask path that did not exist silently turned masking off, and it now raises an
error. A multi record FASTA was silently concatenated, and it now stops. The
contig column of the BED was ignored, and `--mask-contig` was added. The two
library tolerance was fixed at 2.0 degrees, which is 5.5 standard deviations of
the observed distribution and therefore dropped no oligo at all; the tolerance is
now derived from the data as `--tm-cross-k` times the standard deviation. The
pair Tm difference rule accepted 1.5 and now has to be strictly under it.
`--gc-clamp-last 0` applied the rule to the whole oligo instead of turning it off.
Copies of the same locus shifted by 1 or 2 bp filled the output, so
`--min-locus-spacing` was added.

The machine verification of the product passes by definition at this stage,
because the product is cut from the same copy of the template the primers were
derived from. The script now says so plainly and the column is named
`evet_ayni_kalip`. What the check is really for shows up at the specificity
stage.

## The orientation audit: the problem is far more widespread than expected

All 98 consensuses were audited independently with extended variants of five SSU
motifs conserved across all life:

| The group | Forward | Reversed | No motif |
|---|---|---|---|
| A1 | 3 | 17 | 0 |
| A2 | 5 | 15 | 0 |
| B | 6 | 14 | 0 |
| F1 | 0 | 0 | 19 |
| F2 | 7 | 12 | 0 |
| **In total** | **21** | **58** | **19** |

Of the 79 sequences whose motif can be detected, 58 are stored the other way
round. The orientation is random, because `consensus2.sh` picks a seed read for
each taxon and builds the consensus in that read's orientation, and nanopore
reads come in both orientations. The F1 group carries no SSU motif at all,
because that amplicon is the ITS and 18S region; its orientation is settled with
the anchor probes instead.

What this means: orientation normalisation is not a small correction but a
mandatory step touching 59 per cent of the files. Every comparison between
sequences made without it is wrong.

## The universal results per domain

| The domain | Members | Competitors | Conserved oligos | Valid pairs | Product |
|---|---|---|---|---|---|
| Bacteria | 20 | 40 | 137 | 141 | 129-140 bp |
| Archaea | 40 | 20 | 112 | 268 | 76-87 bp |
| Fungi, the long operon | 19 | 15 | 975 | 6,919 | 145-150 bp |
| Fungi, ITS and 18S | 19 | 15 | 28 | **0** | none |

The best candidates:

```
Bacterial universal
F  ACGACAGCCATGCAGCAC     Tm 61.0
R  ACAAGCGGTGGAGCATGTG    Tm 61.0     product 129-140 bp, the reverse primer is the orphan

Archaeal universal
F  CGGCGTTGAGTCCAATTAAAC  Tm 58.1
R  CGCAAGGCTGAAACTTAAAGG  Tm 58.1     product 79 bp, the forward primer is the orphan

Fungal universal, the long operon
F  TAAGAACGGCCATGCACCAC   Tm 61.0
R  AATTGACGGAAGGGCACCAC   Tm 60.9     product 145-147 bp, the reverse primer is the orphan
```

The single reason no pair was found in the ITS and 18S group is
`F1-1-reads_1159556` (*Ustilaginoidea virens*): on its own it blocks every
combination of the 28 conserved oligos that could give a product, 93 drops out of
93. Because that group is the ITS region the separation between species is high
already; either that member has to leave the set, or two separate primer sets
have to be designed instead of one universal.

On the archaeal side `A1-1-reads_2209` is the most restrictive member, 1,163
drops out of 1,163, and on the bacterial side `B-1-reads_1129264`
(*Sphaerochaeta associata*, a spirochaete) is responsible for 1,878. Both are the
expected behaviour: the most distant relative is the member that narrows the
conserved window.

## What still has to be done

The acceptance criterion for these candidates is not that they work on these 98
consensuses. The real coverage of a universal primer has to be measured against
SILVA 138.2 SSU, which the specificity step does. The same step covers the
`mfeprimer` run with `--bind-amp-only`, the `blastn` sweep over the reference
database, judging the competitor ratio with a Wilson lower bound, and confirming
the result on raw reads.

---

# Splitting one group into two sets, and the larger problem that surfaced

A note on naming first. The `F1-*` directories under the consensus folder carry
the naming of the source study and hold the SHORT amplicon (1,155 to 1,649 bp).
On the run sheet the short fungal primer is F2 and corresponds to barcode09 to
barcode12. Below, "F1" means the consensus directory name, that is the short
fungal amplicon.

## The split was made by designability, not by similarity

`split_clusters.py` tries the set with the group engine, and when no valid pair
comes out it removes the member the engine reports as the most blocking and tries
again. As soon as a valid pair appears it counts the remainder as one set and
takes the removed members into a new one.

Clustering by k-mer similarity was tried and abandoned: in the F1 group the
pairwise Jaccard similarity stays between 0.028 and 0.035, so the measurement is
nothing but noise. In the ITS region that is expected.

## The result

| The set | Members | Valid pairs | The best candidate |
|---|---|---|---|
| SET1 | 16 | 12 | F `GTACTTGTTCGCTATCGGTCTC` (Tm 58.6) R `CGAGTCGGGTTGTTTGGG` (Tm 58.4), product 90-91 bp |
| SET2 | 3 | 844 | F `TGTGCGTTCAAAGATTTGATGATTC` (Tm 59.1) R `TTGGTTCTCGCAACGATGAAG` (Tm 59.2), product 90 bp |

The members of SET2 are `F1-1-reads_1159556`, `F1-1-reads_2093779` and
`F1-1-reads_2093780`, all from the same barcode. The split came out barcode based
rather than taxonomic, and because that is suspicious it was investigated.

## The real problem: one taxon's consensuses cover different regions in different years

The consensuses of the same taxon under different barcodes were compared by
12-mer similarity:

| The taxon | The comparison | Similarity | What it means |
|---|---|---|---|
| 44689 *Dictyostelium* | F1-2 against F1-3 | 0.987 | the same region |
| 44689 | F1-1 against F1-2 | 0.065 | a different region |
| 44689 | F1-2 against F1-4 | 0.105 | a different region |
| 2093779 *Podospora* | F1-1 against F1-4 | 0.750 reversed | the same region, the other way round |
| 2093779 | F1-1 against F1-3 | 0.056 | a different region |
| 2093780 *Podospora* | F1-1 against F1-4 | 0.049 | a different region |
| 2545709 *S. osmophilus* | F1-2 against F1-3 | 0.061 | a different region |
| 101201 *Trichoderma* | F1-2 against F1-4 | 0.922 | the same region |

So the consensuses of one taxon across four years can cover three separate pieces
of the rDNA. That is not a strain difference but a difference in which window was
covered.

The cause is `consensus2.sh`: for every taxon it picks a seed read out of the
reads and builds the consensus around it. When the read lengths and the coverage
change, the seed read falls on another part of the operon and the consensus
settles there. So these consensuses are not a defined amplicon but a local
assembly around an arbitrary seed read.

There are two consequences. First, the two sets above are not a taxonomic split
but a REGION split; SET2 gathers three files of one barcode because they sit in a
different window. Second, and more important, a comparison across years is
invalid in this group as things stand.

This problem belongs to that one group alone. In the A1 group all five members
give the same 88 bp product, in the bacterial group 20 members give 129 to 140
bp, and in the F2 group 19 members give 145 to 150 bp. Those are consistent.

## The proposed fix

Those consensuses have to be rebuilt against A SHARED REFERENCE instead of an
arbitrary seed read. `UNITE_ITS.fasta` and `fungi.18SrRNA.fna` are already
present with a BLAST index. The method: align each taxon's reads to that taxon's
UNITE or 18S reference with `minimap2` and call the consensus in the reference
coordinate system. The consensuses of four years then sit in the same window,
which makes both a comparison across years and a single universal primer
possible.

Until that is done the two sets above are usable, as long as it is kept in mind
that SET2 is not a biological group but three files that happen to sit in a
different rDNA window.

---

# The directory names were corrected, and one earlier inference was refuted

## The directory to barcode mapping was proved independently

Trusting the directory names not at all, the taxon set of every consensus and
fastq directory was matched against the top five taxa of that barcode's Bracken
output. Because the extraction step takes the top five anyway, that match has to
be exact, and it was: in all 20 directories the intersection is 5 of 5 or 4 of 4,
with the second best candidate always well behind. The consensus and the fastq
mappings are identical to each other.

| Directory | Barcode | Directory | Barcode |
|---|---|---|---|
| A1-1..A1-4 | barcode01-04 | F1-1..F1-4 | barcode13-16 |
| A2-1..A2-4 | barcode05-08 | F2-1..F2-4 | barcode09-12 |
| B-1..B-4 | barcode17-20 | | |

That is exactly the group to barcode mapping on the run sheet. So the consensus
and fastq directories were right from the start.

## What was renamed

Only the F1 and F2 subdirectories under the Bracken genus, Bracken species and
Kraken results directories were swapped, because those matched the barcodes the
wrong way round. They are all consistent now:

```
bracken results/genus/F1   -> barcode13 14 15 16
bracken results/genus/F2   -> barcode09 10 11 12
bracken results/species/F1 -> barcode13 14 15 16
bracken results/species/F2 -> barcode09 10 11 12
kraken results/F1          -> barcode13 14 15 16
kraken results/F2          -> barcode09 10 11 12
```

A manifest and an undo script were written beside the rename.

## The refuted inference: the Long and Short labels

Earlier, from the consensus lengths alone, the conclusion had been that the run
sheet was right and the source study's directory names were reversed. That was
wrong. Measuring the read lengths straight out of the Kraken output gave this:

| Barcode | Median read | Group | What the run sheet says | What the measurement says |
|---|---|---|---|---|
| barcode01-04 | 1,419-1,428 bp | A1 | Short | Short, correct |
| barcode05-08 | 4,300-4,324 bp | A2 | Long | Long, correct |
| barcode09-12 | 3,697-3,700 bp | F2 | Short | **Long** |
| barcode13-16 | 1,204-1,479 bp | F1 | Long | **Short** |
| barcode17-20 | 1,478-1,483 bp | B | the bacterial panel | 1.5 kb |

The group to barcode mapping on the run sheet is right, but the Long and Short
descriptions of the fungal groups were written the wrong way round. For the
archaeal groups the description is right.

The group descriptions in the workbook were corrected to match: F1 is the short
fungal primer and F2 is the long one. The barcode, year and colour assignments
did not change, because they were already right; only the description column was
corrected, and the measured values were written into the verification and log
sheet.

This correction changes nothing in the earlier findings, only the naming: the
finding that those consensuses cover different regions holds for barcode13 to
barcode16, that is for the short fungal amplicon.

# anchored_reference_consensus.sh

It builds each taxon's consensus in the coordinate system of a shared reference
rather than around an arbitrary seed read.

```bash
bash anchored_reference_consensus.sh --pt <project> --out <project>/reference_consensus --groups F1,F2
```

The flow per taxon: a sample of the reads is taken, the reference giving the
highest total bit score is chosen with `blastn` against a suitable BLAST
database, the reference is pulled out with `blastdbcmd` or `samtools faidx`, the
U's are turned into T for RNA alphabet sources such as SILVA, every read is
aligned to that reference with `minimap2`, and then `samtools consensus` produces
two independent measurements: an IUPAC coded consensus with `-a -A`, and the
depth plus the base string per position with `-a -A -f PILEUP`.

The database is chosen per group: `archaea.16S.fna` for the archaeal groups,
`bacteria.16S.fna` for the bacterial group, and whichever of `fungi.ITS.fna`,
`fungi.18SrRNA.fna` and `fungi.28SrRNA.fna` scores best for the fungal groups.

At the end the script runs an alignment check: it looks at whether the same taxon
sat on the same reference under different barcodes and lists the ones that did
not. Where one did not, the reference has to be pinned by hand for that taxon, or
the years still will not come into one coordinate system.

Because the output is in reference coordinates, every year of one taxon comes out
aligned and orientation normalisation becomes unnecessary as well; the group
engine checks it regardless.

The order: rebuild the consensuses with `anchored_reference_consensus.sh`,
produce the mask with `analyze_ambiguous_bases.sh`, then go back to design with
`design_group_primers.py` and `split_clusters.py`.
