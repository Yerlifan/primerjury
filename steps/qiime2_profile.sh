#!/usr/bin/env bash
# =====================================================================
# qiime2_profile.sh: an independent community profile of every library with
#                    QIIME2 + vsearch, as a second opinion beside the identity chain.
#
#   bash steps/qiime2_profile.sh --root . [--subsample 40000] [--threads 2]
#        [--identity 0.99] [--groups "A2 F1 A1 B F2"] [--out QIIME2_PROFILE]
#
# WHAT IT DOES, per library group (the prefix before "-" of the bin folders):
#   1  all reads of every sample folder are pooled, and a LENGTH WINDOW is taken
#      from that sample's own distribution (10th to 90th percentile, floor 200 bp).
#      No fixed threshold: a 4.3 kb operon library and a 1.5 kb 16S library do
#      not share one. Reads outside the window are fragments or concatemers.
#   2  a RANDOM subsample with a fixed seed (vsearch --fastx_subsample). The first
#      N reads of a nanopore file are the first minutes of the run, not the sample.
#   3  dereplication and de novo clustering at --identity.
#   4  CHIMERA removal (uchime-denovo); both the raw and the cleaned counts are kept.
#   5  TWO thresholds, min 2 and min 10 reads per feature, exported separately, so
#      that a rare taxon lost to the threshold is visible.
#   6  a phylogeny and alpha/beta diversity at the depth of the smallest sample.
#   Every step's read and feature count goes to STEP_COUNTS.tsv.
#
# Needs the qiime2 conda environment active (qiime, vsearch, biom on PATH).
# Memory: the heaviest step is the de novo clustering; keep --threads at 2 on a
# 16 GB machine.
# =====================================================================
set -uo pipefail
ROOT="."; SUB=40000; THREADS=2; IDENT=0.99; GROUPS_GIVEN=""; SEED=42; OUTNAME="QIIME2_PROFILE"
while [ $# -gt 0 ]; do
  case "$1" in
    --root) ROOT="$2"; shift 2;;
    --subsample) SUB="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --identity) IDENT="$2"; shift 2;;
    --groups) GROUPS_GIVEN="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
    --out) OUTNAME="$2"; shift 2;;
    -h|--help) sed -n 2,26p "$0"; exit 0;;
    *) echo "unknown option: $1" >&2; exit 2;;
  esac
done
ROOT="$(cd "$ROOT" && pwd)"
[ -d "$ROOT/fastq files" ] || { echo "ERROR: '$ROOT/fastq files' does not exist" >&2; exit 1; }
KOK="$ROOT/$OUTNAME"; mkdir -p "$KOK"
LOG="$KOK/qiime2_profile.log"; COUNTS="$KOK/STEP_COUNTS.tsv"
log(){ printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" | tee -a "$LOG"; }
err(){ printf '[%s] ERROR: %s\n' "$(date +%H:%M:%S)" "$*" | tee -a "$LOG" >&2; }
count(){ printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" >> "$COUNTS"; }
[ -s "$COUNTS" ] || printf 'library\tsample\tstep\tvalue\n' > "$COUNTS"
for t in qiime vsearch biom; do
  command -v "$t" >/dev/null 2>&1 || { err "$t is not on PATH; activate the qiime2 environment first"; exit 1; }
done
if [ -n "$GROUPS_GIVEN" ]; then GROUPS_LIST="$GROUPS_GIVEN"; else
  GROUPS_LIST="$(for d in "$ROOT/fastq files"/*/; do basename "$d"; done | sed 's/-.*$//' | sort -u | tr '\n' ' ')"
fi
log "=== QIIME2 profile started (subsample $SUB, threads $THREADS, identity $IDENT, groups: $GROUPS_LIST) ==="
log "$(qiime --version 2>&1 | head -1) | $(vsearch --version 2>&1 | head -1)"

for g in $GROUPS_LIST; do
  OUT="$KOK/$g"; mkdir -p "$OUT/raw"
  log "--- library $g ---"
  MAN="$OUT/manifest.tsv"; printf 'sample-id\tabsolute-filepath\n' > "$MAN"; HAVE=0
  for d in "$ROOT/fastq files/${g}-"*; do
    [ -d "$d" ] || continue
    s=$(basename "$d"); READY="$OUT/raw/${s}.fastq"
    if [ ! -s "$READY" ]; then
      ALL="$OUT/raw/${s}_all.fastq"
      cat "$d"/*.fastq > "$ALL" 2>/dev/null
      n0=$(( $(wc -l < "$ALL") / 4 )); count "$g" "$s" "1_raw_reads" "$n0"
      LEN=$(awk 'NR%4==2{print length($0)}' "$ALL" | sort -n)
      P10=$(printf '%s\n' "$LEN" | awk '{a[NR]=$1} END{printf "%d", a[int(NR*0.10)+1]}')
      P90=$(printf '%s\n' "$LEN" | awk '{a[NR]=$1} END{printf "%d", a[int(NR*0.90)+1]}')
      [ "${P10:-0}" -lt 200 ] && P10=200
      count "$g" "$s" "2_length_p10" "$P10"; count "$g" "$s" "2_length_p90" "$P90"
      FILT="$OUT/raw/${s}_filtered.fastq"
      vsearch --fastq_filter "$ALL" --fastq_minlen "$P10" --fastq_maxlen "$P90" --fastq_qmax 93 \
              --fastqout "$FILT" --quiet >>"$LOG" 2>&1
      n1=$(( $(wc -l < "$FILT" 2>/dev/null || echo 0) / 4 )); count "$g" "$s" "3_passed_length_filter" "$n1"
      if [ "$n1" -gt "$SUB" ]; then
        vsearch --fastx_subsample "$FILT" --sample_size "$SUB" --randseed "$SEED" --fastqout "$READY" --quiet >>"$LOG" 2>&1
        count "$g" "$s" "4_random_subsample" "$SUB"
      else
        cp "$FILT" "$READY"; count "$g" "$s" "4_random_subsample" "$n1"
      fi
      rm -f "$ALL" "$FILT"
      log "  $s: $n0 reads -> window ${P10}-${P90} bp -> $n1 -> $(( $(wc -l < "$READY") / 4 )) reads"
    fi
    printf '%s\t%s\n' "$s" "$READY" >> "$MAN"; HAVE=1
  done
  [ "$HAVE" = 1 ] || { err "no sample for $g"; continue; }

  if [ ! -f "$OUT/demux.qza" ]; then
    qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path "$MAN" \
      --input-format SingleEndFastqManifestPhred33V2 --output-path "$OUT/demux.qza" >>"$LOG" 2>&1 \
      || { err "$g could not be imported"; continue; }
  fi
  if [ ! -f "$OUT/feature_table.qza" ]; then
    qiime vsearch dereplicate-sequences --i-sequences "$OUT/demux.qza" \
      --o-dereplicated-table "$OUT/derep_table.qza" --o-dereplicated-sequences "$OUT/derep_seqs.qza" >>"$LOG" 2>&1 \
      || { err "$g could not be dereplicated"; continue; }
    qiime vsearch cluster-features-de-novo --i-table "$OUT/derep_table.qza" --i-sequences "$OUT/derep_seqs.qza" \
      --p-perc-identity "$IDENT" --p-threads "$THREADS" \
      --o-clustered-table "$OUT/feature_table_raw.qza" --o-clustered-sequences "$OUT/rep_seqs_raw.qza" >>"$LOG" 2>&1 \
      || { err "$g could not be clustered"; continue; }
    log "  clustering done ($(awk -v o="$IDENT" 'BEGIN{printf "%.0f", o*100}') per cent identity)"
    if qiime vsearch uchime-denovo --i-table "$OUT/feature_table_raw.qza" --i-sequences "$OUT/rep_seqs_raw.qza" \
         --output-dir "$OUT/chimera" >>"$LOG" 2>&1; then
      qiime feature-table filter-features --i-table "$OUT/feature_table_raw.qza" \
        --m-metadata-file "$OUT/chimera/nonchimeras.qza" --o-filtered-table "$OUT/feature_table.qza" >>"$LOG" 2>&1
      qiime feature-table filter-seqs --i-data "$OUT/rep_seqs_raw.qza" --i-table "$OUT/feature_table.qza" \
        --o-filtered-data "$OUT/rep_seqs.qza" >>"$LOG" 2>&1
      log "  chimera removal applied"
    else
      err "  uchime-denovo failed for $g; chimera removal NOT applied"
      cp "$OUT/feature_table_raw.qza" "$OUT/feature_table.qza"; cp "$OUT/rep_seqs_raw.qza" "$OUT/rep_seqs.qza"
    fi
  else
    log "  feature table already there, skipped"
  fi
  for minf in 2 10; do
    if [ ! -f "$OUT/rep_seqs_min$minf.qza" ]; then
      qiime feature-table filter-features --i-table "$OUT/feature_table.qza" --p-min-frequency "$minf" \
        --o-filtered-table "$OUT/feature_table_min$minf.qza" >>"$LOG" 2>&1
      qiime feature-table filter-seqs --i-data "$OUT/rep_seqs.qza" --i-table "$OUT/feature_table_min$minf.qza" \
        --o-filtered-data "$OUT/rep_seqs_min$minf.qza" >>"$LOG" 2>&1
    fi
    [ -d "$OUT/export_min$minf" ] || {
      qiime tools export --input-path "$OUT/feature_table_min$minf.qza" --output-path "$OUT/export_min$minf" >>"$LOG" 2>&1
      qiime tools export --input-path "$OUT/rep_seqs_min$minf.qza" --output-path "$OUT/export_min$minf" >>"$LOG" 2>&1
      biom convert -i "$OUT/export_min$minf/feature-table.biom" -o "$OUT/export_min$minf/feature-table.tsv" --to-tsv >>"$LOG" 2>&1
      sed -i '1{/^# Constructed from biom file/d}' "$OUT/export_min$minf/feature-table.tsv"
    }
    n=$(grep -c '^>' "$OUT/export_min$minf/dna-sequences.fasta" 2>/dev/null)
    count "$g" "-" "5_features_min$minf" "${n:-0}"; log "  min $minf reads: ${n:-0} features"
  done
  [ -d "$OUT/export_raw" ] || {
    qiime tools export --input-path "$OUT/feature_table_raw.qza" --output-path "$OUT/export_raw" >>"$LOG" 2>&1
    biom convert -i "$OUT/export_raw/feature-table.biom" -o "$OUT/export_raw/feature-table.tsv" --to-tsv >>"$LOG" 2>&1
    sed -i '1{/^# Constructed from biom file/d}' "$OUT/export_raw/feature-table.tsv"
  }
  nr=$(( $(wc -l < "$OUT/export_raw/feature-table.tsv" 2>/dev/null || echo 1) - 1 ))
  count "$g" "-" "5_features_clustered_raw" "$nr"
  if [ ! -f "$OUT/rooted_tree.qza" ]; then
    qiime phylogeny align-to-tree-mafft-fasttree --i-sequences "$OUT/rep_seqs_min10.qza" --p-n-threads "$THREADS" \
      --o-alignment "$OUT/aligned.qza" --o-masked-alignment "$OUT/masked.qza" \
      --o-tree "$OUT/tree.qza" --o-rooted-tree "$OUT/rooted_tree.qza" >>"$LOG" 2>&1 \
      && log "  phylogeny ready" || err "  the tree could not be built for $g"
  fi
  DEPTH=$(awk -F'\t' 'NR==1{next} {for(i=2;i<=NF;i++) t[i]+=$i} END{m=-1; for(i in t) if(m<0||t[i]<m) m=t[i]; printf "%d", m}' \
          "$OUT/export_min10/feature-table.tsv" 2>/dev/null)
  if [ "${DEPTH:-0}" -gt 100 ] && [ ! -f "$OUT/diversity/shannon.tsv" ]; then
    mkdir -p "$OUT/diversity"; log "  diversity (rarefaction depth $DEPTH = the smallest sample)"
    qiime feature-table rarefy --i-table "$OUT/feature_table_min10.qza" --p-sampling-depth "$DEPTH" \
      --o-rarefied-table "$OUT/diversity/rarefied.qza" >>"$LOG" 2>&1
    for m in observed_features shannon simpson; do
      qiime diversity alpha --i-table "$OUT/diversity/rarefied.qza" --p-metric "$m" --o-alpha-diversity "$OUT/diversity/$m.qza" >>"$LOG" 2>&1 \
        && qiime tools export --input-path "$OUT/diversity/$m.qza" --output-path "$OUT/diversity/_$m" >>"$LOG" 2>&1 \
        && mv "$OUT/diversity/_$m/alpha-diversity.tsv" "$OUT/diversity/$m.tsv" 2>/dev/null
      rm -rf "$OUT/diversity/_$m"
    done
    for m in braycurtis jaccard; do
      qiime diversity beta --i-table "$OUT/diversity/rarefied.qza" --p-metric "$m" --o-distance-matrix "$OUT/diversity/beta_$m.qza" >>"$LOG" 2>&1 \
        && qiime tools export --input-path "$OUT/diversity/beta_$m.qza" --output-path "$OUT/diversity/_b$m" >>"$LOG" 2>&1 \
        && mv "$OUT/diversity/_b$m/distance-matrix.tsv" "$OUT/diversity/beta_$m.tsv" 2>/dev/null
      rm -rf "$OUT/diversity/_b$m"
    done
    if [ -f "$OUT/rooted_tree.qza" ]; then
      for m in unweighted_unifrac weighted_unifrac; do
        qiime diversity beta-phylogenetic --i-table "$OUT/diversity/rarefied.qza" --i-phylogeny "$OUT/rooted_tree.qza" \
          --p-metric "$m" --o-distance-matrix "$OUT/diversity/beta_$m.qza" >>"$LOG" 2>&1 \
          && qiime tools export --input-path "$OUT/diversity/beta_$m.qza" --output-path "$OUT/diversity/_u$m" >>"$LOG" 2>&1 \
          && mv "$OUT/diversity/_u$m/distance-matrix.tsv" "$OUT/diversity/beta_$m.tsv" 2>/dev/null
        rm -rf "$OUT/diversity/_u$m"
      done
    fi
    count "$g" "-" "6_rarefaction_depth" "$DEPTH"
  fi
  log "$g done: $OUT"
done
log "=== QIIME2 profile finished; next: bash steps/qiime2_classify.sh --root $ROOT --profile $OUTNAME ==="
