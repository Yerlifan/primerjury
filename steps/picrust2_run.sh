#!/usr/bin/env bash
# =====================================================================
# picrust2_run.sh: functional prediction (PICRUSt2) from the QIIME2 profile, with a
#                  memory guard for a 16 GB machine.
#
#   bash steps/picrust2_run.sh --root . [--profile QIIME2_PROFILE] [--groups "A1 A2 B"]
#        [--max-nsti 2] [--ko yes] [--targets targets_ec_ko.tsv]
#        [--min-free-mb 1500] [--swap-ceiling-mb 9000]
#
# WHAT IT DOES, per library group:
#   1  the domain: a group whose name starts with A is placed on the archaeal
#      reference tree, B on the bacterial one; fungal groups are skipped (no tree).
#   2  a long amplicon (mean feature length > 2,200 bases, the full rDNA operon) is
#      CUT to its SSU part first, by vsearch coordinates against SILVA SSU: PICRUSt2
#      places 16S, not operons.
#   3  place_seqs -> hsp (16S copy number + NSTI, EC) -> metagenome pipeline ->
#      MetaCyc pathways; KO on request (heavy).
#   Every step runs in its own process group under a MEMORY GUARD: it is not
#   started below --min-free-mb of free memory, and it is killed when swap use
#   stays above --swap-ceiling-mb for three checks. Measured on the study machine:
#   the step that eats memory is the trait-table check at start-up, not the
#   calculation; a killed step is reported, never silently skipped.
#   4  --targets: a two-column TSV (code, description) of EC/KO codes to collect
#      into one summary table per group and sample.
#
# Needs the picrust2 conda environment (looked for under ~/miniconda3, ~/anaconda3,
# /opt/conda) and vsearch (the qiime2 environment) for the SSU cut.
# =====================================================================
set -uo pipefail
ROOT="."; PROFILE="QIIME2_PROFILE"; GROUPS_GIVEN=""; NSTI=2; KO="no"; TARGETS=""; MINFREE=1500; SWAPCEIL=9000
while [ $# -gt 0 ]; do
  case "$1" in
    --root) ROOT="$2"; shift 2;;
    --profile) PROFILE="$2"; shift 2;;
    --groups) GROUPS_GIVEN="$2"; shift 2;;
    --max-nsti) NSTI="$2"; shift 2;;
    --ko) KO="$2"; shift 2;;
    --targets) TARGETS="$2"; shift 2;;
    --min-free-mb) MINFREE="$2"; shift 2;;
    --swap-ceiling-mb) SWAPCEIL="$2"; shift 2;;
    -h|--help) sed -n 2,27p "$0"; exit 0;;
    *) echo "unknown option: $1" >&2; exit 2;;
  esac
done
ROOT="$(cd "$ROOT" && pwd)"; KOK="$ROOT/$PROFILE"
[ -d "$KOK" ] || { echo "ERROR: profile directory missing: $KOK" >&2; exit 2; }
REF="$ROOT/REFERENCE_DB"; LOG="$KOK/picrust2_run.log"; DETAIL="$KOK/picrust2_detail.log"
log(){ printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" | tee -a "$LOG"; }
err(){ printf '[%s] ERROR: %s\n' "$(date +%H:%M:%S)" "$*" | tee -a "$LOG" >&2; }
PIC=""
for k in "$HOME/miniconda3" "$HOME/anaconda3" "/opt/conda"; do
  [ -x "$k/envs/picrust2/bin/picrust2_pipeline.py" ] && { PIC="$k/envs/picrust2/bin"; break; }
done
[ -n "$PIC" ] || { err "the picrust2 environment was not found"; exit 1; }
export PATH="$PIC:$PATH"
DEFAULTS="$(ls -d "$(dirname "$PIC")"/lib/python*/site-packages/picrust2/default_files 2>/dev/null | head -1)"
[ -d "$DEFAULTS" ] || { err "picrust2 default_files not found"; exit 1; }
free_mb(){ awk '/^MemAvailable:/{printf "%d", $2/1024}' /proc/meminfo; }
swap_mb(){ awk '/^SwapTotal:/{t=$2} /^SwapFree:/{f=$2} END{printf "%d", (t-f)/1024}' /proc/meminfo; }
guarded(){ # guarded <label> <command...>
  local label="$1"; shift
  local free; free="$(free_mb)"
  if [ "$free" -lt "$MINFREE" ]; then err "$label skipped: $free MB free, floor $MINFREE MB"; return 1; fi
  local t0; t0=$(date +%s)
  setsid nice -n 10 "$@" >>"$DETAIL" 2>&1 &
  local p=$! peak=0 peaks=0 used b s over=0
  while kill -0 "$p" 2>/dev/null; do
    b=$(free_mb); s=$(swap_mb); used=$(( free - b ))
    [ "$used" -gt "$peak" ] && peak=$used; [ "$s" -gt "$peaks" ] && peaks=$s
    if [ "$s" -gt "$SWAPCEIL" ]; then
      over=$(( over + 1 ))
      if [ "$over" -ge 3 ]; then err "  $label: swap $s MB above the ceiling $SWAPCEIL MB; the process group is being killed"; kill -9 -"$p" 2>/dev/null; break; fi
    else over=0; fi
    sleep 5
  done
  wait "$p" 2>/dev/null; local rc=$?; local dur=$(( $(date +%s) - t0 ))
  if [ "$rc" = 0 ]; then log "  $label done (${dur} s, peak memory ~${peak} MB, swap ~${peaks} MB)"
  else err "  $label failed (exit $rc, ${dur} s, peak memory ~${peak} MB, swap ~${peaks} MB)"; fi
  return $rc
}
cut_ssu(){ # cut_ssu <group> <in.fasta> <out.fasta>
  local g="$1" in="$2" out="$3" ssu="$REF/SILVA_138.2_SSURef_NR99.fasta"
  [ -s "$ssu" ] || { err "$g: SILVA SSU missing: $ssu"; return 1; }
  local VS; VS="$(command -v vsearch 2>/dev/null || true)"
  if [ -z "$VS" ]; then
    for k in "$HOME/miniconda3/envs"/*qiime2*/bin/vsearch "$HOME/anaconda3/envs"/*qiime2*/bin/vsearch; do [ -x "$k" ] && { VS="$k"; break; }; done
  fi
  [ -n "$VS" ] || { err "vsearch not found in any environment"; return 1; }
  local u="$KOK/$g/ssu_coordinates.tsv"
  if [ ! -s "$u" ]; then
    log "  $g: measuring SSU coordinates (vsearch)"
    "$VS" --usearch_global "$in" --db "$ssu" --id 0.75 --top_hits_only --maxaccepts 1 --maxrejects 64 \
          --strand both --threads 2 --quiet --userout "$u" --userfields "query+qilo+qihi+qstrand+id" >>"$DETAIL" 2>&1
  fi
  [ -s "$u" ] || { err "$g: no SSU coordinates"; return 1; }
  awk -v coord="$u" '
    BEGIN{ while ((getline line < coord) > 0) { split(line, a, "\t"); if (!(a[1] in lo)) { lo[a[1]]=a[2]; hi[a[1]]=a[3]; st[a[1]]=a[4] } } }
    /^>/ { if (id != "") emit(); id=substr($1, 2); seq=""; next }
    { seq = seq $0 }
    END { if (id != "") emit() }
    function emit(   p, len, s, i, c, t) {
      if (!(id in lo)) return
      p = lo[id]; len = hi[id] - lo[id] + 1
      if (len < 900 || len > 2200) return
      s = substr(seq, p, len)
      if (st[id] == "-") { t = ""; for (i = length(s); i > 0; i--) { c = substr(s, i, 1); t = t (c=="A" ? "T" : c=="T" ? "A" : c=="G" ? "C" : c=="C" ? "G" : "N") }; s = t }
      print ">" id "\n" s
    }' "$in" > "$out"
  local n; n=$(grep -c '^>' "$out" 2>/dev/null)
  [ "${n:-0}" -gt 0 ] || { err "$g: the SSU cut produced nothing"; return 1; }
  log "  $g: SSU part cut from $n features"
}
if [ -n "$GROUPS_GIVEN" ]; then GROUPS_LIST="$GROUPS_GIVEN"; else GROUPS_LIST="$(for d in "$KOK"/*/; do basename "$d"; done | tr '\n' ' ')"; fi
log "=== PICRUSt2 started (groups: $GROUPS_LIST, KO: $KO, swap ceiling ${SWAPCEIL} MB) ==="
log "picrust2: $("$PIC/picrust2_pipeline.py" --version 2>&1 | tail -1)"
for g in $GROUPS_LIST; do
  OUT="$KOK/$g"; PO="$OUT/picrust2"
  S="$OUT/export_min10/dna-sequences.fasta"; T="$OUT/export_min10/feature-table.tsv"
  [ -s "$S" ] && [ -s "$T" ] || { err "$g: inputs missing, skipped"; continue; }
  sed -i '1{/^# Constructed from biom file/d}' "$T"; mkdir -p "$PO"
  case "$g" in
    A*) DOMAIN="archaea"; TREE="$DEFAULTS/archaea/arc_ref";;
    B*) DOMAIN="bacteria"; TREE="$DEFAULTS/bacteria/bac_ref";;
    *) err "$g: no PICRUSt2 domain for this group (fungi have no reference tree), skipped"; continue;;
  esac
  TRAITS="$DEFAULTS/$DOMAIN"
  MEAN=$(awk '/^>/{if(l){s+=l;n++} l=0;next}{l+=length($0)} END{if(l){s+=l;n++} printf "%d", s/n}' "$S")
  log "--- $g: $(grep -c '^>' "$S") features, mean $MEAN bp, domain $DOMAIN, free memory $(free_mb) MB ---"
  IN="$S"
  if [ "$MEAN" -gt 2200 ]; then IN="$OUT/export_min10/ssu-sequences.fasta"; [ -s "$IN" ] || cut_ssu "$g" "$S" "$IN" || continue; fi
  [ -s "$PO/out.tre" ] || guarded "placement ($DOMAIN)" place_seqs.py -s "$IN" -o "$PO/out.tre" -p 1 -r "$TREE" --intermediate "$PO/placement_tmp" || continue
  [ -s "$PO/marker_nsti_predicted.tsv.gz" ] || guarded "16S copy number and NSTI" hsp.py -t "$PO/out.tre" -o "$PO/marker_nsti_predicted.tsv.gz" -p 1 -n --observed_trait_table "$TRAITS/16S.txt.gz" || continue
  [ -s "$PO/EC_predicted.tsv.gz" ] || guarded "EC prediction" hsp.py -t "$PO/out.tre" -o "$PO/EC_predicted.tsv.gz" -p 1 --observed_trait_table "$TRAITS/ec.txt.gz" || continue
  [ -s "$PO/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" ] || guarded "EC metagenome" metagenome_pipeline.py -i "$T" -m "$PO/marker_nsti_predicted.tsv.gz" -f "$PO/EC_predicted.tsv.gz" -o "$PO/EC_metagenome_out" --max_nsti "$NSTI" || continue
  [ -s "$PO/pathways_out/path_abun_unstrat.tsv.gz" ] || guarded "MetaCyc pathways" pathway_pipeline.py -i "$PO/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" -o "$PO/pathways_out" -p 1
  if [ "$KO" = "yes" ]; then
    [ -s "$PO/KO_predicted.tsv.gz" ] || guarded "KO prediction (heavy)" hsp.py -t "$PO/out.tre" -o "$PO/KO_predicted.tsv.gz" -p 1 --observed_trait_table "$TRAITS/ko.txt.gz"
    if [ -s "$PO/KO_predicted.tsv.gz" ] && [ ! -s "$PO/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" ]; then
      guarded "KO metagenome" metagenome_pipeline.py -i "$T" -m "$PO/marker_nsti_predicted.tsv.gz" -f "$PO/KO_predicted.tsv.gz" -o "$PO/KO_metagenome_out" --max_nsti "$NSTI"
    fi
  fi
  for f in "$PO/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" "$PO/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
           "$PO/pathways_out/path_abun_unstrat.tsv.gz" "$PO/EC_metagenome_out/weighted_nsti.tsv.gz"; do
    [ -s "$f" ] && gunzip -kf "$f" 2>/dev/null
  done
  NS="$PO/EC_metagenome_out/weighted_nsti.tsv"; [ -s "$NS" ] && log "  $g weighted NSTI: $(awk 'NR>1{printf "%s=%.2f ", $1, $2}' "$NS")"
  log "$g done: $PO"
done
if [ -n "$TARGETS" ] && [ -s "$TARGETS" ]; then
  SUM="$KOK/PICRUST2_TARGETS.tsv"; printf 'library\tcatalogue\tcode\tdescription\tsample\tpredicted abundance\n' > "$SUM"
  collect(){ # collect <group> <file> <catalogue>
    [ -s "$2" ] || return 0
    awk -F'\t' -v grp="$1" -v cat="$3" -v file="$2" 'NR==FNR{ if ($1 != "") desc[$1]=$2; next }
      END { while ((getline line < file) > 0) { n = split(line, s, "\t"); if (++row == 1) { for (i = 2; i <= n; i++) smp[i] = s[i]; continue }
            if (s[1] in desc) for (i = 2; i <= n; i++) printf "%s\t%s\t%s\t%s\t%s\t%.1f\n", grp, cat, s[1], desc[s[1]], smp[i], s[i] } }' "$TARGETS" /dev/null >> "$SUM"
  }
  for g in $GROUPS_LIST; do
    collect "$g" "$KOK/$g/picrust2/EC_metagenome_out/pred_metagenome_unstrat.tsv" EC
    collect "$g" "$KOK/$g/picrust2/KO_metagenome_out/pred_metagenome_unstrat.tsv" KO
  done
  log "target table: $SUM ($(( $(wc -l < "$SUM") - 1 )) rows)"
fi
log "=== PICRUSt2 finished ==="
