#!/usr/bin/env bash
# =====================================================================
# from_raw.sh: from demultiplexed BAMs (or per-barcode fastq) to bins the
#              identity chain can read, WITHOUT a classifier.
#
#   bash steps/from_raw.sh --map barcodes.tsv [--bam-dir DIR | --fastq-dir DIR]
#                          [--out "fastq files"] [--parallel 2] [--threads 2]
#
#   barcodes.tsv   two columns, no header needed: barcode<TAB>name
#                  e.g.  barcode05    A2-1      (see examples/barcodes_example.tsv)
#
# THE STEPS
#   1  split_barcodes.py   BAM -> barcodeNN.fastq        (skipped with --fastq-dir)
#   2  clean_bins.py       one barcode at a time (v4: one organism per bin), --parallel of them at once,
#                          each with --threads minimap2 threads
#   3  select_bins.py      the bins that enter the chain -> --out
#
# RESUMING: a barcode whose BIN_TABLE.tsv and bin files exist is skipped, so the
# same command continues where it stopped. A barcode another process is binning
# right now (pgrep) is skipped and waited for, so two processes never write the
# same folder.
#
# WHY A POOL AND NOT "ALL AT ONCE": minimap2 with two threads uses two cores;
# the study machine had six and shut itself down once at full load. --parallel 2
# keeps two cores free; raise it on a bigger machine.
# =====================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
MAP=""; BAM_DIR=""; FASTQ_DIR=""; OUT="$ROOT/fastq files"; BINS="$ROOT/fastq files_bins"
PARALLEL=2; THREADS=2; MIN_DISK_GB=8
while [ $# -gt 0 ]; do
  case "$1" in
    --map) MAP="$2"; shift 2 ;;
    --bam-dir) BAM_DIR="$2"; shift 2 ;;
    --fastq-dir) FASTQ_DIR="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --bins) BINS="$2"; shift 2 ;;
    --parallel) PARALLEL="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) sed -n 2,26p "$0"; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done
[ -n "$MAP" ] && [ -f "$MAP" ] || { echo "--map barcodes.tsv is required" >&2; exit 2; }
[ -n "$BAM_DIR" ] || [ -n "$FASTQ_DIR" ] || { echo "give --bam-dir or --fastq-dir" >&2; exit 2; }
LOG="$ROOT/from_raw.log"
log(){ printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG"; }
log "=== from_raw: map=$MAP parallel=$PARALLEL threads=$THREADS ==="

if [ -z "$FASTQ_DIR" ]; then
  FASTQ_DIR="$ROOT/sequences/barcodes"
  log "--- 1 split_barcodes"
  python3 "$HERE/split_barcodes.py" --bam-dir "$BAM_DIR" --out "$FASTQ_DIR" >>"$LOG" 2>&1 || { log "split_barcodes FAILED"; exit 1; }
fi

bin_one(){ # bin_one <barcode> <name>
  local bc="$1" name="$2" t0
  t0=$(date +%s)
  if nice -n 10 python3 "$HERE/clean_bins.py" --fastq "$FASTQ_DIR/$bc.fastq" --name "$name" \
       --out "$BINS" --threads "$THREADS" >>"$LOG" 2>&1; then
    log "    done $name ($(( ($(date +%s)-t0)/60 )) min): $(grep -c '^BIN' "$BINS/$name/BIN_TABLE.tsv") bins"
  else
    log "    FAILED ($name)"
  fi
}
running(){ pgrep -fc "clean_bins.py.*--name " || true; }   # every binning process on the machine

log "--- 2 clean_bins (pool of $PARALLEL)"
while IFS=$'\t' read -r bc name _; do
  [ -n "$bc" ] && [ "${bc#\#}" = "$bc" ] || continue
  if [ -s "$BINS/$name/BIN_TABLE.tsv" ] && ls "$BINS/$name/"*BIN*.fastq >/dev/null 2>&1; then
    log "  $name ($bc): bins exist, skipped"; continue
  fi
  if pgrep -f "clean_bins.py.*--name $name( |$)" >/dev/null; then
    log "  $name ($bc): another process is binning it, skipped and waited for"; continue
  fi
  [ -s "$FASTQ_DIR/$bc.fastq" ] || { log "  $name ($bc): no reads file $FASTQ_DIR/$bc.fastq, skipped"; continue; }
  free=$(df -BG --output=avail "$BINS" 2>/dev/null | tail -1 | tr -dc '0-9')
  if [ "${free:-0}" -lt "$MIN_DISK_GB" ]; then log "STOPPED: disk ${free} GB < $MIN_DISK_GB GB"; wait; exit 3; fi
  while [ "$(running)" -ge "$PARALLEL" ]; do sleep 60; done
  log "  $name ($bc) binning (disk ${free} GB, pool $(( $(running) + 1 ))/$PARALLEL)"
  bin_one "$bc" "$name" &
  sleep 5
done < "$MAP"
wait
while pgrep -f "clean_bins.py.*--name " >/dev/null; do sleep 60; done
missing=0
while IFS=$'\t' read -r bc name _; do
  [ -n "$bc" ] && [ "${bc#\#}" = "$bc" ] || continue
  [ -s "$BINS/$name/BIN_TABLE.tsv" ] || { log "  MISSING: $name"; missing=$((missing+1)); }
done < "$MAP"
[ "$missing" -eq 0 ] || { log "STOPPED: $missing barcode(s) without bins; run the same command again."; exit 5; }
log "bin files: $(find "$BINS" -name '*BIN*.fastq' | wc -l)"

log "--- 3 select_bins -> $OUT"
python3 "$HERE/select_bins.py" --source "$BINS" --target "$OUT" >>"$LOG" 2>&1 || { log "select_bins FAILED"; exit 4; }
log "chosen bins: $(find "$OUT" -name '*BIN*.fastq' | wc -l) (SELECTION_TABLE.tsv beside them)"
log "=== from_raw finished; next: bash steps/anchored_reference_consensus.sh, ./primerjury fungi (fungal bins), ./primerjury consensus, then ./primerjury run ==="
