#!/usr/bin/env bash
# ==== MEETING NOTE ====
# WHAT IT IS FOR   : Reclassifies the raw barcode fastq files against a 106 GB Kraken2 database.
# INPUT            : --in the raw fastq directory, --db the kraken2 database
# OUTPUT           : a .report and a per read output file for every barcode
# HOW TO RUN IT    : bash reclassify_kraken2.sh --db <db> --in <fastq> --out <out> [--confidence 0.1]
# WHY IT IS LIKE THIS : Kraken2 does not abstain under its default setting; if the real organism is not in the database it labels the read with the nearest leaf. The script looks every flag it uses up in the output of kraken2 --help first, it assumes nothing from memory. It measures on the smallest file first and prints an estimate of the time, and does not start the full run without approval.
# =======================
# =====================================================================
# reclassify_kraken2.sh
# The aim: to reclassify the raw barcode fastq files with a large Kraken2
#          database.
#
# Usage:
#   bash reclassify_kraken2.sh \
#        --db   /path/to/the/kraken2/database \
#        --in   /path/to/the/raw/fastq/directory \
#        --out  /path/to/the/output/directory \
#        [--threads N] [--confidence C] [--force-mmap] [--no-mmap]
#        [--only-benchmark] [--yes]
#
# The design decisions:
#   - Every flag that will be used is looked for in the output of
#     `kraken2 --help` before anything runs, and the script stops when one is
#     not there. No flag is assumed from memory.
#   - The memory mapping decision is derived from the RAM and the database
#     size, it is not written by hand.
#   - A measurement is made on the smallest file first, an estimate of the
#     total time is printed, and the full run does not start without a
#     confirmation (--yes skips it).
#   - It does not overwrite an existing output, it skips it; delete the output
#     to produce it again.
#   - When --confidence is given it goes into the output names (_c0.1 for
#     example), so that runs at different confidence thresholds do not
#     overwrite each other. At its default setting (0) Kraken2 does not abstain:
#     when the real organism is not in the database, the read falls onto the
#     highest scoring sibling leaf. Measured: in a real run the reads of four
#     Methanosarcina bins go to the same references and not one of them prefers
#     the species it was assigned to. A threshold around 0.1 gathers those reads
#     at the genus node.
# =====================================================================
set -euo pipefail

DB=""; IN=""; OUT=""; THREADS=""; MMAP="auto"; ONLYBENCH=0; ASSUME_YES=0
CONF=""
while [ $# -gt 0 ]; do
  case "$1" in
    --db) DB="$2"; shift 2;;
    --in) IN="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --confidence) CONF="$2"; shift 2;;
    --force-mmap) MMAP="on"; shift;;
    --no-mmap) MMAP="off"; shift;;
    --only-benchmark) ONLYBENCH=1; shift;;
    --yes) ASSUME_YES=1; shift;;
    *) echo "unknown option: $1" >&2; exit 2;;
  esac
done
[ -n "$DB" ] && [ -n "$IN" ] && [ -n "$OUT" ] || {
  echo "usage: bash $0 --db <kraken2_db> --in <raw_fastq_directory> --out <output_directory>" >&2; exit 2; }

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

command -v kraken2 >/dev/null 2>&1 || die "kraken2 is not on PATH"

# --- 1. the integrity of the database ---------------------------------
for f in hash.k2d opts.k2d taxo.k2d; do
  [ -e "$DB/$f" ] || die "$f is not in the database: $DB"
done
DBBYTES=$(du -sb "$DB" | cut -f1)
DBGB=$(awk -v b="$DBBYTES" 'BEGIN{printf "%.1f", b/1073741824}')
log "the database: $DB  (${DBGB} GB)"

# --- 2. flag verification ---------------------------------------------
HELP=$(kraken2 --help 2>&1 || true)
need_flag() { grep -q -- "$1" <<<"$HELP" || die "this kraken2 version does not know the $1 flag"; }
for fl in --db --threads --report --output --use-names; do need_flag "$fl"; done
HAS_MMAP=0; grep -q -- "--memory-mapping" <<<"$HELP" && HAS_MMAP=1
# --confidence is not assumed from memory, it is looked for in the help output
if [ -n "$CONF" ]; then
  grep -q -- "--confidence" <<<"$HELP" || die "this kraken2 version does not know the --confidence flag"
  case "$CONF" in
    ''|*[!0-9.]*) die "--confidence has to be a number: $CONF";;
  esac
  awk -v c="$CONF" 'BEGIN{ if (c+0 < 0 || c+0 > 1) exit 1 }' || die "--confidence has to be between 0 and 1: $CONF"
  log "the confidence threshold: $CONF (_c$CONF will go into the output names)"
fi
HAS_GZ=0;   grep -q -- "--gzip-compressed" <<<"$HELP" && HAS_GZ=1
HAS_MINHIT=0; grep -q -- "--minimum-hit-groups" <<<"$HELP" && HAS_MINHIT=1
HAS_RMZ=0;  grep -q -- "--report-minimizer-data" <<<"$HELP" && HAS_RMZ=1
log "flag support: memory-mapping=$HAS_MMAP gzip=$HAS_GZ minimum-hit-groups=$HAS_MINHIT report-minimizer-data=$HAS_RMZ"

# --- 3. the thread and memory decision --------------------------------
CORES=$(nproc)
if [ -z "$THREADS" ]; then THREADS=$(( CORES > 2 ? CORES - 2 : 1 )); fi
[ "$THREADS" -gt "$CORES" ] && log "WARNING: more threads were asked for ($THREADS) than there are cores ($CORES)"
AVAILGB=$(awk '/MemAvailable/ {printf "%.1f", $2/1048576}' /proc/meminfo)
log "cores=$CORES  threads used=$THREADS  available RAM=${AVAILGB} GB"

USE_MMAP=0
case "$MMAP" in
  on) USE_MMAP=1; log "memory mapping: turned on by hand";;
  off) USE_MMAP=0; log "memory mapping: turned off by hand";;
  auto)
    NEED=$(awk -v g="$DBGB" 'BEGIN{printf "%.1f", g/0.8}')
    if awk -v a="$AVAILGB" -v n="$NEED" 'BEGIN{exit !(a>=n)}'; then
      USE_MMAP=0
      log "memory mapping: NOT NEEDED (RAM ${AVAILGB} GB >= the needed ${NEED} GB)"
    else
      USE_MMAP=1
      log "memory mapping: NEEDED (RAM ${AVAILGB} GB < the needed ${NEED} GB), the time will depend on the disk"
    fi;;
esac
[ "$USE_MMAP" = 1 ] && [ "$HAS_MMAP" = 0 ] && die "memory mapping is needed but this version does not support it"

# --- 4. the input files -----------------------------------------------
mkdir -p "$OUT"
mapfile -t FQ < <(find "$IN" -maxdepth 2 -type f \
  \( -iname '*.fastq' -o -iname '*.fq' -o -iname '*.fastq.gz' -o -iname '*.fq.gz' \) | sort)
[ ${#FQ[@]} -gt 0 ] || die "no fastq was found in the input directory: $IN"
log "the number of input files: ${#FQ[@]}"

TOTBYTES=0
for f in "${FQ[@]}"; do TOTBYTES=$(( TOTBYTES + $(stat -c%s "$f") )); done
log "input in total: $(awk -v b=$TOTBYTES 'BEGIN{printf "%.2f", b/1073741824}') GB"

# A base name clash check. Two inputs that come down to the same base (say
# barcode03.fastq and barcode03.fastq.gz) produce the same output name and one of
# them is dropped without a word. Duplicate copies in one directory caused exactly
# that kind of silent loss once, so the script stops here and lists them.
declare -A SEEN=()
COLL=0
for f in "${FQ[@]}"; do
  b=$(basename "$f"); b=${b%.gz}; b=${b%.fastq}; b=${b%.fq}
  if [ -n "${SEEN[$b]:-}" ]; then
    printf 'CLASH: "%s" and "%s" both land on the same output name (%s)\n' "${SEEN[$b]}" "$f" "$b" >&2
    COLL=1
  else
    SEEN[$b]="$f"
  fi
done
[ "$COLL" = 1 ] && die "there is a base name clash. Stopped to prevent a silent loss. Take the extra copies out of the input directory."
log "no base name clash; all ${#FQ[@]} inputs will produce a separate output"

run_one() {  # $1 = the path of a fastq
  local fq="$1" base gz=() extra=()
  base=$(basename "$fq"); base=${base%.gz}; base=${base%.fastq}; base=${base%.fq}
  # The confidence threshold goes into the output name; otherwise a run at 0 and a
  # run at 0.1 write into the same file and which result came from which threshold
  # is lost.
  local ek=""; [ -n "$CONF" ] && ek="_c$CONF"
  local rep="$OUT/${base}${ek}_kraken2.report" outp="$OUT/${base}${ek}_output"
  if [ -s "$rep" ] && [ -s "$outp" ]; then echo "SKIPPED $base"; return 0; fi
  case "$fq" in *.gz) [ "$HAS_GZ" = 1 ] && gz=(--gzip-compressed);; esac
  [ "$USE_MMAP" = 1 ] && extra+=(--memory-mapping)
  [ "$HAS_MINHIT" = 1 ] && extra+=(--minimum-hit-groups 3)
  [ -n "$CONF" ] && extra+=(--confidence "$CONF")
  kraken2 --db "$DB" --threads "$THREADS" "${gz[@]}" "${extra[@]}" \
          --use-names --report "$rep" --output "$outp" "$fq"
}

# --- 5. the measurement run (on the smallest file) --------------------
SMALL=$(for f in "${FQ[@]}"; do printf '%s\t%s\n' "$(stat -c%s "$f")" "$f"; done | sort -n | head -1 | cut -f2-)
SMALLB=$(stat -c%s "$SMALL")
log "the file measured: $(basename "$SMALL") ($(awk -v b=$SMALLB 'BEGIN{printf "%.1f", b/1048576}') MB)"
T0=$(date +%s)
run_one "$SMALL"
T1=$(date +%s)
ELAPSED=$(( T1 - T0 )); [ "$ELAPSED" -lt 1 ] && ELAPSED=1
log "the measured time: ${ELAPSED} seconds"
awk -v e="$ELAPSED" -v s="$SMALLB" -v t="$TOTBYTES" 'BEGIN{
  r=e/s; tot=r*t;
  printf "[estimate] rough total time: %.1f minutes (%.2f hours)\n", tot/60, tot/3600;
  printf "[estimate] note: the first call includes loading the database,\n";
  printf "           and later files are faster when memory mapping is off.\n"}'

[ "$ONLYBENCH" = 1 ] && { log "only a measurement was asked for, stopping"; exit 0; }
if [ "$ASSUME_YES" != 1 ]; then
  read -r -p "Start the full run? [y/N] " a
  case "$a" in y|Y) ;; *) log "cancelled"; exit 0;; esac
fi

# --- 6. the full run --------------------------------------------------
i=0
for f in "${FQ[@]}"; do
  i=$(( i + 1 ))
  log "[$i/${#FQ[@]}] $(basename "$f")"
  run_one "$f"
done

# --- 7. the summary ---------------------------------------------------
log "the summary"
python3 - "$OUT" <<'PY'
import sys,glob,os
out=sys.argv[1]
print("%-34s %10s %9s %7s"%("report","total","unclass","unc%"))
for f in sorted(glob.glob(os.path.join(out,"*_kraken2.report"))):
    unc=root=0
    for line in open(f,errors="replace"):
        p=line.rstrip("\n").split("\t")
        if len(p)<6: continue
        if p[3]=="U": unc=int(p[1])
        if p[3]=="R" and p[4]=="1": root=int(p[1])
    tot=unc+root
    if tot: print("%-34s %10d %9d %6.2f"%(os.path.basename(f)[:34],tot,unc,100*unc/tot))
PY
log "finished. The output: $OUT"
echo
echo "The next step is the Bracken re-estimation. Before it, check the reports in"
echo "this directory and the environment report, and pick the kmer_distrib file"
echo "whose length matches your read length; Bracken needs that length to match."
