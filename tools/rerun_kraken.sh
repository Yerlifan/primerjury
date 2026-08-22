#!/bin/bash
# ---------------------------------------------------------------------------
# RECLASSIFYING EVERY TAXON WITH KRAKEN2
#
# WHY ALL OF THEM
# Seven eukaryotes are known to have been mislabelled, and their identities were
# found by an indirect route. Reclassifying only those seven amounts to confirming
# what was already suspected, and it cannot find a NEW mislabel. A criterion lets
# through the mistakes in the places it does not look, and that was the common
# cause of the ten bugs found in this project so far. So EVERY taxon in the fastq
# directory is scanned and none is picked out.
#
# WHY A SINGLE RUN
# kraken2 loads the database from scratch on every call. One call per taxon means
# loading the same database dozens of times, and almost all the work goes into
# loading. Instead every read is merged into one file, each read name carries its
# source taxid, kraken2 runs ONCE, and the result is split by taxon afterwards.
# Since kraken2 classifies every read independently, merging does not change the
# result.
#
# Run:
#   bash tools/rerun_kraken.sh
#   bash tools/rerun_kraken.sh /another/k2db/path
#   IPLIK=8 bash tools/rerun_kraken.sh
#   PROJE=/full/path/to/project bash tools/rerun_kraken.sh
# ---------------------------------------------------------------------------
set -euo pipefail

# THE PROJECT ROOT IS FOUND, NOT HARD CODED.
# This script sits in tools/, so the root is one directory up. BASH_SOURCE is the
# real path of the running file, which makes this a measurement and not a guess.
# If the checkout was moved it can be overridden:  PROJE=/full/path bash <script>
_BETIK_DIZIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJE="${PROJE:-$(cd "$_BETIK_DIZIN/.." && pwd)}"
if [ ! -d "$PROJE/tools" ] || [ ! -d "$PROJE/verification" ]; then
  echo "ERROR: the project root could not be verified ('$PROJE' holds no tools/ and verification/)."
  echo "  Script location: $_BETIK_DIZIN"
  echo "  To give it by hand:  PROJE=/full/path/to/project bash $0"
  exit 1
fi
# THE KRAKEN2 DATABASE PATH IS NOT HARD CODED EITHER.
# It is found by looking for the marker file hash.k2d. A first argument wins.
DB="${1:-}"
if [ -z "$DB" ]; then
  for _k in "$HOME"/*/hash.k2d "$HOME"/*/*/hash.k2d /home/*/*/hash.k2d \
            /opt/*/hash.k2d /mnt/c/*/hash.k2d; do
    [ -f "$_k" ] && DB="$(dirname "$_k")" && break
  done
fi
if [ -z "$DB" ] || [ ! -f "$DB/hash.k2d" ]; then
  echo "ERROR: no Kraken2 database was found (hash.k2d was searched for)."
  echo "  To give it by hand:  bash $0 /full/path/to/k2db"
  exit 1
fi
echo "Kraken2 database found: $DB  (not hard coded, hash.k2d was searched for)"
IS="$PROJE/RESULTS/kraken_yeniden"
KAYNAK="$PROJE/RESULTS/fastq files"
IPLIK="${IPLIK:-12}"
BURASI="$_BETIK_DIZIN"

echo "database : $DB"
echo "source   : $KAYNAK"
echo "output   : $IS"
echo "threads  : $IPLIK"
echo

# --- environment. kraken2 lives in the micromamba environment "mikro", not on PATH
# install.sh tools installs it with micromamba, and in a new shell the environment
# is not active, so kraken2 is invisible. The same thing is done here. A different
# environment name can be given with the ORTAM variable.
ORTAM="${ORTAM:-mikro}"
if ! command -v kraken2 >/dev/null 2>&1; then
  echo "kraken2 is not on PATH, activating the micromamba environment '$ORTAM'"
  export PATH="$HOME/bin:$PATH"
  export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
  if command -v micromamba >/dev/null 2>&1; then
    eval "$(micromamba shell hook -s bash)" || true
    micromamba activate "$ORTAM" || true
  elif command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)" || true
    conda activate "$ORTAM" || true
  fi
fi

# --- pre-checks. Stopping here is preferred over running silently incomplete ---
if ! command -v kraken2 >/dev/null 2>&1; then
  echo "ERROR: kraken2 still could not be found."
  echo "  Activate the environment by hand and try again:"
  echo "    export PATH=\"\$HOME/bin:\$PATH\""
  echo "    export MAMBA_ROOT_PREFIX=\"\$HOME/micromamba\""
  echo "    eval \"\$(micromamba shell hook -s bash)\""
  echo "    micromamba activate $ORTAM"
  echo "  If the environment has a different name: micromamba env list"
  exit 1
fi
echo "kraken2: $(command -v kraken2)"
[ -f "$DB/hash.k2d" ] || { echo "ERROR: $DB/hash.k2d is missing, the database path is wrong"; exit 1; }
[ -d "$KAYNAK" ]      || { echo "ERROR: no fastq directory: $KAYNAK"; exit 1; }
[ -f "$BURASI/kraken_summary.py" ] || { echo "ERROR: tools/kraken_summary.py is missing"; exit 1; }

mkdir -p "$IS"
cd "$IS"

DOSYALAR=$(ls "$KAYNAK"/*/*reads_*.fastq 2>/dev/null || true)
DOSYA_SAYISI=$(printf '%s\n' "$DOSYALAR" | grep -c . || true)
[ "$DOSYA_SAYISI" -gt 0 ] || { echo "ERROR: no fastq file was found"; exit 1; }

# --- memory. If the database does not fit in RAM it is read from disk ---------
# By default kraken2 loads the WHOLE hash table into RAM. A full database such as
# PlusPF passes 100 GB, which is impossible on a 16 GB machine, and kraken2 stops
# with "unable to allocate hash table memory". With --memory-mapping the table is
# read from disk, the RAM requirement goes away, and the work gets slower in
# return. The result is THE SAME; the classification does not change, only the
# time does.
DB_BAYT=$(du -sb "$DB" 2>/dev/null | cut -f1 || echo 0)
RAM_KB=$(awk '/MemAvailable/{print $2}' /proc/meminfo 2>/dev/null || echo 0)
RAM_BAYT=$((RAM_KB * 1024))
echo "database size  : $(numfmt --to=iec "$DB_BAYT" 2>/dev/null || echo "$DB_BAYT") "
echo "available RAM  : $(numfmt --to=iec "$RAM_BAYT" 2>/dev/null || echo "$RAM_BAYT")"
BELLEK_BAYRAGI=""
if [ "$DB_BAYT" -gt 0 ] && [ "$DB_BAYT" -gt "$RAM_BAYT" ]; then
  BELLEK_BAYRAGI="--memory-mapping"
  echo "the database DOES NOT FIT in RAM, --memory-mapping will be used (from disk, slow)"
else
  echo "the database fits in RAM, normal mode"
fi

# --- 1. every read in one file, with its source taxid written into the name ---
# If the KAP variable is set, at most that many reads are taken from each taxon.
# This is NOT a selection of taxa: all 44 taxa are still scanned, only the reads
# within each bin are sampled. A few thousand reads are enough to establish a
# bin's identity, and since the time under --memory-mapping is proportional to the
# read count, this saves hours on the first pass. With no KAP, every read is used.
KAP="${KAP:-0}"
[ "$KAP" -gt 0 ] && echo "at most $KAP reads will be taken from each taxon (every taxon is still scanned)"
# The read name becomes "@tx<taxid>_<old name>". The second field of the kraken2
# output is that name, and it is what the result is split by.
echo "merging $DOSYA_SAYISI fastq files"
: > tum.fastq
: > kaynak_sayim.tsv
: > /tmp/kraken_alinan.tsv
while IFS= read -r f; do
  [ -n "$f" ] || continue
  TX=$(basename "$f" | sed -E 's/.*reads_([0-9]+)\.fastq$/\1/')
  N=$(awk 'END{print int(NR/4)}' "$f")
  # KAP is a total per TAXON, not per file; the reads of one taxon can be spread
  # over more than one sample file.
  ALINAN=$(awk -F'\t' -v t="$TX" '$1==t{s+=$2} END{print s+0}' /tmp/kraken_alinan.tsv)
  if [ "$KAP" -gt 0 ]; then
    KALAN=$((KAP - ALINAN)); [ "$KALAN" -le 0 ] && KALAN=0
    BU=$(( N < KALAN ? N : KALAN ))
  else
    BU=$N
  fi
  printf '%s\t%s\t%s\t%s\n' "$TX" "$N" "$BU" "$(basename "$f")" >> kaynak_sayim.tsv
  printf '%s\t%s\n' "$TX" "$BU" >> /tmp/kraken_alinan.tsv
  [ "$BU" -gt 0 ] || continue
  awk -v tx="$TX" -v lim="$BU" \
    'NR%4==1 {k++; if (k>lim) exit; sub(/^@/, "@tx" tx "_")} {print}' "$f" >> tum.fastq
done <<< "$DOSYALAR"

TOPLAM=$(awk 'END{print int(NR/4)}' tum.fastq)
TAKSON=$(cut -f1 kaynak_sayim.tsv | sort -u | grep -c . || true)
echo "$TAKSON taxa, $TOPLAM reads merged"

# That the merge was lossless is verified separately. A silent loss is the classic
# way "not measured" gets read as "clean" in this project.
BEKLENEN=$(awk -F'\t' '{s+=$3} END{print s+0}' kaynak_sayim.tsv)
KAYNAKTA=$(awk -F'\t' '{s+=$2} END{print s+0}' kaynak_sayim.tsv)
[ "$TOPLAM" = "$BEKLENEN" ] || { echo "ERROR: the read count does not add up ($TOPLAM / $BEKLENEN)"; exit 1; }
echo "read count verified: $TOPLAM (of $KAYNAKTA in the source)"
# Which taxa were sampled is written out explicitly. "Every read was scanned" and
# "a sample was scanned" must not end in the same sentence.
awk -F'\t' '$2!=$3 {print "    " $1 ": " $3 " taken out of " $2 " reads"}' kaynak_sayim.tsv | sort -u | head -50

# --- 2. kraken2, a single run -------------------------------------------------
echo
echo "kraken2 is running. The database is loaded once, and that is the first wait."
[ -n "$BELLEK_BAYRAGI" ] && echo "  --memory-mapping is on, so this can take a long time depending on disk speed"
kraken2 --db "$DB" --threads "$IPLIK" $BELLEK_BAYRAGI \
        --output tum.out --report tum.report \
        --use-names tum.fastq
echo "kraken2 finished"
[ -s tum.out ] || { echo "ERROR: tum.out is empty, kraken2 produced no result"; exit 1; }

# --- 3. splitting and summarising, both on the python side --------------------
# The split is done in python rather than awk: the number of files awk can hold
# open at once is limited in some versions, and past that limit lines are lost
# without an error.
echo
python3 "$BURASI/kraken_summary.py" --job "$IS" --toolkit "$BURASI"

echo
echo "done. Files: $IS"
echo "  tum.report        the kraken report for the whole sample"
echo "  <taxid>.out       the read by read classification of each taxon"
echo "  kraken_ozet.csv   identities per taxon and whether they agree with the label"
