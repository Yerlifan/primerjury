#!/usr/bin/env bash
# =====================================================================
# qiime2_classify.sh: classify every library's features against the reference
#                     that fits it, with vsearch, and build a consensus taxonomy.
#
#   bash steps/qiime2_classify.sh --root . [--profile QIIME2_PROFILE] [--threads 3]
#        [--identity 0.80] [--refs "A1:ssu B:ssu A2:ssu,lsu F1:unite,lsu F2:unite,lsu"]
#
# WHY NOT THE READY-MADE CLASSIFIER: it is trained for 16S and gave meaningless
# labels on the fungal libraries and on the long archaeal operon. Each library is
# aligned to SILVA SSU, SILVA LSU and/or UNITE (REFERENCE_DB), and
# qiime2_reference_taxonomy.py takes the taxonomy the best hits agree on, with
# the alignment-length rule of the identity chain (one source of thresholds).
#
# --refs maps a library group to its reference set(s); when a group is not
# listed: F* -> unite,lsu; a group whose features average over 2,200 bases
# -> ssu,lsu; otherwise ssu.
# =====================================================================
set -uo pipefail
ROOT="."; PROFILE="QIIME2_PROFILE"; THREADS=3; IDENT=0.80; REFS_GIVEN=""
while [ $# -gt 0 ]; do
  case "$1" in
    --root) ROOT="$2"; shift 2;;
    --profile) PROFILE="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --identity) IDENT="$2"; shift 2;;
    --refs) REFS_GIVEN="$2"; shift 2;;
    -h|--help) sed -n 2,18p "$0"; exit 0;;
    *) echo "unknown option: $1" >&2; exit 2;;
  esac
done
ROOT="$(cd "$ROOT" && pwd)"; HERE="$(cd "$(dirname "$0")" && pwd)"
KOK="$ROOT/$PROFILE"; [ -d "$KOK" ] || { echo "ERROR: profile directory missing: $KOK" >&2; exit 1; }
REF="$ROOT/REFERENCE_DB"; LOG="$KOK/qiime2_classify.log"
log(){ printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" | tee -a "$LOG"; }
err(){ printf '[%s] ERROR: %s\n' "$(date +%H:%M:%S)" "$*" | tee -a "$LOG" >&2; }
command -v vsearch >/dev/null 2>&1 || { err "vsearch not found; activate the qiime2 environment"; exit 1; }
log "=== reference classification started ($(vsearch --version 2>&1 | head -1)) ==="

dna_copy(){ # $1 source (may carry RNA letters), $2 target
  if [ -s "$2" ]; then return 0; fi
  [ -s "$1" ] || { err "reference missing: $1"; return 1; }
  awk '/^>/{print;next}{gsub(/U/,"T");gsub(/u/,"t");print}' "$1" > "$2" && [ -s "$2" ]
}
SSU="$REF/SILVA_138.2_SSURef_NR99.fasta"; [ -s "$SSU" ] || SSU=""
LSU="$REF/SILVA_138.2_LSURef_NR99_DNA.fasta"; dna_copy "$REF/SILVA_138.2_LSURef_NR99.fasta" "$LSU" || LSU=""
UNI="$REF/UNITE_ITS.fasta"; [ -s "$UNI" ] || { err "UNITE not found"; UNI=""; }

ref_path(){ case "$1" in ssu) echo "$SSU";; lsu) echo "$LSU";; unite) echo "$UNI";; *) echo "";; esac; }
refs_for(){ # $1 group, $2 query fasta
  local g="$1" q="$2" hit=""
  for kv in $REFS_GIVEN; do [ "${kv%%:*}" = "$g" ] && hit="${kv#*:}"; done
  if [ -n "$hit" ]; then echo "$hit" | tr ',' ' '; return; fi
  case "$g" in F*) echo "unite lsu"; return;; esac
  local mean; mean=$(awk '/^>/{if(l){s+=l;n++} l=0;next}{l+=length($0)} END{if(l){s+=l;n++} printf "%d", (n?s/n:0)}' "$q")
  if [ "${mean:-0}" -gt 2200 ]; then echo "ssu lsu"; else echo "ssu"; fi
}
align(){ # $1 group, $2 query, $3 reference, $4 label
  local g="$1" q="$2" r="$3" et="$4" h="$KOK/$1/hits_$4.tsv"
  [ -n "$r" ] && [ -s "$r" ] || { err "$g: no $et reference, skipped"; return 1; }
  if [ -s "$h" ]; then log "$g $et alignment exists, skipped"; return 0; fi
  log "$g aligning to $et"
  vsearch --usearch_global "$q" --db "$r" --id "$IDENT" --maxaccepts 20 --maxrejects 100 \
          --strand both --threads "$THREADS" --blast6out "$h" --quiet >>"$LOG" 2>&1
  if [ ! -s "$h" ]; then
    log "$g $et: no hit at $IDENT, retrying at 0.70"
    vsearch --usearch_global "$q" --db "$r" --id 0.70 --maxaccepts 20 --maxrejects 100 \
            --strand both --threads "$THREADS" --blast6out "$h" --quiet >>"$LOG" 2>&1
  fi
  [ -s "$h" ] && log "$g $et hit rows: $(wc -l < "$h")" || err "$g $et alignment is empty"
}
for gd in "$KOK"/*/; do
  g=$(basename "$gd"); Q="$gd/export_min10/dna-sequences.fasta"
  [ -s "$Q" ] || continue
  log "--- library $g, $(grep -c '^>' "$Q") features ---"
  PAIRS=()
  for et in $(refs_for "$g" "$Q"); do
    r=$(ref_path "$et"); align "$g" "$Q" "$r" "$et" && PAIRS+=("$gd/hits_$et.tsv:$r")
  done
  [ ${#PAIRS[@]} -gt 0 ] || { err "$g: no alignment could be made"; continue; }
  python3 "$HERE/qiime2_reference_taxonomy.py" "$gd/taxonomy_reference.tsv" "${PAIRS[@]}" "fasta:$Q" >>"$LOG" 2>&1 \
    && log "$g taxonomy ready: $(( $(wc -l < "$gd/taxonomy_reference.tsv") - 1 )) features" \
    || err "$g taxonomy could not be built"
done
log "=== reference classification finished ==="
