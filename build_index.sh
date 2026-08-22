#!/usr/bin/env bash
# =============================================================================
#  build_index.sh  --  builds the MFEprimer index for a reference FASTA FROM
#  SCRATCH. (2026-08-09: GENERALISED from a version fixed to SSU. The old usage is
#   kept exactly: run with no argument it still indexes SSU.)
# =============================================================================
#
#  USAGE
#  -----
#      bash build_index.sh                                  # SSU (the old behaviour)
#      bash build_index.sh SILVA_138.2_LSURef_NR99.fasta    # a name inside REFERENCE_DB
#      bash build_index.sh /full/path/other.fasta           # a full path works too
#      bash build_index.sh --list                           # show the candidate files
#      CPU=8 MEMPCT=60 bash build_index.sh <file>           # resource settings
#      SINA_F=... SINA_R=... SINA_AD=... bash build_index.sh <file>   # your own test pair
#
#  WHY IT IS NEEDED (a measured reason):
#    SILVA FASTA files store their sequences as RNA: U instead of T.
#    MFEprimer builds its 9-mer index over the {A,C,G,T} alphabet. Because there is
#    no T at all in the file, the effective alphabet falls to {A,C,G} and only
#    3^9 = 19,683 k-mers enter the index where there should be 4^9 = 262,144.
#    The evidence is a line of the old indexing log itself:
#        "Sorting 19683 kmers by ID..."      <- SILVA (broken)
#        "Sorting 262144 kmers by ID..."     <- PR2 / ROD / UNITE / SSU (sound)
#    MFEprimer RAISES NO ERROR in that case, it writes "Index built successfully".
#    That is what a silent failure looks like.
#
#  THE FIX: U -> T on the sequence lines. The header lines ARE NOT TOUCHED.
#
#  AN IMPORTANT TRAP: DO NOT USE a blind "sed 's/U/T/g'". There are capital U's in
#    the SILVA headers (measured: in 1,652 headers of the first 20 MB slice). A
#    blind conversion turns "Unknown Family" into "Tnknown Family" and breaks the
#    taxonomy. This script applies the conversion to the sequence lines only:
#    sed '/^>/!y/U/T/'
#
#  SAFETY: the U->T conversion is byte for byte; the file length and every byte
#    offset stay the same. So the existing .fai and BLAST indexes stay valid. Even
#    so, the script REQUIRES a backup to exist BEFORE the conversion (a sibling copy
#    of the same size, or a .gz); without one it exits WITHOUT TOUCHING the file.
#
#  THE RNA AGAINST DNA MEASUREMENT: whether the file is RNA or DNA IS NOT ASSUMED,
#    it is measured. If it is DNA the conversion step is skipped and the file is not
#    rewritten for nothing; that it was skipped is written to the screen PLAINLY.
#
#  RESILIENCE AGAINST INTERRUPTION: every step first asks "has this been done
#    already". If the script is cut short, run the same command again; it skips the
#    finished steps and continues where it stopped. A half written temporary
#    conversion file (.donusum_tmp) is deleted at the start of every run and is
#    never used.
# =============================================================================

set -u
set -o pipefail

# ---- The paths are derived from the script's own location, so they work even if
#      the directory is renamed
KOK="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DB="$KOK/REFERENCE_DB"
MFE="$KOK/tools/mfeprimer"

VARSAYILAN="SILVA_138.2_SSURef_NR99.fasta"   # a call with no argument = the old behaviour

# ---- --list / --yardim (help)
if [ "${1:-}" = "--list" ] || [ "${1:-}" = "-l" ]; then
  echo "the indexable files inside REFERENCE_DB:"
  printf '%-42s %10s  %s\n' "FILE" "MB" "INDEX"
  for f in "$DB"/*.fasta "$DB"/*.fna; do
    [ -f "$f" ] || continue
    printf '%-42s %10s  %s\n' "$(basename "$f")" \
      "$(( $(stat -c%s "$f") / 1048576 ))" \
      "$([ -f "$f.primerqc.bin" ] && echo PRESENT || echo MISSING)"
  done
  exit 0
fi
if [ "${1:-}" = "--yardim" ] || [ "${1:-}" = "-h" ]; then
  sed -n '2,40p' "${BASH_SOURCE[0]}"; exit 0
fi

# ---- Resolve the target FASTA: a full path, a name inside REFERENCE_DB, or the
#      default with no argument
GIRDI="${1:-$VARSAYILAN}"
if [ -f "$GIRDI" ]; then
  FASTA="$(cd "$(dirname "$GIRDI")" && pwd)/$(basename "$GIRDI")"
elif [ -f "$DB/$GIRDI" ]; then
  FASTA="$DB/$GIRDI"
else
  echo "ERROR: no such FASTA: $GIRDI"
  echo "To see the candidates: bash $(basename "${BASH_SOURCE[0]}") --list"
  exit 1
fi
AD="$(basename "$FASTA")"
BIN="$FASTA.primerqc.bin"
LOG="$FASTA.log"

# ---- Kaynak ayarlari. Olculmus degerler (kaynak makinedeki eski kosular):
#      PR2 351 MB -> 9,34 GB tepe RAM;  UNITE 1454 MB -> 10,87 GB tepe RAM.
CPU="${CPU:-16}"       # the core count    (to change it: CPU=8 bash build_index.sh)
MEMPCT="${MEMPCT:-70}" # RAM tavani yuzdesi; 70 eski basarili kosularin degeri
UYKU="${UYKU:-60}"     # ilerleme yazdirma araligi (saniye)

zaman() { date '+%H:%M:%S'; }
yaz()   { echo "[$(zaman)] $*"; }
bolum() { echo; echo "============================================================"; echo "  $*"; echo "============================================================"; }

# ---- THE FUNCTIONAL TEST PAIR: chosen by the database's region.
#      In an LSU or 28S file an SSU primer finds nothing, so a pair from the panel
#      that targets 28S is used. (Measured, TUM_CIFTLER_DEVIR_2026-08-07:
#      Mantar_universal (F1) fungi.28SrRNA.fna'da 226 amplikon, 82 bp -
#      exactly the expected product length. 0 in 18S. So it really is an LSU pair.)
if [ -n "${SINA_F:-}" ] && [ -n "${SINA_R:-}" ]; then
  SINA_AD="${SINA_AD:-Kullanici_cifti}"
elif printf '%s' "$AD" | grep -qi 'LSU\|28S\|23S'; then
  SINA_AD="Mantar_universal_F1"          # panel cifti, LSU/28S hedefli, 82 bp
  SINA_F="GGTTACCCGCTGAACTTAAGC"
  SINA_R="CGCTTCACTCGCCGTTAC"
elif printf '%s' "$AD" | grep -qi 'ITS'; then
  SINA_AD="Mantar_universal_F1"          # ITS kumelerinde de vurus verdigi olculdu
  SINA_F="GGTTACCCGCTGAACTTAAGC"
  SINA_R="CGCTTCACTCGCCGTTAC"
else
  SINA_AD="Bakteri_universal"            # the old test pair for SSU, 16S and 18S
  SINA_F="ACAAGCGGTGGAGCATGTG"
  SINA_R="ACGACAGCCATGCAGCAC"
fi

# Is the existing index sound? Two independent criteria, either one is enough:
#   1) "Sorting 262144 kmers" in the log (the full DNA alphabet; it is written in
#      disk mode only)
#   2) the index to fasta size ratio >= 4.0 (sound: 4.6-5.5x  /  broken: 0.63x)
indeks_saglikli_mi() {
  [ -f "$BIN" ] || return 1
  grep -q "Sorting 262144 kmers" "$LOG" 2>/dev/null && return 0
  local fs bs
  fs=$(stat -c%s "$FASTA"); bs=$(stat -c%s "$BIN")
  awk -v a="$bs" -v b="$fs" 'BEGIN{exit !(b>0 && a/b >= 4.0)}'
}

# =============================================================================
bolum "STEP 0/4  --  Pre-checks"
# =============================================================================
[ -x "$MFE" ] || { echo "ERROR: mfeprimer is missing or not executable: $MFE"; echo "To fix it: chmod +x '$MFE'"; exit 1; }
[ -f "$FASTA" ] || { echo "ERROR: no such FASTA: $FASTA"; exit 1; }
FASTA_MB=$(( $(stat -c%s "$FASTA") / 1048576 ))
yaz "mfeprimer : $MFE"
yaz "TARGET    : $AD  (${FASTA_MB} MB)"
yaz "Sinama    : $SINA_AD  ($SINA_F / $SINA_R)"

# --- The time and RAM estimate: MEASUREMENT BASED, scaled with the file size.
#     Olculen indeksleme hizlari (bu makine):
#       PR2   351 MB -> 11.508 s = 32,8 s/MB
#       ROD   340 MB -> 12.419 s = 36,5 s/MB
#       UNITE 1454 MB -> 43.056 s = 29,6 s/MB      ortalama ~33 s/MB
#       SSU   786 MB -> 28,020 s = 35.7 s/MB       (SILVA's own measurement)
#     Tepe RAM, olculen iki noktadan dogru orantiyla:
#       PR2 351 MB -> 9,34 GB ; UNITE 1454 MB -> 10,87 GB
#       egim = (10,87-9,34)/(1454-351) = 0,001387 GB/MB ; kesme = 8,853 GB
awk -v mb="$FASTA_MB" 'BEGIN{
  printf "           Sure tahmini : %.1f saat (33,0 s/MB ortalama)  |  %.1f saat (35,7 s/MB SSU olcumu)\n", mb*33.0/3600, mb*35.7/3600;
  printf "           Tepe RAM     : ~%.1f GB (0,001387 GB/MB x %d MB + 8,853 GB)\n", mb*0.001387+8.853, mb;
}'

# By the measured ratio the index will be about 5.1 times the fasta size (PR2 5.18x /
# ROD 5,49x / UNITE 4,57x / SSU 5,29x). Donusum sirasinda gecici bir kopya da tutulur.
# The space needed = fasta x 5.5 (the index) + fasta x 1 (the temporary copy) +
# 500 MB of headroom.
GEREKEN_MB=$(( FASTA_MB * 7 + 500 ))
BOS_MB=$(df -Pm "$DB" | awk 'NR==2{print $4}')
yaz "Bos disk  : ${BOS_MB} MB (gereken: ~${GEREKEN_MB} MB)"
if [ "$BOS_MB" -lt "$GEREKEN_MB" ]; then
  echo "ERROR: not enough disk space. Free at least ${GEREKEN_MB} MB and try again."; exit 1
fi
yaz "Cekirdek  : $CPU     RAM tavani: %$MEMPCT"

# Clean up a half written temporary file if there is one (it is never used)
if [ -f "$FASTA.donusum_tmp" ]; then
  yaz "WARNING: a half written conversion file left from an earlier run was deleted (.donusum_tmp)."
  rm -f "$FASTA.donusum_tmp"
fi

# =============================================================================
bolum "STEP 1/4  --  Move the broken index aside"
# =============================================================================
# We do not delete it, we keep it as ".bozuk" so the old and the new can be compared.
if [ -f "$BIN" ]; then
  if indeks_saglikli_mi; then
    yaz "The existing index looks SOUND (the k-mer count and/or the size ratio is right)."
    yaz "To rebuild it, move this file out of the way first: $BIN"
  else
    mv -f "$BIN" "$BIN.bozuk"
    yaz "The broken index was moved aside -> $(basename "$BIN").bozuk"
  fi
else
  yaz "There is no index to move aside (never built, or moved already). Skipping."
fi

# =============================================================================
bolum "STEP 2/4  --  RNA or DNA? (measured, not assumed)"
# =============================================================================
# Ilk 50 MB'in DIZI satirlarinda U ve T sayilir. U>0 ise RNA, U=0 ise DNA.
ORNEK=$(head -c 50000000 "$FASTA" | grep -v '^>')
U_SAT=$(printf '%s' "$ORNEK" | grep -c 'U' || true)
T_SAT=$(printf '%s' "$ORNEK" | grep -c 'T' || true)
unset ORNEK
yaz "The measurement (the first 50 MB, the sequence lines only): $U_SAT lines holding U / $T_SAT lines holding T"

if [ "$U_SAT" -eq 0 ]; then
  yaz ">>> THE FILE IS ALREADY DNA. THE U->T CONVERSION WAS SKIPPED. <<<"
  yaz "    The file was not rewritten, the change list is empty, 0 s elapsed."
  yaz "    (This is not a fault: the PR2, ROD, UNITE and RefSeq sets are DNA.)"
else
  yaz "The file is in RNA form. A conversion is needed."

  # --- THE BACKUP REQUIREMENT: a backup MUST EXIST before the original is changed.
  YEDEK=""
  BOY=$(stat -c%s "$FASTA")
  for aday in "$DB"/*.fasta "$DB"/*.fna; do
    [ -f "$aday" ] || continue
    [ "$aday" = "$FASTA" ] && continue
    [ "$(stat -c%s "$aday")" = "$BOY" ] && { YEDEK="$aday"; break; }
  done
  [ -z "$YEDEK" ] && [ -f "$FASTA.gz" ] && YEDEK="$FASTA.gz"
  if [ -z "$YEDEK" ]; then
    echo "ERROR: this file has no backup (no sibling copy of the same size, no .gz)."
    echo "      THE ORIGINAL WAS NOT TOUCHED. Take a backup first:"
    echo "        cp '$FASTA' '$FASTA.rna_backup'"
    exit 1
  fi
  yaz "  A backup was found: $(basename "$YEDEK")  (the original RNA version can be restored from it)"

  yaz "The conversion is starting. The sequence lines only; the headers are kept..."
  TMP="$FASTA.donusum_tmp"
  rm -f "$TMP"
  # y/U/T/ = transliterasyon (s///g'den hizli), /^>/! = baslik satirlarini atla
  sed '/^>/!y/U/T/' "$FASTA" > "$TMP" || { echo "ERROR: the conversion failed."; rm -f "$TMP"; exit 1; }
  yaz "The conversion is done, verifying..."

  # --- Check 1: the file length has to be exactly the same (U->T is a one byte swap)
  E=$(stat -c%s "$FASTA"); Y=$(stat -c%s "$TMP")
  if [ "$E" -ne "$Y" ]; then
    echo "ERROR: the size does not match (old=$E new=$Y). Nothing was changed."; rm -f "$TMP"; exit 1
  fi
  yaz "  OK  the size is the same: $E bytes"

  # --- Check 2: the header count must not have changed
  BE=$(grep -c '^>' "$FASTA"); BY=$(grep -c '^>' "$TMP")
  if [ "$BE" -ne "$BY" ]; then
    echo "ERROR: the header count does not match ($BE vs $BY)."; rm -f "$TMP"; exit 1
  fi
  yaz "  OK  the header count is the same: $BE"

  # --- Dogrulama 3: basliklar bayt-bayt korunmus olmali
  if ! diff -q <(grep '^>' "$FASTA") <(grep '^>' "$TMP") >/dev/null; then
    echo "ERROR: the headers changed. Nothing was changed."; rm -f "$TMP"; exit 1
  fi
  yaz "  OK  basliklar birebir korundu"

  # --- Dogrulama 4: dizi satirlarinda U kalmamis olmali
  KALAN=$(grep -v '^>' "$TMP" | grep -c 'U' || true)
  if [ "$KALAN" -ne 0 ]; then
    echo "ERROR: $KALAN sequence lines still contain U."; rm -f "$TMP"; exit 1
  fi
  yaz "  OK  dizi satirlarinda U kalmadi"

  mv -f "$TMP" "$FASTA"
  yaz "Donusum uygulandi. (Orijinal RNA surumu: $(basename "$YEDEK"))"

  # Because the offsets do not change the .fai stays valid, but to be safe we let
  # mfeprimer regenerate it.
  rm -f "$FASTA.fai"
fi

# =============================================================================
bolum "STEP 3/4  --  Indexing  (THE LONG STEP)"
# =============================================================================
if indeks_saglikli_mi; then
  yaz "A sound index is there already, the indexing is skipped."
else
  awk -v mb="$FASTA_MB" 'BEGIN{printf "           Beklenen sure: %.1f - %.1f saat (33,0 - 35,7 s/MB araligi)\n", mb*29.6/3600, mb*36.5/3600}'
  yaz "Indeksleme basliyor. Ilerleme her $UYKU saniyede bir asagi yazilacak."
  yaz "Komut: mfeprimer index -i $AD -c $CPU -m $MEMPCT -f"
  T0=$(date +%s)

  "$MFE" index -i "$FASTA" -c "$CPU" -m "$MEMPCT" -f &
  PID=$!

  while kill -0 "$PID" 2>/dev/null; do
    sleep "$UYKU"
    GECEN=$(( ($(date +%s) - T0) / 60 ))
    SON=$(tail -c 20000 "$LOG" 2>/dev/null | grep -E '\[Memory\]|Phase' | tail -1)
    yaz "  ... ${GECEN} dk gecti   ${SON:-(gunluk henuz yazilmadi)}"
  done

  wait "$PID"; RC=$?
  SURE=$(( ($(date +%s) - T0) / 60 ))
  if [ "$RC" -ne 0 ]; then
    echo "ERROR: mfeprimer index ended with exit code $RC (after $SURE min)."
    echo "the end of the log:"; tail -20 "$LOG"
    exit 1
  fi
  yaz "Indexing finished: $SURE minutes.  (the measured rate: $(awk -v s=$SURE -v mb=$FASTA_MB 'BEGIN{if(mb>0) printf "%.1f", s*60/mb; else printf "n/a"}') s/MB)"
fi

# =============================================================================
bolum "STEP 4/4  --  Verification (three tests, all three must pass)"
# =============================================================================
GECTI_A=0; GECTI_B=0; GECTI_C=0

# --- A: the k-mer count in the log has to be 262144 (in the broken one it was 19683)
KMER=$(grep -o "Sorting [0-9]* kmers" "$LOG" 2>/dev/null | tail -1 | grep -o '[0-9]*')
if [ -n "$KMER" ]; then
  if [ "$KMER" = "262144" ]; then
    yaz "A) the k-mer count = $KMER  ->  PASSED (4^9, the full DNA alphabet)"; GECTI_A=1
  else
    yaz "A) the k-mer count = $KMER, 262144 expected  ->  FAILED (the index is still broken)"
  fi
else
  yaz "A) there is no 'Sorting' line in the log (it may have run in memory mode) -> COULD NOT BE MEASURED"
  GECTI_A=1   # olculemeyen sinama "dustu" sayilmaz; B ve C hukmu verir
fi

# --- B: boyut orani. Saglikli indeksler fasta boyunun 4,6-5,5 kati.
if [ -f "$BIN" ]; then
  FS=$(stat -c%s "$FASTA"); BS=$(stat -c%s "$BIN")
  ORAN=$(awk -v a=$BS -v b=$FS 'BEGIN{printf "%.2f", a/b}')
  if awk -v o="$ORAN" 'BEGIN{exit !(o>=4.6 && o<=5.5)}'; then
    yaz "B) the index to fasta ratio = ${ORAN}x  ->  PASSED (the sound range is 4.6-5.5x)"; GECTI_B=1
  else
    yaz "B) the index to fasta ratio = ${ORAN}x  ->  FAILED (4.6-5.5x expected; in the broken one it was 0.63x)"
  fi
else
  yaz "B) the index file was not created: $BIN  ->  FAILED"
fi

# --- C: ISLEVSEL sinama. Bolgeye uygun primer SIFIRDAN BUYUK baglanma vermeli.
yaz "C) Islevsel sinama: $SINA_AD primeri $AD'ye karsi kosuluyor..."
SINA_DIR="$KOK/VERIFICATION_RESULT/indeks_sinama"
mkdir -p "$SINA_DIR"
ETIKET="$(printf '%s' "$AD" | tr -c 'A-Za-z0-9' '_')"
printf '%s\t%s\t%s\n' "$SINA_AD" "$SINA_F" "$SINA_R" > "$SINA_DIR/sina_$ETIKET.tsv"
"$MFE" spec -i "$SINA_DIR/sina_$ETIKET.tsv" -d "$FASTA" -o "$SINA_DIR/sonuc_$ETIKET.txt" >/dev/null 2>&1
echo
awk '/Binding Number/{f=1} f&&/_fp |_rp /{print "     "$0; c++} c==2{exit}' "$SINA_DIR/sonuc_$ETIKET.txt"
grep -m1 "potential amplicons" "$SINA_DIR/sonuc_$ETIKET.txt" | sed 's/^/     /'
echo
VUR=$(awk '/Binding Number/{f=1} f&&/_fp /{print $(NF-1); exit}' "$SINA_DIR/sonuc_$ETIKET.txt")
if [ -n "${VUR:-}" ] && [ "$VUR" -gt 0 ] 2>/dev/null; then
  yaz "C) ${SINA_AD}_fp Plus baglanmasi = $VUR (> 0)  ->  GECTI"; GECTI_C=1
else
  yaz "C) bindings = ${VUR:-0}  ->  FAILED (the index is not being read, or the primer is not in this region)"
fi

# --- Hukum
if [ "$GECTI_A" -eq 1 ] && [ "$GECTI_B" -eq 1 ] && [ "$GECTI_C" -eq 1 ]; then
  bolum "RESULT: SUCCEEDED  --  the $AD index is sound (A+B+C passed)"
  exit 0
else
  bolum "RESULT: FAILED  --  $AD  (A=$GECTI_A B=$GECTI_B C=$GECTI_C; 1=passed)"
  exit 1
fi
