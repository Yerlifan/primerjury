#!/usr/bin/env bash
# =====================================================================
# heavy_jobs.sh
# It runs, in order, the jobs that do not run safely in a small container. Each
# step can also be run on its own; the only reason they go in order is that the
# outputs feed one another.
#
#   bash heavy_jobs.sh                # steps A, B, C, D, E
#   bash heavy_jobs.sh --kraken       # adds step F as well (minutes)
#   SECILEN_ESIK=0.02 bash heavy_jobs.sh --kraken   # apply the threshold directly
#   bash heavy_jobs.sh --only C     # one step only
#
# The steps:
#   A  bin recovery            builds the B-1_2233851 consensus from the reads
#   B  target identity         puts each target's name and what the sequence
#                              shows side by side
#   C  mfeprimer               the second independent measurement of external
#                              specificity
#   D  community trend         builds the abundance workbook on a genus level base
#   E  the self audit          the regression and delivery audits
#   F  the Kraken2 confidence  from the output files; no database is needed
#      threshold               (the threshold scan first, then the chosen threshold)
#   G  the Excel delivery      with the measured identity columns, after B
#   H  the broad external      SILVA, UNITE, ROD, PR2; it takes HOURS and runs
#      database scan           only with "--only H"
#
# Every step writes its own log under $PT/agir_log and none of them deletes
# another's output. If a step fails the script DOES NOT STOP, it marks the failure
# and carries on; at the end it lists which steps failed.
# =====================================================================
set -uo pipefail

# The project directory. It used to default to one person's desktop path, so on
# any other machine the defaults below pointed at nothing. It is derived from
# this script's own location now; $PT still overrides it.
_BETIK_DIZIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PT="${PT:-$(cd "$_BETIK_DIZIN/.." && pwd)}"
HERE="$(cd "$(dirname "$0")" && pwd)"
KONS="${KONS:-$PT/referans_konsensus/baskin/konsensus}"
ADAY="${ADAY:-$PT/primer_adaylari}"
FINAL="${FINAL:-$PT/primer_final}"
LOGD="$PT/agir_log"
MFE="${MFE:-$PT/tools/mfeprimer}"
K2DB="${K2DB:-$HOME/k2db}"
HAM_FASTQ="${HAM_FASTQ:-}"
CONF="${CONF:-0.1}"
IS="${IS:-4}"

KRAKEN=0; YALNIZ=""
while [ $# -gt 0 ]; do
  case "$1" in
    --kraken) KRAKEN=1; shift;;
    --only) YALNIZ="$2"; shift 2;;
    --pt) PT="$2"; shift 2;;
    *) echo "unknown option: $1" >&2; exit 2;;
  esac
done

mkdir -p "$LOGD"
ANA="$LOGD/AGIR_ISLER.log"
ts(){ date '+%Y-%m-%d %H:%M:%S'; }
say(){ printf '[%s] %s\n' "$(ts)" "$*" | tee -a "$ANA"; }
BASARISIZ=()

calistir() {   # $1 = the step letter, $2 = a description, the rest is the command
  local h="$1" ac="$2"; shift 2
  if [ -n "$YALNIZ" ] && [ "$YALNIZ" != "$h" ]; then return 0; fi
  say "----------------------------------------------------------------"
  say "STEP $h  $ac"
  local t0=$(date +%s)
  "$@" 2>&1 | tee -a "$LOGD/adim_$h.log" | tail -n 40
  local rc=${PIPESTATUS[0]}
  say "STEP $h finished, exit=$rc, time=$(( ($(date +%s)-t0)/60 )) minutes"
  if [ "$rc" -ne 0 ]; then
    BASARISIZ+=("$h ($ac), cikis=$rc")
    say "  CAREFUL: this step failed, the full log: $LOGD/adim_$h.log"
  fi
}

say "================================================================"
say "AGIR ISLER BASLANGIC"
say "PT     : $PT"
say "KONS   : $KONS"
say "loglar : $LOGD"

# --- bagimlilik denetimi ---------------------------------------------
say "bagimlilik denetimi"
EKSIK=()
for k in blastn python3; do command -v "$k" >/dev/null 2>&1 || EKSIK+=("$k"); done
[ -x "$MFE" ] || { chmod +x "$MFE" 2>/dev/null || true; }
[ -x "$MFE" ] || EKSIK+=("mfeprimer (not executable: $MFE)")
python3 - <<'PYX' 2>&1 | tee -a "$ANA"
for m in ("primer3", "Bio", "mappy", "openpyxl"):
    try:
        __import__(m); print("   %-10s TAMAM" % m)
    except ImportError:
        print("   %-10s EKSIK  -> pip install --break-system-packages %s"
              % (m, {"Bio": "biopython", "primer3": "primer3-py"}.get(m, m)))
PYX
if [ ${#EKSIK[@]} -gt 0 ]; then
  say "EKSIK: ${EKSIK[*]}"
  say "  for blastn: sudo apt-get install -y ncbi-blast+"
fi

# --- A. bin recovery --------------------------------------------------
# Because the B-1_2233851 self reference is 19 percent IUPAC, minimap2 can pull only
# 2 minimizers out of that sequence and the consensus comes out zero length. The
# reads are there (5914 of them); the template is built from the reads themselves.
KURT_FQ="$PT/fastq files/B-1/B-1-reads_2233851.fastq"
if [ -f "$KURT_FQ" ] && [ ! -s "$KONS/B-1_2233851_baskin_konsensus.fasta" ]; then
  calistir A "bin recovery: B-1_2233851" \
    python3 "$HERE/recover_bins.py" --fastq "$KURT_FQ" \
      --label B-1_2233851 --out "$KONS"
else
  say "STEP A was skipped (there is no fastq, or the consensus is there already)"
fi

# --- B. target identity -----------------------------------------------
calistir B "target identity: the name against what the sequence shows" \
  python3 "$HERE/target_identity.py" --consensus "$KONS" --db "$PT/REFERANS_DB" \
    --targets "$HERE/hedefler.tsv" --names "$HERE/taxid_adlari.tsv" \
    --threads "$IS" --out "$FINAL/hedef_kimlik.tsv"

# --- C. mfeprimer -----------------------------------------------------
# In a small container the bacteria.16S step ran out of memory. There is 16 GB
# here, but a timeout is still set per database and a database that goes over it is
# marked "olculemedi" rather than being counted clean silently.
calistir C "mfeprimer: the second measurement of external specificity" \
  python3 "$HERE/mfeprimer_layer.py" --final "$FINAL" --db "$PT/REFERANS_DB" \
    --mfe "$MFE" --cpu "$IS" --timeout 7200 \
    --out "$FINAL/mfeprimer.tsv"

# --- D. topluluk trendi -----------------------------------------------
KIMLIK=()
for f in "$PT"/t_kimlik/kimlik_*.tsv "$HERE"/kimlik_*.tsv \
         "$HERE"/kimlik/kimlik_*.tsv; do
  [ -f "$f" ] && KIMLIK+=("$f")
done
if [ ${#KIMLIK[@]} -eq 0 ]; then
  say "WARNING: no kimlik_*.tsv was found. The species level reliability mark"
  say "  will be derived from ayirt_edilemez.tsv alone, and the bin identity"
  say "  measurement will be off. The expected location: $HERE/kimlik/"
fi
# The rank coverage is measured first: at which rank can the abundance be read.
# Bracken IS NOT RUN; the reasoning is at the head of abundance_rank.py and in
# bolluk_rutbe_kaniti.md.
RUTBEARG=""
if [ -d "$PT/kraken_c${SECILEN_ESIK:-0.02}" ]; then
  calistir D0 "rank coverage: at which rank can the abundance be read" \
    python3 "$HERE/abundance_rank.py" \
      --kraken "$PT/kraken_c${SECILEN_ESIK:-0.02}" \
      --out "$PT/bolluk_rutbe"
  [ -s "$PT/bolluk_rutbe/ozet.tsv" ] && RUTBEARG="--rank $PT/bolluk_rutbe"
else
  say "  NOTE: there is no kraken_c${SECILEN_ESIK:-0.02}, the rank coverage was skipped."
  say "       First: SECILEN_ESIK=0.02 bash heavy_jobs.sh --only F"
fi
calistir D "the community trend: a rank aware abundance workbook" \
  python3 "$HERE/community_trends.py" --bracken "$PT/bracken results" \
    --distinguishable "$ADAY/ayirt_edilemez.tsv" \
    ${KIMLIK[@]+--identity "${KIMLIK[@]}"} \
    --names "$HERE/taxid_adlari.tsv" $RUTBEARG \
    --out "$PT/PrimerJury_Community_Trends.xlsx"

# --- E. the self audit ------------------------------------------------
calistir E "the self audit: the regression suite" \
  python3 "$HERE/regression_test.py" --real-data --candidates "$ADAY" --consensus "$KONS"
# check_deliverables.py returns exit code 1 when it finds a CRITICAL item. That IS
# NOT A CRASH, it is A FINDING: the raw primer_final.tsv holds mixed domain rows and
# the Excel already leaves them out. Without separating the two, every run says "a
# failed step" and the real crashes go unnoticed.
say "----------------------------------------------------------------"
say "STEP E2  the self audit: the delivery audit"
T0E2=$(date +%s)
python3 "$HERE/check_deliverables.py" --final "$FINAL" --consensus "$KONS" \
  --targets "$HERE/hedefler.tsv" --out "$FINAL/teslim_denetimi.tsv" \
  2>&1 | tee -a "$LOGD/adim_E2.log" | tail -n 40
RCE2=${PIPESTATUS[0]}
say "STEP E2 finished, exit=$RCE2, time=$(( ($(date +%s)-T0E2)/60 )) minutes"
if [ "$RCE2" = 1 ]; then
  say "  E2: there is a CRITICAL finding (not a crash). The findings:"
  awk -F'\t' 'NR>1 && $1=="KRITIK"{say[$4]++} END{for(k in say) printf "     %-24s %d\n", k, say[k]}' \
    "$FINAL/teslim_denetimi.tsv" 2>/dev/null | tee -a "$ANA"
  say "  These findings are already left out of the workbook; they stay in the raw table."
elif [ "$RCE2" -ne 0 ]; then
  BASARISIZ+=("E2 (teslim denetimi cokmesi), cikis=$RCE2")
  say "  CAREFUL: this step CRASHED, the full log: $LOGD/adim_E2.log"
fi

# --- F. the Kraken2 confidence threshold -----------------------------
# IMPORTANT: reclassifying IS NOT NEEDED. Kraken2's --output files already carry the
# k-mer LCA sequence of every read, and the confidence score is computed from
# exactly that sequence. Neither the 106 GB database nor the raw barcode fastq files
# are required. Measured: 99.84 percent of the hit taxa are in the tree built from
# the report files.
#
# THE THRESHOLD IS NOT CHOSEN FROM MEMORY. The 0.1 value often recommended for short
# read data leaves 69 percent of the reads unclassified on this ONT data, because
# most of the k-mers of an ONT read find no counterpart in the database and those
# k-mers enter the denominator of the score. The scan runs first and the threshold is
# chosen from the table.
if [ "$KRAKEN" = 1 ] || [ "$YALNIZ" = "F" ]; then
  if [ ! -d "$PT/kraken results" ]; then
    say "STEP F was skipped: there is no 'kraken results' directory"
  else
    calistir F "the Kraken2 confidence threshold scan" \
      python3 "$HERE/reassign_confidence.py" --kraken "$PT/kraken results" \
        --scan 0,0.002,0.005,0.01,0.02,0.05,0.1 --scan-reads 20000 \
        --out "$PT/kraken_guven"
    say "  The scan table: $PT/kraken_guven/esik_taramasi.tsv"
    say "  Choose a threshold and run this command:"
    say "    python3 $HERE/reassign_confidence.py \\"
    say "        --kraken '$PT/kraken results' --confidence <ESIK> \\"
    say "        --out '$PT/kraken_c<ESIK>'"
    if [ -n "${SECILEN_ESIK:-}" ]; then
      calistir F2 "applying the Kraken2 confidence threshold, threshold=$SECILEN_ESIK" \
        python3 "$HERE/reassign_confidence.py" --kraken "$PT/kraken results" \
          --confidence "$SECILEN_ESIK" --out "$PT/kraken_c$SECILEN_ESIK"
    fi
  fi
fi

# --- G. refresh the Excel delivery ------------------------------------
# The hedef_kimlik.tsv produced in step B enters the workbook as the "olculen
# kimlik" column. That is why the Excel is reproduced AFTER B.
if [ -z "$YALNIZ" ] || [ "$YALNIZ" = "G" ]; then
  REFC="$PT/primer_referans"
  REFARG=""
  [ -s "$REFC/primer_referans.tsv" ] && REFARG="--reference $REFC/primer_referans.tsv"
  calistir G "the Excel delivery, with the measured identity columns" \
    python3 "$HERE/export_excel.py" \
      --candidates "$ADAY" --final "$FINAL" --splits "$ADAY/kume_setleri" \
      --names "$HERE/taxid_adlari.tsv" --targets "$HERE/hedefler.tsv" \
      --consensus "$KONS" --identity "$FINAL/hedef_kimlik.tsv" $REFARG \
      --out "$PT/PrimerJury_Primer_Tasarimi.xlsx"
fi

# --- H. the broad external database scan ------------------------------
# The narrow set (NCBI RefSeq 16S and ITS) holds only type strain and representative
# sequences; uncultured environmental lineages ARE NOT THERE. In this sample not
# even a relative above 90 percent could be found in RefSeq for three of the
# bacterial targets, so much of the community is not represented in the narrow set.
# SILVA SSU NR99, UNITE, ROD and PR2 carry environmental sequences too.
# IT TAKES A LONG TIME (hours), which is why it is a separate step started by hand.
if [ "$YALNIZ" = "H" ]; then
  calistir H "the broad external database scan (external_databases.py --wide)" \
    python3 "$HERE/external_databases.py" --final "$FINAL" \
      --db "$PT/REFERANS_DB" --wide --threads "$IS" --consensus "$KONS" \
      --targets "$HERE/hedefler.tsv" --names "$HERE/taxid_adlari.tsv" \
      --identity "$FINAL/hedef_kimlik.tsv" \
      --timeout 21600 --out "$FINAL/dis_veritabani_genis.tsv"
  calistir H2 "the broad set, the second measurement (mfeprimer_layer.py --wide)" \
    python3 "$HERE/mfeprimer_layer.py" --final "$FINAL" \
      --db "$PT/REFERANS_DB" --mfe "$MFE" --wide --cpu "$IS" \
      --timeout 21600 --blast "$FINAL/dis_veritabani_genis.tsv" \
      --out "$FINAL/mfeprimer_genis.tsv"
fi

# --- the summary ------------------------------------------------------
say "----------------------------------------------------------------"
say "SUMMARY"
for f in "$FINAL/hedef_kimlik.tsv" "$FINAL/mfeprimer.tsv" \
         "$FINAL/teslim_denetimi.tsv" \
         "$PT/PrimerJury_Community_Trends.xlsx" \
         "$KONS/B-1_2233851_baskin_konsensus.fasta"; do
  if [ -s "$f" ]; then
    say "  PRESENT  $(basename "$f")  ($(stat -c%s "$f") bytes)"
  else
    say "  MISSING  $(basename "$f")"
  fi
done
if [ ${#BASARISIZ[@]} -gt 0 ]; then
  say "THE STEPS THAT FAILED:"
  for x in "${BASARISIZ[@]}"; do say "   $x"; done
  say "The full logs: $LOGD"
  exit 1
fi
say "every step finished"
say "The logs: $LOGD"
