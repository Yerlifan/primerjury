#!/usr/bin/env bash
# =====================================================================
# run.sh
# It runs every meeting decision end to end and logs each step with a date and
# time. If it is cut short, the same command continues where it stopped, because
# batch_design.py and specificity.py keep checkpoints.
#
#   bash run.sh
#
# The paths can be changed with environment variables:
#   PT=/baska/yol bash run.sh
#
# The steps:
#   0  the dependency check
#   1  the dominant allele consensus (dominant_allele_consensus.py), from the raw
#      reads, with no ambiguity
#   2  the bulk design (batch_design.py), primer candidates per target
#   3  specificity and the raw read verification (specificity.py)
#   4  the external database scan (external_databases.py)
#   5  the reference based design (design_from_reference.py)
#   6  the Excel delivery (export_excel.py)
#   7  the self audit: the regression test (regression_test.py) plus the delivery
#      audit (check_deliverables.py)
# =====================================================================
set -uo pipefail

# The project directory. It used to default to one person's desktop path, so on
# any other machine the defaults below pointed at nothing. It is derived from
# this script's own location now; $PT still overrides it.
_BETIK_DIZIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PT="${PT:-$(cd "$_BETIK_DIZIN/.." && pwd)}"
KONS_IUPAC="${KONS_IUPAC:-$PT/referans_konsensus/self/konsensus}"
BASKIN="${BASKIN:-$PT/referans_konsensus/baskin}"
KONS="${KONS:-$BASKIN/konsensus}"
ADAY="${ADAY:-$PT/primer_adaylari}"
FINAL="${FINAL:-$PT/primer_final}"
HERE="$(cd "$(dirname "$0")" && pwd)"
ANA_LOG="$PT/CALISTIR.log"
MAX_OKUMA="${MAX_OKUMA:-20000}"
TOP="${TOP:-10}"

ts(){ date '+%Y-%m-%d %H:%M:%S'; }
say(){ printf '[%s] %s\n' "$(ts)" "$*" | tee -a "$ANA_LOG"; }

say "================================================================"
say "BASLANGIC"
say "PT          : $PT"
say "IUPAC kons  : $KONS_IUPAC"
say "baskin kons : $KONS"
say "ADAY        : $ADAY"
say "FINAL       : $FINAL"

if [ ! -d "$KONS_IUPAC" ]; then
  say "ERROR: there is no consensus directory: $KONS_IUPAC"
  say "Once 06 ve 07 --mode self calistirilmali."
  exit 1
fi

# --- 0. bagimlilik denetimi -------------------------------------------
say "----------------------------------------------------------------"
say "STEP 0/7  the dependency check"
python3 - <<'PYX' 2>&1 | tee -a "$ANA_LOG"
eksik = []
for m in ("primer3", "Bio", "mappy", "openpyxl"):
    try:
        __import__(m)
        print("   %-10s TAMAM" % m)
    except ImportError:
        eksik.append(m)
        print("   %-10s EKSIK" % m)
if eksik:
    kur = {"Bio": "biopython", "primer3": "primer3-py"}
    paket = " ".join(kur.get(m, m) for m in eksik)
    # On Debian and Ubuntu based distributions the system Python is protected by
    # PEP 668; a plain 'pip install' gives an externally-managed-environment error.
    print("   To install them: pip install --break-system-packages " + paket)
    print("   Alternatif (sistem Python'una dokunmadan):")
    print("     python3 -m venv ~/primer_venv && source ~/primer_venv/bin/activate")
    print("     pip install " + paket + " primer3-py biopython openpyxl mappy")
    if "mappy" in eksik:
        print("   without mappy the indistinguishable bin measurement and the cross")
        print("   contamination measurement CANNOT BE MADE; the run does not stop but")
PYX

# --- 1. baskin alel konsensusu ---------------------------------------
# samtools consensus -A ciktisi degisken pozisyonlari IUPAC koduyla yaziyor.
# Because the binding rule counts a code as matching every base it stands for, if
# the competitor consensus is rich in ambiguity then whichever primer you pick looks
# as though it binds the competitor too. This step produces a second and
# belirsizliksiz bir konsensus uretir.
say "----------------------------------------------------------------"
say "STEP 1/7  the dominant allele consensus (dominant_allele_consensus.py)"
T0=$(date +%s)
MEVCUT=$(ls "$KONS"/*_konsensus.fasta 2>/dev/null | wc -l)
if [ "$MEVCUT" -gt 0 ] && [ "${YENIDEN_KONS:-0}" != "1" ]; then
  say "  it is there already ($MEVCUT files), skipped. To reproduce it: YENIDEN_KONS=1"
else
  python3 "$HERE/dominant_allele_consensus.py" \
      --consensus "$KONS_IUPAC" \
      --fastq "$PT/fastq files" \
      --out "$BASKIN" 2>&1 | tee -a "$ANA_LOG"
fi
say "STEP 1 finished, time=$(( ($(date +%s)-T0)/60 )) minutes"
N=$(ls "$KONS"/*_konsensus.fasta 2>/dev/null | wc -l)
say "baskin alel konsensus dosyasi: $N"
[ "$N" -gt 0 ] || { say "ERROR: no dominant consensus could be produced"; exit 1; }

# --- 2. toplu tasarim -------------------------------------------------
say "----------------------------------------------------------------"
say "STEP 2/7  the bulk design (batch_design.py)"
T0=$(date +%s)
python3 "$HERE/batch_design.py" \
    --consensus "$KONS" \
    --targets "$HERE/hedefler.tsv" \
    --out "$ADAY" 2>&1 | tee -a "$ANA_LOG"
RC=${PIPESTATUS[0]}
say "STEP 2 finished, exit=$RC, time=$(( ($(date +%s)-T0)/60 )) minutes"
[ "$RC" -eq 0 ] || { say "STEP 2 ended with an error, stopped"; exit "$RC"; }

A=$(ls "$ADAY"/*__*.tsv 2>/dev/null | wc -l)
say "uretilen aday dosyasi: $A"

# --- 3. specificity and the raw read verification ---------------------
say "----------------------------------------------------------------"
say "STEP 3/7  specificity and the raw read verification (specificity.py)"
T0=$(date +%s)
python3 "$HERE/specificity.py" \
    --candidates "$ADAY" \
    --pt "$PT" \
    --consensus "$KONS" \
    --targets "$HERE/hedefler.tsv" \
    --out "$FINAL" \
    --top "$TOP" --max-reads "$MAX_OKUMA" 2>&1 | tee -a "$ANA_LOG"
RC=${PIPESTATUS[0]}
say "STEP 3 finished, exit=$RC, time=$(( ($(date +%s)-T0)/60 )) minutes"

# --- 4. the external database scan ------------------------------------
say "----------------------------------------------------------------"
say "STEP 4/7  the external database scan (external_databases.py)"
T0=$(date +%s)
if command -v blastn >/dev/null 2>&1; then
  python3 "$HERE/external_databases.py" \
      --final "$FINAL" --db "$PT/REFERANS_DB" \
      --out "$FINAL/dis_veritabani.tsv" 2>&1 | tee -a "$ANA_LOG"
else
  say "  blastn was not found, the external database scan was skipped."
  say "  To install it: sudo apt-get install -y ncbi-blast+"
fi
say "STEP 4 finished, time=$(( ($(date +%s)-T0)/60 )) minutes"

# --- 5. referans tabanli tasarim (kapsanamayan hedefler) --------------
# For the targets where no pair could be found because the sample's bins cannot be
# separated, a design is made from a reference database. Those pairs ARE NOT
# CONFIRMED against the sample; only whether the sample supports them is measured,
# ayri bir sayfada, ayri etiketle sunulur.
say "----------------------------------------------------------------"
say "STEP 5/7  the reference based design (design_from_reference.py)"
T0=$(date +%s)
REFC="$PT/primer_referans"
if [ -f "$HERE/hedefler_referans.tsv" ]; then
  python3 "$HERE/design_from_reference.py" \
      --db "$PT/REFERANS_DB" --pt "$PT" \
      --reference-targets "$HERE/hedefler_referans.tsv" \
      --out "$REFC" --max-reads "$MAX_OKUMA" 2>&1 | tee -a "$ANA_LOG"
else
  say "  there is no hedefler_referans.tsv, the step was skipped"
fi
say "STEP 5 finished, time=$(( ($(date +%s)-T0)/60 )) minutes"

# --- 6. Excel teslimati -----------------------------------------------
say "----------------------------------------------------------------"
say "STEP 6/7  the Excel delivery (export_excel.py)"
REFARG=""
[ -s "$REFC/primer_referans.tsv" ] && REFARG="--reference $REFC/primer_referans.tsv"
python3 "$HERE/export_excel.py" \
    --candidates "$ADAY" --final "$FINAL" \
    --splits "$ADAY/kume_setleri" \
    --names "$HERE/taxid_adlari.tsv" \
    --targets "$HERE/hedefler.tsv" --consensus "$KONS" $REFARG \
    --out "$PT/PrimerJury_Primer_Tasarimi.xlsx" 2>&1 | tee -a "$ANA_LOG"

# --- 7. the self audit ------------------------------------------------
# Two separate audits run:
#   17  kod kurallarini bagimsiz referans uygulamalarla karsilastirir
#   check_deliverables.py  measures the delivered table again from scratch (it
#       ice aktarmaz), yazili Tm/GC/urun boyu degerlerini de dogrular
# 18 KRITIK bulgu bulursa cikis kodu 1 doner ve burada acikca raporlanir.
say "----------------------------------------------------------------"
say "STEP 7/7  the self audit (regression_test.py plus check_deliverables.py)"
T0=$(date +%s)
python3 "$HERE/regression_test.py" \
    --real-data --candidates "$ADAY" --consensus "$KONS" 2>&1 | tee -a "$ANA_LOG"
RC17=${PIPESTATUS[0]}
if [ -s "$FINAL/primer_final.tsv" ]; then
  python3 "$HERE/check_deliverables.py" \
      --final "$FINAL" --consensus "$KONS" --targets "$HERE/hedefler.tsv" \
      --out "$FINAL/teslim_denetimi.tsv" 2>&1 | tee -a "$ANA_LOG"
  RC18=${PIPESTATUS[0]}
else
  RC18=2
  say "  there is no primer_final.tsv, the delivery audit could not be run"
fi
say "STEP 7 finished, regression=$RC17, delivery=$RC18, time=$(( ($(date +%s)-T0)/60 )) minutes"
[ "$RC17" -eq 0 ] || say "CAREFUL: an item failed in the regression test, look at the log"
[ "$RC18" -eq 0 ] || say "CAREFUL: there is a CRITICAL finding in the delivery audit, look at the log"

# --- the summary ------------------------------------------------------
say "----------------------------------------------------------------"
if [ -s "$FINAL/primer_final.tsv" ]; then
  G=$(awk -F'\t' 'NR>1 && $4=="GECTI"' "$FINAL/primer_final.tsv" | wc -l)
  H=$(awk -F'\t' 'NR>1 && $4=="GECTI" {print $2}' "$FINAL/primer_final.tsv" | sort -u | wc -l)
  say "RESULT: candidates passing every rule = $G, targets covered = $H"
  say "The file: $FINAL/primer_final.tsv"
  say "Excel: $PT/PrimerJury_Primer_Tasarimi.xlsx"
else
  say "UYARI: primer_final.tsv olusmadi"
fi
if [ -s "$ADAY/ayirt_edilemez.tsv" ]; then
  AE=$(( $(wc -l < "$ADAY/ayirt_edilemez.tsv") - 1 ))
  say "INDISTINGUISHABLE TAXON PAIRS: $AE  (the file: $ADAY/ayirt_edilemez.tsv)"
  say "  Bu ciftler dizi duzeyinde ayrilmiyor; rakip listesinden cikarildilar."
fi
say "Loglar:"
say "  $ANA_LOG"
say "  $ADAY/toplu_tasarim.log"
say "  $FINAL/ozgulluk.log"
say "The checkpoints (if it is cut short the same command continues where it stopped):"
say "  $ADAY/checkpoint.json"
say "  $FINAL/checkpoint.json"
say "FINISHED"
