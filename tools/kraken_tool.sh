#!/bin/bash
# ---------------------------------------------------------------------------
# THE KRAKEN2 RERUN TOOL, A SET OF KEYS
#
# WHY THIS TOOL EXISTS
# The request was that "the results should not differ much from Kraken's". We want
# to answer that in Kraken's own language. But rerunning with the same database
# gives the same answer, so two directions were opened:
#   1. A RAISED CONFIDENCE THRESHOLD. Kraken2 assigns a read only if a certain
#      fraction of its k-mers go to the same clade. At threshold 0 a single k-mer is
#      enough. Raising the threshold and looking at what stays standing is Kraken's
#      own measure of how weak those assignments were to begin with.
#   2. A DATABASE OF BROADER COVERAGE (PlusPFP). Our diagnosis was "the problem was
#      coverage". If that is right, the assignments that collapse under a threshold
#      on the narrow database have to stay standing on the broad one. If they do not,
#      the diagnosis is wrong and it is written down as wrong.
#
# WHAT WAS TAKEN FROM THE SOURCE STUDY'S SCRIPTS
#   troubleshooting tools/kraken2_driver.sh  the call form, the --memory-mapping idea
#   troubleshooting tools/bracken_species.sh the kmer_distrib and -l S usage
#   the source study's Kraken/Bracken script  the database integrity check
#                                            (hash/opts/taxo), the memory check, the
#                                            log format
#   tools/rerun_kraken.sh                    finding the database automatically, the
#                                            micromamba environment, merging into one
#                                            run, the losslessness check
#   tools/kraken_summary.py                  parsing the report, summing at species
#                                            level
# This tool does not rewrite them, it calls them and builds on top.
#
# THE KEYS
#   bash kraken_tool.sh memory-config     suggests a .wslconfig (THE FREEZE FIX,
#                                         DO THIS FIRST)
#   bash kraken_tool.sh verify-sample     does the sample represent the full data
#   bash kraken_tool.sh kraken-path       prints kraken2's resolved full path
#                                         (machine readable)
#   bash kraken_tool.sh status            the environment and database check; it
#                                         runs nothing
#   bash kraken_tool.sh find-db           searches the disk for a kraken2 database
#                                         and lists what it finds
#   bash kraken_tool.sh db-identity       which version the database is, a deep
#                                         detection
#   bash kraken_tool.sh original-db       which database the original run used,
#                                         inferred from evidence
#   bash kraken_tool.sh selftest          every self test
#   bash kraken_tool.sh time              a REAL speed measurement from a small
#                                         trial, plus a time estimate
#   bash kraken_tool.sh threshold-a       the confidence threshold scan, on VT_A
#   bash kraken_tool.sh threshold-b       the confidence threshold scan, on VT_B
#   bash kraken_tool.sh threshold         both, then the curves side by side
#   bash kraken_tool.sh table             the four column comparison table
#   bash kraken_tool.sh all               the tests plus the threshold plus the
#                                         table
#   bash kraken_tool.sh custom-db-build   building a custom database (THE LAST
#                                         RESORT)
#   bash kraken_tool.sh custom-db-run     a run with the custom database
#
# THE VARIABLES
#   VT_A=~/k2db       the first database (the default; PlusPFP is expected here)
#   VT_B=/path        the second database (with one, two curves are drawn side by
#                     side)
#   IPLIK=12          the thread count
#   KAP=3000          at most this many reads per taxon (0 = all of them)
#   ESIKLER="0 0.02 0.05 0.1 0.2 0.5"
#   TABLO_ESIK=0.1    the high threshold used in the table
# ---------------------------------------------------------------------------
set -euo pipefail

TUS="${1:-yardim}"

# --- the project directory is measured, not guessed (the same route as
#     rerun_kraken.sh) ---
_BETIK_DIZIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJE="${PROJE:-$(cd "$_BETIK_DIZIN/.." && pwd)}"
# MEASURED: this used to look for a delivery directory that existed in the source
# project. That directory IS NOT in this repository, so whichever key the script was
# called with it stopped at its first check. The marker is now two directories that
# REALLY ARE in the repository.
if [ ! -d "$PROJE/tools" ] || [ ! -d "$PROJE/verification" ]; then
  echo "ERROR: the project root could not be verified ('$PROJE' holds no tools/ and verification/)."
  echo "  Script location: $_BETIK_DIZIN"
  echo "  To give it by hand:  PROJE=/full/path/to/project bash $0 $TUS"
  exit 1
fi

VT_A="${VT_A:-$HOME/k2db}"
VT_B="${VT_B:-}"
# IPLIK: it used to be nproc, that is, every core. Kraken2's PEAK memory use grows
# with the thread count (each thread keeps its own buffer) and under WSL2 that peak
# chokes Windows directly. The default was pulled down to 3; the measurement result
# DOES NOT CHANGE, it only takes a little longer. On a strong machine it can be
# raised with IPLIK=8.
IPLIK="${IPLIK:-3}"
# ORNEK: the TOTAL read count to be used in the threshold scan. 0 = all of them.
# Gerekce asagida esik_tara icinde yazili.
ORNEK="${ORNEK:-100000}"
KAP="${KAP:-3000}"
ESIKLER="${ESIKLER:-0 0.02 0.05 0.1 0.2 0.5}"
TABLO_ESIK="${TABLO_ESIK:-0.1}"
OZELVT="${OZELVT:-$HOME/k2_ozel}"
KAYNAK="$PROJE/RESULTS/fastq files"
IS_A="$PROJE/RESULTS/kraken_esik_A"
IS_B="$PROJE/RESULTS/kraken_esik_B"
OZEL_IS="$PROJE/RESULTS/kraken_ozelvt"

log_ac() {
  local ad="$1"
  mkdir -p "$PROJE/install_logs"
  LOG="$PROJE/install_logs/kraken_${ad}_$(date '+%Y%m%d_%H%M%S').log"
  exec > >(tee -a "$LOG") 2>&1
  echo "=============================================================="
  echo "kraken_tool  key: $ad  start $(date '+%Y-%m-%d %H:%M:%S %Z')"
  echo "  machine: $(hostname)  $(uname -sr)  $IPLIK threads"
  echo "  log   : $LOG"
  echo "=============================================================="
  trap 'echo; echo "finished $(date "+%Y-%m-%d %H:%M:%S %Z")"; echo "log: $LOG"' EXIT
}

# --- the environment. The same as rerun_kraken.sh: kraken2 lives inside the
#     micromamba "mikro" environment --
# ---------------------------------------------------------------------------
# ortam_ac - kraken2'yi BULUNABILIR hale getirir.
#
# NEDEN GENISLETILDI (2026-08-04)
# Bu fonksiyon kaynak calismanin Kraken betiklerinden
# was inherited: kraken2 is not on PATH, it is in micromamba's "mikro"
# environment. But only the "activate the environment" route was being tried. On the
# user's run the micromamba shell hook did not work and, because kraken2 could not be
# found, every Kraken step was skipped. The binary was sitting on the disk all
# along.
#
# Four routes are now tried IN ORDER and it stops at the first that holds:
#   0) if KRAKEN2_BIN was given from outside it is used directly (the highest
#      authority).
#   1) Zaten PATH'te mi.
#   2) THE SOURCE STUDY'S ROUTE, unchanged and kept exactly: $HOME/bin is added to
#      PATH, MAMBA_ROOT_PREFIX is set, and "$ORTAM" is opened with micromamba or
#      conda.
#   3) The environment DIRECTORIES are looked at directly. Even when the shell hook
#      does not work the binary is in place; there is no need to activate anything,
#      adding the bin directory to PATH is enough.
#   4) Son care: $HOME altinda sinirli derinlikte kraken2 adli calistirilabilir
#      file is searched for.
#
# The path found is written into the KRAKEN2_BIN variable and exported, so the
# later steps do not depend on PATH.
# ---------------------------------------------------------------------------
ortam_ac() {
  ORTAM="${ORTAM:-mikro}"
  export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
  KRAKEN_YONTEM=""
  BAKILAN_YERLER=""

  _kayit() { BAKILAN_YERLER="${BAKILAN_YERLER}
    $1"; }

  # Is kraken2 really RUNNING. The file existing is not enough: conda binaries can
  # depend on the environment's lib directory and, called from outside the
  # environment, they can fail to find a library and blow up. So --version is
  # actually tried.
  _calisiyor_mu() { "$1" --version >/dev/null 2>&1; }

  _benimse() {
    KRAKEN2_BIN="$1"; export KRAKEN2_BIN
    export PATH="$(dirname "$1"):$PATH"
    KRAKEN_YONTEM="$2"
    return 0
  }

  # --- 0) disaridan verilen tam yol her seyin ustundedir -------------------
  if [ -n "${KRAKEN2_BIN:-}" ] && [ -x "${KRAKEN2_BIN}" ]; then
    _benimse "$KRAKEN2_BIN" "disaridan verilen tam yol (KRAKEN2_BIN)"; return 0
  fi

  # --- 1) PATH ------------------------------------------------------------
  _kayit "PATH"
  if command -v kraken2 >/dev/null 2>&1; then
    _benimse "$(command -v kraken2)" "PATH uzerinde bulundu"; return 0
  fi

  # --- 2) THE ENVIRONMENT DIRECTORIES, directly ---------------------------
  # This is looked at before running the shell hook. The reason: in this
  # installation there is a single environment (mikro) and the binary is there;
  # running the hook is both slow and, in some shells, silently unsuccessful.
  # Looking at the file directly is faster and more certain.
  local k
  for k in \
      "$HOME/micromamba/envs/$ORTAM/bin/kraken2" \
      "$MAMBA_ROOT_PREFIX/envs/$ORTAM/bin/kraken2" \
      "$HOME/micromamba/envs"/*/bin/kraken2 \
      "$HOME/miniconda3/envs/$ORTAM/bin/kraken2" \
      "$HOME/miniconda3/envs"/*/bin/kraken2 \
      "$HOME/anaconda3/envs"/*/bin/kraken2 \
      "$HOME/mambaforge/envs"/*/bin/kraken2 \
      "$HOME/miniforge3/envs"/*/bin/kraken2 \
      "$HOME/conda/envs"/*/bin/kraken2 \
      "$HOME/bin/kraken2" \
      /opt/*/envs/*/bin/kraken2 \
      /usr/local/bin/kraken2 ; do
    case "$k" in *'*'*) continue ;; esac
    _kayit "$k"
    [ -x "$k" ] || continue
    if _calisiyor_mu "$k"; then
      _benimse "$k" "found in the environment directory, it runs directly"
      return 0
    fi
    # The binary is there but does not run directly: it needs the environment's
    # libraries. Calling it through micromamba run solves that WITHOUT DISTURBING
    # THE ENVIRONMENT
    # (kaynak calismanin kraken2_driver.sh betigi kraken2'yi cipl ak cagiriyordu; biz
    # ortami etkinlestirmek yerine tek komutluk kabuk aciyoruz).
    local envad; envad="$(basename "$(dirname "$(dirname "$k")")")"
    if command -v micromamba >/dev/null 2>&1 && \
       micromamba run -n "$envad" kraken2 --version >/dev/null 2>&1; then
      local shim; shim="$(mktemp -d)"
      local a
      for a in kraken2 kraken2-build kraken2-inspect bracken; do
        printf '#!/bin/bash\nexec micromamba run -n %s %s "$@"\n' "$envad" "$a" > "$shim/$a"
        chmod +x "$shim/$a"
      done
      export PATH="$shim:$PATH"
      KRAKEN2_BIN="$shim/kraken2"; export KRAKEN2_BIN
      KRAKEN_YONTEM="micromamba run -n $envad (the binary did not run directly)"
      return 0
    fi
  done

  # --- 3) THE SOURCE STUDY'S ROUTE: activate the environment ---------------
  # Kaynak calismanin Kraken betiklerinden aynen devralindi.
  _kayit "the micromamba or conda shell hook plus activating the '$ORTAM' environment"
  export PATH="$HOME/bin:$PATH"
  if command -v micromamba >/dev/null 2>&1; then
    eval "$(micromamba shell hook -s bash)" || true
    micromamba activate "$ORTAM" || true
  elif command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)" || true
    conda activate "$ORTAM" || true
  fi
  if command -v kraken2 >/dev/null 2>&1; then
    _benimse "$(command -v kraken2)" "found by activating the environment ($ORTAM)"
    return 0
  fi

  # --- 4) micromamba run, without activating the environment --------------
  _kayit "micromamba run -n $ORTAM kraken2"
  if command -v micromamba >/dev/null 2>&1 && \
     micromamba run -n "$ORTAM" kraken2 --version >/dev/null 2>&1; then
    local shim; shim="$(mktemp -d)"
    local a
    for a in kraken2 kraken2-build kraken2-inspect bracken; do
      printf '#!/bin/bash\nexec micromamba run -n %s %s "$@"\n' "$ORTAM" "$a" > "$shim/$a"
      chmod +x "$shim/$a"
    done
    export PATH="$shim:$PATH"
    KRAKEN2_BIN="$shim/kraken2"; export KRAKEN2_BIN
    KRAKEN_YONTEM="micromamba run -n $ORTAM"
    return 0
  fi

  # --- 5) son care: ev dizininde arama -----------------------------------
  _kayit "an executable named 'kraken2' under \$HOME, down to a depth of 6"
  k="$(find "$HOME" -maxdepth 6 -type f -name kraken2 -perm -u+x 2>/dev/null | head -1)"
  if [ -n "$k" ] && _calisiyor_mu "$k"; then
    _benimse "$k" "found by searching the home directory"; return 0
  fi
  export BAKILAN_YERLER
  return 1
}

# ---------------------------------------------------------------------------
# tus_kraken_yol prints kraken2's RESOLVED full path, machine readable.
#
# Why it is a key of its own: full_chain.py (key A) used to look for kraken2 on
# PATH by itself and, failing to find it, skipped every Kraken step. Now
# aramayi yapmiyor, bu tusu cagirip cevabi buradan aliyor. Arama mantigi TEK
# ONE PLACE; there are no longer two different searches in two places.
#
# The output: on success one line, KRAKEN2_BIN=<full path>, and exit 0.
#             on failure exit 1 and the install instructions on stderr.
# ---------------------------------------------------------------------------
tus_kraken_yol() {
  if ortam_ac; then
    echo "KRAKEN2_BIN=${KRAKEN2_BIN}"
    echo "KRAKEN_YONTEM=${KRAKEN_YONTEM}"
    echo "KRAKEN2_SURUM=$(kraken2 --version 2>&1 | head -1)"
    return 0
  fi
  {
    echo "kraken2 was not found. These are EXACTLY the places that were checked:"
    printf '%s\n' "$BAKILAN_YERLER" | sed '/^$/d'
    echo
    echo "The place this project expects:  \$HOME/micromamba/envs/${ORTAM:-mikro}/bin/kraken2"
    echo "  (with \$HOME expanded that is: $HOME/micromamba/envs/${ORTAM:-mikro}/bin/kraken2)"
    echo
    echo "INSTALL - the route this project uses:"
    echo "    bash install.sh tools"
    echo "To install it by hand:"
    echo "    micromamba create -n mikro -c bioconda -c conda-forge kraken2 bracken"
    echo
    echo "IF IT IS ALREADY INSTALLED the environment may have a different name:"
    echo "    micromamba env list        (or: conda env list)"
    echo "    ORTAM=<environment_name> bash $0 kraken-path"
    echo "If you know the full path of the binary:"
    echo "    KRAKEN2_BIN=/full/path/kraken2 bash $0 kraken-path"
  } >&2
  return 1
}

kraken_sart() {
  ortam_ac
  if ! command -v kraken2 >/dev/null 2>&1; then
    echo
    echo "ERROR: kraken2 was not found. This is not skipped silently, the work stops here."
    echo
    echo "  Needed: kraken2 version 2.1 or above."
    echo "  The installation route this project uses:"
    echo "      bash install.sh tools"
    echo "  By hand:"
    echo "      micromamba create -n mikro -c bioconda -c conda-forge kraken2 bracken"
    echo "  If it is already installed, activate the environment:"
    echo "      export PATH=\"\$HOME/bin:\$PATH\""
    echo "      export MAMBA_ROOT_PREFIX=\"\$HOME/micromamba\""
    echo "      eval \"\$(micromamba shell hook -s bash)\""
    echo "      micromamba activate ${ORTAM:-mikro}"
    echo "  If you do not know the environment name:  micromamba env list"
    exit 1
  fi
  echo "kraken2: ${KRAKEN2_BIN:-$(command -v kraken2)}   $(kraken2 --version 2>&1 | head -1)"
  [ -n "${KRAKEN_YONTEM:-}" ] && echo "  how it was found: $KRAKEN_YONTEM"
  for _a in kraken2-inspect bracken; do
    if command -v "$_a" >/dev/null 2>&1; then echo "  $_a: $(command -v $_a)"
    else echo "  $_a: MISSING"; fi
  done
}

# =========================================================================
# THE DATABASE READINESS CHECK
# It assumes nothing. All three files are looked for and their size and date are
# written down. If only a .tar.gz is there, the unpack command is printed but
# NOTHING IS UNPACKED BY ITSELF; unpacking tens of gigabytes without the user
# knowing is not acceptable.
# Returns: 0 ready, 1 only the archive is there, 2 nothing at all
# =========================================================================
vt_hazir_mi() {
  local d="$1"
  if [ ! -d "$d" ]; then
    echo "  no such directory: $d"
    return 2
  fi
  local eksik=0
  for g in hash.k2d opts.k2d taxo.k2d; do
    if [ -s "$d/$g" ]; then
      printf "    %-10s present   %8s   %s\n" "$g" \
        "$(du -h "$d/$g" | cut -f1)" \
        "$(date -r "$d/$g" '+%Y-%m-%d %H:%M' 2>/dev/null || echo '?')"
    else
      printf "    %-10s ABSENT\n" "$g"
      eksik=1
    fi
  done
  if [ "$eksik" -eq 0 ]; then
    echo "    -> the index is ready"
    return 0
  fi
  local ars
  ars=$(ls "$d"/*.tar.gz 2>/dev/null | head -3 || true)
  if [ -n "$ars" ]; then
    echo
    echo "    THE INDEX IS NOT READY. There is an archive in the directory but it is not unpacked:"
    while IFS= read -r t; do
      [ -n "$t" ] && printf "      %-52s %8s  %s\n" "$(basename "$t")" \
        "$(du -h "$t" | cut -f1)" "$(date -r "$t" '+%Y-%m-%d' 2>/dev/null || echo '?')"
    done <<< "$ars"
    echo
    echo "    To unpack it (this does NOT happen on its own, it takes tens of GB):"
    echo "        cd $d"
    while IFS= read -r t; do
      [ -n "$t" ] && echo "        tar -xzvf $(basename "$t")"
    done <<< "$ars"
    echo "    First measure the free space:  df -h $d"
    echo "    After unpacking:     bash $0 db-identity"
    return 1
  fi
  echo "    NO INDEX and no archive either: $d"
  return 2
}

# =========================================================================
# DATABASE VERSION DETECTION
# This detection matters: our whole argument rests on the database's COVERAGE. Run
# with the wrong version, we would read the result wrongly.
#
# TWO INDEPENDENT MEASUREMENTS (project rule 1)
#   measurement 1: the taxa INSIDE the database, through kraken2-inspect. Definite.
#   measurement 2: the size of hash.k2d. Rough, but independent.
# If the two disagree, a DISAGREEMENT is written and it does not fall back on one.
#
# What tells them apart: plants (Viridiplantae, taxid 33090).
#   Standard = archaea, bacteria, viruses, plasmids, human, UniVec
#   PlusPF   = Standard + protozoa + fungi
#   PlusPFP  = PlusPF + PLANTS
# =========================================================================
vt_surum() {
  local d="$1"
  echo "  version detection: $d"
  local hb hgb
  hb=$(stat -c%s "$d/hash.k2d" 2>/dev/null || echo 0)
  hgb=$((hb / 1073741824))
  echo "    the size of hash.k2d: ${hgb} GB"

  # --- measurement 2: the size ---
  local o2
  if   [ "$hgb" -ge 125 ]; then o2="PlusPFP (full)"
  elif [ "$hgb" -ge 95  ]; then o2="PlusPF (full)"
  elif [ "$hgb" -ge 60  ]; then o2="Standard (full)"
  elif [ "$hgb" -ge 12  ]; then o2="capped version (16 GB class)"
  elif [ "$hgb" -ge 5   ]; then o2="capped version (8 GB class)"
  else                          o2="a small or custom database"
  fi
  echo "    measurement 2 (size)    : $o2"

  # --- measurement 1: inspect ---
  local o1="not measured"
  local ins="$PROJE/RESULTS/vt_inspect_$(basename "$d").txt"
  if command -v kraken2-inspect >/dev/null 2>&1; then
    # ONBELLEK, AMA KAYNAK DAMGASIYLA (2026-08-04 duzeltmesi)
    # It used to look only at "is the file there" and, if it was, it did not run
    # again. That invited this project's classic kind of fault:
    # baska bir veritabanina ait ya da bayat bir inspect ciktisi sessizce
    # would be reused and the tool would report the WRONG version without any
    # error. The source's size and date stamp is now kept beside the output; if the
    # stamp does not hold, the cache counts as invalid and inspect RUNS AGAIN.
    local damga_dosya="$ins.kaynak"
    local damga_simdi
    damga_simdi="$(stat -c '%s %Y' "$d/hash.k2d" 2>/dev/null || echo 'none')"
    local damga_eski=""
    [ -f "$damga_dosya" ] && damga_eski="$(cat "$damga_dosya" 2>/dev/null)"
    if [ ! -s "$ins" ] || [ "$damga_simdi" != "$damga_eski" ]; then
      if [ -s "$ins" ]; then
        echo "    the cached inspect output belongs to a DIFFERENT database"
        echo "      (the stamp changed: '$damga_eski' -> '$damga_simdi'), re-running"
      fi
      echo "    kraken2-inspect is running (this can take 1 to 5 minutes, once)"
      mkdir -p "$(dirname "$ins")"
      if kraken2-inspect --db "$d" --threads "$IPLIK" > "$ins" 2>/dev/null; then
        printf '%s' "$damga_simdi" > "$damga_dosya"
      else
        echo "    kraken2-inspect DID NOT RUN (its exit code was not zero)."
        : > "$ins"
        rm -f "$damga_dosya"
      fi
    else
      echo "    the kraken2-inspect output is cached and the stamp matches: $ins"
    fi
    if [ -s "$ins" ]; then
      local bitki mantar protozoa insan virus arke
      bitki=$(awk -F'\t' '$5=="33090"{print $2+0}' "$ins" | head -1); bitki=${bitki:-0}
      mantar=$(awk -F'\t' '$5=="4751"{print $2+0}' "$ins" | head -1); mantar=${mantar:-0}
      protozoa=$(awk -F'\t' '$5=="5794"||$5=="33682"{s+=$2} END{print s+0}' "$ins")
      insan=$(awk -F'\t' '$5=="9606"{print $2+0}' "$ins" | head -1); insan=${insan:-0}
      virus=$(awk -F'\t' '$5=="10239"{print $2+0}' "$ins" | head -1); virus=${virus:-0}
      arke=$(awk -F'\t' '$5=="2157"{print $2+0}' "$ins" | head -1); arke=${arke:-0}
      printf "    content: archaea %s, viruses %s, human %s, fungi %s, protozoa %s, PLANTS %s\n" \
             "$arke" "$virus" "$insan" "$mantar" "$protozoa" "$bitki"
      if   [ "$bitki" -gt 0 ] && [ "$mantar" -gt 0 ]; then o1="PlusPFP"
      elif [ "$mantar" -gt 0 ]; then o1="PlusPF (NO plants)"
      elif [ "$arke" -gt 0 ]; then o1="Standard (NO fungi and NO plants)"
      else o1="a custom or unidentifiable database"
      fi
    fi
  else
    echo "    kraken2-inspect is missing, so measurement 1 could not be made."
      echo "      To install it: micromamba install -n ${ORTAM:-mikro} -c bioconda kraken2"
  fi
  echo "    measurement 1 (content) : $o1"

  # --- iki olcumun uzlasisi ---
  echo
  case "$o1:$o2" in
    "not measured:"*)
      echo "    RESULT: UNCERTAIN. Only the size was looked at; the content could not be measured."
      echo "    Size ON ITS OWN is not enough to tell PlusPF from PlusPFP."
      echo "    this question does not close until kraken2-inspect is installed." ;;
    "PlusPFP:PlusPFP (full)")
      echo "    RESULT: PlusPFP. Both measurements say the same thing." ;;
    "PlusPF (NO plants):PlusPF (full)")
      echo "    RESULT: PlusPF. Both measurements say the same thing."
      echo "    CAUTION: this is NOT PlusPFP. It holds no plants. If you expected PlusPFP,"
      echo "    the downloaded archive may be the wrong one. The PlusPFP address:"
      echo "      https://genome-idx.s3.amazonaws.com/kraken/  (k2_pluspfp_*.tar.gz)" ;;
    *)
      echo "    DISAGREEMENT: the content measurement says '$o1', the size measurement says '$o2'."
      echo "    The two measurements disagree, and neither is silently preferred. The content"
      echo "    measurement is the more reliable one, but first confirm the archive unpacked fully." ;;
  esac
  VT_SURUM_SONUC="$o1 / $o2"

  # --- the fingerprint. Which database an old run used is traced from here. ---
  local pi="$PROJE/RESULTS/vt_parmak_izi.tsv"
  mkdir -p "$(dirname "$pi")"
  [ -f "$pi" ] || printf 'tarih\tyol\thash_bayt\thash_tarih\ticerik\tboyut\n' > "$pi"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$(date '+%Y-%m-%d %H:%M')" "$d" "$hb" \
    "$(date -r "$d/hash.k2d" '+%Y-%m-%d' 2>/dev/null || echo '?')" "$o1" "$o2" >> "$pi"
}

# =========================================================================
# WHICH DATABASE THE OLD RUN USED
# The question is critical: if the old run was made with ~/k2db too, then "the old
# database" and "PlusPFP" may be the same thing and the comparison means something
# else.
# =========================================================================
eski_kosu_tespit() {
  echo "WHICH DATABASE THE OLD RUNS USED"
  echo
  echo "  a) A direct record from the logs:"
  local bulundu=0
  for l in "$PROJE"/install_logs/*.log; do
    [ -f "$l" ] || continue
    local yol boy
    # This grep reads OUR OWN logs. The screen text was translated into English,
    # but the old logs are still in Turkish; the pattern takes both words so that
    # either matches.
    yol=$(grep -m1 -oE "(veritabani|database) *: */[^ ]*" "$l" 2>/dev/null | sed 's/.*: *//' || true)
    boy=$(grep -m1 -oE "hash\.k2d *(var|present) *[0-9.]+[KMGT]" "$l" 2>/dev/null | grep -oE "[0-9.]+[KMGT]$" || true)
    if [ -n "$yol" ]; then
      printf "     %-46s  %-28s  hash %s\n" "$(basename "$l")" "$yol" "${boy:-?}"
      bulundu=1
    fi
  done
  [ "$bulundu" -eq 1 ] || echo "     no log record was found"
  echo
  echo "  b) Indirect evidence from the outputs (which domains appear in the reports):"
  for r in "$PROJE/RESULTS/kraken_yeniden/tum.report" \
           "$PROJE/RESULTS/kraken results"/*/*_kraken2.report; do
    [ -f "$r" ] || continue
    local m p b
    m=$(awk -F'\t' '$5=="4751"{print "mantar VAR"}' "$r" | head -1)
    p=$(awk -F'\t' '$5=="5794"||$5=="33682"{print "protozoa VAR"}' "$r" | head -1)
    b=$(awk -F'\t' '$5=="33090"{print "BITKI VAR"}' "$r" | head -1)
    printf "     %-52s %s %s %s\n" "$(basename "$(dirname "$r")")/$(basename "$r")" \
           "${m:-mantar yok}" "${p:-protozoa yok}" "${b:-bitki yok}"
    break
  done
  local ilk="$PROJE/RESULTS/kraken_yeniden/tum.report"
  if [ -f "$ilk" ]; then
    local m p b
    m=$(awk -F'\t' '$5=="4751"{print "VAR"}' "$ilk" | head -1)
    p=$(awk -F'\t' '$5=="5794"{print "VAR"}' "$ilk" | head -1)
    b=$(awk -F'\t' '$5=="33090"{print "VAR"}' "$ilk" | head -1)
    echo "     kraken_yeniden/tum.report: fungi ${m:-none}, protozoa ${p:-none}, plants ${b:-none}"
  fi
  echo
  echo "  THE LIMIT OF THIS EVIDENCE, stated openly:"
  echo "  Fungi and protozoa APPEARING in the report PROVES that the database"
  echo "  contained them. So the old runs used at least PlusPF, not Standard."
  echo "  But plants NOT APPEARING proves nothing: if no plant read was assigned in"
  echo "  a digester sample, no line is written even when the database holds plants."
  echo "  ABSENCE IS NOT EVIDENCE. PlusPF and PlusPFP are separated only by"
  echo "  kraken2-inspect:  bash $0 db-identity"
}

# =========================================================================
# THE KEY: vt-ara   SEARCHING FOR A DATABASE
# If the path is wrong it does not give up, it searches. The whole disk is not
# scanned; particular root directories and a limited depth are used, and the time
# taken is measured and written down. If it finds more than one it DOES NOT CHOOSE
# ON ITS OWN, it lets the user choose.
# =========================================================================
vt_ara() {
  local derin="${DERINLIK:-5}"
  local kokler=("$HOME" /opt /srv /data /media /run/media /mnt)
  echo "DATABASE SEARCH"
  echo "  aranan   : hash.k2d"
  echo "  kokler   : ${kokler[*]}"
  echo "  depth    : $derin (the WHOLE disk is NOT scanned, so that it finishes in reasonable time)"
  echo
  local t0 t1
  t0=$(date +%s)
  local bulunanlar=()
  for k in "${kokler[@]}"; do
    [ -d "$k" ] || continue
    while IFS= read -r h; do
      [ -n "$h" ] && bulunanlar+=("$(dirname "$h")")
    done < <(timeout "${ARAMA_SN:-120}" find "$k" -maxdepth "$derin" -type f \
              -name hash.k2d 2>/dev/null || true)
  done
  t1=$(date +%s)
  # under /mnt/c the Windows side can be very slow, so it is searched separately
  # and less deeply.
  if [ -d /mnt/c/Users ] && [ "${WINDOWS_ARA:-1}" = "1" ]; then
    while IFS= read -r h; do
      [ -n "$h" ] && bulunanlar+=("$(dirname "$h")")
    done < <(timeout "${ARAMA_SN:-120}" find /mnt/c/Users -maxdepth 4 -type f \
              -name hash.k2d 2>/dev/null || true)
    t1=$(date +%s)
  fi
  # yinelenenleri at
  local benzersiz=()
  while IFS= read -r y; do [ -n "$y" ] && benzersiz+=("$y"); done \
    < <(printf '%s\n' "${bulunanlar[@]+"${bulunanlar[@]}"}" | sort -u)
  echo "  search time: $((t1 - t0)) seconds"
  echo
  if [ "${#benzersiz[@]}" -eq 0 ]; then
    echo "  NO KRAKEN2 DATABASE WAS FOUND AT ALL."
    echo
    echo "  There may be only a .tar.gz, an archive that was never unpacked:"
    for k in "${kokler[@]}"; do
      [ -d "$k" ] || continue
      timeout 60 find "$k" -maxdepth "$derin" -name "k2_*.tar.gz" 2>/dev/null | head -5 \
        | while IFS= read -r t; do printf "    %s  %s\n" "$(du -h "$t" | cut -f1)" "$t"; done
    done
    echo
    echo "  If there is an archive, unpack it (it does NOT unpack itself):"
    echo "      cd <directory> && tar -xzvf k2_*.tar.gz"
    echo "  If there is no archive either, download one:"
    echo "      https://benlangmead.github.io/aws-indexes/k2"
    echo "      The full PlusPFP is about 150 GB, and the downloaded archive about 90 GB."
    echo "  If it is somewhere else, give it directly:  VT_A=/full/path bash $0 $TUS"
    return 1
  fi
  echo "  ${#benzersiz[@]} databases found:"
  echo
  local i=0
  for y in "${benzersiz[@]}"; do
    i=$((i+1))
    echo "  [$i] $y"
    printf "      hash.k2d %8s   %s\n" "$(du -h "$y/hash.k2d" 2>/dev/null | cut -f1)" \
           "$(date -r "$y/hash.k2d" '+%Y-%m-%d %H:%M' 2>/dev/null || echo '?')"
    local eksik=""
    for g in opts.k2d taxo.k2d; do [ -s "$y/$g" ] || eksik="$eksik $g"; done
    if [ -n "$eksik" ]; then
      echo "      MISSING FILES:$eksik  (the index cannot be used)"
    else
      echo "      total $(du -sh "$y" 2>/dev/null | cut -f1)"
      vt_surum "$y" | sed 's/^/      /'
    fi
    echo
  done
  if [ "${#benzersiz[@]}" -eq 1 ]; then
    echo "  One database was found. To use it:"
    echo "      VT_A=${benzersiz[0]} bash $0 threshold"
  else
    echo "  THERE IS MORE THAN ONE DATABASE. The script does NOT CHOOSE which one to run,"
    echo "  because running the wrong version means reading the result wrongly."
    echo "  Look at the version lines above, choose, and give the path:"
    for y in "${benzersiz[@]}"; do echo "      VT_A=$y bash $0 threshold"; done
    echo "  If you want to compare two databases:"
    echo "      VT_A=<wide> VT_B=<narrow> bash $0 threshold"
  fi
  return 0
}

# =========================================================================
# THE KEY: ozgun-vt   WHICH DATABASE THE ORIGINAL RUN USED
# The user does not know either and there is nobody to ask. It is inferred from
# evidence. If it cannot be inferred, "belirlenemedi" is written; it is not
# invented.
# =========================================================================
tus_ali_vt() {
  # This key only READS and writes a report, it changes nothing. A file it looks
  # for being absent is normal; missing evidence must not stop the work, it has to
  # be written down as "bulunamadi". That is why errexit is turned off here.
  set +e
  echo "======================================================================"
  echo "WHICH DATABASE THE ORIGINAL RUN USED"
  echo "======================================================================"
  echo
  echo "EVIDENCE 1: the --db argument in the source study's scripts"
  local k1=""
  for f in "$PROJE/troubleshooting tools"/*.sh; do
    [ -f "$f" ] || continue
    grep -HnoE "(KRAKEN[A-Z0-9_]*(DB|PATH|FOLDER)|--db)[= ]*\"?[^\" ]+" "$f" 2>/dev/null \
      | sed 's|^.*/||; s/^/    /' || true
  done
  k1=$(grep -rhoE "/[^\" ]*(pluspf|pluspfp|standard|k2[a-z0-9_]*)[^\" ]*" \
        "$PROJE/troubleshooting tools"/*.sh 2>/dev/null | sort -u | head -3)
  if [ -n "$k1" ]; then
    echo
    echo "    -> the path(s) the scripts point at:"
    printf '       %s\n' $k1
    if printf '%s' "$k1" | grep -qi "pluspf16\|pluspf_16\|_16gb"; then
      echo "       The directory name says PlusPF-16: the version of PlusPF CAPPED to 16 GB."
      echo "       A capped version holds THE SAME taxa as the full one, but the minimizer"
      echo "       count is thinned. The same read can land on a higher node, or even on"
      echo "       a different clade. So it is not the coverage that is lower, it is the RESOLUTION."
    fi
  else
    echo "    no database path was found in the scripts"
  fi

  echo
  echo "EVIDENCE 2: the run logs"
  local k2=0
  for l in "$PROJE"/install_logs/*.log; do
    [ -f "$l" ] || continue
    local yol boy
    # This grep reads OUR OWN logs. The screen text was translated into English,
    # but the old logs are still in Turkish; the pattern takes both words so that
    # either matches.
    yol=$(grep -m1 -oE "(veritabani|database) *: */[^ ]*" "$l" 2>/dev/null | sed 's/.*: *//' || true)
    boy=$(grep -m1 -oE "hash\.k2d *(var|present) *[0-9.]+[KMGT]" "$l" 2>/dev/null | grep -oE "[0-9.]+[KMGT]$" || true)
    [ -n "$yol" ] && { printf "    %-42s %-26s hash %s\n" "$(basename "$l")" "$yol" "${boy:-?}"; k2=1; }
  done
  [ "$k2" -eq 1 ] || echo "    there is no log from the original run (these logs are from our runs)"

  echo
  echo "EVIDENCE 3: report content (which domains WERE in the database)"
  local r="$PROJE/RESULTS/kraken results/A1/edited_barcode01_kraken2.report"
  if [ -f "$r" ]; then
    local m p b v
    m=$(awk -F'\t' '$5=="4751"{print "VAR"}' "$r" | head -1)
    p=$(awk -F'\t' '$5=="5794"{print "VAR"}' "$r" | head -1)
    b=$(awk -F'\t' '$5=="33090"{print "VAR"}' "$r" | head -1)
    v=$(awk -F'\t' '$5=="10239"{print "VAR"}' "$r" | head -1)
    echo "    fungi (4751)       : ${m:-no line}"
    echo "    Apicomplexa (5794) : ${p:-no line}"
    echo "    plants (33090)     : ${b:-no line}"
    echo "    viruses (10239)    : ${v:-no line}"
    echo
    echo "    If fungi and protozoa APPEAR, the database DID CONTAIN them."
    echo "    So the original run's database covers at least PlusPF, not Standard."
    echo "    Plants not appearing proves nothing: if no plant read was assigned in a"
    echo "    digester sample, no line is written even when the database holds plants."
    echo "    ABSENCE IS NOT EVIDENCE."
  else
    echo "    the source study's report was not found: $r"
  fi

  echo
  echo "EVIDENCE 4: did our re-run reproduce the source study's result"
  local oz="$PROJE/RESULTS/kraken_yeniden/kraken_ozet.csv"
  if [ -f "$oz" ]; then
    local top uy ay
    top=$(awk -F',' 'NR>1{n++} END{print n+0}' "$oz")
    uy=$(awk -F',' 'NR>1 && $9=="uyusuyor"{n++} END{print n+0}' "$oz")
    ay=$((top - uy))
    echo "    $uy of $top bins gave the same answer, $ay gave a DIFFERENT one."
    echo
    echo "    THIS PRODUCES AN INFERENCE. Kraken2 is deterministic: the same read, the"
    echo "    same database and the same threshold always give the same assignment. Since"
    echo "    the bins were split by the source study's own output, every bin would repeat"
    echo "    its own label if the database were THE SAME."
    if [ "$ay" -gt 0 ]; then
      echo "    $ay bins did not repeat, so the database behind our re-run is"
      echo "    NOT THE SAME as the database behind the original run."
      echo
      echo "    The limits of that inference, stated honestly:"
      echo "      it ran with extract_kraken_reads.py --include-children, so a bin also"
      echo "      holds reads belonging to the taxon's CHILD nodes. Species level"
      echo "      aggregation closes most of that but not all of it. Against that,"
      echo "      splits that jump genus and family (a Bacteroides bin coming out as"
      echo "      Candidatus Azobacteroides, for instance) cannot be explained by child nodes."
    else
      echo "    Every bin repeated. That points to the two runs using the same or a very"
      echo "    similar database."
    fi
  else
    echo "    there is no kraken_ozet.csv, so this evidence could not be measured"
  fi

  echo
  echo "======================================================================"
  echo "CONCLUSION"
  echo "======================================================================"
  if printf '%s' "$k1" | grep -qi "pluspf16\|pluspf_16\|_16gb"; then
    echo "  the original run's database: PlusPF-16 (the 16 GB capped version of PlusPF),"
    echo "  according to the path written in its script. This is NOT CERTAIN; the path in"
    echo "  a script and what actually ran can differ. But it is the most direct evidence"
    echo "  we have and it does not contradict the report content (fungi and protozoa present)."
  else
    echo "  NOT DETERMINED. The evidence we have says only this much:"
    echo "  the database contained fungi and protozoa, so it covered at least PlusPF."
    echo "  Whether it was capped or full, PlusPFP or not, could not be inferred."
  fi
  echo
  echo "  THIS UNCERTAINTY AFFECTS THE READING, and must be stated openly:"
  echo "  If the old run already used a wide database, the diagnosis that the problem"
  echo "  was coverage weakens, and the nanopore error rate becomes the likelier cause."
  echo "  The way to tell the two explanations apart is still a PlusPFP run:"
  echo "    if it is a coverage problem -> assignments strengthen under PlusPFP and hold at the threshold"
  echo "    if it is a read error       -> assignments stay weak under PlusPFP too and collapse at the threshold"
  echo "  The reading of the table will carry that uncertainty, not hide it."
  set -e
}

# --- memory. The decision from the source Kraken and Bracken script and from
#     rerun_kraken.sh --------
bellek_bayragi() {
  local d="$1"
  local hb rb
  hb=$(stat -c%s "$d/hash.k2d" 2>/dev/null || echo 0)
  rb=$(( $(awk '/MemAvailable/{print $2}' /proc/meminfo 2>/dev/null || echo 0) * 1024 ))
  {
    echo "  hash.k2d       : $(numfmt --to=iec "$hb" 2>/dev/null || echo "$hb")"
    echo "  available      : $(numfmt --to=iec "$rb" 2>/dev/null || echo "$rb")"
  } >&2
  if [ "${ZORLA_MMAP:-1}" = "1" ] || [ "$hb" -gt $((rb - 2147483648)) ]; then
    {
      echo "  --memory-mapping WILL BE USED. The database is not loaded into RAM, it is read"
      echo "  from disk. The RAM requirement drops and the run SLOWS DOWN. The"
      echo "  classification result is THE SAME; the only thing that changes is the time."
    } >&2
    echo "--memory-mapping"
  else
    echo "  RAM is sufficient, the database will be loaded into memory (fast)." >&2
    echo ""
  fi
}

# --- okumalari tek dosyada birlestir (rerun_kraken.sh'ten aynen) --------
birlestir() {
  local hedef="$1"
  local eski="$PROJE/RESULTS/kraken_yeniden/tum.fastq"
  if [ -s "$hedef" ]; then
    echo "the merged file already exists: $(awk 'END{print int(NR/4)}' "$hedef") reads"
    return
  fi
  mkdir -p "$(dirname "$hedef")"
  if [ -s "$eski" ]; then
    echo "the tum.fastq produced by rerun_kraken.sh was found, so THE SAME read set is used."
    echo "  (That is the condition that makes thresholds and databases comparable.)"
    cp "$eski" "$hedef"
    cp "$PROJE/RESULTS/kraken_yeniden/kaynak_sayim.tsv" "$(dirname "$hedef")/" 2>/dev/null || true
    echo "  $(awk 'END{print int(NR/4)}' "$hedef") reads"
    return
  fi
  [ -d "$KAYNAK" ] || { echo "ERROR: no fastq directory: $KAYNAK"; exit 1; }
  local dosyalar sayi
  dosyalar=$(ls "$KAYNAK"/*/*reads_*.fastq 2>/dev/null || true)
  sayi=$(printf '%s\n' "$dosyalar" | grep -c . || true)
  [ "$sayi" -gt 0 ] || { echo "ERROR: no fastq was found: $KAYNAK"; exit 1; }
  echo "merging $sayi fastq files (at most ${KAP} reads per taxon; 0 = all of them)"
  local ks="$(dirname "$hedef")/kaynak_sayim.tsv"
  : > "$hedef"; : > "$ks"; : > /tmp/kraken_alinan_130.tsv
  while IFS= read -r f; do
    [ -n "$f" ] || continue
    local TX N ALINAN KALAN BU
    TX=$(basename "$f" | sed -E 's/.*reads_([0-9]+)\.fastq$/\1/')
    N=$(awk 'END{print int(NR/4)}' "$f")
    ALINAN=$(awk -F'\t' -v t="$TX" '$1==t{s+=$2} END{print s+0}' /tmp/kraken_alinan_130.tsv)
    if [ "$KAP" -gt 0 ]; then
      KALAN=$((KAP - ALINAN)); [ "$KALAN" -le 0 ] && KALAN=0
      BU=$(( N < KALAN ? N : KALAN ))
    else BU=$N; fi
    printf '%s\t%s\t%s\t%s\n' "$TX" "$N" "$BU" "$(basename "$f")" >> "$ks"
    printf '%s\t%s\n' "$TX" "$BU" >> /tmp/kraken_alinan_130.tsv
    [ "$BU" -gt 0 ] || continue
    awk -v tx="$TX" -v lim="$BU" \
      'NR%4==1 {k++; if (k>lim) exit; sub(/^@/, "@tx" tx "_")} {print}' "$f" >> "$hedef"
  done <<< "$dosyalar"
  local toplam beklenen
  toplam=$(awk 'END{print int(NR/4)}' "$hedef")
  beklenen=$(awk -F'\t' '{s+=$3} END{print s+0}' "$ks")
  [ "$toplam" = "$beklenen" ] || { echo "ERROR: the read count does not add up ($toplam / $beklenen)"; exit 1; }
  echo "  $toplam reads merged, and the count verified"
}


# ---------------------------------------------------------------------------
# ornekle builds a REPRESENTATIVE subset for the threshold scan.
#
# WHY
# The point of the threshold scan is to see THE SHAPE OF THE CURVE: which
# assignments collapse as the threshold rises. Reclassifying every read at every
# threshold is not needed for that. Six thresholds x 330 thousand reads, against a
# 110 GB database, was locking Windows up through the page cache under WSL2.
#
# The subset is taken EQUALLY PER BIN, so that abundant taxa do not gain a
# disproportionate weight. The percentages come out of ratios; the ratios stay the
# same in a representative sample and only the absolute numbers shrink.
#
# CONFIRMING THE REPRESENTATIVENESS is a separate job and not this function's:
# the 'dogrula-ornek' key compares the sample against the full data at a single
# threshold. If the ratios hold, the curve is valid. If they do not, ORNEK is
# raised.
# ---------------------------------------------------------------------------
ornekle() {
  local kaynak="$1" hedef="$2" toplam="$3"
  if [ -s "$hedef" ]; then
    echo "  the sample file already exists: $(awk 'END{print int(NR/4)}' "$hedef") reads"
    return
  fi
  local n_kaynak; n_kaynak=$(awk 'END{print int(NR/4)}' "$kaynak")
  if [ "$toplam" -le 0 ] || [ "$n_kaynak" -le "$toplam" ]; then
    echo "  sampling WAS NOT NEEDED (source $n_kaynak reads, target $toplam)"
    cp "$kaynak" "$hedef"
    return
  fi

  # --- ONEK DENETIMI (guard) --------------------------------------------
  # Sharing out per bin rests on the 'tx<taxid>_' prefix birlestir() puts on the
  # read headers. If birlestir() reused an old tum.fastq and that file has no
  # prefix, then EVERY read counts as one bin and the sample is not
  # representative. That would be a silent failure, so it is measured first.
  local onekli
  onekli=$(awk 'NR%4==1{n++; if ($1 ~ /^@tx[0-9]+_/) k++} NR>40000{exit} END{print (n>0? int(100*k/n) : 0)}' "$kaynak")
  if [ "${onekli:-0}" -lt 90 ]; then
    echo "  ERROR: the read headers carry no 'tx<taxid>_' prefix (${onekli:-0}%)."
    echo "  Equal sampling per bin is NOT POSSIBLE, and rather than quietly producing a"
    echo "  wrong sample the work stops. The fix: delete $(dirname "$kaynak")/tum.fastq and"
    echo "  have it regenerated, or run on the full data with ORNEK=0."
    exit 1
  fi

  local kutu_sayisi
  kutu_sayisi=$(awk 'NR%4==1{split(substr($1,2),a,"_"); print a[1]}' "$kaynak" | sort -u | wc -l)
  [ "$kutu_sayisi" -gt 0 ] || kutu_sayisi=1
  local pay=$(( toplam / kutu_sayisi ))
  [ "$pay" -lt 1 ] && pay=1

  # --- RASTGELE ORNEKLEME, SABIT TOHUMLA --------------------------------
  # WHY NOT THE FIRST N (the 2026-08-05 fix)
  # The previous version took the FIRST 'pay' reads of each bin. Over a nanopore
  # run the pore yield falls with time and the read quality changes systematically
  # towards the end; the head of the file differs from its tail. Taking the first N
  # produces a BIASED sample that over-represents the early part of the run, and
  # that bias can shift the threshold curve.
  # Every read is now scored with a hash value INDEPENDENT of its position inside
  # the bin, and the smallest 'pay' of them are chosen. The hash is produced from
  # the read name and a fixed seed, so the choice is RANDOM but REPRODUCIBLE: the
  # same input and the same seed always give the same sample. It can be changed
  # with ORNEK_TOHUM.
  local TOHUM="${ORNEK_TOHUM:-20260805}"
  echo "  sampling: $kutu_sayisi bins x $pay reads = target ~$toplam"
  echo "  method  : RANDOM within a bin, fixed seed $TOHUM (reproducible)"

  # Two passes: the hash of every read first, then a threshold per bin.
  # A name based hash is used instead of awk's own rand(), so that the same reads
  # are chosen even if the read order changes.
  awk -v tohum="$TOHUM" '
    function karma(s,   i,h,c) { h = tohum % 2147483647; for (i = length(s); i >= 1; i--) { c = index("_@0123456789abcdefghijklmnopqrstuvwxyzACGTN", substr(s,i,1)); h = (h * 16777619) % 2147483647; h = h + c * 2654435761; h = h % 2147483647; h = int(h / 7) + (h % 7) * 306783378 } return h % 1000003 }
    NR%4==1 { split(substr($1,2),a,"_"); print a[1] "\t" karma($1) }
  ' "$kaynak" | sort -k1,1 -k2,2n > /tmp/kraken_karma_130.tsv

  awk -v pay="$pay" '{ c[$1]++; if (c[$1]==pay) print $1 "\t" $2 }' /tmp/kraken_karma_130.tsv \
    > /tmp/kraken_esik_130.tsv

  awk -v tohum="$TOHUM" '
    function karma(s,   i,h,c) { h = tohum % 2147483647; for (i = length(s); i >= 1; i--) { c = index("_@0123456789abcdefghijklmnopqrstuvwxyzACGTN", substr(s,i,1)); h = (h * 16777619) % 2147483647; h = h + c * 2654435761; h = h % 2147483647; h = int(h / 7) + (h % 7) * 306783378 } return h % 1000003 }
    NR==FNR { esik[$1] = $2; next }
    FNR%4==1 { split(substr($1,2),a,"_"); tx=a[1]
               al = ((tx in esik) ? (karma($1) <= esik[tx]) : 1) }
    al { print }
  ' /tmp/kraken_esik_130.tsv "$kaynak" > "$hedef"

  local n; n=$(awk 'END{print int(NR/4)}' "$hedef")
  echo "  sample ready: $n reads (out of $n_kaynak)"
  printf 'kaynak_okuma\t%s\nornek_okuma\t%s\nkutu\t%s\npay\t%s\nyontem\trastgele\ntohum\t%s\n' \
    "$n_kaynak" "$n" "$kutu_sayisi" "$pay" "$TOHUM" > "$(dirname "$hedef")/ornek_bilgi.tsv"
  rm -f /tmp/kraken_karma_130.tsv /tmp/kraken_esik_130.tsv
}

# =========================================================================
# THE KEY: sure   A TIME ESTIMATE FROM A REAL MEASUREMENT
# No estimate is invented. A small sample is actually run, the time is measured and
# scaled to the full run. Under --memory-mapping the first reads are the slowest
# (the pages come from disk), so the measurement behaves like an UPPER bound, and
# that this is so is written on the screen.
# =========================================================================
tus_sure() {
  kraken_sart
  local d="${1:-$VT_A}"
  echo
  echo "TIMING MEASUREMENT. A small sample is actually run; no estimate is invented."
  vt_hazir_mi "$d" || { echo "The database is not ready, so no measurement can be made."; exit 1; }
  local BAYRAK; BAYRAK=$(bellek_bayragi "$d")
  local tmp; tmp=$(mktemp -d)
  birlestir "$tmp/tum.fastq" >/dev/null
  local TOPLAM; TOPLAM=$(awk 'END{print int(NR/4)}' "$tmp/tum.fastq")
  local N="${DENEME:-2000}"
  head -n $((N*4)) "$tmp/tum.fastq" > "$tmp/deneme.fastq"
  local gercek; gercek=$(awk 'END{print int(NR/4)}' "$tmp/deneme.fastq")
  echo
  echo "  trial: $gercek reads (the full set is $TOPLAM reads)"
  local t0 t1 sn
  t0=$(date +%s)
  kraken2 --db "$d" --threads "$IPLIK" $BAYRAK --confidence 0 \
          --report "$tmp/d.report" --output "$tmp/d.out" "$tmp/deneme.fastq" \
          >/dev/null 2>"$tmp/d.err" || {
    echo "  kraken2 ERROR:"; sed 's/^/    /' "$tmp/d.err" | tail -5; rm -rf "$tmp"; exit 1; }
  t1=$(date +%s); sn=$((t1 - t0)); [ "$sn" -lt 1 ] && sn=1
  local tam esik_adet
  tam=$(( sn * TOPLAM / gercek ))
  esik_adet=$(echo $ESIKLER | wc -w)
  echo
  echo "  measured : $sn seconds / $gercek reads"
  echo "  full set : about $(( tam / 60 )) minutes (a single threshold)"
  echo "  $esik_adet thresholds: about $(( tam * esik_adet / 60 )) minutes = $(( tam * esik_adet / 3600 )) hours"
  echo
  echo "  THE LIMIT OF THIS NUMBER: under --memory-mapping the first reads are the slowest,"
  echo "  because the pages come from disk. It speeds up as the operating system cache warms."
  echo "  So this estimate is an UPPER BOUND; the real time is equal to it or shorter."
  echo "  The one risk in the other direction: if the machine does other work and frees the"
  echo "  RAM, the cache cools and the time can grow. Keep the machine idle during the run."
  rm -rf "$tmp"
}

# =========================================================================
# TUS: durum
# =========================================================================
tus_durum() {
  echo "ENVIRONMENT CHECK. Nothing is run, only looked at."
  echo
  ortam_ac
  for a in kraken2 kraken2-build kraken2-inspect bracken; do
    if command -v "$a" >/dev/null 2>&1; then printf "  %-16s present  %s\n" "$a" "$(command -v $a)"
    else printf "  %-16s MISSING\n" "$a"; fi
  done
  echo
  echo "VT_A = $VT_A"
  if ! vt_hazir_mi "$VT_A"; then
    echo
    echo "VT_A cannot be used. Searching the disk (this is not given up on):"
    echo
    vt_ara || true
  fi
  if [ -n "$VT_B" ]; then
    echo
    echo "VT_B = $VT_B"
    vt_hazir_mi "$VT_B" || true
    if [ "$(readlink -f "$VT_A" 2>/dev/null)" = "$(readlink -f "$VT_B" 2>/dev/null)" ]; then
      echo
      echo "  WARNING: VT_A and VT_B are the SAME directory. Two curves cannot be drawn, only one."
    fi
  else
    echo
    echo "VT_B = (not given). Work will proceed with a single database."
    echo "  To compare two databases:  VT_B=/path/old_db bash $0 threshold"
  fi
  echo
  local n; n=$(ls "$KAYNAK"/*/*reads_*.fastq 2>/dev/null | wc -l || echo 0)
  echo "reads: $n fastq files ($KAYNAK)"
  echo
  free -g 2>/dev/null | awk '/Mem:/{print "RAM: "$2" GB total, "$7" GB available"}'
  df -h "$VT_A" 2>/dev/null | tail -1 | awk '{print "disk: "$4" bos ("$6")"}'
  echo
  eski_kosu_tespit
  echo
  echo "Next step:  bash $0 db-identity    (which version the database is)"
}

tus_vt_kimlik() {
  log_ac vt_kimlik
  kraken_sart
  echo
  echo "DATABASE VERSION DETECTION"
  echo "The whole argument rests on the COVERAGE of the database. If we run with the"
  echo "wrong version we read the result wrongly, so this step is not skipped."
  echo
  echo "VT_A = $VT_A"
  if vt_hazir_mi "$VT_A"; then
    echo
    vt_surum "$VT_A"
  fi
  if [ -n "$VT_B" ]; then
    echo; echo "VT_B = $VT_B"
    if vt_hazir_mi "$VT_B"; then echo; vt_surum "$VT_B"; fi
  fi
  echo
  eski_kosu_tespit
}

tus_sinav() {
  echo "SELF TESTS. No main work starts before the test passes (project rule 2)."
  local hata=0
  for p in threshold_summary.py comparison_table.py custom_taxonomy.py kraken_summary.py; do
    echo; echo "--- $p"
    if [ -f "$_BETIK_DIZIN/$p" ]; then
      python3 "$_BETIK_DIZIN/$p" --selftest || hata=1
    else echo "  NO SUCH FILE: $p"; hata=1; fi
  done
  echo
  [ "$hata" -eq 0 ] && echo "EVERY SELF TEST PASSED" || { echo "A SELF TEST FAILED"; return 1; }
}

# =========================================================================
# THE KEY: esik-a / esik-b   the confidence threshold scan
# =========================================================================
esik_tara() {
  local d="$1" is="$2" etiket="$3"
  echo
  echo "======================================================================"
  echo "CONFIDENCE THRESHOLD SCAN, $etiket"
  echo "  database: $d"
  echo "  thresholds: $ESIKLER"
  echo "  threads : $IPLIK   (kept low to hold the peak memory down)"
  echo "======================================================================"
  vt_hazir_mi "$d" || {
    echo
    echo "The database is not ready. This is not given up on; the disk is being searched."
    echo
    vt_ara || true
    echo
    echo "ERROR: '$d' could not be used, the work stops. Pick one of the paths above and give it."
    exit 1; }
  echo
  vt_surum "$d"
  local BAYRAK; BAYRAK=$(bellek_bayragi "$d")
  mkdir -p "$is"
  birlestir "$is/tum.fastq"
  # THE THRESHOLD SCAN RUNS ON THE SAMPLE. The reasoning is inside ornekle().
  # To run it on the full data: ORNEK=0
  local GIRDI="$is/tum.fastq"
  if [ "${ORNEK:-0}" -gt 0 ]; then
    echo
    echo "SAMPLING IS ON (ORNEK=$ORNEK). For the full data: ORNEK=0 bash $0 threshold"
    ornekle "$is/tum.fastq" "$is/ornek.fastq" "$ORNEK"
    GIRDI="$is/ornek.fastq"
  fi
  local OKUMA; OKUMA=$(awk 'END{print int(NR/4)}' "$GIRDI")
  echo "  reads to be used in the scan: $OKUMA"
  echo
  for C in $ESIKLER; do
    local ad="esik_${C}"
    if [ -s "$is/${ad}.report" ]; then
      echo "[$C] already exists, skipping"; continue
    fi
    local t0 t1
    t0=$(date +%s)
    echo "[$C] kraken2 --confidence $C  start $(date '+%H:%M:%S')"
    kraken2 --db "$d" --threads "$IPLIK" $BAYRAK --confidence "$C" \
            --use-names --report "$is/${ad}.report" --output "$is/${ad}.out" \
            "$GIRDI" 2>"$is/${ad}.err" >/dev/null || {
      echo "    kraken2 ERROR:"; tail -5 "$is/${ad}.err" | sed 's/^/      /'
      echo "    This threshold is not skipped, the work stops. A half scan misleads."; exit 1; }
    rm -f "$is/${ad}.err"
    t1=$(date +%s)
    local n; n=$(wc -l < "$is/${ad}.out")
    [ "$n" = "$OKUMA" ] || { echo "ERROR: the output has $n lines but $OKUMA reads, they do not match"; exit 1; }
    echo "    done, $n reads, $(( (t1-t0)/60 )) minutes $(( (t1-t0)%60 )) seconds"
  done
}

tus_esik() {
  log_ac esik
  kraken_sart
  local ikili=0
  if [ -n "$VT_B" ] && \
     [ "$(readlink -f "$VT_A" 2>/dev/null)" != "$(readlink -f "$VT_B" 2>/dev/null)" ]; then
    ikili=1
  fi
  if [ -n "$VT_B" ] && [ "$ikili" -eq 0 ]; then
    echo "WARNING: VT_A and VT_B are the same directory. Drawing two curves is"
    echo "meaningless, so a single scan will be made. If you really want to compare"
    echo "two different databases, give the path of the second: VT_B=/path/other bash $0 threshold"
  fi
  esik_tara "$VT_A" "$IS_A" "VT_A"
  if [ "$ikili" -eq 1 ]; then
    esik_tara "$VT_B" "$IS_B" "VT_B"
    echo
    python3 "$_BETIK_DIZIN/threshold_summary.py" --root "$PROJE" \
      --job "$IS_A" --name "$(basename "$VT_A")" \
      --is2 "$IS_B" --ad2 "$(basename "$VT_B")"
  else
    echo
    python3 "$_BETIK_DIZIN/threshold_summary.py" --root "$PROJE" \
      --job "$IS_A" --name "$(basename "$VT_A")"
  fi
  echo
  echo "done. Files: $IS_A"
  echo "  esik_<C>.report          the kraken report for each threshold"
  echo "  esik_egrisi.csv / .txt   assignment percentages per domain, against the threshold"
  # NOTE: "[ ... ] && echo" is never left as the last command. When the condition
  # is false the test returns 1 and that becomes the function's exit code; the work
  # finished successfully and the script still looks as though it FAILED. That trap
  # was fallen into once.
  if [ "$ikili" -eq 1 ]; then
    echo "  esik_iki_veritabani.txt  the two databases side by side, plus what survives"
  fi
  return 0
}

# ---------------------------------------------------------------------------
# tus_tablo - dort sutunlu karsilastirma tablosunu uretir.
#
# TABLO PARCA PARCA TAMAMLANABILIR (2026-08-04 notu)
# comparison_table.py eksik girdiye dayaniklidir: PlusPFP kosusu yoksa 3.
# column empty when there is no PlusPFP run, the 2nd when there is no threshold
# scan, and the 4th when there is no alignment result.
# birakir ve hangilerinin eksik oldugunu basar. Bu yuzden Kraken olcumu
# Even if the measurement is made weeks later, pressing this key again COMPLETES
# the table;
# bastan kosuya gerek yoktur.
#
# But there is a risk in that: with a COMPLETE table in hand, pressing this key at
# a moment when the database is temporarily unreachable would overwrite the table
# with an EMPTIER version.
# degistirebilirdi. Sessiz veri kaybi olurdu. Onlem: mevcut tablo her seferinde
# once yedeklenir. Yedek adi zaman damgalidir, ustune yazilmaz.
# ---------------------------------------------------------------------------
tus_tablo() {
  log_ac tablo
  local hedef="$PROJE/tools/0_TESLIM_RAPOR/KRAKEN_KARSILASTIRMA.md"
  mkdir -p "$(dirname "$hedef")"
  if [ -f "$hedef" ]; then
    local yed="$hedef.yedek_$(date +%Y%m%d_%H%M%S)"
    cp -p "$hedef" "$yed" 2>/dev/null && \
      echo "  the existing table was backed up: $(basename "$yed")"
  fi
  python3 "$_BETIK_DIZIN/comparison_table.py" --root "$PROJE" \
          --job-a "$IS_A" --name-a "$(basename "$VT_A")" \
          --job-b "$IS_B" --name-b "$([ -n "$VT_B" ] && basename "$VT_B" || echo '')" \
          --threshold "$TABLO_ESIK"
}

# =========================================================================
# THE KEY: ozelvt-kur   THE LAST RESORT
# =========================================================================
tus_ozelvt_kur() {
  echo "======================================================================"
  echo "BUILDING A CUSTOM KRAKEN2 DATABASE  (THE LAST RESORT)"
  echo "======================================================================"
  echo
  echo "READ THIS FIRST: if PlusPFP is installed, THIS STEP IS NOT NEEDED."
  echo "PlusPFP already adds protozoa, fungi and plants to Standard, which is"
  echo "every group we measured as missing. This key exists for the case where"
  echo "PlusPFP cannot be installed at all, or where a second opinion at marker"
  echo "gene level (16S/ITS) is wanted. The build takes hours."
  echo
  echo "Check:  bash $0 db-identity    (if PlusPFP is installed, do not come here at all)"
  echo
  local toplam=0 var=0
  for k in ${KUMELER:-silva_ssu unite pr2}; do
    local f; f=$(kume_dosya "$k") || continue
    if [ -s "$f" ]; then
      toplam=$((toplam + $(stat -c%s "$f"))); var=$((var+1))
      printf "  %-12s %-42s %s\n" "$k" "$(basename "$f")" "$(du -h "$f" | cut -f1)"
    else
      printf "  %-12s %-42s ABSENT\n" "$k" "$(basename "$f")"
    fi
  done
  [ "$var" -gt 0 ] || { echo; echo "ERROR: not one reference set was found ($PROJE/REFERANS_DB)"; exit 1; }
  local gb=$((toplam / 1073741824)); [ "$gb" -lt 1 ] && gb=1
  echo
  echo "  total sequence : about ${gb} GB"
  echo "  ESTIMATED TIME : $(( gb*40 + 30 )) to $(( gb*120 + 60 )) minutes"
  echo "  ESTIMATED RAM  : at least $(( gb*6 + 8 )) GB"
  echo "  ESTIMATED DISK : about $(( gb*5 + 10 )) GB"
  echo "  target         : $OZELVT"
  echo
  echo "  This is a MARKER GENE database (16S/18S/ITS); it holds no whole genomes."
  echo "  The taxonomy comes from the lineage strings in the files themselves, not"
  echo "  from NCBI. The taxids are synthetic and cannot be compared with NCBI"
  echo "  taxids, which is why the table compares at name level."
  local mem; mem=$(free -g 2>/dev/null | awk '/Mem:/{print $2}' || echo 0)
  echo "  total RAM on this machine: ${mem} GB"
  [ "${mem:-0}" -lt $(( gb*6 + 8 )) ] && echo "  WARNING: the RAM is below the estimated requirement."
  echo
  if [ "${ONAY:-}" != "evet" ]; then
    read -r -p "Kuruluma baslansin mi? (evet yazin, baska her sey iptal): " c
    [ "$c" = "evet" ] || { echo "cancelled. Nothing was done."; exit 0; }
  fi
  log_ac ozelvt_kur
  kraken_sart
  command -v kraken2-build >/dev/null 2>&1 || {
    echo "ERROR: kraken2-build was not found."
    echo "  micromamba install -n ${ORTAM:-mikro} -c bioconda kraken2"; exit 1; }
  mkdir -p "$OZELVT/library" "$OZELVT/taxonomy"
  echo; echo "1/2  building the taxonomy and library from the lineage strings"
  local args=()
  for k in ${KUMELER:-silva_ssu unite pr2}; do
    local f; f=$(kume_dosya "$k") || continue
    [ -s "$f" ] && args+=(--kume "$k=$f")
  done
  python3 "$_BETIK_DIZIN/custom_taxonomy.py" --output "$OZELVT" "${args[@]}"
  echo; echo "2/2  kraken2-build --build  ($(date '+%H:%M:%S'))"
  kraken2-build --build --db "$OZELVT" --threads "$IPLIK"
  echo; echo "build finished: $OZELVT"
  du -sh "$OZELVT" | awk '{print "  size: "$1}'
}

kume_dosya() {
  case "$1" in
    silva_ssu) echo "$PROJE/REFERANS_DB/SILVA_138.2_SSURef_NR99.fasta" ;;
    silva_lsu) echo "$PROJE/REFERANS_DB/SILVA_138.2_LSURef_NR99.fasta" ;;
    unite)     echo "$PROJE/REFERANS_DB/UNITE_ITS.fasta" ;;
    pr2)       echo "$PROJE/REFERANS_DB/PR2_SSU_taxo_long.fasta" ;;
    rod)       echo "$PROJE/REFERANS_DB/ROD_v1.2_operon_variants.fasta" ;;
    *) return 1 ;;
  esac
}

tus_ozelvt_kos() {
  log_ac ozelvt_kos
  kraken_sart
  [ -f "$OZELVT/hash.k2d" ] || {
    echo "ERROR: the custom database is not built ($OZELVT/hash.k2d is missing)."
    echo "  First:  bash $0 custom-db-build"; exit 1; }
  mkdir -p "$OZEL_IS"
  birlestir "$OZEL_IS/tum.fastq"
  local BAYRAK; BAYRAK=$(bellek_bayragi "$OZELVT")
  local C="${OZEL_ESIK:-0}"
  kraken2 --db "$OZELVT" --threads "$IPLIK" $BAYRAK --confidence "$C" \
          --use-names --report "$OZEL_IS/tum.report" --output "$OZEL_IS/tum.out" \
          "$OZEL_IS/tum.fastq" 2>"$OZEL_IS/hata.txt" >/dev/null || {
    echo "kraken2 ERROR:"; tail -5 "$OZEL_IS/hata.txt" | sed 's/^/  /'; exit 1; }
  rm -f "$OZEL_IS/hata.txt"
  python3 "$_BETIK_DIZIN/kraken_summary.py" --job "$OZEL_IS" --toolkit "$_BETIK_DIZIN" || true
}


# ---------------------------------------------------------------------------
# tus_bellek_ayari produces the .wslconfig that stops WSL2 from choking Windows.
#
# THE PROBLEM
# WSL2 is a separate virtual machine. By default it can take a large part of the
# Windows RAM, AND when a 110 GB database is read with --memory-mapping the page
# cache counts towards the virtual machine's memory too. The result: WSL grows
# without bound, Windows starts swapping and the machine freezes.
#
# THE FIX
# An UPPER BOUND is put ON THE VIRTUAL MACHINE with .wslconfig. Windows is left
# some room to breathe. WSL slows down but the machine stays usable. That is not a
# trade-off but the right configuration: an unbounded WSL is not running fast for
# anyone, it is locking up.
# ---------------------------------------------------------------------------
tus_bellek_ayari() {
  echo "======================================================================"
  echo "WSL MEMORY SETTING - it stops Windows from freezing"
  echo "======================================================================"
  echo
  local TOPLAM_MB=0 KAYNAK="could not be measured"
  if command -v powershell.exe >/dev/null 2>&1; then
    TOPLAM_MB=$(powershell.exe -NoProfile -Command \
      "[math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory/1MB)" \
      2>/dev/null | tr -d '\r' | tr -d ' ')
    [ -n "$TOPLAM_MB" ] && KAYNAK="measured from Windows"
  fi
  case "$TOPLAM_MB" in ''|*[!0-9]*) TOPLAM_MB=0 ;; esac
  if [ "$TOPLAM_MB" -le 0 ]; then
    TOPLAM_MB=$(( $(awk '/MemTotal/{print $2}' /proc/meminfo) / 1024 ))
    KAYNAK="estimated from inside WSL (the Windows total can be LARGER)"
  fi
  local TOPLAM_GB=$(( TOPLAM_MB / 1024 ))
  # 60 percent: enough for WSL and it leaves Windows room to breathe. The floor is 4 GB,
  # because a 110 GB database needs a floor of memory even with mmap.
  local WSL_GB=$(( TOPLAM_GB * 60 / 100 ))
  [ "$WSL_GB" -lt 4 ] && WSL_GB=4
  local SWAP_GB=$(( WSL_GB / 2 ))
  [ "$SWAP_GB" -lt 4 ] && SWAP_GB=4

  echo "  total RAM       : ${TOPLAM_GB} GB   ($KAYNAK)"
  echo "  suggested for WSL: ${WSL_GB} GB    (60 percent of the total)"
  echo "  swap            : ${SWAP_GB} GB"
  echo

  local ICERIK
  ICERIK="[wsl2]
memory=${WSL_GB}GB
swap=${SWAP_GB}GB
processors=${IPLIK}
# Caps the page cache, so that reading a 110 GB database through mmap does
# not let the cache inflate the virtual machine.
pageReporting=true
# Gives idle memory back to Windows.
autoMemoryReclaim=gradual"

  local CIKTI="$PROJE/wslconfig_ONERILEN.txt"
  printf '%s\n' "$ICERIK" > "$CIKTI"

  echo "----------------------------------------------------------------------"
  printf '%s\n' "$ICERIK"
  echo "----------------------------------------------------------------------"
  echo
  echo "WHAT TO DO - three steps, in order"
  echo
  echo "  1) The text above was written to this file:"
  echo "       $CIKTI"
  echo "     Its name on the Windows side is that same file seen through \\\\wsl$."
  echo
  echo "  2) Copy that file to THIS NAME and THIS PLACE:"
  echo "       %USERPROFILE%\\.wslconfig"
  echo "     (the file name starts with a dot and has NO extension)"
  echo "     A shortcut: open Run on Windows (Win+R) and paste this:"
  echo "       notepad %USERPROFILE%\\.wslconfig"
  echo "     If an empty file opens, paste the text and save."
  echo
  echo "  3) Shut WSL down and start it again. In PowerShell or Command Prompt:"
  echo "       wsl --shutdown"
  echo "     Then open WSL again. The setting takes effect only after that."
  echo
  echo "VERIFICATION: after reopening it, run this command inside WSL"
  echo "       free -g"
  echo "  the 'total' column should read about ${WSL_GB} GB."
  echo
  echo "THIS SETTING ON ITS OWN FIXES THE FREEZE. WSL gets slower and Windows stays"
  echo "usable. The threshold scan was also lightened with sampling (see below)."
}

# ---------------------------------------------------------------------------
# tus_dogrula_ornek measures THAT THE SAMPLE IS REPRESENTATIVE.
#
# At a single threshold it classifies both the sample and the FULL data and
# compares THE PERCENTAGES. The absolute numbers will of course differ; what is
# being looked at is whether the ratios hold. If they do, the threshold curve is
# valid.
#
# CAREFUL: this key runs once on the FULL data, so it is heavy. Leave it overnight.
# The threshold curve can be used without running it, but then the report says
# "the representativeness WAS NOT CONFIRMED", and it has to be presented that
# way.
# ---------------------------------------------------------------------------
tus_dogrula_ornek() {
  log_ac dogrula_ornek
  kraken_sart
  local d="$VT_A" is="$IS_A" C="${DOGRULAMA_ESIGI:-0.05}"
  vt_hazir_mi "$d" || { echo "ERROR: the database is not ready"; exit 1; }
  mkdir -p "$is"
  birlestir "$is/tum.fastq"
  ornekle "$is/tum.fastq" "$is/ornek.fastq" "$ORNEK"
  local BAYRAK; BAYRAK=$(bellek_bayragi "$d")
  local a="$is/dogrulama_ornek" b="$is/dogrulama_tam"
  for cift in "ornek.fastq:$a" "tum.fastq:$b"; do
    local gir="${cift%%:*}" cik="${cift##*:}"
    if [ -s "${cik}.report" ]; then echo "$(basename $cik) already exists, skipping"; continue; fi
    echo "[$C] classifying $(basename $gir)  $(date '+%H:%M:%S')"
    kraken2 --db "$d" --threads "$IPLIK" $BAYRAK --confidence "$C" \
            --use-names --report "${cik}.report" --output /dev/null \
            "$is/$gir" 2>"${cik}.err" >/dev/null || {
      echo "  kraken2 ERROR:"; tail -3 "${cik}.err" | sed 's/^/    /'; exit 1; }
    rm -f "${cik}.err"
  done
  echo
  echo "REPRESENTATION COMPARISON  (threshold $C)"
  echo "  What is looked at: domain level percentages. The absolute counts will differ, and that DOES NOT MATTER."
  echo
  printf '  %-24s %10s %10s %10s\n' "domain" "sample %" "full %" "difference"
  local sapma_max=0
  for tx in 2157:Archaea 2:Bacteria 4751:Fungi 33090:Viridiplantae 0:sinifsiz; do
    local id="${tx%%:*}" ad="${tx##*:}"
    local yo yt
    if [ "$id" = "0" ]; then
      yo=$(awk -F'\t' '$5=="0"{print $1+0; exit}' "${a}.report"); yo=${yo:-0}
      yt=$(awk -F'\t' '$5=="0"{print $1+0; exit}' "${b}.report"); yt=${yt:-0}
    else
      yo=$(awk -F'\t' -v i="$id" '$5==i{print $1+0; exit}' "${a}.report"); yo=${yo:-0}
      yt=$(awk -F'\t' -v i="$id" '$5==i{print $1+0; exit}' "${b}.report"); yt=${yt:-0}
    fi
    local fark; fark=$(awk -v x="$yo" -v y="$yt" 'BEGIN{printf "%.2f", (x>y?x-y:y-x)}')
    printf '  %-24s %10s %10s %10s\n' "$ad" "$yo" "$yt" "$fark"
    sapma_max=$(awk -v a="$sapma_max" -v b="$fark" 'BEGIN{print (b>a?b:a)}')
  done
  echo
  echo "  en buyuk sapma: $sapma_max puan"
  awk -v s="$sapma_max" 'BEGIN{
    if (s <= 2.0) {
      print "  SONUC: ORNEK TEMSILIDIR. Sapma 2 puanin altinda, esik egrisi gecerli."
    } else {
      print "  SONUC: SAPMA BUYUK. Ornek temsili degil.";
      print "  ORNEK degerini buyutup tekrarlayin, orn: ORNEK=300000 bash kraken_tool.sh threshold"
    }
  }'
  printf 'esik\t%s\nen_buyuk_sapma_puan\t%s\n' "$C" "$sapma_max" \
    > "$is/temsil_dogrulamasi.tsv"
  echo
  echo "  written: $is/temsil_dogrulamasi.tsv"
}

# =========================================================================
case "$TUS" in
  memory-config) tus_bellek_ayari ;;
  verify-sample) tus_dogrula_ornek ;;
  kraken-path)  tus_kraken_yol ;;
  status)       tus_durum ;;
  find-db)      ortam_ac; vt_ara || true ;;
  original-db)      tus_ali_vt ;;
  db-identity)   tus_vt_kimlik ;;
  selftest)       tus_sinav ;;
  time)        tus_sure "${2:-$VT_A}" ;;
  threshold-a)      tus_sinav >/dev/null || { echo "A SELF TEST FAILED. Detail: bash $0 selftest"; exit 2; }
               log_ac esik_a; kraken_sart; esik_tara "$VT_A" "$IS_A" "VT_A"
               python3 "$_BETIK_DIZIN/threshold_summary.py" --root "$PROJE" --job "$IS_A" --name "$(basename "$VT_A")" ;;
  threshold-b)      [ -n "$VT_B" ] || { echo "ERROR: VT_B was not given.  VT_B=/path bash $0 threshold-b"; exit 1; }
               tus_sinav >/dev/null || { echo "SINAV KALDI"; exit 2; }
               log_ac esik_b; kraken_sart; esik_tara "$VT_B" "$IS_B" "VT_B"
               python3 "$_BETIK_DIZIN/threshold_summary.py" --root "$PROJE" --job "$IS_B" --name "$(basename "$VT_B")" ;;
  threshold)        tus_sinav >/dev/null || { echo "A SELF TEST FAILED. Detail: bash $0 selftest"; exit 2; }
               tus_esik ;;
  table)       tus_tablo ;;
  all)       tus_sinav >/dev/null || { echo "SINAV KALDI"; exit 2; }
               tus_esik; tus_tablo ;;
  custom-db-build)  tus_ozelvt_kur ;;
  custom-db-run)  tus_ozelvt_kos ;;
  *)           sed -n '3,48p' "$0" | sed 's/^# \{0,1\}//'
               echo; echo "For detail: $PROJE/docs/GUIDE.md" ;;
esac
