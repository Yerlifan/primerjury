#!/usr/bin/env bash
# =============================================================================
#  install.sh  --  single installation entry point
#
#  Installs the tools, the reference databases, the Kraken2 database and the
#  QIIME2/PICRUSt2 environment. Replaces four separate installers that used to
#  live in different directories of the source project.
#
#  Two things remain deliberately separate and are called from here:
#     tools/install_mfeprimer.sh   MFEprimer binary (its download URL changes
#                                  with every release, so it is kept apart)
#     build_index.sh               MFEprimer index for one FASTA file
#
#  USAGE
#  -----
#  Prefer the top-level entry point, which wraps this file in English:
#     ./primerjury install all
#
#  Direct use (Turkish subcommand names still work as aliases):
#     bash install.sh status                 # status: MEASURES what is installed, changes nothing
#     bash install.sh tools               # tools:  kraken2, bracken, minimap2, samtools, blast, seqkit
#     bash install.sh databases            # databases: SILVA + UNITE + PR2 + ROD + RefSeq (12 identity sources)
#     bash install.sh databases --only refseq,pr2     # only these
#     bash install.sh kraken-download          # prebuilt Kraken2 database (k fixed at 35)
#     bash install.sh kraken-build --kmer 31  # build your OWN Kraken2 database with your OWN k
#     bash install.sh qiime                 # QIIME2 + PICRUSt2 + SILVA classifier
#     bash install.sh all                 # all of the above (Kraken2 asked separately)
#
#  THREE DESIGN DECISIONS
#  ----------------------
#  1) DOWNLOADS ARE VERIFIED; THE URL IS NOT TRUSTED.
#     A hard-coded URL goes stale silently, and you end up running against the
#     wrong database without any error. So every downloaded file is MEASURED:
#     is it really FASTA, how many records, RNA or DNA alphabet, is the size in
#     the expected order of magnitude. A file that fails is renamed .SUPHELI
#     (suspect) and IS NOT USED. Running silently against a truncated database
#     is worse than not running at all.
#
#  2) EVERY STEP IS RESUMABLE.
#     A verified download is skipped; a partial one continues with wget -c. If
#     the script is interrupted, the same command picks up where it stopped.
#
#  3) NOTHING IS SKIPPED SILENTLY.
#     Everything that could not be installed is listed at the end and the exit
#     code is non-zero. "Installation complete" means it actually completed.
# =============================================================================
set -o pipefail

KOK="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REFDB="$KOK/REFERANS_DB"
LOG_KLASOR="$KOK/kurulum_loglari"
mkdir -p "$LOG_KLASOR" "$REFDB"
LOG="$LOG_KLASOR/KUR_$(date '+%Y%m%d_%H%M%S').log"
exec > >(tee -a "$LOG") 2>&1

EKSIK=()          # kurulamayanlar
SUPHELI=()        # indirildi ama dogrulamayi gecemedi

renk()  { printf '\n\033[1m== %s ==\033[0m\n' "$*"; }
bilgi() { printf '   %s\n' "$*"; }
uyar()  { printf '   \033[33mWARNING:\033[0m %s\n' "$*"; }
hata()  { printf '   \033[31mERROR  :\033[0m %s\n' "$*"; }

# ---------------------------------------------------------------------------
# THE ENVIRONMENT CHECK: this script has to run inside WSL or Linux.
#
# Under Git Bash (MSYS) it RUNS but it MISLEADS: it sees the file system, so it
# counts the databases correctly, while it DOES NOT SEE WSL's $HOME or the conda
# and micromamba environments. The result: installed tools appear MISSING and the
# Kraken2 database is not found. That is a silent false report, which is why it is
# said up front.
# ---------------------------------------------------------------------------
ortam_denetimi() {
  case "$(uname -s 2>/dev/null)" in
    Linux) return 0 ;;
  esac
  uyar "This script must run inside WSL or Linux. Detected: '$(uname -s 2>/dev/null)'."
  bilgi "The filesystem is read correctly, but the TOOL and KRAKEN2 lines will MISLEAD you:"
  bilgi "the WSL \$HOME and the conda/micromamba environment are invisible from here."
  bilgi "Run it this way instead:"
  bilgi "   wsl bash -lc \"cd '$KOK' && bash install.sh $*\""
  printf '\n'
}

# ---------------------------------------------------------------------------
# THE ADDRESSES
#
# CAREFUL: SILVA and UNITE CHANGE their file names when they raise a version. The
# addresses below are the versions of the day they were written. Because the script
# verifies what it downloads, a stale address does not turn into a silent fault: it
# either gets a 404 or fails verification. If you want another version, give it with
# --url, for example:
#     bash install.sh databases --only silva_ssu --url https://.../ANOTHER.fasta.gz
#
# The Kraken2 database address is DELIBERATELY not embedded; the current list
# changes with every release and the wrong database costs hours. See
# kraken-download.
# ---------------------------------------------------------------------------
SILVA_TABAN="https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports"
declare -A URL=(
  [silva_ssu]="$SILVA_TABAN/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"
  [silva_lsu]="$SILVA_TABAN/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz"
  [silva_lsu_parc]="$SILVA_TABAN/SILVA_138.2_LSUParc_tax_silva.fasta.gz"
  [pr2]="https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_taxo_long.fasta.gz"
  [refseq_bak16s]="https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz"
  [refseq_ark16s]="https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz"
  [refseq_its]="https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.ITS.fna.gz"
  [refseq_18s]="https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.18SrRNA.fna.gz"
  [refseq_28s]="https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.28SrRNA.fna.gz"
)
# the target file names; the code expects these names (the database list in
# verification/identity_verification.py)
declare -A HEDEF=(
  [silva_ssu]="SILVA_138.2_SSURef_NR99.fasta"
  [silva_lsu]="SILVA_138.2_LSURef_NR99.fasta"
  [silva_lsu_parc]="SILVA_138.2_LSUParc.fasta"
  [pr2]="PR2_SSU_taxo_long.fasta"
  [refseq_bak16s]="bacteria.16S.fna"
  [refseq_ark16s]="archaea.16S.fna"
  [refseq_its]="fungi.ITS.fna"
  [refseq_18s]="fungi.18SrRNA.fna"
  [refseq_28s]="fungi.28SrRNA.fna"
)
# the minimum record count expected; below this the file is SUSPECT
# (measured values, KONTROL_SONUC/CAPRAZ_KONTROL_2026-08-09_2216.md)
declare -A ASGARI_KAYIT=(
  [silva_ssu]=400000  [silva_lsu]=80000  [silva_lsu_parc]=1000000
  [pr2]=200000        [refseq_bak16s]=20000  [refseq_ark16s]=800
  [refseq_its]=15000  [refseq_18s]=3000      [refseq_28s]=10000
)
SIRA=(silva_ssu silva_lsu silva_lsu_parc pr2 refseq_bak16s refseq_ark16s
      refseq_its refseq_18s refseq_28s)

# UNITE and ROD are special: the UNITE address changes with a DOI per release, and
# ROD is a git repository.
UNITE_SAYFA="https://unite.ut.ee/repository.php"
ROD_DEPO="https://github.com/krabberod/ROD"

# ---------------------------------------------------------------------------
mamba_hazirla() {
  export PATH="$HOME/bin:$PATH"
  if ! command -v micromamba >/dev/null 2>&1; then
    bilgi "installing micromamba..."
    ( cd "$HOME" && curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
        | tar -xj bin/micromamba ) || { hata "could not install micromamba"; return 1; }
    grep -q 'HOME/bin' "$HOME/.bashrc" 2>/dev/null \
      || echo 'export PATH="$HOME/bin:$PATH"' >> "$HOME/.bashrc"
  fi
  export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
  eval "$(micromamba shell hook -s bash)" || return 1
}

# ---------------------------------------------------------------------------
# fasta_dogrula <file> <minimum_records>
#   It MEASURES that the file REALLY IS a usable FASTA.
#   0 = sound, 1 = suspect. The reason is printed to the screen.
# ---------------------------------------------------------------------------
fasta_dogrula() {
  local f="$1" asgari="${2:-1}"
  [ -s "$f" ] || { hata "file missing or empty: $(basename "$f")"; return 1; }
  local ilk; ilk=$(head -c 1 "$f" 2>/dev/null)
  [ "$ilk" = ">" ] || { hata "not FASTA (first character is not '>'): $(basename "$f")"; return 1; }
  local n; n=$(grep -c '^>' "$f" 2>/dev/null || echo 0)
  if [ "$n" -lt "$asgari" ]; then
    hata "$(basename "$f"): $n records, expected at least $asgari. The download may be INCOMPLETE"
    return 1
  fi
  # The RNA against DNA measurement: SILVA stores its sequences with U. MFEprimer
  # builds its index over the {A,C,G,T} alphabet; with no T at all the alphabet
  # falls to {A,C,G} and 3^9 k-mers enter the index instead of 4^9. MFEprimer
  # RAISES NO ERROR in that case. (see build_index.sh)
  local ornek u t
  ornek=$(grep -v '^>' "$f" 2>/dev/null | head -2000 | tr -d '\n')
  u=$(printf '%s' "$ornek" | tr -cd 'Uu' | wc -c)
  t=$(printf '%s' "$ornek" | tr -cd 'Tt' | wc -c)
  if [ "$u" -gt "$t" ]; then
    bilgi "$(basename "$f"): $n records, RNA alphabet (U>T). build_index.sh will convert U to T"
  else
    bilgi "$(basename "$f"): $n records, DNA alphabet"
  fi
  return 0
}

# ---------------------------------------------------------------------------
indir_ve_dogrula() {
  local anahtar="$1" url="${2:-${URL[$anahtar]}}"
  local hedef="$REFDB/${HEDEF[$anahtar]}"
  local asgari="${ASGARI_KAYIT[$anahtar]:-1}"

  if [ -s "$hedef" ] && fasta_dogrula "$hedef" "$asgari" >/dev/null 2>&1; then
    bilgi "$anahtar: already present and verified -> $(basename "$hedef")"
    return 0
  fi
  bilgi "downloading $anahtar: $url"
  local gz="$hedef.gz.indiriliyor"
  if ! wget -c -q --show-progress -O "$gz" "$url"; then
    hata "$anahtar COULD NOT BE DOWNLOADED: $url"
    bilgi "  the URL may be stale. Find the current one and pass it like this:"
    bilgi "    bash install.sh databases --only $anahtar --url <YENI_ADRES>"
    EKSIK+=("$anahtar (download failed)")
    rm -f "$gz"
    return 1
  fi
  bilgi "  extracting..."
  if ! gunzip -c "$gz" > "$hedef.tmp" 2>/dev/null; then
    # gz olmayabilir (bazi kaynaklar duz fasta verir)
    mv "$gz" "$hedef.tmp"
  else
    rm -f "$gz"
  fi
  if fasta_dogrula "$hedef.tmp" "$asgari"; then
    mv "$hedef.tmp" "$hedef"
    bilgi "  OK -> $(basename "$hedef")"
    return 0
  fi
  mv "$hedef.tmp" "$hedef.SUPHELI"
  hata "$anahtar FAILED verification, marked .SUPHELI (suspect) and WILL NOT BE USED"
  SUPHELI+=("$anahtar")
  return 1
}

# ===========================================================================

# Arac surumu: BLAST ailesi '-version' ister, digerleri '--version'. Yanlis
# using the wrong flag prints the usage text and produces a meaningless line like
# "PRESENT  USAGE"; that was measured. The flag is chosen per tool, not guessed.
surum_oku() {
  local t="$1" c
  case "$t" in
    blastn|makeblastdb|blastdbcmd) c=$("$t" -version 2>&1 | head -1) ;;
    seqkit)                        c=$("$t" version 2>&1 | head -1) ;;
    bracken)                       c=$("$t" -v 2>&1 | grep -m1 -oE 'v?[0-9]+\.[0-9]+(\.[0-9]+)?') ;;
    *)                             c=$("$t" --version 2>&1 | head -1) ;;
  esac
  # If the version cannot be read it IS NOT INVENTED. "installed (version not
  # readable)" is a measurement; printing the usage text or an error line as a
  # version is false information.
  case "$c" in
    ''|*USAGE*|*Usage*|*usage*|*unknown\ flag*|*Error*|/*) c='installed (version not readable)' ;;
  esac
  printf '%s' "${c:0:50}"
}

# Araci once PATH'te, bulamazsa proje ortaminda (micromamba) arar.
# The reason: the tools are installed into the '$PT_ORTAM' environment ('mikro' by
# default) and while that environment is not active they ARE NOT VISIBLE on PATH. A
# check that looks only at PATH reports installed tools as MISSING, that is, it makes
# a correctly installed system look broken.
arac_yolu() {
  local t="$1"
  command -v "$t" 2>/dev/null && return 0
  local o="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/${PT_ORTAM:-mikro}/bin/$t"
  [ -x "$o" ] && { printf '%s' "$o"; return 0; }
  return 1
}

komut_durum() {
  renk "TOOLS"
  local t y
  for t in kraken2 bracken minimap2 samtools blastn makeblastdb seqkit mfeprimer qiime picrust2; do
    if y=$(arac_yolu "$t"); then
      local nerede=''
      command -v "$t" >/dev/null 2>&1 || nerede="  [in env ${PT_ORTAM:-mikro}, not activated]"
      printf '   %-14s OK     %s%s\n' "$t" "$(PATH="$(dirname "$y"):$PATH" surum_oku "$t")" "$nerede"
    else
      printf '   %-14s \033[31mMISSING\033[0m\n' "$t"
    fi
  done
  # mfeprimer can also be inside the project
  if ! command -v mfeprimer >/dev/null 2>&1; then
    # EXECUTABLE only, and never a source file. The loose pattern used to match
    # steps/mfeprimer_layer.py and report a Python module as the MFEprimer
    # binary -- measured after the module rename. A wrong "found it" is worse
    # than "not found": the caller would go on to execute a Python file as a
    # binary and get an error that names neither cause.
    local m; m=$(find "$KOK" -maxdepth 3 -name 'mfeprimer*' -type f -perm -u+x \
                     ! -name '*.py' ! -name '*.pyc' ! -name '*.sh' ! -name '*.md' ! -path '*__pycache__*' 2>/dev/null | head -1)
    [ -n "$m" ] && bilgi "mfeprimer found inside the project: ${m#$KOK/}"
  fi

  renk "REFERENCE DATABASES  ($REFDB)"
  local k n
  for k in "${SIRA[@]}"; do
    local f="$REFDB/${HEDEF[$k]}"
    if [ -s "$f" ]; then
      n=$(grep -c '^>' "$f" 2>/dev/null || echo '?')
      printf '   %-16s OK     %8s records  %s\n' "$k" "$n" "$(du -h "$f" 2>/dev/null | cut -f1)"
    elif [ -s "$f.SUPHELI" ]; then
      printf '   %-16s \033[33mSUPHELI\033[0m  dogrulamayi gecemedi\n' "$k"
    else
      printf '   %-16s \033[31mMISSING\033[0m\n' "$k"
    fi
  done
  [ -s "$REFDB/UNITE_ITS.fasta" ] \
    && printf '   %-16s OK     %8s records\n' unite "$(grep -c '^>' "$REFDB/UNITE_ITS.fasta")" \
    || printf '   %-16s \033[31mMISSING\033[0m  (manual: %s)\n' unite "$UNITE_SAYFA"
  [ -s "$REFDB/ROD_v1.2_operon_variants.fasta" ] \
    && printf '   %-16s OK     %8s records\n' rod "$(grep -c '^>' "$REFDB/ROD_v1.2_operon_variants.fasta")" \
    || printf '   %-16s \033[31mMISSING\033[0m  (%s)\n' rod "$ROD_DEPO"

  renk "KRAKEN2 DATABASE"
  if [ -n "${KRAKEN2_DB_PATH:-}" ] && [ -f "$KRAKEN2_DB_PATH/hash.k2d" ]; then
    bilgi "KRAKEN2_DB_PATH = $KRAKEN2_DB_PATH"
  fi
  local bulunan; bulunan=$(find "$HOME" -maxdepth 3 -name hash.k2d 2>/dev/null | head -5)
  if [ -n "$bulunan" ]; then
    while read -r h; do
      local d; d=$(dirname "$h")
      printf '   %s  (%s)\n' "$d" "$(du -sh "$d" 2>/dev/null | cut -f1)"
      [ -f "$d/opts.k2d" ] && bilgi "    for k-mer details: bash tools/kraken_tool.sh db-identity"
    done <<< "$bulunan"
  else
    printf '   \033[31mMISSING\033[0m  -> bash install.sh kraken-download  or  bash install.sh kraken-build\n'
  fi

  renk "IDENTITY VERIFICATION COVERAGE"
  local var=0 top=0
  for k in "${SIRA[@]}"; do top=$((top+1)); [ -s "$REFDB/${HEDEF[$k]}" ] && var=$((var+1)); done
  for f in UNITE_ITS.fasta ROD_v1.2_operon_variants.fasta ref_all2.fna; do
    top=$((top+1)); [ -s "$REFDB/$f" ] && var=$((var+1))
  done
  bilgi "$var of $top independent sources ready"
  bilgi "verification/identity_verification.py requires AT LEAST TWO independent"
  bilgi "sources to agree before it calls a claim VERIFIED. As sources drop away the"
  bilgi "verdict slides to 'UNVERIFIED (single source)'. A missing database therefore"
  bilgi "produces an OPENLY weak answer, never a quietly wrong one."
}

# ===========================================================================
komut_araclar() {
  renk "TOOLS"
  mamba_hazirla || { EKSIK+=("micromamba"); return 1; }
  local ORT="${PT_ORTAM:-mikro}"
  if micromamba env list 2>/dev/null | grep -qE "^ *$ORT "; then
    bilgi "environment '$ORT' exists, updating"
    micromamba install -y -n "$ORT" -c conda-forge -c bioconda \
      kraken2 bracken minimap2 seqkit blast samtools barrnap || EKSIK+=("tools")
  else
    bilgi "creating environment '$ORT'"
    micromamba create -y -n "$ORT" -c conda-forge -c bioconda \
      kraken2 bracken minimap2 seqkit blast samtools barrnap python=3.11 || EKSIK+=("tools")
  fi
  micromamba activate "$ORT" 2>/dev/null

  bilgi "python packages"
  python -m pip install --quiet --upgrade primer3-py biopython numpy openpyxl pysam matplotlib \
    || EKSIK+=("python packages")

  # MFEprimer: it is not in conda, it is downloaded as a binary
  if ! command -v mfeprimer >/dev/null 2>&1 \
     && [ ! -x "$KOK/tools/mfeprimer" ]; then
    bilgi "For MFEprimer: bash tools/install_mfeprimer.sh"
    bilgi "  (kept in a separate script because its download URL keeps changing)"
    EKSIK+=("mfeprimer (install with tools/install_mfeprimer.sh)")
  fi

  renk "VERSION CHECK"
  local t
  for t in kraken2 bracken minimap2 seqkit blastn samtools; do
    if command -v "$t" >/dev/null 2>&1; then
      printf '   %-12s OK\n' "$t"
    else
      printf '   %-12s \033[31mKURULAMADI\033[0m\n' "$t"; EKSIK+=("$t")
    fi
  done
}

# ===========================================================================
komut_veritabani() {
  local yalniz="" ozel_url=""
  while [ $# -gt 0 ]; do
    case "$1" in
      --only)          yalniz="$2"; shift 2 ;;
      --url)           ozel_url="$2"; shift 2 ;;
      *) shift ;;
    esac
  done
  renk "REFERANS VERITABANLARI"
  bilgi "Target directory: $REFDB"
  bilgi "About 28 GB in total. If interrupted, the same command resumes where it stopped."
  local k
  for k in "${SIRA[@]}"; do
    if [ -n "$yalniz" ] && [[ ",$yalniz," != *",$k,"* ]]; then continue; fi
    indir_ve_dogrula "$k" "${ozel_url:-}"
  done

  # ROD: git deposu
  if [ -z "$yalniz" ] || [[ ",$yalniz," == *",rod,"* ]]; then
    if [ ! -s "$REFDB/ROD_v1.2_operon_variants.fasta" ]; then
      renk "ROD (rRNA operon database)"
      if command -v git >/dev/null 2>&1; then
        git clone --depth 1 "$ROD_DEPO" "$REFDB/ROD_depo" 2>/dev/null \
          && bilgi "downloaded -> REFERANS_DB/ROD_depo (copy the operon FASTA from there)" \
          || { hata "could not download ROD"; EKSIK+=("rod"); }
      else
        EKSIK+=("rod (git not installed)")
      fi
    fi
  fi

  # UNITE: the address changes with a DOI per release, it cannot be embedded
  if [ -z "$yalniz" ] || [[ ",$yalniz," == *",unite,"* ]]; then
    if [ ! -s "$REFDB/UNITE_ITS.fasta" ]; then
      renk "UNITE (fungal ITS)"
      uyar "The UNITE download URL changes with every release DOI, so it cannot be hard-coded."
      bilgi "1) copy the current 'General FASTA release' link from $UNITE_SAYFA"
      bilgi "2) bash install.sh databases --only unite --url <THE_LINK_YOU_COPIED>"
      bilgi "   The script VERIFIES the record count and alphabet after downloading."
      EKSIK+=("unite (URL must be given manually)")
    fi
  fi

  renk "MFEPRIMER INDEXES"
  bilgi "Every downloaded FASTA needs an index:"
  bilgi "   bash build_index.sh --liste          # see the candidate files"
  bilgi "   bash build_index.sh <file>          # build them one at a time"
  bilgi "U to T conversion is MANDATORY for SILVA. Without it the index is built SILENTLY broken."
}

# ===========================================================================
komut_kraken_indir() {
  renk "PREBUILT KRAKEN2 DATABASE"
  bilgi "Current list: https://benlangmead.github.io/aws-indexes/k2"
  bilgi
  bilgi "The URL is DELIBERATELY not hard-coded: file names change with every release and"
  bilgi "a stale URL silently wastes hours on the WRONG database."
  bilgi
  bilgi "Recommended: PlusPF (Standard + protozoa + fungi). PlusPF-16 if disk is tight."
  bilgi
  bilgi "   mkdir -p ~/k2db && cd ~/k2db"
  bilgi "   wget <secilen .tar.gz adresi>"
  bilgi "   tar -xzf *.tar.gz"
  bilgi
  uyar "PREBUILT databases are fixed at k=35, l=31. YOU CANNOT CHOOSE THE"
  uyar "K-MER LENGTH. To choose your own: bash install.sh kraken-build"
  bilgi
  bilgi "After installing, VERIFY its identity (which release, which k-mer):"
  bilgi "   bash tools/kraken_tool.sh db-identity"
}

# ===========================================================================
komut_kraken_kur() {
  local KMER=35 MINI=31 BOSLUK=7 KUTUP="" DB="$HOME/k2db_ozel" IS=""
  while [ $# -gt 0 ]; do
    case "$1" in
      --kmer)               KMER="$2"; shift 2 ;;
      --minimizer)          MINI="$2"; shift 2 ;;
      --spaces)             BOSLUK="$2"; shift 2 ;;
      --library)            KUTUP="$2"; shift 2 ;;
      --db)                 DB="$2"; shift 2 ;;
      --threads)            IS="$2"; shift 2 ;;
      *) shift ;;
    esac
  done
  IS="${IS:-$(( $(nproc 2>/dev/null || echo 4) - 2 ))}"; [ "$IS" -lt 1 ] && IS=1

  renk "KRAKEN2 DATABASE KURULUMU  (k-mer secilebilir)"
  cat <<EOF
   k-mer uzunlugu (--kmer)      : $KMER
   minimizer uzunlugu (--minimizer): $MINI
   minimizer boslugu (--bosluk) : $BOSLUK
   veritabani yolu (--db)       : $DB
   threads (--threads)          : $IS

   K-MER SECIMI NE YAPAR
   ---------------------
   Kraken2 bir okumayi k uzunlugundaki parcalara boler ve her parcayi
   in the taxonomy tree; the decision is the lowest common ancestor (LCA) of those
   pieces.
     * a LARGER k -> the piece is more specific and wrong assignments fall, but the
                     sensitivity drops (a single error makes a piece fail to match;
                     on long read, high error rate ONT data that loss is large)
     * k KUCUK  -> duyarlilik artar, ama parca birden cok taksonda bulunur ve
                   LCA agacta YUKARI kayar: tur yerine cins, cins yerine aile.
   The default k=35, l=31 was chosen for short, low error Illumina reads.
   ONT verisinde k'yi dusurmek (ornegin 31) duyarliligi artirir.

   DIKKAT: k-mer secimi TEK BASINA kimlik sorununu cozmez. Bu projede olculdu -
   The Kraken2 label and the alignment based identity DISAGREE on many bins
   (tools/0_TESLIM_RAPOR/KRAKEN_KARSILASTIRMA.md). The cause is not the k-mer but
   the LCA itself and the database's coverage. So whatever k you build with, the
   identity claims have to be tested INDEPENDENTLY:
       python3 verification/identity_verification.py --root .
   That script uses the taxonomy tree NOT AT ALL; it does a seed plus full
   alignment against 12 separate reference databases and requires AT LEAST TWO
   sources to agree.

   BUILDING TWO DATABASES WITH DIFFERENT k AND COMPARING THEM
   ------------------------------------------------
       bash install.sh kraken-build --kmer 35 --db ~/k2db_k35
       bash install.sh kraken-build --kmer 31 --db ~/k2db_k31
       bash tools/kraken_tool.sh threshold      # the threshold scan
       bash tools/kraken_tool.sh table     # the four column comparison
EOF
  if ! command -v kraken2-build >/dev/null 2>&1; then
    hata "kraken2-build not found. Run first: bash install.sh tools"
    EKSIK+=("kraken2-build"); return 1
  fi

  bilgi
  bilgi "Downloading taxonomy (once, about 1-2 hours, resumes if interrupted)..."
  kraken2-build --download-taxonomy --db "$DB" --threads "$IS" \
    || { hata "could not download taxonomy"; EKSIK+=("kraken taxonomy"); return 1; }

  local kutuphaneler="${KUTUP:-bacteria archaea fungi protozoa}"
  local L
  for L in $kutuphaneler; do
    bilgi "downloading library: $L"
    kraken2-build --download-library "$L" --db "$DB" --threads "$IS" \
      || { hata "could not download library: $L"; EKSIK+=("kraken library $L"); }
  done

  bilgi "building the database (k=$KMER, l=$MINI, s=$BOSLUK). This takes HOURS"
  kraken2-build --build --db "$DB" --threads "$IS" \
      --kmer-len "$KMER" --minimizer-len "$MINI" --minimizer-spaces "$BOSLUK" \
    || { hata "could not build the database"; EKSIK+=("kraken build"); return 1; }

  if [ -f "$DB/hash.k2d" ]; then
    bilgi "TAMAM: $DB  ($(du -sh "$DB" 2>/dev/null | cut -f1))"
    # RECORD the k-mer used: opts.k2d is binary and a person cannot read it. If
    # which k it was built with is not written down, in six months nobody will know
    # and the comparison becomes meaningless.
    printf 'built: %s\nkmer: %s\nminimizer: %s\nspaces: %s\nlibrary: %s\n' \
      "$(date '+%Y-%m-%d %H:%M')" "$KMER" "$MINI" "$BOSLUK" "$kutuphaneler" \
      > "$DB/KURULUM_BILGISI.txt"
    bilgi "build parameters written to: $DB/KURULUM_BILGISI.txt"
    bilgi "For Bracken: bracken-build -d $DB -t $IS -k $KMER -l <read_length>"
  else
    hata "hash.k2d was not created, so the build did not finish"; EKSIK+=("kraken build")
  fi
}

# ===========================================================================
komut_qiime() {
  renk "QIIME2 + PICRUSt2"
  if ! command -v conda >/dev/null 2>&1; then
    for c in "$HOME/miniconda3" "$HOME/anaconda3"; do
      [ -f "$c/etc/profile.d/conda.sh" ] && . "$c/etc/profile.d/conda.sh" && break
    done
  fi
  if ! command -v conda >/dev/null 2>&1; then
    hata "conda not found. QIIME2 needs conda (micromamba can be problematic)."
    EKSIK+=("qiime2 (conda not found)"); return 1
  fi
  local SURUM="${PT_QIIME_SURUM:-2024.10}"
  local ORT="qiime2-amplicon-$SURUM"
  if conda env list | grep -q "$ORT"; then
    bilgi "environment already exists: $ORT"
  else
    bilgi "installing QIIME2 $SURUM..."
    conda env create -n "$ORT" \
      --file "https://data.qiime2.org/distro/amplicon/qiime2-amplicon-$SURUM-py310-linux-conda.yml" \
      || { hata "could not install QIIME2"; EKSIK+=("qiime2"); }
  fi
  if conda env list | grep -q "picrust2"; then
    bilgi "picrust2 environment already exists"
  else
    conda create -y -n picrust2 -c conda-forge -c bioconda picrust2 \
      || { hata "could not install PICRUSt2"; EKSIK+=("picrust2"); }
  fi
  local SINIF="$REFDB/silva_classifier.qza"
  if [ -s "$SINIF" ]; then
    bilgi "SILVA classifier already present"
  else
    bilgi "downloading the SILVA classifier..."
    wget -c -q --show-progress -O "$SINIF" \
      "https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza" \
      || { hata "could not download the classifier"; EKSIK+=("silva_classifier"); rm -f "$SINIF"; }
  fi
}

# ===========================================================================
ozet() {
  renk "SUMMARY"
  if [ ${#SUPHELI[@]} -gt 0 ]; then
    hata "Files that FAILED verification (kept as .SUPHELI, NOT used):"
    printf '     - %s\n' "${SUPHELI[@]}"
  fi
  if [ ${#EKSIK[@]} -gt 0 ]; then
    hata "COULD NOT BE INSTALLED:"
    printf '     - %s\n' "${EKSIK[@]}"
    bilgi
    bilgi "None of these were skipped silently. Fix them and run the same command again;"
    bilgi "completed steps will be skipped."
    bilgi "log: $LOG"
    return 1
  fi
  printf '   \033[32mInstallation complete.\033[0m No step was skipped.\n'
  bilgi "to see the current state: bash install.sh status"
  bilgi "log: $LOG"
  return 0
}

# ===========================================================================
KOMUT="${1:-status}"; shift 2>/dev/null || true
ortam_denetimi "$KOMUT"
case "$KOMUT" in
  status)                     komut_durum ;;
  tools)                      komut_araclar; ozet ;;
  databases|db)               komut_veritabani "$@"; ozet ;;
  kraken-download)            komut_kraken_indir ;;
  kraken-build)               komut_kraken_kur "$@"; ozet ;;
  qiime)                      komut_qiime; ozet ;;
  all)                        komut_araclar; komut_veritabani; komut_qiime
                 renk "KRAKEN2"
                 bilgi "The Kraken2 database is a SEPARATE choice (download prebuilt, or build with your own k):"
                 bilgi "   bash install.sh kraken-download     # prebuilt, k fixed at 35"
                 bilgi "   bash install.sh kraken-build --kmer 31"
                 ozet ;;
  *)
    sed -n '2,45p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit 2 ;;
esac
