#!/usr/bin/env bash
# =============================================================================
#  install.sh  --  TEK KURULUM KAPISI
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
#  Direct use (subcommands are still Turkish; being translated):
#     bash install.sh durum                 # status: MEASURES what is installed, changes nothing
#     bash install.sh araclar               # tools:  kraken2, bracken, minimap2, samtools, blast, seqkit
#     bash install.sh veritabani            # databases: SILVA + UNITE + PR2 + ROD + RefSeq (12 identity sources)
#     bash install.sh veritabani --yalniz refseq,pr2     # only these
#     bash install.sh kraken-indir          # prebuilt Kraken2 database (k fixed at 35)
#     bash install.sh kraken-kur --kmer 31  # build your OWN Kraken2 database with your OWN k
#     bash install.sh qiime                 # QIIME2 + PICRUSt2 + SILVA classifier
#     bash install.sh hepsi                 # all of the above (Kraken2 asked separately)
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
uyar()  { printf '   \033[33mUYARI:\033[0m %s\n' "$*"; }
hata()  { printf '   \033[31mHATA :\033[0m %s\n' "$*"; }

# ---------------------------------------------------------------------------
# ORTAM DENETIMI: bu betik WSL/Linux icinde kosmalidir.
#
# Git Bash (MSYS) altinda da CALISIR ama YANILTIR: dosya sistemini gorur, yani
# veritabanlarini dogru sayar; buna karsilik WSL'in $HOME'unu ve conda/micromamba
# ortamini GORMEZ. Sonuc: kurulu araclar "YOK" gorunur, Kraken2 veritabani
# bulunamaz. Bu sessiz bir yanlis rapordur - o yuzden basta soyleniyor.
# ---------------------------------------------------------------------------
ortam_denetimi() {
  case "$(uname -s 2>/dev/null)" in
    Linux) return 0 ;;
  esac
  uyar "Bu betik WSL/Linux icinde kosmalidir; su an '$(uname -s 2>/dev/null)' altindasiniz."
  bilgi "Dosya sistemi dogru okunur ama ARAC ve KRAKEN2 satirlari YANILTIR:"
  bilgi "WSL'in \$HOME'u ve conda/micromamba ortami buradan gorunmez."
  bilgi "Dogrusu:"
  bilgi "   wsl bash -lc \"cd '$KOK' && bash install.sh $*\""
  printf '\n'
}

# ---------------------------------------------------------------------------
# ADRESLER
#
# DIKKAT: SILVA ve UNITE surum yukseltince dosya adlarini DEGISTIRIR. Asagidaki
# adresler yazildigi gunku surumlerdir. Betik indirdigi dosyayi dogruladigi icin
# eskimis bir adres sessiz hataya donusmez: ya 404 alir ya dogrulamayi geceremez.
# Baska bir surum istiyorsaniz --url ile verin, ornek:
#     bash install.sh veritabani --yalniz silva_ssu --url https://.../BASKA.fasta.gz
#
# Kraken2 veritabani adresi BILEREK gomulu degildir; guncel liste her surumde
# degisir ve yanlis veritabani saatler harcatir. Bkz. kraken-indir.
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
# hedef dosya adlari - kod bu adlari bekliyor (verification/identity_verification.py VTB listesi)
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
# asgari beklenen kayit sayisi - bunun altindaysa dosya SUPHELI
# (olculen degerler, KONTROL_SONUC/CAPRAZ_KONTROL_2026-08-09_2216.md)
declare -A ASGARI_KAYIT=(
  [silva_ssu]=400000  [silva_lsu]=80000  [silva_lsu_parc]=1000000
  [pr2]=200000        [refseq_bak16s]=20000  [refseq_ark16s]=800
  [refseq_its]=15000  [refseq_18s]=3000      [refseq_28s]=10000
)
SIRA=(silva_ssu silva_lsu silva_lsu_parc pr2 refseq_bak16s refseq_ark16s
      refseq_its refseq_18s refseq_28s)

# UNITE ve ROD ozel: UNITE adresi surum basina DOI ile degisir, ROD bir git deposudur.
UNITE_SAYFA="https://unite.ut.ee/repository.php"
ROD_DEPO="https://github.com/krabberod/ROD"

# ---------------------------------------------------------------------------
mamba_hazirla() {
  export PATH="$HOME/bin:$PATH"
  if ! command -v micromamba >/dev/null 2>&1; then
    bilgi "micromamba kuruluyor..."
    ( cd "$HOME" && curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
        | tar -xj bin/micromamba ) || { hata "micromamba kurulamadi"; return 1; }
    grep -q 'HOME/bin' "$HOME/.bashrc" 2>/dev/null \
      || echo 'export PATH="$HOME/bin:$PATH"' >> "$HOME/.bashrc"
  fi
  export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
  eval "$(micromamba shell hook -s bash)" || return 1
}

# ---------------------------------------------------------------------------
# fasta_dogrula <dosya> <asgari_kayit>
#   Dosyanin GERCEKTEN kullanilabilir bir FASTA oldugunu OLCER.
#   0 = saglikli, 1 = supheli. Sebep ekrana yazilir.
# ---------------------------------------------------------------------------
fasta_dogrula() {
  local f="$1" asgari="${2:-1}"
  [ -s "$f" ] || { hata "dosya yok ya da bos: $(basename "$f")"; return 1; }
  local ilk; ilk=$(head -c 1 "$f" 2>/dev/null)
  [ "$ilk" = ">" ] || { hata "FASTA degil (ilk karakter '>' degil): $(basename "$f")"; return 1; }
  local n; n=$(grep -c '^>' "$f" 2>/dev/null || echo 0)
  if [ "$n" -lt "$asgari" ]; then
    hata "$(basename "$f"): $n kayit, beklenen en az $asgari - dosya EKSIK indirilmis olabilir"
    return 1
  fi
  # RNA/DNA olcumu: SILVA dizileri U ile saklar. MFEprimer indeksi {A,C,G,T}
  # alfabesi kurar; hic T yoksa alfabe {A,C,G}'ye duser ve indekse 4^9 yerine
  # 3^9 k-mer girer. MFEprimer bu durumda HATA VERMEZ. (bkz. build_index.sh)
  local ornek u t
  ornek=$(grep -v '^>' "$f" 2>/dev/null | head -2000 | tr -d '\n')
  u=$(printf '%s' "$ornek" | tr -cd 'Uu' | wc -c)
  t=$(printf '%s' "$ornek" | tr -cd 'Tt' | wc -c)
  if [ "$u" -gt "$t" ]; then
    bilgi "$(basename "$f"): $n kayit, RNA alfabesi (U>T) - build_index.sh U->T donusumu yapacak"
  else
    bilgi "$(basename "$f"): $n kayit, DNA alfabesi"
  fi
  return 0
}

# ---------------------------------------------------------------------------
indir_ve_dogrula() {
  local anahtar="$1" url="${2:-${URL[$anahtar]}}"
  local hedef="$REFDB/${HEDEF[$anahtar]}"
  local asgari="${ASGARI_KAYIT[$anahtar]:-1}"

  if [ -s "$hedef" ] && fasta_dogrula "$hedef" "$asgari" >/dev/null 2>&1; then
    bilgi "$anahtar: zaten var ve dogrulandi -> $(basename "$hedef")"
    return 0
  fi
  bilgi "$anahtar indiriliyor: $url"
  local gz="$hedef.gz.indiriliyor"
  if ! wget -c -q --show-progress -O "$gz" "$url"; then
    hata "$anahtar INDIRILEMEDI: $url"
    bilgi "  adres eskimis olabilir. Guncelini bulup su sekilde verin:"
    bilgi "    bash install.sh veritabani --yalniz $anahtar --url <YENI_ADRES>"
    EKSIK+=("$anahtar (indirilemedi)")
    rm -f "$gz"
    return 1
  fi
  bilgi "  aciliyor..."
  if ! gunzip -c "$gz" > "$hedef.tmp" 2>/dev/null; then
    # gz olmayabilir (bazi kaynaklar duz fasta verir)
    mv "$gz" "$hedef.tmp"
  else
    rm -f "$gz"
  fi
  if fasta_dogrula "$hedef.tmp" "$asgari"; then
    mv "$hedef.tmp" "$hedef"
    bilgi "  TAMAM -> $(basename "$hedef")"
    return 0
  fi
  mv "$hedef.tmp" "$hedef.SUPHELI"
  hata "$anahtar dogrulamayi GECEMEDI, .SUPHELI olarak isaretlendi - KULLANILMAYACAK"
  SUPHELI+=("$anahtar")
  return 1
}

# ===========================================================================

# Arac surumu: BLAST ailesi '-version' ister, digerleri '--version'. Yanlis
# bayrak kullanmak kullanim metnini basar ve "VAR  USAGE" gibi anlamsiz bir
# satir uretir - olculdu. Bayrak araca gore secilir, tahmin edilmez.
surum_oku() {
  local t="$1" c
  case "$t" in
    blastn|makeblastdb|blastdbcmd) c=$("$t" -version 2>&1 | head -1) ;;
    seqkit)                        c=$("$t" version 2>&1 | head -1) ;;
    bracken)                       c=$("$t" -v 2>&1 | grep -m1 -oE 'v?[0-9]+\.[0-9]+(\.[0-9]+)?') ;;
    *)                             c=$("$t" --version 2>&1 | head -1) ;;
  esac
  # Surum okunamadiysa UYDURULMAZ. "kurulu (surum okunamadi)" bir olcumdur;
  # kullanim metnini ya da hata satirini surum diye basmak yanlis bilgidir.
  case "$c" in
    ''|*USAGE*|*Usage*|*usage*|*unknown\ flag*|*Error*|/*) c='kurulu (surum okunamadi)' ;;
  esac
  printf '%s' "${c:0:50}"
}

# Araci once PATH'te, bulamazsa proje ortaminda (micromamba) arar.
# Sebep: araclar '$PT_ORTAM' (varsayilan 'mikro') ortamina kuruluyor ve o ortam
# etkin degilken PATH'te GORUNMEZ. Yalniz PATH'e bakan bir denetim, kurulu
# araclari "YOK" diye raporlar - yani dogru kurulmus bir sistemi bozuk gosterir.
arac_yolu() {
  local t="$1"
  command -v "$t" 2>/dev/null && return 0
  local o="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/${PT_ORTAM:-mikro}/bin/$t"
  [ -x "$o" ] && { printf '%s' "$o"; return 0; }
  return 1
}

komut_durum() {
  renk "tools"
  local t y
  for t in kraken2 bracken minimap2 samtools blastn makeblastdb seqkit mfeprimer qiime picrust2; do
    if y=$(arac_yolu "$t"); then
      local nerede=''
      command -v "$t" >/dev/null 2>&1 || nerede="  [ortam: ${PT_ORTAM:-mikro}, etkin degil]"
      printf '   %-14s VAR    %s%s\n' "$t" "$(PATH="$(dirname "$y"):$PATH" surum_oku "$t")" "$nerede"
    else
      printf '   %-14s \033[31mYOK\033[0m\n' "$t"
    fi
  done
  # mfeprimer proje icinde de olabilir
  if ! command -v mfeprimer >/dev/null 2>&1; then
    # EXECUTABLE only, and never a source file. The loose pattern used to match
    # steps/mfeprimer_layer.py and report a Python module as the MFEprimer
    # binary -- measured after the module rename. A wrong "found it" is worse
    # than "not found": the caller would go on to execute a Python file as a
    # binary and get an error that names neither cause.
    local m; m=$(find "$KOK" -maxdepth 3 -name 'mfeprimer*' -type f -perm -u+x \
                     ! -name '*.py' ! -name '*.pyc' ! -name '*.sh' ! -name '*.md' ! -path '*__pycache__*' 2>/dev/null | head -1)
    [ -n "$m" ] && bilgi "mfeprimer proje icinde bulundu: ${m#$KOK/}"
  fi

  renk "REFERANS VERITABANLARI  ($REFDB)"
  local k n
  for k in "${SIRA[@]}"; do
    local f="$REFDB/${HEDEF[$k]}"
    if [ -s "$f" ]; then
      n=$(grep -c '^>' "$f" 2>/dev/null || echo '?')
      printf '   %-16s VAR    %8s kayit  %s\n' "$k" "$n" "$(du -h "$f" 2>/dev/null | cut -f1)"
    elif [ -s "$f.SUPHELI" ]; then
      printf '   %-16s \033[33mSUPHELI\033[0m  dogrulamayi gecemedi\n' "$k"
    else
      printf '   %-16s \033[31mYOK\033[0m\n' "$k"
    fi
  done
  [ -s "$REFDB/UNITE_ITS.fasta" ] \
    && printf '   %-16s VAR    %8s kayit\n' unite "$(grep -c '^>' "$REFDB/UNITE_ITS.fasta")" \
    || printf '   %-16s \033[31mYOK\033[0m  (elle: %s)\n' unite "$UNITE_SAYFA"
  [ -s "$REFDB/ROD_v1.2_operon_variants.fasta" ] \
    && printf '   %-16s VAR    %8s kayit\n' rod "$(grep -c '^>' "$REFDB/ROD_v1.2_operon_variants.fasta")" \
    || printf '   %-16s \033[31mYOK\033[0m  (%s)\n' rod "$ROD_DEPO"

  renk "KRAKEN2 VERITABANI"
  if [ -n "${KRAKEN2_DB_PATH:-}" ] && [ -f "$KRAKEN2_DB_PATH/hash.k2d" ]; then
    bilgi "KRAKEN2_DB_PATH = $KRAKEN2_DB_PATH"
  fi
  local bulunan; bulunan=$(find "$HOME" -maxdepth 3 -name hash.k2d 2>/dev/null | head -5)
  if [ -n "$bulunan" ]; then
    while read -r h; do
      local d; d=$(dirname "$h")
      printf '   %s  (%s)\n' "$d" "$(du -sh "$d" 2>/dev/null | cut -f1)"
      [ -f "$d/opts.k2d" ] && bilgi "    k-mer bilgisi icin: bash tools/kraken_tool.sh vt-kimlik"
    done <<< "$bulunan"
  else
    printf '   \033[31mYOK\033[0m  -> bash install.sh kraken-indir   ya da   bash install.sh kraken-kur\n'
  fi

  renk "KIMLIK DOGRULAMA KAPSAMI"
  local var=0 top=0
  for k in "${SIRA[@]}"; do top=$((top+1)); [ -s "$REFDB/${HEDEF[$k]}" ] && var=$((var+1)); done
  for f in UNITE_ITS.fasta ROD_v1.2_operon_variants.fasta ref_all2.fna; do
    top=$((top+1)); [ -s "$REFDB/$f" ] && var=$((var+1))
  done
  bilgi "$var / $top bagimsiz kaynak hazir"
  bilgi "verification/identity_verification.py bir iddiayi DOGRULANDI saymak icin EN AZ IKI"
  bilgi "bagimsiz kaynagin uyusmasini sart kosar. Kaynak sayisi dustukce hukum"
  bilgi "'DOGRULANAMADI (tek kaynak)'a kayar - yani eksik veritabani sessizce"
  bilgi "yanlis cevap degil, ACIKCA zayif cevap uretir."
}

# ===========================================================================
komut_araclar() {
  renk "tools"
  mamba_hazirla || { EKSIK+=("micromamba"); return 1; }
  local ORT="${PT_ORTAM:-mikro}"
  if micromamba env list 2>/dev/null | grep -qE "^ *$ORT "; then
    bilgi "ortam '$ORT' var, guncelleniyor"
    micromamba install -y -n "$ORT" -c conda-forge -c bioconda \
      kraken2 bracken minimap2 seqkit blast samtools barrnap || EKSIK+=("araclar")
  else
    bilgi "ortam '$ORT' kuruluyor"
    micromamba create -y -n "$ORT" -c conda-forge -c bioconda \
      kraken2 bracken minimap2 seqkit blast samtools barrnap python=3.11 || EKSIK+=("araclar")
  fi
  micromamba activate "$ORT" 2>/dev/null

  bilgi "python paketleri"
  python -m pip install --quiet --upgrade primer3-py biopython numpy openpyxl pysam matplotlib \
    || EKSIK+=("python paketleri")

  # MFEprimer: conda'da yok, ikili olarak indirilir
  if ! command -v mfeprimer >/dev/null 2>&1 \
     && [ ! -x "$KOK/tools/mfeprimer" ]; then
    bilgi "MFEprimer icin: bash tools/install_mfeprimer.sh"
    bilgi "  (surum adresi surekli degistigi icin ayri betikte tutuluyor)"
    EKSIK+=("mfeprimer (tools/install_mfeprimer.sh ile kurun)")
  fi

  renk "SURUM DOGRULAMASI"
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
      --yalniz) yalniz="$2"; shift 2 ;;
      --url)    ozel_url="$2"; shift 2 ;;
      *) shift ;;
    esac
  done
  renk "REFERANS VERITABANLARI"
  bilgi "Hedef klasor: $REFDB"
  bilgi "Toplam ~28 GB. Kesilirse ayni komut kaldigi yerden devam eder."
  local k
  for k in "${SIRA[@]}"; do
    if [ -n "$yalniz" ] && [[ ",$yalniz," != *",$k,"* ]]; then continue; fi
    indir_ve_dogrula "$k" "${ozel_url:-}"
  done

  # ROD: git deposu
  if [ -z "$yalniz" ] || [[ ",$yalniz," == *",rod,"* ]]; then
    if [ ! -s "$REFDB/ROD_v1.2_operon_variants.fasta" ]; then
      renk "ROD (rRNA operon veritabani)"
      if command -v git >/dev/null 2>&1; then
        git clone --depth 1 "$ROD_DEPO" "$REFDB/ROD_depo" 2>/dev/null \
          && bilgi "indirildi -> REFERANS_DB/ROD_depo (operon fasta'sini oradan kopyalayin)" \
          || { hata "ROD indirilemedi"; EKSIK+=("rod"); }
      else
        EKSIK+=("rod (git yok)")
      fi
    fi
  fi

  # UNITE: adres surum basina DOI ile degisir, gomulemez
  if [ -z "$yalniz" ] || [[ ",$yalniz," == *",unite,"* ]]; then
    if [ ! -s "$REFDB/UNITE_ITS.fasta" ]; then
      renk "UNITE (mantar ITS)"
      uyar "UNITE indirme adresi her surumde DOI ile degisir; koda gomulemez."
      bilgi "1) $UNITE_SAYFA adresinden guncel 'General FASTA release' baglantisini kopyalayin"
      bilgi "2) bash install.sh veritabani --yalniz unite --url <KOPYALADIGINIZ_ADRES>"
      bilgi "   Betik indirdikten sonra kayit sayisini ve alfabesini DOGRULAR."
      EKSIK+=("unite (adres elle verilmeli)")
    fi
  fi

  renk "MFEPRIMER INDEKSLERI"
  bilgi "Indirilen her FASTA icin indeks kurulmali:"
  bilgi "   bash build_index.sh --liste          # aday dosyalari gorun"
  bilgi "   bash build_index.sh <dosya_adi>      # tek tek kurun"
  bilgi "SILVA dosyalarinda U->T donusumu SART - indeks aksi halde SESSIZCE bozuk kurulur."
}

# ===========================================================================
komut_kraken_indir() {
  renk "KRAKEN2 HAZIR VERITABANI"
  bilgi "Guncel liste: https://benlangmead.github.io/aws-indexes/k2"
  bilgi
  bilgi "Adres BILEREK koda gomulu degildir: dosya adlari her surumde degisir ve"
  bilgi "gomulu bir adres sessizce eskiyerek YANLIS veritabaniyla saatler harcatir."
  bilgi
  bilgi "Onerilen: PlusPF (Standard + protozoa + fungi). Disk yetmezse PlusPF-16."
  bilgi
  bilgi "   mkdir -p ~/k2db && cd ~/k2db"
  bilgi "   wget <secilen .tar.gz adresi>"
  bilgi "   tar -xzf *.tar.gz"
  bilgi
  uyar "HAZIR veritabanlari k=35, l=31 ile kurulmustur. K-MER UZUNLUGUNU"
  uyar "SECEMEZSINIZ. Kendi k-mer'inizi istiyorsaniz: bash install.sh kraken-kur"
  bilgi
  bilgi "Kurduktan sonra kimligini DOGRULAYIN (hangi surum, hangi k-mer):"
  bilgi "   bash tools/kraken_tool.sh vt-kimlik"
}

# ===========================================================================
komut_kraken_kur() {
  local KMER=35 MINI=31 BOSLUK=7 KUTUP="" DB="$HOME/k2db_ozel" IS=""
  while [ $# -gt 0 ]; do
    case "$1" in
      --kmer)      KMER="$2"; shift 2 ;;
      --minimizer) MINI="$2"; shift 2 ;;
      --bosluk)    BOSLUK="$2"; shift 2 ;;
      --kutuphane) KUTUP="$2"; shift 2 ;;
      --db)        DB="$2"; shift 2 ;;
      --is)        IS="$2"; shift 2 ;;
      *) shift ;;
    esac
  done
  IS="${IS:-$(( $(nproc 2>/dev/null || echo 4) - 2 ))}"; [ "$IS" -lt 1 ] && IS=1

  renk "KRAKEN2 VERITABANI KURULUMU  (k-mer secilebilir)"
  cat <<EOF
   k-mer uzunlugu (--kmer)      : $KMER
   minimizer uzunlugu (--minimizer): $MINI
   minimizer boslugu (--bosluk) : $BOSLUK
   veritabani yolu (--db)       : $DB
   is parcacigi (--is)          : $IS

   K-MER SECIMI NE YAPAR
   ---------------------
   Kraken2 bir okumayi k uzunlugundaki parcalara boler ve her parcayi
   taksonomi agacinda arar; karar bu parcalarin en kucuk ortak atasidir (LCA).
     * k BUYUK  -> parca daha ozgul, yanlis atama azalir, ama duyarlilik duser
                   (tek bir hata parcayi eslesmez kilar; uzun okuma / yuksek
                   hata oranli ONT verisinde bu kayip buyuktur)
     * k KUCUK  -> duyarlilik artar, ama parca birden cok taksonda bulunur ve
                   LCA agacta YUKARI kayar: tur yerine cins, cins yerine aile.
   Varsayilan k=35, l=31 kisa ve dusuk hatali Illumina okumasi icin secilmistir.
   ONT verisinde k'yi dusurmek (ornegin 31) duyarliligi artirir.

   DIKKAT: k-mer secimi TEK BASINA kimlik sorununu cozmez. Bu projede olculdu -
   Kraken2 etiketi ile hizalama tabanli kimlik bircok kutuda AYRISIYOR
   (tools/0_TESLIM_RAPOR/KRAKEN_KARSILASTIRMA.md). Sebep k-mer degil, LCA'nin
   kendisi ve veritabani kapsamidir. Bu yuzden hangi k ile kurarsaniz kurun,
   kimlik iddialari BAGIMSIZ olarak sinanmalidir:
       python3 verification/identity_verification.py --kok .
   O betik taksonomi agacini HIC kullanmaz, 12 ayri referans veritabaninda
   tohum + tam hizalama yapar ve EN AZ IKI kaynagin uyusmasini sart kosar.

   FARKLI k ILE IKI VERITABANI KURUP KARSILASTIRMAK
   ------------------------------------------------
       bash install.sh kraken-kur --kmer 35 --db ~/k2db_k35
       bash install.sh kraken-kur --kmer 31 --db ~/k2db_k31
       bash tools/kraken_tool.sh esik      # esik taramasi
       bash tools/kraken_tool.sh tablo     # dort sutunlu karsilastirma
EOF
  if ! command -v kraken2-build >/dev/null 2>&1; then
    hata "kraken2-build bulunamadi. Once: bash install.sh araclar"
    EKSIK+=("kraken2-build"); return 1
  fi

  bilgi
  bilgi "Taksonomi indiriliyor (bir kez, ~1-2 saat, kesilirse devam eder)..."
  kraken2-build --download-taxonomy --db "$DB" --threads "$IS" \
    || { hata "taksonomi indirilemedi"; EKSIK+=("kraken taksonomi"); return 1; }

  local kutuphaneler="${KUTUP:-bacteria archaea fungi protozoa}"
  local L
  for L in $kutuphaneler; do
    bilgi "kutuphane indiriliyor: $L"
    kraken2-build --download-library "$L" --db "$DB" --threads "$IS" \
      || { hata "kutuphane indirilemedi: $L"; EKSIK+=("kraken kutuphane $L"); }
  done

  bilgi "veritabani kuruluyor (k=$KMER, l=$MINI, s=$BOSLUK) - SAATLER surer"
  kraken2-build --build --db "$DB" --threads "$IS" \
      --kmer-len "$KMER" --minimizer-len "$MINI" --minimizer-spaces "$BOSLUK" \
    || { hata "veritabani kurulamadi"; EKSIK+=("kraken build"); return 1; }

  if [ -f "$DB/hash.k2d" ]; then
    bilgi "TAMAM: $DB  ($(du -sh "$DB" 2>/dev/null | cut -f1))"
    # Kurulan k-mer'i KAYDET: opts.k2d ikilidir, insan okuyamaz. Hangi k ile
    # kuruldugu yazilmazsa alti ay sonra bilinmez ve karsilastirma anlamsizlasir.
    printf 'kurulum: %s\nkmer: %s\nminimizer: %s\nbosluk: %s\nkutuphane: %s\n' \
      "$(date '+%Y-%m-%d %H:%M')" "$KMER" "$MINI" "$BOSLUK" "$kutuphaneler" \
      > "$DB/KURULUM_BILGISI.txt"
    bilgi "kurulum bilgisi yazildi: $DB/KURULUM_BILGISI.txt"
    bilgi "Bracken icin: bracken-build -d $DB -t $IS -k $KMER -l <okuma_uzunlugu>"
  else
    hata "hash.k2d olusmadi - kurulum tamamlanmamis"; EKSIK+=("kraken build")
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
    hata "conda bulunamadi. QIIME2 conda ister (micromamba ile sorun cikarabilir)."
    EKSIK+=("qiime2 (conda yok)"); return 1
  fi
  local SURUM="${PT_QIIME_SURUM:-2024.10}"
  local ORT="qiime2-amplicon-$SURUM"
  if conda env list | grep -q "$ORT"; then
    bilgi "ortam zaten var: $ORT"
  else
    bilgi "QIIME2 $SURUM kuruluyor..."
    conda env create -n "$ORT" \
      --file "https://data.qiime2.org/distro/amplicon/qiime2-amplicon-$SURUM-py310-linux-conda.yml" \
      || { hata "QIIME2 kurulamadi"; EKSIK+=("qiime2"); }
  fi
  if conda env list | grep -q "picrust2"; then
    bilgi "picrust2 ortami zaten var"
  else
    conda create -y -n picrust2 -c conda-forge -c bioconda picrust2 \
      || { hata "PICRUSt2 kurulamadi"; EKSIK+=("picrust2"); }
  fi
  local SINIF="$REFDB/silva_classifier.qza"
  if [ -s "$SINIF" ]; then
    bilgi "SILVA siniflandirici zaten var"
  else
    bilgi "SILVA siniflandirici indiriliyor..."
    wget -c -q --show-progress -O "$SINIF" \
      "https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza" \
      || { hata "siniflandirici indirilemedi"; EKSIK+=("silva_classifier"); rm -f "$SINIF"; }
  fi
}

# ===========================================================================
ozet() {
  renk "OZET"
  if [ ${#SUPHELI[@]} -gt 0 ]; then
    hata "DOGRULAMAYI GECEMEYEN dosyalar (.SUPHELI olarak duruyor, KULLANILMIYOR):"
    printf '     - %s\n' "${SUPHELI[@]}"
  fi
  if [ ${#EKSIK[@]} -gt 0 ]; then
    hata "KURULAMAYANLAR:"
    printf '     - %s\n' "${EKSIK[@]}"
    bilgi
    bilgi "Bunlar SESSIZCE atlanmadi. Duzeltip ayni komutu tekrar calistirin;"
    bilgi "tamamlanmis adimlar atlanacaktir."
    bilgi "log: $LOG"
    return 1
  fi
  printf '   \033[32mKurulum tamam.\033[0m Hicbir adim atlanmadi.\n'
  bilgi "durumu gormek icin: bash install.sh durum"
  bilgi "log: $LOG"
  return 0
}

# ===========================================================================
KOMUT="${1:-durum}"; shift 2>/dev/null || true
ortam_denetimi "$KOMUT"
case "$KOMUT" in
  durum)         komut_durum ;;
  araclar)       komut_araclar; ozet ;;
  veritabani)    komut_veritabani "$@"; ozet ;;
  kraken-indir)  komut_kraken_indir ;;
  kraken-kur)    komut_kraken_kur "$@"; ozet ;;
  qiime)         komut_qiime; ozet ;;
  hepsi)         komut_araclar; komut_veritabani; komut_qiime
                 renk "KRAKEN2"
                 bilgi "Kraken2 veritabani AYRI secim ister (hazir indir mi, k-mer secip kur mu):"
                 bilgi "   bash install.sh kraken-indir     # hazir, k=35 sabit"
                 bilgi "   bash install.sh kraken-kur --kmer 31"
                 ozet ;;
  *)
    sed -n '2,45p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit 2 ;;
esac
