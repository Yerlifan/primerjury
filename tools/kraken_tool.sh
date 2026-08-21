#!/bin/bash
# ---------------------------------------------------------------------------
# KRAKEN2 YENIDEN KOSU ARACI, TUS TAKIMI
#
# NEDEN BU ARAC VAR
# Hoca "sonuclar Kraken'den cok farkli olmamali" dedi. Cevabi Kraken'in kendi
# diliyle vermek istiyoruz. Ama ayni veritabaniyla yeniden kosmak ayni cevabi
# verir. Bu yuzden iki yon acildi:
#   1. YUKSELTILMIS GUVEN ESIGI. Kraken2 bir okumayi ancak k-mer'lerinin belli
#      bir bolumu ayni klada gidiyorsa atar. Esik 0 iken tek k-mer bile yeter.
#      Esigi yukseltip neyin ayakta kaldigina bakmak, o atamalarin bastan ne
#      kadar zayif oldugunun Kraken'in kendi agziyla olcusudur.
#   2. DAHA GENIS KAPSAMLI VERITABANI (PlusPFP). Teshisimiz "sorun kapsamdi"
#      idi. Dogruysa, dar veritabaninda esikle cokan atamalar genis
#      veritabaninda ayakta kalmalidir. Kalmazsa teshis yanlistir ve oyle yazilir.
#
# KAYNAK CALISMANIN BETIKLERINDEN NE ALINDI
#   troubleshooting tools/kraken2_driver.sh  cagri bicimi, --memory-mapping fikri
#   troubleshooting tools/bracken_species.sh kmer_distrib ve -l S kullanimi
#   the source study's Kraken/Bracken script  VT butunluk denetimi (hash/opts/taxo),
#                                            bellek denetimi, log bicimi
#   WSL/rerun_kraken.sh                 VT otomatik bulma, micromamba ortami,
#                                            tek kosu birlestirme, kayipsizlik denetimi
#   tools/kraken_summary.py                  rapor ayristirma, tur duzeyi toplama
# Bu arac onlari yeniden yazmaz, cagirir ve uzerine ekler.
#
# TUSLAR
#   bash kraken_tool.sh bellek-ayari  .wslconfig onerir (KILITLENME COZUMU, ONCE BUNU)
#   bash kraken_tool.sh dogrula-ornek ornek tam veriyi temsil ediyor mu
#   bash kraken_tool.sh kraken-yol   kraken2'nin cozulmus tam yolunu basar (makine okunur)
#   bash kraken_tool.sh durum        ortam ve veritabani denetimi, hicbir sey kosmaz
#   bash kraken_tool.sh vt-ara       diskte kraken2 veritabani arar ve listeler
#   bash kraken_tool.sh vt-kimlik    veritabani hangi surum, derin tespit
#   bash kraken_tool.sh ozgun-vt       ozgun kosu hangi veritabanini kullandi, kanittan cikarim
#   bash kraken_tool.sh sinav        butun selftestler
#   bash kraken_tool.sh sure         kucuk denemeyle GERCEK hiz olcumu ve sure tahmini
#   bash kraken_tool.sh esik-a       guven esigi taramasi, VT_A uzerinde
#   bash kraken_tool.sh esik-b       guven esigi taramasi, VT_B uzerinde
#   bash kraken_tool.sh esik         ikisi de, sonra yan yana egri
#   bash kraken_tool.sh tablo        dort sutunlu karsilastirma tablosu
#   bash kraken_tool.sh hepsi        sinav + esik + tablo
#   bash kraken_tool.sh ozelvt-kur   ozel veritabani kurulumu (EN SON CARE)
#   bash kraken_tool.sh ozelvt-kos   ozel veritabaniyla kosu
#
# DEGISKENLER
#   VT_A=~/k2db       birinci veritabani (varsayilan, PlusPFP burada bekleniyor)
#   VT_B=/yol         ikinci veritabani (varsa iki egri yan yana cizilir)
#   IPLIK=12          is parcacigi sayisi
#   KAP=3000          takson basina en fazla okuma (0 = hepsi)
#   ESIKLER="0 0.02 0.05 0.1 0.2 0.5"
#   TABLO_ESIK=0.1    tabloda kullanilacak yuksek esik
# ---------------------------------------------------------------------------
set -euo pipefail

TUS="${1:-yardim}"

# --- proje klasoru olculur, tahmin edilmez (rerun_kraken.sh ile ayni yol) ---
_BETIK_DIZIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJE="${PROJE:-$(cd "$_BETIK_DIZIN/.." && pwd)}"
# OLCULDU: burasi bir zamanlar kaynak projede duran TESLIM klasorunu ariyordu.
# O klasor bu depoda YOK, dolayisiyla betik hangi tusla cagrilirsa cagrilsin ilk
# denetimde duruyordu. Isaret artik depoda GERCEKTEN duran iki klasor.
if [ ! -d "$PROJE/tools" ] || [ ! -d "$PROJE/verification" ]; then
  echo "ERROR: the project root could not be verified ('$PROJE' holds no tools/ and verification/)."
  echo "  Script location: $_BETIK_DIZIN"
  echo "  To give it by hand:  PROJE=/full/path/to/project bash $0 $TUS"
  exit 1
fi

VT_A="${VT_A:-$HOME/k2db}"
VT_B="${VT_B:-}"
# IPLIK: eskiden nproc idi, yani butun cekirdekler. Kraken2'nin TEPE bellek
# tuketimi is parcacigi sayisiyla artar (her parcacik kendi tamponunu tutar) ve
# WSL2'de bu tepe dogrudan Windows'u bogar. Varsayilan 3'e cekildi; olcum sonucu
# DEGISMEZ, yalnizca sure biraz uzar. Guclu makinede IPLIK=8 ile artirilabilir.
IPLIK="${IPLIK:-3}"
# ORNEK: esik taramasinda kullanilacak TOPLAM okuma sayisi. 0 = hepsi.
# Gerekce asagida esik_tara icinde yazili.
ORNEK="${ORNEK:-100000}"
KAP="${KAP:-3000}"
ESIKLER="${ESIKLER:-0 0.02 0.05 0.1 0.2 0.5}"
TABLO_ESIK="${TABLO_ESIK:-0.1}"
OZELVT="${OZELVT:-$HOME/k2_ozel}"
KAYNAK="$PROJE/SONUCLAR/fastq files"
IS_A="$PROJE/SONUCLAR/kraken_esik_A"
IS_B="$PROJE/SONUCLAR/kraken_esik_B"
OZEL_IS="$PROJE/SONUCLAR/kraken_ozelvt"

log_ac() {
  local ad="$1"
  mkdir -p "$PROJE/kurulum_loglari"
  LOG="$PROJE/kurulum_loglari/kraken_${ad}_$(date '+%Y%m%d_%H%M%S').log"
  exec > >(tee -a "$LOG") 2>&1
  echo "=============================================================="
  echo "kraken_tool  key: $ad  start $(date '+%Y-%m-%d %H:%M:%S %Z')"
  echo "  machine: $(hostname)  $(uname -sr)  $IPLIK threads"
  echo "  log   : $LOG"
  echo "=============================================================="
  trap 'echo; echo "finished $(date "+%Y-%m-%d %H:%M:%S %Z")"; echo "log: $LOG"' EXIT
}

# --- ortam. rerun_kraken.sh ile ayni: kraken2 micromamba "mikro" icinde --
# ---------------------------------------------------------------------------
# ortam_ac - kraken2'yi BULUNABILIR hale getirir.
#
# NEDEN GENISLETILDI (2026-08-04)
# Bu fonksiyon kaynak calismanin Kraken betiklerinden
# devralinmisti: kraken2 PATH'te degildir, micromamba'nin "mikro" ortamindadir.
# Ancak yalnizca "ortami etkinlestir" yolu deneniyordu. Kullanicinin kosusunda
# micromamba kabuk kancasi calismadi ve kraken2 bulunamadi diye butun Kraken
# adimlari atlandi. Oysa ikili diskte duruyordu.
#
# Artik dort yol SIRAYLA denenir ve ilki tutunca durulur:
#   0) KRAKEN2_BIN disaridan verilmisse dogrudan o kullanilir (en ust yetki).
#   1) Zaten PATH'te mi.
#   2) KAYNAK CALISMANIN YOLU - degistirilmedi, aynen korundu: $HOME/bin PATH'e eklenir,
#      MAMBA_ROOT_PREFIX kurulur, micromamba ya da conda ile "$ORTAM" acilir.
#   3) Ortam KLASORLERINE dogrudan bakilir. Kabuk kancasi calismasa bile ikili
#      yerindedir; etkinlestirmeye gerek yok, bin klasorunu PATH'e eklemek yeter.
#   4) Son care: $HOME altinda sinirli derinlikte kraken2 adli calistirilabilir
#      dosya aranir.
#
# Bulunan yol KRAKEN2_BIN degiskenine yazilir ve disari verilir; boylece sonraki
# adimlar PATH'e bagimli kalmaz.
# ---------------------------------------------------------------------------
ortam_ac() {
  ORTAM="${ORTAM:-mikro}"
  export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
  KRAKEN_YONTEM=""
  BAKILAN_YERLER=""

  _kayit() { BAKILAN_YERLER="${BAKILAN_YERLER}
    $1"; }

  # kraken2 gercekten CALISIYOR mu. Dosyanin var olmasi yetmez: conda ikilileri
  # ortamin lib klasorune bagimli olabilir ve ortam disindan cagrilinca
  # kutuphane bulamayip patlayabilir. O yuzden --version fiilen denenir.
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

  # --- 2) ORTAM KLASORLERI, dogrudan -------------------------------------
  # Kabuk kancasini calistirmadan once buraya bakilir. Sebep: kullanicinin
  # kurulumunda tek bir ortam var (mikro) ve ikili orada duruyor; kancayi
  # calistirmak hem yavas hem de bazi kabuklarda sessizce basarisiz oluyor.
  # Dosyaya dogrudan bakmak daha hizli ve daha kesin.
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
      _benimse "$k" "ortam klasorunde bulundu, dogrudan calistirilabiliyor"
      return 0
    fi
    # Ikili var ama dogrudan calismiyor: ortam kutuphanelerine ihtiyaci var.
    # Bu durumda micromamba run ile cagirmak ORTAMI BOZMADAN sorunu cozer
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

  # --- 3) KAYNAK CALISMANIN YOLU: ortami etkinlestir -------------------------------
  # Kaynak calismanin Kraken betiklerinden aynen devralindi.
  _kayit "micromamba/conda kabuk kancasi + '$ORTAM' ortaminin etkinlestirilmesi"
  export PATH="$HOME/bin:$PATH"
  if command -v micromamba >/dev/null 2>&1; then
    eval "$(micromamba shell hook -s bash)" || true
    micromamba activate "$ORTAM" || true
  elif command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)" || true
    conda activate "$ORTAM" || true
  fi
  if command -v kraken2 >/dev/null 2>&1; then
    _benimse "$(command -v kraken2)" "ortam etkinlestirilerek bulundu ($ORTAM)"
    return 0
  fi

  # --- 4) micromamba run, ortam etkinlestirmeden --------------------------
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
  _kayit "\$HOME altinda 6 derinlige kadar 'kraken2' adli calistirilabilir dosya"
  k="$(find "$HOME" -maxdepth 6 -type f -name kraken2 -perm -u+x 2>/dev/null | head -1)"
  if [ -n "$k" ] && _calisiyor_mu "$k"; then
    _benimse "$k" "ev dizininde arama ile bulundu"; return 0
  fi
  export BAKILAN_YERLER
  return 1
}

# ---------------------------------------------------------------------------
# tus_kraken_yol - kraken2'nin COZULMUS tam yolunu makine okunur basar.
#
# Neden ayri bir tus: full_chain.py (A tusu) eskiden kraken2'yi kendi basina
# PATH'te ariyordu ve bulamayinca butun Kraken adimlarini atliyordu. Artik
# aramayi yapmiyor, bu tusu cagirip cevabi buradan aliyor. Arama mantigi TEK
# YERDE durur; iki ayri yerde iki farkli arama olmaz.
#
# Cikti: bulursa tek satir  KRAKEN2_BIN=<tam yol>  ve cikis 0.
#        bulamazsa cikis 1 ve stderr'e kurulum yonergesi.
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
    echo "  (for a WSL user named alice that is: /home/alice/micromamba/envs/mikro/bin/kraken2)"
    echo
    echo "INSTALL - the route this project uses:"
    echo "    bash install.sh tools"
    echo "To install it by hand:"
    echo "    micromamba create -n mikro -c bioconda -c conda-forge kraken2 bracken"
    echo
    echo "IF IT IS ALREADY INSTALLED the environment may have a different name:"
    echo "    micromamba env list        (or: conda env list)"
    echo "    ORTAM=<environment_name> bash $0 kraken-yol"
    echo "If you know the full path of the binary:"
    echo "    KRAKEN2_BIN=/full/path/kraken2 bash $0 kraken-yol"
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
  [ -n "${KRAKEN_YONTEM:-}" ] && echo "  bulunma yontemi: $KRAKEN_YONTEM"
  for _a in kraken2-inspect bracken; do
    if command -v "$_a" >/dev/null 2>&1; then echo "  $_a: $(command -v $_a)"
    else echo "  $_a: YOK"; fi
  done
}

# =========================================================================
# VERITABANI HAZIRLIK DENETIMI
# Varsaymaz. Uc dosya da aranir, boyut ve tarih yazilir. Yalnizca .tar.gz
# duruyorsa acma komutu ekrana yazilir ama KENDILIGINDEN ACILMAZ; gigabaytlarca
# dosyayi kullanicinin haberi olmadan acmak kabul edilebilir degil.
# Doner: 0 hazir, 1 sadece arsiv var, 2 hicbir sey yok
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
    echo "    After unpacking:     bash $0 vt-kimlik"
    return 1
  fi
  echo "    NO INDEX and no archive either: $d"
  return 2
}

# =========================================================================
# VERITABANI SURUM TESPITI
# Bu tespit onemli: butun argumanimiz veritabaninin KAPSAMI uzerine kurulu.
# Yanlis surumle kosarsak sonucu yanlis yorumlariz.
#
# IKI BAGIMSIZ OLCUM (proje kurali 1)
#   olcum 1: kraken2-inspect ile veritabaninin ICINDEKI taksonlar. Kesindir.
#   olcum 2: hash.k2d boyutu. Kaba ama bagimsizdir.
# Ikisi ayrilirsa AYRILIK yazilir ve tek bir olcume gecilmez.
#
# Ayirt edici: bitki (Viridiplantae, taxid 33090).
#   Standard = arke, bakteri, virus, plazmit, insan, UniVec
#   PlusPF   = Standard + protozoa + mantar
#   PlusPFP  = PlusPF + BITKI
# =========================================================================
vt_surum() {
  local d="$1"
  echo "  version detection: $d"
  local hb hgb
  hb=$(stat -c%s "$d/hash.k2d" 2>/dev/null || echo 0)
  hgb=$((hb / 1073741824))
  echo "    hash.k2d boyutu: ${hgb} GB"

  # --- olcum 2: boyut ---
  local o2
  if   [ "$hgb" -ge 125 ]; then o2="PlusPFP (full)"
  elif [ "$hgb" -ge 95  ]; then o2="PlusPF (full)"
  elif [ "$hgb" -ge 60  ]; then o2="Standard (full)"
  elif [ "$hgb" -ge 12  ]; then o2="capped version (16 GB class)"
  elif [ "$hgb" -ge 5   ]; then o2="capped version (8 GB class)"
  else                          o2="a small or custom database"
  fi
  echo "    measurement 2 (size)    : $o2"

  # --- olcum 1: inspect ---
  local o1="not measured"
  local ins="$PROJE/SONUCLAR/vt_inspect_$(basename "$d").txt"
  if command -v kraken2-inspect >/dev/null 2>&1; then
    # ONBELLEK, AMA KAYNAK DAMGASIYLA (2026-08-04 duzeltmesi)
    # Eskiden yalnizca "dosya var mi" diye bakiliyordu ve varsa yeniden
    # kosulmuyordu. Bu, projenin klasik hata turunu davet eden bir tasarimdi:
    # baska bir veritabanina ait ya da bayat bir inspect ciktisi sessizce
    # yeniden kullanilir, arac hata vermeden YANLIS surum bildirirdi.
    # Artik ciktinin yaninda kaynagin boyut ve tarih damgasi tutulur; damga
    # tutmuyorsa onbellek gecersiz sayilir ve inspect YENIDEN kosar.
    local damga_dosya="$ins.kaynak"
    local damga_simdi
    damga_simdi="$(stat -c '%s %Y' "$d/hash.k2d" 2>/dev/null || echo 'yok')"
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

  # --- parmak izi. Eski kosunun hangi VT ile yapildigi buradan izlenir. -----
  local pi="$PROJE/SONUCLAR/vt_parmak_izi.tsv"
  mkdir -p "$(dirname "$pi")"
  [ -f "$pi" ] || printf 'tarih\tyol\thash_bayt\thash_tarih\ticerik\tboyut\n' > "$pi"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$(date '+%Y-%m-%d %H:%M')" "$d" "$hb" \
    "$(date -r "$d/hash.k2d" '+%Y-%m-%d' 2>/dev/null || echo '?')" "$o1" "$o2" >> "$pi"
}

# =========================================================================
# ESKI KOSU HANGI VERITABANIYLA YAPILDI
# Bu soru kritik: eski kosu da ~/k2db ile yapildiysa "eski VT" ile "PlusPFP"
# ayni sey olabilir ve karsilastirmanin anlami degisir.
# =========================================================================
eski_kosu_tespit() {
  echo "WHICH DATABASE THE OLD RUNS USED"
  echo
  echo "  a) A direct record from the logs:"
  local bulundu=0
  for l in "$PROJE"/kurulum_loglari/*.log; do
    [ -f "$l" ] || continue
    local yol boy
    # Bu grep KENDI loglarimizi okur. Ekran metni Ingilizceye cevrildi, ama
    # eski loglar hala Turkce; ikisi de tutsun diye desen iki kelimeyi de alir.
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
  for r in "$PROJE/SONUCLAR/kraken_yeniden/tum.report" \
           "$PROJE/SONUCLAR/kraken results"/*/*_kraken2.report; do
    [ -f "$r" ] || continue
    local m p b
    m=$(awk -F'\t' '$5=="4751"{print "mantar VAR"}' "$r" | head -1)
    p=$(awk -F'\t' '$5=="5794"||$5=="33682"{print "protozoa VAR"}' "$r" | head -1)
    b=$(awk -F'\t' '$5=="33090"{print "BITKI VAR"}' "$r" | head -1)
    printf "     %-52s %s %s %s\n" "$(basename "$(dirname "$r")")/$(basename "$r")" \
           "${m:-mantar yok}" "${p:-protozoa yok}" "${b:-bitki yok}"
    break
  done
  local ilk="$PROJE/SONUCLAR/kraken_yeniden/tum.report"
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
  echo "  kraken2-inspect:  bash $0 vt-kimlik"
}

# =========================================================================
# TUS: vt-ara   VERITABANI ARAMA
# Yol yanlissa pes etmez, arar. Butun disk taranmaz; belirli kok dizinler ve
# sinirli derinlik kullanilir, sure olculur ve yazilir. Birden fazla bulursa
# KENDI KAFASINA GORE SECMEZ, kullaniciya sectirir.
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
  # /mnt/c altinda Windows tarafi cok yavas olabilir, ayri ve daha sig aranir.
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
    echo "      VT_A=${benzersiz[0]} bash $0 esik"
  else
    echo "  THERE IS MORE THAN ONE DATABASE. The script does NOT CHOOSE which one to run,"
    echo "  because running the wrong version means reading the result wrongly."
    echo "  Look at the version lines above, choose, and give the path:"
    for y in "${benzersiz[@]}"; do echo "      VT_A=$y bash $0 esik"; done
    echo "  If you want to compare two databases:"
    echo "      VT_A=<wide> VT_B=<narrow> bash $0 esik"
  fi
  return 0
}

# =========================================================================
# TUS: ozgun-vt   OZGUN KOSU HANGI VERITABANINI KULLANDI
# Kullanici da bilmiyor, sorulacak kimse yok. Kanittan cikarilir.
# Cikarilamazsa "belirlenemedi" yazilir, uydurulmaz.
# =========================================================================
tus_ali_vt() {
  # Bu tus yalnizca OKUR ve rapor yazar, hicbir sey degistirmez. Aranan dosyanin
  # bulunmamasi normaldir; eksik kanit isi durdurmamali, "bulunamadi" diye
  # yazilmali. Bu yuzden errexit burada kapatilir.
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
  for l in "$PROJE"/kurulum_loglari/*.log; do
    [ -f "$l" ] || continue
    local yol boy
    # Bu grep KENDI loglarimizi okur. Ekran metni Ingilizceye cevrildi, ama
    # eski loglar hala Turkce; ikisi de tutsun diye desen iki kelimeyi de alir.
    yol=$(grep -m1 -oE "(veritabani|database) *: */[^ ]*" "$l" 2>/dev/null | sed 's/.*: *//' || true)
    boy=$(grep -m1 -oE "hash\.k2d *(var|present) *[0-9.]+[KMGT]" "$l" 2>/dev/null | grep -oE "[0-9.]+[KMGT]$" || true)
    [ -n "$yol" ] && { printf "    %-42s %-26s hash %s\n" "$(basename "$l")" "$yol" "${boy:-?}"; k2=1; }
  done
  [ "$k2" -eq 1 ] || echo "    there is no log from the original run (these logs are from our runs)"

  echo
  echo "EVIDENCE 3: report content (which domains WERE in the database)"
  local r="$PROJE/SONUCLAR/kraken results/A1/edited_barcode01_kraken2.report"
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
  local oz="$PROJE/SONUCLAR/kraken_yeniden/kraken_ozet.csv"
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

# --- bellek. Kaynak Kraken/Bracken betigi ve rerun_kraken.sh'teki karar --------
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
  local eski="$PROJE/SONUCLAR/kraken_yeniden/tum.fastq"
  if [ -s "$hedef" ]; then
    echo "the merged file already exists: $(awk 'END{print int(NR/4)}' "$hedef") reads"
    return
  fi
  mkdir -p "$(dirname "$hedef")"
  if [ -s "$eski" ]; then
    echo "the tum.fastq produced by rerun_kraken.sh was found, so THE SAME read set is used."
    echo "  (That is the condition that makes thresholds and databases comparable.)"
    cp "$eski" "$hedef"
    cp "$PROJE/SONUCLAR/kraken_yeniden/kaynak_sayim.tsv" "$(dirname "$hedef")/" 2>/dev/null || true
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
# ornekle - esik taramasi icin TEMSILI alt kume kurar.
#
# NEDEN
# Esik taramasinin amaci EGRININ SEKLINI gormektir: esik yukseldikce hangi
# atamalar cokuyor. Bunun icin butun okumalari her esikte yeniden
# siniflandirmak gereksizdir. Alti esik x 330 bin okuma, 110 GB'lik bir
# veritabaniyla, WSL2'de sayfa onbellegi uzerinden Windows'u kilitliyordu.
#
# Alt kume KUTU BASINA ESIT alinir, bollugu buyuk taksonlar orantisiz
# agirlik kazanmasin diye. Yuzdeler oranlardan cikar; oranlar temsili bir
# ornekte de ayni kalir, mutlak sayilar kucululur.
#
# TEMSIL DOGRULAMASI ayri bir istir ve bu fonksiyonun isi degildir:
# 'dogrula-ornek' tusu tek bir esikte tam veriyle ornegi karsilastirir.
# Oranlar tutuyorsa egri gecerlidir. Tutmuyorsa ORNEK buyutulur.
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
  # Kutu paylastirmasi, birlestir()'in okuma basliklarina koydugu 'tx<taxid>_'
  # onekine dayanir. birlestir() eski bir tum.fastq'yu yeniden kullanmis ve o
  # dosyada onek yoksa, BUTUN okumalar tek kutu sayilir ve ornek temsili olmaz.
  # Bu sessiz bir basarisizlik olurdu; onun icin once olculur.
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
  # NEDEN ILK N DEGIL (2026-08-05 duzeltmesi)
  # Onceki surum her kutudan ILK 'pay' okumayi aliyordu. Nanopore kosusunda
  # gozenek verimi zamanla duser ve okuma kalitesi kosunun sonuna dogru
  # sistematik olarak degisir; dosyanin basi sonundan farklidir. Ilk N almak,
  # kosunun erken donemini asiri temsil eden YANLI bir ornek uretir ve bu
  # yanlilik esik egrisini kaydirabilir.
  # Artik her okuma, kutu icindeki sirasindan BAGIMSIZ bir karma degerle
  # puanlanir ve en kucuk 'pay' tanesi secilir. Karma, okuma adindan ve sabit
  # bir tohumdan uretilir; yani secim RASTGELE ama TEKRARLANABILIRDIR - ayni
  # girdi ve ayni tohum her zaman ayni ornegi verir. ORNEK_TOHUM ile
  # degistirilebilir.
  local TOHUM="${ORNEK_TOHUM:-20260805}"
  echo "  sampling: $kutu_sayisi bins x $pay reads = target ~$toplam"
  echo "  method  : RANDOM within a bin, fixed seed $TOHUM (reproducible)"

  # Iki gecis: once her okumanin karmasi, sonra kutu basina esik.
  # awk'in kendi rand()'i yerine ad tabanli karma kullanilir; boylece okuma
  # sirasi degisse bile ayni okumalar secilir.
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
# TUS: sure   GERCEK OLCUMLE SURE TAHMINI
# Tahmin uydurulmaz. Kucuk bir ornek gercekten kosulur, sure olculur ve
# tam kosuya oranlanir. --memory-mapping'te ilk okumalar en yavas olandir
# (sayfalar diskten gelir), bu yuzden olcum bir UST sinir gibi davranir ve
# bunun boyle oldugu ekrana yazilir.
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
    if command -v "$a" >/dev/null 2>&1; then printf "  %-16s var   %s\n" "$a" "$(command -v $a)"
    else printf "  %-16s YOK\n" "$a"; fi
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
    echo "  To compare two databases:  VT_B=/path/old_db bash $0 esik"
  fi
  echo
  local n; n=$(ls "$KAYNAK"/*/*reads_*.fastq 2>/dev/null | wc -l || echo 0)
  echo "reads: $n fastq files ($KAYNAK)"
  echo
  free -g 2>/dev/null | awk '/Mem:/{print "RAM: "$2" GB toplam, "$7" GB kullanilabilir"}'
  df -h "$VT_A" 2>/dev/null | tail -1 | awk '{print "disk: "$4" bos ("$6")"}'
  echo
  eski_kosu_tespit
  echo
  echo "Next step:  bash $0 vt-kimlik    (which version the database is)"
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
# TUS: esik-a / esik-b   guven esigi taramasi
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
  # ESIK TARAMASI ORNEK UZERINDE KOSAR. Gerekce ornekle() icinde yazili.
  # Tam veriyle kosmak icin: ORNEK=0
  local GIRDI="$is/tum.fastq"
  if [ "${ORNEK:-0}" -gt 0 ]; then
    echo
    echo "SAMPLING IS ON (ORNEK=$ORNEK). For the full data: ORNEK=0 bash $0 esik"
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
    echo "two different databases, give the path of the second: VT_B=/path/other bash $0 esik"
  fi
  esik_tara "$VT_A" "$IS_A" "VT_A"
  if [ "$ikili" -eq 1 ]; then
    esik_tara "$VT_B" "$IS_B" "VT_B"
    echo
    python3 "$_BETIK_DIZIN/threshold_summary.py" --kok "$PROJE" \
      --is "$IS_A" --ad "$(basename "$VT_A")" \
      --is2 "$IS_B" --ad2 "$(basename "$VT_B")"
  else
    echo
    python3 "$_BETIK_DIZIN/threshold_summary.py" --kok "$PROJE" \
      --is "$IS_A" --ad "$(basename "$VT_A")"
  fi
  echo
  echo "done. Files: $IS_A"
  echo "  esik_<C>.report          the kraken report for each threshold"
  echo "  esik_egrisi.csv / .txt   assignment percentages per domain, against the threshold"
  # NOT: son komut olarak "[ ... ] && echo" birakilmaz. Kosul yanlis oldugunda
  # test 1 doner ve fonksiyonun cikis kodu olur; is basariyla bittigi halde
  # betik HATA vermis gibi gorunur. Bu tuzaga bir kez dusuldu.
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
# sutunu, esik taramasi yoksa 2. sutunu, hizalama sonucu yoksa 4. sutunu bos
# birakir ve hangilerinin eksik oldugunu basar. Bu yuzden Kraken olcumu
# haftalar sonra yapilsa bile bu tus yeniden cagrilarak tablo TAMAMLANIR;
# bastan kosuya gerek yoktur.
#
# Ama bunun bir riski var: elde TAM bir tablo varken, veritabani gecici olarak
# erisilemez oldugu bir anda bu tusa basmak tabloyu DAHA BOS bir surumle
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
  python3 "$_BETIK_DIZIN/comparison_table.py" --kok "$PROJE" \
          --is-a "$IS_A" --ad-a "$(basename "$VT_A")" \
          --is-b "$IS_B" --ad-b "$([ -n "$VT_B" ] && basename "$VT_B" || echo '')" \
          --esik "$TABLO_ESIK"
}

# =========================================================================
# TUS: ozelvt-kur   EN SON CARE
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
  echo "Check:  bash $0 vt-kimlik    (if PlusPFP is installed, do not come here at all)"
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
    [ "$c" = "evet" ] || { echo "iptal edildi. Hicbir sey yapilmadi."; exit 0; }
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
  python3 "$_BETIK_DIZIN/custom_taxonomy.py" --cikti "$OZELVT" "${args[@]}"
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
    echo "  First:  bash $0 ozelvt-kur"; exit 1; }
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
# tus_bellek_ayari - WSL2'nin Windows'u bogmasini engelleyen .wslconfig uretir.
#
# SORUN
# WSL2 ayri bir sanal makinedir. Varsayilan olarak Windows RAM'inin buyuk bir
# kismini alabilir VE 110 GB'lik bir veritabani --memory-mapping ile okundugunda
# sayfa onbellegi de sanal makinenin bellegine sayilir. Sonuc: WSL sinirsizca
# buyur, Windows takas etmeye baslar ve makine kilitlenir.
#
# COZUM
# .wslconfig ile SANAL MAKINEYE UST SINIR konur. Windows'a nefes payi kalir.
# WSL yavaslar, ama makine kullanilabilir kalir. Bu bir takas degil, dogru
# yapilandirmadir: sinirsiz WSL zaten kimseye hizli calismiyor, kilitliyor.
# ---------------------------------------------------------------------------
tus_bellek_ayari() {
  echo "======================================================================"
  echo "WSL MEMORY SETTING - it stops Windows from freezing"
  echo "======================================================================"
  echo
  local TOPLAM_MB=0 KAYNAK="olculemedi"
  if command -v powershell.exe >/dev/null 2>&1; then
    TOPLAM_MB=$(powershell.exe -NoProfile -Command \
      "[math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory/1MB)" \
      2>/dev/null | tr -d '\r' | tr -d ' ')
    [ -n "$TOPLAM_MB" ] && KAYNAK="Windows'tan olculdu"
  fi
  case "$TOPLAM_MB" in ''|*[!0-9]*) TOPLAM_MB=0 ;; esac
  if [ "$TOPLAM_MB" -le 0 ]; then
    TOPLAM_MB=$(( $(awk '/MemTotal/{print $2}' /proc/meminfo) / 1024 ))
    KAYNAK="estimated from inside WSL (the Windows total can be LARGER)"
  fi
  local TOPLAM_GB=$(( TOPLAM_MB / 1024 ))
  # Yuzde 60: WSL'e yeter, Windows'a nefes payi birakir. Alt sinir 4 GB,
  # cunku 110 GB'lik veritabani mmap ile bile bir taban bellek ister.
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
# tus_dogrula_ornek - ORNEGIN TEMSIL ETTIGINI olcer.
#
# Tek bir esikte hem ornegi hem TAM veriyi siniflandirir ve YUZDELERI
# karsilastirir. Mutlak sayilar elbette farkli olacak; bakilan sey oranlarin
# tutup tutmadigidir. Tutuyorsa esik egrisi gecerlidir.
#
# DIKKAT: bu tus TAM veriyle bir kez kosar, yani agirdir. Gece birakin.
# Kosmadan da esik egrisi kullanilabilir, ama o zaman raporda "temsil
# dogrulamasi YAPILMADI" yazar - ve oyle sunulmalidir.
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
      print "  ORNEK degerini buyutup tekrarlayin, orn: ORNEK=300000 bash kraken_tool.sh esik"
    }
  }'
  printf 'esik\t%s\nen_buyuk_sapma_puan\t%s\n' "$C" "$sapma_max" \
    > "$is/temsil_dogrulamasi.tsv"
  echo
  echo "  written: $is/temsil_dogrulamasi.tsv"
}

# =========================================================================
case "$TUS" in
  bellek-ayari) tus_bellek_ayari ;;
  dogrula-ornek) tus_dogrula_ornek ;;
  kraken-yol)  tus_kraken_yol ;;
  durum)       tus_durum ;;
  vt-ara)      ortam_ac; vt_ara || true ;;
  ozgun-vt)      tus_ali_vt ;;
  vt-kimlik)   tus_vt_kimlik ;;
  sinav)       tus_sinav ;;
  sure)        tus_sure "${2:-$VT_A}" ;;
  esik-a)      tus_sinav >/dev/null || { echo "A SELF TEST FAILED. Detail: bash $0 sinav"; exit 2; }
               log_ac esik_a; kraken_sart; esik_tara "$VT_A" "$IS_A" "VT_A"
               python3 "$_BETIK_DIZIN/threshold_summary.py" --kok "$PROJE" --is "$IS_A" --ad "$(basename "$VT_A")" ;;
  esik-b)      [ -n "$VT_B" ] || { echo "ERROR: VT_B was not given.  VT_B=/path bash $0 esik-b"; exit 1; }
               tus_sinav >/dev/null || { echo "SINAV KALDI"; exit 2; }
               log_ac esik_b; kraken_sart; esik_tara "$VT_B" "$IS_B" "VT_B"
               python3 "$_BETIK_DIZIN/threshold_summary.py" --kok "$PROJE" --is "$IS_B" --ad "$(basename "$VT_B")" ;;
  esik)        tus_sinav >/dev/null || { echo "A SELF TEST FAILED. Detail: bash $0 sinav"; exit 2; }
               tus_esik ;;
  tablo)       tus_tablo ;;
  hepsi)       tus_sinav >/dev/null || { echo "SINAV KALDI"; exit 2; }
               tus_esik; tus_tablo ;;
  ozelvt-kur)  tus_ozelvt_kur ;;
  ozelvt-kos)  tus_ozelvt_kos ;;
  *)           sed -n '3,48p' "$0" | sed 's/^# \{0,1\}//'
               echo; echo "For detail: $PROJE/docs/GUIDE.md" ;;
esac
