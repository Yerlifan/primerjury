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
#   ozgun Kraken/Bracken betigi                 VT butunluk denetimi (hash/opts/taxo),
#                                            bellek denetimi, log bicimi
#   WSL/rerun_kraken.sh                 VT otomatik bulma, micromamba ortami,
#                                            tek kosu birlestirme, kayipsizlik denetimi
#   WSL/ozgun Kraken ozet betigi                    rapor ayristirma, tur duzeyi toplama
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
if [ ! -d "$PROJE/TESLIM" ]; then
  echo "HATA: proje klasoru dogrulanamadi ('$PROJE' icinde TESLIM yok)."
  echo "  Betik konumu: $_BETIK_DIZIN"
  echo "  Elle vermek icin:  PROJE=/tam/yol/proje bash $0 $TUS"
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
  mkdir -p "$PROJE/WSL/loglar"
  LOG="$PROJE/WSL/loglar/130_${ad}_$(date '+%Y%m%d_%H%M%S').log"
  exec > >(tee -a "$LOG") 2>&1
  echo "=============================================================="
  echo "kraken_tool  tus: $ad  baslangic $(date '+%Y-%m-%d %H:%M:%S %Z')"
  echo "  makine: $(hostname)  $(uname -sr)  $IPLIK iplik"
  echo "  log   : $LOG"
  echo "=============================================================="
  trap 'echo; echo "bitis $(date "+%Y-%m-%d %H:%M:%S %Z")"; echo "log: $LOG"' EXIT
}

# --- ortam. rerun_kraken.sh ile ayni: kraken2 micromamba "mikro" icinde --
# ---------------------------------------------------------------------------
# ortam_ac - kraken2'yi BULUNABILIR hale getirir.
#
# NEDEN GENISLETILDI (2026-08-04)
# Bu fonksiyon kaynak calismanin ozgun Kraken betiklerinden
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
      KRAKEN_YONTEM="micromamba run -n $envad (ikili dogrudan calismadi)"
      return 0
    fi
  done

  # --- 3) KAYNAK CALISMANIN YOLU: ortami etkinlestir -------------------------------
  # ozgun Kraken betiklerinden aynen devralindi.
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
    echo "kraken2 bulunamadi. TAM OLARAK su yerlere bakildi:"
    printf '%s\n' "$BAKILAN_YERLER" | sed '/^$/d'
    echo
    echo "Bu projedeki beklenen yer:  \$HOME/micromamba/envs/${ORTAM:-mikro}/bin/kraken2"
    echo "  (WSL kullanicisi yerlifan ise: /home/yerlifan/micromamba/envs/mikro/bin/kraken2)"
    echo
    echo "KURULUM - bu projede kullanilan yol:"
    echo "    bash $install.sh araclar"
    echo "Elle kurmak icin:"
    echo "    micromamba create -n mikro -c bioconda -c conda-forge kraken2 bracken"
    echo
    echo "ZATEN KURULUYSA ortam adi farkli olabilir:"
    echo "    micromamba env list        (ya da: conda env list)"
    echo "    ORTAM=<ortam_adi> bash $0 kraken-yol"
    echo "Ikilinin tam yolunu biliyorsaniz:"
    echo "    KRAKEN2_BIN=/tam/yol/kraken2 bash $0 kraken-yol"
  } >&2
  return 1
}

kraken_sart() {
  ortam_ac
  if ! command -v kraken2 >/dev/null 2>&1; then
    echo
    echo "HATA: kraken2 bulunamadi. Sessizce atlanmiyor, is burada duruyor."
    echo
    echo "  Gereken: kraken2 surum 2.1 ve ustu."
    echo "  Bu projede kullanilan kurulum yolu:"
    echo "      bash $install.sh araclar"
    echo "  Elle:"
    echo "      micromamba create -n mikro -c bioconda -c conda-forge kraken2 bracken"
    echo "  Zaten kuruluysa ortami etkinlestirin:"
    echo "      export PATH=\"\$HOME/bin:\$PATH\""
    echo "      export MAMBA_ROOT_PREFIX=\"\$HOME/micromamba\""
    echo "      eval \"\$(micromamba shell hook -s bash)\""
    echo "      micromamba activate ${ORTAM:-mikro}"
    echo "  Ortam adini bilmiyorsaniz:  micromamba env list"
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
    echo "  klasor yok: $d"
    return 2
  fi
  local eksik=0
  for g in hash.k2d opts.k2d taxo.k2d; do
    if [ -s "$d/$g" ]; then
      printf "    %-10s var   %8s   %s\n" "$g" \
        "$(du -h "$d/$g" | cut -f1)" \
        "$(date -r "$d/$g" '+%Y-%m-%d %H:%M' 2>/dev/null || echo '?')"
    else
      printf "    %-10s YOK\n" "$g"
      eksik=1
    fi
  done
  if [ "$eksik" -eq 0 ]; then
    echo "    -> indeks hazir"
    return 0
  fi
  local ars
  ars=$(ls "$d"/*.tar.gz 2>/dev/null | head -3 || true)
  if [ -n "$ars" ]; then
    echo
    echo "    INDEKS HAZIR DEGIL. Klasorde arsiv var ama acilmamis:"
    while IFS= read -r t; do
      [ -n "$t" ] && printf "      %-52s %8s  %s\n" "$(basename "$t")" \
        "$(du -h "$t" | cut -f1)" "$(date -r "$t" '+%Y-%m-%d' 2>/dev/null || echo '?')"
    done <<< "$ars"
    echo
    echo "    Acmak icin (kendiliginden ACILMADI, onlarca GB yer kaplar):"
    echo "        cd $d"
    while IFS= read -r t; do
      [ -n "$t" ] && echo "        tar -xzvf $(basename "$t")"
    done <<< "$ars"
    echo "    Once bos yer olcun:  df -h $d"
    echo "    Actiktan sonra:      bash $0 vt-kimlik"
    return 1
  fi
  echo "    INDEKS YOK, arsiv da yok: $d"
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
  echo "  surum tespiti: $d"
  local hb hgb
  hb=$(stat -c%s "$d/hash.k2d" 2>/dev/null || echo 0)
  hgb=$((hb / 1073741824))
  echo "    hash.k2d boyutu: ${hgb} GB"

  # --- olcum 2: boyut ---
  local o2
  if   [ "$hgb" -ge 125 ]; then o2="PlusPFP (tam)"
  elif [ "$hgb" -ge 95  ]; then o2="PlusPF (tam)"
  elif [ "$hgb" -ge 60  ]; then o2="Standard (tam)"
  elif [ "$hgb" -ge 12  ]; then o2="kapakli surum (16 GB sinifi)"
  elif [ "$hgb" -ge 5   ]; then o2="kapakli surum (8 GB sinifi)"
  else                          o2="kucuk ya da ozel veritabani"
  fi
  echo "    olcum 2 (boyut)  : $o2"

  # --- olcum 1: inspect ---
  local o1="olculemedi"
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
        echo "    onbellekteki inspect ciktisi BASKA bir veritabanina ait"
        echo "      (damga degisti: '$damga_eski' -> '$damga_simdi'), yeniden kosuluyor"
      fi
      echo "    kraken2-inspect calisiyor (1 ile 5 dakika surebilir, bir kez)"
      mkdir -p "$(dirname "$ins")"
      if kraken2-inspect --db "$d" --threads "$IPLIK" > "$ins" 2>/dev/null; then
        printf '%s' "$damga_simdi" > "$damga_dosya"
      else
        echo "    kraken2-inspect KOSMADI (cikis kodu sifir degil)."
        : > "$ins"
        rm -f "$damga_dosya"
      fi
    else
      echo "    kraken2-inspect ciktisi onbellekte ve damga tutuyor: $ins"
    fi
    if [ -s "$ins" ]; then
      local bitki mantar protozoa insan virus arke
      bitki=$(awk -F'\t' '$5=="33090"{print $2+0}' "$ins" | head -1); bitki=${bitki:-0}
      mantar=$(awk -F'\t' '$5=="4751"{print $2+0}' "$ins" | head -1); mantar=${mantar:-0}
      protozoa=$(awk -F'\t' '$5=="5794"||$5=="33682"{s+=$2} END{print s+0}' "$ins")
      insan=$(awk -F'\t' '$5=="9606"{print $2+0}' "$ins" | head -1); insan=${insan:-0}
      virus=$(awk -F'\t' '$5=="10239"{print $2+0}' "$ins" | head -1); virus=${virus:-0}
      arke=$(awk -F'\t' '$5=="2157"{print $2+0}' "$ins" | head -1); arke=${arke:-0}
      printf "    icerik: arke %s, virus %s, insan %s, mantar %s, protozoa %s, BITKI %s\n" \
             "$arke" "$virus" "$insan" "$mantar" "$protozoa" "$bitki"
      if   [ "$bitki" -gt 0 ] && [ "$mantar" -gt 0 ]; then o1="PlusPFP"
      elif [ "$mantar" -gt 0 ]; then o1="PlusPF (bitki YOK)"
      elif [ "$arke" -gt 0 ]; then o1="Standard (mantar ve bitki YOK)"
      else o1="ozel ya da tanimlanamayan veritabani"
      fi
    fi
  else
    echo "    kraken2-inspect yok, olcum 1 yapilamadi."
    echo "      Kurmak icin: micromamba install -n ${ORTAM:-mikro} -c bioconda kraken2"
  fi
  echo "    olcum 1 (icerik) : $o1"

  # --- iki olcumun uzlasisi ---
  echo
  case "$o1:$o2" in
    "olculemedi:"*)
      echo "    SONUC: BELIRSIZ. Yalnizca boyuta bakildi, icerik olculemedi."
      echo "    Boyut PlusPF ile PlusPFP'yi ayirmaya TEK BASINA yetmez."
      echo "    kraken2-inspect kurulmadan bu soru kapanmaz." ;;
    "PlusPFP:PlusPFP (tam)")
      echo "    SONUC: PlusPFP. Iki olcum de ayni seyi soyluyor." ;;
    "PlusPF (bitki YOK):PlusPF (tam)")
      echo "    SONUC: PlusPF. Iki olcum de ayni seyi soyluyor."
      echo "    DIKKAT: bu PlusPFP DEGIL. Bitki icermez. Beklediginiz PlusPFP ise"
      echo "    indirilen arsiv yanlis olabilir. PlusPFP adresi:"
      echo "      https://genome-idx.s3.amazonaws.com/kraken/  (k2_pluspfp_*.tar.gz)" ;;
    *)
      echo "    AYRILIK: icerik olcumu '$o1', boyut olcumu '$o2' diyor."
      echo "    Iki olcum ayrildi, birine sessizce gecilmiyor. Icerik olcumu daha"
      echo "    guvenilirdir ama once arsivin tam acildigini dogrulayin." ;;
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
  echo "ESKI KOSULAR HANGI VERITABANIYLA YAPILDI"
  echo
  echo "  a) Loglardan dogrudan kayit:"
  local bulundu=0
  for l in "$PROJE"/WSL/loglar/*.log; do
    [ -f "$l" ] || continue
    local yol boy
    yol=$(grep -m1 -oE "veritabani *: */[^ ]*" "$l" 2>/dev/null | sed 's/.*: *//' || true)
    boy=$(grep -m1 -oE "hash\.k2d *var *[0-9.]+[KMGT]" "$l" 2>/dev/null | grep -oE "[0-9.]+[KMGT]$" || true)
    if [ -n "$yol" ]; then
      printf "     %-46s  %-28s  hash %s\n" "$(basename "$l")" "$yol" "${boy:-?}"
      bulundu=1
    fi
  done
  [ "$bulundu" -eq 1 ] || echo "     log kaydi bulunamadi"
  echo
  echo "  b) Ciktilardan dolayli kanit (raporlarda hangi alanlar gorunuyor):"
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
    echo "     kraken_yeniden/tum.report: mantar ${m:-yok}, protozoa ${p:-yok}, bitki ${b:-yok}"
  fi
  echo
  echo "  BU KANITIN SINIRI, acikca yazilmali:"
  echo "  Raporda mantar ve protozoa GORUNMESI, veritabaninin onlari icerdigini"
  echo "  KANITLAR. Yani eski kosular en az PlusPF ile yapilmistir, Standard degil."
  echo "  Ama bitkinin GORUNMEMESI hicbir sey kanitlamaz: cururucu numunesinde"
  echo "  bitki okumasi atanmadiysa, veritabani bitki icerse bile rapora satir"
  echo "  yazilmaz. Yoklugun kaniti degildir. PlusPF ile PlusPFP ancak"
  echo "  kraken2-inspect ile ayrilir:  bash $0 vt-kimlik"
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
  echo "VERITABANI ARAMASI"
  echo "  aranan   : hash.k2d"
  echo "  kokler   : ${kokler[*]}"
  echo "  derinlik : $derin (butun disk TARANMAZ, makul surede bitmesi icin)"
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
  echo "  arama suresi: $((t1 - t0)) saniye"
  echo
  if [ "${#benzersiz[@]}" -eq 0 ]; then
    echo "  HICBIR KRAKEN2 VERITABANI BULUNAMADI."
    echo
    echo "  Yalnizca .tar.gz duruyor olabilir, acilmamis arsivler:"
    for k in "${kokler[@]}"; do
      [ -d "$k" ] || continue
      timeout 60 find "$k" -maxdepth "$derin" -name "k2_*.tar.gz" 2>/dev/null | head -5 \
        | while IFS= read -r t; do printf "    %s  %s\n" "$(du -h "$t" | cut -f1)" "$t"; done
    done
    echo
    echo "  Arsiv varsa acin (kendiliginden ACILMAZ):"
    echo "      cd <klasor> && tar -xzvf k2_*.tar.gz"
    echo "  Arsiv da yoksa indirin:"
    echo "      https://benlangmead.github.io/aws-indexes/k2"
    echo "      PlusPFP tam surum yaklasik 150 GB, indirilmis arsiv yaklasik 90 GB."
    echo "  Baska bir yerdeyse dogrudan verin:  VT_A=/tam/yol bash $0 $TUS"
    return 1
  fi
  echo "  ${#benzersiz[@]} veritabani bulundu:"
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
      echo "      EKSIK DOSYA:$eksik  (indeks kullanilamaz)"
    else
      echo "      toplam $(du -sh "$y" 2>/dev/null | cut -f1)"
      vt_surum "$y" | sed 's/^/      /'
    fi
    echo
  done
  if [ "${#benzersiz[@]}" -eq 1 ]; then
    echo "  Tek veritabani bulundu. Kullanmak icin:"
    echo "      VT_A=${benzersiz[0]} bash $0 esik"
  else
    echo "  BIRDEN FAZLA VERITABANI VAR. Hangisiyle kosulacagini betik SECMEZ,"
    echo "  cunku yanlis surumle kosarsak sonucu yanlis yorumlariz."
    echo "  Yukaridaki surum satirlarina bakip secin ve yolu verin:"
    for y in "${benzersiz[@]}"; do echo "      VT_A=$y bash $0 esik"; done
    echo "  Iki veritabanini karsilastirmak isterseniz:"
    echo "      VT_A=<genis> VT_B=<dar> bash $0 esik"
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
  echo "OZGUN KOSU HANGI VERITABANINI KULLANDI"
  echo "======================================================================"
  echo
  echo "KANIT 1: kaynak calismanin betiklerindeki --db argumani"
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
    echo "    -> betiklerin gosterdigi yol(lar):"
    printf '       %s\n' $k1
    if printf '%s' "$k1" | grep -qi "pluspf16\|pluspf_16\|_16gb"; then
      echo "       Klasor adi PlusPF-16 diyor: PlusPF'in 16 GB'a KAPAKLI surumu."
      echo "       Kapakli surum, tam surumle AYNI taksonlari icerir ama minimizer"
      echo "       sayisi seyreltilmistir. Ayni okuma daha ust bir dugume, hatta"
      echo "       baska bir klada dusebilir. Yani kapsam degil, COZUNURLUK dusuktur."
    fi
  else
    echo "    betiklerde veritabani yolu bulunamadi"
  fi

  echo
  echo "KANIT 2: kosu loglari"
  local k2=0
  for l in "$PROJE"/WSL/loglar/*.log; do
    [ -f "$l" ] || continue
    local yol boy
    yol=$(grep -m1 -oE "veritabani *: */[^ ]*" "$l" 2>/dev/null | sed 's/.*: *//' || true)
    boy=$(grep -m1 -oE "hash\.k2d *var *[0-9.]+[KMGT]" "$l" 2>/dev/null | grep -oE "[0-9.]+[KMGT]$" || true)
    [ -n "$yol" ] && { printf "    %-42s %-26s hash %s\n" "$(basename "$l")" "$yol" "${boy:-?}"; k2=1; }
  done
  [ "$k2" -eq 1 ] || echo "    ozgun kosuya ait log yok (loglar bizim kosularimiza ait)"

  echo
  echo "KANIT 3: rapor icerigi (hangi alanlar veritabaninda VARDI)"
  local r="$PROJE/SONUCLAR/kraken results/A1/edited_barcode01_kraken2.report"
  if [ -f "$r" ]; then
    local m p b v
    m=$(awk -F'\t' '$5=="4751"{print "VAR"}' "$r" | head -1)
    p=$(awk -F'\t' '$5=="5794"{print "VAR"}' "$r" | head -1)
    b=$(awk -F'\t' '$5=="33090"{print "VAR"}' "$r" | head -1)
    v=$(awk -F'\t' '$5=="10239"{print "VAR"}' "$r" | head -1)
    echo "    mantar (4751)      : ${m:-satir yok}"
    echo "    Apicomplexa (5794) : ${p:-satir yok}"
    echo "    bitki (33090)      : ${b:-satir yok}"
    echo "    virus (10239)      : ${v:-satir yok}"
    echo
    echo "    Mantar ve protozoa GORUNUYORSA veritabani onlari ICERIYORDU."
    echo "    Yani ozgun kosunun veritabani en az PlusPF kapsamindadir, Standard degil."
    echo "    Bitkinin gorunmemesi hicbir sey kanitlamaz: cururucu numunesinde"
    echo "    bitki okumasi atanmadiysa, veritabani bitki icerse bile satir yazilmaz."
    echo "    YOKLUK, KANIT DEGILDIR."
  else
    echo "    kaynak calismanin raporu bulunamadi: $r"
  fi

  echo
  echo "KANIT 4: bizim yeniden kosumuz kaynak calismanin sonucunu tekrarladi mi"
  local oz="$PROJE/SONUCLAR/kraken_yeniden/kraken_ozet.csv"
  if [ -f "$oz" ]; then
    local top uy ay
    top=$(awk -F',' 'NR>1{n++} END{print n+0}' "$oz")
    uy=$(awk -F',' 'NR>1 && $9=="uyusuyor"{n++} END{print n+0}' "$oz")
    ay=$((top - uy))
    echo "    $top kutunun $uy tanesinde ayni sonuc, $ay tanesinde FARKLI sonuc."
    echo
    echo "    BU BIR CIKARIM URETIR. Kraken2 belirlenimcidir: ayni okuma, ayni"
    echo "    veritabani, ayni esik her zaman ayni atamayi verir. Kutular zaten"
    echo "    kaynak calismanin ciktisina gore ayrildigina gore, veritabani AYNI olsaydi"
    echo "    butun kutular kendi etiketini tekrarlardi."
    if [ "$ay" -gt 0 ]; then
      echo "    $ay kutu tekrarlamadi, dolayisiyla bizim yeniden kosumuzun"
      echo "    veritabani ozgun kosunun veritabaniyla AYNI DEGILDIR."
      echo
      echo "    Bu cikarimin sinirlari, durustce:"
      echo "      extract_kraken_reads.py --include-children ile calisti, yani bir"
      echo "      kutuda o taksonun ALT dugumlerine ait okumalar da var. Tur duzeyi"
      echo "      toplama bunun cogunu kapatir ama tamamini degil. Buna karsilik"
      echo "      cins ve aile atlayan ayrismalar (ornegin Bacteroides kutusunun"
      echo "      Candidatus Azobacteroides cikmasi) alt dugumle aciklanamaz."
    else
      echo "    Butun kutular tekrarlandi. Bu, iki kosunun ayni ya da cok benzer"
      echo "    veritabanini kullandigina isaret eder."
    fi
  else
    echo "    kraken_ozet.csv yok, bu kanit olculemedi"
  fi

  echo
  echo "======================================================================"
  echo "SONUC"
  echo "======================================================================"
  if printf '%s' "$k1" | grep -qi "pluspf16\|pluspf_16\|_16gb"; then
    echo "  ozgun kosunun veritabani: PlusPF-16 (PlusPF'in 16 GB'a kapakli surumu),"
    echo "  betiginde yazili yola gore. Bu KESIN degil, betikteki yol ile fiilen"
    echo "  kosulan sey ayrilabilir; ama elimizdeki en dogrudan kanit budur ve"
    echo "  rapor icerigiyle (mantar ve protozoa var) celismiyor."
  else
    echo "  BELIRLENEMEDI. Elimizdeki kanit yalnizca su kadarini soyluyor:"
    echo "  veritabani mantar ve protozoa iceriyordu, yani en az PlusPF"
    echo "  kapsamindaydi. Kapakli mi tam mi, PlusPFP mi degil mi, cikarilamadi."
  fi
  echo
  echo "  BU BELIRSIZLIK YORUMU ETKILER, acikca yazilmali:"
  echo "  Eski kosu zaten genis bir veritabaniyla yapildiysa, 'sorun kapsamdi'"
  echo "  teshisi zayiflar ve sebebin nanopore hata orani olma ihtimali one cikar."
  echo "  Iki aciklamayi ayirt etmenin yolu yine de PlusPFP kosusudur:"
  echo "    kapsam sorunuysa  -> PlusPFP'de atamalar guclenir ve esige dayanir"
  echo "    okuma hatasiysa   -> PlusPFP'de de atamalar zayif kalir ve esikle coker"
  echo "  Tablonun yorumu bu belirsizligi tasiyacak, gizlemeyecek."
  set -e
}

# --- bellek. ozgun Kraken/Bracken betigi ve rerun_kraken.sh'teki karar ---------
bellek_bayragi() {
  local d="$1"
  local hb rb
  hb=$(stat -c%s "$d/hash.k2d" 2>/dev/null || echo 0)
  rb=$(( $(awk '/MemAvailable/{print $2}' /proc/meminfo 2>/dev/null || echo 0) * 1024 ))
  {
    echo "  hash.k2d       : $(numfmt --to=iec "$hb" 2>/dev/null || echo "$hb")"
    echo "  kullanilabilir : $(numfmt --to=iec "$rb" 2>/dev/null || echo "$rb")"
  } >&2
  if [ "${ZORLA_MMAP:-1}" = "1" ] || [ "$hb" -gt $((rb - 2147483648)) ]; then
    {
      echo "  --memory-mapping KULLANILACAK. Veritabani RAM'e yuklenmez, diskten"
      echo "  okunur. RAM ihtiyaci duser, kosu YAVASLAR. Siniflandirma sonucu"
      echo "  AYNIDIR, degisen tek sey suredir."
    } >&2
    echo "--memory-mapping"
  else
    echo "  RAM yeterli, veritabani bellege yuklenecek (hizli)." >&2
    echo ""
  fi
}

# --- okumalari tek dosyada birlestir (rerun_kraken.sh'ten aynen) --------
birlestir() {
  local hedef="$1"
  local eski="$PROJE/SONUCLAR/kraken_yeniden/tum.fastq"
  if [ -s "$hedef" ]; then
    echo "birlestirilmis dosya zaten var: $(awk 'END{print int(NR/4)}' "$hedef") okuma"
    return
  fi
  mkdir -p "$(dirname "$hedef")"
  if [ -s "$eski" ]; then
    echo "rerun_kraken.sh'in urettigi tum.fastq bulundu, AYNI okuma kumesi kullanilacak."
    echo "  (Esikleri ve veritabanlarini karsilastirilabilir kilan sart budur.)"
    cp "$eski" "$hedef"
    cp "$PROJE/SONUCLAR/kraken_yeniden/kaynak_sayim.tsv" "$(dirname "$hedef")/" 2>/dev/null || true
    echo "  $(awk 'END{print int(NR/4)}' "$hedef") okuma"
    return
  fi
  [ -d "$KAYNAK" ] || { echo "HATA: fastq klasoru yok: $KAYNAK"; exit 1; }
  local dosyalar sayi
  dosyalar=$(ls "$KAYNAK"/*/*reads_*.fastq 2>/dev/null || true)
  sayi=$(printf '%s\n' "$dosyalar" | grep -c . || true)
  [ "$sayi" -gt 0 ] || { echo "HATA: hic fastq bulunamadi: $KAYNAK"; exit 1; }
  echo "$sayi fastq dosyasi birlestiriliyor (takson basina en fazla ${KAP} okuma; 0 = hepsi)"
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
  [ "$toplam" = "$beklenen" ] || { echo "HATA: okuma sayisi tutmuyor ($toplam / $beklenen)"; exit 1; }
  echo "  $toplam okuma birlestirildi ve sayisi dogrulandi"
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
    echo "  ornek dosyasi zaten var: $(awk 'END{print int(NR/4)}' "$hedef") okuma"
    return
  fi
  local n_kaynak; n_kaynak=$(awk 'END{print int(NR/4)}' "$kaynak")
  if [ "$toplam" -le 0 ] || [ "$n_kaynak" -le "$toplam" ]; then
    echo "  ornekleme GEREKMEDI (kaynak $n_kaynak okuma, hedef $toplam)"
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
    echo "  HATA: okuma basliklarinda 'tx<taxid>_' oneki yok (%${onekli:-0})."
    echo "  Kutu basina esit ornekleme YAPILAMAZ; sessizce yanli bir ornek uretmektense"
    echo "  is durduruluyor. Cozum: $(dirname "$kaynak")/tum.fastq dosyasini silip"
    echo "  yeniden urettirin, ya da ORNEK=0 ile tam veriyle kosun."
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
  echo "  ornekleme: $kutu_sayisi kutu x $pay okuma = hedef ~$toplam"
  echo "  yontem   : kutu ici RASTGELE, sabit tohum $TOHUM (tekrarlanabilir)"

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
  echo "  ornek hazir: $n okuma ($n_kaynak icinden)"
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
  echo "SURE OLCUMU. Kucuk bir ornek gercekten kosulur, tahmin uydurulmaz."
  vt_hazir_mi "$d" || { echo "Veritabani hazir degil, olcum yapilamaz."; exit 1; }
  local BAYRAK; BAYRAK=$(bellek_bayragi "$d")
  local tmp; tmp=$(mktemp -d)
  birlestir "$tmp/tum.fastq" >/dev/null
  local TOPLAM; TOPLAM=$(awk 'END{print int(NR/4)}' "$tmp/tum.fastq")
  local N="${DENEME:-2000}"
  head -n $((N*4)) "$tmp/tum.fastq" > "$tmp/deneme.fastq"
  local gercek; gercek=$(awk 'END{print int(NR/4)}' "$tmp/deneme.fastq")
  echo
  echo "  deneme: $gercek okuma (tam kume $TOPLAM okuma)"
  local t0 t1 sn
  t0=$(date +%s)
  kraken2 --db "$d" --threads "$IPLIK" $BAYRAK --confidence 0 \
          --report "$tmp/d.report" --output "$tmp/d.out" "$tmp/deneme.fastq" \
          >/dev/null 2>"$tmp/d.err" || {
    echo "  kraken2 HATA:"; sed 's/^/    /' "$tmp/d.err" | tail -5; rm -rf "$tmp"; exit 1; }
  t1=$(date +%s); sn=$((t1 - t0)); [ "$sn" -lt 1 ] && sn=1
  local tam esik_adet
  tam=$(( sn * TOPLAM / gercek ))
  esik_adet=$(echo $ESIKLER | wc -w)
  echo
  echo "  olculen  : $sn saniye / $gercek okuma"
  echo "  tam kume : yaklasik $(( tam / 60 )) dakika (tek esik)"
  echo "  $esik_adet esik: yaklasik $(( tam * esik_adet / 60 )) dakika = $(( tam * esik_adet / 3600 )) saat"
  echo
  echo "  BU SAYININ SINIRI: --memory-mapping ile ilk okumalar en yavas olandir,"
  echo "  sayfalar diskten gelir. Isletim sistemi onbellegi isindikca hizlanir."
  echo "  Yani bu tahmin bir UST SINIRDIR, gercek sure buna esit ya da daha kisa."
  echo "  Ters yonde tek risk: makine baska is yapip RAM'i bosaltirsa onbellek"
  echo "  soguyup sure uzayabilir. Kosu sirasinda makineyi mesgul etmeyin."
  rm -rf "$tmp"
}

# =========================================================================
# TUS: durum
# =========================================================================
tus_durum() {
  echo "ORTAM DENETIMI. Hicbir sey kosulmaz, yalnizca bakilir."
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
    echo "VT_A kullanilamiyor. Diskte araniyor (pes edilmiyor):"
    echo
    vt_ara || true
  fi
  if [ -n "$VT_B" ]; then
    echo
    echo "VT_B = $VT_B"
    vt_hazir_mi "$VT_B" || true
    if [ "$(readlink -f "$VT_A" 2>/dev/null)" = "$(readlink -f "$VT_B" 2>/dev/null)" ]; then
      echo
      echo "  UYARI: VT_A ile VT_B AYNI klasor. Iki egri cizilemez, tek egri cizilir."
    fi
  else
    echo
    echo "VT_B = (verilmedi). Tek veritabaniyla calisilacak."
    echo "  Iki veritabanini karsilastirmak icin:  VT_B=/yol/eski_db bash $0 esik"
  fi
  echo
  local n; n=$(ls "$KAYNAK"/*/*reads_*.fastq 2>/dev/null | wc -l || echo 0)
  echo "okumalar: $n fastq dosyasi ($KAYNAK)"
  echo
  free -g 2>/dev/null | awk '/Mem:/{print "RAM: "$2" GB toplam, "$7" GB kullanilabilir"}'
  df -h "$VT_A" 2>/dev/null | tail -1 | awk '{print "disk: "$4" bos ("$6")"}'
  echo
  eski_kosu_tespit
  echo
  echo "Sonraki adim:  bash $0 vt-kimlik    (veritabani hangi surum)"
}

tus_vt_kimlik() {
  log_ac vt_kimlik
  kraken_sart
  echo
  echo "VERITABANI SURUM TESPITI"
  echo "Butun argumanimiz veritabaninin KAPSAMI uzerine kurulu. Yanlis surumle"
  echo "kosarsak sonucu yanlis yorumlariz, o yuzden bu adim atlanmaz."
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
  echo "SELFTESTLER. Sinav gecmeden ana is baslamaz (proje kurali 2)."
  local hata=0
  for p in threshold_summary.py comparison_table.py custom_taxonomy.py ozgun Kraken ozet betigi; do
    echo; echo "--- $p"
    if [ -f "$_BETIK_DIZIN/$p" ]; then
      python3 "$_BETIK_DIZIN/$p" --selftest || hata=1
    else echo "  DOSYA YOK: $p"; hata=1; fi
  done
  echo
  [ "$hata" -eq 0 ] && echo "BUTUN SINAVLAR GECTI" || { echo "SINAV KALDI"; return 1; }
}

# =========================================================================
# TUS: esik-a / esik-b   guven esigi taramasi
# =========================================================================
esik_tara() {
  local d="$1" is="$2" etiket="$3"
  echo
  echo "======================================================================"
  echo "GUVEN ESIGI TARAMASI, $etiket"
  echo "  veritabani: $d"
  echo "  esikler   : $ESIKLER"
  echo "  is parcacigi: $IPLIK   (tepe bellegi dusurmek icin dusuk tutuluyor)"
  echo "======================================================================"
  vt_hazir_mi "$d" || {
    echo
    echo "Veritabani hazir degil. Pes edilmiyor, diskte araniyor."
    echo
    vt_ara || true
    echo
    echo "HATA: '$d' kullanilamadi, is duruyor. Yukaridan bir yol secip verin."
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
    echo "ORNEKLEME ACIK (ORNEK=$ORNEK). Tam veri icin: ORNEK=0 bash $0 esik"
    ornekle "$is/tum.fastq" "$is/ornek.fastq" "$ORNEK"
    GIRDI="$is/ornek.fastq"
  fi
  local OKUMA; OKUMA=$(awk 'END{print int(NR/4)}' "$GIRDI")
  echo "  taramada kullanilacak okuma: $OKUMA"
  echo
  for C in $ESIKLER; do
    local ad="esik_${C}"
    if [ -s "$is/${ad}.report" ]; then
      echo "[$C] zaten var, atlaniyor"; continue
    fi
    local t0 t1
    t0=$(date +%s)
    echo "[$C] kraken2 --confidence $C  baslangic $(date '+%H:%M:%S')"
    kraken2 --db "$d" --threads "$IPLIK" $BAYRAK --confidence "$C" \
            --use-names --report "$is/${ad}.report" --output "$is/${ad}.out" \
            "$GIRDI" 2>"$is/${ad}.err" >/dev/null || {
      echo "    kraken2 HATA:"; tail -5 "$is/${ad}.err" | sed 's/^/      /'
      echo "    Bu esik atlanmiyor, is duruyor. Yarim tarama yaniltir."; exit 1; }
    rm -f "$is/${ad}.err"
    t1=$(date +%s)
    local n; n=$(wc -l < "$is/${ad}.out")
    [ "$n" = "$OKUMA" ] || { echo "HATA: cikti satiri $n, okuma $OKUMA, tutmuyor"; exit 1; }
    echo "    bitti, $n okuma, $(( (t1-t0)/60 )) dakika $(( (t1-t0)%60 )) saniye"
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
    echo "UYARI: VT_A ile VT_B ayni klasor. Iki egri cizmenin anlami yok,"
    echo "tek tarama yapilacak. Gercekten iki farkli veritabani karsilastirmak"
    echo "istiyorsaniz ikincisinin yolunu verin: VT_B=/yol/oteki bash $0 esik"
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
  echo "bitti. Dosyalar: $IS_A"
  echo "  esik_<C>.report          her esigin kraken raporu"
  echo "  esik_egrisi.csv / .txt   alan bazinda atama yuzdeleri, esige gore"
  # NOT: son komut olarak "[ ... ] && echo" birakilmaz. Kosul yanlis oldugunda
  # test 1 doner ve fonksiyonun cikis kodu olur; is basariyla bittigi halde
  # betik HATA vermis gibi gorunur. Bu tuzaga bir kez dusuldu.
  if [ "$ikili" -eq 1 ]; then
    echo "  esik_iki_veritabani.txt  iki veritabani yan yana + ayakta kalma"
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
  local hedef="$PROJE/TESLIM/KRAKEN_KARSILASTIRMA.md"
  if [ -f "$hedef" ]; then
    local yed="$hedef.yedek_$(date +%Y%m%d_%H%M%S)"
    cp -p "$hedef" "$yed" 2>/dev/null && \
      echo "  mevcut tablo yedeklendi: $(basename "$yed")"
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
  echo "OZEL KRAKEN2 VERITABANI KURULUMU  (EN SON CARE)"
  echo "======================================================================"
  echo
  echo "ONCE BUNU OKUYUN: PlusPFP kuruluysa BU ADIMA GEREK YOKTUR."
  echo "PlusPFP zaten Standard'a protozoa, mantar ve bitki ekler; yani eksik"
  echo "oldugunu olctugumuz gruplarin tamami. Bu tus, PlusPFP hicbir sekilde"
  echo "kurulamiyorsa ya da marker gen (16S/ITS) duzeyinde ikinci bir gorus"
  echo "isteniyorsa vardir. Kurulum saatler surer."
  echo
  echo "Kontrol:  bash $0 vt-kimlik    (PlusPFP kuruluysa buraya hic gelmeyin)"
  echo
  local toplam=0 var=0
  for k in ${KUMELER:-silva_ssu unite pr2}; do
    local f; f=$(kume_dosya "$k") || continue
    if [ -s "$f" ]; then
      toplam=$((toplam + $(stat -c%s "$f"))); var=$((var+1))
      printf "  %-12s %-42s %s\n" "$k" "$(basename "$f")" "$(du -h "$f" | cut -f1)"
    else
      printf "  %-12s %-42s YOK\n" "$k" "$(basename "$f")"
    fi
  done
  [ "$var" -gt 0 ] || { echo; echo "HATA: hicbir referans kumesi bulunamadi ($PROJE/REFERANS_DB)"; exit 1; }
  local gb=$((toplam / 1073741824)); [ "$gb" -lt 1 ] && gb=1
  echo
  echo "  toplam dizi  : yaklasik ${gb} GB"
  echo "  TAHMINI SURE : $(( gb*40 + 30 )) ile $(( gb*120 + 60 )) dakika"
  echo "  TAHMINI RAM  : en az $(( gb*6 + 8 )) GB"
  echo "  TAHMINI DISK : yaklasik $(( gb*5 + 10 )) GB"
  echo "  Hedef        : $OZELVT"
  echo
  echo "  Bu bir MARKER GEN veritabanidir (16S/18S/ITS), tam genom icermez."
  echo "  Taksonomi NCBI'dan degil dosyalarin kendi soy dizgilerinden uretilir;"
  echo "  taxid'ler sentetiktir ve NCBI taxid'leriyle karsilastirilamaz, tablo"
  echo "  bu yuzden isim duzeyinde karsilastirir."
  local mem; mem=$(free -g 2>/dev/null | awk '/Mem:/{print $2}' || echo 0)
  echo "  Bu makinede toplam RAM: ${mem} GB"
  [ "${mem:-0}" -lt $(( gb*6 + 8 )) ] && echo "  UYARI: RAM tahmini ihtiyacin altinda."
  echo
  if [ "${ONAY:-}" != "evet" ]; then
    read -r -p "Kuruluma baslansin mi? (evet yazin, baska her sey iptal): " c
    [ "$c" = "evet" ] || { echo "iptal edildi. Hicbir sey yapilmadi."; exit 0; }
  fi
  log_ac ozelvt_kur
  kraken_sart
  command -v kraken2-build >/dev/null 2>&1 || {
    echo "HATA: kraken2-build bulunamadi."
    echo "  micromamba install -n ${ORTAM:-mikro} -c bioconda kraken2"; exit 1; }
  mkdir -p "$OZELVT/library" "$OZELVT/taxonomy"
  echo; echo "1/2  soy dizgilerinden taksonomi ve kutuphane uretiliyor"
  local args=()
  for k in ${KUMELER:-silva_ssu unite pr2}; do
    local f; f=$(kume_dosya "$k") || continue
    [ -s "$f" ] && args+=(--kume "$k=$f")
  done
  python3 "$_BETIK_DIZIN/custom_taxonomy.py" --cikti "$OZELVT" "${args[@]}"
  echo; echo "2/2  kraken2-build --build  ($(date '+%H:%M:%S'))"
  kraken2-build --build --db "$OZELVT" --threads "$IPLIK"
  echo; echo "kurulum bitti: $OZELVT"
  du -sh "$OZELVT" | awk '{print "  boyut: "$1}'
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
    echo "HATA: ozel veritabani kurulu degil ($OZELVT/hash.k2d yok)."
    echo "  Once:  bash $0 ozelvt-kur"; exit 1; }
  mkdir -p "$OZEL_IS"
  birlestir "$OZEL_IS/tum.fastq"
  local BAYRAK; BAYRAK=$(bellek_bayragi "$OZELVT")
  local C="${OZEL_ESIK:-0}"
  kraken2 --db "$OZELVT" --threads "$IPLIK" $BAYRAK --confidence "$C" \
          --use-names --report "$OZEL_IS/tum.report" --output "$OZEL_IS/tum.out" \
          "$OZEL_IS/tum.fastq" 2>"$OZEL_IS/hata.txt" >/dev/null || {
    echo "kraken2 HATA:"; tail -5 "$OZEL_IS/hata.txt" | sed 's/^/  /'; exit 1; }
  rm -f "$OZEL_IS/hata.txt"
  python3 "$_BETIK_DIZIN/ozgun Kraken ozet betigi" --is "$OZEL_IS" --kok "$PROJE" || true
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
  echo "WSL BELLEK AYARI - Windows'un kilitlenmesini onler"
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
    KAYNAK="WSL icinden tahmin (Windows toplami DAHA BUYUK olabilir)"
  fi
  local TOPLAM_GB=$(( TOPLAM_MB / 1024 ))
  # Yuzde 60: WSL'e yeter, Windows'a nefes payi birakir. Alt sinir 4 GB,
  # cunku 110 GB'lik veritabani mmap ile bile bir taban bellek ister.
  local WSL_GB=$(( TOPLAM_GB * 60 / 100 ))
  [ "$WSL_GB" -lt 4 ] && WSL_GB=4
  local SWAP_GB=$(( WSL_GB / 2 ))
  [ "$SWAP_GB" -lt 4 ] && SWAP_GB=4

  echo "  Toplam RAM   : ${TOPLAM_GB} GB   ($KAYNAK)"
  echo "  WSL'e onerilen: ${WSL_GB} GB     (toplamin yuzde 60'i)"
  echo "  Swap         : ${SWAP_GB} GB"
  echo

  local ICERIK
  ICERIK="[wsl2]
memory=${WSL_GB}GB
swap=${SWAP_GB}GB
processors=${IPLIK}
# Sayfa onbellegini sinirlar; 110 GB'lik veritabani mmap ile okunurken
# onbellegin sanal makineyi sisirmesini engeller.
pageReporting=true
# Bosta kalan bellegi Windows'a geri verir.
autoMemoryReclaim=gradual"

  local CIKTI="$PROJE/wslconfig_ONERILEN.txt"
  printf '%s\n' "$ICERIK" > "$CIKTI"

  echo "----------------------------------------------------------------------"
  printf '%s\n' "$ICERIK"
  echo "----------------------------------------------------------------------"
  echo
  echo "NE YAPACAKSINIZ - sirayla, uc adim"
  echo
  echo "  1) Yukaridaki metin su dosyaya yazildi:"
  echo "       $CIKTI"
  echo "     Windows tarafinda gorunen adi:"
  echo "       C:\\Users\\yerli\\Masaustu\\PROJE\\wslconfig_ONERILEN.txt"
  echo
  echo "  2) O dosyayi SU ADA ve SU YERE kopyalayin:"
  echo "       C:\\Users\\yerli\\.wslconfig"
  echo "     (dosya adi nokta ile baslar ve uzantisi YOKTUR)"
  echo "     Kisa yol: Windows'ta Calistir (Win+R) acip sunu yapistirin:"
  echo "       notepad C:\\Users\\yerli\\.wslconfig"
  echo "     Bos dosya acilirsa metni yapistirip kaydedin."
  echo
  echo "  3) WSL'i kapatip acin. Windows'ta PowerShell ya da Komut Istemi:"
  echo "       wsl --shutdown"
  echo "     Sonra WSL'i yeniden acin. Ayar ancak bundan sonra gecerlidir."
  echo
  echo "DOGRULAMA: yeniden actiktan sonra WSL icinde su komutu calistirin"
  echo "       free -g"
  echo "  'total' sutunu yaklasik ${WSL_GB} GB gostermelidir."
  echo
  echo "BU AYAR TEK BASINA KILITLENMEYI COZER. WSL yavaslar, Windows kullanilabilir"
  echo "kalir. Esik taramasi ayrica ornekleme ile hafifletildi (asagi bakin)."
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
  vt_hazir_mi "$d" || { echo "HATA: veritabani hazir degil"; exit 1; }
  mkdir -p "$is"
  birlestir "$is/tum.fastq"
  ornekle "$is/tum.fastq" "$is/ornek.fastq" "$ORNEK"
  local BAYRAK; BAYRAK=$(bellek_bayragi "$d")
  local a="$is/dogrulama_ornek" b="$is/dogrulama_tam"
  for cift in "ornek.fastq:$a" "tum.fastq:$b"; do
    local gir="${cift%%:*}" cik="${cift##*:}"
    if [ -s "${cik}.report" ]; then echo "$(basename $cik) zaten var, atlaniyor"; continue; fi
    echo "[$C] $(basename $gir) siniflandiriliyor  $(date '+%H:%M:%S')"
    kraken2 --db "$d" --threads "$IPLIK" $BAYRAK --confidence "$C" \
            --use-names --report "${cik}.report" --output /dev/null \
            "$is/$gir" 2>"${cik}.err" >/dev/null || {
      echo "  kraken2 HATA:"; tail -3 "${cik}.err" | sed 's/^/    /'; exit 1; }
    rm -f "${cik}.err"
  done
  echo
  echo "TEMSIL KARSILASTIRMASI  (esik $C)"
  echo "  Bakilan: alan duzeyi yuzdeleri. Mutlak sayilar farkli olacak, ONEMLI DEGIL."
  echo
  printf '  %-24s %10s %10s %10s\n' "alan" "ornek %" "tam %" "fark"
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
  echo "  yazildi: $is/temsil_dogrulamasi.tsv"
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
  esik-a)      tus_sinav >/dev/null || { echo "SINAV KALDI. Ayrinti: bash $0 sinav"; exit 2; }
               log_ac esik_a; kraken_sart; esik_tara "$VT_A" "$IS_A" "VT_A"
               python3 "$_BETIK_DIZIN/threshold_summary.py" --kok "$PROJE" --is "$IS_A" --ad "$(basename "$VT_A")" ;;
  esik-b)      [ -n "$VT_B" ] || { echo "HATA: VT_B verilmedi.  VT_B=/yol bash $0 esik-b"; exit 1; }
               tus_sinav >/dev/null || { echo "SINAV KALDI"; exit 2; }
               log_ac esik_b; kraken_sart; esik_tara "$VT_B" "$IS_B" "VT_B"
               python3 "$_BETIK_DIZIN/threshold_summary.py" --kok "$PROJE" --is "$IS_B" --ad "$(basename "$VT_B")" ;;
  esik)        tus_sinav >/dev/null || { echo "SINAV KALDI. Ayrinti: bash $0 sinav"; exit 2; }
               tus_esik ;;
  tablo)       tus_tablo ;;
  hepsi)       tus_sinav >/dev/null || { echo "SINAV KALDI"; exit 2; }
               tus_esik; tus_tablo ;;
  ozelvt-kur)  tus_ozelvt_kur ;;
  ozelvt-kos)  tus_ozelvt_kos ;;
  *)           sed -n '3,48p' "$0" | sed 's/^# \{0,1\}//'
               echo; echo "Ayrinti icin: $PROJE/NASIL_DEVAM_EDILIR.md" ;;
esac
