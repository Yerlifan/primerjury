#!/bin/bash
# ---------------------------------------------------------------------------
# BUTUN TAKSONLARIN KRAKEN2 ILE YENIDEN SINIFLANDIRILMASI
#
# NEDEN HEPSI
# Yedi okaryotun yanlis etiketlendigi biliniyor, kimlikleri dolayli yoldan bulundu.
# Yalnizca o yediyi yeniden sinirlandirmak, zaten supheli olanlari teyit etmekten
# ibaret kalir ve YENI bir yanlis etiketi bulamaz. Bir olcut, bakmadigi yerdeki
# hatayi gecirir; bu projede bugune kadarki on hatanin ortak sebebi buydu.
# Bu yuzden fastq klasorundeki BUTUN taksonlar taranir, hicbiri secilmez.
#
# NEDEN TEK KOSU
# kraken2 her cagrida veritabanini bastan yukler. Takson basina ayri cagri, ayni
# veritabanini onlarca kez yuklemek demektir ve is neredeyse tamamen yuklemeye
# gider. Bunun yerine butun okumalar tek dosyada birlestirilir, her okumanin adina
# kaynak taxid yazilir, kraken2 BIR KEZ calisir, sonuc sonradan taksona bolunur.
# Kraken2 her okumayi bagimsiz siniflandirdigi icin birlestirme sonucu degistirmez.
#
# Calistirma:
#   bash rerun_kraken.sh
#   bash rerun_kraken.sh /baska/k2db/yolu
#   IPLIK=8 bash rerun_kraken.sh
# ---------------------------------------------------------------------------
set -euo pipefail

# ALI KLASORU OTOMATIK BULUNUR, SABIT YAZILMAZ.
# Bu betik tools/WSL icinde duruyor; bir ust klasor ALI'dir. BASH_SOURCE
# calisan dosyanin gercek yolu oldugu icin bu tahmin degil olcumdur.
# Kurulum tasindiysa disaridan bastirilabilir:  ALI=/tam/yol/ALI bash <betik>
_BETIK_DIZIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ALI="${ALI:-$(cd "$_BETIK_DIZIN/.." && pwd)}"
if [ ! -d "$tools/0_TESLIM_RAPOR" ]; then
  echo "HATA: ALI klasoru dogrulanamadi ('$ALI' icinde 0_TESLIM_RAPOR yok)."
  echo "  Betik konumu: $_BETIK_DIZIN"
  echo "  Elle vermek icin:  ALI=/tam/yol/ALI bash $0"
  exit 1
fi
# KRAKEN2 VERITABANI YOLU DA SABIT YAZILMAZ.
# Isaret dosyasi hash.k2d aranarak bulunur. Ilk arguman verilirse o kazanir.
DB="${1:-}"
if [ -z "$DB" ]; then
  for _k in "$HOME"/*/hash.k2d "$HOME"/*/*/hash.k2d /home/*/*/hash.k2d \
            /opt/*/hash.k2d /mnt/c/*/hash.k2d; do
    [ -f "$_k" ] && DB="$(dirname "$_k")" && break
  done
fi
if [ -z "$DB" ] || [ ! -f "$DB/hash.k2d" ]; then
  echo "HATA: Kraken2 veritabani bulunamadi (hash.k2d aranmasina ragmen)."
  echo "  Elle vermek icin:  bash $0 /tam/yol/k2db"
  exit 1
fi
echo "Kraken2 veritabani bulundu: $DB  (sabit yazilmadi, hash.k2d arandi)"
IS="$tools/SONUCLAR/kraken_yeniden"
KAYNAK="$tools/SONUCLAR/fastq files"
IPLIK="${IPLIK:-12}"
BURASI="$(cd "$(dirname "$0")" && pwd)"

echo "veritabani : $DB"
echo "kaynak     : $KAYNAK"
echo "cikti      : $IS"
echo "iplik      : $IPLIK"
echo

# --- ortam. kraken2 micromamba'nin "mikro" ortaminda, PATH'te degil ------------
# install.sh araclar onu micromamba ile kuruyor; yeni bir kabukta ortam etkin olmadigi
# icin kraken2 gorunmez. Burada aynisi yapilir. ORTAM degiskeniyle baska bir
# ortam adi verilebilir.
ORTAM="${ORTAM:-mikro}"
if ! command -v kraken2 >/dev/null 2>&1; then
  echo "kraken2 PATH'te yok, micromamba ortami '$ORTAM' etkinlestiriliyor"
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

# --- on kontroller. Sessizce eksik calismaktansa burada durmak yeglenir --------
if ! command -v kraken2 >/dev/null 2>&1; then
  echo "HATA: kraken2 hala bulunamadi."
  echo "  Ortami elle etkinlestirip tekrar deneyin:"
  echo "    export PATH=\"\$HOME/bin:\$PATH\""
  echo "    export MAMBA_ROOT_PREFIX=\"\$HOME/micromamba\""
  echo "    eval \"\$(micromamba shell hook -s bash)\""
  echo "    micromamba activate $ORTAM"
  echo "  Ortam adi farkliysa: micromamba env list"
  exit 1
fi
echo "kraken2: $(command -v kraken2)"
[ -f "$DB/hash.k2d" ] || { echo "HATA: $DB/hash.k2d yok, veritabani yolu yanlis"; exit 1; }
[ -d "$KAYNAK" ]      || { echo "HATA: fastq klasoru yok: $KAYNAK"; exit 1; }
[ -f "$BURASI/ozgun Kraken ozet betigi" ] || { echo "HATA: ozgun Kraken ozet betigi yok"; exit 1; }

mkdir -p "$IS"
cd "$IS"

DOSYALAR=$(ls "$KAYNAK"/*/*reads_*.fastq 2>/dev/null || true)
DOSYA_SAYISI=$(printf '%s\n' "$DOSYALAR" | grep -c . || true)
[ "$DOSYA_SAYISI" -gt 0 ] || { echo "HATA: hic fastq bulunamadi"; exit 1; }

# --- bellek. Veritabani RAM'e sigmiyorsa disk uzerinden okunur ----------------
# kraken2 varsayilan olarak hash tablosunun TAMAMINI RAM'e yukler. PlusPF gibi tam
# veritabanlari 100 GB'i asar; 16 GB'lik bir makinede bu mumkun degildir ve kraken2
# "unable to allocate hash table memory" deyip durur. --memory-mapping ile tablo
# diskten okunur, RAM sarti kalkar, karsiliginda is yavaslar. Sonuc AYNIDIR,
# siniflandirma degismez, yalnizca sure uzar.
DB_BAYT=$(du -sb "$DB" 2>/dev/null | cut -f1 || echo 0)
RAM_KB=$(awk '/MemAvailable/{print $2}' /proc/meminfo 2>/dev/null || echo 0)
RAM_BAYT=$((RAM_KB * 1024))
echo "veritabani boyutu : $(numfmt --to=iec "$DB_BAYT" 2>/dev/null || echo "$DB_BAYT") "
echo "kullanilabilir RAM: $(numfmt --to=iec "$RAM_BAYT" 2>/dev/null || echo "$RAM_BAYT")"
BELLEK_BAYRAGI=""
if [ "$DB_BAYT" -gt 0 ] && [ "$DB_BAYT" -gt "$RAM_BAYT" ]; then
  BELLEK_BAYRAGI="--memory-mapping"
  echo "veritabani RAM'e SIGMIYOR, --memory-mapping kullanilacak (disk uzerinden, yavas)"
else
  echo "veritabani RAM'e siginiyor, normal mod"
fi

# --- 1. butun okumalar tek dosyada, adlarina kaynak taxid yazilarak -----------
# KAP degiskeni verilirse her taksondan en fazla o kadar okuma alinir. Bu bir
# TAKSON secimi DEGILDIR: 44 taksonun hepsi yine taranir, yalnizca her kutudan
# okuma orneklenir. Kutunun kimligini belirlemek icin birkac bin okuma yeter,
# --memory-mapping ile sure okuma sayisiyla dogru orantili oldugu icin bu ilk
# turda saatler kazandirir. KAP verilmezse butun okumalar kullanilir.
KAP="${KAP:-0}"
[ "$KAP" -gt 0 ] && echo "her taksondan en fazla $KAP okuma alinacak (butun taksonlar taranir)"
# Okuma adi "@tx<taxid>_<eski ad>" olur. Kraken2 ciktisinin ikinci alani bu addir,
# sonucu taksona bolmek icin kullanilir.
echo "$DOSYA_SAYISI fastq dosyasi birlestiriliyor"
: > tum.fastq
: > kaynak_sayim.tsv
: > /tmp/kraken_alinan.tsv
while IFS= read -r f; do
  [ -n "$f" ] || continue
  TX=$(basename "$f" | sed -E 's/.*reads_([0-9]+)\.fastq$/\1/')
  N=$(awk 'END{print int(NR/4)}' "$f")
  # KAP her TAKSON icin toplamdir, dosya basina degil; bir taksonun okumalari
  # birden cok numune dosyasina dagilmis olabilir.
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
echo "$TAKSON takson, $TOPLAM okuma birlestirildi"

# Birlestirmenin kayipsiz oldugu ayrica dogrulanir. Sessiz kayip, bu projede
# "olculmedi"nin "temiz" diye okunmasinin tipik yoludur.
BEKLENEN=$(awk -F'\t' '{s+=$3} END{print s+0}' kaynak_sayim.tsv)
KAYNAKTA=$(awk -F'\t' '{s+=$2} END{print s+0}' kaynak_sayim.tsv)
[ "$TOPLAM" = "$BEKLENEN" ] || { echo "HATA: okuma sayisi tutmuyor ($TOPLAM / $BEKLENEN)"; exit 1; }
echo "okuma sayisi dogrulandi: $TOPLAM (kaynakta toplam $KAYNAKTA)"
# Hangi taksonlarda ornekleme yapildigi acikca yazilir. "Butun okumalar tarandi"
# ile "ornek tarandi" ayni cumleyle bitmemeli.
awk -F'\t' '$2!=$3 {print "    " $1 ": " $2 " okumadan " $3 " alindi"}' kaynak_sayim.tsv | sort -u | head -50

# --- 2. kraken2, tek kosu -----------------------------------------------------
echo
echo "kraken2 calisiyor. Veritabani bir kez yuklenecek, ilk bekleme odur."
[ -n "$BELLEK_BAYRAGI" ] && echo "  --memory-mapping etkin, disk hizina bagli olarak uzun surebilir"
kraken2 --db "$DB" --threads "$IPLIK" $BELLEK_BAYRAGI \
        --output tum.out --report tum.report \
        --use-names tum.fastq
echo "kraken2 bitti"
[ -s tum.out ] || { echo "HATA: tum.out bos, kraken2 sonuc uretmedi"; exit 1; }

# --- 3. bolme ve ozet, ikisi de python tarafinda ------------------------------
# Bolme awk yerine python'da yapilir: awk'in ayni anda acabilecegi dosya sayisi
# bazi surumlerde sinirlidir ve sinir asilinca hata vermeden satir kaybedilir.
echo
python3 "$BURASI/ozgun Kraken ozet betigi" --is "$IS" --ali "$ALI"

echo
echo "bitti. Dosyalar: $IS"
echo "  tum.report        butun numunenin kraken raporu"
echo "  <taxid>.out       her taksonun okuma okuma siniflandirmasi"
echo "  kraken_ozet.csv   takson basina kimlikler ve etiketle uyusma durumu"
