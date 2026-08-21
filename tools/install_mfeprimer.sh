#!/bin/bash
# ---------------------------------------------------------------------------
# MFEPRIMER KURULUMU
#
# NEDEN GEREKLI
# tools/linux-x64 dosyasi MFEprimer DEGIL, "claude-science" adli baska bir program.
# dogrulama asamasi onu MFEprimer sanip surumunu bile ekrana yaziyordu, cunku eski
# kimlik kontrolu "cikis kodu 0 VEYA ciktida mfeprimer gecsin" diyordu ve sifir donen
# her ikili bu testi geciyordu. Klasordeki tek gercek MFEprimer mfeprimer.exe, o da
# Windows ikilisi. Yani WSL'de calisan bir MFEprimer YOK.
#
# MFEprimer, projenin IKINCI BAGIMSIZ OLCUMU. O olmadan 1 numarali degismez kural
# ("iki bagimsiz olcum ayrilirsa aday elenir") islemez ve siparis verilemez.
#
# BU BETIK NE YAPAR
# Dosya adini TAHMIN ETMEZ. GitHub surum API'sinden varliklarin gercek adlarini okur,
# Linux amd64 olani secer, indirir, kurar ve KIMLIGINI DOGRULAR: --version ciktisinda
# "mfeprimer" gecmiyorsa kurulumu basarisiz sayar. Bugun tam da bu dogrulamanin
# eksikligi yuzunden yanlis bir ikili MFEprimer sanildi.
#
# Calistirma:
#   bash install_mfeprimer.sh              # v4.4.0
#   bash install_mfeprimer.sh v4.3.0       # baska surum
#   LISTE=1 bash install_mfeprimer.sh      # yalnizca varliklari listele, indirme
# ---------------------------------------------------------------------------
set -euo pipefail

# PROJE KLASORU OTOMATIK BULUNUR, SABIT YAZILMAZ.
# Bu betik PROJE/WSL icinde duruyor; bir ust klasor PROJE'dir. BASH_SOURCE
# calisan dosyanin gercek yolu oldugu icin bu tahmin degil olcumdur.
# Kurulum tasindiysa disaridan bastirilabilir:  PROJE=/tam/yol/proje bash <betik>
_BETIK_DIZIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJE="${PROJE:-$(cd "$_BETIK_DIZIN/.." && pwd)}"
if [ ! -d "$PROJE/0_TESLIM_RAPOR" ]; then
  echo "HATA: proje klasoru dogrulanamadi ('$PROJE' icinde 0_TESLIM_RAPOR yok)."
  echo "  Betik konumu: $_BETIK_DIZIN"
  echo "  Elle vermek icin:  PROJE=/tam/yol/proje bash $0"
  exit 1
fi
SURUM="${1:-v4.4.0}"
DEPO="quwubin/MFEprimer-3.0"
HEDEF="$PROJE/tools/mfeprimer"
GECICI="$(mktemp -d)"
trap 'rm -rf "$GECICI"' EXIT

command -v curl >/dev/null || { echo "HATA: curl yok. sudo apt install curl"; exit 1; }
command -v python3 >/dev/null || { echo "HATA: python3 yok"; exit 1; }

echo "depo   : $DEPO"
echo "surum  : $SURUM"
echo "hedef  : $HEDEF"
echo

API="https://api.github.com/repos/$DEPO/releases/tags/$SURUM"
echo "surum bilgisi aliniyor: $API"
if ! curl -fsSL --max-time 60 -H "Accept: application/vnd.github+json" "$API" -o "$GECICI/surum.json"; then
  echo "HATA: surum bilgisi alinamadi. Internet var mi, surum etiketi dogru mu?"
  echo "  Etiketleri gormek icin: curl -s https://api.github.com/repos/$DEPO/releases | grep tag_name"
  exit 1
fi

# Varliklar listelenir ve Linux amd64 olan SECILIR. Ad tahmin edilmez, API'den okunur.
python3 - "$GECICI/surum.json" > "$GECICI/varliklar.tsv" <<'PY'
import json, sys, re
d = json.load(open(sys.argv[1]))
varlik = d.get("assets", [])
print(f"# surum: {d.get('tag_name')}  yayin: {d.get('published_at')}", file=sys.stderr)
if not varlik:
    print("# bu surumde indirilebilir varlik yok", file=sys.stderr); sys.exit(2)
for a in varlik:
    print(f"{a['name']}\t{a['size']}\t{a['browser_download_url']}")
PY

echo
echo "surumdeki varliklar:"
awk -F'\t' '{printf "  %-46s %8.1f MB\n", $1, $2/1000000}' "$GECICI/varliklar.tsv"
echo

if [ "${LISTE:-0}" = "1" ]; then echo "LISTE=1 verildi, indirme yapilmadi."; exit 0; fi

# SECIM OLCUTU. Mimari acikca aranir: amd64, x86_64 ya da x64. Ciplak "64" YETMEZ,
# cunku "linux-arm64" de icinde 64 gecirir ve x86 makineye ARM ikilisi kurulurdu.
# Bunu sahte bir varlik listesiyle denerken yakaladim, once tam olarak bu oluyordu.
# arm, aarch, i386 ve 386 acikca DISLANIR.
MIMARI="$(uname -m)"
[ "$MIMARI" = "x86_64" ] || echo "UYARI: makine mimarisi $MIMARI, x86_64 varligi araniyor"
SEC=$(awk -F'\t' 'BEGIN{IGNORECASE=1}
  tolower($1) ~ /linux/ \
  && tolower($1) ~ /amd64|x86_64|x64/ \
  && tolower($1) !~ /arm|aarch|i386|[^0-9]386/ {print $2"\t"$1"\t"$3}' \
  "$GECICI/varliklar.tsv" | sort -n | head -1 || true)

if [ -z "$SEC" ]; then
  echo "HATA: Linux amd64 varligi bulunamadi. Yukaridaki listeden dogru dosyayi secip"
  echo "elle indirin, sonra su iki adimi yapin:"
  echo "  chmod +x <dosya> && mv <dosya> \"$HEDEF\""
  echo "  \"$HEDEF\" --version    # ciktida 'mfeprimer' gecmeli"
  exit 1
fi

AD=$(printf '%s' "$SEC" | cut -f2)
URL=$(printf '%s' "$SEC" | cut -f3)
echo "secilen varlik: $AD"
echo "indiriliyor   : $URL"
curl -fsSL --max-time 600 "$URL" -o "$GECICI/$AD"
ls -la "$GECICI/$AD"

# ACMA. Duz .gz de olabilir, tar.gz de. Ilk surumde yalnizca tar.gz ve zip vardi;
# gercek varlik "mfeprimer-4.4.0-linux-amd64.gz" cikti ve dosya SIKISTIRILMIS halde
# kuruldu. Sira onemli: .tar.gz once denenmeli, yoksa .gz dali onu da yakalar.
cd "$GECICI"
case "$AD" in
  *.tar.gz|*.tgz) tar xzf "$AD" ;;
  *.tar.bz2|*.tbz) tar xjf "$AD" ;;
  *.tar.xz)       tar xJf "$AD" ;;
  *.tar)          tar xf  "$AD" ;;
  *.zip)          command -v unzip >/dev/null && unzip -q "$AD" || { echo "HATA: unzip yok"; exit 1; } ;;
  *.gz)           gunzip -kf "$AD" ;;
  *.bz2)          bunzip2 -kf "$AD" ;;
  *.xz)           unxz -kf "$AD" ;;
esac

# KIMLIK ARAMASI, IKI KAPILI.
# 1. Dosya once ELF calistirilabilir mi diye MAGIC BAYTLARINDAN bakilir. Sikistirilmis
#    ya da metin dosyalari daha calistirilmadan elenir.
# 2. Sonra --version calistirilir ve ciktida "mfeprimer" aranir, ama once DOSYA YOLU
#    ciktidan silinir.
# Ikinci kapinin sebebi somut: ilk surumde sikistirilmis dosya calistirilmaya
# calisildi, kabuk "bash: .../MFEprimer-4.4.0-linux-amd64.gz: cannot execute binary
# file" dedi, bu hata mesaji DOSYA ADINI iceriyordu ve grep "mfeprimer" ile esletti.
# Yani kimlik kontrolu, dosyanin ADIYLA kendini kandirdi. Duzelttigim hatanin
# aynisini kurulum betiginde tekrarlamisim.
elf_mi() {
  [ -f "$1" ] || return 1
  head -c 4 "$1" 2>/dev/null | od -An -tx1 2>/dev/null | tr -d ' \n' | grep -qi '^7f454c46'
}
BULUNAN=""
while IFS= read -r f; do
  elf_mi "$f" || continue
  chmod +x "$f" 2>/dev/null || true
  out=$("$f" --version 2>&1 || true)
  # yol ciktidan cikarilir ki ad uzerinden yanlis eslesme olmasin
  temiz=$(printf '%s' "$out" | sed "s|$(printf '%s' "$f" | sed 's/[]\/$*.^[]/\\&/g')||g")
  if printf '%s' "$temiz" | grep -qi mfeprimer; then BULUNAN="$f"; break; fi
done < <(find "$GECICI" -type f -size +1M 2>/dev/null)

if [ -z "$BULUNAN" ]; then
  echo
  echo "HATA: indirilen dosyalarin hicbiri kimligini MFEprimer olarak bildirmedi."
  echo "Bu, bugun yasanan hatanin aynisidir ve sessiz gecmemelidir."
  echo "Indirilenler: $GECICI"
  find "$GECICI" -type f -size +1M -printf "  %p  %s bayt\n" 2>/dev/null
  echo "Dosya turleri:"; find "$GECICI" -type f -size +1M -exec file {} \; 2>/dev/null | head
  exit 1
fi

mkdir -p "$PROJE/tools"
cp "$BULUNAN" "$HEDEF"
chmod +x "$HEDEF"
echo
echo "kuruldu: $HEDEF"
echo "kimlik : $("$HEDEF" --version 2>&1 | head -1)"
echo
echo "verification/mfeprimer_layer.py bu dosyayi tools/mfeprimer adiyla kendisi buluyor."
echo "Not: tools/linux-x64 MFEprimer DEGIL, silinmesi gerekmez ama kullanilmaz."
echo
echo "Sinama:"
echo "  cd $PROJE/WSL"
echo "  python3 verification/specificity_round.py --kok \"$PROJE\""
echo "Beklenen: 'MFEprimer bulundu' satirinda surum MFEprimer olmali ve"
echo "indeks kurulumu bu sefer basarili olmali."
