#!/usr/bin/env bash
# =============================================================================
#  INDEKS_KUR.sh  --  Bir referans FASTA icin MFEprimer indeksini SIFIRDAN kurar
#  (2026-08-09: SSU'ya sabitlenmis surumden GENELLESTIRILDI. Eski kullanim
#   birebir korundu: argumansiz calistirilirsa yine SSU'yu indeksler.)
# =============================================================================
#
#  KULLANIM
#  --------
#      bash INDEKS_KUR.sh                                  # SSU (eski davranis)
#      bash INDEKS_KUR.sh SILVA_138.2_LSURef_NR99.fasta    # REFERANS_DB icinden ad
#      bash INDEKS_KUR.sh /tam/yol/baska.fasta             # tam yol da olur
#      bash INDEKS_KUR.sh --liste                          # aday dosyalari goster
#      CPU=8 MEMPCT=60 bash INDEKS_KUR.sh <dosya>          # kaynak ayari
#      SINA_F=... SINA_R=... SINA_AD=... bash INDEKS_KUR.sh <dosya>   # kendi sinama cifti
#
#  NEDEN GEREKIYOR (olculmus sebep):
#    SILVA FASTA dosyalari dizileri RNA olarak saklar: T yerine U.
#    MFEprimer 9'lu k-mer indeksini {A,C,G,T} alfabesi uzerinde kurar. Dosyada
#    hic T olmadigi icin gecerli alfabe {A,C,G}'ye duser ve indekse yalnizca
#    3^9 = 19.683 k-mer girer; olmasi gereken 4^9 = 262.144.
#    Kanit, eski indeksleme gunlugunun kendi satiri:
#        "Sorting 19683 kmers by ID..."      <- SILVA (bozuk)
#        "Sorting 262144 kmers by ID..."     <- PR2 / ROD / UNITE / SSU (saglikli)
#    MFEprimer bu durumda HATA VERMEZ, "Index built successfully" yazar.
#    Sessiz basarisizlik budur.
#
#  COZUM: dizi satirlarindaki U -> T. Baslik satirlarina DOKUNULMAZ.
#
#  ONEMLI TUZAK: kor "sed 's/U/T/g'" KULLANMAYIN. SILVA basliklarinda buyuk U
#    vardir (olculdu: ilk 20 MB'lik dilimdeki 1.652 baslikta). Kor donusum
#    "Unknown Family" -> "Tnknown Family" yapar, taksonomiyi bozar.
#    Bu betik donusumu yalniz dizi satirlarina uygular: sed '/^>/!y/U/T/'
#
#  GUVENLIK: U->T donusumu bayt-bayt birebirdir; dosya boyu ve butun bayt
#    ofsetleri degismez. Bu yuzden mevcut .fai ve BLAST indeksleri gecerli
#    kalir. Betik yine de donusumden ONCE bir yedegin varligini SART kosar
#    (ayni boyda kardes kopya ya da .gz); yoksa dosyaya DOKUNMADAN cikar.
#
#  RNA/DNA OLCUMU: dosyanin RNA mi DNA mi oldugu VARSAYILMAZ, olculur. DNA ise
#    donusum adimi atlanir ve dosya bosuna yeniden yazilmaz; atlandigi ekrana
#    ACIKCA yazilir.
#
#  KESINTIYE DAYANIKLILIK: her adim once "zaten yapilmis mi" diye bakar.
#    Betik yarida kalirsa ayni komutu tekrar calistirin; tamamlanmis adimlari
#    atlar, kaldigi yerden devam eder. Yarim kalan gecici donusum dosyasi
#    (.donusum_tmp) her koşunun basinda silinir, asla kullanilmaz.
# =============================================================================

set -u
set -o pipefail

# ---- Yollar: betigin kendi konumundan turetilir, klasor adi degisse de calisir
KOK="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DB="$KOK/REFERANS_DB"
MFE="$KOK/ARACLAR/mfeprimer"

VARSAYILAN="SILVA_138.2_SSURef_NR99.fasta"   # argumansiz cagri = eski davranis

# ---- --liste / --yardim
if [ "${1:-}" = "--liste" ] || [ "${1:-}" = "-l" ]; then
  echo "REFERANS_DB icindeki indekslenebilir dosyalar:"
  printf '%-42s %10s  %s\n' "DOSYA" "MB" "INDEKS"
  for f in "$DB"/*.fasta "$DB"/*.fna; do
    [ -f "$f" ] || continue
    printf '%-42s %10s  %s\n' "$(basename "$f")" \
      "$(( $(stat -c%s "$f") / 1048576 ))" \
      "$([ -f "$f.primerqc.bin" ] && echo VAR || echo YOK)"
  done
  exit 0
fi
if [ "${1:-}" = "--yardim" ] || [ "${1:-}" = "-h" ]; then
  sed -n '2,40p' "${BASH_SOURCE[0]}"; exit 0
fi

# ---- Hedef FASTA'yi coz: tam yol / REFERANS_DB icinde ad / argumansiz varsayilan
GIRDI="${1:-$VARSAYILAN}"
if [ -f "$GIRDI" ]; then
  FASTA="$(cd "$(dirname "$GIRDI")" && pwd)/$(basename "$GIRDI")"
elif [ -f "$DB/$GIRDI" ]; then
  FASTA="$DB/$GIRDI"
else
  echo "HATA: FASTA bulunamadi: $GIRDI"
  echo "Adaylari gormek icin: bash $(basename "${BASH_SOURCE[0]}") --liste"
  exit 1
fi
AD="$(basename "$FASTA")"
BIN="$FASTA.primerqc.bin"
LOG="$FASTA.log"

# ---- Kaynak ayarlari. Olculmus degerler (Burak'in makinesindeki eski kosular):
#      PR2 351 MB -> 9,34 GB tepe RAM;  UNITE 1454 MB -> 10,87 GB tepe RAM.
CPU="${CPU:-16}"       # cekirdek sayisi   (degistirmek icin: CPU=8 bash INDEKS_KUR.sh)
MEMPCT="${MEMPCT:-70}" # RAM tavani yuzdesi; 70 eski basarili kosularin degeri
UYKU="${UYKU:-60}"     # ilerleme yazdirma araligi (saniye)

zaman() { date '+%H:%M:%S'; }
yaz()   { echo "[$(zaman)] $*"; }
bolum() { echo; echo "============================================================"; echo "  $*"; echo "============================================================"; }

# ---- ISLEVSEL SINAMA CIFTI: veritabaninin bolgesine gore secilir.
#      LSU/28S dosyalarinda SSU primeri hicbir sey bulmaz; o yuzden panelden
#      28S hedefleyen bir cift kullanilir. (Olculdu, TUM_CIFTLER_DEVIR_2026-08-07:
#      Mantar_universal (F1) fungi.28SrRNA.fna'da 226 amplikon, 82 bp -
#      beklenen urun boyuyla birebir. 18S'te 0. Yani gercekten LSU primeri.)
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
  SINA_AD="Bakteri_universal"            # SSU/16S/18S icin eski sinama cifti
  SINA_F="ACAAGCGGTGGAGCATGTG"
  SINA_R="ACGACAGCCATGCAGCAC"
fi

# Mevcut indeks saglikli mi? Iki bagimsiz olcut, biri yeterli:
#   1) gunlukte "Sorting 262144 kmers" (tam DNA alfabesi; yalniz disk modunda yazilir)
#   2) indeks/fasta boyut orani >= 4.0 (saglikli: 4,6-5,5x  /  bozuk: 0,63x)
indeks_saglikli_mi() {
  [ -f "$BIN" ] || return 1
  grep -q "Sorting 262144 kmers" "$LOG" 2>/dev/null && return 0
  local fs bs
  fs=$(stat -c%s "$FASTA"); bs=$(stat -c%s "$BIN")
  awk -v a="$bs" -v b="$fs" 'BEGIN{exit !(b>0 && a/b >= 4.0)}'
}

# =============================================================================
bolum "ADIM 0/4  --  On kontroller"
# =============================================================================
[ -x "$MFE" ] || { echo "HATA: mfeprimer bulunamadi/calistirilamaz: $MFE"; echo "Duzeltmek icin: chmod +x '$MFE'"; exit 1; }
[ -f "$FASTA" ] || { echo "HATA: FASTA yok: $FASTA"; exit 1; }
FASTA_MB=$(( $(stat -c%s "$FASTA") / 1048576 ))
yaz "mfeprimer : $MFE"
yaz "HEDEF     : $AD  (${FASTA_MB} MB)"
yaz "Sinama    : $SINA_AD  ($SINA_F / $SINA_R)"

# --- Sure ve RAM tahmini: OLCUME DAYALI, dosya boyuyla oranlanir.
#     Olculen indeksleme hizlari (bu makine):
#       PR2   351 MB -> 11.508 s = 32,8 s/MB
#       ROD   340 MB -> 12.419 s = 36,5 s/MB
#       UNITE 1454 MB -> 43.056 s = 29,6 s/MB      ortalama ~33 s/MB
#       SSU   786 MB -> 28.020 s = 35,7 s/MB       (SILVA'nin kendi olcumu)
#     Tepe RAM, olculen iki noktadan dogru orantiyla:
#       PR2 351 MB -> 9,34 GB ; UNITE 1454 MB -> 10,87 GB
#       egim = (10,87-9,34)/(1454-351) = 0,001387 GB/MB ; kesme = 8,853 GB
awk -v mb="$FASTA_MB" 'BEGIN{
  printf "           Sure tahmini : %.1f saat (33,0 s/MB ortalama)  |  %.1f saat (35,7 s/MB SSU olcumu)\n", mb*33.0/3600, mb*35.7/3600;
  printf "           Tepe RAM     : ~%.1f GB (0,001387 GB/MB x %d MB + 8,853 GB)\n", mb*0.001387+8.853, mb;
}'

# Indeks, olculmus orana gore fasta boyutunun ~5,1 kati olacak (PR2 5,18x /
# ROD 5,49x / UNITE 4,57x / SSU 5,29x). Donusum sirasinda gecici bir kopya da tutulur.
# Gereken alan = fasta x 5,5 (indeks) + fasta x 1 (gecici kopya) + 500 MB pay.
GEREKEN_MB=$(( FASTA_MB * 7 + 500 ))
BOS_MB=$(df -Pm "$DB" | awk 'NR==2{print $4}')
yaz "Bos disk  : ${BOS_MB} MB (gereken: ~${GEREKEN_MB} MB)"
if [ "$BOS_MB" -lt "$GEREKEN_MB" ]; then
  echo "HATA: yeterli disk alani yok. En az ${GEREKEN_MB} MB bosaltip tekrar deneyin."; exit 1
fi
yaz "Cekirdek  : $CPU     RAM tavani: %$MEMPCT"

# Yarim kalmis gecici dosya varsa temizle (asla kullanilmaz)
if [ -f "$FASTA.donusum_tmp" ]; then
  yaz "UYARI: onceki kosudan kalan yarim donusum dosyasi silindi (.donusum_tmp)."
  rm -f "$FASTA.donusum_tmp"
fi

# =============================================================================
bolum "ADIM 1/4  --  Bozuk indeksi kenara al"
# =============================================================================
# Silmiyoruz, ".bozuk" olarak sakliyoruz ki eski/yeni karsilastirilabilsin.
if [ -f "$BIN" ]; then
  if indeks_saglikli_mi; then
    yaz "Mevcut indeks SAGLIKLI gorunuyor (k-mer sayisi ve/veya boyut orani dogru)."
    yaz "Yeniden kurmak isterseniz once su dosyayi tasiyin: $BIN"
  else
    mv -f "$BIN" "$BIN.bozuk"
    yaz "Bozuk indeks kenara alindi -> $(basename "$BIN").bozuk"
  fi
else
  yaz "Kenara alinacak indeks yok (hic kurulmamis ya da zaten alinmis). Atlaniyor."
fi

# =============================================================================
bolum "ADIM 2/4  --  RNA mi DNA mi? (olculur, varsayilmaz)"
# =============================================================================
# Ilk 50 MB'in DIZI satirlarinda U ve T sayilir. U>0 ise RNA, U=0 ise DNA.
ORNEK=$(head -c 50000000 "$FASTA" | grep -v '^>')
U_SAT=$(printf '%s' "$ORNEK" | grep -c 'U' || true)
T_SAT=$(printf '%s' "$ORNEK" | grep -c 'T' || true)
unset ORNEK
yaz "Olcum (ilk 50 MB, yalniz dizi satirlari): U iceren $U_SAT satir / T iceren $T_SAT satir"

if [ "$U_SAT" -eq 0 ]; then
  yaz ">>> DOSYA ZATEN DNA. U->T DONUSUMU ATLANDI. <<<"
  yaz "    Dosya yeniden yazilmadi, degisiklik listesi bos, gecen sure 0 sn."
  yaz "    (Bu bir hata degil: PR2 / ROD / UNITE / RefSeq kumeleri DNA'dir.)"
else
  yaz "Dosya RNA formatinda. Donusum gerekiyor."

  # --- YEDEK SARTI: orijinali degistirmeden once bir yedek VAR OLMALI.
  YEDEK=""
  BOY=$(stat -c%s "$FASTA")
  for aday in "$DB"/*.fasta "$DB"/*.fna; do
    [ -f "$aday" ] || continue
    [ "$aday" = "$FASTA" ] && continue
    [ "$(stat -c%s "$aday")" = "$BOY" ] && { YEDEK="$aday"; break; }
  done
  [ -z "$YEDEK" ] && [ -f "$FASTA.gz" ] && YEDEK="$FASTA.gz"
  if [ -z "$YEDEK" ]; then
    echo "HATA: bu dosyanin yedegi yok (ne ayni boyda kardes kopya, ne .gz)."
    echo "      Orijinale DOKUNULMADI. Once yedek alin:"
    echo "        cp '$FASTA' '$FASTA.rna_yedek'"
    exit 1
  fi
  yaz "  Yedek bulundu: $(basename "$YEDEK")  (orijinal RNA surumu buradan geri alinabilir)"

  yaz "Donusum basliyor. Sadece dizi satirlari; basliklar korunuyor..."
  TMP="$FASTA.donusum_tmp"
  rm -f "$TMP"
  # y/U/T/ = transliterasyon (s///g'den hizli), /^>/! = baslik satirlarini atla
  sed '/^>/!y/U/T/' "$FASTA" > "$TMP" || { echo "HATA: donusum basarisiz."; rm -f "$TMP"; exit 1; }
  yaz "Donusum bitti, dogrulaniyor..."

  # --- Dogrulama 1: dosya boyu birebir ayni olmali (U->T tek bayt degisimi)
  E=$(stat -c%s "$FASTA"); Y=$(stat -c%s "$TMP")
  if [ "$E" -ne "$Y" ]; then
    echo "HATA: boyut uyusmuyor (eski=$E yeni=$Y). Degisiklik yapilmadi."; rm -f "$TMP"; exit 1
  fi
  yaz "  OK  boyut ayni: $E bayt"

  # --- Dogrulama 2: baslik sayisi degismemis olmali
  BE=$(grep -c '^>' "$FASTA"); BY=$(grep -c '^>' "$TMP")
  if [ "$BE" -ne "$BY" ]; then
    echo "HATA: baslik sayisi uyusmuyor ($BE vs $BY)."; rm -f "$TMP"; exit 1
  fi
  yaz "  OK  baslik sayisi ayni: $BE"

  # --- Dogrulama 3: basliklar bayt-bayt korunmus olmali
  if ! diff -q <(grep '^>' "$FASTA") <(grep '^>' "$TMP") >/dev/null; then
    echo "HATA: basliklar degismis! Degisiklik yapilmadi."; rm -f "$TMP"; exit 1
  fi
  yaz "  OK  basliklar birebir korundu"

  # --- Dogrulama 4: dizi satirlarinda U kalmamis olmali
  KALAN=$(grep -v '^>' "$TMP" | grep -c 'U' || true)
  if [ "$KALAN" -ne 0 ]; then
    echo "HATA: hala U iceren $KALAN dizi satiri var."; rm -f "$TMP"; exit 1
  fi
  yaz "  OK  dizi satirlarinda U kalmadi"

  mv -f "$TMP" "$FASTA"
  yaz "Donusum uygulandi. (Orijinal RNA surumu: $(basename "$YEDEK"))"

  # Ofsetler degismedigi icin .fai gecerli kalir, ama garanti olsun diye
  # mfeprimer'in yeniden uretmesine izin veriyoruz.
  rm -f "$FASTA.fai"
fi

# =============================================================================
bolum "ADIM 3/4  --  Indeksleme  (UZUN ADIM)"
# =============================================================================
if indeks_saglikli_mi; then
  yaz "Saglikli indeks zaten var, indeksleme atlaniyor."
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
    echo "HATA: mfeprimer index cikis kodu $RC ile bitti ($SURE dk sonra)."
    echo "Gunluk sonu:"; tail -20 "$LOG"
    exit 1
  fi
  yaz "Indeksleme bitti: $SURE dakika.  (olculen hiz: $(awk -v s=$SURE -v mb=$FASTA_MB 'BEGIN{if(mb>0) printf "%.1f", s*60/mb; else printf "n/d"}') s/MB)"
fi

# =============================================================================
bolum "ADIM 4/4  --  Dogrulama (uc sinama, ucu de gecmeli)"
# =============================================================================
GECTI_A=0; GECTI_B=0; GECTI_C=0

# --- A: gunlukteki k-mer sayisi 262144 olmali (bozukta 19683 idi)
KMER=$(grep -o "Sorting [0-9]* kmers" "$LOG" 2>/dev/null | tail -1 | grep -o '[0-9]*')
if [ -n "$KMER" ]; then
  if [ "$KMER" = "262144" ]; then
    yaz "A) k-mer sayisi = $KMER  ->  GECTI (4^9, tam DNA alfabesi)"; GECTI_A=1
  else
    yaz "A) k-mer sayisi = $KMER, beklenen 262144  ->  DUSTU (indeks hala bozuk)"
  fi
else
  yaz "A) Gunlukte 'Sorting' satiri yok (bellek modunda kosmus olabilir) -> OLCULEMEDI"
  GECTI_A=1   # olculemeyen sinama "dustu" sayilmaz; B ve C hukmu verir
fi

# --- B: boyut orani. Saglikli indeksler fasta boyunun 4,6-5,5 kati.
if [ -f "$BIN" ]; then
  FS=$(stat -c%s "$FASTA"); BS=$(stat -c%s "$BIN")
  ORAN=$(awk -v a=$BS -v b=$FS 'BEGIN{printf "%.2f", a/b}')
  if awk -v o="$ORAN" 'BEGIN{exit !(o>=4.6 && o<=5.5)}'; then
    yaz "B) indeks/fasta orani = ${ORAN}x  ->  GECTI (saglikli aralik 4,6-5,5x)"; GECTI_B=1
  else
    yaz "B) indeks/fasta orani = ${ORAN}x  ->  DUSTU (beklenen 4,6-5,5x; bozukta 0,63x idi)"
  fi
else
  yaz "B) indeks dosyasi olusmamis: $BIN  ->  DUSTU"
fi

# --- C: ISLEVSEL sinama. Bolgeye uygun primer SIFIRDAN BUYUK baglanma vermeli.
yaz "C) Islevsel sinama: $SINA_AD primeri $AD'ye karsi kosuluyor..."
SINA_DIR="$KOK/DOGRULAMA_SONUC/indeks_sinama"
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
  yaz "C) baglanma = ${VUR:-0}  ->  DUSTU (indeks okunmuyor ya da primer bu bolgede yok)"
fi

# --- Hukum
if [ "$GECTI_A" -eq 1 ] && [ "$GECTI_B" -eq 1 ] && [ "$GECTI_C" -eq 1 ]; then
  bolum "SONUC: BASARILI  --  $AD indeksi saglam (A+B+C gecti)"
  exit 0
else
  bolum "SONUC: BASARISIZ  --  $AD  (A=$GECTI_A B=$GECTI_B C=$GECTI_C; 1=gecti)"
  exit 1
fi
