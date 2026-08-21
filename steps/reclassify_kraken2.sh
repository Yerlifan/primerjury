#!/usr/bin/env bash
# =====================================================================
# reclassify_kraken2.sh
# Amaç: ham barkod fastq dosyalarını 106 GB Kraken2 veritabanıyla
#       yeniden sınıflandırmak.
#
# Kullanım:
#   bash reclassify_kraken2.sh \
#        --db   /yol/kraken2_106gb \
#        --in   /yol/ham_fastq_klasoru \
#        --out  /yol/kraken_yeni \
#        [--threads N] [--confidence C] [--force-mmap] [--no-mmap]
#        [--only-benchmark] [--yes]
#
# Tasarım kararları:
#   - Kullanılacak her bayrak, çalıştırmadan önce `kraken2 --help` çıktısında
#     aranır; bulunmazsa betik durur. Hiçbir bayrak hafızadan varsayılmaz.
#   - Bellek eşlemesi kararı RAM ile veritabanı boyutundan türetilir,
#     elle yazılmaz.
#   - Önce en küçük dosyada ölçüm yapılır, toplam süre tahmini basılır,
#     onay alınmadan tam koşu başlamaz (--yes ile atlanır).
#   - Var olan çıktının üzerine yazmaz, atlar; yeniden üretmek için silin.
#   - --confidence verilirse çıktı adlarına eklenir (ör. _c0.1), böylece
#     farklı güven eşikli koşular birbirinin üzerine yazmaz. Kraken2
#     varsayılan ayarında (0) çekimser kalmaz: gerçek organizma
#     veritabanında yoksa okuma en yüksek puanlı kardeş yaprağa düşer.
#     Ölçüldü: A2-4'teki dört Methanosarcina kutusunun okumaları aynı
#     referanslara gidiyor, hiçbiri kendi atandığı türü tercih etmiyor.
#     0.1 dolayında bir eşik bu okumaları cins düğümünde toplar.
# =====================================================================
set -euo pipefail

DB=""; IN=""; OUT=""; THREADS=""; MMAP="auto"; ONLYBENCH=0; ASSUME_YES=0
CONF=""
while [ $# -gt 0 ]; do
  case "$1" in
    --db) DB="$2"; shift 2;;
    --in) IN="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --confidence) CONF="$2"; shift 2;;
    --force-mmap) MMAP="on"; shift;;
    --no-mmap) MMAP="off"; shift;;
    --only-benchmark) ONLYBENCH=1; shift;;
    --yes) ASSUME_YES=1; shift;;
    *) echo "bilinmeyen secenek: $1" >&2; exit 2;;
  esac
done
[ -n "$DB" ] && [ -n "$IN" ] && [ -n "$OUT" ] || {
  echo "kullanim: bash $0 --db <kraken2_db> --in <ham_fastq_klasoru> --out <cikti_klasoru>" >&2; exit 2; }

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die() { printf 'HATA: %s\n' "$*" >&2; exit 1; }

command -v kraken2 >/dev/null 2>&1 || die "kraken2 PATH'te yok"

# --- 1. veritabanı bütünlüğü -----------------------------------------
for f in hash.k2d opts.k2d taxo.k2d; do
  [ -e "$DB/$f" ] || die "veritabaninda $f yok: $DB"
done
DBBYTES=$(du -sb "$DB" | cut -f1)
DBGB=$(awk -v b="$DBBYTES" 'BEGIN{printf "%.1f", b/1073741824}')
log "veritabani: $DB  (${DBGB} GB)"

# --- 2. bayrak doğrulaması -------------------------------------------
HELP=$(kraken2 --help 2>&1 || true)
need_flag() { grep -q -- "$1" <<<"$HELP" || die "bu kraken2 surumu $1 bayragini tanimiyor"; }
for fl in --db --threads --report --output --use-names; do need_flag "$fl"; done
HAS_MMAP=0; grep -q -- "--memory-mapping" <<<"$HELP" && HAS_MMAP=1
# --confidence hafizadan varsayilmaz, yardim ciktisinda aranir
if [ -n "$CONF" ]; then
  grep -q -- "--confidence" <<<"$HELP" || die "bu kraken2 surumu --confidence bayragini tanimiyor"
  case "$CONF" in
    ''|*[!0-9.]*) die "--confidence sayisal olmali: $CONF";;
  esac
  awk -v c="$CONF" 'BEGIN{ if (c+0 < 0 || c+0 > 1) exit 1 }' || die "--confidence 0 ile 1 arasinda olmali: $CONF"
  log "guven esigi: $CONF (cikti adlarina _c$CONF eklenecek)"
fi
HAS_GZ=0;   grep -q -- "--gzip-compressed" <<<"$HELP" && HAS_GZ=1
HAS_MINHIT=0; grep -q -- "--minimum-hit-groups" <<<"$HELP" && HAS_MINHIT=1
HAS_RMZ=0;  grep -q -- "--report-minimizer-data" <<<"$HELP" && HAS_RMZ=1
log "bayrak destegi: memory-mapping=$HAS_MMAP gzip=$HAS_GZ minimum-hit-groups=$HAS_MINHIT report-minimizer-data=$HAS_RMZ"

# --- 3. iş parçacığı ve bellek kararı --------------------------------
CORES=$(nproc)
if [ -z "$THREADS" ]; then THREADS=$(( CORES > 2 ? CORES - 2 : 1 )); fi
[ "$THREADS" -gt "$CORES" ] && log "UYARI: istenen is parcacigi ($THREADS) cekirdek sayisindan ($CORES) fazla"
AVAILGB=$(awk '/MemAvailable/ {printf "%.1f", $2/1048576}' /proc/meminfo)
log "cekirdek=$CORES  kullanilan is parcacigi=$THREADS  kullanilabilir RAM=${AVAILGB} GB"

USE_MMAP=0
case "$MMAP" in
  on) USE_MMAP=1; log "bellek eslemesi: elle acildi";;
  off) USE_MMAP=0; log "bellek eslemesi: elle kapatildi";;
  auto)
    NEED=$(awk -v g="$DBGB" 'BEGIN{printf "%.1f", g/0.8}')
    if awk -v a="$AVAILGB" -v n="$NEED" 'BEGIN{exit !(a>=n)}'; then
      USE_MMAP=0
      log "bellek eslemesi: GEREKMEZ (RAM ${AVAILGB} GB >= gerekli ${NEED} GB)"
    else
      USE_MMAP=1
      log "bellek eslemesi: GEREKIR (RAM ${AVAILGB} GB < gerekli ${NEED} GB), sure diske bagli olacak"
    fi;;
esac
[ "$USE_MMAP" = 1 ] && [ "$HAS_MMAP" = 0 ] && die "bellek eslemesi gerekli ama bu surum desteklemiyor"

# --- 4. girdi dosyaları ----------------------------------------------
mkdir -p "$OUT"
mapfile -t FQ < <(find "$IN" -maxdepth 2 -type f \
  \( -iname '*.fastq' -o -iname '*.fq' -o -iname '*.fastq.gz' -o -iname '*.fq.gz' \) | sort)
[ ${#FQ[@]} -gt 0 ] || die "girdi klasorunde fastq bulunamadi: $IN"
log "girdi dosyasi sayisi: ${#FQ[@]}"

TOTBYTES=0
for f in "${FQ[@]}"; do TOTBYTES=$(( TOTBYTES + $(stat -c%s "$f") )); done
log "toplam girdi: $(awk -v b=$TOTBYTES 'BEGIN{printf "%.2f", b/1073741824}') GB"

# Taban ad cakismasi kontrolu. Ayni tabana inen iki girdi (ornegin barcode03.fastq
# ve barcode03.fastq.gz) ayni cikti adini uretir ve biri sessizce dusurulur.
# Bu projede F2-1 klasorundeki kopyalar tam bu tur bir sessiz kayba yol acmisti,
# bu yuzden burada durup listeliyoruz.
declare -A SEEN=()
COLL=0
for f in "${FQ[@]}"; do
  b=$(basename "$f"); b=${b%.gz}; b=${b%.fastq}; b=${b%.fq}
  if [ -n "${SEEN[$b]:-}" ]; then
    printf 'CLASH: "%s" and "%s" both land on the same output name (%s)\n' "${SEEN[$b]}" "$f" "$b" >&2
    COLL=1
  else
    SEEN[$b]="$f"
  fi
done
[ "$COLL" = 1 ] && die "taban ad cakismasi var, sessiz kaybi onlemek icin duruldu. Fazla kopyalari girdi klasorunden cikarin."
log "taban ad cakismasi yok, ${#FQ[@]} girdinin tamami ayri cikti uretecek"

run_one() {  # $1 = fastq yolu
  local fq="$1" base gz=() extra=()
  base=$(basename "$fq"); base=${base%.gz}; base=${base%.fastq}; base=${base%.fq}
  # Guven esigi cikti adina girer; aksi halde 0 ile 0.1 kosulari ayni
  # dosyaya yazar ve hangi sonucun hangi esikten geldigi kaybolur.
  local ek=""; [ -n "$CONF" ] && ek="_c$CONF"
  local rep="$OUT/${base}${ek}_kraken2.report" outp="$OUT/${base}${ek}_output"
  if [ -s "$rep" ] && [ -s "$outp" ]; then echo "ATLANDI $base"; return 0; fi
  case "$fq" in *.gz) [ "$HAS_GZ" = 1 ] && gz=(--gzip-compressed);; esac
  [ "$USE_MMAP" = 1 ] && extra+=(--memory-mapping)
  [ "$HAS_MINHIT" = 1 ] && extra+=(--minimum-hit-groups 3)
  [ -n "$CONF" ] && extra+=(--confidence "$CONF")
  kraken2 --db "$DB" --threads "$THREADS" "${gz[@]}" "${extra[@]}" \
          --use-names --report "$rep" --output "$outp" "$fq"
}

# --- 5. ölçüm koşusu (en küçük dosya) --------------------------------
SMALL=$(for f in "${FQ[@]}"; do printf '%s\t%s\n' "$(stat -c%s "$f")" "$f"; done | sort -n | head -1 | cut -f2-)
SMALLB=$(stat -c%s "$SMALL")
log "olcum dosyasi: $(basename "$SMALL") ($(awk -v b=$SMALLB 'BEGIN{printf "%.1f", b/1048576}') MB)"
T0=$(date +%s)
run_one "$SMALL"
T1=$(date +%s)
ELAPSED=$(( T1 - T0 )); [ "$ELAPSED" -lt 1 ] && ELAPSED=1
log "olcum suresi: ${ELAPSED} saniye"
awk -v e="$ELAPSED" -v s="$SMALLB" -v t="$TOTBYTES" 'BEGIN{
  r=e/s; tot=r*t;
  printf "[estimate] rough total time: %.1f minutes (%.2f hours)\n", tot/60, tot/3600;
  printf "[estimate] note: the first call includes loading the database,\n";
  printf "           and later files are faster when memory mapping is off.\n"}'

[ "$ONLYBENCH" = 1 ] && { log "yalnizca olcum istendi, duruluyor"; exit 0; }
if [ "$ASSUME_YES" != 1 ]; then
  read -r -p "Tam kosuyu baslatalim mi? [e/H] " a
  case "$a" in e|E|y|Y) ;; *) log "iptal edildi"; exit 0;; esac
fi

# --- 6. tam koşu -----------------------------------------------------
i=0
for f in "${FQ[@]}"; do
  i=$(( i + 1 ))
  log "[$i/${#FQ[@]}] $(basename "$f")"
  run_one "$f"
done

# --- 7. özet ---------------------------------------------------------
log "ozet"
python3 - "$OUT" <<'PY'
import sys,glob,os
out=sys.argv[1]
print("%-34s %10s %9s %7s"%("rapor","toplam","unclass","unc%"))
for f in sorted(glob.glob(os.path.join(out,"*_kraken2.report"))):
    unc=root=0
    for line in open(f,errors="replace"):
        p=line.rstrip("\n").split("\t")
        if len(p)<6: continue
        if p[3]=="U": unc=int(p[1])
        if p[3]=="R" and p[4]=="1": root=int(p[1])
    tot=unc+root
    if tot: print("%-34s %10d %9d %6.2f"%(os.path.basename(f)[:34],tot,unc,100*unc/tot))
PY
log "bitti. cikti: $OUT"
echo
echo "The next step is the Bracken re-estimation, but first share the reports in"
echo "this directory and ortam_raporu.txt; the kmer_distrib file that suits"
echo "kmer_distrib uzunlugunu veriden secip betigi kesinlestirecegim."
