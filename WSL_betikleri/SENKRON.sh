#!/usr/bin/env bash
# =====================================================================
# SENKRON.sh
# İki taraftaki dosyaların aynı olduğunu SÖYLEMEK yerine DOĞRULAR.
#
#   bash SENKRON.sh              # manifesto üret ve ekrana bas
#   bash SENKRON.sh --dogrula    # kayıtlı manifestoyla karşılaştır
#
# Neden gerekli: "gönderdim" ile "sende var" aynı şey değil. Dosya
# aktarımı sessizce yarım kalabilir, eski sürüm üzerine yazılmayabilir,
# bir adım eski betikle çalışıp doğru görünen yanlış sonuç üretebilir.
# Bu betik her dosyanın sha256 özetini çıkarır; iki taraf aynı özeti
# görüyorsa senkron gerçekten sağlanmıştır.
#
# Manifesto iki bölümdür:
#   BETIK   sürüm kontrolü gereken kod ve tablolar (iki tarafta AYNI olmalı)
#   CIKTI   koşu çıktıları (sizin tarafınızda üretilir, bende kopyası olur)
# =====================================================================
set -uo pipefail

PT="${PT:-/mnt/c/Users/yerli/Masaüstü/PrimerTasarlama}"
HERE="$(cd "$(dirname "$0")" && pwd)"
MANIFEST="$HERE/MANIFEST.sha256"
DOGRULA=0
[ "${1:-}" = "--dogrula" ] && DOGRULA=1

ozet() {   # $1 = kok klasor, gerisi = desenler
  local kok="$1"; shift
  ( cd "$kok" 2>/dev/null || return 0
    for d in "$@"; do
      for f in $d; do
        [ -f "$f" ] && sha256sum "$f"
      done
    done ) 2>/dev/null
}

uret() {
  echo "# MicRhoBooster dosya manifestosu"
  echo "# uretildi: $(date '+%Y-%m-%d %H:%M:%S')"
  echo "# BETIK bolumu iki tarafta BIREBIR AYNI olmali."
  echo "[BETIK]"
  ozet "$HERE" '[0-9][0-9]_*.py' '[0-9][0-9]_*.sh' 'alan_denetimi.py' \
       'hizalama.py' 'CALISTIR.sh' 'AGIR_ISLER.sh' 'SENKRON.sh' \
       'hedefler.tsv' 'hedefler_referans.tsv' 'taxid_adlari.tsv' \
       'bol_*.sh' 'bol2_*.sh' 'kimlik/kimlik_*.tsv' \
       'rapor/rapor_uret.js' 'rapor/veri_uret.py' 'rapor/veri.json' | sort -k2
  echo "[CIKTI]"
  ozet "$PT" 'primer_final/primer_final.tsv' 'primer_final/dis_veritabani.tsv' \
       'primer_final/mfeprimer.tsv' 'primer_final/hedef_kimlik.tsv' \
       'primer_final/teslim_denetimi.tsv' \
       'primer_referans/primer_referans.tsv' \
       'primer_adaylari/ozet.tsv' 'primer_adaylari/ayirt_edilemez.tsv' \
       'primer_adaylari/dislanan_takson.tsv' \
       'kraken_guven/esik_taramasi.tsv' \
       'MicRhoBooster_Primer_Tasarimi.xlsx' \
       'MicRhoBooster_Primer_Raporu.docx' \
       'Microbooster_Topluluk_Trend_Analizi.xlsx' | sort -k2
}

if [ "$DOGRULA" = 1 ]; then
  [ -f "$MANIFEST" ] || { echo "manifesto yok: $MANIFEST"; exit 2; }
  YENI=$(mktemp); uret > "$YENI"
  # yalniz BETIK bolumu karsilastirilir; CIKTI her kosuda degisir
  kes() { awk '/^\[BETIK\]/{f=1;next} /^\[CIKTI\]/{f=0} f' "$1"; }
  ESKI_T=$(mktemp); YENI_T=$(mktemp)
  kes "$MANIFEST" | sort -k2 > "$ESKI_T"
  kes "$YENI"     | sort -k2 > "$YENI_T"
  echo "=== BETIK bolumu karsilastirmasi ==="
  FARK=0
  # yalniz birinde olanlar
  comm -23 <(cut -d' ' -f3- "$ESKI_T" | sort) <(cut -d' ' -f3- "$YENI_T" | sort) \
    | while read -r x; do echo "  KAYIP    $x (manifestoda var, diskte yok)"; FARK=1; done
  comm -13 <(cut -d' ' -f3- "$ESKI_T" | sort) <(cut -d' ' -f3- "$YENI_T" | sort) \
    | while read -r x; do echo "  FAZLA    $x (diskte var, manifestoda yok)"; done
  # ikisinde de olup ozeti farkli olanlar
  join -j2 -o 0,1.1,2.1 "$ESKI_T" "$YENI_T" 2>/dev/null \
    | awk '$2!=$3{printf "  FARKLI   %s\n     manifesto: %s\n     disk     : %s\n",$1,substr($2,1,16),substr($3,1,16); n++} END{if(!n) print "  butun dosyalar manifestoyla ayni"}'
  rm -f "$YENI" "$ESKI_T" "$YENI_T"
  exit 0
fi

uret | tee "$MANIFEST"
echo
echo "yazildi: $MANIFEST"
echo "Bu dosyayi bana gonderin; ben kendi tarafimdaki ozetlerle karsilastirir,"
echo "ayrilan dosya varsa hangisi oldugunu ve hangi yonde guncel oldugunu"
echo "soylerim. Sonraki kontrol icin: bash SENKRON.sh --dogrula"
