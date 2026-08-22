#!/usr/bin/env bash
# =====================================================================
# sync.sh
# It VERIFIES that the files on the two sides are the same, rather than SAYING so.
#
#   bash sync.sh              # produce the manifest and print it
#   bash sync.sh --dogrula    # compare against the recorded manifest
#
# Why it is needed: "I sent it" and "you have it" are not the same thing. A file
# transfer can stop half way in silence, an old version may not get overwritten, and
# a step can run with an old script and produce a wrong result that looks right.
# This script takes the sha256 digest of every file; if the two sides see the same
# digest then they really are in sync.
#
# The manifest has two sections:
#   BETIK   the code and tables that need version control (they must be THE SAME on
#           both sides)
#   CIKTI   the run outputs (produced on your side, a copy of them on mine)
# =====================================================================
set -uo pipefail

# The project directory. It used to default to one person's desktop path, so on
# any other machine the defaults below pointed at nothing. It is derived from
# this script's own location now; $PT still overrides it.
_BETIK_DIZIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PT="${PT:-$(cd "$_BETIK_DIZIN/.." && pwd)}"
HERE="$(cd "$(dirname "$0")" && pwd)"
MANIFEST="$HERE/MANIFEST.sha256"
DOGRULA=0
[ "${1:-}" = "--dogrula" ] && DOGRULA=1

ozet() {   # $1 = the root directory, the rest are patterns
  local kok="$1"; shift
  ( cd "$kok" 2>/dev/null || return 0
    for d in "$@"; do
      for f in $d; do
        [ -f "$f" ] && sha256sum "$f"
      done
    done ) 2>/dev/null
}

uret() {
  echo "# PrimerJury file manifest"
  echo "# uretildi: $(date '+%Y-%m-%d %H:%M:%S')"
  echo "# The SCRIPT section must be EXACTLY THE SAME on both sides."
  echo "[BETIK]"
  ozet "$HERE" '[0-9][0-9]_*.py' '[0-9][0-9]_*.sh' 'field_audit.py' \
       'alignment.py' 'run.sh' 'heavy_jobs.sh' 'sync.sh' \
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
       'PrimerJury_Primer_Tasarimi.xlsx' \
       'PrimerJury_Primer_Raporu.docx' \
       'PrimerJury_Community_Trends.xlsx' | sort -k2
}

if [ "$DOGRULA" = 1 ]; then
  [ -f "$MANIFEST" ] || { echo "no manifest: $MANIFEST"; exit 2; }
  YENI=$(mktemp); uret > "$YENI"
  # only the BETIK section is compared; CIKTI changes on every run
  kes() { awk '/^\[BETIK\]/{f=1;next} /^\[CIKTI\]/{f=0} f' "$1"; }
  ESKI_T=$(mktemp); YENI_T=$(mktemp)
  kes "$MANIFEST" | sort -k2 > "$ESKI_T"
  kes "$YENI"     | sort -k2 > "$YENI_T"
  echo "=== comparing the [BETIK] section ==="
  FARK=0
  # yalniz birinde olanlar
  comm -23 <(cut -d' ' -f3- "$ESKI_T" | sort) <(cut -d' ' -f3- "$YENI_T" | sort) \
    | while read -r x; do echo "  MISSING  $x (in the manifest, not on disk)"; FARK=1; done
  comm -13 <(cut -d' ' -f3- "$ESKI_T" | sort) <(cut -d' ' -f3- "$YENI_T" | sort) \
    | while read -r x; do echo "  EXTRA    $x (on disk, not in the manifest)"; done
  # the ones present in both whose digest differs
  join -j2 -o 0,1.1,2.1 "$ESKI_T" "$YENI_T" 2>/dev/null \
    | awk '$2!=$3{printf "  DIFFERENT %s\n     manifest: %s\n     disk    : %s\n",$1,substr($2,1,16),substr($3,1,16); n++} END{if(!n) print "  butun dosyalar manifestoyla ayni"}'
  rm -f "$YENI" "$ESKI_T" "$YENI_T"
  exit 0
fi

uret | tee "$MANIFEST"
echo
echo "yazildi: $MANIFEST"
echo "Send me this file; I will compare it against the digests on my side,"
echo "If a file has diverged, this says which one and which side is newer."
echo "For the next check: bash sync.sh --dogrula"
