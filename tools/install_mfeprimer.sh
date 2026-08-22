#!/bin/bash
# ---------------------------------------------------------------------------
# INSTALLING MFEPRIMER
#
# WHY IT IS NEEDED
# The tools/linux-x64 file is NOT MFEprimer, it is another program called
# "claude-science". The verification stage took it for MFEprimer and even printed
# its version, because the old identity check said "exit code 0 OR the output
# mentions mfeprimer" and every binary returning zero passed it. The only real
# MFEprimer in the directory is mfeprimer.exe, and that is a Windows binary, so
# there is NO MFEprimer that runs under WSL.
#
# MFEprimer is the project's SECOND INDEPENDENT MEASUREMENT. Without it, invariant
# rule number 1 ("if two independent measurements disagree the candidate is
# eliminated") does not work and no order can be placed.
#
# WHAT THIS SCRIPT DOES
# It DOES NOT GUESS the file name. It reads the real asset names from the GitHub
# release API, picks the Linux amd64 one, downloads it, installs it and CONFIRMS
# ITS IDENTITY: if the --version output does not mention "mfeprimer" the install
# counts as failed. It was exactly the absence of that check that made a wrong
# binary be taken for MFEprimer today.
#
# To run:
#   bash install_mfeprimer.sh              # v4.4.0
#   bash install_mfeprimer.sh v4.3.0       # another version
#   LISTE=1 bash install_mfeprimer.sh      # list the assets only, no download
# ---------------------------------------------------------------------------
set -euo pipefail

# THE PROJECT DIRECTORY IS FOUND AUTOMATICALLY, IT IS NOT HARDCODED.
# This script sits inside PROJECT/tools; one directory up is the PROJECT. Because
# BASH_SOURCE is the real path of the running file, this is a measurement rather
# than a guess. If the install was moved it can be overridden from outside:
#   PROJE=/full/path/to/project bash <script>
_BETIK_DIZIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJE="${PROJE:-$(cd "$_BETIK_DIZIN/.." && pwd)}"
# MEASURED: 0_TESLIM_RAPOR IS NOT in this repository; the script used to stop at
# its very first check.
if [ ! -d "$PROJE/tools" ] || [ ! -d "$PROJE/verification" ]; then
  echo "ERROR: the project root could not be verified ('$PROJE' holds no tools/ and verification/)."
  echo "  Script location: $_BETIK_DIZIN"
  echo "  To give it by hand:  PROJE=/full/path/to/project bash $0"
  exit 1
fi
SURUM="${1:-v4.4.0}"
DEPO="quwubin/MFEprimer-3.0"
HEDEF="$PROJE/tools/mfeprimer"
GECICI="$(mktemp -d)"
trap 'rm -rf "$GECICI"' EXIT

command -v curl >/dev/null || { echo "ERROR: curl is missing. sudo apt install curl"; exit 1; }
command -v python3 >/dev/null || { echo "ERROR: python3 is missing"; exit 1; }

echo "depo   : $DEPO"
echo "version: $SURUM"
echo "target : $HEDEF"
echo

API="https://api.github.com/repos/$DEPO/releases/tags/$SURUM"
echo "fetching the release information: $API"
if ! curl -fsSL --max-time 60 -H "Accept: application/vnd.github+json" "$API" -o "$GECICI/surum.json"; then
  echo "ERROR: the release information could not be fetched. Is there internet, and is the tag right?"
  echo "  To see the tags: curl -s https://api.github.com/repos/$DEPO/releases | grep tag_name"
  exit 1
fi

# The assets are listed and the Linux amd64 one is CHOSEN. The name is not
# guessed, it is read from the API.
python3 - "$GECICI/surum.json" > "$GECICI/varliklar.tsv" <<'PY'
import json, sys, re
d = json.load(open(sys.argv[1]))
varlik = d.get("assets", [])
print(f"# version: {d.get('tag_name')}  published: {d.get('published_at')}", file=sys.stderr)
if not varlik:
    print("# this release has no downloadable asset", file=sys.stderr); sys.exit(2)
for a in varlik:
    print(f"{a['name']}\t{a['size']}\t{a['browser_download_url']}")
PY

echo
echo "the assets in this release:"
awk -F'\t' '{printf "  %-46s %8.1f MB\n", $1, $2/1000000}' "$GECICI/varliklar.tsv"
echo

if [ "${LISTE:-0}" = "1" ]; then echo "LISTE=1 was given, nothing was downloaded."; exit 0; fi

# THE SELECTION CRITERION. The architecture is searched for explicitly: amd64,
# x86_64 or x64. A bare "64" IS NOT ENOUGH, because "linux-arm64" holds a 64 too
# and an ARM binary would be installed on an x86 machine. I caught this while
# trying it against a fake asset list; that is exactly what was happening at
# first. arm, aarch, i386 and 386 are explicitly EXCLUDED.
MIMARI="$(uname -m)"
[ "$MIMARI" = "x86_64" ] || echo "WARNING: the machine architecture is $MIMARI, an x86_64 asset is being looked for"
SEC=$(awk -F'\t' 'BEGIN{IGNORECASE=1}
  tolower($1) ~ /linux/ \
  && tolower($1) ~ /amd64|x86_64|x64/ \
  && tolower($1) !~ /arm|aarch|i386|[^0-9]386/ {print $2"\t"$1"\t"$3}' \
  "$GECICI/varliklar.tsv" | sort -n | head -1 || true)

if [ -z "$SEC" ]; then
  echo "ERROR: no Linux amd64 asset was found. Pick the right file from the list above,"
  echo "download it by hand, then do these two steps:"
  echo "  chmod +x <file> && mv <file> \"$HEDEF\""
  echo "  \"$HEDEF\" --version    # the output must mention 'mfeprimer'"
  exit 1
fi

AD=$(printf '%s' "$SEC" | cut -f2)
URL=$(printf '%s' "$SEC" | cut -f3)
echo "the asset chosen: $AD"
echo "downloading    : $URL"
curl -fsSL --max-time 600 "$URL" -o "$GECICI/$AD"
ls -la "$GECICI/$AD"

# UNPACKING. It can be a plain .gz as well as a tar.gz. The first version handled
# only tar.gz and zip; the real asset turned out to be
# "mfeprimer-4.4.0-linux-amd64.gz" and the file was installed STILL COMPRESSED.
# The order matters: .tar.gz has to be tried first, otherwise the .gz branch
# catches it too.
cd "$GECICI"
case "$AD" in
  *.tar.gz|*.tgz) tar xzf "$AD" ;;
  *.tar.bz2|*.tbz) tar xjf "$AD" ;;
  *.tar.xz)       tar xJf "$AD" ;;
  *.tar)          tar xf  "$AD" ;;
  *.zip)          command -v unzip >/dev/null && unzip -q "$AD" || { echo "ERROR: unzip is missing"; exit 1; } ;;
  *.gz)           gunzip -kf "$AD" ;;
  *.bz2)          bunzip2 -kf "$AD" ;;
  *.xz)           unxz -kf "$AD" ;;
esac

# THE IDENTITY CHECK, WITH TWO GATES.
# 1. The file is first checked FROM ITS MAGIC BYTES for being an ELF executable.
#    Compressed or text files are eliminated before ever being run.
# 2. Then --version is run and "mfeprimer" is looked for in the output, but THE
#    FILE PATH is stripped from the output first.
# The reason for the second gate is concrete: in the first version the compressed
# file was tried, the shell said "bash: .../MFEprimer-4.4.0-linux-amd64.gz: cannot
# execute binary file", that error message HELD THE FILE NAME, and it matched
# grep "mfeprimer". So the identity check fooled itself with the file's NAME. I
# had repeated, in the install script, the very fault I was fixing.
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
  echo "ERROR: not one of the downloaded files identified itself as MFEprimer."
  echo "This is the same failure seen before, and it must not pass silently."
  echo "Downloaded: $GECICI"
  find "$GECICI" -type f -size +1M -printf "  %p  %s bytes\n" 2>/dev/null
  echo "File types:"; find "$GECICI" -type f -size +1M -exec file {} \; 2>/dev/null | head
  exit 1
fi

mkdir -p "$PROJE/tools"
cp "$BULUNAN" "$HEDEF"
chmod +x "$HEDEF"
echo
echo "installed: $HEDEF"
echo "identity : $("$HEDEF" --version 2>&1 | head -1)"
echo
echo "verification/mfeprimer_layer.py finds this file on its own, as tools/mfeprimer."
echo "Note: tools/linux-x64 is NOT MFEprimer. It need not be deleted, but it is not used."
echo
echo "Sinama:"
echo "  cd $PROJE/WSL"
echo "  python3 verification/specificity_round.py --kok \"$PROJE\""
echo "Expected: on the 'MFEprimer bulundu' line the version must read MFEprimer, and"
echo "the index build must succeed this time."
