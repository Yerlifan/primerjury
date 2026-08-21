#!/usr/bin/env bash
# =====================================================================
# check_environment.sh
# Amaç: yeniden sınıflandırma öncesi ortamı ölçmek. Hiçbir şeyi
# değiştirmez, yalnızca okur ve rapor eder.
# Kullanım (WSL içinde):
#     bash check_environment.sh 2>&1 | tee ortam_raporu.txt
# Çıktı: ortam_raporu.txt ve yardim_ciktilari/ klasörü.
# =====================================================================
set -uo pipefail
OUT="yardim_ciktilari"; mkdir -p "$OUT"
bar() { printf '\n%s\n%s\n' "=== $* ===" "-----------------------------------------------"; }

bar "1. Makine"
echo "tarih          : $(date -Iseconds)"
echo "cekirdek       : $(nproc)"
if [ -r /proc/meminfo ]; then
  awk '/MemTotal|MemAvailable|SwapTotal/ {printf "%-15s: %.1f GB\n", $1, $2/1048576}' /proc/meminfo
fi
echo "kernel         : $(uname -r)"
grep -qi microsoft /proc/version && echo "environment    : WSL" || echo "environment    : not WSL"

bar "2. Araç sürümleri"
for t in kraken2 kraken2-build bracken bracken-build est_abundance.py samtools minimap2 seqkit mash python3 blastn makeblastdb primer3_core; do
  p=$(command -v "$t" 2>/dev/null)
  if [ -n "$p" ]; then
    v=$("$t" --version 2>&1 | head -2 | tr '\n' ' ' | cut -c1-70)
    [ -z "${v// }" ] && v=$("$t" -version 2>&1 | head -1 | cut -c1-70)
    printf "%-18s VAR   %-46s %s\n" "$t" "$v" "$p"
  else
    printf "%-18s YOK\n" "$t"
  fi
done
echo
echo "Biopython:"; python3 -c "import Bio,sys;print('  Bio',Bio.__version__)" 2>&1 | head -1
echo "primer3-py:"; python3 -c "import primer3;print('  primer3-py',primer3.__version__)" 2>&1 | head -1

bar "3. Yardım çıktıları dosyaya yazılıyor (bayrakları hafızadan varsaymamak için)"
for pair in "kraken2 --help" "bracken --help" "samtools --version" "minimap2 --help" "seqkit version" "blastn -help"; do
  set -- $pair; tool=$1; shift
  if command -v "$tool" >/dev/null 2>&1; then
    f="$OUT/${tool}_help.txt"
    { echo "### $tool $*"; "$tool" "$@"; } > "$f" 2>&1
    echo "written: $f ($(wc -l < "$f") lines)"
  fi
done
if command -v samtools >/dev/null 2>&1; then
  { echo "### samtools consensus --help"; samtools consensus --help; } > "$OUT/samtools_consensus_help.txt" 2>&1
  echo "yazildi: $OUT/samtools_consensus_help.txt"
fi
if command -v est_abundance.py >/dev/null 2>&1; then
  { echo "### est_abundance.py -h"; est_abundance.py -h; } > "$OUT/est_abundance_help.txt" 2>&1
  echo "yazildi: $OUT/est_abundance_help.txt"
fi

bar "4. Kraken2 veritabanları (hash.k2d dosyasına göre otomatik aranır)"
# Aramayı makul köklerle sınırlıyoruz, tüm diski taramıyoruz.
ROOTS=("$HOME" /mnt/c /mnt/d /mnt/e /mnt/f /opt /srv /data)
FOUND=()
for r in "${ROOTS[@]}"; do
  [ -d "$r" ] || continue
  while IFS= read -r h; do FOUND+=("$(dirname "$h")"); done < <(
    timeout 240 find "$r" -maxdepth 7 -name hash.k2d -not -path '*/\.*' 2>/dev/null)
done
if [ ${#FOUND[@]} -eq 0 ]; then
  echo "hash.k2d was not found. Give the database path by hand."
else
  printf '%s\n' "${FOUND[@]}" | sort -u | while read -r db; do
    echo
    echo "DB: $db"
    tot=$(du -sb "$db" 2>/dev/null | cut -f1)
    printf "  total size   : %.1f GB\n" "$(echo "$tot" | awk '{print $1/1073741824}')"
    for f in hash.k2d opts.k2d taxo.k2d seqid2taxid.map; do
      if [ -e "$db/$f" ]; then
        printf "  %-16s %10.2f GB\n" "$f" "$(stat -c%s "$db/$f" | awk '{print $1/1073741824}')"
      else
        printf "  %-16s YOK\n" "$f"
      fi
    done
    echo "  Bracken kmer_distrib dosyalari:"
    ls -1 "$db"/database*mers.kmer_distrib 2>/dev/null | sed 's/^/    /' || echo "    YOK (bracken-build gerekir)"
    echo "  inspect.txt / library:"
    ls -d "$db"/inspect.txt "$db"/library "$db"/taxonomy 2>/dev/null | sed 's/^/    /'
  done
fi

bar "5. RAM ile veritabanı boyutu karşılaştırması"
python3 - <<'PY'
import os,subprocess,shutil
mem=0
for l in open('/proc/meminfo'):
    if l.startswith('MemAvailable'): mem=int(l.split()[1])/1048576
print("kullanilabilir RAM: %.1f GB"%mem)
print("Karar kurali: veritabani boyutu kullanilabilir RAM'in yuzde 80'inden")
print("kucukse --memory-mapping GEREKMEZ, buyukse GEREKIR ve sure diske baglidir.")
PY

bar "6. HAM okuma dosyaları (en kritik nokta)"
cat <<'NOT'
PrimerTasarlama\fastq files klasorundeki dosyalar HAM OKUMA DEGIL.
kaynak calismanin sequence_extraction.sh betigi her ornekten yalnizca en bol ILK BES
taksonu cikarmis. Kraken raporlarindan olculen tutulma oranlari:
  A1 grubu (bc01-04) : yuzde 71 ile 85
  A2 grubu (bc05-08) : yuzde 90 ile 98
  F2 grubu (bc09-12) : yuzde 12 ile 63
  F1 grubu (bc13-16) : yuzde 17 ile 41
  B  grubu (bc17-20) : yuzde 12 ile 18
Siniflandirilamayan okuma orani yuzde 0,04 ile 1,68 arasinda, yani kayip
siniflandirmadan degil ilk bes filtresinden geliyor.
Bu yuzden 106 GB veritabaniyla yeniden siniflandirma icin ORIJINAL barkod
fastq dosyalari gerekiyor (kraken2_driver.sh betigindeki ./inputs/*.fastq).
Asagidaki arama onlari bulmayi dener.
NOT
echo
for r in "$HOME" /mnt/c /mnt/d /mnt/e; do
  [ -d "$r" ] || continue
  timeout 180 find "$r" -maxdepth 7 \( -iname 'barcode*.fastq' -o -iname 'barcode*.fastq.gz' \
     -o -iname '*barcode??.fastq*' \) -not -path '*/\.*' 2>/dev/null | head -40
done | sort -u | while read -r f; do printf "  %10s  %s\n" "$(du -h "$f" | cut -f1)" "$f"; done
echo "(no line here means no raw fastq was found)"

bar "7. Okuma uzunlugu (Bracken kmer_distrib secimi icin veriden olculur)"
PT=""
for c in "/mnt/c/Users/yerli/Masaüstü/PrimerTasarlama" "/mnt/c/Users/yerli/Masaustu/PrimerTasarlama"; do
  [ -d "$c" ] && PT="$c" && break
done
if [ -n "$PT" ]; then
  echo "PrimerTasarlama: $PT"
  python3 - "$PT" <<'PY'
import sys,glob,os,statistics
pt=sys.argv[1]
lens=[]
for f in sorted(glob.glob(os.path.join(pt,"kraken results","*","*_output")))[:20]:
    n=0
    for line in open(f,errors="replace"):
        p=line.split("\t")
        if len(p)>=4 and p[3].strip().isdigit():
            lens.append(int(p[3])); n+=1
        if n>=20000: break
if lens:
    lens.sort()
    print("okunan okuma sayisi: %d"%len(lens))
    print("medyan uzunluk : %d bp"%statistics.median(lens))
    print("ortalama       : %d bp"%(sum(lens)/len(lens)))
    print("yuzde 10 / 90  : %d / %d bp"%(lens[len(lens)//10],lens[9*len(lens)//10]))
    print("NOTE: the source study used database300mers.kmer_distrib. The median is far")
    print("     uzerindeyse Bracken icin uygun uzunlukta kmer_distrib")
    print("     uretilmeli (bracken-build -l <uzunluk>).")
else:
    print("kraken _output dosyalari okunamadi")
PY
else
  echo "The project directory was not found under /mnt/c, give the path by hand."
fi

bar "8. Disk boslugu"
df -h "$HOME" /mnt/c 2>/dev/null | sed 's/^/  /'

bar "BITTI"
echo "Share this file and the yardim_ciktilari/ directory; the scripts are adjusted"
echo "kesinlestirecegim. Hicbir bayragi hafizadan varsaymayacagim."
