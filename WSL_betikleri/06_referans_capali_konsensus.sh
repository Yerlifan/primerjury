#!/usr/bin/env bash
# =====================================================================
# 06_referans_capali_konsensus.sh
# Amaç: her taksonun konsensüsünü keyfî bir çekirdek okuma yerine ORTAK BİR
#       REFERANS dizisinin koordinat sisteminde kurmak.
#
# Neden gerekli:
#   consensus2.sh her takson için okumalardan bir çekirdek okuma seçip
#   konsensüsü onun etrafında kuruyor. Okuma uzunluğu ve kapsama değişince
#   çekirdek okuma operonun başka bir parçasına düşüyor. Ölçtüm: F1 grubunda
#   taxid 44689'un dört yıla ait konsensüsü üç ayrı bölgeyi kapsıyor
#   (F1-2 ile F1-3 arasında 12-mer benzerliği 0,987 ama F1-1 ile F1-2
#   arasında 0,065). Yani konsensüsler tanımlı bir amplikon değil, keyfî bir
#   çekirdek okumanın etrafındaki yerel birleştirme.
#   Referans çapalı kurulumda dört yılın konsensüsü aynı pencereye oturur,
#   hem yıllar arası karşılaştırma hem tek bir universal primer mümkün olur.
#
# Yöntem, takson başına:
#   1. Okumalardan bir örneklem alınıp uygun BLAST veritabanına karşı blastn
#      ile en iyi referans seçilir (en yüksek toplam bit skoru).
#   2. Referans dizi blastdbcmd ya da samtools faidx ile çıkarılır.
#      SILVA gibi RNA alfabeli veritabanlarında U harfleri T'ye çevrilir.
#   3. Bütün okumalar minimap2 ile o referansa hizalanır.
#   4. samtools consensus iki bağımsız ölçüm üretir: -a -A ile IUPAC kodlu
#      konsensüs, -a -A -f PILEUP ile pozisyon başına derinlik ve baz dizisi.
#   5. Çıktı referans koordinatındadır, yani bütün yıllar hizalıdır.
#
# Kullanım:
#   bash 06_referans_capali_konsensus.sh \
#        --pt  "/mnt/c/Users/yerli/Masaüstü/PrimerTasarlama" \
#        --out "/mnt/c/Users/yerli/Masaüstü/PrimerTasarlama/referans_konsensus" \
#        [--groups F1,F2] [--threads N] [--sample 50] [--min-depth N]
#        [--db ROD_v1.2_operon_variants.fasta]   REFERANS_DB icinde, virgulle
# =====================================================================
set -euo pipefail

# DIKKAT: GROUPS bash'in yerlesik dizisidir (kullanicinin grup kimlikleri),
# ona atama sessizce yok sayilir. Bu yuzden GRUP_SEC adi kullanilir.
PT=""; OUT=""; THREADS=""; SAMPLE=50; GRUP_SEC=""; MINDEPTH=""; DB_OVERRIDE=""
while [ $# -gt 0 ]; do
  case "$1" in
    --pt) PT="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --sample) SAMPLE="$2"; shift 2;;
    --groups) GRUP_SEC="$2"; shift 2;;
    --min-depth) MINDEPTH="$2"; shift 2;;
    --db) DB_OVERRIDE="$2"; shift 2;;
    *) echo "bilinmeyen secenek: $1" >&2; exit 2;;
  esac
done
[ -n "$PT" ] && [ -n "$OUT" ] || {
  echo "kullanim: bash $0 --pt <PrimerTasarlama> --out <cikti>" >&2; exit 2; }
[ -d "$PT/fastq files" ] || { echo "HATA: '$PT/fastq files' yok" >&2; exit 1; }
[ -d "$PT/REFERANS_DB" ] || { echo "HATA: '$PT/REFERANS_DB' yok" >&2; exit 1; }

log(){ printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die(){ printf 'HATA: %s\n' "$*" >&2; exit 1; }
for t in minimap2 samtools blastn python3; do
  command -v "$t" >/dev/null 2>&1 || die "$t PATH'te yok"
done
HAS_BDC=0; command -v blastdbcmd >/dev/null 2>&1 && HAS_BDC=1
[ -z "$THREADS" ] && THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))

# bayrak dogrulamasi, hafizadan varsayim yok
CH=$(samtools consensus 2>&1 || true)
for fl in "--ambig" "--format"; do
  grep -q -- "$fl" <<<"$CH" || die "samtools consensus $fl tanimiyor"
done
log "araclar hazir, is parcacigi=$THREADS, blastdbcmd=$HAS_BDC"

mkdir -p "$OUT"/{ref,bam,konsensus,pileup,maske,blast,log}
DB="$PT/REFERANS_DB"

# --- grup basina veritabani secimi ------------------------------------
# F gruplari icin once ROD (tam boy rDNA operonlari, yalnizca okaryot) denenir.
# Sebep: NCBI'nin fungi.ITS ve fungi.28S kayitlari tek bir alt bolgeyi kapsar
# (yaklasik 750-1450 bp). F1 okumalari 1204-1479 bp, F2 okumalari 3697-3700 bp
# oldugu icin amplikon referanstan uzundur ve referansin kenarlari okumalarla
# ortusmez; ilk kosumda F1 grubunda 900 bp referansin ancak 600 bp'si kapsandi.
# ROD tam operon tuttugu icin amplikonun tamamini kapsar.
# A ve B gruplarinda ROD kullanilamaz, cunku ROD bakteri ve arkeyi disliyor.
db_for() {
  case "$1" in
    A*) echo "$DB/archaea.16S.fna $DB/SILVA_138.2_SSURef_NR99.fasta";;
    B*) echo "$DB/bacteria.16S.fna $DB/SILVA_138.2_SSURef_NR99.fasta";;
    F*) echo "$DB/ROD_v1.2_operon_variants.fasta $DB/fungi.ITS.fna $DB/fungi.28SrRNA.fna $DB/fungi.18SrRNA.fna";;
    *)  echo "";;
  esac
}
# --db ile elle secim, virgulle ayrilmis dosya adlari (REFERANS_DB icinde)
if [ -n "$DB_OVERRIDE" ]; then
  db_for() { echo "$DB_OVERRIDE" | tr ',' '\n' | sed "s#^#$DB/#" | tr '\n' ' '; }
fi

MAP="$OUT/referans_secimi.tsv"
printf 'grup\ttaxid\tfastq\tveritabani\treferans_id\tref_uzunluk\tbit_toplam\thit_okuma\n' > "$MAP"

shopt -s nullglob
TOTAL=0; DONE=0; SKIPPED=0
for fq in "$PT/fastq files"/*/*.fastq; do TOTAL=$((TOTAL+1)); done
log "toplam fastq: $TOTAL"
[ "$TOTAL" -gt 0 ] || die "hic fastq bulunamadi: $PT/fastq files"
if [ -n "$GRUP_SEC" ]; then log "grup suzgeci: $GRUP_SEC"; fi

for fq in "$PT/fastq files"/*/*.fastq; do
  base=$(basename "$fq" .fastq)
  grp=$(basename "$(dirname "$fq")")
  if [ -n "$GRUP_SEC" ]; then
    case ",$GRUP_SEC," in
      *",${grp%%-*},"*) ;;
      *) SKIPPED=$((SKIPPED+1)); continue;;
    esac
  fi
  tid=$(echo "$base" | sed -n 's/.*reads[-_]\([0-9]\+\).*/\1/p')
  [ -n "$tid" ] || { echo "taxid cikarilamadi, atlandi: $fq" >&2; continue; }
  tag="${grp}_${tid}"
  DONE=$((DONE+1))
  log "[$DONE] $tag"

  reffa="$OUT/ref/${tag}_ref.fasta"
  if [ ! -s "$reffa" ]; then
    # 1. okuma ornegi -> fasta
    q="$OUT/blast/${tag}_sorgu.fasta"
    awk -v n="$SAMPLE" 'NR%4==1{h=substr($1,2)} NR%4==2{print ">"h"\n"$0; c++; if(c>=n) exit}' \
        "$fq" > "$q"
    [ -s "$q" ] || { echo "  okuma yok, atlandi" >&2; continue; }
    # 2. aday veritabanlarina blastn, en yuksek toplam bit skoru kazanir
    best_db=""; best_id=""; best_bit=0; best_n=0
    for d in $(db_for "$grp"); do
      [ -e "$d.nin" ] || { echo "  indeks yok, atlandi: $d" >&2; continue; }
      bo="$OUT/blast/${tag}_$(basename "$d").tsv"
      blastn -query "$q" -db "$d" -outfmt '6 qseqid sseqid bitscore pident length' \
             -max_target_seqs 5 -evalue 1e-20 -num_threads "$THREADS" \
             > "$bo" 2> "$OUT/log/${tag}_blast.log" || true
      [ -s "$bo" ] || continue
      read -r sid bit nq < <(awk -F'\t' '{b[$2]+=$3; q[$2"|"$1]=1}
          END{for(s in b){n=0; for(k in q){split(k,a,"|"); if(a[1]==s) n++}
              if(b[s]>m){m=b[s]; ms=s; mn=n}} print ms, m+0, mn+0}' "$bo")
      if [ -n "${sid:-}" ] && awk -v a="$bit" -v b="$best_bit" 'BEGIN{exit !(a>b)}'; then
        best_db="$d"; best_id="$sid"; best_bit="$bit"; best_n="$nq"
      fi
    done
    [ -n "$best_id" ] || { echo "  referans bulunamadi, atlandi" >&2; continue; }
    # 3. referans diziyi cikar, RNA alfabesi varsa T'ye cevir
    if [ "$HAS_BDC" = 1 ] && blastdbcmd -db "$best_db" -entry "$best_id" \
         > "$reffa.tmp" 2>/dev/null && [ -s "$reffa.tmp" ]; then :
    else
      [ -e "$best_db.fai" ] || samtools faidx "$best_db"
      samtools faidx "$best_db" "$best_id" > "$reffa.tmp" || {
        echo "  referans cikarilamadi: $best_id" >&2; rm -f "$reffa.tmp"; continue; }
    fi
    python3 - "$reffa.tmp" "$reffa" <<'PY'
import sys
inp,outp=sys.argv[1],sys.argv[2]
name=None; seq=[]
for l in open(inp):
    if l.startswith(">"):
        if name is None: name=l.strip()
        else: break            # yalnizca ilk kayit
    else: seq.append(l.strip())
s="".join(seq).upper().replace("U","T")     # SILVA gibi RNA alfabeleri icin
open(outp,"w").write((name or ">ref")+"\n"+"\n".join(s[i:i+70] for i in range(0,len(s),70))+"\n")
print("  referans uzunlugu: %d bp"%len(s))
PY
    rm -f "$reffa.tmp"
    rl=$(python3 -c "import sys;print(sum(len(l.strip()) for l in open(sys.argv[1]) if not l.startswith('>')))" "$reffa")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$grp" "$tid" "$(basename "$fq")" "$(basename "$best_db")" "$best_id" \
      "$rl" "$best_bit" "$best_n" >> "$MAP"
  fi

  # 4. butun okumalari referansa hizala
  bam="$OUT/bam/${tag}.bam"
  if [ ! -s "$bam" ]; then
    minimap2 -ax map-ont -t "$THREADS" "$reffa" "$fq" 2> "$OUT/log/${tag}_minimap2.log" \
      | samtools sort -@ "$THREADS" -o "$bam" 2>> "$OUT/log/${tag}_minimap2.log"
    samtools index "$bam"
  fi

  # 5. iki bagimsiz olcum, KATI referans koordinatinda
  #    --show-ins no  : eklemeler gosterilmez, aksi halde cikti referanstan
  #                     uzar ve koordinatlar kayar
  #    --show-del yes : silmeler * olarak korunur, aksi halde cikti kisalir
  #    Ikisi birlikte cikti uzunlugunu referansa esitler. Sentetik BAM ile
  #    dogrulandi: 200 bp referans, 6 baz ekleme ve 5 baz silme tasiyan
  #    okumalar; bayraksiz 201, yalniz --show-ins no ile 195, yalniz
  #    --show-del yes ile 206, ikisi birlikte tam 200.
  samtools consensus -a -A --show-ins no --show-del yes \
    -o "$OUT/konsensus/${tag}_ref_konsensus.fasta.raw" "$bam" \
    2> "$OUT/log/${tag}_cons.log"
  # silme isareti * primer tasarimi icin bir baz degildir; N'e cevrilir ki
  # 02'nin maskesi orayi dusuk_derinlik gibi yasaklasin
  python3 - "$OUT/konsensus/${tag}_ref_konsensus.fasta.raw" \
            "$OUT/konsensus/${tag}_ref_konsensus.fasta" <<'PY2'
import sys
inp,outp=sys.argv[1],sys.argv[2]
hdr=None; seq=[]
for l in open(inp):
    if l.startswith(">"):
        if hdr is None: hdr=l.rstrip("\n")
    else: seq.append(l.strip())
s="".join(seq).upper().replace("*","N")
open(outp,"w").write((hdr or ">ref")+"\n"+"\n".join(s[i:i+70] for i in range(0,len(s),70))+"\n")
PY2
  rm -f "$OUT/konsensus/${tag}_ref_konsensus.fasta.raw"
  samtools consensus -a -A --show-ins no --show-del yes -f PILEUP \
    -o "$OUT/pileup/${tag}_pileup.txt" "$bam" 2> "$OUT/log/${tag}_pileup.log"
done

# --- islenen dosya yoksa sessizce bitme --------------------------------
log "islenen: $DONE, grup suzgeciyle atlanan: $SKIPPED, toplam: $TOTAL"
if [ "$DONE" -eq 0 ]; then
  echo >&2
  echo "HATA: hicbir dosya islenmedi." >&2
  if [ "$SKIPPED" -eq "$TOTAL" ] && [ -n "$GRUP_SEC" ]; then
    echo "Sebep: --groups '$GRUP_SEC' hicbir klasore uymadi." >&2
    echo "Mevcut grup onekleri:" >&2
    for d in "$PT/fastq files"/*/; do basename "$d"; done | sed 's/-[0-9]*$//' \
      | sort -u | sed 's/^/   /' >&2
  else
    echo "Sebep: fastq adlarindan taxid cikarilamadi ya da referans bulunamadi." >&2
    echo "Ayrinti icin $OUT/log/ altindaki dosyalara bakin." >&2
  fi
  exit 1
fi

# --- ozet ve hizalilik denetimi ---------------------------------------
log "ozet"
python3 - "$OUT" "$MAP" "${MINDEPTH:-AUTO}" <<'PY'
import sys,os,glob,csv,re,statistics,collections
out,mapf,mind=sys.argv[1],sys.argv[2],sys.argv[3]
rows=list(csv.DictReader(open(mapf),delimiter="\t"))
print("referans secilen takson-barkod: %d"%len(rows))
byref=collections.defaultdict(list)
for r in rows: byref[(r["taxid"],r["referans_id"])].append(r["grup"])
# Ayni taksonun farkli barkodlarda AYNI referansa oturmasi, hizaliligin sarti
bytax=collections.defaultdict(set)
for r in rows: bytax[r["taxid"]].add(r["referans_id"])
coklu={t:v for t,v in bytax.items() if len(v)>1}
print("ayni taksonda farkli referans secilen: %d"%len(coklu))
for t,v in sorted(coklu.items()):
    print("   taxid %-9s -> %s"%(t,", ".join(sorted(v))))
print("   (bu taksonlarda referansi elle sabitlemek gerekir, aksi halde")
print("    yillar yine ayni koordinat sistemine oturmaz; 07 betigi bunu yapar)")
print()
print("### REFERANS TAKSONOMISI UYARISI")
print("ROD genomlardan turetilmis bir veritabani ve yalnizca 11.935 genom")
print("kapsiyor. Hedef taksonun kendisi yoksa blastn en yakin akrabayi secer.")
print("Asagida her hedefin oturdugu referansin taksonomisi var; hedefle ayni")
print("cins ya da aile degilse konsensus yabanci bir iskelet uzerine kuruluyor")
print("demektir. Dizinin kendisi okumalardan cagriliyor, ama indeller ve")
print("kenarlar referansa gore yorumlanir.")
for r in sorted(rows,key=lambda x:(x["grup"],x["taxid"])):
    rid=r["referans_id"]
    tax=rid.split("|")[2] if rid.count("|")>=2 else ""
    parts=[p for p in tax.split(";") if p][-2:] if tax else []
    print("   %-6s taxid %-9s -> %s"%(r["grup"],r["taxid"],
          " / ".join(parts) if parts else rid[:60]))
# kapsama ozeti
print()
# KOORDINAT DENETIMI: konsensus uzunlugu referans uzunluguna esit olmali
reflen={}
for r in rows:
    p=os.path.join(out,"ref","%s_%s_ref.fasta"%(r["grup"],r["taxid"]))
    if os.path.exists(p):
        reflen[(r["grup"],r["taxid"])]=sum(len(l.strip()) for l in open(p) if not l.startswith(">"))
print()
print("%-26s %8s %8s %8s %8s %s"%("hedef","ref_uz","kapsanan","N","medyan_der","koordinat"))
for p in sorted(glob.glob(os.path.join(out,"konsensus","*_ref_konsensus.fasta"))):
    tag=os.path.basename(p).replace("_ref_konsensus.fasta","")
    s="".join(l.strip() for l in open(p) if not l.startswith(">")).upper()
    pil=os.path.join(out,"pileup","%s_pileup.txt"%tag)
    dep=[]
    if os.path.exists(pil):
        for line in open(pil,errors="replace"):
            f=line.rstrip("\n").split("\t")
            if len(f)>=4 and f[3].isdigit(): dep.append(int(f[3]))
    cov=sum(1 for c in s if c not in "N")
    key=tuple(tag.rsplit("_",1))
    rl=reflen.get(key)
    koord = "TAMAM" if (rl is not None and rl==len(s)) else ("KAYIK(%s)"%rl if rl else "?")
    print("%-26s %8d %8d %8d %8s %s"%(tag,len(s),cov,s.count("N"),
          int(statistics.median([d for d in dep if d>0])) if any(dep) else "-", koord))
PY

log "bitti"
echo
echo "Cikti:"
echo "  $OUT/referans_secimi.tsv        hangi takson hangi referansa oturdu"
echo "  $OUT/konsensus/*_ref_konsensus.fasta  referans koordinatinda konsensus"
echo "  $OUT/pileup/*_pileup.txt        pozisyon basina derinlik ve IUPAC cagri"
echo
echo "Sonraki adim: bu konsensusleri 02_N_analizi.sh ile isleyip maske uretin,"
echo "sonra 04_grup_primer.py'yi bunlarla calistirin. Ayni taksonun butun"
echo "yillari artik ayni koordinat sisteminde oldugu icin yon normalizasyonu"
echo "da gereksizlesir, ancak betik yine de kontrol eder."
