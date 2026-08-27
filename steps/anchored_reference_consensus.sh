#!/usr/bin/env bash
# =====================================================================
# anchored_reference_consensus.sh
# The aim: to build each taxon's consensus in the coordinate system of A SHARED
#          REFERENCE sequence rather than around an arbitrary seed read.
#
# Why it is needed:
#   consensus2.sh picks a seed read out of the reads for each taxon and builds the
#   consensus around it. When the read length and the coverage change, the seed read
#   falls on another part of the operon. Measured: in group F1 the consensuses of
#   taxid 44689 across four years cover three separate regions (the 12-mer
#   similarity between F1-2 and F1-3 is 0.987, and between F1-1 and F1-2 it is
#   0.065). So the consensuses are not a defined amplicon but a local assembly
#   around an arbitrary seed read.
#   Anchored to a reference, the consensuses of four years sit in the same window,
#   which makes both a year to year comparison and a single universal primer
#   possible.
#
# The method, per taxon:
#   1. A sample of the reads is taken and the best reference is chosen with blastn
#      against a suitable BLAST database (the highest total bit score).
#   2. The reference sequence is pulled out with blastdbcmd or samtools faidx. In
#      databases with an RNA alphabet such as SILVA the U's are turned into T.
#   3. Every read is aligned to that reference with minimap2.
#   4. samtools consensus produces two independent measurements: an IUPAC coded
#      consensus with -a -A, and the depth plus base string per position with
#      -a -A -f PILEUP.
#   5. The output is in reference coordinates, so every year is aligned.
#
# Usage:
#   bash anchored_reference_consensus.sh \
#        --pt  /path/to/project \
#        --out /path/to/project/referans_konsensus \
#        [--groups F1,F2] [--threads N] [--sample 50] [--min-depth N]
#        [--db ROD_v1.2_operon_variants.fasta]   inside REFERENCE_DB, comma separated
# =====================================================================
set -euo pipefail

# CAREFUL: GROUPS is a bash builtin array (the user's group ids) and assigning to
# it is silently ignored. That is why the name GRUP_SEC is used.
PT=""; OUT=""; THREADS=""; SAMPLE=50; GRUP_SEC=""; MINDEPTH=""; DB_OVERRIDE=""
DBKOK=""
while [ $# -gt 0 ]; do
  case "$1" in
    --pt) PT="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --sample) SAMPLE="$2"; shift 2;;
    --groups) GRUP_SEC="$2"; shift 2;;
    --min-depth) MINDEPTH="$2"; shift 2;;
    --db) DB_OVERRIDE="$2"; shift 2;;
    # --db-root: the directory the reference databases sit in.
    # MEASURED on 2026-08-25: with the reference databases on a Windows drive
    # mounted into WSL, the SILVA SSU BLAST index of about 302 MB is read from
    # scratch for every bin, and a bin took about 11 minutes. Copying the same
    # files onto the WSL disk and pointing this option at them brings the step
    # down to about 1 minute.
    --db-root) DBKOK="$2"; shift 2;;
    *) echo "unknown option: $1" >&2; exit 2;;
  esac
done
[ -n "$PT" ] && [ -n "$OUT" ] || {
  echo "usage: bash $0 --pt <project directory> --out <output>" >&2; exit 2; }
[ -d "$PT/fastq files" ] || { echo "ERROR: no such directory: '$PT/fastq files'" >&2; exit 1; }
[ -d "$PT/REFERENCE_DB" ] || { echo "ERROR: no such directory: '$PT/REFERENCE_DB'" >&2; exit 1; }

log(){ printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die(){ printf 'ERROR: %s\n' "$*" >&2; exit 1; }
for t in minimap2 samtools blastn python3; do
  command -v "$t" >/dev/null 2>&1 || die "$t is not on PATH"
done
HAS_BDC=0; command -v blastdbcmd >/dev/null 2>&1 && HAS_BDC=1
[ -z "$THREADS" ] && THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))

# flag verification; nothing is assumed from memory
CH=$(samtools consensus 2>&1 || true)
for fl in "--ambig" "--format"; do
  grep -q -- "$fl" <<<"$CH" || die "samtools consensus $fl tanimiyor"
done
log "araclar hazir, is parcacigi=$THREADS, blastdbcmd=$HAS_BDC"

mkdir -p "$OUT"/{ref,bam,konsensus,pileup,maske,blast,log}
DB="${DBKOK:-$PT/REFERENCE_DB}"

# --- choosing a database per group -------------------------------------
# For the F groups ROD (full length rDNA operons, eukaryotes only) is tried first.
# The reason: NCBI's fungi.ITS and fungi.28S records cover a single subregion (about
# 750-1450 bp). Because the F1 reads are 1204-1479 bp and the F2 reads 3697-3700 bp,
# the amplicon is longer than the reference and the reference's edges do not overlap
# the reads; on the first run only 600 bp of a 900 bp reference was covered in group
# F1. Because ROD holds the full operon it covers the whole amplicon.
# ROD cannot be used for the A and B groups, because ROD excludes bacteria and
# archaea.
db_for() {
  case "$1" in
    A*) echo "$DB/archaea.16S.fna $DB/SILVA_138.2_SSURef_NR99.fasta";;
    B*) echo "$DB/bacteria.16S.fna $DB/SILVA_138.2_SSURef_NR99.fasta";;
    F*) echo "$DB/ROD_v1.2_operon_variants.fasta $DB/fungi.ITS.fna $DB/fungi.28SrRNA.fna $DB/fungi.18SrRNA.fna";;
    *)  echo "";;
  esac
}
# a hand selection with --db, comma separated file names (inside REFERENCE_DB)
if [ -n "$DB_OVERRIDE" ]; then
  db_for() { echo "$DB_OVERRIDE" | tr ',' '\n' | sed "s#^#$DB/#" | tr '\n' ' '; }
fi

MAP="$OUT/referans_secimi.tsv"
printf 'grup\ttaxid\tfastq\tveritabani\treferans_id\tref_uzunluk\tbit_toplam\thit_okuma\n' > "$MAP"

shopt -s nullglob
TOTAL=0; DONE=0; SKIPPED=0
for fq in "$PT/fastq files"/*/*.fastq; do TOTAL=$((TOTAL+1)); done
log "fastq files in total: $TOTAL"
[ "$TOTAL" -gt 0 ] || die "no fastq was found at all: $PT/fastq files"
if [ -n "$GRUP_SEC" ]; then log "the group filter: $GRUP_SEC"; fi

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
  [ -n "$tid" ] || { echo "the taxid could not be extracted, skipped: $fq" >&2; continue; }
  tag="${grp}_${tid}"
  DONE=$((DONE+1))
  # Resuming: a bin whose consensus has already been produced is not processed
  # again, so a run that was interrupted carries on from where it stopped.
  if [ -s "$OUT/konsensus/${tag}_ref_konsensus.fasta" ]; then
    log "[$DONE] $tag is already there, skipped"
    continue
  fi
  log "[$DONE] $tag"

  reffa="$OUT/ref/${tag}_ref.fasta"
  if [ ! -s "$reffa" ]; then
    # 1. a sample of the reads -> fasta
    q="$OUT/blast/${tag}_sorgu.fasta"
    awk -v n="$SAMPLE" 'NR%4==1{h=substr($1,2)} NR%4==2{print ">"h"\n"$0; c++; if(c>=n) exit}' \
        "$fq" > "$q"
    [ -s "$q" ] || { echo "  no reads, skipped" >&2; continue; }
    # 2. blastn against the candidate databases; the highest total bit score wins
    best_db=""; best_id=""; best_bit=0; best_n=0
    for d in $(db_for "$grp"); do
      [ -e "$d.nin" ] || { echo "  no index, skipped: $d" >&2; continue; }
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
    [ -n "$best_id" ] || { echo "  no reference found, skipped" >&2; continue; }
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
s="".join(seq).upper().replace("U","T")     # for RNA alphabets such as SILVA
open(outp,"w").write((name or ">ref")+"\n"+"\n".join(s[i:i+70] for i in range(0,len(s),70))+"\n")
print("  referans uzunlugu: %d bp"%len(s))
PY
    rm -f "$reffa.tmp"
    rl=$(python3 -c "import sys;print(sum(len(l.strip()) for l in open(sys.argv[1]) if not l.startswith('>')))" "$reffa")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$grp" "$tid" "$(basename "$fq")" "$(basename "$best_db")" "$best_id" \
      "$rl" "$best_bit" "$best_n" >> "$MAP"
  fi

  # 4. align every read to the reference
  bam="$OUT/bam/${tag}.bam"
  if [ ! -s "$bam" ]; then
    minimap2 -ax map-ont -t "$THREADS" "$reffa" "$fq" 2> "$OUT/log/${tag}_minimap2.log" \
      | samtools sort -@ "$THREADS" -o "$bam" 2>> "$OUT/log/${tag}_minimap2.log"
    samtools index "$bam"
  fi

  # 5. two independent measurements, in STRICT reference coordinates
  #    --show-ins no  : insertions are not shown; otherwise the output grows longer
  #                     than the reference and the coordinates shift
  #    --show-del yes : deletions are kept as *; otherwise the output shortens
  #    Together they make the output length equal the reference's. Confirmed with a
  #    synthetic BAM: a 200 bp reference with reads carrying a 6 base insertion and
  #    a 5 base deletion; with no flags 201, with --show-ins no alone 195, with
  #    --show-del yes alone 206, and with both exactly 200.
  samtools consensus -a -A --show-ins no --show-del yes \
    -o "$OUT/konsensus/${tag}_ref_konsensus.fasta.raw" "$bam" \
    2> "$OUT/log/${tag}_cons.log"
  # the deletion mark * is not a base for primer design; it is turned into an N so
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

# --- do not finish quietly when no file was processed -------------------
log "processed: $DONE, skipped by the group filter: $SKIPPED, in total: $TOTAL"
if [ "$DONE" -eq 0 ]; then
  echo >&2
  echo "ERROR: not one file was processed." >&2
  if [ "$SKIPPED" -eq "$TOTAL" ] && [ -n "$GRUP_SEC" ]; then
    echo "Cause: --groups '$GRUP_SEC' matched no directory." >&2
    echo "Mevcut grup onekleri:" >&2
    for d in "$PT/fastq files"/*/; do basename "$d"; done | sed 's/-[0-9]*$//' \
      | sort -u | sed 's/^/   /' >&2
  else
    echo "Cause: no taxid could be read from the fastq names, or no reference was found." >&2
    echo "For detail see the files under $OUT/log/." >&2
  fi
  exit 1
fi

# --- the summary and the alignment check --------------------------------
log "the summary"
python3 - "$OUT" "$MAP" "${MINDEPTH:-AUTO}" <<'PY'
import sys,os,glob,csv,re,statistics,collections
out,mapf,mind=sys.argv[1],sys.argv[2],sys.argv[3]
rows=list(csv.DictReader(open(mapf),delimiter="\t"))
print("referans secilen takson-barkod: %d"%len(rows))
byref=collections.defaultdict(list)
for r in rows: byref[(r["taxid"],r["referans_id"])].append(r["grup"])
# The same taxon sitting on THE SAME reference across different barcodes is the
# condition for them to be aligned
bytax=collections.defaultdict(set)
for r in rows: bytax[r["taxid"]].add(r["referans_id"])
coklu={t:v for t,v in bytax.items() if len(v)>1}
print("taxa where a different reference was chosen: %d"%len(coklu))
for t,v in sorted(coklu.items()):
    print("   taxid %-9s -> %s"%(t,", ".join(sorted(v))))
print("   (on those taxa the reference has to be fixed by hand, otherwise")
print("    the years still do not sit in one coordinate system)")
print()
print("### REFERANS TAKSONOMISI UYARISI")
print("ROD is a database derived from genomes and covers only 11,935 genomes.")
print("If the target taxon itself is absent, blastn picks the nearest relative.")
print("Below is the taxonomy of the reference each target sits on; where it is not")
print("cins ya da aile degilse konsensus yabanci bir iskelet uzerine kuruluyor")
print("demektir. Dizinin kendisi okumalardan cagriliyor, ama indeller ve")
print("the same as the target, the edges are read against the reference.")
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
print("%-26s %8s %8s %8s %8s %s"%("target","ref_len","covered","N","median_depth","coordinate"))
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

log "finished"
echo
echo "Output:"
echo "  $OUT/referans_secimi.tsv        which taxon settled on which reference"
echo "  $OUT/konsensus/*_ref_konsensus.fasta  referans koordinatinda konsensus"
echo "  $OUT/pileup/*_pileup.txt        depth and IUPAC call per position"
echo
echo "Next step: run these consensuses through analyze_ambiguous_bases.sh to build a mask,"
echo "then run design_group_primers.py with them. Since every year of the same taxon"
echo "now shares one coordinate system, orientation normalisation becomes"
echo "unnecessary, though the script still checks it."
