#!/usr/bin/env bash
# =====================================================================
# freeze_reference.sh
# The aim: to rebuild the consensuses the anchored step produced so that every
#          year of one taxon sits in ONE AND THE SAME coordinate system.
#
# There are two modes. They write into separate folders, so the two can be
# compared side by side.
#
#   --mode pin   (the default)
#     Pins an outside reference. The anchored step chooses the best reference
#     independently for every taxon and barcode pair, and on some taxa the years
#     land on different references. This mode declares the reference with the
#     highest total bit score per taxon the winner and aligns every year of that
#     taxon to it.
#     Its advantage: the coordinates stay comparable with the outside world.
#
#   --mode self
#     Takes your own data as the reference. For every taxon the best covered
#     consensus in the anchored output is chosen, its uncovered ends are clipped,
#     and every year of that taxon is aligned to that sequence.
#     The reasoning: ROD covers only 11,935 genomes, and when the target taxon
#     itself is absent blastn picks the nearest relative. In a real run taxid 4896
#     (Schizosaccharomyces pombe, an Ascomycota yeast) landed on Earliella
#     scabrosa, that is on a Basidiomycota fungus, and taxid 44689 (Dictyostelium,
#     Amoebozoa) landed on Pleurotus giganteus. The sequence itself is called from
#     the reads, but a foreign scaffold spoils how the indels and the edges are
#     read. Taking your own data as the reference removes that problem entirely.
#     Its drawback: the coordinates are specific to this study alone.
#
# In both modes the consensus is called in STRICT reference coordinates:
#   --show-ins no --show-del yes makes the output length equal the reference's
#   (confirmed with a synthetic BAM: a 200 bp reference with reads carrying a 6
#   base insertion and a 5 base deletion; with no flags 201, with ins-no alone
#   195, with del-yes alone 206, and with both exactly 200). The deletion mark *
#   is turned into an N afterwards, because a deleted position is not a base and
#   no primer footprint can sit there.
#
# Usage:
#   bash freeze_reference.sh --pt <project> --out <anchored output> \
#        [--mode pin|self] [--threads N] [--dry-run]
# =====================================================================
set -euo pipefail

PT=""; OUT=""; THREADS=""; DRY=0; MODE="pin"
while [ $# -gt 0 ]; do
  case "$1" in
    --pt) PT="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --mode) MODE="$2"; shift 2;;
    --self-reference) MODE="self"; shift;;
    --dry-run) DRY=1; shift;;
    *) echo "unknown option: $1" >&2; exit 2;;
  esac
done
[ -n "$PT" ] && [ -n "$OUT" ] || {
  echo "usage: bash $0 --pt <project directory> --out <freeze output directory> [--mode pin|self]" >&2
  exit 2; }
case "$MODE" in pin|self) ;; *) echo "ERROR: --mode has to be pin or self" >&2; exit 2;; esac
case "$PT$OUT" in *'...'*) echo "ERROR: the path contains '...'. Give a full path, not a placeholder." >&2; exit 2;; esac
[ -d "$PT" ] || { echo "ERROR: no such --pt directory: $PT" >&2; exit 1; }
MAP="$OUT/referans_secimi.tsv"
[ -s "$MAP" ] || { echo "ERROR: $MAP is missing. Run the mapping step first." >&2; exit 1; }

log(){ printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die(){ printf 'ERROR: %s\n' "$*" >&2; exit 1; }
for t in minimap2 samtools python3; do
  command -v "$t" >/dev/null 2>&1 || die "$t is not on PATH"
done
[ -z "$THREADS" ] && THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))
log "the mode: $MODE, threads: $THREADS"

konsensus_cagir() {   # $1=bam  $2=output fasta  $3=output pileup  $4=log prefix
  samtools consensus -a -A --show-ins no --show-del yes -o "$2.raw" "$1" 2> "$4_cons.log"
  python3 - "$2.raw" "$2" <<'PY2'
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
  rm -f "$2.raw"
  samtools consensus -a -A --show-ins no --show-del yes -f PILEUP -o "$3" "$1" 2> "$4_pileup.log"
}

fastq_bul() { ls "$PT/fastq files/$1"/*reads[-_]"$2".fastq 2>/dev/null | head -1; }

# =====================================================================
if [ "$MODE" = "pin" ]; then
# =====================================================================
PLAN="$OUT/sabitleme_plani.tsv"
python3 - "$MAP" "$PLAN" <<'PY'
import sys,csv,collections
mapf,planf=sys.argv[1],sys.argv[2]
rows=list(csv.DictReader(open(mapf),delimiter="\t"))
bytax=collections.defaultdict(list)
for r in rows: bytax[r["taxid"]].append(r)
plan=[]; degisen=0
print("%-9s %-6s %-16s %-16s %s"%("taxid","barcode","old","new","state"))
for tid,rs in sorted(bytax.items(),key=lambda x:int(x[0])):
    agg=collections.defaultdict(lambda:[0.0,0])
    for r in rs:
        try: b=float(r["bit_toplam"])
        except ValueError: b=0.0
        try: L=int(r["ref_uzunluk"])
        except ValueError: L=0
        agg[r["referans_id"]][0]+=b
        agg[r["referans_id"]][1]=max(agg[r["referans_id"]][1],L)
    win=sorted(agg.items(),key=lambda kv:(-kv[1][0],-kv[1][1],kv[0]))[0][0]
    if len(agg)==1: continue
    for r in rs:
        st="same" if r["referans_id"]==win else "WILL_CHANGE"
        if st=="WILL_CHANGE": degisen+=1
        print("%-9s %-6s %-16s %-16s %s"%(tid,r["grup"],
              r["referans_id"].split("|")[0][:16],win.split("|")[0][:16],st))
        plan.append((r["grup"],tid,r["referans_id"],win,st))
with open(planf,"w",newline="",encoding="utf-8") as fh:
    w=csv.writer(fh,delimiter="\t",lineterminator="\n")
    w.writerow(["grup","taxid","eski_referans","yeni_referans","durum"]); w.writerows(plan)
print()
print("taxa with more than one reference: %d"%sum(1 for t,rs in bytax.items()
      if len(set(r["referans_id"] for r in rs))>1))
print("targets to be realigned: %d"%degisen)
PY
[ "$DRY" = 1 ] && { log "a dry run. The plan: $PLAN"; exit 0; }

n=0
while IFS=$'\t' read -r grp tid eski yeni durum; do
  [ "$durum" = "WILL_CHANGE" ] || continue
  tag="${grp}_${tid}"; src=""
  for cand in "$OUT/ref"/*_ref.fasta; do
    head -1 "$cand" | grep -qF "$yeni" && { src="$cand"; break; }
  done
  [ -n "$src" ] || { echo "  the winning reference file was not found ($tag)" >&2; continue; }
  fq=$(fastq_bul "$grp" "$tid")
  [ -n "$fq" ] || { echo "  no fastq found: $grp $tid" >&2; continue; }
  n=$((n+1)); log "[$n] realigning $tag"
  cp -f "$src" "$OUT/ref/${tag}_ref.fasta"
  rm -f "$OUT/bam/${tag}.bam" "$OUT/bam/${tag}.bam.bai"
  minimap2 -ax map-ont -t "$THREADS" "$OUT/ref/${tag}_ref.fasta" "$fq" \
      2> "$OUT/log/${tag}_mm2_sabit.log" \
    | samtools sort -@ "$THREADS" -o "$OUT/bam/${tag}.bam" 2>> "$OUT/log/${tag}_mm2_sabit.log"
  samtools index "$OUT/bam/${tag}.bam"
  konsensus_cagir "$OUT/bam/${tag}.bam" "$OUT/konsensus/${tag}_ref_konsensus.fasta" \
                  "$OUT/pileup/${tag}_pileup.txt" "$OUT/log/${tag}_sabit"
done < <(tail -n +2 "$PLAN")
log "targets realigned: $n"
BEKLENEN=$(awk -F'\t' 'NR>1 && $5=="WILL_CHANGE"' "$PLAN" | wc -l)
[ "$n" -eq "$BEKLENEN" ] || echo "WARNING: the plan said $BEKLENEN, $n were processed." >&2
DENETIM_DIZIN="$OUT/konsensus"

# =====================================================================
else   # MODE = self
# =====================================================================
SELF="$OUT/self"
mkdir -p "$SELF"/{ref,bam,konsensus,pileup,log}
SECIM="$SELF/self_referans_secimi.tsv"

log "choosing the best covered consensus for every taxon"
python3 - "$OUT" "$SECIM" <<'PY'
import sys,glob,os,re,csv,collections,statistics
out,secim=sys.argv[1],sys.argv[2]
cand=collections.defaultdict(list)
for p in sorted(glob.glob(os.path.join(out,"konsensus","*_ref_konsensus.fasta"))):
    tag=os.path.basename(p).replace("_ref_konsensus.fasta","")
    m=re.match(r"(.+)_(\d+)$",tag)
    if not m: continue
    grp,tid=m.group(1),m.group(2)
    s="".join(l.strip() for l in open(p) if not l.startswith(">")).upper()
    pil=os.path.join(out,"pileup","%s_pileup.txt"%tag)
    dep=[]
    if os.path.exists(pil):
        for line in open(pil,errors="replace"):
            f=line.rstrip("\n").split("\t")
            if len(f)>=4 and f[3].isdigit() and int(f[3])>0: dep.append(int(f[3]))
    cand[tid].append(dict(tag=tag,grup=grp,path=p,
                          kapsanan=sum(1 for c in s if c!="N"),
                          medyan=int(statistics.median(dep)) if dep else 0,
                          uzunluk=len(s)))
rows=[]
print("%-9s %-4s %-26s %10s %11s"%("taxid","members","chosen","covered","median_depth"))
for tid,items in sorted(cand.items(),key=lambda x:int(x[0])):
    # the most covered wins; on a tie the deeper one, and on a tie there the name order
    best=sorted(items,key=lambda d:(-d["kapsanan"],-d["medyan"],d["tag"]))[0]
    print("%-9s %-4d %-26s %10d %11d"%(tid,len(items),best["tag"],
          best["kapsanan"],best["medyan"]))
    rows.append((tid,best["tag"],best["path"],best["kapsanan"],best["medyan"],
                 ";".join(d["tag"] for d in items)))
with open(secim,"w",newline="",encoding="utf-8") as fh:
    w=csv.writer(fh,delimiter="\t",lineterminator="\n")
    w.writerow(["taxid","secilen_tag","kaynak_dosya","kapsanan","medyan_derinlik","uyeler"])
    w.writerows(rows)
print()
print("the number of taxa: %d"%len(rows))
PY
[ "$DRY" = 1 ] && { log "a dry run. The choice: $SECIM"; exit 0; }

log "clipping the self references and realigning every year"
n=0
while IFS=$'\t' read -r tid stag spath kap med uyeler; do
  ref="$SELF/ref/${tid}_self_ref.fasta"
  # the uncovered ends are clipped and the inner N kept (we invent no information)
  python3 - "$spath" "$ref" "$tid" "$stag" <<'PY2'
import sys
src,dst,tid,stag=sys.argv[1:5]
s="".join(l.strip() for l in open(src) if not l.startswith(">")).upper()
i=0
while i<len(s) and s[i]=="N": i+=1
j=len(s)
while j>i and s[j-1]=="N": j-=1
core=s[i:j]
if len(core)<200:
    sys.exit("  taxid %s: the clipped reference is too short (%d bp), skipped"%(tid,len(core)))
open(dst,"w").write(">self_%s kaynak=%s kirpma=%d-%d uzunluk=%d\n"
                    %(tid,stag,i+1,j,len(core))
                    +"\n".join(core[k:k+70] for k in range(0,len(core),70))+"\n")
print("  taxid %-9s %d bp -> clipped to %d bp (inner N: %d)"%(tid,len(s),len(core),core.count("N")))
PY2
  [ -s "$ref" ] || { echo "  no reference could be produced, taxon skipped: $tid" >&2; continue; }
  for tag in $(echo "$uyeler" | tr ';' ' '); do
    grp="${tag%_*}"; t2="${tag##*_}"
    fq=$(fastq_bul "$grp" "$t2")
    [ -n "$fq" ] || { echo "  no fastq found: $grp $t2" >&2; continue; }
    n=$((n+1))
    minimap2 -ax map-ont -t "$THREADS" "$ref" "$fq" 2> "$SELF/log/${tag}_mm2.log" \
      | samtools sort -@ "$THREADS" -o "$SELF/bam/${tag}.bam" 2>> "$SELF/log/${tag}_mm2.log"
    samtools index "$SELF/bam/${tag}.bam"
    konsensus_cagir "$SELF/bam/${tag}.bam" "$SELF/konsensus/${tag}_self_konsensus.fasta" \
                    "$SELF/pileup/${tag}_pileup.txt" "$SELF/log/${tag}"
  done
done < <(tail -n +2 "$SECIM")
log "targets realigned: $n"
DENETIM_DIZIN="$SELF/konsensus"
fi

# =====================================================================
# the final check, shared by both modes
# =====================================================================
log "the final check: $DENETIM_DIZIN"
python3 - "$DENETIM_DIZIN" "$OUT" <<'PY'
import sys,glob,os,re,collections
d,out=sys.argv[1],sys.argv[2]
seqs={}
for p in sorted(glob.glob(os.path.join(d,"*_konsensus.fasta"))):
    tag=re.sub(r"_(ref|self)_konsensus\.fasta$","",os.path.basename(p))
    hdr=open(p).readline().strip()
    s="".join(l.strip() for l in open(p) if not l.startswith(">")).upper()
    seqs[tag]=(hdr,s)
bytax=collections.defaultdict(list)
for tag,(h,s) in seqs.items():
    m=re.search(r"_(\d+)$",tag)
    if m: bytax[m.group(1)].append((tag,h,s))
bad=0
print("%-9s %-4s %-16s %-24s %s"%("taxid","years","length","covered","aligned"))
for tid,items in sorted(bytax.items(),key=lambda x:int(x[0])):
    L=set(len(s) for _,_,s in items)
    H=set(h.split()[0] for _,h,_ in items)
    ok=len(L)==1 and len(H)==1
    if not ok: bad+=1
    print("%-9s %-4d %-16s %-24s %s"%(tid,len(items),
          ",".join(str(x) for x in sorted(L)),
          ",".join(str(sum(1 for c in s if c!="N")) for _,_,s in items),
          "YES" if ok else "NO"))
print()
print("taxa that are not aligned: %d"%bad)
print("In every taxon the years are in the same coordinate system." if bad==0 else
      "On the taxa marked NO the coordinates still diverge.")
base=os.path.join(out,"konsensus")
if os.path.abspath(base)!=os.path.abspath(d):
    print()
    print("### this mode compared with the anchored step (an outside reference)")
    print("%-26s %12s %12s %8s %10s"%("target","old_covered","new_covered","diff","new_length"))
    tot0=tot1=0
    for tag,(h,s) in sorted(seqs.items()):
        p0=os.path.join(base,"%s_ref_konsensus.fasta"%tag)
        if not os.path.exists(p0): continue
        s0="".join(l.strip() for l in open(p0) if not l.startswith(">")).upper()
        c0=sum(1 for c in s0 if c!="N"); c1=sum(1 for c in s if c!="N")
        tot0+=c0; tot1+=c1
        print("%-26s %12d %12d %+8d %10d"%(tag,c0,c1,c1-c0,len(s)))
    print("%-26s %12d %12d %+8d"%("TOTAL",tot0,tot1,tot1-tot0))
    print()
    print("The number of covered bases measures the area open to primer design.")
    print("In self mode the reference is our own covered region already, so a high")
    print("covered fraction is expected; the real gain is that the scaffold is")
    print("taxonomically right.")
PY
log "finished"
echo
echo "Next step: analyze_ambiguous_bases.sh processes these consensuses into a mask,"
echo "after which design_group_primers.py and split_clusters.py take design back up."
