#!/usr/bin/env bash
# =====================================================================
# analyze_ambiguous_bases.sh
# The aim: to classify every N position in the consensus sequences with two
#          independent measurements, and to produce a mask for primer placement.
#          The consensus files ARE NOT TOUCHED and no N is replaced by a base.
#
# Why it is needed (confirmed with samtools 1.19.2):
#   The source study's consensus2.sh calls samtools consensus twice:
#       samtools consensus -a --use-qual              ("degenerate" in the comment)
#       samtools consensus -a --use-qual -c 0.9       ("strict" in the comment)
#   But samtools' default mode is "bayesian", and -c/--call-fract and -q/--use-qual
#   are options of the "simple" mode only. It was tested with a synthetic BAM: the
#   output of the two calls is IDENTICAL BYTE FOR BYTE, so -c 0.9 has no effect at
#   all. On top of that, because -A/--ambig is not given, samtools cannot print an
#   IUPAC code under any condition; even at a real biallelic position splitting 60
#   to 40 it writes an N. The same position comes back as M with -A.
#   The result: the N's in this project merge two different causes indistinguishably,
#   low depth and real biallelism. This script separates the two.
#
# Usage:
#   bash analyze_ambiguous_bases.sh --pt /path/to/project \
#        --out /path/to/project/N_analizi \
#        [--threads N] [--min-depth N] [--het-fract 0.15] [--config r10.4_sup]
# =====================================================================
set -euo pipefail

PT=""; OUT=""; THREADS=""; MINDEPTH=""; HET="0.15"; SMCONF=""
while [ $# -gt 0 ]; do
  case "$1" in
    --pt) PT="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --min-depth) MINDEPTH="$2"; shift 2;;
    --het-fract) HET="$2"; shift 2;;
    --config) SMCONF="$2"; shift 2;;
    *) echo "unknown option: $1" >&2; exit 2;;
  esac
done
[ -n "$PT" ] && [ -n "$OUT" ] || { echo "usage: bash $0 --pt <project> --out <output>" >&2; exit 2; }
[ -d "$PT/consensus sequences" ] || { echo "ERROR: no such directory: '$PT/consensus sequences'" >&2; exit 1; }

log(){ printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die(){ printf 'ERROR: %s\n' "$*" >&2; exit 1; }

for t in minimap2 samtools python3; do command -v "$t" >/dev/null 2>&1 || die "$t is not on PATH"; done
[ -z "$THREADS" ] && THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))

# --- flag verification; nothing is assumed from memory -----------------
CH=$(samtools consensus 2>&1 || true)
for fl in "--ambig" "--format" "--min-depth" "--mode"; do
  grep -q -- "$fl" <<<"$CH" || die "bu samtools surumunun consensus alt komutu $fl tanimiyor"
done
HAS_CONF=0; grep -q -- "--config" <<<"$CH" && HAS_CONF=1
if [ -n "$SMCONF" ]; then
  [ "$HAS_CONF" = 1 ] || die "--config istendi ama samtools desteklemiyor"
  grep -q "$SMCONF" <<<"$CH" || die "--config $SMCONF is not defined in this version"
fi
log "samtools consensus bayraklari dogrulandi (ambig, format, min-depth, mode; config=$HAS_CONF)"
samtools --version | head -1
minimap2 --version

mkdir -p "$OUT"/{bam,ambig,pileup,maske,log}
MAP_TSV="$OUT/hedef_fastq_eslemesi.tsv"; : > "$MAP_TSV"

# --- 1. matching a consensus to a fastq (by the barcode and taxid in the name) ---
log "konsensus ve fastq eslemesi kuruluyor"
python3 - "$PT" "$MAP_TSV" <<'PY'
import sys,glob,os,re
pt,outtsv=sys.argv[1],sys.argv[2]
cons=sorted(glob.glob(os.path.join(pt,"consensus sequences","*","*_consensus_strict.fasta")))
fq=glob.glob(os.path.join(pt,"fastq files","*","*.fastq"))
# the key: (the normalised group, the taxid). A1_1 and A1-1 count as the same.
def key(p):
    b=os.path.basename(p)
    m=re.search(r"reads[-_](\d+)",b)
    if not m: return None
    grp=os.path.basename(os.path.dirname(p)).replace("_","-")
    return (grp,m.group(1))
fmap={}
for p in fq:
    k=key(p)
    if k: fmap.setdefault(k,[]).append(p)
rows=[];miss=[]
for c in cons:
    k=key(c)
    got=fmap.get(k,[])
    if len(got)==1: rows.append((k[0],k[1],c,got[0]))
    elif len(got)>1: rows.append((k[0],k[1],c,got[0])); print("WARNING: %d fastq files for %s, the first will be used"%(len(got),str(k)))
    else: miss.append((k,c))
with open(outtsv,"w") as fh:
    fh.write("grup\ttaxid\tkonsensus\tfastq\n")
    for r in rows: fh.write("\t".join(r)+"\n")
print("targets matched: %d"%len(rows))
print("fastq'i olmayan konsensus: %d"%len(miss))
for k,c in miss: print("   EKSIK FASTQ:",k,os.path.basename(c))
PY

N=$(( $(wc -l < "$MAP_TSV") - 1 ))
[ "$N" -gt 0 ] || die "no target matched"
log "targets to process: $N"

# --- 2. an alignment plus two independent measurements per target ------
i=0
tail -n +2 "$MAP_TSV" | while IFS=$'\t' read -r grp taxid cons fastq; do
  i=$(( i + 1 )); tag="${grp}_${taxid}"
  bam="$OUT/bam/${tag}.bam"
  log "[$i/$N] $tag"
  if [ ! -s "$bam" ]; then
    minimap2 -ax map-ont -t "$THREADS" "$cons" "$fastq" 2> "$OUT/log/${tag}_minimap2.log" \
      | samtools sort -@ "$THREADS" -o "$bam" 2>> "$OUT/log/${tag}_minimap2.log"
    samtools index "$bam"
  fi
  # measurement 1: a consensus that allows IUPAC ambiguity codes
  cfg=(); [ -n "$SMCONF" ] && cfg=(-X "$SMCONF")
  samtools consensus -a -A "${cfg[@]}" -o "$OUT/ambig/${tag}_ambig.fa" "$bam" \
    2> "$OUT/log/${tag}_ambig.log"
  # measurement 2: the depth per position, the IUPAC call and the base string.
  # -A burada da verilir; PILEUP bicimindeki 5. sutun boylece IUPAC kodunu tasir
  # and the classification rests on the explicit position number, not on the FASTA length.
  samtools consensus -a -A -f PILEUP "${cfg[@]}" -o "$OUT/pileup/${tag}_pileup.txt" "$bam" \
    2> "$OUT/log/${tag}_pileup.log"
done

# --- 3. N pozisyonlarinin siniflandirilmasi ---------------------------
log "N pozisyonlari siniflandiriliyor"
python3 - "$PT" "$OUT" "$MAP_TSV" "${MINDEPTH:-AUTO}" "$HET" <<'PY'
import sys,os,csv,re,statistics,collections
pt,out,maptsv,mindepth_arg,het=sys.argv[1:6]
het=float(het)
IUPAC2={"R":"AG","Y":"CT","S":"CG","W":"AT","K":"GT","M":"AC"}

def readfa(p):
    s=[]
    for l in open(p):
        if not l.startswith(">"): s.append(l.strip())
    return "".join(s).upper()

targets=list(csv.DictReader(open(maptsv),delimiter="\t"))

# --- deriving the threshold: measure the depth distribution of every target first ---
alldepth=[]
pil={}
for t in targets:
    tag="%s_%s"%(t["grup"],t["taxid"])
    p=os.path.join(out,"pileup","%s_pileup.txt"%tag)
    if not os.path.exists(p): continue
    d={}
    for line in open(p,errors="replace"):
        f=line.rstrip("\n").split("\t")
        if len(f)<7: continue
        try: pos=int(f[1]); dep=int(f[3])
        except ValueError: continue
        d[pos]=(dep,f[4],f[6])   # derinlik, cagri, baz dizisi
        alldepth.append(dep)
    pil[tag]=d
alldepth.sort()
if not alldepth: sys.exit("pileup okunamadi")
# Because zero depth positions pull the distribution down, the threshold is derived
# from the median of the per target median depth. The floor is 5, because below 5
# reads the notion of a majority means nothing.
permed=[]
for tag,d in pil.items():
    dl=[v[0] for v in d.values() if v[0]>0]
    if dl: permed.append(statistics.median(dl))
typ=statistics.median(permed) if permed else 0
med=statistics.median(alldepth)
p10=alldepth[len(alldepth)//10]
if mindepth_arg=="AUTO":
    MIND=max(5,int(round(0.10*typ)))
    derive=("derived from the data: max(5, half the median of the per target "
            "medyaninin yuzde 10'u = %.1f)"%(0.10*typ))
else:
    MIND=int(mindepth_arg); derive="given by hand"
print("derinlik dagilimi (tum pozisyonlar): medyan=%d  yuzde10=%d  min=%d  maks=%d"
      %(med,p10,alldepth[0],alldepth[-1]))
print("the median of the per target median depth: %.1f"%typ)
print("kullanilan min-depth esigi: %d  (%s)"%(MIND,derive))
print("kullanilan het-fract esigi: %.2f"%het)

rows=[]; summary=[]
for t in targets:
    tag="%s_%s"%(t["grup"],t["taxid"])
    strict=readfa(t["konsensus"])
    apath=os.path.join(out,"ambig","%s_ambig.fa"%tag)
    if not os.path.exists(apath) or tag not in pil: continue
    ambig=readfa(apath)
    d=pil[tag]
    # A cross check: do the FASTA and PILEUP calls agree. If the lengths are equal
    # every position is compared and the number of disagreements is logged. The
    # classification rests on PILEUP, because the position number is explicit there.
    xchk_len_ok = (len(ambig)==len(strict))
    xchk_diff=0
    if xchk_len_ok:
        for j in range(len(strict)):
            pcall=d.get(j+1,(0,"?",""))[1]
            if pcall!="?" and ambig[j]!=pcall: xchk_diff+=1
        if xchk_diff: print("WARNING: %s, FASTA and PILEUP disagree at %d positions"%(tag,xchk_diff))
    else:
        print("WARNING: %s, the FASTA length (%d) does not equal strict (%d), the cross "
              "check was skipped; the classification rests on the PILEUP positions"
              %(tag,len(ambig),len(strict)))
    cnt=collections.Counter(); bed=[]
    for i,ch in enumerate(strict):
        if ch!="N": continue
        pos=i+1
        dep,call,bases=d.get(pos,(0,"?",""))
        ab=call                      # PILEUP column 5; with -A it carries an IUPAC code
        b=collections.Counter(x.upper() for x in bases if x.upper() in "ACGT")
        tot=sum(b.values()); maj=b.most_common(1)[0] if tot else ("-",0)
        minf=(b.most_common(2)[1][1]/tot) if len(b)>1 and tot else 0.0
        if dep<MIND:
            cls="dusuk_derinlik"
        elif ab in IUPAC2 and minf>=het:
            cls="iki_allelli"
        elif ab in "ACGT":
            cls="kurtarilabilir"
        else:
            cls="belirsiz"
        cnt[cls]+=1
        rows.append(dict(grup=t["grup"],taxid=t["taxid"],pozisyon=pos,derinlik=dep,
                         ambig_cagri=ab,baskin_baz=maj[0],baskin_oran=round(maj[1]/tot,3) if tot else 0,
                         ikinci_oran=round(minf,3),sinif=cls))
        # maskeye giren siniflar: iki_allelli ve dusuk_derinlik kalici yasak,
        # even a recoverable one stays forbidden for a primer footprint (we do not change the sequence).
        bed.append((pos-1,pos,cls))
    # BED birlestirme
    bedp=os.path.join(out,"maske","%s_maske.bed"%tag)
    with open(bedp,"w") as fh:
        cur=None
        for s,e,c in bed:
            if cur and s==cur[1] and c==cur[2]: cur=(cur[0],e,c)
            else:
                if cur: fh.write("%s\t%d\t%d\t%s\n"%(tag,cur[0],cur[1],cur[2]))
                cur=(s,e,c)
        if cur: fh.write("%s\t%d\t%d\t%s\n"%(tag,cur[0],cur[1],cur[2]))
    dl=[d[p][0] for p in d]
    summary.append(dict(grup=t["grup"],taxid=t["taxid"],uzunluk=len(strict),
                        N=strict.count("N"),
                        medyan_derinlik=int(statistics.median(dl)) if dl else 0,
                        dusuk_derinlik=cnt["dusuk_derinlik"],iki_allelli=cnt["iki_allelli"],
                        kurtarilabilir=cnt["kurtarilabilir"],belirsiz=cnt["belirsiz"],
                        fasta_pileup_ayrilik=(xchk_diff if xchk_len_ok else -1)))

with open(os.path.join(out,"N_pozisyonlari.tsv"),"w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys()) if rows else
                     ["grup","taxid","pozisyon","derinlik","ambig_cagri","baskin_baz",
                      "baskin_oran","ikinci_oran","sinif"],delimiter="\t")
    w.writeheader(); w.writerows(rows)
with open(os.path.join(out,"hedef_ozeti.tsv"),"w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(summary[0].keys()),delimiter="\t")
    w.writeheader(); w.writerows(summary)

tot=collections.Counter(r["sinif"] for r in rows)
print()
print("=== THE TOTAL N CLASSIFICATION ===")
for k in ("dusuk_derinlik","iki_allelli","kurtarilabilir","belirsiz"):
    print("  %-16s %6d"%(k,tot[k]))
print("  %-16s %6d"%("TOTAL",sum(tot.values())))
print()
print("iki_allelli  : real within strain variability, PERMANENTLY forbidden for a primer footprint")
print("dusuk_derinlik: no read support, forbidden under the mask")
print("kurtarilabilir: tek baz acik cogunlukta, dizi DEGISTIRILMEDI ama not edildi")
print()
print("The 10 targets with the most biallelic positions:")
for s in sorted(summary,key=lambda x:-x["iki_allelli"])[:10]:
    print("  %-6s %-8s uzunluk=%5d N=%4d medyan_derinlik=%5d iki_allelli=%4d dusuk=%4d kurtarilabilir=%4d"
          %(s["grup"],s["taxid"],s["uzunluk"],s["N"],s["medyan_derinlik"],
            s["iki_allelli"],s["dusuk_derinlik"],s["kurtarilabilir"]))
PY

log "finished"
echo
echo "Output:"
echo "  $OUT/N_pozisyonlari.tsv   every N position, its depth, the IUPAC call and the class"
echo "  $OUT/hedef_ozeti.tsv      a summary per target"
echo "  $OUT/maske/*.bed          regions where a primer may not be placed"
echo "  $OUT/ambig/*.fa           IUPAC kodlarina izin veren konsensus (referans, kullanilmayacak)"
echo "  $OUT/pileup/*.txt         depth and base string per position"
echo
echo "Share these outputs. The count of biallelic positions says which targets"
echo "primer tasarimina uygun oldugunu birlikte karara baglayacagiz."
