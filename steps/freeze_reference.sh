#!/usr/bin/env bash
# =====================================================================
# freeze_reference.sh
# Amaç: 06'nın ürettiği konsensüsleri, aynı taksonun bütün yılları TEK ve
#       AYNI koordinat sisteminde olacak şekilde yeniden kurmak.
#
# İki kip var, ikisi ayrı klasöre yazar, böylece yan yana karşılaştırılabilir.
#
#   --mode pin   (öntanımlı)
#     Dış referansı sabitler. 06 her takson-barkod için bağımsız olarak en
#     iyi referansı seçiyor ve bazı taksonlarda yıllar farklı referansa
#     oturuyor. Bu kip, takson başına toplam bit skoru en yüksek referansı
#     kazanan ilan edip o taksonun bütün yıllarını ona hizalar.
#     Avantajı: koordinatlar dış dünyayla karşılaştırılabilir kalır.
#
#   --mode self
#     Kendi verinizi referans alır. Her takson için 06 çıktısındaki en iyi
#     kapsanan konsensüs seçilir, kapsanmayan uçları kırpılır, o taksonun
#     bütün yılları bu diziye hizalanır.
#     Gerekçe: ROD yalnızca 11.935 genom kapsıyor ve hedef taksonun kendisi
#     yoksa blastn en yakın akrabayı seçiyor. Gerçek koşuda taxid 4896
#     (Schizosaccharomyces pombe, bir Ascomycota mayası) Earliella scabrosa'ya
#     yani bir Basidiomycota mantarına, taxid 44689 (Dictyostelium, Amoebozoa)
#     Pleurotus giganteus'a oturdu. Dizinin kendisi okumalardan çağrılıyor,
#     ama yabancı bir iskelet indellerin ve kenarların yorumunu bozar.
#     Kendi verinizi referans almak bu sorunu tümüyle ortadan kaldırır.
#     Dezavantajı: koordinatlar yalnızca bu çalışmaya özgüdür.
#
# Her iki kipte de konsensüs KATI referans koordinatında çağrılır:
#   --show-ins no --show-del yes  ile çıktı uzunluğu referansa eşitlenir
#   (sentetik BAM ile doğrulandı: 200 bp referans, 6 baz ekleme ve 5 baz
#   silme taşıyan okumalar; bayraksız 201, yalnız ins-no ile 195, yalnız
#   del-yes ile 206, ikisi birlikte tam 200), silme işareti * sonradan N'e
#   çevrilir çünkü silinmiş bir pozisyon baz değildir ve orada primer ayağı
#   olamaz.
#
# Kullanım:
#   bash freeze_reference.sh --pt <PrimerTasarlama> --out <06 cikti> \
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
    *) echo "bilinmeyen secenek: $1" >&2; exit 2;;
  esac
done
[ -n "$PT" ] && [ -n "$OUT" ] || {
  echo "kullanim: bash $0 --pt <PrimerTasarlama> --out <06 cikti klasoru> [--mode pin|self]" >&2
  exit 2; }
case "$MODE" in pin|self) ;; *) echo "HATA: --mode pin ya da self olmali" >&2; exit 2;; esac
case "$PT$OUT" in *'...'*) echo "ERROR: the path contains '...'. Give a full path, not a placeholder." >&2; exit 2;; esac
[ -d "$PT" ] || { echo "ERROR: no such --pt directory: $PT" >&2; exit 1; }
MAP="$OUT/referans_secimi.tsv"
[ -s "$MAP" ] || { echo "ERROR: $MAP is missing. Run the mapping step first." >&2; exit 1; }

log(){ printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die(){ printf 'HATA: %s\n' "$*" >&2; exit 1; }
for t in minimap2 samtools python3; do
  command -v "$t" >/dev/null 2>&1 || die "$t PATH'te yok"
done
[ -z "$THREADS" ] && THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))
log "kip: $MODE, is parcacigi: $THREADS"

konsensus_cagir() {   # $1=bam  $2=cikti fasta  $3=cikti pileup  $4=log oneki
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
print("%-9s %-6s %-16s %-16s %s"%("taxid","barkod","eski","yeni","durum"))
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
        st="ayni" if r["referans_id"]==win else "DEGISECEK"
        if st=="DEGISECEK": degisen+=1
        print("%-9s %-6s %-16s %-16s %s"%(tid,r["grup"],
              r["referans_id"].split("|")[0][:16],win.split("|")[0][:16],st))
        plan.append((r["grup"],tid,r["referans_id"],win,st))
with open(planf,"w",newline="",encoding="utf-8") as fh:
    w=csv.writer(fh,delimiter="\t",lineterminator="\n")
    w.writerow(["grup","taxid","eski_referans","yeni_referans","durum"]); w.writerows(plan)
print()
print("cok referansli takson: %d"%sum(1 for t,rs in bytax.items()
      if len(set(r["referans_id"] for r in rs))>1))
print("yeniden hizalanacak hedef: %d"%degisen)
PY
[ "$DRY" = 1 ] && { log "kuru kosu. Plan: $PLAN"; exit 0; }

n=0
while IFS=$'\t' read -r grp tid eski yeni durum; do
  [ "$durum" = "DEGISECEK" ] || continue
  tag="${grp}_${tid}"; src=""
  for cand in "$OUT/ref"/*_ref.fasta; do
    head -1 "$cand" | grep -qF "$yeni" && { src="$cand"; break; }
  done
  [ -n "$src" ] || { echo "  the winning reference file was not found ($tag)" >&2; continue; }
  fq=$(fastq_bul "$grp" "$tid")
  [ -n "$fq" ] || { echo "  no fastq found: $grp $tid" >&2; continue; }
  n=$((n+1)); log "[$n] $tag yeniden hizalaniyor"
  cp -f "$src" "$OUT/ref/${tag}_ref.fasta"
  rm -f "$OUT/bam/${tag}.bam" "$OUT/bam/${tag}.bam.bai"
  minimap2 -ax map-ont -t "$THREADS" "$OUT/ref/${tag}_ref.fasta" "$fq" \
      2> "$OUT/log/${tag}_mm2_sabit.log" \
    | samtools sort -@ "$THREADS" -o "$OUT/bam/${tag}.bam" 2>> "$OUT/log/${tag}_mm2_sabit.log"
  samtools index "$OUT/bam/${tag}.bam"
  konsensus_cagir "$OUT/bam/${tag}.bam" "$OUT/konsensus/${tag}_ref_konsensus.fasta" \
                  "$OUT/pileup/${tag}_pileup.txt" "$OUT/log/${tag}_sabit"
done < <(tail -n +2 "$PLAN")
log "yeniden hizalanan hedef: $n"
BEKLENEN=$(awk -F'\t' 'NR>1 && $5=="DEGISECEK"' "$PLAN" | wc -l)
[ "$n" -eq "$BEKLENEN" ] || echo "UYARI: plan $BEKLENEN diyordu, $n islendi." >&2
DENETIM_DIZIN="$OUT/konsensus"

# =====================================================================
else   # MODE = self
# =====================================================================
SELF="$OUT/self"
mkdir -p "$SELF"/{ref,bam,konsensus,pileup,log}
SECIM="$SELF/self_referans_secimi.tsv"

log "her takson icin en iyi kapsanan konsensus seciliyor"
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
print("%-9s %-4s %-26s %10s %11s"%("taxid","uye","secilen","kapsanan","medyan_der"))
for tid,items in sorted(cand.items(),key=lambda x:int(x[0])):
    # en cok kapsanan kazanir; esitlikte daha derin, o da esitse ad sirasi
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
print("takson sayisi: %d"%len(rows))
PY
[ "$DRY" = 1 ] && { log "kuru kosu. Secim: $SECIM"; exit 0; }

log "self referanslar kirpiliyor ve butun yillar yeniden hizalaniyor"
n=0
while IFS=$'\t' read -r tid stag spath kap med uyeler; do
  ref="$SELF/ref/${tid}_self_ref.fasta"
  # kapsanmayan uclar kirpilir, ic N korunur (bilgi uydurmuyoruz)
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
    sys.exit("  taxid %s: kirpilmis referans cok kisa (%d bp), atlandi"%(tid,len(core)))
open(dst,"w").write(">self_%s kaynak=%s kirpma=%d-%d uzunluk=%d\n"
                    %(tid,stag,i+1,j,len(core))
                    +"\n".join(core[k:k+70] for k in range(0,len(core),70))+"\n")
print("  taxid %-9s %d bp -> kirpilmis %d bp (ic N: %d)"%(tid,len(s),len(core),core.count("N")))
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
log "yeniden hizalanan hedef: $n"
DENETIM_DIZIN="$SELF/konsensus"
fi

# =====================================================================
# son denetim, iki kip icin ortak
# =====================================================================
log "son denetim: $DENETIM_DIZIN"
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
print("%-9s %-4s %-16s %-24s %s"%("taxid","yil","uzunluk","kapsanan","hizali mi"))
for tid,items in sorted(bytax.items(),key=lambda x:int(x[0])):
    L=set(len(s) for _,_,s in items)
    H=set(h.split()[0] for _,h,_ in items)
    ok=len(L)==1 and len(H)==1
    if not ok: bad+=1
    print("%-9s %-4d %-16s %-24s %s"%(tid,len(items),
          ",".join(str(x) for x in sorted(L)),
          ",".join(str(sum(1 for c in s if c!="N")) for _,_,s in items),
          "EVET" if ok else "HAYIR"))
print()
print("hizalanmamis takson: %d"%bad)
print("Butun taksonlarda yillar ayni koordinat sisteminde." if bad==0 else
      "HAYIR yazan taksonlarda koordinat hala ayrisiyor.")
base=os.path.join(out,"konsensus")
if os.path.abspath(base)!=os.path.abspath(d):
    print()
    print("### 06 (dis referans) ile bu kipin karsilastirmasi")
    print("%-26s %12s %12s %8s %10s"%("hedef","06_kapsanan","yeni_kapsanan","fark","yeni_uzunluk"))
    tot0=tot1=0
    for tag,(h,s) in sorted(seqs.items()):
        p0=os.path.join(base,"%s_ref_konsensus.fasta"%tag)
        if not os.path.exists(p0): continue
        s0="".join(l.strip() for l in open(p0) if not l.startswith(">")).upper()
        c0=sum(1 for c in s0 if c!="N"); c1=sum(1 for c in s if c!="N")
        tot0+=c0; tot1+=c1
        print("%-26s %12d %12d %+8d %10d"%(tag,c0,c1,c1-c0,len(s)))
    print("%-26s %12d %12d %+8d"%("TOPLAM",tot0,tot1,tot1-tot0))
    print()
    print("Kapsanan baz sayisi primer tasarimina acik alanin olcusudur.")
    print("Self kipte referans zaten kendi kapsanan bolgemiz oldugu icin")
    print("kapsanan oranin yuksek olmasi beklenir; asil kazanc ise iskeletin")
    print("taksonomik olarak dogru olmasidir.")
PY
log "bitti"
echo
echo "Next step: analyze_ambiguous_bases.sh processes these consensuses into a mask,"
echo "after which design_group_primers.py and split_clusters.py take design back up."
