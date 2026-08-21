"""A fast pair (F+R) product scan over the pool; the criterion is exactly ispcr's."""
# -------------------------------------------------------------------------
# urunler - finds in which reads of the pool the F+R pair gives a product.
#
# THE CRITERION (exactly ispcr.amplify's):
#   - F binds in the forward direction (hv.bul; the last 2 bases at the 3' end EXACT,
#     mismatches <= mm)
#   - not R itself but its REVERSE COMPLEMENT is searched for (hv.bul5). Since the
#     engine scans only the plus strand, finding rc(R) on the plus strand means R
#     really binds the opposite strand. bul5 requires the FIRST 2 bases of rc(R),
#     which correspond to R's last 2 bases at the 3' end, to match exactly.
#   - The product length window lo..hi. The default (70, 250) comes from the qPCR
#     constraint: in the configuration URUN_IDEAL is 60-150 and URUN_KABUL 150-250 bp.
#   - seg >= f + len(F): R's binding site must lie to the RIGHT of F and NOT OVERLAP
#     it; the two primers must face one another. Otherwise two bindings sitting on the
#     same region would produce a false product.
# -------------------------------------------------------------------------
import sys,pickle,numpy as np
# Oturuma bagli mutlak yollar - bu betikler bir WSL calisma dizininde yazildi.
# engine_gateway.py ice aktarmadan once bu yollari gecerli hale getirir.
sys.path.insert(0,'/tmp/wk2/engine'); sys.path.insert(0,'/tmp/wk2/is')
import ispcr

# -------------------------------------------------------------------------
# urunler - finds in which reads of the pool the F+R pair gives a product.
#
# THE CRITERION (exactly ispcr.amplify's):
#   - F binds in the forward direction (hv.bul; the last 2 bases at the 3' end EXACT,
#     mismatches <= mm)
#   - not R itself but its REVERSE COMPLEMENT is searched for (hv.bul5). Since the
#     engine scans only the plus strand, finding rc(R) on the plus strand means R
#     really binds the opposite strand. bul5 requires the FIRST 2 bases of rc(R),
#     which correspond to R's last 2 bases at the 3' end, to match exactly.
#   - The product length window lo..hi. The default (70, 250) comes from the qPCR
#     constraint: in the configuration URUN_IDEAL is 60-150 and URUN_KABUL 150-250 bp.
#   - seg >= f + len(F): R's binding site must lie to the RIGHT of F and NOT OVERLAP
#     it; the two primers must face one another. Otherwise two bindings sitting on the
#     same region would produce a false product.
# -------------------------------------------------------------------------
def urunler(hv,F,R,lo=70,hi=250,mm=1):
    fs=hv.bul(F,mm)
    # If F binds nowhere there can be no product - an early exit with an empty mask.
    if fs.size==0: return np.zeros(hv.n,bool),{}
    rs=hv.bul5(ispcr.rc(R),mm)
    if rs.size==0: return np.zeros(hv.n,bool),{}
    # sid: WHICH read each global position in the pool belongs to. In the 'N' padding
    # regions between reads, sid = -1.
    fid=hv.sid[fs]; rid=hv.sid[rs]
    # sid < 0 olan vuruslar dolgu bolgesine dusmustur; atilir. Dolgu hem iki
    # okumanin ucunun yanlislikla birlesip sahte bir motif olusturmasini, hem
    # de bir urunun iki ayri okumaya yayilmasini engeller.
    okf=fid>=0; fs=fs[okf]; fid=fid[okf]
    okr=rid>=0; rs=rs[okr]; rid=rid[okr]
    m=np.zeros(hv.n,bool); boy={}
    # The R hits are sorted by read identity and each read's range is extracted in
    # advance with searchsorted. That way the loop below looks, for one F hit, ONLY at
    # the R hits IN THE SAME READ, which is how the product is guaranteed to stay inside
    # a single read. It also takes an O(1) slice instead of scanning the whole R list for
    # every F; the speed comes from there.
    order=np.argsort(rid,kind='stable'); rs=rs[order]; rid=rid[order]
    starts=np.searchsorted(rid,np.arange(hv.n),'left')
    ends=np.searchsorted(rid,np.arange(hv.n),'right')
    for i,f in zip(fid,fs):
        # If the read is already marked as giving a product it is not looked at again: this
        # is a PRESENCE measurement (in how many reads is there a product), not depth.
        if m[i]: continue
        seg=rs[starts[i]:ends[i]]
        if seg.size==0: continue
        # Urun boyu: F'nin basindan rc(R)'nin sonuna; primerler urune dahildir.
        sz=seg+len(R)-f
        sel=(sz>=lo)&(sz<=hi)&(seg>=f+len(F))
        if sel.any():
            m[i]=True
            # ONE value per read goes into the length distribution (the first matching product).
            v=int(sz[sel][0]); boy[v]=boy.get(v,0)+1
    return m,boy
