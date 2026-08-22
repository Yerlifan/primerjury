"""A fast STRICT scanner: a 3' terminal 10-mer seed (the last 2 bases exact) plus a
full primer verification.
The criterion is the same as ispcr.find_sites: mismatches <=1, the last 2 bases at
the 3' end exact.

"""
import sys, numpy as np
sys.path.insert(0,'/tmp/wk2/engine')
import ispcr
SEED=10

class Havuz:
    def __init__(self, seqs, sep=40):
        self.n=len(seqs)
        parts=[]; sid=[]; off=[]
        cur=0
        for i,s in enumerate(seqs):
            parts.append(s); sid.append(np.full(len(s),i,dtype=np.int32))
            off.append(cur); cur+=len(s)
            parts.append('N'*sep); sid.append(np.full(sep,-1,dtype=np.int32)); cur+=sep
        self.big=''.join(parts)
        self.sid=np.concatenate(sid)
        self.enc=ispcr.encode(self.big)
        e=self.enc.astype(np.int64)
        L=len(e); m=L-SEED+1
        code=np.zeros(m,dtype=np.int64); bad=np.zeros(m,dtype=bool)
        for k in range(SEED):
            c=e[k:k+m]
            code=code*4+np.where(c<0,0,c)
            bad|=(c<0)
        code[bad]=-1
        ok=np.nonzero(code>=0)[0]
        o=np.argsort(code[ok],kind='stable')
        self.pos=ok[o].astype(np.int32)
        self.key=code[ok][o]

    def seed_pos(self,kmer):
        c=0
        for ch in kmer:
            v='ACGT'.find(ch)
            if v<0: return np.empty(0,dtype=np.int32)
            c=c*4+v
        a=np.searchsorted(self.key,c,'left'); b=np.searchsorted(self.key,c,'right')
        return self.pos[a:b]

    def bul(self,primer,max_mm=1):
        """primerin ileri yonde bagladigi (global baslangic) dizisi."""
        L=len(primer); tail=primer[-SEED:]
        vars=set([tail])
        for i in range(SEED-2):           # son 2 baz TAM kalmali
            for b in 'ACGT':
                if b!=tail[i]: vars.add(tail[:i]+b+tail[i+1:])
        cand=[]
        for v in vars:
            p=self.seed_pos(v)
            if p.size: cand.append(p)
        if not cand: return np.empty(0,dtype=np.int32)
        st=np.unique(np.concatenate(cand))-(L-SEED)
        st=st[(st>=0)&(st+L<=len(self.enc))]
        if st.size==0: return st
        mm=np.zeros(st.size,dtype=np.int16)
        pe=np.array(['ACGT'.find(c) for c in primer],dtype=np.int8)
        for k in range(L):
            col=self.enc[st+k]
            mm+=((col!=pe[k])|(col<0))
        return st[mm<=max_mm]

    def dizi_seti(self,primer,max_mm=1):
        st=self.bul(primer,max_mm)
        if st.size==0: return np.zeros(self.n,dtype=bool), {}
        ids=self.sid[st]
        m=np.zeros(self.n,dtype=bool); m[ids[ids>=0]]=True
        return m, dict(zip(st.tolist(), ids.tolist()))

    def bul5(self,rcrev,max_mm=1):
        """For rc(R): the last 2 bases at the 5' end (index 0, 1) EXACT, mismatches <=1."""
        L=len(rcrev); head=rcrev[:SEED]
        vars=set([head])
        for i in range(2,SEED):
            for b in 'ACGT':
                if b!=head[i]: vars.add(head[:i]+b+head[i+1:])
        cand=[]
        for v in vars:
            p=self.seed_pos(v)
            if p.size: cand.append(p)
        if not cand: return np.empty(0,dtype=np.int32)
        st=np.unique(np.concatenate(cand)).astype(np.int64)
        st=st[(st>=0)&(st+L<=len(self.enc))]
        if st.size==0: return st.astype(np.int32)
        mm=np.zeros(st.size,dtype=np.int16)
        pe=np.array(['ACGT'.find(c) for c in rcrev],dtype=np.int8)
        for k in range(L):
            col=self.enc[st+k]
            mm+=((col!=pe[k])|(col<0))
        return st[mm<=max_mm].astype(np.int32)
