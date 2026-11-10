# Direct (t,i) simulator of second-order ECA F(x,y)_i=(y_i, x_i ^ f(y_{i-1},y_i,y_{i+1})).
# Batched over cases. Rows t=-1,0 are given by a placed pair with corner [0,0] (t0=0,i0=0):
#   columns i<0 carry pi(t-i) (Phi+ point), columns i>=0 carry xi(t+i) (Phi- point).
import numpy as np
def ftab(rule): return np.array([(rule>>k)&1 for k in range(8)],np.uint8)
NSUP={94:[2],91:[0],95:[0,2]}  # support of N where f=a^c^N
def ext(seq,tail,n):
    """seq: list for indices 0..len-1 ; tail: tuple (v_even,v_odd) by absolute index parity; return array length n"""
    out=np.empty(n,np.uint8)
    for k in range(n): out[k]= seq[k] if k<len(seq) else tail[k%2]
    return out
def run(rule,PIs,XIs,T,neg_bg,ray_q=None,ray_p=None):
    """PIs: list of arrays pi(q) for q in [0,Lq) (q<0 -> pi uses neg_bg), XIs: arrays xi(p) for p in [0,Lp).
       neg_bg: (b0,b1) background for negative indices by parity (front background beta).
       Returns gamma (cases x P), o (cases x Q), firing list per case."""
    nc=len(PIs); W=4*T+8; off=2*T+4   # column index i -> array idx i+off
    tbl=ftab(rule); Nt=np.zeros(8,np.uint8)
    for k in NSUP[rule]: Nt[k]=1
    I=np.arange(W)-off
    def val(arrs,idx):
        out=np.empty((nc,W),np.uint8)
        for c in range(nc):
            a=arrs[c]; L=len(a)
            ii=np.clip(idx,0,L-1)
            v=a[ii].copy()
            neg=idx<0; v[neg]=np.array(neg_bg,np.uint8)[np.mod(idx[neg],2)]
            out[c]=v
        return out
    rows=[]
    for t in (-1,0):
        q=t-I; p=t+I
        Pv=val(PIs,q); Xv=val(XIs,p)
        rows.append(np.where(I<0,Pv,Xv).astype(np.uint8))
    prev,cur=rows
    fires=[[] for _ in range(nc)]
    for t in range(0,T):
        l=np.roll(cur,1,axis=1); r=np.roll(cur,-1,axis=1)
        idx=4*l+2*cur+r
        nxt=prev^tbl[idx]
        n=Nt[idx]
        cs,js=np.nonzero(n[:,2:-2]); js=js+2
        for c,j in zip(cs,js):
            i=j-off
            if abs(i)<=T: fires[c].append((t-i,t+i))
        prev,cur=cur,nxt
    # cur is row T
    return cur,off,fires
def outputs(cur,off,T,P,Q):
    # gamma(p)=Y(T,p-T), o(q)=Y(T,T-q)
    g=np.stack([cur[:,p-T+off] for p in range(P)],1)
    o=np.stack([cur[:,T-q+off] for q in range(Q)],1)
    return g,o
