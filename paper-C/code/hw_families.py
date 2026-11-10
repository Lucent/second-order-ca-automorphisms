"""The families with a hidden involution at the middle position of a window of length 5 (as in hw_search.py),
with the involution itself. Run as `hw_families.py RULES` it checks that the flip of all walls, as a table
of radius 2, commutes with F on all 4^7 windows. Used by hw_flipall.py.
"""
import sys, numpy as np

sys.argv=[sys.argv[0]]+sys.argv[1:]
from centralizer import make_rule, index_word
from classify_r4 import check_commutes_np_R
import hw_search as S
k=5
def families(rule,cap=200):
    f,Fl,Fil=make_rule(rule)
    # reuse analyse but return iotas
    Ft=np.array([[[Fl(a,b,c) for c in range(4)] for b in range(4)] for a in range(4)])
    Fit=np.array([[[Fil(a,b,c) for c in range(4)] for b in range(4)] for a in range(4)])
    n=4**k; idx=np.arange(n)
    W=np.stack([(idx>>(2*(k-1-i)))&3 for i in range(k)],1)
    succ=[]
    for T in (Ft,Fit):
        for a in range(4):
            for b in range(4):
                U=np.concatenate([np.full((n,1),a),W,np.full((n,1),b)],1)
                V=T[U[:,:-2],U[:,1:-1],U[:,2:]]
                succ.append(sum(V[:,i]<<(2*(k-1-i)) for i in range(k)))
    succ=np.stack(succ,1)
    M=np.zeros((n,n),np.float32); M[np.repeat(idx,32),succ.ravel()]=1; M[idx,idx]=1
    for _ in range(12): M=((M@M)>0).astype(np.float32)
    R=M>0; sizes=R.sum(1); seen=set(); res=[]
    j=2; sh=2*(k-1-j)
    for w0 in np.nonzero(sizes<=cap)[0]:
        fam=frozenset(np.nonzero(R[w0])[0].tolist())
        if fam in seen: continue
        seen.add(fam)
        for d in (1,2,3):
            p0=w0^(d<<sh)
            if p0 not in fam: continue
            iota={w0:p0}; stack=[w0]; ok=True
            while stack and ok:
                w=stack.pop(); v=iota[w]
                for e in range(32):
                    a,b=succ[w,e],succ[v,e]
                    if a in iota:
                        if iota[a]!=b: ok=False;break
                    else:
                        if ((a^b)&~(3<<sh)) or b not in fam: ok=False;break
                        iota[a]=b; stack.append(a)
            if ok and len(iota)==len(fam) and all(iota[iota[w]]==w and iota[w]!=w for w in iota):
                res.append(iota); break
    return res, Fl
for rule in (map(int,sys.argv[1].split(",")) if __name__=="__main__" else []):
    res,Fl=families(rule)
    for iota in res:
        tab=[index_word(e,5)[2] for e in range(1024)]
        for w,v in iota.items(): tab[w]=index_word(v,5)[2]
        ok=check_commutes_np_R(tab,Fl,2)
        print(rule,len(iota),'commutes:',ok, sorted(''.join(map(str,index_word(w,5))) for w in iota)[:8])
