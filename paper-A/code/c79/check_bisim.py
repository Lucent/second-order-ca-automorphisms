# Confirmation (not used in the proof): the greatest bisimulation G_L on the level-L states of the move graph
# Gamma_L with respect to the edge values. Usage: python3 check_bisim.py L   (L=3: G_3 is the diagonal; L=2: it is not).
import numpy as np, sys
from graphs import build
if __name__=='__main__':
    L=int(sys.argv[1])
    valid,edges,Z,allc,K=build(L,-L,L)
    n=6*K
    print('L',L,'states',int(valid.sum()),'move batches',len(edges),flush=True)
    # move function per (tau,u): dest[id] ; pairs only within same tau with same moves
    mv={}  # u -> dest array over all ids (-1 if none)
    for u,kind,src,dst in edges:
        d=mv.setdefault(u,-np.ones(n,np.int64)); d[src]=dst
    # initial relation per tau: same z key
    key=np.zeros(n,np.int64)
    for j in range(Z.shape[1]): key=key*2+Z[:,j]
    ids=np.nonzero(valid)[0]
    tau_of=ids//K
    R=np.zeros((n,n),bool) if n<=20000 else None
    I,J=np.meshgrid(ids,ids,indexing='ij'); I=I.ravel(); J=J.ravel()
    keep=(I//K==J//K)&(key[I]==key[J])
    I=I[keep]; J=J[keep]
    print('initial pairs',len(I),flush=True)
    while True:
        R=set()
        code=I*n+J
        Rs=np.sort(code)
        ok=np.ones(len(I),bool)
        for u,d in mv.items():
            di=d[I]; dj=d[J]
            has=(di>=0)|(dj>=0)
            both=(di>=0)&(dj>=0)
            assert (has==both).all()
            c=di*n+dj
            pos=np.searchsorted(Rs,c); pos[pos>=len(Rs)]=0
            inR=Rs[pos]==c
            ok&=~both|inR
        rem=int((~ok).sum())
        print('removed',rem,'left',int(ok.sum()),flush=True)
        if rem==0: break
        I=I[ok]; J=J[ok]
    print('diagonal:',bool((I==J).all()),'pairs',len(I))
    if L==3: assert (I==J).all() and len(I)==1152
    if L==2: assert len(I)==1660
    off=I!=J
    if off.any():
        # what differs
        ci=allc[I[off]%K]; cj=allc[J[off]%K]
        diffcols=sorted(set(int(L-c) for a,b in zip(ci,cj) for c in np.nonzero(a!=b)[0]))
        print('offdiag pairs',int(off.sum()),'differ at distances from edge',diffcols)
