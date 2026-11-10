import numpy as np, random
from sim import step, mk
def place_np(X,Y,kk,st,start,p):
    L=2*kk['b']; j=kk['j']
    while (start-(j+L))%p: start+=1
    i=start; Wd=len(X)
    for c in kk['path']:
        if i>=Wd: break
        X[i]=c>>1;Y[i]=c&1;i+=1
    if kk['cyc'] is None:
        cells=[2*st[0][q]+st[1][q] for q in range(p)]; o=(kk['gend']+L)%p
        while i<Wd: c=cells[o%p];X[i]=c>>1;Y[i]=c&1;i+=1;o+=1
    else:
        jj=0
        while i<Wd: c=kk['cyc'][jj%len(kk['cyc'])];X[i]=c>>1;Y[i]=c&1;i+=1;jj+=1
    return start
def legit_gap(rule,kk,st,p,gmax=10,nback=36,trials=6,W=500,r=200):
    """smallest gap g such that d|K(front at r+g) = F^n(F^-n d | F^-n K) for n<=nback (numerically)"""
    ft=mk(rule); rng=random.Random(7)
    for g in range(2,gmax+1):
        X=np.zeros(W,np.int8);Y=np.zeros(W,np.int8)
        for c in range(W): X[c]=st[0][c%p];Y[c]=st[1][c%p]
        place_np(X,Y,kk,st,r+g,p)
        ok=True
        for trial in range(trials):
            Xd=X.copy();Yd=Y.copy()
            for c2 in range(r-30,r+1): Xd[c2]=rng.randrange(2);Yd[c2]=rng.randrange(2)
            Xg,Yg=Xd.copy(),Yd.copy(); Xk,Yk=X.copy(),Y.copy(); Xo,Yo=Xd.copy(),Yd.copy()
            for c2 in range(r+1,W): Xo[c2]=st[0][c2%p];Yo[c2]=st[1][c2%p]
            for jj in range(1,nback):
                Xg,Yg=step(ft,Xg,Yg,-1);Xk,Yk=step(ft,Xk,Yk,-1);Xo,Yo=step(ft,Xo,Yo,-1)
                lo=r-30+jj+1
                if not (np.array_equal(Xg[lo:r+1],Xo[lo:r+1]) and np.array_equal(Yg[lo:r+1],Yo[lo:r+1]) and np.array_equal(Xg[r+1:W-jj-1],Xk[r+1:W-jj-1]) and np.array_equal(Yg[r+1:W-jj-1],Yk[r+1:W-jj-1])):
                    ok=False; break
            if not ok: break
        if ok: return g
    return None
