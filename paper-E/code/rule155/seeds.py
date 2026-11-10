import numpy as np, itertools, sys
from strip import ext
def allcycles(beta,K):
    cyc=[[(0,)],[(1,)]]   # line 0 free: both values
    for k in range(1,K):
        new=[]
        for c in cyc: new+=[x[1] for x in ext(c,beta,k)]
        cyc=new
    return cyc
def base(R):
    W=2*R+1; S=set()
    for u in itertools.product((0,1),repeat=W+1):
        s=''.join(map(str,u))
        if '11' in s or '000' in s: continue
        # Phi+: y_i=u[i], x_i=u[i+1]  (i=0..W-1) ; Phi-: y_i=u[i+1], x_i=u[i]
        for (ys,xs) in ((u[:W],u[1:]),(u[1:],u[:W])):
            X=int(''.join(map(str,xs)),2); Y=int(''.join(map(str,ys)),2); S.add((X<<W)|Y)
    return S
def fronts(R,smin,smax,betas,rev,genuine_only=False):
    W=2*R+1; S=set()
    for beta in betas:
        K=max(1,smax+R+1)
        for c in allcycles(beta,K):
            if genuine_only and c[0][0]!=beta: continue
            L=len(c)
            for t in range(L):
                Zt=c[t]; Zp=c[(t-1)%L]
                for s in range(smin,smax+1):
                    xs=[];ys=[]
                    for j in range(s-R,s+R+1):
                        if j<0:
                            d=-j; ys.append(beta if d%2 else 1-beta); xs.append(1-beta if d%2 else beta)
                        else:
                            ys.append(Zt[j]); xs.append(beta if j==0 else Zp[j-1])
                    X=int(''.join(map(str,xs)),2); Y=int(''.join(map(str,ys)),2)
                    S.add((X<<W)|Y)
                    if rev: S.add((Y<<W)|X)
    return S
def load(fn,R):
    W=2*R+1; a=np.unpackbits(np.fromfile(fn,dtype=np.uint8),bitorder='little')[:4**W]
    return set(np.nonzero(a)[0].tolist())
def save(S,R,fn):
    W=2*R+1; Z=np.zeros(max(4**W,64),bool); Z[list(S)]=True
    np.packbits(Z,bitorder='little').tofile(fn)
