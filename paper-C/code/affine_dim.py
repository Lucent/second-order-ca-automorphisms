"""Dimension of the affine F2-space of radius-R block maps phi with phi∘F = F∘phi, for an AFFINE
second-order rule f (cells c = 2x+y, F(x,y)_i = (y_i, x_i ^ f(y_{i-1},y_i,y_{i+1}))).
Unknowns: X[w], Y[w] = bits of phi on window w (4^(2R+1) windows).
For each window c of length 2R+3:  X[F(c)] = Y[c_mid]   (2-term, merged by union-find)
                                    Y[F(c)] = X[c_mid] + a*Y[c_left] + b*Y[c_mid] + g*Y[c_right] + eps
Independent of centralizer.py and of the SAT encoding. Usage: affine_dim.py RULE R1 [R2 ...]."""
import sys, time, numpy as np
def run(rule,R):
    t0=time.time()
    f=lambda a,b,c:(rule>>((a<<2)|(b<<1)|c))&1
    e0=f(0,0,0); al=f(1,0,0)^e0; be=f(0,1,0)^e0; ga=f(0,0,1)^e0
    for a in (0,1):
        for b in (0,1):
            for c in (0,1): assert f(a,b,c)==(al&a)^(be&b)^(ga&c)^e0
    n=2*R+1; m=n+2; N=4**n; W=4**m
    w=np.arange(W,dtype=np.int64)
    cell=[(w>>(2*(m-1-i)))&3 for i in range(m)]
    # F on interior cells 1..m-2
    def Floc(l,c,r):
        x=c>>1; y=c&1
        fy=(al&(l&1))^(be&y)^(ga&(r&1))^e0
        return (y<<1)|(x^fy)
    Fc=[Floc(cell[i-1],cell[i],cell[i+1]) for i in range(1,m-1)]
    imF=np.zeros(W,dtype=np.int64)
    for v in Fc: imF=(imF<<2)|v
    def sub(lo):
        r=np.zeros(W,dtype=np.int64)
        for i in range(lo,lo+n): r=(r<<2)|cell[i]
        return r
    L,Cm,Rt=sub(0),sub(1),sub(2)
    X=lambda e:2*e; Y=lambda e:2*e+1
    nv=2*N; parent=np.arange(nv)
    def find(a):
        r=a
        while parent[r]!=r: r=parent[r]
        while parent[a]!=r: parent[a],a=r,parent[a]
        return r
    for a,b in set(zip(X(imF).tolist(),Y(Cm).tolist())):
        ra,rb=find(a),find(b)
        if ra!=rb: parent[max(ra,rb)]=min(ra,rb)
    root=np.array([find(i) for i in range(nv)])
    roots,cls=np.unique(root,return_inverse=True); C=len(roots)
    cols=[cls[Y(imF)],cls[X(Cm)]]
    if al: cols.append(cls[Y(L)])
    if be: cols.append(cls[Y(Cm)])
    if ga: cols.append(cls[Y(Rt)])
    words=(C+1+63)//64   # last column C = constant
    rows=np.zeros((W,words),dtype=np.uint64)
    ar=np.arange(W)
    for col in cols:
        np.bitwise_xor.at(rows,(ar,col//64),np.uint64(1)<<(col%64).astype(np.uint64))
    if e0:
        rows[:,C//64]^=np.uint64(1)<<np.uint64(C%64)
    rows=np.unique(rows,axis=0); rows=rows[rows.any(axis=1)]
    rank=0; used=np.zeros(len(rows),bool); inconsistent=False
    for col in range(C+1):
        wd,b=col//64,np.uint64(col%64)
        has=((rows[:,wd]>>b)&np.uint64(1)).astype(bool)
        cand=np.flatnonzero(has&~used)
        if len(cand)==0: continue
        if col==C: inconsistent=True; break
        p=cand[0]; used[p]=True; rank+=1; has[p]=False
        t=np.flatnonzero(has)
        if len(t): rows[t]^=rows[p]
    print(f"rule {rule} R={R}: {C} classes, rank {rank}, consistent {not inconsistent}, solution-space dim {C-rank} [{time.time()-t0:.0f}s]",flush=True)
if __name__=="__main__":
    for R in map(int,sys.argv[2:]): run(int(sys.argv[1]),R)
