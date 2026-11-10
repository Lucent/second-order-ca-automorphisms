# Theorem 2.4, parts (1) and (2), by exhaustion for all 256 rules and n <= NMAX (argument, default 4).
import numpy as np, sys, itertools
def fb(r,a,b,c): return (r>>(4*a+2*b+c))&1
def reads_c(r): return any(fb(r,a,b,0)!=fb(r,a,b,1) for a in (0,1) for b in (0,1))
def reads_a(r): return any(fb(r,0,b,c)!=fb(r,1,b,c) for b in (0,1) for c in (0,1))
def evolve(r,X,Yr,n):
    """X=row -1, Yr=row 0 (arrays shape (M,L)); returns rows dict t-> array, shrinking by 1 each side"""
    ft=np.array([(r>>k)&1 for k in range(8)],dtype=np.int8)
    rows={-1:X,0:Yr}
    for t in range(1,n+1):
        a=rows[t-1]; pp=rows[t-2]
        # align: rows[t-1] has width L-2(t-1) centred; rows[t-2] width L-2(t-2) (row -1 has width L)
        off=(pp.shape[1]-a.shape[1])//2
        new=pp[:,off+1:off+a.shape[1]-1]^ft[a[:,:-2]*4+a[:,1:-1]*2+a[:,2:]]
        rows[t]=new
    return rows
# (1) exhaustive reach of F^n at cell 0 (y and x output) for n<=4, all 256 rules
def reach(r,n):
    L=2*n+1; N=4**L
    idx=np.arange(N,dtype=np.int64)
    cells=np.stack([(idx>>(2*(L-1-i)))&3 for i in range(L)],axis=1).astype(np.int8)
    X=(cells>>1)&1; Yr=cells&1
    rows=evolve(r,X,Yr,n)
    outy=rows[n][:,0]; outx=rows[n-1][:,(rows[n-1].shape[1])//2]
    out=outx*2+outy
    # dependence on cell j (index j+n)
    dep=[]
    for j in range(-n,n+1):
        sh=2*(L-1-(j+n))
        d=False
        for m in (1,2,3):
            if np.any(out!=out[idx^(m<<sh)]): d=True;break
        dep.append(d)
    js=[j for j,d in zip(range(-n,n+1),dep) if d]
    return (max(js),-min(js)) if js else (None,None)
if __name__=="__main__":
    bad=[]
    nmax=int(sys.argv[1]) if len(sys.argv)>1 else 4
    cnt=0
    for r in range(256):
        rc,ra=reads_c(r),reads_a(r)
        for n in range(1,nmax+1):
            rr,ll=reach(r,n)
            if rc and rr!=n: bad.append(('right',r,n,rr))
            if ra and ll!=n: bad.append(('left',r,n,ll))
        cnt+= rc and ra
    print('rules reading a and c:',cnt,'violations:',bad)
