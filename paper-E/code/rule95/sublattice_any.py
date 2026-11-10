"""Block maps of sublattice type commuting with F_f, for an elementary rule f that does not read its middle argument.
A second encoder for the search of sat_sublattice.py; all solutions of both types are enumerated for R = 1..Rmax,
classified as sigma^a F^b, constant or other, and rechecked against hF = Fh on random periodic configurations.
type par=0: y'_i=TY[b], b_m=y_{i+m} (m even), x_{i+m} (m odd); x'_i=TX[b], b_m=x_{i+m} (m even), y_{i+m} (m odd)
type par=1: even/odd exchanged.   usage: sublattice_any.py RULE Rmax [cap [SOLVER]]     (default solver cadical153)"""
import sys, itertools, time, numpy as np
from pysat.solvers import Solver
rule=int(sys.argv[1]); Rmax=int(sys.argv[2]); cap=int(sys.argv[3]) if len(sys.argv)>3 else 20000
solver_name=sys.argv[4] if len(sys.argv)>4 else 'cadical153'
assert all(((rule>>(4*a+c))&1)==((rule>>(4*a+2+c))&1) for a in (0,1) for c in (0,1)), "rule reads middle"
FT=np.array([[(rule>>(4*a+c))&1 for c in (0,1)] for a in (0,1)],dtype=np.int64)
def Fw(X,Y,lo,hi):   # dicts cell->array on lo..hi ; returns on lo+1..hi-1
    return {k:Y[k] for k in range(lo+1,hi)}, {k:X[k]^FT[Y[k-1],Y[k+1]] for k in range(lo+1,hi)}
def Finvw(X,Y,lo,hi):
    return {k:Y[k]^FT[X[k-1],X[k+1]] for k in range(lo+1,hi)}, {k:X[k] for k in range(lo+1,hi)}
def idx(R,par,X,Y,i,out):
    w=0
    for m in range(-R,R+1):
        same=(m-par)%2==0
        lay=(Y if same else X) if out=='y' else (X if same else Y)
        w=w|(lay[i+m]<<(m+R))
    return w
def assignments(cl,lo,hi,other):
    n=hi-lo+1; a=np.arange(1<<n,dtype=np.int64); X={};Y={}
    for j,k in enumerate(range(lo,hi+1)):
        b=(a>>j)&1; o=np.full(1<<n,other,dtype=np.int64)
        if (k%2==0)==(cl==0): Y[k],X[k]=b,o
        else: X[k],Y[k]=b,o
    return X,Y
def named(R,par):
    """tables (TX,TY) of sigma^a F^b, |a|+|b|<=R, that are of this type"""
    out={}
    N=1<<(2*R+1)
    for b in range(-R,R+1):
        for a in range(-(R-abs(b)),R-abs(b)+1):
            tabs={}
            okk=True
            for ob in ('x','y'):
                res=[]
                for other in (0,1):
                    # build config on -R..R from window index of output bit ob at cell 0
                    w=np.arange(N,dtype=np.int64); X={};Y={}
                    for m in range(-R,R+1):
                        bit=(w>>(m+R))&1; o=np.full(N,other,dtype=np.int64)
                        same=(m-par)%2==0
                        inY = same if ob=='y' else (not same)
                        if inY: Y[m],X[m]=bit,o
                        else: X[m],Y[m]=bit,o
                    lo,hi=-R,R
                    for _ in range(abs(b)):
                        X,Y=(Fw if b>0 else Finvw)(X,Y,lo,hi); lo+=1;hi-=1
                    res.append((X if ob=='x' else Y)[a])
                if not (res[0]==res[1]).all(): okk=False;break
                tabs[ob]=res[0].astype(np.uint8)
            if okk: out[(a,b)]=(tabs['x'].tobytes(),tabs['y'].tobytes())
    return out
def run(R,par):
    N=1<<(2*R+1); vx=lambda w:1+w; vy=lambda w:1+N+w
    eqs=set(); rels=set()
    for cl in (0,1):
        for other in (0,1):
            X,Y=assignments(cl,-R-1,R+1,other)
            FX,FY=Fw(X,Y,-R-1,R+1)
            a=idx(R,par,FX,FY,0,'x'); b=idx(R,par,X,Y,0,'y')
            eqs|=set(zip(vx(a).tolist(),vy(b).tolist()))
            a=idx(R,par,FX,FY,0,'y'); b=idx(R,par,X,Y,0,'x'); c=idx(R,par,X,Y,-1,'y'); d=idx(R,par,X,Y,1,'y')
            rels|=set(zip(vy(a).tolist(),vx(b).tolist(),vy(c).tolist(),vy(d).tolist()))
    cnf=[]
    for a,b in eqs:
        if a!=b: cnf+=[[a,-b],[-a,b]]
    for t in rels:
        vs=sorted(set(t))
        for vals in itertools.product((0,1),repeat=len(vs)):
            d=dict(zip(vs,vals))
            if d[t[0]]!=(d[t[1]]^int(FT[d[t[2]],d[t[3]]])):
                cnf.append([-v if d[v] else v for v in vs])
    s=Solver(name=solver_name,bootstrap_with=cnf)
    sols=[]
    while len(sols)<cap and s.solve():
        m=s.get_model(); bits=np.array([1 if v>0 else 0 for v in m[:2*N]],dtype=np.uint8)
        sols.append((bits[:N].tobytes(),bits[N:].tobytes()))
        s.add_clause([-v for v in m[:2*N]])
    capped = len(sols)>=cap
    s.delete()
    return sols,len(cnf),capped
def apply_map(R,par,TX,TY,c):   # c: list of (x,y) cyclic, returns image
    n=len(c); TX=np.frombuffer(TX,dtype=np.uint8); TY=np.frombuffer(TY,dtype=np.uint8); out=[]
    for i in range(n):
        wy=wx=0
        for m in range(-R,R+1):
            x,y=c[(i+m)%n]; same=(m-par)%2==0
            wy|=((y if same else x)<<(m+R)); wx|=((x if same else y)<<(m+R))
        out.append((int(TX[wx]),int(TY[wy])))
    return tuple(out)
def Fc(c):
    n=len(c); return tuple((c[i][1], c[i][0]^int(FT[c[i-1][1],c[(i+1)%n][1]])) for i in range(n))
import random
random.seed(1)
for R in range(1,Rmax+1):
    for par in (0,1):
        t0=time.time()
        sols,ncl,capped=run(R,par)
        nm=named(R,par); nmset=set(nm.values())
        N=1<<(2*R+1)
        const=[s for s in sols if len(set(s[0]))==1 and len(set(s[1]))==1]
        others=[s for s in sols if s not in nmset and s not in const]
        missing=[k for k,v in nm.items() if v not in set(sols)]
        # random commutation test of every 'other' (and a few named)
        bad=0
        for s in others[:200]+[v for v in list(nmset)[:5]]:
            for _ in range(20):
                n=2*R+6; c=tuple((random.randint(0,1),random.randint(0,1)) for _ in range(n))
                if apply_map(R,par,s[0],s[1],Fc(c))!=Fc(apply_map(R,par,s[0],s[1],c)): bad+=1
        # injectivity necessary test for others: periodic configs of period n<=7
        inj_cand=0
        for s in others[:3000]:
            ok=True
            for n in (1,2,3,4,5,6,7):
                seen=set()
                for cc in itertools.product(range(4),repeat=n):
                    c=tuple((v>>1,v&1) for v in cc)
                    im=apply_map(R,par,s[0],s[1],c)
                    if im in seen: ok=False;break
                    seen.add(im)
                if not ok: break
            if ok: inj_cand+=1
        print(f"rule {rule} R={R} type={par}: clauses={ncl} solutions={len(sols)}{'(CAP)' if capped else ''} named_of_type={len(nm)} named_missing={len(missing) if not capped else '?'} constants={len(const)} others={len(others)} others_inj_on_periodic<=7={inj_cand} (of {min(len(others),3000)} tested) commut_fail={bad} time={time.time()-t0:.0f}s",flush=True)
