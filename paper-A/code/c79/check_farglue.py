# Confirmation (not used in the proof): far-glue probe by brute-force simulation; compare debris with S_u^{tau'} rule (forward and rho-backward).
import numpy as np, random, sys
from base import *
from kinks import vcell, Kminus0, Kplus0
import sys as _s; _a=_s.argv; _s.argv=[_a[0],"0","1","3"]
from glue import km as _km, kp as _kp, KM, KP, SPAN
_s.argv=_a
def km(j,sh,i):
    k=i-sh
    if k>SPAN[1]-60: return KM[j][-30]
    return _km(j,sh,i)
def kp(n,i):
    m,j=divmod(n,3)
    k=i-m
    if k>SPAN[1]-60: return KP[j][-30]
    if k<SPAN[0]+60: return KP[j][k-SPAN[0]] if k-SPAN[0]>=0 else vcell(j,i-m) if False else None
    return KP[j][k-SPAN[0]]
RHO=[0,2,1,3]
rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 1)
LO=-1500; HI=1500
def evolve(st,lo,n):
    """st list of cells starting at column lo; exact shrinking"""
    w=st
    for _ in range(abs(n)): w=Fstep_word(w) if n>0 else Finv_word(w)
    return w,lo+abs(n)
def tail_phase(st,lo,i0,i1):
    for ph in range(6):
        if all(st[i-lo]==vcell(ph,i) for i in range(i0,i1)): return ph
    return None
def probe(cells,lo,tau,a,N,S):
    """cells: time-0 list from lo, tail tau beyond r. glue F^? K- at time a with first nonvac col near N."""
    st,l=evolve(cells,lo,a) if a!=0 else (cells,lo)
    ph=(tau+a)%6
    # kink of phase ph: F^ph K-_0 shifted by even shift so that first nonvac col >= N
    sh=0
    while True:
        p=next(i for i in range(-100,400) if km(ph,sh,i)!=vcell(ph,i)) if True else None
        if p>=N: break
        sh+=2
    hi=l+len(st)-1
    g=[st[i-l] if i<N-5 else km(ph,sh,i) for i in range(l,hi+1)]
    assert all(st[i-l]==vcell(ph,i) for i in range(1,N-5))
    st2,l2=evolve(g,l,S); t=a+S
    return st2,l2,t
def S_rule(cells,lo,tau,u,tn,Lk=60):
    st,l=evolve(cells,lo,u)
    st=[st[i-l] if i<=0 else vcell(tn+u,i) for i in range(l,l+len(st))]
    st,l=evolve(st,l,-u)
    return st,l
def zrow(cells,lo,t):
    if t>=0:
        st,l=evolve(cells,lo,t); return st[0-l]&1
    st,l=evolve(cells,lo,t+1); return st[0-l]>>1
bad=0; n=0; Ns=set()
for trial in range(int(sys.argv[2]) if len(sys.argv)>2 else 40):
    tau=rng.randrange(6)
    cells=[rng.randrange(4) if i<=0 else vcell(tau,i) for i in range(LO,HI)]
    while cells[0-LO]==vcell(tau,0): cells[0-LO]=rng.randrange(4)
    a=rng.randrange(-12,13); N=rng.randrange(25,60); S=3*N+60
    st,l,t=probe(cells,LO,tau,a,N,S)
    # right part: find K+ trajectory n with F^n K+_0 on cols [5, 3N] at time t
    found=None
    for nn in range(-150,560):
        if all(kp(nn,i)==st[i-l] for i in range(5,min(200,l+len(st)-1))):
            found=nn; break
    assert found is not None,'no K+'
    tE=(found)%6  # phase at time t of left vacuum = found mod 6; tail at time 0 = found - t
    tE0=(found-t)%6
    # debris: cols <=0 at time t, plus vacuum tE
    d=[st[i-l] if i<=0 else vcell(tE0+t,i) for i in range(l,l+len(st))]
    assert all(st[i-l]==vcell(tE0+t,i) for i in range(1,5))
    D,lD=evolve(d,l,-t)
    # rule: forward u = the one with (u+tau)%6==2 ; search all candidates u in plausible window and report which match
    lam=tau
    cands=[u for u in range(a+3*N-40,a+3*N+20) if (u+lam)%6==2]
    match=[]
    for u in cands:
        z=zrow(cells,LO,u-1); tn=(tau+2)%6 if z==0 else (tau-1)%6
        if tn!=tE0: continue
        E,lE=S_rule(cells,LO,tau,u,tn)
        lo2=max(lD,lE)
        if all(E[i-lE]==D[i-lD] for i in range(lo2+5,1)): match.append(u)
    n+=1
    if not match: bad+=1; print('NO MATCH',tau,a,N)
    Ns.add((found-t) - (-3*0))
print('forward far-glue tests',n,'bad',bad)
assert bad==0
# --- n0 formula & rho-probe backward rule
bad=0;n=0
for trial in range(int(sys.argv[2]) if len(sys.argv)>2 else 40):
    tau=rng.randrange(6)
    cells=[rng.randrange(4) if i<=0 else vcell(tau,i) for i in range(LO,HI)]
    while cells[0-LO]==vcell(tau,0): cells[0-LO]=rng.randrange(4)
    rc=[RHO[x] for x in cells]
    rtau=tail_phase(rc,LO,1,40)
    a=rng.randrange(-12,13); N=rng.randrange(25,50); S=3*N+60
    st,l,t=probe(rc,LO,rtau,a,N,S)
    found=next(nn for nn in range(-150,560) if all(kp(nn,i)==st[i-l] for i in range(5,200)))
    tE0=(found-t)%6
    d=[st[i-l] if i<=0 else vcell(tE0+t,i) for i in range(l,l+len(st))]
    D,lD=evolve(d,l,-t)
    rD=[RHO[x] for x in D]; tD=tail_phase(rD,lD,1,40)
    # backward rule on c: v with (v+tau)%6==0, z_v, tau' = tau-2 if z=0 else tau+1
    match=[]
    for v in range(-a-3*N-40,-a-3*N+40):
        if (v+tau)%6: continue
        z=zrow(cells,LO,v); tn=(tau-2)%6 if z==0 else (tau+1)%6
        if tn!=tD: continue
        E,lE=S_rule(cells,LO,tau,v,tn)
        lo2=max(lD,lE)
        if all(E[i-lE]==rD[i-lD] for i in range(lo2+5,1)): match.append(v)
    n+=1
    if not match: bad+=1; print('NO MATCH rho',tau,a,N)
print('rho far-glue tests',n,'bad',bad)
assert bad==0
