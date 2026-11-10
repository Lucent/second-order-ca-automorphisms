# Exhaustive near-glue simulation (library for check_gluing.py). Column r-1 = arbitrary input sequence w_t.
import sys
from base import *
from kinks import vcell, Kminus0, Kplus0, evolve_list
def family(K,n,lo,hi):
    m=abs(n)+2
    w=[K(i) for i in range(lo-m,hi+m+1)]
    w=evolve_list(w,n,0)
    return w[m-abs(n):m-abs(n)+hi-lo+1]
SPAN=(-200,200)
KM=[family(Kminus0,j,*SPAN) for j in range(6)]
KP=[family(Kplus0,j,*SPAN) for j in range(3)]
def km(j,shift,i): return KM[j][i-shift-SPAN[0]]
def kp(n,i):
    m,j=divmod(n,3); return KP[j][i-m-SPAN[0]]
def const_rows(c0,T):
    rows={-1:c0>>1,0:c0&1}
    for t in range(0,T+2):
        y=rows[t]; rows[t+1]=rows[t-1]^f(y,y,y)
    return rows
def match_kplus(r,t,Pa,Pb,W,cright):
    cells=[2*Pa[i]+Pb[i] for i in range(1,W+1)]
    # first column (>=r+1) where cells differ from every vacuum phase consistent with r+1,r+2
    for nm in range(6):
        if vcell(nm,r+1)!=cells[0] or vcell(nm,r+2)!=cells[1]: continue
        q=next((i for i in range(1,W+1) if cells[i-1]!=vcell(nm,r+i)),None)
        if q is None: return None
        for n in range(3*(r+q)-15,3*(r+q)+15):
            if n%6!=nm: continue
            if all(kp(n,r+i)==cells[i-1] for i in range(1,W+1)) and \
               all(kp(n+d,j)==vcell(n+d,j) for d in range(0,4) for j in (r+1,r+2,r+3)):
                # beyond the window: K^+_n is constant, equal to the right constant of the
                # simulated configuration at time t, on [r+W+1, oo)
                assert all(kp(n,i)==cright for i in range(r+W+1,n//3+SPAN[1]-1)), (r,t,n)
                return n
    return None
def run(r,tau,edge,shift,Tmax=None):
    for i in range(r+1,r+3):
        if km(tau,shift,i)!=vcell(tau,i): return None
    p=next(i for i in range(r+1,SPAN[1]) if km(tau,shift,i)!=vcell(tau,i))
    if Tmax is None: Tmax=3*(p-r)+40
    # cs: the kink is constant on [cs, oo).  The columns r..r+W are simulated, with the evolution of
    # the right constant as the value of column r+W+1.  By the light cone, column r+W+1 of the true
    # configuration at time t depends only on its initial cells on [r+W+1-t, r+W+1+t], which are
    # constant for t <= r+W+1-cs; with W = Tmax+3+(cs-r) this covers every t <= Tmax+4, so the
    # simulated columns r..r+W are exact up to time Tmax, and the true columns > r+W are constant.
    cs=max(i for i in range(r+1,SPAN[1]-40) if km(tau,shift,i)!=km(tau,shift,i-1))
    cR=km(tau,shift,cs)
    assert all(km(tau,shift,i)==cR for i in range(cs,SPAN[1]-40))
    W=Tmax+3+(cs-r)
    assert r+W+1-cs>=Tmax+4
    CR=const_rows(cR,Tmax+3)
    P0=tuple(edge>>1 if i==0 else km(tau,shift,r+i)>>1 for i in range(W+1))
    P1=tuple(edge&1 if i==0 else km(tau,shift,r+i)&1 for i in range(W+1))
    lam=(tau+3*(r&1))%6
    ucands=[u for u in range(max(0,3*(p-r)-14),Tmax) if (u+lam)%6==2]
    states={(P0,P1,(edge>>1,edge&1),())}
    t=0; results=[]; unresolved=0
    while states:
        if t>=Tmax: unresolved=len(states); break
        new=set()
        for (Pa,Pb,ds,Ec) in states:
            if t in ucands and ds is not None:
                Ec=Ec+tuple(((t,Pa[0],Pb[0]),tE,Pa[0],Pb[0],True) for tE in range(6))
            if t>=1:
                n=match_kplus(r,t,Pa,Pb,W,2*CR[t-1]+CR[t])
                if n is not None:
                    alive=tuple(sorted((u,tE) for (u,tE,ea,eb,al) in Ec if al and (tE+t-n)%6==0))
                    results.append((alive,n-t,t))
                    continue
            m = (Pb[0]==1 and Pb[1]==1) or (ds is not None and ds[1]==1 and Vac(t+tau,r+1)==1) or \
                any(al and eb==1 and Vac(t+tE,r+1)==1 for (u,tE,ea,eb,al) in Ec)
            for w in ((0,1) if m else (0,)):
                nb=[0]*(W+1)
                nb[0]=Pa[0]^f(w,Pb[0],Pb[1])
                for i in range(1,W+1):
                    right=Pb[i+1] if i<W else CR[t]
                    nb[i]=Pa[i]^f(Pb[i-1],Pb[i],right)
                nds=None
                if ds is not None:
                    nd=ds[0]^f(w,ds[1],Vac(t+tau,r+1))
                    if nd==nb[0]: nds=(ds[1],nd)
                nE=tuple((u,tE,eb,ea^f(w,eb,Vac(t+tE,r+1)),al and (ea^f(w,eb,Vac(t+tE,r+1)))==nb[0]) for (u,tE,ea,eb,al) in Ec)
                new.add((Pb,tuple(nb),nds,nE))
        states=new; t+=1
    return p,ucands,results,unresolved
if __name__=='__main__':
    r=int(sys.argv[1])
    allok=True
    for tau in range(6):
        for edge in range(4):
            if edge==vcell(tau,r): continue
            for shift in range(0,2*int(sys.argv[2]),2):
                out=run(r,tau,edge,shift)
                if out is None: continue
                p,uc,res,unres=out
                if p-r>int(sys.argv[3]): continue
                alives=set(a for a,_,_ in res)
                n0s=sorted(set((a,n0) for a,n0,_ in res))
                tmax=max(t for _,_,t in res)
                ok = unres==0 and all(len(a)>=1 for a in alives)
                if not ok: allok=False
                print('r',r,'tau',tau,'edge',edge,'p',p,'unres',unres,'branches',len(res),'tres<=',tmax,'(alive (u,tE), n0):',n0s,flush=True)
    print('ALLOK',allok)
