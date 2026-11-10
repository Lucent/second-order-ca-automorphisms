"""Exact verification for rule 133R (and 161R via the complement):
 (1) W = 0330 is fixed by F in every context (all 16 contexts);
 (2) D: the radius-4 block map that applies the car rotation s -> F(s) (= 0->1->3->2->0) to every cell s in
     an occurrence of W s W and is the identity elsewhere, commutes with F on all 4^11 windows (numpy check
     from classify_r4.check_commutes_np_R);
 (3) D^4 = id and D is injective on all periodic configurations of period <= 9 (sanity);
 (4) D is not sigma^a F^b (numerically: D((0330 s)^inf) = F((0330 s)^inf), D(0^Z)=0^Z, F(0^Z)=1^Z)."""
import sys, itertools, numpy as np

from classify_r4 import check_commutes_np_R
from centralizer import make_rule
for rule,W in ((133,(0,3,3,0)),(161,(3,0,0,3))):
    f,F_loc,Finv=make_rule(rule)
    ok=all(tuple(F_loc(*ctx[i:i+3]) for i in range(4))==W for ctx in (tuple([a])+W+tuple([b]) for a in range(4) for b in range(4)))
    car={s:F_loc(W[-1],s,W[0]) for s in range(4)}
    print(f"rule {rule}: wall {''.join(map(str,W))} fixed in every context: {ok}; single-car map {car}")
    R=4; n=9; idx=np.arange(4**n,dtype=np.int64); cells=[(idx>>(2*(n-1-i)))&3 for i in range(n)]
    isW=lambda off: np.all([cells[off+i]==W[i] for i in range(4)],axis=0)
    m=isW(0)&isW(5); ct=np.array([car[s] for s in range(4)])
    tab=np.where(m,ct[cells[4]],cells[4])
    print(f"  D: radius {R}, {int(m.sum())} windows changed; commutes with F on all 4^{2*R+3} windows:", check_commutes_np_R(tab,F_loc,R))
    def apply(tab,c):
        L=len(c); return [int(tab[sum(c[(i+j-R)%L]<<(2*(n-1-j)) for j in range(n))]) for i in range(L)]
    inj=True; ord4=True
    for p in range(1,10):
        seen={}
        for e in range(4**p):
            c=[(e>>(2*(p-1-i)))&3 for i in range(p)]; d=apply(tab,c); t=tuple(d)
            if t in seen: inj=False
            seen[t]=1
            dd=d
            for _ in range(3): dd=apply(tab,dd)
            if dd!=c: ord4=False
    print("  injective on periodic configs of period<=9:",inj,"; D^4=id there:",ord4)
    Fc=lambda c:[F_loc(c[i-1],c[i],c[(i+1)%len(c)]) for i in range(len(c))]
    for s in range(4):
        c=list(W)+[s]; assert apply(tab,c)==Fc(c)
    z=[0]*5; print("  D(0^Z)==0^Z:",apply(tab,z)==z,"; F(0^Z) =",Fc(z)[0],"^Z (const orbit:", [Fc([s]*5)[0] for s in range(4)],")")
