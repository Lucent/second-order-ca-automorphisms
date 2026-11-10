"""Rules 127R, 179R, 223R, Proposition 6.19: Psi on 150 rings of length 1500 for each rule and each gap
d = 5, 6, 8 (fixed seed), seeded with trains of walls of one class with about one wall in ten of another
class mixed in, all m = 2, ..., 5 permuted at once: commutation with F along 2P+3 steps and with F^-1, and
the homomorphism law.
"""
import numpy as np, itertools, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    import hwtp_walls as t5_tp
from hw_ca import *
k=5
for rule,(W,orb,iot,Ft,Fi) in t5_tp.res.items():
    P=len(orb); h=P//2; pos={w:i for i,w in enumerate(orb)}
    codes={}
    for w in W: codes[sum(x<<(2*(k-1-j)) for j,x in enumerate(w))]=w
    def occ(c):
        e=np.zeros(len(c),dtype=np.int64)
        for j in range(k): e=(e<<2)|np.roll(c,-j)
        return {int(i):codes[int(e[i])] for i in np.nonzero(np.isin(e,list(codes)))[0]}
    def Psi(c,perm,d):
        n=len(c); O=occ(c); out=c.copy()
        for i in O:
            if (i-d)%n in O: continue
            ch=[i]
            while (ch[-1]+d)%n in O and len(ch)<=n: ch.append((ch[-1]+d)%n)
            m=len(ch)
            if m not in perm: continue
            st=[pos[O[x]] for x in ch]
            if len({s%h for s in st})!=1: continue
            b=[s//h for s in st]; nb=[b[perm[m].index(j)] for j in range(m)]
            for x,s,bo,bn in zip(ch,st,b,nb):
                if bo!=bn:
                    nw=iot[O[x]]
                    for j in range(k): out[(x+j)%n]=nw[j]
        return out
    rng=np.random.default_rng(rule)
    def seed(n,d):
        c=rng.integers(0,4,n); s=0
        while s<n-d*8-10:
            if rng.random()<0.5:
                m=int(rng.integers(1,6)); s0=int(rng.integers(P))
                for j in range(m):
                    w=orb[(s0+h*int(rng.integers(2)))%P] if rng.random()<0.9 else orb[int(rng.integers(P))]
                    c[s+j*d:s+j*d+k]=w
                s+=m*d+int(rng.integers(1,d+3))
            else: s+=int(rng.integers(1,8))
        return c
    for d in (5,6,8):
        bad=0; hb=0; nt=0; ntr=0
        for t in range(150):
            c=seed(1500,d); perm={m:list(rng.permutation(m)) for m in range(2,6)}; perm2={m:list(rng.permutation(m)) for m in range(2,6)}
            g=Psi(c,perm,d); nt+=(g!=c).any()
            x,y=c.copy(),g.copy()
            for s in range(2*P+3):
                x=step(Ft,x); y=step(Ft,y)
                if not (Psi(x,perm,d)==y).all(): bad+=1; break
            if not (Psi(step(Fi,c),perm,d)==step(Fi,g)).all(): bad+=1
            comp={m:[perm[m][perm2[m][j]] for j in range(m)] for m in perm}
            if not (Psi(Psi(c,perm2,d),perm,d)==Psi(c,comp,d)).all(): hb+=1
        print(rule,'gap',d,'commute failures',bad,'hom failures',hb,'nontrivial',nt,'/150')
