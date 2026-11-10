"""Rule 59R, Theorem 6.14 on rings (random data, fixed seed).
(a) 300 rings of length 3000: F and F^-1 commute with the flip of a random set of walls of both families, and
    the sets of walls are invariant.
(b) 200 rings of length 2000 seeded with trains of 1 to 6 walls, all m = 2, ..., 6 permuted at once:
    F Psi = Psi F along 60 steps, F^-1 Psi = Psi F^-1, the homomorphism law, and the inverse.
(c) The witness of injectivity has walls of family A exactly at the m centres and no wall of family B, m <= 7.
(d) The value of Psi at a cell does not change when all cells outside radius 2m+2 are re-randomized.
"""
import numpy as np, itertools
from hw_ca import *
Ft,Fi=tables(59)
F=lambda c: step(Ft,c); Fb=lambda c: step(Fi,c)
def wallmask(c,dl):
    a=np.roll(c,2); L=np.roll(c,1); R=np.roll(c,-1)
    return (L<3)&((a==L)|(a==3))&(R==(L+dl)%3)&((c==(L+dl)%3)|(c==3))
def flip(c,I,dl):
    c=c.copy(); L=np.roll(c,1)
    for i in I:
        p=L[i]; c[i]= 3 if c[i]==(p+dl)%3 else (p+dl)%3
    return c
def bit(c,i): return int(c[i]==3)^int(c[i-1]==2)
def trains(c):
    n=len(c); A=set(np.nonzero(wallmask(c,1))[0].tolist()); out=[]
    for i in A:
        if (i-2)%n in A: continue
        ch=[i]
        while (ch[-1]+2)%n in A: ch.append((ch[-1]+2)%n)
        out.append(ch)
    return out   # whole-ring chains (no start) are ignored = not finite maximal
def Psi(c,perm):   # perm: dict m -> list, new bit j = old bit perm^-1(j)
    I=[]
    for ch in trains(c):
        m=len(ch)
        if m not in perm: continue
        P=perm[m]; b=[bit(c,i) for i in ch]
        nb=[b[P.index(j)] for j in range(m)]
        I+=[ch[j] for j in range(m) if nb[j]!=b[j]]
    return flip(c,I,1)
rng=np.random.default_rng(12345)
def seed(n,density=0.3):
    c=rng.integers(0,4,n)
    s=0
    while s<n-30:
        if rng.random()<density:
            m=int(rng.integers(1,7)); p=int(rng.integers(3))
            c[s]=rng.choice([p,3]); c[s+1]=p
            for j in range(m):
                c[s+2+2*j]=rng.choice([(p+1+j)%3,3]); c[s+3+2*j]=(p+1+j)%3
            s+=2*m+4
        else: s+=int(rng.integers(1,6))
    return c
if __name__=="__main__":
    # (a) commutation of F with iota_I for random subsets I of all walls (A and B), many walls
    bad=0; nw=0
    for t in range(300):
        c=seed(3000) if t%2 else rng.integers(0,4,3000)
        for T in (Ft,Fi):
            A=np.nonzero(wallmask(c,1))[0]; B=np.nonzero(wallmask(c,2))[0]; nw+=len(A)+len(B)
            IA=A[rng.random(len(A))<0.5]; IB=B[rng.random(len(B))<0.5]
            g=flip(flip(c,IA,1),IB,2)
            Tc=step(T,c)
            h=flip(flip(Tc,IA,1),IB,2)
            if not (step(T,g)==h).all(): bad+=1
            # wall sets invariant
            if not ((wallmask(Tc,1)==wallmask(c,1)).all() and (wallmask(g,1)==wallmask(c,1)).all() and (wallmask(g,2)==wallmask(c,2)).all()): bad+=1
            c=Tc
    print('iota_I commutation / wall-set invariance failures:',bad,' total walls seen',nw)
    # (b) Theorem 6.14 on long rings, several m simultaneously
    bad=0; hb=0; nontriv=0; ntr=0
    for t in range(200):
        c=seed(2000)
        perm={m:list(rng.permutation(m)) for m in range(2,7)}
        perm2={m:list(rng.permutation(m)) for m in range(2,7)}
        g=Psi(c,perm); nontriv+=(g!=c).any(); ntr+=sum(1 for ch in trains(c) if len(ch)>1)
        x=c.copy(); y=g.copy()
        for s in range(60):
            x=F(x); y=F(y)
            if not (y==Psi(x,perm)).all(): bad+=1; break
        if not (Fb(g)==Psi(Fb(c),perm)).all(): bad+=1
        comp={m:[perm[m][perm2[m][j]] for j in range(m)] for m in perm}
        if not (Psi(Psi(c,perm2),perm)==Psi(c,comp)).all(): hb+=1
        inv={m:[perm[m].index(j) for j in range(m)] for m in perm}
        if not (Psi(Psi(c,perm),inv)==c).all(): hb+=1
    print('Theorem 6.14: commute failures',bad,'hom/inverse failures',hb,'nontrivial',nontriv,'of 200; trains(len>=2) seen',ntr)
    # (c) injectivity witness: ...33 p C1 (p+1) C2 ... Cm (p+m) 33...
    for m in range(1,8):
        for p in range(3):
            for bits in itertools.product([0,1],repeat=m):
                w=[3,3,3,p]
                for j in range(1,m+1): w+=[3 if bits[j-1] else (p+j)%3,(p+j)%3]
                w+=[3,3,3]
                c=np.array(w); A=np.nonzero(wallmask(c,1))[0].tolist(); B=np.nonzero(wallmask(c,2))[0].tolist()
                exp=[4+2*(j-1) for j in range(1,m+1)]
                assert A==exp and B==[], (m,p,bits,A,B)
    print('injectivity witness: A-walls exactly at the m centres, no B walls: OK')
    # (d) radius: block map test. Psi output at 0 compared with Psi applied to modified config outside radius 2m+2
    bad=0
    for t in range(3000):
        m=int(rng.integers(1,5)); c=seed(200); P={m:list(rng.permutation(m))}
        g=Psi(c,P); i=int(rng.integers(40,160)); d=c.copy()
        Rr=2*m+2
        mask=np.ones(200,bool); mask[i-Rr:i+Rr+1]=False
        d[mask]=rng.integers(0,4,mask.sum())
        if Psi(d,P)[i]!=g[i]: bad+=1
    print('radius 2m+2 violations',bad)
