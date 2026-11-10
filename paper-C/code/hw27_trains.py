"""Rule 27R, Theorem 6.20: (H4) on all words of length 9 to 11 (a flip does not change the bit beta of another
wall), and Psi_m on 150 random rings of length 400 for each gap d = 3, 6, 7 (fixed seed): commutation with F
and F^-1 and the homomorphism law.
"""
import itertools, random, numpy as np
from hw_common import tables, step
T,Ti=tables(27)
W=[(a,p,p,C,(p+1)%3,(p+1)%3) for p in range(3) for a in (p,3) for C in ((p+1)%3,3)]
S=set(W)
def io(w): w=list(w); w[3]=3 if w[3]!=3 else (w[1]+1)%3; return tuple(w)
def beta(w): return (int(w[3]==3)+int(w[1]==2))%2
def occ(c):
    n=len(c); return {i for i in range(n) if tuple(int(c[(i+t)%n]) for t in range(6)) in S}
# H4 exhaustive (linear words)
bad=0
def occl(c): return [i for i in range(len(c)-5) if tuple(c[i:i+6]) in S]
for n in range(9,12):
  for c in itertools.product(range(4),repeat=n):
    o=occl(c)
    if len(o)<2: continue
    for i in o:
      ci=list(c); ci[i:i+6]=io(tuple(c[i:i+6]))
      for i2 in o:
        if i2!=i and beta(tuple(ci[i2:i2+6]))!=beta(tuple(c[i2:i2+6])): bad+=1
print("H4 violations",bad)
def Psi(c,m,sig,d):
    n=len(c); o=occ(c); flips=[]
    for i in o:
        if (i-d)%n in o: continue
        L=1
        while (i+L*d)%n in o and L<=n: L+=1
        if L!=m: continue
        mem=[(i+r*d)%n for r in range(m)]
        b=[beta(tuple(int(c[(q+t)%n]) for t in range(6))) for q in mem]
        inv=[0]*m
        for r in range(m): inv[sig[r]]=r
        flips+=[mem[r] for r in range(m) if b[inv[r]]!=b[r]]
    c2=c.copy()
    for q in flips: c2[(q+3)%n]=io(tuple(int(c[(q+t)%n]) for t in range(6)))[3]
    return c2
random.seed(7)
for d in (3,6,7):
  fails=0; nt=0
  for _ in range(150):
    n=400; c=[random.randrange(4) for _ in range(n)]; pos=0
    while pos<n-60:
        pos+=random.randrange(1,10); m=random.randrange(1,5); p=random.randrange(3)
        if d==3:
            seg=[random.choice([p,3]),p,p]
            for t in range(1,m+1): seg+=[random.choice([(p+t)%3,3]),(p+t)%3,(p+t)%3]
            c[pos:pos+len(seg)]=seg; pos+=len(seg)
        else:
            for r in range(m):
                p=random.randrange(3); c[pos:pos+6]=[random.choice([p,3]),p,p,random.choice([(p+1)%3,3]),(p+1)%3,(p+1)%3]; pos+=d
            pos+=6
    c=np.array(c[:n],dtype=np.int8)
    for m in range(1,5):
        sig=tuple(random.sample(range(m),m)); rho=tuple(random.sample(range(m),m))
        if not (Psi(step(T,c),m,sig,d)==step(T,Psi(c,m,sig,d))).all(): fails+=1
        if not (Psi(step(Ti,c),m,sig,d)==step(Ti,Psi(c,m,sig,d))).all(): fails+=1
        comp=tuple(sig[rho[r]] for r in range(m))
        if not (Psi(Psi(c,m,rho,d),m,sig,d)==Psi(c,m,comp,d)).all(): fails+=1
        if (Psi(c,m,sig,d)!=c).any(): nt+=1
  print("rule27 d",d,"fails",fails,"nontrivial",nt)
