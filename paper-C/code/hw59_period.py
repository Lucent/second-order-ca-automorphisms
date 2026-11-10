"""Rule 59R: for every word w of X_A and X_B and p in {6, 12, 30, 60, 120}, the number of 20000 random contexts
(fixed seed) in which the window of F^p equals w again. A count below 20000 exhibits a context in which
F^p changes the word (its protector), as stated after Corollary 6.9 in Section 6.2.
"""
import numpy as np
from hw_ca import *
Ft,Fi=tables(59)
rng=np.random.default_rng(5)
def fam(dl): return sorted({(a,p,C,(p+dl)%3) for p in range(3) for a in (p,3) for C in ((p+dl)%3,3)})
N=20000
for w in fam(1)+fam(2):
    res=[]
    for p in (6,12,30,60,120):
        L=p+3; c=rng.integers(0,4,(N,2*L+4)); c[:,L:L+4]=w
        x=c
        for _ in range(p): x=Ft[x[:,:-2],x[:,1:-1],x[:,2:]]
        off=L-p; win=x[:,off:off+4]
        res.append(int((win==np.array(w)).all(1).sum()))
    print(w,'#contexts (of %d) with F^p window == w for p=6,12,30,60,120:'%N,res)
