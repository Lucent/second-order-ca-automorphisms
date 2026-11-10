import numpy as np, random
from core import *
from strip import tree
K=16
treesets={}
for beta in (0,1):
    cyc=tree(beta,K,verbose=False)[-1]
    S=set()
    for c in cyc:
        L=len(c)
        for t in range(L): S.add(tuple(tuple(c[(t+d)%L]) for d in range(32)))  # 32-step windows of the periodic sequence
    treesets[beta]=S
rng=random.Random(5); badrec=0; notin=0; tot=0
for tr in range(300):
    T=400; L=2*T+200; s=T+100
    ph=rng.randint(0,1); i=np.arange(L)
    Y0=((i+ph)%2).astype(np.uint8); Ym1=((i+ph+1)%2).astype(np.uint8)
    w=rng.randint(1,15)
    for j in range(s,s+w): Y0[j]=rng.randint(0,1); Ym1[j]=rng.randint(0,1)
    rows=diagram(Ym1,Y0,T)
    beta=int(Y0[s-1])   # Z_{-1}
    assert Y0[s-2]==1-beta
    Z=lambda t,j: int(rows[t][s-t+j])
    for t in range(1,T-1):
        for j in range(0,K):
            zm2=lambda tt: Z(tt,j-2)
            if Z(t+1,j)!=(Z(t,j)^Z(t,j-2)^Z(t-1,j-2)^N(Z(t,j-2),Z(t,j-1),Z(t,j))): badrec+=1
    # late strip: lines 0..K-1 at t=T-40..T-9 ; genuine fronts only (Z_0==beta)
    if Z(0,0)==beta:
        tot+=1
        seq=tuple(tuple(Z(t,j) for j in range(K)) for t in range(T-40,T-8))
        if seq not in treesets[beta]: notin+=1
print('recursion violations',badrec,'genuine fronts',tot,'late strip not on a tree cycle',notin)
