# Right transparency of 155R, parts (a)-(d). usage: thmW.py SEED TRIALS
import numpy as np, random, sys
from core import *
rng=random.Random(int(sys.argv[1])); trials=int(sys.argv[2])
viol={'a':0,'b':0,'c':0,'d':0,'phi':0}; checks=0
for tr in range(trials):
    T=rng.randint(5,60); w0=rng.randint(0,10)
    L=4*T+60; W=L//2      # medium agrees on [W, L)
    phi=phi_extend_right(rng.choice([(1,0),(0,1),(0,0)]),L+3,rng)  # phi(k) for k=0..; phi in Phi_155 (no 11, no 000)
    s=''.join(map(str,phi))
    if '11' in s or '000' in s: continue
    # Phi+ point: Y(t,i)=phi(i-t) ; cell i = (Y(-1,i),Y(0,i)) = (phi(i+1),phi(i))
    Y0=np.array([phi[i] for i in range(L)],np.uint8); Ym1=np.array([phi[i+1] for i in range(L)],np.uint8)
    # arbitrary left part on [0,W)
    mode=rng.random()
    for i in range(W):
        if mode<0.5 or i>=W-w0: Y0[i]=rng.randint(0,1); Ym1[i]=rng.randint(0,1)
    rows=diagram(Ym1,Y0,T)
    lo=W; hi=L-2*T-4   # keep away from right seam
    for t in range(0,T):
        checks+=1
        fr=fires(rows[t])
        if fr[W+1:hi].any(): viol['b']+=1
        if t>=0 and (rows[t][W+1:hi]!=rows[t-1][W:hi-1]).any(): viol['a']+=1
        # (c): S(t)=Y(t,W) for t>=1
        if t>=1:
            S1,S2=rows[t-1][W],rows[t-2][W]
            if S1==1: exp=0
            elif S2==0: exp=1
            else: exp=rows[t-1][W-1]
            if rows[t][W]!=exp: viol['c']+=1
        # profile u(k)=Y(t,k+t) for k+t>=W : check in Phi
        u=''.join(str(int(v)) for v in rows[t][W:hi])
        if '11' in u or '000' in u: viol['phi']+=1
    # (d): windows centred i>=W+R are windows of a Phi+ point: rows t-1,t on [i-R,i+R] satisfy Y(t-1,j)=Y(t,j+1) and profile in Phi -> implied by (a)+(phi); check directly the cell relation x_j=y_{j+1}
    for t in range(0,T):
        if (rows[t-1][W:hi-1]!=rows[t][W+1:hi]).any(): viol['d']+=1
print('checks',checks,'violations',viol)
