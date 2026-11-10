# Right transparency of 155R: c agrees with a Phi+ point on [w,inf). Claim: for t>=0, no firing at centres i>=w+1,
# and Y(t,i)=Y(t-1,i-1) for i>=w+1.  Exhaustive over junk on [0,w) (both rows), random Phi+ right media and random left parts.
import numpy as np, random, itertools, sys
RULE=155
ftab=np.array([(RULE>>k)&1 for k in range(8)],np.uint8)
Ntab=np.array([((RULE>>k)&1)^((k>>2)&1)^(k&1) for k in range(8)],np.uint8)
def phiword(n,rng,start):
    w=list(start)
    while len(w)<n:
        if w[-1]==1: w.append(0)
        elif w[-2]==0: w.append(1)
        else: w.append(rng.randint(0,1))
    return w
T=int(sys.argv[1]); wmax=int(sys.argv[2]); nmed=int(sys.argv[3])
bad=0; tot=0
rng=random.Random(1)
for w in range(0,wmax+1):
  for jx in range(4**w):
    for m in range(nmed):
        L=4*T+2*w+20; base=2*T+10   # junk at [base,base+w)
        u=phiword(L+2,rng,(1,0))     # Phi+ medium: y_i=u[i+1], x_i=u[i+2]
        y=np.array(u[1:L+1],np.uint8); x=np.array(u[2:L+2],np.uint8)
        # random left part (left of base) arbitrary
        for i in range(base): x[i]=rng.randint(0,1); y[i]=rng.randint(0,1)
        for k in range(w):
            x[base+k]=(jx>>(2*k))&1; y[base+k]=(jx>>(2*k+1))&1
        prev=None
        for t in range(T):
            k=4*np.roll(y,1)+2*y+np.roll(y,-1); fr=Ntab[k]
            seg=slice(base+w+1,base+w+1+T)   # stay away from ring wrap
            tot+=1
            if fr[seg].any(): bad+=1; break
            if prev is not None:
                # Y(t,i)=Y(t-1,i-1) for i>=base+w+1
                if (y[base+w+1:base+w+1+T]!=prev[base+w:base+w+T]).any(): bad+=1; break
            prev=y.copy()
            x,y=y,x^ftab[k]
print('W-check T',T,'wmax',wmax,'violations',bad,'checks',tot)
