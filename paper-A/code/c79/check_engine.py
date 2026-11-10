# Confirmation (not used in the proof): the three-valued engine against brute-force full simulation of S_u on random wide configs
import numpy as np, random
from base import *
from engine import *
rng=random.Random(5)
def full_switch(cells_left,tau,u,tn,Wl=80,Wr=40):
    # cells on [-Wl,0] arbitrary, tail tau on (0,Wr]; evolve exactly with shrinking windows; pad big
    pad=abs(u)+5
    cols=list(range(-Wl-2*pad,Wr+2*pad+1))
    w=[cells_left[i] if i<=0 else vcell(tau,i) for i in cols]
    # pad left with random
    for n in range(abs(u)):
        w=Fstep_word(w) if u>0 else Finv_word(w)
    lo=cols[0]+abs(u)
    st={lo+j:w[j] for j in range(len(w))}
    st2=[st[i] if i<=0 else vcell(tn+u,i) for i in range(lo,lo+len(w))]
    w=st2
    for n in range(abs(u)):
        w=Finv_word(w) if u>0 else Fstep_word(w)
    lo2=lo+abs(u)
    return {lo2+j:w[j] for j in range(len(w))}
L=4; bad=0; tests=0
for trial in range(300):
    tau=rng.randrange(6)
    cl={i:rng.randrange(4) for i in range(-300,1)}
    if cl[0]==vcell(tau,0): continue
    cells=np.array([[cl[i] for i in range(-L,1)]],np.uint8)
    D=Diag(tau,cells,L,L+1)
    for u in range(-L,L+1):
        m=move(D,u)
        if m is None: continue
        kind,nt,nc=m
        # true z
        res=full_switch(cl,tau,u,int(nt[0]))
        tests+=1
        exp=[res[i] for i in range(-L,1)]
        got=list(nc[0])
        if any(g!=10 and g!=e for g,e in zip(got,exp)) or any(g==10 for g in got): bad+=1; print('mismatch',tau,u,kind,got,exp)
        # also z check from full sim
print('tests',tests,'bad',bad)
assert bad==0
# z check
def fullY(cl,tau,t,i,W=200):
    cols=list(range(-W,W)); w=[cl.get(j,0) if j<=0 else vcell(tau,j) for j in cols]
    if t>=0:
        for n in range(t): w=Fstep_word(w)
        return w[i-cols[0]-t]&1
    else:
        for n in range(-t-1): w=Finv_word(w)
        return w[i-cols[0]-(-t-1)]>>1
bad=0
for trial in range(100):
    tau=rng.randrange(6); cl={i:rng.randrange(4) for i in range(-150,1)}
    if cl[0]==vcell(tau,0): continue
    cells=np.array([[cl[i] for i in range(-L,1)]],np.uint8); D=Diag(tau,cells,L,L+1)
    for t in range(-L-1,L+1):
        if D.z(t)[0]!=fullY(cl,tau,t,0): bad+=1
print('z bad',bad)
assert bad==0
