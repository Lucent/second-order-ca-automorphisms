# Exhaustive separation check for finite perturbations of valid backgrounds. usage: separ.py RULE W T R
import numpy as np, sys, itertools
rule=int(sys.argv[1]); w=int(sys.argv[2]); T=int(sys.argv[3]); R=int(sys.argv[4])
ft=np.array([(rule>>k)&1 for k in range(8)],dtype=np.int8)
Nt=np.array([ft[k]^((k>>2)&1)^(k&1) for k in range(8)],dtype=np.int8)
bgs=[b for b in itertools.product((0,1),repeat=2) if Nt[4*b[0]+2*b[1]+b[0]]==0 and Nt[4*b[1]+2*b[0]+b[1]]==0]
M=2*T+w+40; M+=M%2; s0=M//2
idx=np.arange(4**w)
P=np.stack([(idx>>(2*k))&3 for k in range(w)],axis=1).astype(np.int8)
def isS(x,y):
    # x,y: (B,M). returns bool (B,) : every (2R+1)-window (inside safe zone) is Phi+ or Phi- window
    trip=Nt[4*np.roll(y,1,1)+2*y+np.roll(y,-1,1)]   # forbidden triple centred at j in the y-row
    plus=(x==np.roll(y,-1,1)); minus=(x==np.roll(y,1,1))
    ok=np.ones(x.shape[0],bool)
    B=x.shape[0]
    # sliding all via cumulative sums
    def allwin(a,lo,hi):  # all(a[j] for j in [i+lo, i+hi]) for each i
        c=np.concatenate([np.zeros((B,1),np.int32),np.cumsum(1-a.astype(np.int32),1)],1)
        n=a.shape[1]; out=np.zeros_like(a,dtype=bool)
        i=np.arange(n); L=(i+lo); H=(i+hi)
        good=(L>=0)&(H<n)
        out[:,good]=(c[:,H[good]+1]-c[:,L[good]])==0
        out[:,~good]=True   # ignore ring boundary windows (background there anyway)
        return out
    # Phi+ word: y_{i-R..i+R}, x_{i+R}=y_{i+R+1}: plus on [i-R,i+R], no forbidden triple centred in [i-R+1, i+R]
    pw=allwin(plus,-R,R)&allwin(trip==0,-R+1,R)
    mw=allwin(minus,-R,R)&allwin(trip==0,-R,R-1)
    return (pw|mw).all(1)
for b in bgs:
    y=np.tile(np.array([b[i%2] for i in range(M)],np.int8),(4**w,1))
    x=np.tile(np.array([b[(i-1)%2] for i in range(M)],np.int8),(4**w,1))
    x[:,s0:s0+w]=P>>1; y[:,s0:s0+w]=P&1
    late=0
    first_sep=np.full(4**w,-1); stays=np.ones(4**w,bool); lastfire=np.full(4**w,-1)

    for t in range(0,T):
        trip=Nt[4*np.roll(y,1,1)+2*y+np.roll(y,-1,1)]
        late+= int(trip.sum()) if t>=T//2 else 0
        any_f=trip.any(1); lastfire[any_f]=t
        s=isS(x,y)
        newly=(first_sep<0)&s; first_sep[newly]=t
        stays&= ~((first_sep>=0)&~s)
        x,y=y,x^ft[4*np.roll(y,1,1)+2*y+np.roll(y,-1,1)]

    # d'Alembert check on a sample: Y(q,p) via formula vs simulation, using frame data only (future t>=1)
    print(f'rule {rule} bg {b} w={w}: separated by T={T}: {np.sum(first_sep>=0)}/{4**w}; max sep time {first_sep.max()}; '
          f'left S_R after entering: {np.sum(~stays)}; last firing time max {lastfire.max()} (horizon {T-1}); '
          f'firings at t>=T/2: {late}')
