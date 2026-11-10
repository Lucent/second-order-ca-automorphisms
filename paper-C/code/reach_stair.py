# The construction of Lemma 2.3 (free staircases) and the step of Theorem 2.4, with random staircase data,
# for all rules that read c and n up to 31.
import random
def fb(r,a,b,c): return (r>>(4*a+2*b+c))&1
def build(r,n,S):
    """S: dict (t,i)->bit on the staircase. Returns frame-0 config dict i->(x,y) on [n-K,n-1] via the construction of the lemma."""
    K=T=2*n+2; Y={}
    tau=lambda m: T-(n-1-m)
    for m in (n-1,n-2):
        for t in range(-1,tau(m)+1):
            Y[(t,m-t)]=S.get((t,m-t),0)
    for m in range(n-3,n-2-K,-1):
        Y[(tau(m),m-tau(m))]=0
        for t in range(tau(m),-1,-1):
            i=m-t+1
            Y[(t-1,i)]=Y[(t+1,i)]^fb(r,Y[(t,i-1)],Y[(t,i)],Y[(t,i+1)])
    return {i:(Y[(-1,i)],Y[(0,i)]) for i in range(n-K,n)}
def run(r,cfg,lo,hi,n):
    x={i:cfg[i][0] for i in range(lo,hi+1)}; y={i:cfg[i][1] for i in range(lo,hi+1)}
    rows={-1:x,0:y}
    for t in range(1,n+1):
        rows[t]={i:rows[t-2][i]^fb(r,rows[t-1][i-1],rows[t-1][i],rows[t-1][i+1]) for i in range(lo+t,hi-t+1)}
    return rows
random.seed(1)
bad=0; badfree=0; tests=0
full=[r for r in range(256) if any(fb(r,a,b,0)!=fb(r,a,b,1) for a in (0,1) for b in (0,1))]
for r in full:
    PR=[(a,b) for a in (0,1) for b in (0,1) if fb(r,a,b,0)!=fb(r,a,b,1)]
    for n in (1,2,3,7,15,31):
        # (a) random staircase values: check realised (Lemma 2.3)
        S={}
        for k in range(n):
            for j in (n-k-2,n-k-1): S[(k,j)]=random.randint(0,1)
        cfg=build(r,n,S); lo=min(cfg); hi=n-1
        rows=run(r,cfg,lo,hi,n)
        if any(rows[k][j]!=v for (k,j),v in S.items()): badfree+=1
        # (b) Theorem 2.4: staircase pairs = (a*,b*) in P_R, flip y_n
        ab=random.choice(PR); S={}
        for k in range(n): S[(k,n-k-2)]=ab[0]; S[(k,n-k-1)]=ab[1]
        cfg=build(r,n,S); cfg[n]=(random.randint(0,1),random.randint(0,1))
        cfg2=dict(cfg); cfg2[n]=(cfg[n][0],1-cfg[n][1])
        lo=min(cfg)
        # extend left/right enough so cell 0 at time n is computed: need [-n,n]
        assert lo<=-n
        o1=run(r,cfg,lo,n,n)[n][0]; o2=run(r,cfg2,lo,n,n)[n][0]
        tests+=1
        if o1==o2: bad+=1
print('rules',len(full),'tests',tests,'Lemma 2.3 failures',badfree,'Theorem 2.4 failures',bad)
