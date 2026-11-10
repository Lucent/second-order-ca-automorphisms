# Runs the proof as an algorithm: for an object X, produce and check the chain of scatterings that
# certifies X in the closure, with window length LW.  Every scattering identity is checked by simulation.
import verify1, sys, itertools
# usage: python3 verify2.py LW NMAX [lib|direct|goursat_cc]
LW=int(sys.argv[1]); NMAX=int(sys.argv[2])
verify1.use(sys.argv[3] if len(sys.argv)>3 else 'lib')
def S(P,X): return verify1.S(P,X)
stats={'passes':0,'maxpass':0}
def shift(X,s): return (X[0],[k+s for k in X[1]])
def norm(X):
    # translate so that first kink is at 8 or 9 keeping parity (translation by 2 is harmless)
    if not X[1]: return X
    s=X[1][0]-8; s-=s%2
    return shift(X,-s)
def prepend_kink(b1,U):
    """U certified, tail c, first kink u1>=b1+2. Checks the growth chain for X=[b1;U]; returns number of passes."""
    lam=1-U[0]; n=0
    while True:
        u1=U[1][0]
        assert u1-b1>=2
        if u1-b1>=LW: return n          # (C0)
        X=(lam,[b1]+U[1])
        if lam==1:   # A-pass: probe barK_a, a !≡ b1
            a=2+((b1+1)%2)
            G,O,F=S((1,[a]),X); G2,O2,F2=S((0,[a]),U)
            assert G==(0,[b1]+G2[1]) and O==O2 and G2[0]==1 and G2[1][0]>=u1
        else:        # B-pass: probe K_a, a ≡ u1
            a=2+(u1%2)
            G,O,F=S((0,[a]),X); G2,O2,F2=S((1,[a]),U)
            assert G==(1,[b1]+G2[1]) and O==O2 and G2[0]==0 and G2[1][0]>=u1+1
        U=G2; lam=1-lam; n+=1
        assert n<10*LW
memo={}
def cert(X):
    X=norm(X); key=(X[0],tuple(X[1]))
    if key in memo: return True
    ks=X[1]; lam=X[0]
    if len(ks)<=1: memo[key]=1; return True
    if lam==1 and ks[1]-ks[0]==1:
        if len(ks)==2: memo[key]=1; return True      # Z, by the seed lemma
        d=ks[0]; U=(1,ks[2:]); u1=U[1][0]; y=u1-d-1
        assert y>=2 and cert(U)
        if y==2:
            a=2+((d+1)%2)                             # Z-probe, a !≡ d
            G,O,F=S((1,[a,a+1]),X); G2,O2,F2=S((1,[a,a+1]),U)
            assert G2[0]==1 and G==(1,[d,d+1]+G2[1]) and O==O2 and G2[1][0]>=u1+1
            U=G2; u1=U[1][0]; y=u1-d-1; X=(1,[d,d+1]+U[1])
        a=2+(d%2)                                     # kink probe a ≡ d
        G,O,F=S((1,[a]),X); G2,O2,F2=S((1,[a+1,a+2]),U)
        assert G2[0]==1 and G==(0,[d+2]+G2[1]) and O==O2 and G2[1][0]>=u1
        n=prepend_kink(d+2,G2)
    else:
        assert ks[1]-ks[0]>=2
        U=(1-lam,ks[1:]); assert cert(U)
        n=prepend_kink(ks[0],U)
    stats['passes']+=n; stats['maxpass']=max(stats['maxpass'],n)
    memo[key]=1; return True
cnt=0
for n in range(1,NMAX+1):
    for lam in (0,1):
        for bits in itertools.product('01',repeat=n):
            c=''.join(bits)
            if int(c[0])==lam: continue
            if n>1 and c[-1]==c[-2]: continue
            w=str(lam)*8+c+c[-1]*4
            if '010' in w: continue
            ks=[k for k in range(1,len(w)) if w[k]!=w[k-1]]
            assert cert((lam,ks)); cnt+=1
    print('core length',n,'objects certified so far',cnt,stats,flush=True)
