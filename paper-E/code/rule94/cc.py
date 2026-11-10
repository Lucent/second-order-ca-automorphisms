# Rule 94 second-order (Fredkin) ECA: scattering by direct (t,i) simulation and by a light-cone (Goursat) solver.
# Object X in Phi^-: (lam, kinks) : xi(p) = lam ^ #{b in kinks : b<=p} mod 2.
# Probe P in Phi^+: (lam, kinks in q).
import numpy as np
def val(pt,k):
    lam,ks=pt; return (lam+sum(1 for x in ks if x<=k))&1
def rep(S,off=0):
    S=list(S); return (int(S[0]),[k+off for k in range(1,len(S)) if S[k]!=S[k-1]])
RULE=94
TBL=np.array([(RULE>>k)&1 for k in range(8)],np.uint8)
def direct(P,X,T=None):
    """Direct (t,i) simulation of F(x,y)_i=(y_i, x_i ^ f(y_{i-1},y_i,y_{i+1})).
    Corner [0,0]: rows t=-1,0; columns i<0 carry P (pi(t-i)), i>=0 carry X (xi(t+i)).
    Requires P,X equal to lam on q<=1, p<=1 resp.  Returns (Gamma, O) kink reps."""
    mx=max([2]+P[1]+X[1])
    if T is None: T=3*mx+30
    W=T+mx+10
    I=np.arange(-W,W+1)
    def row(t):
        r=np.empty(len(I),np.uint8)
        for n,i in enumerate(I):
            r[n]=val(P,t-i) if i<0 else val(X,t+i)
        return r
    prev=row(-1); cur=row(0)
    for t in range(T):
        idx=4*cur[:-2]+2*cur[1:-1]+cur[2:]
        nb=prev.copy(); nb[1:-1]=prev[1:-1]^TBL[idx]
        prev,cur=cur,nb
    # gamma(p)=Y(T,p-T) for p in [-5, 2T-q1]; o(q)=Y(T,T-q)
    L=mx+15
    g=[cur[(p-T)+W] for p in range(0,L)]
    o=[cur[(T-q)+W] for q in range(0,L)]
    G=rep(g); O=rep(o)
    return G,O
def goursat(P,X,M=None,fires=False):
    """Light-cone solver on quadrant q,p>=0 (corner [0,0]), pure python dict."""
    mx=max([2]+P[1]+X[1])
    if M is None: M=3*mx+30
    Y={}
    for k in range(M):
        Y[(k,k%2)]=val(P,k); Y[(k%2,k)]=val(X,k)
    F=[]
    for s in range(4,2*M,2):
        for q in range(2,s-1):
            p=s-q
            if p<2 or p>=M or q>=M: continue
            W_=Y[(q,p-2)];S=Y[(q-2,p-2)];E=Y[(q-2,p)];C=Y[(q-1,p-1)]
            n=1 if (W_,C,E)==(0,1,0) else 0
            Y[(q,p)]=W_^S^E^n
            if n: F.append((q-1,p-1))
    L=mx+15
    def at(q,p): return Y[(q,p)]
    g=[at(M-1 if (M-1-p)%2==0 else M-2,p) for p in range(0,L)]
    o=[at(q,M-1 if (M-1-q)%2==0 else M-2) for q in range(0,L)]
    r=(rep(g),rep(o))
    return (r,F) if fires else r
def word(pt,lo,hi): return ''.join(str(val(pt,k)) for k in range(lo,hi))
def inphi(pt):
    lam,ks=pt; hi=(max(ks) if ks else 0)+4
    return '010' not in word(pt,-3,hi)
def T(X,e,impl=direct,start=None):
    """single bare kink probe, relative parity e w.r.t. first kink of X. X kinks must be >=2+..; returns Gamma"""
    lam,ks=X
    a=2+((ks[0]+e-2)%2)   # a ≡ b1+e, a>=2
    G,O=impl((lam,[a]),X)
    return G
def TO(X,e,impl=direct):
    lam,ks=X; a=2+((ks[0]+e-2)%2)
    return impl((lam,[a]),X)
def mu(X,N=None):
    lam,ks=X; rt=lam^(len(ks)&1)
    if N is None: N=max(ks)+ min(ks)
    # xi'(p)=xi(N-1-p): kink at b (xi(b-1)!=xi(b)) -> xi'(N-b)!=xi'(N-1-b)... kink at N-b
    return (rt,sorted(N-b for b in ks))
def norm(X,start=10):
    lam,ks=X; s=ks[0]-start; return (lam,[k-s for k in ks])
def fromword(w,lam,start=10):
    # object: lam^inf w (last letter)^inf ; w[0] placed at start
    ks=[];v=lam
    for k,ch in enumerate(w):
        if int(ch)!=v: ks.append(start+k); v=int(ch)
    return (lam,ks)
def show(X):
    lam,ks=X
    if not ks: return str(lam)+'^inf'
    w=word(X,ks[0],ks[-1]+1)
    return f"{lam}^ {w} {val(X,ks[-1])}^"
