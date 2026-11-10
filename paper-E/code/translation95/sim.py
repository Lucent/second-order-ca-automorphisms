import ctypes, numpy as np, os
_l=ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)),'libgs.so'))
u8=np.ctypeslib.ndpointer(np.uint8,flags='C'); i32=np.ctypeslib.ndpointer(np.int32,flags='C')
_l.gsolve.argtypes=[ctypes.c_int,u8,u8,ctypes.c_int,u8,u8,i32,i32,ctypes.c_int]; _l.gsolve.restype=ctypes.c_int
MASK={94:1<<2,95:(1<<0)|(1<<2),91:1}
# a "point" is given by kink representation: (bg, kinks) meaning value(k)=bg ^ #{kinks<=k} mod 2 ; for q/p coordinate
def val(pt,k):
    bg,ks=pt; return (bg+sum(1 for x in ks if x<=k))&1
def gs(rule,pi,xi,M=None,retfire=False):
    """pi,xi: kink reps in coords q,p (corner [0,0]); require val(pi,0..1)==val(xi,0..1) const bg and kinks>=2.
    returns (gamma kinks rep, o kinks rep) as reps in p and q coordinates"""
    ka=pi[1]; kb=xi[1]
    if M is None: M=2*(max([0]+ka)+max([0]+kb))+40
    A=np.array([val(pi,k) for k in range(M)],np.uint8); B=np.array([val(xi,k) for k in range(M)],np.uint8)
    assert A[0]==B[0] and A[1]==B[1]
    G=np.zeros(M,np.uint8);O=np.zeros(M,np.uint8); mf=100000
    fq=np.zeros(mf,np.int32);fp=np.zeros(mf,np.int32)
    n=_l.gsolve(MASK[rule],A,B,M,G,O,fq,fp,mf)
    def rep(S):
        S=list(S[:M-8]); return (S[0],[k for k in range(1,len(S)) if S[k]!=S[k-1]])
    r=(rep(G),rep(O))
    if retfire: return r,list(zip(fq[:n].tolist(),fp[:n].tolist()))
    return r
def direct(rule,pi,xi,M=None):
    """direct (t,i) simulation of the second-order CA; same interface."""
    ka=pi[1]; kb=xi[1]
    if M is None: M=2*(max([0]+ka)+max([0]+kb))+40
    T=M  # run long enough: firings have q,p< ~M/2
    tbl=[(rule>>k)&1 for k in range(8)]
    Lc=-(T+M+4); Rc=T+M+4
    def init(t,i):
        q=t-i; p=t+i
        return val(pi,q) if i<0 else val(xi,p)
    I=range(Lc,Rc+1)
    prev=[init(-1,i) for i in I]; cur=[init(0,i) for i in I]
    n=len(cur)
    for t in range(T):
        nb=prev[:]
        for j in range(1,n-1): nb[j]=prev[j]^tbl[4*cur[j-1]+2*cur[j]+cur[j+1]]
        prev,cur=cur,nb
    # gamma(p)=Y(T,p-T) for p small rel M ; o(q)=Y(T,T-q)
    G=[cur[(p-T)-Lc] for p in range(M-8)]; O=[cur[(T-q)-Lc] for q in range(M-8)]
    def rep(S): return (S[0],[k for k in range(1,len(S)) if S[k]!=S[k-1]])
    return rep(G),rep(O)
