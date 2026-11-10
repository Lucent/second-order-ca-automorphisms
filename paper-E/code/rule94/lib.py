# Rule 94R scattering tools based on the Goursat (light-cone) recurrence.
# A word-object: (lam, core) meaning lam^inf core tail^inf, tail=core[-1] (or lam if core empty); core[0]!=lam.
# Goursat solver on quadrant q,p>=0 with boundary A(q), B(p) (lists), returns dict Y and firings.
def solve(A,B):
    M=len(A); Mp=len(B)
    Y={}
    for q in range(M): Y[(q,q%2)]=A[q]
    for p in range(Mp): Y[(p%2,p)]=B[p]
    assert A[0]==B[0] and A[1]==B[1]
    F=[]
    for q in range(2,M):
        for p in range(2+(q%2),Mp,2):
            W=Y[(q,p-2)];S=Y[(q-2,p-2)];E=Y[(q-2,p)];C=Y[(q-1,p-1)]
            n=1 if (W==0 and C==1 and E==0) else 0
            Y[(q,p)]=W^S^E^n
            if n: F.append((q-1,p-1))
    return Y,F
def scatter(pw,xw,pad=None,full=False):
    """pw: probe word (string) with pw[0]=pw[1]=lam, xw: object word with xw[0]=xw[1]=lam.
    Both are extended by their last letter. Returns (gamma word, omega word, firings)."""
    if pad is None: pad=3*(len(pw)+len(xw))+20
    A=[int(c) for c in pw]; B=[int(c) for c in xw]
    la=len(A); lb=len(B)
    A=A+[A[-1]]*pad; B=B+[B[-1]]*pad
    Y,F=solve(A,B)
    M=len(A); Mp=len(B)
    g=[Y[(M-1 if (M-1-p)%2==0 else M-2,p)] for p in range(Mp)]
    o=[Y[(q,Mp-1 if (Mp-1-q)%2==0 else Mp-2)] for q in range(M)]
    # check stabilization: firings far from the far edges
    mq=max([f[0] for f in F],default=0); mp=max([f[1] for f in F],default=0)
    assert mq<M-8 and mp<Mp-8,(mq,mp,M,Mp)
    gs=''.join(map(str,g)); os=''.join(map(str,o))
    if full: return gs,os,F,Y
    return gs,os,F
def trim(s):
    """strip constant tail, keep 2 lead letters"""
    t=s.rstrip(s[-1]); return t+s[-1]*2
def core(s):
    """(lam, core) of an eventually constant word given with lam prefix"""
    lam=s[0]; t=s.lstrip(lam); 
    if not t: return lam,''
    tail=t[-1]; t=t.rstrip(tail)+tail
    return lam,t
def inphi(s): return '010' not in s
def kinks(s): return [i for i in range(1,len(s)) if s[i]!=s[i-1]]
