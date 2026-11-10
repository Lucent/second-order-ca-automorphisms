# Second three-valued move engine for rule 93 at the edge r=0, written independently of engine.py, with states
# (theta, tuple of cells (x,y) on [-L,0]); U=2 is the unknown value. Used by the confirmations only.
import itertools
U=2
def f2(a,b,c): return (93>>(4*a+2*b+c))&1
def f(a,b,c):
    vals={f2(x,y,z) for x in ((0,1) if a==U else (a,)) for y in ((0,1) if b==U else (b,)) for z in ((0,1) if c==U else (c,))}
    return vals.pop() if len(vals)==1 else U
def xor(a,b): return U if U in (a,b) else a^b
def V(t,i): return [i&1,i&1,0,1-(i&1),1-(i&1),0][t%6]
def vcell(w,i): return (V(w-1,i),V(w,i))
W=14
class Diag:
    def __init__(s,th,rows,TR,base=None,cone=None):
        # rows: dict {(t0-1,i),(t0,i)} values for i in [-L-W..0]; build all rows in [-TR-1,TR]
        s.th=th; s.Y=dict(rows); s.lo=-W-12
        s.base=base; s.cone=cone
        ts=sorted({t for t,_ in rows}); t0=ts[1]
        for t in range(t0,TR):     # forward
            for i in range(s.lo,1):
                if (t+1,i) in s.Y: continue
                s.Y[(t+1,i)]=s._calc(t+1,i,t,t-1)
        for t in range(t0-1,-TR-1,-1):
            for i in range(s.lo,1):
                if (t-1,i) in s.Y: continue
                s.Y[(t-1,i)]=s._calc(t-1,i,t,t+1)
    def _calc(s,tn,i,t,tp):
        if s.cone is not None and not s.cone(tn,i): return s.base.g(tn,i)
        return xor(s.g(tp,i),f(s.g(t,i-1),s.g(t,i),s.g(t,i+1)))
    def g(s,t,i):
        if i>0: return V(t+s.th,i)
        if i<s.lo: return U
        return s.Y.get((t,i),U)
    def z(s,t): return s.g(t,0)
def diag(state,TR):
    th,cells=state; L=len(cells)-1
    rows={}
    for i in range(-W-12,1):
        if i>=-L: x,y=cells[i+L]
        else: x=y=U
        rows[(-1,i)]=x; rows[(0,i)]=y
    return Diag(th,rows,TR)
def move(state,m,TR=None,D=None):
    th,cells=state; L=len(cells)-1
    if TR is None: TR=max(L,abs(m))+4
    if D is None: D=diag(state,TR+abs(m)+2)
    n=(m+th)%6
    if n==2:
        z=D.z(m-1); assert z!=U; th2=(th+2)%6 if z==0 else (th-1)%6
    elif n==0:
        z=D.z(m); assert z!=U; th2=(th-2)%6 if z==0 else (th+1)%6
    else: return None
    rows={(t,i):D.g(t,i) for t in (m-1,m) for i in range(D.lo,1)}
    if m>=1: cone=lambda t,i: not (t<=m-2 and i<=-(m-1-t))
    else: cone=lambda t,i: not (t>=m+1 and i<=-(t-m))
    D2=Diag(th2,rows,TR+abs(m)+2,base=D,cone=cone)
    new=tuple((D2.g(-1,i),D2.g(0,i)) for i in range(-L,1))
    assert all(U not in c for c in new), (state,m,new)
    assert new[-1]!=vcell(th2,0)
    return (th2,new)
def states(L):
    A=[(0,0),(0,1),(1,0),(1,1)]
    for th in range(6):
        for cs in itertools.product(A,repeat=L+1):
            if cs[-1]!=vcell(th,0): yield (th,cs)
