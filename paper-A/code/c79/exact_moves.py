# Exact implementation of the moves T_m for rule 93 (edge r=0) by direct simulation of F^{\pm1} on finite words
# with an explicit filling to the left of the window (library for check_graphs_exact.py).
import random
def f(a,b,c): return (93>>(4*a+2*b+c))&1
def V(t,i): return [i&1,i&1,0,1-(i&1),1-(i&1),0][t%6]
def vcell(w,i): return (V(w-1,i),V(w,i))
def Fw(w):  # list of (x,y) -> F, length-2
    return [(w[i][1], w[i][0]^f(w[i-1][1],w[i][1],w[i+1][1])) for i in range(1,len(w)-1)]
def Bw(w):
    return [(w[i][1]^f(w[i-1][0],w[i][0],w[i+1][0]), w[i][0]) for i in range(1,len(w)-1)]
def power(w,m):
    for _ in range(abs(m)): w=Fw(w) if m>0 else Bw(w)
    return w
PAD=40; RIGHT=40
def full(state,pad):
    """state=(om,cells tuple on [-L,0]); pad: list of cells for columns < -L (length PAD).
    returns word on [-L-PAD, RIGHT] and offset"""
    om,cells=state
    L=len(cells)-1
    w=list(pad)+list(cells)+[vcell(om,i) for i in range(1,RIGHT+1)]
    return w, -L-PAD
def ycol(state,pad,t):
    w,off=full(state,pad)
    if t>=0: w2=power(w,t); return w2[-off-t][1]
    else: w2=power(w,t+1); return w2[-off+t+1][0]
def T(state,m,pad):
    """returns new state (om', cells on [-L,0]) or None if no move"""
    om,cells=state; L=len(cells)-1
    th=om%6
    if (m+th)%6==2:
        z=ycol(state,pad,m-1); om2=(om+2)%6 if z==0 else (om-1)%6
    elif (m+th)%6==0:
        z=ycol(state,pad,m); om2=(om-2)%6 if z==0 else (om+1)%6
    else: return None
    w,off=full(state,pad)
    a=abs(m)
    wm=power(w,m)            # columns off+a .. RIGHT-a
    o2=off+a
    wm=[wm[i-o2] if i<=0 else vcell(om2+m,i) for i in range(o2,RIGHT-a+1)]
    wb=power(wm,-m)          # columns off+2a ..
    o3=off+2*a
    new=tuple(wb[i-o3] for i in range(-L,1))
    # checks: right part is vacuum om2, edge not vacuum
    assert all(wb[i-o3]==vcell(om2,i) for i in range(1,RIGHT-3*a))
    assert new[-1]!=vcell(om2,0)
    return (om2,new)
def allstates(L):
    import itertools
    A=[(0,0),(0,1),(1,0),(1,1)]
    for om in range(6):
        for cs in itertools.product(A,repeat=L+1):
            if cs[-1]!=vcell(om,0): yield (om,cs)
def randpad(rng): return [(rng.randint(0,1),rng.randint(0,1)) for _ in range(PAD)]
