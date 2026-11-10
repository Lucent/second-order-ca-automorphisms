import numpy as np
def f(a,b,c): return (155>>(4*a+2*b+c))&1
def N(a,b,c): return f(a,b,c)^a^c
NT=np.array([N((k>>2)&1,(k>>1)&1,k&1) for k in range(8)],np.uint8)
FT=np.array([f((k>>2)&1,(k>>1)&1,k&1) for k in range(8)],np.uint8)
assert [k for k in range(8) if NT[k]]==[0,6,7]
def diagram(Ym1,Y0,T):
    """rows -1..T on a ring (caller keeps away from seam). returns dict t->row"""
    rows={-1:Ym1.astype(np.uint8).copy(),0:Y0.astype(np.uint8).copy()}
    for t in range(0,T):
        y=rows[t]; k=4*np.roll(y,1)+2*y+np.roll(y,-1)
        rows[t+1]=rows[t-1]^FT[k]
    return rows
def back(Y0,Y1,T):
    rows={1:Y1.copy(),0:Y0.copy()}
    for t in range(0,-T,-1):
        y=rows[t]; k=4*np.roll(y,1)+2*y+np.roll(y,-1)
        rows[t-1]=rows[t+1]^FT[k]
    return rows
def fires(y):
    k=4*np.roll(y,1)+2*y+np.roll(y,-1); return NT[k]
def phi_extend_right(u,n,rng):
    u=list(u)
    while len(u)<n:
        if u[-1]==1: u.append(0)
        elif u[-2]==0: u.append(1)
        else: u.append(rng.randint(0,1))
    return u
