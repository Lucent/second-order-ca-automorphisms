# Exact three-valued move engine at the right edge r=0: cells left of the window are unknown (U); results are certified known.
import numpy as np, itertools
from base import f, Vac
U=2  # unknown
FT3=np.zeros(27,np.uint8)
for a,b,c in itertools.product(range(3),repeat=3):
    vals=set(f(x,y,z) for x in ((0,1) if a==U else (a,)) for y in ((0,1) if b==U else (b,)) for z in ((0,1) if c==U else (c,)))
    FT3[a*9+b*3+c]=vals.pop() if len(vals)==1 else U
XOR3=np.array([[0,1,U],[1,0,U],[U,U,U]],np.uint8)
def vcell(phi,i): return 2*Vac(phi-1,i)+Vac(phi,i)
class Diag:
    """old diagrams for a batch of states with tail tau, time-0 cells on [-L,0]; unknown left of -L"""
    def __init__(s,tau,cells,L,trange):
        s.tau=tau; s.L=L; s.n=cells.shape[0]; s.M=L+2*trange+4
        tmin,tmax=-trange-2,trange+2
        s.Y={}
        s.unk=np.full(s.n,U,np.uint8)
        for i in range(-s.M,1):
            if i>=-L:
                s.Y[(-1,i)]=(cells[:,i+L]>>1).astype(np.uint8); s.Y[(0,i)]=(cells[:,i+L]&1).astype(np.uint8)
            else:
                s.Y[(-1,i)]=s.unk; s.Y[(0,i)]=s.unk
        for t in range(0,tmax):
            for i in range(-s.M,1):
                s.Y[(t+1,i)]=XOR3[s.g(t-1,i),FT3[s.g(t,i-1)*9+s.g(t,i)*3+s.g(t,i+1)]]
        for t in range(-1,tmin,-1):
            for i in range(-s.M,1):
                s.Y[(t-1,i)]=XOR3[s.g(t+1,i),FT3[s.g(t,i-1)*9+s.g(t,i)*3+s.g(t,i+1)]]
    def g(s,t,i):
        if i>0: return np.full(s.n,Vac(t+s.tau,i),np.uint8)
        if i<-s.M: return s.unk
        return s.Y[(t,i)]
    def z(s,t): return s.Y[(t,0)]
def move(D,u):
    """returns (kind, newtau array, newcells (n,L+1) with 3-valued entries) or None if no move at time u"""
    tau=D.tau; L=D.L
    if (u+tau)%6==2:
        z=D.z(u-1); kind='F'
        nt=np.where(z==0,(tau+2)%6,(tau-1)%6); nt[z==U]=99
    elif (u+tau)%6==0:
        z=D.z(u); kind='B'
        nt=np.where(z==0,(tau-2)%6,(tau+1)%6); nt[z==U]=99
    else: return None
    out_cells=np.zeros((D.n,L+1),np.uint8); out_tau=nt.copy()
    for tn in np.unique(nt):
        sel=np.nonzero(nt==tn)[0]
        if tn==99:
            out_cells[sel]=U*5; continue
        new={}
        def G(t,i):
            if i>0: return np.full(len(sel),Vac(t+tn,i),np.uint8)
            if (t,i) in new: return new[(t,i)]
            return D.g(t,i)[sel] if i>=-D.M else np.full(len(sel),U,np.uint8)
        if u>=1:
            for t in range(u-2,-2,-1):
                for i in range(-(u-2-t),1):
                    new[(t,i)]=XOR3[G(t+2,i),FT3[G(t+1,i-1)*9+G(t+1,i)*3+G(t+1,i+1)]]
        else:
            for t in range(u+1,1):
                for i in range(-(t-u-1),1):
                    new[(t,i)]=XOR3[G(t-2,i),FT3[G(t-1,i-1)*9+G(t-1,i)*3+G(t-1,i+1)]]
        for i in range(-L,1):
            a=G(-1,i); b=G(0,i)
            c=np.where((a==U)|(b==U),U*5,2*a+b)
            out_cells[sel,i+L]=c
    return kind,out_tau,out_cells
def all_states(L):
    allc=np.array(list(itertools.product(range(4),repeat=L+1)),np.uint8)
    return allc
