"""Light-speed window graph Lambda_M built directly from the local rule:
vertices = words w_0..w_M over A (cell = 2x+y), edge w -> F(w r1 r2) restricted to sites 1..M+1.
Reports the number of strongly connected components (Tarjan, iterative) for M <= Mmax.
Also: for the linear part, the controllability rank of the window system (dimension of the space of
windows reachable from the zero window), for M up to Mlin."""
import sys, numpy as np
sys.setrecursionlimit(10000)
def f_of(rule): return lambda a,b,c: (rule >> (4*a+2*b+c)) & 1
def Floc(f, s0, s1, s2):
    x, y = s1 >> 1, s1 & 1
    return (y << 1) | (x ^ f(s0 & 1, y, s2 & 1))
def scc_count(n, succ):
    index=[-1]*n; low=[0]*n; onst=[False]*n; st=[]; idx=0; comps=0
    for v0 in range(n):
        if index[v0]!=-1: continue
        work=[(v0,0)]
        while work:
            v,i=work.pop()
            if i==0:
                index[v]=low[v]=idx; idx+=1; st.append(v); onst[v]=True
            recurse=False
            for j in range(i,len(succ[v])):
                w=succ[v][j]
                if index[w]==-1:
                    work.append((v,j+1)); work.append((w,0)); recurse=True; break
                elif onst[w]: low[v]=min(low[v],index[w])
            if recurse: continue
            if low[v]==index[v]:
                comps+=1
                while True:
                    w=st.pop(); onst[w]=False
                    if w==v: break
            if work:
                u,_=work[-1]; low[u]=min(low[u],low[v])
    return comps
def lam(rule, M):
    f=f_of(rule); L=M+1; n=4**L
    succ=[]
    for e in range(n):
        w=[(e>>(2*(L-1-k)))&3 for k in range(L)]
        out=set()
        for r1 in range(4):
            for r2 in range(4):
                ww=w+[r1,r2]
                new=[Floc(f,ww[k],ww[k+1],ww[k+2]) for k in range(L)]  # sites 1..M+1
                out.add(sum(s<<(2*(L-1-k)) for k,s in enumerate(new)))
        succ.append(sorted(out))
    return scc_count(n,succ)
def ctrl_rank(rule, M, T=None):
    """linear part: state = window c_0..c_M (bits x_k,y_k), zero initial window, free cells right of M.
    Reachable subspace after T steps = span of images of unit free inputs; computed by propagating unit
    perturbations of the free tail with the linear rule (constant term dropped)."""
    f=f_of(rule); k0=f(0,0,0); al=f(1,0,0)^k0; be=f(0,1,0)^k0; ga=f(0,0,1)^k0
    if T is None: T=2*(M+1)+2
    W=M+1+2*T+2
    rows=[]
    for j in range(M+1, W):
        for bit in (0,1):
            x=np.zeros(W,dtype=np.uint8); y=np.zeros(W,dtype=np.uint8)
            (x if bit==0 else y)[j]=1
            for t in range(T):
                ym=np.roll(y,1); ym[0]=0; yp=np.roll(y,-1); yp[-1]=0
                x,y=y,(x^(al*ym)^(be*y)^(ga*yp))&1
            rows.append(np.concatenate([x[T:T+M+1],y[T:T+M+1]]))
    A=np.array(rows,dtype=np.uint8)
    # rank over F2
    A=A.copy(); r=0
    for c in range(A.shape[1]):
        piv=[i for i in range(r,A.shape[0]) if A[i,c]]
        if not piv: continue
        A[[r,piv[0]]]=A[[piv[0],r]]
        for i in range(A.shape[0]):
            if i!=r and A[i,c]: A[i]^=A[r]
        r+=1
    return r, 2*(M+1)
if __name__=="__main__":
    if len(sys.argv)==3:      # window_graph.py RULE MMAX: the numbers of components only
        rule,Mmax=int(sys.argv[1]),int(sys.argv[2])
        print(rule, "SCCs of Lambda_M, M=0..%d:"%Mmax, [lam(rule,M) for M in range(Mmax+1)], flush=True); sys.exit()
    rules=[240,15,60,195,90,165,150,105]
    for rule in rules:
        print(rule, "SCCs of Lambda_M, M=0..5:", [lam(rule,M) for M in range(6)], flush=True)
    for rule in rules:
        print(rule, "controllability rank/dim, M=0..30:", [ "%d/%d"%ctrl_rank(rule,M) for M in (0,1,2,3,5,8,12,20,30)], flush=True)
