# Independent exact engine for rule 93 second-order CA at a right edge r=0.
# cell = (x,y); config in Omega_0 = (theta, list of cells on sites -W..0), vacuum of phase theta to the right.
import itertools, random
def f(a,b,c): return (1-c) | ((1-a)&b)
assert [f(a,b,c) for a in (0,1) for b in (0,1) for c in (0,1)]==[1,0,1,1,1,0,1,0]
def V(t,i):
    t%=6
    if t in (0,1): return i%2
    if t in (2,5): return 0
    return 1-(i%2)
def vcell(w,i): return (V(w-1,i),V(w,i))
def Fw(cells,lo,w):
    """cells on sites lo..0 (list of (x,y)), vacuum phase w to the right. returns F(c) on sites lo+1..0."""
    ext=cells+[vcell(w,1)]
    out=[]
    for k in range(1,len(cells)):
        x,y=ext[k]
        out.append((y, x^f(ext[k-1][1],y,ext[k+1][1])))
    return out
def Finvw(cells,lo,w):
    ext=cells+[vcell(w,1)]
    out=[]
    for k in range(1,len(cells)):
        x,y=ext[k]
        out.append((y^f(ext[k-1][0],x,ext[k+1][0]), x))
    return out
def Fpow(cells,w,m):
    """returns (cells of F^m c on sites lo+|m|..0, phase w+m)"""
    c=list(cells)
    for _ in range(abs(m)):
        if m>0: c=Fw(c,None,w); w=(w+1)%6
        else: c=Finvw(c,None,w); w=(w-1)%6
    return c,w
def zcol(cells,w,ts):
    """edge column values z_t for t in ts"""
    res={}
    for t in ts:
        if t>=0:
            c,_=Fpow(cells,w,t); res[t]=c[-1][1]
        else:
            c,_=Fpow(cells,w,t+1); res[t]=c[-1][0]
    return res
def move(cells,w,m):
    """T_m on (cells,w). returns (newcells,neww,kind) or None. cells list on sites -W..0; result same length."""
    n=(m+w)%6
    if n not in (0,2): return None
    c,wm=Fpow(cells,w,m)
    zm1,zm=c[-1]
    if n==2:
        d=2 if zm1==0 else -1; kind='F+' if zm1==0 else 'F-'
    else:
        d=-2 if zm==0 else 1; kind='B-' if zm==0 else 'B+'
    w2=(wm+d)%6
    assert c[-1]!=vcell(w2,0)
    c2,w3=Fpow(c,w2,-m)
    assert w3==(w+d)%6
    lost=len(cells)-len(c2)
    new=list(cells[:lost])+c2
    return new,w3,kind
def states(L):
    for w in range(6):
        for cs in itertools.product([(0,0),(0,1),(1,0),(1,1)],repeat=L+1):
            if cs[-1]!=vcell(w,0): yield (w,cs)
def graph(L,times=None,pad=14,seeds=(1,2)):
    """edges[(state)] = dict m -> (newstate,kind)."""
    if times is None: times=range(-L,L+1)
    G={}
    for (w,cs) in states(L):
        outs=None
        for sd in seeds:
            rnd=random.Random(hash((sd,w,cs))&0xffffffff)
            left=[(rnd.randint(0,1),rnd.randint(0,1)) for _ in range(pad)]
            full=left+list(cs)
            o={}
            for m in times:
                r=move(full,w,m)
                if r is None: continue
                nc,nw,kind=r
                assert nc[:pad+1-0]==full[:pad+1] if False else True
                assert nc[:pad]==left and nc[pad]==cs[0], 'exterior changed'
                o[m]=((nw,tuple(nc[pad:])),kind)
            if outs is None: outs=o
            else: assert outs==o, ('depends on left part',w,cs)
        G[(w,cs)]=outs
    return G
def scc(G):
    # Tarjan iterative; G: node -> iterable of successors
    index={}; low={}; onst=set(); st=[]; comps=[]; idx=[0]
    for root in G:
        if root in index: continue
        work=[(root,iter(G[root]))]
        index[root]=low[root]=idx[0]; idx[0]+=1; st.append(root); onst.add(root)
        while work:
            v,it=work[-1]
            adv=False
            for u in it:
                if u not in index:
                    index[u]=low[u]=idx[0]; idx[0]+=1; st.append(u); onst.add(u)
                    work.append((u,iter(G[u]))); adv=True; break
                elif u in onst: low[v]=min(low[v],index[u])
            if adv: continue
            work.pop()
            if work: low[work[-1][0]]=min(low[work[-1][0]],low[v])
            if low[v]==index[v]:
                comp=[]
                while True:
                    u=st.pop(); onst.discard(u); comp.append(u)
                    if u==v: break
                comps.append(comp)
    return comps
def cstr(c): return '%d%d'%c
def sstr(s): return '%d|'%s[0]+' '.join(cstr(c) for c in s[1])
