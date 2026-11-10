import itertools, numpy as np
from collections import defaultdict
from common import *
def ftab(rule): return np.array([(rule>>t)&1 for t in range(8)],np.int8)
def Fword(ft,W):  # W: array (...,n) of cells 0..3 -> (...,n-2)
    X=W>>1; Y=W&1
    NY=X[...,1:-1]^ft[(Y[...,:-2]<<2)|(Y[...,1:-1]<<1)|Y[...,2:]]
    return (Y[...,1:-1]<<1)|NY
def fix_graph(rule,a,b):
    """valid (2b+1)-windows w with (F^b w)_center == w[b+a]"""
    ft=ftab(rule); L=2*b+1
    allw=np.array(list(itertools.product(range(4),repeat=L)),np.int8)
    W=allw.copy()
    for _ in range(b): W=Fword(ft,W)
    ok=(W[:,0]==allw[:,b+a])
    valid=[tuple(w) for w in allw[ok]]
    succ=defaultdict(set); pred=defaultdict(set)
    for w in valid: succ[w[:-1]].add(w[1:]); pred[w[1:]].add(w[:-1])
    V=set(succ)|set(pred)
    ch=True
    while ch:
        ch=False
        for v in list(V):
            if not (succ[v]&V) or not (pred[v]&V): V.discard(v); ch=True
    return V,succ,pred
def sccs(V,succ):
    idx={};low={};st=[];on=set();out=[];c=[0]
    import sys; sys.setrecursionlimit(1000000)
    def sc(v):
        idx[v]=low[v]=c[0];c[0]+=1;st.append(v);on.add(v)
        for w in succ[v]&V:
            if w not in idx: sc(w); low[v]=min(low[v],low[w])
            elif w in on: low[v]=min(low[v],idx[w])
        if low[v]==idx[v]:
            comp=[]
            while True:
                w=st.pop();on.discard(w);comp.append(w)
                if w==v: break
            out.append(comp)
    for v in V:
        if v not in idx: sc(v)
    return out
def sccs(V,succ):
    index={};low={};onst=set();st=[];out=[];cnt=0
    for root in V:
        if root in index: continue
        work=[(root,iter(sorted(succ[root]&V)))]; index[root]=low[root]=cnt;cnt+=1;st.append(root);onst.add(root)
        while work:
            v,it=work[-1]; adv=False
            for w in it:
                if w not in index:
                    index[w]=low[w]=cnt;cnt+=1;st.append(w);onst.add(w);work.append((w,iter(sorted(succ[w]&V))));adv=True;break
                elif w in onst: low[v]=min(low[v],index[w])
            if adv: continue
            work.pop()
            if work: low[work[-1][0]]=min(low[work[-1][0]],low[v])
            if low[v]==index[v]:
                comp=[]
                while True:
                    w=st.pop();onst.discard(w);comp.append(w)
                    if w==v: break
                out.append(comp)
    return out
def kinks(rule,a,b,Bstate,maxpaths=20):
    """Bstate: (xs,ys) periodic tuple of B phase. returns list of (word, rightcycle) : word = finite cells after B prefix"""
    V,succ,pred=fix_graph(rule,a,b)
    xs,ys=Bstate; p=len(xs); cells=[2*xs[i]+ys[i] for i in range(p)]
    L=2*b
    Bv=[tuple(cells[(j+i)%p] for i in range(L)) for j in range(p)]
    if not all(v in V for v in Bv): return None
    comps=sccs(V,succ); cid={}
    for i,c in enumerate(comps):
        for v in c: cid[v]=i
    Bc=cid[Bv[0]]
    info={'nV':len(V),'nscc':len(comps),'Bscc_size':len(comps[Bc])}
    # simple-cycle SCC test: every vertex has exactly one successor inside
    def simple(ci): return all(len([w for w in succ[v]&V if cid[w]==ci])==1 for v in comps[ci])
    info['Bsimple']=simple(Bc)
    # DFS from B vertices out of B SCC, record paths until reaching another SCC, then that SCC's cycle
    res=[]; Bset=set(Bv)
    from collections import deque
    for j,v0 in enumerate(Bv):
        stack=[(v0,[])]
        while stack and len(res)<maxpaths:
            v,path=stack.pop()
            for w in sorted(succ[v]&V):
                if cid[w]==Bc:
                    if w in Bset:
                        if path: res.append({'Bphase':j,'path':path+[w[-1]],'cycle':None,'glider_end':Bv.index(w),'b':b,'scc':Bc,'scc_size':len(comps[Bc]),'scc_simple':False})
                        continue
                    np_=path+[w[-1]]
                    if len(np_)<=24: stack.append((w,np_))
                    continue
                np_=path+[w[-1]]
                if len(np_)>40: continue
                ci=cid[w]
                if len(comps[ci])>1 or w in succ[w]:   # nontrivial SCC reached
                    # extend around a cycle in ci to give periodic right tail
                    seq=[w];u=w
                    while True:
                        nx=[z for z in sorted(succ[u]&V) if cid[z]==ci][0]
                        if nx in seq: loopstart=seq.index(nx);break
                        seq.append(nx);u=nx
                    # cells after w: seq[1..] then loop seq[loopstart..] repeated (loop begins again at seq[loopstart])
                    trans=[z[-1] for z in seq[1:]]
                    loop=[z[-1] for z in seq[loopstart:]]
                    np_=np_+trans
                    cyc=None
                    res.append({'Bphase':j,'path':np_,'cycle':loop,'scc':ci,'scc_size':len(comps[ci]),'scc_simple':simple(ci)})
                else:
                    stack.append((w,np_))
    return info,res
