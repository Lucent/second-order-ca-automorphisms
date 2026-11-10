"""base bitmap = all (2R+1)-windows of bi-infinite points of Fix(sigma^a F^b), (a,b) in a list.
usage: mkbase.py rule R "a,b;a,b;..." out"""
import sys, itertools, numpy as np
sys.setrecursionlimit(10000)
r=int(sys.argv[1]); R=int(sys.argv[2]); pairs=[tuple(map(int,p.split(','))) for p in sys.argv[3].split(';')]; out=sys.argv[4]
f=lambda a,b,c:(r>>((a<<2)|(b<<1)|c))&1
def Fw(w):  # word of symbols s=(x<<1)|y -> F(word), shorter by 2
    return tuple(((w[i]&1)<<1)|((w[i]>>1)^f(w[i-1]&1,w[i]&1,w[i+1]&1)) for i in range(1,len(w)-1))
def Bw(w):
    return tuple((((w[i]&1)^f(w[i-1]>>1,w[i]>>1,w[i+1]>>1))<<1)|(w[i]>>1) for i in range(1,len(w)-1))
def windows(a,b,n):
    M=max(abs(b),abs(a)); L=2*M+1
    def ok(w):  # (F^b w)_centre == w_{centre - a}   [sigma^a F^b c = c  <=>  (F^b c)_{i+a} = c_i ... use (F^b c)_i = c_{i-a}]
        u=w
        for _ in range(abs(b)): u=Fw(u) if b>0 else Bw(u)
        k=len(u)//2
        return u[k]==w[M-a]
    edges=[w for w in itertools.product(range(4),repeat=L) if ok(w)]
    # vertices: words of length L-1; trim to bi-infinite part
    out_={};in_={}
    for w in edges: out_.setdefault(w[:-1],[]).append(w[1:]); in_.setdefault(w[1:],[]).append(w[:-1])
    V=set(out_)|set(in_)
    while True:
        V2={v for v in V if any(u in V for u in out_.get(v,[])) and any(u in V for u in in_.get(v,[]))}
        if V2==V: break
        V=V2
    S=set()
    if not V: return S
    if n<=L-1:
        for v in V:
            for s in range(L-1-n+1): S.add(v[s:s+n])
        return S
    paths={v:None for v in V}; cur=[v for v in V]
    for _ in range(n-(L-1)):
        cur=[p+(u[-1],) for p in cur for u in out_.get(p[-(L-1):],[]) if u in V]
        cur=list(set(cur))
    return set(cur)
n=2*R+1; NW=4**n; Z=np.zeros(max(NW,64),bool); tot=set()
for (a,b) in pairs:
    S=windows(a,b,n); tot|=S
for w in tot:
    X=Y=0
    for s in w: X=(X<<1)|(s>>1); Y=(Y<<1)|(s&1)
    Z[(X<<n)|Y]=True
np.packbits(Z,bitorder='little').tofile(out); print(r,R,pairs,'base',int(Z.sum()))
