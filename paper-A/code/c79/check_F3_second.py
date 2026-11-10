# Confirmation (not used in the proof): the search of c_F3.py (Proposition (F*_3)) repeated on Gamma_3 built twice more,
# by the second engine engine2.py and by exact simulation (exact2.py) of whole configurations with three random
# left parts per state and move, which are required to give the same result; and the two graphs are compared.
import random, itertools, sys
from graphs2 import build
from exact2 import Conf
from engine2 import vcell
L=3
def exact_graph(L):
    rng=random.Random(7); A=[(0,0),(0,1),(1,0),(1,1)]; G={}
    for th in range(6):
        for cs in itertools.product(A,repeat=L+1):
            if cs[-1]==vcell(th,0): continue
            s=(th,cs); G[s]={}
            for m in range(-L,L+1):
                res=set()
                for k in range(3):
                    left=[(rng.randint(0,1),rng.randint(0,1)) for _ in range(12)]
                    d=Conf(th,left+list(cs),PADR=30).move(m)
                    res.add(None if d is None else (d.th,tuple(d.w[d.n+i] for i in range(-L,1))))
                assert len(res)==1,(s,m,res)
                v=res.pop()
                if v is not None: G[s][m]=v
    return G
EXC={(0,(1,1)),(1,(1,1)),(2,(1,1)),(3,(1,0)),(5,(0,1))}
def inE(s):
    th,cs=s; return cs[1]==(0,0) and cs[2]==(1,1) and (th,cs[3]) in EXC
def key(s):
    th,cs=s; return (cs[0],(cs[1]==(1,1)) if cs[0]==(0,0) else None)
def reach(G,s):
    seen={s}; st=[s]
    while st:
        x=st.pop()
        for y in G[x].values():
            if y not in seen: seen.add(y); st.append(y)
    return seen
S2,Z2,M2=build(L); X=exact_graph(L)
same=set(X)==set(S2) and all(X[s]==M2[s] for s in S2)
print('Gamma_3 by engine2.py == Gamma_3 by exact simulation:',same,'edges',sum(len(v) for v in X.values()))
ok=same
for name,G in (('engine2',M2),('exact',X)):
    cls={}
    for s in G: cls.setdefault(key(s),set()).add(s)
    fs=all({x for x in cls[key(s)] if not inE(x)}<=reach(G,s)<=cls[key(s)] for s in G)
    sharp=all(not any(inE(x) for x in reach(G,s)) for s in G if not inE(s))
    print(name,'states',len(G),'(F*_3):',fs,' no state outside E reaches E:',sharp)
    ok&=fs and sharp
if not ok: print('FAILURE'); sys.exit(1)
print('(F*_3) holds on both graphs: True')
