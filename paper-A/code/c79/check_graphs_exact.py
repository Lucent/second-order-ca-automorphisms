# Confirmation (not used in the proof): the move graphs Gamma_3, Gamma_4 built by exact simulation
# (exact_moves.py: two random fillings to the left of the window, required to give the same result)
# are identical to the graphs computed by the three-valued engine (engine.py, via graphs.build).
import random
import exact_moves as X
from graphs import build
def exact_graph(L):
    rng=random.Random(L)
    pads=[X.randpad(rng) for _ in range(2)]
    G={}
    for s in X.allstates(L):
        out={}
        for m in range(-L,L+1):
            rs=[X.T(s,m,p) for p in pads]
            assert rs[0]==rs[1]
            if rs[0] is not None:
                assert rs[0][1][0]==s[1][0]
                out[m]=rs[0]
        G[s]=out
    return G
for L in (3,4):
    G=exact_graph(L)
    valid,edges,Z,allc,K=build(L,-L,L)
    def dec(i):
        tau,code=divmod(int(i),K); cs=allc[code]
        return (tau,tuple(((int(c)>>1),int(c)&1) for c in cs))
    n=0
    for u,kind,src,dst in edges:
        for s,d in zip(src,dst):
            a=dec(s); b=dec(d)
            assert G[a][u]==b,(L,u,a,b,G[a].get(u)); n+=1
    assert n==sum(len(x) for x in G.values())
    print('L',L,'engine graph == exact-simulation graph,',n,'edges')
