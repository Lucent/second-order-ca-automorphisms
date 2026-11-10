# Confirmation (not used in the proof): Proposition (F*_3), the orbits of the moves at level 3, by breadth-first
# search on Gamma_3.
# Gamma_3 is the move graph of graphs.py / engine.py: 1152 level-3 states, the moves at the times |m|<=3.
# A state is (theta, cells on [r-3, r]), a cell being 2x+y; index 0 is the exterior cell, 1 the tip.
# Checked: for every state xi, S(xi) \ E  is contained in R(xi), and R(xi) is contained in S(xi), where R(xi) is
# the set of states reachable from xi and S(xi) the set of states with the exterior cell of xi satisfying the
# condition on (1,1) (if the exterior cell is (0,0): tip (1,1) iff the tip of xi is (1,1)).
# Also checked: no state outside E reaches a state in E.
import sys
from graphs import build
L=3
valid,edges,Z,allc,K=build(L,-L,L)
ids=[i for i in range(6*K) if valid[i]]
adj={i:set() for i in ids}
for u,kind,src,dst in edges:
    for s,d in zip(src.tolist(),dst.tolist()): adj[s].add(d)
nedges=sum(len(src) for _,_,src,_ in edges)
def theta(i): return i//K
def cells(i): return [int(c) for c in allc[i%K]]
EXC={(0,3),(1,3),(2,3),(3,2),(5,1)}        # (theta, b_r) with b_{r-2}=(0,0), b_{r-1}=(1,1)
def inE(i):
    c=cells(i); return c[1]==0 and c[2]==3 and (theta(i),c[3]) in EXC
def key(i):
    c=cells(i); return (c[0], c[1]==3 if c[0]==0 else None)
def reach(s):
    seen={s}; st=[s]
    while st:
        x=st.pop()
        for y in adj[x]:
            if y not in seen: seen.add(y); st.append(y)
    return seen
S={}
for i in ids: S.setdefault(key(i),set()).add(i)
E=[i for i in ids if inE(i)]
print('Gamma_3: states',len(ids),'edges',nedges,'; states of E',len(E),
      'per exterior cell',sorted({e:sum(1 for i in E if cells(i)[0]==e) for e in range(4)}.items()))
assert len(ids)==1152 and nedges==2688 and len(E)==20
ok=True; sharp=True
for i in ids:
    R=reach(i); Si=S[key(i)]
    if not ({x for x in Si if not inE(x)}<=R<=Si): ok=False; print('FAILURE at state',theta(i),cells(i))
    if not inE(i) and any(inE(x) for x in R): sharp=False
print('S(xi)\\E <= R(xi) <= S(xi) for all states:',ok)
print('not used: no state outside E reaches E:',sharp)
if not ok: print('FAILURE'); sys.exit(1)
print('(F*_3) holds on Gamma_3:',ok)
