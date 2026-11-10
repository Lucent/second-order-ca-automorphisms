"""Rule 59R, the walls X_A and X_B (Lemmas 6.11, 6.12 and 6.13), exhaustively: each family is invariant under
F and F^-1 in all 16 contexts; the flip iota commutes with F and F^-1 in all contexts; the cells L, C, R of
the image do not depend on the context; the map on cores; the parity bit beta; and, on all words of length
5 to 9, the overlaps of walls (distances and types), with every flip preserving the set of walls and every
two flips commuting.
"""
import numpy as np, itertools
from hw_ca import *
from centralizer import make_rule
Ft,Fi=tables(59)
f,Fl,Fil=make_rule(59)
assert all(Ft[a,b,c]==Fl(a,b,c) and Fi[a,b,c]==Fil(a,b,c) for a,b,c in itertools.product(range(4),repeat=3))
# check Fi inverts Ft on random ring
rng=np.random.default_rng(1); c=rng.integers(0,4,1000); assert (step(Fi,step(Ft,c))==c).all()
def fam(delta): return {(a,p,C,(p+delta)%3) for p in range(3) for a in (p,3) for C in ((p+delta)%3,3)}
XA,XB=fam(1),fam(2); assert len(XA)==12 and len(XB)==12 and not (XA&XB)
def iota(w,delta):
    p=w[1]; C=w[2]; alt=[(p+delta)%3,3]; return (w[0],w[1],alt[1-alt.index(C)],w[3])
for name,W,dl in [('A',XA,1),('B',XB,2)]:
    inv=True; flip=True; core=set(); Lctx=True
    for w in W:
        for l,r in itertools.product(range(4),repeat=2):
            u=(l,)+w+(r,)
            for T in (Ft,Fi):
                Fw=step_word(T,u)
                if Fw not in W: inv=False
                v=(l,)+iota(w,dl)+(r,)
                # flip lemma: need outputs at all cells; changed input cell 0 (index 3 in u) affects outputs idx 1..3 of w -> we compare whole Fw window plus outside cells
                Fv=step_word(T,v)
                if Fw in W and Fv!=iota(Fw,dl): flip=False
            # also outside cells: cell right of R depends on R only; cell left of a depends on a only -> unchanged by flip. ok.
            core.add(((w[1],w[2]),step_word(Ft,u)[1:3]))
    # context-freeness of L,C,R
    cf=all(len({step_word(Ft,(l,)+w+(r,))[1:] for l in range(4) for r in range(4)})==1 for w in W)
    print(name,'invariant F,F^-1:',inv,' flip lemma F,F^-1:',flip,' L,C,R context-free:',cf)
    print('  core map:',sorted(core))
# parity
cycB=[(0,2),(1,0),(2,1),(0,3),(1,3),(2,3)]
def beta(w,dl):
    if dl==1: return int(w[2]==3)^int(w[1]==2)
    return cycB.index((w[1],w[2]))%2
for W,dl in [(XA,1),(XB,2)]:
    ok=all(beta(iota(w,dl),dl)!=beta(w,dl) for w in W)
    ok2=all(beta(step_word(Ft,(l,)+w+(r,)),dl)!=beta(w,dl) for w in W for l in range(4) for r in range(4))
    print('parity', dl, ok, ok2)
# overlaps: exhaustive over all words of length n, union family, check that each flip preserves full labelled occurrence set and flips commute
def occs(w):
    o={}
    for i in range(len(w)-3):
        s=tuple(w[i:i+4])
        if s in XA: o[i]=1
        elif s in XB: o[i]=2
    return o
def flipat(w,i,dl):
    w=list(w); s=iota(tuple(w[i:i+4]),dl); w[i:i+4]=s; return tuple(w)
for n in (5,6,7,8,9):
    bad=0; pairs=0; dists=set()
    for w in itertools.product(range(4),repeat=n):
        o=occs(w)
        if len(o)<2: continue
        for i,di in o.items():
            wi=flipat(w,i,di)
            if occs(wi)!=o: bad+=1
            for j,dj in o.items():
                if j<=i: continue
                pairs+=1; dists.add((j-i,di,dj))
                if flipat(flipat(w,i,di),j,dj)!=flipat(flipat(w,j,dj),i,di): bad+=1
    print('len',n,'overlap failures',bad,'pairs',pairs,'distances(type)',sorted(d for d in dists if d[0]<4))
