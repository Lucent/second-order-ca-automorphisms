"""For each rule given (usage: hw_flipall.py 127,179,223) and each family found by hw_families.py: the flip of
all walls as a block map of radius 2 commutes with F, its square, composed exactly at radius 4 on all 4^9
windows, is the identity, and it is injective by the pair-graph test (Proposition 6.19 and Corollary 6.3).
Also prints the number of values of each column of the window along the orbit.
"""
import sys, numpy as np, itertools

from centralizer import make_rule, index_word, BlockMap, injective_pairgraph_np
from classify_r4 import check_commutes_np_R
import hw_families as TT
for rule in map(int,sys.argv[1].split(',')):
    res,Fl=TT.families(rule); f,Fl,Fil=make_rule(rule)
    print('RULE',rule, 'f truth', [f(a,b,c) for a in (0,1) for b in (0,1) for c in (0,1)])
    for iota in res:
        W=sorted(iota); ws={w:index_word(w,5) for w in W}
        tab=[index_word(e,5)[2] for e in range(1024)]
        for w,v in iota.items(): tab[w]=index_word(v,5)[2]
        g=BlockMap(2,tab)
        sq=g.compose(g,4); idok=all(sq.table[e]==index_word(e,9)[4] for e in range(4**9))
        print('  family',[''.join(map(str,ws[w]))+'<->'+str(index_word(iota[w],5)[2]) for w in W])
        print('  commutes',check_commutes_np_R(tab,Fl,2),'tau^2=id',idok,'pairgraph injective',injective_pairgraph_np(g))
        # set evolution per column
        S={tuple(ws[W[0]])}; cols=[]
        for t in range(10):
            cols.append(''.join(str(len({u[k] for u in S})) for k in range(5)))
            S={tuple(int(Fl(*x[i-1:i+2])) for i in range(1,6)) for u in S for a in range(4) for b in range(4) for x in [(a,)+u+(b,)]}
        print('  #values per column t=0..9:',cols)
        # orbit of the deterministic columns
        S={tuple(ws[W[0]])}; seq=[]
        for t in range(13):
            seq.append(''.join(str(next(iter({u[k] for u in S}))) if len({u[k] for u in S})==1 else '*' for k in range(5)))
            S={tuple(int(Fl(*x[i-1:i+2])) for i in range(1,6)) for u in S for a in range(4) for b in range(4) for x in [(a,)+u+(b,)]}
        print('  orbit',seq)
