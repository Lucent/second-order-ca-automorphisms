"""Rule 59R, the classes of walls used in Proposition 6.16: on all words of length 9, a flip of any wall
preserves, for the wall at position 3, its family and its two class bits (level: the cell after the wall
repeats its last cell; pp: the protector side), and all 8 classes (family, level, pp) occur.
"""
import itertools, sys
sys.argv=['x']
from hw59_table import XA,XB,iota
SA=set(XA()); SB=set(XB())
def walls(c):
    out={}
    for i in range(len(c)-3):
        w=tuple(c[i:i+4])
        if w in SA: out[i]='A'
        elif w in SB: out[i]='B'
    return out
def cls(c,i):
    lev = (i+4<len(c)) and c[i+4]==c[i+3]
    pp = i>=1 and c[i]==c[i+1] and (c[i-1]==c[i+1] or c[i-1]==3)
    return lev,pp
def flip(c,i,X):
    c=list(c); c[i:i+4]=iota(tuple(c[i:i+4]),1 if X=='A' else 2); return tuple(c)
bad=0; combos=set()
for c in itertools.product(range(4),repeat=9):
    ws=walls(c)
    if ws.get(3) is None: continue
    combos.add((ws[3],)+cls(c,3))
    for j,X in ws.items():
        c2=flip(c,j,X)
        if walls(c2).get(3)!=ws[3] or cls(c2,3)!=cls(c,3): bad+=1; print(c,j)
print("class changes",bad); print(sorted(combos))
