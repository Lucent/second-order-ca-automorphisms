"""Rule 59R, Lemma 6.13: on all words of length 4 to 9, the distances at which two walls of the families A, B
occur (only 2 and 3 among the distances below 4); on all words of length at most 8, a flip of one wall
preserves the set of walls with their families, leaves the bit beta of every other wall unchanged (H4) and
commutes with the flip of every other wall (H2). Prints the distances and the number of failures.
"""
import itertools
from hw59_table import XA,XB,iota,beta
SA=set(XA()); SB=set(XB())
def walls(c,fams=('A','B')):
    out={}
    for i in range(len(c)-3):
        w=tuple(c[i:i+4])
        if 'A' in fams and w in SA: out[i]='A'
        if 'B' in fams and w in SB: out[i]='B'
    return out
def flip(c,i,X):
    c=list(c); w=tuple(c[i:i+4]); c[i:i+4]=iota(w,1 if X=='A' else 2); return tuple(c)
dists={}
bad=0
for n in range(4,10):
  for c in itertools.product(range(4),repeat=n):
    ws=walls(c)
    if len(ws)<2: continue
    pos=sorted(ws)
    for i in pos:
      for i2 in pos:
        if i2<=i: continue
        dists.setdefault((ws[i],ws[i2]),set()).add(i2-i)
    if n<=8:
     for i in pos:
      ci=flip(c,i,ws[i]); wsi=walls(ci)
      if wsi!=ws: bad+=1; print("walls changed",c,i)
      for i2 in pos:
        if i2==i: continue
        # bit and core unchanged
        if beta(ws[i2],tuple(ci[i2:i2+4]))!=beta(ws[i2],tuple(c[i2:i2+4])): print("H4 bit fail",c,i,i2); bad+=1
        if flip(flip(c,i,ws[i]),i2,ws[i2])!=flip(flip(c,i2,ws[i2]),i,ws[i]): print("H2 comm fail",c,i,i2); bad+=1
print({k:sorted(v if max(v)<4 else [x for x in v if x<4]) for k,v in dists.items()})
print("bad",bad); import sys
