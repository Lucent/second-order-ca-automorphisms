"""Rule 27R, Theorem 6.20: (W1) for F and F^-1 and (H1) for the 12 words of W_27 in all 16 contexts; the images
of 111222; the flip of all walls, as a table of radius 3, commutes with F on all 4^9 windows; the offsets
at which two words of W_27 overlap; and, on all words of length 7 to 11, every flip preserves the set of
occurrences and every two flips commute (H2).
"""
import itertools, numpy as np, sys

from classify_r4 import check_commutes_np_R
from hw_common import Floc, Finvloc, wimg
F=Floc(27); Fi=Finvloc(27)
W=[(a,p,p,C,(p+1)%3,(p+1)%3) for p in range(3) for a in (p,3) for C in ((p+1)%3,3)]
S=set(W); assert len(S)==12
def io(w): w=list(w); w[3]=3 if w[3]!=3 else (w[1]+1)%3; return tuple(w)
ok=True
for G in (F,Fi):
  for w in W:
    for l in range(4):
      for r in range(4):
        im=wimg(G,(l,)+w+(r,))
        if im not in S: ok=False; print("W1 fail",w,l,r,im)
        if G is F and wimg(G,(l,)+io(w)+(r,))!=io(im): ok=False; print("H1 fail",w,l,r)
print("W1,H1",ok)
print("111222 images",{wimg(F,(l,1,1,1,2,2,2,r)) for l in range(4) for r in range(4)})
# flip-all radius 3 table: cell n flipped if occurrence at n-3 (reads n-3..n+2)
idx=np.arange(4**7); cells=[(idx>>(2*(6-i)))&3 for i in range(7)]
tab=cells[3].copy()
for e in range(4**7):
    w=tuple(int(cells[i][e]) for i in range(6))
    if w in S: tab[e]=io(w)[3]
print("flipall commutes",check_commutes_np_R(tab,F,3))
# check H2/H4-type: overlaps of W27, flips preserving occurrences, commuting
def occ(c): return [i for i in range(len(c)-5) if tuple(c[i:i+6]) in S]
def fl(c,i): c=list(c); c[i:i+6]=io(tuple(c[i:i+6])); return tuple(c)
offs=set(); bad=0
for a in W:
  for b in W:
    for d in range(1,6):
      if a[d:]==b[:6-d]: offs.add(d)
print("overlap offsets",sorted(offs))
for n in range(7,12):
  for c in itertools.product(range(4),repeat=n):
    o=occ(c)
    if len(o)<2: continue
    for i in o:
      ci=fl(c,i)
      if occ(ci)!=o: bad+=1
      for i2 in o:
        if i2!=i and fl(fl(c,i),i2)!=fl(fl(c,i2),i): bad+=1
  print(n,"bad",bad)
