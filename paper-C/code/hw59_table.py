"""Rule 59R: the clock formula for the local rule on all 64 triples; the table of Lemma 6.11 (for each word of
X_A and X_B, in all 16 contexts, for F and F^-1: the image lies in the same family, its cells L, C, R do not
depend on the context, the cores follow the two 6-cycles; the possible protectors are printed); and Lemma
6.12 (the flip iota is an involution without fixed points, commutes with F and F^-1 in every context, and
beta o iota = beta + 1, beta o F = beta + 1). Defines XA, XB, iota, beta for the other hw59_*.py programs.
"""
from hw_common import *
import itertools
F=Floc(59); Fi=Finvloc(59)
def XA(): return [(a,p,C,(p+1)%3) for p in range(3) for a in (p,3) for C in ((p+1)%3,3)]
def XB(): return [(a,p,C,(p+2)%3) for p in range(3) for a in (p,3) for C in ((p+2)%3,3)]
fam={'A':(XA(),1),'B':(XB(),2)}
assert not set(XA())&set(XB())
# the clock formula (equation (5) of the paper)
def trig(l,r): return (l&1)==0 and (r&1)==1
for l in range(4):
  for s in range(4):
    for r in range(4):
      exp={0:1,2:0,1:3 if trig(l,r) else 2,3:2 if trig(l,r) else 3}[s]
      assert F(l,s,r)==exp
print("clock ok")
cyc={'A':[(0,1),(1,3),(2,3),(0,3),(1,2),(2,0)],'B':[(0,2),(1,0),(2,1),(0,3),(1,3),(2,3)]}
for X,(W,dl) in fam.items():
  Ws=set(W)
  for G,name in ((F,'F'),(Fi,'Finv')):
    for w in W:
      imgs=set(); prot=set()
      for l in range(4):
        for r in range(4):
          im=wimg(G,(l,)+w+(r,))
          assert im in Ws,(X,name,w,l,r,im)
          if name=='F': assert im[1]==(w[1]+1)%3
          else: assert im[1]==(w[1]-1)%3
          imgs.add(im[1:]); prot.add(im[0])
      assert len(imgs)==1
      if name=='F':
        core=(w[1],w[2]); i=cyc[X].index(core); nc=(list(imgs)[0][0],list(imgs)[0][1])
        assert nc==cyc[X][(i+1)%6],(X,w,nc)
        print(X,w,'->',sorted(prot),list(imgs)[0])
print("walls ok")
# H1 flip, beta
def iota(w,dl):
  a,L,C,R=map(int,w); C2=3 if C!=3 else (L+dl)%3; return (a,L,C2,R)
def beta(X,w):
  a,L,C,R=map(int,w); return (int(C==3)+int(L==(2 if X=="A" else 1)))%2
for X,(W,dl) in fam.items():
  for w in W:
    assert iota(iota(w,dl),dl)==w and iota(w,dl)!=w and beta(X,iota(w,dl))==(beta(X,w)+1)%2
    for G in (F,Fi):
     for l in range(4):
      for r in range(4):
        im=wimg(G,(l,)+w+(r,)); im2=wimg(G,(l,)+iota(w,dl)+(r,))
        assert im2==iota(im,dl),(X,w,l,r)
        if G is F: assert beta(X,im)==(beta(X,w)+1)%2
print("flip ok")
