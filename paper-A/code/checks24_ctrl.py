# Direct check of the pulled-back controllability for 66: for tailed b with edge r and k<=3,
# some composite of 66 moves at times |m|<=k-1 gives Y(0,r-k+1) != Y(0,r-k+2).
import random
from eng24 import *

def tailflip(b, be, r, s, R=66):
    u = b.Fn(R, s)
    u2 = mk(u.lo, u.hi, lambda i: u.cell(i) if i <= r else (1 - be, 1 - be))
    return u2.Fn(R, -s)

def M66(b, be, r, m, bar=False):
    x, y = b.Fn(66, m).cell(r)
    want = (1 - be, be) if bar else (be, 1 - be)
    return tailflip(b, be, r, m) if (x, y) == want else b

rng=random.Random(11)
LO,HI=-60,80
def key(c,r): return tuple(c.cell(i) for i in range(r-6,r+3))
for trial in range(60):
    be=rng.randint(0,1); r=0
    b=rand_tailed(be,r,LO,HI,rng)
    for k in (1,2,3):
        seen={key(b,r)}; frontier=[(b,be)]; hit=False
        while frontier and not hit:
            nf=[]
            for c,bc in frontier:
                if c.cell(r-k+1)[1]!=c.cell(r-k+2)[1]: hit=True;break
                for m in range(-(k-1),k):
                    for bar in (False,True):
                        c2=M66(c,bc,r,m,bar)
                        bc2=c2.cell(r+1)[0]
                        # re-pad to the original window (moves change only cells >= r-|m|+1)
                        c2=mk(LO,HI,lambda i: c2.cell(i) if c2.lo<=i<c2.hi else (b.cell(i) if i<r else (bc2,bc2)))
                        kk=key(c2,r)
                        if kk not in seen: seen.add(kk); nf.append((c2,bc2))
            frontier=nf
        assert hit,(trial,k)
print('direct 66 controllability k<=3 ok')
print('ALL CHECKS PASSED')
