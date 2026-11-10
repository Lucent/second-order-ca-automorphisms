"""Rules 127R, 179R, 223R, Proposition 6.19, exhaustively: the classes of 127, 179, 223, 59, 27 under reflection
and complementation; each word of the listed cycle is mapped to the next by F and to the previous by F^-1
in all 16 contexts; the words are distinct; words half a period apart differ exactly at position 2; the
offsets at which two words of the cycle overlap.
"""
import itertools
from hw_common import Floc, Finvloc, wimg, frule
cyc={127:"03112 12230 30021 21103 03312 12330 30221 21003",
     179:"01110 12321 20302 01310 12221 20002",
     223:"00122 11300 33211 22133 00322 11200 33011 22033"}
def conj(r):
    f=frule(r); out=set()
    for refl in (0,1):
      for comp in (0,1):
        g=0
        for a,b,c in itertools.product((0,1),repeat=3):
          A,B,C=(c,b,a) if refl else (a,b,c)
          if comp: A,B,C=1-A,1-B,1-C
          g|=f(A,B,C)<<(4*a+2*b+c)
        out.add(g)
    return sorted(out)
for r in (127,179,223,59,27): print(r,conj(r))
for r,s in cyc.items():
    W=[tuple(int(ch) for ch in w) for w in s.split()]; P=len(W); h=P//2; k=5; j=2
    F=Floc(r); Fi=Finvloc(r)
    assert len(set(W))==P
    for s_,w in enumerate(W):
        for l in range(4):
            for rr in range(4):
                assert wimg(F,(l,)+w+(rr,))==W[(s_+1)%P],(r,w,l,rr)
                assert wimg(Fi,(l,)+w+(rr,))==W[(s_-1)%P]
        diff=[q for q in range(k) if w[q]!=W[(s_+h)%P][q]]
        assert diff==[j],(r,w,diff)
    # overlaps
    offs=set()
    for a in W:
        for b in W:
            for dlt in range(1,k):
                if a[dlt:]==b[:k-dlt]: offs.add(dlt)
    print(r,"P",P,"offsets",sorted(offs))
    # check (H2),(H4) directly and flip-all/trains by brute force on words len<=10 not needed given (c)
