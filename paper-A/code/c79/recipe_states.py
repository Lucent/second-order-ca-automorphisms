# Configurations for the recipes recipe_H2.py and recipe_H3.py.
# A(state,L): every configuration with a given level-L state, in three-valued logic (engine2.py); a value is used
#   only if it is known, that is the same for every completion of the cells left of r-L; otherwise Unknown is raised.
# C(conf): one whole configuration, simulated exactly (exact2.py).
# Both check the one-move lemma (phase and blocks of T_m c) at every move used.
import sys,itertools
from engine2 import diag,move,states,U,vcell
from exact2 import Conf
class Unknown(Exception): pass
class A:
    def __init__(s,state,L): s.st=state; s.th=state[0]; s.L=L; s.D=diag(state,L+8); s.cache={}
    def col(s,t,j):
        v=s.D.g(t,-j)
        if v==U: raise Unknown((s.st,t,j))
        return v
    z=lambda s,t: s.col(t,0)
    y=lambda s,t: s.col(t,1)
    w=lambda s,t: s.col(t,2)
    def n(s,t): return (t+s.th)%6
    def blk(s,t):
        n=s.n(t); pos=n if n else 6; t1=t-pos+1
        p,q=s.z(t1),s.z(t1+1)
        assert [s.z(t1+i) for i in range(6)]==[p,q,1-p,1-q,p,q] and (p,q)!=(0,0)
        return (p,q),pos
    def move(s,m):
        if abs(m)>s.L: raise Unknown(('m',m))
        try: st=move(s.st,m)
        except AssertionError: raise Unknown(('mv',m))
        assert st is not None
        D=A(st,s.L)
        n=s.n(m)
        if n==2 and s.z(m-1)==1:   assert D.th==(s.th-1)%6 and D.blk(m)==((s.z(m),1^(s.z(m)&(1-s.y(m)))),1)
        elif n==2:                 assert D.th==(s.th+2)%6 and D.blk(m)==((1,1-s.z(m)),4)
        elif n==0 and s.z(m)==0:   assert D.th==(s.th-2)%6 and D.blk(m)==((1-s.z(m-1),1),4)
        elif n==0:                 assert D.th==(s.th+1)%6 and D.blk(m)==((1,s.z(m+1)),1) and D.blk(m-1)==((1^(s.z(m-1)&(1-s.y(m-1))),s.z(m-1)),6)
        else: raise Exception
        return D,None
class C:
    def __init__(s,conf): s.c=conf; s.th=conf.th; s.cache={}
    def col(s,t,j):
        k=(t,j)
        if k not in s.cache: s.cache[k]=s.c.col(t,-j)
        return s.cache[k]
    z=lambda s,t: s.col(t,0)
    y=lambda s,t: s.col(t,1)
    w=lambda s,t: s.col(t,2)
    def n(s,t): return (t+s.th)%6
    def blk(s,t):
        n=s.n(t); pos=n if n else 6; t1=t-pos+1
        p,q=s.z(t1),s.z(t1+1)
        assert [s.z(t1+i) for i in range(6)]==[p,q,1-p,1-q,p,q] and (p,q)!=(0,0)
        return (p,q),pos
    def move(s,m):
        d=s.c.move(m); assert d is not None; D=C(d)
        # one-move lemma
        n=s.n(m)
        if n==2 and s.z(m-1)==1:   kind='F-'; assert D.th==(s.th-1)%6; assert D.blk(m)==((s.z(m),1^(s.z(m)&(1-s.y(m)))),1)
        elif n==2:                 kind='F+'; assert D.th==(s.th+2)%6; assert D.blk(m)==((1,1-s.z(m)),4)
        elif n==0 and s.z(m)==0:   kind='B-'; assert D.th==(s.th-2)%6; assert D.blk(m)==((1-s.z(m-1),1),4)
        elif n==0:                 kind='B+'; assert D.th==(s.th+1)%6; assert D.blk(m)==((1,s.z(m+1)),1); assert D.blk(m-1)==((1^(s.z(m-1)&(1-s.y(m-1))),s.z(m-1)),6)
        else: raise Exception('no move')
        for t in (m-1,m): assert D.z(t)==s.z(t)
        for t in range(m-2,m+2): assert D.y(t)==s.y(t)
        for t in range(m-3,m+3): assert D.w(t)==s.w(t)
        return D,kind
