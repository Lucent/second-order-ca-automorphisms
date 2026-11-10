"""Rule 59R, Proposition 6.16: each of the eight generators (the flips of the walls of one class) commutes with
F on all 4^9 windows (radius 3); the number of entries of its radius-3 table that differ from the identity
is printed. The flip of all walls of a family is the same map whether composed from the generators of T_R
or from those of T_L.
"""
import sys

from classify_r4 import check_commutes_np_R
from hw59_census_lib import *
from hw_common import Floc
F=Floc(59)
def tab3(cls):
    n=7; W=allwords(n); w=flipgen(Wd(W,-3),cls)  # out start 0, length 2
    return w.a[:,0-w.start].astype(np.int64)
for g in GR+GL:
    t=tab3([g]); print(g,check_commutes_np_R(t,F,3), (t!=allwords(7)[:,3]).sum())
# also flip-all of A and B
for X in 'AB':
    print('tau',X,check_commutes_np_R(tab3([(X,'='),(X,'!=')]),F,3), (tab3([(X,'='),(X,'!=')])==tab3([(X,'pp'),(X,'npp')])).all())
