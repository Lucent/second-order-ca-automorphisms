# Confirmation (not used in the proof): exhaustive near-glue check of the move lemma for rule 93.
# For r in {0,1}, every tail phase w (written tau here), every edge cell != vacuum cell, and every
# kink K^-_m = F^m K^- with m = w (mod 6) whose vacuum part covers (r, r+2] and whose first
# non-vacuum cell p satisfies 3 <= p-r <= 14, the glued configuration x |_r K^-_m is simulated
# with column r-1 an arbitrary input sequence (all sequences, branching only where f reads it).
# Asserted, for every branch:
#   * the outgoing kink K^+_n appears and F^t y agrees with K^+_{t+U} on (r,oo);
#   * the debris equals S_u^{w'} x with u = 2 - m - 3r, w' = w+2 (z_{u-1}=0) or w-1 (z_{u-1}=1);
#   * U = m + 6r + 3 z_{u-1} - 4;  the new edge cell is not a vacuum cell.
import sys
sys.argv=[sys.argv[0],'0','12','14']
from base import *
from kinks import vcell
import glue
from glue import run, kp
def qplus(n): return next(i for i in range(-150,150) if kp(n,i)!=vcell(n%6,i))
tot=0; bad=0; branches=0
for r in (0,1):
    for tau in range(6):
        for edge in range(4):
            if edge==vcell(tau,r): continue
            for shift in range(0,24,2):
                out=run(r,tau,edge,shift)
                if out is None: continue
                p,uc,res,unres=out
                if p-r>14: continue
                m=tau-3*shift          # glue.py: the kink km(tau,shift) is sigma^{-shift} F^tau K^- = F^{tau-3 shift} K^-
                u=2-m-3*r
                assert unres==0
                tot+=1
                for alive,n0,t in res:
                    branches+=1
                    hit=[a for a in alive if a[0][0]==u]
                    if len(hit)!=1: bad+=1; print('no unique alive candidate at u',r,tau,edge,m,alive); continue
                    (uu,z,zu),tE=hit[0]
                    ok = tE==(tau+2-3*z)%6 and n0==m+6*r+3*z-4 and 2*z+zu!=vcell(tE+u,r) \
                         and qplus(n0+t)<=r+3*(p-r)+30
                    if not ok: bad+=1; print('BAD',r,tau,edge,m,hit,n0)
print('near-glue cases',tot,'branches',branches,'failures',bad)
assert bad==0 and tot==216
