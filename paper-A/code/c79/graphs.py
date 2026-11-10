# The move graphs Gamma_L on the level-L states, from the three-valued engine (engine.py).
# build(L,umin,umax) returns the valid states, the move edges at the times umin..umax, and the edge values.
import numpy as np
from engine import *
def build(L,umin,umax,ztr=None):
    """states: global id = tau*K + code(cells -L..0). returns valid, list of (u,src_ids,dst_ids), zkeys"""
    allc=all_states(L); K=len(allc)
    pw=4**np.arange(L,-1,-1)
    valid=np.zeros(6*K,bool); edges=[]; Z=np.zeros((6*K,0),np.uint8)
    zt=list(range(-L-1,L+1)) if ztr is None else ztr
    Z=np.zeros((6*K,len(zt)),np.uint8)
    for tau in range(6):
        ok=np.nonzero(allc[:,L]!=vcell(tau,0))[0]
        valid[tau*K+ok]=True
        cells=allc[ok]
        D=Diag(tau,cells,L,max(abs(umin),abs(umax),L)+1)
        for j,t in enumerate(zt):
            zz=D.z(t); assert (zz!=U).all(),('z unknown',t); Z[tau*K+ok,j]=zz
        for u in range(umin,umax+1):
            m=move(D,u)
            if m is None: continue
            kind,nt,nc=m
            assert (nt!=99).all() and (nc!=10).all(), ('inexact move',L,tau,u)
            assert (nc[:,0]==cells[:,0]).all(), 'exterior changed'
            for t2 in range(6):
                s=nt==t2
                assert (nc[s,L]!=vcell(t2,0)).all(),'edge became vacuum'
            dst=nt.astype(np.int64)*K+(nc.astype(np.int64)*pw).sum(1)
            edges.append((u,kind,tau*K+ok,dst))
    return valid,edges,Z,allc,K
