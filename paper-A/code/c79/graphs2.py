# The move graphs Gamma_L built with the second engine (engine2.py): states, edge values Z and moves M.
from functools import lru_cache
from engine2 import *
@lru_cache(None)
def build(L):
    S=list(states(L)); Z={}; M={}
    for s in S:
        D=diag(s,L+3); Z[s]={t:D.z(t) for t in range(-L-1,L+1)}
        assert U not in Z[s].values()
        M[s]={m:move(s,m) for m in range(-L,L+1)}
        M[s]={m:v for m,v in M[s].items() if v is not None}
    return S,Z,M
