"""Rule 59R, Proposition 6.16: tools for comparing the census with the maps (i)-(iii). Reads the census
census_59_R12.json (written by `python3 centralizer.py --rules 59 --R 1,2 --out census_59_R12.json`) and
asserts that it has 133 invertible maps of radius 2. Words are NumPy arrays of all words of a given length;
Fmap, Sig, Sinv apply F, the shift and its inverse to them; flipgen applies the flip of all walls of given
classes; table_of composes a list of such operations exactly, on all windows of the given radius, and
returns the radius-2 table of the result (None if the result is not a function of the central 5 cells).
"""
import json, itertools, numpy as np
from hw_common import tables
T,Ti=tables(59)
d=json.load(open('census_59_R12.json'))
sols=[s for s in d[1]['solutions'] if s['invertible']]
soltabs={s['table']:s['label'] for s in sols}
assert len(soltabs)==133
# word arrays: (N, len) with absolute start position
def allwords(n):
    idx=np.arange(4**n,dtype=np.int64)
    return np.stack([(idx>>(2*(n-1-i)))&3 for i in range(n)],axis=1).astype(np.int8)
class Wd:
    def __init__(s,a,start): s.a=a; s.start=start
def Fmap(w,inv=False):
    TT=Ti if inv else T; a=w.a
    return Wd(TT[a[:,:-2],a[:,1:-1],a[:,2:]],w.start+1)
def Sig(w): # (sigma c)_n = c_{n+1}
    return Wd(w.a[:,1:],w.start)
def Sinv(w):
    return Wd(w.a[:,:-1],w.start+1)
FAM={'A':1,'B':2}
def wallat(a,i,dl):  # columns i..i+3 (array index); returns bool
    L=a[:,i+1]; C=a[:,i+2]; R=a[:,i+3]; p=a[:,i]
    ok=(L<3)&(R==(L+dl)%3)&((C==R)|(C==3))&((p==L)|(p==3))
    return ok
def flipgen(w,cls):
    # cls: set of (X,prop) with prop in {'=','!=','pp','npp'}; flip cell n if wall at n-2 with class in cls
    a=w.a; N=a.shape[1]
    lo,hi=(3,2)  # read [n-3,n+2]
    out=a[:,lo:N-hi].copy()
    for k in range(out.shape[1]):
        n=k+lo
        for X,prop in cls:
            dl=FAM[X]; i=n-2
            m=wallat(a,i,dl)
            if prop=='=': m&=(a[:,n+2]==a[:,n+1])
            elif prop=='!=': m&=(a[:,n+2]!=a[:,n+1])
            elif prop=='pp': m&=(a[:,i]==a[:,i+1])&((a[:,i-1]==a[:,i+1])|(a[:,i-1]==3))
            elif prop=='npp': m&=~((a[:,i]==a[:,i+1])&((a[:,i-1]==a[:,i+1])|(a[:,i-1]==3)))
            L=a[:,n-1]
            newv=np.where(a[:,n]==3,(L+dl)%3,3).astype(np.int8)
            out[:,k]=np.where(m,newv,out[:,k])
    return Wd(out,w.start+lo)
def table_of(ops,Rn):
    """ops: list of functions applied in order to word; return radius-2 table string or None if not radius 2"""
    n=2*Rn+1; W=allwords(n); w=Wd(W,-Rn)
    for op in ops: w=op(w)
    col=0-w.start
    assert 0<=col<w.a.shape[1]
    val=w.a[:,col]
    # radius-2 key
    c=Rn
    key=sum(W[:,c-2+i].astype(np.int64)<<(2*(4-i)) for i in range(5))
    tab=np.full(1024,-1,dtype=np.int64)
    tab[key]=val
    # check consistency
    if not (tab[key]==val).all(): return None
    return ''.join(map(str,tab))
GR=[('A','='),('A','!='),('B','='),('B','!=')]
GL=[('A','pp'),('A','npp'),('B','pp'),('B','npp')]
