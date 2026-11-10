# exhaustive: rule 95 links with core length <=K (zeros on bg 1, no two zeros at distance 2).
# check (1) full Goursat (gs) and (2) direct (t,i) simulation: Gamma and O are sublattice-wise translates.
from sim import *
import itertools, sys
K=int(sys.argv[1]); usedirect=int(sys.argv[2])
def cores(K):
    out=[]
    for n in range(1,K+1):
        for c in itertools.product((0,1),repeat=n):
            if c[0]==1 or c[-1]==1: continue
            z=[i for i,x in enumerate(c) if x==0]
            if any(b-a==2 for a in z for b in z): continue
            out.append(z)
    return out
def rep(z):
    zs=set(z); ks=[k for k in range(1,max(z)+3) if ((k in zs)!=((k-1) in zs))]
    return (1,ks)
def zeros(r,n): return [k for k in range(n) if val(r,k)==0]
def transl(a,b):
    for s in (0,1):
        x=[z for z in a if z%2==s]; y=[z for z in b if z%2==s]
        if len(x)!=len(y): return False
        if x and len(set(v-u for u,v in zip(x,y)))!=1: return False
    return True
C=cores(K); n=0; bad=0
for zp in C:
    for zx in C:
        for sh in (0,1):
            P=rep([z+4 for z in zp]); X=rep([z+4+sh for z in zx])
            M=2*(len(zp)+len(zx)+K)+40
            for impl in ([gs,direct] if usedirect else [gs]):
                G,O=impl(95,P,X,M=M)
                if not (transl([z+4+sh for z in zx],zeros(G,M-10)) and transl([z+4 for z in zp],zeros(O,M-10))): bad+=1
            n+=1
print('K',K,'pairs',n,'failures',bad)
