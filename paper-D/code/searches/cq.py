"""Exact (rational) conserved quantities h: A^r -> Q of F, reduced modulo constants and coboundaries; prints a basis
of the quotient as functions on A^r (normal form: h with h(w)=0 on a complement chosen by elimination)."""
import itertools,sys
from fractions import Fraction as Fr
rule=int(sys.argv[1]); r=int(sys.argv[2])
f=lambda a,b,c:(rule>>((a<<2)|(b<<1)|c))&1
Floc=lambda s0,s1,s2: ((s1&1)<<1)|((s1>>1)^f(s0&1,s1&1,s2&1))
nh=4**r; ng=4**(r+1)
idx=lambda w: sum(s<<(2*(len(w)-1-i)) for i,s in enumerate(w))
rows=[]
for w in itertools.product(range(4),repeat=r+2):
    row={}
    def add(j,v): row[j]=row.get(j,0)+v
    Fw=tuple(Floc(*w[i:i+3]) for i in range(r)); add(idx(Fw),1); add(idx(w[1:r+1]),-1)
    add(nh+idx(w[1:]),-1); add(nh+idx(w[:-1]),1)
    rows.append({k:Fr(v) for k,v in row.items() if v})
# add trivial h's as extra "quotient" : we want solutions modulo trivial. Do: compute nullspace basis, then reduce mod trivial.
# Gaussian elimination (RREF) on dict rows
ncol=nh+ng; piv={}
for row in rows:
    row=dict(row)
    for c in sorted(row):
        pass
    # reduce
    while True:
        ks=[k for k in row if k in piv]
        if not ks: break
        c=ks[0]; pr=piv[c]; fc=row[c]
        for k,v in pr.items():
            nv=row.get(k,0)-fc*v
            if nv: row[k]=nv
            else: row.pop(k,None)
    if not row: continue
    c=min(row); fc=row[c]; row={k:v/fc for k,v in row.items()}
    for pc,pr in piv.items():
        if c in pr:
            fct=pr[c]
            for k,v in row.items():
                nv=pr.get(k,0)-fct*v
                if nv: pr[k]=nv
                else: pr.pop(k,None)
    piv[c]=row
free=[c for c in range(ncol) if c not in piv]
basis=[]
for fc in free:
    v=[Fr(0)]*ncol; v[fc]=Fr(1)
    for pc,pr in piv.items(): v[pc]=-pr.get(fc,0)
    basis.append(v[:nh])
# trivial: constants and coboundaries g'(w[1:])-g'(w[:-1]) with g' on A^{r-1}
triv=[[Fr(1)]*nh]
for e in range(4**(r-1)):
    v=[Fr(0)]*nh
    for w in itertools.product(range(4),repeat=r):
        v[idx(w)]=Fr(int(idx(w[1:])==e))-Fr(int(idx(w[:-1])==e))
    triv.append(v)
def rank(vs):
    M=[list(v) for v in vs]; rk=0; n=len(M[0]) if M else 0
    for c in range(n):
        p=next((i for i in range(rk,len(M)) if M[i][c]),None)
        if p is None: continue
        M[rk],M[p]=M[p],M[rk]
        for i in range(len(M)):
            if i!=rk and M[i][c]:
                fct=M[i][c]/M[rk][c]; M[i]=[a-fct*b for a,b in zip(M[i],M[rk])]
        rk+=1
    return rk
rt=rank(triv); cur=list(triv); out=[]
for v in basis:
    if rank(cur+[v])>rank(cur):
        cur.append(v); out.append(v)
print(f"rule {rule} r {r}: solution dim {len(basis)}, trivial rank {rt}, nontrivial {len(out)}")
for v in out:
    print("  h:", {''.join(map(str,w)):str(v[idx(w)]) for w in itertools.product(range(4),repeat=r) if v[idx(w)]})
