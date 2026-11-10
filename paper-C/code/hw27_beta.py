"""Rule 27R, Theorem 6.20: the images of the 12 words of W_27 over all 16 contexts, and all maps
beta: W_27 -> Z/2 with beta o iota = beta + 1 that satisfy (H3) for some epsilon. Expected: two solutions,
beta and beta + 1, both with epsilon = 1.
"""
import itertools
from hw_common import Floc, wimg
F=Floc(27)
W=[(a,p,p,C,(p+1)%3,(p+1)%3) for p in range(3) for a in (p,3) for C in ((p+1)%3,3)]
S=set(W)
def io(w): w=list(w); w[3]=3 if w[3]!=3 else (w[1]+1)%3; return tuple(w)
trans={}
for w in W:
    trans[w]={wimg(F,(l,)+w+(r,)) for l in range(4) for r in range(4)}
    print(w,'->',sorted(trans[w]))
reps=[]
for w in W:
    if io(w) not in reps: reps.append(w)
sols=[]
for bits in itertools.product([0,1],repeat=6):
    beta={}
    for b,w in zip(bits,reps): beta[w]=b; beta[io(w)]=1-b
    for eps in (0,1):
        if all(beta[v]==(beta[w]+eps)%2 for w in W for v in trans[w]): sols.append((bits,eps))
print("beta solutions |K|=1:",sols)
