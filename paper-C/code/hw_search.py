"""Search for invariant families with a hidden involution, conditions (W1) and (H1) of Section 6.1, among the
words of length k (usage: hw_search.py K [RULES]; one rule per class under reflection and complementation).
For each word the closure under F and F^-1 in all 16 contexts is computed; for each closure of at most `cap`
words and each interior position j, the program looks for an involution iota of the family, without fixed
points, that changes only position j and commutes with F and F^-1 in every context. Writes search_kK.json.
"""
import sys, numpy as np, itertools, json

from centralizer import make_rule
def refl(r): return sum(((r>>((a<<2)|(b<<1)|c))&1)<<((c<<2)|(b<<1)|a) for a in (0,1) for b in (0,1) for c in (0,1))
def comp(r): return sum(((r>>(((1-a)<<2)|((1-b)<<1)|(1-c)))&1)<<((a<<2)|(b<<1)|c) for a in (0,1) for b in (0,1) for c in (0,1))
def cls(r): return {r,refl(r),comp(r),refl(comp(r))}
reps=sorted({min(cls(r)) for r in range(256)})
def analyse(rule,k,cap=200):
    f,Fl,Fil=make_rule(rule)
    Ft=np.array([[[Fl(a,b,c) for c in range(4)] for b in range(4)] for a in range(4)])
    Fit=np.array([[[Fil(a,b,c) for c in range(4)] for b in range(4)] for a in range(4)])
    n=4**k; idx=np.arange(n)
    W=np.stack([(idx>>(2*(k-1-i)))&3 for i in range(k)],1)
    succ=[]
    for T in (Ft,Fit):
        for a in range(4):
            for b in range(4):
                U=np.concatenate([np.full((n,1),a),W,np.full((n,1),b)],1)
                V=T[U[:,:-2],U[:,1:-1],U[:,2:]]
                succ.append(sum(V[:,i]<<(2*(k-1-i)) for i in range(k)))
    succ=np.stack(succ,1)   # n x 32 ; first 16 forward, last 16 backward
    # transitive closure by boolean matmul
    M=np.zeros((n,n),np.float32); M[np.repeat(idx,32),succ.ravel()]=1; M[idx,idx]=1
    for _ in range(int(np.ceil(np.log2(n)))+1):
        M=((M@M)>0).astype(np.float32)
    R=M>0
    sizes=R.sum(1)
    found=[]
    seen=set()
    for w0 in np.nonzero(sizes<=cap)[0]:
        fam=frozenset(np.nonzero(R[w0])[0].tolist())
        if fam in seen: continue
        seen.add(fam)
        famset=fam
        for j in range(1,k-1):
            sh=2*(k-1-j)
            partners=[w0 ^ (d<<sh) for d in (1,2,3) if (w0 ^ (d<<sh)) in famset]
            for p0 in partners:
                iota={w0:p0}; stack=[w0]; ok=True
                while stack and ok:
                    w=stack.pop(); v=iota[w]
                    for e in range(32):
                        a,b=succ[w,e],succ[v,e]
                        if a in iota:
                            if iota[a]!=b: ok=False;break
                        else:
                            if (a ^ b) & ~(3<<sh): ok=False;break
                            if b not in famset: ok=False;break
                            iota[a]=b; stack.append(a)
                if ok and len(iota)==len(famset) and all(iota[iota[w]]==w for w in iota) and all(iota[w]!=w for w in iota):
                    found.append((len(famset),j,sorted(famset)[:3]))
                    break
    return found, sizes
def fixed_word_exists(rule,k):
    pass
if __name__=='__main__':
    k=int(sys.argv[1]); out={}
    only=None if len(sys.argv)<3 else set(map(int,sys.argv[2].split(',')))
    for r in reps:
        if only and r not in only: continue
        fnd,sizes=analyse(r,k)
        if fnd:
            out[r]=[(int(a),int(b)) for a,b,_ in fnd]
            print(r, sorted(cls(r)), 'fixedword' if sizes.min()==1 else 'NOfixed', 'families with hidden involution (size,j):', sorted(set(out[r]))[:12], 'n=',len(fnd), 'min closure size', int(sizes.min()), flush=True)
    json.dump(out,open(f'search_k{k}.json','w'))
