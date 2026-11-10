# All-length search for words fixed by F_f in every context.
# A fixed word W has x(W_j)=y(W_j) (so W in {0,3}^k), and with w=y(W):
#   f(w_{j-1},w_j,w_{j+1})=0 for interior j, f(u,w_0,w_1)=0 and f(w_{k-2},w_{k-1},v)=0 for u,v in {0,1}
#   (k=1: f(u,w_0,v)=0 for all u,v).  Fixed words of length k>=2 are walks of k-1 states in the
#   graph on the 4 pairs (w_j,w_{j+1}); so if one exists, one of length <=5 exists.
import itertools
def f(rule,a,b,c): return (rule>>((a<<2)|(b<<1)|c))&1
def fixed_bin(rule,w):
    k=len(w)
    if k==1: return all(f(rule,u,w[0],v)==0 for u in (0,1) for v in (0,1))
    if any(f(rule,u,w[0],w[1]) for u in (0,1)): return False
    if any(f(rule,w[-2],w[-1],v) for v in (0,1)): return False
    return all(f(rule,w[j-1],w[j],w[j+1])==0 for j in range(1,k-1))
def lengths_dp(rule,K=40):
    # DP over states (pair) : reachable set of end-pairs for each length
    L=[]
    if any(fixed_bin(rule,(b,)) for b in (0,1)): L.append(1)
    start={(a,b) for a in (0,1) for b in (0,1) if all(f(rule,u,a,b)==0 for u in (0,1))}
    end={(a,b) for a in (0,1) for b in (0,1) if all(f(rule,a,b,v)==0 for v in (0,1))}
    cur=start
    for k in range(2,K+1):
        if cur & end: L.append(k)
        cur={(b,c) for (a,b) in cur for c in (0,1) if f(rule,a,b,c)==0}
    return L
# brute-force cross-check over the full alphabet for k<=4 (uses the F_f formula directly)
def Floc(rule,s0,s1,s2):
    x,y=s1>>1,s1&1; return (y<<1)|(x^f(rule,s0&1,y,s2&1))
def fixed_full(rule,W):
    k=len(W)
    return all(tuple(Floc(rule,*ctx[i:i+3]) for i in range(k))==W
               for ctx in ((a,)+W+(b,) for a in range(4) for b in range(4)))
mism=0
for rule in range(256):
    for k in (1,2,3,4):
        full={W for W in itertools.product(range(4),repeat=k) if fixed_full(rule,W)}
        binw={tuple(3*b for b in w) for w in itertools.product((0,1),repeat=k) if fixed_bin(rule,w)}
        if full!=binw: mism+=1
print("brute force over A^k (k<=4) agrees with the {0,3} reduction for all 256 rules:", mism==0)
# classes under reflection and complement kappa
def refl(r): return sum(f(r,c,b,a)<<((a<<2)|(b<<1)|c) for a in (0,1) for b in (0,1) for c in (0,1))
def kap(r):  return sum(f(r,1-a,1-b,1-c)<<((a<<2)|(b<<1)|c) for a in (0,1) for b in (0,1) for c in (0,1))
classes={}
for r in range(256):
    cl=frozenset({r,refl(r),kap(r),refl(kap(r))}); classes[min(cl)]=sorted(cl)
print("number of classes:",len(classes))
short=[rep for rep,cl in classes.items() if any(k<=3 for k in lengths_dp(rep))]
nshort_rules=sum(len(classes[r]) for r in short)
print("classes with a fixed word of length <=3:",len(short),"rules:",nshort_rules)
free=sorted(rep for rep in classes if rep not in short)
print("classes with none of length <=3:",len(free))
print(free)
for rep in free:
    for r in classes[rep]:
        L=lengths_dp(r,40)
        if L: print(f"  rule {r} (class {rep}): fixed-word lengths <=40: {L[:12]}{'...' if len(L)>12 else ''}")
