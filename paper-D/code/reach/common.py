import itertools
def fv(rule,a,b,c): return (rule>>(4*a+2*b+c))&1
def mirror(r): return sum(fv(r,c,b,a)<<(4*a+2*b+c) for a in (0,1) for b in (0,1) for c in (0,1))
def comp(r): return sum(fv(r,1-a,1-b,1-c)<<(4*a+2*b+c) for a in (0,1) for b in (0,1) for c in (0,1))
def cls(r): return sorted({r,mirror(r),comp(r),mirror(comp(r))})
def Fstep(rule,x,y):  # periodic rows (tuples)
    p=len(y); return y, tuple(x[i]^fv(rule,y[i-1],y[i],y[(i+1)%p]) for i in range(p))
def canon_rot(state):
    x,y=state; p=len(x)
    return min((x[i:]+x[:i],y[i:]+y[:i]) for i in range(p))
def backgrounds(rule,pmax=4):
    """spatially periodic F-cycles (minimal spatial period p), left-blind ones"""
    seen=set(); out=[]
    for p in range(1,pmax+1):
        for xs in itertools.product((0,1),repeat=p):
            for ys in itertools.product((0,1),repeat=p):
                s=(xs,ys)
                # minimal period p
                if any(xs==xs[d:]+xs[:d] and ys==ys[d:]+ys[:d] for d in range(1,p) if p%d==0): continue
                if canon_rot(s) in seen: continue
                cyc=[s]; t=Fstep(rule,*s)
                while t!=s: cyc.append(t); t=Fstep(rule,*t)
                for c in cyc: seen.add(canon_rot(c))
                blind=all(fv(rule,0,y[i],y[(i+1)%p])==fv(rule,1,y[i],y[(i+1)%p]) for (x,y) in cyc for i in range(p))
                out.append((p,len(cyc),cyc,blind))
    return out
