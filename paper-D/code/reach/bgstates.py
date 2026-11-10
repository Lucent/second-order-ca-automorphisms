from sim import Fst, bg_orbit
def shift(s,d):
    xs,ys=s;p=len(xs); return (tuple(xs[(i+d)%p] for i in range(p)),tuple(ys[(i+d)%p] for i in range(p)))
def all_states(rule,s):
    """all B states: orbit x spatial shifts, as list of (xs,ys) aligned to global index"""
    out=[]
    for t in bg_orbit(rule,s):
        for d in range(len(s[0])):
            u=shift(t,d)
            if u not in out: out.append(u)
    return out
def rho(s): return (s[1],s[0])
