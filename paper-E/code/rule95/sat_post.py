"""Post-processing of the output sat_R<R>_<mode>.pkl of sat_sublattice.py (read from the current directory): are the restrictions to Phi^+ / Phi^- of all solutions particle-conserving block
   maps of G (or the zero map of the constant solution)?  For mode 'aut': describe the solutions that are not of the
   form sigma^a F^b and write their tables to sat_others_R<R>_<type>.npz.   usage: sat_post.py R all|aut"""
import pickle, sys, numpy as np
from enum_golden import gwords, injective
def conserving(w, f):
    mask = (1 << (w - 1)) - 1; lam = {0: 0}; stack = [0]; E = gwords(w); adj = {}
    for e in E:
        s, t = e >> 1, e & mask; wt = f[e] - ((e >> (w - 1)) & 1)
        adj.setdefault(s, []).append((t, -wt)); adj.setdefault(t, []).append((s, wt))
    while stack:
        v = stack.pop()
        for u, d in adj[v]:
            if u not in lam: lam[u] = lam[v] + d; stack.append(u)
            elif lam[u] != lam[v] + d: return False
    for e in E:
        for b in (0, 1):
            e2 = ((e & mask) << 1) | b
            if e2 in f and f[e] and f[e2]: return False
    return True
def minwin(w, t):
    """minimal window [lo,hi] (positions 0..w-1, first letter = position 0) on which the table depends"""
    W = gwords(w); f = dict(zip(W, t)); dep = []
    for p in range(w):
        bit = 1 << (w - 1 - p)
        if any((e ^ bit) in f and f[e ^ bit] != f[e] for e in W): dep.append(p)
    return (min(dep), max(dep)) if dep else None
R = int(sys.argv[1]); mode = sys.argv[2]
d = pickle.load(open("sat_R%d_%s.pkl" % (R, mode), "rb")); W = gwords(R + 1)
for par in (0, 1):
    r = d[par]
    if mode == "all":
        S = set()
        for lab, key, iP, iM in r["rows"]:
            for t in iP + iM: S.add(t)
        cons = sum(conserving(R + 1, dict(zip(W, t))) for t in S); zero = sum(1 for t in S if not any(t))
        print("R=%d type %s: %d solutions; distinct class maps occurring in restrictions: %d, conserving %d, zero map %d; minimal windows: %s" % (
            R, "PS"[par], len(r["rows"]), len(S), cons, zero, sorted(set(str(minwin(R + 1, t)) for t in S))))
    else:
        for sgn in ("+1", "-1"):
            S = set(t for info in r["restr" + sgn] for t in info)
            cons = sum(conserving(R + 1, dict(zip(W, t))) for t in S); zero = sum(1 for t in S if not any(t))
            print("R=%d type %s Phi^%s: %d restrictions; distinct class maps %d, conserving %d, zero %d" % (R, "PS"[par], sgn[0], r["nrestr" + sgn], len(S), cons, zero))
        others = [x for x in r["found"] if x[0] is None]
        print("   found %d, named %d, others %d" % (len(r["found"]), len(r["found"]) - len(others), len(others)))
        for x in others:
            print("   other: Phi+ class maps (x-bit,y-bit) windows", [minwin(R + 1, t) for t in x[2]], "| Phi- windows", [minwin(R + 1, t) for t in x[3]],
                  "| shift-type Phi+:", [sum(t) == sum(1 for e in W if (e >> k) & 1) and all(t[i] == ((W[i] >> k) & 1) for i in range(len(W))) for t in x[2] for k in range(R + 1) if all(t[i] == ((W[i] >> k) & 1) for i in range(len(W)))],
                  "Phi-:", [k for t in x[3] for k in range(R + 1) if all(t[i] == ((W[i] >> k) & 1) for i in range(len(W)))])
        # one row per solution: the 2^(2R+1) entries of TX followed by those of TY
        np.savez_compressed("sat_others_R%d_%s.npz" % (R, "PS"[par]),
                            tables=np.array([np.frombuffer(x[1], dtype=np.uint8) for x in others], dtype=np.uint8).reshape(len(others), 4 << (2 * R)))
