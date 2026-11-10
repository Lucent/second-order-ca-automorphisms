"""Exhaustive test of (T_P) at radius 2:
for every g:{0,1}^4->{0,1} (65536 functions) test (T_1), then (T_2), then (T_3);
for those failing all, build Lambda_M (M=2,3) directly. Report by affine / non-affine."""
import itertools, sys
from checks3 import lambda_sc

def TPtab(tab, P):
    # tab: 16-bit int; g(xi0..xi3) = bit (xi0 + 2 xi1 + 4 xi2 + 8 xi3)
    # state: (s1,s2,s3,s4) each a P-bit int (bit t = value at time t). xi_k = s_{4-k}
    mask = (1 << P) - 1
    nst = 1 << (4*P)
    phis = [0]*nst; par = [0]*nst
    for s in range(nst):
        s1 = s & mask; s2 = (s >> P) & mask; s3 = (s >> 2*P) & mask; s4 = (s >> 3*P) & mask
        ph = 0
        for t in range(P):
            xi = ((s4 >> t) & 1) | (((s3 >> t) & 1) << 1) | (((s2 >> t) & 1) << 2) | (((s1 >> t) & 1) << 3)
            b = ((s1 >> ((t-1) % P)) & 1) ^ ((tab >> xi) & 1)
            ph |= b << t
        phis[s] = ph; par[s] = bin(ph).count("1") & 1
    def succs(s):
        if par[s]: return ()
        ph = phis[s]; res = []
        for z0 in (0, 1):
            z = z0; cur = z0
            for t in range(P-1):
                cur ^= (ph >> t) & 1
                z |= cur << (t+1)
            res.append((s >> P) | (z << 3*P))
        return res
    cur = set(range(nst))
    while True:
        nxt = set()
        for s in cur: nxt.update(succs(s))
        if nxt == cur: break
        cur = nxt
    return any(par[s] for s in cur)

def is_affine(tab):
    f = lambda x: (tab >> x) & 1
    c = f(0); co = [f(1 << k) ^ c for k in range(4)]
    return all(f(x) == c ^ (sum(co[k] for k in range(4) if (x >> k) & 1) & 1) for x in range(16))

if __name__ == "__main__":
    stats = {}
    fails = []
    for tab in range(65536):
        aff = is_affine(tab)
        lvl = None
        for P in (1, 2, 3):
            if TPtab(tab, P): lvl = P; break
        stats[(aff, lvl)] = stats.get((aff, lvl), 0) + 1
        if lvl is None: fails.append((tab, aff))
    print("counts by (affine, least P<=3 with T_P):", stats)
    print("functions with no T_P, P<=3:", len(fails))
    for tab, aff in fails:
        g = lambda x, tab=tab: (tab >> (x[0] | (x[1] << 1) | (x[2] << 2) | (x[3] << 3))) & 1
        sc = [lambda_sc(g, 2, M) for M in (2, 3)]
        print("  tab=%d (0x%04x) affine=%s Lambda_M SC for M=2,3: %s" % (tab, tab, aff, sc))
