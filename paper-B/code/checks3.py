"""General radius r. g takes a tuple xi=(xi_0..xi_{2r-1}), xi_k = y_{i-r+1+k};
f = y_{i-r} xor g.  Provides F on words and on periodic configurations, strong connectivity
of Lambda_M, the certificate (T_P), J(g), self-duality, and a search for (I)_q witnesses.
Run as a script: prints these data for the radius-1 rules and the named radius-2 rules."""
import itertools, math
from ca import scc_count

def Fword_g(g, r, w):
    out = []
    for i in range(r, len(w)-r):
        x, y = w[i] >> 1, w[i] & 1
        xi = tuple(w[i-r+1+k] & 1 for k in range(2*r))
        out.append(2*y + (x ^ (w[i-r] & 1) ^ g(xi)))
    return tuple(out)

def Fper_g(g, r, c):
    n = len(c)
    return tuple(2*(c[i] & 1) + ((c[i] >> 1) ^ (c[(i-r) % n] & 1) ^ g(tuple(c[(i-r+1+k) % n] & 1 for k in range(2*r)))) for i in range(n))

def lambda_sc(g, r, M):
    words = list(itertools.product(range(4), repeat=M+1))
    idx = {w: i for i, w in enumerate(words)}
    succ = []
    for w in words:
        s = set()
        for iota in itertools.product(range(4), repeat=2*r):
            s.add(idx[Fword_g(g, r, w + iota)])
        succ.append(sorted(s))
    return scc_count(len(words), lambda v: succ[v]) == 1

def TP(g, r, P):
    """certificate (T_P). column state = tuple of 2r P-bit tuples (s_1..s_2r)."""
    series = list(itertools.product((0, 1), repeat=P))
    states = list(itertools.product(series, repeat=2*r))
    def phi(s):
        return tuple(s[0][(t-1) % P] ^ g(tuple(s[2*r-1-k][t] for k in range(2*r))) for t in range(P))
    def succs(s):
        ph = phi(s)
        if sum(ph) % 2: return []
        res = []
        for z0 in (0, 1):
            z = [z0]
            for t in range(P-1): z.append(z[-1] ^ ph[t])
            res.append(s[1:] + (tuple(z),))
        return res
    cur = set(states)
    while True:
        nxt = set()
        for s in cur: nxt.update(succs(s))
        if nxt == cur: break
        cur = nxt
    return any(sum(phi(s)) % 2 == 1 for s in cur), len(cur)

def J_of(g, r):
    n = 2*r
    def d(fn, j):
        return lambda xi: fn(xi) ^ fn(tuple(v ^ (k == j) for k, v in enumerate(xi)))
    pts = list(itertools.product((0, 1), repeat=n))
    d0 = d(g, 0)
    const = len(set(d0(p) for p in pts)) == 1
    js = [j for j in range(1, n) if any(d(d0, j)(p) for p in pts)]
    return (max(js) if js else None), const

def selfdual(g, r):
    return all(g(tuple(1-v for v in p)) == 1 - g(p) for p in itertools.product((0, 1), repeat=2*r))

def witnesses(g, r, q, nmax):
    """search (I)_q witnesses in Per_n, n<=nmax. Return list of (n,b,s,k)."""
    aq = 0 if q % 2 else r % 2
    found = []
    for n in range(1, nmax+1):
        seen = set(); stabs = {}
        for c in itertools.product(range(4), repeat=n):
            if c in seen: continue
            # least period n?
            if any(c[d:] + c[:d] == c for d in range(1, n)):
                continue
            # Z^2 orbit: F-orbit and shifts
            forb = [c]; d_ = Fper_g(g, r, c)
            while d_ != c: forb.append(d_); d_ = Fper_g(g, r, d_)
            fset = {x: i for i, x in enumerate(forb)}
            # smallest s>0 with sigma^s c in F-orbit
            orbit = set()
            for sft in range(n):
                for x in forb: orbit.add(x[sft:] + x[:sft])
            seen |= orbit
            s0 = next(s for s in range(1, n+1) if (c[s % n:] + c[:s % n]) in fset)
            # sigma^{s0} c = F^{j} c  => (s0, -j) in Stab ; Stab = <(s0,-j),(0,L)>
            j = fset[c[s0 % n:] + c[:s0 % n]]; L = len(forb)
            key = (s0, (-j) % L, L)
            stabs[key] = stabs.get(key, 0) + 1
        for (s0, b0, L), k in stabs.items():
            # lattice gens (s0,b0),(0,L). (aq,1) in qZ^2 + lattice ?
            inl = False
            for m1 in range(q):
                for m2 in range(q):
                    v = ((m1*s0) % q, (m1*b0 + m2*L) % q)
                    if v == (aq % q, 1 % q): inl = True
            if k % q and not inl:
                found.append((n, (s0, b0), (0, L), k))
    return found

if __name__ == "__main__":
    print("== r=1: (T_P) ==")
    from ca import fbit
    for rule in (30, 45, 75, 120, 180, 210, 135, 225, 90, 150, 165, 240, 15, 60, 105, 195):
        g = (lambda rule: lambda xi: fbit(rule, 0, xi[0], xi[1]))(rule)
        print(rule, [TP(g, 1, P)[0] for P in range(1, 6)], "J,const d0:", J_of(g, 1), "selfdual", selfdual(g, 1), "I_2 witnesses n<=3:", witnesses(g, 1, 2, 3)[:3])
    print("== r=2 ==")
    maj = lambda a, b, c: (a & b) | (a & c) | (b & c)
    rules2 = {
        "xi0xi2": lambda x: x[0] & x[2],
        "xi0xi3": lambda x: x[0] & x[3],
        "xi3+xi0xi1xi2": lambda x: x[3] ^ (x[0] & x[1] & x[2]),
        "xi0xi1": lambda x: x[0] & x[1],
        "xi1xi2": lambda x: x[1] & x[2],
        "maj(0,2,3)": lambda x: maj(x[0], x[2], x[3]),
        "xi0+xi3": lambda x: x[0] ^ x[3],
        "xi1": lambda x: x[1],
        "xi3": lambda x: x[3],
    }
    for name, g in rules2.items():
        tp = [TP(g, 2, P)[0] for P in range(1, 4)]
        lam = [lambda_sc(g, 2, M) for M in range(0, 4)]
        print(name, "J", J_of(g, 2), "sd", selfdual(g, 2), "T_P P=1..3", tp, "Lambda_M SC M=0..3", lam, "I_2 wit n<=4:", witnesses(g, 2, 2, 4)[:3])
