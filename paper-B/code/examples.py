"""Data of the examples: J(g), self-duality, least P with (T_P), witnesses of (I)_q with
their exact k(L), strong connectivity of Lambda_M, and explicit (T_P) certificates.

Rule: f(y_{i-r},...,y_{i+r}) = y_{i-r} XOR g(xi_0,...,xi_{2r-1}), xi_k = y_{i-r+1+k};
F(c)_i = (y_i, x_i XOR f(y_{i-r..i+r})); sigma(c)_i = c_{i+1}.
g is an integer truth table: g(xi) = bit number (xi_0 + 2 xi_1 + 4 xi_2 + ...) of tab.

Stabilisers.  For c of least spatial period n, Stab(c) = {(a,b): sigma^a F^b c = c} is
computed from the definition on the torus [0,n) x [0,l), l the F-period of c, and written
(n,b,s), meaning <(n,0),(s,b)>.  k(L) is the number of Z^2-orbits with stabiliser exactly L,
obtained by enumerating all points of least period n.  L is a witness of (I)_q iff q does
not divide k(L) and (a_q,1) is not in qZ^2 + L, with a_q = 0 for odd q and a_2 = r mod 2.

Lambda_M is built from the layer equations
    x'[e] = y[e],  y'[d] = v_d (d in B_y),
    y'[d] = y[d] XOR x[d-2r] XOR g(y[d-1],...,y[d-2r])  (d >= -M+r),  y[k] = u_k (k in I_u),
and tested for strong connectivity by forward and backward closure from one vertex.

Run:  python examples.py            prints all tables
"""
import itertools
from collections import Counter
import numpy as np
import tp_direct as tp


def gval(tab, xi):
    return (tab >> sum(b << k for k, b in enumerate(xi))) & 1


# ---------------------------------------------------------------- periodic points

def F_periodic(tab, r, c):
    """c = (x, y), two tuples of length n (one spatial period)."""
    x, y = c
    n = len(y)
    return (y, tuple(x[i] ^ y[(i - r) % n] ^ gval(tab, [y[(i - r + 1 + k) % n] for k in range(2 * r)])
                     for i in range(n)))


def shift(c, a):
    x, y = c
    n = len(y)
    return (tuple(x[(i + a) % n] for i in range(n)), tuple(y[(i + a) % n] for i in range(n)))


def stabilisers(tab, r, n):
    """Counter {(n,b,s): number of Z^2-orbits of points of least period n with this stabiliser},
    and one representative point for each."""
    seen = set(); count = Counter(); rep = {}
    for bits in itertools.product((0, 1), repeat=2 * n):
        c = (bits[:n], bits[n:])
        if c in seen or any(shift(c, a) == c for a in range(1, n)):
            continue
        forb = [c]; d = F_periodic(tab, r, c)
        while d != c:
            forb.append(d); d = F_periodic(tab, r, d)
        l = len(forb); pos = {d: j for j, d in enumerate(forb)}
        # elements (a, b) of Stab(c) modulo <(n,0),(0,l)>:  sigma^a c = F^j c  <=>  (a, -j) in Stab
        elems = [(a, (-pos[shift(c, a)]) % l) for a in range(n) if shift(c, a) in pos]
        pos_b = [(b, a) for a, b in elems if b > 0]
        b, s = min(pos_b) if pos_b else (l, 0)
        if b == l:
            s = 0
        # consistency: the stabiliser has index n*b, which is the size of the orbit
        orbit = {shift(d, a) for d in forb for a in range(n)}
        assert len(orbit) == n * b and len(elems) * b == l
        assert shift(forb[b % l], s) == c
        seen |= orbit
        count[(n, b, s)] += 1; rep.setdefault((n, b, s), c)
    return count, rep


def in_qZ2_plus_L(q, L, vec):
    n, b, s = L
    return any(((i * n + j * s - vec[0]) % q, (j * b - vec[1]) % q) == (0, 0) for i in range(q) for j in range(q))


def witness(tab, r, q, nmax):
    """first witness of (I)_q: least n, then least index n*b, then least s.
    Returns ((n,b,s), k, all witnesses with that n as a list of ((n,b,s),k)) or None."""
    aq = 0 if q % 2 else r % 2
    for n in range(1, nmax + 1):
        count, _ = stabilisers(tab, r, n)
        ws = sorted((L, k) for L, k in count.items() if k % q and not in_qZ2_plus_L(q, L, (aq, 1)))
        if ws:
            return ws[0][0], ws[0][1], ws
    return None


# ---------------------------------------------------------------- Lambda_M from the layer equations

def lambda_strongly_connected(tab, r, M):
    ny = M + 1                                   # layers y[d], d = -M..0   -> bit d+M
    xe = list(range(-M - r, -r + 1))             # layers x[e]              -> bit ny + index
    N = 1 << (2 * (M + 1))
    v = np.arange(N, dtype=np.int64)
    Y = {d: (v >> (d + M)) & 1 for d in range(-M, 1)}
    X = {e: (v >> (ny + i)) & 1 for i, e in enumerate(xe)}
    Iu = list(range(-M - r, -M))
    By = list(range(-M, min(0, -M + r - 1) + 1))
    T = np.array([(tab >> a) & 1 for a in range(1 << (2 * r))], dtype=np.int64)
    succ = []
    for u in itertools.product((0, 1), repeat=len(Iu)):
        y = dict(Y)
        for k, uk in zip(Iu, u):
            y[k] = np.full(N, uk, dtype=np.int64)
        fixed = np.zeros(N, dtype=np.int64)
        for i, e in enumerate(xe):               # x'[e] = y[e]
            fixed |= y[e] << (ny + i)
        for d in range(-M + r, 1):               # y'[d]
            arg = np.zeros(N, dtype=np.int64)
            for k in range(2 * r):               # xi_k = y[d-1-k]
                arg |= y[d - 1 - k] << k
            fixed |= (y[d] ^ X[d - 2 * r] ^ T[arg]) << (d + M)
        for vv in itertools.product((0, 1), repeat=len(By)):
            w = fixed.copy()
            for d, b in zip(By, vv):
                w |= b << (d + M)
            succ.append(w)
    succ = np.stack(succ, axis=1)
    fw = np.zeros(N, dtype=bool); fw[0] = True
    while True:
        new = fw.copy(); new[succ[fw].ravel()] = True
        if new.sum() == fw.sum(): break
        fw = new
    bw = np.zeros(N, dtype=bool); bw[0] = True
    while True:
        new = bw | bw[succ].any(axis=1)
        if new.sum() == bw.sum(): break
        bw = new
    return bool(fw.all() and bw.all())


# ---------------------------------------------------------------- (T_P), written out

def column_graph(tab, r, P):
    """dict state -> (phi, pi, successors) with states as tuples (s_1..s_2r) of P-tuples."""
    series = list(itertools.product((0, 1), repeat=P))
    G = {}
    for s in itertools.product(series, repeat=2 * r):
        phi = tuple(s[0][(t - 1) % P] ^ gval(tab, [s[2 * r - 1 - k][t] for k in range(2 * r)]) for t in range(P))
        pi = sum(phi) % 2
        succ = []
        if pi == 0:
            for z0 in (0, 1):
                z = [z0]
                for t in range(P - 1):
                    z.append(z[-1] ^ phi[t])
                succ.append(s[1:] + (tuple(z),))
        G[s] = (phi, pi, succ)
    return G


def images(G):
    """list Im_0, Im_1, ... up to and including the first repetition (the last one is S_inf)."""
    ims = [set(G)]
    while True:
        nxt = {t for s in ims[-1] for t in G[s][2]}
        if nxt == ims[-1]:
            return ims
        ims.append(nxt)


def shortest_certificate(G):
    """a cycle of the successor relation followed by a chain of successors to a state with
    pi = 1, with the least total number of states (cycle states + chain states after the cycle).
    Returns (cycle, chain) where chain starts at a state of the cycle."""
    best = None
    states = sorted(G)
    def dist_from(src):
        dist = {src: [src]}; queue = [src]
        while queue:
            a = queue.pop(0)
            for b in G[a][2]:
                if b not in dist:
                    dist[b] = dist[a] + [b]; queue.append(b)
        return dist
    paths = {s: dist_from(s) for s in states}
    for s in states:
        # shortest cycle through s
        cyc = None
        for t in G[s][2]:
            if s in paths[t]:
                c = [s] + paths[t][s][:-1] if t != s else [s]
                if cyc is None or len(c) < len(cyc): cyc = c
        if cyc is None: continue
        for e in states:
            if G[e][1] == 1 and e in paths[s]:
                chain = paths[s][e]
                cost = len(cyc) + len(chain) - 1
                if best is None or cost < best[0]:
                    best = (cost, cyc, chain)
    return best


def fmt(s):
    return "(" + ",".join("".join(map(str, sj)) for sj in s) + ")"


def print_certificate(name, tab, r, P, full=False):
    G = column_graph(tab, r, P)
    ims = images(G)
    S = ims[-1]
    print("%s, P=%d: |Im_m| = %s, |S_inf| = %d, states of S_inf with pi=1: %s" % (
        name, P, [len(i) for i in ims], len(S), " ".join(fmt(s) for s in sorted(S) if G[s][1])))
    if full:
        for s in sorted(G):
            phi, pi, succ = G[s]
            print("    s=%s phi=%s pi=%d successors: %s" % (fmt(s), "".join(map(str, phi)), pi,
                                                         " ".join(fmt(t) for t in succ) or "none"))
        for m, im in enumerate(ims):
            print("    Im_%d = {%s}" % (m, " ".join(fmt(s) for s in sorted(im))))
    best = shortest_certificate(G)
    if best is None:
        print("    no certificate: (T_%d) fails" % P); return
    cost, cyc, chain = best
    print("    shortest certificate (%d states): cycle %s -> back to %s ; chain %s" % (
        cost, " -> ".join(fmt(s) for s in cyc), fmt(cyc[0]), " -> ".join(fmt(s) for s in chain)))
    for s in cyc + chain[1:]:
        phi, pi, succ = G[s]
        print("      s=%s phi=%s pi=%d successors: %s" % (fmt(s), "".join(map(str, phi)), pi,
                                                       " ".join(fmt(t) for t in succ) or "none"))


# ---------------------------------------------------------------- tables

def from_anf(n, monomials):
    tab = 0
    for x in range(1 << n):
        v = 0
        for m in monomials:
            v ^= all((x >> k) & 1 for k in m)
        tab |= v << x
    return tab


def eca_g(rule):
    """g(b,c) = f(0,b,c) for the elementary rule f = a XOR g(b,c); xi_0 = b, xi_1 = c."""
    assert all(((rule >> (2 * b + c)) ^ (rule >> (4 + 2 * b + c))) & 1 for b in (0, 1) for c in (0, 1))
    return sum(((rule >> (2 * b + c)) & 1) << (b + 2 * c) for b in (0, 1) for c in (0, 1))


RANDOM_R2 = [0xad44, 0x34ff, 0x3f7b, 0xbaec, 0xfa9b, 0x5ff8, 0x5980, 0x4f98, 0xe98e, 0x6152]
RANDOM_R3 = [0x02f8aca0c96b3825, 0x9207dbf3d68f59b4, 0x6241805bcece3963,
             0xdfc519edc1a14412, 0x822cd8c8f65243a3, 0x8557d3586049b4bb]
maj = [(0, 2), (0, 3), (2, 3)]
NAMED_R2 = [("xi0xi2", [(0, 2)]), ("xi0xi3", [(0, 3)]), ("xi3+xi0xi1xi2", [(3,), (0, 1, 2)]),
            ("xi0xi1", [(0, 1)]), ("xi1xi2", [(1, 2)]), ("maj(xi0,xi2,xi3)", maj),
            ("xi0+xi3", [(0,), (3,)]), ("xi1", [(1,)]), ("xi3", [(3,)])]


def row(name, tab, r, Pmax, nmax, digits):
    n = 2 * r
    lp = tp.least_P([tab], r, Pmax)[tab]
    out = "%-18s tab=0x%0*x J=%s self-dual=%s least P<=%d: %s" % (
        name, digits, tab, tp.J(tab, n), "yes" if tp.self_dual(tab, n) else "no", Pmax, lp)
    for q in sorted({p for p in (2, 3, 5, 7) if (2 * r) % p == 0}):
        w = witness(tab, r, q, nmax[q] if isinstance(nmax, dict) else nmax)
        nm = nmax[q] if isinstance(nmax, dict) else nmax
        out += " | (I)_%d: %s" % (q, "none with n<=%d" % nm if w is None else
                                  "(n,b,s)=%s k=%d [all witnesses with this n: %s]" % (
                                      w[0], w[1], " ".join("%s:k=%d" % (L, k) for L, k in w[2])))
    return out


if __name__ == "__main__":
    print("== radius 1 ==")
    for rule in (30, 45, 75, 135, 120, 180, 210, 225, 90, 150):
        tab = eca_g(rule)
        print(row("rule %d" % rule, tab, 1, 5, 4, 1))
    print("rule 90: (T_P), P=1..8:", [bool(tp.Columns(1, P).holds(tp.tables_from_ints([eca_g(90)], 2))[0]) for P in range(1, 9)])
    print("rule 90: Lambda_M strongly connected, M=0..8:", [lambda_strongly_connected(eca_g(90), 1, M) for M in range(9)])
    print("rule 150: Lambda_M strongly connected, M=0..8:", [lambda_strongly_connected(eca_g(150), 1, M) for M in range(9)])
    for n in range(1, 9):
        c, _ = stabilisers(eca_g(150), 1, n)
        ws = [(L, k) for L, k in sorted(c.items()) if k % 2 and not in_qZ2_plus_L(2, L, (1, 1))]
        print("rule 150: n=%d stabilisers (n,b,s):k = %s ; (I)_2 witnesses: %s" % (n, dict(sorted(c.items())), ws or "none"))
    print("== radius 2, named rules ==")
    for name, ms in NAMED_R2:
        tab = from_anf(4, ms)
        print(row(name, tab, 2, 4, 5, 4), "| ANF:", tp.anf_str(tab, 4))
    t = from_anf(4, [(0,), (3,)])
    for n in range(1, 8):
        c, _ = stabilisers(t, 2, n)
        ws = [(L, k) for L, k in sorted(c.items()) if k % 2 and not in_qZ2_plus_L(2, L, (0, 1))]
        print("xi0+xi3: n=%d stabilisers (n,b,s):k = %s ; (I)_2 witnesses: %s" % (n, dict(sorted(c.items())), ws or "none"))
    t = from_anf(4, [(3,)])
    print("xi3: (T_P), P=1..5:", [bool(tp.Columns(2, P).holds(tp.tables_from_ints([t], 4))[0]) for P in range(1, 6)])
    print("xi3: Lambda_M strongly connected, M=0..5:", [lambda_strongly_connected(t, 2, M) for M in range(6)])
    for name, ms in NAMED_R2[:6] + NAMED_R2[7:8]:
        print("%s: Lambda_M strongly connected, M=0..5:" % name, [lambda_strongly_connected(from_anf(4, ms), 2, M) for M in range(6)])
    print("== radius 2, ten further rules ==")
    for i, tab in enumerate(RANDOM_R2):
        print(row("#%d" % i, tab, 2, 4, 5, 4), "| ANF:", tp.anf_str(tab, 4))
    print("== radius 3, six rules ==")
    for i, tab in enumerate(RANDOM_R3):
        print(row("#%d" % i, tab, 3, 3, {2: 6, 3: 4}, 16))
    print("== certificates, radius 1 ==")
    for rule in (30, 45, 75, 135):
        print_certificate("rule %d (g tab 0x%x)" % (rule, eca_g(rule)), eca_g(rule), 1, 1, full=True)
    print_certificate("rule 135", eca_g(135), 1, 2, full=True)
