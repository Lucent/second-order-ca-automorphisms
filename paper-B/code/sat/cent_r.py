"""SAT centralizer for a radius-r second-order rule F(x,y)=(y, x^f(y)) among ALL block maps of radius R
on A={0,1}^2 (centralizer.py is the same computation for r=1).  Unknown: table h[e] = (X_e, Y_e) for
e in A^{2R+1}.

Commutation h o F = F o h on every window w in A^{2R+2r+1}:
   X_{eL} = Y_{e_r}
   Y_{eL} ^ X_{e_r} ^ Y_{e_0} = g(Y_{e_1}, ..., Y_{e_{2r}})
with eL = index of F(w) (a (2R+1)-word), e_k = index of w[k:k+2R+1].  g is Tseitin-encoded with one aux
variable per distinct tuple (e_1..e_{2r}).  All solutions are enumerated; each is re-verified directly,
compared with the maps sigma^a F^b (built by composing local rules; the comparison is exact), and, if it
is not one of them, tested for injectivity (pair graph).  The program asserts that every sigma^a F^b with
|a| + r|b| <= R occurs among the solutions.

Usage:  python cent_r.py ECA<n> R        (elementary rule n, left-permutive)
        python cent_r.py <gtab> r R      (gtab = the 2^(2r) values of g as a string of 0/1, xi_0 most
                                          significant; for example 0000000000110011 2 1 for g = xi_0 xi_2)
Output: one JSON record: nsol = number of commuting block maps of radius <= R, named_found = the pairs
(a,b) found, extra_automorphisms = number of injective solutions not of the form sigma^a F^b,
non_injective = number of non-injective solutions, constant = number of constant maps among them.
"""
import sys, time, itertools, json
import numpy as np
from pysat.solvers import Solver
from common import Rule, elementary

def word_index(w):
    e = 0
    for s in w: e = (e << 2) | s
    return e

def index_word(e, n):
    return tuple((e >> (2 * (n - 1 - i))) & 3 for i in range(n))

def apply_loc(loc, word, r):
    k = 2 * r + 1
    return tuple(loc(tuple(word[i:i + k])) for i in range(len(word) - k + 1))

def add_truth(solver, roles, ok):
    """add clauses forbidding assignments of the distinct variables in `roles` for which ok(values) is False."""
    vs = sorted(set(roles)); n = 0
    for bits in itertools.product((0, 1), repeat=len(vs)):
        val = dict(zip(vs, bits))
        if not ok([val[v] for v in roles]):
            solver.add_clause([-v if val[v] else v for v in vs]); n += 1
    return n

def encode(rule, R, solver):
    r = rule.r; n = 2 * R + 1; m = 2 * R + 2 * r + 1
    X = lambda e: 2 * e + 1; Y = lambda e: 2 * e + 2
    base = 2 * 4 ** n
    aux = {}; ncl = 0
    gl = lambda vals: rule.gval(vals)
    for ew in range(4 ** m):
        w = index_word(ew, m)
        eL = word_index(apply_loc(rule.F_cell, w, r))
        es = [word_index(w[k:k + n]) for k in range(2 * r + 1)]
        er = es[r]
        solver.add_clause([-X(eL), Y(er)]); solver.add_clause([X(eL), -Y(er)]); ncl += 2
        key = tuple(es[1:])
        if key not in aux:
            gv = base + 1 + len(aux); aux[key] = gv
            roles = [gv] + [Y(e) for e in key]
            ncl += add_truth(solver, roles, lambda a: a[0] == gl(a[1:]))
        roles = [Y(eL), X(er), Y(es[0]), aux[key]]
        ncl += add_truth(solver, roles, lambda a: (a[0] ^ a[1] ^ a[2] ^ a[3]) == 0)
    return base, ncl

def solve_all(rule, R, maxsol=10000):
    n = 2 * R + 1; nvar = 2 * 4 ** n
    with Solver(name="cadical153") as s:
        t0 = time.time(); _, ncl = encode(rule, R, s); te = time.time() - t0
        sols = []
        while s.solve():
            model = s.get_model()
            val = [1 if lit > 0 else 0 for lit in model[:nvar]]
            sols.append(tuple((val[2 * e] << 1) | val[2 * e + 1] for e in range(4 ** n)))
            s.add_clause([-lit for lit in model[:nvar]])
            if len(sols) >= maxsol: break
    return sols, ncl, te

def check_commutes(rule, R, tab):
    r = rule.r; n = 2 * R + 1; m = 2 * R + 2 * r + 1
    h = lambda word: tuple(tab[word_index(word[i:i + n])] for i in range(len(word) - n + 1))
    for ew in range(4 ** m):
        w = index_word(ew, m)
        if h(apply_loc(rule.F_cell, w, r)) != apply_loc(rule.F_cell, h(w), r): return False
    return True

def injective(R, tab):
    """pair-graph injectivity test on A^Z (as in centralizer.py), numpy version."""
    L = 2 * R; N = 4 ** L; mask = N - 1
    T = np.array(tab, dtype=np.int64)
    if L == 0:
        return len(set(tab)) == 4
    u = np.arange(N, dtype=np.int64)
    # successor lists: for nodes (u,v), symbols (a,b) with T[u a] == T[v b]
    src = []; dst = []
    ua = (u[:, None] << 2) | np.arange(4)[None, :]           # N x 4
    Tua = T[ua]
    for a in range(4):
        for b in range(4):
            eq = Tua[:, a][:, None] == Tua[:, b][None, :]       # N x N
            iu, iv = np.nonzero(eq)
            src.append(iu * N + iv)
            dst.append(((ua[iu, a] & mask) * N + (ua[iv, b] & mask)))
    src = np.concatenate(src); dst = np.concatenate(dst)
    alive = np.ones(N * N, dtype=bool)
    while True:
        e = alive[src] & alive[dst]
        has_out = np.zeros(N * N, dtype=bool); has_out[src[e]] = True
        has_in = np.zeros(N * N, dtype=bool); has_in[dst[e]] = True
        new = alive & has_out & has_in
        if new.sum() == alive.sum(): break
        alive = new
    diag = (np.arange(N * N) // N) == (np.arange(N * N) % N)
    return not np.any(alive & ~diag)

def named(rule, R):
    """tables of sigma^a F^b of radius <= R (radius of F^b is r|b|)."""
    r = rule.r; n = 2 * R + 1; out = {}
    for b in range(-(R // r), R // r + 1):
        rb = r * abs(b)
        for a in range(-(R - rb), R - rb + 1):
            tab = []
            for e in range(4 ** n):
                w = index_word(e, n)
                for _ in range(abs(b)):
                    w = apply_loc(rule.F_cell if b > 0 else rule.Finv_cell, w, r)
                # w now has length n - 2 rb, centre index R - rb ; sigma^a: output c_{i+a}
                tab.append(w[R - rb + a])
            out[tuple(tab)] = (a, b)
    return out

def identify(rule, R, tab, trials=6, n=48, seed=0):
    """return (a,b) if h = sigma^a F^b on random rings (|a|<=R+2r, |b|<=3), else None."""
    import random
    rng = random.Random(seed); r = rule.r
    rings = [tuple(rng.randrange(4) for _ in range(n)) for _ in range(trials)]
    def ring_apply(loc, c, rad):
        return tuple(loc(tuple(c[(i + k) % n] for k in range(-rad, rad + 1))) for i in range(n))
    hc = [tuple(tab[word_index(tuple(c[(i + k) % n] for k in range(-R, R + 1)))] for i in range(n)) for c in rings]
    for b in range(-3, 4):
        Fb = []
        for c in rings:
            d = c
            for _ in range(abs(b)): d = ring_apply(rule.F_cell if b > 0 else rule.Finv_cell, d, r)
            Fb.append(d)
        for a in range(-(R + 2 * r), R + 2 * r + 1):
            if all(tuple(d[(i + a) % n] for i in range(n)) == h for d, h in zip(Fb, hc)): return (a, b)
    return None

def equals_named(rule, R, tab, a, b):
    """exact test: the block map with table tab equals sigma^a F^b (comparison on all words of length 2 rho + 1,
    rho = max(R, |a| + r|b|))."""
    r = rule.r; rho = max(R, abs(a) + r * abs(b)); n = 2 * R + 1
    for w in itertools.product(range(4), repeat=2 * rho + 1):
        v = w
        for _ in range(abs(b)):
            v = apply_loc(rule.F_cell if b > 0 else rule.Finv_cell, v, r)
        if v[rho - r * abs(b) + a] != tab[word_index(w[rho - R:rho + R + 1])]: return False
    return True

def classify(rule, R, maxsol=10000):
    sols, ncl, te = solve_all(rule, R, maxsol)
    nm = dict(named(rule, R))
    if len(sols) < maxsol:
        missing = [ab for tab, ab in nm.items() if tab not in set(sols)]
        assert not missing, f"sigma^a F^b missing from the solutions: {missing}"
    # sigma^a F^b can have radius <= R without |a| + r|b| <= R (when g does not depend on xi_{2r-1}):
    # candidates are found on random rings and then compared exactly, as local rules.
    for tab in sols:
        if tab not in nm:
            ab = identify(rule, R, tab)
            if ab is not None and equals_named(rule, R, tab, *ab): nm[tab] = ab
    rep = dict(rule=rule.name, r=rule.r, R=R, clauses=ncl, nsol=len(sols), capped=len(sols) >= maxsol)
    autos = []; known = []; endo = 0; bad = 0
    for tab in sols:
        if not check_commutes(rule, R, tab): bad += 1; continue
        if tab in nm: known.append(nm[tab]); continue
        if injective(R, tab): autos.append(tab)
        else: endo += 1
    rep.update(verified_fail=bad, named_found=sorted(known), extra_automorphisms=len(autos), non_injective=endo,
               constant=sum(1 for tab in sols if len(set(tab)) == 1))
    # describe extra automorphisms of radius 0 (cellwise permutations)
    rep["extra_examples"] = []
    for tab in autos[:4]:
        cellwise = all(tab[e] == tab[(e & 3) | ((e >> 2) << 2)] for e in range(len(tab))) if R == 0 else None
        rep["extra_examples"].append("".join(map(str, tab[:64])) + ("..." if len(tab) > 64 else ""))
    return rep, autos

if __name__ == "__main__":
    spec = sys.argv[1]
    if spec.startswith("ECA"): rule = elementary(int(spec[3:])); R = int(sys.argv[2])
    else: rule = Rule(int(sys.argv[2]), [int(c) for c in spec]); R = int(sys.argv[3])
    rep, _ = classify(rule, R)
    print(json.dumps(rep))
