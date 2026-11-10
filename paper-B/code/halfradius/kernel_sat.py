"""Right-sided block maps commuting with F, for a rule f(a,b,c) = g(b,c) that does not read its left
neighbour (Section 6 of the paper).

Usage:  python kernel_sat.py RULE A1,A2,...        for example  python kernel_sat.py 187 1,2,3,4,5,6

Cells are z = 2x + y, and F(c)_i = phi_{y(c_{i+1})}(c_i) with phi_nu(x,y) = (y, x XOR g(y,nu)).
For each anticipation A the unknown is the table K : cells^(A+1) -> cells of a map k(c)_i = K(c_i, ..., c_{i+A}).
The clauses say that
  (1) k F = F k, on every word of length A+2;
  (2) K is bijective in its first argument when the other A arguments are fixed.
Every right-sided automorphism whose inverse is right-sided satisfies (2).  All solutions are enumerated
with blocking clauses (pysat, CaDiCaL), each is compared with the tables of F^b, |b| <= A, and each
solution that is not a power of F is tested for injectivity on the configurations of period n <= 6.

Output, for each A: the number of solutions, the exponents b of the powers F^b found, and the numbers of
other solutions that are injective / not injective on the configurations of period n <= 6.
"""
import sys, itertools, time
from pysat.solvers import Solver


def g_of(rule):
    f = lambda a, b, c: (rule >> (4 * a + 2 * b + c)) & 1
    assert all(f(0, b, c) == f(1, b, c) for b in (0, 1) for c in (0, 1)), "the rule reads its left neighbour"
    return lambda b, c: f(0, b, c)


def phis(rule):
    g = g_of(rule)
    return {nu: tuple((y << 1) | (x ^ g(y, nu)) for z in range(4) for x, y in [(z >> 1, z & 1)]) for nu in (0, 1)}


def idx(w):
    e = 0
    for s in w:
        e = (e << 2) | s
    return e


def solutions(rule, A, maxsol=100000):
    g = g_of(rule); P = phis(rule); n = A + 1; N = 4 ** n
    X = lambda e: 2 * e + 1
    Y = lambda e: 2 * e + 2
    S = Solver(name='cadical153')
    for W in range(4 ** (n + 1)):
        w = [(W >> (2 * (n - i))) & 3 for i in range(n + 1)]
        Fw = [P[w[j + 1] & 1][w[j]] for j in range(n)]
        eL = idx(Fw); e0 = idx(w[:n]); e1 = idx(w[1:])
        # x(K(F w)) = y(K(w_0..w_A))
        S.add_clause([-X(eL), Y(e0)]); S.add_clause([X(eL), -Y(e0)])
        # y(K(F w)) = x(K(w_0..)) XOR g(y(K(w_0..)), y(K(w_1..)))
        vs = sorted(set([Y(eL), X(e0), Y(e0), Y(e1)]))
        for vals in itertools.product((0, 1), repeat=len(vs)):
            a = dict(zip(vs, vals))
            if a[Y(eL)] != (a[X(e0)] ^ g(a[Y(e0)], a[Y(e1)])):
                S.add_clause([-v if a[v] else v for v in vs])
    for u in range(4 ** A):
        es = [(s << (2 * A)) | u for s in range(4)]
        for e1, e2 in itertools.combinations(es, 2):
            for xv in (0, 1):
                for yv in (0, 1):
                    S.add_clause([(-1 if xv else 1) * X(e1), (-1 if xv else 1) * X(e2),
                                  (-1 if yv else 1) * Y(e1), (-1 if yv else 1) * Y(e2)])
    out = []
    while S.solve() and len(out) < maxsol:
        m = S.get_model()
        tab = tuple(((m[X(e) - 1] > 0) << 1) | (m[Y(e) - 1] > 0) for e in range(N))
        out.append(tab)
        S.add_clause([-X(e) if tab[e] >> 1 else X(e) for e in range(N)] + [-Y(e) if tab[e] & 1 else Y(e) for e in range(N)])
    S.delete()
    return out


def commutes(rule, A, tab):
    """direct check of k F = F k on all words of length A+2."""
    P = phis(rule); n = A + 1
    for w in itertools.product(range(4), repeat=n + 1):
        Fw = [P[w[j + 1] & 1][w[j]] for j in range(n)]
        k0 = tab[idx(w[:n])]; k1 = tab[idx(w[1:])]
        if tab[idx(Fw)] != P[k1 & 1][k0]:
            return False
    return True


def power_tables(rule, A):
    """tables of F^b, |b| <= A, as maps of anticipation A."""
    g = g_of(rule); P = phis(rule); n = A + 1; out = {}
    for b in range(-A, A + 1):
        tab = []
        for w in itertools.product(range(4), repeat=n):
            w = list(w)
            for _ in range(abs(b)):
                if b > 0:
                    w = [P[w[j + 1] & 1][w[j]] for j in range(len(w) - 1)]
                else:   # F^{-1}(d)_i = (y XOR g(x, x(d_{i+1})), x)
                    w = [(((w[j] & 1) ^ g(w[j] >> 1, w[j + 1] >> 1)) << 1) | (w[j] >> 1) for j in range(len(w) - 1)]
            tab.append(w[0])
        out[tuple(tab)] = b
    return out


def injective_periodic(tab, A, nmax=6):
    for n in range(1, nmax + 1):
        seen = set()
        for c in itertools.product(range(4), repeat=n):
            img = tuple(tab[idx([c[(i + j) % n] for j in range(A + 1)])] for i in range(n))
            if img in seen:
                return False
            seen.add(img)
    return True


if __name__ == '__main__':
    rule = int(sys.argv[1])
    for A in map(int, sys.argv[2].split(',')):
        t = time.time()
        sols = solutions(rule, A)
        known = power_tables(rule, A)
        assert all(commutes(rule, A, tab) for tab in (sols if A <= 4 else sols[:50]))
        assert all(tab in set(sols) for tab in known), "a power of F is missing from the solutions"
        powers = sorted(known[tab] for tab in sols if tab in known)
        others = [tab for tab in sols if tab not in known]
        inj = sum(1 for tab in others if injective_periodic(tab, A))
        print(f"rule {rule}, anticipation <= {A}: {len(sols)} solutions; powers F^b: b = {powers}; "
              f"other solutions: {len(others)} ({inj} injective on periods <= 6, {len(others) - inj} not injective); "
              f"{time.time() - t:.1f}s", flush=True)
