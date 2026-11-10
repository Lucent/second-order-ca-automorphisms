"""SAT search for nonzero solutions gamma of the stripping equation (*) (general r):
  gamma(G^ th) + gamma(H^ th) + gamma(th) = g(xi_0..xi_{2r-1}) + g(xi_0+gamma(s th), ..., xi_{2r-1}+gamma(s^{2r} th))
  G^(th)_k = (zeta_k+xi_k+g(xi_{k+1..k+2r}), xi_{k+2r}),  H^(th)_k = (zeta_k, xi_{k+2r}+zeta_k+g(zeta_{k+1..k+2r}))
gamma depends on th_0..th_m, th_k = (xi_k, zeta_k); its truth table is indexed by the base-4 number with
digits 2*xi_k + zeta_k, k = 0..m (th_0 most significant).  s is the shift of th.
The clause "gamma is not identically 0" is added and the enumeration stops after 20 solutions.
Usage: python strip.py <gtab> r m     (gtab as in cent_r.py)
Output: J(g), self-duality, and the number of solutions, each marked const1 (gamma = 1) or nonconst.
"""
import sys, itertools
from pysat.solvers import Solver
from common import Rule

def run(rule, m, want="nonzero"):
    r = rule.r; n2 = 2 * r; L = m + n2 + 1
    var = lambda idx: idx + 1          # gamma table index: sum over k<=m of (2*xi_k+zeta_k) << 2(m-k)
    def gidx(th):                      # th: list of (xi,zeta), use first m+1
        e = 0
        for k in range(m + 1): e = (e << 2) | (2 * th[k][0] + th[k][1])
        return e
    s = Solver(name="cadical153"); ncl = 0
    for code in range(4 ** L):
        th = [((code >> (2 * (L - 1 - k) + 1)) & 1, (code >> (2 * (L - 1 - k))) & 1) for k in range(L)]
        xi = [a for a, _ in th]; ze = [b for _, b in th]
        G = [(ze[k] ^ xi[k] ^ rule.gval(xi[k + 1:k + 1 + n2]), xi[k + n2]) for k in range(m + 1)]
        H = [(ze[k], xi[k + n2] ^ ze[k] ^ rule.gval(ze[k + 1:k + 1 + n2])) for k in range(m + 1)]
        vG, vH, v0 = var(gidx(G)), var(gidx(H)), var(gidx(th))
        vs_shift = [var(gidx(th[j + 1:])) for j in range(n2)]
        roles = [vG, vH, v0] + vs_shift
        base = rule.gval(xi[:n2])
        vs = sorted(set(roles))
        for bits in itertools.product((0, 1), repeat=len(vs)):
            val = dict(zip(vs, bits)); a = [val[v] for v in roles]
            lhs = a[0] ^ a[1] ^ a[2]
            rhs = base ^ rule.gval([xi[j] ^ a[3 + j] for j in range(n2)])
            if lhs != rhs:
                s.add_clause([-v if val[v] else v for v in vs]); ncl += 1
    N = 4 ** (m + 1)
    s.add_clause([var(i) for i in range(N)])            # gamma not identically 0
    sols = []
    while s.solve() and len(sols) < 20:
        mod = s.get_model()[:N]
        sols.append(tuple(1 if x > 0 else 0 for x in mod))
        s.add_clause([-x for x in mod])
    return sols

if __name__ == "__main__":
    rule = Rule(int(sys.argv[2]), [int(c) for c in sys.argv[1]]); m = int(sys.argv[3])
    sols = run(rule, m)
    print(rule.name, rule.anf_str(), "J=", rule.J(), "selfdual=", rule.self_dual(), "m=", m,
          "nonzero solutions:", len(sols), ["const1" if all(sols_i) else "nonconst" for sols_i in sols])
