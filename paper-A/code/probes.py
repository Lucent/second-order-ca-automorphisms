"""probes.py -- the moves computed from real probes, and the bounce formula.

For f in {110, 42, 46}: exact second-order dynamics on configurations that are
constant outside a finite interval.  For random 0-tailed configurations b with edge r
and times |m| <= 6 it computes the probe P_n(b) = F^{-n}(F^n b |_r rho eta_{r+2}),
the separation of the outgoing kink and the debris D_n(b) (Lemma cl:lem:probe), and the
moves M_m(b) = D_{m-2}(b), Mbar_m(b) = rho D_{-m-2}(rho b), and checks:
  (P) at s = n+2+eps, F^s P_n(b) = F^s D |_r eta_{r+s-n-1-eps}, with D 0-tailed with edge r,
      and the same D is obtained at s+5;
  (T) F^m M_m(b), F^m Mbar_m(b) differ from F^m b at most in the edge cell, changed as in
      the move table (Lemma cl:lem:edgemove for 110, Lemma cl:lem:moves1446 for 42, 46);
  (C) M_m(b), Mbar_m(b) agree with b outside [r-|m|, r] (Lemma cl:lem:movelocal);
  (B) the bounce formula U(d |_e rho eta_N) = 2e - N + Y_d(N-e-1, e) (Lemma cl:lem:bounce(c)).
Cells are integers 2x+y.  Usage: python3 probes.py [ntests]
"""
import random, sys, time

def rulefn(n):
    T = [(n >> j) & 1 for j in range(8)]
    return lambda a, b, c: T[4 * a + 2 * b + c]

class Conf:
    """cells a[0..] at positions lo.., constant cell lt to the left and rt to the right"""
    __slots__ = ('lo', 'a', 'lt', 'rt')
    def __init__(s, lo, a, lt=0, rt=0):
        s.lo, s.a, s.lt, s.rt = lo, list(a), lt, rt
    def at(s, i):
        j = i - s.lo
        if j < 0: return s.lt
        if j >= len(s.a): return s.rt
        return s.a[j]
    def hi(s): return s.lo + len(s.a) - 1
    def eq(s, o, lo, hi): return all(s.at(i) == o.at(i) for i in range(lo, hi + 1))

def step(f, C, inv=False):
    lo, hi = C.lo - 1, C.hi() + 1
    def img(c, l, r):
        x, y = c >> 1, c & 1
        if not inv: return 2 * y + (x ^ f(l & 1, y, r & 1))
        return 2 * (y ^ f(l >> 1, x, r >> 1)) + x
    out = [img(C.at(i), C.at(i - 1), C.at(i + 1)) for i in range(lo, hi + 1)]
    return Conf(lo, out, img(C.lt, C.lt, C.lt), img(C.rt, C.rt, C.rt))

def Fp(f, C, n):
    for _ in range(abs(n)): C = step(f, C, n < 0)
    return C

def rho(C):
    sw = lambda c: 2 * (c & 1) + (c >> 1)
    return Conf(C.lo, [sw(c) for c in C.a], sw(C.lt), sw(C.rt))

def edge(C):
    assert C.rt == 0
    for i in range(C.hi(), C.lo - 1, -1):
        if C.at(i): return i
    return None

def ycol(f, d, t, i):
    """Y_d(t,i): y-bit of F^t d at i (t >= -1); row -1 is the x-row of d"""
    if t == -1: return d.at(i) >> 1
    return Fp(f, d, t).at(i) & 1

def glue(C, r, K):
    """C |_r K : C on (-inf, r], K on (r, inf)"""
    hi = max(C.hi(), K.hi()) + 2
    return Conf(C.lo, [C.at(i) if i <= r else K.at(i) for i in range(C.lo, hi + 1)], C.lt, K.rt)

def eta(p):  return Conf(p, [2], 0, 3)
def rhoeta(p): return Conf(p, [1], 0, 3)
def cut(C, r):
    return Conf(C.lo, [C.at(i) for i in range(C.lo, r + 1)], C.lt, 0)

def debris(f, b, n):
    """D_n(b) and eps; checks the separation statement (P)"""
    r = edge(b)
    P = Fp(f, glue(Fp(f, b, n), r, rhoeta(r + 2)), -n)
    eps = 1 if (Fp(f, b, n + 1).at(r) & 1) == 0 else 0
    res = []
    for s in (n + 2 + eps, n + 7 + eps):
        G = Fp(f, P, s)
        N = r + s - n - 1 - eps
        assert G.eq(eta(N), r + 1, G.hi() + 3), 'outgoing kink not separated'
        D = Fp(f, cut(G, r), -s)
        assert edge(D) == r, 'debris edge'
        res.append(D)
    lo = min(res[0].lo, res[1].lo) - 2
    assert res[0].eq(res[1], lo, r + 3), 'debris depends on the separation time'
    return res[0]

def M(f, b, m): return debris(f, b, m - 2)
def Mbar(f, b, m): return rho(debris(f, rho(b), -m - 2))

def table(rule, xr, yr, xl, yl, bar):
    if not bar:
        if (xr, yr) == (0, 1): return (1, 0)
        if xr == 1 and rule == 110: return (1, yr ^ xl)
        if xr == 1 and rule == 42: return (1, yr ^ 1 ^ xl)
    else:
        if (xr, yr) == (1, 0): return (0, 1)
        if yr == 1 and rule == 110: return (xr ^ yl, 1)
        if yr == 1 and rule == 42: return (xr ^ 1 ^ yl, 1)
    return (xr, yr)

def random_tailed(rng, r, width):
    cells = [rng.randrange(4) for _ in range(width - 1)] + [rng.randrange(1, 4)]
    return Conf(r - width + 1, cells, 0, 0)

def U(f, c, t0, span=6):
    """eventual value of mu_t - t, mu_t the largest i with Y_c(t,i)=0; checked constant"""
    vals = []
    G = Fp(f, c, t0 - 1)
    for t in range(t0, t0 + span):
        G = step(f, G)
        mu = max(i for i in range(G.lo - 1, G.hi() + 1) if (G.at(i) & 1) == 0)
        vals.append(mu - t)
    assert len(set(vals)) == 1, ('no stable front', vals)
    return vals[0]

def main(ntests):
    rng = random.Random(2026)
    t0 = time.time(); ok = True
    for rule in (110, 42, 46):
        f = rulefn(rule); nt = nc = nb = 0
        for _ in range(ntests):
            r = 0; b = random_tailed(rng, r, rng.randrange(2, 12))
            m = rng.randrange(-6, 7)
            S = Fp(f, b, m)
            for bar in (0, 1):
                b2 = (Mbar if bar else M)(f, b, m)
                assert edge(b2) == r
                S2 = Fp(f, b2, m)
                lo = min(S.lo, S2.lo) - 2
                # (T)
                if not S2.eq(S, lo, r - 1): ok = False; print('FAIL T (cells)', rule, m)
                if not S2.eq(S, r + 1, r + 5): ok = False; print('FAIL T (right)', rule, m)
                xr, yr = S.at(r) >> 1, S.at(r) & 1
                xl, yl = S.at(r - 1) >> 1, S.at(r - 1) & 1
                want = table(rule, xr, yr, xl, yl, bar)
                if (S2.at(r) >> 1, S2.at(r) & 1) != want: ok = False; print('FAIL T (edge)', rule, m)
                nt += 1
                # (C)
                if not b2.eq(b, min(b.lo, b2.lo) - 3 * abs(m) - 5, r - abs(m) - 1): ok = False; print('FAIL C', rule, m)
                nc += 1
        for _ in range(max(ntests // 4, 1)):
            e = 0; d = random_tailed(rng, e, rng.randrange(1, 10))
            N = e + 1 + rng.randrange(0, 12)
            c = glue(d, e, rhoeta(N))
            want = 2 * e - N + ycol(f, d, N - e - 1, e)
            got = U(f, c, N - e + 4)
            if got != want: ok = False; print('FAIL B', rule, N, got, want)
            nb += 1
        print('rule %d: %d move tests (table and edge), %d cone tests, %d bounce tests' % (rule, nt, nc, nb), flush=True)
    print('time %.1f s' % (time.time() - t0))
    if not ok: raise SystemExit('FAILURE')
    print('ALL CHECKS PASSED')

if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 1500)
