"""Common utilities: left-permutive second-order rules of radius r over GF(2).

f(y_{i-r},...,y_{i+r}) = y_{i-r} XOR g(xi_0,...,xi_{2r-1}),   xi_k = y_{i-r+1+k}.
gtab[idx] with idx = sum_k xi_k << (2r-1-k)   (xi_0 most significant).
Cells are encoded s = 2x + y (x = previous row, y = current row), as in centralizer.py.
F(c)_i = (y_i, x_i ^ f(y_{i-r..i+r})).
"""
import itertools, random

class Rule:
    def __init__(self, r, gtab, name=None):
        self.r = r; self.n = 2 * r
        assert len(gtab) == 1 << self.n
        self.g = tuple(int(v) & 1 for v in gtab)
        self.name = name or ("g=" + "".join(map(str, self.g)))
        # f table on 2r+1 bits, index = sum y_{i-r+k} << (2r-k)
        self.ftab = [((e >> self.n) & 1) ^ self.g[e & ((1 << self.n) - 1)] for e in range(1 << (self.n + 1))]

    def gval(self, xs):
        e = 0
        for b in xs: e = (e << 1) | b
        return self.g[e]

    def f(self, ys):          # ys = (y_{i-r}, ..., y_{i+r})
        return ys[0] ^ self.gval(ys[1:])

    # ---------- algebraic invariants
    def anf(self):
        """ANF coefficients: dict monomial(frozenset of indices k of xi_k) -> 1."""
        n = self.n; a = list(self.g)
        # Moebius transform over index bits (bit (n-1-k) <-> xi_k)
        for i in range(n):
            for e in range(1 << n):
                if e >> i & 1: a[e] ^= a[e ^ (1 << i)]
        out = set()
        for e in range(1 << n):
            if a[e]: out.add(frozenset(k for k in range(n) if e >> (n - 1 - k) & 1))
        return out

    def is_affine(self):
        return all(len(m) <= 1 for m in self.anf())

    def deriv2_nonzero(self, i, j):
        """is d_i d_j g not identically 0 ?"""
        n = self.n; bi = 1 << (n - 1 - i); bj = 1 << (n - 1 - j)
        return any(self.g[e] ^ self.g[e ^ bi] ^ self.g[e ^ bj] ^ self.g[e ^ bi ^ bj] for e in range(1 << n))

    def J(self):
        """max j>=1 with d_0 d_j g != 0, or None if d_0 g is constant."""
        js = [j for j in range(1, self.n) if self.deriv2_nonzero(0, j)]
        return max(js) if js else None

    def self_dual(self):
        m = (1 << self.n) - 1
        return all(self.g[e ^ m] == 1 ^ self.g[e] for e in range(1 << self.n))

    def anf_str(self):
        ms = sorted(self.anf(), key=lambda m: (len(m), sorted(m)))
        if not ms: return "0"
        def mon(m):
            if not m: return "1"
            return "".join(f"y{k - self.r + 1:+d}".replace("+0", "0") for k in sorted(m))
        return " + ".join(mon(m) for m in ms)

    # ---------- cell dynamics
    def F_cell(self, cells):
        """cells: tuple of 2r+1 cell symbols (2x+y); returns new centre cell."""
        r = self.r
        x, y = cells[r] >> 1, cells[r] & 1
        ys = [c & 1 for c in cells]
        return (y << 1) | (x ^ self.f(ys))

    def Finv_cell(self, cells):
        r = self.r
        x, y = cells[r] >> 1, cells[r] & 1
        xs = [c >> 1 for c in cells]
        return ((y ^ self.f(xs)) << 1) | x

    def F_periodic(self, c):
        """F on a spatially periodic configuration given as a tuple of cells (period len(c))."""
        n = len(c); r = self.r
        return tuple(self.F_cell(tuple(c[(i + k) % n] for k in range(-r, r + 1))) for i in range(n))


def rule_from_anf(r, monomials, name=None):
    """monomials: iterable of tuples of xi-indices (k in 0..2r-1); () = constant 1."""
    n = 2 * r; tab = []
    for e in range(1 << n):
        xs = [(e >> (n - 1 - k)) & 1 for k in range(n)]
        v = 0
        for m in monomials:
            t = 1
            for k in m: t &= xs[k]
            v ^= t
        tab.append(v)
    return Rule(r, tab, name)

def elementary(rule_no):
    """elementary rule number -> Rule with r=1 if left-permutive (f = a ^ g(b,c))."""
    f = lambda a, b, c: (rule_no >> ((a << 2) | (b << 1) | c)) & 1
    assert all(f(0, b, c) ^ f(1, b, c) for b in (0, 1) for c in (0, 1)), "not left-permutive"
    g = [f(0, (e >> 1) & 1, e & 1) for e in range(4)]
    return Rule(1, g, f"ECA{rule_no}")

def random_rule(r, rng, nonlinear=True):
    n = 2 * r
    while True:
        R = Rule(r, [rng.randrange(2) for _ in range(1 << n)])
        if not nonlinear or not R.is_affine(): return R
