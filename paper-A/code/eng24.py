# Engine for second-order ECA on finite windows.
# A configuration is (lo, X, Y): cells i in [lo, lo+len) with c_i=(X[i-lo],Y[i-lo]).
# Each step of F or F^{-1} shrinks the trusted window by one cell on each side.
import random
import sys, traceback

def _fail(t, v, tb):
    traceback.print_exception(t, v, tb)
    print('FAILURE')
    sys.exit(1)

sys.excepthook = _fail

def fval(R, a, b, c):
    return (R >> (4 * a + 2 * b + c)) & 1

class C:
    def __init__(self, lo, X, Y):
        self.lo, self.X, self.Y = lo, list(X), list(Y)
    @property
    def hi(self):
        return self.lo + len(self.X)  # exclusive
    def cell(self, i):
        return (self.X[i - self.lo], self.Y[i - self.lo])
    def F(self, R):
        X, Y = self.X, self.Y
        nX = Y[1:-1]
        nY = [X[i] ^ fval(R, Y[i - 1], Y[i], Y[i + 1]) for i in range(1, len(X) - 1)]
        return C(self.lo + 1, nX, nY)
    def Finv(self, R):
        X, Y = self.X, self.Y
        nY = X[1:-1]
        nX = [Y[i] ^ fval(R, X[i - 1], X[i], X[i + 1]) for i in range(1, len(X) - 1)]
        return C(self.lo + 1, nX, nY)
    def Fn(self, R, n):
        c = self
        for _ in range(abs(n)):
            c = c.F(R) if n > 0 else c.Finv(R)
        return c
    def rho(self):
        return C(self.lo, self.Y, self.X)
    def kappa(self):
        return C(self.lo, [1 - v for v in self.X], [1 - v for v in self.Y])
    def agree(self, other, a, b):
        """agree on [a,b] (inclusive); both must cover it"""
        assert self.lo <= a and b < self.hi and other.lo <= a and b < other.hi, (self.lo, self.hi, other.lo, other.hi, a, b)
        return all(self.cell(i) == other.cell(i) for i in range(a, b + 1))

def mk(lo, hi, fn):
    cells = [fn(i) for i in range(lo, hi)]
    return C(lo, [c[0] for c in cells], [c[1] for c in cells])

def phi(p, be, i):
    return be if i <= p else be ^ ((i - p) & 1)

def K(p, be, lo, hi):
    return mk(lo, hi, lambda i: (phi(p, be, i + 1), phi(p, be, i)))

def L(p, be, lo, hi):  # L_{p,be} = rho K_{p,be}
    return mk(lo, hi, lambda i: (phi(p, be, i), phi(p, be, i + 1)))

def const(s, lo, hi):
    return mk(lo, hi, lambda i: s)

def eta(p, lo, hi):
    return mk(lo, hi, lambda i: (0, 0) if i < p else ((1, 0) if i == p else (1, 1)))

def glue(d, lam, r):
    assert d.lo == lam.lo and d.hi == lam.hi
    return mk(d.lo, d.hi, lambda i: d.cell(i) if i <= r else lam.cell(i))

def delta(c):
    # (delta c)_i = c_i xor c_{i+1}; window loses the last cell
    return mk(c.lo, c.hi - 1, lambda i: (c.cell(i)[0] ^ c.cell(i + 1)[0], c.cell(i)[1] ^ c.cell(i + 1)[1]))

def rand_tailed(be, e, lo, hi, rng):
    def g(i):
        if i > e:
            return (be, be)
        if i == e:
            while True:
                s = (rng.randint(0, 1), rng.randint(0, 1))
                if s != (be, be):
                    return s
        return (rng.randint(0, 1), rng.randint(0, 1))
    return mk(lo, hi, g)
