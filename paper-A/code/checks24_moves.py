# (1) 66 moves (conditional tail flips) versus 46 moves computed from 46 probes, through delta.
# (2) Fix(sigma^n F^n), left-asymptotic to 00, for n=1..5 (and F189^4).
# (3) Light cone of the tail flip.
import random, itertools, sys
from eng24 import *
rng = random.Random(7)
LO, HI = -200, 300

def tailflip(b, be, r, s, R=66):
    u = b.Fn(R, s)
    u2 = mk(u.lo, u.hi, lambda i: u.cell(i) if i <= r else (1 - be, 1 - be))
    return u2.Fn(R, -s)

def M66(b, be, r, m, bar=False):
    x, y = b.Fn(66, m).cell(r)
    want = (1 - be, be) if bar else (be, 1 - be)
    return tailflip(b, be, r, m) if (x, y) == want else b

def debris46(c, r, n):
    """D_n(c) for rule 46, c 0-tailed with edge r, from the probe."""
    u = c.Fn(46, n)
    P = glue(u, eta(r + 2, u.lo, u.hi).rho(), r)   # state at time n of the probe
    s = 6
    v = P.Fn(46, s)  # state at time n+s
    # check it is D-part + eta
    w = mk(v.lo, v.hi, lambda i: v.cell(i) if i <= r else (0, 0))
    return w.Fn(46, -(n + s))

def M46(c, r, m, bar=False):
    if not bar:
        return debris46(c, r, m - 2)
    return debris46(c.rho(), r, -m - 2).rho()

def agree_common(a, b, lo=None, hi=None):
    lo_ = max(a.lo, b.lo) if lo is None else lo
    hi_ = min(a.hi, b.hi) if hi is None else hi
    return all(a.cell(i) == b.cell(i) for i in range(lo_, hi_))

ok = 0
for trial in range(400):
    be = rng.randint(0, 1); r = rng.randint(-2, 2)
    b = rand_tailed(be, r, LO, HI, rng)
    m = rng.randint(-6, 6); bar = rng.random() < 0.5
    b2 = M66(b, be, r, m, bar)
    db = delta(b)
    lhs = delta(b2)
    rhs = M46(db, r, m, bar)
    # compare on a central window
    assert agree_common(lhs, rhs, -60, 80), trial
    # light cone: b2 agrees with b on (-inf, r-|m|]
    assert agree_common(b, b2, -60, r - abs(m) + 1)
    ok += 1
print('66 moves = delta-lifts of 46 moves:', ok, 'cases')

# (2) Fix(sigma^n F^n) left-asymptotic to 00
def Fn_local(R, n, w):
    """y,x of (F^n c) at centre of window w (len 2n+1 cells)"""
    X = [c[0] for c in w]; Y = [c[1] for c in w]
    for _ in range(n):
        X, Y = Y[1:-1], [X[i] ^ fval(R, Y[i - 1], Y[i], Y[i + 1]) for i in range(1, len(X) - 1)]
    return (X[0], Y[0])

def fixsearch(R, n, Lw, margin):
    Z = (0, 0)
    res = []
    def ok_at(word, i):
        # constraint (F^n c)_i = c_{i-n}; window [i-n, i+n]
        w = [word[j] if j >= 0 else Z for j in range(i - n, i + n + 1)]
        tgt = word[i - n] if i - n >= 0 else Z
        return Fn_local(R, n, w) == tgt
    def dfs(word):
        L = len(word)
        i = L - 1 - n  # newly checkable constraint
        if i >= -n and not ok_at(word, i):
            return
        if L == Lw:
            res.append(tuple(word)); return
        for s in itertools.product((0, 1), repeat=2):
            dfs(word + [s])
    dfs([])
    return res

def Kword(p, Lw):
    return tuple((phi(p, 0, i + 1), phi(p, 0, i)) for i in range(Lw))

for R, n in [(66, 1), (66, 2), (66, 3), (66, 4), (66, 5), (189, 4)]:
    Lw = 14 + 2 * n; L0 = 14
    words = fixsearch(R, n, Lw, 0)
    pref = set(w[:L0] for w in words)
    expect = set([tuple([(0, 0)] * L0)]) | set(Kword(p, L0) for p in range(0, L0))
    assert pref == expect, (R, n)
    print('R', R, 'n', n, 'survivors', len(words), 'prefixes', len(pref), 'expected', len(expect), 'equal', pref == expect)
print('ALL CHECKS PASSED')
