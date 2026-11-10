"""Direct (t,i) simulation of F_f(x,y)_i=(y_i, x_i ^ f(y_{i-1},y_i,y_{i+1})), using the local rule only.
Rows are Python ints (bit k <-> column i = iL+k).
A sequence (diagram of a point of Phi^+ / Phi^-) with constant tails is
   Seq = (lt, start, word, rt): value lt for x<start, word[x-start], rt after;  normalised."""
import itertools

def norm(lt, start, word, rt):
    word = list(word)
    while word and word[0] == lt:
        word.pop(0); start += 1
    while word and word[-1] == rt:
        word.pop()
    if not word and lt == rt:
        start = 0
    return (lt, start, tuple(word), rt)

def val(s, x):
    lt, st, w, rt = s
    if x < st: return lt
    if x < st + len(w): return w[x - st]
    return rt

def fromfun(fn, lo, hi):
    """sequence equal to fn(lo) left of lo, fn(hi) right of hi"""
    return norm(fn(lo), lo, [fn(x) for x in range(lo, hi + 1)], fn(hi))

def kinks(s):
    lt, st, w, rt = s
    return [x for x in range(st, st + len(w) + 1) if val(s, x - 1) != val(s, x)]

def is_bg(s): return not s[2] and s[0] == s[3]
def rev(s):   # mu: xi(p) -> xi(-1-p)
    lt, st, w, rt = s
    return norm(rt, -1 - (st + len(w) - 1), w[::-1], lt)
def shift(s, k): return norm(s[0], s[1] + k, s[2], s[3]) if not is_bg(s) else s
def canon(s):
    s = norm(*s); return (s[0], s[2], s[3])
def words(s, L):
    """all words of length L of the bi-infinite sequence"""
    lt, st, w, rt = s
    return {tuple(val(s, x + j) for j in range(L)) for x in range(st - L, st + len(w) + 1)}

def theta(a):    return (0, a, (), 1)
def K(a):        return theta(a)            # K^+_a / K^-_a  (diagram theta_a)
def Kb(a):       return (1, a, (), 0)       # bar K
def Z(a):        return (1, a, (0,), 1)     # isolated zero / bare defect
BG0 = (0, 0, (), 0); BG1 = (1, 0, (), 1)

def make_rule(num):
    tab = [(num >> k) & 1 for k in range(8)]
    def f(a, b, c, mask):
        na, nb, nc = a ^ mask, b ^ mask, c ^ mask
        r = 0
        for k in range(8):
            if tab[k]:
                r |= (a if k & 4 else na) & (b if k & 2 else nb) & (c if k & 1 else nc)
        return r
    return f
def f94(a, b, c, mask): return a ^ c ^ ((a ^ mask) & b & (c ^ mask))
def f95(a, b, c, mask): return a ^ c ^ ((a ^ mask) & (c ^ mask))
def f90(a, b, c, mask): return a ^ c
RULES = {94: f94, 95: f95, 90: f90}
def getrule(n):
    return RULES.get(n) or make_rule(n)
# sanity: fast forms agree with table forms
for _n in (94, 95, 90):
    _g = make_rule(_n)
    for _a, _b, _c in itertools.product((0, 1), repeat=3):
        assert RULES[_n](_a, _b, _c, 1) == _g(_a, _b, _c, 1) == (_n >> (4*_a + 2*_b + _c)) & 1

class NotSettled(Exception): pass

def scatter(P, X, rule=94, T=None, record=False, extra=0, corner=None):
    """Direct simulation of the configuration 'P (Phi^+, diagram pi(q), q=t-i) on the left,
    X (Phi^-, diagram xi(p), p=t+i) on the right' exactly as in the definition of a placed pair.
    Returns (Gamma, Omega[, info]).  Gamma: sequence in p, Omega: sequence in q."""
    f = getrule(rule)
    assert P[0] == X[0], "not placed: backgrounds differ"
    if corner is None:
        q0 = P[1] - 2; p0 = X[1] - 2
        if (q0 - p0) % 2: q0 -= 1
    else:
        q0, p0 = corner
        assert (q0 - p0) % 2 == 0 and all(val(P, q) == P[0] for q in (q0, q0+1)) and all(val(X, p) == X[0] for p in (p0, p0+1))
    t0 = (q0 + p0) // 2; i0 = (p0 - q0) // 2
    Pend = P[1] + len(P[2]); Xend = X[1] + len(X[2])
    bL, bR = P[3], X[3]
    if T is None: T = 2 * (len(P[2]) + len(X[2])) + 12
    while True:
        iL = t0 - Pend - T - 6; iR = Xend - t0 + T + 6
        W = iR - iL + 1; mask = (1 << W) - 1
        def row(t):
            r = 0
            for i in range(iL, iR + 1):
                v = val(P, t - i) if i < i0 else val(X, t + i)
                r |= v << (i - iL)
            return r
        prev, cur = row(t0 - 1), row(t0)
        rows = {t0 - 1: prev, t0: cur} if record else None
        fires = []
        hi = bR << (W - 1)
        for t in range(t0, t0 + T):
            a = ((cur << 1) | bL) & mask          # Y(t,i-1)
            c = (cur >> 1) | hi                   # Y(t,i+1)
            fv = f(a, cur, c, mask)
            n = fv ^ a ^ c                        # N-term: firing at centre (t,i)
            if n:
                k = 0; nn = n
                while nn:
                    low = (nn & -nn).bit_length() - 1
                    fires.append((t - (iL + low), t + (iL + low)))   # [q,p]
                    nn &= nn - 1
            prev, cur = cur, prev ^ fv
            if record: rows[t + 1] = cur
        t1 = t0 + T
        # edges must be untouched
        assert (cur & 7) == (7 if bL else 0) and (cur >> (W - 3)) == (7 if bR else 0)
        # D(i)=Y(t1,i)^Y(t1-1,i+1) (0 iff left-moving), E(i)=Y(t1,i)^Y(t1-1,i-1)
        D = cur ^ ((prev >> 1) | hi)
        E = cur ^ (((prev << 1) | bL) & mask)
        dmin = (D & -D).bit_length() - 1 if D else W
        emax = E.bit_length() - 1 if E else -1
        # need m with D=0 for k<=m+1, E=0 for k>=m-1, both rows constant on [m-2,m+2]
        m = None
        for k in range(max(emax + 2, 2), min(dmin - 1, W - 2)):
            seg = (cur >> (k - 2)) & 31; seg2 = (prev >> (k - 2)) & 31
            if seg in (0, 31) and seg2 == seg:
                m = k; break
        if m is not None: break
        T *= 2
        if T > 20000: raise NotSettled((P, X))
    bp = (cur >> m) & 1
    # gamma(p) = Y(t1, p-t1) for column k<=m ; o(q)=Y(t1,t1-q) for k>=m
    def g(p):
        k = p - t1 - iL
        if k < 0: return bL
        if k > m: return bp
        return (cur >> k) & 1
    def o(q):
        k = t1 - q - iL
        if k >= W: return bR
        if k < m: return bp
        return (cur >> k) & 1
    Gam = fromfun(g, t1 + iL, t1 + iL + m)
    Om = fromfun(o, t1 - iL - W + 1, t1 - iL - m)
    if record:
        def Y(q, p):
            assert (q - p) % 2 == 0
            t = (q + p) // 2; i = (p - q) // 2
            if t not in rows or not (iL <= i <= iR): return None
            return (rows[t] >> (i - iL)) & 1
        return Gam, Om, dict(fires=sorted(fires), Y=Y, q0=q0, p0=p0, t0=t0, t1=t1, iL=iL, iR=iR, T=T)
    return Gam, Om, sorted(fires)

def in_phi(s, sign, rule=94, T=12):
    """check by direct simulation that the configuration with rows s(t-i) (sign +) or s(t+i) (sign -)
    is a point of Phi^+ / Phi^-: i.e. the diagram stays a function of q (resp. p)."""
    f = getrule(rule)
    lo = s[1] - T - 6; hi_ = s[1] + len(s[2]) + T + 6
    W = hi_ - lo + 1; mask = (1 << W) - 1
    def row(t):
        r = 0
        for i in range(lo, hi_ + 1):
            r |= val(s, t - i if sign > 0 else t + i) << (i - lo)
        return r
    bl = s[3] if sign > 0 else s[0]; br = s[0] if sign > 0 else s[3]
    prev, cur = row(-1), row(0)
    for t in range(0, T):
        a = ((cur << 1) | bl) & mask; c = (cur >> 1) | (br << (W - 1))
        prev, cur = cur, prev ^ f(a, cur, c, mask)
        # compare interior (away from array edges, which the moving pattern does not respect)
        exp = row(t + 1)
        inner = mask >> (2 * (T + 3)) << (T + 3)
        if (cur ^ exp) & inner: return False
    return True

def objects(maxcore, rule=94, tails=(0, 1)):
    """all sequences with constant tails, core length <= maxcore, in Phi_f (no triple with N=1)"""
    f = make_rule(rule)
    def ok(s):
        lt, st, w, rt = s
        for x in range(st - 2, st + len(w) + 2):
            a, b, c = val(s, x - 1), val(s, x), val(s, x + 1)
            if f(a, b, c, 1) ^ a ^ c: return False
        return True
    out = []
    for lt in tails:
        for rt in tails:
            for n in range(0, maxcore + 1):
                for w in itertools.product((0, 1), repeat=n):
                    if n and (w[0] == lt or w[-1] == rt): continue
                    s = (lt, 0, w, rt)
                    if ok(s): out.append(s)
    return out
