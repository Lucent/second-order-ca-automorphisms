"""Rule 95R: the map F, block maps of the golden mean shift G, and the maps hat H_{phi,psi} sigma^n.
   F(x,y)_i = (y_i, x_i XOR f(y_{i-1},y_i,y_{i+1})), f = Wolfram rule 95.
   Finite perturbations of the background (1,1) are stored as (lo, X, Y): cells lo..lo+n-1, (1,1) outside.
   hat H_{phi,psi} sigma^n is computed from its characterisation: it commutes with F, and on a configuration that is a
   point of Phi^+ on the left and a point of Phi^- on the right (far apart) it applies phi to the two classes of the
   right-movers and psi to the two classes of the left-movers."""
import numpy as np, itertools

RULE = 95
def f95(a, b, c): return (RULE >> (4 * a + 2 * b + c)) & 1
FT = np.array([[[f95(a, b, c) for c in (0, 1)] for b in (0, 1)] for a in (0, 1)], dtype=np.uint8)

def Fstep(X, Y):
    """one step on arrays, background 1 beyond both ends; arrays must have >=1 background cell of margin (asserted)."""
    assert X[0] == 1 and Y[0] == 1 and X[-1] == 1 and Y[-1] == 1
    Yl = np.concatenate(([1], Y[:-1])); Yr = np.concatenate((Y[1:], [1]))
    return Y.copy(), X ^ FT[Yl, Y, Yr]

def Finv(X, Y):
    assert X[0] == 1 and Y[0] == 1 and X[-1] == 1 and Y[-1] == 1
    Xl = np.concatenate(([1], X[:-1])); Xr = np.concatenate((X[1:], [1]))
    return Y ^ FT[Xl, X, Xr], X.copy()

def pad(X, Y, m):
    o = np.ones(m, dtype=np.uint8)
    return np.concatenate((o, X, o)), np.concatenate((o, Y, o))

# ---------- periodic configurations (tuples of (x,y)) ----------
def Fper(c):
    n = len(c)
    return tuple((c[i][1], c[i][0] ^ f95(c[i - 1][1], c[i][1], c[(i + 1) % n][1])) for i in range(n))
def Sper(c, a=1):
    n = len(c); return tuple(c[(i + a) % n] for i in range(n))

# ---------- block maps of the golden mean shift ----------
class GMap:
    """(phi x)_i = tab[word x_{i+a} ... x_{i+a+w-1}] (word read as a binary number, first letter most significant)."""
    def __init__(self, w, a, tab):
        self.w, self.a, self.tab = w, a, tab            # tab: dict word-int -> 0/1 on words of G; missing -> error
    def apply_set(self, S):
        """S: finite set of positions of 1s (no two adjacent). returns the set of 1s of the image."""
        if not S: return set()
        lo, hi = min(S) - self.w - 1, max(S) + self.w + 1
        x = np.zeros(hi - lo + 1, dtype=np.int64)
        for s in S: x[s - lo] = 1
        n = len(x) - self.w + 1
        idx = np.zeros(n, dtype=np.int64)
        for k in range(self.w): idx = (idx << 1) | x[k:k + n]
        out = set()
        for m in np.unique(idx):
            if self.tab[int(m)]:
                for p in np.nonzero(idx == m)[0]: out.add(int(p) + lo - self.a)   # window starts at i+a = p+lo
        return out

def gm_shift(k):
    """sigma^k as a GMap: (sigma^k x)_i = x_{i+k}"""
    return GMap(1, k, {0: 0, 1: 1})
def gm_swap(u, v):
    """involution exchanging the words u and v (equal length, as strings) wherever they occur"""
    L = len(u); d = [i for i in range(L) if u[i] != v[i]]
    lo, hi = min(d), max(d)
    # output at position i: is i at index k (lo<=k<=hi) of an occurrence of u or v?  window x_{i-hi .. i+L-1-lo}
    w = hi + (L - 1 - lo) + 1; a = -hi
    tab = {}
    for m in range(1 << w):
        s = format(m, '0%db' % w)
        if '11' in s: continue
        val = int(s[hi])
        for k in range(lo, hi + 1):
            st = hi - k
            if s[st:st + L] == u: val = int(v[k])
            elif s[st:st + L] == v: val = int(u[k])
        tab[m] = val
    return GMap(w, a, tab)
ID = gm_shift(0)

# ---------- hat H_{phi,psi} sigma^n on finite perturbations ----------
def separate(lo, X, Y, margin, tmax=10**6, step=None):
    """evolve backwards until the configuration is [Phi^+ part][>=margin cells (1,1)][Phi^- part]; returns t, lo, X, Y, (u,v)"""
    n = len(X); t = 0
    step = step or max(8, n // 8)
    while True:
        bp = np.nonzero(X[:-1] != Y[1:])[0]            # cells i (index) violating x_i = y_{i+1}
        bm = np.nonzero(X[1:] != Y[:-1])[0] + 1        # cells violating x_i = y_{i-1}
        a = (bm.max() + 1) if len(bm) else 0           # allowed split region: indices in [a, b)
        b = bp.min() if len(bp) else len(X)
        if b - a >= margin:
            bg = (X[a:b] & Y[a:b]).astype(np.int8)
            # longest run of ones in bg
            d = np.diff(np.concatenate(([0], bg, [0])))
            st = np.nonzero(d == 1)[0]; en = np.nonzero(d == -1)[0]
            if len(st):
                k = np.argmax(en - st)
                if en[k] - st[k] >= margin: return t, lo, X, Y, (a + st[k], a + en[k] - 1)
        assert t < tmax
        X, Y = pad(X, Y, step); lo -= step
        for _ in range(step): X, Y = Finv(X, Y)
        t += step

def trim(lo, X, Y, keep=2):
    nz = np.nonzero((X & Y) == 0)[0]
    if len(nz) == 0: return 0, np.ones(2 * keep, np.uint8), np.ones(2 * keep, np.uint8)
    a, b = nz.min(), nz.max()
    o = np.ones(keep, dtype=np.uint8)
    return lo + a - keep, np.concatenate((o, X[a:b + 1], o)), np.concatenate((o, Y[a:b + 1], o))

def hatH(lo, X, Y, phi, psi, nsig=0, margin=64, extra=0):
    """returns (lo', X', Y') = hat H_{phi,psi} sigma^nsig (c).  sigma^n c has cells c_{i+n}: lo -> lo - n."""
    lo -= nsig
    lo, X, Y = trim(lo, X, Y, 2)
    t, lo, X, Y, (u, v) = separate(lo, X, Y, margin)
    if extra:
        X, Y = pad(X, Y, extra); lo -= extra; u += extra; v += extra
        for _ in range(extra): X, Y = Finv(X, Y)
        t += extra
        # after `extra` more steps the gap has grown by 2*extra: [u-extra, v+extra]
        u -= extra; v += extra
        assert (X[u:v + 1] & Y[u:v + 1]).all()
    s = (u + v) // 2                                    # split index; cells <= s belong to P, > s to X
    # checks: left part is a point of Phi^+ (x_i = y_{i+1}; pi(q) = y_{-q} has no two zeros at distance 2), same for right
    assert (X[:s] == Y[1:s + 1]).all() and (X[s + 1:] == Y[s:-1]).all()
    zl = np.nonzero(Y[:s + 1] == 0)[0]; zr = np.nonzero(Y[s:] == 0)[0] + s
    assert 2 not in set(np.diff(zl)) and 2 not in set(np.diff(zr))
    # particles.  absolute cell i = lo + index.  right-movers: pi(q)=y_{-q}: zero at cell i <-> q=-i, class e=q mod 2, j=(q-e)/2
    P = [set(), set()]; Q = [set(), set()]
    for k in zl:
        q = -(lo + int(k)); P[q & 1].add((q - (q & 1)) // 2)
    for k in zr:
        p = lo + int(k); Q[p & 1].add((p - (p & 1)) // 2)
    P2 = [phi.apply_set(P[0]), phi.apply_set(P[1])]; Q2 = [psi.apply_set(Q[0]), psi.apply_set(Q[1])]
    zq = sorted(2 * j + e for e in (0, 1) for j in P2[e]); zp = sorted(2 * k + e for e in (0, 1) for k in Q2[e])
    cs = lo + s                                         # absolute split cell
    # new configuration: cells i <= cs: (pi'(-1-i), pi'(-i)); cells i > cs: (xi'(i-1), xi'(i))
    # require the images to stay on their sides with room
    assert all(-q <= cs - 3 for q in zq) and all(p >= cs + 3 for p in zp), "margin too small"
    lo2 = min([-q for q in zq] + [cs]) - 4; hi2 = max(zp + [cs]) + 4
    n2 = hi2 - lo2 + 1
    X2 = np.ones(n2, np.uint8); Y2 = np.ones(n2, np.uint8)
    for q in zq:
        Y2[-q - lo2] = 0; X2[-q - 1 - lo2] = 0          # y_i = pi(-i): i=-q ; x_i = pi(-1-i): i=-1-q
    for p in zp:
        Y2[p - lo2] = 0; X2[p + 1 - lo2] = 0            # y_i = xi(i): i=p ; x_i = xi(i-1): i=p+1
    X2, Y2 = pad(X2, Y2, t + 2); lo2 -= t + 2
    for _ in range(t): X2, Y2 = Fstep(X2, Y2)
    return trim(lo2, X2, Y2, 2)

def to_dict(lo, X, Y):
    return {lo + i: (int(X[i]), int(Y[i])) for i in range(len(X)) if not (X[i] and Y[i])}
def from_cells(lo, cells):
    X = np.array([c[0] for c in cells], np.uint8); Y = np.array([c[1] for c in cells], np.uint8)
    X, Y = pad(X, Y, 2); return lo - 2, X, Y

# ---------- self-test ----------
if __name__ == "__main__":
    # hat H_{phi,psi} for (phi,psi) = (sigma,id), (id,sigma), (sigma^-1,sigma) is sigma^-1 F, sigma F, sigma^2,
    # on pseudo-random finite perturbations of the background
    import random
    random.seed(1)
    def shifted(c, n): lo, X, Y = c; return trim(lo - n, X, Y, 2)        # sigma^n
    def Fc(c):
        lo, X, Y = c; X, Y = pad(X, Y, 2); X, Y = Fstep(X, Y); return trim(lo - 2, X, Y, 2)
    cases = [("(sigma, id) = sigma^-1 F", gm_shift(1), ID, lambda c: shifted(Fc(c), -1)),
             ("(id, sigma) = sigma F", ID, gm_shift(1), lambda c: shifted(Fc(c), 1)),
             ("(sigma^-1, sigma) = sigma^2", gm_shift(-1), gm_shift(1), lambda c: shifted(c, 2))]
    for name, phi, psi, ref in cases:
        bad = 0; n = 300
        for _ in range(n):
            w = random.randint(1, 12)
            c = from_cells(0, [(random.randint(0, 1), random.randint(0, 1)) for _ in range(w)])
            if to_dict(*hatH(*c, phi, psi)) != to_dict(*ref(c)): bad += 1
        print(f"hat H for (phi, psi) = {name}: configurations {n} mismatches {bad}")
