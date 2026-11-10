"""Properties of the radius-7 table map h = h_* (tx7.npy, ty7.npy):
   exact dependence sets; restriction to Phi^+ and Phi^- (all 3025 words of length 16 of Phi_95);
   explicit periodic witness for h not in <sigma,F>; direct tests on all periodic configurations of period <= 9.
   usage: props7.py [tx7.npy ty7.npy]"""
import numpy as np, itertools, os, sys
from core import f95, Fper, Sper, gm_swap
HERE = os.path.dirname(os.path.abspath(__file__))
TX = np.load(sys.argv[1] if len(sys.argv) > 2 else os.path.join(HERE, "tx7.npy"))
TY = np.load(sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "ty7.npy")); R = 7
idx = np.arange(1 << 15)
print("exact dependence (offset m: number of windows in which flipping the bit at m changes the output):")
for name, T, lay in (("y'", TY, lambda m: "y" if m % 2 else "x"), ("x'", TX, lambda m: "x" if m % 2 else "y")):
    dep = {m: int((T[idx] != T[idx ^ (1 << (m + R))]).sum()) // 2 for m in range(-R, R + 1)}
    print("  %s:" % name, {"%s%+d" % (lay(m), m): v for m, v in dep.items()})

def g_per(c):
    n = len(c); out = []
    for i in range(n):
        wy = wx = 0
        for m in range(-R, R + 1):
            x, y = c[(i + m) % n]
            wy |= (y if m % 2 else x) << (m + R); wx |= (x if m % 2 else y) << (m + R)
        out.append((int(TX[wx]), int(TY[wy])))
    return tuple(out)

# restrictions
def gwords(n): return [m for m in range(1 << n) if not (m & (m >> 1))]
phi1 = gm_swap('0001001', '0000101')
def phi1_at0(xw):     # xw: dict j -> bit, j in -4..3 ; (phi1 x)_0 = tab[x_{-4..3}]
    e = 0
    for j in range(-4, 4): e = (e << 1) | xw[j]
    return phi1.tab[e]
badp = badm = n = 0
G8 = gwords(8)
for a in G8:
    for b in G8:
        # pi on [-8,7]: class 0 at q=-8,-6,...,6 (x^0_j, j=-4..3), class 1 at q=-7,...,7 (x^1_j, j=-4..3); pi = 1 - x
        x0 = {j: (a >> (3 - j)) & 1 for j in range(-4, 4)}; x1 = {j: (b >> (3 - j)) & 1 for j in range(-4, 4)}
        pi = {}
        for j in range(-4, 4): pi[2 * j] = 1 - x0[j]; pi[2 * j + 1] = 1 - x1[j]
        n += 1
        # point of Phi^+: c_i = (pi(-1-i), pi(-i)), cells -7..7
        c = {i: (pi[-1 - i], pi[-i]) for i in range(-7, 8)}
        wy = sum((c[m][1] if m % 2 else c[m][0]) << (m + R) for m in range(-R, R + 1))
        wx = sum((c[m][0] if m % 2 else c[m][1]) << (m + R) for m in range(-R, R + 1))
        got = (int(TX[wx]), int(TY[wy]))
        # expected: hatH c has diagram pi'(2j+e) = 1 - (phi_* x^e)_j ; h(c)_0 = hatH(c)_{-1} = (pi'(0), pi'(1))
        exp = (1 - phi1_at0(x0), 1 - phi1_at0(x1))
        badp += got != exp
        # point of Phi^-: c_i = (xi(i-1), xi(i)), xi on [-8,7] ; h = sigma^{-1}: h(c)_0 = c_{-1} = (xi(-2), xi(-1))
        xi = pi
        c = {i: (xi[i - 1], xi[i]) for i in range(-7, 8)}
        wy = sum((c[m][1] if m % 2 else c[m][0]) << (m + R) for m in range(-R, R + 1))
        wx = sum((c[m][0] if m % 2 else c[m][1]) << (m + R) for m in range(-R, R + 1))
        badm += (int(TX[wx]), int(TY[wy])) != (xi[-2], xi[-1])
print("restriction to Phi^+ = sigma^-1 o (phi_* on each class): %d words, %d failures; restriction to Phi^- = sigma^-1: %d failures" % (n, badp, badm))

# periodic configurations
def orbit_sF(c):
    S = set(); d = c
    while True:
        for a in range(len(c)): S.add(Sper(d, a))
        d = Fper(d)
        if d == c: break
    return S
first = None; tot = {}
for n in range(1, 10):
    moved = bad = out = 0
    for cells in itertools.product(range(4), repeat=n):
        c = tuple((s >> 1, s & 1) for s in cells)
        gc = g_per(c)
        if g_per(Fper(c)) != Fper(gc): bad += 1
        sg = Sper(gc, 1)
        if Sper(g_per(sg), 1) != c: bad += 1
        if sg != c:
            moved += 1
            if n <= 7 and gc not in orbit_sF(c):
                out += 1
                if first is None: first = c
    print("period %d: %d configurations, failures of hF=Fh or (sigma h)^2=id: %d; sigma h moves %d; (n<=7) h(c) outside the <sigma,F>-orbit of c: %d" % (n, 4 ** n, bad, moved, out), flush=True)
print("witness c =", first, " h(c) =", g_per(first), " F-period", next(k for k in range(1, 10**4) if (lambda d: d == first)(__import__('functools').reduce(lambda d, _: Fper(d), range(k), first))))
