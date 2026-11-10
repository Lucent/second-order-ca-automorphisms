"""Second implementation of check (1) of t91b.py: direct simulation of F_91 on cells.
The configuration is Delta^+ (one zero, right-moving) on the left and an object X of Phi^- on the right, background
(1,1).  Transparency  <=>  for large T the configuration F^T(c) equals the free superposition
Y0(t,i) = pi(t-i) xor xi(t+i) xor 1 on the rows T-1, T."""
import sys
from scat import phi_words
RULE = 91
f = lambda a, b, c: (RULE >> (4*a + 2*b + c)) & 1
def step(X, Y):
    n = len(X)
    return Y, [X[i] ^ f(Y[i-1] if i else 1, Y[i], Y[i+1] if i < n-1 else 1) for i in range(n)]
def transparent(core, off):
    L = len(core); T = 2*L + 24; n = 4*T + L + 40
    a = T + 8                      # pi(q) = 0 iff q = -a   (so at t=0 the zero is at i = a ... moving right)
    s = a + 6 + off                # xi(p) = core[p-s] for s <= p < s+L
    pi = lambda q: 0 if q == -a else 1
    xi = lambda p: core[p-s] if s <= p < s+L else 1
    Y0 = lambda t, i: pi(t-i) ^ xi(t+i) ^ 1
    X = [Y0(-1, i) for i in range(n)]; Y = [Y0(0, i) for i in range(n)]
    for _ in range(T): X, Y = step(X, Y)
    return all(X[i] == Y0(T-1, i) and Y[i] == Y0(T, i) for i in range(n))
def sparse3(c):
    z = [i for i, v in enumerate(c) if v == 0]
    return all(z[i+1] - z[i] >= 3 for i in range(len(z)-1))
L = int(sys.argv[1]) if len(sys.argv) > 1 else 11
tot = tr = bad = 0
for n in range(1, L+1):
    for w in phi_words(RULE, n):
        if w[0] or w[-1]: continue
        tot += 1
        t = transparent(w, 0) and transparent(w, 1)
        tr += t; bad += (t != sparse3(w))
print("cores of length <=", L, ":", tot, "transparent:", tr, "disagreements with 3-sparse:", bad)
