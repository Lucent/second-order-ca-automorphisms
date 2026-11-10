"""Rule 91, background 11, exact finite checks (Goursat solver scat.solve; second implementation: direct cell
simulation in t91c.py).
(1) Objects X in Phi^- (cores = words of Phi_91 starting and ending with 0, length <= L) that are transparent to the
    bare defect as a probe in both parities, S(Delta^+,X_b) = (X_b,Delta^+), b = 0,1.  Claim: exactly the 3-sparse ones.
(2) First-ray lemma and second-ray lemma on random data."""
import sys, itertools, random
from scat import solve, phi_words
from t91a import transparent
RULE = 91
def sparse3(c):
    z = [i for i, v in enumerate(c) if v == 0]
    return all(z[i+1] - z[i] >= 3 for i in range(len(z)-1))
L = int(sys.argv[1]) if len(sys.argv) > 1 else 11
tot = 0; tr = 0; bad = []
for n in range(1, L+1):
    for w in phi_words(RULE, n):
        if w[0] or w[-1]: continue
        tot += 1
        t = transparent((0,), w, 0) and transparent((0,), w, 1)
        tr += t
        if t != sparse3(w): bad.append(w)
print("(1) cores of Phi_91 of length <=", L, ":", tot, " transparent to the bare defect:", tr, " disagreements with 3-sparse:", len(bad), bad[:5])
# (2)
random.seed(4)
def rand_word(n, isolated):
    w = [1, 1]
    while len(w) < n:
        b = random.random() < 0.5
        if b == 0 and (w[-1] == 0 and (isolated or w[-2] == 0)): b = 1
        w.append(int(b))
    return w + [1]*8
viol1 = viol2 = 0; n1 = n2 = 0
for _ in range(4000):
    M = 44
    iso = random.random() < 0.5
    B = rand_word(M-8, iso)[:M]
    a1 = random.randint(2, 7)
    A = [1]*M; A[a1] = 0
    for q in range(a1+1+iso, M-8): A[q] = random.randint(0, 1)   # arbitrary continuation (need not be in Phi)
    Y = solve(RULE, A, B)
    N = lambda a, b, c: int(a == 0 and b == 0 and c == 0)
    def fires(q, p): return N(Y[q+1][p-1], Y[q][p], Y[q-1][p+1])
    # ray a1-1: centres [a1-1,p'], p' = a1-1 mod 2, p' >= 1
    f1 = [p for p in range(1, M-1) if (p - (a1-1)) % 2 == 0 and fires(a1-1, p)]
    xs = [x for x in range(3, M) if B[x] == 0 and B[x-1] == 0 and (x - a1) % 2 == 0]
    exp = [xs[0]-1] if xs else []
    n1 += 1; viol1 += (f1 != exp)
    if iso:
        f2 = [p for p in range(1, M-1) if (p - a1) % 2 == 0 and fires(a1, p)]
        zs = [x for x in range(2, M) if B[x] == 0 and (x - (a1+1)) % 2 == 0]
        par = (B[zs[-1]-2] == 0) if zs else False
        n2 += 1; viol2 += ((len(f2) % 2 == 1) != par)
print("(2) first-ray lemma: cases", n1, "violations", viol1, "; second-ray parity lemma: cases", n2, "violations", viol2)
