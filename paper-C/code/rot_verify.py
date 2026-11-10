# Second implementation of the checks for the rotating walls of 219R (Theorem 6.23, Proposition 6.24) and for
# the hidden-bit data of 126R (Proposition 6.25). Uses make_rule of centralizer.py.
import sys, itertools, random

from centralizer import make_rule

A = range(4)
pi = lambda s: ((s & 1) << 1) | (1 - (s >> 1))      # (x,y) -> (y, not x): 0->1->3->2->0
kap = lambda s: 3 - s

def Fring(Fl, c):
    n = len(c); return [Fl(c[i-1], c[i], c[(i+1) % n]) for i in range(n)]

# ---------------- 219R
f, Fl, Fil = make_rule(219)
print('219: f = 0 exactly on', [t for t in itertools.product((0,1), repeat=3) if f(*t) == 0])
# (1) bonded cells evolve by pi in every context (bond on the left or on the right)
bad = sum(1 for l, a, r in itertools.product(A, repeat=3) if (l == a or a == r) and Fl(l, a, r) != pi(a))
print('219 (1) bonded cell not mapped by pi, over all 64 triples:', bad)
# (2) bond invariance under F and F^-1 on all 4-windows
b2 = sum(1 for l, a, b, r in itertools.product(A, repeat=4) if (Fl(l, a, b) == Fl(a, b, r)) != (a == b))
b3 = sum(1 for l, a, b, r in itertools.product(A, repeat=4) if (Fil(l, a, b) == Fil(a, b, r)) != (a == b))
print('219 (2) bond not F-invariant:', b2, '| not F^-1-invariant:', b3)
# (3) car maps phi_s(a) = F(s a s), a != s
phi = {s: {a: Fl(s, a, s) for a in A if a != s} for s in A}
for s in A:
    img = set(phi[s].values())
    print('219 (3) phi_%d =' % s, phi[s], ' bijection onto A\\{pi(s)}:', img == set(A) - {pi(s)})
def Rret(s, a):
    for _ in range(4): a = phi[s][a]; s = pi(s)
    return a
for s in A:
    R = {a: Rret(s, a) for a in A if a != s}
    print('219 (4) R_%d =' % s, R, ' 3-cycle:', all(R[a] != a for a in R) and all(Rret(s, Rret(s, Rret(s, a))) == a for a in R))

# trains: flanked cars p with c[p-2]=c[p-1]=c[p+1]=c[p+2]=s != c[p]; train = maximal step-3 chain of same s
def flanked(c, p):
    n = len(c); s = c[(p+1) % n]
    return c[(p-2) % n] == c[(p-1) % n] == s == c[(p+2) % n] and c[p] != s
def trains(c):
    n = len(c); fl = {p for p in range(n) if flanked(c, p)}
    out = []
    for p in fl:
        s = c[(p+1) % n]
        if (p-3) % n in fl and c[(p-2) % n] == s: continue      # not the first member
        T = [p]
        while (T[-1]+3) % n in fl and c[(T[-1]+4) % n] == s and len(T) <= n: T.append((T[-1]+3) % n)
        if len(T) <= n // 3 - 1: out.append(T)
    return out
def Psi(c, m, perm):           # car r receives the old value of car perm[r]
    d = list(c)
    for T in trains(c):
        if len(T) == m:
            for r, q in enumerate(T): d[q] = c[T[perm[r]]]
    return d
def Dmap(c, k=1):
    d = list(c)
    for p in range(len(c)):
        if flanked(c, p):
            a = c[p]
            for _ in range(k): a = Rret(c[(p+1) % len(c)], a)
            d[p] = a
    return d
random.seed(219)
def plant(n):
    c = [random.randrange(4) for _ in range(n)]; p = 0
    while p < n - 20:
        if random.random() < 0.6:
            s = random.randrange(4); m = random.randint(1, 4)
            w = [s, s, s]
            for _ in range(m): w += [random.choice([v for v in A if v != s]), s, s]
            w += [s]
            for k, v in enumerate(w):
                if p + k < n: c[p+k] = v
            p += len(w) + random.randint(0, 3)
        else: p += random.randint(1, 4)
    return c
fails = dict(comm=0, hom=0, D=0, D3=0); nontriv = 0; noncomm = 0
for trial in range(1500):
    c = plant(120); Fc = Fring(Fl, c)
    for m in (1, 2, 3):
        for perm in itertools.permutations(range(m)):
            d = Psi(c, m, perm)
            nontriv += d != c
            fails['comm'] += Fring(Fl, d) != Psi(Fc, m, perm)
    a = Psi(Psi(c, 3, (1, 0, 2)), 3, (0, 2, 1)); b = Psi(c, 3, (1, 2, 0))   # (01) after (12) as maps on positions
    fails['hom'] += Psi(Psi(c, 3, (0, 2, 1)), 3, (1, 0, 2)) != Psi(c, 3, tuple((0, 2, 1)[(1, 0, 2)[r]] for r in range(3)))
    noncomm += a != Psi(Psi(c, 3, (0, 2, 1)), 3, (1, 0, 2))
    fails['D'] += Fring(Fl, Dmap(c)) != Dmap(Fc)
    fails['D3'] += Dmap(Dmap(Dmap(c))) != c
print('219 (5) planted rings: failures', fails, '| nontrivial Psi applications', nontriv, '| rings where (01),(12) do not commute', noncomm)
# (6) D on all windows: F o D = D o F as radius-3 local rules on all 4^7 windows
def Dloc(w):  # w of length 5, centre w[2]
    if w[0] == w[1] == w[3] == w[4] != w[2]: return Rret(w[1], w[2])
    return w[2]
bad = 0
for w in itertools.product(A, repeat=7):
    Fw = [Fl(*w[i:i+3]) for i in range(5)]              # F on cells 1..5
    lhs = Dloc(Fw[0:5])                                 # (D F)(c) at cell 3
    Dw = [Dloc(w[i:i+5]) for i in range(3)]             # D on cells 2..4
    if lhs != Fl(*Dw): bad += 1
print('219 (6) D F = F D on all 4^7 windows: failures', bad)

# ---------------- 126R: hidden-bit data for Theorem 6.6
f, Fl, Fil = make_rule(126)
W = [(b, kap(b), q, b, kap(b)) for b in A for q in A if q not in (b, kap(b))]
Ws = set(W)
img = lambda w, l, r: tuple(Fl(*((l,) + w + (r,))[i:i+3]) for i in range(5))
W1 = all(img(w, l, r) in Ws for w in W for l in A for r in A)
ctxfree = all(len({img(w, l, r) for l in A for r in A}) == 1 for w in W)
iota = lambda w: (w[0], w[1], kap(w[2]), w[3], w[4])
H1 = all(img(iota(w), l, r) == iota(img(w, l, r)) for w in W for l in A for r in A)
ov = sorted({d for d in range(1, 5) for a in W for b in W if all(a[d+k] == b[k] for k in range(5-d))})
cls = all(img(w, 0, 0)[0] == pi(w[0]) for w in W)
beta = lambda w: int(w[2] == pi(w[0]))
H3 = all(beta(img(w, 0, 0)) == beta(w) for w in W) and all(beta(iota(w)) == 1 - beta(w) for w in W)
bb = sum(1 for l, a, r in itertools.product(A, repeat=3) if (kap(l) == a or kap(a) == r) and Fl(l, a, r) != pi(a))
print("126: anti-bonded cells evolve by pi (failures %d); |W| = %d; (W1) %s; words context-free %s; (H1) %s; overlap offsets %s;"
      " class w_0 -> pi(w_0) %s; beta=[q=pi(b)] invariant, beta o iota = beta+1: %s" % (bb, len(W), W1, ctxfree, H1, ov, cls, H3))
