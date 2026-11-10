"""Rules 219R and 126R: which theorem of Section 6 applies.
(a) Corollary 6.9: is there a cycle of 5-words w_s with w_{s+h} differing from w_s only at position 2?
(b) Theorem 6.6 for 126R with W = {(b, kb, q, b, kb) : q not in {b,kb}}, j = 2, iota = kappa on the car:
    (W1), (H1) on all contexts; overlaps (H2), (H4); (H3') with K the phase of the wall, eps = 0.
"""
from centralizer import make_rule
import itertools
pi = lambda s: ((s & 1) << 1) | (1 - (s >> 1)); kap = lambda s: 3 - s
def Fword(Fl, l, w, r):            # image of word w in context l,r
    ext = (l,) + tuple(w) + (r,)
    return tuple(Fl(*ext[i:i + 3]) for i in range(len(w)))
# (a)
for rule, words in ((219, [(s, s, a, s, s) for s in range(4) for a in range(4) if a != s]),
                    (126, [(b, kap(b), q, b, kap(b)) for b in range(4) for q in range(4) if q not in (b, kap(b))])):
    f, Fl, _ = make_rule(rule)
    ctxfree = all(len({Fword(Fl, l, w, r) for l in range(4) for r in range(4)}) == 1 for w in words)
    orbit_ok = all(Fword(Fl, 0, w, 0) in words for w in words)
    # orbit lengths and whether F^h w differs from w only at position 2
    res = []
    for w in words:
        orb = [w]
        while True:
            nw = Fword(Fl, 0, orb[-1], 0)
            if nw == w: break
            orb.append(nw)
        P = len(orb); h = P // 2
        only2 = P % 2 == 0 and [k for k in range(5) if orb[h][k] != w[k]] == [2]
        res.append((P, only2))
    print(rule, 'car words context-free:', ctxfree, 'family F-invariant:', orbit_ok,
          'orbit periods:', sorted({p for p, _ in res}), 'Corollary 6.9 (b) holds for some word:', any(o for _, o in res))
# (b) 126R hidden-bit
f, Fl, _ = make_rule(126)
W = [(b, kap(b), q, b, kap(b)) for b in range(4) for q in range(4) if q not in (b, kap(b))]
Wset = set(W)
iota = lambda w: (w[0], w[1], kap(w[2]), w[3], w[4])
W1 = all(Fword(Fl, l, w, r) in Wset for w in W for l in range(4) for r in range(4))
# (H1): F(iota_0 c) = iota_0 F(c): only cells 1..3 can differ; check on 7-windows (l, w, r) + note cells outside unaffected
H1 = all(Fword(Fl, l, iota(w), r) == iota(Fword(Fl, l, w, r)) for w in W for l in range(4) for r in range(4))
# overlaps: offsets delta with occurrences at 0 and delta in some word of length 5+delta
ov = set()
for dlt in range(1, 5):
    for a in W:
        for b in W:
            if all(a[dlt + k] == b[k] for k in range(5 - dlt)): ov.add(dlt)
# (H3'): class = wall symbol b = w[0]; theta = pi; beta(w) = [car == canon(b)], canon chosen along the pi-orbit
canon = {0: 1}
for t in range(3):
    b = [0, 1, 3, 2][t]; canon[pi(b)] = pi(canon[b])
assert pi(canon[2]) == canon[0]
beta = lambda w: int(w[2] == canon[w[0]])
H3 = True; epsv = {}
for w in W:
    Fw = Fword(Fl, 0, w, 0)
    if Fw[0] != pi(w[0]): H3 = False
    epsv.setdefault(w[0], beta(Fw) ^ beta(w))
    if epsv[w[0]] != (beta(Fw) ^ beta(w)): H3 = False
print('126R hidden-bit: (W1)', W1, '(H1)', H1, 'overlap offsets', sorted(ov), '(max(j,k-1-j)=2, so (H2),(H4) hold if 1,2 absent)',
      "(H3') with K=Z/4, theta=pi", H3, 'eps', epsv, 'beta o iota = beta+1:', all(beta(iota(w)) == 1 - beta(w) for w in W))
