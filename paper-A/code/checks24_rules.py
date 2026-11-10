# Rule identities, delta factor, kinks, gluing and the bounce lemma for F_66.
import itertools, random
from eng24 import *

B3 = list(itertools.product((0, 1), repeat=3))
f66 = lambda a, b, c: fval(66, a, b, c)
f46 = lambda a, b, c: fval(46, a, b, c)
# formulas
assert all(f66(a, b, c) == ((a == b) & (b ^ c)) for a, b, c in B3)
assert all(fval(24, a, b, c) == f66(c, b, a) for a, b, c in B3)
assert all(fval(189, a, b, c) == 1 ^ f66(a, b, c) for a, b, c in B3)
assert all(fval(231, a, b, c) == 1 ^ fval(24, a, b, c) for a, b, c in B3)
assert all(f66(a, b, c) == f66(1 - a, 1 - b, 1 - c) for a, b, c in B3)
# allowed triples of 66 and 189
print('66 forbidden', [t for t in B3 if (t[0] ^ t[2]) != f66(*t)])
print('189 allowed', [t for t in B3 if (t[0] ^ t[2]) == 1 ^ f66(*t)])
# delta identity
assert all((f66(a, b, c) ^ f66(b, c, d)) == f46(a ^ b, b ^ c, c ^ d) for a, b, c, d in itertools.product((0, 1), repeat=4))
# f46(a,0,c)=c, f46(a,1,c)=not a
assert all(f46(a, 0, c) == c and f46(a, 1, c) == 1 - a for a, c in itertools.product((0, 1), repeat=2))
print('rule identities ok')

rng = random.Random(1)
LO, HI = -150, 250
# operator identities on random windows
for _ in range(200):
    c = mk(LO, HI, lambda i: (rng.randint(0, 1), rng.randint(0, 1)))
    # delta F66 = F46 delta, delta F189 = F46 delta
    for R in (66, 189):
        a = delta(c.F(R)); b = delta(c).F(46)
        assert a.agree(b, b.lo, b.hi - 1)
    # F189 = nu F66 with nu(x,y)=(x,not y); F189^2 = kappa F66^2; F189^4 = F66^4
    a, b = c.F(189), c.F(66)
    assert all(a.cell(i) == (b.cell(i)[0], 1 - b.cell(i)[1]) for i in range(a.lo, a.hi))
    a, b = c.Fn(189, 2), c.Fn(66, 2).kappa()
    assert a.agree(b, a.lo, a.hi - 1)
    a, b = c.Fn(189, 4), c.Fn(66, 4)
    assert a.agree(b, a.lo, a.hi - 1)
    # kappa commutes with F66, F189
    for R in (66, 189):
        assert c.kappa().F(R).agree(c.F(R).kappa(), c.lo + 1, c.hi - 2)
print('operator identities ok')

# kinks
for be in (0, 1):
    for p in range(-5, 6):
        k = K(p, be, LO, HI); l = L(p, be, LO, HI)
        assert k.F(66).agree(K(p + 1, be, LO, HI), LO + 1, HI - 2)
        assert l.F(66).agree(L(p - 1, be, LO, HI), LO + 1, HI - 2)
        assert k.kappa().agree(K(p, 1 - be, LO, HI), LO, HI - 1)
        assert delta(k).agree(eta(p - 1, LO, HI), LO, HI - 2)
        assert delta(l).agree(eta(p - 1, LO, HI).rho(), LO, HI - 2)
        # right tail of K_{p,0} is sigma^{-p}A, of L_{p,0} is sigma^{1-p}A ; A_i=((i+1)&1,i&1)
        if be == 0:
            assert all(k.cell(i) == (((i - p) + 1) & 1, (i - p) & 1) for i in range(p + 1, HI))
            assert all(l.cell(i) == (((i + 1 - p) + 1) & 1, (i + 1 - p) & 1) for i in range(p + 1, HI))
print('kinks ok')

def tailflip(b, be, r, s, R=66):
    u = b.Fn(R, s)
    assert u.cell(r) != (1 - be, 1 - be)
    u2 = mk(u.lo, u.hi, lambda i: u.cell(i) if i <= r else (1 - be, 1 - be))
    return u2.Fn(R, -s)

# gluing facts (G1),(G2) on arbitrary left parts
for _ in range(300):
    be = rng.randint(0, 1); q = rng.randint(-5, 5)
    junk = mk(LO, HI, lambda i: (rng.randint(0, 1), rng.randint(0, 1)))
    p = q + rng.randint(1, 6)
    u = glue(junk, K(p, be, LO, HI), q - 1)
    assert u.F(66).agree(K(p + 1, be, LO, HI), q, HI - 2)
    p = q + rng.randint(2, 6)
    u = glue(junk, L(p, be, LO, HI), q - 1)
    assert u.F(66).agree(L(p - 1, be, LO, HI), q, HI - 2)
    # sharpness: p=q (K) or p=q+1 (L) can fail
print('gluing ok')

# bounce lemma
cnt = {}
for trial in range(500):
    be = rng.randint(0, 1); e = rng.randint(-3, 3)
    d = rand_tailed(be, e, LO, HI, rng)
    p = e + 2 + rng.randint(0, 14)
    c = glue(d, L(p, be, LO, HI), e)
    s0 = p - e - 2
    eps = int(d.Fn(66, s0).cell(e)[1] == be)
    U = 2 * e - p + 3 - eps
    a = c; dt = d
    for t in range(0, s0 + 1):
        lt = L(p - t, be, dt.lo, dt.hi)
        assert a.agree(glue(dt, lt, e), a.lo, a.hi - 1)
        a = a.F(66); dt = dt.F(66)
    D = d if eps == 0 else tailflip(d, be, e, s0 + 1)
    # D is (be^eps)-tailed with edge e
    nb = be ^ eps
    assert all(D.cell(i) == (nb, nb) for i in range(e + 1, D.hi)) and D.cell(e) != (nb, nb)
    a = c.Fn(66, s0 + 1 + eps); Dt = D.Fn(66, s0 + 1 + eps)
    for t in range(s0 + 1 + eps, s0 + 25):
        lo = max(a.lo, Dt.lo); hi = min(a.hi, Dt.hi)
        Dt2 = mk(lo, hi, Dt.cell); kk = K(t + U, nb, lo, hi)
        assert mk(lo, hi, a.cell).agree(glue(Dt2, kk, e), lo, hi - 1), (trial, t)
        a = a.F(66); Dt = Dt.F(66)
    # when eps=1 the state at time s0+1 is not yet of the form D' + K_{.,.} on [e+1,inf)
    if eps == 1:
        a = c.Fn(66, s0 + 1)
        assert a.cell(e + 1) == (be, 1 - be)
    cnt[eps] = cnt.get(eps, 0) + 1
print('bounce lemma ok', cnt)
print('ALL CHECKS PASSED')
