"""The induction of the scattering-closure theorem (every point with finitely many kinks lies in C_L) run
mechanically: every scattering of the proof is simulated and compared with the asserted form.
argv: L maxcore [RULE [NRANDOM]]  (default rule 94): all objects with core <= maxcore; if NRANDOM is given,
NRANDOM random objects with cores of length 14..30 instead."""
import sys
from rc import *
L = int(sys.argv[1]); NC = int(sys.argv[2]); RULE = int(sys.argv[3]) if len(sys.argv) > 3 else 94
nsc = 0; maxcore = 0
class Fail(Exception): pass
def S(P, X):
    global nsc, maxcore
    nsc += 1
    G, O, _ = scatter(P, X, RULE)
    maxcore = max(maxcore, len(G[2]), len(O[2]))
    return G, O
def need(c, *info):
    if not c: raise Fail(info)
def prepend(b, U):
    c = U[0]
    return fromfun(lambda p: (1 ^ c) if p < b else val(U, p), b - 1, U[1] + len(U[2]) + 1)
def prependZ(d, U):
    return fromfun(lambda p: 0 if p == d else val(U, p), d - 1, U[1] + len(U[2]) + 1)

def chain(b, U, tag):
    """Lemma prepend: U object (assumed in C_L), b <= u-2.  Returns number of steps."""
    need(not is_bg(U), tag, 'U background')
    u0 = kinks(U)[0]; need(b <= u0 - 2, tag, 'distance')
    Ui = U; hist = []
    for i in range(2 * L):
        u = kinks(Ui)[0]; c = Ui[0]
        Xi = prepend(b, Ui)
        # alternate the two admissible probe positions to exercise both
        if c == 0:
            a = b + 1 + 2 * (i % 3); Pi, Pp = Kb(a), K(a)
        else:
            a = u + 2 * (i % 3); Pi, Pp = K(a), Kb(a)
        G, O = S(Pi, Xi)
        Gp, Op = S(Pp, Ui)
        need(not is_bg(Gp), tag, 'step', i, 'Gamma\' background')
        up = kinks(Gp)[0]
        need(Gp[0] == 1 ^ c, tag, 'step', i, 'value left of first kink does not flip')
        need(up >= u + (1 if c == 1 else 0), tag, 'step', i, 'c=%d' % c, 'first kink bound', u, up)
        need(G == prepend(b, Gp) and O == Op, tag, 'step', i, 'c=%d' % c, 'S(P_i,X_i) != ([b;U_{i+1}],Omega_i)')
        hist.append(u)
        Ui = Gp
    un = kinks(Ui)[0]
    need(un >= u0 + L and un >= b + L + 2, tag, 'final gap', u0, un)
    Xn = prepend(b, Ui)
    bare = (1 ^ Ui[0], b, (), Ui[0])
    need(words(Xn, L) <= (words(Ui, L) | words(bare, L)), tag, '(C0) window condition')
    # stronger, as in the proof: windows starting at s>=b are windows of U_n at the same place, s<b of the bare kink
    for s in range(b - L - 2, Ui[1] + len(Ui[2]) + 3):
        w = tuple(val(Xn, s + j) for j in range(L))
        ref = Ui if s >= b else bare
        need(w == tuple(val(ref, s + j) for j in range(L)), tag, '(C0) positional', s)

def step_Z(d, U, tag):
    """Lemma prependZ: U in C_L with value 1 left of first kink u, d<=u-3: <d;U> in C_L."""
    u = kinks(U)[0]; need(U[0] == 1 and d <= u - 3, tag, 'hyp')
    X = prependZ(d, U)
    if u >= d + 4:
        a = d
        G, O = S(Kb(a), X); Gp, Op = S(Z(a + 1), U)
        need(not is_bg(Gp) and Gp[0] == 1 and kinks(Gp)[0] >= u, tag, 'Z case1 gamma\'')
        need(G == prepend(d + 2, Gp) and O == Op, tag, 'Z case1 identity')
        chain(d + 2, Gp, tag + ('Zc1',))
    else:
        a = d + 1
        G, O = S(Z(a), X); Gp, Op = S(Z(a), U)
        need(not is_bg(Gp) and Gp[0] == 1 and kinks(Gp)[0] >= u + 1, tag, 'Z case2 gamma\'')
        need(G == prependZ(d, Gp) and O == Op, tag, 'Z case2 identity')
        step_Z(d, Gp, tag + ('Zc2',))

def top(X):
    ks = kinks(X); k = len(ks)
    if k <= 1: return 'bare'
    b1, b2 = ks[0], ks[1]
    if b2 - b1 >= 2:
        U = fromfun(lambda p: val(X, p) if p >= b1 else val(X, b1), b1 - 1, X[1] + len(X[2]) + 1)
        need(kinks(U) == ks[1:], 'U kinks'); need(prepend(b1, U) == X, 'X=[b1;U]')
        chain(b1, U, (X,)); return 'prepend'
    need(val(X, b1) == 0 and val(X, b1 - 1) == 1 and val(X, b1 + 1) == 1, 'isolated zero form')
    if k == 2:
        need(X == Z(b1), 'Z'); return 'zero'
    need(ks[2] - b2 >= 2, 'b3-b2>=2')
    U = fromfun(lambda p: val(X, p) if p > b1 else 1, b1 - 1, X[1] + len(X[2]) + 1)
    need(kinks(U) == ks[2:] and U[0] == 1 and ks[2] >= b1 + 3 and prependZ(b1, U) == X, 'U for Z')
    step_Z(b1, U, (X,)); return 'prependZ'

# Corollary zero (Lemma B with A=2), for both parities
for b in (0, 1):
    for a in (0, 1):
        if (a - b) % 2: continue
        X = prepend(b, prepend(b + 2, K(b + 4)))
        G, O = S(K(a), X)
        if O != Z(a + 1): print('rule', RULE, 'Corollary zero / Lemma B (A=2) FAILS:', a, b, 'response', O)
objs = objects(NC, 94)
if len(sys.argv) > 4:
    import random; random.seed(5)
    def rnd(n):
        while True:
            lt = random.randint(0, 1); w = [1 ^ lt]
            while len(w) < n: w.append(random.randint(0, 1))
            rt = 1 ^ w[-1]; full = [lt, lt] + w + [rt, rt]
            if all(tuple(full[i:i+3]) != (0, 1, 0) for i in range(len(full) - 2)): return (lt, 0, tuple(w), rt)
    objs = [rnd(random.randint(14, 30)) for _ in range(int(sys.argv[4]))]
kinds = {}; fl = []
for X in objs:
    try:
        r = top(X); kinds[r] = kinds.get(r, 0) + 1
    except Fail as e:
        fl.append((X, e.args))
print('rule', RULE, 'L', L, 'maxcore', NC, 'objects', len(objs), kinds, 'scatterings', nsc, 'largest core seen', maxcore, 'FAILURES', len(fl))
for x in fl[:5]: print('   ', str(x)[:400])
if fl:
    from collections import Counter
    print(Counter(str([z for z in e if isinstance(z, (str,int))]) for _, e in fl).most_common(8))
