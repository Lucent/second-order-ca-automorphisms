"""The identities and inequalities of the lemmas on pairs of bare kinks, the seed event, the exchange with an
isolated zero, and the prepending of a kink and of an isolated zero, by direct simulation.
argv: [RULE [NC [NRANDOM]]]  (defaults 94, 9): all objects with core <= NC; if NRANDOM is given, NRANDOM random
objects with cores of length 14..40 instead."""
import sys, random
from rc import *
RULE = int(sys.argv[1]) if len(sys.argv) > 1 else 94
NC = int(sys.argv[2]) if len(sys.argv) > 2 else 9
fails = {}
count = {}
def chk(name, cond, info=None):
    count[name] = count.get(name, 0) + 1
    if not cond:
        fails.setdefault(name, []).append(info)
def S(P, X, **kw): return scatter(P, X, RULE, **kw)

# step 0: the bare points and isolated zeros are points of Phi^+-
for a in range(-2, 3):
    for nm, s in (('K', K(a)), ('Kbar', Kb(a)), ('Z', Z(a))):
        for sg in (1, -1):
            chk('0 bare points in Phi: ' + nm, in_phi(s, sg, RULE), (a, sg))

if fails:
    print('rule', RULE, ': STEP 0 FAILS -- bare points are not points of Phi^+-; nothing later is defined')
    for k in count: print('%-55s n=%-7d fail=%d' % (k, count[k], len(fails.get(k, []))))
    sys.exit()
def diag_check(name, P, X, formula, region, firesexp=None, T=60):
    G, O, info = S(P, X, record=True, T=T)
    Y = info['Y']; q0, p0 = info['q0'], info['p0']
    bad = 0; n = 0
    for q in range(q0, q0 + 40):
        for p in range(p0, p0 + 40):
            if (q - p) % 2: continue
            if not region(q, p): continue
            y = Y(q, p)
            if y is None: continue
            n += 1
            if y != formula(q, p): bad += 1
    chk(name + ' [diagram]', bad == 0 and n > 0, (P, X, bad, n))
    if firesexp is not None:
        chk(name + ' [firings]', info['fires'] == sorted(firesexp), (P, X, info['fires'], firesexp))
    return G, O

R = range(-4, 5)
# Lemma pair and (cross94)
for a in R:
    for b in R:
        th = lambda a_: (lambda x: 1 if x >= a_ else 0)
        G, O = diag_check('pair(i) m=0', K(a), K(b), lambda q, p: (q >= a) ^ (p >= b), lambda q, p: True, [])
        chk('cross94 first', (G, O) == (Kb(b), Kb(a)), (a, b, G, O))
        if (a - b) % 2:
            G, O = diag_check('pair(i) m=1 a!=b', Kb(a), Kb(b), lambda q, p: 1 ^ (q >= a) ^ (p >= b), lambda q, p: True, [])
            chk('cross94 second, a!=b', (G, O) == (K(b), K(a)), (a, b, G, O))
        else:
            G, O = diag_check('pair(ii) m=1 a==b', Kb(a), Kb(b),
                              lambda q, p: 1 if (q < a and p < b) or (q > a and p > b) else 0, lambda q, p: True,
                              [(a - 1, b - 1), (a + 1, b + 1)])
            chk('cross94 second, a==b (delay)', (G, O) == (K(b + 1), K(a + 1)), (a, b, G, O))

# Lemma B (seed event)
for A in range(2, 13):
    for a in range(-3, 4):
        for b in range(-3, 4):
            r = b + A
            X = fromfun(lambda p: 0 if p < b else 1 if p < r else 0 if p in (r, r + 1) else 1, b - 1, r + 3)
            if (a - b - A) % 2: continue
            for corner in (None, (a - 2 - ((a - b) % 2), b - 2), (a - 4 - ((a - b) % 2), b - 6)):
                G, O, fires = S(K(a), X, corner=corner)
                Gexp = fromfun(lambda p: 1 if p < b else 0 if p <= r else 1 if p in (r + 1, r + 2) else 0 if p == r + 3 else 1, b - 1, r + 5)
                chk('lem B output', (G, O) == (Gexp, Z(a + 1)), (A, a, b, G, O))
                chk('lem B firings', fires == sorted([(a - 1, r - 1), (a + 2, r + 4)]), (A, a, b, fires))
            # Cor zero: X = [b;[r;K^-_{r+2}]]
            if A == 2:
                chk('cor zero: X=[b;[r;K_{r+2}]]', X == norm(0, b, [1, 1, 0, 0], 1), X)

# Lemma exchange
for a in R:
    for d in R:
        if (a - d) % 2 == 0:
            G, O = diag_check('exch(i)', Kb(a), Z(d), lambda q, p: 1 ^ (q == a + 1), lambda q, p: p >= d + 2,
                              [(a - 1, d - 1), (a + 2, d + 2)])
            chk('exch(i) output', (G, O) == (K(d + 2), Z(a + 1)), (a, d, G, O))
            # full closed form of Y_Sigma
            def YS(q, p):
                x, y = q - a, p - d
                return 1 ^ (x >= 0) ^ (y == 0) ^ (x % 2 == 0 and x >= 0 and y >= 0) ^ (x % 2 == 1 and x >= 3 and y >= 3)
            diag_check('exch(i) Y_Sigma closed form', Kb(a), Z(d), YS, lambda q, p: True)
            u = d
            diag_check('exch(ii)', Z(a), Kb(u), lambda q, p: 1 ^ (p == u + 1), lambda q, p: q >= a + 2)
        else:
            G, O = diag_check('exch(iii)', Z(a), Z(d), lambda q, p: 1 ^ (q == a) ^ (p == d), lambda q, p: True, [])
            chk('exch(iii) output', (G, O) == (Z(d), Z(a)), (a, d, G, O))

# Lemma B closed form Y_Sigma on G' is internal to the proof; checked through output+firings.

def first_kink(U): return kinks(U)[0]
def prepend(b, U):
    u = first_kink(U); c = U[0]; assert b <= u - 2
    return fromfun(lambda p: (1 ^ c) if p < b else val(U, p), b - 1, U[1] + len(U[2]) + 1)
def prependZ(d, U):
    u = first_kink(U); assert U[0] == 1 and d <= u - 3
    return fromfun(lambda p: 0 if p == d else val(U, p), d - 1, U[1] + len(U[2]) + 1)
def glue(X, X2, p2):
    return fromfun(lambda p: val(X, p) if p <= p2 + 1 else val(X2, p), min(X[1], X2[1]) - 2, max(X[1] + len(X[2]), X2[1] + len(X2[2])) + 2)

objs = [U for U in objects(NC, RULE if RULE != 95 else 94) if not is_bg(U)]
if len(sys.argv) > 3:
    random.seed(1)
    # random larger objects
    def rnd(n):
        while True:
            lt = random.randint(0, 1); w = [1 ^ lt]
            while len(w) < n:
                nx = random.randint(0, 1)
                if len(w) >= 2 and (w[-2], w[-1], nx) == (0, 1, 0): continue
                if len(w) == 1 and (lt, w[-1], nx) == (0, 1, 0): continue
                w.append(nx)
            rt = 1 ^ w[-1]
            if (w[-2] if len(w) > 1 else lt, w[-1], rt) == (0, 1, 0): continue
            return (lt, 0, tuple(w), rt)
    objs = [rnd(random.randint(14, 40)) for _ in range(int(sys.argv[3]))]
for U in objs:
    u = first_kink(U); c = U[0]
    for dist in (2, 3, 4, 5, 8):
        b = u - dist
        Xi = prepend(b, U)
        chk('prepend: [b;U] kinks = b + kinks(U)', kinks(Xi) == [b] + kinks(U), (b, U))
        for a in range(-1, 3):
            if c == 0:
                if (a - b) % 2 == 0: continue
                Pi, Pp, Kout = Kb(a), K(a), K(b)
            else:
                if (a - u) % 2: continue
                Pi, Pp, Kout = K(a), Kb(a), Kb(b)
            G, O, _ = S(Pi, Xi)
            Gp, Op, _ = S(Pp, U)
            name = 'prepend c=%d' % c
            chk(name + ': S(P_i,X_i)=(K*Gamma\',Omega\')', G == glue(Kout, Gp, b) and O == Op, (b, a, U, G, O, Gp, Op))
            chk(name + ': Gamma\' not a background', not is_bg(Gp), (b, a, U))
            if not is_bg(Gp):
                up = first_kink(Gp)
                if c == 0:
                    chk(name + ': gamma\'=1 on p<u (u\'>=u, new c=1)', Gp[0] == 1 and up >= u, (b, a, U, Gp))
                else:
                    chk(name + ': gamma\'=0 on p<=u (u\'>=u+1, new c=0)', Gp[0] == 0 and up >= u + 1, (b, a, U, Gp))
                chk(name + ': K*Gamma\' = [b;Gamma\']', up >= b + 2 and G == prepend(b, Gp), (b, a, U, G, Gp))
    if c == 1:
        for dist in (3, 4, 5, 6, 9):
            d = u - dist
            X = prependZ(d, U)
            chk('prependZ: <d;U> kinks', kinks(X) == [d, d + 1] + kinks(U), (d, U))
            for a in range(-1, 3):
                if dist >= 4:
                    if (a - d) % 2: continue
                    G, O, _ = S(Kb(a), X)
                    Gp, Op, _ = S(Z(a + 1), U)
                    nm = 'prependZ u>=d+4'
                    chk(nm + ': S=(K_{d+2}*Gamma\',Omega\')', G == glue(K(d + 2), Gp, d + 2) and O == Op, (d, a, U, G, O, Gp, Op))
                    ok = (not is_bg(Gp)) and Gp[0] == 1 and first_kink(Gp) >= u
                    chk(nm + ': gamma\'=1 on p<u, not bg', ok, (d, a, U, Gp))
                    if ok: chk(nm + ': result=[d+2;Gamma\']', G == prepend(d + 2, Gp), (d, a, U, G, Gp))
                else:
                    if (a - d) % 2 == 0: continue
                    G, O, _ = S(Z(a), X)
                    Gp, Op, _ = S(Z(a), U)
                    nm = 'prependZ u=d+3'
                    chk(nm + ': S=(Z_d*Gamma\',Omega\')', G == glue(Z(d), Gp, d + 1) and O == Op, (d, a, U, G, O, Gp, Op))
                    ok = (not is_bg(Gp)) and Gp[0] == 1 and first_kink(Gp) >= u + 1
                    chk(nm + ': gamma\'=1 on p<=u (u\'>=d+4)', ok, (d, a, U, Gp))
                    if ok: chk(nm + ': result=<d;Gamma\'>', G == prependZ(d, Gp), (d, a, U, G, Gp))

print('rule', RULE, 'objects', len(objs))
for k in count:
    fl = fails.get(k, [])
    print('%-55s n=%-7d fail=%d' % (k, count[k], len(fl)), ('  e.g. ' + str(fl[0])[:150]) if fl else '')
