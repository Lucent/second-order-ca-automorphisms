"""Brute-force closure of the bare points under (C0)-(C3) inside the universe of all
sequences with core <= NC (mod translation; Phi^+ and Phi^- identified through rev, as C_L is rev-invariant).
argv: rule NC L1,L2,... [props]   (props: the closure properties used, default 0123)   Points whose core exceeds NC are discarded (so the result is a LOWER bound for C_L)."""
import sys, time
from rc import *
RULE = int(sys.argv[1]); NC = int(sys.argv[2]); Ls = [int(x) for x in sys.argv[3].split(',')]
USE = sys.argv[4] if len(sys.argv) > 4 else '0123'      # which closure properties to use
tails = (1,) if RULE == 95 else (0, 1)
objs = objects(NC, RULE, tails)
idx = {canon(s): i for i, s in enumerate(objs)}
n = len(objs)
print('rule', RULE, 'universe: core <=', NC, 'size', n, flush=True)
t0 = time.time()
table = {}      # (i,j,par) -> (g,o) indices or None
nsim = 0
for i, P in enumerate(objs):
    for j, X in enumerate(objs):
        if P[0] != X[0]: continue
        for par in (0, 1):
            G, O, _ = scatter(P, shift(X, par), RULE); nsim += 1
            table[(i, j, par)] = (idx.get(canon(G)), idx.get(canon(O)))
print('scattering table:', nsim, 'simulations', round(time.time() - t0), 's', flush=True)
# duality check (Prop. dual) on the table: S(mu Om, mu Gam) = (mu X, mu P) whenever all are in the universe
rv = [idx[canon(rev(s))] for s in objs]
# parity bookkeeping for the duality check is done by direct simulation on a sample instead
import random; random.seed(7); bad = 0; nd = 0
for _ in range(3000):
    P = random.choice(objs); X = shift(random.choice(objs), random.randint(0, 1))
    if P[0] != X[0]: continue
    G, O, _ = scatter(P, X, RULE)
    G2, O2, _ = scatter(rev(O), rev(G), RULE); nd += 1
    if (G2, O2) != (rev(X), rev(P)): bad += 1
print('duality S(mu Omega, mu Gamma)=(mu X, mu P): sample', nd, 'failures', bad, flush=True)

bare = [s for s in objs if len(kinks(s)) <= 1] if RULE != 95 else [s for s in objs if len(s[2]) <= 1]
for L in Ls:
    wd = [frozenset(words(s, L)) for s in objs]
    C = set(idx[canon(s)] for s in bare)
    rounds = 0; via = {}
    while True:
        rounds += 1; before = len(C)
        if '1' in USE: C |= {rv[i] for i in C}
        if '2' in USE or '3' in USE:
            for (i, j, par), (g, o) in table.items():
                if '2' in USE and i in C and j in C:
                    if g is not None and g not in C: C.add(g); via[g] = 'C2'
                    if o is not None and o not in C: C.add(o); via[o] = 'C2'
                if '3' in USE and g is not None and o is not None and g in C and o in C:
                    if i not in C: C.add(i); via[i] = 'C3'
                    if j not in C: C.add(j); via[j] = 'C3'
        if '0' in USE:
            W = set()
            for i in C: W |= wd[i]
            for i in range(n):
                if i not in C and wd[i] <= W: C.add(i); via[i] = 'C0'
        if len(C) == before: break
    miss = [objs[i] for i in range(n) if i not in C]
    from collections import Counter
    print('rule', RULE, 'L', L, 'props', USE, ': closure has', len(C), 'of', n, 'rounds', rounds, 'first-reached-by', dict(Counter(via.values())), flush=True)
    if miss:
        miss.sort(key=lambda s: len(s[2]))
        print('   missing (shortest):', [''.join(map(str, (s[0],) + s[2] + (s[3],))) for s in miss[:12]])
        if RULE == 95:
            def d4(s): 
                z = [x for x in range(len(s[2])) if s[2][x] == 0]
                return any(b - a == 4 for a in z for b in z)
            print('   missing with two zeros at distance 4:', sum(d4(s) for s in miss), ' present with distance 4:', sum(d4(objs[i]) for i in C),
                  ' missing without:', sum(not d4(s) for s in miss))
            def evd(s):
                z = [x for x in range(len(s[2])) if s[2][x] == 0]
                return sorted({b - a for a in z for b in z if b > a and (b - a) % 2 == 0})
            print('   same-class zero distances occurring in closure:', sorted({d for i in C for d in evd(objs[i])}))
