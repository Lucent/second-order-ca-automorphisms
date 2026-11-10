"""Rules 219R and 126R (Section 6.6, Lemmas on bonds and cars), on all contexts: bonded cells evolve by
pi: (x,y) -> (y, 1-x); the bond bits are invariant under F and F^-1; the maps phi_s on cars.
Cells s = 2x+y as in centralizer.py.
"""
import sys, itertools

from centralizer import make_rule
pi  = lambda s: ((s & 1) << 1) | (1 - (s >> 1))     # (x,y) -> (y, 1-x)
kap = lambda s: 3 - s
A = range(4)
for rule in (219, 126):
    f, Floc, Finv = make_rule(rule)
    bond = (lambda a, b: a == b) if rule == 219 else (lambda a, b: b == kap(a))
    # (1) bonded pair evolves by pi in every context: window (l, a, b, r) with bond(a,b)
    bad = 0
    for l, a, b, r in itertools.product(A, repeat=4):
        if bond(a, b):
            if Floc(l, a, b) != pi(a) or Floc(a, b, r) != pi(b): bad += 1
    # (2) bond invariance both ways: bond(F c_i, F c_{i+1}) <=> bond(c_i,c_{i+1}), all 4-windows
    bad2 = sum(1 for l, a, b, r in itertools.product(A, repeat=4)
               if bond(Floc(l, a, b), Floc(a, b, r)) != bond(a, b))
    bad3 = sum(1 for l, a, b, r in itertools.product(A, repeat=4)
               if bond(Finv(l, a, b), Finv(a, b, r)) != bond(a, b))
    # pi commutes with kappa, order of pi
    print(rule, 'pi-evolution failures', bad, '| bond not F-invariant', bad2, '| not F^-1-invariant', bad3)
print('pi cycle on cells:', [0, pi(0), pi(pi(0)), pi(pi(pi(0)))], 'pi^2 == kappa:', all(pi(pi(s)) == kap(s) for s in A))
# (3) car dynamics. 219: car q with both neighbours s (wall cells): new car = Floc(s, a, s) for a != s.
f, Floc, _ = make_rule(219)
for s in A:
    phi = {a: Floc(s, a, s) for a in A if a != s}
    print('219 car map between s=%d walls:' % s, phi, ' image avoids pi(s)=%d:' % pi(s), pi(s) not in phi.values())
f, Floc, _ = make_rule(126)
for b in A:
    a = kap(b)   # 126 train: left neighbour a = kappa(b), right neighbour b
    phi = {q: Floc(a, q, b) for q in A if q not in (a, b)}
    print('126 car map, left nb %d right nb %d:' % (a, b), phi, ' == pi:', all(phi[q] == pi(q) for q in phi),
          ' image avoids {pi a, pi b}:', not set(phi.values()) & {pi(a), pi(b)})
