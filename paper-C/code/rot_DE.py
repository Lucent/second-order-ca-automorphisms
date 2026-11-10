"""Rules 219R and 126R, Propositions 6.24 and 6.25: the automorphism D of radius 2 of 219R and the involution E
of radius 2 of 126R, verified exactly.
219R:  D(c)_q = R_s(c_q)  if (c_{q-2},c_{q-1},c_q,c_{q+1},c_{q+2}) = (s,s,a,s,s) with a != s, else c_q,
       where R_s = phi_{pi^3 s} o phi_{pi^2 s} o phi_{pi s} o phi_s  and phi_s(a) = F_loc(s,a,s)
       (R_s is four steps of F applied to an isolated car, a permutation of the cells other than s).
126R:  E(c)_q = kappa(c_q) if the window is (b, kb, q, b, kb) with q not in {b, kb}, else c_q.
Checked: commutation with F and F^-1 on all 4^7 windows, injectivity by the pair graph, the inverse, and that
the map is not sigma^a F^b and not affine.
"""
from centralizer import *
import itertools
pi = lambda s: ((s & 1) << 1) | (1 - (s >> 1)); kap = lambda s: 3 - s

def check_all(rule, g, ginv, label):
    f, Fl, Fil = make_rule(rule)
    R = g.R
    comm = check_commutes(g, Fl, R)                      # all 4^(2R+3) windows
    commi = check_commutes(g, Fil, R)                    # also with F^{-1} (redundant, cheap)
    inj = injective_pairgraph_np(g)
    idR = BlockMap.shift(2 * R, 0)
    left = ginv.compose(g, 2 * R).key() == idR.key()     # exact on all 4^(4R+1) windows
    right = g.compose(ginv, 2 * R).key() == idR.key()
    nchanged = sum(1 for e in range(4 ** (2 * R + 1)) if g.table[e] != index_word(e, 2 * R + 1)[R])
    names = named_maps(R, Fl, Fil)
    kapR = BlockMap.from_local(R, lambda s: 3 - s, 0)
    in_named = g.key() in names or any(kapR.compose(BlockMap(R, list(k)), 2 * R).key() ==
                                       BlockMap(R, g.table).compose(BlockMap.shift(R, 0), 2 * R).key() for k in [])
    print(f'{label}: radius {R}, entries differing from id {nchanged}, commutes with F on all 4^{2*R+3} windows: {comm}, '
          f'with F^-1: {commi}, pair-graph injective: {inj}, ginv o g = id: {left}, g o ginv = id: {right}, '
          f'equals some sigma^aF^b (|a|+|b|<=R): {g.key() in names}, affine: {is_affine(g)}')
    return g

# ---------------- 219
f, Fl, Fil = make_rule(219)
phi = lambda s, a: Fl(s, a, s)
def Rret(s, a):
    for t in range(4):
        a = phi(s, a); s = pi(s)
    return a
for s in range(4):
    print('219 return map R_%d on cars:' % s, {a: Rret(s, a) for a in range(4) if a != s})
def Dloc(power):
    def loc(w0, w1, w2, w3, w4):
        if w0 == w1 == w3 == w4 and w2 != w1:
            a = w2
            for _ in range(power): a = Rret(w1, a)
            return a
        return w2
    return loc
D = BlockMap.from_local(2, Dloc(1), 2); D2 = BlockMap.from_local(2, Dloc(2), 2)
check_all(219, D, D2, '219R  D')
# order 3: D o D = D2 as radius-4 tables
print('219R  D o D == D^2-table:', D.compose(D, 4).key() == BlockMap(2, D2.table).compose(BlockMap.shift(2, 0), 4).key())
# D commutes with kappa?  and D is not kappa^e sigma^a F^b: it has order 3 (shown), the group Z^2 x Z/2 has no 3-torsion.
K = BlockMap.from_local(2, kap, 0)
print('219R  D kappa == kappa D:', D.compose(K, 4).key() == K.compose(D, 4).key())

# ---------------- 126
def Eloc(w0, w1, w2, w3, w4):
    b = w0
    if w1 == kap(b) and w3 == b and w4 == kap(b) and w2 not in (b, kap(b)):
        return kap(w2)
    return w2
E = BlockMap.from_local(2, Eloc, 2)
check_all(126, E, E, '126R  E')
