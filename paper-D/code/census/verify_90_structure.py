"""90R: verify that the full solution space (all block maps of radius R commuting with F) is spanned by
  (A) affine maps  c -> (aI + bM) c + k,  k in {(0,0),(1,1)}                                   [dim 4R+1]
  (U) g_h(c) = ( h(u), sigma^{-1} h(u) ),  u = x + sigma^{-1} y  (u_i = x_i xor y_{i-1}),  h any boolean
      block map on u with neighbourhood [-(R-2), R]                                               [dim 2^{2R-1}]
  (W) g'_h(c) = ( h(w), sigma h(w) ),  w = x + sigma y  (w_i = x_i xor y_{i+1}),  h with neighbourhood [-R, R-2]
Reason: u o F = sigma^{-1} o u and F(v, sigma^{-1} v) = (sigma^{-1} v, sigma^{-2} v), i.e. u is a factor map onto
the binary shift and v -> (v, sigma^{-1} v) embeds it back (same for w with sigma). Prediction: dim = 2^{2R} + 1."""
import sys, itertools
from centralizer import *
from linear_90 import build_rows, kernel_basis
from analyze_90 import affine_table_basis, rank_of

def family_tables(R, kind):
    n = 2 * R + 1; N = 4 ** n; vecs = []
    lo, hi = ((-(R - 2), R) if kind == "U" else (-R, R - 2))
    cells = list(range(lo, hi + 1))
    if R == 1: cells = [c for c in cells if -R <= c <= R - 1]        # R=1: u-neighbourhood must fit; handled below
    def u_at(w, j):            # u_{j} (kind U) or w_{j} (kind W) of the window w, centre index R
        s = w[R + j]; x = s >> 1
        if kind == "U": t = w[R + j - 1] & 1
        else: t = w[R + j + 1] & 1
        return x ^ t
    # X_out = h(u)_0 uses u_j for j in cells; Y_out = h(u)_{-1} (U) or h(u)_{+1} (W)
    sh = -1 if kind == "U" else 1
    for mono in range(1 << len(cells)):           # monomial h = product of u_j over j in mono
        vec = 0
        for e in range(N):
            w = index_word(e, n)
            try:
                X = 1
                for i, j in enumerate(cells):
                    if mono >> i & 1: X &= u_at(w, j)
                Y = 1
                for i, j in enumerate(cells):
                    if mono >> i & 1: Y &= u_at(w, j + sh)
            except IndexError:
                X = Y = None; break
            if X: vec |= 1 << (2 * e)
            if Y: vec |= 1 << (2 * e + 1)
        if X is not None: vecs.append(vec)
    return vecs

for R in (1, 2):
    rows, nvar = build_rows(R); rank, basis = kernel_basis(rows, nvar)
    A = affine_table_basis(R); U = family_tables(R, "U"); W = family_tables(R, "W")
    f, F_loc, _ = make_rule(90)
    # every U/W table commutes (direct check), and lies in the solution space
    piv = {}
    for r in basis:
        while r:
            p = r & -r
            if p in piv: r ^= piv[p]
            else: piv[p] = r; break
    def in_space(v):
        while v:
            p = v & -v
            if p in piv: v ^= piv[p]
            else: return False
        return True
    from linear_90 import vec_to_table
    okU = all(in_space(v) and check_commutes(BlockMap(R, vec_to_table(v, R)), F_loc, R) for v in U)
    okW = all(in_space(v) and check_commutes(BlockMap(R, vec_to_table(v, R)), F_loc, R) for v in W)
    Asol = [v for v in A]  # affine table basis is not itself inside the solution space; use intersection dims
    dA = rank_of(A) + len(basis) - rank_of(A + basis)
    span = rank_of(A + U + W)                     # A spans ALL affine tables (dim 22 at R=2), so intersect with solutions:
    # span of (affine solutions) + U + W  =  (A+U+W) ∩ solutions ; compute via dim formula
    d_AUW_cap_sol = rank_of(A + U + W) + len(basis) - rank_of(A + U + W + basis)
    print(f"R={R}: solution dim {len(basis)}; affine solutions dim {dA}; |U|={len(U)} tables (all commute & in space: {okU}); "
          f"|W|={len(W)} ({okW}); dim of (affine solutions + U + W) = {d_AUW_cap_sol}; prediction 2^(2R)+1 = {4**R+1}; "
          f"EXACT SPAN: {d_AUW_cap_sol == len(basis)}")
