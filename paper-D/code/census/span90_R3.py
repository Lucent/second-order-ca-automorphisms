"""90R: the three families A (affine commutant), U = iota∘h∘u, W = iota'∘h∘w span a space of dimension 2^{2R}+1
inside the centralizer, for R = 1..3 (every family table re-checked to commute by direct evaluation)."""
import sys, itertools
from centralizer import make_rule, index_word, BlockMap, check_commutes
from classify_r4 import check_commutes_np_R
f, F_loc, Finv_loc = make_rule(90)
def rank_of(vecs):
    piv = {}
    for r in vecs:
        while r:
            p = r & -r
            if p in piv: r ^= piv[p]
            else: piv[p] = r; break
    return len(piv)
def tab2vec(tab):
    v = 0
    for e, s in enumerate(tab): v |= ((s >> 1) << (2 * e)) | ((s & 1) << (2 * e + 1))
    return v
def family(R, kind):
    n = 2 * R + 1; N = 4 ** n
    lo, hi = ((-(R - 2), R) if kind == "U" else (-R, R - 2))
    cells = [j for j in range(lo, hi + 1)]
    sh = -1 if kind == "U" else 1
    def lin(w, j):   # u_j = x_j + y_{j-1}  /  w_j = x_j + y_{j+1}, window centred at R
        s = w[R + j]; t = w[R + j + sh]
        return (s >> 1) ^ (t & 1)
    tabs = []
    for mono in range(1 << len(cells)):
        tab = []
        for e in range(N):
            w = index_word(e, n)
            X = 1; Y = 1
            for i, j in enumerate(cells):
                if mono >> i & 1:
                    X &= lin(w, j)
                    jj = j + sh
                    Y &= lin(w, jj) if (0 <= R + jj + sh < n and 0 <= R + jj < n) else 0
            tab.append((X << 1) | Y)
        tabs.append(tab)
    return tabs
def affine(R):
    n = 2 * R + 1; N = 4 ** n; tabs = []
    for j in range(-R, R + 1):                       # a = s^j : (x,y) -> (x_j, y_j)
        tabs.append([index_word(e, n)[R + j] for e in range(N)])
    Fm = BlockMap.from_local(1, F_loc, 1)
    for j in range(-(R - 1), R):                     # b = s^j : sigma^j F
        g = BlockMap.shift(R - 1, j).compose(Fm, R); tabs.append(list(g.table))
    tabs.append([3] * N)                             # constant (1,1)
    return tabs
for R in (1, 2, 3):
    A = affine(R); U = family(R, "U"); W = family(R, "W")
    allt = A + U + W
    ok = all(check_commutes_np_R(t, F_loc, R) for t in allt)
    print(f"R={R}: |A|={len(A)} |U|={len(U)} |W|={len(W)}; all commute: {ok}; rank(A)={rank_of([tab2vec(t) for t in A])}, "
          f"rank(U)={rank_of([tab2vec(t) for t in U])}, rank(W)={rank_of([tab2vec(t) for t in W])}, rank(A+U+W)={rank_of([tab2vec(t) for t in allt])}, 2^(2R)+1={4**R+1}")
