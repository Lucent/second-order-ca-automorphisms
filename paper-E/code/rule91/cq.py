"""Additive conserved quantities of F_f with density h on r consecutive cells: sum_i h(c_i..c_{i+r-1}) invariant under F.
h(F(w)) - h(w[1..r]) = g(w[1..r+1]) - g(w[0..r]) for all w in A^{r+2}.  Exact rational nullspace."""
import sys, itertools
from fractions import Fraction

def rref(rows):
    rows = [[Fraction(a) for a in r] for r in rows]
    piv = []; rk = 0; ncol = len(rows[0])
    for c in range(ncol):
        p = None
        for i in range(rk, len(rows)):
            if rows[i][c] != 0: p = i; break
        if p is None: continue
        rows[rk], rows[p] = rows[p], rows[rk]
        inv = 1 / rows[rk][c]
        rows[rk] = [a*inv for a in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] != 0:
                fct = rows[i][c]; rr = rows[rk]
                rows[i] = [a - fct*b for a, b in zip(rows[i], rr)]
        piv.append(c); rk += 1
        if rk == len(rows): break
    return rows[:rk], piv
def nullspace(rows):
    R, piv = rref(rows); ncol = len(rows[0]); free = [c for c in range(ncol) if c not in piv]
    out = []
    for fc in free:
        v = [Fraction(0)]*ncol; v[fc] = Fraction(1)
        for i, pc in enumerate(piv): v[pc] = -R[i][fc]
        out.append(v)
    return out
def rank(rows): return len(rref(rows)[1])
rule = int(sys.argv[1]); r = int(sys.argv[2])
f = lambda a, b, c: (rule >> (4*a + 2*b + c)) & 1
def Floc(s0, s1, s2):
    x, y = s1 >> 1, s1 & 1
    return (y << 1) | (x ^ f(s0 & 1, y, s2 & 1))
nh = 4**r; ng = 4**(r+1)
rows = []
for w in itertools.product(range(4), repeat=r+2):
    row = [0]*(nh+ng)
    fw = 0
    for i in range(r): fw = fw*4 + Floc(w[i], w[i+1], w[i+2])
    mid = 0
    for i in range(1, r+1): mid = mid*4 + w[i]
    gr = 0
    for i in range(1, r+2): gr = gr*4 + w[i]
    gl = 0
    for i in range(0, r+1): gl = gl*4 + w[i]
    row[fw] += 1; row[mid] -= 1; row[nh+gr] -= 1; row[nh+gl] += 1
    rows.append(row)
ns = nullspace(rows)
print("nullity", len(ns))
# project to h-part, quotient by trivial: constants and coboundaries k(w[0..r-2]) - k(w[1..r-1])
H = [list(v[:nh]) for v in ns]
triv = [[1]*nh]
for kk in itertools.product(range(4), repeat=r-1):
    row = [0]*nh
    for w in itertools.product(range(4), repeat=r):
        idx = 0
        for a in w: idx = idx*4 + a
        if tuple(w[:r-1]) == kk: row[idx] += 1
        if tuple(w[1:]) == kk: row[idx] -= 1
    triv.append(row)
rt = rank(triv)
ra = rank(triv + H)
print("rank trivial", rt, "rank all h", ra, "nontrivial", ra - rt)
# print a basis of nontrivial h modulo trivial, in terms of cell words
basis = []
cur = list(triv); currank = rt
for v in H:
    if rank(cur + [v]) > currank:
        cur.append(v); currank += 1; basis.append(v)
for v in basis:
    print("h:")
    for w in itertools.product(range(4), repeat=r):
        idx = 0
        for a in w: idx = idx*4 + a
        if v[idx] != 0: print("   ", w, v[idx])
