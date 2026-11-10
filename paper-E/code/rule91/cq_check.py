"""Exact check (rational linear algebra): the densities of Q1 and Q2 on two adjacent cells are additive conserved
densities of F_f, i.e. there is g on 3 cells with h(F(w)) - h(w[1..2]) = g(w[1..3]) - g(w[0..2]) for all w in A^4."""
import sys, itertools
from fractions import Fraction
rule = int(sys.argv[1])
f = lambda a, b, c: (rule >> (4*a + 2*b + c)) & 1
def Floc(s0, s1, s2):
    x, y = s1 >> 1, s1 & 1
    return (y << 1) | (x ^ f(s0 & 1, y, s2 & 1))
def rank(rows):
    rows = [[Fraction(a) for a in r] for r in rows]; rk = 0
    for c in range(len(rows[0])):
        p = next((i for i in range(rk, len(rows)) if rows[i][c] != 0), None)
        if p is None: continue
        rows[rk], rows[p] = rows[p], rows[rk]
        for i in range(rk+1, len(rows)):
            if rows[i][c] != 0:
                fct = rows[i][c] / rows[rk][c]
                rows[i] = [a - fct*b for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk
def solvable(h):
    # unknown g on A^3 (64): equations for each w in A^4:  g(w1w2w3) - g(w0w1w2) = h(F(w)) - h(w1,w2)
    rows = []
    for w in itertools.product(range(4), repeat=4):
        row = [0]*65
        row[w[1]*16 + w[2]*4 + w[3]] += 1; row[w[0]*16 + w[1]*4 + w[2]] -= 1
        row[64] = h(Floc(w[0], w[1], w[2]), Floc(w[1], w[2], w[3])) - h(w[1], w[2])
        rows.append(row)
    return rank(rows) == rank([r[:64] for r in rows])
xy = lambda s: (s >> 1, s & 1)
Q1 = lambda s, t: int(s != 3) - int(s != t and s != 3 and t != 3)
Q2 = lambda s, t: int(xy(s)[1] != xy(t)[0]) - int(xy(s)[0] != xy(t)[1])
print("rule", rule, "Q1 conserved:", solvable(Q1), " Q2 conserved:", solvable(Q2))
