"""Generic scattering for f in {91,94,95}: Goursat solver in light-cone coordinates.
Y[q][p] for q,p >= 0 (only q=p mod 2 are lattice points). A(q) on p in {0,1}, B(p) on q in {0,1}."""
import itertools, sys, random
def N_of(rule):
    f = lambda a, b, c: (rule >> (4*a + 2*b + c)) & 1
    return lambda a, b, c: f(a, b, c) ^ a ^ c
def solve(rule, A, B):
    """A, B lists (same parity conventions): A[q] = Y[q, p in {0,1}, p=q mod 2], B[p] = Y[q in {0,1}, q = p mod 2]."""
    N = N_of(rule)
    Q, P = len(A), len(B)
    Y = [[None]*P for _ in range(Q)]
    for q in range(Q): Y[q][q % 2] = A[q]
    for p in range(P): Y[p % 2][p] = B[p]
    for q in range(2, Q):
        for p in range(2 + (q % 2), P, 2):
            Y[q][p] = Y[q][p-2] ^ Y[q-2][p-2] ^ Y[q-2][p] ^ N(Y[q][p-2], Y[q-1][p-1], Y[q-2][p])
    return Y
def phi_words(rule, n):
    N = N_of(rule)
    out = []
    def rec(w):
        if len(w) == n: out.append(tuple(w)); return
        for b in (0, 1):
            if len(w) >= 2 and N(w[-2], w[-1], b): continue
            rec(w + [b])
    rec([])
    return out
