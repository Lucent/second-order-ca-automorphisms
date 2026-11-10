"""Rule 91, background 11.  Which pairs (E,D) of objects with isolated zeros are mutually transparent in both
parities, S(E,D_b) = (D_b,E) for b = 0,1 ?   (For h in Aut, E = h Delta^+, D = h Delta^- must be such a pair.)
Also: test of the first-ray lemma and of the second-ray parity rule."""
import sys, itertools, random
from scat import solve
RULE = 91
def cores(maxlen):
    """cores = words starting and ending with 0, no two adjacent zeros (isolated zeros), length <= maxlen"""
    out = []
    for n in range(1, maxlen+1):
        for w in itertools.product((0, 1), repeat=n):
            if w[0] or w[-1]: continue
            if any(w[i] == 0 and w[i+1] == 0 for i in range(n-1)): continue
            out.append(w)
    return out
def transparent(pc, xc, off):
    M = 3*(len(pc)+len(xc)) + 24
    M += M % 2
    A = [1]*M; B = [1]*M
    for i, v in enumerate(pc): A[4+i] = v
    for i, v in enumerate(xc): B[4+off+i] = v
    Y = solve(RULE, A, B)
    g = [Y[M-1 - ((M-1-p) % 2)][p] for p in range(M)]
    o = [Y[q][M-1 - ((M-1-q) % 2)] for q in range(M)]
    return g[:M-6] == B[:M-6] and o[:M-6] == A[:M-6]
if __name__ == "__main__":
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    cs = cores(L)
    print(len(cs), "cores with isolated zeros, length <=", L)
    good = []
    for pc in cs:
        for xc in cs:
            if transparent(pc, xc, 0) and transparent(pc, xc, 1): good.append((pc, xc))
    print("mutually transparent pairs:", len(good))
    for g in good[:60]: print("  ", "".join(map(str, g[0])), "".join(map(str, g[1])))
    # objects transparent to themselves-mirror? and to the bare defect
    tb = [c for c in cs if transparent(c, (0,), 0) and transparent(c, (0,), 1)]
    print("cores transparent to the bare defect (as probe, both parities):", len(tb), ["".join(map(str, c)) for c in tb][:40])
