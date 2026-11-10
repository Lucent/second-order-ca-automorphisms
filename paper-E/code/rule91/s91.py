"""Rule 91 on the background 11: scattering of finite objects (u = 1 xor Y has finite support)."""
import sys, random
from scat import solve, N_of, phi_words
rule = int(sys.argv[1]) if len(sys.argv) > 1 else 91
def run(pc, xc, M=None):
    M = M or (4*(len(pc)+len(xc)) + 40)
    A = [1]*M; B = [1]*M
    for i, v in enumerate(pc): A[4+i] = v
    for i, v in enumerate(xc): B[4+i] = v
    Y = solve(rule, A, B)
    # outgoing: gamma(p) = Y[q][p] for large q ; o(q) for large p
    g = [Y[M-1 - ((M-1-p) % 2)][p] for p in range(M)]
    o = [Y[q][M-1 - ((M-1-q) % 2)] for q in range(M)]
    g2 = [Y[M-3 - ((M-3-p) % 2)][p] for p in range(M)]
    o2 = [Y[q][M-3 - ((M-3-q) % 2)] for q in range(M)]
    return A, B, g, o, (g[:M-12] == g2[:M-12] and o[:M-12] == o2[:M-12])
def core(s):
    z = [i for i, v in enumerate(s) if v == 0]
    return (z[0], tuple(s[z[0]:z[-1]+1])) if z else (None, ())
words = [w for n in range(1, 8) for w in phi_words(rule, n+4) if w[:2] == (1, 1) and w[-2:] == (1, 1) and w[2] == 0 and w[-3] == 0]
print(len(words), "cores")
stats = {}
ex = {}
for pw in words:
    for xw in words:
        for off in (0, 1):
            A, B, g, o, stable = run(pw[2:-2], (1,)*off + xw[2:-2])
            cp, cx, cg, co = core(A), core(B), core(g), core(o)
            if not stable: key = "unstable"
            elif g[-14:] != [1]*14 or o[-14:] != [1]*14: key = "tail changed"
            else:
                key = ("shape kept" if (cg[1] == cx[1] and co[1] == cp[1]) else "shape changed")
            stats[key] = stats.get(key, 0) + 1
            ex.setdefault(key, []).append((pw[2:-2], xw[2:-2], off, cg, co, cx[0], cp[0]))
print(stats)
for k in ex:
    for e in ex[k][:6]: print(k, e)
