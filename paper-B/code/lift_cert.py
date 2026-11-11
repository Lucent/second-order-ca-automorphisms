"""The lift certificate (the supplement to the paper, the section "The lift certificate").

Usage:  python lift_cert.py RULE NMAX TMAX D

RULE is an elementary rule f(a,b,c) = a XOR g(b,c); NMAX, TMAX and D are the parameters
n_max, T_max and the largest depth D.

A node of depth k is a pair (X, V): X is the time series, up to rotation and over its minimal
period p, of the four bits (Z^a_t(k-1), Z^a_t(k), Z^b_t(k-1), Z^b_t(k)), and V is the 2-adic
valuation of the minimal period of a periodic history realizing the node.

R_0 has one node for each F-orbit of a configuration of spatial period n <= NMAX: copy a is the
comoving diagram Z_t(d) = y_t(t - d) of the orbit and copy b is the same diagram with layer 0
complemented.

test(R_k) is true if some node of R_k has v_2(p) = V and different increment sums, over the
period p, in the two copies at the step y[k+1].

R_{k+1} is the set of nodes obtained from the nodes of R_k by adding a layer k+1 in each copy
(two series if the increment sum over p is even, the series with initial value 0 if it is odd)
and keeping the joint series of layers (k, k+1) of length at most TMAX.

The program prints |R_k| and test(R_k) for k = 0, ..., D-1 and compares each R_{k+1},
k+1 <= D, with all R_j, j <= k.  Sets of at most KEEP nodes are stored and compared as sets.
Larger sets are compared through the SHA-256 digest of their sorted list of nodes; different
digests prove that the sets are different.  The program stops at the first repetition
R_{k0+rho} = R_{k0}.
"""
import sys, hashlib, itertools
from math import lcm

KEEP = 20000


def f_of(rule):
    return lambda a, b, c: (rule >> (4 * a + 2 * b + c)) & 1


def g_of(rule):
    f = f_of(rule)
    assert all(f(0, b, c) != f(1, b, c) for b in (0, 1) for c in (0, 1)), "rule is not left-permutive"
    return {(b, c): f(0, b, c) for b in (0, 1) for c in (0, 1)}


def orbits(rule, nmax):
    """The F-orbits of the configurations of spatial period n <= nmax, F(x,y) = (y, x + f(y))."""
    f = f_of(rule)
    out = []
    for n in range(1, nmax + 1):
        seen = set()
        for x in itertools.product((0, 1), repeat=n):
            for y in itertools.product((0, 1), repeat=n):
                c = (x, y)
                if c in seen:
                    continue
                cyc = [c]
                z = c
                while True:
                    X, Y = z
                    z = (Y, tuple(X[i] ^ f(Y[(i - 1) % n], Y[i], Y[(i + 1) % n]) for i in range(n)))
                    if z == c:
                        break
                    cyc.append(z)
                seen |= set(cyc)
                out.append((n, cyc))
    return out


def comoving(cyc, n, d, t):
    """Z_t(d) = y_t(t - d) for the orbit cyc."""
    return cyc[t % len(cyc)][1][(t - d) % n]


def v2(n):
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


def minper(ser):
    P = len(ser)
    for p in range(1, P + 1):
        if P % p == 0 and all(ser[t] == ser[t % p] for t in range(P)):
            return p


def canon(ser):
    p = minper(ser)
    ser = list(ser[:p])
    return min(tuple(ser[k:] + ser[:k]) for k in range(p))


def digest(R):
    return hashlib.sha256(repr(sorted(R)).encode()).hexdigest()


def run(rule, nmax, Tmax, D, verbose=True):
    g = g_of(rule)

    def incs(lo, hi):
        # phi_t for the step to the next layer: Z_{t+1}(k+1) = Z_t(k+1) + Z_{t-1}(k-1) + g(Z_t(k), Z_t(k-1))
        P = len(lo)
        return [lo[(t - 1) % P] ^ g[hi[t], lo[t]] for t in range(P)]

    def lift1(lo, hi):
        P = len(lo)
        ph = incs(lo, hi)
        if sum(ph) % 2 == 0:
            outs = []
            for z0 in (0, 1):
                z = [z0]
                for t in range(P - 1):
                    z.append(z[-1] ^ ph[t])
                outs.append(z)
            return outs
        z = [0]
        for t in range(2 * P - 1):
            z.append(z[-1] ^ ph[t % P])
        return [z]

    start = set()
    for n, cyc in orbits(rule, nmax):
        P = lcm(len(cyc), n)
        Pt = next(p for p in range(1, P + 1) if P % p == 0 and all(
            comoving(cyc, n, d, t + p) == comoving(cyc, n, d, t) for d in range(n) for t in range(P)))
        lo = [comoving(cyc, n, -1, t) for t in range(P)]
        hi = [comoving(cyc, n, 0, t) for t in range(P)]
        start.add((canon([(lo[t], hi[t], lo[t], 1 - hi[t]) for t in range(P)]), v2(Pt)))

    R = frozenset(start)
    seen = {digest(R): 0}
    kept = {0: R} if len(R) <= KEEP else {}
    tests = []
    for k in range(D):
        ok = False
        for (X, V) in R:
            p = len(X)
            if v2(p) != V:
                continue
            sa = sum(incs([x[0] for x in X], [x[1] for x in X])) % 2
            sb = sum(incs([x[2] for x in X], [x[3] for x in X])) % 2
            if sa != sb:
                ok = True
                break
        tests.append(ok)
        if verbose:
            print(f"k={k}: |R_k|={len(R)}, test(R_k)={ok}", flush=True)
        new = set()
        for (X, V) in R:
            P = len(X)
            la = lift1([x[0] for x in X], [x[1] for x in X])
            lb = lift1([x[2] for x in X], [x[3] for x in X])
            for za in la:
                for zb in lb:
                    Q = max(len(za), len(zb))
                    if Q > Tmax:
                        continue
                    ser = [(X[t % P][1], za[t % len(za)], X[t % P][3], zb[t % len(zb)]) for t in range(Q)]
                    new.add((canon(ser), max(V, v2(minper(za)), v2(minper(zb)))))
        R = frozenset(new)
        h = digest(R)
        if h in seen:
            k0 = seen[h]
            how = "as sets" if (k0 in kept and kept[k0] == R) else "by SHA-256 digest"
            if verbose:
                print(f"R_{k+1} = R_{k0} ({how}): k0={k0}, rho={k+1-k0}")
            return tests, (k0, k + 1 - k0), k + 1
        seen[h] = k + 1
        if len(R) <= KEEP:
            kept[k + 1] = R
    if verbose:
        print(f"|R_{D}|={len(R)}")
    return tests, None, D


if __name__ == '__main__':
    rule, nmax, Tmax, D = map(int, sys.argv[1:5])
    tests, rep, last = run(rule, nmax, Tmax, D)
    print(f"rule {rule}, n_max={nmax}, T_max={Tmax}: test(R_k) holds for all 0 <= k <= {len(tests)-1}: {all(tests)}")
    if rep:
        print(f"repetition R_{rep[0]+rep[1]} = R_{rep[0]}; the sets R_0, ..., R_{last-1} are pairwise distinct")
    else:
        print(f"no repetition: the sets R_0, ..., R_{last} are pairwise distinct")
