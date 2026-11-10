"""The certificate (T_P) for second-order automata of radius r over {0,1},
computed from its definition, for many rules g at once.

Rule: f(y_{i-r},...,y_{i+r}) = y_{i-r} XOR g(xi_0,...,xi_{2r-1}), xi_k = y_{i-r+1+k}.
A function g is given by its truth table, an integer `tab` with
    g(xi_0,...,xi_{2r-1}) = bit number (xi_0 + 2 xi_1 + 4 xi_2 + ...) of tab.

Definition of (T_P).  A column state is a 2r-tuple s = (s_1,...,s_2r) of
functions Z/P -> {0,1}.
    phi(s)_t = s_1(t-1) XOR g(s_2r(t), s_{2r-1}(t), ..., s_1(t)),
    pi(s)    = sum_t phi(s)_t mod 2.
If pi(s) = 0 the successors of s are the two states (s_2,...,s_2r,z) with
z(t+1) = z(t) XOR phi(s)_t; if pi(s) = 1 there is none.  Im_0 is the set of all
states, Im_{m+1} the set of successors of elements of Im_m, S_inf their
intersection.  (T_P) holds iff some s in S_inf has pi(s) = 1.

Data representation.  A state is stored time-major: with n = 2r,
    code(s) = sum_t c_t * 2^(n t),   c_t = sum_{j=1..n} s_j(t) * 2^(j-1),
so c_t is the snapshot of all n layers at time t.  All rules of a batch are
processed together as rows of numpy arrays; sets of states are boolean vectors.

Run as a script (radius 2, all 65536 functions g):
    python tp_direct.py            writes out_leastP_r2.txt and prints a summary
"""
import sys
import numpy as np


def tables_from_ints(tabs, n):
    """array B x 2^n of truth-table bits, T[b, x] = bit x of tabs[b]."""
    return np.array([[(int(tab) >> x) & 1 for x in range(1 << n)] for tab in tabs], dtype=np.uint8)


class Columns:
    """state tables for n = 2r layers and period P (independent of g)."""

    def __init__(self, r, P):
        n = 2 * r
        self.r, self.n, self.P = r, n, P
        self.NS = 1 << (n * P)
        code = np.arange(self.NS, dtype=np.int64)
        snap = [(code >> (n * t)) & ((1 << n) - 1) for t in range(P)]
        # argument of g at time t: xi_k = s_{n-k}(t), so the bit index
        # sum_k xi_k 2^k is the n-bit reversal of the snapshot c_t.
        rev = np.array([sum(((v >> j) & 1) << (n - 1 - j) for j in range(n)) for v in range(1 << n)])
        self.arg = np.stack([rev[snap[t]] for t in range(P)], axis=1)            # NS x P
        self.old = np.stack([snap[(t - 1) % P] & 1 for t in range(P)], axis=1).astype(np.uint8)  # s_1(t-1)
        self.dropped = sum((snap[t] >> 1) << (n * t) for t in range(P))          # (s_2,...,s_n, 0)
        self.top = [1 << (n * t + n - 1) for t in range(P)]                      # position of z(t)

    def phi(self, T):
        """B x NS x P array of phi(s)_t for the truth tables T (B x 2^n)."""
        return self.old[None, :, :] ^ T[:, self.arg]

    def successors(self, ph):
        """the two successor codes (z(0) = 0 and z(0) = 1); valid where pi = 0."""
        P = self.P
        z = np.zeros(ph.shape[:2], dtype=np.uint8)
        s0 = np.broadcast_to(self.dropped, ph.shape[:2]).copy()
        for t in range(1, P):
            z = z ^ ph[:, :, t - 1]                       # z(t) = z(t-1) XOR phi_{t-1}
            s0 += z.astype(np.int64) * self.top[t]
        return s0, s0 ^ sum(self.top)

    def s_inf(self, T):
        """returns (pi, S_inf, list of |Im_m|) ; pi and S_inf are B x NS arrays."""
        ph = self.phi(T)
        pi = (ph.sum(axis=2) & 1).astype(bool)
        s0, s1 = self.successors(ph)
        B, NS = pi.shape
        row = np.arange(B, dtype=np.int64)[:, None] * NS
        f0 = (row + s0).ravel(); f1 = (row + s1).ravel()
        has_succ = (~pi).ravel()
        cur = np.ones(B * NS, dtype=bool)
        sizes = [cur.reshape(B, NS).sum(axis=1)]
        while True:
            src = cur & has_succ
            nxt = np.zeros(B * NS, dtype=bool)
            nxt[f0[src]] = True
            nxt[f1[src]] = True
            if np.array_equal(nxt, cur):
                break
            cur = nxt
            sizes.append(cur.reshape(B, NS).sum(axis=1))
        return pi, cur.reshape(B, NS), sizes

    def holds(self, T):
        """boolean vector: (T_P) for each row of T."""
        pi, S, _ = self.s_inf(T)
        return (pi & S).any(axis=1)

    def decode(self, code):
        """state code -> tuple (s_1,...,s_n), each s_j a tuple (s_j(0),...,s_j(P-1))."""
        n, P = self.n, self.P
        return tuple(tuple((code >> (n * t + j)) & 1 for t in range(P)) for j in range(n))


def least_P(tabs, r, Pmax, batch=2048):
    """dict tab -> least P <= Pmax with (T_P), or None."""
    n = 2 * r
    res = {int(t): None for t in tabs}
    todo = [int(t) for t in tabs]
    for P in range(1, Pmax + 1):
        if not todo:
            break
        C = Columns(r, P)
        b = max(1, min(batch, (1 << 24) // (C.NS * P)))
        rest = []
        for i in range(0, len(todo), b):
            chunk = todo[i:i + b]
            ok = C.holds(tables_from_ints(chunk, n))
            for tab, o in zip(chunk, ok):
                if o:
                    res[tab] = P
                else:
                    rest.append(tab)
        todo = rest
    return res


# ---- algebraic normal form and the conditions (N), (D)

def anf(tab, n):
    """set of monomials; a monomial is the bit mask of the variables xi_k it contains."""
    a = [(tab >> x) & 1 for x in range(1 << n)]
    for i in range(n):
        for x in range(1 << n):
            if x >> i & 1:
                a[x] ^= a[x ^ (1 << i)]
    return {x for x in range(1 << n) if a[x]}


def anf_str(tab, n):
    ms = sorted(anf(tab, n), key=lambda m: (bin(m).count("1"), [k for k in range(n) if m >> k & 1]))
    if not ms:
        return "0"
    return " + ".join("1" if m == 0 else "".join("xi%d" % k for k in range(n) if m >> k & 1) for m in ms)


def is_affine(tab, n):
    return all(bin(m).count("1") <= 1 for m in anf(tab, n))


def J(tab, n):
    """max j in [1,n-1] with d_0 d_j g not identically 0 (None if d_0 g is constant).
    d_0 d_j g is not 0 iff some monomial of g contains both xi_0 and xi_j."""
    js = [j for j in range(1, n) if any((m & 1) and (m >> j & 1) for m in anf(tab, n))]
    return max(js) if js else None


def self_dual(tab, n):
    full = (1 << n) - 1
    return all(((tab >> x) ^ (tab >> (x ^ full))) & 1 for x in range(1 << n))


if __name__ == "__main__":
    r, n, Pmax = 2, 4, 4
    tabs = list(range(1 << (1 << n)))
    lp = least_P(tabs, r, Pmax)
    with open("out_leastP_r2.txt", "w") as fh:
        fh.write("# truth table of g (bit index xi0+2xi1+4xi2+8xi3), least P<=4 with (T_P) ('-' if none)\n")
        for tab in tabs:
            fh.write("0x%04x %s\n" % (tab, lp[tab] if lp[tab] else "-"))
    counts = {}
    nonaff = cN = cND = sd_nonaff = 0
    for tab in tabs:
        aff = is_affine(tab, n)
        key = ("affine" if aff else "non-affine", lp[tab])
        counts[key] = counts.get(key, 0) + 1
        if aff:
            continue
        nonaff += 1
        j = J(tab, n); sd = self_dual(tab, n)
        sd_nonaff += sd
        if j is not None and j >= r:
            cN += 1
            cND += not sd
    print("radius 2, all %d functions g; least P <= %d with (T_P):" % (len(tabs), Pmax))
    for key in sorted(counts, key=lambda k: (k[0], k[1] or 99)):
        print("  %-10s least P = %-4s : %d" % (key[0], key[1], counts[key]))
    print("non-affine g:", nonaff)
    print("non-affine g with (N):", cN)
    print("non-affine g with (N) and (D):", cND)
    print("non-affine self-dual g:", sd_nonaff)
    print("affine g with no (T_P), P<=%d:" % Pmax, ["0x%04x" % t for t in tabs if is_affine(t, n) and lp[t] is None])
    four = ["0x%04x" % t for t in tabs if lp[t] == 4]
    print("g with least P = 4 (%d):" % len(four), " ".join(four))
