"""Identification of the solutions of the radius-4 census for the rules 24, 189, 63 and 187.

Input: the file raw_solutions_<rule>_R4.txt written by enum_r4.py (one table per line, 4^9 digits).

Rules 24 and 189 (kappa f = f): every solution is a constant, a map sigma^a F^b with |a|+|b| <= 4,
or kappa sigma^a F^b with |a|+|b| <= 4, where kappa is the cellwise complement s -> 3 - s.

Rules 63 and 187 (f ignores an outer variable): F reads the cells i-1, i (rule 63) or i, i+1
(rule 187), and so does its inverse, so sigma^a F^b reads an interval of |b|+1 cells and is a block
map of radius 4 whenever this interval lies in [-4, 4]. There are 81 such pairs (a, b). Every
solution is one of these 81 maps or a constant.

The tables of the maps sigma^a F^b are computed here directly from the local rule on all 4^9
windows; cells outside the window are filled with 0 and are never read by the entries used (this
is the statement about the interval, which the program asserts by recomputing with the outside
cells filled with 3).

Usage: python identify_r4.py RULE [directory of the raw solution files]
"""
import sys, os
import numpy as np

R = 4; n = 2 * R + 1

def fbit(rule, a, b, c): return (rule >> ((a << 2) | (b << 1) | c)) & 1

def step(rule, W, inverse):
    x = W >> 1; y = W & 1
    if not inverse:
        v = (rule >> (4 * y[:, :-2] + 2 * y[:, 1:-1] + y[:, 2:])) & 1
        return (y[:, 1:-1] << 1) | (x[:, 1:-1] ^ v)
    v = (rule >> (4 * x[:, :-2] + 2 * x[:, 1:-1] + x[:, 2:])) & 1
    return ((y[:, 1:-1] ^ v) << 1) | x[:, 1:-1]

def table(rule, a, b, fill):
    """table of sigma^a F^b on all windows [-4, 4], outside cells = fill"""
    N = 4 ** n; idx = np.arange(N, dtype=np.int64); pad = abs(b) + abs(a) + 1
    W = np.full((N, n + 2 * pad), fill, dtype=np.int64)
    for i in range(n): W[:, pad + i] = (idx >> (2 * (n - 1 - i))) & 3
    for _ in range(abs(b)): W = step(rule, W, b < 0)
    centre = pad + R - abs(b)                      # column of cell 0 after |b| steps
    return W[:, centre + a].astype(np.uint8)

def main():
    rule = int(sys.argv[1]); d = sys.argv[2] if len(sys.argv) > 2 else "."
    sols = [np.frombuffer(line.strip().encode(), dtype=np.uint8) - 48
            for line in open(os.path.join(d, f"raw_solutions_{rule}_R4.txt")) if line.strip()]
    left_blind = all(fbit(rule, 0, b, c) == fbit(rule, 1, b, c) for b in (0, 1) for c in (0, 1))
    right_blind = all(fbit(rule, a, b, 0) == fbit(rule, a, b, 1) for a in (0, 1) for b in (0, 1))
    self_dual = all(fbit(rule, a, b, c) == fbit(rule, 1 - a, 1 - b, 1 - c)
                    for a in (0, 1) for b in (0, 1) for c in (0, 1))
    named = {}
    if left_blind or right_blind:
        for b in range(-8, 9):
            for a in range(-12, 13):
                lo, hi = (a - abs(b), a) if right_blind else (a, a + abs(b))
                if lo < -R or hi > R: continue
                t0 = table(rule, a, b, 0); t3 = table(rule, a, b, 3)
                assert (t0 == t3).all(), "entry depends on cells outside the window"
                named[t0.tobytes()] = f"sigma^{a} F^{b}"
    else:
        for b in range(-R, R + 1):
            for a in range(-(R - abs(b)), R - abs(b) + 1):
                t0 = table(rule, a, b, 0); assert (t0 == table(rule, a, b, 3)).all()
                named[t0.tobytes()] = f"sigma^{a} F^{b}"
                if self_dual: named[(3 - t0).astype(np.uint8).tobytes()] = f"kappa sigma^{a} F^{b}"
    counts = {"sigma^a F^b": 0, "kappa sigma^a F^b": 0, "constant": 0, "unidentified": 0}
    seen = set()
    for t in sols:
        key = t.astype(np.uint8).tobytes()
        if key in named:
            seen.add(key); counts["kappa sigma^a F^b" if named[key].startswith("kappa") else "sigma^a F^b"] += 1
        elif (t == t[0]).all(): counts["constant"] += 1
        else: counts["unidentified"] += 1
    print(f"rule {rule}R, radius 4: {len(sols)} solutions; {len(named)} candidate maps, "
          f"{len(seen)} of them among the solutions; {counts}")

if __name__ == "__main__":
    main()
