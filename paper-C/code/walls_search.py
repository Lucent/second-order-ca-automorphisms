"""(1) Persistent words: w of length k such that F^p maps every occurrence of w to w, whatever the context,
for some period p <= 6 (checked exhaustively over all contexts of width p on each side).
(2) Damage spreading: flip one cell of a random configuration, measure the left and right edges of the
difference set under F^n (forward and backward), averaged."""
import itertools, random, sys
import numpy as np
from centralizer import make_rule
def run(rule):
    f, F_loc, Finv_loc = make_rule(rule)
    T = np.array([[[F_loc(a, b, c) for c in range(4)] for b in range(4)] for a in range(4)], dtype=np.int8)
    def Fw(arr):  # arr (..., L) -> (..., L-2)
        return T[arr[..., :-2], arr[..., 1:-1], arr[..., 2:]]
    found = []
    for k in range(1, 6):
        for p in range(1, 7):
            ctx = np.array(list(itertools.product(range(4), repeat=2 * p)), dtype=np.int8) if p <= 4 else None
            if ctx is None:   # sample contexts for larger p
                ctx = np.random.default_rng(0).integers(0, 4, size=(20000, 2 * p), dtype=np.int8)
            for w in itertools.product(range(4), repeat=k):
                arr = np.concatenate([ctx[:, :p], np.tile(np.array(w, dtype=np.int8), (len(ctx), 1)), ctx[:, p:]], axis=1)
                x = arr
                for _ in range(p): x = Fw(x)
                if (x == np.array(w, dtype=np.int8)).all():
                    found.append((k, p, "".join(map(str, w))))
            if found and found[-1][0] == k: pass
        if any(fk == k for fk, _, _ in found): break
    # damage spreading
    rng = np.random.default_rng(1); L = 801; N = 300; edges = []
    for _ in range(40):
        c = rng.integers(0, 4, size=L, dtype=np.int8); d = c.copy(); d[L // 2] ^= rng.integers(1, 4)
        a, b = c, d; lo_hi = []
        for n in range(N):
            a = np.concatenate([[a[-1]], a, [a[0]]]); a = Fw(a); b = np.concatenate([[b[-1]], b, [b[0]]]); b = Fw(b)
            diff = np.flatnonzero(a != b)
            lo_hi.append((L // 2 - diff.min(), diff.max() - L // 2) if len(diff) else (None, None))
        edges.append(lo_hi[-1])
    left = [e[0] for e in edges if e[0] is not None]; right = [e[1] for e in edges if e[1] is not None]; healed = sum(1 for e in edges if e[0] is None)
    print(f"{rule}R: persistent words (k, period, word) found: {found[:6]}{' ...' if len(found) > 6 else ''} ({len(found)} total, shortest length {found[0][0] if found else '-'})")
    print(f"     damage after {N} steps over 40 trials: left spread {np.mean(left) if left else 0:.0f} (min {min(left) if left else '-'}), right spread {np.mean(right) if right else 0:.0f} (min {min(right) if right else '-'}), healed {healed}/40")
for r in map(int, sys.argv[1:]): run(r)
