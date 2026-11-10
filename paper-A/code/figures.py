"""Spacetime figures for the paper, computed from exact simulations.

Run from this folder (needs numpy and matplotlib):
    python3 figures.py
Writes the figures to ../tex/figs/*.pdf (to ../figs/*.pdf if ../tex/figs does
not exist) and the log out_figures.txt.

Conventions: a configuration c has cells c_i=(x_i,y_i); its diagram Y_c has
row -1 = x and row 0 = y, and Y(t+1,i) = Y(t-1,i) xor f(Y(t,i-1),Y(t,i),Y(t,i+1)).
Configurations are eventually constant (0,0) or (1,1) on each side; these
constants are fixed by F for the rules used here.  Simulations run on a padded
window whose boundary cells are held at the constant; the padding exceeds the
number of time steps, so every displayed cell is exact.
"""
import os, random
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.colors import ListedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "..", "tex", "figs")
if not os.path.isdir(FIGS):
    FIGS = os.path.join(HERE, "..", "figs")
os.makedirs(FIGS, exist_ok=True)
LOG = []
def log(s):
    print(s); LOG.append(s)

ACC = "#d95f02"          # accent colour
ACC2 = "#1b9e77"

def rule(n):
    return lambda a, b, c: (n >> (4 * a + 2 * b + c)) & 1

class Diagram:
    """Exact diagram of a configuration that is constant Lc left of lo and Rc right of hi."""
    def __init__(self, n, cells, lo, Lc, Rc, t0=0, tmin=-5, tmax=30, pad=None):
        # cells: list of (x,y) for columns lo..lo+len-1, the state at time t0
        f = rule(n)
        T = max(abs(tmax - t0), abs(tmin - t0)) + 3
        pad = pad or T + 4
        self.lo = lo - pad
        W = len(cells) + 2 * pad
        xs = [Lc] * pad + [c[0] for c in cells] + [Rc] * pad
        ys = [Lc] * pad + [c[1] for c in cells] + [Rc] * pad
        rows = {t0 - 1: np.array(xs, np.int8), t0: np.array(ys, np.int8)}
        def nxt(prev, cur):
            out = prev.copy()
            for i in range(W):
                a = cur[i - 1] if i > 0 else (Lc)
                c = cur[i + 1] if i < W - 1 else (Rc)
                out[i] = prev[i] ^ f(a, cur[i], c)
            return out
        for t in range(t0, tmax):
            rows[t + 1] = nxt(rows[t - 1], rows[t])
        for t in range(t0 - 1, tmin, -1):
            rows[t - 1] = nxt(rows[t + 1], rows[t])
        self.rows = rows
        self.Lc, self.Rc = Lc, Rc
        self.tmin, self.tmax = tmin, tmax
    def __call__(self, t, i):
        j = i - self.lo
        row = self.rows[t]
        if j < 0: return self.Lc
        if j >= len(row): return self.Rc
        return int(row[j])
    def state(self, t, a, b):
        return [(self(t - 1, i), self(t, i)) for i in range(a, b + 1)]
    def array(self, t1, t2, a, b):
        return np.array([[self(t, i) for i in range(a, b + 1)] for t in range(t1, t2 + 1)])

def show(ax, A, t1, a, title=None, cmap=ListedColormap(["white", "black"])):
    T, C = A.shape
    ax.imshow(A, origin="lower", cmap=cmap, vmin=0, vmax=1, interpolation="nearest",
              extent=(a - 0.5, a + C - 0.5, t1 - 0.5, t1 + T - 0.5), aspect="equal")
    ax.set_xlabel("$i$", fontsize=8); ax.set_ylabel("$t$", fontsize=8)
    ax.tick_params(labelsize=6)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    if title: ax.set_title(title, fontsize=9)

def overlay(ax, mask, t1, a, color=ACC, alpha=0.85):
    T, C = mask.shape
    rgba = np.zeros((T, C, 4))
    rgba[mask] = matplotlib.colors.to_rgba(color, alpha)
    ax.imshow(rgba, origin="lower", interpolation="nearest",
              extent=(a - 0.5, a + C - 0.5, t1 - 0.5, t1 + T - 0.5), aspect="equal")

def edge_line(ax, r, t1, t2, label="$r$"):
    ax.axvline(r + 0.5, color=ACC2, lw=0.8, ls="--")
    ax.text(r + 0.5, t1 - 0.7, label, color=ACC2, fontsize=7, ha="center", va="top",
            bbox=dict(fc="white", ec="none", pad=0.5))

def save(fig, name):
    fig.savefig(os.path.join(FIGS, name), bbox_inches="tight")
    plt.close(fig)
    log("wrote figs/" + name)

def last_zero(D, t, a, b):
    z = [i for i in range(a, b + 1) if D(t, i) == 0]
    return max(z) if z else None

def zero_tailed(rng, lo, r):
    cells = [(rng.randint(0, 1), rng.randint(0, 1)) for _ in range(lo, r + 1)]
    while cells[-1] == (0, 0):
        cells[-1] = (rng.randint(0, 1), rng.randint(0, 1))
    return cells

# ---------------------------------------------------------------- (1) front
def fig_front():
    rng = random.Random(3)
    for trial in range(200):
        L = 14
        cells = [(rng.randint(0, 1), rng.randint(0, 1)) for _ in range(L)]
        D = Diagram(110, cells, 0, 0, 1, tmin=-1, tmax=34)
        # find a front at (t0,q): Y(t0,q)=0, Y(t0,i)=1 for i>q, Y(t0-1,i)=1 for i>=q
        found = None
        for t0 in range(0, 12):
            q = last_zero(D, t0, -40, 80)
            if q is None: continue
            if all(D(t0 - 1, i) == 1 for i in range(q, 80)):
                found = (t0, q); break
        if found and found[0] >= 4: break
    t0, q = found
    for t in range(t0, 34):
        s = q + t - t0
        assert all(D(t, i) == 0 for i in range(q, s + 1)), t
        assert all(D(t, i) == 1 for i in range(s + 1, 80)), t
    log(f"front: rule 110, front at (t0,q)=({t0},{q}); zero wake and speed-1 front verified for t<=33; U={q-t0}")
    a, b, t1, t2 = -6, 45, -1, 33
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    show(ax, D.array(t1, t2, a, b), t1, a)
    ax.plot([q - 0.5, q + (t2 - t0) + 0.5], [t0, t2 + 0.5], color=ACC, lw=1.2)
    ax.plot([q], [t0], "o", color=ACC, ms=3)
    ax.text(q - 1, t0, "$(t_0,q)$", color=ACC, fontsize=7, ha="right", va="center",
            bbox=dict(fc="white", ec="none", pad=0.6))
    save(fig, "fig_front.pdf")

# ---------------------------------------------------------------- (2) bounce
def fig_bounce():
    rng = random.Random(11)
    e, N = 0, 14
    d = zero_tailed(rng, -10, e)
    Dd = Diagram(110, d, -10, 0, 0, tmin=-1, tmax=40)
    cells = d + [(0, 0)] * (N - e - 1) + [(0, 1)] + [(1, 1)] * 30
    Dc = Diagram(110, cells, -10, 0, 1, tmin=-1, tmax=40)
    phi = Dd(N - e - 1, e)
    U = None
    lam = [last_zero(Dc, t, -60, 120) - t for t in range(25, 40)]
    assert len(set(lam)) == 1
    U = lam[0]
    assert U == 2 * e - N + phi, (U, e, N, phi)
    # incoming kink rows 0..N-e-1: rightmost zero at N-1-t
    for t in range(0, N - e - 1):
        assert last_zero(Dc, t, -60, 120) == N - 1 - t
    log(f"bounce: rule 110, e={e}, N={N}, phi={phi}, U(c)={U} = 2e-N+phi verified")
    for seed in range(40):
        r2 = random.Random(seed)
        for N2 in range(1, 12):
            d2 = zero_tailed(r2, -8, 0)
            D2 = Diagram(110, d2, -8, 0, 0, tmin=-1, tmax=N2 + 4)
            c2 = d2 + [(0, 0)] * (N2 - 1) + [(0, 1)] + [(1, 1)] * 20
            C2 = Diagram(110, c2, -8, 0, 1, tmin=-1, tmax=N2 + 14)
            ph = D2(N2 - 1, 0)
            l = [last_zero(C2, t, -40, 60) - t for t in range(N2 + 4, N2 + 14)]
            assert len(set(l)) == 1 and l[0] == -N2 + ph
    log("bounce: formula U=2e-N+phi also verified on 440 random cases (N-e=1..11)")
    a, b, t1, t2 = -10, 24, -1, 30
    fig, ax = plt.subplots(figsize=(4.6, 4.0))
    show(ax, Dc.array(t1, t2, a, b), t1, a)
    edge_line(ax, e, t1, t2, "$e$")
    ax.annotate("$N$", xy=(N, -1), xytext=(N, -4.5), fontsize=7, ha="center", color=ACC,
                arrowprops=dict(arrowstyle="-", color=ACC, lw=0.6), annotation_clip=False)
    ts = np.arange(30 - 6, t2 + 1)
    ax.plot([U + t for t in (N - e, t2)], [N - e, t2], color=ACC, lw=1.0, ls=":")
    ax.text(U + t2 - 1, t2 - 3, "$\\eta$", color=ACC, fontsize=8, ha="right")
    ax.text(N - 4, 3, "$\\rho\\eta_N$", color=ACC, fontsize=8)
    save(fig, "fig_bounce.pdf")

# ---------------------------------------------------------------- (3) probe and debris
def probe_debris(c, lo, r, n, tmin, tmax):
    """Return diagrams of c, P_n(c), D_n(c), eps, s."""
    Dc = Diagram(110, c, lo, 0, 0, tmin=min(tmin, n - 3), tmax=max(tmax, n + 3))
    st = Dc.state(n, lo - 20, r)   # state at time n on [lo-20, r]
    glued = st + [(0, 0), (0, 1)] + [(1, 1)] * 30
    P = Diagram(110, glued, lo - 20, 0, 1, t0=n, tmin=tmin, tmax=tmax)
    eps = 1 if Dc(n + 1, r) == 0 else 0
    s = n + 2 + eps
    stP = P.state(s, lo - 40, r + 40)
    left = stP[: (r - (lo - 40)) + 1]
    right = stP[(r - (lo - 40)) + 1:]
    Npos = r + s - n - 1 - eps
    # right part equals eta_N on columns > r
    for j, cell in enumerate(right):
        i = r + 1 + j
        exp = (0, 0) if i < Npos else ((1, 0) if i == Npos else (1, 1))
        assert cell == exp, (i, cell, exp)
    assert left[-1] != (0, 0)
    Dd = Diagram(110, left, lo - 40, 0, 0, t0=s, tmin=tmin, tmax=tmax)
    for t in range(n + 1 + eps, tmax):
        for i in range(lo - 10, r + 1):
            assert Dd(t, i) == P(t, i), (t, i)
        for i in range(r + 1, r + 10):
            assert Dd(t, i) == 0
    return Dc, P, Dd, eps, s

def fig_probe():
    rng = random.Random(5)
    lo, r = -12, 0
    for trial in range(50):
        c = zero_tailed(rng, lo, r)
        n = 2
        Dc = Diagram(110, c, lo, 0, 0, tmin=-3, tmax=8)
        if Dc(n + 1, r) == 0: break    # take a case with eps = 1
    Dc, P, Dd, eps, s = probe_debris(c, lo, r, n, -4, 20)
    cnt = 0
    for seed in range(100):
        r3 = random.Random(1000 + seed)
        c3 = zero_tailed(r3, -10, 0)
        probe_debris(c3, -10, 0, r3.randint(-4, 4), -6, 14); cnt += 1
    log(f"probe: rule 110, n={n}, eps={eps}, separation time n+2+eps={s}; "
        f"F^sP = F^sD |_| eta_(r+s-n-1-eps) and agreement of diagrams on columns <= r at rows >= n+1+eps verified; "
        f"also on {cnt} random cases")
    a, b, t1, t2 = -12, 16, -4, 20
    A = P.array(t1, t2, a, b); B = Dd.array(t1, t2, a, b)
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 3.7), sharey=True)
    show(axs[0], A, t1, a, "$P_n(c)$")
    show(axs[1], B, t1, a, "$D_n(c)$")
    diff = (A != B)
    diff[:, (r - a) + 1:] = False          # columns > r: kink versus (0,0), not marked
    overlay(axs[1], diff, t1, a, alpha=0.8)
    for ax in axs:
        edge_line(ax, r, t1, t2)
        for tt, lab in ((n, "$n$"), (s, "$n{+}2{+}\\varepsilon$")):
            ax.axhline(tt, color=ACC, lw=0.6, ls=":")
        axs[1].text(b + 0.8, n, "$n$", fontsize=7, va="center", color=ACC)
        axs[1].text(b + 0.8, s, "$n{+}2{+}\\varepsilon$", fontsize=7, va="center", color=ACC)
    axs[1].set_ylabel("")
    save(fig, "fig_probe.pdf")

# ---------------------------------------------------------------- (4) cone of a move
def fig_cone():
    rng = random.Random(8)
    lo, r, m = -16, 0, 3
    for trial in range(200):
        b = zero_tailed(rng, lo, r)
        Db, P, Dm, eps, s = probe_debris(b, lo, r, m - 2, -12, 16)
        if Dm.state(0, lo, r) != b: break
    # move table at time m
    x = [Db(m - 1, i) for i in (r - 1, r)]; y = [Db(m, i) for i in (r - 1, r)]
    nx, ny = Dm(m - 1, r), Dm(m, r)
    if (x[1], y[1]) == (0, 1): exp = (1, 0)
    elif x[1] == 1: exp = (1, y[1] ^ x[0])
    else: exp = (x[1], y[1])
    assert (nx, ny) == exp
    assert Dm(m - 1, r - 1) == x[0] and Dm(m, r - 1) == y[0]
    for i in range(lo - 10, r - 1):
        assert Dm(m - 1, i) == Db(m - 1, i) and Dm(m, i) == Db(m, i)
    bad = 0
    for t in range(-12, 16):
        for i in range(lo - 10, r + 1):
            if Dm(t, i) != Db(t, i):
                assert (t - m >= r - i) or (m - 1 - t >= r - i), (t, i)
    log(f"cone: rule 110, m={m}; edge cell at time m {(x[1],y[1])} -> {(nx,ny)} (table of the move lemma); "
        f"differences of the diagrams lie in the cone t-m>=r-i or m-1-t>=r-i")
    a, bb, t1, t2 = -16, 3, -12, 16
    A = Db.array(t1, t2, a, bb); B = Dm.array(t1, t2, a, bb)
    fig, ax = plt.subplots(figsize=(3.6, 4.6))
    show(ax, B, t1, a)
    overlay(ax, A != B, t1, a)
    edge_line(ax, r, t1, t2)
    # cone boundaries: i = r - (t-m) and i = r - (m-1-t)
    ax.plot([r + 0.5, r - (t2 - m) - 0.5], [m - 0.5, t2 + 0.5], color=ACC2, lw=0.8)
    ax.plot([r + 0.5, r - (m - 1 - t1) - 0.5], [m - 0.5, t1 - 0.5], color=ACC2, lw=0.8)
    ax.axhspan(-1.5, 0.5, color=ACC2, alpha=0.12, lw=0)
    ax.text(a - 0.6, -0.5, "$c$", fontsize=7, ha="right", va="center", color=ACC2)
    ax.axhline(m - 0.5, color=ACC, lw=0.5, ls=":")
    ax.text(r + 1.2, m - 0.5, "$m$", fontsize=7, va="center", color=ACC)
    save(fig, "fig_cone.pdf")

# ---------------------------------------------------------------- (5) three windows (schematic)
def fig_windows():
    k = 6
    fig, ax = plt.subplots(figsize=(5.6, 2.4))
    cols = {"lab": [-1, 0, 1]}
    for j, l in enumerate((-1, 0, 1)):
        y = 2 - j
        lo, hi = l - (k - 2), l + (k - 2)
        ax.add_patch(plt.Rectangle((lo - 0.45, y - 0.3), hi - lo + 0.9, 0.6,
                                   color=[ACC, "0.55", ACC2][j], alpha=0.35, lw=0))
        for m in range(lo, hi + 1):
            ax.plot(m, y, "o", color="black", ms=3)
        ax.text(-k - 1.2, y, f"$\\ell={l}$: $|m-\\ell|\\leq k-2$", fontsize=7, ha="right", va="center")
    for m in range(-(k - 1), k):
        ax.plot(m, -0.85, "s", color="black", ms=3)
    ax.text(-k - 1.2, -0.85, "$|m|\\leq k-1$", fontsize=7, ha="right", va="center")
    ax.set_yticks([]); ax.set_xticks(range(-(k - 1), k))
    ax.set_xticklabels([f"${m}$" for m in range(-(k - 1), k)], fontsize=6)
    ax.set_xlabel("time $m$ of a move ($k=6$)", fontsize=8)
    ax.set_xlim(-k - 9, k + 0.5); ax.set_ylim(-1.3, 2.5)
    for s in ("left", "right", "top"): ax.spines[s].set_visible(False)
    save(fig, "fig_windows.pdf")

# ---------------------------------------------------------------- (6) collision for 14
def fig_collision14():
    N = 12
    cells = [(1, 0)] + [(1, 1)] * (N - 1)
    D = Diagram(14, cells, 0, 0, 0, tmin=-6, tmax=2 * N + 4)
    for t in range(-6, 2 * N + 4):
        for i in range(-20, N + 20):
            assert D(t, i) == int(N - 1 - abs(t - N + 1) < i < N), (t, i)
    for N2 in range(2, 40):
        cl = [(1, 0)] + [(1, 1)] * (N2 - 1)
        D2 = Diagram(14, cl, 0, 0, 0, tmin=-3, tmax=2 * N2 + 3)
        for t in range(-3, 2 * N2 + 3):
            for i in range(-6, N2 + 6):
                assert D2(t, i) == int(N2 - 1 - abs(t - N2 + 1) < i < N2)
    log(f"collision: rule 14, diagram of E(0,N) = [N-1-|t-N+1| < i < N] verified for N=2..39")
    a, b, t1, t2 = -4, N + 4, -4, 2 * N + 2
    fig, ax = plt.subplots(figsize=(3.4, 4.4))
    show(ax, D.array(t1, t2, a, b), t1, a)
    ax.axvline(N - 0.5, color=ACC2, lw=0.8, ls="--")
    ax.text(N - 0.5, t2 + 0.8, "$N$", color=ACC2, fontsize=7, ha="center", va="bottom")
    ax.text(-2.5, 5, "$\\eta_t$", color=ACC, fontsize=9)
    ax.text(-2.5, 2 * N - 8, "$\\zeta_{2N-2-t}$", color=ACC, fontsize=9)
    ax.text(N + 1.2, N - 1, "$\\lambda_N$", color=ACC, fontsize=8, va="center")
    save(fig, "fig_collision14.pdf")

if __name__ == "__main__":
    for fn in (fig_front, fig_bounce, fig_probe, fig_cone, fig_windows, fig_collision14):
        fn()
    log("all figure checks passed")
    with open(os.path.join(HERE, "out_figures.txt"), "w") as fh:
        fh.write("\n".join(LOG) + "\n")
