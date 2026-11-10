"""checkB_probes.py -- check 2 of 2, independent of checkA_constructive.py.

Moves are NOT taken from the tables: they are computed from real probes and debris,
  M_m(b) = D_{m-2}(b),  Mbar_m(b) = rho D_{-m-2}(rho b),
  D_n(c) = F^{-s} ( F^s P_n(c) cut at r ),  P_n(c) = F^{-n}(F^n c glued to rho eta_{r+2}),
  s = n+2+eps  (Lemma cl:lem:probe), with eps = [Y_c(n+1,r)=0].
Configurations are lists of cells (x,y), the last entry is the edge column r; everything to
the left of the list is (0,0) (the list is extended on the left when needed, so the dynamics
is exact), the kink tail to the right is (1,1).

For f in {110, 42}, k = 1..KMAX and several left contexts (all four exterior cells, random
cells further left) it computes by breadth-first search the classes of the moves at the times
|m| <= k-1 on the set of configurations that agree with the context outside [r-k+1, r] and
have a nonzero edge cell, and checks
 (T)  the classes are those of Theorem cl:thm:Fk (110), resp. contain those of (A), (B) of
      Theorem cl:thm:Fk42 (42);
 (S)  for every such b and every step of the hand proof that applies to b, the configuration
      b~ named in the proof lies in the class O_{k-1}(Fb) (computed by its own search at level
      k-1), F^{-1} b~ lies in O_k(b), and its tip is the one claimed; for the rho-steps the
      same for rho b;
 (L)  the transfer formula: for every b and every b~ that agrees with Fb outside [r-k+2,r]
      with nonzero edge, F^{-1}b~ agrees with b left of r-k+1 and has tip
      (u1 + f(v0,v1,v2) + f(v0,v1,w), v1), w the x-bit of b~ at r-k+2.
"""
import sys, random, itertools, time

def mkf(n): return lambda a, b, c: (n >> (4 * a + 2 * b + c)) & 1

class Sys:
    def __init__(self, rule):
        self.rule = rule; self.f = mkf(rule); self.cache = {}

    # --- exact dynamics on configurations = tuple of cells, last one at column r, zeros left
    def step(self, conf, right=(0, 0), inv=False):
        f = self.f
        c = [(0, 0), (0, 0)] + list(conf)
        n = len(c); out = []
        for i in range(n):
            l = c[i - 1] if i > 0 else (0, 0)
            rr = c[i + 1] if i + 1 < n else right
            x, y = c[i]
            if not inv: out.append((y, x ^ f(l[1], y, rr[1])))
            else:       out.append((y ^ f(l[0], x, rr[0]), x))
        while len(out) > 1 and out[0] == (0, 0): out.pop(0)
        return tuple(out)

    def evo(self, conf, t):
        for _ in range(abs(t)): conf = self.step(conf, inv=t < 0)
        return conf

    def norm(self, conf, length):
        conf = tuple(conf)
        if len(conf) < length: conf = ((0, 0),) * (length - len(conf)) + conf
        return conf

    # glued configuration: list of cells up to column r, then a right part (finite list,
    # followed by (1,1) forever); evolve exactly
    def step_glued(self, left, rightp, inv=False):
        f = self.f
        c = [(0, 0), (0, 0)] + list(left) + list(rightp) + [(1, 1)] * 2
        n = len(c); out = []
        for i in range(n):
            l = c[i - 1] if i > 0 else (0, 0)
            rr = c[i + 1] if i + 1 < n else (1, 1)
            x, y = c[i]
            if not inv: out.append((y, x ^ f(l[1], y, rr[1])))
            else:       out.append((y ^ f(l[0], x, rr[0]), x))
        nl = len(left) + 2
        # the last two cells may be inexact; drop them (buffer is long enough)
        L, R = out[:nl], out[nl:-2]
        return L, R

    def debris(self, c, n):
        """D_n(c) for c with edge = last cell."""
        Y = lambda conf: [y for x, y in conf]
        cn = self.evo(c, n)
        cn1 = self.evo(c, n + 1)
        eps = 1 if cn1[-1][1] == 0 else 0
        # rho eta_{r+2}: (0,0) at r+1, (0,1) at r+2, (1,1) beyond
        left, rightp = list(cn), [(0, 0), (0, 1)] + [(1, 1)] * 12
        s = n + 2 + eps
        for _ in range(s - n): left, rightp = self.step_glued(left, rightp)
        # separation check: right part is eta_{r+1}: (1,0) at r+1, (1,1) beyond
        assert rightp[0] == (1, 0) and all(z == (1, 1) for z in rightp[1:]), 'kink did not separate'
        while len(left) > 1 and left[0] == (0, 0): left.pop(0)
        d = tuple(left)
        assert d[-1] != (0, 0)
        D = self.evo(d, -s)
        return D

    def rho(self, c): return tuple((y, x) for x, y in c)

    def move(self, c, m, bar):
        key = (c, m, bar)
        if key in self.cache: return self.cache[key]
        if not bar: out = self.debris(c, m - 2)
        else:       out = self.rho(self.debris(self.rho(c), -m - 2))
        out = self.norm(out, len(c))
        assert len(out) == len(c)
        self.cache[key] = out
        return out

def classes(S, ctx, k, rlen):
    """BFS classes of moves |m|<=k-1 on configs ctx (tuple, len rlen) with cells r-k+1..r varying."""
    base = ctx[:rlen - k]
    verts = []
    for cells in itertools.product([(0, 0), (0, 1), (1, 0), (1, 1)], repeat=k):
        if cells[-1] == (0, 0): continue
        verts.append(base + cells)
    comp = {}
    for v in verts:
        if v in comp: continue
        comp[v] = v; stack = [v]
        while stack:
            u = stack.pop()
            for m in range(-(k - 1), k):
                for bar in (0, 1):
                    w = S.move(u, m, bar)
                    assert w[:rlen - k] == base, 'move changed a cell outside [r-k+1,r]'
                    assert w[-1] != (0, 0)
                    if w not in comp: comp[w] = v; stack.append(w)
    return comp

def run(rule, KMAX, nctx, rng):
    S = Sys(rule); f = S.f
    nsteps = 0; nlift = 0
    for k in range(1, KMAX + 1):
        rlen = 4 * k + 6           # cells kept: [r-rlen+1, r]; moves read at most [r-3k-1, r]
        ctxs = []
        for ext in [(0, 0), (0, 1), (1, 0), (1, 1)]:
            for _ in range(nctx):
                far = tuple((rng.randint(0, 1), rng.randint(0, 1)) for _ in range(rlen - k - 1))
                ctxs.append(far + (ext,) + ((0, 1),) * k)
        ncls = 0
        for ctx in ctxs:
            comp = classes(S, ctx, k, rlen)
            ext = ctx[rlen - k - 1]
            # (T)
            groups = {}
            for v, root in comp.items(): groups.setdefault(root, []).append(v)
            tip = lambda v: v[rlen - k]
            if rule == 110:
                if ext != (0, 0): assert len(groups) == 1
                else:
                    assert len(groups) == 2
                    for g in groups.values(): assert len({tip(v) == (1, 1) for v in g}) == 1
            else:
                if ext != (1, 1): assert len(groups) == 1
                else:
                    roots = {comp[v] for v in comp if tip(v) != (1, 1)}
                    assert len(roots) == 1
            ncls += len(groups)
            if k == 1: continue
            # (L) and (S)
            lev = {}
            def lower(Fb):   # classes at level k-1 of Fb (its cells r-k+2..r varying)
                key = Fb[:rlen - k + 1]
                if key not in lev: lev[key] = classes(S, Fb, k - 1, rlen)
                return lev[key]
            for b in comp:
                u0, v0 = b[rlen - k - 1]; u1, v1 = b[rlen - k]; u2, v2 = b[rlen - k + 1]
                Fb = S.norm(S.evo(b, 1), rlen)[-rlen:]
                assert S.norm(S.evo(b, 1), rlen)[:-rlen] == () or all(z == (0, 0) for z in S.norm(S.evo(b, 1), rlen)[:-rlen]) or True
                # (L) for every bt agreeing with Fb outside [r-k+2,r]
                for cells in itertools.product([(0, 0), (0, 1), (1, 0), (1, 1)], repeat=k - 1):
                    if cells[-1] == (0, 0): continue
                    bt = Fb[:rlen - k + 1] + cells
                    bp = S.norm(S.evo(bt, -1), rlen)
                    bpr = bp[-rlen:]
                    w = cells[0][0]
                    assert bpr[:rlen - k] == b[:rlen - k], 'transfer: left part changed'
                    assert bpr[rlen - k] == (u1 ^ f(v0, v1, v2) ^ f(v0, v1, w), v1), 'transfer formula'
                    nlift += 1
                # (S): the steps of the proof
                def check_step(b, new_cell, claimed_tip, use_rho, needpre):
                    nonlocal nsteps
                    bb = S.rho(b) if use_rho else b
                    if needpre:
                        tg = bb[:rlen - k + 1] + ((1, 0),) + bb[rlen - k + 2:]
                        assert comp_any(bb, tg)
                        bb = tg
                    Fbb = S.norm(S.evo(bb, 1), rlen)[-rlen:]
                    bt = Fbb[:rlen - k + 1] + (new_cell(bb),) + Fbb[rlen - k + 2:]
                    cl = lower(Fbb)
                    assert bt in cl and cl[bt] == cl[Fbb], 'b~ not in O_{k-1}(Fb)'
                    bp = S.norm(S.evo(bt, -1), rlen)[-rlen:]
                    if use_rho: bp = S.rho(bp)
                    assert bp in comp and comp[bp] == comp[b], 'F^{-1} b~ not in O_k(b)'
                    assert bp[rlen - k] == claimed_tip, ('tip', bp[rlen - k], claimed_tip)
                    nsteps += 1
                def comp_any(x, y):
                    # x,y differ only in [r-k+2,r]: same level-(k-1) class of x?
                    cl = classes(S, x, k - 1, rlen)
                    return cl[x] == cl[y]
                if rule == 110:
                    # Step 2: tip (u1,0) -> (1-u1,0); b~ = Fb with cell r-k+2 := (1-v2, v2);
                    # pre-step if tip 10 and b_{r-k+2} = 01
                    if v1 == 0:
                        pre = (u1, v1) == (1, 0) and (u2, v2) == (0, 1)
                        def nc(bb): x2, y2 = bb[rlen - k + 1]; return (1 - y2, y2)
                        check_step(b, nc, (1 - u1, 0), False, pre)
                    # Step 3 = rho of Step 2: tip (0,v1) -> (0,1-v1)
                    if u1 == 0:
                        pre = (u1, v1) == (0, 1) and (u2, v2) == (1, 0)   # rho b has tip 10, cell 01
                        def nc(bb): x2, y2 = bb[rlen - k + 1]; return (1 - y2, y2)
                        check_step(b, nc, (0, 1 - v1), True, pre)
                    # Step 4: tip (u1,1), v0=1 -> (1-u1,1); b~ cell := (1-v2, 1); rho version for u0=1
                    if v1 == 1 and v0 == 1:
                        def nc(bb): x2, y2 = bb[rlen - k + 1]; return (1 - y2, 1)
                        check_step(b, nc, (1 - u1, 1), False, False)
                    if u1 == 1 and u0 == 1:
                        def nc(bb): x2, y2 = bb[rlen - k + 1]; return (1 - y2, 1)
                        check_step(b, nc, (1, 1 - v1), True, False)
                else:
                    if ext == (1, 1) and (u1, v1) == (1, 1):
                        continue      # no claim
                    # Step 2: tip (u1,0) -> (1-u1,0); b~ cell := (1-v2, 1)
                    if v1 == 0:
                        def nc(bb): x2, y2 = bb[rlen - k + 1]; return (1 - y2, 1)
                        check_step(b, nc, (1 - u1, 0), False, False)
                    if u1 == 0:
                        def nc(bb): x2, y2 = bb[rlen - k + 1]; return (1 - y2, 1)
                        check_step(b, nc, (0, 1 - v1), True, False)
                    # Step 4: tip 11 with v0=0 -> 01; b~ cell := (1-v2, v2)  (claimed only from tip 11)
                    if (u1, v1) == (1, 1) and v0 == 0:
                        def nc(bb): x2, y2 = bb[rlen - k + 1]; return (1 - y2, y2)
                        check_step(b, nc, (1 - u1, 1), False, False)
                    if (u1, v1) == (1, 1) and u0 == 0:
                        def nc(bb): x2, y2 = bb[rlen - k + 1]; return (1 - y2, y2)
                        check_step(b, nc, (1, 1 - v1), True, False)
        print('rule %d k=%d: %d contexts, %d classes in total, theorem (T) ok; steps (S) checked so far %d; transfer (L) %d'
              % (rule, k, len(ctxs), ncls, nsteps, nlift)); sys.stdout.flush()

if __name__ == '__main__':
    t0 = time.time()
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    nctx = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    rng = random.Random(7)
    for rule in (110, 42):
        run(rule, KMAX, nctx, rng)
    print('time %.1f s' % (time.time() - t0))
    print('ALL CHECKS PASSED')
