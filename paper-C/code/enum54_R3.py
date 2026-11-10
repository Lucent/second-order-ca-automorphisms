"""Exact invertible centralizer of 54R at radius 3, in parallel.

Phase 1 (`run`): enumerate block maps g of radius 3 with
   (a) g o F = F o g                      (centralizer.encode)
   (b) Hedlund balance                    (necessary for surjectivity; minicard native cardinality)
   (c) g injective on all spatially periodic configurations of period <= P   (necessary for injectivity)
 split into 4^k cubes by fixing k table entries, one solver per cube, W worker processes.
 Every solution is written to results/cube_<id>.txt as soon as it is found.
Phase 2 (`classify`): dedupe, re-verify commutation, match against the certified-invertible candidates
 (sigma^a F^b, F^b o Phi(phi1, phi2) with phi_L in the centralizer C_L of F_L, sigma^{+-1} F^b Phi_1),
 and decide the rest with periodic witnesses (period <= 8) and then the exact pair-graph test (in parallel).
"""
import sys, os, time, json, argparse, itertools
from collections import Counter
import numpy as np
from pysat.solvers import Solver
from pysat.card import CardEnc, EncType
SOLVER = "minicard"
NOBAL1 = False
from centralizer import make_rule, encode, BlockMap, index_word, word_index, named_maps, injective_periodic, injective_pairgraph_np
from classify_r4 import check_commutes_np_R
from propsearch import build

RULE, R = 54, 3
n = 2 * R + 1; N = 4 ** n; NV = 2 * N
f, F_loc, Finv_loc = make_rule(RULE)
PI = {0: 1, 1: 2, 2: 0, 3: 3}

# ---------------------------------------------------------------- wall family (certified invertible candidates)
def F_seg(seg):
    w = (3,) + tuple(seg) + (3,); return tuple(F_loc(*w[i:i+3]) for i in range(len(seg)))
def S_L(L):
    return [s for s in itertools.product(range(4), repeat=L)
            if s[0] != 3 and s[-1] != 3 and not any(s[t] == 3 and s[t+1] == 3 for t in range(L-1))]
def centralizer_of(F_map, S):
    """all permutations of S commuting with F_map (brute force; |S| <= 9)."""
    idx = {s: i for i, s in enumerate(S)}; Fp = [idx[F_map(s)] for s in S]; out = []
    for perm in itertools.permutations(range(len(S))):
        if all(perm[Fp[i]] == Fp[perm[i]] for i in range(len(S))):
            out.append({S[i]: S[perm[i]] for i in range(len(S))})
    return out
def seg_at_centre(w):
    """w: 7 cells, centre 3. Returns (L, pos, content) if the centre lies in a walled segment of length 1 or 2."""
    if w[1] == w[2] == 3 and w[4] == w[5] == 3 and w[3] != 3: return (1, 0, (w[3],))
    if w[1] == w[2] == 3 and w[3] != 3 and w[4] != 3 and w[5] == w[6] == 3: return (2, 0, (w[3], w[4]))
    if w[0] == w[1] == 3 and w[2] != 3 and w[3] != 3 and w[4] == w[5] == 3: return (2, 1, (w[2], w[3]))
    return None
def candidates():
    names = named_maps(R, F_loc, Finv_loc)          # table -> (a, b)
    by_ab = {v: k for k, v in names.items()}
    S1, S2 = S_L(1), S_L(2)
    C1 = centralizer_of(lambda s: F_seg(s), S1); C2 = centralizer_of(lambda s: F_seg(s), S2)
    assert len(C1) == 3 and len(C2) == 162, (len(C1), len(C2))
    def FL_pow(s, b):
        for _ in range(b % 3): s = F_seg(s)      # F_1 and F_2 have order 3
        return s
    cands = {}
    for tab, (a, b) in names.items(): cands[tab] = f"sigma^{a} F^{b}"
    def build(b, phi1, phi2):
        base = by_ab[(0, b)]; tab = list(base)
        for e in range(N):
            w = index_word(e, n); si = seg_at_centre(w)
            if si is None: continue
            L, pos, s = si
            t = FL_pow((phi1 if L == 1 else phi2)[s], b)
            tab[e] = t[pos]
        return tuple(tab)
    id1 = next(p for p in C1 if all(p[s] == s for s in S1)); id2 = next(p for p in C2 if all(p[s] == s for s in S2))
    for b in range(-3, 4):
        for i1, phi1 in enumerate(C1):
            for i2, phi2 in enumerate(C2):
                if phi1 is id1 and phi2 is id2: continue
                cands.setdefault(build(b, phi1, phi2), f"F^{b} Phi(C1[{i1}],C2[{i2}])")
    # sigma^{+-1} F^b Phi_1(phi1), |b| <= 2 (radius 2 maps, shifted)
    for b in range(-2, 3):
        for i1, phi1 in enumerate(C1):
            if phi1 is id1: continue
            G = build(b, phi1, id2)
            for a in (-1, 1):
                tab = []
                for e in range(N):
                    w = index_word(e, n)
                    w2 = (w[2:7] + (0,) * 0) if a == 1 else w[0:5]
                    # radius-2 table lookup of G: embed the 5-window into a 7-window with dummy outer cells
                    tab.append(G[word_index((0,) + w2 + (0,))])
                cands.setdefault(tuple(tab), f"sigma^{a} F^{b} Phi_1(C1[{i1}])")
    return cands

# ---------------------------------------------------------------- phase 1

PMAX = 10
PAIRS = [(3, 3)]
PW = {}
def periodic_cells(p):
    idx = np.arange(4 ** p, dtype=np.int64)
    return [(idx >> (2 * (p - 1 - j))) & 3 for j in range(p)]
def periodic_collision(tab, pmax):
    """smallest p <= pmax and two distinct period-p configurations with g(c) == g(c2), or None."""
    tab = np.asarray(tab, dtype=np.int64)
    for p in range(1, pmax + 1):
        if p not in PW:
            cells = periodic_cells(p)
            PW[p] = np.stack([sum(cells[(i + d) % p] << (2 * (n - 1 - (d + R))) for d in range(-R, R + 1)) for i in range(p)], 1)
        W = PW[p]
        outs = tab[W]                                   # (4^p, p)
        code = np.zeros(4 ** p, dtype=np.int64)
        for i in range(p): code = (code << 2) | outs[:, i]
        order = np.argsort(code, kind="stable"); sc = code[order]
        dup = np.flatnonzero(sc[1:] == sc[:-1])
        if len(dup): return p, int(order[dup[0]]), int(order[dup[0] + 1])
    return None
CANDS = None

def periodic_collisions_all(tab, pmax, cap):
    """all pairs of distinct period-p configurations (p <= pmax) with equal images, up to cap per period."""
    tab = np.asarray(tab, dtype=np.int64); found = []
    for p in range(1, pmax + 1):
        if p not in PW:
            cells = periodic_cells(p)
            PW[p] = np.stack([sum(cells[(i + d) % p] << (2 * (n - 1 - (d + R))) for d in range(-R, R + 1)) for i in range(p)], 1)
        W = PW[p]; outs = tab[W]
        code = np.zeros(4 ** p, dtype=np.int64)
        for i in range(p): code = (code << 2) | outs[:, i]
        order = np.argsort(code, kind="stable"); sc = code[order]
        dup = np.flatnonzero(sc[1:] == sc[:-1])
        if len(dup):
            step = max(1, len(dup) // cap)
            found += [(p, int(order[d]), int(order[d + 1])) for d in dup[::step][:cap]]
    return found

KDIAM = 300
MEM_LIMIT_MB = 3500
CTX = np.arange(4 ** 12, dtype=np.int32); CTX_X = CTX >> 12; CTX_Y = CTX & (4 ** 6 - 1); del CTX
def m1_diamonds(tab, kmax):
    """m = 1 diamonds: contexts (x, y) with |x| = |y| = 6 and centres s != t such that the 7-letter images of
    x s y and x t y coincide. Returns up to kmax tuples (s, t, ctx), spread over the collision set."""
    tab = np.asarray(tab, dtype=np.int32); outs = []
    for s in range(4):
        w13 = (CTX_X << 14) | np.int32(s << 12) | CTX_Y; o = np.zeros(4 ** 12, dtype=np.int32)
        for j in range(7): o = (o << 2) | tab[(w13 >> (2 * (6 - j))) & (N - 1)]
        del w13; outs.append(o)
    found = []
    for s in range(4):
        for tt in range(s + 1, 4):
            eq = np.flatnonzero(outs[s] == outs[tt])
            if len(eq):
                step = max(1, len(eq) // max(1, kmax // 6))
                found += [(s, tt, int(c)) for c in eq[::step][:max(1, kmax // 6)]]
    return found
def diamond_windows(s, ctx):
    x, y = ctx >> 12, ctx & (4 ** 6 - 1); w13 = (x << 14) | (s << 12) | y
    return [(w13 >> (2 * (6 - j))) & (N - 1) for j in range(7)]

def add_exactly(s, inds, k, top):
    """sum(inds) == k; native for minicard, cardinality network CNF otherwise. Returns new top var."""
    if SOLVER == "minicard":
        s.add_atmost(inds, k); s.add_atmost([-i for i in inds], len(inds) - k); return top
    cnf = CardEnc.equals(lits=inds, bound=k, top_id=top, encoding=EncType.cardnetwrk)
    s.append_formula(cnf.clauses); return max(top, cnf.nv)
def periodic_windows(p):
    """for each period-p configuration (int 0..4^p-1, cell j = (c >> 2(p-1-j)) & 3) and position i, the window index."""
    out = np.zeros((4 ** p, p), dtype=np.int64)
    for c in range(4 ** p):
        cells = [(c >> (2 * (p - 1 - j))) & 3 for j in range(p)]
        for i in range(p):
            out[c, i] = word_index(tuple(cells[(i + d) % p] for d in range(-R, R + 1)))
    return out

def worker(args):
    global CANDS
    cube_id, cube_entries, P, outdir = args
    t0 = time.time()
    cj = json.load(open(f"{outdir}/candidates.json"))
    CANDS = dict(zip(cj["tables"], cj["labels"]))          # keyed by the table string (memory)
    S = build(RULE, R); cls = S["cls"]
    reps = {}
    for bit in range(NV): reps.setdefault(int(cls[bit]), bit + 1)
    rep_vars = sorted(reps.values())
    # resume: solutions already found for this cube are reloaded and blocked again; a finished cube is skipped
    prev = []
    if os.path.exists(f"{outdir}/cube_{cube_id:03d}.log") and "DONE" in open(f"{outdir}/cube_{cube_id:03d}.log").read():
        n_prev = sum(1 for l in open(f"{outdir}/cube_{cube_id:03d}.txt") if l.strip())
        return cube_id, n_prev, 0.0
    if os.path.exists(f"{outdir}/cube_{cube_id:03d}.txt"):
        prev = [l.split()[-1] for l in open(f"{outdir}/cube_{cube_id:03d}.txt") if l.strip()]
    out = open(f"{outdir}/cube_{cube_id:03d}.txt", "a"); log = open(f"{outdir}/cube_{cube_id:03d}.log", "a")
    def rss_mb(): return int(open("/proc/self/status").read().split("VmRSS:")[1].split()[0]) // 1024
    state = {"top": NV, "Tvar": {}, "s": None, "blocks": []}
    def T(e, e2):
        key = (min(e, e2), max(e, e2)); Tvar = state["Tvar"]; s = state["s"]
        if key in Tvar: return Tvar[key]
        state["top"] += 3; tX, tY, Tv = state["top"] - 2, state["top"] - 1, state["top"]
        Xe, Ye, Xf, Yf = 2 * key[0] + 1, 2 * key[0] + 2, 2 * key[1] + 1, 2 * key[1] + 2
        s.add_clause([-tX, Xe, Xf]); s.add_clause([-tX, -Xe, -Xf])
        s.add_clause([-tY, Ye, Yf]); s.add_clause([-tY, -Ye, -Yf])
        s.add_clause([-Tv, tX, tY]); Tvar[key] = Tv
        return Tv
    def build_solver():
        """(re)build the solver: commutation, optional balance, static periodic injectivity, cube, and the
        blocking clauses of solutions already found. Learned clauses are re-added by sync_shared from offset 0."""
        if state["s"] is not None: state["s"].delete()
        s = Solver(name=SOLVER); state["s"] = s; state["Tvar"] = {}
        nv, nc = encode(R, f, F_loc, s); top = NV
        for sym in ([] if NOBAL1 else range(4)):
            xb, yb = sym >> 1, sym & 1; inds = []
            for e in range(N):
                top += 1; I = top; X, Y = 2 * e + 1, 2 * e + 2
                lx = X if xb else -X; ly = Y if yb else -Y
                s.add_clause([-I, lx]); s.add_clause([-I, ly]); s.add_clause([I, -lx, -ly]); inds.append(I)
            top = add_exactly(s, inds, 4 ** (2 * R), top)
        state["top"] = top; ninj = 0
        for p in range(1, P + 1):
            W = periodic_windows(p); M = 4 ** p
            for c in range(M):
                for c2 in range(c + 1, M):
                    s.add_clause([T(int(W[c, i]), int(W[c2, i])) for i in range(p) if W[c, i] != W[c2, i]]); ninj += 1
        for e, val in cube_entries:
            X, Y = 2 * e + 1, 2 * e + 2
            s.add_clause([X if (val >> 1) & 1 else -X]); s.add_clause([Y if val & 1 else -Y])
        for blk in state["blocks"]: s.add_clause(blk)
        return ninj
    shared = open(f"{outdir}/learned.txt", "a"); shared_pos = [0]; nsynced = [0]
    def sync_shared():
        """add every clause any worker has learned (each is a valid necessary condition for every cube)."""
        with open(f"{outdir}/learned.txt") as fh:
            fh.seek(shared_pos[0])
            while True:
                line = fh.readline()
                if not line.endswith("\n"): break
                shared_pos[0] = fh.tell()
                prs = [tuple(int(v) for v in x.split(",")) for x in line.split()]
                state["s"].add_clause([T(e, e2) for e, e2 in prs]); nsynced[0] += 1
        return nsynced[0]
    def block_of(tab):
        """blocking clause (over the class representatives) excluding the table tab; variable v is the X bit
        (v odd) or Y bit (v even) of entry (v-1)//2."""
        lits = []
        for v in rep_vars:
            e = (v - 1) // 2; bit = (tab[e] >> 1) & 1 if (v - 1) % 2 == 0 else tab[e] & 1
            lits.append(-v if bit else v)
        return lits
    for tabs in prev:                                   # blocking clauses of solutions found before a restart
        state["blocks"].append(block_of([int(ch) for ch in tabs]))
    ninj = build_solver(); t_enc = time.time() - t0
    if prev: print(f"cube {cube_id}: resumed with {len(prev)} previously found solutions blocked", file=log, flush=True)
    print(f"cube {cube_id}: encoded {state['top']} vars, {ninj} static injectivity clauses, {len(state['Tvar'])} T-vars [{t_enc:.0f}s]", file=log, flush=True)
    cands = CANDS
    nsol = len(prev); nwit = 0; npg = 0; ndia = 0; nclause = 0; nround = 0; nrebuild = 0; t1 = time.time(); wit_hist = Counter()
    sync_shared()
    while state["s"].solve():
        s = state["s"]; model = s.get_model()
        val = [1 if lit > 0 else 0 for lit in model[:NV]]
        tab = [(val[2 * e] << 1) | val[2 * e + 1] for e in range(N)]
        learned = []
        for (sa, sb, ctx) in m1_diamonds(tab, KDIAM):
            wa, wb = diamond_windows(sa, ctx), diamond_windows(sb, ctx)
            learned.append([(int(wa[j]), int(wb[j])) for j in range(7)])
        ndia += bool(learned)
        if not learned:
            cols = periodic_collisions_all(tab, PMAX, 100)
            for (p, c, c2) in cols:
                W = PW[p]; wit_hist[p] += 1
                learned.append([(int(W[c, i]), int(W[c2, i])) for i in range(p) if W[c, i] != W[c2, i]])
            nwit += bool(cols)
        if learned:
            for prs in learned:
                s.add_clause([T(e, e2) for e, e2 in prs])
                shared.write(" ".join(f"{e},{e2}" for e, e2 in prs) + "\n")
            shared.flush(); nclause += len(learned); nround += 1
            if rss_mb() > MEM_LIMIT_MB:
                t2 = time.time(); shared_pos[0] = 0; nsynced[0] = 0; build_solver(); nsync = sync_shared(); nrebuild += 1
                print(f"    rebuilt solver ({nrebuild}): {rss_mb()} MB after rebuild, {nsync} clauses re-added [{time.time()-t2:.0f}s]", file=log, flush=True)
            if nround % 5 == 0:
                nsync = sync_shared()
                print(f"    round {nround}: {nclause} own clauses ({ndia} diamond, {nwit} periodic rounds {dict(wit_hist)}), {nsync} synced, {nsol} solutions, {npg} pair-graph, {rss_mb()} MB [{time.time()-t1:.0f}s]", file=log, flush=True)
            continue
        key = "".join(map(str, tab))
        if key in cands: label, inv = cands[key], True
        else:
            npg += 1; inv = injective_pairgraph_np(BlockMap(R, tab)); label = "OTHER"
        out.write(f"{int(inv)} {label} " + "".join(map(str, tab)) + "\n"); out.flush()
        nsol += 1
        blk = [-model[v - 1] for v in rep_vars]; state["blocks"].append(blk); s.add_clause(blk)
    state["s"].delete()
    print(f"cube {cube_id}: DONE {nsol} solutions passing periods<={PMAX} ({npg} needed the pair graph), {nround} learning rounds / {nclause} own clauses ({ndia} diamond, {nwit} periodic {dict(wit_hist)}), {nsynced[0]} synced, {nrebuild} rebuilds, encode {t_enc:.0f}s, enumerate {time.time()-t1:.0f}s", file=log, flush=True)
    out.close(); log.close(); shared.close()
    return cube_id, nsol, time.time() - t0

def choose_cube_entries(cands, k):
    A = np.array(list(cands.keys()), dtype=np.int8)
    H = []
    for e in range(N):
        c = np.bincount(A[:, e], minlength=4) / len(A)
        H.append(-(c[c > 0] * np.log2(c[c > 0])).sum())
    return [int(e) for e in np.argsort(-np.array(H))[:k]]

def run(args):
    os.makedirs(args.outdir, exist_ok=True)
    cands = candidates()
    print(f"{len(cands)} certified-invertible candidates built", flush=True)
    assert all(check_commutes_np_R(list(t), F_loc, R) for t in cands), "a candidate does not commute"
    json.dump({"labels": list(cands.values()), "tables": ["".join(map(str, t)) for t in cands]}, open(f"{args.outdir}/candidates.json", "w"))
    entries = choose_cube_entries(cands, args.cube_bits)
    cubes = [(i, list(zip(entries, vals)), args.maxperiod, args.outdir) for i, vals in enumerate(itertools.product(range(4), repeat=args.cube_bits))]
    if args.only is not None: cubes = [cubes[i] for i in args.only]
    print(f"{len(cubes)} cubes on entries {entries} ({[''.join(map(str, index_word(e, n))) for e in entries]}), {args.workers} workers, P={args.maxperiod}", flush=True)
    from multiprocessing import Pool
    summary = {}; t0 = time.time()
    with Pool(args.workers) as pool:
        for cube_id, nsol, dt in pool.imap_unordered(worker, cubes):
            summary[cube_id] = nsol
            print(f"  cube {cube_id:3d}: {nsol} solutions in {dt:.0f}s  (total so far {sum(summary.values())}, {len(summary)}/{len(cubes)} cubes, {time.time()-t0:.0f}s)", flush=True)
            json.dump(summary, open(f"{args.outdir}/summary.json", "w"))
    print(f"PHASE1 DONE: {sum(summary.values())} solutions in {time.time()-t0:.0f}s", flush=True)

# ---------------------------------------------------------------- phase 2
def pg_worker(tab):
    g = BlockMap(R, list(tab))
    wit = next((p for p in range(1, 9) if not injective_periodic(g, p)), None)
    if wit is not None: return tab, False, f"periodic witness p={wit}"
    return tab, injective_pairgraph_np(g), "pair-graph"

def classify(args):
    import glob
    cj = json.load(open(f"{args.outdir}/candidates.json"))
    cands = {tuple(int(ch) for ch in t): lab for t, lab in zip(cj["tables"], cj["labels"])}
    sols = set()
    for fn in sorted(glob.glob(f"{args.outdir}/cube_*.txt")):
        for l in open(fn):
            l = l.strip()
            if l: sols.add(tuple(int(ch) for ch in l))
    print(f"{len(sols)} distinct solutions from phase 1; {len(cands)} candidates", flush=True)
    t0 = time.time()
    bad = [t for t in sols if not check_commutes_np_R(list(t), F_loc, R)]
    print(f"commutation re-check: {len(bad)} failures [{time.time()-t0:.0f}s]", flush=True)
    missing = [lab for t, lab in cands.items() if t not in sols]
    print(f"candidates missing from the solution set: {len(missing)} {missing[:5]}", flush=True)
    matched = {t: cands[t] for t in sols if t in cands}
    left = [t for t in sols if t not in cands]
    print(f"matched (certified invertible): {len(matched)}; unmatched: {len(left)}", flush=True)
    from multiprocessing import Pool
    res = Counter(); new_inv = []
    with Pool(args.workers) as pool:
        for tab, inv, how in pool.imap_unordered(pg_worker, left, chunksize=4):
            res[(inv, how if not inv and how.startswith("periodic") else how)] += 1
            if inv: new_inv.append("".join(map(str, tab))); print("  NEW INVERTIBLE MAP outside the wall family", flush=True)
    print(f"unmatched classification: {dict(res)} [{time.time()-t0:.0f}s]", flush=True)
    total_inv = len(matched) + len(new_inv)
    print(f"RESULT 54R radius 3: invertible commuting maps = {total_inv} (named: {sum(1 for l in matched.values() if l.startswith('sigma') and 'Phi' not in l)}, "
          f"wall family: {sum(1 for l in matched.values() if 'Phi' in l)}, outside wall family: {len(new_inv)}); balanced+periodic-injective but not injective: {len(left) - len(new_inv)}", flush=True)
    json.dump({"n_solutions": len(sols), "n_matched": len(matched), "n_unmatched": len(left), "n_invertible": total_inv,
               "n_new_invertible": len(new_inv), "new_invertible_tables": new_inv, "missing_candidates": missing,
               "unmatched_classification": {str(k): v for k, v in res.items()}}, open(f"{args.outdir}/classify_54_R3.json", "w"), indent=1)
    open(f"{args.outdir}/classify.done", "w").close()

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["run", "classify"])
    ap.add_argument("--workers", type=int, default=24); ap.add_argument("--cube-bits", type=int, default=3)
    ap.add_argument("--maxperiod", type=int, default=5); ap.add_argument("--outdir", default="results")
    ap.add_argument("--only", type=int, nargs="*", default=None); ap.add_argument("--pmax", type=int, default=10); ap.add_argument("--pairs", default="33"); ap.add_argument("--solver", default="minicard"); ap.add_argument("--memlimit", type=int, default=3500); ap.add_argument("--nobalance1", action="store_true")
    a = ap.parse_args(); PMAX = a.pmax; SOLVER = a.solver; MEM_LIMIT_MB = a.memlimit; NOBAL1 = a.nobalance1; PAIRS = [(int(x[0]), int(x[1])) for x in a.pairs.split(",") if x]; run(a) if a.cmd == "run" else classify(a)
