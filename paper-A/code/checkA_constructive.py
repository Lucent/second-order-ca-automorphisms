"""checkA_constructive.py -- check 1 of 2 of the proofs of (F_k) (110) and (A),(B) (42).

The hand proof is run as an algorithm.  For a 0-tailed configuration c with edge r and a
level k it produces an explicit word of moves (m, bar) with |m| <= k-1 that takes c to the
canonical configuration K_k(c): c with the cells r-k+1..r replaced by (1,0), except that for
110 a frozen tip (exterior 00, tip 11) is kept.  The word is built exactly as in the proof:
  * tip phase: Steps 2, 3, 4 of the proof, each realised as F^{-1} O' F with O' a word at
    level k-1 joining Fb to the configuration b~ named in the proof (b~ = Fb with the cell
    r-k+2 replaced by the cell given in the proof), and its rho-conjugates;
  * Step 1: the word for level k-1 applied to b itself.
The words at level k-1 joining Fb and b~ are canon(Fb) followed by the reverse of canon(b~);
the script ASSERTS that K_{k-1}(Fb) = K_{k-1}(b~), i.e. that b~ lies in the class the proof
claims (this is where the hypotheses at level k-1 are used), and that the tip after each step
is the one claimed in the proof.  Finally every word is replayed from c with the moves of the
tables (Lemma cl:lem:edgemove, Lemma cl:lem:moves1446(a)) and the end point compared with K_k(c).

Configurations: finitely supported, as pairs of Python integers (X, Y) = (row -1, row 0),
bit d = column r-d; unbounded integers make the dynamics exact (no truncation).
"""
import sys, random, itertools, time
sys.setrecursionlimit(100000)

RULE = None
def fw(Y):
    a, b, c = Y >> 1, Y, Y << 1          # left neighbour, cell, right neighbour
    if RULE == 110: return (b ^ c) | (b & ~a)
    if RULE == 42:  return c & ~(a & b)
    raise ValueError

def Fs(s):  X, Y = s; return (Y, X ^ fw(Y))
def Fi(s):  X, Y = s; return (Y ^ fw(X), X)
def ev(s, t):
    while t > 0: s = Fs(s); t -= 1
    while t < 0: s = Fi(s); t += 1
    return s
def rho(s): return (s[1], s[0])
def cell(s, d): return ((s[0] >> d) & 1, (s[1] >> d) & 1)
def setcell(s, d, xy):
    X, Y = s; m = 1 << d
    X = (X & ~m) | (xy[0] << d); Y = (Y & ~m) | (xy[1] << d)
    return (X, Y)

def table(xr, yr, xl, yl, bar):
    if not bar:
        if (xr, yr) == (0, 1): return (1, 0)
        if xr == 1:
            if RULE == 110: return (1, yr ^ xl)
            return (1, yr ^ 1 ^ xl)
        return (xr, yr)
    else:
        if (xr, yr) == (1, 0): return (0, 1)
        if yr == 1:
            if RULE == 110: return (xr ^ yl, 1)
            return (xr ^ 1 ^ yl, 1)
        return (xr, yr)

def move(s, m, bar):
    t = ev(s, m)
    xr, yr = cell(t, 0); xl, yl = cell(t, 1)
    t = setcell(t, 0, table(xr, yr, xl, yl, bar))
    return ev(t, -m)

def apply(s, word):
    for m, bar in word: s = move(s, m, bar)
    return s

def reverse_word(s, word):
    """word takes s to s'; return a word taking s' back to s (each move undone at its time)."""
    states = [s]
    for m, bar in word: states.append(move(states[-1], m, bar))
    out = []
    for i in range(len(word) - 1, -1, -1):
        m = word[i][0]; ok = None
        for bb in (0, 1):
            if move(states[i + 1], m, bb) == states[i]: ok = (m, bb); break
        assert ok is not None, 'move not undone at the same time'
        out.append(ok)
    return out

def frozen110(s, k):  # exterior 00 and tip 11
    return cell(s, k) == (0, 0) and cell(s, k - 1) == (1, 1)

def K(s, k):
    if RULE == 42:
        assert not (cell(s, k) == (1, 1) and cell(s, k - 1) == (1, 1)), 'canonical form undefined'
    t = s
    for d in range(k):
        t = setcell(t, d, (1, 0))
    if RULE == 110 and frozen110(s, k):
        t = setcell(t, k - 1, (1, 1))
    return t

STATS = {}
memo = {}
def canon(s, k):
    key = (s, k)
    if key in memo: return memo[key]
    assert cell(s, 0) != (0, 0)
    if k == 1:
        # (F_1) / k=1 of (A),(B): moves at time 0 on the edge cell, at most two steps
        tgt = K(s, 1); seen = {s: []}; frontier = [s]
        while frontier:
            nf = []
            for u in frontier:
                for bar in (0, 1):
                    v = move(u, 0, bar)
                    if v not in seen: seen[v] = seen[u] + [(0, bar)]; nf.append(v)
            frontier = nf
        assert tgt in seen, 'k=1 failed'
        memo[key] = seen[tgt]; return seen[tgt]
    word = []; cur = s
    for _ in range(4):                       # tip phase: at most 3 steps
        tip = cell(cur, k - 1); ext = cell(cur, k)
        tau = (1, 1) if (RULE == 110 and frozen110(s, k)) else (1, 0)
        if tip == tau: break
        u0, v0 = ext; u1, v1 = tip
        if RULE == 110:
            if tip == (0, 0):   op = ('step2', False)          # 00 -> 10
            elif tip == (0, 1): op = ('step2', True)           # 01 -> 00 (rho of 10 -> 00)
            elif tip == (1, 1):
                assert ext != (0, 0)
                op = ('step4', False) if v0 == 1 else ('step4', True)   # 11->01 / 11->10
            else: raise AssertionError
        else:
            if tip == (0, 0):   op = ('step2', False)
            elif tip == (0, 1): op = ('step2', True)
            elif tip == (1, 1):
                assert ext != (1, 1)
                op = ('step4', False) if v0 == 0 else ('step4', True)
            else: raise AssertionError
        w = lift(cur, k, op[0], op[1])
        word += w; cur = apply(cur, w)
    assert cell(cur, k - 1) == ((1, 1) if (RULE == 110 and frozen110(s, k)) else (1, 0))
    w = canon(cur, k - 1)                     # Step 1
    word += w; cur = apply(cur, w)
    assert cur == K(s, k), 'canonical form not reached'
    memo[key] = word
    return word

def lift(b, k, step, use_rho):
    """one application of Step 2 or Step 4 (or its rho-conjugate, Step 3)."""
    if use_rho:
        w = lift(rho(b), k, step, False)
        return [(-m, 1 - bar) for m, bar in w]
    pre = []
    if RULE == 110 and step == 'step2' and cell(b, k - 1) == (1, 0) and cell(b, k - 2) == (0, 1):
        # pre-step (Step 1): replace the cell r-k+2 by (1,0); the tip 10 is not 00
        tg = setcell(b, k - 2, (1, 0))
        assert K(b, k - 1) == K(tg, k - 1)
        pre = canon(b, k - 1) + reverse_word(tg, canon(tg, k - 1))
        assert apply(b, pre) == tg
        b0 = b; b = tg
        STATS[(k, 'pre')] = STATS.get((k, 'pre'), 0) + 1
    u0, v0 = cell(b, k); u1, v1 = cell(b, k - 1); u2, v2 = cell(b, k - 2)
    Fb = Fs(b)
    if RULE == 110:
        if step == 'step2': assert v1 == 0; new = (1 - v2, v2)
        else:               assert v1 == 1 and v0 == 1; new = (1 - v2, 1)
    else:
        if step == 'step2': assert v1 == 0; new = (1 - v2, 1)
        else:               assert v1 == 1 and v0 == 0; new = (1 - v2, v2)
    bt = setcell(Fb, k - 2, new)
    assert cell(bt, 0) != (0, 0)
    # b~ must lie in O_{k-1}(Fb): same canonical form at level k-1
    assert K(Fb, k - 1) == K(bt, k - 1), ('b~ not in the claimed class', step, k)
    w1 = canon(Fb, k - 1)
    w2 = reverse_word(bt, canon(bt, k - 1))
    wl = w1 + w2
    assert apply(Fb, wl) == bt
    wb = [(m + 1, bar) for m, bar in wl]      # F^{-1} M_m F = M_{m+1}
    assert all(abs(m) <= k - 1 for m, _ in wb)
    b2 = apply(b, wb)
    assert b2 == Fi(bt)
    # claimed tip: (u1 ^ 1, v1) for Step 2 / Step 4; left part unchanged
    assert cell(b2, k - 1) == (1 - u1, v1), ('tip claim', step)
    for d in range(k, k + 3 * k + 5): assert cell(b2, d) == cell(b, d)
    STATS[(k, step)] = STATS.get((k, step), 0) + 1
    return pre + wb

def random_conf(k, extra, rng):
    """cells r-k..r random (edge != 00), then `extra` random cells further left"""
    while True:
        X = rng.getrandbits(k + 1 + extra); Y = rng.getrandbits(k + 1 + extra)
        if (X | Y) & 1: return (X, Y)

def run(rule, kmax_exh, kmax_rand, nrand, seed=1):
    global RULE, memo
    RULE = rule; memo = {}
    rng = random.Random(seed)
    # sanity: rho M_m rho = Mbar_{-m} on random configurations
    for _ in range(300):
        s = random_conf(5, 8, rng); m = rng.randint(-4, 4); bar = rng.randint(0, 1)
        assert rho(move(rho(s), m, bar)) == move(s, -m, 1 - bar)
    total = 0
    for k in range(1, kmax_rand + 1):
        confs = []
        if k <= kmax_exh:
            # all values of the cells r-k..r (edge != 00) and all values of the 2 cells
            # further left, plus 3 random fillings of 3k+3 cells further left
            for bits in range(1 << (2 * (k + 3))):
                X = bits & ((1 << (k + 3)) - 1); Y = bits >> (k + 3)
                if (X | Y) & 1: confs.append((X, Y))
            for _ in range(3):
                lx = rng.getrandbits(3 * k + 3) << (k + 1); ly = rng.getrandbits(3 * k + 3) << (k + 1)
                for bits in range(1 << (2 * (k + 1))):
                    X = bits & ((1 << (k + 1)) - 1); Y = bits >> (k + 1)
                    if (X | Y) & 1: confs.append((X | lx, Y | ly))
        else:
            confs = [random_conf(k, 3 * k + 3, rng) for _ in range(nrand)]
        nconf = 0; maxlen = 0
        for s in confs:
            if rule == 42 and cell(s, k) == (1, 1) and cell(s, k - 1) == (1, 1): continue
            w = canon(s, k)
            assert all(abs(m) <= k - 1 for m, _ in w)
            assert apply(s, w) == K(s, k)          # independent replay
            nconf += 1; maxlen = max(maxlen, len(w)); total += 1
        print('rule %d k=%d: %d configurations taken to the canonical form, longest word %d, %s'
              % (rule, k, nconf, maxlen, 'exhaustive on r-k-2..r (+3 random fillings)' if k <= kmax_exh else 'random'))
        sys.stdout.flush()
    print('rule %d: lift steps used:' % rule, dict(sorted(STATS.items())))
    STATS.clear()
    return total

if __name__ == '__main__':
    t0 = time.time()
    ke = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    kr = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    nr = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    for rule in (110, 42):
        run(rule, ke, kr, nr)
    print('time %.1f s' % (time.time() - t0))
    print('ALL CHECKS PASSED')
