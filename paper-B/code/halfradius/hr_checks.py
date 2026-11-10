"""Finite checks of Section 6 of the paper (rules that read two neighbours), for f = 187 and f = 119.

Usage:  python hr_checks.py [MW [MS]]      (defaults MW = 6, MS = 10)

Cells are z = 2x + y.  phi_nu(x,y) = (y, x XOR g(y,nu)) with f(a,b,c) = g(b,c), and
s_nu(x w) = phi_nu(x) s_{y(x)}(w) on words over the four cells (the first letter is the cell next to the
incoming bit).  Words in s_0^{+-1}, s_1^{+-1} are composed from right to left.

The program prints, for each of the two rules:
 1. the permutations phi_0, phi_1 (Table hr:tab:phi);
 2. the permutations pi of the cells with pi(phi_{y(c1)}(c0)) = phi_{y(pi c1)}(pi c0) for all cells c0, c1
    (Lemma hr:lem:P: only the identity);
 3. for the elements t_1, t_2 of Lemma hr:lem:sections: the image of x_*, the section at x_* as a word
    (product rule and free reduction), and a test of t(x_* w) = x_* t|_{x_*}(w) on all words w of length 7;
 4. the order of <phi_0^2, phi_1> and its 2-transitivity (Lemma hr:lem:level1);
 5. for M <= MW: strong connectivity of the window graph Lambda^0_M and of the two-copy window graph W_M,
    both built from Definition hr:def:graphs;
 6. for M <= MS: the order of the level-1 image of { t|_u : t in G, t(u) = u }, u = x_*^M, G = <s_0, s_1>,
    computed from Schreier generators of the stabilizer of u.
"""
import sys, itertools

XSTAR = {187: 1, 119: 3}
# the elements of Lemma hr:lem:sections, leftmost generator first; (nu, e) stands for s_nu^e
T = {187: {'t1': ([(1, 1), (0, -1), (0, -1), (1, 1), (1, 1)], [(1, 1)]),
           't2': ([(0, 1), (0, 1), (1, 1), (0, -1), (0, -1), (1, 1)], [(0, 1), (0, 1)])},
     119: {'t1': ([(1, 1)], [(1, 1)]),
           't2': ([(0, -1), (0, -1), (1, 1), (1, 1), (1, 1), (0, 1), (0, 1)], [(0, 1), (0, 1), (1, 1)])}}


def phis(rule):
    f = lambda a, b, c: (rule >> (4 * a + 2 * b + c)) & 1
    assert all(f(0, b, c) == f(1, b, c) for b in (0, 1) for c in (0, 1))
    return [tuple(((z & 1) << 1) | ((z >> 1) ^ f(0, z & 1, nu)) for z in range(4)) for nu in (0, 1)]


def inv(p):
    r = [0] * len(p)
    for i, v in enumerate(p):
        r[v] = i
    return tuple(r)


def comp(p, q):          # p after q
    return tuple(p[q[i]] for i in range(len(q)))


def cycles(p):
    seen = set(); out = []
    for i in range(len(p)):
        if i in seen:
            continue
        c = [i]; seen.add(i); j = p[i]
        while j != i:
            c.append(j); seen.add(j); j = p[j]
        if len(c) > 1:
            out.append("(" + " ".join(map(str, c)) + ")")
    return "".join(out) or "id"


def closure(gens):
    e = tuple(range(4)); G = {e}; fr = [e]
    while fr:
        nf = []
        for a in fr:
            for g in gens:
                c = comp(g, a)
                if c not in G:
                    G.add(c); nf.append(c)
        fr = nf
    return G


def s(ph, nu, w):
    out = []
    for x in w:
        out.append(ph[nu][x]); nu = x & 1
    return tuple(out)


def s_inv(ph, nu, w):
    out = []
    for x in w:
        xp = inv(ph[nu])[x]; out.append(xp); nu = xp & 1
    return tuple(out)


def act(ph, word, w):    # word: leftmost generator first, applied from right to left
    for nu, e in reversed(word):
        w = s(ph, nu, w) if e == 1 else s_inv(ph, nu, w)
    return w


def section(ph, word, x):
    """image of the letter x and the section at x, as a freely reduced word (leftmost generator first)."""
    sec = []
    for nu, e in reversed(word):
        if e == 1:
            sec.insert(0, (x & 1, 1)); x = ph[nu][x]
        else:
            x = inv(ph[nu])[x]; sec.insert(0, (x & 1, -1))
    red = []
    for gen in sec:
        if red and red[-1][0] == gen[0] and red[-1][1] == -gen[1]:
            red.pop()
        else:
            red.append(gen)
    return x, red


def show(word):
    return " ".join(f"s{nu}" if e == 1 else f"s{nu}^-1" for nu, e in word) or "1"


def strongly_connected(n, succ):
    """succ: list of successor lists on vertices 0..n-1."""
    pred = [[] for _ in range(n)]
    for v in range(n):
        for u in succ[v]:
            pred[u].append(v)
    for nb in (succ, pred):
        seen = [False] * n; seen[0] = True; st = [0]; cnt = 1
        while st:
            v = st.pop()
            for u in nb[v]:
                if not seen[u]:
                    seen[u] = True; cnt += 1; st.append(u)
        if cnt != n:
            return False
    return True


def window_graphs(ph, M):
    """Lambda^0_M and W_M of Definition hr:def:graphs; words omega_0 ... omega_M."""
    def Fnu(nu, om):
        return tuple(ph[om[k + 1] & 1][om[k]] for k in range(M)) + (ph[nu][om[M]],)
    words = list(itertools.product(range(4), repeat=M + 1)); wi = {w: i for i, w in enumerate(words)}
    lam = [[wi[Fnu(nu, w)] for nu in (0, 1)] for w in words]
    verts = [(w, z) for w in words for z in range(4) if z != w[0]]; vi = {v: i for i, v in enumerate(verts)}
    W = [[vi[(Fnu(nu, w), ph[(w[1] & 1) if M > 0 else nu][z])] for nu in (0, 1)] for (w, z) in verts]
    return strongly_connected(len(words), lam), strongly_connected(len(verts), W)


def stabilizer_image(ph, xs, M):
    """order of the level-1 image of the sections at u = xs^M of the stabilizer of u in G."""
    u = (xs,) * M
    P = {u: tuple(range(4))}; order = [u]; gens = set(); i = 0
    while i < len(order):
        v = order[i]; i += 1
        for nu in (0, 1):
            out = []; st = nu
            for x in v:
                out.append(ph[st][x]); st = x & 1
            v2 = tuple(out)
            q = comp(ph[st], P[v])     # level-1 action below v of the path u -> v followed by s_nu
            if v2 not in P:
                P[v2] = q; order.append(v2)
            else:
                gens.add(comp(inv(P[v2]), q))
    assert len(P) == 4 ** M               # G is transitive on level M
    return len(closure(gens))


if __name__ == '__main__':
    MW = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    MS = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    for rule in (187, 119):
        ph = phis(rule); xs = XSTAR[rule]
        print(f"rule {rule}, x_* = {xs}")
        for nu in (0, 1):
            print(f"  1. phi_{nu} = {ph[nu]} = {cycles(ph[nu])}")
        good = [p for p in itertools.permutations(range(4))
                if all(p[ph[c1 & 1][c0]] == ph[p[c1] & 1][p[c0]] for c0 in range(4) for c1 in range(4))]
        print(f"  2. cell permutations commuting with F: {[cycles(p) for p in good]}")
        for name, (word, target) in T[rule].items():
            x, sec = section(ph, word, xs)
            ok = all(act(ph, word, (xs,) + w) == (xs,) + act(ph, target, w) for w in itertools.product(range(4), repeat=7))
            print(f"  3. {name} = {show(word)}: x_* -> {x}; section = {show(sec)}; stated section {show(target)}: "
                  f"{sec == target}; on all words of length 7: {ok}")
        a = comp(ph[0], ph[0]); H1 = closure([a, ph[1]])
        two = len({(p[0], p[1]) for p in H1}) == 12
        print(f"  4. phi_0^2 = {cycles(a)}, phi_1 = {cycles(ph[1])}, |<phi_0^2, phi_1>| = {len(H1)}, 2-transitive: {two}")
        for M in range(MW + 1):
            l, w = window_graphs(ph, M)
            print(f"  5. M = {M}: Lambda^0_M strongly connected: {l}; W_M strongly connected: {w}", flush=True)
        orders = [stabilizer_image(ph, xs, M) for M in range(MS + 1)]
        print(f"  6. order of the level-1 image of the sections at x_*^M of the stabilizer, M = 0..{MS}: {orders}", flush=True)
