# Check of the marker-permutation construction for 90R and 165R in the paper (Theorem 7.10):
# segment words U_L (first/last letter 0, no 11, no 000, even- and odd-position sums 0), frames 111 u 111
# and complementary frames 000 ubar 000; psi_rho permutes the segment words; J is the integral;
# Phi_a(c) = c + iota(J(d c)), Phi_b(c) = c + iota'(J(e c)).
import itertools, random
random.seed(90)
def U(L):
    out = []
    for w in itertools.product((0, 1), repeat=L):
        s = ''.join(map(str, w))
        if w[0] or w[-1] or '11' in s or '000' in s: continue
        if sum(w[0::2]) % 2 or sum(w[1::2]) % 2: continue
        out.append(w)
    return out
print('|U_L| for L=1..18:', [len(U(L)) for L in range(1, 19)])
def fam(N):
    L = 4 * N + 11
    ws = [tuple(map(int, '0' + '10' * (2 * i + 1) + '100' + '10' + '100' + '10' * (2 * (N - i)))) for i in range(N + 1)]
    assert all(len(w) == L for w in ws)
    return L, ws
for N in range(0, 2):
    L, ws = fam(N); UL = set(U(L))
    print('N=%d L=%d: explicit family of %d words lies in U_L: %s, distinct: %s' % (N, L, len(ws), all(w in UL for w in ws), len(set(ws)) == len(ws)))

L, Uw = fam(2)            # three words of length 19 -> an S_3 action
Lens = {L: Uw}
def frames(z):
    """list of (p, L, sign, index) of frames of the cyclic binary word z."""
    n = len(z); out = []
    for p in range(n):
        for L, ws in Lens.items():
            seg = tuple(z[(p + 3 + k) % n] for k in range(L))
            for sg in (0, 1):
                mk = all(z[(p + k) % n] == 1 - sg for k in range(3)) and all(z[(p + L + 3 + k) % n] == 1 - sg for k in range(3))
                if not mk: continue
                u = tuple(b ^ sg for b in seg)
                if u in ws: out.append((p, L, sg, ws.index(u)))
    return out
def psi_J(z, rho, comp=True):
    n = len(z); y = list(z); J = [0] * n
    for p, L, sg, i in frames(z):
        if sg and not comp: continue
        u = Lens[L][i]; v = Lens[L][rho[i]]
        delta = [a ^ b for a, b in zip(u, v)]
        for k in range(L):
            y[(p + 3 + k) % n] = v[k] ^ sg
            J[(p + 3 + k) % n] = sum(delta[k::2]) % 2
    return y, J
def F(rule, x, y):
    n = len(x); eps = 1 if rule == 165 else 0
    return list(y), [x[i] ^ y[i - 1] ^ y[(i + 1) % n] ^ eps for i in range(n)]
d_of = lambda x, y: [x[(i + 1) % len(x)] ^ y[(i + 2) % len(x)] for i in range(len(x))]
e_of = lambda x, y: [y[i] ^ x[(i + 1) % len(x)] for i in range(len(x))]
def Phi(x, y, rho, mover, comp=True):
    n = len(x)
    if mover == 'a':
        _, J = psi_J(d_of(x, y), rho, comp); return [x[i] ^ J[i - 1] for i in range(n)], [y[i] ^ J[i] for i in range(n)]
    _, J = psi_J(e_of(x, y), rho, comp); return [x[i] ^ J[(i + 1) % n] for i in range(n)], [y[i] ^ J[i] for i in range(n)]
def plant(n):
    z = [random.randrange(2) for _ in range(n)]; p = 0
    while p + L + 8 < n:
        if random.random() < 0.6:
            sg = random.randrange(2); w = [1, 1, 1] + list(random.choice(Uw)) + [1, 1, 1]
            for k, b in enumerate(w): z[p + k] = b ^ sg
            p += L + 3 + random.choice((0, 0, 1, 2, 5))   # sometimes sharing a marker
        else: p += random.randint(1, 4)
    return z
def ring(n):                 # (x,y) on an odd ring with planted d and e
    while True:
        d = plant(n); e = plant(n)
        if sum(d) % 2 == sum(e) % 2: break
    # x_{j} = d_{j-1} + y_{j+1};  y_j = e_j + x_{j+1} = e_j + d_j + y_{j+2}
    y = [None] * n; y[0] = random.randrange(2); j = 0
    for _ in range(n - 1):
        jj = (j - 2) % n; y[jj] = e[jj] ^ d[jj] ^ y[j]; j = jj
    x = [d[(j - 1) % n] ^ y[(j + 1) % n] for j in range(n)]
    assert d_of(x, y) == d and e_of(x, y) == e
    return x, y
perms = list(itertools.permutations(range(3)))
comp_ = lambda s, r: tuple(s[r[i]] for i in range(3))
for rule in (90, 165):
    bad = dict(integral=0, frameset=0, dpsi=0, comm=0, hom=0, ab=0, inv=0); nontriv = 0; noncomm = 0
    for t in range(150):
        x, y = ring(181)
        z = d_of(x, y); Fr = sorted((p, Lq, sg) for p, Lq, sg, i in frames(z))
        for rho in perms:
            pz, J = psi_J(z, rho)
            n = len(z)
            bad['integral'] += any(z[i] ^ J[i] ^ J[(i + 2) % n] != pz[i] for i in range(n))
            bad['frameset'] += sorted((p, Lq, sg) for p, Lq, sg, i in frames(pz)) != Fr
            for mover in 'ab':
                x1, y1 = Phi(x, y, rho, mover)
                nontriv += (x1, y1) != (x, y)
                if mover == 'a': bad['dpsi'] += d_of(x1, y1) != pz or e_of(x1, y1) != e_of(x, y)
                bad['comm'] += F(rule, x1, y1) != Phi(*F(rule, x, y), rho, mover)
                inv = tuple(sorted(range(3), key=lambda i: rho[i]))
                bad['inv'] += Phi(x1, y1, inv, mover) != (x, y)
        for s, r in [(perms[1], perms[2]), (perms[3], perms[4])]:
            for mover in 'ab':
                bad['hom'] += Phi(*Phi(x, y, r, mover), s, mover) != Phi(x, y, comp_(s, r), mover)
                noncomm += Phi(*Phi(x, y, r, mover), s, mover) != Phi(*Phi(x, y, s, mover), r, mover)
        bad['ab'] += Phi(*Phi(x, y, perms[4], 'a'), perms[1], 'b') != Phi(*Phi(x, y, perms[1], 'b'), perms[4], 'a')
    print(rule, 'failures', bad, '| nontrivial', nontriv, '| noncommuting pairs', noncomm)
# negative control: 165 with only the 111-frames (J not complement-invariant)
bad = 0
for t in range(60):
    x, y = ring(181)
    for rho in perms[1:]:
        bad += F(165, *Phi(x, y, rho, 'a', comp=False)) != Phi(*F(165, x, y), rho, 'a', comp=False)
print('negative control 165, no complement frames: commutation failures', bad, 'of', 60 * 5)
# direct membership check of the explicit family for larger N (no enumeration of U_L)
for N in range(8):
    ws = ['0' + '10' * (2 * i + 1) + '100' + '10' + '100' + '10' * (2 * (N - i)) for i in range(N + 1)]
    ok = all(len(w) == 4 * N + 11 and w[0] == '0' == w[-1] and '11' not in w and '000' not in w
             and sum(map(int, w[0::2])) % 2 == 0 and sum(map(int, w[1::2])) % 2 == 0 for w in ws)
    print(N, ok, len(set(ws)))
