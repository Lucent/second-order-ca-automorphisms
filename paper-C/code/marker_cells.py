"""Markers for 90R and 165R at the level of cells (Theorem 7.9 and Example 7.11).
(1) exact identities on all small windows (these make commutation an exact consequence of J o sigma = sigma o J
    and J(not z) = J(z));
(2) the S_3 instance Phi^a_rho(c) = c + iota(J_rho(d(c))), Phi^b_rho(c) = c + iota'(J_rho(e(c))) on ALL rings of
    period 10 and 11 (cyclic configurations, exact for those) and on planted random rings.
Conventions: cell (x,y); (sigma v)_i = v_{i+1};  d_n = x_{n+1}^y_{n+2};  e_m = y_m^x_{m+1};
 iota(v):  x_i = v_{i-1}, y_i = v_i;     iota'(v): x_i = v_{i+1}, y_i = v_i."""
import numpy as np, itertools
U = ['0010110', '0101010', '0110100']
perms = {'id': (0, 1, 2), '(01)': (1, 0, 2), '(12)': (0, 2, 1), '(02)': (2, 1, 0), '(012)': (1, 2, 0), '(021)': (2, 0, 1)}
def compose(s, r): return tuple(s[r[i]] for i in range(3))          # (s o r)(i) = s(r(i))
def inverse(r): inv = [0]*3; [inv.__setitem__(r[i], i) for i in range(3)]; return tuple(inv)
Ub = [int(u, 2) for u in U]
def build_tables(comp=True):
    fid = np.zeros(8192, dtype=np.int8)
    for i, u in enumerate(Ub):
        pat = (7 << 10) | (u << 3) | 7; fid[pat] = 1 + i
        if comp: fid[(~pat) & 0x1FFF] = -(1 + i)
    return fid
def Jtab(rho):
    t = np.zeros((3, 7), dtype=np.int8)
    for i in range(3):
        delta = [((Ub[i] ^ Ub[rho[i]]) >> (6 - k)) & 1 for k in range(7)]
        for k in range(7): t[i, k] = sum(delta[k::2]) % 2
    return t
def J_of(z, rho, fid):
    """z: (M,n) binary, cyclic.  returns J (M,n)."""
    M, n = z.shape
    idx = np.zeros((M, n), dtype=np.int64)
    for k in range(13): idx = (idx << 1) | np.roll(z, -k, axis=1)       # idx[:,p] = z[p..p+12]
    f = fid[idx]; i = np.abs(f).astype(np.int64) - 1; has = f != 0
    t = Jtab(rho); J = np.zeros_like(z)
    for k in range(7):
        contrib = np.where(has, t[np.clip(i, 0, 2), k], 0).astype(z.dtype)
        J ^= np.roll(contrib, 3 + k, axis=1)                              # position p+3+k
    return J
R_ = lambda a, k: np.roll(a, -k, axis=-1)                                 # R_(a,k)_i = a_{i+k}
def F(rule, x, y):
    f = (rule >> (4 * R_(y, -1) + 2 * y + R_(y, 1))) & 1
    return y.copy(), x ^ f
d_of = lambda x, y: R_(x, 1) ^ R_(y, 2)
e_of = lambda x, y: y ^ R_(x, 1)
def Phi(x, y, rho, mover, fid):
    if mover == 'a': J = J_of(d_of(x, y), rho, fid); return x ^ R_(J, -1), y ^ J
    else:            J = J_of(e_of(x, y), rho, fid); return x ^ R_(J, 1), y ^ J

# ------------------------------------------------------------ (1) exact window identities
def exact_identities():
    out = {}
    for rule, eps in ((90, 0), (165, 1)):
        f = lambda a, b, c: (rule >> (4*a + 2*b + c)) & 1
        lin = lambda a, b, c: a ^ c
        # F affine with linear part F_90: on all pairs of 3-windows
        ok = all(((x1 ^ x2) == (x1 ^ x2)) and ((f(*y1) ^ f(*y2)) == lin(*(p ^ q for p, q in zip(y1, y2))))
                 for y1 in itertools.product((0, 1), repeat=3) for y2 in itertools.product((0, 1), repeat=3) for x1 in (0, 1) for x2 in (0, 1))
        out[(rule, 'F affine, linear part F90')] = ok
        # d(Fc)_n = d(c)_{n+1} + eps  : reads c[n..n+3]
        ok = True; ok2 = True
        for w in itertools.product(range(4), repeat=4):
            X = [s >> 1 for s in w]; Y = [s & 1 for s in w]
            Fx = lambda i: Y[i]; Fy = lambda i: X[i] ^ f(Y[i-1], Y[i], Y[i+1])
            # positions 0..3 = n..n+3 ; d_n(Fc) = Fx(n+1) ^ Fy(n+2)
            if (Fx(1) ^ Fy(2)) != (X[2] ^ Y[3] ^ eps): ok = False
            # e_m(Fc) = Fy(m) ^ Fx(m+1) ; with m = n+2 (needs Y[m-1..m+1]) -> e_{m-1}(c) + eps = Y[m-1]^X[m]
            m = 2
            if (Fy(m) ^ Fx(m + 1)) != (Y[m - 1] ^ X[m] ^ eps): ok2 = False
        out[(rule, 'd(Fc) = sigma d(c) + eps')] = ok; out[(rule, 'e(Fc) = sigma^-1 e(c) + eps')] = ok2
    # linear-part identities with binary v (window of 6 bits, positions -2..3)
    ok = {k: True for k in ['F90(iota v) = iota(sigma v)', "F90(iota' v) = iota'(sigma^-1 v)", 'd(iota v) = v + sigma^2 v',
                            'e(iota v) = 0', "d(iota' v) = 0", "e(iota' v) = v + sigma^2 v"]}
    for vv in itertools.product((0, 1), repeat=9):
        v = lambda i: vv[i + 4]                      # i in -4..4
        ix = lambda i: v(i - 1); iy = lambda i: v(i)          # iota
        jx = lambda i: v(i + 1); jy = lambda i: v(i)          # iota'
        i = 0
        # F90 linear: (x,y) -> (y_i, x_i ^ y_{i-1} ^ y_{i+1})
        if not (iy(i) == v(i - 1 + 1) and (ix(i) ^ iy(i - 1) ^ iy(i + 1)) == v(i + 1)): ok['F90(iota v) = iota(sigma v)'] = False
        if not (jy(i) == v(i + 1 - 1) and (jx(i) ^ jy(i - 1) ^ jy(i + 1)) == v(i - 1)): ok["F90(iota' v) = iota'(sigma^-1 v)"] = False
        n = -1
        if (ix(n + 1) ^ iy(n + 2)) != (v(n) ^ v(n + 2)): ok['d(iota v) = v + sigma^2 v'] = False
        if (iy(n) ^ ix(n + 1)) != 0: ok['e(iota v) = 0'] = False
        if (jx(n + 1) ^ jy(n + 2)) != 0: ok["d(iota' v) = 0"] = False
        if (jy(n) ^ jx(n + 1)) != (v(n) ^ v(n + 2)): ok["e(iota' v) = v + sigma^2 v"] = False
    out.update(ok)
    return out

for k, v in exact_identities().items(): print('exact identity', k, ':', v)

# ------------------------------------------------------------ (2) rings
def ring_checks(rule, x, y, fid, label):
    bad = dict(comm=0, inv=0, hom=0, ab=0); nontriv = 0; noncomm = 0
    for mover in ('a', 'b'):
        for name, rho in perms.items():
            x1, y1 = Phi(x, y, rho, mover, fid)
            nontriv += int(((x1 != x).any(1) | (y1 != y).any(1)).sum()) if name != 'id' else 0
            a1 = Phi(*F(rule, x, y), rho, mover, fid); a2 = F(rule, x1, y1)
            bad['comm'] += int(((a1[0] != a2[0]).any(1) | (a1[1] != a2[1]).any(1)).sum())
            z = Phi(x1, y1, inverse(rho), mover, fid)
            bad['inv'] += int(((z[0] != x).any(1) | (z[1] != y).any(1)).sum())
        for s, r in (('(01)', '(12)'), ('(012)', '(01)'), ('(02)', '(021)')):
            z1 = Phi(*Phi(x, y, perms[r], mover, fid), perms[s], mover, fid)
            z2 = Phi(x, y, compose(perms[s], perms[r]), mover, fid)
            bad['hom'] += int(((z1[0] != z2[0]).any(1) | (z1[1] != z2[1]).any(1)).sum())
        z1 = Phi(*Phi(x, y, perms['(01)'], mover, fid), perms['(12)'], mover, fid)
        z2 = Phi(*Phi(x, y, perms['(12)'], mover, fid), perms['(01)'], mover, fid)
        noncomm += int(((z1[0] != z2[0]).any(1) | (z1[1] != z2[1]).any(1)).sum())
    # a- and b-movers commute
    z1 = Phi(*Phi(x, y, perms['(012)'], 'a', fid), perms['(01)'], 'b', fid)
    z2 = Phi(*Phi(x, y, perms['(01)'], 'b', fid), perms['(012)'], 'a', fid)
    bad['ab'] += int(((z1[0] != z2[0]).any(1) | (z1[1] != z2[1]).any(1)).sum())
    print(f'  {rule}R {label}: rings {x.shape[0]} x {x.shape[1]}, failures {bad}, nontrivial applications {nontriv}, '
          f'rings where (01),(12) do not commute {noncomm}')

def all_rings(n):
    e = np.arange(4 ** n, dtype=np.int64)
    c = np.stack([(e >> (2 * k)) & 3 for k in range(n)], axis=1).astype(np.int64)
    return c >> 1, c & 1
rng = np.random.default_rng(16)
def planted(M, n, comp):
    x = rng.integers(0, 2, (M, n)); y = rng.integers(0, 2, (M, n))
    for r in range(M):
        p = 0
        while p + 14 < n - 2:
            if rng.random() < 0.6:
                sgn = rng.integers(0, 2) if comp else 0
                w = [1, 1, 1] + [int(ch) for ch in U[rng.integers(0, 3)]] + [1, 1, 1]
                for k, b in enumerate(w): x[r, (p + k + 1) % n] = (b ^ sgn) ^ y[r, (p + k + 2) % n]   # d_n = x_{n+1}^y_{n+2}
                if rng.random() < 0.5:   # also plant a b-mover frame on e_m = y_m ^ x_{m+1} at a shifted place
                    q = (p + 3) % n
                p += 10 + rng.integers(0, 4)
            else: p += rng.integers(1, 5)
    # plant b-frames in separate rings: e_m = y_m ^ x_{m+1}  -> set y_m
    return x, y
def planted_b(M, n, comp):
    x = rng.integers(0, 2, (M, n)); y = rng.integers(0, 2, (M, n))
    for r in range(M):
        p = 0
        while p + 14 < n - 2:
            if rng.random() < 0.6:
                sgn = rng.integers(0, 2) if comp else 0
                w = [1, 1, 1] + [int(ch) for ch in U[rng.integers(0, 3)]] + [1, 1, 1]
                for k, b in enumerate(w): y[r, (p + k) % n] = (b ^ sgn) ^ x[r, (p + k + 1) % n]
                p += 10 + rng.integers(0, 4)
            else: p += rng.integers(1, 5)
    return x, y

fid_sym = build_tables(True); fid_plain = build_tables(False)
for rule in (90, 165):
    for n in (10, 11):
        x, y = all_rings(n)
        ring_checks(rule, x, y, fid_sym, f'ALL rings of period {n}, complement-symmetric frames')
    x, y = planted(3000, 240, True); ring_checks(rule, x, y, fid_sym, 'planted a-frames (both signs)')
    x, y = planted_b(3000, 240, True); ring_checks(rule, x, y, fid_sym, 'planted b-frames (both signs)')
print('negative control: 165R with 111-frames only (J not complement-symmetric)')
x, y = planted(3000, 240, True); ring_checks(165, x, y, fid_plain, 'planted, plain frames')
x, y = all_rings(10); ring_checks(165, x, y, fid_plain, 'ALL rings period 10, plain frames')
