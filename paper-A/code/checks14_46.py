"""checks14_46.py -- confirmations for the classes of 14 and 46 (all proved by hand).

(1) Psi = kappa o rev conjugates F_14 to F_42 (Psi F_14 = F_42 Psi on random
    configurations), and Psi eta_p = rho eta_{-p}, Psi zeta_p = eta_{-p-1},
    Psi lambda_p = lambda_{1-p}  (equation cl:eq:Psi).
(2) The collision formula of Lemma cl:lem:collision for F_14, 1 <= N <= 30, for all
    times -5 <= t <= 3N+5.
(3) Rule 46: the bond bits chi_i(c) = [c_i = 00] xor [c_{i+1} = 11] satisfy
    chi_i(Fc) = chi_i(c) and chi_i(rho c) = chi_i(c), on all 4^4 values of c_{i-1..i+2}
    (Lemma cl:lem:bond(1)).
"""
import random, time
from probes import rulefn, Conf, step, Fp, rho

def Psi(C):
    return Conf(-C.hi(), [3 - c for c in reversed(C.a)], 3 - C.rt, 3 - C.lt)

def same(A, B, lo, hi): return all(A.at(i) == B.at(i) for i in range(lo, hi + 1))

t0 = time.time(); ok = True
f14, f42, f46 = rulefn(14), rulefn(42), rulefn(46)
rng = random.Random(1)
bad = 0
for _ in range(2000):
    C = Conf(-8, [rng.randrange(4) for _ in range(16)], rng.randrange(4), rng.randrange(4))
    if not same(Psi(step(f14, C)), step(f42, Psi(C)), -30, 30): bad += 1
print('(1) Psi F14 = F42 Psi: 2000 random configurations, %d mismatches' % bad); ok &= bad == 0
eta = lambda p: Conf(p, [2], 0, 3)
rhoeta = lambda p: Conf(p, [1], 0, 3)
zeta = lambda p: Conf(p + 1, [1], 0, 3)
lam = lambda p: Conf(p, [0], 3, 0)
bad = 0
for p in range(-10, 11):
    bad += not same(Psi(eta(p)), rhoeta(-p), -30, 30)
    bad += not same(Psi(zeta(p)), eta(-p - 1), -30, 30)
    bad += not same(Psi(lam(p)), lam(1 - p), -30, 30)
print('(1) images of eta_p, zeta_p, lambda_p, |p| <= 10: %d mismatches' % bad); ok &= bad == 0
bad = 0
for N in range(1, 31):
    E = Conf(0, [2] + [3] * (N - 1), 0, 0)
    for t in range(-5, 3 * N + 6):
        G = Fp(f14, E, t)
        for i in range(-N - 12, N + 12):
            if (G.at(i) & 1) != int(N - 1 - abs(t - N + 1) < i < N): bad += 1
print('(2) collision formula, 1 <= N <= 30: %d mismatches' % bad); ok &= bad == 0
chi = lambda C, i: int(C.at(i) == 0) ^ int(C.at(i + 1) == 3)
bad = 0
for w in range(256):
    C = Conf(-1, [(w >> (2 * j)) & 3 for j in range(4)], 0, 0)   # cells -1,0,1,2
    bad += chi(step(f46, C), 0) != chi(C, 0)
    bad += chi(rho(C), 0) != chi(C, 0)
print('(3) bond bits of 46 under F and rho, 256 cases: %d failures' % bad); ok &= bad == 0
print('time %.1f s' % (time.time() - t0))
if not ok: raise SystemExit('FAILURE')
print('ALL CHECKS PASSED')
