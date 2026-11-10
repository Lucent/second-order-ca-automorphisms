# h = h_* = sigma^{-1} hatH_{phi_*,id} for 95R is not in <sigma,F>.   usage: python3 standalone7.py tx7.npy ty7.npy
import numpy as np, sys
TX, TY = np.load(sys.argv[1]), np.load(sys.argv[2])
f = lambda a, b, c: (95 >> (4*a + 2*b + c)) & 1                                          # Wolfram rule 95
F = lambda c: tuple((c[i][1], c[i][0] ^ f(c[i-1][1], c[i][1], c[(i+1) % len(c)][1])) for i in range(len(c)))
S = lambda c, a=1: tuple(c[(i+a) % len(c)] for i in range(len(c)))                       # shift
w = lambda c, i, s: sum(c[(i+m) % len(c)][(s+m) % 2] << (m+7) for m in range(-7, 8))     # s=0: x at even m, y at odd m
g = lambda c: tuple((int(TX[w(c, i, 1)]), int(TY[w(c, i, 0)])) for i in range(len(c)))   # cells (x_i, y_i)
c = ((0,0), (0,0), (0,0), (0,1), (1,0)); orbit = [c]                                     # spatial period 5
while F(orbit[-1]) != c: orbit.append(F(orbit[-1]))
print("c    =", c, "\nh(c) =", g(c), "\nF-period of c:", len(orbit))
print("h(F(d)) == F(h(d))      on the F-orbit:", all(g(F(d)) == F(g(d)) for d in orbit))
print("(sigma h)^2 (d) == d    on the F-orbit:", all(S(g(S(g(d)))) == d for d in orbit))
print("sigma h (d) != d        on the F-orbit:", all(S(g(d)) != d for d in orbit))
print("h(c) = sigma^a F^b (c) for some a, b  :", any(S(d, a) == g(c) for d in orbit for a in range(5)))
