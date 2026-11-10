"""Rule 91: direct check on all periodic configurations of period <= 7 and random ones that
   Q1(c) = #{i: c_i != (1,1)} - #{i: c_i != c_{i+1}, both != (1,1)}      (number of maximal blocks without 11)
   Q2(c) = #{i: y_i != x_{i+1}} - #{i: x_i != y_{i+1}}                    (left-movers minus right-movers)
are invariant under F_91 (c_i = (x_i,y_i))."""
import itertools, random, sys
RULE = int(sys.argv[1]) if len(sys.argv) > 1 else 91
f = lambda a, b, c: (RULE >> (4*a + 2*b + c)) & 1
def F(x, y):
    n = len(x)
    return y, tuple(x[i] ^ f(y[i-1], y[i], y[(i+1) % n]) for i in range(n))
def Q1(x, y):
    n = len(x); c = [(x[i], y[i]) for i in range(n)]
    return sum(c[i] != (1, 1) for i in range(n)) - sum(c[i] != c[(i+1) % n] and c[i] != (1, 1) and c[(i+1) % n] != (1, 1) for i in range(n))
def Q2(x, y):
    n = len(x)
    return sum(y[i] != x[(i+1) % n] for i in range(n)) - sum(x[i] != y[(i+1) % n] for i in range(n))
bad1 = bad2 = cnt = 0
for n in range(1, 8):
    for bits in itertools.product((0, 1), repeat=2*n):
        x, y = bits[:n], bits[n:]; x2, y2 = F(x, y); cnt += 1
        bad1 += Q1(x, y) != Q1(x2, y2); bad2 += Q2(x, y) != Q2(x2, y2)
random.seed(1)
for _ in range(3000):
    n = random.randint(8, 60)
    x = tuple(random.randint(0, 1) for _ in range(n)); y = tuple(random.randint(0, 1) for _ in range(n))
    x2, y2 = F(x, y); cnt += 1
    bad1 += Q1(x, y) != Q1(x2, y2); bad2 += Q2(x, y) != Q2(x2, y2)
print("rule", RULE, "configs", cnt, "Q1 violations", bad1, "Q2 violations", bad2)
