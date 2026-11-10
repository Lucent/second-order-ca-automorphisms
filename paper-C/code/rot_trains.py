"""Moving (time-periodic) walls for 219R and 126R.
219: bond b_i=[c_i==c_{i+1}] is frozen; a run of >=2 equal cells evolves by pi=(x,y)->(y,1-x) in every context.
126: bond b_i=[c_{i+1}==kappa(c_i)] is frozen; an alternating run s,ks,s,... evolves by pi in every context.
Train: cars (single cells, not in any bonded run) separated by bonded runs of length exactly 2, all the
bonded cells adjacent to cars in the train equal to the same symbol pattern; bounded by runs of length >=3.
Psi(sigma) permutes the cars of every m-car train."""
import numpy as np, itertools
rng = np.random.default_rng(7)
kap = lambda s: 3 - s
pi = lambda s: ((s & 1) << 1) | (1 - (s >> 1))
def F(rule, c):
    x = c >> 1; y = c & 1
    f = (rule >> (4*np.roll(y,1) + 2*y + np.roll(y,-1))) & 1
    return (y << 1) | (x ^ f)
def bond(rule, a, b):
    return a == b if rule == 219 else b == kap(a)
def partner(rule, s):          # the symbol following s inside a bonded run
    return s if rule == 219 else kap(s)

def trains(rule, c, m):
    """list of (positions of cars) for every train with exactly m cars."""
    N = len(c); B = [bond(rule, c[i], c[(i+1)%N]) for i in range(N)]
    # runs: maximal bonded intervals; a cell is in a run if B[i-1] or B[i]
    inrun = [B[(i-1)%N] or B[i] for i in range(N)]
    # run lengths
    out = []
    for i in range(N):
        # a train starts at a car i whose left run has length >=3
        if inrun[i]: continue
        # left run length
        L = 0; j = (i-1) % N
        while inrun[j] and B[(j) % N] if False else inrun[j]:
            L += 1; j = (j-1) % N
            if L > 6: break
        if L < 3: continue
        cars = [i]; p = i
        ok = False
        while True:
            # right run starting at p+1
            R = 0; j = (p+1) % N
            while inrun[j] and R < 7:
                R += 1; j = (j+1) % N
            if R == 2:
                q = (p+3) % N
                if inrun[q]: break
                cars.append(q); p = q
                if len(cars) > m: break
            elif R >= 3:
                ok = True; break
            else: break
        if ok and len(cars) == m:
            # require all wall cells adjacent to cars equal: left neighbours and right neighbours
            nb = set()
            for q in cars:
                nb.add(c[(q-1)%N]); nb.add(c[(q+1)%N])
            # for 219 neighbours left/right are both the run symbol; for 126 left neighbour s', right t
            lefts = {c[(q-1)%N] for q in cars}; rights = {c[(q+1)%N] for q in cars}
            if len(lefts) == 1 and len(rights) == 1:
                out.append(cars)
    return out

def Psi(rule, c, m, perm):
    d = c.copy()
    for cars in trains(rule, c, m):
        vals = [c[q] for q in cars]
        for j, q in enumerate(cars): d[q] = vals[perm[j]]   # car j receives old car perm[j]
    return d

def plant(rule, N):
    c = rng.integers(0, 4, N)
    p = 0
    while p < N - 40:
        if rng.random() < 0.5:
            s = rng.integers(0, 4); m = rng.integers(1, 5)
            seq = [s, partner(rule, s), s] if rule == 219 else [kap(s), s, kap(s)]  # left run length 3 ending ...
            # build: run(>=3) car run(2) car ... run(>=3); keep wall cells adjacent to cars equal
            if rule == 219:
                w3 = [s]*3; w2 = [s]*2
                cars = [rng.choice([v for v in range(4) if v != s]) for _ in range(m)]
                word = w3[:]
                for k, cv in enumerate(cars):
                    word.append(cv); word += (w2 if k < m-1 else w3)
            else:
                # left wall ends with t, right wall starts with u; inner walls of length 2 are (t',?) need
                # left neighbour of each car = a, right neighbour = b, inner wall = (b, kap(b)) must end in a => a = kap(b)
                b = rng.integers(0, 4); a = kap(b)
                cars = [rng.choice([v for v in range(4) if v != kap(a) and v != kap(b)]) for _ in range(m)]
                word = [kap(a), a, kap(a), a][-3:]   # run ending with a, length 3: a?,  make (a, ka, a)
                word = [a, kap(a), a]
                for k, cv in enumerate(cars):
                    word.append(cv)
                    word += ([b, kap(b)] if k < m-1 else [b, kap(b), b])
            for k, v in enumerate(word): c[p+k] = v
            p += len(word) + rng.integers(1, 4)
        else:
            p += rng.integers(1, 5)
    return c

for rule in (219, 126):
    bad = 0; nontriv = 0; noncomm = 0
    for trial in range(400):
        c = plant(rule, 400)
        # wall invariance check
        for m in (2, 3, 4):
            for perm in itertools.permutations(range(m)):
                d = Psi(rule, c, m, perm)
                if (d != c).any(): nontriv += 1
                if (F(rule, d) != Psi(rule, F(rule, c), m, perm)).any(): bad += 1
        m = 3
        a = Psi(rule, Psi(rule, c, 3, (1,0,2)), 3, (0,2,1)); b = Psi(rule, Psi(rule, c, 3, (0,2,1)), 3, (1,0,2))
        if (a != b).any(): noncomm += 1
    print(rule, 'commutation failures', bad, 'nontrivial', nontriv, 'noncommuting S3 pairs', noncomm)
