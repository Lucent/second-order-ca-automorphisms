"""Radius 1: the second-order automaton
F(x,y)_i = (y_i, x_i XOR f(y_{i-1},y_i,y_{i+1})), f an elementary rule in Wolfram numbering.
Cells are integers 2x+y.  Provides F on words and on periodic configurations, a count of
strongly connected components (Tarjan), and the light-speed window graph Lambda_M.
Run as a script: prints whether Lambda_M is strongly connected, M = 0..5, for the 16
elementary rules permutive in the left variable."""
import sys, itertools
sys.setrecursionlimit(1000000)

def fbit(rule, a, b, c):
    return (rule >> (4*a + 2*b + c)) & 1

def Fword(rule, w):
    """apply F to a word of cells; returns word shorter by 2 (radius 1)."""
    out = []
    for i in range(1, len(w)-1):
        x, y = w[i] >> 1, w[i] & 1
        out.append(2*y + (x ^ fbit(rule, w[i-1] & 1, y, w[i+1] & 1)))
    return tuple(out)

def Fper(rule, c):
    """F on a spatially periodic configuration given by one period (tuple of cells)."""
    n = len(c)
    return tuple(2*(c[i] & 1) + ((c[i] >> 1) ^ fbit(rule, c[(i-1) % n] & 1, c[i] & 1, c[(i+1) % n] & 1)) for i in range(n))

def scc_count(nv, succ):
    """number of SCCs (iterative Tarjan). succ(v)->iterable"""
    index = [-1]*nv; low = [0]*nv; onst = [False]*nv; st = []; cnt = 0; idx = 0
    for s in range(nv):
        if index[s] != -1: continue
        work = [(s, iter(succ(s)))]
        index[s] = low[s] = idx; idx += 1; st.append(s); onst[s] = True
        while work:
            v, it = work[-1]
            adv = False
            for w in it:
                if index[w] == -1:
                    index[w] = low[w] = idx; idx += 1; st.append(w); onst[w] = True
                    work.append((w, iter(succ(w)))); adv = True; break
                elif onst[w]:
                    low[v] = min(low[v], index[w])
            if adv: continue
            work.pop()
            if work:
                u = work[-1][0]; low[u] = min(low[u], low[v])
            if low[v] == index[v]:
                cnt += 1
                while True:
                    w = st.pop(); onst[w] = False
                    if w == v: break
    return cnt

def lambda_graph(rule, M):
    """light-speed window graph Lambda_M (r=1): vertices words of length M+1,
    edges omega -> F(omega r1 r2) restricted to sites p+1..p+M+1."""
    words = list(itertools.product(range(4), repeat=M+1))
    idx = {w: i for i, w in enumerate(words)}
    succ = []
    for w in words:
        s = set()
        for r1 in range(4):
            for r2 in range(4):
                s.add(idx[Fword(rule, w + (r1, r2))])
        succ.append(sorted(s))
    return succ

def is_sc(succ):
    return scc_count(len(succ), lambda v: succ[v]) == 1

if __name__ == "__main__":
    print("== Lambda_M strongly connected? (r=1) ==")
    nonaff = [30, 120, 45, 180, 75, 210, 135, 225]
    aff = [15, 60, 90, 105, 150, 165, 195, 240]
    for rule in nonaff + aff:
        res = []
        for M in range(0, 6):
            res.append(is_sc(lambda_graph(rule, M)))
        print(rule, res)
