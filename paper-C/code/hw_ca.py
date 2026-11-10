"""Local rule of the second-order automaton F_f and of its inverse as 4x4x4 NumPy tables (cells 2x+y),
a step on a cyclic configuration, and the image of a word. Used by hw59_walls.py, hw59_trains.py,
hw59_period.py, hwtp_walls.py and hwtp_trains.py.
"""
import numpy as np, itertools
def rule_f(rule):
    return lambda a,b,c: (rule>>((a<<2)|(b<<1)|c))&1
def tables(rule):
    f=rule_f(rule)
    Ft=np.zeros((4,4,4),dtype=np.int64); Fi=np.zeros((4,4,4),dtype=np.int64)
    for l,s,r in itertools.product(range(4),repeat=3):
        x,y=s>>1,s&1
        Ft[l,s,r]=(y<<1)|(x^f(l&1,y,r&1))          # x'=y, y'=x^f(y_l,y,y_r)
        # inverse: given c' (x',y'): y=x', x=y'^f(x'_l,x',x'_r)
        xp,yp=x,y
        Fi[l,s,r]=((yp^f(l>>1,xp,r>>1))<<1)|xp
    return Ft,Fi
def step(T,c): return T[np.roll(c,1),c,np.roll(c,-1)]
def step_word(T,w):  # word -> word shorter by 2
    return tuple(int(T[w[i-1],w[i],w[i+1]]) for i in range(1,len(w)-1))
