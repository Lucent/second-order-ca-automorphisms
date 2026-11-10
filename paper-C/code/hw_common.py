"""Local rule of the second-order automaton F_f and of its inverse (cells 2x+y), as functions and as tables;
a step on a cyclic configuration; the image of the interior of a word. Used by hw59_table.py,
hw59_overlap.py, hw59_classes.py, hw59_census_lib.py, hw59_generators.py, hw27_*.py and hwtp_cycles.py.
"""
import itertools, numpy as np
def frule(rule):
    return lambda a,b,c: (rule>>((a<<2)|(b<<1)|c))&1
def Floc(rule):
    f=frule(rule)
    def g(s0,s1,s2):
        x,y=s1>>1,s1&1
        return (y<<1)|(x^f(s0&1,y,s2&1))
    return g
def Finvloc(rule):
    f=frule(rule)
    def g(s0,s1,s2):
        x,y=s1>>1,s1&1
        return ((y^f(s0>>1,x,s2>>1))<<1)|x
    return g
def tables(rule):
    F=Floc(rule); Fi=Finvloc(rule)
    T=np.array([[[F(a,b,c) for c in range(4)] for b in range(4)] for a in range(4)],dtype=np.int8)
    Ti=np.array([[[Fi(a,b,c) for c in range(4)] for b in range(4)] for a in range(4)],dtype=np.int8)
    return T,Ti
def step(T,c):  # cyclic
    return T[np.roll(c,1),c,np.roll(c,-1)]
def wimg(loc,w): # image of interior of word, length len-2
    return tuple(loc(*w[i:i+3]) for i in range(len(w)-2))
