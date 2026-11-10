# The identity behind the conservation of Q_1 for rule 91, checked exhaustively.
# Cells are written (u, m): u the first bit, m the second; F replaces u by
# u XOR s, s_i = f(m_{i-1}, m_i, m_{i+1}), and exchanges the two bits.
# (1) the density of Q_1 and the regrouped density of the proof differ by a
#     telescoping term; (2) the description of s in the proof, and the
#     invariance of every term of the regrouped density; (3) Q_1 on all
#     configurations of spatial period at most 8.
# usage: python3 q1_identity.py
import itertools
def f(rule,a,b,c): return (rule>>(4*a+2*b+c))&1
def h(c0,c1):   # original density on two adjacent cells, cells are (x,y)
    return int(c0!=(1,1)) - int(c0!=c1 and c0!=(1,1) and c1!=(1,1))
def k(um0,um1): # k_i = [m_i=0, m_{i+1}=1](1-u_{i+1}); cells written (u,m)
    return int(um0[1]==0 and um1[1]==1)*(1-um1[0])
def hp(l,c,r):  # regrouped density h'_i, cells (u,m) at i-1,i,i+1
    (ul,ml),(u,m),(ur,mr)=l,c,r
    if m==1: return (1-int(ml==0)-int(mr==0))*(1-u)
    return 1-int(mr==0)*(u^ur)
def run(rule):
    bad_regroup=bad_inv=bad_s=0
    # (1) h_i - h'_i = k_i - k_{i-1} for all triples of cells (u,m)
    cells=[(u,m) for u in (0,1) for m in (0,1)]
    for l,c,r in itertools.product(cells,repeat=3):
        # h uses cells as unordered pairs (x,y); (u,m) or (m,u) give the same value
        if hp(l,c,r)-h(c,r) != k(c,r)-k(l,c): bad_regroup+=1
        if h(c[::-1],r[::-1])!=h(c,r): bad_regroup+=1
    # (2) h'_i(u+s(m), m) = h'_i(u, m): m on cells i-2..i+2, u on i-1..i+1
    for m in itertools.product((0,1),repeat=5):
        s=[f(rule,m[j-1],m[j],m[j+1]) for j in (1,2,3)]
        # description of s used in the proof
        for j in (1,2,3):
            if m[j]==1: exp=int((m[j-1]==0)!=(m[j+1]==0))
            else: exp=int(m[j-1]==0 or m[j+1]==0)
            if exp!=s[j-1]: bad_s+=1
        for u in itertools.product((0,1),repeat=3):
            v=[u[j]^s[j] for j in range(3)]
            A=hp((u[0],m[1]),(u[1],m[2]),(u[2],m[3])); B=hp((v[0],m[1]),(v[1],m[2]),(v[2],m[3]))
            if A!=B: bad_inv+=1
    # (3) brute force on periodic configurations
    badQ=0; n_cfg=0
    for n in range(1,9):
        for x in itertools.product((0,1),repeat=n):
            for y in itertools.product((0,1),repeat=n):
                c=[(x[i],y[i]) for i in range(n)]
                d=[(y[i], x[i]^f(rule,y[i-1],y[i],y[(i+1)%n])) for i in range(n)]
                Q=lambda c: sum(h(c[i],c[(i+1)%n]) for i in range(n))
                n_cfg+=1
                if Q(c)!=Q(d): badQ+=1
    print("rule",rule,"regroup violations",bad_regroup,"description of s violations",bad_s,"termwise invariance violations",bad_inv,"| periodic configs",n_cfg,"Q1 violations",badQ)
for r in (91,94,95): run(r)
