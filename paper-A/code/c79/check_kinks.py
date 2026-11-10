# Confirmation (not used in the proof): vacuum, constants, kinks, and the shift of finite type Fix(sigma^{-s} F^3) for rule 93.
# All objects are periodic outside a finite window, so every identity below is checked exactly
# on a window wide enough to contain the non-periodic part plus the radius of the maps involved.
import itertools
from collections import defaultdict
from base import f, Vac, Fstep_word, Finv_word
def vcell(w,i): return 2*Vac(w-1,i)+Vac(w,i)          # cell i of the vacuum V_w (state at time w)
RHO=[0,2,1,3]
# --- vacuum -----------------------------------------------------------------------------------
for t in range(6):
    for i in range(2):
        assert Vac(t+1,i)^Vac(t-1,i)==f(Vac(t,i-1),Vac(t,i),Vac(t,i+1))      # spacetime relation
        assert Vac(t+3,i)==Vac(t,i+1)                                         # F^3 V = sigma V
        assert Vac(1-t,i)==Vac(t,i)                                           # rho V_w = V_{2-w}
        b,c=Vac(t,i),Vac(t,i+1); assert (b,c)!=(1,1) and f(0,b,c)==f(1,b,c)  # left-blindness
words=[(vcell(w,0),vcell(w,1)) for w in range(6)]
assert len(set(words))==6
print('vacuum: relation, F^3V=sigma V, rho V_w=V_{2-w}, left-blind; (even,odd) cells of V_0..V_5:',words)
print('constants: F maps 0,1,2,3 to',[Fstep_word([c]*3)[0] for c in range(4)])
# --- F-orbits on Fix(sigma^2) -------------------------------------------------------------------
def Fper(p):
    r=Fstep_word([p[i%2] for i in range(6)]); return (r[1],r[2])
seen=set(); orbs=[]
for p in itertools.product(range(4),repeat=2):
    if p in seen: continue
    o=[p]
    while Fper(o[-1])!=p: o.append(Fper(o[-1]))
    seen|=set(o); seen|=set((b,a) for a,b in o); orbs.append(o)
print('F-orbits on Fix(sigma^2) (pairs = cells at even, odd sites):',orbs)
six=[set(o) for o in orbs if len(o)==6]
assert len(six)==2 and {(0,1),(1,0)} <= six[0]|six[1] and any({(1,2),(2,1)}<=o for o in six), six
assert sorted(len(o) for o in orbs)==[1,3,6,6]
print('period-2 configurations: up to sigma, the 6-orbits of F are exactly the vacuum and W')
Wset=[o for o in orbs if len(o)==6 and (0,1) not in o and (1,0) not in o][0]
# --- kinks ----------------------------------------------------------------------------------------
LO,HI=-60,60
def Km0(i): return vcell(0,i) if i<=1 else 0
def Kp0(i): return vcell(0,i) if i<=0 else 1
def ev(fn,n,lo=LO,hi=HI):
    m=abs(n); w=[fn(i) for i in range(lo-m,hi+m+1)]
    for _ in range(m): w=Fstep_word(w) if n>0 else Finv_word(w)
    return {lo+j:w[j] for j in range(len(w))}
K3=ev(Km0,3); assert all(K3[i]==Km0(i+1) for i in range(LO,HI))
K3=ev(Kp0,3); assert all(K3[i]==Kp0(i-1) for i in range(LO,HI))
K2=ev(Kp0,2); assert all(K2[i]==RHO[Km0(i)] for i in range(LO,HI))
print('F^3 K^- = sigma K^-, F^3 K^+ = sigma^{-1} K^+, rho K^- = F^2 K^+ : ok')
def first_nonvac(cfg,w): return next(i for i in range(LO,HI) if cfg[i]!=vcell(w,i))
pm=[first_nonvac(ev(Km0,m),m%6) for m in range(-12,13)]
qp=[first_nonvac(ev(Kp0,n),n%6) for n in range(-12,13)]
assert all(pm[j]==2-((m-1)//3) for j,m in enumerate(range(-12,13)))
assert all(qp[j]==2+((n+1)//3) for j,n in enumerate(range(-12,13)))
assert all(ev(Km0,m)[50]==m%3 for m in range(6)) and all(ev(Kp0,n)[50]==(n+1)%3 for n in range(6))
print('first non-vacuum cell: p^-(m) = 2-floor((m-1)/3), q^+(n) = 2+floor((n+1)/3); right constants m mod 3, n+1 mod 3')
def const_from(cfg):  # least s such that cfg is constant on [s, HI]
    c=cfg[HI-1]; return max(i for i in range(LO,HI) if cfg[i]!=c)+1
assert all(const_from(ev(Km0,m))==2-((m+1)//3) for m in range(-12,13))
print('K^-_m is constant exactly on [2-floor((m+1)/3), oo), -12 <= m <= 12')
# --- Fix(sigma^s F^3) as a one-step SFT on 7-windows ---------------------------------------------
def analyse(s):
    E=set(w for w in itertools.product(range(4),repeat=7)
          if Fstep_word(Fstep_word(Fstep_word(list(w))))[0]==w[3+s])
    while True:   # prune to windows lying on a bi-infinite path
        heads=set(e[1:] for e in E); tails=set(e[:-1] for e in E)
        E2=set(e for e in E if e[:-1] in heads and e[1:] in tails)
        if E2==E: break
        E=E2
    succ=defaultdict(set)
    for e in E: succ[e[:-1]].add(e[1:])
    V=set(e[:-1] for e in E)|set(e[1:] for e in E)
    # strongly connected components (Kosaraju)
    order=[]; vis=set()
    for v0 in V:
        if v0 in vis: continue
        st=[(v0,iter(succ[v0]))]; vis.add(v0)
        while st:
            v,it=st[-1]
            for w in it:
                if w not in vis: vis.add(w); st.append((w,iter(succ[w]))); break
            else: st.pop(); order.append(v)
    pred=defaultdict(set)
    for v in succ:
        for w in succ[v]: pred[w].add(v)
    cid={}; comps=[]
    for v0 in reversed(order):
        if v0 in cid: continue
        c=[v0]; cid[v0]=len(comps); st=[v0]
        while st:
            v=st.pop()
            for w in pred[v]:
                if w not in cid: cid[w]=len(comps); c.append(w); st.append(w)
        comps.append(c)
    rec=[i for i,c in enumerate(comps) if any(cid[w]==i for v in c for w in succ[v])]
    def kind(i):
        c=comps[i]
        if len(c)==1 and len(set(c[0]))==1: return 'const%d'%c[0][0]
        if len(c)==2 and all(len(set(v))==2 and v[0]==v[2] for v in c):
            pr=(c[0][0],c[0][1]); return 'vac' if (pr in words or pr[::-1] in words) else 'per2'
        return 'big%d'%len(c)
    kinds={i:kind(i) for i in rec}
    # transient paths between recurrent components (the transient part is acyclic)
    paths=defaultdict(list)
    for A in rec:
        for v in comps[A]:
            for w in succ[v]:
                if cid[w]==A: continue
                st=[(w,(w[-1],))]
                while st:
                    x,p=st.pop()
                    if cid[x] in rec: paths[(A,cid[x])].append((v,p)); continue
                    for y in succ[x]: st.append((y,p+(y[-1],)))
    return V,succ,cid,comps,rec,kinds,paths
for s,name in ((1,'Fix(sigma^{-1}F^3)'),(-1,'Fix(sigma F^3)')):
    V,succ,cid,comps,rec,kinds,paths=analyse(s)
    print(name,': vertices',len(V),'recurrent components',sorted(kinds.values()))
    big=[i for i in rec if kinds[i].startswith('big')]; assert len(big)==1; big=big[0]
    assert sorted(kinds.values())==sorted(['const0','const1','const2','const3','vac','vac','vac',kinds[big]])
    # W lies in the big component
    for (a,b) in Wset: assert cid[(a,b)*3]==big
    # the constant loops have no exits: every edge from a vertex of a constant loop stays in it
    for A in rec:
        if kinds[A].startswith('const'):
            assert all(cid[w]==A for v in comps[A] for w in succ[v]), kinds[A]
            assert not any(A2==A for (A2,B) in paths)
    # no transient path leaves the big component; no path from a vacuum to another vacuum or to W
    assert not any(A==big for (A,B) in paths)
    for (A,B),lst in paths.items():
        if kinds[A]=='vac':
            assert kinds[B]!='vac'
            if kinds[B].startswith('const'): assert len(lst)==1, (A,B,lst)
    for A in rec:
        if kinds[A]=='vac':
            outs=[B for (A2,B) in paths if A2==A and kinds[B].startswith('const')]
            assert len(outs)==1
            v,p=paths[(A,outs[0])][0]
            print('   vacuum',''.join(map(str,comps[A][0][:2])),'-> unique kink to',kinds[outs[0]],':',''.join(map(str,v))+'|'+''.join(map(str,p)))
print('kinks ok: W lies in the big component, which has no exit; the constant loops have no exits; each vacuum pair has exactly one kink to a constant')
