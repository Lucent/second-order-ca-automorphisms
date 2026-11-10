# Confirmation of the hand proof of Proposition (F*_3) (class of 79, rule 93).  Not used in the proof.
#
# Checked, on the move graph Gamma_3 (1152 level-3 states, moves at the times |m|<=3):
#   L. Lemma "the moves at the times 0,...,3": the table of part (a) and the formulas of part (b), for every
#      state and every m in {0,1,2,3}, and the mirror rule T_{-m} = rho T_m rho for m in {1,2,3}.
#   1. Step 1: the two chains of arrows of the table of arrows; the 64 sectors of 12 configurations are strongly connected.
#   2. Step 2: its arrows in the table of arrows and the case distinction: every irregular configuration reaches a sector.
#   3. Step 3: its arrows in the table of arrows and the case distinction: every irregular configuration outside E is
#      reached from a sector, in the order of the text; the 20 states of E are not.
#   4. Step 4: the formula for (tau',mu'), the table of J, the closed walks.
#   5. Conclusion: the pairs (c,c') joined by the proof are exactly those of the Proposition; comparison with
#      breadth-first search; a single move into E starts in E.
#
# In the comments, (1)-(14) number the arrows of the proof in the order of the paper: (1), (2) the two groups of
# Step 1, (3)-(7) those of Step 2 and (8)-(12) those of Step 3 in the table of arrows, (13), (14) the two arrows of
# Step 4.  The functions P, Pbar are the relations (J), (J-bar) of the text.
#
# Usage:  python3 check_F3_hand.py            (Gamma_3 from graphs.py / engine.py of this folder)
#         python3 check_F3_hand.py core       (Gamma_3 from core.py: exact simulation with random left parts)
import sys
C4=[(0,0),(0,1),(1,0),(1,1)]
def load_graph(which):
    """G[state][m] = T_m(state) for the moves that act; a state is (theta,(a,tau,mu,eps)), cells (x,y)"""
    if which=='core':
        from core import graph
        G0=graph(3,seeds=(21,22,23))
        return {s:{m:t for m,(t,k) in G0[s].items()} for s in G0}
    from graphs import build
    L=3
    valid,edges,Z,allc,K=build(L,-L,L)
    def dec(i):
        return (i//K, tuple((int(c)>>1,int(c)&1) for c in allc[i%K]))
    G={dec(i):{} for i in range(6*K) if valid[i]}
    for u,kind,src,dst in edges:
        for s,d in zip(src.tolist(),dst.tolist()): G[dec(s)][u]=dec(d)
    return G
G=load_graph(sys.argv[1] if len(sys.argv)>1 else 'paper')
assert len(G)==1152 and sum(len(v) for v in G.values())==2688
def T(s,m): return G[s].get(m)
def f(a,b,c): return (1-c) | ((1-a)&b)
br=lambda u,v: int((u,v)!=(0,1))                         # <uv>
sw=lambda c:(c[1],c[0])
def st(th,a,tau,mu,eps): return (th,(a,tau,mu,eps))
def rho(s): return ((2-s[0])%6, tuple(sw(c) for c in s[1]))
def cs(c): return '%d%d'%c
def ss(s): return '%d|'%s[0]+' '.join(cs(c) for c in s[1])
VAC=[(0,0),(0,0),(0,0),(0,1),(1,1),(1,0)]
assert all(s[1][3]!=VAC[s[0]] for s in G) and len({(s[0],s[1][3]) for s in G})==18

# ---------------- L. the lemma on the moves at the times 0..3, literally ----------------
def lemma(s,m):
    """(theta', D) from the table of (a) and the formulas of (b); None if the move does not act"""
    th,(a,tau,mu,eps)=s
    (am,a0),(wm,w0),(ym,y0),(zm,z0)=a,tau,mu,eps
    if (m+th)%6 not in (0,2): return None
    y1=ym^f(w0,y0,z0); w1=wm^f(a0,w0,y0)
    pi=1^((1-z0)&y0)
    z1=None; z1n=None
    if th==0: z1=zm^(z0&(1-y0))
    if th in (3,5): z1=1-zm
    if th==5: z2=z0^(z1&(1-y1))
    row=None
    # table: (theta', D_r, z'_1)  ; D_r may refer to y'_0, filled in below
    if (m,th)==(0,2): row=(4 if zm==0 else 1, eps, None)
    if (m,th)==(0,0): row=(4 if z0==0 else 1, eps, None)
    if (m,th)==(1,1): row=(3,(0,0),None) if z0==0 else (0,eps,None)
    if (m,th)==(1,5): row=(3,eps,None) if zm==1 else (0,(1^(z0&(1-y0)),z0),None)
    if (m,th)==(2,0): row=(2,(1,0),None) if z1==0 else (5,(0,z0),None)
    if (m,th)==(2,4): row=(2,eps,None) if z0==1 else (5,(zm,1^((1-zm)&ym)),None)
    if (m,th)==(3,5):
        if z2==0: assert z1==1; row=(1,'F+',0)
        else: row=(4,(zm,0),z1)
    if (m,th)==(3,3):
        if zm==0: row=(1,(1,0),z1)
        else: assert z1==0; row=(4,(1-pi,z0),pi)
    thn,Dr,z1n=row
    # part (b)
    taun,mun=tau,mu
    if m==2 or (m==3 and z1n==z1):
        z0n=Dr[1]
        mun=(ym^(br(w0,y0)&(z0^z0n)),y0)
    elif m==3:
        dl=br(w1,y1); y0n=y0^dl
        if Dr=='F+': Dr=(1-y0n,1)
        z0n=Dr[1]
        taun=(wm^(dl&br(a0,w0)),w0)
        mun=(y1^f(w0,y0n,z0n),y0n)
    return st(thn,a,taun,mun,Dr)
nL=0
for s in G:
    for m in range(0,4):
        assert lemma(s,m)==T(s,m), ('lemma',ss(s),m)
        nL+=1
    for m in range(1,4):
        r=lemma(rho(s),m)
        assert (None if r is None else rho(r))==T(s,-m), ('mirror',ss(s),m)
        nL+=1
print('L ok: the table and part (b) of the lemma, and the mirror rule, on all states:',nL,'cases')

# ---------------- notation of the proof ----------------
def mx(tau,mu): return (mu[0]^br(tau[1],mu[1]), mu[1])   # mu^x
def my(tau,mu): return (mu[0], mu[1]^br(tau[0],mu[0]))   # mu^y
count=[0]
def arrow(s,m,t):
    assert T(s,m)==t, ('arrow fails',ss(s),m,ss(t))
    assert T(rho(s),-m)==rho(t), 'mirror image fails'
    count[0]+=1
E=set()
for a in C4:
    for th,eps in ((0,(1,1)),(1,(1,1)),(2,(1,1)),(3,(1,0)),(5,(0,1))):
        E.add(st(th,a,(0,0),(1,1),eps))
A_=[(4,(1,0)),(0,(1,0)),(5,(0,0)),(1,(0,1)),(0,(0,1))]
D_=[(4,(0,1)),(2,(0,1)),(3,(0,0)),(1,(1,0)),(2,(1,0))]
def Aset(a,tau,nu): return [st(th,a,tau,nu,e) for th,e in A_]
def Dset(a,tau,nu): return [st(th,a,tau,nu,e) for th,e in D_]
def Hset(a,tau,mu): return [st(3,a,tau,mu,(1,1)),st(5,a,tau,mu,(1,1))]
def Sigma(a,tau,mu): return Hset(a,tau,mu)+Aset(a,tau,mx(tau,mu))+Dset(a,tau,my(tau,mu))
sector={}
for a in C4:
    for tau in C4:
        for mu in C4:
            assert mx(tau,mx(tau,mu))==mu and my(tau,my(tau,mu))==mu
            assert sw(mx(tau,mu))==my(sw(tau),sw(mu))
            assert {rho(s) for s in Aset(a,tau,mu)}==set(Dset(sw(a),sw(tau),sw(mu)))
            assert {rho(s) for s in Sigma(a,tau,mu)}==set(Sigma(sw(a),sw(tau),sw(mu)))
            for s in Sigma(a,tau,mu):
                assert s in G and s not in sector
                sector[s]=(tau,mu)
assert len(sector)==4*16*12
irregular=[s for s in G if s not in sector]
assert len(irregular)==4*16*6
assert {(s[0],s[1][3]) for s in irregular}=={(0,(1,1)),(1,(1,1)),(2,(1,1)),(5,(0,1)),(3,(1,0)),(4,(0,0))}
assert E<=set(irregular)
def strongly_connected(nodes,succ):
    nodes=list(nodes)
    def reach(x,nb):
        seen={x}; stk=[x]
        while stk:
            u=stk.pop()
            for v in nb[u]:
                if v not in seen: seen.add(v); stk.append(v)
        return seen
    pred={u:[] for u in nodes}
    for u in nodes:
        for v in succ[u]: pred[v].append(u)
    return len(reach(nodes[0],succ))==len(nodes)==len(reach(nodes[0],pred))

# ---------------- Step 1 ----------------
for a in C4:
    for tau in C4:
        for nu in C4:
            c=lambda th,e: st(th,a,tau,nu,e)
            arrow(c(4,(1,0)),-2,c(0,(1,0))); arrow(c(0,(1,0)),2,c(5,(0,0))); arrow(c(5,(0,0)),-3,c(1,(0,1)))
            arrow(c(1,(0,1)),1,c(0,(0,1))); arrow(c(0,(0,1)),0,c(1,(0,1))); arrow(c(1,(0,1)),-1,c(5,(0,0)))
            arrow(c(5,(0,0)),1,c(0,(1,0))); arrow(c(0,(1,0)),0,c(4,(1,0)))
            arrow(c(4,(1,0)),2,st(5,a,tau,mx(tau,nu),(1,1)))
            arrow(st(5,a,tau,nu,(1,1)),3,st(4,a,tau,mx(tau,nu),(1,0)))
            arrow(st(5,a,tau,nu,(1,1)),1,st(3,a,tau,nu,(1,1)))
            arrow(st(3,a,tau,nu,(1,1)),-1,st(5,a,tau,nu,(1,1)))
for a in C4:
    for tau in C4:
        for mu in C4:
            S=set(Sigma(a,tau,mu))
            assert strongly_connected(S,{s:[t for t in G[s].values() if t in S] for s in S})
print('1 ok: 64 sectors of 12 configurations, each strongly connected')

# ---------------- Step 2 ----------------
def fwd(s):
    th,(a,tau,nu,eps)=s
    ym,y0=nu
    c=lambda th,e,nu2=None: st(th,a,tau,nu if nu2 is None else nu2,e)
    out=[]
    def go(m,t):
        cur=out[-1][1] if out else s
        assert T(cur,m)==t,('step 2',ss(cur),m,ss(t)); out.append((m,t)); count[0]+=1
    if (th,eps)==(4,(0,0)):
        if ym==1: go(2,c(5,(0,0)))                                     # (7)
        elif y0==1: go(-2,c(3,(0,0)))                                  # mirror image of (7)
        else: go(2,c(5,(0,1),(1,0))); go(1,c(0,(0,1),(1,0)))           # (7), (5)
        return out
    if (th,eps)==(5,(0,1)):
        if y0==0: go(1,c(0,(0,1)))                                     # (5)
        elif ym==1: go(-3,c(4,(0,1)))                                  # mirror image of (6)
        else: go(1,c(0,(1,1)))                                         # (5): into C(tau,01)
        if out[-1][1] in sector: return out
        s2=out[-1][1]; return out+fwd(s2)
    if (th,eps)==(3,(1,0)):
        if ym==0: go(-1,c(2,(1,0)))
        elif y0==1: go(3,c(4,(1,0)))                                   # (6)
        else: go(-1,c(2,(1,1)))
        if out[-1][1] in sector: return out
        s2=out[-1][1]; return out+fwd(s2)
    assert eps==(1,1) and th in (0,1,2)                                # C(tau,nu)
    if y0==0:
        if th==2: go(0,c(1,(1,1)))
        if th in (1,2): go(1,c(0,(1,1)))
        go(2,c(2,(1,0),mx(tau,nu))); return out                       # (3), (4)
    if ym==0:
        if th==0: go(0,c(1,(1,1)))
        if th in (0,1): go(-1,c(2,(1,1)))
        go(-2,c(0,(0,1),my(tau,nu))); return out                      # (3), mirror image of (4)
    if th==2: go(0,c(1,(1,1)))                                        # nu=11
    if th in (1,2): go(1,c(0,(1,1)))
    go(2,c(5,(0,1))); go(-3,c(4,(0,1)))                                # (4), mirror image of (6)
    return out
FW={}
for s in irregular:
    p=fwd(s); assert p[-1][1] in sector; FW[s]=p
print('2 ok: every irregular configuration reaches a sector (at most',max(len(p) for p in FW.values()),'moves)')

# ---------------- Step 3 ----------------
BW={}
def reach(pre,m,t):
    assert T(pre,m)==t,('step 3',ss(pre),m,ss(t)); count[0]+=1
    assert pre in sector or pre in BW, ('predecessor not yet reached',ss(pre))
    if t not in sector and t not in BW: BW[t]=(pre,m)
def mreach(pre,m,t): reach(rho(pre),-m,rho(t))
for R in (reach,mreach):
    for a in C4:
        for tau in C4:
            for nu in C4:
                c=lambda th,e: st(th,a,tau,nu,e)
                if nu[0]==1: R(c(5,(0,0)),3,c(4,(0,0)))                # (8)
                if nu[1]==0: R(c(0,(0,1)),2,c(5,(0,1)))                # (10)
    for a in C4:
        for tau in C4:
            R(st(5,a,tau,(1,0),(0,1)),3,st(4,a,tau,(0,0),(0,0)))       # (9)
assert all(st(4,a,tau,nu,(0,0)) in BW for a in C4 for tau in C4 for nu in C4)
def whole_class(a,tau,nu):
    c=lambda th,e: st(th,a,tau,nu,e)
    known=lambda x: x in BW or x in sector
    assert any(known(x) for x in (c(0,(1,1)),c(1,(1,1)),c(2,(1,1)))+((c(5,(0,1)),) if nu[1]==1 else ())+((c(3,(1,0)),) if nu[0]==1 else ())), ('class not reached',a,tau,nu)
    for _ in range(3):
        if nu[1]==1 and known(c(5,(0,1))): reach(c(5,(0,1)),1,c(0,(1,1)))        # (5)
        if nu[0]==1 and known(c(3,(1,0))): reach(c(3,(1,0)),-1,c(2,(1,1)))
        for (p,m,t) in ((c(0,(1,1)),0,c(1,(1,1))),(c(1,(1,1)),-1,c(2,(1,1))),(c(2,(1,1)),0,c(1,(1,1))),(c(1,(1,1)),1,c(0,(1,1)))):
            if known(p): reach(p,m,t)                                              # (3)
        if nu[1]==1 and known(c(0,(1,1))): reach(c(0,(1,1)),2,c(5,(0,1)))        # (4)
        if nu[0]==1 and known(c(2,(1,1))): reach(c(2,(1,1)),-2,c(3,(1,0)))
for a in C4:
    for tau in C4:
        wm,w0=tau
        reach(st(4,a,tau,(0,1),(0,0)),2,st(5,a,tau,mx(tau,(0,1)),(0,1)))          # (7): nu=11 (w0=1) or nu=01 (w0=0)
        assert mx(tau,(0,1))==((1,1) if w0==1 else (0,1))
        reach(st(4,a,tau,(1,0),(0,0)),-2,st(3,a,tau,my(tau,(1,0)),(1,0)))         # mirror image
        if wm==1: reach(st(5,a,tau,(0,0),(0,0)),3,st(1,a,tau,(1,0),(1,1)))        # (11)
        if w0==1: reach(st(3,a,tau,(0,0),(0,0)),-3,st(1,a,tau,(0,1),(1,1)))       # mirror image
        for nu in ((0,1),(1,0),(1,1)):
            if (tau,nu)!=((0,0),(1,1)): whole_class(a,tau,nu)
for a in C4:
    for u in (0,1):
        for w0 in (0,1):
            if (u,w0)==(0,0): continue                                             # source would lie in E
            reach(st(5,a,(u,w0),(1-w0,1),(0,1)),3,st(1,a,(u^br(a[1],w0),w0),(0,0),(1,1)))     # (12)
            reach(rho(st(5,sw(a),(u,w0),(1-w0,1),(0,1))),-3,rho(st(1,sw(a),(u^br(a[0],w0),w0),(0,0),(1,1))))
    for tau in C4: whole_class(a,tau,(0,0))
assert {s for s in irregular if s not in BW}==E
print('3 ok: every irregular configuration outside E is reached from a sector')

# ---------------- Step 4 ----------------
def P(a0,tau,mu):
    """the pair (tau',mu') of (14), in the two-case form of the text, and P_{a0}(tau,mu)=(tau',(mu')^y)"""
    wm,w0=tau; ym,y0=mu
    if w0==0: tp,mp=(tau,mu) if (wm==ym and ym!=y0) else ((1-wm,0),(1-ym,1-y0))
    else: tp,mp=(tau,mu) if (ym==1 and wm==1-(a0&y0)) else ((wm^a0,1),(ym,1-y0))
    return tp,mp,(tp,my(tp,mp))
PT={}
for a in C4:
    for tau in C4:
        for mu in C4:
            tp,mp,img=P(a[1],tau,mu)
            arrow(st(3,a,tau,mu,(1,1)),3,st(4,a,tp,mp,(0,1)))                      # (14)
            assert sector[st(4,a,tp,mp,(0,1))]==img
            PT[(a[1],tau,mu)]=img
            # general formula with delta
            wm,w0=tau; ym,y0=mu
            y1=ym^((1-w0)&y0); w1=wm^f(a[1],w0,y0); dl=br(w1,y1)
            assert (tp,mp)==((wm^(dl&br(a[1],w0)),w0),(ym^(dl&(1-w0)),y0^dl))
            if mu[1]==1:                                                            # (13) and (S)
                nu=mx(tau,mu)
                arrow(st(0,a,tau,nu,(0,1)),2,st(2,a,tau,mu,(1,0)))
                assert sector[st(0,a,tau,nu,(0,1))]==(tau,mu) and sector[st(2,a,tau,mu,(1,0))]==(tau,my(tau,mu))
TABLE={ # as printed in the text: rows tau (and a0), columns mu=00,01,10,11
 ('00',None):'10.10 00.00 10.00 10.01', ('10',None):'00.11 00.10 10.11 00.01',
 ('01',0):'01.00 01.01 01.11 01.10', ('11',0):'11.00 11.01 11.11 11.10',
 ('01',1):'11.00 11.01 11.10 01.11', ('11',1):'01.00 01.01 11.11 01.10'}
pc=lambda s:(int(s[0]),int(s[1]))
ps=lambda s:(pc(s[:2]),pc(s[3:]))
for (t,a0),row in TABLE.items():
    for mu,ent in zip(C4,row.split()):
        for a00 in ((0,1) if a0 is None else (a0,)):
            assert PT[(a00,pc(t),mu)]==ps(ent), ('table',t,a0,mu)
def Pmv(a,A): return PT[(a[1],)+A]
def Pbar(a,A):
    r=PT[(a[0],sw(A[0]),sw(A[1]))]; return (sw(r[0]),sw(r[1]))
def Smv(a,A): return (A[0],my(*A)) if A[1][1]==1 else None
def Sbar(a,A): return (A[0],mx(*A)) if A[1][0]==1 else None
MV={'P':Pmv,'Pb':Pbar,'S':Smv,'Sb':Sbar}
def walk(a,desc):
    toks=desc.split(); cur=ps(toks[0]); seen=[cur]
    for lab,nxt in zip(toks[1::2],toks[2::2]):
        assert MV[lab](a,cur)==ps(nxt),('walk',a,desc,lab,nxt)
        cur=ps(nxt); seen.append(cur)
    return seen
Wd='00.00 P 10.10 P 10.11 P 00.01 Pb 01.00 Pb 00.11 P 10.01 P 00.10 Pb 00.00'
Wm='00.00 Pb 01.01 Pb 01.11 Pb 00.10 P 10.00 P 00.11 Pb 01.10 Pb 00.01 P 00.00'
for a in C4:
    v=walk(a,Wd)+walk(a,Wm)
    assert set(v)=={(t,m) for t in C4 if t!=(1,1) for m in C4}
    if a[1]==1: walk(a,'01.00 P 11.00 P 01.00'); walk(a,'01.01 P 11.01 P 01.01'); walk(a,'01.10 P 11.10 P 11.11 P 01.10')
    if a[0]==1: walk(a,'10.00 Pb 11.00 Pb 10.00'); walk(a,'10.10 Pb 11.10 Pb 10.10'); walk(a,'10.01 Pb 11.01 Pb 11.11 Pb 10.01')
a=(0,0); tau=(1,1)
arrow(st(3,a,tau,(0,1),(0,0)),-3,st(4,a,tau,(0,1),(0,0)))
arrow(st(4,a,tau,(0,1),(0,0)),2,st(5,a,tau,(1,1),(0,1)))
arrow(st(5,a,tau,(1,1),(0,1)),-3,st(4,a,tau,(1,1),(0,1)))
assert sector[st(3,a,tau,(0,1),(0,0))]==ps('11.00') and sector[st(4,a,tau,(1,1),(0,1))]==ps('11.10')
walk(a,'11.10 P 11.11 Sb 11.01 S 11.00')
print('4 ok: the formula for (tau\',mu\'), the table of J and the closed walks')

# ---------------- Conclusion ----------------
def secgraph(a):   # only the relations established in Step 4
    g={(t,m):set() for t in C4 for m in C4}
    def add(desc):
        v=walk(a,desc)
        for x,y in zip(v,v[1:]): g[x].add(y)
    add(Wd); add(Wm)
    if a[1]==1: add('01.00 P 11.00 P 01.00'); add('01.01 P 11.01 P 01.01'); add('01.10 P 11.10 P 11.11 P 01.10')
    if a[0]==1: add('10.00 Pb 11.00 Pb 10.00'); add('10.10 Pb 11.10 Pb 10.10'); add('10.01 Pb 11.01 Pb 11.11 Pb 10.01')
    if a==(0,0): g[ps('11.00')].add(ps('11.10')); add('11.10 P 11.11 Sb 11.01 S 11.00')
    return g
def closure(g,A):
    seen={A}; stk=[A]
    while stk:
        x=stk.pop()
        for y in g[x]:
            if y not in seen: seen.add(y); stk.append(y)
    return seen
pairs=0
REACH={}
for a in C4:
    g=secgraph(a); R={A:closure(g,A) for A in g}
    for c in (s for s in G if s[1][0]==a):
        A=sector[c] if c in sector else sector[FW[c][-1][1]]
        for c2 in (s for s in G if s[1][0]==a):
            if c2 in E: continue
            if a==(0,0) and ((c2[1][1]==(1,1))!=(c[1][1]==(1,1))): continue
            x=c2
            while x not in sector: x=BW[x][0]
            assert sector[x] in R[A], ('not joined',ss(c),ss(c2))
            pairs+=1
def bfs(s):
    seen={s}; stk=[s]
    while stk:
        x=stk.pop()
        for t in G[x].values():
            if t not in seen: seen.add(t); stk.append(t)
    return seen
npairs=0
for c in G:
    Rc=bfs(c); a=c[1][0]
    for c2 in G:
        if c2[1][0]!=a: assert c2 not in Rc; continue
        allowed=(a!=(0,0)) or ((c2[1][1]==(1,1))==(c[1][1]==(1,1)))
        if allowed and c2 not in E: assert c2 in Rc; npairs+=1
        if not allowed: assert c2 not in Rc
        if c2 in E and c not in E: assert c2 not in Rc
assert npairs==pairs
for s in G:
    for m,t in G[s].items():
        if t in E: assert s in E
print('5 ok: the proof joins',pairs,'pairs (c,c\'), the same as breadth-first search; no move leads into E from outside')
print('hand proof of (F*_3) confirmed on Gamma_3: True')
