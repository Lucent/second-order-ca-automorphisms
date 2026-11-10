# Verification of the scattering identities used in the proof of the scattering-closure theorem (rule 94R).
# Points are kink lists: (lam, [kinks]) with value(p) = lam xor #{b<=p}; kinks >= 2 (corner [0,0]).
# usage: python3 verify1.py [lib|direct|goursat_cc] [ncases]
#   lib: Goursat solver of lib.py (default); direct: (t,i) simulation of cc.py; goursat_cc: Goursat solver of cc.py
import random, sys
from lib import scatter
def word(lam,ks,n):
    out=[];v=lam;s=set(ks)
    for p in range(n):
        if p in s: v^=1
        out.append(str(v))
    return ''.join(out)
def rep(s): return (int(s[0]),[k for k in range(1,len(s)) if s[k]!=s[k-1]])
def S_lib(P,X):
    n=max(P[1]+X[1]+[2])+3
    g,o,F=scatter(word(P[0],P[1],n),word(X[0],X[1],n))
    return rep(g),rep(o),F
S=S_lib
def use(impl):
    global S
    if impl=='lib': S=S_lib
    else:
        import cc
        if impl=='direct':
            def S2(P,X):
                G,O=cc.direct((P[0],list(P[1])),(X[0],list(X[1])))
                return (int(G[0]),list(G[1])),(int(O[0]),list(O[1])),None
        else:
            def S2(P,X):
                (G,O),F=cc.goursat((P[0],list(P[1])),(X[0],list(X[1])),fires=True)
                return (int(G[0]),list(G[1])),(int(O[0]),list(O[1])),F
        S=S2
def valid(X):
    n=max(X[1]+[2])+4
    return '010' not in word(X[0],X[1],n) and all(k>=2 for k in X[1])
def randobj(lam,first,nk,maxgap=5):
    while True:
        ks=[first]
        for i in range(nk-1): ks.append(ks[-1]+random.randint(1,maxgap))
        if valid((lam,ks)): return (lam,ks)
def main(N):
    random.seed(1)
    bad=0; n=0
    for t in range(N):   # first-kink claims
        lam=random.randint(0,1); U=randobj(lam,random.randint(4,7),random.randint(1,7))
        a=random.randint(2,5); G,O,F=S((lam,[a]),U); u1=U[1][0]
        if lam==0: ok=(G[0]==1 and G[1][0]==u1)
        elif (a-u1)%2==0: ok=(G[0]==0 and G[1][0]>=u1+1)
        else: ok=(G[0]==0 and G[1][0]==u1)
        n+=1; bad+=not ok
    print('first-kink claims',n,'bad',bad)
    bad=0;n=0
    for t in range(N):   # A-pass
        U=randobj(0,random.randint(6,9),random.randint(1,7)); u1=U[1][0]
        g=random.randint(2,u1-3); b1=u1-g; a=random.randint(2,5)
        if (a-b1)%2==0: a+=1
        G,O,F=S((1,[a]),(1,[b1]+U[1])); G2,O2,F2=S((0,[a]),U)
        ok=(G==(0,[b1]+G2[1]) and G2[0]==1 and O==O2 and G2[1][0]>=u1)
        n+=1; bad+=not ok
    print('A-pass',n,'bad',bad)
    bad=0;n=0
    for t in range(N):   # B-pass
        U=randobj(1,random.randint(6,9),random.randint(1,7)); u1=U[1][0]
        g=random.randint(2,u1-3); b1=u1-g; a=random.randint(2,5)
        G,O,F=S((0,[a]),(0,[b1]+U[1])); G2,O2,F2=S((1,[a]),U)
        ok=(G==(1,[b1]+G2[1]) and G2[0]==0 and O==O2)
        if (a-u1)%2==0: ok=ok and G2[1][0]>=u1+1
        n+=1; bad+=not ok
    print('B-pass',n,'bad',bad)
    bad=0;n=0
    for a in range(2,12):
        for d in range(2,12):
            G,O,F=S((1,[a]),(1,[d,d+1])); n+=1     # exchange
            if (a-d)%2==0:
                ok=(G==(0,[d+2]) and O==(1,[a+1,a+2]) and (F is None or sorted(F)==[(a-1,d-1),(a+2,d+2)]))
                bad+=not ok
            G,O,F=S((1,[d,d+1]),(1,[a])); n+=1     # mirrored exchange
            if (a-d)%2==0:
                ok=(G==(1,[a+1,a+2]) and O==(0,[d+2]) and (F is None or sorted(F)==[(d-1,a-1),(d+2,a+2)]))
                bad+=not ok
            G,O,F=S((1,[a,a+1]),(1,[d,d+1])); n+=1 # two zeros of different parity
            if (a-d)%2==1:
                ok=(G==(1,[d,d+1]) and O==(1,[a,a+1]) and (F is None or F==[]))
                bad+=not ok
    print('exchange / mirrored / zero-zero',n,'bad',bad)
    bad=0;n=0
    for A in range(2,9):                           # seed lemma
        for b in (4,5):
            r=b+A
            for a in (2,3,4,5):
                if (a-r)%2: continue
                G,O,F=S((0,[a]),(0,[b,r,r+2])); n+=1
                ok=(G==(1,[b,r+1,r+3,r+4]) and O==(1,[a+1,a+2]) and (F is None or sorted(F)==[(a-1,r-1),(a+2,r+4)]))
                bad+=not ok
    print('seed',n,'bad',bad)
    bad=0;n=0
    for t in range(N):   # Z-prepend, y>=3
        U=randobj(1,random.randint(8,11),random.randint(1,7)); u1=U[1][0]
        y=random.randint(3,u1-4); d=u1-1-y; a=random.choice([2,3,4,5])
        if (a-d)%2: a+=1
        G,O,F=S((1,[a]),(1,[d,d+1]+U[1])); G2,O2,F2=S((1,[a+1,a+2]),U)
        ok=(G2[0]==1 and G==(0,[d+2]+G2[1]) and O==O2 and G2[1][0]>=u1)
        n+=1; bad+=not ok
    print('zero prepend y>=3',n,'bad',bad)
    bad=0;n=0
    for t in range(N):   # Z-probe pass
        U=randobj(1,random.randint(7,9),random.randint(1,7)); u1=U[1][0]
        y=random.choice([2,2,2,3,4]) ; d=u1-1-y
        if d<2: continue
        a=random.choice([2,3,4,5])
        if (a-d)%2==0: a+=1
        G,O,F=S((1,[a,a+1]),(1,[d,d+1]+U[1])); G2,O2,F2=S((1,[a,a+1]),U)
        ok=(G2[0]==1 and G==(1,[d,d+1]+G2[1]) and O==O2)
        if (a-u1)%2==0: ok=ok and G2[1][0]>=u1+1
        n+=1; bad+=not ok
    print('zero-probe pass',n,'bad',bad)
if __name__=='__main__':
    impl=sys.argv[1] if len(sys.argv)>1 else 'lib'
    use(impl); print('implementation:',impl)
    main(int(sys.argv[2]) if len(sys.argv)>2 else 3000)
