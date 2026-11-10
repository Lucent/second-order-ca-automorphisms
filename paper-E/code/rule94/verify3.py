# Closed forms of the exchange lemma and of the crossing of two isolated zeros, checked on the whole quadrant.
# usage: python3 verify3.py
from lib import solve
bad=0;n=0
M=40
for a in range(2,10):
    for d in range(2,10):
        if (a-d)%2==0:
            A=[1^(1 if q>=a else 0) for q in range(M)]; B=[1^(1 if p==d else 0) for p in range(M)]
            Y,F=solve(A,B)
            assert sorted(F)==[(a-1,d-1),(a+2,d+2)]
            for (q,p),v in Y.items():
                x=q-a;y=p-d
                if x<0: e=0 if y==0 else 1
                elif x%2==0: e=1 if y>=2 else 0
                else: e=1 if (x>=3 and y>=3) else 0
                n+=1; bad+=(e!=v)
                if p>=d+2:   # equals the outgoing probe: 1 except at q=a+1
                    bad+=(v!=(0 if q==a+1 else 1))
                if q>=a+3:   # equals the outgoing object: theta_{d+2}
                    bad+=(v!=(1 if p>=d+2 else 0))
        else:
            A=[1^(1 if q==a else 0) for q in range(M)]; B=[1^(1 if p==d else 0) for p in range(M)]
            Y,F=solve(A,B)
            assert F==[]
            for (q,p),v in Y.items():
                n+=1; bad+=(v!=(1^(q==a)^(p==d)))
print('points checked',n,'bad',bad)
# Lemma (isolated zeros)(ii): pair (Z^+_a, barK^-_u), a ≡ u: diagram is 1 xor [p=u+1] wherever q>=a+2
bad=0;n=0
for a in range(2,10):
    for u in range(2,10):
        if (a-u)%2: continue
        A=[1^(1 if q==a else 0) for q in range(M)]; B=[1^(1 if p>=u else 0) for p in range(M)]
        Y,F=solve(A,B)
        for (q,p),v in Y.items():
            if q>=a+2:
                n+=1; bad+=(v!=(1^(p==u+1)))
print('part (ii): points checked',n,'bad',bad)
