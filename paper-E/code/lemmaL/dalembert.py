import random
random.seed(0)
def check(rule,trials=200):
    f=lambda a,b,c:(rule>>(4*a+2*b+c))&1
    N=lambda a,b,c: f(a,b,c)^a^c
    bad=0; tested=0
    for _ in range(trials):
        W=60; T=25
        Y={}
        for i in range(W): Y[(-1,i)]=random.randint(0,1); Y[(0,i)]=random.randint(0,1)
        for t in range(0,T):
            for i in range(t+1,W-t-1): Y[(t+1,i)]=Y[(t-1,i)]^f(Y[(t,i-1)],Y[(t,i)],Y[(t,i+1)])
        g=lambda q,p: Y.get(((q+p)//2,(p-q)//2))
        # corner
        q0=random.randint(-32,-26); p0=-q0+2*random.randint(0,2)
        if g(q0,p0) is None: continue
        A=lambda q: g(q,p0) if (q-p0)%2==0 else g(q,p0+1)
        B=lambda p: g(q0,p) if (p-q0)%2==0 else g(q0+1,p)
        kap=lambda q: g(q0,p0) if (q-q0)%2==0 else g(q0+1,p0+1)
        for q in range(q0,q0+12):
            for p in range(p0,p0+12):
                if (q-p)%2: continue
                if g(q,p) is None: continue
                cnt=0; okk=True
                for qc in range(q0+1,q):
                    for pc in range(p0+1,p):
                        if (qc-pc)%2 or (qc-q)%2==0: continue
                        W_,C,E=g(qc+1,pc-1),g(qc,pc),g(qc-1,pc+1)
                        if None in (W_,C,E): okk=False;break
                        cnt^=N(W_,C,E)
                    if not okk: break
                if not okk or None in (A(q),B(p),kap(q)): continue
                tested+=1
                if g(q,p)!=A(q)^B(p)^kap(q)^cnt: bad+=1
    return bad,tested
for r in (91,94,95,155): print(r,'formula violations:',check(r))
