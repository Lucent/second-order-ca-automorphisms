# The edge bit of 155R: Y(t,s-t)=Y(0,s) behind a left edge into a valid background, and the 11-front fires at every step.
import random
RULE=155
f=lambda a,b,c:(RULE>>(4*a+2*b+c))&1
N=lambda a,b,c:f(a,b,c)^a^c
print("N support",[ (a,b,c) for a in (0,1) for b in (0,1) for c in (0,1) if N(a,b,c)])
for beta in (0,1):
    print("beta",beta,"edge map c->(1-b)^f(1-b,b,c):",[(1-beta)^f(1-beta,beta,c) for c in (0,1)])
# simulation on a finite line with light-cone shrinking (no wrap): Y rows as lists
def run(x,y,T):
    rows=[x,y]
    for t in range(T):
        a=rows[-2]; b=rows[-1]; n=len(b)
        new=[None]*n
        for i in range(1,n-1):
            new[i]=a[i]^f(b[i-1],b[i],b[i+1])
        # boundaries unknown: we treat ends as invalid by trimming
        rows.append(new)
    return rows
viol=0; forever=0; trials=0
for seed in range(3000):
    rng=random.Random(seed)
    T=60; L=4*T+40; s=2*T+10
    phi=rng.randint(0,1)
    bg=lambda t,i:(t+i+phi)%2
    x=[bg(-1,i) for i in range(L)]; y=[bg(0,i) for i in range(L)]
    # perturbation: y_s arbitrary, x_s arbitrary, cells >s arbitrary
    for j in range(s,L):
        x[j]=rng.randint(0,1); y[j]=rng.randint(0,1)
    # exact computation with free right boundary: compute full Y in light cone of known cells: use extended right randomness, trust only i<=L-1-t
    Y={(-1,i):x[i] for i in range(L)}; Y.update({(0,i):y[i] for i in range(L)})
    for t in range(0,T):
        for i in range(1+t, L-1-t):
            Y[(t+1,i)]=Y[(t-1,i)]^f(Y[(t,i-1)],Y[(t,i)],Y[(t,i+1)])
    trials+=1
    ok=True
    for t in range(0,T):
        # left of edge is background
        for i in range(1+t, s-t):
            if Y[(t,i)]!=bg(t,i): ok=False
        if Y[(t,s-t)]!=y[s]: ok=False
    if not ok: viol+=1
    beta=bg(0,s-1)
    if y[s]==beta and beta==1:
        if all(N(Y[(t,s-t-1)],Y[(t,s-t)],Y[(t,s-t+1)])==1 for t in range(T)): forever+=1
        else: viol+=1000
print("trials",trials,"violations",viol,"11-front cases firing every step",forever)
