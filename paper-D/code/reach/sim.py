import numpy as np
W=720
def mk(rule): return np.array([(rule>>t)&1 for t in range(8)],np.int8)
def step(ft,X,Y,n):
    for _ in range(abs(n)):
        if n>0:
            X,Y=Y,X^ft[(np.roll(Y,1)<<2)|(Y<<1)|np.roll(Y,-1)]
        else:
            X,Y=Y^ft[(np.roll(X,1)<<2)|(X<<1)|np.roll(X,-1)],X
    return X,Y
def Fst(rule,s):  # periodic state tuple
    xs,ys=s;p=len(xs)
    return ys,tuple(xs[i]^((rule>>(4*ys[i-1]+2*ys[i]+ys[(i+1)%p]))&1) for i in range(p))
def bg_orbit(rule,s):
    orb=[s];t=Fst(rule,s)
    while t!=s: orb.append(t);t=Fst(rule,t)
    return orb
def fill_bg(X,Y,lo,hi,s):
    xs,ys=s;p=len(xs)
    for i in range(lo,hi): X[i]=xs[i%p];Y[i]=ys[i%p]
def is_bg(X,Y,i,s):
    xs,ys=s;p=len(xs); return X[i]==xs[i%p] and Y[i]==ys[i%p]
