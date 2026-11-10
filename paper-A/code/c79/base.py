import numpy as np
RULE=93
def f(a,b,c): return (RULE>>(4*a+2*b+c))&1
FT=np.array([(RULE>>t)&1 for t in range(8)],np.uint8)
# vacuum diagram V(t,i), period 6 in t and 2 in i (checked in check_kinks.py)
def Vac(t,i):
    m=t%6
    return [i&1, i&1, 0, 1-(i&1), 1-(i&1), 0][m]
def Fstep_word(w):
    """w: tuple of cells (2x+y), returns F(w) of length len-2"""
    out=[]
    for i in range(1,len(w)-1):
        x,y=w[i]>>1,w[i]&1
        yl,yr=w[i-1]&1,w[i+1]&1
        out.append(2*y+(x^f(yl,y,yr)))
    return out
def Finv_word(w):
    out=[]
    for i in range(1,len(w)-1):
        x,y=w[i]>>1,w[i]&1
        xl,xr=w[i-1]>>1,w[i+1]>>1
        out.append(2*(y^f(xl,x,xr))+x)
    return out
def stepF(X,Y):  # numpy rows along axis0 columns; non periodic: returns shrunk
    return Y[1:-1], X[1:-1]^FT[(Y[:-2]<<2)|(Y[1:-1]<<1)|Y[2:]]
def stepB(X,Y):
    return Y[1:-1]^FT[(X[:-2]<<2)|(X[1:-1]<<1)|X[2:]], X[1:-1]
