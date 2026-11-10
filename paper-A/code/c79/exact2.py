# Exact (two-valued) simulation of rule 93 on finite words with a vacuum to the right; used by the confirmations only.
import random
from engine2 import f2,V,vcell
def Fw(w): return [(w[i][1], w[i][0]^f2(w[i-1][1],w[i][1],w[i+1][1])) for i in range(1,len(w)-1)]
def Bw(w): return [(w[i][1]^f2(w[i-1][0],w[i][0],w[i+1][0]), w[i][0]) for i in range(1,len(w)-1)]
def pw(w,m):
    for _ in range(abs(m)): w=Fw(w) if m>0 else Bw(w)
    return w
class Conf:
    def __init__(s,th,left,PADR=60):  # left: list of cells on [-n..0]
        s.th=th; s.n=len(left)-1; s.w=list(left)+[vcell(th,i) for i in range(1,PADR+1)]
    def col(s,t,i):  # Y(t,i)
        off=-s.n
        if t>=0:
            w=pw(s.w,t); k=i-off-t; assert 0<=k<len(w) and k>=0; return w[k][1]
        w=pw(s.w,t+1); k=i-off+t+1; assert 0<=k<len(w); return w[k][0]
    def move(s,m):
        n=(m+s.th)%6
        if n==2: z=s.col(m-1,0); th2=(s.th+2)%6 if z==0 else (s.th-1)%6
        elif n==0: z=s.col(m,0); th2=(s.th-2)%6 if z==0 else (s.th+1)%6
        else: return None
        a=abs(m); off=-s.n
        wm=pw(s.w,m); o2=off+a
        wm=[wm[i-o2] if i<=0 else vcell(th2+m,i) for i in range(o2,o2+len(wm))]
        wb=pw(wm,-m); o3=off+2*a
        new=[wb[i-o3] for i in range(o3,1)]
        c=Conf(th2,new); 
        assert all(wb[i-o3]==vcell(th2,i) for i in range(1,len(wb)+o3-3*a))
        return c
