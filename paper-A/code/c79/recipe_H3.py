# The proof of Proposition (H_3) as a checked recipe: prove(c,t) follows Steps 1-4 of the proof (Y1 is Step 1,
# Y0 Steps 2-4) and asserts every rule (K1)-(K3), identity (Y), (Z), (W) and block type it uses.
import sys
sys.setrecursionlimit(10000)
cnt={}
def tick(k): cnt[k]=cnt.get(k,0)+1
def R3(c,m,t):
    assert m-3<=t<=m+2; return c.move(m)[0]
def Y1(c,t):
    """y_t=1  ->  g_t=0 : z-control by moves at m in [t-1,t+2]"""
    assert c.y(t)==1
    steps=0
    while c.z(t)==0:
        steps+=1; assert steps<=4
        (p,q),pos=c.blk(t); tick(('Y1',p,q,pos))
        if pos==3: m=t-1          # F- (p=1)
        elif pos==2: m=t          # [P10,2]: F-
        elif pos==1: m=t+1        # [P01,1]: F+
        elif pos==6: m=t          # [P10,6]: B-
        elif pos==4: m=t+2        # B+
        elif pos==5: m=t+1        # [P01,5]: B+
        assert t-1<=m<=t+2
        c=R3(c,m,t); assert c.y(t)==1
    tick('R1')   # y_t=z_t=1
def Y0(c,t,ref):
    assert c.y(t)==0
    if c.z(t)==0:                      # Case A
        assert c.y(t+1)!=c.y(t-1)
        d=1 if c.y(t+1)==0 else -1     # (rho)
        tick(('A',c.z(t+d)))
        if c.z(t+d)==0:
            assert c.y(t+2*d)==1; Y1(c,t+2*d)           # R2 at t+d
        else:
            assert c.y(t+2*d)==0 and c.z(t+2*d)==1 and c.z(t+3*d)==0 and c.y(t+3*d)==0 and c.y(t+4*d)==1
            Y1(c,t+4*d)                                  # R2 at t+d and t+3d
        return
    assert c.y(t+1)==c.y(t-1) and c.z(t+1)!=c.z(t-1) and c.w(t+1)!=c.w(t-1)
    if c.y(t+1)==0:                    # Case B
        d=1 if c.z(t+1)==0 else -1; tick('B')
        assert c.y(t+2*d)==1; Y1(c,t+2*d); return
    # Case C
    (p,q),pos=c.blk(t)
    if c.z(t-1)==0:   # orientation z_{t-1}=0, z_{t+1}=1
        key=((p,q),pos); tick(('C',key))
        if key==((1,1),1): D=R3(c,t-1,t); assert D.z(t)==0
        elif key==((0,1),2): D=R3(c,t-2,t); assert D.z(t)==0
        elif key==((1,1),5): D=R3(c,t+3,t); assert D.z(t)==0
        elif key==((0,1),6): D=R3(R3(c,t,t),t-1,t); assert D.z(t)==0
        elif key==((1,0),4):
            D=R3(c,t-2,t)
            if c.w(t-1)==1: assert D.y(t)==1
            else: D=R3(D,t+3,t); assert D.y(t)==1
        else: raise Exception(key)
    else:             # mirror orientation (rho images): z_{t-1}=1, z_{t+1}=0
        key=((p,q),pos); tick(('Cm',key))
        if key==((1,1),6): D=R3(c,t+2,t); assert D.z(t)==0
        elif key==((1,0),5): D=R3(c,t+3,t); assert D.z(t)==0
        elif key==((1,1),2): D=R3(c,t-2,t); assert D.z(t)==0
        elif key==((1,0),1): D=R3(R3(c,t+1,t),t+2,t); assert D.z(t)==0
        elif key==((0,1),3):
            D=R3(c,t+3,t)
            if c.w(t+1)==1: assert D.y(t)==1
            else: D=R3(D,t-2,t); assert D.y(t)==1
        else: raise Exception(key)
    prove(D,t)
def prove(c,t):
    if c.y(t)==1: Y1(c,t)
    else: Y0(c,t,c)

