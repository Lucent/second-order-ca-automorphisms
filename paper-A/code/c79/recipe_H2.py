# The proof of Proposition (H_2) as a checked recipe: h2(c,t) follows the steps of the proof for the time t and
# asserts every rule (G1)-(G3) and every block type it uses. c provides th, z, y, n, blk, move (recipe_states.py).
def G1(c,t): assert c.n(t) in (0,1) and c.z(t)==1
def G2(c,mid,a,b): assert c.z(mid)==0 and {a,b}=={mid-1,mid+1}
def G3(c,m,t): assert m-2<=t<=m+1; return c.move(m)[0]
def h2(c,t):
    (p,q),j=c.blk(t); T=lambda k: t-j+k
    if (p,q)==(1,1):
        if j in (1,6): G1(c,t); return
        if j==2: D=G3(c,T(2),t); assert D.n(t)==1 and D.z(t)==1; G1(D,t); return     # F- at t_2
        if j==4: G2(c,T(3),t,T(2)); h2(c,T(2)); return
        if j==5: D=G3(c,T(6),t); assert D.n(t)==0 and D.z(t)==1; G1(D,t); return     # B+ at t_6
        if j==3: G2(c,T(4),t,T(5)); h2(c,T(5)); return
    if (p,q)==(0,1):
        if j==6: G1(c,t); return
        if j==4: G2(c,T(5),t,T(6)); G1(c,T(6)); return
        if j==1:   # F+ at t_2: D-block of type (1,0) on [t_2-3,t_2+2], t_1 at position 3
            D=G3(c,T(2),t); assert D.blk(t)==((1,0),3); G2(D,t-1,t,t-2); G1(D,t-2); return
        if j==2: G2(c,T(1),t,T(0)); assert c.blk(T(0))[0][1]==1; G1(c,T(0)); return
        if j==3: G2(c,T(4),t,T(5)); h2(c,T(5)); return
        if j==5:   # B+ at t_6: t_5 at position 6 of a D-block of type (1,0)
            D=G3(c,T(6),t); assert D.blk(t)==((1,0),6); h2(D,t); return
    if (p,q)==(1,0):
        if j==1: G1(c,t); return
        if j==3: G2(c,T(2),t,T(1)); G1(c,T(1)); return
        if j==6:   # B- at t_6: t_6 at position 4 of a D-block of type (0,1)
            D=G3(c,T(6),t); assert D.blk(t)==((0,1),4); h2(D,t); return
        if j==2:   # F- at t_2: t_2 at position 1 of a D-block of type (0,1)
            D=G3(c,T(2),t); assert D.blk(t)==((0,1),1); h2(D,t); return
        if j==4: G2(c,T(3),t,T(2)); h2(c,T(2)); return
        if j==5: G2(c,T(6),t,T(7)); assert c.blk(T(7))==((1,0) if False else c.blk(T(7))[0],1) and c.z(T(7))==1; G1(c,T(7)); return
    raise Exception((p,q,j))
