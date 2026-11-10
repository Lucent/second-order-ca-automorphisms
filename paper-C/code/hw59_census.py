"""Rule 59R, Proposition 6.16: builds, by exact composition on all windows of radius up to 4, the tables of
(i) the 13 maps sigma^a F^b with |a|+|b| <= 2, (ii) the 75 maps F^b t with |b| <= 2 and t a nontrivial
element of T_R, (iii) the 45 maps sigma F^b t with |b| <= 1 and t a nontrivial element of T_L; checks that
each has radius 2, that the 133 tables are pairwise distinct, and that they are exactly the 133 invertible
maps of the census.
Expected output: distinct tables found 133, matched solutions 133, unmatched solutions 0, candidates not in
census 0.
"""
from hw59_census_lib import *
import sys
found={}
def add(name,tab):
    if tab is None: print("NOT radius 2:",name); return
    found.setdefault(tab,[]).append(name)
def Fb(b):
    return [ (lambda w: Fmap(w)) ]*b if b>0 else [ (lambda w: Fmap(w,True)) ]*(-b)
def Sa(a):
    return [Sig]*a if a>0 else [Sinv]*(-a)
# (i)
for a in range(-2,3):
    for b in range(-2,3):
        if abs(a)+abs(b)<=2:
            add(('S',a,'F',b),table_of(Sa(a)+Fb(b),max(2,abs(a)+abs(b))))
# (ii) F^b t, t in T_R nontrivial, t as union flip
for b in range(-2,3):
    for sub in itertools.product([0,1],repeat=4):
        if not any(sub): continue
        cls=[g for g,s in zip(GR,sub) if s]
        add(('F',b,'TR',sub),table_of([lambda w,cls=cls: flipgen(w,cls)]+Fb(b),3+abs(b)))
# (iii) sigma F^b t, t in T_L nontrivial
for b in range(-1,2):
    for sub in itertools.product([0,1],repeat=4):
        if not any(sub): continue
        cls=[g for g,s in zip(GL,sub) if s]
        add(('S',1,'F',b,'TL',sub),table_of([lambda w,cls=cls: flipgen(w,cls)]+Fb(b)+[Sig],4+abs(b)))
print("distinct tables found",len(found), "multi:",[v for v in found.values() if len(v)>1])
inS=[k for k in found if k in soltabs]; print("matched solutions",len(inS))
miss=[s for s in soltabs if s not in found]; print("unmatched solutions",len(miss),[soltabs[s] for s in miss])
notsol=[found[k] for k in found if k not in soltabs]; print("candidates not in census",len(notsol),notsol[:5])
