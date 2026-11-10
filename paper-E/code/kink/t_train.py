import sim, numpy as np
T=90; L=2*T+20
def kinks(s): return [x for x in range(1,len(s)) if s[x]!=s[x-1]]
res={}
for beta in (0,1):
  cases=[];meta=[]
  for a in (8,9):
    for b in (8,9):
      for g in list(range(1,9))+[15,22,31]:
        PI=[beta if q<a else 1-beta for q in range(L)]
        XI=[beta if (p<b or p>=b+g) else 1-beta for p in range(L)]
        cases.append((np.array(PI,np.uint8),np.array(XI,np.uint8))); meta.append((a,b,g))
  cur,off,f=sim.run(94,[c[0] for c in cases],[c[1] for c in cases],T,(beta,beta))
  G,O=sim.outputs(cur,off,T,60,60)
  for c,(a,b,g) in enumerate(meta):
    d = int(a%2==(b+g)%2) if beta==0 else int(a%2==b%2)
    ok = kinks(O[c])==[a+d] and O[c][0]==beta
    if (not ok and not (beta==0 and g==1)) or g in (3,31): print('beta',beta,'a',a,'b',b,'g',g,'pred delta',d,'O kinks',kinks(O[c]),'O front',O[c][0],'Gamma kinks',kinks(G[c]),'ok',ok)
