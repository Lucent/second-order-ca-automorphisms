from base import *
def vcell(phi,i): return 2*Vac(phi-1,i)+Vac(phi,i)
def evolve_list(cells,n,pad):
    """cells dict-free: list over columns lo..hi; return list shrunk by n each side (exact)"""
    w=list(cells)
    for _ in range(abs(n)): w=Fstep_word(w) if n>0 else Finv_word(w)
    return w
def Kminus0(i): return vcell(0,i) if i<=1 else 0
def Kplus0(i):  return vcell(0,i) if i<=0 else 1
if __name__=='__main__':
    lo,hi=-40,40
    for name,K,s in (('K-',Kminus0,1),('K+',Kplus0,-1)):
        w=[K(i) for i in range(lo,hi+1)]
        w3=evolve_list(w,3,0)  # columns lo+3..hi-3
        ok=all(w3[j]==K(lo+3+j+s) for j in range(len(w3)))
        print(name,'F^3 = sigma^%d'%s,ok)
    # rho images
    RHO=[0,2,1,3]
    # rho K- should be F^n K+ for some n (up to the family); test by matching on window
    def fam(K,n,lo,hi):  # F^n K restricted, computed by evolving a wide window
        m=abs(n)+2
        w=[K(i) for i in range(lo-m,hi+m+1)]
        w=evolve_list(w,n,0)
        return w[m-abs(n):m-abs(n)+hi-lo+1]
    for name,K,T in (('rho K-',Kminus0,Kplus0),('rho K+',Kplus0,Kminus0)):
        target=[RHO[K(i)] for i in range(-20,21)]
        hits=[n for n in range(-12,13) if fam(T,n,-20,20)==target]
        print(name,'= F^n of other kink for n in',hits)
