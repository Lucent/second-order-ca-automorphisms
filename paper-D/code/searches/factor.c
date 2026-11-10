/* Exact search for factor maps u (radius r, any alphabet) with u o F^p = sigma^k o u  (mode 0)
   or binary u with u o F^p = NOT sigma^k o u (mode 1, parity union-find).
   Table entries = (2r+1)-words. For every window w of length 2r+1+2p: u(F^p(w)) = [NOT] u(w[p+k .. p+k+2r+1)).
   Output: number of equivalence classes (mode 0) [>=2 iff a nonconstant factor exists], consistency (mode 1). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static int rule;
static inline int f(int a,int b,int c){return (rule>>((a<<2)|(b<<1)|c))&1;}
static inline int Floc(int s0,int s1,int s2){int x=s1>>1,y=s1&1;return (y<<1)|(x^f(s0&1,y,s2&1));}
static int32_t *par; static uint8_t *pp;
static int32_t find(int32_t a, int *p){int q=0; int32_t r=a; while(par[r]!=r){q^=pp[r]; r=par[r];}
  /* path compression */ int32_t x=a; int qx=q; while(par[x]!=x){int32_t nx=par[x]; int nq=qx^pp[x]; par[x]=r; pp[x]=qx; x=nx; qx=nq;} *p=q; return r;}
int main(int argc,char**argv){
  rule=atoi(argv[1]); int r=atoi(argv[2]); int p=atoi(argv[3]); int k=atoi(argv[4]); int mode=atoi(argv[5]);
  int n=2*r+1, m=n+2*p; int64_t N=1LL<<(2*n), M=1LL<<(2*m);
  par=malloc(N*4); pp=calloc(N,1); for(int64_t i=0;i<N;i++)par[i]=i;
  int buf[64]; int bad=0;
  for(int64_t w=0; w<M && !bad; w++){
    for(int i=0;i<m;i++) buf[i]=(w>>(2*(m-1-i)))&3;
    int len=m; for(int t=0;t<p;t++){ for(int i=0;i<len-2;i++) buf[i]=Floc(buf[i],buf[i+1],buf[i+2]); len-=2; }
    int64_t a=0; for(int i=0;i<n;i++) a=(a<<2)|buf[i];
    int64_t b=(w>>(2*(m-(p+k)-n)))&(N-1);
    int pa,pb; int32_t ra=find(a,&pa), rb=find(b,&pb);
    int want = mode; /* parity between a and b */
    if(ra==rb){ if((pa^pb)!=want) bad=1; }
    else { par[ra]=rb; pp[ra]=pa^pb^want; }
  }
  if(mode==1){ printf("rule %d r %d p %d k %d NOT-mode: %s\n",rule,r,p,k,bad?"inconsistent (no factor)":"CONSISTENT (binary factor exists)"); return 0;}
  int64_t cl=0; for(int64_t i=0;i<N;i++) if(par[i]==i) cl++;
  printf("rule %d r %d p %d k %d classes %lld%s\n",rule,r,p,k,(long long)cl, cl>1?"  <-- NONTRIVIAL FACTOR":"");
  return 0;}
