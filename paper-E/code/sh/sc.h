// Scattering for second-order ECA F(x,y)=(y, x^f(y)), direct (t,i) simulation.
// Phi+ : Y(t,i)=phi(i-t) (moves right), Phi- : Y(t,i)=psi(i+t) (moves left).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int RULE;
static inline int fr(int a,int b,int c){return (RULE>>(4*a+2*b+c))&1;}
static inline int Nf(int a,int b,int c){return fr(a,b,c)^a^c;}
#define WMAX 4096
#define OFS 2048
// tails: 2-periodic, t[par] where par = j&1 (j may be negative -> use ((j%2)+2)%2)
static inline int par(int j){return j&1;}
typedef struct { int c[2], d[2]; } tails_t;
// phiP(j): probe on background c, core p[0..lp) at positions sP..sP+lp-1
// psiQ(j): object, core q[0..lq) at sQ.., left tail c, right tail d.
static unsigned char R0[WMAX],R1[WMAX],R2[WMAX];
// result arrays indexed by j+OFS
static unsigned char GO[WMAX], OO[WMAX]; static int MID[2];
static int scatter(const int*c,const int*d,const unsigned char*p,int lp,int sP,const unsigned char*q,int lq,int sQ,int T,int W){
  // rows over i in [-W,W]
  unsigned char *ym=R0,*y0=R1,*y1=R2; int n=2*W+1;
  for(int k=0;k<n;k++){int i=k-W; int vm,v0;
    if(i<0){ int j=i; v0=(j>=sP&&j<sP+lp)?p[j-sP]:c[par(j)]; j=i+1; vm=(j>=sP&&j<sP+lp)?p[j-sP]:c[par(j)]; }
    else { int j=i; v0=(j<sQ)?c[par(j)]:(j<sQ+lq?q[j-sQ]:d[par(j)]); j=i-1; vm=(j<sQ)?c[par(j)]:(j<sQ+lq?q[j-sQ]:d[par(j)]); }
    ym[k]=vm; y0[k]=v0; }
  int fires=0;
  for(int t=0;t<T;t++){
    // boundary k=0 and k=n-1 at time t+1: tails c(i+t+1), d(i+t+1)
    y1[0]=c[par(-W+t+1)]; y1[n-1]=d[par(W+t+1)];
    for(int k=1;k<n-1;k++){ y1[k]=ym[k]^fr(y0[k-1],y0[k],y0[k+1]); fires+=Nf(y0[k-1],y0[k],y0[k+1]); }
    unsigned char*tmp=ym; ym=y0; y0=y1; y1=tmp; }
  // now y0=Y(T), ym=Y(T-1). verify free state
  int ok=1;
  // mid tail
  int mi=W; // index of i=0
  MID[par(0+T)]=y0[mi]; MID[par(-1+T)]=y0[mi-1];
  for(int i=-4;i<4;i++) if(y0[mi+i]!=MID[par(i+T)]) ok=0;
  for(int i=-W+2;i<=-1;i++){ if(ym[mi+i]!=y0[mi+i-1]) ok=0; }
  for(int i=0;i<=W-2;i++){ if(ym[mi+i]!=y0[mi+i+1]) ok=0; }
  for(int k=1;k<n-1;k++) if(Nf(y0[k-1],y0[k],y0[k+1])) ok=0;
  // boundaries far: left must be c, right d
  for(int i=-W;i<-W+4;i++) if(y0[mi+i]!=c[par(i+T)]) ok=0;
  for(int i=W-3;i<=W;i++) if(y0[mi+i]!=d[par(i+T)]) ok=0;
  // extract psi_out(j)=Y(T,j-T) for j-T<=-1 ; phi_out(j)=Y(T,j+T) for j+T>=0; store for j in [-OFS,OFS)
  for(int j=-OFS;j<OFS;j++){ int i=j-T; GO[j+OFS]= (i<=-1 && i>=-W)? y0[mi+i] : (i>-1? MID[par(j)] : c[par(j)]);
                             i=j+T; OO[j+OFS]= (i>=0 && i<=W)? y0[mi+i] : (i<0? MID[par(j)] : d[par(j)]); }
  return ok?fires:-1;
}
// support of a sequence g (indexed j+OFS) wrt left tail L[] and right tail Rt[]
static int suppof(const unsigned char*g,const int*L,const int*Rt){
  int a=-OFS; while(a<OFS && g[a+OFS]==L[par(a)]) a++;
  int e=OFS-1; while(e>=-OFS && g[e+OFS]==Rt[par(e)]) e--;
  if(a>e) return 0; return e-a+1; }
static int okseq(const unsigned char*s,int n){ for(int k=0;k+2<n;k++) if(Nf(s[k],s[k+1],s[k+2])) return 0; return 1;}
