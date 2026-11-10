/* Rule 219R, after Proposition 6.24: exact verification of the 2-car train swap Psi_2 for 219R.
   Cells s = 2x+y.  F(c)_i = (y_i, x_i ^ f(y_{i-1},y_i,y_{i+1})), f = rule 219.
   Occurrence at first car p:  c[p-3..p+6] = s s s a s s b s s s  with a,b != s.
   Psi(c)_p = c_{p+3}, Psi(c)_{p+3} = c_p for every occurrence p; identity elsewhere.  Radius 6.
   T1: Psi o F = F o Psi on ALL 4^15 windows (exact).
   T2: occurrence lemma  occ_q(Psi c) == occ_q(c)  on all 22-cell windows having an occurrence with a car in
       [q-3,q+6] (all other windows are trivial since Psi is the identity on [q-3,q+6] there).  This gives
       Psi o Psi = id exactly.
   T3: on all periodic configurations of period n<=11: Psi o Psi = id, F Psi = Psi F; count nontrivial. */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#define RULE 219
static inline int Floc(int l,int c,int r){ int x=c>>1,y=c&1; int f=(RULE>>(((l&1)<<2)|(y<<1)|(r&1)))&1; return (y<<1)|(x^f); }
static inline int occ(const uint8_t*w,int p){ int s=w[p-1];
  return w[p-3]==s && w[p-2]==s && w[p]!=s && w[p+1]==s && w[p+2]==s && w[p+3]!=s && w[p+4]==s && w[p+5]==s && w[p+6]==s; }
/* psi on a 13-window starting at w (centre w[6]) */
static inline int psi13(const uint8_t*w){ if(occ(w,6)) return w[9]; if(occ(w,3)) return w[3]; return w[6]; }
int main(){
  /* ---------- T1 */
  long long bad1=0, nontriv1=0; const long long N15=1LL<<30;
  #pragma omp parallel for reduction(+:bad1,nontriv1) schedule(static)
  for(long long e=0;e<N15;e++){
    uint8_t v[15],Fv[13]; long long t=e; for(int k=14;k>=0;k--){v[k]=t&3;t>>=2;}
    for(int k=0;k<13;k++) Fv[k]=Floc(v[k],v[k+1],v[k+2]);
    int lhs=psi13(Fv);
    int g0=psi13(v),g1=psi13(v+1),g2=psi13(v+2);
    int rhs=Floc(g0,g1,g2);
    if(lhs!=rhs) bad1++;
    if(g1!=v[7]) nontriv1++;
  }
  printf("T1 commutation on all 4^15 windows: failures %lld (windows where Psi changes the centre: %lld)\n",bad1,nontriv1);
  /* ---------- T2 */
  long long bad2=0, tested2=0, changed2=0;
  for(int p=3;p<=15;p++) for(int s=0;s<4;s++) for(int a=0;a<4;a++) for(int b=0;b<4;b++){
    if(a==s||b==s) continue;
    int fixedpos[10]={p-3,p-2,p-1,p,p+1,p+2,p+3,p+4,p+5,p+6}; int fixedval[10]={s,s,s,a,s,s,b,s,s,s};
    int isfixed[22]; memset(isfixed,0,sizeof isfixed); int val[22];
    for(int k=0;k<10;k++){isfixed[fixedpos[k]]=1; val[fixedpos[k]]=fixedval[k];}
    int freepos[12],nf=0; for(int k=0;k<22;k++) if(!isfixed[k]) freepos[nf++]=k;
    #pragma omp parallel for reduction(+:bad2,tested2,changed2) schedule(static)
    for(long long e=0;e<(1LL<<24);e++){
      uint8_t u[22],g[22]; long long t=e;
      for(int k=0;k<22;k++) u[k]=val[k]*isfixed[k];
      for(int k=0;k<nf;k++){u[freepos[k]]=t&3;t>>=2;}
      memcpy(g,u,22);
      for(int j=6;j<=15;j++) g[j]=psi13(u+j-6);
      int o1=occ(u,9), o2=occ(g,9);
      if(o1!=o2) bad2++;
      tested2++; if(memcmp(g+6,u+6,10)) changed2++;
    }
  }
  printf("T2 occurrence lemma on structured 22-windows: tested %lld, Psi nontrivial on [q-3,q+6] in %lld, failures %lld\n",tested2,changed2,bad2);
  /* ---------- T3 */
  for(int n=1;n<=11;n++){
    long long M=1LL<<(2*n), bad=0, nt=0;
    #pragma omp parallel for reduction(+:bad,nt) schedule(static)
    for(long long e=0;e<M;e++){
      uint8_t c[16],g[16],gg[16],Fc[16],gF[16],Fg[16],w[13]; long long t=e;
      for(int k=0;k<n;k++){c[k]=t&3;t>>=2;}
      #define APPLY(src,dst) for(int i=0;i<n;i++){ for(int k=0;k<13;k++) w[k]=src[((i+k-6)%n+n)%n]; dst[i]=psi13(w);}
      #define FAPP(src,dst) for(int i=0;i<n;i++) dst[i]=Floc(src[(i+n-1)%n],src[i],src[(i+1)%n]);
      APPLY(c,g); APPLY(g,gg); FAPP(c,Fc); APPLY(Fc,gF); FAPP(g,Fg);
      if(memcmp(gg,c,n)||memcmp(gF,Fg,n)) bad++;
      if(memcmp(g,c,n)) nt++;
    }
    printf("T3 period %2d: rings %lld, Psi nontrivial on %lld, failures (Psi^2!=id or F Psi != Psi F) %lld\n",n,M,nt,bad);
  }
  return 0;
}
