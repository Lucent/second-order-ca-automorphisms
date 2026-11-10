// Exhaustive check of the defect-climbing lemma for N=[000] (91) and N=[010] (94). usage: climb L RULE
// Word Z on positions 0..L-1 (p0=0). Step with updated class u in {0,1}: positions p≡u.
// D(p_b)=0 at lowest p_b of class; D(p)=D(p-2)^N(Z'(p-2),Z(p-1),Z(p)); Z'(p)=Z(p)^D(p) on class, Z' = Z elsewhere.
// Check (i) Z'=Z below h, h'>=h ; (iii) firings at p>=h ; (ii) h''>=h+3 (two steps, classes u then 1-u).
#include <stdio.h>
#include <stdlib.h>
int L, Nmask;
static int Nf(int a,int b,int c){return (Nmask>>(4*a+2*b+c))&1;}
static int low(const int*Z){for(int p=2;p<L;p++) if(Nf(Z[p-2],Z[p-1],Z[p])) return p; return 1<<20;}
static int step(const int*Z,int*Zn,int u,int*minfire){ int D=0; *minfire=1<<20;
  for(int p=0;p<L;p++) Zn[p]=Z[p];
  int pb=u; D=0; // D_bot=0
  for(int p=pb+2;p<L;p+=2){ int n=Nf(Zn[p-2],Z[p-1],Z[p]); if(n && p<*minfire) *minfire=p; D^=n; Zn[p]=Z[p]^D; }
  return 0;}
int main(int argc,char**argv){ L=atoi(argv[1]); int rule=atoi(argv[2]); Nmask = rule==91? 1 : 4; // [000] -> bit0 ; [010] -> bit 2
  long viol_i=0,viol_iii=0,viol_ii=0,cnt=0; int Z[64],Z1[64],Z2[64];
  for(long w=0;w<(1L<<L);w++){ for(int p=0;p<L;p++) Z[p]=(w>>p)&1;
    int h=low(Z);
    for(int u=0;u<2;u++){ int mf,mf2; step(Z,Z1,u,&mf); int h1=low(Z1);
      for(int p=0;p<L&&p<h;p++) if(Z1[p]!=Z[p]) {viol_i++;break;}
      if(h1<h) viol_i++;
      if(mf<h) viol_iii++;
      step(Z1,Z2,1-u,&mf2); int h2=low(Z2);
      if(h+3<=L-1 && h2<h+3) viol_ii++;
      cnt++; }}
  printf("rule %d L=%d steps %ld  viol(i)=%ld viol(iii)=%ld viol(ii)=%ld\n",rule,L,cnt,viol_i,viol_iii,viol_ii);}
