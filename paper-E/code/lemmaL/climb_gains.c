// Two-step gains of the climbing lemmas, computed with the Goursat recurrence. usage: climb_gains RULE W G
// Y[q][p], q in {0..3}, p in [0,W). Column word Z_k(p)=Y[k+((p-k)&1)][p].
// Start: Z_0 arbitrary (Y[0][even p], Y[1][odd p]); bottom D=0: Y[2][0]=Y[0][0], Y[3][1]=Y[1][1].
// Option e0: first step updates even class (k=0). For e0=1, we shift: use k=1 start (Z_1 given).
#include <stdio.h>
#include <stdlib.h>
int N[8];
int W;
int lowdef(int *Z){ for(int p=2;p<W;p++) if(N[4*Z[p-2]+2*Z[p-1]+Z[p]]) return p; return W; }
int main(int argc,char**argv){
  int r=atoi(argv[1]); W=atoi(argv[2]); int g=atoi(argv[3]);
  for(int k=0;k<8;k++){int a=k>>2,c=k&1; N[k]=((r>>k)&1)^a^c;}
  int worst2=1000, worst1=1000; long viol=0, cnt=0;
  for(int k0=0;k0<2;k0++)
  for(long w=0;w<(1L<<W);w++){
    static int Y[6][64]; 
    // rows q=k0, k0+1 hold Z_{k0}
    for(int p=0;p<W;p++){ int b=(w>>p)&1; int q=k0+((p-k0)&1); Y[q][p]=b; }
    // need Y[k0+((p-k0)&1)] only. compute rows q=k0+2, k0+3
    for(int q=k0+2;q<=k0+3;q++){
      int pb = (q&1); // least p with p≡q, p>=0
      Y[q][pb]=Y[q-2][pb];
      for(int p=pb+2;p<W;p+=2){ int n=N[4*Y[q][p-2]+2*Y[q-1][p-1]+Y[q-2][p]]; Y[q][p]=Y[q][p-2]^Y[q-2][p-2]^Y[q-2][p]^n; }
    }
    int Z0[64],Z1[64],Z2[64];
    for(int p=0;p<W;p++){ Z0[p]=Y[k0+((p-k0)&1)][p]; Z1[p]=Y[k0+1+((p-k0-1)&1)][p]; Z2[p]=Y[k0+2+((p-k0)&1)][p]; }
    int m0=lowdef(Z0), m1=lowdef(Z1), m2=lowdef(Z2);
    if(m0>=W-g-1) continue; cnt++;
    int g1=(m1<W? m1-m0 : 999), g2=(m2<W? m2-m0: 999);
    if(m1<m0) viol++;
    if(g1<worst1) worst1=g1; if(g2<worst2) worst2=g2;
  }
  printf("rule %d W=%d tested %ld: min 1-step gain %d, min 2-step gain %d, monotone violations %ld\n",r,W,cnt,worst1,worst2,viol);
}
