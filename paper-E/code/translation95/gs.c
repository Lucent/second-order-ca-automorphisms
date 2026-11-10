// Goursat light-cone solver. Y[q][p] on 0<=q,p<M, q==p mod 2.
// boundary: Y[q][q%2]=A[q], Y[p%2][p]=B[p]. rule given by N-mask over (W,C,E) index 4W+2C+E.
#include <stdlib.h>
int gsolve(int nmask,const unsigned char*A,const unsigned char*B,int M,unsigned char*G,unsigned char*O,int*fq,int*fp,int maxf){
  unsigned char *Y=malloc((size_t)M*M);
  #define Y_(q,p) Y[(size_t)(q)*M+(p)]
  for(int k=0;k<M;k++){Y_(k,k%2)=A[k]; Y_(k%2,k)=B[k];}
  int cnt=0;
  for(int s=4;s<2*M;s+=2) for(int q=2;q<M;q++){ int p=s-q; if(p<2||p>=M) continue;
     int W=Y_(q,p-2),S=Y_(q-2,p-2),E=Y_(q-2,p),C=Y_(q-1,p-1);
     int n=(nmask>>(4*W+2*C+E))&1; Y_(q,p)=W^S^E^n;
     if(n){ if(cnt<maxf){fq[cnt]=q-1;fp[cnt]=p-1;} cnt++; }
  }
  for(int p=0;p<M;p++){int q=M-1; if((q^p)&1) q--; G[p]=Y_(q,p);}
  for(int q=0;q<M;q++){int p=M-1; if((q^p)&1) p--; O[q]=Y_(q,p);}
  free(Y); return cnt;
}
