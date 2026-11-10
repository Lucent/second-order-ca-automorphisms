/* Window closure (diamond D=2, five-point stencil).
   Base: windows of Fix(sigma F) U Fix(sigma^-1 F) built from Phi_f (no triple in supp N, N=f^a^c).
   Window: cells j=0..W-1, x_j=Y(t-1,j), y_j=Y(t,j). Pass p uses forced set from passes < p (Jacobi).
   usage: dc RULE R  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int RULE,R,W,M; static int fb(int a,int b,int c){return (RULE>>(4*a+2*b+c))&1;}
static int Nf(int a,int b,int c){return fb(a,b,c)^a^c;}
static int *stamp;
#define OFF 8
static int Yd[8][64]; /* row t+4, col i+M */
#define Y(t,i) Yd[(t)+4][(i)+M]
static uint64_t wincode(int t,int i){ uint64_t X=0,Yv=0; for(int k=i-R;k<=i+R;k++){X=(X<<1)|Y(t-1,k);Yv=(Yv<<1)|Y(t,k);} return (X<<W)|Yv; }
int main(int argc,char**argv){ RULE=atoi(argv[1]); R=atoi(argv[2]); W=2*R+1; M=R+2; int D=2;
 uint64_t NW=1ULL<<(2*W); stamp=malloc(NW*sizeof(int)); for(uint64_t w=0;w<NW;w++) stamp[w]=-1;
 /* base */
 int n=W+1; long nb=0; int K=6, nn=n+2*K;
 for(uint64_t v=0;v<(1ULL<<nn);v++){ int ok=1; int bb[80]; for(int k=0;k<nn;k++) bb[k]=(v>>(nn-1-k))&1;
   for(int k=1;k+1<nn&&ok;k++) if(Nf(bb[k-1],bb[k],bb[k+1])) ok=0; if(!ok) continue;
   int *b=bb+K;
   uint64_t Xp=0,Yp=0,Xm=0,Ym=0; for(int j=0;j<W;j++){ Yp=(Yp<<1)|b[j]; Xp=(Xp<<1)|b[j+1]; Ym=(Ym<<1)|b[j+1]; Xm=(Xm<<1)|b[j]; }
   uint64_t a1=(Xp<<W)|Yp, a2=(Xm<<W)|Ym; if(stamp[a1]<0){stamp[a1]=0;nb++;} if(stamp[a2]<0){stamp[a2]=0;nb++;} }
 printf("rule %d R %d base %ld of %llu\n",RULE,R,nb,(unsigned long long)NW);
 int pts[13][2]; int np=0; for(int t=-D;t<=D;t++) for(int i=-D;i<=D;i++) if(abs(t)+abs(i)<=D){pts[np][0]=t;pts[np][1]=i;np++;}
 int cen[5][2]={{0,0},{1,0},{-1,0},{0,1},{0,-1}};
 for(int pass=1;;pass++){ long added=0, remaining=0;
  for(uint64_t w=0;w<NW;w++){ if(stamp[w]>=0) continue; remaining++;
   uint64_t Xw=w>>W, Yw=w&((1ULL<<W)-1); int forced=0;
   for(int e=0;e<256 && !forced;e++){
     for(int k=0;k<W;k++){ Y(-1,k-R)=(Xw>>(W-1-k))&1; Y(0,k-R)=(Yw>>(W-1-k))&1; }
     Y(-1,-M)=e&1; Y(-1,-M+1)=(e>>1)&1; Y(-1,M-1)=(e>>2)&1; Y(-1,M)=(e>>3)&1;
     Y(0,-M)=(e>>4)&1; Y(0,-M+1)=(e>>5)&1; Y(0,M-1)=(e>>6)&1; Y(0,M)=(e>>7)&1;
     for(int t=1;t<=2;t++) for(int i=-M+t;i<=M-t;i++) Y(t,i)=Y(t-2,i)^fb(Y(t-1,i-1),Y(t-1,i),Y(t-1,i+1));
     for(int t=-2;t>=-3;t--){ int s=-1-t; for(int i=-M+s;i<=M-s;i++) Y(t,i)=Y(t+2,i)^fb(Y(t+1,i-1),Y(t+1,i),Y(t+1,i+1)); }
     int z[5][5]; /* known-zero flag at (t,i) indexed t+2,i+2 */
     for(int q=0;q<np;q++){ int t=pts[q][0],i=pts[q][1]; if(t==0&&i==0) continue; int s=stamp[wincode(t,i)]; z[t+2][i+2]=(s>=0 && s<pass); }
     for(int c=0;c<5 && !forced;c++){ int ct=cen[c][0],ci=cen[c][1];
       int sp[5][2]={{ct+1,ci},{ct-1,ci},{ct,ci-1},{ct,ci},{ct,ci+1}}; int unk[5],nu=0,tg=-1;
       for(int k=0;k<5;k++){ if(sp[k][0]==0&&sp[k][1]==0) tg=k; else if(!z[sp[k][0]+2][sp[k][1]+2]) unk[nu++]=k; }
       int sat=0;
       for(int a=0;a<(1<<nu)&&!sat;a++){ int g[5]={0,0,0,0,0}; g[tg]=1; for(int u=0;u<nu;u++) g[unk[u]]=(a>>u)&1;
         int yw=Y(sp[2][0],sp[2][1]),yc=Y(sp[3][0],sp[3][1]),ye=Y(sp[4][0],sp[4][1]);
         if((g[0]^g[1])==(fb(yw^g[2],yc^g[3],ye^g[4])^fb(yw,yc,ye))) sat=1; }
       if(!sat) forced=1; }
   }
   if(forced){ stamp[w]=pass; added++; }
  }
  printf("pass %d added %ld remaining %ld\n",pass,added,remaining-added); fflush(stdout);
  if(remaining-added==0){ printf("COMPLETE after %d passes\n",pass); break; }
  if(added==0){ printf("STALLED\n"); break; }
 }
}
