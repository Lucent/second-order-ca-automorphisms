// (Sh) check by direct simulation. args: rule c0 c1 d0 d1 nmin nmax LP mode(0: kink images supp0, 1: kink images pure c)
#include "sc.h"
static int run(const int*c,const int*d,const unsigned char*p,int lp,int sP,const unsigned char*q,int lq,int sQ){
  int T=2*(lp+lq)+40; for(;;){ int W=T+lp+lq+40; if(W>=WMAX/2-2){fprintf(stderr,"W too big\n");exit(2);} int r=scatter(c,d,p,lp,sP,q,lq,sQ,T,W); if(r>=0) return r; T*=2; } }
typedef struct{unsigned char w[64]; int l,s;} probe_t;
int main(int argc,char**argv){ RULE=atoi(argv[1]); int c[2]={atoi(argv[2]),atoi(argv[3])}, d[2]={atoi(argv[4]),atoi(argv[5])};
 int nmin=atoi(argv[6]),nmax=atoi(argv[7]),LP=atoi(argv[8]),mode=atoi(argv[9]);
 int kink=(c[0]!=d[0]||c[1]!=d[1]);
 int cap=1<<22; probe_t*P=malloc(sizeof(probe_t)*cap); int np=0; int ntot=0;
 int sQ=2;
 for(int L=1;L<=LP;L++) for(long m=0;m<(1L<<L);m++) for(int e=0;e<2;e++){
   int sP=-(L+3)-e; unsigned char w[64]; for(int k=0;k<L;k++) w[k]=(m>>k)&1;
   if(w[0]==c[par(sP)]||w[L-1]==c[par(sP+L-1)]) continue;
   unsigned char s[80]; for(int k=0;k<L+4;k++){int j=sP-2+k; s[k]=(j>=sP&&j<sP+L)?w[j-sP]:c[par(j)];} if(!okseq(s,L+4)) continue;
   ntot++;
   if(kink){ int bad=0; for(int a=sQ;a<=sQ+1&&!bad;a++){ run(c,d,w,L,sP,w,0,a);
       if(mode==0){ if(suppof(GO,c,MID)) bad=1; } else { for(int j=-OFS;j<OFS;j++) if(GO[j+OFS]!=c[par(j)]){bad=1;break;} } }
     if(bad) continue; }
   memcpy(P[np].w,w,L); P[np].l=L; P[np].s=sP; np++; }
 printf("rule %d tails %d%d|%d%d mode %d: probes %d of %d admissible\n",RULE,c[0],c[1],d[0],d[1],mode,np,ntot); fflush(stdout);
 int*ord=malloc(sizeof(int)*np); for(int i=0;i<np;i++)ord[i]=i;
 for(int n=nmin;n<=nmax;n++){ long cnt=0,fail=0; int maxneed=0;
  for(long m=0;m<(1L<<n);m++){ unsigned char q[64]; for(int k=0;k<n;k++) q[k]=(m>>k)&1;
   if(q[0]==c[par(sQ)]||q[n-1]==d[par(sQ+n-1)]) continue;
   unsigned char s[80]; for(int k=0;k<n+4;k++){int j=sQ-2+k; s[k]=(j<sQ)?c[par(j)]:(j<sQ+n?q[j-sQ]:d[par(j)]);} if(!okseq(s,n+4)) continue;
   cnt++; int found=-1;
   for(int oi=0;oi<np;oi++){ int b=ord[oi]; run(c,d,P[b].w,P[b].l,P[b].s,q,n,sQ); int L[2]={c[0],c[1]}; if(suppof(GO,L,MID)<n){found=oi;break;} }
   if(found<0){ fail++; printf("FAIL n=%d q=",n); for(int k=0;k<n;k++)printf("%d",q[k]); printf("\n"); fflush(stdout); continue; }
   int b=ord[found]; if(P[b].l>maxneed) maxneed=P[b].l; for(int k=found;k>0;k--) ord[k]=ord[k-1]; ord[0]=b; }
  printf("n=%d objects %ld fail %ld maxneed %d\n",n,cnt,fail,maxneed); fflush(stdout); }
}
