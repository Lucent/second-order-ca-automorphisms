/* Example 7.11: exact checks on the binary d-shift for an explicit S_3 marker construction.
   Frames: 111 u 111 with u in U = {0010110, 0101010, 0110100} (all of parity key (even-sum,odd-sum)=(0,1)),
   and the complementary frames 000 ubar 000.  rho in S_3 acts by u_i -> u_{rho(i)} (same on complements).
   psi_rho(d): replace every frame segment u_i by u_{rho(i)}.
   J_rho(d)_n = XOR_{k'>=k, k'=k mod 2} delta[k'] if n = p+3+k lies in the segment of a frame at p
                (delta = u_i xor u_{rho(i)}), else 0.
   B0: distinct frames never share segment cells (all pairs of offsets 1..12).
   B1: d + J + sigma^2 J = psi(d)  on all 2^21 windows d[n-9..n+11]            (exact)
   B2: J(not d) = J(d)             on all 2^19 windows d[n-9..n+9]             (exact)
   B3: frame set preserved: frame_p(psi_rho d) == frame_p(d), with the segment of psi_rho d equal to rho(segment),
       on ALL 2^31 windows d[p-9..p+21]                                         (exact)  */
#include <stdio.h>
#include <stdint.h>
static const char *Us[3]={"0010110","0101010","0110100"};
static int U[3]; static int16_t frameid[8192];  /* 13-bit pattern (bit 12 = first cell) -> 1+i (111 frames), -(1+i) (000 frames), 0 */
static int perms[6][3]={{0,1,2},{1,0,2},{0,2,1},{2,1,0},{1,2,0},{2,0,1}};
static inline int bitat(uint64_t w,int len,int pos){ return (w>>(len-1-pos))&1; }   /* pos 0 = leftmost */
static inline int frame_at(uint64_t w,int len,int p){ if(p<0||p+13>len) return 99; return frameid[(w>>(len-13-p))&0x1FFF]; }
/* J and psi at position n of window w (needs frames at p=n-3-k, k=0..6) */
static inline int Jat(uint64_t w,int len,int n,const int*rho){
  for(int k=0;k<7;k++){ int p=n-3-k; int f=frame_at(w,len,p); if(f==99) return -1; if(f){ int i=(f>0?f:-f)-1;
      int delta=U[i]^U[rho[i]]; int J=0; for(int kk=k;kk<7;kk+=2) J^=(delta>>(6-kk))&1; return J; } }
  return 0;
}
static inline int psiat(uint64_t w,int len,int n,const int*rho){
  for(int k=0;k<7;k++){ int p=n-3-k; int f=frame_at(w,len,p); if(f==99) return -1; if(f){ int i=(f>0?f:-f)-1;
      int v=U[rho[i]]; int b=(v>>(6-k))&1; return f>0?b:1-b; } }
  return bitat(w,len,n);
}
int main(){
  for(int i=0;i<3;i++){U[i]=0; for(int k=0;k<7;k++) U[i]=(U[i]<<1)|(Us[i][k]-'0');}
  for(int i=0;i<3;i++){ int pat=(7<<10)|(U[i]<<3)|7; frameid[pat]=1+i; int comp=(~pat)&0x1FFF; frameid[comp]=-(1+i); }
  /* B0 */
  long long overl=0;
  for(int d=1;d<=12;d++) for(uint64_t w=0;w<(1ULL<<(13+d));w++){ int a=frame_at(w,13+d,0),b=frame_at(w,13+d,d);
      if(a&&b){ if(d!=10) overl++; } }
  printf("B0 pairs of frames at offsets 1..12 other than 10 (shared marker): %lld\n",overl);
  for(int r=1;r<6;r++){ const int*rho=perms[r];
    long long b1=0,b2=0,nt=0;
    for(uint64_t w=0;w<(1ULL<<21);w++){ int n=9; int J0=Jat(w,21,n,rho),J2=Jat(w,21,n+2,rho),P=psiat(w,21,n,rho);
      if((bitat(w,21,n)^J0^J2)!=P) b1++; if(J0) nt++; }
    for(uint64_t w=0;w<(1ULL<<19);w++){ uint64_t c=(~w)&((1ULL<<19)-1); if(Jat(w,19,9,rho)!=Jat(c,19,9,rho)) b2++; }
    long long b3=0, b3seg=0;
    #pragma omp parallel for reduction(+:b3,b3seg) schedule(static)
    for(long long e=0;e<(1LL<<31);e++){ uint64_t w=(uint64_t)e; const int len=31, p=9;
      int f0=frame_at(w,len,p);
      uint64_t img=0; for(int j=p;j<p+13;j++) img=(img<<1)|psiat(w,len,j,rho);
      int f1=frameid[img];
      if((f0!=0)!=(f1!=0) || (f0>0)!=(f1>0)) b3++;
      else if(f0){ int i=(f0>0?f0:-f0)-1, i1=(f1>0?f1:-f1)-1; if(i1!=rho[i]) b3seg++; }
    }
    printf("rho=(%d%d%d): B1 failures %lld (J!=0 at %lld windows) | B2 failures %lld | B3 frame-set failures %lld, segment-label failures %lld\n",
           rho[0],rho[1],rho[2],b1,nt,b2,b3,b3seg);
  }
  return 0;
}
