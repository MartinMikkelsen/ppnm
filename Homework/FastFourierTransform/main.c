#include<complex.h>
#include<tgmath.h>
#include<stdlib.h>
#include<stdio.h>
#define RND (double)rand()/RAND_MAX
#define FOR(i) for(int i=0;i<N;i++)
void fft(int N,complex*x,complex*c);
void ift(int N,complex*c,complex*x);
int main(int argc,char** argv){
if(argc>1){
	int N=atoi(argv[1]);
	printf("N=%i\n",N);
	complex xs[N],cs[N];
	FOR(i) xs[i]=I;
	fft(N,xs,cs);
	}
else	{
	//int N=128;
	int N=3*64;
	double noise=2,Afilter=0.5,Ffilter=0.02;
	complex ts[N],xs[N];
	FOR(i){
		ts[i]=(double)i/N;
		double omega=2*M_PI;
		xs[i]=sin(omega*ts[i])+cos(2*omega*ts[i]-2);
		}
	FOR(i) printf("%g %g\n",creal(ts[i]),creal(xs[i]));
	printf("\n\n");
	FOR(i) xs[i] += noise*RND*(RND-0.5)*2;
	FOR(i) printf("%g %g\n",creal(ts[i]),creal(xs[i]));
	printf("\n\n");
	complex cs[N];
	fft(N,xs,cs);
	double cmax=0;
	FOR(i) cmax=fabs(cs[i])>cmax?fabs(cs[i]):cmax;
	complex cA[N];
	FOR(i) if(fabs(cs[i]) < cmax*Afilter) cA[i]=0; else cA[i]=cs[i];
	complex ys[N];
	ift(N,cA,ys);
	FOR(i) printf("%g %g\n",creal(ts[i]),creal(ys[i]));
	printf("\n\n");
	complex cF[N];
	FOR(i) if(i>N*Ffilter && i<N-N*Ffilter) cF[i]=0; else cF[i]=cs[i];
	ift(N,cF,ys);
	FOR(i) printf("%g %g\n",creal(ts[i]),creal(ys[i]));
	}
}
