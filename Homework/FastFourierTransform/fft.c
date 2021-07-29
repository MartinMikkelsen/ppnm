#include<complex.h>
#include<tgmath.h>

void dfts(int sign, int N, complex* x, complex* c, int stride){
	for(int k=0;k<N;k++){
		c[k]=0;
		for(int n=0;n<N;n++) c[k]+=x[n*stride]*exp(sign*2*M_PI*I*n*k/N);
	}
}

void fft2s(int sign,int N,complex*x,complex*c,int stride){
if(N==1) c[0]=x[0];
else if(N%2==0){
	fft2s(sign,N/2,x       ,c    ,2*stride);
	fft2s(sign,N/2,x+stride,c+N/2,2*stride);
	for(int k=0;k<N/2;k++){
		complex p=c[k];
		complex q=exp(sign*2*M_PI*I*k/N)*c[k+N/2];
		c[k    ]=p+q;
		c[k+N/2]=p-q;
		}
	}
else dfts(sign,N,x,c,stride);
}

void fft(int N,complex*x,complex*c){
	fft2s(-1,N,x,c,1);
	}
void ift(int N,complex*c,complex*x){
	fft2s(+1,N,c,x,1);
	for(int i=0;i<N;i++)x[i]/=N;
	}
