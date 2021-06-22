#include<math.h>
#include<stdlib.h>
#include<complex.h>
#define RND (double)rand()/RAND_MAX
//Implement plain Monte Carlo multi-dimensional integration. 

void randomx(int dim, double* a, double* b, double* x)
{
    for(int i=0;i<dim;i++)
    {
        x[i] = a[i]+RND*(b[i]-a[i]);
    }
}

complex plainmc(int dim,double f(int dim,double* x),double* a,double* b,int N){
	double V=1; for(int i=0;i<dim;i++)V*=b[i]-a[i];
	double sum=0,sum2=0,x[dim];
	for(int i=0;i<N;i++){
		for(int i=0;i<dim;i++)x[i]=a[i]+RND*(b[i]-a[i]);
		double fx=f(dim,x);
		sum+=fx;
		sum2+=fx*fx;
		}
	double mean=sum/N, sigma=sqrt(sum2/N-mean*mean);
	complex result=mean*V+I*sigma*V/sqrt(N);
	return result;
}

 

