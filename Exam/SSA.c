#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#define FOR(k) for (int k = 0; k<dim;k++)
#define RND (double)rand()/(RAND_MAX)

double strata(int dim,double f(int dim,double*x),double*a,double*b,double acc,double eps,int n_reuse,double mean_reuse)
{
int N=16*dim;
double V=1; FOR(k)V*=b[k]-a[k];
int n_left[dim],n_right[dim];
double x[dim],mean_left[dim],mean_right[dim],mean=0;
FOR(k){ mean_left[k]=0; mean_right[k]=0; n_left[k]=0; n_right[k]=0; }
for(int i=0;i<N;i++){
	FOR(k) x[k]=a[k]+RND*(b[k]-a[k]);
	double fx=f(dim,x);
	FOR(k){
		if(x[k]>(a[k]+b[k])/2){ n_right[k]++; mean_right[k]+=fx; }
		else                  { n_left[k]++; mean_left[k]+=fx; }
		}
	mean+=fx;
	}
mean/=N;
FOR(k){ mean_left[k]/=n_left[k]; mean_right[k]/=n_right[k]; }

int kdiv=0; double maxvar=0;
FOR(k){
	double var=fabs(mean_right[k]-mean_left[k]);
	if(var>maxvar){ maxvar=var; kdiv=k; }
	}

double integ=(mean*N+mean_reuse*n_reuse)/(N+n_reuse)*V;
double error=fabs(mean_reuse-mean)*V;
double toler=acc+fabs(integ)*eps;
if(error<toler)return integ;

double a2[dim],b2[dim];
FOR(k)a2[k]=a[k];
FOR(k)b2[k]=b[k];
a2[kdiv]=(a[kdiv]+b[kdiv])/2;
b2[kdiv]=(a[kdiv]+b[kdiv])/2;
double integ_left=strata(dim,f,a,b2,acc/sqrt(2.),eps,n_left[kdiv],mean_left[kdiv]);
double integ_right=strata(dim,f,a2,b,acc/sqrt(2.),eps,n_right[kdiv],mean_right[kdiv]);
return integ_left+integ_right;
}
