#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#define FOR(k) for (int k = 0; k<dim;k++)
#define RND (double)rand()/(RAND_MAX)

double strata(
	int dim,
	double f(int dim,double*x),
	double*a,double*b,
	double acc,double eps,
	int n_reuse,double mean_reuse)
{
int N=16*dim;
double V=1; FOR(k)V*=b[k]-a[k];
int n_l[dim],n_r[dim];
double x[dim],mean_l[dim],mean_r[dim],mean=0;
FOR(k){ mean_l[k]=0; mean_r[k]=0; n_l[k]=0; n_r[k]=0; }
for(int i=0;i<N;i++){
	FOR(k) x[k]=a[k]+RND*(b[k]-a[k]);
	double fx=f(dim,x);
	FOR(k){
		if(x[k]>(a[k]+b[k])/2){ n_r[k]++; mean_r[k]+=fx; }
		else                  { n_l[k]++; mean_l[k]+=fx; }
		}
	mean+=fx;
	}
mean/=N;
FOR(k){ mean_l[k]/=n_l[k]; mean_r[k]/=n_r[k]; }

int kdiv=0; double maxvar=0;
FOR(k){
	double var=fabs(mean_r[k]-mean_l[k]);
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
double integ_l=strata(dim,f,a,b2,acc/sqrt(2.),eps,n_l[kdiv],mean_l[kdiv]);
double integ_r=strata(dim,f,a2,b,acc/sqrt(2.),eps,n_r[kdiv],mean_r[kdiv]);
return integ_l+integ_r;
}
