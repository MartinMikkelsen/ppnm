#include<complex.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#define SQR2 1.41421356237309504880
#define R 0.8


void randomx(int dim, double* a, double* b, double* x);

void plainmc(int dim,double f(int dim,double* x),double* a,double* b,int N, double* result, double* error);

void Halton(int n, int dim, double *a, double *b, double *x);

void pMonteCarlo(int dim, double f(int dim, double* x), double* a, double* b, double N, double* result2, double* error2);

double strata(int dim, double f(int dim, double* x),double* a, double* b, double acc, double eps, int n_reuse,double mean_reuse);

//From quadratures

double Fa1(int dim,double x)
{
    return sqrt(x);
};
double Fa2(int dim,double x)
{
    return 4.*sqrt(1-(x)*(x));
};
double Fa3(double x)
{
    return 1/(sqrt(x));
}
double Fa4(double x)
{
    return log(x)/sqrt(x);
}

double func(int dim, double* x)
{
   return 1.0/(1 - cos(x[0])*cos(x[1])*cos(x[2]))/(M_PI*M_PI*M_PI);
}
int main()
{

	double a[]={0,0},b[]={1,1},acc=1e-3,eps=1e-3;

	int dim=sizeof(a)/sizeof(a[0]);
	double integ=strata(dim,Fa1,a,b,acc,eps,1e6,1e6);
	double exact=2/3.;
	printf("∫ dx √(x) from 0 to 1 using strata =%g\n",R);
	printf("a = {%g,%g} b = {%g,%g}\n",a[0],a[1],b[0],b[1]);
	printf("acc = %g eps = %g\n\n",acc,eps);
	printf("integ = %g error-estimate = %g\n",integ,acc+fabs(integ)*eps);
	printf("exact = %g actual-error   = %g\n",exact,fabs(integ-exact));

    
    
    return 0;
}
