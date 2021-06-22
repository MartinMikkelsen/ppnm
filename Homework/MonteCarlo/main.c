#include<complex.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#define SQR2 1.41421356237309504880


void randomx(int dim, double* a, double* b, double* x);

complex plainmc(int dim,double f(int dim,double*x),double*a,double*b,int N);
#define R 0.9

/*
double f(int dim, double*p )
{   assert(dim==2);
    double x = p[0];
    double y = p[1];
    if(x*x+y*y<R*R)return 1;
    else return 0;
};
*/
double Fa1(int dim, double*p)
{
    assert(dim==2);
    double x = p[0];
    if(sqrt(x)<R*R)return 1;
    else return 0;
}

int main(int argc, char** argv)
{   
    double a[]={0,0},b[]={1,1};
    int dim=sizeof(a)/sizeof(a[0]);
    if(argc>1){
	srand(42);
	int N=(int)atof(argv[1]);	
	complex result_p=plainmc(dim,Fa1,a,b,N);
	double integ_p=creal(result_p);
	double exact=M_PI/4*R*R;
	double error_p=fabs(integ_p-exact);
	printf("%i %g\n",N,error_p);
	}
else {
	int N=(int)1e7;
	complex result_p=plainmc(dim,Fa1,a,b,N);
	double integ_p=creal(result_p), esterr_p=cimag(result_p);
	double exact=2/3.*R*R;
	double error_p=fabs(integ_p-exact);
	printf("integrating sqrt(x)<%g*%g?1:0 with N=%i\n",R,R,N);
	printf("pseudo-random:\n");
	printf("integral = %f error estimate = %f\n",integ_p,esterr_p);
	printf("exact    = %f actual error   = %f\n",exact,error_p);
	}
    
return 0;
}
