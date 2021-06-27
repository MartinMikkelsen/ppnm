#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#define SQR2 1.41421356237309504880
#define R 0.8


void randomx(int dim, double* a, double* b, double* x);

void plainmc(int dim,double f(int dim,double* x),double* a,double* b,int N, double* result, double* error);

void Halton(int n, int dim, double *a, double *b, double *x);

void pMonteCarlo(int dim, double f(int dim, double* x), double* a, double* b, double N, double* result2, double* error2);

double strata(int dim, double f(int dim, double* x),double* a, double* b, double acc, double eps, int n_reuse,double mean_reuse);

double intrun(double f(double), double a, double b, double acc, double eps, double f2, double f3, double nrec);

double integrate(double f(double), double a, double b,double acc,double eps);

double clenshaw(double f(double),double a,double b, double acc,double eps);

//These functions below are from the quadratures homework. 

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
double exact = M_PI;
double a[]={0,0},b[]={1,1},acc=1e-3,eps=1e-3;
    {
    double result;
    double error;
    double result2;
    double error2;
    double N=1e6;
	int dim=sizeof(a)/sizeof(a[0]);
	double integ=strata(dim,Fa2,a,b,acc,eps,1e6,1e6);
    plainmc(dim, Fa2, a, b, N, &result, &error);
    pMonteCarlo(dim, Fa2, a, b, N, &result2, &error2);
    
    printf("------------------------------------------------------------\n");
	printf("∫ dx 4√(1-x²) from 0 to 1 using different methods\n");
	printf("The exact value of the integral is π\n");
    printf("Absolute accuracy goal    : %.1e\n"   , acc);
    printf("Relative accuracy goal    : %.1e\n"   , eps);
    printf("------------------------------------------------------------\n");
	printf("Using adaptive 1D integrator with random nodes\n");
	printf("Value       = %.10f\n",integ);
    printf("Error est.  = %.10f\n", acc+fabs(integ)*eps);
    printf("Diff        = %.10f\n",fabs(exact-integ));
    printf("Using n_reuse       = %g\n",1e6);
    printf("Using mean_reuse    = %g\n",1e6);
    printf("------------------------------------------------------------\n");
    printf("Using plain Monte Carlo \n");
    printf("Value       = %.10f \n", result);
    printf("Error       = %.10f \n", error);
    printf("Diff        = %.10f \n", fabs(exact-result));
    printf("Using  N    = %g\n",N);
    printf("------------------------------------------------------------\n");
    printf("Using pseudo Monte Carlo \n");
    printf("Value       = %.10f \n", result2);
    printf("Error       = %.10f \n", error2);
    printf("Diff        = %.10f \n", exact-result2);
    printf("Using  N    = %g\n",N);
    }
    {
    double a = 0.;
    double b = 1.;
    double I1 = integrate(Fa2,a,b,acc,eps);
    double diff = exact-I1;
    printf("------------------------------------------------------------\n");
    printf("Using recursive adaptive integrator \n");
    printf("Value       = %.10f\n"  , I1);
    printf("Diff        = %.10f\n"  , fabs(diff));
    printf("------------------------------------------------------------\n");
    }
    {
    double a = 0.;
    double b = 1.;
    double I3 = clenshaw(Fa2,a,b,acc,eps);
    double diff = exact-I3;
    printf("Using Clenshaw–Curtis transformation\n");
    printf("Value       = %.10f\n"  , I3);
    printf("Diff        = %.10f\n"  , fabs(-diff));
    printf("------------------------------------------------------------\n");
    }
    {
    double res, err;

    double xl[1] = { 0 };
    double xu[1] = { 1 };

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G = { &Fa2, 1, 0 };

    size_t calls = 500000;

    gsl_rng_env_setup ();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    gsl_monte_miser_state *s = gsl_monte_miser_alloc(1);
    gsl_monte_miser_integrate (&G, xl, xu, 1, calls, r, s,
                               &res, &err);
    gsl_monte_miser_free (s);    
    printf("Using GSL MISER\n");
    printf("Value       = %.10f\n"  , res);
    printf("Error       = %.10f\n"  , err);
    }

return 0;
}
