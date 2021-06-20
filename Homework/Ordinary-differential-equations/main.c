
#include <assert.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <time.h>
#include <math.h>

#include "utilities.h"
#include "Runge-Kutta.h"
#include <stdio.h>


void rkstep12(void f(double x,gsl_vector* y, gsl_vector* dydx), double x, gsl_vector* yx, double h, gsl_vector* yh, gsl_vector* err);

void rkstep23(void (*f)(int n, double x, double* y, double* dydx), int n, double x, double* y_curr, double h, double* y_next, double* err);

int driver(
	void (*f)(int n, double x, double* y, double* dydx), 
    int n,          // size of vectors
	double  a,      // the start-point a
	double  b,      // the end-point of the integration
	double* ya,     // y(a)
	double* yb,     // y(b) to be calculated
	double  h,      // initial step-size
	double  acc,    // absolute accuracy goal
	double  eps,    // relative accuracy goal
    char*   outfile // target file
    );


void f(int n,double x,double*y,double*dydt){
	dydt[0]=+y[1];
	dydt[1]=2*((x)*(x)/2-0.5)*y[0];
	}

double N;
double Tr;
double Tc;

void SIR_model(int n, double x, double* y, double*dydx) {
   dydx[0] = -y[0]*y[1]/(N*Tc);
   dydx[1] =  y[0]*y[1]/(N*Tc) - y[1]/Tr;
   dydx[2] =  y[1]/Tr;
}

int main() {
{
    int n = 2;
    double a = 0;
    double b = 10;
    double ya[n];
    double yb[n];
    double h = 0.001;
    double acc = 1e-4;
    double eps = 1e-4;
    char* outfile = "solu1.txt";
    ya[0] = 1;

    driver(&f, n, a, b, ya, yb, h, acc, eps, outfile);
}


{
    int n = 3;
    int a = 0;
    int b = 60;
    double ya[n];
    double yb[n];
    double h = 0.05;
    double acc = 1;
    double eps = 1e-6;
    char* outfile;

    outfile = "sir_model_1.txt";
    N = 5.806e6;
    Tr = 2;
    Tc = 1.0;

    ya[0] = N;
    ya[1] = 15;
    ya[2] = 0;

    driver(&SIR_model, n, a, b, ya, yb, h, acc, eps, outfile);
}
{
    int n = 3;
    int a = 0;
    int b = 60;
    double ya[n];
    double yb[n];
    double h= 0.05;
    double acc = 1;
    double eps = 1e-6;
    char* outfile;

    outfile = "sir_model_2.txt";
    N = 5.806e6;
    Tr = 6;
    Tc = 2;

    ya[0] = N;
    ya[1] = 15;
    ya[2] = 0;

    driver(&SIR_model, n, a, b, ya, yb, h, acc, eps, outfile);
}
{
  int n = 3;
  int a = 0;
  int b = 60;
  double ya[n];
  double yb[n];
  double h= 0.05;
  double acc = 1;
  double eps = 1e-6;
  char* outfile;

  outfile = "sir_model_3.txt";
  N = 5.806e6;
  Tr = 8;
  Tc = 2;

  ya[0] = N;
  ya[1] = 15;
  ya[2] = 0;

  driver(&SIR_model, n, a, b, ya, yb, h, acc, eps, outfile);
}
   return 0;
}
