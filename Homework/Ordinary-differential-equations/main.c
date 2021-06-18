 
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
	void (*f)(int n, double x, double* y, double* dydx), // right-hand-side of dy/dx = f(x, y) 
   int n,          // size of vectors 
	double  a,      // the start-point a 
	double  b,      // the end-point of the integration 
	double* ya,     // y(a) 
	double* yb,     // y(b) to be calculated 
	double  h,      // initial step-size 
	double  acc,    // absolute accuracy goal 
	double  eps,    // relative accuracy goal
    char*   outfile // trajectory file
    );


void f(int n,double x,double*y,double*dydt){
	dydt[0]=+y[1];
	dydt[1]=2*(x*x/2-0.5)*y[0];
	}

double N;
double Tr;
double Tc;

void SIR_model(int n, double x, double* y, double*dydx) {
   dydx[0] = -y[0]*y[1]/(N*Tc);
   dydx[1] =  y[0]*y[1]/(N*Tc) - y[1]/Tr;
   dydx[2] =  y[1]/Tr;
}



int main(){

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

   outfile = "sir1.txt";
   N = 5.8e6;
   Tr = 7;
   Tc = 1.0;

   ya[0] = N;
   ya[1] = 50;
   ya[2] = 0;
   
   driver(&f_SIR, n, a, b, ya, yb, h, acc, eps, outfile);




return 0;    
}
