
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

void threebody(int n, double x, double*y, double* dydx){

    //Inital values from the article: https://arxiv.org/abs/math/0011268

    dydx[0] = y[6];
    dydx[1] = y[7];
    dydx[2] = y[8];
    dydx[3] = y[9];
    dydx[4] = y[10];
    dydx[5] = y[11];

    double r12 = pow(y[2] - y[0], 2) + pow(y[3] - y[1], 2);
    double r13 = pow(y[4] - y[0], 2) + pow(y[5] - y[1], 2);
    double r23 = pow(y[4] - y[2], 2) + pow(y[5] - y[3], 2);

    dydx[6]  = (y[2] - y[0])*pow(r12, -3.0/2) + (y[4] - y[0])*pow(r13, -3.0/2);
    dydx[7]  = (y[3] - y[1])*pow(r12, -3.0/2) + (y[5] - y[1])*pow(r13, -3.0/2);
    dydx[8]  = (y[0] - y[2])*pow(r12, -3.0/2) + (y[4] - y[2])*pow(r23, -3.0/2);
    dydx[9]  = (y[1] - y[3])*pow(r12, -3.0/2) + (y[5] - y[3])*pow(r23, -3.0/2);
    dydx[10] = (y[0] - y[4])*pow(r13, -3.0/2) + (y[2] - y[4])*pow(r23, -3.0/2);
    dydx[11] = (y[1] - y[5])*pow(r13, -3.0/2) + (y[3] - y[5])*pow(r23, -3.0/2);
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
{
    int n = 12;
    double a = 0;
    double b = 6; //time
    double ya[n];
    double yb[n];
    double h = 0.001;
    double acc = 1e-6;
    double eps = 1e-6;
    char* outfile = "threebody.txt";

    //Initial values from https://arxiv.org/abs/math/0011268
    ya[0]  =  0.97000436;
    ya[1]  = -0.24308753;
    ya[2]  = -0.97000436;
    ya[3]  =  0.24308753;
    ya[4]  =  0;
    ya[5]  =  0;

    ya[6]  =  0.93240737/2;
    ya[7]  =  0.86473146/2;
    ya[8]  =  0.93240737/2;
    ya[9]  =  0.86473146/2;
    ya[10] = -0.93240737;
    ya[11] = -0.86473146;

   driver(&threebody, n, a, b, ya, yb, h, acc, eps, outfile);

}
return 0;
}
