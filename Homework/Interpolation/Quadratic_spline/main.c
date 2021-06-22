#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_interp.h>
#include <assert.h>

typedef struct {
   int n;      // number of data points
   double* x;  // arrays
   double* y;  // .
   double* b;  // .
   double* c;  // .
} qspline;

qspline* qspline_alloc(int n, double* x, double* y);
double qspline_integ(qspline *s, double input);
double qspline_eval(qspline *s, double input);
double qspline_deriv(qspline* s, double x_new);
void qspline_free(qspline *s);

int main() 
{
    int n=20;
	double x[n],y[n];

	for(int i=0;i<n;i++)
    {
		x[i]=i;
		y[i]=sin(i*i);
	}

	printf("# index 0: Data points\n");
	for(int i=0;i<n;i++) printf("%g %g\n",x[i],y[i]);
	printf("\n\n");

    printf("# index 1: Interpolations\n");
    qspline* s = qspline_alloc(n, x, y);
        
	double dz=0.1;
	for(double z=x[0];z<=x[n-1];z+=dz){
		double qz=qspline_eval(s,z);
		printf("%g %g\n",z,qz);       
    }     
	printf("\n\n");
    printf("# index 2: Integral\n");
	for(double z=x[0];z<=x[n-1];z+=dz){
		double integ=qspline_integ(s,z);
		printf("%g %g\n",z,integ);
	}
	printf("\n\n");
    printf("# index 2: Derivative\n"); 
	for(double z=x[0];z<=x[n-1];z+=dz){
		double deriv=qspline_deriv(s,z);
		printf("%g %g\n",z,deriv);
    }

    qspline_free(s);
return 0;
}
