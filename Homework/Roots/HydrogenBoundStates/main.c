#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<stdio.h>
#include<float.h>
#include<math.h>

void Jacobian(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x);

int Newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps);

void vector_print(gsl_vector* vec);

static int ncalls;

void f(gsl_vector* p,gsl_vector* fx){
	ncalls++;
	double x=gsl_vector_get(p,0), y=gsl_vector_get(p,1);
	gsl_vector_set(fx,0, 2*(1-x)*(-1)+100*2*(y-x*x)*(-1)*2*x);
	gsl_vector_set(fx,1, 100*2*(y-x*x));
	}

void Rosenbrock(gsl_vector* x, gsl_vector* fx){
    double X = gsl_vector_get(x,0);
    double Y = gsl_vector_get(x,1);
    gsl_vector_set(fx,0,2*(1-X)-4*X*100*(Y-X*X));
    gsl_vector_set(fx,1,2*100*(Y-X*X));
}


int main()
{
    gsl_vector* x = gsl_vector_alloc(2);
    gsl_vector_set(x,0,1.3);
    gsl_vector_set(x,1,1.3);
    double eps = 1e-5;
    Newton(Rosenbrock,x,eps);
    printf("The extremum of the Rosenbrock valley function is located at (x= %g,y=%g) \n", gsl_vector_get(x,0),gsl_vector_get(x,1));
return 0;    
}

