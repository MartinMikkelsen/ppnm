#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include<stdio.h>


void numeric_gradient(double F(gsl_vector*), gsl_vector*x, gsl_vector*grad);
void Broyden(int dim, double F(int, double*), double* x, double acc);
int qnewton(double F(gsl_vector*), gsl_vector*x, double acc);
void vector_print(const char* s,gsl_vector* v){
	printf("%s",s);
	for(int i=0;i<v->size;i++)printf("%9.3g ",gsl_vector_get(v,i));
	printf("\n");
	}
double rosen(gsl_vector*v){
	double x=gsl_vector_get(v,0);
	double y=gsl_vector_get(v,1);
	return pow(1-x,2)+100*pow(y-x*x,2);
}

double himmel(gsl_vector*v){
	double x=gsl_vector_get(v,0);
	double y=gsl_vector_get(v,1);
	return pow(x*x+y-11,2)+pow(x+y*y-7,2);
}

int main(){
    int n=2,nsteps;
	gsl_vector* x=gsl_vector_alloc(n);
	gsl_vector* g=gsl_vector_alloc(n);
	double acc=1e-3;

	printf("\nMINIMIZATION OF ROSENBROCK'S FUNCTION\n");
	gsl_vector_set(x,0,18);
	gsl_vector_set(x,1,15);
	vector_print("start point   :",x);
	printf      ("value(start)  : %g\n",rosen(x));
	nsteps=qnewton(rosen,x,acc);
	vector_print("minimum found :",x);
	printf      ("steps         : %i\n",nsteps);
	printf      ("value(min)    : %g\n",rosen(x));
	numeric_gradient(rosen,x,g);
	printf      ("acc           : %.1e\n",acc);
	vector_print("gradient(min) :",g);

    printf("\nMINIMIZATION OF HIMMELBLAU'S FUNCTION\n");
	gsl_vector_set(x,0,18);
	gsl_vector_set(x,1,15);
	vector_print("start point   :",x);
	printf      ("value(start)  : %g\n",himmel(x));
	nsteps=qnewton(himmel,x,acc);
	printf      ("steps         : %i\n",nsteps);
	vector_print("minimum       :",x);
	printf      ("value(min)    : %g\n",himmel(x));
	numeric_gradient(himmel,x,g);
	printf      ("acc           : %.1e\n",acc);
	vector_print("gradient(min) :",g);


return 0;
}

