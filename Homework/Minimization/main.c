#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include<stdio.h>

void numeric_gradient(double F(gsl_vector*), gsl_vector*x, gsl_vector*grad);

void Broyden(int dim, double F(int, double*), double* x, double acc);

int qnewton(double F(gsl_vector* x), gsl_vector*x, double acc);

void vector_print(const char* s,gsl_vector* v){
	printf("%s",s);
	for(int i=0;i<v->size;i++)printf("%9.3g ",gsl_vector_get(v,i));
	printf("\n");
	}

double Rosenbrock(gsl_vector*v){
	double x=gsl_vector_get(v,0);
	double y=gsl_vector_get(v,1);
	return pow(1-x,2)+100*pow(y-x*x,2);
}

double Himmelblau(gsl_vector*v){
	double x=gsl_vector_get(v,0);
	double y=gsl_vector_get(v,1);
	return pow(x*x+y-11,2)+pow(x+y*y-7,2);
}

int main(){

    int dim=2;
    int nsteps;
	gsl_vector* x=gsl_vector_alloc(dim);
	gsl_vector* g=gsl_vector_alloc(dim);
	double acc=1e-3;

	printf("\n===============================================\n");
	printf("\nMinimization of Rosenbrock's valley function\n");
	printf("\n          f(x,y)=(1-x)²+100(y-x²)²\n");
	printf("\n-----------------------------------------------\n");
	gsl_vector_set(x,0,M_PI);
	gsl_vector_set(x,1,6);
	vector_print("Starting point, x0    :",x);
	printf      ("Initial value, f(x0)  : %g\n",Rosenbrock(x));
	nsteps=qnewton(Rosenbrock,x,acc);
	vector_print("Minimum found, xmin   :",x);
	printf      ("Done in %i steps \n", nsteps);            
	printf      ("Value at min, f(xmin) : %g\n",Rosenbrock(x));
	numeric_gradient(Rosenbrock,x,g);
	printf      ("Accuracy goal         : %.1e\n",acc);
	vector_print("Gradient(min)         :",g);

	printf("\n===============================================\n");
	printf("\nMinimization of Himmelblau's function\n");
	printf("\n          f(x,y)=(x²+y-11)²+(x+y²-7)²\n");
	printf("\n-----------------------------------------------\n");
	gsl_vector_set(x,0,M_PI*M_PI);
	gsl_vector_set(x,1,12);
	vector_print("Starting point, x0    :",x);
	printf      ("Initial value, f(x0)  : %g\n",Himmelblau(x));
	nsteps=qnewton(Himmelblau,x,acc);
	vector_print("Minimum found, xmin   :",x);
	printf      ("Done in %i steps \n", nsteps);            
	printf      ("Value at min, f(xmin) : %g\n",Himmelblau(x));
	numeric_gradient(Himmelblau,x,g);
	vector_print("Gradient(min)         :",g);
	printf("\n-----------------------------------------------\n");


return 0;
}
