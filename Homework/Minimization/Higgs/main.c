#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include<stdio.h>

//functions from other files

void numeric_gradient(double F(gsl_vector*), gsl_vector*x, gsl_vector*grad);

void Broyden(int dim, double F(int, double*), double* x, double acc);

int qnewton(double F(gsl_vector* x), gsl_vector*x, double acc);

void vector_print(const char* s,gsl_vector* v){
	printf("%s",s);
	for(int i=0;i<v->size;i++)printf("%9.3g ",gsl_vector_get(v,i));
	printf("\n");
	}

int n;

//define BreitWigner function
double BreitWigner(double E,gsl_vector*v){
	double m = gsl_vector_get(v,0);
	double Gamma = gsl_vector_get(v,1);
    double A = gsl_vector_get(v,2);
	return A/(pow((E-m),2)+pow(Gamma/2,2));
}
//define deviation function (as a function of BreitWigner). Also load data.
double deviation(gsl_vector* v) {
    double s = 0;
    double Energy[] =  {101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127, 129, 131, 133, 135, 137, 139, 141, 143, 145, 147, 149, 151, 153, 155, 157, 159};
	double Cross[] = {-0.25, -0.3, -0.15, -1.71, 0.81, 0.65, -0.91, 0.91, 0.96, -2.52, -1.01, 2.01, 4.83, 4.58, 1.26, 0.45, 0.15, -0.9,1 -0.81, -1.41, 1.36, 0.5, -0.45, 1.61, -2.21, -1.86, 1.76, -0.5};
	double Error[] = {2, 2, 1.9, 1.9, 1.9, 1.9, 1.9, 1.9, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 0.9, 0.9, 0.9};
    for (int i=0;i<n;i++){
	double breitwigner = BreitWigner(Energy[i],v);
	double cross = Cross[i];
	double error= Error[i];
	s+=(breitwigner-cross)*(breitwigner-cross)/(error*error);
	}
	return s;
}

int main(){

    n= 20;
    gsl_vector* x0     = gsl_vector_alloc(n);
    gsl_vector_set(x0, 0, 125); // m
    gsl_vector_set(x0, 1, 5); // gamma
    gsl_vector_set(x0, 2, 5); // A
    double eps = 1e-6;
    qnewton(deviation,x0,eps);
    printf("================ Higgs ===============\n");
    printf("Mass = %f GeV\nWidth = %f\nA = %f \n",fabs(gsl_vector_get(x0,0)),fabs(gsl_vector_get(x0,1)), gsl_vector_get(x0,2));
    printf("--------------------------------------\n");

    gsl_vector_free(x0);

return 0;
}
