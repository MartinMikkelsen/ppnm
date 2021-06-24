#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include"ann.h"

double activationfunc2(double x){return exp(-x*x);}
double activationfunc3(double x){return cos(5*x)*exp(-x*x);}
double function_to_fit(double x){return cos(5*x-1)*exp(-x*x);}

double activationfunc1(double x);
int main()
{
    int n=5;
	ann* network=ann_alloc(n,activationfunc1);
	double a=-1,b=1;
	int nx=20;
	gsl_vector* vx=gsl_vector_alloc(nx);
	gsl_vector* vy=gsl_vector_alloc(nx);

	for(int i=0;i<nx;i++){
		double x=a+(b-a)*i/(nx-1);
		double f=function_to_fit(x);
		gsl_vector_set(vx,i,x);
		gsl_vector_set(vy,i,f);
	}

	for(int i=0;i<network->n;i++){
		gsl_vector_set(network->params,3*i+0,a+(b-a)*i/(network->n-1));
		gsl_vector_set(network->params,3*i+1,1);
		gsl_vector_set(network->params,3*i+2,1);
	}

	ann_train(network,vx,vy);

	for(int i=0;i<vx->size;i++){
		double x=gsl_vector_get(vx,i);
		double f=gsl_vector_get(vy,i);
		printf("%g %g\n",x,f);
	}
	printf("\n\n");

	double dz=1.0/64;
	for(double z=a;z<=b;z+=dz){
	double y  = ann_response(network,z);
    double yprime  = ann_feedDeriv(network,z);
    double yprimeprime  = ann_feedInt(network,z);
		printf(" %.12e %.12e %.12e %.12e\n",z,y,yprime, yprimeprime);
    }

	for(int i=0;i<network->n;i++){
		double ai=gsl_vector_get(network->params,3*i+0);
		double bi=gsl_vector_get(network->params,3*i+1);
		double wi=gsl_vector_get(network->params,3*i+2);
		fprintf(stderr,"i=%i ai,bi,wi = %g %g %g\n",i,ai,bi,wi);
	}




gsl_vector_free(vx);
gsl_vector_free(vy);
ann_free(network);

return 0;
}
