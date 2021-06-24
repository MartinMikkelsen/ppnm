#include<gsl/gsl_vector.h>
#include<stdio.h>
#include<math.h>
#include"ann.h"

double activationfunc1(double x){return x*exp(-x*x);}
double f_deriv(double x){return exp(-x*x)*(1-2*pow(x,2));}
double f_int(double x){return (-1)*exp(-x*x)/2;}

int n, startPoint;

ann* ann_alloc(int n,double(*f)(double)){
	ann* network = malloc(sizeof(ann));
	network->n=n;
	network->f=f;
	network->params=gsl_vector_alloc(3*n);
	return network;
}
void ann_free(ann* network){
	gsl_vector_free(network->params);
	free(network);
}

double ann_response(ann* network,double x){
	double s=0;
	for(int i=0;i<network->n;i++){
		double a=gsl_vector_get(network->params,3*i+0);
		double b=gsl_vector_get(network->params,3*i+1);
		double w=gsl_vector_get(network->params,3*i+2);
		s+=network->f((x-a)/b)*w;
	}
	return s;
}

double ann_feedDeriv(ann* network,double x)
{
	double sum = 0;
	for(int i=0; i<network->n; i++){
		double a = gsl_vector_get(network->params,3*i+0);
		double b = gsl_vector_get(network->params,3*i+1);
		double w = gsl_vector_get(network->params,3*i+2);
		sum += f_deriv((x-a)/b)*w/b;
		}
	return sum;
}
double ann_feedInt(ann* network,double x)
{
	double sum = 0;
	for(int i=0; i<network->n; i++){
		double a = gsl_vector_get(network->params,3*i+0);
		double b = gsl_vector_get(network->params,3*i+1);
		double w = gsl_vector_get(network->params,3*i+2);
		sum += f_int((x-a)/b)*b*w;
		sum -= f_int((startPoint-a)/b)*b*w;
		}
	return sum;
}


int qnewton(double F(gsl_vector* x), gsl_vector*x, double acc);

void ann_train(ann* network,gsl_vector* xs,gsl_vector* ys){

	double cost_function(gsl_vector* p){
		gsl_vector_memcpy(network->params,p);
		double sum=0;
		for(int i=0;i<xs->size;i++){
			double xi=gsl_vector_get(xs,i);
			double yi=gsl_vector_get(ys,i);
			double fi=ann_response(network,xi);
			sum+=fabs(fi-yi);
		}
		return sum/xs->size;
	}

	gsl_vector* p=gsl_vector_alloc(network->params->size);
	gsl_vector_memcpy(p,network->params);
	qnewton(cost_function,p,1e-6);
	gsl_vector_memcpy(network->params,p);
	gsl_vector_free(p);
}
