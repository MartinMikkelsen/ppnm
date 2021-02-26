#include<stdio.h>
#include<math.h>
#include <gsl/gsl_integration.h>

double f (double t, void* params) {
  double t = *(double*)params;
  double f = exp(-pow(t,2));
  return f;
}

double myerror(double z){
	gsl_function F;
	F.function=&f;
	F.params =(void*)&t;
	int limit = 999;
	gsl_integration_workspace* w;
	w = gsl_integration_workspace_alloc (limit);
	double a=0, acc=1e-6, eps=1e-6, result, error;
	gsl_integration_qagiu(&F,a,z,acc,eps,limit,w,&result);
	gsl_integration_workspace_free(w);
	return result;	
}
int main(){
	for(double z=1.0/8;z<=5;z+=1.0/8)
		printf("%10g %10g \n",z,myerror(z));
}
