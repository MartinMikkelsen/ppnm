#include<stdio.h>
#include<math.h>
#include <gsl/gsl_integration.h>

double f (double x, void * params) {
  double z = *(double *) params;
  double f = pow(x,z-1)*exp(-x);
  return f;
}

double mygamma(double z){
	gsl_function F;
	F.function=&f;
	F.params=(void*)&z;
	int limit=999;
	gsl_integration_workspace* w;
	w = gsl_integration_workspace_alloc (limit);
	double a=0,acc=1e-6,eps=1e-6,result,error;
	gsl_integration_qagiu(&F,a,acc,eps,limit,w,&result,&error);
	gsl_integration_workspace_free(w);
	return result;
}

int main(){
	for(double x=1.0/8;x<=5;x+=1.0/8)
		printf("%10g %10g\n",x,mygamma(x));
return 0;
}
