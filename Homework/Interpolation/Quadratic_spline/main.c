#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_interp.h>
#include <assert.h>
#include"utilities.h"


double linterp(int n, double x[], double y[], double z);

double linterp_int(int n, double x[], double y[], double z);

int main() 
{
	int n=15;
	double a=-2,b=2,x[n],y[n];

	for(int i=0;i<n;i++){
		x[i]=a+(b-a)*i/(n-1);
		y[i]=1/(1+pow(x[i],2));
	}

	printf("# index 0: data\n");
	for(int i=0;i<n;i++) printf("%g %g\n",x[i],y[i]);
	printf("\n\n");
	
	printf("# index 1: exact integral\n");
	for(int i=0;i<n;i++) printf("%g %g\n",x[i],atan(x[i])-atan(x[0]));
	printf("\n\n");

	gsl_interp* linear = gsl_interp_alloc(gsl_interp_linear    ,n);
	gsl_interp_init(linear,x,y,n);
	
	double z=1.0/5;
	printf("# index 2: interpolations\n");

	for(double i=x[0];i<=x[n-1];i+=z){
		printf("%g %g\n",i,linterp(n,x,y,i));
      	}	
  	
	printf("\n\n");
	
	printf("# index 3: integral\n");
	for(double i=x[0];i<=x[n-1];i+=z){
		printf("%g %g\n",i,linterp_int(n,x,y,i));
	}

	printf("\n\n");

//CODE TO COMPARE
//
	printf("# index 4: interpolations\n");
	for(double z=x[0];z<=x[n-1];z+=1./16){
		double interp_l=gsl_interp_eval(linear    ,x,y,z,NULL);
		printf("%g %g\n",z,interp_l);
	}
	printf("\n\n");

	printf("# index 5: integrals\n");
	for(double z=x[0];z<=x[n-1];z+=1./16){
		double integ_int=gsl_interp_eval_integ(linear    ,x,y,x[0],z,NULL);
		printf("%g %g\n",z,integ_int);
	}

	gsl_interp_free(linear);
return 0;
}
