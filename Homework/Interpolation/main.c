#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_interp.h>
#include <assert.h>

int binsearch(int n, double* x, double z){/* locates the interval for z by bisection */ 
	assert(n>1 && x[0]<=z && z<=x[n-1]);
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
		}
	return i;
	}

double linterp(int n, double x[], double y[], double z){
		int i=binsearch(n,x,z);
                double pi=(y[i+1]-y[i])/(x[i+1]-x[i]);
                
		double s=y[i]+pi*(z-x[i]);
		return s; 
              	}

double linterp_int(int n, double x[], double y[], double z){
	int i=binsearch(n,x,z);
	double s_intg=0;
	for(int p=1; p<=i;p++){
		s_intg+=y[p-1]*(x[p]-x[p-1])+0.5*(y[p]-y[p-1])*(x[p]-x[p-1]);
	}
	
	double sz=linterp(n,x,y,z);

	s_intg+=y[i]*(z-x[i])+0.5*(sz-y[i])*(z-x[i]);
	return s_intg;
}

int main (void) {
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

	gsl_interp* linear     = gsl_interp_alloc(gsl_interp_linear    ,n);
	gsl_interp_init(linear    ,x,y,n);
	
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
