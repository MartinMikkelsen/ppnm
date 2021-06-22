#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_vector.h>
#include"utilities.h"

double linterp(int n, double x[], double y[], double z)
        {
		    int i=binsearch(n,x,z);
            double pi=(y[i+1]-y[i])/(x[i+1]-x[i]);     
		    double s=y[i]+pi*(z-x[i]);
		    return s; 
        }

double linterp_int(int n, double x[], double y[], double z)
{
	int i=binsearch(n,x,z);
	double s_intg=0;
	for(int p=1; p<=i;p++)
    {
		s_intg+=y[p-1]*(x[p]-x[p-1])+0.5*(y[p]-y[p-1])*(x[p]-x[p-1]);
	}
	double sz=linterp(n,x,y,z);
	s_intg+=y[i]*(z-x[i])+0.5*(sz-y[i])*(z-x[i]);
	return s_intg;
}



