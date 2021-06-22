#include <stdlib.h>
#include <assert.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>

typedef struct {int n; double* x, *y, *b, *c;} qspline;

qspline* qspline_alloc(int n, double* x, double* y)
{
    qspline *s = (qspline*)malloc(sizeof(qspline));
	s->b = (double*)malloc((n-1)*sizeof(double));
	s->c = (double*)malloc((n-1)*sizeof(double));
	s->x = (double*)malloc(n*sizeof(double));
	s->y = (double*)malloc(n*sizeof(double));
	s->n = n;

	int i=0;
	for(i=0;i<n;i++){
		s->x[i]=x[i];
		s->y[i]=y[i];
	}
	double dy[n-1], dx[n-1];
	for (i=0;i<n-1;i++){
		dx[i]=x[i+1]-x[i];
		dy[i]=(y[i+1]-y[i])/dx[i];
	}
	s->c[0]=0;
	for (i=0;i<n-2;i++){
		s->c[i+1]=(dy[i+1]-dy[i]-s->c[i]*dx[i])/dx[i+1];
	}
	s->c[n-2]/=2;
	for (i=n-3;i>=0;i--){
		s->c[i]=(dy[i+1]-dy[i]-s->c[i+1]*dx[i+1])/dx[i];
	}
	for (i=0;i<n-1;i++){
		s->b[i]=dy[i]-s->c[i]*dx[i];

	}
return s;
}

double qspline_eval(qspline *s, double input){
	int i = binsearch(s->n,s->x,input);
	double dx=input-s->x[i];
	double output=s->y[i]+dx*(s->b[i]+dx*s->c[i]);
return output;
}


double qspline_integ(qspline *s, double input)
{
	int i = binsearch(s->n,s->x,input);
	double output = 0;
	double dx=0;
	for (int p=1;p<=i;p++){
		dx=s->x[p]-s->x[p-1];
		output+=dx*(s->y[p-1]+dx*(s->b[p-1]/2+dx*s->c[p-1]/3));
	}
	dx=input-s->x[i];
	output+=dx*(s->y[i]+dx*(s->b[i]/2+dx*s->c[i]/3));

return output;
}

double qspline_deriv(qspline* s, double x_new)
{
   int i = binsearch(s->n, s->x, x_new);

   double deriv = s->b[i] + 2*(s->c[i])*(x_new - s->x[i]);
   return deriv;
}

void qspline_free(qspline *s)
{
	free(s->x);
	free(s->y);
	free(s->b);
	free(s->c);
	free(s);
}
