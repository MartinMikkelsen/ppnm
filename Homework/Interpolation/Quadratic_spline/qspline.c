#include <stdlib.h>
#include <assert.h>
#include "qspline.h"

typedef struct {int n; double *x, *y, *b, *c;} qspline;

qspline* qspline_alloc(int n,double* x,double* y){ //builds qspline
	qspline* s = malloc(sizeof(qspline));//spline
	s->b = malloc((n-1)*sizeof(double));  // b_i
	s->c = malloc((n-1)*sizeof(double));  // c_i
	s->x = malloc(n*sizeof(double));      // x_i
	s->y = malloc(n*sizeof(double));      // y_i
	s->n = n;
	for(int i=0;i<n;i++){
		s->x[i]=x[i];
		s->y[i]=y[i];
	}
	int i;
	double p[n-1], h[n-1];                  //VLA from C99
	for(i=0;i<n-1;i++){
		h[i]=x[i+1]-x[i];
		p[i]=(y[i+1]-y[i])/h[i];
	}
	s->c[0]=0;                                     //recursion up:
	for(i=0;i<n-2;i++)
		s->c[i+1]=(p[i+1]-p[i]-s->c[i]*h[i])/h[i+1];
	s->c[n-2]/=2;                                 //recursion down:
	for(i=n-3;i>=0;i--)
		s->c[i]=(p[i+1]-p[i]-s->c[i+1]*h[i+1])/h[i];
	for(i=0;i<n-1;i++)
		s->b[i]=p[i]-s->c[i]*h[i];
	return s;
}

double qspline_eval(qspline *s, double z){     //evaluates s(z)
	assert(z>=s->x[0] && z<=s->x[s->n-1]);
	int i=0, j=s->n-1;                     //binary search:
	while(j-i>1){
		int m=(i+j)/2;
		if(z>s->x[m]) i=m;
		else j=m;
	}
	double h=z-s->x[i];
	return s->y[i]+h*(s->b[i]+h*s->c[i]);
}//inerpolating polynomial

void qspline_free(qspline *s){ //free the allocated memory
	free(s->x); free(s->y); free(s->b); free(s->c); free(s);
}
