#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_vector.h>


int binsearch(int n, double* x, double z){/* locates the interval for z by bisection */
	assert(n>1 && x[0]<=z && z<=x[n-1]);
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
		}
	return i;
	}

typedef struct {
   int N;      // number of data points
   double* x;  // arrays
   double* y;  // .
   double* b;  // .
   double* c;  // .
   double* d;  // .
} cspline;


cspline* cspline_init(int N, double* x, double* y) {
   cspline* s = (cspline*)malloc(sizeof(cspline));
   s->N = N;
   s->x = (double*)malloc(N*sizeof(double));
   s->y = (double*)malloc(N*sizeof(double));
   s->b = (double*)malloc(N*sizeof(double));
   s->c = (double*)malloc((N - 1)*sizeof(double));
   s->d = (double*)malloc((N - 1)*sizeof(double));
   for (int i = 0; i < N; ++i) {
      s->x[i] = x[i];
      s->y[i] = y[i];
   }
   int i;
   double dx[N-1];
   double  p[N-1];
   for (i = 0; i < N - 1; ++i) {
      dx[i] = x[i+1] - x[i];
      assert(dx[i] > 0);
      p[i]  = (y[i+1] - y[i])/dx[i];
   }
   double D[N];
   double B[N];
   double Q[N-1];
   D[0]   = 2;
   D[N-1] = 2;
   B[0]   = 3*p[0];
   B[N-1] = 3*p[N-2];
   Q[0]   = 1;

   for (i = 0; i < N - 2; ++i) {
      D[i+1] = 2*dx[i]/dx[i+1] + 2;
      B[i+1] = 3*(p[i] + p[i+1]*dx[i]/dx[i+1]);
      Q[i+1] = dx[i]/dx[i+1];
   }
   for (i = 1; i < N; ++i) {
      D[i] -= Q[i-1]/D[i-1];
      B[i] -= B[i-1]/D[i-1];
   }
   s->b[N-1] = B[N-1]/D[N-1];

   for (i = N - 2; i >= 0; --i) {
      s->b[i] = (B[i] - Q[i]*s->b[i+1])/D[i];
   }

   for (i = 0; i < N - 1; ++i) {
      s->c[i] = (-2*s->b[i] - s->b[i+1] + 3*p[i])/dx[i];
      s->d[i] = (s->b[i] + s->b[i+1] - 2*p[i])/dx[i]/dx[i];
   }

   return s;
}

double cspline_eval(cspline* s, double x_new) {
   int i = binsearch(s->N, s->x, x_new);
   double dx = x_new - s->x[i];
   double y_new = s->y[i] + s->b[i]*dx + s->c[i]*pow(dx, 2) + s->d[i]*pow(dx, 3);
   return y_new;
}

double cspline_eval_deriv(cspline* s, double x_new) {

   int i = binsearch(s->N, s->x, x_new);
   double dx = x_new - s->x[i];
   double deriv = s->b[i] + 2*s->c[i]*dx + 3*s->d[i]*pow(dx, 2);
   return deriv;

}

double cspline_eval_integral(cspline* s, double x_new) {
   int i = binsearch(s->N, s->x, x_new);
   double integral = 0;
   for (int j = 0; j < i; ++j) {
      double dx = s->x[j+1] - s->x[j];
      integral += s->y[j]*dx + s->b[j]*pow(dx, 2)/2 + s->c[j]*pow(dx, 3)/3 + s->d[j]*pow(dx, 4)/4;
   }
   double dx = x_new - s->x[i];
   integral += s->y[i]*dx + s->b[i]*pow(dx, 2)/2 + s->c[i]*pow(dx, 3)/3 + s->d[i]*pow(dx, 4)/4;
   return integral;

}
void cspline_free(cspline* s) {
   free(s->x);
   free(s->y);
   free(s->b);
   free(s->c);
   free(s->d);
   free(s);
}
