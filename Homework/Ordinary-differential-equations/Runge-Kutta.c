#include<stdio.h>
#include<assert.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include"utilities.h"
#include"Runge-Kutta.h"
#define FOR(i) for(int i=0;i<n;i++)


//Implement an embedded Runge-Kutta stepper, rkstepXY -- where XY is the imbedded method. The function must estimate the values y_k(t+h) and store them in a vector* yh, estimate the errors delta-y_k and store in a vector err. Using the same approach as in the book but with gsl. 


void rkstep12(void f(double x, gsl_vector* y, gsl_vector* dydx)
              , double x, gsl_vector* yx, double h,
              gsl_vector* yh, gsl_vector* err) {
    int n = yx->size;
    gsl_vector *k0 = gsl_vector_alloc(n);
    gsl_vector *k1 = gsl_vector_alloc(n);
    gsl_vector *yt = gsl_vector_alloc(n);
    //Calculate first order
    f(x, yx, k0);
    for (int i = 0; i < n; i++) {
        double yxi = gsl_vector_get(yx, i);
        double k0i = gsl_vector_get(k0, i);
        double yti = yxi + h / 2 * k0i;
        gsl_vector_set(yt, i, yti);
    }
    //Calculate second order
    f(x + h / 2, yt, k1);
    for (int i = 0; i < n; i++) {
        double yxi = gsl_vector_get(yx, i);
        double k1i = gsl_vector_get(k1, i);
        double yhi = yxi + h * k1i;
        gsl_vector_set(yh, i, yhi);
    }
    //Error estimate
    for (int i = 0; i < n; i++) {
        double k0i = gsl_vector_get(k0, i);
        double k1i = gsl_vector_get(k1, i);
        double erri = (k0i - k1i) * h / 2;
        gsl_vector_set(err,i,erri);
    }

    gsl_vector_free(k0);
    gsl_vector_free(k1);
    gsl_vector_free(yt);
}

void rkstep23(void (*f)(int n, double x, double* y, double* dydx), int n, double x, double* y_curr, double h, double* y_next, double* err) {
   
   // declare variables
   double k1[n];
   double k2[n];
   double k3[n];
   double k4[n];
   double y_temp[n];

   // first point
   f(n, x, y_curr, k1);
   for (int i = 0; i < n; ++i) {
      y_temp[i] = y_curr[i] + (1.0/2)*k1[i]*h;
   }

   // second point
   f(n, x + (1.0/2)*h, y_temp, k2);
   for (int i = 0; i < n; ++i) {
      y_temp[i] = y_curr[i] + (3.0/4)*k2[i]*h;
   }

   // third point and second order estimate of y(x + h)
   f(n, x + (3.0/4)*h, y_temp, k3);
   for (int i = 0; i < n; ++i) {
      y_temp[i] = y_curr[i] + ((2.0/9)*k1[i] + (1.0/3)*k2[i] + (4.0/9)*k3[i])*h;
   }

   // fourth point, third order estimate of y(x + h), and error estimate
   f(n, x + h, y_temp, k4);
   for (int i = 0; i < n; ++i) {
      y_next[i] = y_curr[i] + ((7.0/24)*k1[i] + (1./4)*k2[i] + (1.0/3)*k3[i] + (1.0/8)*k4[i])*h;
      err[i] = y_next[i] - y_temp[i];
   }
  

}

//Implement an adaptive-step-size driver routine.
//
int driver(
	void (*f)(int n, double x, double* y, double* dydx), // right-hand-side of dy/dx = f(x, y) 
    int n,          // size of vectors 
	double  a,      // the start-point a 
	double  b,      // the end-point of the integration 
	double* ya,     // y(a) 
	double* yb,     // y(b) to be calculated 
	double  h,      // initial step-size 
	double  acc,    // absolute accuracy goal 
	double  eps,    // relative accuracy goal 
   char*   outfile // trajectory file
) { 

   // declare variables
   int k = 0;     // step counter
   double x;      // current x
   double y[n];   // current y vector
   double yh[n];  // estimate of y(x + h)
   double err[n]; // error estimate vector
   double sum;    // intermediate variable
   double tau;    // local tolerance
   FILE* file = fopen(outfile, "w");

   // prepare first step and write it to file
   x = a;
   fprintf(file, "%20.12e ", x);
   for (int i = 0; i < n; ++i) {
      y[i] = ya[i];
      fprintf(file, "%20.12e ", y[i]);
   }
   fprintf(file, "\n");


   // loop until endpoint is reached
   while (x < b) {
      //printf("===============\n");
      //printf("Step no. %3d\n", k);
      //printf("x = %.4f\ny =\n", x);
      //cvector_print(y, n);

      // avoid overstepping
      if (x + h > b) {
         h = b - x;
      }

      // make step
      rkstep23(f, n, x, y, h, yh, err);
      
      // calculate local error as norm of error estimate vector
      sum = 0;
      for (int i = 0; i < n; ++i) {
         sum += err[i]*err[i];
      }  
      double e = sqrt(sum);
      //printf("Error = %.5f\n", e);

      // calculate norm of estimated y(x + h)
      sum = 0;
      for (int i = 0; i < n; ++i) {
         sum += yh[i]*yh[i];
      }
      double norm = sqrt(sum);
      //printf("Norm = %.5f\n", e);

      // calculate local tolerance
      tau = (norm*eps + acc)*sqrt(h/(b-a));

      // accept or reject step?
      if (e < tau) {
         k++;     // increment counter
         x += h;  // increment x

         // increment elements of y vector and print to file
         fprintf(file, "%20.12e ", x);
         for (int i = 0; i < n; ++i) {
            y[i] = yh[i];
            fprintf(file, "%20.12e ", y[i]);
         }
         fprintf(file, "\n");

      }
      
      // adjust step size
      if (e > 0) { h *= pow(tau/e, 0.25)*0.95; }
      else       { h *= 2; } 
     
   } /* end while */

   // set final y vector
   for (int i = 0; i < n; ++i) {
      y[i] = yh[i];
   }

   //printf("===============\n");
   //printf("Step no. %3d\n", k);
   //printf("x = %.4f\ny =\n", x);
   //cvector_print(y, n);

   // close output file
   fclose(file);

   // return number of steps
   return k;

}









