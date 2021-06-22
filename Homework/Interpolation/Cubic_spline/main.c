#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_interp.h>


typedef struct {
   int N;      // number of data points
   double* x;  // arrays
   double* y;  // .
   double* b;  // .
   double* c;  // .
   double* d;  // .
} cspline;

cspline* cspline_init(int N, double* x, double* y);

double cspline_eval(cspline* s, double x_new);

double cspline_eval_deriv(cspline* s, double x_new);

double cspline_eval_integral(cspline* s, double x_new);

void cspline_free(cspline* s);

int main() {

   int N = 10;
   gsl_vector* x = gsl_vector_alloc(N);
   gsl_vector* y = gsl_vector_alloc(N);

   FILE* x_file = fopen("xs.txt", "r");
   FILE* y_file = fopen("ys.txt", "r");
   gsl_vector_fscanf(x_file, x);
   gsl_vector_fscanf(y_file, y);

   double xa[x->size];
   double ya[y->size];
   for (int i = 0; i < x->size; ++i) {
      xa[i] = gsl_vector_get(x, i);
      ya[i] = gsl_vector_get(y, i);
   }

   gsl_interp* gsl_interp3 = gsl_interp_alloc(gsl_interp_cspline, x->size);
   gsl_interp_init(gsl_interp3, xa, ya, x->size);
   gsl_interp_accel* acc = gsl_interp_accel_alloc();

   cspline* s = cspline_init(N, xa, ya);
   double x_min = gsl_vector_get(x, 0);
   double x_max = gsl_vector_get(x, x->size - 1);
   double step = 0.05;
   int n = (x_max - x_min)/step;
   FILE* cspline_file = fopen("cspline.txt", "w");
   for (int i = 0; i <= n; ++i) {
      double x_new = i*step;
      double y_new = cspline_eval(s, x_new);
      double d_new = cspline_eval_deriv(s, x_new);
      double i_new = cspline_eval_integral(s, x_new);
      double y_gsl = gsl_interp_eval(gsl_interp3, xa, ya, x_new, acc);
      double d_gsl = gsl_interp_eval_deriv(gsl_interp3, xa, ya, x_new, acc);
      double i_gsl = gsl_interp_eval_integ(gsl_interp3, xa, ya, x_min, x_new, acc);
      fprintf(cspline_file, "%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n", x_new, y_gsl, y_new, d_gsl, d_new, i_gsl, i_new);
   }

   gsl_interp_accel_free(acc);
   gsl_interp_free(gsl_interp3);
   cspline_free(s);
   fclose(cspline_file);
   fclose(x_file);
   fclose(y_file);
   return 0;
}
