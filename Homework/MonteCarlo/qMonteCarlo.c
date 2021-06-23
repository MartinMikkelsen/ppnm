#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<complex.h>

//Implement a multidimensional Monte-Carlo integrator that uses low-discrepancy (quasi-random) sequences. The error could be estimated by using two different sequences. Compare the scaling of the error with pseudo-random Monte-Carlo integrator.

//Van der Corput and Halton sequences

double corput(int n, int base)
{
    double q = 0;
    double bk = (double)1/base;
    while(n>0)
    {
        q += (n % base)*bk;
        n /= base;
        bk /= base;
    }
    return q;
}
double Halton(int type, int n, int dim, int i) {
   if (type == 0) {
      int base[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67};
      int maxdim = sizeof(base)/sizeof(int);
      assert(dim <= maxdim);
      return corput(n, base[i]);
   }
   else {
      int base[] = {71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163};
      int maxdim = sizeof(base)/sizeof(int);
      assert(dim <= maxdim);
      return corput(n, base[i]);
   }
}

void pMonteCarlo(int dim, double f(int dim, double* x), double* a, double* b, double N, double* result2, double* error2)
{
   double V = 1;
   for (int i = 0; i < dim; ++i) V *= b[i] - a[i];

   double sum;
   double x[dim];
   double fx;

   sum = 0;
   for (int i = 0; i < N; ++i) {
      for(int j = 0; j < dim; ++j) x[j] = a[j] + Halton(0, i, dim, j)*(b[j] - a[j]);
      fx = f(dim, x);
      sum += fx;
   }
   double mean_a = sum/N;

   sum = 0;
   for (int i = 0; i < N; ++i) {
      for(int j = 0; j < dim; ++j) x[j] = a[j] + Halton(1, i, dim, j)*(b[j] - a[j]);
      fx = f(dim, x);
      sum += fx;
   }
   double mean_b = sum/N;

   *result2 = V*(mean_a + mean_b)/2;
   *error2  = fabs(V*mean_a - V*mean_b);
}
