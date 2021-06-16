#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort_vector.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_linalg.h>
#include"Jacobi.h"
#include"utilities.h"
#define RND ((double)rand()/RAND_MAX-0.5)*2


// A <- A*J
void timesJ(gsl_matrix* A, int p, int q, double theta);

// A <- J*A
void Jtimes(gsl_matrix* A, int p, int q, double theta);

// Jacobi rotation
void Jacobi_diag(gsl_matrix* A, gsl_matrix* V);

int main(){
   
    
   FILE* eigenfunctions_file = fopen("eigenfunctions.txt", "w");
   double x;
   int nplots = 3;

   //Now make new file containing eigenfunctions (numerical and analytical) and plot
    
    // Hamiltonian matrix
    int n = 100;
    double s = 1.0/(n+1);
    double k = -1/(s*s);
    double scale = sqrt(2.0/(n + 1));
    gsl_matrix* H = gsl_matrix_alloc(n, n);

    gsl_matrix_set(H, 0, 0, -2*k);
    gsl_matrix_set(H, 0, 1,  1*k);

    for (int i = 1; i < n-1; ++i) {
    gsl_matrix_set(H, i, i-1,  1*k);
      gsl_matrix_set(H, i, i  , -2*k);
      gsl_matrix_set(H, i, i+1,  1*k);
   }  

    gsl_matrix_set(H, n-1, n-2,  1*k);
    gsl_matrix_set(H, n-1, n-1, -2*k);

    // Use the Jacobi-implementation

    gsl_matrix* V = gsl_matrix_alloc(n, n);
    gsl_matrix_set_identity(V);
    Jacobi_diag(H, V);

   printf("Particle in a box energies:\n");
   for (int k = 0; k < 10; ++k) {
      double exact = M_PI*M_PI*(k + 1)*(k + 1);
      double calc  = gsl_matrix_get(H, k, k);
      printf("n = %2d: exact = %10.6f, calc = %10.6f\n", k+1, exact, calc);
   }
   
   for (int i = 0; i < n; ++i) {
      gsl_vector_view view = gsl_matrix_column(V, i);
      gsl_vector* vi = &view.vector;
      double first = gsl_vector_get(vi,0);
      if (first < 0) {
         gsl_vector_scale(vi, -1);
      }
   }

   fprintf(eigenfunctions_file, "0 ");
   for (int k = 0; k < nplots; ++k) {
      fprintf(eigenfunctions_file, "0 0 ");
   }
   fprintf(eigenfunctions_file, "\n");

   for(int i = 0; i < n; i++) {
      x = (i + 1.0)/(n + 1);
      fprintf(eigenfunctions_file, "%g ", x);
      for (int k = 0; k < nplots; ++k) {
         fprintf(eigenfunctions_file, "%g %g ", sin((k+1)*M_PI*x), gsl_matrix_get(V, i, k)/scale);
      }
      fprintf(eigenfunctions_file, "\n");
   }

   fprintf(eigenfunctions_file, "1 ");
   for (int k = 0; k < nplots; ++k) {
      fprintf(eigenfunctions_file, "0 0 ");
   }
   fprintf(eigenfunctions_file, "\n");
   

   fclose(eigenfunctions_file);
   gsl_matrix_free(H);
   gsl_matrix_free(V);

return 0;    
}

