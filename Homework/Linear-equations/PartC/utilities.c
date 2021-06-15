#include<stdio.h>
#include<assert.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_vector.h>
#include"utilities.h"


double dot(gsl_vector* x, gsl_vector* y){
	double x_dot_y;
	gsl_blas_ddot(x, y, &x_dot_y); //GSL function that calculates dot product
	return x_dot_y;
}

double norm(gsl_vector* x){
	return sqrt(dot(x,x));
}

// print gsl vector
void vector_print(gsl_vector* vec){

	for(int i = 0; i < vec->size; i++) {
      printf("%10.4f\n", gsl_vector_get(vec,i));
   }
	printf("\n");
}

// print gsl matrix
void matrix_print(gsl_matrix* mat){

	for(int i = 0; i < mat->size1; i++) {
	   for(int j = 0; j < mat->size2; j++) {
         printf("%10.4f ", gsl_matrix_get(mat, i, j));
      }
	   printf("\n");
   }
	printf("\n");
}

// binary search for gsl_vector
int binsearch_vector(gsl_vector* x, double x_new) {
   int N = x->size;
	assert(gsl_vector_get(x, 0) <= x_new && x_new <= gsl_vector_get(x, N-1));
	int i = 0;
   int j = N-1;

	while (j-i > 1) {
		int mid = (i + j)/2;
		if (x_new > gsl_vector_get(x, mid)) {
         i = mid;
      } else {
         j = mid;
      }
   }

	return i;
}

// binary search for plain c array
int binsearch_array(int N, double* x, double x_new) {
	assert(x[0] <= x_new && x_new <= x[N-1]);
	int i = 0;
   int j = N-1;

	while (j-i > 1) {
		int mid = (i + j)/2;
		if (x_new > x[mid]) {
         i = mid;
      } else {
         j = mid;
      }
   }

	return i;
}
