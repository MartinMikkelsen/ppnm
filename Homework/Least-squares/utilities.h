#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


double dot(gsl_vector* x, gsl_vector* y);

double norm(gsl_vector* x);

// print gsl vector
void vector_print(gsl_vector* vec);

// print gsl matrix
void matrix_print(gsl_matrix* mat);

// binary search for gsl_vector
int binsearch_vector(gsl_vector* x, double x_new);

// binary search for plain c array
int binsearch_array(int N, double* x, double x_new);

//Pass a set of functions
double funs(int i, double x);
