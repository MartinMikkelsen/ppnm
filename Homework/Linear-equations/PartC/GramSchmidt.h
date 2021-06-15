#include <gsl/gsl_vector.h> 
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h> 

void GramSchmidt(gsl_matrix* A, gsl_matrix* R);

void GramSchmidtSolver(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

void GramSchmidtInverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);
