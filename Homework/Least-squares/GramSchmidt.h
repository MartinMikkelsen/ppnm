#include <gsl/gsl_vector.h> 
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h> 

void backsub(gsl_matrix* M, gsl_vector* v);

void forwardsub(gsl_matrix* L, gsl_vector* c);

void GS_decomp(gsl_matrix* A, gsl_matrix* R);

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

void GramSchmidtInverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);

