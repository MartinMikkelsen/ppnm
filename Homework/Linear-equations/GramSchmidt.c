#include <assert.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include "gramSchmidt.h"
#include "backsub.h"
#include "forwardsub.h"

//Need to implement normalization



void GramSchmidt_decomp(gsl_matrix* A, gsl_matrix* Q)

    int numOfRows   =  (int) A -> size1;
    int numOfCols   =  (int) A -> size2;
    assert( numOfRows >= numOfCols );

    gsl_vector* col   =   gsl_vector_alloc(numOfRows);
    *col              =   (gsl_matrix_column( matrixToQR, colId )).vector;
    double colNorm    =   gsl_blas_dnrm2(col);                                 
