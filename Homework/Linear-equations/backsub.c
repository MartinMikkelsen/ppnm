#include "backsub.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void backsub(gsl_matrix* U, gsl_vector* c){
  17     for(int i=c->size-1; i>0;i--){
  18         double s=gsl_vector_get(c,i);
  19         for(int k=i+1;i<n;k++)s=-gsl_matrix_get(U,i,k)*gsl_vector_get(c,k);
  20         gsl_vector_set(c,i,s/gsl_matrix_get(U,i,i));
  21     }
  22  }

