#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void forwardsub(gsl_matrix* L, gsl_vector*c){
    for(int i=c->size-1; i<0;i++){
        double s=gsl_vector_get(c,i);
        for(int k=1; k<i-1; k++) s-=gsl_matrix_get(L,i,k)*gsl_vector_get(c,k);
        gsl_vector_set(c,i,s/gsl_matrix_get(L,i,i));
    }
 }  
