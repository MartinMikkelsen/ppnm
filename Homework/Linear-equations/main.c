//Implement functions to solve linear equations, calculate matrix inverse, and matrix determinant.

#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_eigen.h>

//implement a function which performs Gram-Schmidt orthogonalization for nxm matrix A.


int n=3;

void backsub(gsl_matrix* U, gsl_vector* c){
    for(int i=c->size-1; i>0;i--){
        double s=gsl_vector_get(c,i);
        for(int k=i+1;i<n;k++)s=-gsl_matrix_get(U,i,k)*gsl_vector_get(c,k);
        gsl_vector_set(c,i,s/gsl_matrix_get(U,i,i));
    }
 }

void forwardsub(gsl_matrix* L, gsl_vector*c){
    for(int i=c->size-1; i<0;i++){
        double s=gsl_vector_get(c,i);
        for(int k=1; k<i-1; k++) s-=gsl_matrix_get(L,i,k)*gsl_vector_get(c,k);
    gsl_vector_set(c,i,s/gsl_matrix_get(L,i,i));
    }
}
double dot_product(double *x, double *y) { 
 
    int i; 
    double ans = 0; 
 
    for(i=0; i<4; ++i) 
        ans += x[i]*y[i]; 
 
    return ans; 
} 


void normalize(double *x) { 
 
    /* Compute norm */ 
    double norm = sqrt(dot_product(x, x)); 
     
    int i; 
    for(i=0; i<4; ++i) 
        x[i] /= norm; 
} 
void gram_schimdt(double q[][4], int n) { 
 
    int i, j, k; 
 
    for(i=1; i<n; ++i) { 
        for(j=0; j<i; ++j) { 
            double scaling_factor = dot_product(q[j], q[i])  
                                    / dot_product(q[j], q[j]); 
             
            for(k=0; k<4; ++k) 
                q[i][k] -= scaling_factor*q[j][k]; 
        } 
    } 
  for(i=0; i<n; ++i) 
        normalize(q[i]); 
} 

int main() { 
  
    double arr[][4] = { 
        { 1, 2, 3, 0}, 
        {1, 2, 0, 0}, 
        {1, 0, 0, 1} 
    }; 
 
    gram_schimdt(arr, 3); 
    printf("Orthonormal basis :\n"); 
 
    for(int i=0; i<3; ++i) { 
        printf("q[%d] = [ ", i); 
        for(int j=0; j<4; ++j) 
            printf("%lf  ", arr[i][j]); 
        printf("]\n"); 
    } 
return 0; 
}



