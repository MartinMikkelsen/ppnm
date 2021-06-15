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
    int n = 4;
    unsigned int SEED = time(NULL);
    srandom(SEED);

    gsl_matrix* A =gsl_matrix_alloc(n,n);
    gsl_matrix* V =gsl_matrix_alloc(n,n);
    gsl_matrix_set_identity(V);
 
    for(int i=0;i<n;i++)
    for(int j=i;j<n;j++){
	double aij=RND;
	gsl_matrix_set(A,i,j,aij);
	gsl_matrix_set(A,j,i,aij);
	}   
    printf("A random real symmetric matrix, A:\n");
    matrix_print(A);

    Jacobi_diag(A, V);
    
    printf("D is the diagonal matrix with the corresponding eigenvalues:\n");
   matrix_print(A);

   printf("Check that V^t*V=1 :\n");
   gsl_matrix* VtV = gsl_matrix_alloc(n, n);
   gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, V, V, 0, VtV);
   matrix_print(VtV);

    gsl_matrix* Acopy = gsl_matrix_alloc(n, n);
    gsl_matrix_memcpy(Acopy, A);

   printf("Check that V^t*A*V=D:\n");
   gsl_matrix* VtAV = gsl_matrix_alloc(n, n);
   gsl_matrix* temp = gsl_matrix_alloc(n, n);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Acopy, V, 0, temp);
   gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, V, temp, 0, VtAV);
   matrix_print(VtAV);

   printf("Check that V*D*V^t=A :\n");
   gsl_matrix* VDVt = gsl_matrix_alloc(n, n);
   gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, A, V, 0, temp);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, V, temp, 0, VDVt);
   matrix_print(VDVt);

    gsl_matrix_free(A);
    gsl_matrix_free(V);
    gsl_matrix_free(temp);
    gsl_matrix_free(VtAV);
    gsl_matrix_free(VDVt);

return 0;    
}

