//Implement functions to solve linear equations, calculate matrix inverse, and matrix determinant.

#include<stdio.h>
#include<stdlib.h>
#include <time.h>

#include"utilities.h"
#include"GramSchmidt.h"


#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>

// Part A:  Solving linear equations using QR-decomposition by modified Gram-Schmidt orthogonalization
//implement a function which performs Gram-Schmidt orthogonalization for nxm matrix A. 

void GS_decomp(gsl_matrix* A, gsl_matrix* R);

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);


int main() { 
   

   int N = 4;   
   unsigned int SEED = time(NULL);
   srandom(SEED);
   gsl_matrix* A = gsl_matrix_alloc(N, N);
   for (int i = 0; i < A->size1; ++i) {
      for (int j = 0; j < A->size2; ++j) {
         gsl_matrix_set(A, i, j, 5*((double)random())/RAND_MAX);
      }
   }

   printf("This is the original A matrix:\n");
   matrix_print(A);

   gsl_matrix* R = gsl_matrix_alloc(N, N);
   GS_decomp(A, R);

   gsl_matrix* QR = gsl_matrix_alloc(N, N);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, R, 0, QR);
   printf("Q*R=A :\n");
   matrix_print(QR);

   gsl_matrix* B = gsl_matrix_alloc(N, N);
   GS_inverse(A, R, B);

   printf("B:\n");
   matrix_print(B);

   gsl_matrix* AB = gsl_matrix_alloc(N, N);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, QR, B, 0, AB);
   printf("A*B:\n");
   matrix_print(AB);

   gsl_matrix* BA = gsl_matrix_alloc(N, N);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, B, QR, 0, BA);
   printf("B*A:\n");
   matrix_print(BA);  

   gsl_matrix_free(A);
   gsl_matrix_free(R);
   gsl_matrix_free(QR);
   gsl_matrix_free(B);


 
return 0;
}


