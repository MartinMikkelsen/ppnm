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

//implement a function which performs Gram-Schmidt orthogonalization for nxm matrix A. 

// QR decomposition Gram-Schmidt algorithm (A <- Q)
void GS_decomp(gsl_matrix* A, gsl_matrix* R);

// Solve QRx = b by back-substitution
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

int main() { 
   
   // print task
   printf("====  TASK A.1  ====\n");

   // random gsl matrix A (N >= M)
   int N = 6;
   int M = 4;
   unsigned int SEED = time(NULL);
   srandom(SEED);
   gsl_matrix* A = gsl_matrix_alloc(N, M);
   for (int i = 0; i < A->size1; ++i) {
      for (int j = 0; j < A->size2; ++j) {
         gsl_matrix_set(A, i, j, 5*((double)random())/RAND_MAX);
      }
   }

   // print original A matrix
   printf("This is the original A matrix:\n");
   matrix_print(A);

   // initialize empty gsl matrix R
   gsl_matrix* R = gsl_matrix_alloc(M, M);

   // do QR decomposition and check results
   GS_decomp(A, R);

   printf("This is Q (in place of A):\n");
   matrix_print(A);

   gsl_matrix* QtQ = gsl_matrix_alloc(M, M);
   gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, A, A, 0, QtQ);
   printf("This is Qt*Q (should be identity):\n");
   matrix_print(QtQ);
   
   printf("This is R (should be upper triangular):\n");
   matrix_print(R);

   gsl_matrix* QR = gsl_matrix_alloc(N, M);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, R, 0, QR);
   printf("This is Q*R (should be same as original A):\n");
   matrix_print(QR);


   // print task
   printf("====================\n");
   printf("====  TASK A.2  ====\n");
   printf("====================\n\n");

   // initialize random square matrix a
   int n = 4;
   unsigned int seed = time(NULL);
   srandom(seed);
   gsl_matrix* a = gsl_matrix_alloc(n, n);
   for (int i = 0; i < a->size1; ++i) {
      for (int j = 0; j < a->size2; ++j) {
         gsl_matrix_set(a, i, j, 5*((double)random())/RAND_MAX);
      }
   }

   // print original A matrix
   printf("This is the matrix A in A*x = Q*R*x = b:\n");
   matrix_print(a);

   // also initialize empty matrix R and do QR decomposition
   gsl_matrix* r = gsl_matrix_alloc(n, n);
   GS_decomp(a, r);
   gsl_matrix* qr = gsl_matrix_alloc(n, n);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, a, r, 0, qr);

   // initialize random b and empty x for A*x = Q*R*x = b
   gsl_vector* x = gsl_vector_alloc(n);
   gsl_vector* b = gsl_vector_alloc(n);
   for (int i = 0; i < b->size; ++i) {
      gsl_vector_set(b, i, 5*((double)random())/RAND_MAX);
   }

   // solve system of equations and print results
   printf("This is the righ-hand side b:\n");
   vector_print(b);

   GS_solve(a, r, b, x);
   
   printf("This is the solution x:\n");
   vector_print(x);

   gsl_vector* qrx = gsl_vector_alloc(n);
   gsl_blas_dgemv(CblasNoTrans, 1, qr, x, 0, qrx);
   printf("This is Q*R*x (should be same as b):\n");
   vector_print(qrx);

   // free memory  
   gsl_matrix_free(A);
   gsl_matrix_free(R);
   gsl_matrix_free(QtQ);
   gsl_matrix_free(QR);
   gsl_matrix_free(a);
   gsl_matrix_free(r);
   gsl_matrix_free(qr);
   gsl_vector_free(qrx);
   gsl_vector_free(x);
   gsl_vector_free(b);

   return 0;
}


