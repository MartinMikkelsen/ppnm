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

// Question 1) QR decomposition Gram-Schmidt algorithm (A <- Q)
void GS_decomp(gsl_matrix* A, gsl_matrix* R);

// Question 2) Solve QRx = b by back-substitution
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

//Part B: Matrix inverse by Gram-Schmidt QR factorization

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);



//Given Q and R, calculates the inverse of the matrix A into the matrix B.
void GS_inverse(gsl_matrix* Q,gsl_matrix* R, gsl_matrix* B);

//Part C: Operations count for QR-decomposition and comparison with GSL


int main() { 
   
   // print task
   printf("====================\n");
   printf("====  TASK A.1  ====\n");
   printf("====================\n\n");

   // random gsl matrix A (N > M)
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
    //Here we want to test that the implementation actually works. We decompose a random matrix and check that QR=A.

   // Make a random matrix, A
   printf("This is the original A matrix:\n");
   matrix_print(A);

   // Empty gsl matrix R
   gsl_matrix* R = gsl_matrix_alloc(M, M);

   // Do QR 
   GS_decomp(A, R);

   printf("This is the matrix Q:\n");
   matrix_print(A);

   gsl_matrix* QtQ = gsl_matrix_alloc(M, M);
   gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, A, A, 0, QtQ);
   printf("This is Qt*Q:\n");
   matrix_print(QtQ);
   
   printf("This is R which is upper triangular:\n");
   matrix_print(R);

   gsl_matrix* QR = gsl_matrix_alloc(N, M);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, R, 0, QR);
   printf("This is Q*R (should be same as original A):\n");
   matrix_print(QR);

    
   // free memory  
   gsl_matrix_free(A);
   gsl_matrix_free(R);
   gsl_matrix_free(QtQ);
   gsl_matrix_free(QR);
 
return 0;
}


