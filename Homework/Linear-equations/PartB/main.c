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


int main() { 
   
   printf("====================\n");
   printf("====   TASK B   ====\n");
   printf("====================\n\n");


   // initialize random square matrix A and empty matrix R
   int N = 4;   
   unsigned int SEED = time(NULL);
   srandom(SEED);
   gsl_matrix* A = gsl_matrix_alloc(N, N);
   for (int i = 0; i < A->size1; ++i) {
      for (int j = 0; j < A->size2; ++j) {
         gsl_matrix_set(A, i, j, 5*((double)random())/RAND_MAX);
      }
   }

   // print original A matrix
   printf("This is the original A matrix:\n");
   matrix_print(A);

   // do QR decomposition and check result
   gsl_matrix* R = gsl_matrix_alloc(N, N);
   GS_decomp(A, R);

   gsl_matrix* QR = gsl_matrix_alloc(N, N);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, R, 0, QR);
   printf("Q*R=A :\n");
   matrix_print(QR);

   // Initialize empty matrix B and find inverse
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

   // free memory
   gsl_matrix_free(A);
   gsl_matrix_free(R);
   gsl_matrix_free(QR);
   gsl_matrix_free(B);


 
return 0;
}


