#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_eigen.h>

int hh(){
  int n=3;
  
  gsl_matrix* A=gsl_matrix_alloc(n,n);
  gsl_matrix* Acopy=gsl_matrix_alloc(n,n);
  gsl_vector* b=gsl_vector_alloc(n);
  gsl_vector* x=gsl_vector_alloc(n);
  gsl_vector* y=gsl_vector_calloc(n);

  double A_data[3][3]={{6.13, -2.90, 5.86},{8.08, -6.31, -3.89},{-4.36, 1.00, 0.19}};
  for (int i=0;i<n;i++){
    for (int j=0;j<n;j++){
      gsl_matrix_set(A,i,j,A_data[i][j]);
    }
  }
  gsl_matrix_memcpy(Acopy,A);
  double b_data[]={6.23, 5.37, 2.29};
  for (int i=0;i<n;i++){
    gsl_vector_set(b,i,b_data[i]);
  }

  gsl_linalg_HH_solve(Acopy,b,x);
  gsl_blas_dgemv(CblasNoTrans,1,A,x,0,y);

  printf("Right hand side x is x=\n");

  gsl_vector_fprintf(stdout, x, "%g");
  printf("Check: A*x=\n");
  gsl_vector_fprintf(stdout, y, "%g");

  gsl_matrix_free(A);
  gsl_matrix_free(Acopy);
  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_vector_free(y);
  return 0;
}

int hilbert(){
  int n=4;

  gsl_matrix* A=gsl_matrix_alloc(n,n);
  gsl_matrix* Acopy=gsl_matrix_alloc(n,n);
  gsl_vector* b=gsl_vector_alloc(n);
  gsl_vector* x=gsl_vector_alloc(n);
  gsl_vector* y=gsl_vector_calloc(n);

  for (int i=0;i<n;i++){
    for (int j=0; j<n;j++){
      double Aij=1.0/(i+j+1);
      gsl_matrix_set(A,i,j,Aij);
    }
  }
  gsl_eigen_symmv_workspace * w= gsl_eigen_symmv_alloc(n);
  gsl_eigen_symmv (A,x,Acopy,w);
  gsl_eigen_symmv_free(w);
  gsl_eigen_symmv_sort(x,Acopy,GSL_EIGEN_SORT_ABS_ASC);

  for (int i=0;i<n;i++){
    double x_i=gsl_vector_get (x,i);
    gsl_vector_view evec_i = gsl_matrix_column (Acopy, i);

    printf("Eigenval = %g\n", x_i);
    printf("Eigenvector = \n");
    gsl_vector_fprintf(stdout, &evec_i.vector,"%g");
  }

  gsl_matrix_free(A);
  gsl_matrix_free(Acopy);
  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_vector_free(y);
  return 0;
}
int main(){

  printf("######################### Linear System Solver #########################\n");

  hh();

  printf("############################ Hilbert matrix ############################\n");

  hilbert();

  return 0;
}
