#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<stdarg.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#ifndef DBL_EPSILON
#define DBL_EPSILON 2.22045e-16
#endif

double DELTA=sqrt(DBL_EPSILON);
// First introduce the functions from other routines. From utilities.c

void GS_decomp(gsl_matrix* A, gsl_matrix* R);
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
double norm(gsl_vector* x);
void vector_print(gsl_vector* vec);
void matrix_print(gsl_matrix* mat);


// Newton's method with numerical Jacobian
void Newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps){
    //Allocate
    int n = x->size;
    gsl_matrix* J = gsl_matrix_alloc(n,n);
    gsl_matrix* R = gsl_matrix_alloc(n,n);
    gsl_vector* Dx = gsl_vector_alloc(n);
    gsl_vector* xLamb = gsl_vector_alloc(n);
    gsl_vector* fx = gsl_vector_alloc(n);
    gsl_vector* Dfx = gsl_vector_alloc(n);
    f(x,fx);
    int N_MAX =0;
    //Jacobian
    while(gsl_blas_dnrm2(fx)>eps && N_MAX<10000) {
        N_MAX++;
        gsl_vector_memcpy(Dx, x);
        for (int j = 0; j < n; j++) {
            double xj = gsl_vector_get(x, j);
            gsl_vector_set(Dx, j, xj + sqrt(DBL_EPSILON));
            f(Dx, Dfx);
            for (int i = 0; i < n; i++) {
                double J_ij = (gsl_vector_get(Dfx, i) - gsl_vector_get(fx, i))/sqrt(DBL_EPSILON);
                gsl_matrix_set(J, i, j, -J_ij);
            }
        }
        double lambda = 1.;
        // GS from Gramschimidt.c
        GS_decomp(J,R);
        GS_solve(J, R, fx, Dx);
        for (int i = 0; i < n; i++) {
            double xi = gsl_vector_get(x, i);
            double Dxi = gsl_vector_get(Dx, i);
            gsl_vector_set(xLamb, i, xi + lambda * Dxi);
        }
        f(xLamb, Dfx);
        while (gsl_blas_dnrm2(Dfx) > (1 - lambda / 2) * gsl_blas_dnrm2(fx) && lambda > 1. / 64) {
            lambda /= 2;
            for (int i = 0; i < n; i++) {
                double xi = gsl_vector_get(x, i);
                double Dxi = gsl_vector_get(Dx, i);
                gsl_vector_set(xLamb, i, xi + lambda * Dxi);
            }
            f(xLamb, Dfx);
        }
        gsl_vector_memcpy(x, xLamb);
        f(x, fx);
    }
    //free
    gsl_vector_free(Dx);
    gsl_vector_free(xLamb);
    gsl_vector_free(fx);
    gsl_vector_free(Dfx);
    gsl_matrix_free(J);
}
