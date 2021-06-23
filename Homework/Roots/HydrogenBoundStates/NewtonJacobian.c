#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<float.h>
#ifndef DBL_EPSILON
#define DBL_EPSILON 2.22045e-16
#endif
double DELTA=sqrt(DBL_EPSILON);
int LIMIT=10000;

// Functions from different folders. 

void GS_decomp(gsl_matrix* A, gsl_matrix* R);

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);

double dot(gsl_vector* x, gsl_vector* y);

double norm(gsl_vector* x);

//Now Implement the Newton's method with simple backtracking linesearch algorithm where the derivatives in the Jacobian matrix are calculated numerically using finite differences. 

int Newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps)
{
	int n=x->size, steps=0;
	gsl_matrix* J  = gsl_matrix_alloc(n,n);
	gsl_matrix* R  = gsl_matrix_alloc(n,n);
	gsl_vector* fx = gsl_vector_alloc(n);
	gsl_vector* z  = gsl_vector_alloc(n);
	gsl_vector* fz = gsl_vector_alloc(n);
	gsl_vector* df = gsl_vector_alloc(n);
	gsl_vector* Dx = gsl_vector_alloc(n);
	while(steps<LIMIT){
		f(x,fx);
		for (int j=0;j<n;j++){ /* Jacobian */
			double xj=gsl_vector_get(x,j);
			gsl_vector_set(x,j,xj+DELTA);
			f(x,df);
			gsl_vector_sub(df,fx); /* df=f(x+DELTA)-f(x) */
			for(int i=0;i<n;i++)
				gsl_matrix_set(J,i,j,gsl_vector_get(df,i)/DELTA);
			gsl_vector_set(x,j,xj);
			}
		GS_decomp(J,R);
		GS_solve(J,R,fx,Dx);
		gsl_vector_scale(Dx,-1);
		double s=1;
		while(s>1./64){
			gsl_vector_memcpy(z,x);
			gsl_vector_add(z,Dx);
			f(z,fz);
			if( gsl_blas_dnrm2(fz)<(1-s/2)*gsl_blas_dnrm2(fx)) break;
			s*=0.5;
			gsl_vector_scale(Dx,0.5);
			}
		gsl_vector_memcpy(x,z);
		gsl_vector_memcpy(fx,fz);
		if( gsl_blas_dnrm2(Dx)<DELTA || gsl_blas_dnrm2(fx)<eps ) break;
		steps++; } //while
	gsl_matrix_free(J);
	gsl_vector_free(fx);
	gsl_vector_free(z);
	gsl_vector_free(fz);
	gsl_vector_free(df);
	gsl_vector_free(Dx);
return steps;
}
