#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<float.h>
#ifndef DBL_EPSILON
#define DBL_EPSILON 2.22045e-16
#endif
static const double DELTA=sqrt(DBL_EPSILON);


//Implement the quasi-Newton minimization method with numerical gradient, back-tracking linesearch, and a rank-1 update.

//First define the gradient, grad(f(x)) = df/dxi. 

void numeric_gradient(double F(gsl_vector*), gsl_vector*x, gsl_vector*grad)
{
	int n = x->size;
    double fx=F(x);
	for(int i=0;i<n;i++){
		double dx;
        double xi=gsl_vector_get(x,i);
		if(fabs(xi)<sqrt(DELTA)) dx=DELTA;
		else dx=fabs(xi)*DELTA;
		gsl_vector_set(x,i,xi+dx);
		gsl_vector_set(grad,i,(F(x)-fx)/dx);
		gsl_vector_set(x,i,xi);
	}
}

int qnewton(double F(gsl_vector* x), gsl_vector*x, double acc) 
{
	int dim=x->size;
    int nsteps=0;
    int nbad=0;
    int ngood=0;

	gsl_matrix* B=gsl_matrix_alloc(dim,dim);
	gsl_vector* gradient=gsl_vector_alloc(dim);
	gsl_vector* Dx=gsl_vector_alloc(dim);
	gsl_vector* z=gsl_vector_alloc(dim); 
	gsl_vector* gz=gsl_vector_alloc(dim);
	gsl_vector* y=gsl_vector_alloc(dim);
	gsl_vector* u=gsl_vector_alloc(dim);
	gsl_vector* a=gsl_vector_alloc(dim);
	gsl_matrix_set_identity(B);
	numeric_gradient(F,x,gradient);
	double fx=F(x);
    double fz; //fxpdx
	while(nsteps<1000){
		nsteps++;
		gsl_blas_dgemv(CblasNoTrans,-1,B,gradient,0,Dx);
		if(gsl_blas_dnrm2(Dx)<DELTA*gsl_blas_dnrm2(x))
			{fprintf(stderr,"qnewton: |Dx|<DELTA*|x|\n"); break;} 
		if(gsl_blas_dnrm2(gradient)<acc)
			{fprintf(stderr,"qnewton: |grad|<acc\n"); break;}
		double lambda=1;
		while(1){
			gsl_vector_memcpy(z,x);
			gsl_vector_add(z,Dx);
			fz=F(z);
			double sTg; gsl_blas_ddot(Dx,gradient,&sTg);
			if(fz<fx+0.01*sTg){ ngood++; break; } 
			if(lambda<DELTA){
				nbad++;
				gsl_matrix_set_identity(B);
				break;
				}
			lambda*=0.5;
			gsl_vector_scale(Dx,0.5);
		}
		numeric_gradient(F,z,gz);
		gsl_vector_memcpy(y,gz);
		gsl_blas_daxpy(-1,gradient,y);
		gsl_vector_memcpy(u,Dx);
		gsl_blas_dgemv(CblasNoTrans,-1,B,y,1,u); 
		double sTy,uTy;
		gsl_blas_ddot(Dx,y,&sTy);
		if(fabs(sTy)>1e-12){ 
			gsl_blas_ddot(u,y,&uTy);
			double gamma=uTy/2/sTy;
			gsl_blas_daxpy(-gamma,Dx,u); // u=u-gamma*s
			gsl_blas_dger(1.0/sTy,u,Dx,B); // B= B + u*s^T/s^T*y
			gsl_blas_dger(1.0/sTy,Dx,u,B); // B= B + s*u^T/s^T*y
		}
		gsl_vector_memcpy(x,z);
		gsl_vector_memcpy(gradient,gz);
		fx=fz;
	}
gsl_matrix_free(B);
gsl_vector_free(gradient);
gsl_vector_free(Dx);
gsl_vector_free(z);
gsl_vector_free(gz);
gsl_vector_free(y);
gsl_vector_free(u);
gsl_vector_free(a);
fprintf(stderr,"qnewton: nsteps=%i ngood=%i nbad=%i fx=%.1e\n",nsteps,ngood,nbad,fx);
return nsteps;
}
