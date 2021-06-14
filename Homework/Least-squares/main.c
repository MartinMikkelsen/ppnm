#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include<assert.h>


void GS_decomp(gsl_matrix* A, gsl_matrix* R){
  int m=A->size2;

  double DotP;

  for(int i=0; i<m; i++){
    gsl_vector_view col=gsl_matrix_column(A,i);
    gsl_matrix_set(R, i, i, gsl_blas_dnrm2(&col.vector));
    gsl_vector_scale(&col.vector, 1/gsl_matrix_get(R, i, i));

    for(int j=i+1; j<m; j++){
      gsl_vector_view col2=gsl_matrix_column(A,j);
      gsl_blas_ddot(&col.vector, &col2.vector, &DotP);
      gsl_matrix_set(R,i,j,DotP);
      gsl_blas_daxpy(-gsl_matrix_get(R,i,j),&col.vector, &col2.vector);
    }
  }
}

void backsub(gsl_matrix* M, gsl_vector* v){
  for(int i=v->size-1;i>=0;i--){
    double Sum=gsl_vector_get(v,i);
    for(int j=i+1; j<v->size; j++){
      Sum-=gsl_matrix_get(M,i,j)*gsl_vector_get(v,j);
    }
    gsl_vector_set(v,i,Sum/gsl_matrix_get(M,i,i));
  }
}

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* v, gsl_vector* w){
  gsl_blas_dgemv(CblasTrans, 1, Q, v, 0, w);
  backsub(R,w);
}

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B){
  gsl_vector* basevec=gsl_vector_alloc(Q->size2);
  for(int i=0; i<Q-> size2; i++){
    gsl_vector_set_basis(basevec,i);
    gsl_vector_view col=gsl_matrix_column(B,i);
    GS_solve(Q,R,basevec, &col.vector);
  }
  gsl_vector_free(basevec);
}

/* Implement a function that makes a least-squares fit—using your QR-decomposition routines—of a given data-set */

void lsfit(
        int m, double f(int i, double x),
        gsl_vector* x, gsl_vector* y, gsl_vector* dy,
        gsl_vector* c, gsl_matrix* S)
{
int n = x->size;

gsl_matrix *A = gsl_matrix_alloc(n,m);
gsl_vector *b = gsl_vector_alloc(n);
gsl_matrix *R = gsl_matrix_alloc(m,m);
gsl_matrix *invR = gsl_matrix_alloc(m,m);
gsl_matrix *I = gsl_matrix_alloc(m,m);

for(int i =0; i<n; i++){
        double xi = gsl_vector_get(x,i);
        double yi = gsl_vector_get(y,i);
        double dyi = gsl_vector_get(dy,i); assert(dyi>0);
        gsl_vector_set(b,i,yi/dyi);
        for(int k=0;k<m;k++)gsl_matrix_set(A,i,k,f(k,xi)/dyi);
        }
GS_decomp(A,R);
GS_solve(A,R,b,c);

gsl_matrix_set_identity(I);
GS_inverse(I,R,invR);
gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,invR,invR,0,S); /*matrix multiplication */

gsl_matrix_free(A);
gsl_vector_free(b);
gsl_matrix_free(R);
gsl_matrix_free(invR);
gsl_matrix_free(I);
}

double funs(int i, double x){
	switch(i){
		case  0: return 1  ; break;
		case  1: return x  ; break;
		case  2: return x*x; break;
		default: return NAN;
		}
	}

double fit(double x){
	double s=0;
	for(int k=0;k<m;k++)s+=gsl_vector_get(c,k)*funs(k,x);
	return s;
	}

double fit_plus(int i, double x){
	return fit(x)+gsl_vector_get(dc,i)*funs(i,x);
	}

double fit_minus(int i, double x){
	return fit(x)-gsl_vector_get(dc,i)*funs(i,x);
	}
double c1 =gsl_vector_get(c, 1);
double dc1=gsl_vector_get(dc,1);
double T=-1/c1*log(2.0);
double dT=dc1/c1/c1;

int main(){
	double x[]  = {1,   2,   3,  4,  6,  9,    10,   13,   15};
	double y[]  = {117, 100, 88, 72, 53, 29.5, 25.2, 15.2, 11.1};
	int n=sizeof(x)/sizeof(x[0]);
	double dy[n];
	for(int i=0;i<n;i++)dy[i]=0.05*y[i];

	for(int i=0;i<n;i++){
		dy[i]/=y[i];
		y[i]=log(y[i]);
		}

	gsl_vector* vx = gsl_vector_alloc(n);
	gsl_vector* vy = gsl_vector_alloc(n);
	gsl_vector* vdy = gsl_vector_alloc(n);
	for(int i=0;i<n;i++){
		gsl_vector_set(vx,i,x[i]);
		gsl_vector_set(vy,i,y[i]);
		gsl_vector_set(vdy,i,dy[i]);
		}
    int m=2;
	gsl_vector* c = gsl_vector_alloc(m);
	gsl_matrix* S = gsl_matrix_alloc(m,m);
	lsfit(m,funs,vx,vy,vdy,c,S);

    fprintf(stderr,"c=\n");
    for(int i=0;i<c->size;i++)fprintf(stderr,"%10g ",gsl_vector_get(c,i));
    fprintf(stderr,"\n");

    fprintf(stderr,"S=\n");
    for(int i=0;i<S->size1;i++){
    for(int j=0;j<S->size2;j++)
		    fprintf(stderr,"%10g ",gsl_matrix_get(S,i,j));
	    fprintf(stderr,"\n");
	    }

	    gsl_vector* dc = gsl_vector_alloc(m);
	    for(int k=0;k<m;k++){
		    double skk=gsl_matrix_get(S,k,k);
		    gsl_vector_set(dc,k,sqrt(skk));
		    }
        printf("# half-life = %.3g +- %.2g days\n",T,dT);

	    printf("# time log(activity) delta(log(activity))\n");
	    for(int i=0;i<n;i++)printf("%g %g %g\n",x[i],y[i],dy[i]);
	    printf("\n\n");

	    for(int i=0;i<m;i++){
		    printf("# time fit fit_plus fit_minus; k=%i\n",i);
		    for(double z=x[0],dz=(x[n-1]-x[0])/64;z<=x[n-1];z+=dz)
		    printf("%g %g %g %g\n",z,fit(z),fit_plus(i,z),fit_minus(i,z));
	    printf("\n\n");
	    }

return 0;
}
