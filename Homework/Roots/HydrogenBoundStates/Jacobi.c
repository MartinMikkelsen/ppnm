//Implement the Jacobi algorithm for matrix diagonalizaion (also called matrix eigenvalue decomposition or eigendecomposition).
//In this algorithm we want to implement 3 functions, timesJ, Jtimes and jacobi_diag. 
//The goal is to generate a random symmetric matrix A, apply the routine to perform the eigenvalue-decomposition, A=VDVT (where V is the orthogonal matrix of eigenvectors and D is the diagonal matrix with the corresponding eigenvalues), and check that VTAV==D, VDVT==A, VTV==1,

#include<stdio.h>
#include<assert.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>


//multiplies the given matrix A and the Jacobi matrix J(p,q,θ) in O(n) operations, A←AJ
void timesJ(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
	for(int i=0;i<A->size1;i++){
		double new_aip=c*gsl_matrix_get(A,i,p)-s*gsl_matrix_get(A,i,q);
		double new_aiq=s*gsl_matrix_get(A,i,p)+c*gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,p,new_aip);
		gsl_matrix_set(A,i,q,new_aiq);
		}
}
//multiplies the Jacobi matrix J(p,q,θ) and the given matrix A in O(n) operations, A←JA
void Jtimes(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
	for(int j=0;j<A->size2;j++){
		double new_apj= c*gsl_matrix_get(A,p,j)+s*gsl_matrix_get(A,q,j);
		double new_aqj=-s*gsl_matrix_get(A,p,j)+c*gsl_matrix_get(A,q,j);
		gsl_matrix_set(A,p,j,new_apj);
		gsl_matrix_set(A,q,j,new_aqj);
		}
}
void Jacobi_diag(gsl_matrix*A, gsl_matrix* V){
    int n = A->size1;
    int changed;
do{
	changed=0;
	for(int p=0;p<n-1;p++)
	for(int q=p+1;q<n;q++){
		double apq=gsl_matrix_get(A,p,q);
		double app=gsl_matrix_get(A,p,p);
		double aqq=gsl_matrix_get(A,q,q);
		double theta=0.5*atan2(2*apq,aqq-app);
		double c=cos(theta),s=sin(theta);
		double new_app=c*c*app-2*s*c*apq+s*s*aqq;
		double new_aqq=s*s*app+2*s*c*apq+c*c*aqq;
		if(new_app!=app || new_aqq!=aqq) // do rotation
			{
			changed=1;
			timesJ(A,p,q, theta);
			Jtimes(A,p,q,-theta); // A←J^T*A*J
			timesJ(V,p,q, theta); // V←V*J
			}
	}
}while(changed!=0);
}
