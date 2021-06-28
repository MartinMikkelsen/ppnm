#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <assert.h>
#include"utilities.h"
#include"GramSchmidt.h"
#include"lsfit.h"

void lsfit(
        int m, double f(int i, double x),
        gsl_vector* x, gsl_vector* y, gsl_vector* dy,
        gsl_vector* c, gsl_matrix* S);

void matrix_print(gsl_matrix* mat);

int main(){
    FILE* outfile = fopen("output.txt","w");
    double t[] = {1,  2,  3, 4, 6, 9,   10,  13,  15};
    double y[]=  {117,100,88,72,53,29.5,25.2,15.2,11.1};
    int n=sizeof(t)/sizeof(t[0]);
    double dy[n];
    double ylog[n];
    double dlogy[n];
    for(int i=0;i<n;i++){
		dy[i] = 0.05 * y[i];
		dlogy[i] = dy[i]/y[i];
		ylog[i] = log(y[i]);
	}
    gsl_vector* lt = gsl_vector_alloc(n);
	gsl_vector* ly = gsl_vector_alloc(n);
	gsl_vector* ldy = gsl_vector_alloc(n);
	for(int i=0;i<n;i++){
		gsl_vector_set(lt,i,t[i]);
		gsl_vector_set(ly,i,ylog[i]);
		gsl_vector_set(ldy,i,dlogy[i]);
	}
    double m = 2;
    gsl_matrix* cov = gsl_matrix_alloc(m,m);
	gsl_vector* c = gsl_vector_alloc(m);
	lsfit(m, funs, lt, ly, ldy, c, cov);
	gsl_vector* dfit  = gsl_vector_alloc(m);
	for(int k = 0; k < m; k++){
		gsl_vector_set(dfit, k, sqrt(gsl_matrix_get(cov, k, k)));
	}

    fprintf(outfile, "---------------------------\n");
    fprintf(outfile,"The fit is given by:\n");
	for(int i = 0; i < m; i++) {
		fprintf(outfile, "%g*x^%i", gsl_vector_get(c, i), i);
	}
	fprintf(outfile, "\n\nThe value of c1 is %g.\n \n",gsl_vector_get(c, 1));
	double t0 = -log(2)/gsl_vector_get(c, 1);
	fprintf(outfile, "The estimation of the value for the half-life is %g days +/- %g days =  [%g %g] days\n\n",t0,t0*0.05,t0+t0*0.05,t0-t0*0.05);
    fprintf(outfile, "The correct value is %g days. The two half-lives do not agree.\n\n", 3.6319);
	fprintf(outfile, "The covariance matrix is\n\n");
    fprintf(outfile, "    [");
	for (int i = 0; i < m ; i++){
		for(int j = 0; j < m; j++){
		fprintf(outfile, "\t%3.g", gsl_matrix_get(cov, i, j));
		}
		if (i < m - 1) {
			fprintf(outfile, "\n");
		} else {
			fprintf(outfile, "\t]\n");
		}
	}
    fclose(outfile);
    double fit(double x){
		double sum =0;
		for(int i=0;i<m;i++){
			sum += gsl_vector_get(c, i)*funs(i, x);
		}
		return sum;
	}

    FILE* resultfile = fopen("fitparams.txt", "w");
	// Data points
	fprintf(resultfile,"# Data points\n");
	for(int i = 0; i < n; i++) {
		fprintf(resultfile,"%g %g %g\n", t[i], ylog[i], dlogy[i]);
	}

	fprintf(resultfile,"\n\n# Fit\n");
	for(int i = 0; i < n; i++) {
		fprintf(resultfile,"%g %g \n", t[i], fit(t[i]));
	}

	fprintf(resultfile,"\n\n# Fit with uncertainties= +delta(c0)+delta(c1)\n");
	for(int i = 0; i < n; i++) {
		fprintf(resultfile,"%g %g \n", t[i], fit(t[i]) + gsl_vector_get(dfit, 0)*funs(0, t[i]) + gsl_vector_get(dfit, 1)*funs(1, t[i]));
	}

	fprintf(resultfile,"\n\n# Fit with uncertainties = -delta(c0) - delta(c1)\n");
	for(int i = 0; i < n; i++) {
		fprintf(resultfile,"%g %g \n", t[i], fit(t[i]) - gsl_vector_get(dfit, 0)*funs(0, t[i]) - gsl_vector_get(dfit, 1)*funs(1, t[i]));
	}

	fclose(resultfile);

	gsl_matrix_free(cov);
	gsl_vector_free(c);
	gsl_vector_free(dfit);
return 0;
}
