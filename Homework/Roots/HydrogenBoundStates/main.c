#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<stdio.h>
#include<float.h>
#include<math.h>

double E, R_max;

void vector_print(gsl_vector* vec);

void matrix_print(gsl_matrix* mat);

int Newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps);

int driver(void f(double x, gsl_vector* y, gsl_vector* dydx),double a, gsl_vector* ya, double b, double h,double acc, double eps);

int driver2(void f(double x, gsl_vector* y, gsl_vector* dydx),double a, gsl_vector* ya, double **Yal,int steps, double **Xal, double b, double h, double acc, double eps);

void schrodinger_eq(double x, gsl_vector* y, gsl_vector* dydx){
    gsl_vector_set(dydx,0,gsl_vector_get(y,1));
    double fx = gsl_vector_get(y,0);
    double ddy = -2.*(E+1./x)*fx;
    gsl_vector_set(dydx,1,ddy);
}

void ground_state(gsl_vector* eps, gsl_vector* fx)
{
    double a = 0.01;
    double acc = 0.001;
    double rel_acc = 0.001;
    E=gsl_vector_get(eps,0);
    gsl_vector* ya = gsl_vector_alloc(2);
    gsl_vector_set(ya,0,(a-a*a));
    gsl_vector_set(ya,1,(1-2.*a));
    driver(schrodinger_eq,a,ya,R_max,0.1,acc,rel_acc);
    gsl_vector_set(fx,0,gsl_vector_get(ya,0));
}

int main()
{
    gsl_vector* x = gsl_vector_alloc(2);
    gsl_vector* y = gsl_vector_alloc(1);
    gsl_vector_set(x,0,1.3);
    gsl_vector_set(x,1,1.3);
    gsl_vector_set(y,0,-3);

    R_max = 8.;
    Newton(ground_state,y,0.001);
    printf("# index 1:Lowest energy: E0 = %g\n",E);
    FILE* hydrogen = fopen("hydrogen.txt","w");
    
    double eps = 0.001;
    double abs = 0.001;
    double a=0.01, b=8.;

    gsl_vector* ya = gsl_vector_calloc(2);
    gsl_vector_set(ya,0,(a-a*a));
    gsl_vector_set(ya,1,(1-2.*a));
    int n=2, m=200;
    double* Yal = malloc(sizeof(double)*((m*n)));
    double* Xal = malloc(sizeof(double)*m);
    driver2(schrodinger_eq,a,ya,&Yal,m,&Xal,b,(b-a)/10,abs,eps);
    gsl_matrix_view Y = gsl_matrix_view_array(Yal, m, n);
    gsl_vector_view X = gsl_vector_view_array(Xal, m);
    for(int i = 0; i<m; i++){
        double xi = gsl_vector_get(&X.vector,i);
        if(xi!=0){
            double ana = xi*exp(-xi);
            double num = gsl_matrix_get(&Y.matrix,i,0);
            fprintf(hydrogen,"%g %g %g \n",xi, num, ana);
        }
    }
    fclose(hydrogen);
    free(Yal);
    free(Xal);
    gsl_vector_free(x);
    gsl_vector_free(y);
    gsl_vector_free(ya);

return 0;
}
