#include<complex.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#define SQR2 1.41421356237309504880


void randomx(int dim, double* a, double* b, double* x);

void plainmc(int dim,double f(int dim,double* x),double* a,double* b,int N, double* result, double* error);

void Halton(int n, int dim, double *a, double *b, double *x);

void pMonteCarlo(int dim, double f(int dim, double* x), double* a, double* b, double N, double* result2, double* error2);

//From quadratures

double Fa1(double x)
{
    return sqrt(x);
};
double Fa2(double x)
{
    return 4.*sqrt(1-(x)*(x));
};
double Fa3(double x)
{
    return 1/(sqrt(x));
}
double Fa4(double x)
{
    return log(x)/sqrt(x);
}

double func(int dim, double* x)
{
   return 1.0/(1 - cos(x[0])*cos(x[1])*cos(x[2]))/(M_PI*M_PI*M_PI);
}

int main()
{
    {
        double result;
        double error;
        double result2;
        double error2;
        int dim = 1;
        int N = 1e3;
        double a[dim];
        double b[dim];
        //int dim=sizeof(a)/sizeof(a[0]);


        for (int i = 0; i < dim; ++i)
        {
            a[i] = 0;
            b[i] = 1;
        }
        plainmc(dim, Fa1, a, b, N, &result, &error);
        pMonteCarlo(dim, Fa1, a, b, N, &result2, &error2);
        double exact = 2/3.;
        printf("-------------------------------\n");
        printf("∫ dx √(x) from 0 to 1 using plain Monte Carlo \n");
        printf("Value = %.10f \n", result);
        printf("Error = %.10f \n", error);
        printf("Exact = %.10f \n", exact);
        printf("Diff  = %.10f \n", exact-result);
        printf("-------------------------------\n");
        printf("∫ dx √(x) from 0 to 1 using pseudo Monte Carlo \n");
        printf("Value = %.10f \n", result2);
        printf("Error = %.10f \n", error2);
        printf("Exact = %.10f \n", exact);
        printf("Diff  = %.10f \n", exact-result2);

    }
    {
        double result;
        double error;
        double result2;
        double error2;
        int dim = 1;
        int N = 1e3;
        double a[dim];
        double b[dim];
        for (int i = 0; i < dim; ++i)
        {
            a[i] = 0;
            b[i] = 1;
        }
        plainmc(dim, Fa2, a, b, N, &result, &error);
        pMonteCarlo(dim, Fa2, a, b, N, &result2, &error2);
        double exact = M_PI;
        printf("-------------------------------\n");
        printf("∫ dx 4√(1-x²) from 0 to 1 using plain Monte Carlo \n");
        printf("Value = %.10f \n", result);
        printf("Error = %.10f \n", error);
        printf("Exact = %.10f \n", exact);
        printf("Diff  = %.10f \n", exact-result);
        printf("-------------------------------\n");
        printf("∫ dx 4√(1-x²) from 0 to 1 using pseudo Monte Carlo \n");
        printf("Value = %.10f \n", result2);
        printf("Error = %.10f \n", error2);
        printf("Exact = %.10f \n", exact);
        printf("Diff  = %.10f \n", exact-result2);
        printf("-------------------------------\n");

        //also compare the scaling using this example

        FILE* errorscaling = fopen("errorscaling.txt", "w");
        for (int N = 1e5; N <= 2e6; N += 1e5) {
        plainmc (dim, Fa2, a, b, N, &result , &error);
        pMonteCarlo(dim, Fa2, a, b, N, &result2, &error2);
        fprintf(errorscaling, "%10.2e %20.15e %20.15e\n", (double)N, fabs(result - exact), fabs(result2 - exact));
        }


    }
    {
        double result;
        double error;
        double result2;
        double error2;
        int dim = 1;
        int N = 1e3;
        double a[dim];
        double b[dim];
        for (int i = 0; i < dim; ++i)
        {
            a[i] = 0;
            b[i] = 1;
        }
        plainmc(dim, Fa3, a, b, N, &result, &error);
        pMonteCarlo(dim, Fa3, a, b, N, &result2, &error2);
        double exact = 2.;
        printf("-------------------------------\n");
        printf("∫ dx/√(x) from 0 to 1 using plain Monte Carlo \n");
        printf("Value = %.10f \n", result);
        printf("Error = %.10f \n", error);
        printf("Exact = %.10f \n", exact);
        printf("Diff  = %.10f \n", exact-result);
        printf("-------------------------------\n");


    }
    {
        double result;
        double error;
        double result2;
        double error2;
        int dim = 1;
        int N = 1e3;
        double a[dim];
        double b[dim];
        for (int i = 0; i < dim; ++i)
        {
            a[i] = 0;
            b[i] = 1;
        }
        plainmc(dim, Fa4, a, b, N, &result, &error);
        pMonteCarlo(dim, Fa4, a, b, N, &result2, &error2);
        double exact = -4.;
        printf("-------------------------------\n");
        printf("∫ ln(x)/√(x) from 0 to 1 using plain Monte Carlo \n");
        printf("Value = %.10f \n", result);
        printf("Error = %.10f \n", error);
        printf("Exact = %.10f \n", exact);
        printf("Diff  = %.10f \n", exact-result);
        printf("-------------------------------\n");
    }
    {
        double result;
        double error;
        double result2;
        double error2;
        int dim = 3;
        int N = 1e3;
        double a[dim];
        double b[dim];
        for (int i = 0; i < dim; ++i)
        {
            a[i] = 0;
            b[i] = M_PI;
        }
        plainmc(dim, func, a, b, N, &result, &error);
        pMonteCarlo(dim, func, a, b, N, &result2, &error2);
        double exact = 1.3932039296856768591842462603255;
        printf("-------------------------------\n");
        printf("1/π³ ∫ dxdydz ([1-cos(x)cos(y)cos(z)]-1) from 0 to π using plain Monte Carlo \n");
        printf("Value = %.10f \n", result);
        printf("Error = %.10f \n", error);
        printf("Exact = %.10f \n", exact);
        printf("Diff  = %.10f \n", exact-result);
        printf("-------------------------------\n");
    }
return 0;
}
