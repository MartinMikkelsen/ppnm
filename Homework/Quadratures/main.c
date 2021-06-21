#include<math.h>
#include<assert.h>
#include<stdio.h>
#include"Recursive_Adaptive_Integrator.h"

#define SQR2 1.41421356237309504880

double intrun(double f(double), double a, double b, double acc, double eps, double f2, double f3, double nrec);


double integrate(double f(double), double a, double b,double acc,double eps);

int calls;

double Fa1(double x)
{
    return sqrt(x);
}
double Fa2(double x)
{
    return 4.*sqrt(1-(x)*(x));
}
int main()
{
    {

        double a = 0.;
        double b = 1.;
        double acc = 0.0001;
        double eps = 0.0001;

        double exact = 2/3.;
        double I1 = integrate(Fa1,a,b,acc,eps);
        double diff = exact-I1;
        printf("∫ dx √(x) from 0 to 1 \n");
        printf("-------------------------------\n");
        printf("  Value (exact = 2/3)       : %.10f\n"  , I1);
        printf("  Absolute accurary goal    : %.1e\n"   , acc);
        printf("  Relative accurary goal    : %.1e\n"   , eps);
        printf("  Difference (value-exact)  : %.10f\n"  , -diff);
        printf("-------------------------------\n");
    }
    {
        double a = 0.;
        double b = 1.;
        double acc = 0.0001;
        double eps = 0.0001;

        double exact = M_PI;
        double I2 = integrate(Fa2,a,b,acc,eps);
        double diff = exact-I2;
        printf("∫ dx 4√(1-x²) from 0 to 1 \n");
        printf("-------------------------------\n");
        printf("  Value (exact = π)         : %.10f\n"  , I2);
        printf("  Absolute accurary goal    : %.1e\n"   , acc);
        printf("  Relative accurary goal    : %.1e\n"   , eps);
        printf("  Difference (value-exact)  : %.10f\n"  ,-diff);
        printf("-------------------------------\n");
    }

return 0;
}
