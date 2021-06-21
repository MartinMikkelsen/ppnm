#include<math.h>
#include<assert.h>
#include<gsl/gsl_integration.h>
#include<stdio.h>
#include"Recursive_Adaptive_Integrator.h"
#include"Clenshaw–Curtis.h"
#define SQR2 1.41421356237309504880

double intrun(double f(double), double a, double b, double acc, double eps, double f2, double f3, double nrec);

double integrate(double f(double), double a, double b,double acc,double eps);

double clenshaw(double f(double),double a,double b, double acc,double eps);

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
        printf("  Absolute accurary goal    : %.1e\n"   , acc);
        printf("  Relative accurary goal    : %.1e\n"   , eps);
        printf("-------------------------------\n");
        printf("∫ dx √(x) from 0 to 1 \n");
        printf("-------------------------------\n");
        printf("  Value (exact = 2/3)       : %.18f\n"  , I1);
        printf("  Difference (value-exact)  : %.18f\n"  , -diff);
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
        printf("  Value (exact = π)         : %.18f\n"  , I2);
        printf("  Difference (value-exact)  : %.18f\n"  ,-diff);
        printf("-------------------------------\n");
    }
    {
        double a = 0.;
        double b = 1.;
        double acc = 0.0001;
        double eps = 0.0001;
        double exact = 2.;

        double I3 = clenshaw(Fa3,a,b,acc,eps);
        double diff = exact-I3;
        printf("∫ dx 1/√(x) from 0 to 1. With Clenshaw–Curtis transformation\n");
        printf("------------------------------------------------------------\n");
        printf("  Value (exact = 2)                         : %.18f\n"  , I3);
        printf("  Difference (value-exact)                  : %.18f\n"  , -diff);
        printf("------------------------------------------------------------\n");
    }
    {
        double a = 0.;
        double b = 1.;
        double acc = 0.0001;
        double eps = 0.0001;
        double exact = -4.;

        double I4 = clenshaw(Fa4,a,b,acc,eps);
        double diff = exact-I4;
        printf("∫ dx ln(x)/√(x) from 0 to 1. With Clenshaw–Curtis transformation\n");
        printf("----------------------------------------------------------------\n");
        printf("  Value (exact = -4)                        : %.18f\n"  , I4);
        printf("  Difference (value-exact)                  : %.18f\n"  , diff);
        printf("----------------------------------------------------------------\n");
    }
    {
        double a = 0.;
        double b = 1.;
        double acc = 0.0001;
        double eps = 0.0001;
        double exact = M_PI;

        double I5 = clenshaw(Fa2,a,b,acc,eps);
        double diff = exact-I5;
        printf("∫ dx 4√(1-x²) from 0 to 1. With Clenshaw-Curtis transformation \n");
        printf("--------------------------------------------------------------\n");
        printf("  Value (exact = π)                           : %.18f\n"  , I5);
        printf("  Difference (value-exact)                    : %.18f\n"  , diff);
        printf("---------------------------------------------------------------\n");
    }
    {
        printf("Now compare to GSL \n");
        printf("=================================================================\n");
        gsl_integration_workspace * w
        = gsl_integration_workspace_alloc (1000);

        double a = 0.;
        double b = 1.;
        double acc = 0.0001;
        double eps = 0.0001;
        double result, error;
        double exact = 2/3.;
        gsl_function F;
        F.function = &Fa1;

        gsl_integration_qags (&F, a, b,acc , eps, 1000, w, &result, &error);
        double diff = exact-result;
        printf("∫ dx √(x)) from 0 to 1. With GSL \n");
        printf("------------------------------------\n");
        printf("  Value (exact = 2/3)       : %.18f\n"  , result);
        printf("  Estimated error           :% .18f\n", error);
        printf("-------------------------------\n");

        gsl_integration_workspace_free (w);
    }
    {
        printf("=================================================================\n");
        gsl_integration_workspace * w
        = gsl_integration_workspace_alloc (1000);

        double a = 0.;
        double b = 1.;
        double acc = 0.0001;
        double eps = 0.0001;
        double result, error;
        double exact = M_PI;
        gsl_function F;
        F.function = &Fa2;

        gsl_integration_qags (&F, a, b,acc , eps, 1000, w, &result, &error);
        double diff = exact-result;
        printf("∫ dx 4√(1-x²) from 0 to 1. With GSL \n");
        printf("------------------------------------\n");
        printf("  Value (exact = π)         : %.18f\n"  , result);
        printf("  Difference (value-exact)  : %.18f\n"  , -diff);
        printf("  Estimated error           :% .18f\n", error);
        printf("-------------------------------\n");

        gsl_integration_workspace_free (w);
    }
    {
        printf("=================================================================\n");
        gsl_integration_workspace * w
        = gsl_integration_workspace_alloc (1000);

        double a = 0.;
        double b = 1.;
        double acc = 0.0001;
        double eps = 0.0001;
        double result, error;
        double exact = 2;
        gsl_function F;
        F.function = &Fa3;

        gsl_integration_qags (&F, a, b,acc , eps, 1000, w, &result, &error);
        double diff = exact-result;
        printf("∫ dx 1/√(x) from 0 to 1. With GSL \n");
        printf("------------------------------------\n");
        printf("  Value (exact = 2)         : %.18f\n"  , result);
        printf("  Difference (value-exact)  : %.18f\n"  , -diff);
        printf("  Estimated error           : %.18f\n", error);
        printf("-------------------------------\n");

        gsl_integration_workspace_free (w);
    }
    {
        printf("=================================================================\n");
        gsl_integration_workspace * w
        = gsl_integration_workspace_alloc (1000);

        double a = 0.;
        double b = 1.;
        double acc = 0.0001;
        double eps = 0.0001;
        double result, error;
        double exact = -4;
        gsl_function F;
        F.function = &Fa4;

        gsl_integration_qags (&F, a, b,acc , eps, 1000, w, &result, &error);
        double diff = exact-result;
        printf("∫ dx ln(x)/√(x) from 0 to 1. With GSL \n");
        printf("------------------------------------\n");
        printf("  Value (exact = -4)        : %.18f\n"  , result);
        printf("  Difference (value-exact)  : %.18f\n"  , -diff);
        printf("  Estimated error           :% .18f\n", error);
        printf("-------------------------------\n");

        gsl_integration_workspace_free (w);
    }
return 0;
}
