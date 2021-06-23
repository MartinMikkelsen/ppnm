#include<math.h>
#include<assert.h>
#include<stdio.h>
#define SQR2 1.41421356237309504880

//Implement a recursive adaptive integrator that estimates the integral of a given function f(x) on a given interval [a,b] with the required absolute, δ, or relative, ε, accuracy goals.

double intrun(double f(double),double f2, double f3, double a, double b, double acc, double eps, double nrec)
{
    assert(nrec <1000000);
    double f1 = f(a+1.*(b-a)/6.);
    double f4 = f(a+5.*(b-a)/6.);
    double Q = (b-a)*(2.*f1+f2+f3+2*f4)/6.;
    double q = (b-a)*(f1+f4+f2+f3)/4.;
    double toll = acc+eps*fabs(Q);
    double err = fabs(Q-q);

    if(err < toll) return Q;
    else
    {
        double Q1 = intrun(f,f1,f2,a,(a+b)/2,acc/sqrt(2),eps,nrec+1);
        double Q2 = intrun(f,f3,f4,(a+b)/2,b,acc/sqrt(2),eps,nrec+1);
        return Q1+Q2;
    }
}
double integrate(double f(double), double a, double b, double acc, double eps)
    {
        double f2 = f(a+2.*(b-a)/6.);
        double f3 = f(a+4.*(b-a)/6.);
        int nrec = 0;
        double value = intrun(f,f2,f3,a,b,acc,eps,nrec);

        return value;
    }
