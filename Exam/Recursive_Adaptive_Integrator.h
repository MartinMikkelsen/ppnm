#ifndef HAVE_ADAPT_H
#define HAVE_ADAPT_H

double intrun(double f(double), double a, double b, double acc, double eps, double f2, double f3, double nrec);

double integrate(double f(double), double a, double b,double acc,double eps);

#endif
