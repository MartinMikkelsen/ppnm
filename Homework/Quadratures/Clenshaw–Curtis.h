#ifndef HAVE_CC_H
#define HAVE_CC_H

double CC24( double f(double),double a, double b,double acc, double eps, double f2, double f3, int nrec);

double clenshaw(
	double f(double),double a,double b,double acc,double eps);
#endif
