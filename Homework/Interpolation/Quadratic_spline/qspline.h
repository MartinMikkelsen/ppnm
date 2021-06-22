#ifndef HAVE_QSPLINE_H
#define HAVE_QSPLINE_H
typedef struct {int n; double *x, *y, *b, *c;} qspline;
qspline* qspline_alloc(int n,double* x,double* y);
double qspline_eval(qspline *s, double z);
void qspline_free(qspline *s);
#endif
