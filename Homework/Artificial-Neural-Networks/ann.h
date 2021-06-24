#include"gsl/gsl_vector.h"
#ifndef HAVE_ANN_H
#define HAVE_ANN_H
typedef struct { int n; double(*f)(double); gsl_vector* params; } ann;
ann*   ann_alloc   (int n,double(*f)(double));
void   ann_free    (ann* network);
double ann_response(ann* network,double x);
void   ann_train   (ann* network,gsl_vector* xs,gsl_vector* ys);
#endif
