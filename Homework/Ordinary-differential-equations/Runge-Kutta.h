
void rkstep12(void f(double x,gsl_vector* y, gsl_vector* dydx), double x, gsl_vector* yx, double h, gsl_vector* yh, gsl_vector* err);

void rkstep23(void (*f)(int n, double x, double* y, double* dydx), int n, double x, double* y_curr, double h, double* y_next, double* err);

int driver(
	void (*f)(int n, double x, double* y, double* dydx), // right-hand-side of dy/dx = f(x, y) 
    int n,          // size of vectors 
	double  a,      // the start-point a 
	double  b,      // the end-point of the integration 
	double* ya,     // y(a) 
	double* yb,     // y(b) to be calculated 
	double  h,      // initial step-size 
	double  acc,    // absolute accuracy goal 
	double  eps,    // relative accuracy goal 
    char*   outfile // trajectory file
    );



