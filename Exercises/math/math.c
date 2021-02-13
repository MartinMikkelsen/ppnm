#include<math.h>
#include<complex.h>
#include<stdio.h>

int main(){
  printf("tgamma(5)= %g\n",tgamma(5));
  printf("j1(0.5)=%g \n",j1(0.05));
  double complex z=csqrt(-2);
  printf("sqrt(-2)=%g + %g*I\n", creal(z), cimag(z));
  printf("exp(I*PI)=%g + %g*I\n", creal(cexp(I*M_PI)), cimag(cexp(I*M_PI)));
  printf("exp(I)=%g + %g*I\n", creal(cexp(I)), cimag(cexp(I)));
  printf("cpow(I,M_E)=%g + %g*I\n", creal(cpow(I,M_E)), cimag(cpow(I,M_E)));
  printf("cpow(I,I)=%g+%g*I\n", creal(cpow(I,I)), cimag(cpow(I,I)));
  float x_float = 1.f/9;
  double x_double = 1./9;
  long double x_long_double = 1.L/9;
  printf("x_float=%.25g\n", x_float);
  printf("x_double=%.25lg\n", x_double);
  printf("x_long_double=%.25Lg\n", x_long_double);
return 0;
}
