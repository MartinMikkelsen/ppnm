#include<stdio.h>
#include<math.h>

double square(double x){return x*x;}

void print_table(double (*f)(double),double a, double b, double dx){
	for(double x=a;x<=b;x+=dx)
		printf("%10g %10g \n",x,(*f)(x));
}

int main() {
	double (*g)(double);
	g=&square;
	g=&cos;
	printf("g(2)=%g \n",(*g)(2));
	print_table(sin,0,M_PI,M_PI/8);
return 0;
}
