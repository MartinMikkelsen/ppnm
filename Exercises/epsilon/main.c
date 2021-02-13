#include<limits.h>
#include<float.h>
#include<stdio.h>

/* INT_MAX */

int equal(double a, double b, double tau, double epsilon);

int main(){
  int i=1; while(i+1>i) {i++;}
  printf("max int while = %i\n",i);
  printf("The maximum value of INT_MAX = %d\n",INT_MAX);
  int j=1;
  for(j=1; j+1>j; j++){}
  printf("max int for = %i\n",j);
  int k=1;
	do
	{
		k++;
	}while (k+1>k);
  printf("max int for do = %i\n", k);
  printf("Now for the minimum value\n");
  while(i-1<i) {i++;}
  printf("min int while = %i\n",i);
  printf("The minimum value of INT_MAX = %d\n",INT_MIN);
  for(j=1; j-1<j; j++){}
  printf("min int for = %i\n",j);
  do
	{
		k++;
	}while (k-1<k);
  printf("min int for do = %i\n", k);
  printf("Epsilon machine\n");
  double x1=1; while(1+x1!=1){x1/=2;} x1*=2;
  printf("The number is: %g\n",x1);
  double z1=1;
  for(z1=1; 1+z1!=1; z1/=2){}
  z1*=2;
  printf("The number is: %g\n",z1);
  double y1=1;
  do
	{
    y1/=2;
	}while (1+y1!=1); y1*=2;
  printf("The number is: %g\n",y1);

  float x2=1; while(1+x2!=1){x2/=2;} x2*=2;
  printf("The number is: %g\n",x2);
  float z2=1;
  for(z2=1; 1+z2!=1; z2/=2){}
  z2*=2;
  printf("The number is: %g\n",z2);
  float y2=1;
  do
  {
    y2/=2;
  }while (1+y2!=1); y2*=2;
  printf("The number is: %g\n",y2);

  long double x3=1; while(1+x3!=1){x3/=2;} x3*=2;
  printf("The number is: %Lg\n",x3);
  long double z3=1;
  for(z3=1; 1+z3!=1; z3/=2){}
  z3*=2;
  printf("The number is: %Lg\n",z3);
  long double y3=1;
  do
  {
    y3/=2;
  }while (1+y3!=1); y3*=2;
  printf("The number is: %Lg\n",y3);
  printf("The value for FLT_EPSILON = %g\n", FLT_EPSILON);
  printf("The value for DBL_EPSILON= %g\n", DBL_EPSILON);
  printf("The value for LDBL_EPSILON= %Lg\n", LDBL_EPSILON);

  float f = INT_MAX/300;
  float sum_up_float = 0;
  for(float i=1; i<f; i++){
  sum_up_float += 1.0*f/(i);
  }
  printf("sum up = %g\n",sum_up_float);

  float sum_down_float = 0;
  for(float i=1; i<f; i++){
  sum_down_float += 1.0*f/(f-i);
  }
  printf("sum down = %g\n",sum_down_float);

/*The difference is summing from 1 and up or from up to 1. The sums converges but takes more time for larger f. */

double g = INT_MAX/300;
double sum_up_double = 0;
for(double i=1; i<g; i++){
sum_up_double += 1.0*g/(i);
}
printf("sum up double = %g\n",sum_up_double);

double sum_down_double = 0;
for(double i=1; i<g; i++){
sum_down_double += 1.0*f/(g-i);
}
printf("sum down double = %g\n",sum_down_double);
/*Float up float yields a differenet answer for our sum. */

equal(2,2,2,2);
printf("%i\n",equal(2,2,2,2));
return 0;
}
