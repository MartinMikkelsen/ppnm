#include <stdio.h>
#include <limits.h>
#include <float.h>

int E1(){
  //Exercise i
  printf("Exercise 1.i:\n\n");
  int i,j,k=1;
  while(i+1>i){
     i++;
  }
  printf("My value of i for the \"while\" loop is %i:\n",i);
  for (j=1; j+1>j; j=j+1){
  }
  printf("My value of i for the \"for\" loop is %i:\n",j);
  do {
    k=k+1;
  } while (k+1>k);
  printf("My value of i for the \"do-while\" loop is %i:\n",k);

  printf("The value of INT_MAX in limits.h is %i:\n", INT_MAX);

  i=j=k=1;

  //Exercise ii
  printf("\n\nExercise 1.ii\n\n");

  while (i-1<i){
    i--;
  }
  printf("My value of i for the \"while\" loop is %i:\n",i);
  for (j=1; j-1<j; j=j-1){
  }
  printf("My value of i for the \"for\" loop is %i:\n",j);
  do {
    k=k-1;
  } while (k-1<k);
  printf("My value of i for the \"do-while\" loop is %i:\n",k);

  printf("The value of INT_MAX in limits.h is %i:\n", INT_MIN);

  //Exercise iii
  printf("\n\nExercise 1.iii\n\n");

  double x = 1;

  while (1+x!=1){
    x/=2;
  }
  x*=2;
  
  printf("The double x= %g\n",x);
  float y = 1;

  while (1+y!=1){
    y/=2;
  }
  y*=2;

  printf("The float x= %g\n",y);

  long double z = 1;

  while (1+z!=1){
    z/=2;
  }
  z*=2;

  printf("The long double x= %Lg\n",z);

  
  printf("FLT_EPSILON: %g \nDBL_EPSILON: %g \nLDBL_EPSILON: %Lg \n",FLT_EPSILON, DBL_EPSILON, LDBL_EPSILON);


  for (x=1;1+x!=1; x/=2){
  }
  x*=2;

  printf("for loop for double e= %g\n",x);


  for (y=1;1+y!=1; y/=2){
  }
  y*=2;

  printf("for loop for float e= %g\n",y);


  for (z=1;1+z!=1; z/=2){
  }
  z*=2;

  printf("for loop for long double e= %Lg\n",z);

  x = 1;

  do {
    x/=2;
  } while (1+x!=1);
  x*=2;
  printf("do while loop for double ex = %g\n",x);

  y =1;

  do {
    y/=2;
  } while (1+y!=1);
  y*=2;
  printf("do while loop for float ex = %g\n",y);

  z = 1;

  do {
    z/=2;
  } while (1+z!=1);
  z*=2;
  printf("do while loop for long double ex = %Lg\n",z);

  return 0;
}

int E2(){
  int max = INT_MAX/3e4;
  printf("max=%i \n", max);
  int i=1;
  float sum_up_float, sum_down_float=0.0;
  while(i<max){
    sum_up_float = 1.f/i+sum_up_float;
    i++;
  }
  
  printf("The sum_up_float becomes: %f \n",sum_up_float);

  while(i>0){
    sum_down_float = sum_down_float+1.f/i;
    i--;
  }
  printf("The sum_down_float becomes: %f \n", sum_down_float);


  double sum_up_double, sum_down_double=0.0;

  i=1;

  while(i<max){
    sum_up_double=sum_up_double+1.f/i;
    i++;
  }
  while(i>0){
    sum_down_double=sum_down_double+1.f/i;
    i--;
  }
  printf("The sum_up_double is: %f \nThe sum_down_double is: %f \n",sum_up_double, sum_down_double);

  sum_up_float=0;
  sum_down_float=0;
  

  for (i=1;i<max;i++){
    sum_up_float=1.0f/i+sum_up_float;
  }
  for (i=max;i>0;i--){
    sum_down_float=1.0f/i+sum_down_float;
  }
  printf("The sum_up_float for the for loop is: %f \nThe sum_down_float for the for loop is: %f \n",sum_up_float, sum_down_float);

  sum_up_double=0;
  sum_down_double=0;

  for (i=1;i<max;i++){
    sum_up_double=1.0f/i+sum_up_double;
  }
  for (i=max;i>0;i--){
    sum_down_double=1.0f/i+sum_down_double;
  }
  printf("The sum_up_double for the for loop is: %f \nThe sum_down_double for the for loop is: %f \n",sum_up_double, sum_down_double);
  i=1;
  
  sum_up_double=0;
  sum_down_double=0;

  do {
    sum_up_double=sum_up_double+1.0f/i;
    i++;
  } while(i<max);

  do {
    sum_down_double=sum_down_double+1.0f/i;
    i--;
  } while(i>0);

  printf("The sum_up_double for the do while loop is: %f \nThe sum_down_double for the do while loop is: %f \n",sum_up_double, sum_down_double);

  sum_up_float=0;
  sum_down_float=0;

  i=1;

  do {
    sum_up_float=sum_up_float+1.0f/i;
    i++;
  } while(i<max);

  do {
    sum_down_float=sum_down_float+1.0f/i;
    i--;
  } while(i>0);

  printf("The sum_up_float for the do while loop is: %f \nThe sum_down_float for the do while loop is: %f \n",sum_up_float, sum_down_float);
  return 0;
}

int equal(double a, double b, double epsilon, double tau);

void name_digit(int i);

int main(){
  printf("\n\n\n################## Exercise 1 #################\n");
  E1();
  printf("\n\n\n################## Exercise 2 #################\n");
  E2();
  printf("\n\n\n################## Exercise 3 #################\n");
  printf("The return value of the function is: %i\n", equal(2.3,2.5,0,0));
  printf("\n\n\n################## Exercise 4 #################\n");
  name_digit(12);
  return 0;
}
