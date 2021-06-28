#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <math.h>

int equal(double a, double b, double epsilon, double tau){
  if (fabs(a-b)<tau || fabs (a-b)/ fabs (a+b)<epsilon/2){
    return 1;
  }
  else{
    return 0;
  }
}

void name_digit(int i){
  switch (i){
  case 1 :
    printf("One \n");
    break;
  case 2 :
    printf("Two \n");
    break;
  case 3 :
    printf("Three \n");
    break;
  case 4 :
    printf("Four \n");
    break;
  case 5 :
    printf("Five \n");
    break;
  case 6 :
    printf("Six \n");
    break;
  case 7 :
    printf("Seven \n");
    break;
  case 8 :
    printf("Eight \n");
    break;
  case 9 :
    printf("Nine \n");
    break;
  case 0 :
    printf("Zero \n");
  default :
    printf("Not a digit. \n");
  }
}
