#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

double f (double x, void * params){
  double z = *(double *) params;
  double f = log(z*x) / sqrt(x);
  return f;
}

double erfg (double x, void * param){

  double erfg = (2/sqrt(M_PI))*exp(-pow(x,2));
  return erfg;
}
double myG (double z){
  gsl_function G;
  G.function = &erfg;

  int limit=999;

  gsl_integration_workspace* w = gsl_integration_workspace_alloc (limit);
  double a=0,result,error;

  gsl_integration_qags(&G,a,z,0,1e-8,limit,w,&result,&error);

  gsl_integration_workspace_free (w);
  return result;
}
double bessel(double x, void * param){
  double t = *(double*)param;
  double bessel=1.0/M_PI *cos(x-t*sin(x));
  return bessel;
}
double mybessel(){
  int limit = 1000;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (limit);
  FILE* DataBessel=fopen("DataBessel.txt","w");
  for(double t=0;t<=20;t+=1.0/100){
    gsl_function H;
    H.function = &bessel;
    H.params = &t;

    double result,error;

    gsl_integration_qags (&H, 0, M_PI, 1e-8, 1e-8, limit, w, &result, &error);
    fprintf(DataBessel,"%10g %10g\n", t, result);
  }
  gsl_integration_workspace_free(w);
  fclose(DataBessel);
  return 0;
}

int main (){
  int limit = 1000;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (limit);

  double result, error;
  double z=1;

  gsl_function F;
  F.function=&f;
  F.params = &z;

  gsl_integration_qags (&F, 0, 1, 0, 1e-8, limit,w, &result, &error);

  printf("Result= %g\n", result);

  FILE * datab=fopen("DataB.txt","w");
  for (double x= -3; x<=3;x+=1.0/100){
    fprintf(datab,"%10g %10g\n",x,myG(x));
  }

  mybessel();

  gsl_integration_workspace_free (w);
  fclose(datab);

  return 0;
}
