#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_erf.h>
#include<gsl/gsl_sf_gamma.h>

double Erf(double);
double Gamma(double);

int main(){
        double xmin=-2,xmax=2;
        for(double x=xmin;x<=xmax;x+=1.0/8){
                printf("%10g %10g %10g %10g\n",x,erf(x),gsl_sf_erf(x),Erf(x));
                }
        double ymin=-5.05,ymax=5;
        FILE* data2=fopen("data2.txt","w");
        for(double y=ymin;y<=ymax;y+=1.0/8){
                fprintf(data2,"%10g %10g %10g %10f \n",y,tgamma(y),gsl_sf_gamma(y),Gamma(y));
                }

return 0;
}
