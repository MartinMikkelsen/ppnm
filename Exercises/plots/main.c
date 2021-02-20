#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_erf.h>

double Erf(double);
double tgamma(double);

int main(){
	double xmin=-2,xmax=2;
	for(double x=xmin;x<=xmax;x+=1.0/8){
		printf("%10g %10g %10g %10g\n",x,erf(x),gsl_sf_erf(x),Erf(x));
		}

	FILE* data2=fopen("data2.txt","w");
	for(double x=xmin;x<=xmax;x+=1.0/8){
	fprintf(data2,"%10g %10g \n",x,tgamma(x));
		}
return 0;
}
