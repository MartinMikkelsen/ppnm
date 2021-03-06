#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double ex(double x){
	if(x<0)return 1/ex(-x);
	if(x>1./8)return pow(ex(x/2),2);
	return
1+x*(1+x/2*(1+x/3*(1+x/4*(1+x/5*(1+x/6*(1+x/7*(1+x/8*(1+x/9*(1+x/10)))))))));
}

int main(int argc,char** argv){
	double x=1;
	if(argc>1)x=atof(argv[1]);
	/*printf("x=%g\nex =%.20g\nexp=%.20g\n",x,ex(x),exp(x));
	*/
	for(double x=0;x<10;x+=1./2)
		printf("%10g %10g %10g\n",x,ex(x),exp(x));

return 0;
}
