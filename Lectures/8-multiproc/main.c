#include <stdio.h>
#include <pthread.h>
#include <math.h>

void* bar(void* arg){
	double* x = (double*)arg;
	for(int i=0;i<1e8;i++) *x=cos(*x);
return NULL;
}

int main(){
	double x=0,y=100,z=-100;
	bar((void*)&x);
	bar((void*)&y);
	bar((void*)&z);
	printf("x=%g cos(x)=%g\n",x,cos(x));
	printf("y=%g cos(y)=%g\n",y,cos(y));
	printf("z=%g cos(z)=%g\n",z,cos(z));
return 0;
}
