#include<stdio.h>
#include<math.h>
int main(){
	double x;
	int items;
	FILE* my_out_stream=fopen("out.file.txt","w");
	do{
		//items=scanf("%lg",&x); // from stdin
		items=fscanf(stdin,"%lg",&x); // from stdin
		printf("x=%g sin(x)=%g\n",x,sin(x)); // to stdout
		fprintf(stderr,"x=%g sin(x)=%g\n",x,sin(x)); // to stderr
		fprintf(my_out_stream,"x=%g sin(x)=%g\n",x,sin(x));
	}while(items!=EOF);
fclose(my_out_stream);
return 0;
}
