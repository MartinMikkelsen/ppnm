#include<stdio.h>
int main(){
	double x=9;
	printf("this goes to stdout: x=%g\n",x);
	fprintf(stdout,"this goes to stdout too: x=%g\n",x);
	fprintf(stderr,"this goes to stderr: x=%g\n",x);
	FILE* my_output_stream = fopen("outfile.txt","w");
	fprintf(my_output_stream, "this goes to my stream: x=%g\n",x);
fclose(my_output_stream);
return 0;
}
