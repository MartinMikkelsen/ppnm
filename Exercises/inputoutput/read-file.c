#include<stdio.h>
#include<math.h>

int main() {
	double x;
	int items;
	FILE* fp=fopen("input.txt","r");
	FILE* outputfile=fopen("out.txt","w+");
	do{
		items=fscanf(fp,"%lg",&x);
		fprintf(outputfile,"x=%g,cos(x)=%g\n",x,cos(x));
	}while(items!=EOF);
fclose(fp);
fclose(outputfile);
return 0;
}
