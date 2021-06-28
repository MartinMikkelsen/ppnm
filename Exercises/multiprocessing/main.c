#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#define N 1000

void* monte_c(void* arg1){
  FILE* data=fopen("data.txt","a");
  double x,y,z;
  int *count= (int*)arg1;

  unsigned int SEED;
  for (int i=0; i<=N; i++){
    x= (double)rand_r(&SEED)/RAND_MAX;
    y= (double)rand_r(&SEED)/RAND_MAX;
    z=pow((x),2)+pow((y),2);
    fprintf(data,"%10g %10g \n",x,y);
    if (z<=1){
      *count=*count+1;
    }
  }
  fclose(data);
  return NULL;
  
}

int main() {
  double pi;
  int count1=0, count2=0;
  pthread_t thread1;
  FILE* data1=fopen("data.txt","w");

  pthread_create(&thread1, NULL, monte_c, (void*) &count1);

  monte_c((void*)&count2);

  void* returnval = NULL;

  pthread_join(thread1,returnval);

  int tot_count=count1+count2;
  
  pi=(double)(tot_count)/(N*2) *4;
  printf("# of interation=%i, estimate of pi=%g\n",N*2,pi);
  FILE* circle=fopen("circle.txt","w");
  double dj=1e-2;
  for(double j=0;j<=M_PI/2;j=j+dj){
    fprintf(circle,"%10g %10g\n",cos(j),sin(j));
  }

  fclose(data1);
  fclose(circle);
  return 0;
}
