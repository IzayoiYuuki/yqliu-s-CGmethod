#include<stdio.h>

int cgmethod(float *,float *,float *,float*,int,int);

int main(){
	float A[4]={1,2,2,1},x[2],b[2]={5,4},x0[2]={0,0};
	int N=2,STEP=8;
	cgmethod(A,x0,b,x,N,STEP);
	
	int i,j;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			printf("%f ",A[i+N*j]);}
		printf("%f ",x[i]);
		printf("%f ",b[i]);
		printf("\n");}
	return 1;}
		
