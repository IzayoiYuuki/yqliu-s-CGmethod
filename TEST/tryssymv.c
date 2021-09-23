#include<stdio.h>

float ssymv_(char *,int *,float *,float *,int *,float *,int *,float *,float *,int *);

int main(){
	float A[4]={1,3,2,4},x[2]={5,6},y[2];
	float alpha=1,beta=0;
	int incx=1,lda=2,M=2,N=2,i;
	char trans='N',uplo='L';

	ssymv_(&uplo,&N,&alpha,A,&lda,x,&incx,&beta,y,    &incx);
	
	for(i=0;i<M;i++){
		printf("%f ",y[i]);}

	return 1;}

