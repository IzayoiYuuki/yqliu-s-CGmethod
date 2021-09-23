#include<stdlib.h>

//some function from BLAS
//use gfortran add -lblas when you compile
float saxpy_(int *,float *,float *,int *,float *,int *);
float ssymv_(char *,int *,float *,float *,int *,float *,int *,float *,float *,int *);
float sdot_(int *,float *,int *,float *,int *);
float scopy_(int *,float *,int *,float *,int *);
float snrm2_(int *,float *,int *);
float sscal_(int *,float *,float *,int *);

int cgmethod(float *A,float *x0,float *b,float *x,int N,int maxSTEP){
	// slove Ax=b with CG Algorithm
	// initially slove x0
	// N is the length of x0
	// Algorithm stop when ||r||<eps or maxSTEP was reached
	// and the final solve will be saved in x
	
	// p is the gradient for each loop
	float *r,*Ap,*p;
	//some parm for blas function
	float scaalpha=-1,scabeta=1,scagamma=0;
	int incx=1,stepnow;
	char uplo='U';
	//some parm for CG
	float dotrr,dotrr1,alpha,beta;
	
	r=(float*)malloc(N*sizeof(float));
	p=(float*)malloc(N*sizeof(float));
	Ap=(float*)malloc(N*sizeof(float));
	
	// step No. 0
	scopy_(&N,b,&incx,r,&incx);
	scopy_(&N,x0,&incx,x,&incx);
	// compute r0=b-Ax0 p0=r0 dotrr=(r0,r0) 
	ssymv_(&uplo,&N,&scaalpha,A,&N,x,&incx,&scabeta,r,&incx);
	scopy_(&N,r,&incx,p,&incx);
	dotrr=snrm2_(&N,r,&incx);
	dotrr=dotrr*dotrr;	

	//loop
	for(stepnow=0;stepnow<maxSTEP;stepnow++){
		//compute A*p
		ssymv_(&uplo,&N,&scabeta,A,&N,p,&incx,&scagamma,Ap,&incx);
		//compute alpha=(rj,rj)/(Ap,p)
		alpha=dotrr/sdot_(&N,Ap,&incx,p,&incx);
		//compute x for each step 
		saxpy_(&N,&alpha,p,&incx,x,&incx);
		//compute r=r-alpha*A*p
		alpha=-alpha;
		saxpy_(&N,&alpha,Ap,&incx,r,&incx);
		alpha=-alpha;
		//compute dotrr1=(rj+1,rj+1)
		dotrr1=snrm2_(&N,r,&incx);
		//if ||r||<eps stop
		if(dotrr1<1./1000000){
			break;}
		dotrr1=dotrr1*dotrr1;
		//compute beta=(rj+1,rj+1)/(rj,rj)
		beta=dotrr1/dotrr;
		dotrr=dotrr1;
		//compute p=r+beta*A*p
		scaalpha=1/beta;
		saxpy_(&N,&scaalpha,r,&incx,p,&incx);
		sscal_(&N,&beta,p,&incx);
	}
	return 0;}
		



