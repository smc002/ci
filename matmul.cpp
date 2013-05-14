
#include "mex.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#define THREAD_NUMBER 4
#define M 9

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
	double *u, *v;
    mwIndex *irs,*jcs;
    double *sr;
	mwSize N = mxGetN(prhs[0]);

    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
    u = (double*)mxGetPr(prhs[0]);
	//printf("wow\n");

    v = (double*)mxGetPr(plhs[0]);

    sr  = mxGetPr(prhs[1]);
    irs = mxGetIr(prhs[1]);
    jcs = mxGetJc(prhs[1]);
//	for (int i=0;i<6;i++)
//		printf("sr = %f, irs = %d, jcs = %d\n",sr[i], irs[i], jcs[i]);
//	double A[M*M]={0};
	//double vs[N][THREAD_NUMBER]={0};//many v vector
	double *vs = (double*)malloc(sizeof(double)*N*THREAD_NUMBER);//many v vector
	memset(vs,0,N*THREAD_NUMBER);

	#pragma omp parallel shared(u,sr,irs,jcs) num_threads(THREAD_NUMBER)
	{
	int tid = omp_get_thread_num();
	int begini, endi;
	begini = N/THREAD_NUMBER*tid;
	endi   = N/THREAD_NUMBER + begini;
	if (tid == THREAD_NUMBER-1)
		endi = N;
	//printf("tid = %d, begin at %d, end at %d\n",tid, begini, endi);

	for (int i=begini;i<endi;i++)
	{
		for (int j=jcs[i];j<jcs[i+1];j++)
		{
			vs[i+tid*N] += u[irs[j]] * sr[j];
		}
	}
	}//omp ends
	for (int k=0;k<N;k++)
	{
		v[k]=0;
		for (int p=0;p<THREAD_NUMBER;p++)
			v[k] += vs[k+p*N];
	}

}
