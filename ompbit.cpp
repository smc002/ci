//guji size of sparse
//step = thread_number -done. inefficient at line 170: if (j%THREAD_NUMBER != tid)
//crash on sorted -done.
//put in the para -done.
//
//eigs
#define __STDC_FORMAT_MACROS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mat.h"
#include <inttypes.h> 
#include <math.h>
#include <omp.h>

#define CONFIG_SIZE 42378 
//#define CONFIG_SIZE 40000
#define SIZEJ 20
#define M 10
#define MIN_VALUE 1e-8
#define MAX_UINT 0xffffffffffffffff
#define INDJ(i,j,k,l) (J[(i)+(j)*SIZEJ+(k)*SIZEJ*SIZEJ+(l)*SIZEJ*SIZEJ*SIZEJ])

#define THREAD_NUMBER 4


int main()
{
	/* Variables to read mat file. */
	MATFile *pmat;
	const char *name;
	int	  ndir;
	int	  i,j;
	mxArray *pa;

	/* Variables to get config and its dimensions. */
    unsigned char *config;              

	uint64_t tempNumebr;
	uint64_t u64config[CONFIG_SIZE]={0};
	mwSize nzmax=0.18e8;


	/*
 	* Open file to get directory
 	*/
	pmat = matOpen("./config_unsorted.mat", "r");
	if (pmat == NULL) {
  		printf("Error opening file.\n");
  		return(1);
	}
	pa = matGetNextVariable(pmat, &name);
	if ((pa==NULL)||(strcmp(name, "config")!=0)) {
		printf("Error reading in file.\n");
		return(1);
	} 

    config = (unsigned char*)mxGetPr(pa);

	/* convert the config(10*42378) to u64config(42378) */
	for (i=0;i<CONFIG_SIZE;i++)
		for (j=0;j<M;j++)
		{
			tempNumebr = (uint64_t)config[i*M+j];
			asm("btsq %1,%0" : "+r" (u64config[i]) : "g" (tempNumebr));
		}

	/* Read E and J in cidat_1119_forBit.mat */
	double *E;
	double *J;
	pmat = matOpen("./cidat_1119_forBit.mat", "r");
	if (pmat == NULL) {
  		printf("Error opening file.\n");
  		return(1);
	}
	pa = matGetNextVariable(pmat, &name);
	if ((pa==NULL)||(strcmp(name, "E")!=0)) {
		printf("Error reading in file.\n");
		return(1);
	} 
    E = (double*)mxGetPr(pa);

	pa = matGetNextVariable(pmat, &name);
	if ((pa==NULL)||(strcmp(name, "J")!=0)) {
		printf("Error reading in file.\n");
		return(1);
	} 
	matClose(pmat);
    J = (double*)mxGetPr(pa);

	
	printf("para begin!\n");	

	//todo: use #define to do the decalaration
//	mxArray *ResSparse0, *ResSparse1, *ResSparse2, *ResSparse3;
	mxArray *ResSparse[THREAD_NUMBER];
	MATFile *pmatpara;
	for (i=0;i<THREAD_NUMBER;i++)
	{
    ResSparse[i] = mxCreateSparse(CONFIG_SIZE,CONFIG_SIZE,nzmax/THREAD_NUMBER,mxREAL);
    //ResSparse0 = mxCreateSparse(CONFIG_SIZE,CONFIG_SIZE,nzmax/THREAD_NUMBER,mxREAL);
	printf("Try to create sparse! p = %lld\n", (long long int)ResSparse[i]);	
	}
//    ResSparse1 = mxCreateSparse(CONFIG_SIZE,CONFIG_SIZE,nzmax/THREAD_NUMBER,mxREAL);
//	printf("Try to create sparse! p = %lld\n", (long long int)ResSparse1);	
//    ResSparse2 = mxCreateSparse(CONFIG_SIZE,CONFIG_SIZE,nzmax/THREAD_NUMBER,mxREAL);
//	printf("Try to create sparse! p = %lld\n", (long long int)ResSparse2);	
//    ResSparse3 = mxCreateSparse(CONFIG_SIZE,CONFIG_SIZE,nzmax/THREAD_NUMBER,mxREAL);
//	printf("Try to create sparse! p = %lld\n", (long long int)ResSparse3);	
	
	#pragma omp parallel shared(E,J,ResSparse) private(i,j) num_threads(THREAD_NUMBER)
{
	/* Create a sparse matrix using matlab functions */
	int tid = omp_get_thread_num();
    mwIndex *irs,*jcs;
    double *sr;
    sr  = mxGetPr(ResSparse[tid]);
    irs = mxGetIr(ResSparse[tid]);
    jcs = mxGetJc(ResSparse[tid]);

//this is not needed now.
//	int begin, end;
//	begin = tid;//tid 0 start with j=0; tid 1 start with j=1; and so on.
//	switch (tid)
//	{
//		case 0:
//	begin = 0, end = CONFIG_SIZE/THREAD_NUMBER;
//	printf("Calculate the begin/end number. tid=%d. begin=%d, end=%d\n",tid,begin,end);	
//	break;
//
//		case 1:
//    sr  = mxGetPr(ResSparse[tid]);
//    irs = mxGetIr(ResSparse[tid]);
//    jcs = mxGetJc(ResSparse[tid]);
//	begin = CONFIG_SIZE/THREAD_NUMBER, end = CONFIG_SIZE/THREAD_NUMBER*2;
//	printf("Calculate the begin/end number. tid=%d. begin=%d, end=%d\n",tid,begin,end);	
//	break;
//		case 2:
//    sr  = mxGetPr(ResSparse[tid]);
//    irs = mxGetIr(ResSparse[tid]);
//    jcs = mxGetJc(ResSparse[tid]);
//	begin = CONFIG_SIZE/THREAD_NUMBER*2, end = CONFIG_SIZE/THREAD_NUMBER*3;
//	printf("Calculate the begin/end number. tid=%d. begin=%d, end=%d\n",tid,begin,end);	
//	break;
//		case 3:
//    sr  = mxGetPr(ResSparse[tid]);
//    irs = mxGetIr(ResSparse[tid]);
//    jcs = mxGetJc(ResSparse[tid]);
//	begin = CONFIG_SIZE/THREAD_NUMBER*3, end = CONFIG_SIZE;
//	printf("Calculate the begin/end number. tid=%d. begin=%d, end=%d\n",tid,begin,end);	
//	break;
//	}

	/* Real Calculation */
	double pr;
    int k = 0;
	uint64_t DiffNumber;
	int a, b;
	int DiffInd[4];
	//unsigned char* BArray = config + tid*M;//*sizeof(unsigned char);
	unsigned char* BArray;//*sizeof(unsigned char);
	printf("hello from tid %d! Cal begin now\n", tid);	
//		for (a=0; a<M; a++)
//	printf("hello from tid %d! The 1st config is %d\n", tid, (int)BArray[a]);	

	/* cal diag */
//    for (j=tid, BArray = config + tid*M; (j<CONFIG_SIZE); j+=THREAD_NUMBER,BArray+=M*THREAD_NUMBER ) {
    for (j=0,BArray = config; (j<CONFIG_SIZE); j+=1,BArray+=M) {
        jcs[j] = k;
		if (j%THREAD_NUMBER != tid)
			continue;
		
		pr = 0;
		for (a=0; a<M; a++)
			pr += E[BArray[a]];
		for (a=0; a<M; a++)
			for (b=0; b<a; b++)
				pr += INDJ(BArray[a],BArray[b],BArray[a],BArray[b]);
		sr[k] = pr/2;
		irs[k] = j;
		k++;

		/* cal others */
        for (i=j+1; (i<CONFIG_SIZE); i++) {
			DiffNumber = u64config[i] ^ u64config[j];
			a =  __builtin_popcountll (DiffNumber);
			//if(i==1452&&j==107)
				//printf("a = %d\n", a);
			if (a>4)
				continue;
			else if (a==4){
				uint64_t u64i = u64config[i] & DiffNumber;
				DiffInd[0] = __builtin_ctzll(u64i);
				u64i &= (u64i - 1);
				DiffInd[2] = __builtin_ctzll(u64i);

				uint64_t u64j = u64config[j] & DiffNumber;
				DiffInd[1] = __builtin_ctzll(u64j);
				u64j &= (u64j - 1);
				DiffInd[3] = __builtin_ctzll(u64j);
				pr = INDJ(DiffInd[0], DiffInd[2], DiffInd[1], DiffInd[3]);

			//if(i==1754&&j==22)
				//printf("%d,%d,%d,%d,pr=%f\n", DiffInd[0],DiffInd[1],DiffInd[2],DiffInd[3],pr);

				/*pr = __builtin_parityll(
					((MAX_UINT<<DiffInd[0])^(MAX_UINT<<DiffInd[1]))&
					u64config[i]&
					u64config[j]
					)
					? (-pr) : (pr);
				pr = __builtin_parityll(
					((MAX_UINT<<DiffInd[2])^(MAX_UINT<<DiffInd[3]))&
					u64config[i]&
					u64config[j]
					)
					? (-pr) : (pr);*/
				pr =__builtin_parityll(
					(((MAX_UINT<<DiffInd[0])^(MAX_UINT<<DiffInd[1]))^
					((MAX_UINT<<DiffInd[2])^(MAX_UINT<<DiffInd[3])))&
					u64config[i]&
					u64config[j]
					)
					? (-pr) : (pr);
				/*pr =( __builtin_parityll(
					((MAX_UINT<<DiffInd[0])^(MAX_UINT<<DiffInd[1]))&
					u64config[i]&
					u64config[j]
					))
					==
					(__builtin_parityll(
					((MAX_UINT<<DiffInd[2])^(MAX_UINT<<DiffInd[3]))&
					u64config[i]&
					u64config[j]
					))
					? (-pr) : (pr);*/
				/* printf("Diff Con Is 4. i = %d, j= %d; DiffInd is %d, %d, %d and %d. pr is %f\n", i, j, DiffInd[0], DiffInd[1], DiffInd[2], DiffInd[3], pr);*/

				goto NOTZERO;
			}
			else{
				uint64_t u64i = u64config[i] & DiffNumber;
				DiffInd[0] = __builtin_ctzll(u64i);
				uint64_t u64j = u64config[j] & DiffNumber;
				DiffInd[1] = __builtin_ctzll(u64j);

				pr = 0;
				for (a=0; a<SIZEJ; a++)
				{
					pr += (((1<<a)&u64config[i])?(INDJ(DiffInd[0], a, DiffInd[1], a)):0);/*todo. what if a>32*/
						/*printf("wow, %d\n",a);
						printf("wow, %d\n",(1<<a));
						printf("wow, %d\n",((1<<a)&a));*/
				}
				pr = __builtin_parityll(
					((MAX_UINT<<DiffInd[0])^(MAX_UINT<<DiffInd[1]))&
					u64config[i]&
					u64config[j]
					)
					? (-pr) : (pr);
			/*if(i==1452&&j==107)
				printf("%d,%d,pr=%f\n", DiffInd[0],DiffInd[1],pr);*/
			}
		NOTZERO:
			if (fabs(pr) > MIN_VALUE)
			{
				sr[k] = pr;
				irs[k] = i;
				k++;
	//printf("hello from tid %d! k = %d, i = %d, j = %d, pr = %f\n", tid,k,i,j,pr);	
//				if ((j==12012))//&&(i%(1<<4)== 0))
//	printf("hello from tid %d! k = %d, i = %d, j = %d\n", tid,k,i,j);	
			}
		}
	}
    jcs[j] = k;

	printf("tid %d cal end! k = %d\n", tid,k);	
}
//end of OpenMP para cal

	/* Save the sparse matrix. */
	//#pragma omp critical
	char filename[30];
	char varname[30];
	for(i=0;i<THREAD_NUMBER;i++)
	{
	sprintf(filename, "./Res_forBit%d.mat",i);
	sprintf(varname, "res%d",i);
	pmatpara = matOpen(filename, "w");
	printf("put res%d returns %d\n", i, matPutVariable(pmatpara, varname, ResSparse[i]));
	mxDestroyArray(ResSparse[i]);
	matClose(pmatpara);
	printf("File %d is saved!\n",i);	
	}

//	pmatpara = matOpen("./Res_forBit1.mat", "w");
//	printf("put res1 returns %d\n", matPutVariable(pmatpara, "res1", ResSparse[1]));
//	mxDestroyArray(ResSparse[1]);
//	matClose(pmatpara);
//	printf("File 1 is saved!\n");	
//
//	pmatpara = matOpen("./Res_forBit2.mat", "w");
//	printf("put res2 returns %d\n", matPutVariable(pmatpara, "res2", ResSparse[2]));
//	mxDestroyArray(ResSparse[2]);
//	matClose(pmatpara);
//	printf("File 2 is saved!\n");	
//
//	pmatpara = matOpen("./Res_forBit3.mat", "w");
//	printf("put res3 returns %d\n", matPutVariable(pmatpara, "res3", ResSparse[3]));
//	mxDestroyArray(ResSparse[3]);
//	matClose(pmatpara);
//	printf("File 3 is saved!\n");	

//	printf("hello from tid %d! j= %d\n", tid,j);	
//	pmatpara = matOpen("./Res_forBit2.mat", "w");
//	matPutVariable(pmatpara, "res_forBit2", ResSparse0);
//	mxDestroyArray(ResSparse0);
//	matClose(pmatpara);
//
//	printf("hello from tid %d! j= %d\n", tid,j);	
//	pmatpara = matOpen("./Res_forBit3.mat", "w");
//	matPutVariable(pmatpara, "res_forBit3", ResSparse0);
//	mxDestroyArray(ResSparse0);
//	matClose(pmatpara);
}
