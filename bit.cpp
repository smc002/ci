#define __STDC_FORMAT_MACROS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mat.h"
#include <inttypes.h> 
#include <math.h>

#define CONFIG_SIZE 42378 
#define SIZEJ 20
#define M 10
#define MIN_VALUE 1e-8
#define MAX_UINT 0xffffffffffffffff
#define INDJ(i,j,k,l) (J[(i)+(j)*SIZEJ+(k)*SIZEJ*SIZEJ+(l)*SIZEJ*SIZEJ*SIZEJ])

int main()
{
	/* Variables to read mat file. */
	MATFile *pmat;
	const char **dir;
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
    J = (double*)mxGetPr(pa);

	/* Create a sparse matrix using matlab functions */
	mxArray* ResSparse;
    mwIndex *irs,*jcs;
    double *sr;

    ResSparse = mxCreateSparse(CONFIG_SIZE,CONFIG_SIZE,nzmax,mxREAL);
    sr  = mxGetPr(ResSparse);
    irs = mxGetIr(ResSparse);
    jcs = mxGetJc(ResSparse);

	/* Real Calculation */
	double pr;
    int k = 0;
	uint64_t DiffNumber;
	int a, b;
	int DiffInd[4];
	unsigned char* BArray = config;

	/* cal diag */
    for (j=0; (j<CONFIG_SIZE); j++, BArray+=M) {
        jcs[j] = k;
		
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
			if(i==1452&&j==107)
				printf("a = %d\n", a);
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
			}
		}
	}
    jcs[CONFIG_SIZE] = k;
	printf("wow \n" );

	/* Save the sparse matrix. */
	pmat = matOpen("./Res_forBit.mat", "w");
	matPutVariable(pmat, "res_forBit", ResSparse);
	mxDestroyArray(ResSparse);
	matClose(pmat);
	
}
