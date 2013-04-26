
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include "engine.h"
#include "mat.h"
#define  BUFSIZE 256
#define THREAD_NUMBER 4

int main()

{
	#pragma omp parallel num_threads(THREAD_NUMBER)
	{
	int tid = omp_get_thread_num();

	/* open matlab engine */
	Engine *ep;
	char buffer[BUFSIZE+1];
	if (!(ep = engOpen("matlab -nojvm"))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
	}
	buffer[BUFSIZE] = '\0';
	engOutputBuffer(ep, buffer, BUFSIZE);

	/* read the res from .mat file */
	MATFile *pmat;
	mxArray *pa;
	const char *name;
	char filename[30];
	sprintf(filename, "./Res_forBit%d.mat",tid);
	pmat = matOpen(filename, "r");
	if (pmat == NULL) {
  		printf("Error opening file.\n");
	}
	pa = matGetNextVariable(pmat, &name);
	if (pa==NULL) {
		printf("Error reading in file.\n");
	} 
	engPutVariable(ep, "res", pa);
	engEvalString(ep, "nnz(res)");


	printf("%s\n", buffer);
	}//#pragma omp parallel




	/*
	MATFile *pmat;
	mxArray *pa;
	double *E;
	double *J;
	const char *name;
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
	*/

}
