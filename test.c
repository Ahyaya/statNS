#include <stdio.h>
#include <stdlib.h>
#include "src/statNS.h"

int main(int argc, char *argv[]){
	int pf;
	double cdenArray[]={6.5e17,7e17,7.2e17,8e17,8.2e17,8.4e17,8.6e17,8.8e17,9e17,9.2e17,9.5e17,9.8e17,1e18,1.05e18,1.08e18,1.1e18,1.12e18,1.14e18,1.2e18};

	CompactStar_t Results[18];
	EoS_t myEoS;

	loadEoS(&myEoS, "APR.txt");

	fmode_mt(Results,&myEoS,cdenArray,18, 5);

	fprintf(stdout,"fmode:\n");
	for(pf=0;pf<18;pf++){
		fprintf(stderr,"rho=%e, M=%.2lf, R=%.2lf, freq=%.2lf, dmpTime=%.4lf\n",
		Results[pf].Rho,
		Results[pf].M,
		Results[pf].r,
		Results[pf].freq,
		Results[pf].dampTime);
	}

	solveTOV_mt(Results,&myEoS,cdenArray,18, 7);

	fprintf(stdout,"\nsolveTOV:\n");
	for(pf=0;pf<18;pf++){
		fprintf(stderr,"rho=%e, M=%.2lf, R=%.2lf, I=%.2lf, Lambda=%.4lf\n",
		Results[pf].Rho,
		Results[pf].M,
		Results[pf].r,
		Results[pf].I,
		Results[pf].Lambda);
	}

	CompactStar_t myResult;

	double rhoc=M2Rhoc(&myEoS, 1.440091);

	solveTOV(&myResult, &myEoS, rhoc);

	fprintf(stdout,"\nTry to find M=1.440091:\nrhoc=%e M=%.8lf\n",rhoc,myResult.M);

	return 0;
}
