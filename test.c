#include <stdio.h>
#include <stdlib.h>
#include "include/statNS.h"

int main(int argc, char *argv[]){
	int pf;
	
//=====================================================================
//Example 0: Compute f-mode frequencies given an array of central density

	//define central density sequence in the SI unit kg/m^3
	double RhocSI[]={6.5e17,7e17,7.2e17,8e17,8.2e17,8.4e17,8.6e17,8.8e17,9e17,9.2e17,9.5e17,9.8e17,1e18,1.05e18,1.08e18,1.1e18,1.12e18,1.14e18,1.2e18};

	//allocate array for storing the output
	//you may register a large enough array to store results when using longer input.
	CompactStar_t Results[32];
	
	//allocate EoS type variable to store your EoS density and pressure info in log10() and ascending format.
	EoS_t myEoS;

	//load file to you EoS variable
	loadEoS(&myEoS, "EoS_lib/APR.txt");

	//call f-mode multi-thread computation, it will return the results to Results[],
	//the last two numbers represent: 18 the length of RhocSI[], 4 the total threads for computation.
	fmode_mt (Results, &myEoS, RhocSI, 18, 4);

	//print results to console
	fprintf(stderr,"\nExample 0 output:\n");
	for(pf=0;pf<18;pf++){
		fprintf(stderr,"Rhoc=%.5e, M=%.4lf, R=%.4lf, freq=%.4lf, dmpTime=%.4lf\n",
		Results[pf].Rho,
		Results[pf].M,
		Results[pf].r,
		Results[pf].freq,
		Results[pf].dampTime);
	}

//=====================================================================
//Example 1: Compute f-mode frequencies given an array of mass

	//define mass array sequence in the unit of solar mass
	double massArray[]={1.25, 1.325, 1.4, 1.44 ,1.46, 1.48, 1.5, 1.55, 1.62, 1.64, 1.67, 1.7123, 1.7456, 1.8257, 1.90123, 1.95};

	double rhocArray[32];

	//This function compute the central density that corresponds to the mass array one by one,
	//The central density will be written into rhocArray[],
	//the last two numbers are, 16 the length fo massArray[], 6 the desired threads to use
	M2Rhoc_Arr_fm(rhocArray, &myEoS, massArray, 16, 6);

	//Call the f-mode multi-thread computation with rhocArray[], see also in Example 0.
	fmode_mt(Results, &myEoS, rhocArray, 16, 6);

	//print results to console
	fprintf(stderr,"\nExample 1 output:\n");
	for(pf=0;pf<16;pf++){
		fprintf(stderr,"Rhoc=%.5e, M=%.4lf, R=%.4lf, freq=%.4lf, dmpTime=%.4lf\n",
		Results[pf].Rho,
		Results[pf].M,
		Results[pf].r,
		Results[pf].freq,
		Results[pf].dampTime);
	}

//=====================================================================
//Example 2: Compute TOV equation given an array of central density

	//Similar to the fmode_mt() in Example 0, 18 is the length of RhocSI[]. 
	//There will be in total 6 parallel worker threads.
	solveTOV_mt(Results, &myEoS, RhocSI, 18, 6);

	//print results to console
	fprintf(stderr,"\nExample 2 output:\n");
	for(pf=0;pf<18;pf++){
		fprintf(stderr,"Rhoc=%.5e, M=%.4lf, R=%.4lf, I=%.4lf, Lambda=%.4lf\n",
		Results[pf].Rho,
		Results[pf].M,
		Results[pf].r,
		Results[pf].I,
		Results[pf].Lambda);
	}

//=====================================================================
//Example 3: Compute TOV equation given an array of mass

	//Similar usage as M2Rhoc_Arr_fm() in Example 1,
	//NOTICE: 
	//integrator for fmode-xx related computation is different from that for solveTOV-xx
	//So I implement the central density determining function M2Rhoc() in two.
	M2Rhoc_Arr_s(rhocArray, &myEoS, massArray, 16, 4);
	
	//Call the multi-thread TOV solver as Example 2
	solveTOV_mt(Results, &myEoS, rhocArray, 16, 4);
	
	//print results to console
	fprintf(stderr,"\nExample 3 output:\n");
	for(pf=0;pf<16;pf++){
		fprintf(stderr,"Rhoc=%.5e, M=%.4lf, R=%.4lf, I=%.4lf, Lambda=%.4lf\n",
		Results[pf].Rho,
		Results[pf].M,
		Results[pf].r,
		Results[pf].I,
		Results[pf].Lambda);
	}

//=====================================================================
//Example 4: Compute TOV equation given an array of central density and EoS parameters

	//We integrate the amEoS (Nai-Bo Zhang and Bao-An Li)
	//As an example here we choose XL=58.7 MeV, Ksym=-105 MeV, Jsym=400 MeV, J0=115 MeV
	if(genAmEoS(&myEoS, 58.7, -105, 400, 115, NULL)<0){
		fprintf(stderr,"AmEoS failed!\n");
		return -1;
	}
	solveTOV_mt(Results, &myEoS, RhocSI, 16, 4);
	//print results to console
	fprintf(stderr,"\nExample 4 output:\n");
	for(pf=0;pf<16;pf++){
		fprintf(stderr,"Rhoc=%.5e, M=%.4lf, R=%.4lf, I=%.4lf, Lambda=%.4lf\n",
		Results[pf].Rho,
		Results[pf].M,
		Results[pf].r,
		Results[pf].I,
		Results[pf].Lambda);
	}

//=====================================================================
//Example 5: Use EoS_opt_t to alter the default density grid of amEoS

	//Firstly, you need to clearify an EoS_opt_t object
	EoS_opt_t myeosopt;
	set_EoS_default_opt(&myeosopt);
	myeosopt.core_length=250;	/*set the core segment's length to 250*/
	myeosopt.maxU=8.0;			/*maximum hadronic density set to 8.0 saturation density*/
	myeosopt.dU=-1;				/*disable the default dU (hadronic density's stepsize)*/

	if(genAmEoS(&myEoS, 58.7, -105, 400, 115, &myeosopt)<0){
		fprintf(stderr,"AmEoS failed!\n");
		return -1;
	}
	solveTOV_mt(Results, &myEoS, RhocSI, 16, 4);
	//print results to console
	fprintf(stderr,"\nExample 5 output:\n");
	for(pf=0;pf<16;pf++){
		fprintf(stderr,"Rhoc=%.5e, M=%.4lf, R=%.4lf, I=%.4lf, Lambda=%.4lf\n",
		Results[pf].Rho,
		Results[pf].M,
		Results[pf].r,
		Results[pf].I,
		Results[pf].Lambda);
	}

	//save the generated EoS to file
	saveEoS(&myEoS, "EoS_lib/test_EoSgen.txt");
	fprintf(stderr,"\nTest EoS saved to EoS_lib/test_EoSgen.txt\n\n");

	return 0;
}
