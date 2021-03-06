#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>

#include "../include/statNS.h"

static const double pi=3.14159265358979, 
			  G=6.673e-11,
			  c=2.99792458e8,
			  Msun=1.989e30,
			  Cri=1e-7,
			  Gc2=7.42471382405e-28,
			  Gc4=8.26110825251e-45,
			  Mscale=2.033931261665867e5;


ComputeStatus_t ComputeStatus={0,0,0,
RungeKutta_RK5L_roll,RungeKutta_RK5L_join,
RungeKutta_RK5L_roll_s,RungeKutta_RK5L_join_s};

double interp_p2rho(EoS_t *EoS, double cp, double *VsOUT, int *ref) {
	if(cp<EoS->Pmin) {
		return -1;
	}
	double logp=log10(cp),logrho,crho;
	int n=EoS->length-1;
	n=((*ref)+8)<n?((*ref)+8):n;
	for(;n>0;n--)
	{
		if(logp>EoS->lgP[n-1]) {
			logrho=(logp-EoS->lgP[n-1])*(EoS->lgRho[n]-EoS->lgRho[n-1])/(EoS->lgP[n]-EoS->lgP[n-1])+EoS->lgRho[n-1];
			crho=pow(10,logrho);
			/*compute local speed of sound*/
			*VsOUT=cp/crho*(EoS->lgP[n]-EoS->lgP[n-1])/(EoS->lgRho[n]-EoS->lgRho[n-1]);
			*ref=n-1;
			return(crho);
		}
	}
	fprintf(stderr,"%s interrupted by interp_p2rho(), reason: small P!\n", EoS->FilePath);
	return 0;
}

double interp_rho2p(EoS_t *EoS, double crho) {
	if(crho<EoS->Rhomin){
		return -1;
	}
	double logp,logrho=log10(crho),cp;
	int n=EoS->length-1;
	for(;n>0;n--)
	{
		if(logrho>EoS->lgRho[n-1]) {
			logp=(logrho-EoS->lgRho[n-1])*(EoS->lgP[n]-EoS->lgP[n-1])/(EoS->lgRho[n]-EoS->lgRho[n-1])+EoS->lgP[n-1];
			cp=pow(10,logp);
			return(cp);
		}
	}
	fprintf(stderr,"%s interrupted by interp_rho2p(), reason: small Rho!\n", EoS->FilePath);
	return 0;
}

double interp_p2rho_SI(EoS_t *EoS, double cp, int *ref) {
    double logp=log10(cp),logrho,crho;
	int n=EoS->length-1;
	n=((*ref)+8)<n?((*ref)+8):n;
    for(;n>0;n--){
        if(logp>EoS->lgP_SI[n-1]){
            logrho=(logp-EoS->lgP_SI[n-1])*(EoS->lgRho_SI[n]-EoS->lgRho_SI[n-1])/(EoS->lgP_SI[n]-EoS->lgP_SI[n-1])+EoS->lgRho_SI[n-1];
            crho=pow(10,logrho);
			*ref=n-1;
            return(crho);
        }
    }
    fprintf(stderr,"%s meet small pressure at interp_p2rho_fm()\n", EoS->FilePath);
    return 0;
}

double interp_rho2p_SI(EoS_t *EoS, double crho) {
    double logp,logrho=log10(crho),cp;
    int n=EoS->length-1;
    for(;n>0;n--){
        if(logrho>EoS->lgRho_SI[n-1]){
            logp=(logrho-EoS->lgRho_SI[n-1])*(EoS->lgP_SI[n]-EoS->lgP_SI[n-1])/(EoS->lgRho_SI[n]-EoS->lgRho_SI[n-1])+EoS->lgP_SI[n-1];
            cp=pow(10,logp);
            return(cp);
        }
    }
    fprintf(stderr,"%s meet small density at interp_rho2p_fm()\n", EoS->FilePath);
    return 0;
}
/*======================================================
*	Array adders are matrix operation function.
*	C language is not good at operating matrix!
*/
int RungeKutta_Array_add (RK_Arr_t *Result, double h, RK_Arr_t *K, RK_Arr_t *X) {
	Result->P=X->P+h*(K->P);
	Result->M=X->M+h*(K->M);
	Result->I=X->I+h*(K->I);
	Result->y=X->y+h*(K->y);
	Result->Ag00=X->Ag00+h*(K->Ag00);
	return 0;
}

int RungeKutta_Array_adds (RK_Arr_t *Result, double *h, RK_Arr_t *K, RK_Arr_t *X, int dim) {
	int pf;
	Result->P=X->P;
	Result->M=X->M;
	Result->I=X->I;
	Result->y=X->y;
	Result->Ag00=X->Ag00;

	for (pf=0;pf<dim;pf++){
		Result->P+=h[pf]*((K+pf)->P);
		Result->M+=h[pf]*((K+pf)->M);
		Result->I+=h[pf]*((K+pf)->I);
		Result->y+=h[pf]*((K+pf)->y);
		Result->Ag00+=h[pf]*((K+pf)->Ag00);
	}
	return 0;
}

/*Core function in statNS*/
int dFunc(EoS_t *EoS, RK_Arr_t *K, double r, RK_Arr_t *X, int *ref) {
	double Vs;
	double Rho=interp_p2rho(EoS, X->P,&Vs, ref);
	if(Rho<0){
		return -1;
	}
	K->M = 4.0*pi*Rho*r*r;
	K->P = (Rho+X->P)*(X->M+4*pi*r*r*r*(X->P))/(2*(X->M)*r-r*r);
	K->I = 8.0/3*pi*r*r*r*r*(Rho+X->P)/sqrt(1-2*(X->M)/r)/sqrt(X->Ag00);
	K->Ag00 = -2*(X->Ag00)*(K->P)/(Rho+X->P);
	K->y=-(X->y)*(X->y)/r-(X->y)/r*(1+4.0*pi*r*r*(X->P-Rho))/(1-2*(X->M)/r)-r*(4*pi*(5*Rho+9*(X->P)+(X->P+Rho)/Vs)/(1-2*(X->M)/r)-6/r/r/(1-2*(X->M)/r)-4/r/r/r/r/(1-2*(X->M)/r)/(1-2*(X->M)/r)*(X->M+4*pi*r*r*r*(X->P))*(X->M+4*pi*r*r*r*(X->P)));
	return 0;
}

int dFunc_s(EoS_t *EoS, RK_Arr_s *K, double r, RK_Arr_s *X, int *ref) {
	double Vs;
	double Rho=interp_p2rho(EoS, X->P,&Vs, ref);
	if(Rho<0){
		return -1;
	}
	K->M = 4.0*pi*Rho*r*r;
	K->P = (Rho+X->P)*(X->M+4*pi*r*r*r*(X->P))/(2*(X->M)*r-r*r);
	return 0;
}

int RungeKutta_RK4_roll(EoS_t *EoS, RK_Arr_t *K, double h, double r, RK_Arr_t *X, int *ref) {
	RK_Arr_t rkVar;
	if(dFunc(EoS,&K[0], r, X, ref)){return -1;}
	RungeKutta_Array_add(&rkVar, h/2, &K[0], X);
	if(dFunc(EoS,&K[1], r+h/2, &rkVar, ref)){return -1;}
	RungeKutta_Array_add(&rkVar, h/2, &K[1], X);
	if(dFunc(EoS,&K[2], r+h/2, &rkVar, ref)){return -1;}
	RungeKutta_Array_add(&rkVar, h, &K[2], X);
	if(dFunc(EoS,&K[3], r+h, &rkVar, ref)){return -1;}
	return 0;
}

int RungeKutta_RK4_join (RK_Arr_t *Result, double h, RK_Arr_t *K) {
	Result->P += h/6*(K[0].P+2*(K[1].P)+2*(K[2].P)+K[3].P);
	Result->M += h/6*(K[0].M+2*(K[1].M)+2*(K[2].M)+K[3].M);
	Result->I += h/6*(K[0].I+2*(K[1].I)+2*(K[2].I)+K[3].I);
	Result->y += h/6*(K[0].y+2*(K[1].y)+2*(K[2].y)+K[3].y);
	Result->Ag00 += h/6*(K[0].Ag00+2*(K[1].Ag00)+2*(K[2].Ag00)+K[3].Ag00);
	return 0;
}

int RungeKutta_RK4M_roll(EoS_t *EoS, RK_Arr_t *K, double h, double r, RK_Arr_t *X, int *ref) {
	RK_Arr_t rkVar;
	double H_1[1]={h/2};
	double H_2[2]={0,h/2};
	double H_2M[2]={-h,2*h};
	double H_3[3]={0,0,h};
	if(dFunc(EoS,&K[0], r, X, ref)){return -1;}
	RungeKutta_Array_adds(&rkVar, H_1, K, X,1);
	if(dFunc(EoS,&K[1], r+h/2, &rkVar, ref)){return -1;}
	RungeKutta_Array_adds(&rkVar, H_2, K, X,2);
	if(dFunc(EoS,&K[2], r+h/2, &rkVar, ref)){return -1;}
	RungeKutta_Array_adds(&rkVar, H_3, K, X,3);
	if(dFunc(EoS,&K[3], r+h, &rkVar, ref)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_2M,K,X,2);
	if(dFunc(EoS,&K[4], r+h, &rkVar, ref)){return -1;}
	return 0;
}

int RungeKutta_RK4M_join (RK_Arr_t *Result, double h, RK_Arr_t *K) {
	RK_Arr_t rkVar;
	double H[4]={h/6,h/3,h/3,h/6};
	RungeKutta_Array_adds(&rkVar, H, K, Result, 4);
	RungeKutta_Array_adds(Result, H, K, &rkVar, 0);/*dim=0 means copy the array*/
	/*ComputeStatus.RkErr=(-2*(K[1].P)+2*(K[2].P)+(K[3].P)-(K[4].P))/(K[0].P+2*(K[1].P)+2*(K[2].P)+K[3].P);
	ComputeStatus.RkErr=ComputeStatus.RkErr>0?ComputeStatus.RkErr:(-ComputeStatus.RkErr);*/
	return 0;
}

int RungeKutta_RK5F_roll(EoS_t *EoS, RK_Arr_t *K, double h, double r, RK_Arr_t *X, int *ref) {
	RK_Arr_t rkVar;
	double H_1[1]={0.25*h};
	double H_2[2]={0.093750*h,0.281250*h};
	double H_3[3]={1932.0*h/2179,-7200.0*h/2179,7296.0*h/2179};
	double H_4[4]={439.0*h/216,-8*h,3680.0*h/513,-845.0*h/4104};
	double H_5[5]={-8.0*h/27,2*h,-3544.0*h/2565,1859.0*h/4104,-11.0*h/40};
	if(dFunc(EoS,&K[0], r, X, ref)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_1,K,X,1);
	if(dFunc(EoS,&K[1], r+0.25*h, &rkVar, ref)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_2,K,X,2);
	if(dFunc(EoS,&K[2], r+0.375*h, &rkVar, ref)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_3,K,X,3);
	if(dFunc(EoS,&K[3], r+12.0*h/13, &rkVar, ref)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_4,K,X,4);
	if(dFunc(EoS,&K[4], r+h, &rkVar, ref)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_5,K,X,5);
	if(dFunc(EoS,&K[5], r+0.5*h, &rkVar, ref)){return -1;}
	return 0;
}

int RungeKutta_RK5F_join (RK_Arr_t *Result, double h, RK_Arr_t *K) {
	RK_Arr_t rkVar;
	double H[6]={16.0/135*h, 0, 6656.0/12825*h, 28561.0/56430*h, - 0.18*h, 2.0/55*h};
	RungeKutta_Array_adds(&rkVar, H, K, Result, 6);
	RungeKutta_Array_adds(Result, H, K, &rkVar, 0);/*dim=0 means copy the array*/
	/*ComputeStatus.RkErr=(0.0027778*(K[0].P) - 0.0299415*(K[2].P)  - 0.0291999*(K[3].P) + 0.0200000*(K[4].P) +  0.0363636*(K[5].P))/(16.0/135*(K[0].P) + 6656.0/12825*(K[2].P) + 28561.0/56430*(K[3].P) - 0.18*(K[4].P) + 2.0/55*(K[5].P) );
	ComputeStatus.RkErr=ComputeStatus.RkErr>0?ComputeStatus.RkErr:(-ComputeStatus.RkErr);*/
	return 0;
}

int RungeKutta_RK5M_roll(EoS_t *EoS, RK_Arr_t *K, double h, double r, RK_Arr_t *X, int *ref){
	RK_Arr_t rkVar;
	double H_1[1]={0.2*h};
	double H_2[2]={0.075000*h,0.225000*h};
	double H_3[3]={0.3*h,-0.9*h,1.2*h};
	double H_4[4]={226.0/729*h,-25.0/27*h,880.0/729*h,55.0/729*h};
	double H_5[5]={-181.0/270*h,2.5*h,-266.0/297*h,-91.0/27*h,189.0/55*h};
	if(dFunc(EoS, &K[0], r, X, ref)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_1,K,X,1);
	if(dFunc(EoS, &K[1], r+0.2*h, &rkVar, ref)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_2,K,X,2);
	if(dFunc(EoS, &K[2], r+0.3*h, &rkVar, ref)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_3,K,X,3);
	if(dFunc(EoS, &K[3], r+0.6*h, &rkVar, ref)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_4,K,X,4);
	if(dFunc(EoS, &K[4], r+2.0*h/3, &rkVar, ref)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_5,K,X,5);
	if(dFunc(EoS, &K[5], r+h, &rkVar, ref)){return -1;}
	return 0;
}

int RungeKutta_RK5M_join(RK_Arr_t *Result, double h, RK_Arr_t *K) {
	RK_Arr_t rkVar;
	double H[6]={31.0/540*h, 0, 190.0/297*h, -145.0/108*h, 351.0/220*h, 0.05*h};
	RungeKutta_Array_adds(&rkVar, H, K, Result, 6);
	RungeKutta_Array_adds(Result, H, K, &rkVar, 0);/*dim=0 means copy the array*/
	/*ComputeStatus.RkErr=(-0.030556*(K[0].P) + 0.158730*(K[2].P)  - 0.763889*(K[3].P) + 0.675000*(K[4].P) - 0.039286*(K[5].P))/(31.0/540*(K[0].P) + 190.0/297*(K[2].P) - 145.0/108*(K[3].P) + 351.0/220*(K[4].P) + 0.05*(K[5].P) );
	ComputeStatus.RkErr=ComputeStatus.RkErr>0?ComputeStatus.RkErr:(-ComputeStatus.RkErr);*/
	return 0;
}

int RungeKutta_RK5L_roll(EoS_t *EoS, RK_Arr_t *K, double h, double r, RK_Arr_t *X, int *ref) {
	RK_Arr_t rkVar;
    double H_1[1]={0.08333333333333333*h};
    double H_2[2]={-0.125*h, 0.375*h};
    double H_3[3]={0.6*h,-0.9*h,0.8*h};
    double H_4[4]={0.4875*h,-0.45*h,0.15*h,0.5625*h};
    double H_5[5]={-1.685714285714286*h,1.885714285714286*h,1.371428571428571*h,-1.714285714285714*h,1.142857142857143*h};
    if(dFunc(EoS, &K[0], r, X, ref)){return -1;}
    RungeKutta_Array_adds(&rkVar,H_1,K,X,1);
    if(dFunc(EoS, &K[1], r+0.083333333333333333*h, &rkVar, ref)){return -1;}
    RungeKutta_Array_adds(&rkVar,H_2,K,X,2);
    if(dFunc(EoS, &K[2], r+0.25*h, &rkVar,ref)){return -1;}
    RungeKutta_Array_adds(&rkVar,H_3,K,X,3);
    if(dFunc(EoS, &K[3], r+0.5*h, &rkVar, ref)){return -1;}
    RungeKutta_Array_adds(&rkVar,H_4,K,X,4);
    if(dFunc(EoS, &K[4], r+0.75*h, &rkVar, ref)){return -1;}
    RungeKutta_Array_adds(&rkVar,H_5,K,X,5);
    if(dFunc(EoS, &K[5], r+h, &rkVar, ref)){return -1;}
    return 0;
}

int RungeKutta_RK5L_join (RK_Arr_t *Result, double h, RK_Arr_t *K) {
	RK_Arr_t rkVar;
    double H[6]={0.077777777777777777778*h,0,0.355555555555555555556*h,0.13333333333333333333*h,0.355555555555555555556*h,0.077777777777777777777778*h};
    RungeKutta_Array_adds(&rkVar, H, K, Result, 6);
    RungeKutta_Array_adds(Result, H, K, &rkVar, 0);/*dim=0 means copy the array*/
    return 0;
}

int loadEoS(EoS_t *EoS, char *Path) {
	FILE *inf;
	int n;

	if((inf=fopen(Path,"r"))==NULL) {
		fprintf(stderr,"cannot open %s\n",Path);
		return -1;
	}
	sprintf(EoS->FilePath,"%s",Path);
	EoS->length=0;
	while(fscanf(inf,"%lf",&EoS->lgRho_SI[EoS->length])==1) {
		fscanf(inf,"%lf%*[^\n]",&EoS->lgP_SI[EoS->length]);
		EoS->length++;
	}
	fclose(inf);
	EoS->RhomaxSI=pow(10,EoS->lgRho_SI[EoS->length-1]);
	EoS->RhominSI=pow(10,EoS->lgRho_SI[0]);

/*Convert to natural unit*/
	for(n=0;n<EoS->length;n++) {
		EoS->lgRho[n]=EoS->lgRho_SI[n]-10.175607290470733;
		EoS->lgP[n]=EoS->lgP_SI[n]-27.129849799910058;
	}
	EoS->Pmax=pow(10,EoS->lgP[EoS->length-1]);
	EoS->Pmin=pow(10,EoS->lgP[0]);
	EoS->Rhomax=pow(10,EoS->lgRho[EoS->length-1]);
	EoS->Rhomin=pow(10,EoS->lgRho[0]);
	return 0;
}

int RungeKutta_Array_adds_s (RK_Arr_s *Result, double *h, RK_Arr_s *K, RK_Arr_s *X, int dim) {
	int pf;
	Result->P=X->P;
	Result->M=X->M;

	for (pf=0;pf<dim;pf++){
		Result->P+=h[pf]*((K+pf)->P);
		Result->M+=h[pf]*((K+pf)->M);
	}
	return 0;
}

int RungeKutta_RK5L_roll_s (EoS_t *EoS, RK_Arr_s *K, double h, double r, RK_Arr_s *X, int *ref) {
	RK_Arr_s rkVar;
    double H_1[1]={0.08333333333333333*h};
    double H_2[2]={-0.125*h, 0.375*h};
    double H_3[3]={0.6*h,-0.9*h,0.8*h};
    double H_4[4]={0.4875*h,-0.45*h,0.15*h,0.5625*h};
    double H_5[5]={-1.685714285714286*h,1.885714285714286*h,1.371428571428571*h,-1.714285714285714*h,1.142857142857143*h};
    if(dFunc_s(EoS, &K[0], r, X, ref)){return -1;}
    RungeKutta_Array_adds_s(&rkVar,H_1,K,X,1);
    if(dFunc_s(EoS, &K[1], r+0.083333333333333333*h, &rkVar, ref)){return -1;}
    RungeKutta_Array_adds_s(&rkVar,H_2,K,X,2);
    if(dFunc_s(EoS, &K[2], r+0.25*h, &rkVar, ref)){return -1;}
    RungeKutta_Array_adds_s(&rkVar,H_3,K,X,3);
    if(dFunc_s(EoS, &K[3], r+0.5*h, &rkVar, ref)){return -1;}
    RungeKutta_Array_adds_s(&rkVar,H_4,K,X,4);
    if(dFunc_s(EoS, &K[4], r+0.75*h, &rkVar, ref)){return -1;}
    RungeKutta_Array_adds_s(&rkVar,H_5,K,X,5);
    if(dFunc_s(EoS, &K[5], r+h, &rkVar, ref)){return -1;}
    return 0;
}

int RungeKutta_RK5L_join_s (RK_Arr_s *Result, double h, RK_Arr_s *K) {
	RK_Arr_s rkVar;
    double H[6]={0.077777777777777777778*h,0,0.355555555555555555556*h,0.13333333333333333333*h,0.355555555555555555556*h,0.077777777777777777777778*h};
    RungeKutta_Array_adds_s(&rkVar, H, K, Result, 6);
    RungeKutta_Array_adds_s(Result, H, K, &rkVar, 0);/*dim=0 means copy the array*/
    return 0;
}



double getM_s(EoS_t *EoS, double RhocSI) {
	
	RK_Arr_s K[8];
	double h;

	double rho=RhocSI*6.6741e-11, r=0.125/c, r_offset;
	double m=4.0/3*r*r*r*rho, p=interp_rho2p(EoS,rho);
	int pf=0, interpRef=EoS->length-1;
	RK_Arr_s X={p,m};

	/*core section with proceeding stepsize*/
	h=0.125/c;
	for(pf=0;pf<64;pf++) {
		if(ComputeStatus.K_roll_s(EoS, K, h, r, &X, &interpRef)) break;
		ComputeStatus.X_join_s(&X,h, K);
		r+=h;
	}
	h=1.0/c;
	for(pf=0;pf<128;pf++) {
		if(ComputeStatus.K_roll_s(EoS, K, h, r, &X, &interpRef)) break;
		ComputeStatus.X_join_s(&X,h, K);
		r+=h;
	}
	h=4.0/c;
	for(pf=0;pf<256;pf++) {
		if(ComputeStatus.K_roll_s(EoS, K, h, r, &X, &interpRef)) break;
		ComputeStatus.X_join_s(&X,h, K);
		r+=h;
	}
	/*end of core section*/

	/*regular computation*/
	h=16.0/c;
	while(r<1e-4) {
		if(ComputeStatus.K_roll_s(EoS, K, h, r, &X, &interpRef)) break;
		ComputeStatus.X_join_s(&X,h, K);
		r+=h;
	}

	/*use a small stepsize to reach the surface*/
	h=2.0/c;
	while(r<1e-4) {
		if(ComputeStatus.K_roll_s(EoS, K, h, r, &X, &interpRef)) break;
		ComputeStatus.X_join_s(&X,h, K);
		r+=h;
	}
	/*interpolate the radius at surface, according to pressure and its derivative*/
	r_offset=K[0].P<0?X.P/(K[0].P):0;
	r_offset=r_offset>-h?r_offset:-h;
	r-=r_offset;

	m=X.M;

	return m*Mscale;
}

int solveTOV(CompactStar_t *Results, EoS_t *EoS, double RhocSI) {
	
	RK_Arr_t K[8];
	double h;

	double rho=RhocSI*6.6741e-11, r=0.125/c, r_offset;
	double m=4.0/3*r*r*r*rho, p=interp_rho2p(EoS,rho), y=2.0, I=0.0, Ag00=1.0;
	int pf=0, interpRef=EoS->length-1;
	RK_Arr_t X={p,m,I,Ag00,y};

	/*core section with proceeding stepsize*/
	h=0.125/c;
	for(pf=0;pf<64;pf++) {
		if(ComputeStatus.K_roll(EoS, K, h, r, &X, &interpRef)) break;
		ComputeStatus.X_join(&X,h, K);
		r+=h;
	}
	h=1.0/c;
	for(pf=0;pf<128;pf++) {
		if(ComputeStatus.K_roll(EoS, K, h, r, &X, &interpRef)) break;
		ComputeStatus.X_join(&X,h, K);
		r+=h;
	}
	h=4.0/c;
	for(pf=0;pf<256;pf++) {
		if(ComputeStatus.K_roll(EoS, K, h, r, &X, &interpRef)) break;
		ComputeStatus.X_join(&X,h, K);
		r+=h;
	}
	/*end of core section*/

	/*regular computation*/
	h=16.0/c;
	while(r<1e-4) {
		if(ComputeStatus.K_roll(EoS, K, h, r, &X, &interpRef)) break;
		ComputeStatus.X_join(&X,h, K);
		r+=h;
	}

	/*use a small stepsize to reach the surface*/
	h=2.0/c;
	while(r<1e-4) {
		if(ComputeStatus.K_roll(EoS, K, h, r, &X, &interpRef)) break;
		ComputeStatus.X_join(&X,h, K);
		r+=h;
	}
	/*interpolate the radius at surface, according to pressure and its derivative*/
	r_offset=K[0].P<0?X.P/(K[0].P):0;
	r_offset=r_offset>-h?r_offset:-h;
	r-=r_offset;

	m=X.M;
	I=X.I;
	y=X.y;
	p=X.P;
	Ag00=X.Ag00;

	double bt=m/r;
	double Bg00=(1-2.0*m/r)/Ag00;
	Results->Rho=RhocSI;
	Results->I=I/(m*r*r*sqrt(Bg00));
	Results->r=r*c*1e-3;
	Results->M=m*Mscale;
	Results->k2=1.6*bt*bt*bt*bt*bt*(1-2*bt)*(1-2*bt)*(2-y+2*bt*(y-1))/(2*bt*(6-3*y+3*bt*(5*y-8))+4*bt*bt*bt*(13-11*y+bt*(3*y-2)+2*bt*bt*(1+y))+3*(1-2*bt)*(1-2*bt)*(2-y+2*bt*(y-1))*log(1-2*bt));
	Results->Lambda=9495*(Results->k2)*pow(Results->r/10,5)/pow(Results->M,5);

	return 0;
}

void *solveTOVGate (void *args) {
	struct solveTOVparams_t *thisParams;
	thisParams=(struct solveTOVparams_t *) args;
	solveTOV(thisParams->Results,thisParams->EoS,thisParams->RhocSI);
	return 0x00;
}

int solveTOV_mt(CompactStar_t *Results, EoS_t *EoS, double *RhocSI, int arrayLen, int threads) {
	threads=threads<arrayLen?threads:arrayLen;
	if(threads<1){
		return 0;
	}
	struct solveTOVparams_t solveTOVparams[threads];
	pthread_t solveTOVthread[threads];
	int pf,section=0;

	while(section+threads<=arrayLen){
		for(pf=0;pf<threads;pf++,section++){
			solveTOVparams[pf].RhocSI=RhocSI[section];
			solveTOVparams[pf].Results=&Results[section];
			solveTOVparams[pf].EoS=EoS;
			while(pthread_create(&solveTOVthread[pf],NULL,solveTOVGate,&solveTOVparams[pf])){usleep(50000);}
		}
		for(pf=0;pf<threads;pf++){
			pthread_join(solveTOVthread[pf],NULL);
		}
	}

	solveTOV_mt(&Results[section],EoS,&RhocSI[section],arrayLen-section,threads);
	return 0;
}

double getMmax(EoS_t *EoS) {
	double ma,mb,mc,md,Ea,Eb,Ec,Ed,dE=5e16;
    double Mmax,Mmin,m0,m1,m2;
    int n, pf;
	double tmpMassM0M1[2], tmpRhocE0E1[2];
    
	tmpRhocE0E1[1]=(EoS->RhomaxSI<4.2e18)?EoS->RhomaxSI:4.2e18;
	tmpRhocE0E1[0]=tmpRhocE0E1[1]-dE;
	getM_s_mt(tmpMassM0M1, EoS, tmpRhocE0E1, 2, 2);

	m1=tmpMassM0M1[1];
	m0=tmpMassM0M1[0];
	if(m0-m1<=0) {
		EoS->Rhoc_MmaxSI=tmpRhocE0E1[1];
		EoS->Mmax=m1;
		return m1;
	}
	for(n=0;n<7;n++) {
		tmpRhocE0E1[1] *= 0.72;
		tmpRhocE0E1[0]=tmpRhocE0E1[1]-dE;
		getM_s_mt(tmpMassM0M1, EoS, tmpRhocE0E1, 2, 2);
		m1=tmpMassM0M1[1];
		m0=tmpMassM0M1[0];
		
		if(m0-m1<=0) break;
	}
	if(m0-m1>0) {
		EoS->Rhoc_MmaxSI=tmpRhocE0E1[0];
		EoS->Mmax=m0;
		return m0;
	}
	Ea=tmpRhocE0E1[1];Eb=tmpRhocE0E1[1]*1.3888888889;
	Ec=Ea+0.382*(Eb-Ea);
	Ed=Eb-Ec+Ea;

	tmpRhocE0E1[0]=Ec;
	tmpRhocE0E1[1]=Ed;
	getM_s_mt(tmpMassM0M1, EoS, tmpRhocE0E1, 2, 2);
	mc=tmpMassM0M1[0];
	md=tmpMassM0M1[1];
	for(n=0;n<15;n++) {
		if(mc-md>0) {
			Eb=Ed;mb=md;
			Ed=Ec;md=mc;
			Ec=Ea+0.382*(Eb-Ea);
			mc=getM_s(EoS,Ec);
		}else{
			Ea=Ec;ma=mc;
			Ec=Ed;mc=md;
			Ed=Ea+0.618*(Eb-Ea);
			md=getM_s(EoS,Ed);
		}
	}
	EoS->Rhoc_MmaxSI=(mc-md>0)?(0.5*(Ea+Ed)):(0.5*(Ec+Eb));
	EoS->Mmax=getM_s(EoS,EoS->Rhoc_MmaxSI);
	return(EoS->Mmax);
}

double M2Rhoc(EoS_t *EoS, double fM) {
	double Mmax,Mmin,E0,E1,E2,m0,m1,m2;
	int n, threads=8;
	double RhocFrame[threads], massFrame[threads], massGuess, RhocGuess;

	ArXiv_t thisArxiv={0};
	setIndex(thisArxiv.index,sizeof(thisArxiv.index));

	Mmax=getMmax(EoS);
	Mmin=getM_s(EoS,5e17);
	if(Mmax-fM<0) return EoS->Rhoc_MmaxSI;
	if(Mmin>fM) {
		fprintf(stderr,"warning: attempt to find a low mass neutron star!\n");
		return 5e17;
	}
	E0=5e17;E2=EoS->Rhoc_MmaxSI;m0=Mmin;m2=Mmax;

	for(n=0;n<threads;n++){
		RhocFrame[n]=5e17+(n+1)*(EoS->Rhoc_MmaxSI-5e17)/(threads+1);
	}
	getM_s_mt(massFrame, EoS, RhocFrame, threads, threads);

	for(n=0;n<threads;n++){
		rec2arxiv(&thisArxiv,RhocFrame[n],massFrame[n]);
	}

	arcSimSort(0, thisArxiv.length-1, thisArxiv.index, thisArxiv.M);

	for(n=0;n<8;n++){
		interp_ArXiv_M2Rhoc_arr(&RhocGuess, &thisArxiv, &fM, 1);
		massGuess=getM_s(EoS, RhocGuess);
		rec2arxiv(&thisArxiv,RhocGuess,massGuess);
		arcSimSort(0, thisArxiv.length-1, thisArxiv.index, thisArxiv.M);
	}

	interp_ArXiv_M2Rhoc_arr(&RhocGuess, &thisArxiv, &fM, 1);

	return(RhocGuess);
}

double gf(EoS_t *EoS, double e){
	double le, p, gamma;
	int i=EoS->length-1;
	e=e/G*c*c;
	le=log10(e);
	for(;i>0;i--){
  		if(le>EoS->lgRho_SI[i-1]){
    		p=EoS->lgP_SI[i-1]*(le-EoS->lgRho_SI[i])/(EoS->lgRho_SI[i-1]-EoS->lgRho_SI[i])+EoS->lgP_SI[i]*(le-EoS->lgRho_SI[i-1])/(EoS->lgRho_SI[i]-EoS->lgRho_SI[i-1]);
    		p=pow(10,p);
    		gamma=(e+p/c/c)/e*(EoS->lgP_SI[i]-EoS->lgP_SI[i-1])/(EoS->lgRho_SI[i]-EoS->lgRho_SI[i-1]);
    		return(gamma);
		}
	}
	return((pow(10,EoS->lgRho_SI[0])+pow(10,EoS->lgP_SI[0])/c/c)/e*(EoS->lgP_SI[i]-EoS->lgP_SI[i-1])/(EoS->lgRho_SI[i]-EoS->lgRho_SI[i-1]));
}

double Ff(double r,double A,double e,double p,double m,double Dp){
	double F;
	F=8*pi/r/r*sqrt(A);
	F=F*(e+3*p+r*Dp-(m/4/pi/r/r/r+p)*(4+A*(m/r-4*pi*r*r*e)));
	return(F);
}

double fp(double r,double p,double e,double m) {
	return(-G*(e+p/c/c)*(m+4*pi*r*r*r*p/c/c)/(r*r-2*r*m*Gc2));
}

double fm(double r,double e) {
	return(4*pi*r*r*e);
}

double Bf(double r,double p,double m,double B) {
	return(2*Gc2/r/r*(m+4*pi*r*r*r*p/c/c)*B/(1-2*Gc2*m/r));
}

double DH1(double r,double m,double A,double p,double e,double H1,double H0,double K,double V) {
	return(-1/r*(2+1+2*m*A/r+4*pi*r*r*A*(p-e))*H1+1/r*A*(H0+K-16*pi*(e+p)*V));
}

double DK(double r,double H0,double H1,double Dv,double K,double e,double p,double A,double W) {
	return(1/r*H0+0.5*2*3.0/r*H1-(3.0/r-0.5*Dv)*K-8*pi*(e+p)*sqrt(A)/r*W);
}

double DW(double r,double W,double A,double gamma,double p,double B,double X,double V,double H0,double K) {
	return(-(2+1)/r*W+r*sqrt(A)*(1/gamma/p/sqrt(B)*X-6.0/r/r*V+0.5*H0+K));
}

double DX(double r,double X,double e,double p,double B,double Dv,double H0,double w,double H1,double K,double V,double A,double F,double W) {
	return(-2/r*X+(e+p)*sqrt(B)*(0.5*(1/r-0.5*Dv)*H0+0.5*(r*w*w/B+3.0/r)*H1+0.5*(1.5*Dv-1/r)*K-3.0*Dv/r/r*V-1/r*(4*pi*(e+p)*sqrt(A)+w*w*sqrt(A)/B-0.5*r*r*F)*W));
}

double H0f(double r,double B,double X,double m,double p,double A,double H1,double K,double w) {
	return (8*pi*r*r*r/sqrt(B)*X-(3*(m+4*pi*r*r*r*p)-w*w*r*r*r/A/B)*H1+(2*r-w*w*r*r*r/B-1/r*A*(m+4*pi*r*r*r*p)*(3*m-r+4*pi*r*r*r*p))*K)/(3*m+2*r+4*pi*r*r*r*p);
}

double Vf(double r,double w,double e,double p,double B,double A,double Dp,double W,double H0,double X) {
	return(1/w/w/(e+p)*B*(1/sqrt(B)*X+1/r*Dp/sqrt(A)*W-0.5*(e+p)*H0));
}

int fmode(CompactStar_t *Results, EoS_t *EoS, double RhocSI, auxSpace_t *auxSpace) {
	double dr=0.5, l=2;
	double r,r0=1,R,RR,drx,rx;
	double ne,p,e,m,mR,A,B=1.0,BR,Bfactor,m1,m2,m3,m4,p1,p2,p3,p4,B1,B2,B3,B4,pressure,I=0,J,DDf,Df=0,f=1;
	double H1,H0,K,W,X,F,V,Dv,gamma,Dp,N;
	double DH11,DK1,DW1,DX1,H01,H02,K1,K2,x,X1,X2,Xp1,Xp2,W1,W2,V1,V2,V01,V02,V0;
	double w,wcheck,o[2],wi,aR,bR,gR,hR,kR,n,Y1,Y2,Z,DZ,DDZ,VZ,Ar1,Ar2,Ai1,Ai2,ar,ai,Ar[2],Ai[2],Br[2],Bi[2];
	int t,q,wpf,rpf,pfEnd,n4,status=0,interpRef=EoS->length-1;

	double *mfile,*pfile,*rhofile,*Bfile,*Bcor,*Wfile,*Wfile1,*Wfile2,*Vfile,*Vfile1,*Vfile2;
	mfile=auxSpace->mfile;
	pfile=auxSpace->pfile;
	rhofile=auxSpace->rhofile;
	Bfile=auxSpace->Bfile;
	Bcor=auxSpace->Bcor;
	Wfile=auxSpace->Wfile;
	Wfile1=auxSpace->Wfile1;
	Wfile2=auxSpace->Wfile2;
	Vfile=auxSpace->Vfile;
	Vfile1=auxSpace->Vfile1;
	Vfile2=auxSpace->Vfile2;
	
	pressure=pow(10,EoS->lgP_SI[0]);
	r=r0;
	ne=RhocSI;
	e=RhocSI;
	p=interp_rho2p_SI(EoS,e);
	m=1.3333333*r*r*pi*e*r;
	wpf=-1;
	for(;p>pressure;r=r+dr){
		if(wpf>49995){
			Results->Rho=RhocSI;
			Results->M=0;
			Results->r=0;
			Results->freq=0;
			Results->dampTime=0;
			status=1;
			goto funcExit;
		}
		wpf++;
		rhofile[wpf]=e*Gc2;   /*-- e,p,m in G=c=1 --*/
		pfile[wpf]=p*Gc4;
		Bfile[wpf]=B;
		mfile[wpf]=m*Gc2;
		A=1/(1-2*m*Gc2/r);
		p1=fp(r,p,e,m);
		m1=fm(r,e);
		if((p+dr*p1/2)>pressure) e=interp_p2rho_SI(EoS,p+dr*p1/2,&interpRef); else break;
		B1=Bf(r,p,m,B);
		p2=fp(r+dr/2,p+dr*p1/2,e,m+dr*m1/2);
		m2=fm(r+dr/2,e);
		if((p+dr*p2/2)>pressure) e=interp_p2rho_SI(EoS,p+dr*p2/2,&interpRef); else break;
		B2=Bf(r+dr/2,p+dr*p1/2,m+dr*m1/2,B+dr*B1/2);
		p3=fp(r+dr/2,p+dr*p2/2,e,m+dr*m2/2);
		m3=fm(r+dr/2,e);
		if((p+dr*p3)>pressure) e=interp_p2rho_SI(EoS,p+dr*p3,&interpRef); else break;
		B3=Bf(r+dr/2,p+dr*p2/2,m+dr*m2/2,B+dr*B2/2);
		p4=fp(r+dr,p+dr*p3,e,m+dr*m3);
		m4=fm(r+dr,e);
		B4=Bf(r+dr,p+dr*p3,m+dr*m3,B+dr*B3);
		J=-4*pi*(e+p/c/c)*Gc2*r*A;
		DDf=-(4/r*Df+J*Df+4/r*J*f);
		I=I-2.0/3*f*J/sqrt(A*B)*r*r*r*dr;
		p=p+dr*(p1+2*p2+2*p3+p4)/6;
		e=interp_p2rho_SI(EoS,p,&interpRef);
		m=m+dr*(m1+2*m2+2*m3+m4)/6;
		B=B+dr*(B1+2*B2+2*B3+B4)/6;
		f=f+Df*dr;
		Df=Df+DDf*dr;
	}
	pfEnd=wpf;
	R=r;
	mR=m*Gc2;
	BR=1-2*Gc2*m/r;
	Bfactor=BR/B;
	gamma=(EoS->lgP_SI[1]-EoS->lgP_SI[0])/(EoS->lgRho_SI[1]-EoS->lgRho_SI[0]);
	N=1/(gamma-1);
	RR=R-(N+1)*(p-dr*(p1+2*p2+2*p3+p4)/6)/(p1+2*p2+2*p3+p4)*6;
	wpf=-1;
	rpf=-1;
	for(n4=0;n4<pfEnd+1;n4++){
		Bcor[n4]=Bfile[n4]*Bfactor;
	}
	I=I/sqrt(Bfactor);
	I=I/(f+2*I/r/r/r)/Gc2;
	I=m*sqrt(m/I)*Gc2;
	o[0]=(-0.0047+0.133*I+0.575*I*I)/mR-0.1e-5;
	o[1]=o[0]+0.2e-5;
	q=1;
	wcheck=0;
	for(t=0;;t++){
		if(t>20){
			Results->Rho=RhocSI;
			Results->M=0;
			Results->r=0;
			Results->freq=0;
			Results->dampTime=0;
			status=2;
			goto funcExit;
		}
		w=t?o[q]:o[t];
		e=rhofile[0];
		p=pfile[0];
		B=Bcor[0];
		W=1.0;
		K=(e+p);
		X=(e+p)*sqrt(B)*((4*pi/3*(e+3*p)-w*w/B/l)*W+0.5*K);
		H1=(2*l*K+16*pi*(e+p)*W)/l/(l+1);
		rpf=-1;
		wpf=-1;
		r=r0;
		while(rpf<pfEnd){
			rpf++;
			p=pfile[rpf];
			e=rhofile[rpf];
			B=Bcor[rpf];
			m=mfile[rpf];
			Dp=-(e+p)*(m+4*pi*r*r*r*p)/r/r/(1-2*m/r);
			Dv=-2*Dp/(e+p);
			A=1/(1-2*m/r);
			gamma=gf(EoS,e);
			H0=H0f(r,B,X,m,p,A,H1,K,w);
			V=Vf(r,w,e,p,B,A,Dp,W,H0,X);
			if(r==r0)V01=V;
			if(fabs((wcheck-w)/w)<Cri){
				wpf++;
				Wfile[wpf]=sqrt(1-2*m/r)*W;
				Vfile[wpf]=V;
			}
			F=Ff(r,A,e,p,m,Dp);
			DH11=DH1(r,m,A,p,e,H1,H0,K,V);
			DK1=DK(r,H0,H1,Dv,K,e,p,A,W);
			DW1=DW(r,W,A,gamma,p,B,X,V,H0,K);
			DX1=DX(r,X,e,p,B,Dv,H0,w,H1,K,V,A,F,W);
			H1=H1+DH11*dr;
			K=K+DK1*dr;
			W=W+DW1*dr;
			X=X+DX1*dr;
			r=r+dr;
		}
		wpf=-1;
		rpf=-1;
		X1=X;Xp1=DX1;K1=K;H01=H0f(r,B,X,m,p,A,H1,K,w);W1=W;V1=V;
		p=pfile[0];
		e=rhofile[0];
		B=Bcor[0];
		W=1.0;K=-(e+p);
		X=(e+p)*sqrt(B)*((4*pi/3*(e+3*p)-w*w/B/l)*W+0.5*K);
		H1=(2*l*K+16*pi*(e+p)*W)/l/(l+1);
		r=r0;
		while(rpf<pfEnd){
			rpf++;
			p=pfile[rpf];
			e=rhofile[rpf];
			B=Bcor[rpf];
			m=mfile[rpf];
			Dp=-(e+p)*(m+4*pi*r*r*r*p)/r/r/(1-2*m/r);
			Dv=-2*Dp/(e+p);
			A=1/(1-2*m/r);
			gamma=gf(EoS,e);
			H0=H0f(r,B,X,m,p,A,H1,K,w);
			V=Vf(r,w,e,p,B,A,Dp,W,H0,X);
			if(r==r0)V02=V;
			if(fabs((wcheck-w)/w)<Cri){
				wpf++;
				Wfile[wpf]=sqrt(1-2*m/r)*W;
				Vfile[wpf]=V;
			}
			F=Ff(r,A,e,p,m,Dp);
			DH11=DH1(r,m,A,p,e,H1,H0,K,V);
			DK1=DK(r,H0,H1,Dv,K,e,p,A,W);
			DW1=DW(r,W,A,gamma,p,B,X,V,H0,K);
			DX1=DX(r,X,e,p,B,Dv,H0,w,H1,K,V,A,F,W);
			H1=H1+DH11*dr;
			K=K+DK1*dr;
			W=W+DW1*dr;
			X=X+DX1*dr;
			r=r+dr;
		}
		wpf=-1;
		rpf=-1;
		X2=X;Xp2=DX1;K2=K;H02=H0f(r,B,X,m,p,A,H1,K,w);W2=W;V2=V;
		x=-(X1-(RR-R)/(N+1)*Xp1)/(X2-(RR-R)/(N+1)*Xp2);
		H0=H01+x*H02; K=K1+x*K2; W=W1+x*W2; V=V1+x*V2;V0=V01+x*V02;
		if(fabs((wcheck-w)/w)<Cri){
			r=r0;
			while(rpf<pfEnd){
				rpf++;
				W1=Wfile1[rpf];
				W2=Wfile2[rpf];
				W=W1+x*W2;
				wpf++;
				Wfile[wpf]=W/(1+x);
				V1=Vfile1[rpf];
				V2=Vfile2[rpf];
				V=V1+x*V2;
				Vfile[wpf]=V/V0;
				r=r+dr;
			}
			break;
		}
		wcheck=w;
		n=0.5*(l-1)*(l+2);
		aR=-(n*R+3*mR)/(w*w*R*R-(n+1)*mR/R);
		bR=(n*R*(R-2*mR)-w*w*R*R*R*R+mR*(R-3*mR));
		bR=bR/(R-2*mR)/(w*w*R*R-(n+1)*mR/R);
		gR=n*(n+1)*R*R+3*n*mR*R+6*mR*mR;
		gR=gR/R/R/(n*R+3*mR);
		hR=-n*R*R+3*n*mR*R+3*mR*mR;
		hR=hR/(R-2*mR)/(n*R+3*mR);
		kR=-R*R/(R-2*mR);
		Y1=K;
		Y2=aR*H0+bR*K;
		Z=(kR*Y1-Y2)/(kR*gR-hR);
		DZ=(gR*Y2-hR*Y1)/(gR*kR-hR);
		if(w<2e-6){
			Results->Rho=RhocSI;
			Results->M=0;
			Results->r=0;
			Results->freq=0;
			Results->dampTime=0;
			status=3;
			goto funcExit;
		}
		for(r=R;r<25.0/w;r=r+dr){
			drx=dr/(1-2*mR/r);
			VZ=(1-2*mR/r)/r/r/r/(n*r+3*mR)/(n*r+3*mR);
			VZ=VZ*(2.0*n*n*(n+1)*r*r*r+6.0*n*n*mR*r*r+18.0*n*mR*mR*r+18*mR*mR*mR);
			DDZ=(VZ-w*w)*Z;
			Z=Z+DZ*drx;
			DZ=DZ+DDZ*drx;
		}
		r=r-dr;
		rx=r+2*mR*log(r/2/mR-1);
		Ar1=2*cos(w*rx)-2*(n+1)/w/r*sin(w*rx)+1/w/w/r/r*(1.5*mR*w*(1+2/n)*sin(w*rx)-n*(n+1)*cos(w*rx));
		Ai1=2*sin(w*rx)+2*(n+1)/w/r*cos(w*rx)-1/w/w/r/r*(1.5*mR*w*(1+2/n)*cos(w*rx)+n*(n+1)*sin(w*rx));
		Ar2=-2*w*sin(w*rx)-2*(n+1)*cos(w*rx)/r+1/w/r/r*(1.5*mR*w*(1+2/n)*cos(w*rx)+n*(n+1)*sin(w*rx))+(1-2*mR/r)*2*(n+1)/w/r/r*sin(w*rx);
		Ai2=2*w*cos(w*rx)-2*(n+1)*sin(w*rx)/r+1/w/r/r*(1.5*mR*w*(1+2/n)*sin(w*rx)-n*(n+1)*cos(w*rx))-(1-2*mR/r)*2*(n+1)/w/r/r*cos(w*rx);
		ar=(Ai2*Z-Ai1*DZ)/(Ar1*Ai2-Ar2*Ai1);
		ai=-(Ar1*DZ-Ar2*Z)/(Ar1*Ai2-Ar2*Ai1);
		if(t==0){
			Ar[t]=ar;
			Ai[t]=ai;
		}else{
			Ar[q]=ar;
			Ai[q]=ai;
			Br[0]=(o[0]*Ar[1]-o[1]*Ar[0])/(o[0]-o[1]);
			Br[1]=(Ar[0]-Ar[1])/(o[0]-o[1]);
			Bi[0]=(o[0]*Ai[1]-o[1]*Ai[0])/(o[0]-o[1]);
			Bi[1]=(Ai[0]-Ai[1])/(o[0]-o[1]);
			w=-(Br[0]*Br[1]+Bi[0]*Bi[1])/(Br[1]*Br[1]+Bi[1]*Bi[1]);
			if (w<=o[0]){
				o[1]=o[0];o[0]=w;Ar[1]=Ar[0];Ai[1]=Ai[0];q=0;
			}else if(w>=o[1]){
				o[0]=o[1];o[1]=w;Ar[0]=Ar[1];Ai[0]=Ai[1];q=1;
			}else if((o[1]-w)>(w-o[0])){
				o[1]=w;q=1;
			}else{
				o[0]=w;q=0;
			}
		}
	}
	wi=(Br[0]*Bi[1]-Bi[0]*Br[1])/(Br[1]*Br[1]+Bi[1]*Bi[1]);
	Results->Rho=RhocSI;
	Results->M=mR/(Msun*Gc2);
	Results->r=RR/1000;
	Results->freq=w*c/(2000*pi);
	Results->dampTime=1/(wi*c);
funcExit:
	return status;
}

void * fmodeGate (void *args){
	struct fmodeParams_t *thisParams;
	thisParams=(struct fmodeParams_t *) args;
	/*start the computation*/
	fmode(thisParams->Results,thisParams->EoS,thisParams->RhocSI,thisParams->auxSpace);
	return 0x00;
}

int fmode_mt(CompactStar_t *Results, EoS_t *EoS, double *RhocSI, int arrayLen, int threads) {
	threads=threads<arrayLen?threads:arrayLen;
	if(threads<1){
		return 0;
	}
	auxSpace_t auxSpace[threads];
	struct fmodeParams_t fmodeParams[threads];
	pthread_t fmodeThread[threads];
	int pf,section=0;
	/*register auxiliary space for each thread*/
	for(pf=0;pf<threads;pf++){
		auxSpace[pf].mfile=(double*)malloc(50000*sizeof(double));
		auxSpace[pf].pfile=(double*)malloc(50000*sizeof(double));
		auxSpace[pf].rhofile=(double*)malloc(50000*sizeof(double));
		auxSpace[pf].Bfile=(double*)malloc(50000*sizeof(double));
		auxSpace[pf].Bcor=(double*)malloc(50000*sizeof(double));
		auxSpace[pf].Wfile=(double*)malloc(50000*sizeof(double));
		auxSpace[pf].Wfile1=(double*)malloc(50000*sizeof(double));
		auxSpace[pf].Wfile2=(double*)malloc(50000*sizeof(double));
		auxSpace[pf].Vfile=(double*)malloc(50000*sizeof(double));
		auxSpace[pf].Vfile1=(double*)malloc(50000*sizeof(double));
		auxSpace[pf].Vfile2=(double*)malloc(50000*sizeof(double));
	}
	
	while(section+threads<=arrayLen){
		for(pf=0;pf<threads;pf++,section++){
			fmodeParams[pf].RhocSI=RhocSI[section];
			fmodeParams[pf].Results=&Results[section];
			fmodeParams[pf].auxSpace=&auxSpace[pf];
			fmodeParams[pf].EoS=EoS;
			while(pthread_create(&fmodeThread[pf],NULL,fmodeGate,&fmodeParams[pf])){usleep(50000);}
		}
		for(pf=0;pf<threads;pf++){
			pthread_join(fmodeThread[pf],NULL);
		}
	}

	/*release the aux space*/
	for(pf=0;pf<threads;pf++){
		free(auxSpace[pf].mfile);
		free(auxSpace[pf].pfile);
		free(auxSpace[pf].rhofile);
		free(auxSpace[pf].Bfile);
		free(auxSpace[pf].Bcor);
		free(auxSpace[pf].Wfile);
		free(auxSpace[pf].Wfile1);
		free(auxSpace[pf].Wfile2);
		free(auxSpace[pf].Vfile);
		free(auxSpace[pf].Vfile1);
		free(auxSpace[pf].Vfile2);
	}
	fmode_mt(&Results[section],EoS,&RhocSI[section],arrayLen-section,threads);
	return 0;
}


/*Function for fm*/

double getM_fm(EoS_t *EoS, double RhocSI){
    double dr=0.5;
    double r,p,e,m,m1,m2,m3,m4,p1,p2,p3,p4,p_surf;
    int pf, interpRef=EoS->length;

    p_surf=pow(10,EoS->lgP_SI[0]);
    r=1.0;
    e=RhocSI;
    p=interp_rho2p_SI(EoS,e);
    m=1.333333333333333*r*r*pi*e*r;

    for(pf=0;p>p_surf;r+=dr,pf++){
        if(pf>49995){
            return 0;
        }

        p1=fp(r,p,e,m);
        m1=fm(r,e);
        if((p+dr*p1/2)>p_surf) e=interp_p2rho_SI(EoS,p+dr*p1/2,&interpRef); else break;
        p2=fp(r+dr/2,p+dr*p1/2,e,m+dr*m1/2);
        m2=fm(r+dr/2,e);
        if((p+dr*p2/2)>p_surf) e=interp_p2rho_SI(EoS,p+dr*p2/2,&interpRef); else break;
        p3=fp(r+dr/2,p+dr*p2/2,e,m+dr*m2/2);
        m3=fm(r+dr/2,e);
        if((p+dr*p3)>p_surf) e=interp_p2rho_SI(EoS,p+dr*p3,&interpRef); else break;
        p4=fp(r+dr,p+dr*p3,e,m+dr*m3);
        m4=fm(r+dr,e);
        p=p+dr*(p1+2*p2+2*p3+p4)/6;
        e=interp_p2rho_SI(EoS,p,&interpRef);
        m=m+dr*(m1+2*m2+2*m3+m4)/6;
    }
    return m/Msun;
}

double getMmax_fm(EoS_t *EoS) {
    int n;
    double m0,m1,ma,mb,mc,md,E0,E1,Ea,Eb,Ec,Ed,Mmax,dE=5e15;

    E1=(EoS->RhomaxSI<4.2e18)?EoS->RhomaxSI:4.2e18;
    E0=E1-dE;
    m1=getM_fm(EoS,E1);
    m0=getM_fm(EoS,E0);
    if(m0-m1<=0) {
        EoS->Rhoc_MmaxSI=E1;
        EoS->Mmax=m1;
        return m1;
    }
    for(n=0;n<7;n++) {
        E1=E1*0.72;
        E0=E1-dE;
        m1=getM_fm(EoS,E1);
        m0=getM_fm(EoS,E0);
        if(m0-m1<=0) break;
    }
    if(m0-m1>0) {
        EoS->Rhoc_MmaxSI=E0;
        EoS->Mmax=m0;
        return m0;
    }
    Ea=E1;Eb=E1*1.3888888888889;
    Ec=Ea+0.382*(Eb-Ea);
    Ed=Eb-Ec+Ea;
    mc=getM_fm(EoS,Ec);
    md=getM_fm(EoS,Ed);
    for(n=0;n<15;n++) {
        if(mc-md>0) {
            Eb=Ed;mb=md;
            Ed=Ec;md=mc;
            Ec=Ea+0.382*(Eb-Ea);
            mc=getM_fm(EoS,Ec);
        }else{
            Ea=Ec;ma=mc;
            Ec=Ed;mc=md;
            Ed=Ea+0.618*(Eb-Ea);
            md=getM_fm(EoS,Ed);
        }
    }
    EoS->Rhoc_MmaxSI=(mc-md>0)?(0.5*(Ea+Ed)):(0.5*(Ec+Eb));
    EoS->Mmax=getM_fm(EoS,EoS->Rhoc_MmaxSI);
    return(EoS->Mmax);
}

double M2Rhoc_fm(EoS_t *EoS, double fM) {
    double Mmax,Mmin,E0,E1,E2,m0,m1,m2;
    int n;
    Mmax=getMmax_fm(EoS);
    Mmin=getM_fm(EoS,5e17);
    if(Mmax-fM<0) return EoS->Rhoc_MmaxSI;
    if(Mmin>fM) {
        fprintf(stderr,"error: attempt to find a low mass neutron star! EoS path: %s Mmax=%lf\n", EoS->FilePath, EoS->Mmax);
        return 5e17;
    }
    E0=5e17;E2=EoS->Rhoc_MmaxSI;m0=Mmin;m2=Mmax;
    for(n=0;n<25;n++) {
        E1=E2-(m2-fM)*(E2-E0)/(m2-m0);
        m1=getM_fm(EoS,E1);
        if(m1-fM<0){
            E0=E1;m0=m1;
        }else{
            E2=E1;m2=m1;
        }
    }
    return(E0+(E2-E0)*(fM-m0)/(m2-m0));
}

void * m2fmodeGate_fm(void *args){
    struct fmodeParams_t *thisParams;
    thisParams = (struct fmodeParams_t * ) args;
    loadEoS(thisParams->EoS,thisParams->EoS->FilePath);
    double rhoc=M2Rhoc_fm(thisParams->EoS,thisParams->fM);
    fmode(thisParams->Results, thisParams->EoS, rhoc, thisParams->auxSpace);
    return 0x00;
}

int m2fmode_mt(CompactStar_t *Results, double fM, Path_t *EoSlist, int EoSlistLen, int threads) {
    threads=threads<EoSlistLen?threads:EoSlistLen;
    if(threads<1){
        return 0;
    }
    auxSpace_t auxSpace[threads];
    struct fmodeParams_t fmodeParams[threads];
    EoS_t EoS[threads];
    pthread_t fmodeThread[threads];
    int pf,section=0;

    /*register auxiliary space for each thread*/
    for(pf=0;pf<threads;pf++){
        auxSpace[pf].mfile=(double*)malloc(50000*sizeof(double));
        auxSpace[pf].pfile=(double*)malloc(50000*sizeof(double));
        auxSpace[pf].rhofile=(double*)malloc(50000*sizeof(double));
        auxSpace[pf].Bfile=(double*)malloc(50000*sizeof(double));
        auxSpace[pf].Bcor=(double*)malloc(50000*sizeof(double));
        auxSpace[pf].Wfile=(double*)malloc(50000*sizeof(double));
        auxSpace[pf].Wfile1=(double*)malloc(50000*sizeof(double));
        auxSpace[pf].Wfile2=(double*)malloc(50000*sizeof(double));
        auxSpace[pf].Vfile=(double*)malloc(50000*sizeof(double));
        auxSpace[pf].Vfile1=(double*)malloc(50000*sizeof(double));
        auxSpace[pf].Vfile2=(double*)malloc(50000*sizeof(double));
    }
    while(section+threads<=EoSlistLen){
        for(pf=0;pf<threads;pf++,section++){
            strcpy(EoS[pf].FilePath,EoSlist[section].FilePath);
            fmodeParams[pf].EoS=&EoS[pf];
            fmodeParams[pf].fM=fM;
            fmodeParams[pf].Results=&Results[section];
            fmodeParams[pf].auxSpace=&auxSpace[pf];
            while(pthread_create(&fmodeThread[pf],NULL,m2fmodeGate_fm,&fmodeParams[pf])){usleep(50000);}
        }
        for(pf=0;pf<threads;pf++){
            pthread_join(fmodeThread[pf],NULL);
        }
    }

    /*release the aux space*/
    for(pf=0;pf<threads;pf++){
        free(auxSpace[pf].mfile);
        free(auxSpace[pf].pfile);
        free(auxSpace[pf].rhofile);
        free(auxSpace[pf].Bfile);
        free(auxSpace[pf].Bcor);
        free(auxSpace[pf].Wfile);
        free(auxSpace[pf].Wfile1);
        free(auxSpace[pf].Wfile2);
        free(auxSpace[pf].Vfile);
        free(auxSpace[pf].Vfile1);
        free(auxSpace[pf].Vfile2);
    }

    m2fmode_mt(&Results[section],fM,&EoSlist[section],EoSlistLen-section,threads);
    return 0;
}

void setIndex(int *index, int len){
	while(len-- > 0){
		index[len]=len;
	}
	return;
}

void rec2arxiv(ArXiv_t *arxiv, double RhocSI, double M){
	arxiv->RhocSI[arxiv->length]=RhocSI;
	arxiv->M[arxiv->length]=M;
	arxiv->length++;
	return;
}

int arcSimSort(int head, int tail, int* index, double* data){
    int pf, spf=0, unsort, swap;
    for(unsort=head;unsort<tail+1;unsort++){
        pf=unsort;spf=unsort;
        while(pf<tail+1){
                spf=data[index[pf]]<data[index[spf]]?pf:spf;
                pf++;
        }
        swap=index[spf];
        index[spf]=index[unsort];
        index[unsort]=swap;
    }
    return 0;
}

int interp_ArXiv_M2Rhoc_arr(double *RhocGuess, ArXiv_t *Arxiv, double *massArr, int arrayLen) {
	int pf, *mid=Arxiv->index;
	while(arrayLen-- > 0){
		for(pf=Arxiv->length-1; pf>0; pf--){
			if(massArr[arrayLen]>Arxiv->M[mid[pf-1]]){
				break;
			}
		}
		RhocGuess[arrayLen]=(massArr[arrayLen]-(Arxiv->M[mid[pf-1]]) )*((Arxiv->RhocSI[mid[pf]])-(Arxiv->RhocSI[mid[pf-1]]))/((Arxiv->M[mid[pf]])-(Arxiv->M[mid[pf-1]]))+(Arxiv->RhocSI[mid[pf-1]]);
	}
	return 0;
}

void *getM_fmGate (void *args) {
	struct getM_params_t *thisParams;
	thisParams=(struct getM_params_t *) args;
	*(thisParams->M)=getM_fm(thisParams->EoS,thisParams->RhocSI);
	return 0x00;
}

int getM_fm_mt(double *massArr, EoS_t *EoS, double *RhocSI, int arrayLen, int threads) {
	threads=threads<arrayLen?threads:arrayLen;
	if(threads<1){
		return 0;
	}
	struct getM_params_t gm_params[8];
	pthread_t getMthread[threads];
	int pf,section=0;

	while(section+threads<=arrayLen){
		for(pf=0;pf<threads;pf++,section++){
			gm_params[pf].RhocSI=RhocSI[section];
			gm_params[pf].M=&massArr[section];
			gm_params[pf].EoS=EoS;
			while(pthread_create(&getMthread[pf],NULL,getM_fmGate,&gm_params[pf])){usleep(50000);}
		}
		for(pf=0;pf<threads;pf++){
			pthread_join(getMthread[pf],NULL);
		}
	}

	getM_fm_mt(&massArr[section],EoS,&RhocSI[section],arrayLen-section,threads);
	return 0;
}

int M2Rhoc_Arr_fm(double *RhocSI, EoS_t *EoS, double *massArr, int arrayLen, int threads) {
	threads=threads<arrayLen?threads:arrayLen;
	threads=threads<2?2:threads;
	double ma,mb,mc,md,Ea,Eb,Ec,Ed,dE=5e16;
    double Mmax,Mmin,E0,E1,E2,m0,m1,m2;
    int n, pf;

	double RhocFrame[threads], massFrame[threads], RhocGuess[arrayLen], massGuess[arrayLen];
	double tmpMassM0M1[2], tmpRhocE0E1[2];

	ArXiv_t thisArxiv={0};
	setIndex(thisArxiv.index,sizeof(thisArxiv.index));
    
	tmpRhocE0E1[1]=(EoS->RhomaxSI<4.2e18)?EoS->RhomaxSI:4.2e18;
	tmpRhocE0E1[0]=tmpRhocE0E1[1]-dE;
	getM_fm_mt(tmpMassM0M1, EoS, tmpRhocE0E1, 2, 2);

    E1=tmpRhocE0E1[1];m1=tmpMassM0M1[1];rec2arxiv(&thisArxiv,E1,m1);
    E0=tmpRhocE0E1[0];m0=tmpMassM0M1[0];rec2arxiv(&thisArxiv,E0,m0);
    if(m0-m1<=0) {
        EoS->Rhoc_MmaxSI=E1;
        EoS->Mmax=m1;
        goto MmaxFound;
    }
    for(n=0;n<7;n++) {
		tmpRhocE0E1[1] *= 0.72;
		tmpRhocE0E1[0] = tmpRhocE0E1[1]-dE;
		getM_fm_mt(tmpMassM0M1, EoS, tmpRhocE0E1, 2, 2);

		E1=tmpRhocE0E1[1];m1=tmpMassM0M1[1];rec2arxiv(&thisArxiv,E1,m1);
		E0=tmpRhocE0E1[0];m0=tmpMassM0M1[0];rec2arxiv(&thisArxiv,E0,m0);
        if(m0-m1<=0) break;
    }
    if(m0-m1>0) {
        EoS->Rhoc_MmaxSI=E0;
        EoS->Mmax=m0;
        goto MmaxFound;
    }
    Ea=E1;Eb=E1*1.3888888888889;
    Ec=Ea+0.382*(Eb-Ea);
    Ed=Eb-Ec+Ea;
	tmpRhocE0E1[0]=Ec;tmpRhocE0E1[1]=Ed;
	getM_fm_mt(tmpMassM0M1, EoS, tmpRhocE0E1, 2, 2);

	mc=tmpMassM0M1[0];rec2arxiv(&thisArxiv,Ec,mc);
	md=tmpMassM0M1[1];rec2arxiv(&thisArxiv,Ed,md);

    for(n=0;n<15;n++) {
        if(mc-md>0) {
            Eb=Ed;mb=md;
            Ed=Ec;md=mc;
            Ec=Ea+0.382*(Eb-Ea);
            mc=getM_fm(EoS,Ec);rec2arxiv(&thisArxiv,Ec,mc);
        }else{
            Ea=Ec;ma=mc;
            Ec=Ed;mc=md;
            Ed=Ea+0.618*(Eb-Ea);
            md=getM_fm(EoS,Ed);rec2arxiv(&thisArxiv,Ed,md);
        }
    }
    EoS->Rhoc_MmaxSI=(mc-md>0)?(0.5*(Ea+Ed)):(0.5*(Ec+Eb));
    EoS->Mmax=getM_fm(EoS,EoS->Rhoc_MmaxSI);rec2arxiv(&thisArxiv,EoS->Rhoc_MmaxSI,EoS->Mmax);
	Mmin=getM_fm(EoS,5e17);rec2arxiv(&thisArxiv,5e17,Mmin);

MmaxFound:

	for(n=0;n<threads;n++){
		RhocFrame[n]=5e17+(n+1)*(EoS->Rhoc_MmaxSI-5e17)/(threads+1);
	}
	getM_fm_mt(massFrame, EoS, RhocFrame, threads, threads);

	for(n=0;n<threads;n++){
		rec2arxiv(&thisArxiv,RhocFrame[n],massFrame[n]);
	}

	arcSimSort(0, thisArxiv.length-1, thisArxiv.index, thisArxiv.M);

	for(n=0;n<8;n++){
		interp_ArXiv_M2Rhoc_arr(RhocGuess, &thisArxiv, massArr, arrayLen);
		getM_fm_mt(massGuess, EoS, RhocGuess, arrayLen, threads);
		for(pf=0;pf<arrayLen;pf++){
			rec2arxiv(&thisArxiv,RhocGuess[pf],massGuess[pf]);
		}
		arcSimSort(0, thisArxiv.length-1, thisArxiv.index, thisArxiv.M);
	}

	interp_ArXiv_M2Rhoc_arr(RhocSI, &thisArxiv, massArr, arrayLen);

	return 0;
}

void *getM_sGate (void *args) {
	struct getM_params_t *thisParams;
	thisParams=(struct getM_params_t *) args;
	*(thisParams->M)=getM_s(thisParams->EoS,thisParams->RhocSI);
	return 0x00;
}


int getM_s_mt(double *massArr, EoS_t *EoS, double *RhocSI, int arrayLen, int threads) {
	threads=threads<arrayLen?threads:arrayLen;
	if(threads<1){
		return 0;
	}
	struct getM_params_t gm_params[8];
	pthread_t getMthread[threads];
	int pf,section=0;

	while(section+threads<=arrayLen){
		for(pf=0;pf<threads;pf++,section++){
			gm_params[pf].RhocSI=RhocSI[section];
			gm_params[pf].M=&massArr[section];
			gm_params[pf].EoS=EoS;
			while(pthread_create(&getMthread[pf],NULL,getM_sGate,&gm_params[pf])){usleep(50000);}
		}
		for(pf=0;pf<threads;pf++){
			pthread_join(getMthread[pf],NULL);
		}
	}

	getM_s_mt(&massArr[section],EoS,&RhocSI[section],arrayLen-section,threads);
	return 0;
}

int M2Rhoc_Arr_s(double *RhocSI, EoS_t *EoS, double *massArr, int arrayLen, int threads) {
	threads=threads<arrayLen?threads:arrayLen;
	threads=threads<2?2:threads;
	double ma,mb,mc,md,Ea,Eb,Ec,Ed,dE=5e16;
    double Mmax,Mmin,E0,E1,E2,m0,m1,m2;
    int n, pf;

	double RhocFrame[threads], massFrame[threads], RhocGuess[arrayLen], massGuess[arrayLen];
	double tmpMassM0M1[2], tmpRhocE0E1[2];

	ArXiv_t thisArxiv={0};
	setIndex(thisArxiv.index,sizeof(thisArxiv.index));
    
	tmpRhocE0E1[1]=(EoS->RhomaxSI<4.2e18)?EoS->RhomaxSI:4.2e18;
	tmpRhocE0E1[0]=tmpRhocE0E1[1]-dE;
	getM_s_mt(tmpMassM0M1, EoS, tmpRhocE0E1, 2, 2);

    E1=tmpRhocE0E1[1];m1=tmpMassM0M1[1];rec2arxiv(&thisArxiv,E1,m1);
    E0=tmpRhocE0E1[0];m0=tmpMassM0M1[0];rec2arxiv(&thisArxiv,E0,m0);
    if(m0-m1<=0) {
        EoS->Rhoc_MmaxSI=E1;
        EoS->Mmax=m1;
        goto MmaxFound;
    }
    for(n=0;n<7;n++) {
		tmpRhocE0E1[1] *= 0.72;
		tmpRhocE0E1[0] = tmpRhocE0E1[1]-dE;
		getM_s_mt(tmpMassM0M1, EoS, tmpRhocE0E1, 2, 2);

		E1=tmpRhocE0E1[1];m1=tmpMassM0M1[1];rec2arxiv(&thisArxiv,E1,m1);
		E0=tmpRhocE0E1[0];m0=tmpMassM0M1[0];rec2arxiv(&thisArxiv,E0,m0);
        if(m0-m1<=0) break;
    }
    if(m0-m1>0) {
        EoS->Rhoc_MmaxSI=E0;
        EoS->Mmax=m0;
        goto MmaxFound;
    }
    Ea=E1;Eb=E1*1.3888888888889;
    Ec=Ea+0.382*(Eb-Ea);
    Ed=Eb-Ec+Ea;
	tmpRhocE0E1[0]=Ec;tmpRhocE0E1[1]=Ed;
	getM_s_mt(tmpMassM0M1, EoS, tmpRhocE0E1, 2, 2);

	mc=tmpMassM0M1[0];rec2arxiv(&thisArxiv,Ec,mc);
	md=tmpMassM0M1[1];rec2arxiv(&thisArxiv,Ed,md);

    for(n=0;n<15;n++) {
        if(mc-md>0) {
            Eb=Ed;mb=md;
            Ed=Ec;md=mc;
            Ec=Ea+0.382*(Eb-Ea);
            mc=getM_fm(EoS,Ec);rec2arxiv(&thisArxiv,Ec,mc);
        }else{
            Ea=Ec;ma=mc;
            Ec=Ed;mc=md;
            Ed=Ea+0.618*(Eb-Ea);
            md=getM_fm(EoS,Ed);rec2arxiv(&thisArxiv,Ed,md);
        }
    }
    EoS->Rhoc_MmaxSI=(mc-md>0)?(0.5*(Ea+Ed)):(0.5*(Ec+Eb));
    EoS->Mmax=getM_fm(EoS,EoS->Rhoc_MmaxSI);rec2arxiv(&thisArxiv,EoS->Rhoc_MmaxSI,EoS->Mmax);
	Mmin=getM_fm(EoS,5e17);rec2arxiv(&thisArxiv,5e17,Mmin);

MmaxFound:

	for(n=0;n<threads;n++){
		RhocFrame[n]=5e17+(n+1)*(EoS->Rhoc_MmaxSI-5e17)/(threads+1);
	}
	getM_s_mt(massFrame, EoS, RhocFrame, threads, threads);

	for(n=0;n<threads;n++){
		rec2arxiv(&thisArxiv,RhocFrame[n],massFrame[n]);
	}

	arcSimSort(0, thisArxiv.length-1, thisArxiv.index, thisArxiv.M);

	for(n=0;n<8;n++){
		interp_ArXiv_M2Rhoc_arr(RhocGuess, &thisArxiv, massArr, arrayLen);
		getM_s_mt(massGuess, EoS, RhocGuess, arrayLen, threads);
		for(pf=0;pf<arrayLen;pf++){
			rec2arxiv(&thisArxiv,RhocGuess[pf],massGuess[pf]);
		}
		arcSimSort(0, thisArxiv.length-1, thisArxiv.index, thisArxiv.M);
	}

	interp_ArXiv_M2Rhoc_arr(RhocSI, &thisArxiv, massArr, arrayLen);

	return 0;
}
