#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>

#include "statNS.h"

static const double pi=3.14159265358979, 
			  G=6.673e-11,
			  c=2.99792458e8,
			  Msun=1.989e30,
			  Cri=1e-7,
			  Gc2=7.42471382405e-28,
			  Gc4=8.26110825251e-45,
			  Mscale=2.033931261665867e5;


struct EoS_t EoS;

struct ComputeStatus_t ComputeStatus={0,0,0,RungeKutta_RK5L_roll,RungeKutta_RK5L_join};

double interp_p2rho(double cp, double *VsOUT) {
	if(cp<EoS.Pmin) {
		return -1;
	}
	double logp=log10(cp),logrho,crho;
	int n=EoS.length-1;
	for(;n>0;n--)
	{
		if(logp>EoS.lgP[n-1]) {
			logrho=(logp-EoS.lgP[n-1])*(EoS.lgRho[n]-EoS.lgRho[n-1])/(EoS.lgP[n]-EoS.lgP[n-1])+EoS.lgRho[n-1];
			crho=pow(10,logrho);
			/*compute local speed of sound*/
			*VsOUT=cp/crho*(EoS.lgP[n]-EoS.lgP[n-1])/(EoS.lgRho[n]-EoS.lgRho[n-1]);
			return(crho);
		}
	}
	fprintf(stderr,"%s interrupted by interp_p2rho(), reason: small P!\n", EoS.FilePath);
	return -1;
}

double interp_rho2p(double crho) {
	if(crho<EoS.Rhomin){
		return -1;
	}
	double logp,logrho=log10(crho),cp;
	int n=EoS.length-1;
	for(;n>0;n--)
	{
		if(logrho>EoS.lgRho[n-1]) {
			logp=(logrho-EoS.lgRho[n-1])*(EoS.lgP[n]-EoS.lgP[n-1])/(EoS.lgRho[n]-EoS.lgRho[n-1])+EoS.lgP[n-1];
			cp=pow(10,logp);
			return(cp);
		}
	}
	fprintf(stderr,"%s interrupted by interp_rho2p(), reason: small Rho!\n", EoS.FilePath);
	return -1;
}

/*======================================================
*	Array adders are matrix operation function.
*	C language is not good at operating matrix!
*/
int RungeKutta_Array_add (struct RungeKutta_Array_t *Result, double h, struct RungeKutta_Array_t *K, struct RungeKutta_Array_t *X) {
	Result->P=X->P+h*(K->P);
	Result->M=X->M+h*(K->M);
	Result->I=X->I+h*(K->I);
	Result->y=X->y+h*(K->y);
	Result->Ag00=X->Ag00+h*(K->Ag00);
	return 0;
}

int RungeKutta_Array_adds (struct RungeKutta_Array_t *Result, double *h, struct RungeKutta_Array_t *K, struct RungeKutta_Array_t *X, int dim) {
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
int dFunc(struct RungeKutta_Array_t *K, double r, struct RungeKutta_Array_t *X){
	double Vs;
	double Rho=interp_p2rho(X->P,&Vs);
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

int RungeKutta_RK4_roll(struct RungeKutta_Array_t *K, double h, double r, struct RungeKutta_Array_t *X){
	struct RungeKutta_Array_t rkVar;
	if(dFunc(&K[0], r, X)){return -1;}
	RungeKutta_Array_add(&rkVar, h/2, &K[0], X);
	if(dFunc(&K[1], r+h/2, &rkVar)){return -1;}
	RungeKutta_Array_add(&rkVar, h/2, &K[1], X);
	if(dFunc(&K[2], r+h/2, &rkVar)){return -1;}
	RungeKutta_Array_add(&rkVar, h, &K[2], X);
	if(dFunc(&K[3], r+h, &rkVar)){return -1;}
	return 0;
}

int RungeKutta_RK4_join (struct RungeKutta_Array_t *Result, double h, struct RungeKutta_Array_t *K) {
	Result->P += h/6*(K[0].P+2*(K[1].P)+2*(K[2].P)+K[3].P);
	Result->M += h/6*(K[0].M+2*(K[1].M)+2*(K[2].M)+K[3].M);
	Result->I += h/6*(K[0].I+2*(K[1].I)+2*(K[2].I)+K[3].I);
	Result->y += h/6*(K[0].y+2*(K[1].y)+2*(K[2].y)+K[3].y);
	Result->Ag00 += h/6*(K[0].Ag00+2*(K[1].Ag00)+2*(K[2].Ag00)+K[3].Ag00);
	return 0;
}

int RungeKutta_RK4M_roll(struct RungeKutta_Array_t *K, double h, double r, struct RungeKutta_Array_t *X){
	struct RungeKutta_Array_t rkVar;
	double H_1[1]={h/2};
	double H_2[2]={0,h/2};
	double H_2M[2]={-h,2*h};
	double H_3[3]={0,0,h};
	if(dFunc(&K[0], r, X)){return -1;}
	RungeKutta_Array_adds(&rkVar, H_1, K, X,1);
	if(dFunc(&K[1], r+h/2, &rkVar)){return -1;}
	RungeKutta_Array_adds(&rkVar, H_2, K, X,2);
	if(dFunc(&K[2], r+h/2, &rkVar)){return -1;}
	RungeKutta_Array_adds(&rkVar, H_3, K, X,3);
	if(dFunc(&K[3], r+h, &rkVar)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_2M,K,X,2);
	if(dFunc(&K[4], r+h, &rkVar)){return -1;}
	return 0;
}

int RungeKutta_RK4M_join (struct RungeKutta_Array_t *Result, double h, struct RungeKutta_Array_t *K) {
	struct RungeKutta_Array_t rkVar;
	double H[4]={h/6,h/3,h/3,h/6};
	RungeKutta_Array_adds(&rkVar, H, K, Result, 4);
	RungeKutta_Array_adds(Result, H, K, &rkVar, 0);/*dim=0 means copy the array*/
	/*ComputeStatus.RkErr=(-2*(K[1].P)+2*(K[2].P)+(K[3].P)-(K[4].P))/(K[0].P+2*(K[1].P)+2*(K[2].P)+K[3].P);
	ComputeStatus.RkErr=ComputeStatus.RkErr>0?ComputeStatus.RkErr:(-ComputeStatus.RkErr);*/
	return 0;
}

int RungeKutta_RK5F_roll(struct RungeKutta_Array_t *K, double h, double r, struct RungeKutta_Array_t *X){
	struct RungeKutta_Array_t rkVar;
	double H_1[1]={0.25*h};
	double H_2[2]={0.093750*h,0.281250*h};
	double H_3[3]={1932.0*h/2179,-7200.0*h/2179,7296.0*h/2179};
	double H_4[4]={439.0*h/216,-8*h,3680.0*h/513,-845.0*h/4104};
	double H_5[5]={-8.0*h/27,2*h,-3544.0*h/2565,1859.0*h/4104,-11.0*h/40};
	if(dFunc(&K[0], r, X)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_1,K,X,1);
	if(dFunc(&K[1], r+0.25*h, &rkVar)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_2,K,X,2);
	if(dFunc(&K[2], r+0.375*h, &rkVar)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_3,K,X,3);
	if(dFunc(&K[3], r+12.0*h/13, &rkVar)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_4,K,X,4);
	if(dFunc(&K[4], r+h, &rkVar)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_5,K,X,5);
	if(dFunc(&K[5], r+0.5*h, &rkVar)){return -1;}
	return 0;
}

int RungeKutta_RK5F_join (struct RungeKutta_Array_t *Result, double h, struct RungeKutta_Array_t *K) {
	struct RungeKutta_Array_t rkVar;
	double H[6]={16.0/135*h, 0, 6656.0/12825*h, 28561.0/56430*h, - 0.18*h, 2.0/55*h};
	RungeKutta_Array_adds(&rkVar, H, K, Result, 6);
	RungeKutta_Array_adds(Result, H, K, &rkVar, 0);/*dim=0 means copy the array*/
	/*ComputeStatus.RkErr=(0.0027778*(K[0].P) - 0.0299415*(K[2].P)  - 0.0291999*(K[3].P) + 0.0200000*(K[4].P) +  0.0363636*(K[5].P))/(16.0/135*(K[0].P) + 6656.0/12825*(K[2].P) + 28561.0/56430*(K[3].P) - 0.18*(K[4].P) + 2.0/55*(K[5].P) );
	ComputeStatus.RkErr=ComputeStatus.RkErr>0?ComputeStatus.RkErr:(-ComputeStatus.RkErr);*/
	return 0;
}

int RungeKutta_RK5M_roll(struct RungeKutta_Array_t *K, double h, double r, struct RungeKutta_Array_t *X){
	struct RungeKutta_Array_t rkVar;
	double H_1[1]={0.2*h};
	double H_2[2]={0.075000*h,0.225000*h};
	double H_3[3]={0.3*h,-0.9*h,1.2*h};
	double H_4[4]={226.0/729*h,-25.0/27*h,880.0/729*h,55.0/729*h};
	double H_5[5]={-181.0/270*h,2.5*h,-266.0/297*h,-91.0/27*h,189.0/55*h};
	if(dFunc(&K[0], r, X)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_1,K,X,1);
	if(dFunc(&K[1], r+0.2*h, &rkVar)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_2,K,X,2);
	if(dFunc(&K[2], r+0.3*h, &rkVar)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_3,K,X,3);
	if(dFunc(&K[3], r+0.6*h, &rkVar)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_4,K,X,4);
	if(dFunc(&K[4], r+2.0*h/3, &rkVar)){return -1;}
	RungeKutta_Array_adds(&rkVar,H_5,K,X,5);
	if(dFunc(&K[5], r+h, &rkVar)){return -1;}
	return 0;
}

int RungeKutta_RK5M_join (struct RungeKutta_Array_t *Result, double h, struct RungeKutta_Array_t *K) {
	struct RungeKutta_Array_t rkVar;
	double H[6]={31.0/540*h, 0, 190.0/297*h, -145.0/108*h, 351.0/220*h, 0.05*h};
	RungeKutta_Array_adds(&rkVar, H, K, Result, 6);
	RungeKutta_Array_adds(Result, H, K, &rkVar, 0);/*dim=0 means copy the array*/
	/*ComputeStatus.RkErr=(-0.030556*(K[0].P) + 0.158730*(K[2].P)  - 0.763889*(K[3].P) + 0.675000*(K[4].P) - 0.039286*(K[5].P))/(31.0/540*(K[0].P) + 190.0/297*(K[2].P) - 145.0/108*(K[3].P) + 351.0/220*(K[4].P) + 0.05*(K[5].P) );
	ComputeStatus.RkErr=ComputeStatus.RkErr>0?ComputeStatus.RkErr:(-ComputeStatus.RkErr);*/
	return 0;
}

int RungeKutta_RK5L_roll(struct RungeKutta_Array_t *K, double h, double r, struct RungeKutta_Array_t *X){
	struct RungeKutta_Array_t rkVar;
    double H_1[1]={0.08333333333333333*h};
    double H_2[2]={-0.125*h, 0.375*h};
    double H_3[3]={0.6*h,-0.9*h,0.8*h};
    double H_4[4]={0.4875*h,-0.45*h,0.15*h,0.5625*h};
    double H_5[5]={-1.685714285714286*h,1.885714285714286*h,1.371428571428571*h,-1.714285714285714*h,1.142857142857143*h};
    if(dFunc(&K[0], r, X)){return -1;}
    RungeKutta_Array_adds(&rkVar,H_1,K,X,1);
    if(dFunc(&K[1], r+0.083333333333333333*h, &rkVar)){return -1;}
    RungeKutta_Array_adds(&rkVar,H_2,K,X,2);
    if(dFunc(&K[2], r+0.25*h, &rkVar)){return -1;}
    RungeKutta_Array_adds(&rkVar,H_3,K,X,3);
    if(dFunc(&K[3], r+0.5*h, &rkVar)){return -1;}
    RungeKutta_Array_adds(&rkVar,H_4,K,X,4);
    if(dFunc(&K[4], r+0.75*h, &rkVar)){return -1;}
    RungeKutta_Array_adds(&rkVar,H_5,K,X,5);
    if(dFunc(&K[5], r+h, &rkVar)){return -1;}
    return 0;
}

int RungeKutta_RK5L_join (struct RungeKutta_Array_t *Result, double h, struct RungeKutta_Array_t *K) {
	struct RungeKutta_Array_t rkVar;
    double H[6]={0.077777777777777777778*h,0,0.355555555555555555556*h,0.13333333333333333333*h,0.355555555555555555556*h,0.077777777777777777777778*h};
    RungeKutta_Array_adds(&rkVar, H, K, Result, 6);
    RungeKutta_Array_adds(Result, H, K, &rkVar, 0);/*dim=0 means copy the array*/
    return 0;
}

int loadEoS(char *Path) {
	FILE *inf;
	int n;

	if((inf=fopen(Path,"r"))==NULL) {
		fprintf(stderr,"cannot open %s\n",Path);
		return -1;
	}
	sprintf(EoS.FilePath,"%s",Path);
	EoS.length=0;
	while(fscanf(inf,"%lf",&EoS.lgRho_SI[EoS.length])==1) {
		fscanf(inf,"%lf%*[^\n]",&EoS.lgP_SI[EoS.length]);
		EoS.length++;
	}
	fclose(inf);
	EoS.RhomaxSI=pow(10,EoS.lgRho_SI[EoS.length-1]);
	EoS.RhominSI=pow(10,EoS.lgRho_SI[0]);

/*Convert to natural unit*/
	for(n=0;n<EoS.length;n++) {
		EoS.lgRho[n]=EoS.lgRho_SI[n]-10.175607290470733;
		EoS.lgP[n]=EoS.lgP_SI[n]-27.129849799910058;
	}
	EoS.Pmax=pow(10,EoS.lgP[EoS.length-1]);
	EoS.Pmin=pow(10,EoS.lgP[0]);
	EoS.Rhomax=pow(10,EoS.lgRho[EoS.length-1]);
	EoS.Rhomin=pow(10,EoS.lgRho[0]);
	return 0;
}

int solveTOV(double RhocSI, struct CompactStar_t *Results) {
	
	struct RungeKutta_Array_t K[8];
	double h;

	double rho=RhocSI*6.6741e-11, r=0.125/c, r_offset;
	double m=4.0/3*r*r*r*rho, p=interp_rho2p(rho), y=2.0, I=0.0, Ag00=1.0;
	int pf=0;
	struct RungeKutta_Array_t X={p,m,I,Ag00,y};

	/*core section with proceeding stepsize*/
	h=0.125/c;
	for(pf=0;pf<64;pf++) {
		if(ComputeStatus.K_roll(K, h, r, &X)) break;
		ComputeStatus.X_join(&X,h, K);
		r+=h;
	}
	h=1.0/c;
	for(pf=0;pf<128;pf++) {
		if(ComputeStatus.K_roll(K, h, r, &X)) break;
		ComputeStatus.X_join(&X,h, K);
		r+=h;
	}
	h=4.0/c;
	for(pf=0;pf<256;pf++) {
		if(ComputeStatus.K_roll(K, h, r, &X)) break;
		ComputeStatus.X_join(&X,h, K);
		r+=h;
	}
	/*end of core section*/

	/*regular computation*/
	h=16.0/c;
	while(r<1e-4) {
		if(ComputeStatus.K_roll(K, h, r, &X)) break;
		ComputeStatus.X_join(&X,h, K);
		r+=h;
	}

	/*use a small stepsize to reach the surface*/
	h=2.0/c;
	while(r<1e-4) {
		if(ComputeStatus.K_roll(K, h, r, &X)) break;
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
	Results->I=I/(m*r*r*sqrt(Bg00));
	Results->r=r*c*1e-3;
	Results->M=m*Mscale;
	Results->k2=1.6*bt*bt*bt*bt*bt*(1-2*bt)*(1-2*bt)*(2-y+2*bt*(y-1))/(2*bt*(6-3*y+3*bt*(5*y-8))+4*bt*bt*bt*(13-11*y+bt*(3*y-2)+2*bt*bt*(1+y))+3*(1-2*bt)*(1-2*bt)*(2-y+2*bt*(y-1))*log(1-2*bt));
	Results->Lambda=9495*(Results->k2)*pow(Results->r/10,5)/pow(Results->M,5);

	return 0;
}

void *solveTOVGate (void *args){
	struct solveTOVparams_t *thisParams;
	thisParams=(struct solveTOVparams_t *) args;
	solveTOV(thisParams->RhocSI,thisParams->Results);
	return 0x00;
}

int solveTOV_mt(int threadNum, struct CompactStar_t *Results, double *cdenArray, int arrayLen) {
	threadNum=threadNum<arrayLen?threadNum:arrayLen;
	if(threadNum<1){
		return 0;
	}
	struct solveTOVparams_t solveTOVparams[threadNum];
	pthread_t solveTOVthread[threadNum];
	int pf,section=0;

	while(section+threadNum<=arrayLen){
		for(pf=0;pf<threadNum;pf++,section++){
			solveTOVparams[pf].RhocSI=cdenArray[section];
			solveTOVparams[pf].Results=&Results[section];
			while(pthread_create(&solveTOVthread[pf],NULL,solveTOVGate,&solveTOVparams[pf])){usleep(50000);}
		}
		for(pf=0;pf<threadNum;pf++){
			pthread_join(solveTOVthread[pf],NULL);
		}
	}

	solveTOV_mt(threadNum,&Results[section],&cdenArray[section],arrayLen-section);
	return 0;
}

double getMmax() {
	int n;
	double m0,m1,ma,mb,mc,md,E0,E1,Ea,Eb,Ec,Ed,Mmax,dE=1e12;
	struct CompactStar_t nsVar;

	E1=(EoS.RhomaxSI<4.2e18)?EoS.RhomaxSI:4.2e18;
	E0=E1-dE;
	solveTOV(E1,&nsVar);
	m1=nsVar.M;
	solveTOV(E0,&nsVar);
	m0=nsVar.M;
	if(m0-m1<=0) {
		EoS.Rhoc_MmaxSI=E1;
		EoS.Mmax=m1;
		goto getMmaxRTB;
	}
	for(n=0;n<7;n++) {
		E1=E1*0.72;
		E0=E1-dE;
		solveTOV(E1,&nsVar);
		m1=nsVar.M;
		solveTOV(E0,&nsVar);
		m0=nsVar.M;
		if(m0-m1<=0) break;
	}
	if(m0-m1>0) {
		EoS.Rhoc_MmaxSI=E0;
		EoS.Mmax=m0;
		goto getMmaxRTB;
	}
	Ea=E1;Eb=E1*1.3888888889;
	Ec=Ea+0.382*(Eb-Ea);
	Ed=Eb-Ec+Ea;
	solveTOV(Ec,&nsVar);
	mc=nsVar.M;
	solveTOV(Ed,&nsVar);
	md=nsVar.M;
	for(n=0;n<15;n++) {
		if(mc-md>0) {
			Eb=Ed;mb=md;
			Ed=Ec;md=mc;
			Ec=Ea+0.382*(Eb-Ea);
			solveTOV(Ec,&nsVar);
			mc=nsVar.M;
		}else{
			Ea=Ec;ma=mc;
			Ec=Ed;mc=md;
			Ed=Ea+0.618*(Eb-Ea);
			solveTOV(Ed,&nsVar);
			md=nsVar.M;
		}
	}
	EoS.Rhoc_MmaxSI=(mc-md>0)?(0.5*(Ea+Ed)):(0.5*(Ec+Eb));
	solveTOV(EoS.Rhoc_MmaxSI,&nsVar);
	EoS.Mmax=nsVar.M;
getMmaxRTB:
	return(nsVar.M);
}

double M2Rhoc(double fM) {
	double Mmax,Mmin,E0,E1,E2,m0,m1,m2;
	int n;
	struct CompactStar_t nsVar;
	Mmax=getMmax();
	solveTOV(5e17,&nsVar);
	Mmin=nsVar.M;
	if(Mmax-fM<0) return EoS.Rhoc_MmaxSI;
	if(Mmin>fM) {
		fprintf(stderr,"warning: attempt to find a low mass neutron star!\n");
		return 5e17;
	}
	E0=5e17;E2=EoS.Rhoc_MmaxSI;m0=Mmin;m2=Mmax;
	for(n=0;n<25;n++) {
		/*E1=0.5*(E0+E2);*/
		/*The binsection method is replaced by false point method here*/
		E1=E2-(m2-fM)*(E2-E0)/(m2-m0);
		solveTOV(E1,&nsVar);
		m1=nsVar.M;
		if(m1-fM<0)
		{
			E0=E1;m0=m1;
		}else{
			E2=E1;m2=m1;
		}
	}
	return(E0+(E2-E0)*(fM-m0)/(m2-m0));
}

double ch(double cp) {
	double logp=log10(cp),logrho,crho;
	int n=EoS.length-1;
	for(;n>0;n--)
	{
		if(logp>EoS.lgP_SI[n-1])
		{
			logrho=(logp-EoS.lgP_SI[n-1])*(EoS.lgRho_SI[n]-EoS.lgRho_SI[n-1])/(EoS.lgP_SI[n]-EoS.lgP_SI[n-1])+EoS.lgRho_SI[n-1];
			crho=pow(10,logrho);
			return(crho);
		}
	}
	fprintf(stderr,"%s meet small pressure at ch()\n", EoS.FilePath);
	return 0;
}

double che(double crho) {
	double logp,logrho=log10(crho),cp;
	int n=EoS.length-1;
	for(;n>0;n--)
	{
		if(logrho>EoS.lgRho_SI[n-1])
		{
			logp=(logrho-EoS.lgRho_SI[n-1])*(EoS.lgP_SI[n]-EoS.lgP_SI[n-1])/(EoS.lgRho_SI[n]-EoS.lgRho_SI[n-1])+EoS.lgP_SI[n-1];
			cp=pow(10,logp);
			return(cp);
		}
	}
	fprintf(stderr,"%s meet small density at che()\n", EoS.FilePath);
	return 0;
}

double gf(double e){
	double le, p, gamma;
	int i=EoS.length-1;
	e=e/G*c*c;
	le=log10(e);
	for(;i>0;i--){
  		if(le>EoS.lgRho_SI[i-1]){
    		p=EoS.lgP_SI[i-1]*(le-EoS.lgRho_SI[i])/(EoS.lgRho_SI[i-1]-EoS.lgRho_SI[i])+EoS.lgP_SI[i]*(le-EoS.lgRho_SI[i-1])/(EoS.lgRho_SI[i]-EoS.lgRho_SI[i-1]);
    		p=pow(10,p);
    		gamma=(e+p/c/c)/e*(EoS.lgP_SI[i]-EoS.lgP_SI[i-1])/(EoS.lgRho_SI[i]-EoS.lgRho_SI[i-1]);
    		return(gamma);
		}
	}
	return((pow(10,EoS.lgRho_SI[0])+pow(10,EoS.lgP_SI[0])/c/c)/e*(EoS.lgP_SI[i]-EoS.lgP_SI[i-1])/(EoS.lgRho_SI[i]-EoS.lgRho_SI[i-1]));
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

int fmode(double RhocSI, struct CompactStar_t *Results, struct auxSpace_t *auxSpace){
	double dr=0.5;
	double r,r0=1,R,RR,drx,rx;
	double ne,p,e,m,mR,A,B=1.0,BR,Bfactor,m1,m2,m3,m4,p1,p2,p3,p4,B1,B2,B3,B4,pressure,I=0,J,DDf,Df=0,f=1;
	double H1,H0,K,W,X,F,V,Dv,gamma,Dp,N;
	double DH11,DK1,DW1,DX1,H01,H02,K1,K2,x,X1,X2,Xp1,Xp2,W1,W2,V1,V2,V01,V02,V0;
	double w,wcheck,o[2],wi,aR,bR,gR,hR,kR,n,Y1,Y2,Z,DZ,DDZ,VZ,Ar1,Ar2,Ai1,Ai2,ar,ai,Ar[2],Ai[2],Br[2],Bi[2];
	int t,q,wpf,rpf,pfEnd,n4,status=0;

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
	
	pressure=pow(10,EoS.lgP_SI[0]);
	r=r0;
	ne=RhocSI;
	e=RhocSI;
	p=che(e);
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
		if((p+dr*p1/2)>pressure) e=ch(p+dr*p1/2); else break;
		B1=Bf(r,p,m,B);
		p2=fp(r+dr/2,p+dr*p1/2,e,m+dr*m1/2);
		m2=fm(r+dr/2,e);
		if((p+dr*p2/2)>pressure) e=ch(p+dr*p2/2); else break;
		B2=Bf(r+dr/2,p+dr*p1/2,m+dr*m1/2,B+dr*B1/2);
		p3=fp(r+dr/2,p+dr*p2/2,e,m+dr*m2/2);
		m3=fm(r+dr/2,e);
		if((p+dr*p3)>pressure) e=ch(p+dr*p3); else break;
		B3=Bf(r+dr/2,p+dr*p2/2,m+dr*m2/2,B+dr*B2/2);
		p4=fp(r+dr,p+dr*p3,e,m+dr*m3);
		m4=fm(r+dr,e);
		B4=Bf(r+dr,p+dr*p3,m+dr*m3,B+dr*B3);
		J=-4*pi*(e+p/c/c)*Gc2*r*A;
		DDf=-(4/r*Df+J*Df+4/r*J*f);
		I=I-2.0/3*f*J/sqrt(A*B)*r*r*r*dr;
		p=p+dr*(p1+2*p2+2*p3+p4)/6;
		e=ch(p);
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
	gamma=(EoS.lgP_SI[1]-EoS.lgP_SI[0])/(EoS.lgRho_SI[1]-EoS.lgRho_SI[0]);
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
		X=(e+p)*sqrt(B)*((4*pi/3*(e+3*p)-w*w/B/2)*W+0.5*K);
		H1=(2*K+16*pi*(e+p)*W)/3;
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
			gamma=gf(e);
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
		X=(e+p)*sqrt(B)*((4*pi/3*(e+3*p)-w*w/B/2.0)*W+0.5*K);
		H1=(4.0*K+16*pi*(e+p)*W)/6.0;
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
			gamma=gf(e);
			H0=H0f(r,B,X,m,p,A,H1,K,w);
			V=Vf(r,w,e,p,B,A,Dp,W,H0,X);
			if(r==r0)V02=V;
			if((wcheck-w)/w>-Cri&&(wcheck-w)/w<Cri){
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
		n=1.5;
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
	Results->Rho=ne;
	Results->M=mR/(Msun*Gc2);
	Results->r=RR/1000;
	Results->freq=w*c/(2000*pi);
	Results->dampTime=1/(wi*c);
funcExit:
	return status;
}

void * fmodeGate (void *args){
	struct fmodeParams_t *thisParams;
	/*safely transfer the parameters from global dashboard fmodeParams to thisParams*/
	thisParams=(struct fmodeParams_t *) args;
	/*start the computation*/
	fmode(thisParams->RhocSI,thisParams->Results,thisParams->auxSpace);
	return 0x00;
}

int fmode_mt(int threadNum, struct CompactStar_t *Results, double *cdenArray, int arrayLen) {
	threadNum=threadNum<arrayLen?threadNum:arrayLen;
	if(threadNum<1){
		return 0;
	}
	struct auxSpace_t auxSpace[threadNum];
	struct fmodeParams_t fmodeParams[threadNum];
	pthread_t fmodeThread[threadNum];
	int pf,section=0;
	/*register auxiliary space for each thread*/
	for(pf=0;pf<threadNum;pf++){
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
	
	while(section+threadNum<=arrayLen){
		for(pf=0;pf<threadNum;pf++,section++){
			fmodeParams[pf].RhocSI=cdenArray[section];
			fmodeParams[pf].Results=&Results[section];
			fmodeParams[pf].auxSpace=&auxSpace[pf];
			while(pthread_create(&fmodeThread[pf],NULL,fmodeGate,&fmodeParams[pf])){usleep(50000);}
		}
		for(pf=0;pf<threadNum;pf++){
			pthread_join(fmodeThread[pf],NULL);
		}
	}

	/*release the aux space*/
	for(pf=0;pf<threadNum;pf++){
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
	fmode_mt(threadNum,&Results[section],&cdenArray[section],arrayLen-section);
	return 0;
}
