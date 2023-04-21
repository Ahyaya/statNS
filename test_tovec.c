#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "include/statNS.h"

static const double pi=3.14159265358979, 
					G=6.673e-11,
					c=2.99792458e8,
					Msun=1.989e30,
					Mscale=2.033931261665867e5;

int solveTOV_db (char dblogPath[], double errTolerance, EoS_t *EoS, double RhocSI){
	
	RK_plan_t RK_plan_db = {RungeKutta_RK5M_roll, NULL, RungeKutta_RK5M_join_ec};

	double rho=RhocSI*6.6741e-11, r=0.125/c, r_offset;
	double m=4.0/3*r*r*r*rho, p=interp_rho2p(EoS,rho), y=2.0, I=0.0, Ag00=1.0;
	int pf, interpRef[16];
	/*maximum ref space for RK method and interp is 16*/
	for(pf=0;pf<16;++pf){
		interpRef[pf] = EoS->length-1;
	}
	int *ref=interpRef;

	RK_Arr_t X={p,m,I,Ag00,y};
	RK_Arr_t nextX = X;
	RK_Arr_t K[8] = {0};
	double err = 0;

	FILE *fp = fopen(dblogPath, "w");
	if(fp==NULL){
		fprintf(stderr, "Error: unable to write %s\n", dblogPath);
		return -1;
	}

	double h = 1.0/c, hmin = h/64.0, hmax = h*64.0;

	/*regular computation*/
	while (r<1e-4) {
		if (RK_plan_db.K_roll (EoS, K, h, r, &X, ref)) break;

		err = RK_plan_db.X_join_ec (&nextX, h, K);

		if (errTolerance < 1e-8){
			X = nextX;
			r+=h;
			fprintf(fp, "%.4e, %.4e, %.4e\n", r*c*1e-3, err, h*c);
		}else{
			if (err > errTolerance && h > hmin){
				h *= 0.5;
				continue;
			}
			X = nextX;
			r += h;
			fprintf(fp, "%.4e, %.4e, %.4e\n", r*c*1e-3, err, h*c);
			h *= (err < 0.25*errTolerance && h < hmax)+1;
		}
	}

	fclose(fp);

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
	
	double Iout = I/(m*r*r*sqrt(Bg00));
	double Rout = r*c*1e-3;
	double Mout = m*Mscale;
	double k2 = 1.6*bt*bt*bt*bt*bt*(1-2*bt)*(1-2*bt)*(2-y+2*bt*(y-1))/(2*bt*(6-3*y+3*bt*(5*y-8))+4*bt*bt*bt*(13-11*y+bt*(3*y-2)+2*bt*bt*(1+y))+3*(1-2*bt)*(1-2*bt)*(2-y+2*bt*(y-1))*log(1-2*bt));
	double Lambda=9495*k2*pow(Rout/10,5)/pow(Mout,5);

	fprintf(stderr, "M=%.4lf, R=%.4lf, I=%.4lf, Lambda=%.4lf\n",
	Mout,
	Rout,
	Iout,
	Lambda);

	return 0;
}

int main(int argc, char *argv[]){

	//allocate EoS type variable to store your EoS density and pressure info in log10() and ascending format.
	EoS_t myEoS;

	//load file to EoS_t variable
	loadEoS(&myEoS, "EoS_lib/APR.txt");

	double rhoc = 1.00183e+18;

	solveTOV_db ("tovErr_fixed.log", 0, &myEoS, rhoc);

	solveTOV_db ("tovErr_adaptive.log", 1e-6, &myEoS, rhoc);

	return 0;
}
