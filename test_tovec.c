#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "include/statNS.h"

#define NONE 0
#define SOFT 1
#define HARD 2

static const double pi=3.14159265358979, 
					G=6.673e-11,
					c=2.99792458e8,
					Msun=1.989e30,
					Mscale=2.033931261665867e5;

int solveTOV_db (char dblogPath[], RK_plan_t *RK_plan_db, int ecflag, double errTolerance, EoS_t *EoS, double RhocSI){

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

	double h0 = 10.0/c;

	double h = h0, hmin = h0/64.0, hmax = h0*64.0, hscale;

	int upflag, downflag;

	while (r<1e-4) {

		if (RK_plan_db->K_roll (EoS, K, h, r, &X, ref)) break;

		err = RK_plan_db->X_join_ec (&nextX, h, K);

		upflag = (err < 0.01*errTolerance && h < hmax);

		downflag = (err > errTolerance && h > hmin);

		if(ecflag == HARD && downflag){
			/*forced roll back when error control strategy is HARD while downflag is set*/
			nextX = X;
		}else{
			X = nextX;
			r+=h;
			fprintf(fp, "%.4e, %.4e, %.4e\n", r*c*1e-3, err, h*c);
		}

		/*  0 .. 1 -> 1 .. a  */
		/* (ecflag != NONE) * (a-1) + 1 */
		hscale = (ecflag != NONE) * ((upflag+1)*(1-0.5*downflag) - 1) + 1;
		
		/*scale the step size according to ecflag and err*/
		h *= hscale;
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

	/*run with RK5M debug mode*/
	RK_plan_t rk5M_db = {RungeKutta_RK5M_roll, NULL, RungeKutta_RK5M_join_ec};
	solveTOV_db ("tovErr_rk5m_fixed.log", &rk5M_db, NONE, 0, &myEoS, rhoc);
	solveTOV_db ("tovErr_rk5m_soft.log", &rk5M_db, SOFT, 1e-8, &myEoS, rhoc);
	solveTOV_db ("tovErr_rk5m_hard.log", &rk5M_db, HARD, 1e-8, &myEoS, rhoc);

	/*run with RK5DP debug mode*/
	RK_plan_t rk5DP_db = {RungeKutta_RK5DP_roll, NULL, RungeKutta_RK5DP_join_ec};
	solveTOV_db ("tovErr_rk5dp_fixed.log", &rk5DP_db, NONE, 0, &myEoS, rhoc);
	solveTOV_db ("tovErr_rk5dp_soft.log", &rk5DP_db, SOFT, 1e-8, &myEoS, rhoc);
	solveTOV_db ("tovErr_rk5dp_hard.log", &rk5DP_db, HARD, 1e-8, &myEoS, rhoc);

	/*run with RK5F debug mode*/
	RK_plan_t rk5F_db = {RungeKutta_RK5F_roll, NULL, RungeKutta_RK5F_join_ec};
	solveTOV_db ("tovErr_rk5f_fixed.log", &rk5F_db, NONE, 0, &myEoS, rhoc);
	solveTOV_db ("tovErr_rk5f_soft.log", &rk5F_db, SOFT, 1e-8, &myEoS, rhoc);
	solveTOV_db ("tovErr_rk5f_hard.log", &rk5F_db, HARD, 1e-8, &myEoS, rhoc);

	/*run with RK5CP debug mode*/
	RK_plan_t rk5CK_db = {RungeKutta_RK5CK_roll, NULL, RungeKutta_RK5CK_join_ec};
	solveTOV_db ("tovErr_rk5ck_fixed.log", &rk5CK_db, NONE, 0, &myEoS, rhoc);
	solveTOV_db ("tovErr_rk5ck_soft.log", &rk5CK_db, SOFT, 1e-8, &myEoS, rhoc);
	solveTOV_db ("tovErr_rk5ck_hard.log", &rk5CK_db, HARD, 1e-8, &myEoS, rhoc);

	/*run with RK4M debug mode*/
	RK_plan_t rk4M_db = {RungeKutta_RK4M_roll, NULL, RungeKutta_RK4M_join_ec};
	solveTOV_db ("tovErr_rk4m_fixed.log", &rk4M_db, NONE, 0, &myEoS, rhoc);
	solveTOV_db ("tovErr_rk4m_soft.log", &rk4M_db, SOFT, 1e-8, &myEoS, rhoc);
	solveTOV_db ("tovErr_rk4m_hard.log", &rk4M_db, HARD, 1e-8, &myEoS, rhoc);

	return 0;
}
