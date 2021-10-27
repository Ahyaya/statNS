/*	Neutron Star computation code to derive:
*       mass
*       radius
*       moment of inertia
*       tidal deformability
*       f-mode frequency
*
*	remix version, please notice that the origin version has been completely out of date.
*
*	this code is tested on CentOS-8 with gcc 10.3.1
*
*	Windows compatible code is NOT included in this version
*
*	Contact: c.houyuan@mail.scut.edu.cn
*
*	or you can leave your comment on my website:
*	https://arxiv.cloud/
*
*/

#ifndef _STATNS_H
#define _STATNS_H
#endif

struct EoS_t {
	int length;
	double RhomaxSI;
	double RhominSI;
	double Rhomax;
	double Rhomin;
	double Pmax;
	double Pmin;
	double Mmax;
	double Rhoc_MmaxSI;
	double lgRho[3000];
	double lgP[3000];
	double lgRho_SI[3000];
	double lgP_SI[3000];
	char FilePath[128];
};

struct CompactStar_t {
	double P;
	double r;
	double Rho;
	double M;
	double Ma;
	double Mp;
	double I;
	double Ag00;
	double y;
	double k2;
	double Lambda;
	double Vs;
	double freq;
	double dampTime;
};

struct RungeKutta_Array_t {
	double P;
	double M;
	double I;
	double Ag00;
	double y;
};

struct solveTOVparams_t{
	double RhocSI;
	struct CompactStar_t *Results;
};

struct ComputeStatus_t {
	int RkStop;
	int AllStop;
	double RkErr;
	int (*K_roll)(struct RungeKutta_Array_t *K, double h, double r, struct RungeKutta_Array_t *X);
	int (*X_join)(struct RungeKutta_Array_t *Result, double h, struct RungeKutta_Array_t *K);
};

struct auxSpace_t {
	double *mfile,*pfile,*rhofile,*Bfile,*Bcor,*Wfile,*Wfile1,*Wfile2,*Vfile,*Vfile1,*Vfile2;
};

struct fmodeParams_t {
	double RhocSI;
	struct CompactStar_t *Results;
	struct auxSpace_t *auxSpace;
};

int dFunc(struct RungeKutta_Array_t *K, double r, struct RungeKutta_Array_t *X);
double interp_p2rho(double cp, double *VsOUT);
double interp_rho2p(double crho);
int RungeKutta_Array_add (struct RungeKutta_Array_t *Result, double h, struct RungeKutta_Array_t *K, struct RungeKutta_Array_t *X);
int RungeKutta_Array_adds (struct RungeKutta_Array_t *Result, double *h, struct RungeKutta_Array_t *K, struct RungeKutta_Array_t *X, int dim);
int loadEoS(char *Path);
int RungeKutta_RK4_roll(struct RungeKutta_Array_t *K, double h, double r, struct RungeKutta_Array_t *X);
int RungeKutta_RK4_join (struct RungeKutta_Array_t *Result, double h, struct RungeKutta_Array_t *K);
int RungeKutta_RK4M_roll(struct RungeKutta_Array_t *K, double h, double r, struct RungeKutta_Array_t *X);
int RungeKutta_RK4M_join (struct RungeKutta_Array_t *Result, double h, struct RungeKutta_Array_t *K);
int RungeKutta_RK5F_roll(struct RungeKutta_Array_t *K, double h, double r, struct RungeKutta_Array_t *X);
int RungeKutta_RK5F_join (struct RungeKutta_Array_t *Result, double h, struct RungeKutta_Array_t *K);
int RungeKutta_RK5M_roll(struct RungeKutta_Array_t *K, double h, double r, struct RungeKutta_Array_t *X);
int RungeKutta_RK5M_join (struct RungeKutta_Array_t *Result, double h, struct RungeKutta_Array_t *K);
int RungeKutta_RK5L_roll(struct RungeKutta_Array_t *K, double h, double r, struct RungeKutta_Array_t *X);
int RungeKutta_RK5L_join (struct RungeKutta_Array_t *Result, double h, struct RungeKutta_Array_t *K);

double getMmax();
double M2Rhoc(double fM);
void *solveTOVGate (void *args);
int solveTOV(double RhocSI, struct CompactStar_t *Results);
int solveTOV_mt(int threadNum, struct CompactStar_t *Results, double *cdenArray, int arrayLen);

double ch(double cp);
double che(double crho);
double gf(double e);
double Ff(double r,double A,double e,double p,double m,double Dp);

void * fmodeGate (void *args);
int fmode(double RhocSI, struct CompactStar_t *Results, struct auxSpace_t *auxSpace);
int fmode_mt(int threadNum, struct CompactStar_t *Results, double *cdenArray, int arrayLen);

double fp(double r,double p,double e,double m);
double fm(double r,double e);
double Bf(double r,double p,double m,double B);
double DH1(double r,double m,double A,double p,double e,double H1,double H0,double K,double V);
double DK(double r,double H0,double H1,double Dv,double K,double e,double p,double A,double W);
double DW(double r,double W,double A,double gamma,double p,double B,double X,double V,double H0,double K);
double DX(double r,double X,double e,double p,double B,double Dv,double H0,double w,double H1,double K,double V,double A,double F,double W);
double H0f(double r,double B,double X,double m,double p,double A,double H1,double K,double w);
double Vf(double r,double w,double e,double p,double B,double A,double Dp,double W,double H0,double X);
