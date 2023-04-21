/*
*	Neutron Star computation code to derive:
*       mass
*       radius
*       moment of inertia
*       tidal deformability
*       f-mode frequency
*
*	version Cleveland
*
*	this code is tested on CentOS-8 with gcc 10.3.1
*
*	Windows compatible code is NOT included in this version
*
*	Please contact: c.houyuan@mail.scut.edu.cn
*
*	or you can leave your comment on my website:
*	https://zoo-404.cn/
*
*/

/*
 * AmEoS decleration:
 * 
 * This is a series of asymmetric nuclear matters EoS
 * Origin code is developed by Nai-Bo Zhang and Bao-An in f90
 *
 * This is a testing-purposed C/C++ lib transplated by Houyuan Chen
 *
 * Given four asymmetric energy parameters as inputs,
 * the transition density is automatically considered. 
 *
 * User of this code should quote the publication:  
 * Combined constraints on the equation of state of dense neutron-rich matter from terrestrial nuclear experiments and observations of neutron stars
 * Nai-Bo Zhang, Bao-An Li, and Jun Xu, APJ (2018) 589, 90, arXiv:1801.06855.
 *
*/


#ifndef _STATNS_H
#define _STATNS_H
#endif

#define EOS_MAX_LENGTH 2048
#define ARXIV_MAX_LENGTH 2048

#ifndef _AMEOS_H
#define _AMEOS_H
#endif
#define PT_LENGTH 5166

struct statNS_EoS_t {
	int length;
	double RhomaxSI;
	double RhominSI;
	double Rhomax;
	double Rhomin;
	double Pmax;
	double Pmin;
	double Mmax;
	double Rhoc_MmaxSI;
	double lgRho[EOS_MAX_LENGTH];
	double lgP[EOS_MAX_LENGTH];
	double lgRho_SI[EOS_MAX_LENGTH];
	double lgP_SI[EOS_MAX_LENGTH];
	char FilePath[128];
	double lgDen[EOS_MAX_LENGTH];
	double lgDen_SI[EOS_MAX_LENGTH];
	double XL;
	double Ksym;
	double Jsym;
	double J0;
};
typedef struct statNS_EoS_t 	EoS_t;

struct Am_EoS_opt_t {
	int secure;
	int causality;
	int total_length;
	int core_length;
	double dU;
	double maxU;
};
typedef struct Am_EoS_opt_t		EoS_opt_t;

struct statNS_Path_t {
	char FilePath[128];
};
typedef struct statNS_Path_t	Path_t;

struct statNS_Asset_t {
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
typedef struct statNS_Asset_t	CompactStar_t;

struct statNS_RK_Arr_t {
	double P;
	double M;
	double I;
	double Ag00;
	double y;
};
typedef struct statNS_RK_Arr_t	RK_Arr_t;

struct solveTOVplan_t{
	int section;
	int length;
	double *RhocSI;
	EoS_t *EoS;
	CompactStar_t *Results;
};
typedef struct solveTOVplan_t solveTOVplan_t;

struct getMplan_t {
	int section;
	int length;
	double *RhocSI;
	EoS_t *EoS;
	double *massArr;
	double (*getM)(EoS_t *EoS, double RhocSI);
};
typedef struct getMplan_t getMplan_t;

struct statNS_auxSpace_t {
	double *mfile,*pfile,*rhofile,*Bfile,*Bcor,*Wfile,*Wfile1,*Wfile2,*Vfile,*Vfile1,*Vfile2;
};
typedef struct statNS_auxSpace_t	auxSpace_t;

struct fmodeplan_t {
	int section;
	int length;
    double *RhocSI;
	EoS_t *EoS;
    CompactStar_t *Results;
    auxSpace_t *auxSpace;
	int workerId;
};
typedef struct fmodeplan_t fmodeplan_t;

struct statNS_RK_plan_t {
	int (*K_roll)(EoS_t *EoS, RK_Arr_t *K, double h, double r, RK_Arr_t *X, int *ref);
	int (*X_join)(RK_Arr_t *Result, double h, RK_Arr_t *K);
	double (*X_join_ec)(RK_Arr_t *Result, double h, RK_Arr_t *K);
};
typedef struct statNS_RK_plan_t	RK_plan_t;

struct statNS_ArXiv_t {
	int length;
	int index[ARXIV_MAX_LENGTH];
	double RhocSI[ARXIV_MAX_LENGTH];
	double M[ARXIV_MAX_LENGTH];
};

typedef struct statNS_ArXiv_t	ArXiv_t;


/*Low level function*/

int dFunc(EoS_t *EoS, RK_Arr_t *K, double r, RK_Arr_t *X, int *ref);
int dFunc_s(EoS_t *EoS, RK_Arr_t *K, double r, RK_Arr_t *X, int *ref);

double interp_p2rho(EoS_t *EoS, double cp, double *VsOUT, int *ref);
double interp_rho2p(EoS_t *EoS, double crho);
double interp_p2rho_SI(EoS_t *EoS, double cp, int *ref);
double interp_rho2p_SI(EoS_t *EoS, double crho);

int RungeKutta_Array_add (RK_Arr_t *Result, double h, RK_Arr_t *K, RK_Arr_t *X);
int RungeKutta_Array_adds (RK_Arr_t *Result, double *H, RK_Arr_t *K, RK_Arr_t *X, int dim);
int RungeKutta_Array_adds_s (RK_Arr_t *Result, double *H, RK_Arr_t *K, RK_Arr_t *X, int dim);

int loadEoS(EoS_t *EoS, char *Path);

int RungeKutta_RK4_roll(EoS_t *EoS, RK_Arr_t *K, double h, double r, RK_Arr_t *X, int *ref);
int RungeKutta_RK4_join (RK_Arr_t *Result, double h, RK_Arr_t *K);
int RungeKutta_RK4M_roll(EoS_t *EoS, RK_Arr_t *K, double h, double r, RK_Arr_t *X, int *ref);
int RungeKutta_RK4M_join (RK_Arr_t *Result, double h, RK_Arr_t *K);
int RungeKutta_RK5F_roll(EoS_t *EoS, RK_Arr_t *K, double h, double r, RK_Arr_t *X, int *ref);
int RungeKutta_RK5F_join (RK_Arr_t *Result, double h, RK_Arr_t *K);
int RungeKutta_RK5M_roll(EoS_t *EoS, RK_Arr_t *K, double h, double r, RK_Arr_t *X, int *ref);
int RungeKutta_RK5M_join (RK_Arr_t *Result, double h, RK_Arr_t *K);

/*Postscript with _ec mean join Runge Kutta array with error control*/
/*Test features only build for the Numeric Computation Course in SYSU*/
double RungeKutta_RK4M_join_ec (RK_Arr_t *Result, double h, RK_Arr_t *K);
double RungeKutta_RK5F_join_ec (RK_Arr_t *Result, double h, RK_Arr_t *K);
double RungeKutta_RK5M_join_ec (RK_Arr_t *Result, double h, RK_Arr_t *K);

/*Default choice for RK_plan_t is the RK5L algorithm*/
/*Notice that it is faster, but it simply has no error control capability*/
int RungeKutta_RK5L_roll(EoS_t *EoS, RK_Arr_t *K, double h, double r, RK_Arr_t *X, int *ref);
int RungeKutta_RK5L_join (RK_Arr_t *Result, double h, RK_Arr_t *K);

/*Postscript with _s mean RK_plan_t call dFunc_s() instead of dFunc()*/
/*also, join_s() calls RK_Array_adds_s() instead of RK_Array_adds()*/
/*dFunc_s() and join_s() would only calculate P and M, thus faster*/
int RungeKutta_RK5L_roll_s (EoS_t *EoS, RK_Arr_t *K, double h, double r, RK_Arr_t *X, int *ref);
int RungeKutta_RK5L_join_s (RK_Arr_t *Result, double h, RK_Arr_t *K);

/*fmode dedicated function*/
/*strongly not recommend to modify them!!*/
double fp(double r,double p,double e,double m);
double fm(double r,double e);
double Bf(double r,double p,double m,double B);
double DH1(double r,double m,double A,double p,double e,double H1,double H0,double K,double V);
double DK(double r,double H0,double H1,double Dv,double K,double e,double p,double A,double W);
double DW(double r,double W,double A,double gamma,double p,double B,double X,double V,double H0,double K);
double DX(double r,double X,double e,double p,double B,double Dv,double H0,double w,double H1,double K,double V,double A,double F,double W);
double H0f(double r,double B,double X,double m,double p,double A,double H1,double K,double w);
double Vf(double r,double w,double e,double p,double B,double A,double Dp,double W,double H0,double X);
double gf(EoS_t *EoS, double e);
double Ff(double r,double A,double e,double p,double m,double Dp);

/*End of low level function*/


/*High level function*/

double getMmax(EoS_t *EoS);
double M2Rhoc(EoS_t *EoS, double fM);
void *solveTOVworker (void *args);
int solveTOV(CompactStar_t *Results, EoS_t *EoS, double RhocSI);
int solveTOV_mt(CompactStar_t *Results, EoS_t *EoS, double *RhocSI, int length, int threads);

void * fmode_worker (void *args);
int fmode(CompactStar_t *Results, EoS_t *EoS, double RhocSI, auxSpace_t *auxSpace);
int fmode_mt(CompactStar_t *Results, EoS_t *EoS, double *RhocSI, int length, int threads);

double getM_fm(EoS_t *EoS, double RhocSI);
double getMmax_fm(EoS_t *EoS);
double M2Rhoc_fm(EoS_t *EoS, double fM);

/*
 * new frame to fast convert multiple-mass (array) to Rhoc, with one EoS loaded.
*/
void setIndex(int *index, int len);
void rec2arxiv(ArXiv_t *arxiv, double RhocSI, double M);
int arcSimSort(int head, int tail, int* index, double* data);
int interp_ArXiv_M2Rhoc_arr(double *RhocGuess, ArXiv_t *Arxiv, double *massArr, int length);
int getM_fm_mt(double *massArr, EoS_t *EoS, double *RhocSI, int length, int threads);
int M2Rhoc_Arr_fm(double *RhocSI, EoS_t *EoS, double *massArr, int length, int threads);

double getM_s(EoS_t *EoS, double RhocSI);
void *getM_worker (void *args);
int getM_s_mt(double *massArr, EoS_t *EoS, double *RhocSI, int length, int threads);
int M2Rhoc_Arr_s(double *RhocSI, EoS_t *EoS, double *massArr, int length, int threads);

/*
 * amEoS functionality
*/
void set_EoS_default_opt(EoS_opt_t *opt);
int genAmEoS(EoS_t *eos, double XL, double Ksym, double Jsym, double J0, EoS_opt_t *opt);
int saveEoS(EoS_t *eos, char path[]);
