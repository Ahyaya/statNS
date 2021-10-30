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
	double lgRho[3000];
	double lgP[3000];
	double lgRho_SI[3000];
	double lgP_SI[3000];
	char FilePath[128];
};
typedef struct statNS_EoS_t 	EoS_t;

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

struct statNS_RK_Arr_s {
	double P;
	double M;
};
typedef struct statNS_RK_Arr_s	RK_Arr_s;

struct statNS_RK_Arr_t {
	double P;
	double M;
	double I;
	double Ag00;
	double y;
};
typedef struct statNS_RK_Arr_t	RK_Arr_t;

struct solveTOVparams_t{
	double RhocSI;
	CompactStar_t *Results;
	EoS_t *EoS;
};

struct statNS_ComputeStatus_t {
	int RkStop;
	int AllStop;
	double RkErr;
	int (*K_roll)(EoS_t *EoS, RK_Arr_t *K, double h, double r, RK_Arr_t *X);
	int (*X_join)(RK_Arr_t *Result, double h, RK_Arr_t *K);
	int (*K_roll_s)(EoS_t *EoS, RK_Arr_s *K, double h, double r, RK_Arr_s *X);
	int (*X_join_s)(RK_Arr_s *Result, double h, RK_Arr_s *K);
};
typedef struct statNS_ComputeStatus_t	ComputeStatus_t;

struct statNS_auxSpace_t {
	double *mfile,*pfile,*rhofile,*Bfile,*Bcor,*Wfile,*Wfile1,*Wfile2,*Vfile,*Vfile1,*Vfile2;
};
typedef struct statNS_auxSpace_t	auxSpace_t;

struct fmodeParams_t {
    double RhocSI;
    double fM;
    CompactStar_t *Results;
    auxSpace_t *auxSpace;
    EoS_t *EoS;
};

struct statNS_ArXiv_t {
	int length;
	int index[2048];
	double RhocSI[2048];
	double M[2048];
};

typedef struct statNS_ArXiv_t	ArXiv_t;

struct getM_params_t {
	EoS_t *EoS;
	double RhocSI;
	double *M;
};

int dFunc(EoS_t *EoS, RK_Arr_t *K, double r, RK_Arr_t *X);
double interp_p2rho(EoS_t *EoS, double cp, double *VsOUT);
double interp_rho2p(EoS_t *EoS, double crho);
double interp_p2rho_SI(EoS_t *EoS, double cp);
double interp_rho2p_SI(EoS_t *EoS, double crho);
int RungeKutta_Array_add (RK_Arr_t *Result, double h, RK_Arr_t *K, RK_Arr_t *X);
int RungeKutta_Array_adds (RK_Arr_t *Result, double *h, RK_Arr_t *K, RK_Arr_t *X, int dim);
int loadEoS(EoS_t *EoS, char *Path);
int RungeKutta_RK4_roll(EoS_t *EoS, RK_Arr_t *K, double h, double r, RK_Arr_t *X);
int RungeKutta_RK4_join (RK_Arr_t *Result, double h, RK_Arr_t *K);
int RungeKutta_RK4M_roll(EoS_t *EoS, RK_Arr_t *K, double h, double r, RK_Arr_t *X);
int RungeKutta_RK4M_join (RK_Arr_t *Result, double h, RK_Arr_t *K);
int RungeKutta_RK5F_roll(EoS_t *EoS, RK_Arr_t *K, double h, double r, RK_Arr_t *X);
int RungeKutta_RK5F_join (RK_Arr_t *Result, double h, RK_Arr_t *K);
int RungeKutta_RK5M_roll(EoS_t *EoS, RK_Arr_t *K, double h, double r, RK_Arr_t *X);
int RungeKutta_RK5M_join (RK_Arr_t *Result, double h, RK_Arr_t *K);
int RungeKutta_RK5L_roll(EoS_t *EoS, RK_Arr_t *K, double h, double r, RK_Arr_t *X);
int RungeKutta_RK5L_join (RK_Arr_t *Result, double h, RK_Arr_t *K);

int dFunc_s(EoS_t *EoS, RK_Arr_s *K, double r, RK_Arr_s *X);
int RungeKutta_Array_adds_s (RK_Arr_s *Result, double *h, RK_Arr_s *K, RK_Arr_s *X, int dim);
int RungeKutta_RK5L_roll_s (EoS_t *EoS, RK_Arr_s *K, double h, double r, RK_Arr_s *X);
int RungeKutta_RK5L_join_s (RK_Arr_s *Result, double h, RK_Arr_s *K);

double getMmax(EoS_t *EoS);
double M2Rhoc(EoS_t *EoS, double fM);
void *solveTOVGate (void *args);
int solveTOV(CompactStar_t *Results, EoS_t *EoS, double RhocSI);
int solveTOV_mt(CompactStar_t *Results, EoS_t *EoS, double *RhocSI, int arrayLen, int threads);

double gf(EoS_t *EoS, double e);
double Ff(double r,double A,double e,double p,double m,double Dp);

void * fmodeGate (void *args);
int fmode(CompactStar_t *Results, EoS_t *EoS, double RhocSI, auxSpace_t *auxSpace);
int fmode_mt(CompactStar_t *Results, EoS_t *EoS, double *RhocSI, int arrayLen, int threads);

double fp(double r,double p,double e,double m);
double fm(double r,double e);
double Bf(double r,double p,double m,double B);
double DH1(double r,double m,double A,double p,double e,double H1,double H0,double K,double V);
double DK(double r,double H0,double H1,double Dv,double K,double e,double p,double A,double W);
double DW(double r,double W,double A,double gamma,double p,double B,double X,double V,double H0,double K);
double DX(double r,double X,double e,double p,double B,double Dv,double H0,double w,double H1,double K,double V,double A,double F,double W);
double H0f(double r,double B,double X,double m,double p,double A,double H1,double K,double w);
double Vf(double r,double w,double e,double p,double B,double A,double Dp,double W,double H0,double X);

/*
 * Function for fm part
 * EoS is no longer given as global var in this part
 * all function need to specify address of EoS instead.
 * This frame is optimized for multi-EoS computation.
 */

double getM_fm(EoS_t *EoS, double RhocSI);
double getMmax_fm(EoS_t *EoS);
double M2Rhoc_fm(EoS_t *EoS, double fM);
void * m2fmodeGate(void *args);
int m2fmode_mt(CompactStar_t *Results, double fM, Path_t *EoSlist, int EoSlistLen, int threads);

/*
 * new frame to fast convert multiple-mass (array) to Rhoc, with one EoS loaded.
*/
void setIndex(int *index, int len);
void rec2arxiv(ArXiv_t *arxiv, double RhocSI, double M);
int arcSimSort(int head, int tail, int* index, double* data);
int interp_ArXiv_M2Rhoc_arr(double *RhocGuess, ArXiv_t *Arxiv, double *massArr, int arrayLen);
void *getM_fmGate (void *args);
int getM_fm_mt(double *massArr, EoS_t *EoS, double *RhocSI, int arrayLen, int threads);
int M2Rhoc_Arr_fm(double *RhocSI, EoS_t *EoS, double *massArr, int arrayLen, int threads);

double getM_s(EoS_t *EoS, double RhocSI);
void *getM_sGate (void *args);
int getM_s_mt(double *massArr, EoS_t *EoS, double *RhocSI, int arrayLen, int threads);
int M2Rhoc_Arr_s(double *RhocSI, EoS_t *EoS, double *massArr, int arrayLen, int threads);