//#define C_VERSION 0          /* for Borland */
  #define C_VERSION 1          /* for Visual */

#include "..\headers\version.h"
#include "..\headers\mylib.h"
#include "..\headers\matlib.h"

#define  PI  3.14159265358979
#define EPS  1E-40

#define can X1  parm X2  x X3  fc X4 f X5  xmin X6  xmax X7  multX X8 \
  multF X9  lam X10  mzamp X11  ko X12  komod X13  teta X14 \
  n X15  m X16  me X17  name X18  extin X19  extout X20
#define DECLARE0(can,parm,x,fc,f,xmin,xmax,multX,multF,lam,mzamp,\
  ko,komod,teta,n,m,me,name,extin,extout) \
  static double *x,*f,*xmin,*xmax,*multX,*multF,*lam,*ko,*komod,*teta; \
  static double fc; \
  static int *mzamp; static int n,m,me, can,parm;\
  char name[15],extin[3],extout[3]
#undef can parm x fc f xmin xmax multX multF lam mzamp ko komod teta \
  n m me name extin extout
/* for recursive optimization */

#define can X1  parm X2 x X3  fc X4 f X5  xmin X6  xmax X7  multX X8 \
  multF X9  lam X10  mzamp X11  ko X12  komod X13  teta X14 \
  n X15  m X16  me X17
#define DECLARE1(can,parm,x,fc,f,xmin,xmax,multX,multF,lam,mzamp,\
  ko,komod,teta,n,m,me) \
  static double *x,fc,*f,*xmin,*xmax,*multX,*multF,*lam,*ko,*komod,*teta;\
  static int *mzamp; static int n,m,me; int can,parm
#undef can parm x fc f xmin xmax multX multF lam mzamp ko komod teta n m me
/* for recursive optimization */

#define  can_ X1  parm_ X2  x_ X3  fc_ X4  f_ X5  xmin_ X6  xmax_ X7 \
  multX_ X8 multF_ X9  lam_ X10  mzamp_ X11  ko_ X12  komod_ X13  teta_ X14 \
  n_ X15  m_ X16  me_ X17  x0_ X18  x__ X19  name_ X20 \
  extin_ X21  extout_ X22  prmod_ X23  iprmod_ X24
  #define DECLARE2(can_,parm_,x_,fc_,f_,xmin_,xmax_,multX_,multF_,lam_,mzamp_,\
  ko_,komod_,teta_,n_,m_,me_,x0_,x__,name_,extin_,extout_,prmod_,iprmod_) \
  static double *x_,*f_,*xmin_,*xmax_,*multX_,*multF_,*lam_,*ko_,\
  *komod_,*teta_,*x0_,*x__; \
  static double fc_; \
  static int *mzamp_; static int n_,m_,me_,can_,parm_=2,prmod_,iprmod_;\
  char name_[15],extin_[3],extout_[3]
#undef can_ parm_ x_ fc_ f_ xmin_ xmax_ multX_ multF_ lam_ mzamp_ ko_ \
  komod_ teta_ n_ m_ me_ x0_ x__ name_ extin_ extout_ \
  prmod_ iprmod_
/* for recursive optimization */

#define can X1  parm X2  parmget X3  x X4  fc X5  f X6  xmin X7  xmax X8 \
  multX X9 \ multF X10  lam X11  mzamp X12  ko X13  komod X14  teta X15 \
  n X16  m X17  me X18  namef X19  name X20  extin X21  extout X22 \
  npar X23  par X24
#define CALLROPUD0(can,parm,parmget,x,fc,f,xmin,xmax,multX,multF,lam,mzamp,\
  ko,komod,teta,n,m,me,namef,name,extin,extout,npar,par) \
  callropud0(can,parm,parmget,&x,&fc,&f,&xmin,&xmax,&multX,&multF,&lam, \
    &mzamp,&ko,&komod,&teta,&n,&m,&me,namef,name,extin,extout,npar,par);
#undef can parm parmget x fc f xmin xmax multX multF lam mzamp ko komod \
  teta n m me namef name extin extout npar par
#define  X1 can  X2 parm  parmget X3  x X4  fc X5  f X6  xmin X7  xmax X8 \
  multX X9 multF X10  lam X11  mzamp X12  ko X13  komod X14  teta X15 \
  n X16  m X17  me X18
#define CALLROPUD1(can,parm,parmget,x,fc,f,xmin,xmax,multX,multF,lam,mzamp,\
  ko,komod,teta,n,m,me) \
  callropud1(can,parm,parmget,&x,&fc,&f,&xmin,&xmax,&multX,&multF,&lam,\
  &mzamp,&ko,&komod,&teta,&n,&m,&me)
#undef can parm parmget x fc f xmin xmax multX multF lam mzamp ko komod\
   teta n m me parmget

#define can_ X1  parm_ X2  x_ X3  fc_ X4  f_ X5  xmin_ X6  xmax_ X7 \
  multX_ X8  multF_ X9  lam_ X10  mzamp_ X11  ko_ X12  komod_ X13  teta_ X14 \
  n_ X15  m_ X16  me_ X17  ko X18  komod X19  teta X20  x0_ X21  x__ X22 \
  name_ X23 extin_ X24  extout_ X25  prmod_ X26  iprmod_ X27  koderpa X28
#define CALLROPUD2(can_,parm_,x_,fc_,f_,xmin_,xmax_,multX_,multF_,lam_, \
  mzamp_,ko_,komod_,teta_,n_,m_,me_,ko,komod,teta,x0_,x__,name_, \
  extin_,extout_, prmod_,iprmod_,kocallrpd) \
  callropud2(can_,&parm_,&x_,&fc_,&f_,&xmin_,&xmax_,&multX_,&multF_,&lam_, \
    &mzamp_,&ko_,&komod_,&teta_,&n_,&m_,&me_,ko,komod,teta,&x0_,&x__,name_, \
    extin_,extout_,prmod_,iprmod_,kocallrpd)
#undef can_ parm_ x_ fc_ f_ xmin_ xmax_ multX_ multF_ lam_ mzamp_ ko_ \
  komod_ teta_ n_ m_ me_ ko komod teta x0_ x__ name_ extin_ extout_ \
  prmod_ iprmod_ koderpa
/* for recursive optimization */

#define g X1  dfdx X2  ko X3  komod X4  teta X5 \
  x X6  x_ X7  x__ X8  n_ X9  prmod2 X10  iprmod2 X11  koderpa X12
#define DERPA(g,dfdx,ko,komod,teta,x,x_,x__,n_,prmod2,iprmod2,koderpa) \
  derpa(g,dfdx,ko,komod,teta,x,x_,x__,n_,prmod2,iprmod2,koderpa)
#undef g dfdx ko komod teta x x_ x__ n_ prmod2 iprmod2 koderpa
/* for recursive optimization */

void callropud0(int can,int parm,int parmget,double **x,double *pfc,
    double **f,double **xmin,double **xmax,double **multX,double **multF,
    double **lam,int **mzamp,double **ko,double **komod,double **teta,
    int *pn,int *pm,int *pme,char namef[],char name[],char extin[],
    char extout[],int npar,char **par);
void callropud1(int can,int parm,int parmget,double **x,double *pfc,
    double **f,double **xmin,double **xmax,double **multX,double **multF,
    double **lam,int **mzamp,double **ko,double **komod,double **teta,
    int *pn,int *pm,int *pme);
void callropud2(int can_,int *parm_,double **x_,double *pfc_,double **f_,
   double **xmin_,double **xmax_,double **multX_,double **multF_,
   double **lam_,int **mzamp_,double **ko_,double **komod_,double **teta_,
   int *pn_,int *pm_,int *pme_,double ko[],double komod[],double teta[],
   double **x0_,double **x__,char name_[],char extin_[],char extout_[],
   int prmod_,int iprmod_,double kocallrpd[]);
void chngvd(double x[],double z[],int n,int zamp,double xmin[],double xmax[],
	    int mzamp[],double multX[]);
void chngvr(double x[],double z[],int n,int zamp,double xmin[],double xmax[],
	    int mzamp[],double multX[]);
void chngdd(double z[],double g[],int n,int zamp,double xmin[],
	    double xmax[],int mzamp[],double multX[]);
void chngddm(double z[],DOUBLE **dfdz,int m,int n,int zamp,double xmin[],
	    double xmax[],int mzamp[],double multX[],double multF[]);
void chngfd(double f[],int m,double multF[]);
void chngfr(double f[],int m,double multF[]);
void derpa(double g[],DOUBLE **dfdx,double ko[],
     double komod[],double teta[],double x[],double x_[],double x__[],
     int l,int prmod,int iprmod,double kocallrpd[]);
void dirs(double z[],double *pfc,double *pfcg, double f[],double g[],
	  DOUBLE **dfdz,double dfidalf,int n,int m,int me,double multX[],
          double multF[],double xmin[],double xmax[],double lam[],int mzamp[],
          double d[],double *palfa,double ko[],double komod[],double teta[]);
void  getname(char name1[],char name2[],char extin[],char extout[],
              int can,int parm);
void getname0(char namef[],char name[],char extin[],char extout[],
	     int can,int npar,char **par);
void getname2(char name[],char extn[],char extout[],int can);
void getropud(double **x,double **xmin,double **xmax,double **f,
     double **multX,double **multF,double **lam,int **mzamp,double **ko,
     double **komod,double **teta,int *pn,int *pm,int *pme,int *pnmod,
     int *pnteta,int can,int parm);
void getmem(double **x,double **xmin,double **xmax,double **f,
	  double **multX,double **multF,double **lam,int **mzamp,double **ko,
	  double **komod,double **teta,int can,int *nAl,int *mAl);
void getmem1(double **x,double **xmin,double **xmax,double **f,
	  double **multX,double **multF,double **lam,int **mzamp,double **ko,
	  double **komod,double **teta,int can,int *nAl,int *mAl);
void getmem2(double ko[],int nAl,int mAl);
void ropud(double x[],double *pfc,double *pfcg,double xmin[],double xmax[],
           double f[],double multX[],double multF[],double lam[],int mzamp[],
           int n,int m,int me,double ko[],double komod[],double teta[]);
void sqpr(double x[],double *pfc,double *pfcg,double xmin[],double xmax[],
             double f[],double multX[],double multF[],double lam[],int mzamp[],
	     int n,int m,int me,double ko[],double komod[],double teta[]);
void newtonc(double x[],double *pfc,double *pfcg,double xmin[],double xmax[],
             double f[],double multX[],double multF[],double lam[],int mzamp[],
             int n,int m,int me,double ko[],double komod[],double teta[]);
void newtonop(double x[],double *pfc,double *pfcg,double xmin[],double xmax[],
             double f[],double multX[],double multF[],double lam[],int mzamp[],
             int n,int m,int me,double ko[],double komod[],double teta[]);
void functru0(double z[],double *pfc,double *pfcg,double f[],double g[],
			 DOUBLE **dfdz,int n,int m,int me,double multX[],double multF[],
             double xmin[],double xmax[],double lam[],int mzamp[],
             double ko[],double komod[],double teta[]);
void functru1(double z[],double *pfc,double f[],double g[],DOUBLE **dfdz,
             int n,int m,int me,double multX[],double multF[],double xmin[],
             double xmax[],int mzamp[],double ko[],double komod[],
             double teta[]);
void functru2(double x[],double *pfc,double f[],double g[],DOUBLE **dfdx,
              int n,int m,int me,double ko[],double komod[],double teta[]);
void functrr(double x[],double *pfc,double f[],double g[],DOUBLE **dfdx,
             int n,int m,int me,double ko[],double komod[],double teta[]);
void functr(double x[],double *pfc,double f[],double g[],DOUBLE **dfdx,
             int n,int m,int me,double ko[],double komod[],double teta[]);
void funct(double x[],double *pfc,double f[],double g[],DOUBLE **dfdx,
             int n,int m,int me,double ko[],double komod[],double teta[]);
//void fcg_cal(double *pfcg,double ko[],int m,int me);
//void dfcg_cal(double g[],double f[],DOUBLE **dfdz,int n,int m,int me);
//void fcg_lagr(double *pfc,double *pfcg,double f[],double lam[],
//             int m,int me,double ko[]);
//void dfcg_lagr(double g[],double gg[],double f[],DOUBLE **dfdz,
//             double lam[],int m,int me,int n,double ko[]);
//void fcg_sqpr(double *pfc,double *pfcg,double f[],double z[],double xmin[],
//   double xmax[],double multX[],double lam[],int n,int m, int me,double ko[]);
//void gLg2_sqpr(double g[],double gL[],double g2[],DOUBLE **dfdz,double f[],
//   double lam[],int n,int m,int me);
//void matrH_Br(DOUBLE **H,int n,double z[],double zold[],double f[],
//             double fold[],int *pifail,double ko[]);
/*void matrH_N(DOUBLE **H,int n,double *pdet,long *varm1,long *varm2,
	     char path[]);*/
//void matrH_N(DOUBLE **H,int n,double *pdet,int *varm1,int *varm2,
//           char path[]);
//void matrH_BFG(DOUBLE **H,int n,double z[],double zold[],double g[],
//             double gold[],int parm,double gammaH,int *ifail,
//             double *teta,double ko[]);
//void matrB_BFG(DOUBLE **B,int n,double s[],double y[],double gamma,
//             int *pifail,double *pteta,double ko[]);
void print(double x[],double z[],int prz,double *pfc,double *pfcg,double f[],
		  double g[],int prg,double gL[],int prgL,double gg[],int prgg,
		  int n, int m,int me,double ko[],double nevf);
double test_max(double x[],int m,int me,int m0);
void mult(double *fc,double x[],double f[],double g[],DOUBLE **dfdx,
	 double multX[],double multF[],int n,int m,int me,double ko[]);
void ext_parm();
void delmem1(double *x,double *xmin,double *xmax,double *f,
	  double *multX,double *multF,double *lam,int *mzamp,double *ko,
	  double *komod,double *teta);
void delmem2(double ko[]);
void delmem3(double kocallrpd[]);
void check1(DOUBLE **H,DOUBLE ***H1a,int n,int can);
void check2(DOUBLE **H,DOUBLE **H1,double det,int n,int can);
void convIzD(double x[],double ko[],int n);
void convIzR(double x[],double ko[],int n);
void putCan(void);