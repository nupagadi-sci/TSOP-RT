
#include "..\headers\ropud.h"
#define _Nj 11
#define _Nt 30
#define _Nq 5
#define _Nd 2
#define _Nz 2
#define _Na 6
#define EPS_MAXFI 1e-07
#define EPS_UL	1e-02
#define EPS_T 1e-02
#define	K	0.0001
#define KA 1000
#define RMIN 0.01

#define EPS_QCR 1e-06

extern int nfunU,nfungU,nfunM,nfungM,nfunL,nfungL;
// раскомментарить для вывода отладочной информации
//#define DEBUG 0

extern int Nit, Nz, Nq, Nd, Nj, Na, Nt, Nap, Ns, is_linearization, is_b_min_inheritance, is_b_split_inheritance;
extern double z_start[_Nz], zmin[_Nz], zmax[_Nz], **b_start, bmin, bmax, d[_Nd], dmin[_Nd], dmax[_Nd], **E, **A, *P, Rabs[_Nq], alpha[_Na];
extern double start_q[_Nq], start_minq[_Nq], start_maxq[_Nq];
extern double gamma, mltf[5];
extern struct T *ROPUD_pT;
extern struct TK *ROPUD_pTK;
//extern FILE *log;
//extern int ROPUD_ConstraintToMaximize;

struct T
{
	int id;
	//	//	Вектор управляющих параметров 
	//double *z;
	//	//	Матрица коэффициентов апроксимаций управляющих параметров
	//	//	b[Nz][Nq+1]; b[i] - коэффициенты i-го управляющего параметра
	//double **b;	//	b[i][Nq] - свободный член
//границы
	double *lbound;
	double *rbound;
//поисковые значения границ
	double *slbound;
	double *srbound;
	int isSoftCons;	//	Мягкое ли ограничение
//вероятность на область (поисковая)
	//double alpha;
//номер ограничения, которому принадлежит эта область
	int NCons;
//	критические точки, макс. значения неравенств в области, флаг активности
	double *Qcr;
	int NQcr;
	//double *cv;
	int isActive;
	//int *ActiveConstraints;
	int isLocked;

	//int *fork;	//	Ветвление
		//	Количество лишних дроблений огр-й. Обл-ти со зн-ем атр-та ">0" при плановом
	//int NumExtraSplits;	//	дроблении по критерию - игнорируются, а атрибут уменьшается на 1
//	размер области и номер макс. ребра	
	double r;
	int qr;
//	коэффициент для вычисляемого в области критерия (части критерия при сложении в общий критерий)
//	у областей для ограничений =0 и ни на что не влияет
	struct TK *relK;	//	Область критерия, которой принадлежит данная область ограничения
	struct T *pnext;
};

struct TK
{
	int id;
		//	Вектор управляющих параметров 
	double *z;
		//	Матрица коэффициентов апроксимаций управляющих параметров
		//	b[Nz][Nq+1]; b[i] - коэффициенты i-го управляющего параметра
	double **b;	//	b[i][Nq] - свободный член
		//	Границы
	double *lbound;
	double *rbound;

	int iE[_Nq];	//	Индекс матожидания в глобальном массиве E
	double E[_Nq];

	int isActive;
		//	Размер области и номер макс. ребра	
	double r;
	int qr;
		//	Коэффициент для вычисляемого в области критерия (части критерия при сложении в общий критерий)
		//	у областей для ограничений =0 и ни на что не влияет
	double a;
	double ai[_Nq];
	struct TK *pnext;
};

struct BinTree
{
	struct T *region;
	struct BinTree *lBranch;
	struct BinTree *rBranch;
};

struct TSplit
{
	double *lbound;
	double *rbound;
	int iteration;
	struct TSplit *pnext;
};

	//	Переназначает всем ячейкам дерева новую соответствующую область pTK по критерию 
void replant(struct BinTree *BT, struct TK *pTK);
	//	Пересчитывает управления z, для апроксимаций b и точек q
int set_z(struct TK *pT, double* q);
	//	Получает ссылку на значение id исходной обласити, присваивает ему значение левой части этой области,
int getNewID(int *lid);		//	возвращает значение правой области

struct TSplit* CreateTS(double *lbound, double *rbound, int Nq, int i);

struct T* CreateT(int Nq, int Nj, double *lbound, double *rbound, double *Rabs);
struct TK* CreateTK(double *lbound, double *rbound);
void DeleteT(struct T **pT);
void CalculateSize(struct T *pT,double *Rabs);
void CalculateSize2(struct T *pT);
struct T* Split(struct T *pT1, double *Rabs);
int CreateNewQcr(int NQcr, double *pQcr, struct T *pT);
void DeleteList(struct T **proot);
int AssignConstraintsValues(int t, int Nj, double *cv, struct T *pT);
int AddQcr(struct T *pT, double *pQcr);
int DeleteQcr(struct T *pT);

int start();

int CreateDatFile(
				  char cdf_sz_name[], 
				  char cdf_sz_number[], 
				  int cdf_n, 
				  int cdf_m, 
				  int cdf_me, 
				  int cdf_nmod, 
				  int cdf_nteta, 
				  int cdf_prTask, 
				  double cdf_x[], 
				  double cdf__xmin[], 
				  double cdf_xmax[], 
				  double cdf_komod[]
				  );

int MaximizeConstraint(								// возвращает ifail ропуда
					   FILE *F,						// файл отчета
					   struct T *pT,				// область, в которой максимизируется ограничение
					   int ConstraintToMaximize,	// номер ограничения для максимизации
					   double *ConstraintValue,		// значение максимизированного ораничения
					   double *tetamax,				// координаты точки максимума
					   double U,
					   double gamma,
					   int iPrintFlag				// печать на экран промежуточной информации
					   );

int MaximizeAllConstraints(
						   FILE *F,						// файл отчета
						   struct T *proot,				// указатель на корень областей
						   double **S,					// массив критических точек (размерность Nq*nkrit (приходит=NULL)
						   int *nCrit,					// суммарное число крит точек по всем областям (=0)
						   double U,					// U
						   int Nd,
						   double gamma,
						   int iPrintFlag				// печать на экран промежуточной информации
						   );

int mvadzo(											// возвращает ifail ропуда
					   FILE *F,						// файл отчета
					   struct T *proot,				// области для поиска
					   struct TK *prootK,			// области критерия
					   double **S,					// возвращаемые критические точки
					   int *nCrit,					// возвращаемое число критических точек
					   double *U,					// возвращаемое значение критерия
					   int Nd,
					   double gamma,
					   int Nap,
					   double *pQap,
					   int Na,
					   double *alpha,
					   int iPrintFlag			// печать на экран промежуточной информации
					   );

void SolveMinUdzo(									// возвращает ifail ропуда
					   FILE *F,						// файл отчета
					   struct T *proot,				// области для поиска
					   struct TK *prootK,			// области критерия
					   double *U,					// 
					   double *pQap,
					   double *alpha,
					   int iPrintFlag,			// печать на экран промежуточной информации
					   int ifail[]
					   );


void SplitT(int Nt, int Ntsplit,int *V, int *qsplit, double *array_teta,
	   double *tetamin, double *tetamax, double *z1,
	   double *Qcr, int *NQcr);



int MaximizeAllCriteria(
						   FILE *F,						// файл отчета
						   struct TK *proot,				// указатель на корень областей
						   double U,					// U
						   double **d,
						   int Nd,
						   double gamma,
						   int iPrintFlag				// печать на экран промежуточной информации
						   );


int MaximizeCriteria(								// возвращает ifail ропуда
					   FILE *F,						// файл отчета
					   struct TK *pT,				// область, в которой максимизируется ограничение
					   int ConstraintToMaximize,	// номер ограничения для максимизации
					   double *ConstraintValue,		// значение максимизированного ораничения
					   double *tetamax,				// координаты точки максимума
					   double U,
					   double gamma,
					   double **d,
					   int Nd,
					   int iPrintFlag				// печать на экран промежуточной информации
					   );

void upperBound(double x[],double *pfc,double f[],double g[],DOUBLE **dfdx,int n,int m,int me,double ko[],double komod[],double teta[]);
void maxCriteria(double x[],double *pfc,double f[],double g[],DOUBLE **dfdx,int n,int m,int me,double ko[],double komod[],double teta[]);
void maximize(double x[],double *pfc,double f[],double g[],DOUBLE **dfdx,int n,int m,int me,double ko[],double komod[],double teta[]);

double calcCriteria();
double calcConstraint(int NCons);
double calcConstraintOnZ(int NCons);

double varSpread(double bound, double spread, double leftBound, double rightBound, int isRight);


