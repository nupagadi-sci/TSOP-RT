#define TRUE  1
#define FALSE 0

//void invert(DOUBLE a[],int *n,double *d,int l[],int m[]);
void invert(DOUBLE a[],int n,double *det,int l[],int m[]);
//void  qp(long *m, long *me, long *n, DOUBLE *c, double d[], DOUBLE *a,
//   double b[], double xl[],double xu[], double x[], double u[],
//   long *ifail);
void  qp(long m, long me, long mPen, long n, long l, DOUBLE **C,
   double c[], DOUBLE **A,double b[], double xl[],double xu[],
   double x[], double lam[],long *ifail, double ko[]);
void  qpm1(long m, long me, long n, long l, DOUBLE **C, double c[],
  DOUBLE **A,double b[], double xl[],double xu[], double x[], double lam[],
  double ko[],long *ifail);
void  qpm2(long m, long me, long n, long l, DOUBLE **C, double c[],
  DOUBLE **A,double b[], double xl[],double xu[], double x[], double lam[],
  double ko[],long *ifail);
int ARITHIF(double number);
int IARITHIF(int number);
double *ADR(double x,double y);
double powi(double x,int n);
int DOCNT(int i1, int i2, int dlti);
/*double RC_INI(struct{long numb;double ini;} arr[]);*/