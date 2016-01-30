#include "..\headers\my.h"
/******************************************************************************/
void functr(double x[],double *fc,double f[],double g[],DOUBLE **dfdx,
	    int n,int m,int me,double ko[],double komod[],double teta[])
/* functr for EXAMPLE1 */
{
char number[15];
int prTask;
extern long nModel; int nLevel,level,prDerpa;
nLevel=(int)ko[39]; level=(int)ko[127]; prDerpa=(int)ko[38];
strcpy(number,(char*)(ko+25));
prTask=(int)ko[28];
if (strcmp(number,"max")==0)  maximize(x,fc,f,g,dfdx,n,m,me,ko,komod,teta);
if (strcmp(number,"XY")==0)    upperBound(x,fc,f,g,dfdx,n,m,me,ko,komod,teta);
if (strcmp(number,"KU")==0)    maxCriteria(x,fc,f,g,dfdx,n,m,me,ko,komod,teta);
//if (strcmp(number,"KL")==0)    k23(x,fc,f,g,dfdx,n,m,me,ko,komod,teta);
} /* functr */
/****************************************************************************/
