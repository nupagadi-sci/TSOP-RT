#include "..\headers\my.h"

//static double V1, V2;	//	Конструктивные
//static double RT1, RT2;	//	Управляющие
//static double E1, E2, K10, K20;	//	Неопределенные
//static double K1, K2, K3, K4, CA1, CA2, CB1, CB2;	//	Прочие
//static double C = 0.5;	//	Параметр

static double ER=5.6, CA0=32.04, DH=-232.6, CPR =1.674, CPRW=4.190,
				V,F1,T2,T1,TW2,ATEP,CONV,CONV1,F0,KR,T0,UU,TW1,FW,QHE;



void calcModel(double *d, double *z, double *Q)
{
	F0 = Q[0];
	T0 = Q[1];
	TW1 = Q[2];
	UU = Q[3];
	KR = Q[4];

	V = d[0];
	ATEP = d[1];

	T1 = z[0];
	TW2 = z[1];

		//	Модель
	CONV1= V*KR*CA0*exp(-ER/T1);
	CONV = CONV1/(10*F0+CONV1);
 	T2   = 2*(-DH)*F0*CONV/ATEP/UU/100;
	T2   = T2-2*F0*CPR*(T1-T0)/ATEP/UU;
	T2   = T2-T1+TW2+TW1;
	QHE  = ATEP*UU*((T1-TW2)+(T2-TW1))/2;
	F1   = 10*QHE/CPR/(T1-T2);
	FW   = 1000*QHE/CPRW/(TW2-TW1);
}

double* calcDerivative(double *d, double *z, double **b, double *Q)
{
	double dF_dQ[_Nq];

	return dF_dQ;
}

double calcConstraint(int NCons)
{
	double f;
	switch(NCons)
	{
			//	Мягкие
	case 0:
		f=(0.9-CONV);
		break;
	case 1:
		f=(T2-3.89);
		break;
	case 2:
		f=(TW1-TW2+0.01);
		break;
	case 3:
		f=(TW1- T2+0.111);
		break;
	case 4:
		f=(3.11-T2);
		break;
	case 5:
		f=((T2-T1)+0.01);	//	0.95
		break;

			//	Жесткие
	case 6:
		f=(TW2-T1+0.111);
		break;
	case 7:
		f=(zmin[0]-T1);
		break;
	case 8:
		f=(T1-zmax[0]);
		break;
	case 9:
		f=(zmin[1]-TW2);
		break;
	case 10:
		f=(TW2-zmax[1]);
		break;
	//case 11:
	//	f=(-F1);
	//	break;
	//case 12:
	//	f=(-FW);
	//	break;
	}
	return f;
}

double calcConstraintOnZ(int NCons)
{
	double f = 0;
	switch(NCons)
	{
	case 0:
		f=(zmin[0]-T1);
		break;
	case 1:
		f=(T1-zmax[0]);
		break;
	case 2:
		f=(zmin[1]-TW2);
		break;
	case 3:
		f=(TW2-zmax[1]);
		break;	}
	return f;
}

double calcCriteria()
{
	return	691.2 * pow(V, 0.7) + 873.0 * pow(ATEP, 0.6) + 1.76 * FW + 7.056 * F1;
}
