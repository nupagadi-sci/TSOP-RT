#include "..\headers\my.h"

double CalculateF(
				  int Nd, 
				  int Nq, 
				  int Nz, 
				  double *d, 
				  struct Q *qi, 
				  double *z,
				  double koef)
{
	double f;
	double V,F1,T2, T1,TW2,ATEP, CONV, CONV1;
	double CA0,F0,KR,ER,DH,CPR,T0,UU, TW1,FW,CPRW, QHE;

	T1 = d[0];
	TW2 = d[1];
	V = d[2];
	ATEP = d[3];
//	T1 = z[0];
//	TW2 = z[1];

	ER  =560.0;
	CA0 =32.04;
	DH  =-23260.0;
	CPR =167.4;
	CPRW=4.190;

	F0=qi->F0;
	T0=qi->T0;
	TW1=qi->TW1;
	UU=qi->UU;
	KR=qi->KR;

	CONV1= V*KR*CA0*exp(-ER/T1);
	CONV = CONV1/(F0+CONV1);
 
	T2   = 2*(-DH)*F0*CONV/ATEP/UU;
	T2   = T2-2*F0*CPR*(T1-T0)/ATEP/UU;
	T2   = T2-T1+TW2+TW1;

	QHE  = ATEP*UU*((T1-TW2)+(T2-TW1))/2;
 
	F1   = QHE/CPR/(T1-T2);
	FW   = QHE/CPRW/(TW2-TW1);  

	f = 691.2 * pow(V, 0.7) + 873.0 * pow(ATEP, 0.6) + 1.76 * FW + 7.056 * F1;
	printf("\nFi = %f",f);
	f = f * koef;
	if (f<0)
		printf("-");
	return f;
}

double CalculateIntegral5(
						  int Nap, 
						  int Nd, 
						  int Nq, 
						  int Nz, 
						  double *Zap, 
						  double *Qap, 
						  double *d)
{
	double f=0;
	int i,j;
	struct Q *pq;
	double *pz;
	double *pkoef;
	pkoef=malloc(Nap*sizeof(double));
	if (pkoef==NULL) exit(2);
	pkoef[0]=0.5;
	pkoef[1]=0.125;
	pkoef[2]=0.125;
	pkoef[3]=0.125;
	pkoef[4]=0.125;
	pq=malloc(sizeof(struct Q));
	if(pq==NULL) exit(2);
	pz=malloc(Nz*sizeof(double));
	if(pz==NULL) exit(2);
	for(i=0;i<Nap;i++)
	{
		for(j=0;j<5;j++)
		{
//проверить: правильный ли порядок точек в строке?
			pq->F0=Qap[i*5+j];
			j++;
			pq->T0=Qap[i*5+j];
			j++;
			pq->TW1=Qap[i*5+j];
			j++;
			pq->UU=Qap[i*5+j];
			j++;
			pq->KR=Qap[i*5+j];
		}
		for(j=0;j<Nz;j++)
		{
			pz[j]=Zap[i*_Nz+j];
		}
		printf("\nApproximation point #%d",i);
		printf("\n\tF0 = %10.5f\n\tT0 = %10.5f\n\tTW1= %10.5f\n\tUU = %10.5f\n\tKR = %10.5f",pq->F0,pq->T0,pq->TW1,pq->UU,pq->KR);
		printf("\n\n\tZ[0] = %10.5f\n\tZ[1] = %10.5f",pz[0],pz[1]);
		printf("\n\n\td[0] = %10.5f\n\td[1] = %10.5f",d[0],d[1]);
//		if(i==4)
//			printf(".");
		f+=CalculateF(Nd,Nq,Nz,d,pq,pz,pkoef[i]);
		printf("\n-==-");
	}
	free(pq);
	free(pz);
	free(pkoef);
	return f;
}
