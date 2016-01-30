#include "..\headers\my.h"


int mvadzo(											// возвращает ifail ропуда
					   FILE *F,						// файл отчета
					   struct T *proot,				// области для ограничения
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
			)
{
	int nIteration, nTotalCrit;
	int bContinue, bFirstTime;
	int nRegion;
	int indexq, indexj;
	double *Sl=NULL;
	struct T *pT, *pTnext;
	int ifail[2],ifailM;
	int isScndMinNeed;

	FILE *fmva;
	int i;

	if(iPrintFlag)
	{
		fmva=fopen("mva.txt","wt");
		if(fmva==NULL) 
		{
			printf("\nERROR: can\'t create mva.txt");
			exit(2);
		}
	}

	nIteration=0;
	bContinue=TRUE;
	bFirstTime=TRUE;
	nTotalCrit=0;
	nTotalCrit=*nCrit;		//	?????
//	nCrit=0;
	while(bContinue)// - МВА цикл
	{
		if(Sl!=NULL) 
		{
			free(Sl);
			Sl=NULL;
		}
	
		ifailM=MaximizeAllConstraints(F,proot,&Sl,nCrit,0,Nd,gamma,iPrintFlag);
		if(ifailM>1)	
			exit(32);

	//if(iPrintFlag)
	//{
		pT=proot;
		nRegion=0;
//		printf("\n======================================================");
		//while(pT->pnext)
		//{
			pT=pT->pnext;
			printf("\nIteration#%3d\tRegion #%3d\tnCrit=%3d",nIteration,nRegion,*nCrit);
			printf("\n\tU=%12.6f\n\teps=%10.6f",*U,EPS_MAXFI);
			//indexj=pT->NCons;
			//printf("\n cv[%2d]=%10.6f\tac[%2d]=%d",indexj,pT->cv[indexj],indexj,pT->ActiveConstraints[indexj]);
			//nRegion++;
	//	}
	//}

		//if(*S==NULL)
		//	*S=(double*)malloc((*nCrit)*Nq*sizeof(double));
		//else
		//	*S=(double*)realloc(*S,(nTotalCrit+(*nCrit))*Nq*sizeof(double));
		//if((*S==NULL)&&(nTotalCrit>0)&&(nCrit>0)) exit(11);
		//for(indexq=0;indexq<(*nCrit)*Nq;indexq++)
		//	(*S)[nTotalCrit*Nq+indexq]=Sl[indexq];
		nTotalCrit+=(*nCrit);
		if(bFirstTime||(*nCrit>0)||(ifail[1]>0))
		{
			//calcDD();
			SolveMinUdzo(F,proot,prootK,U,pQap,alpha,iPrintFlag, ifail);
 			nIteration++;
			bFirstTime=FALSE;
			printf("\n\t%lf\n", *U);

			if(!ifail[0])
			{	
				pT = proot;
				isScndMinNeed = 0;
				while(pT->pnext)
				{
					if(pT->pnext->isSoftCons && isLocked(pT->pnext))
						{
							isScndMinNeed = 1;
							pT->pnext->isLocked = 1;
							pTnext = pT->pnext->pnext;
							free(pT->pnext->Qcr);
							pT->pnext->Qcr = NULL;
							pT->pnext->NQcr = 0;
							//DeleteT(&pT->pnext);
							//Ns--; Nt--;
							//pT->pnext = pTnext;
						};
					pT=pT->pnext;
				}
				
				if(isScndMinNeed)
				{
					SolveMinUdzo(F,proot,prootK,U,pQap,alpha,iPrintFlag, ifail);
 					nIteration++;
					printf("\t%lf", *U);
				}
			}

		}
		else
			bContinue = FALSE;

		if(nIteration > 30) bContinue = FALSE;


		if(iPrintFlag)
		{
			fprintf(fmva,"\nIteration #%d\nV= %10.5f\tATEP= %10.5f",nIteration,(d)[0],(d)[1]);
			pT=proot;
			i=0;
			for(i=0;i<Nap;i++)
			{
				fprintf(fmva,"\nApprox #%d\tF0= %10.5f\tT0= %10.5f",i,(pQap)[0+i*Nap],(pQap)[1+i*Nap]);
			}
			fprintf(fmva,"\nf= %10.5f\n\n",*U);
		}

	}//	while(bContinue) - МВА цикл




// возвращаем все точки, найденные МВА
	*nCrit=nTotalCrit;
	free(Sl);
	if(iPrintFlag)
		fclose(F);

	return ifail[0];
}
