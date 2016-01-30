#include "..\headers\my.h"
/*
struct T
{
//	границы и значени€ параметров области (констр. и упр. переменные)
	int Nd;
	int Nz;
	int Nq;
	int Nj;
	double *lbound;
	double *rbound;
//	критические точки, макс. значени€ неравенств в области, флаг активности
	double *Qcr;
	int NQcr;
	double *cv;
	int isActive;
	int *ActiveConstraints;
//	размер области и номер макс. ребра	
	double r;
	int qr;
	struct T *pnext;
};
*/

struct TSplit* CreateTS(double *lbound, double *rbound, int Nq, int i)
{
	int bOK;
	struct TSplit *pT;
	bOK=1;
	pT=(struct TSplit*)malloc(sizeof(struct TSplit));
	if (pT==NULL) bOK=0;
	if (bOK)	// если структура успешно создана
	{
		pT->lbound=(double*)malloc(Nq*sizeof(double));
		if(pT->lbound==NULL) bOK=0;
		pT->rbound=(double*)malloc(Nq*sizeof(double));
		if(pT->rbound==NULL) bOK=0;
		pT->iteration=i;
		pT->pnext=NULL;
	}
	if (bOK)	//	если все пол€ структуры проинициализированы
	{
		for(i=0;i<Nq;i++)
			pT->lbound[i]=lbound[i];
		for(i=0;i<Nq;i++)
			pT->rbound[i]=rbound[i];
	}
	return pT;
}

struct T* CreateT(int Nq, int Nj, double *lbound, double *rbound, double *Rabs)
{
	int i,bOK;
	struct T *pT;
	bOK=1;
	pT=(struct T*)malloc(sizeof(struct T));
	if (pT==NULL) bOK=0;
	if (bOK)	// если структура успешно создана
	{
		pT->lbound=(double*)malloc(Nq*sizeof(double));
		if(pT->lbound==NULL) bOK=0;
		pT->rbound=(double*)malloc(Nq*sizeof(double));
		if(pT->rbound==NULL) bOK=0;
	}
	if (bOK)	//	если все пол€ структуры проинициализированы
	{
		for(i=0;i<Nq;i++)
		{
			pT->lbound[i]=lbound[i];
			pT->rbound[i]=rbound[i];
		}
		//for(i=0;i<Nj;i++)
		//{
		//	//pT->cv[i]=0;
		//	pT->ActiveConstraints[i]=0;
		//}
		CalculateSize(pT,Rabs);
		pT->isLocked = 0;
		pT->NQcr=0;
		pT->Qcr=NULL;
		pT->isActive=0;
		pT->pnext=NULL;
	}
	else		//	если хот€ бы 1 из массивов не создан
		DeleteT(&pT);
	return pT;
}
struct TK* CreateTK(double *lbound, double *rbound)
{
	int i, bOK;
	struct TK *pT;
	bOK = 1;
	pT = malloc(sizeof(struct TK));
	if (!pT) bOK = 0;

	if (bOK)	//	≈сли структура успешно создана
	{
		pT->z = malloc(Nz*sizeof(double));
		if(!pT->z)			bOK = 0;
		if(!reset_b(&pT->b))	bOK = 0;
		pT->lbound = malloc(Nq*sizeof(double));
		if(!pT->lbound)		bOK = 0;
		pT->rbound = malloc(Nq*sizeof(double));
		if(!pT->rbound)		bOK = 0;

		//pT->relC = malloc(Nj*sizeof(struct T*));
		//if(!pT->relC)	bOK = 0;
		//for(i=0; i<Nj; i++)
		//{
		//	pT->relCons[i] = malloc(sizeof(struct BinTree));
		//	if(!pT->relCons[i])	bOK = 0;
		//	pT->relCons[i]->region = NULL;
		//	pT->relCons[i]->lBranch = NULL;
		//	pT->relCons[i]->rBranch = NULL;
		//}
	}


	if (bOK)	//	если все пол€ структуры проинициализированы
	{
		for(i=0; i<Nq; i++)
		{
			pT->lbound[i] = lbound[i];
			pT->rbound[i] = rbound[i];
		}

		CalculateSize2(pT);
		pT->isActive = 0;
		pT->pnext = NULL;
	}
	else		//	≈сли хот€ бы 1 из массивов не создан
		DeleteT(&pT);
	return pT;
}

void DeleteT(struct T **pT)
{
	if((*pT)->lbound)
	{	free((*pT)->lbound);	
	(*pT)->lbound = NULL;}
	if((*pT)->rbound!=NULL) 
		free((*pT)->rbound);
	if((*pT)->slbound)
	{
		free((*pT)->slbound);
		(*pT)->slbound = NULL;
	}
	if((*pT)->srbound)
	{
		free((*pT)->srbound);
		(*pT)->srbound = NULL;
	}
	if((*pT)->Qcr)
	{
		free((*pT)->Qcr);
		(*pT)->Qcr = NULL;
	}
	free(*pT);
	*pT=NULL;
}

void DeleteTCr(struct TK **pT)
{
	int i;
	if((*pT)->z)
	{
		free((*pT)->z);
		(*pT)->z = NULL;
	}
	if((*pT)->b)
	{
		for(i=0; i<Nz; i++)
		{
			if((*pT)->b[i])
			{
				free((*pT)->b[i]);
				(*pT)->b[i] = NULL;
			}
		}
		free((*pT)->b);
		(*pT)->b = NULL;
	}
	if((*pT)->lbound)
	{
		free((*pT)->lbound);
		(*pT)->lbound = NULL;
	}
	if((*pT)->rbound)
	{
		free((*pT)->rbound);
		(*pT)->lbound = NULL;
	}
	//if((*pT)->ai)
	//{
	//	free((*pT)->ai);
	//	//(*pT)->ai = NULL;
	//}

	free(*pT);
	*pT = NULL;
}

void DeleteList(struct T **proot)
{
	struct T *ptmp;
	ptmp=*proot;
	while(ptmp->pnext!=NULL)
	{
		ptmp=(*proot)->pnext;
		DeleteT(proot);
		*proot=ptmp;
	}
	DeleteT(proot);
	*proot=NULL;
}
void DeleteListCr(struct TK **proot)
{
	struct TK *ptmp;
	ptmp = *proot;
	while(ptmp->pnext)
	{
		ptmp = (*proot)->pnext;
		DeleteTCr(proot);
		*proot = ptmp;
	}
	DeleteTCr(proot);
	*proot = NULL;
}

void CalculateSize(struct T *pT, double *Rabs)
{
	double tempr, r;
	int i, nq;

	r = 0;	nq = 0;
	tempr = 0;

	for(i=0;i<Nq;i++)
	{
		tempr = (pT->rbound[i]-pT->lbound[i])/Rabs[i];

		if(tempr > r)
		{
			r = tempr;	nq = i;
		}
	}
	pT->r = r;
	pT->qr = nq;
}

void CalculateSize2(struct TK *pT)
{
	int tempr, r;
	int i, nq;

	r = 10000000000;	nq = 0;
	tempr = 0;

	for(i=0;i<Nq;i++)
	{
		tempr = round(Rabs[i]/(pT->rbound[i]-pT->lbound[i]));

		if(tempr < r)
		{
			r = tempr;	nq = i;
		}
	}
	pT->r = 1.0/(double)r;
	pT->qr = nq;
}

int CreateNewQcr(int NQcr, double *pQcr, struct T *pT)	//	ј’“”Ќ√!!111 ≈сли кр. точек 0, то добавл€ет (0.51; 0.51) 
{
	double *ptemp;
	int i,j,nq;
	nq=Nq;
	if(NQcr)
	{
		ptemp=(double*)malloc(NQcr*nq*sizeof(double));
		if(ptemp==NULL) return -1;
		for(i=0;i<NQcr;i++)
			for(j=0;j<nq;j++)
				ptemp[i*nq+j]=pQcr[i*nq+j];
		if(pT->Qcr!=NULL)
			free(pT->Qcr);
		pT->Qcr=ptemp;
	}
	//else 
	//{
	//	pT->Qcr = NULL;
	//						//	ј’“”Ќ√!!!11
	//	pT->Qcr = (double*)malloc(Nq*sizeof(double));
	//	for(i=0; i<Nq; i++)
	//		pT->Qcr[i] = 0.51;
	//	pT->NQcr = 1;
	//}
	return 0;
}

int AddQcr(struct T *pT, double *pQcr)
{
	int i, iq, NQcr;
	int flag;

	flag = 0;
	NQcr = pT->NQcr;

	pT->Qcr = (double*)realloc(pT->Qcr, (NQcr+1)*_Nq*sizeof(double));
	if(!(pT->Qcr))	return 0;

		//	ѕроверка на идентичность точек, дл€ избежани€ их дублировани€
	//for(i=0; i<NQcr; i++)
	//{
	//	for(iq=0; iq<_Nq; iq++)
	//	{
	//		if(fabs(pQcr[iq]) < EPS_QCR)	//	≈сли близко к нолю, сравниваем абсолютные значени€
	//		{
	//			if( fabs(pT->Qcr[i*_Nq+iq]-pQcr[iq]) < EPS_QCR )
	//				flag++;
	//			else iq = _Nq;	//	≈сли хоть одна коор-та точки не соответствует, остальные не провер€ем
	//		}
	//		else							//	»наче - относительные
	//			if( fabs(pT->Qcr[i*_Nq+iq]-pQcr[iq]) / pQcr[iq] < EPS_QCR )
	//				flag++;
	//			else iq = _Nq;
	//	}
	//	if(flag == _Nq)	return 0;
	//	else	flag = 0;
	//}


	for(iq=0; iq<_Nq; iq++)
		pT->Qcr[NQcr*_Nq+iq] = pQcr[iq];

	pT->NQcr++;

	return 1;
}

int DeleteQcr(struct T *pT)
{
	if(pT->Qcr!=NULL)
	{
		free(pT->Qcr);
		pT->Qcr=NULL;
	}
	return 0;
}

struct T* Split(struct T *pT1, double *Rabs)
{
	struct T *pT2;
	int indexq, indexn, iq, iz, ij, icrit, iT1, iT2;
	double *pQcrTsplit, *pQcrTsplitIf, *pQcrT1, *pQcrT2; // [NQcr][Nq]
	double tetasplit;	//	т.делени€
	int nq,T2_NQcr;
	double *pT2_lb, *pT2_rb;
	
	icrit=pT1->NQcr;
	nq=Nq;

	pT2_lb=(double*)malloc(nq*sizeof(double));
	if(pT2_lb==NULL) return NULL;
	pT2_rb=(double*)malloc(nq*sizeof(double));
	if(pT2_rb==NULL) return NULL;


	if(icrit)
	{
		pQcrTsplit=(double*)malloc(icrit*nq*sizeof(double));
		if(pQcrTsplit==NULL) return NULL;
		pQcrTsplitIf=(double*)malloc(icrit*nq*sizeof(double));
		if(pQcrTsplitIf==NULL) return NULL;
		pQcrT1=(double*)malloc(icrit*nq*sizeof(double));
		if(pQcrT1==NULL) return NULL;
		pQcrT2=(double*)malloc(icrit*nq*sizeof(double));
		if(pQcrT2==NULL) return NULL;
	}

	for(indexq=0;indexq<nq*icrit;indexq++)
	{
		pQcrTsplit[indexq]=pT1->Qcr[indexq];
		if(pT1->isSoftCons)
			pQcrTsplitIf[indexq]=pT1->slbound[indexq%Nq]+(pT1->srbound[indexq%Nq]-pT1->slbound[indexq%Nq])*pT1->Qcr[indexq];
		else
			pQcrTsplitIf[indexq]=pT1->lbound[indexq%Nq]+(pT1->rbound[indexq%Nq]-pT1->lbound[indexq%Nq])*pT1->Qcr[indexq];
	}

	indexq = pT1->qr;	//	–ебро дроблени€
	tetasplit=(pT1->rbound[indexq]+pT1->lbound[indexq])/2;
	for(iq=0;iq<nq;iq++)
	{
		pT2_lb[iq]=pT1->lbound[iq];
		pT2_rb[iq]=pT1->rbound[iq];
	}
	pT2_lb[indexq]=tetasplit;
	pT1->rbound[indexq]=tetasplit;
	iT1=iT2=0;

	pT2 = CreateT(Nq,Nj,pT2_lb,pT2_rb,Rabs);
	if(!pT2) return NULL;
	pT2->isSoftCons = pT1->isSoftCons;

	if(pT1->isSoftCons)
	{
		pT2->slbound = (double*)malloc(nq*sizeof(double));
		pT2->srbound = (double*)malloc(nq*sizeof(double));
		for (ij=0;ij<nq;ij++)
		{
			pT2->slbound[ij] = pT1->slbound[ij];
			pT2->srbound[ij] = pT1->srbound[ij];
			if (pT1->slbound[ij] > pT1->rbound[ij]) pT1->slbound[ij] = pT1->rbound[ij];
			if (pT1->srbound[ij] > pT1->rbound[ij]) pT1->srbound[ij] = pT1->rbound[ij];
			if (pT2->slbound[ij] < pT2->lbound[ij]) pT2->slbound[ij] = pT2->lbound[ij];
			if (pT2->srbound[ij] < pT2->lbound[ij]) pT2->srbound[ij] = pT2->lbound[ij];
		}
	}

//	раскидываем точки по 2 област€м

	for(ij=0;ij<icrit;ij++)
	{
		if(pQcrTsplitIf[ij*nq+indexq]<=tetasplit)
		{
			for(iq=0; iq<nq; iq++)
				pQcrT1[iT1*nq+iq] = pQcrTsplit[ij*nq+iq];
			if(pT1->isSoftCons && pT1->srbound[indexq]!=pT1->slbound[indexq])	//	„тобы не делить на ноль
				pQcrT1[iT1*nq+indexq] = (pQcrTsplitIf[ij*nq+indexq]-pT1->slbound[indexq])
										/	(pT1->srbound[indexq]-pT1->slbound[indexq]);
			else
				pQcrT1[iT1*nq+indexq] = (pQcrTsplitIf[ij*nq+indexq]-pT1->lbound[indexq])
										/	(pT1->rbound[indexq]-pT1->lbound[indexq]);
			iT1++;
		}
		else
		{
			for(iq=0;iq<nq;iq++)
				pQcrT2[iT2*nq+iq] = pQcrTsplit[ij*nq+iq];
			if(pT2->isSoftCons && pT2->srbound[indexq]!=pT2->slbound[indexq])	//	„тобы не делить на ноль
				pQcrT2[iT2*nq+indexq] = (pQcrTsplitIf[ij*nq+indexq]-pT2->slbound[indexq])
										/	(pT2->srbound[indexq]-pT2->slbound[indexq]);
			else
				pQcrT2[iT2*nq+indexq] = (pQcrTsplitIf[ij*nq+indexq]-pT2->lbound[indexq])
										/	(pT2->rbound[indexq]-pT2->lbound[indexq]);	
			iT2++;
		}
	}

	pT1->NQcr=iT1;
	T2_NQcr=iT2;

	pT2->NQcr=iT2;
	if(icrit)
	{
		CreateNewQcr(iT1,pQcrT1,pT1);
		CreateNewQcr(iT2,pQcrT2,pT2);
	}

	if(!pT1->isSoftCons)
	{
		pT1->slbound = pT1->lbound;
		pT1->srbound = pT1->rbound;
		pT2->slbound = pT2->lbound;
		pT2->srbound = pT2->rbound;
	}
//		нужно заполнить новую область своими крит точками и перераспределить 
//		точки дл€ исходной области


//дублируем номер ограничени€
	pT2->NCons = pT1->NCons;

	pT2->id = getNewID(&pT1->id);

	CalculateSize(pT1, Rabs);
	CalculateSize(pT2, Rabs);
//	чистим пам€ть
	if(icrit)
	{
		free(pQcrTsplit);
		free(pQcrTsplitIf);
		free(pQcrT1);
		free(pQcrT2);
	}
	free(pT2_lb);
	free(pT2_rb);

	return pT2;
}

struct TK* SplitK(struct TK *pT1, double *Rabs)
{
	struct TK *pT2;
	int indexq, iq, iz, i;
	double tetasplit;	//	т.делени€

	//Nsp = pT1->NumSplits;
	
	indexq = pT1->qr;	//	–ебро дроблени€
	tetasplit = (pT1->rbound[indexq] + pT1->lbound[indexq]) /2;

	pT2 = CreateTK(pT1->lbound, pT1->rbound);
	if(!pT2) return NULL;

	pT2->lbound[indexq] = tetasplit;
	pT1->rbound[indexq] = tetasplit;

	pT2->id = getNewID(&pT1->id);

	for(iq=0; iq<Nq; iq++)
	{
		pT2->iE[iq] = pT1->iE[iq];
		pT2->E[iq] = pT1->E[iq];
		pT2->ai[iq] = pT1->ai[iq];
	}
	pT2->iE[indexq] = pT1->iE[indexq]*=2;
	pT2->iE[indexq]++;


	pT1->E[indexq] = E[getSplitDepth(pT1->r/2)][pT1->iE[indexq]];
	pT2->E[indexq] = E[getSplitDepth(pT2->r/2)][pT2->iE[indexq]];
	pT1->ai[indexq] = A[getSplitDepth(pT1->r/2)][pT1->iE[indexq]];
	pT2->ai[indexq] = A[getSplitDepth(pT2->r/2)][pT2->iE[indexq]];

	pT1->a = pT2->a = A[0][0];
	for(iq=0; iq<Nq; iq++)
	{
		pT1->a *= pT1->ai[iq];
		pT2->a *= pT2->ai[iq];
	}


	if(is_b_split_inheritance)
		for(iz=0; iz<Nz; iz++)	//	Ќаследование значений b дл€ дробимых областей
			for(iq=0; iq<Nq+1; iq++) 
				pT2->b[iz][iq] = pT1->b[iz][iq];

	CalculateSize2(pT1);
	CalculateSize2(pT2);

	return pT2;
}

int isLocked(struct T *pT)
{
	int i;
	for(i=0; i<Nq; i++)
		if(fabs(pT->srbound[i] - pT->slbound[i]) < EPS_MAXFI) return 1;
	return 0;
}

int getNewID(int *lid)
{
	*lid = *lid*2 + 1;
	return *lid+1;
}

int getSplitDepth(double r)
{
	int i=0, ar;
	ar = round(1.0/r);
	while(ar > 1)
	{
		ar/=2;	i++;
	}
	return i;
}

int isRelated(int CritID, int ConsID)	//	ƒолжно выполн€тьс€ CritID <= ConsID
{
	while(CritID < ConsID)
	{
		if(ConsID%2)		ConsID--;
		else				ConsID-=2;
		ConsID/=2;
	}
	return CritID==ConsID;
}

int fork(int **f1, int **f2, int Nsp)
{
	int i;
	*f1 = realloc(*f1, Nsp*sizeof(int));
	*f2 = malloc(Nsp*sizeof(int));
	if(!((*f1)&&(*f2)))	return 0;	//	≈сли пам€ть на какой-либо массив не выделилась
	for(i=0; i<Nsp-1; i++)
		(*f2)[i] = (*f1)[i];
	(*f1)[Nsp-1] = 0;	(*f2)[Nsp-1] = 1;
}

int round(double d)
{
	int i;

	if((d - (int)d) > 0.9)
		i = (int)++d;
	else
		i = (int)d;
	return i;
}

