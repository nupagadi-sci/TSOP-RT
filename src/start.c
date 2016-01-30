#include "..\headers\my.h"

int start()
{
	struct TK *prootK, *pTCr, *pTK, *pTK2=NULL, *ptmpK=NULL;

	FILE *foutput,*foutputt, *itroutp, *log;
	int nIteration;
	int indext, indexz, indexd, indexq, iq;
	struct T *pspt, *proot=NULL, *pT=NULL, *pT2=NULL, *ptmp=NULL, *prootU=NULL, *prootU1=NULL, *prootU2=NULL, *prootU3=NULL, *prootU4=NULL,*pTU=NULL;
	struct T **ppT;
	int bContinue, bInT;
	int iStopClause;
	int nIterationUbelowZero,nIterationLaboveZero;
	int ifail, ifailC;
	double U,L,XU,XL,ul;

	double start_d[_Nd],start_mind[_Nd],start_maxd[_Nd];
	double start_u,start_minu,start_maxu;

	double *plb,*prb,*pmyd;
	double *Sl=NULL, *S=NULL;
	int nCrit=0, nTotalCrit=0;
	char ss1[20];
	double sigma;
	int iCriterion; // 0 - min U, 1 - min sum(d0-d)^2+kb

	double *pQap, *pQapL;
	double *pZap;
	double *pQapTemp;

	int splitDepth;
	FILE *fstart;
	char buf[255];
	char szname[255],szvalue[255];
	char symb;
	int i,j,k;

	int __maxIt,__nDirMax,__Task;
	double __epsF,__epsZ,__epsG,__epsFcgFg;
	int __iprint,__iprGetR,__iprRop,__iptTras;
	double __epsM,__epsT,__epsUL,__epsULabs,__rmin,__gamma;
	//double *alpha;

	//double *test[3], **p, t1=1,t2=2,t3=3;

	clock_t t_start, t_end;
	double total_time;
	int /*iPrintFlag,*/ ScreenOutput;
///************************************
//
//		стартовые условия (U,Z,Q)
//
///************************************
	foutput=NULL;
	foutputt=NULL;
	remove("log.dat");


	iCriterion=0;

	Nt=1;

	//test[0] = &t1;
	//test[1] = &t2;
	//test[2] = &t3;
	//test[3] = NULL;
	//p=test;
	//while(*p)
	//	printf("%lf", **p++);

///^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//		задание значений "по-умолчанию" для переменных
__maxIt=1000;
__nDirMax=32;
__iprint=1;
__iprGetR=1;
__iprRop=1;
__iptTras=1;
__epsF=1e-7;
__epsZ=1e-5;
__epsG=1e-7;
__epsFcgFg=1e-7;
__epsM=1e-7;
__epsT=1e-2;
__epsUL=1e-2;
__epsULabs=1e-3;
__rmin=1e-2;
__gamma=1;
__Task=0;
ScreenOutput=0;

//----------
//		чтение файла входных параметров для комплекса (файл должен быть 
//		расположен в одном каталоге с исполняемым файлом )
//	fstart=fopen("c:\\ropud\\dzowstr\\start.dat","r");

	fstart=fopen("start.dat","r");
	if(!fstart)	exit(1007);
	while(!feof(fstart))
	{
		fscanf(fstart, "%s", buf);
		//if(!strcmp(buf, "maxIt"))		fscanf(fstart, "%d", &__maxIt);
		//if(!strcmp(buf, "nDirMax"))		fscanf(fstart, "%d", &__nDirMax);
		if(!strcmp(buf, "ScreenOutput"))		fscanf(fstart, "%d", &ScreenOutput);
		if(!strcmp(buf, "iprRop"))		fscanf(fstart, "%d", &__iprRop);
		//if(!strcmp(buf, "Na"))		fscanf(fstart, "%d", &Na);
		//if(!strcmp(buf, "Nq"))		fscanf(fstart, "%d", &Nq);
		//if(!strcmp(buf, "Nj"))		fscanf(fstart, "%d", &Nj);
		//if(!strcmp(buf, "Nd"))		fscanf(fstart, "%d", &Nd);

		if(!strcmp(buf, "epsZ"))	fscanf(fstart, "%lf", &__epsZ);
		if(!strcmp(buf, "epsG"))	fscanf(fstart, "%lf", &__epsG);
		if(!strcmp(buf, "epsFcgFg"))	fscanf(fstart, "%lf", &__epsFcgFg);
		if(!strcmp(buf, "epsUL"))	fscanf(fstart, "%lf", &__epsUL);
		if(!strcmp(buf, "epsULabs"))	fscanf(fstart, "%lf", &__epsULabs);

		if(!strcmp(buf, "linerize"))	fscanf(fstart, "%d", &is_linearization );
		if(!strcmp(buf, "b_min_inheritance"))	fscanf(fstart, "%d", &is_b_min_inheritance);
		if(!strcmp(buf, "b_split_inheritance"))	fscanf(fstart, "%d", &is_b_split_inheritance);

		if(!strcmp(buf, "gamma"))	fscanf(fstart, "%lf", &gamma);

		if(!strcmp(buf, "alpha"))
		{
			//alpha = malloc(Na*sizeof(double));
			for(i=0; i<Na; i++)	fscanf(fstart, "%lf", &alpha[i]);
		}

		if(!strcmp(buf, "d"))
			for(i=0; i<Nd; i++)	fscanf(fstart, "%lf", &d[i]);
		if(!strcmp(buf, "dmin"))
			for(i=0; i<Nd; i++)	fscanf(fstart, "%lf", &dmin[i]);
		if(!strcmp(buf, "dmax"))
			for(i=0; i<Nd; i++)	fscanf(fstart, "%lf", &dmax[i]);

		if(!strcmp(buf, "zstart"))
			for(i=0; i<Nz; i++)	fscanf(fstart, "%lf", &z_start[i]);
		if(!strcmp(buf, "zmin"))
			for(i=0; i<Nz; i++)	fscanf(fstart, "%lf", &zmin[i]);
		if(!strcmp(buf, "zmax"))
			for(i=0; i<Nz; i++)	fscanf(fstart, "%lf", &zmax[i]);

		if(!strcmp(buf, "Q"))
			for(i=0; i<Nq; i++)	fscanf(fstart, "%lf", &start_q[i]);
		if(!strcmp(buf, "Qmin"))
			for(i=0; i<Nq; i++)	fscanf(fstart, "%lf", &start_minq[i]);
		if(!strcmp(buf, "Qmax"))
			for(i=0; i<Nq; i++)	fscanf(fstart, "%lf", &start_maxq[i]);

		if(!strcmp(buf, "bmin"))	fscanf(fstart, "%lf", &bmin);
		if(!strcmp(buf, "bmax"))	fscanf(fstart, "%lf", &bmax);

		if(!strcmp(buf, "splitDepth"))		fscanf(fstart, "%d", &splitDepth);
		if(!strcmp(buf, "E"))
		{
			E = malloc((splitDepth+1)*sizeof(double));
			for(i=0; i<splitDepth+1; i++)	
			{
				E[i] = malloc(pow(2,i)*sizeof(double));
				for(k=0; k<pow(2,i); k++)	
					fscanf(fstart, "%lf", &E[i][k]);
			}
		}
		if(!strcmp(buf, "A"))
		{
			A = malloc((splitDepth+1)*sizeof(double));
			for(i=0; i<splitDepth+1; i++)	
			{
				A[i] = malloc(pow(2,i)*sizeof(double));
				for(k=0; k<pow(2,i); k++)	
					fscanf(fstart, "%lf", &A[i][k]);
			}
		}
		if(!strcmp(buf, "P"))
		{
			P = malloc(Nq*sizeof(double));
			for(i=0; i<Nq; i++)	fscanf(fstart, "%lf", &P[i]);
		}

		fgets(buf, 255, fstart);
	}


//		вывод на экран стартовых значений задачи
	printf("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
	printf("\nmaxIt =    %d",__maxIt);
	printf("\nnDirMax =  %d",__nDirMax);
	printf("\niprint =   %d",__iprint);
	printf("\niprGetR =  %d",__iprGetR);
	printf("\niprRop =   %d",__iprRop);
	printf("\niptTras =  %d",__iptTras);
	printf("\nepsF =     %.10f",__epsF);
	printf("\nepsZ =     %.10f",__epsZ);
	printf("\nepsG =     %.10f",__epsG);
	printf("\nepsFcgFg = %.10f",__epsFcgFg);
	printf("\nepsM =     %.10f",__epsM);
	printf("\nepsT =     %.10f",__epsT);
	printf("\nepsUL =    %.10f",__epsUL);
	printf("\nepsULabs = %.10f",__epsULabs);
	//printf("\nNj =       %d",__Nj);
	//printf("\nNq =       %d",__Nq);
	//printf("\nNd =       %d",__Nd);
	printf("\nNa =       %d",Na);
	//printf("\nNt =       %d",__Nt);
	//printf("\ngamma =    %.10f",__gamma);

	for(k=0;k<Na;k++)
		printf("\nalpha[%d] =    %.10f",k,alpha[k]); 

	if(ScreenOutput)
		printf("\nPrint working process to screen");
	else
		printf("\nPrint results only");
	printf("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
	fclose(fstart);

	reset_b(&b_start);


//		определение номинального значения и размера области
//		неопределенности по каждому неопределенному параметру
//		в зависимости от заданного gamma

//		Q
	for(iq=0; iq<Nq; iq++)
	{
		//start_q[iq] = E[0][0];
		start_minq[iq] = start_q[iq] - (start_q[iq]-start_minq[iq])*gamma;
		start_maxq[iq] = start_q[iq] + (start_maxq[iq]-start_q[iq])*gamma;
		Rabs[iq] = start_maxq[iq] - start_minq[iq];
	}

	Nap = 1;
	Ns = Na;
	pQap = malloc(Nap*Nq*sizeof(double));
	if(!pQap)	exit(1012);
	for(iq=0; iq<Nq; iq++)
		pQap[iq] = start_q[iq];

		//	Создание корня областей критерия
	prootK = malloc(sizeof(struct TK));
	prootK->z = NULL;
	prootK->b = NULL;
	prootK->lbound = NULL;
	prootK->rbound = NULL;
	pTK = prootK;

		//	Инициализация области критерия
	pTK->pnext = CreateTK(start_minq, start_maxq);
	pTK = pTK->pnext;
	//pTK->NumSplits = 0;
	pTK->id = 0;
	for(i=0; i<Nq; i++)
	{
		pTK->iE[i] = 0;
		pTK->E[i] = E[0][0];
	}
	pTK->a = 1;
	for(i=0; i<Nq; i++)
		pTK->a *= pTK->ai[i] = A[0][0];
	pTK->pnext = NULL;

	Nt = Nj;


		//	Создание корня областей для верхней оценки
	prootU=malloc(sizeof(struct T));
	prootU->lbound=NULL;
	prootU->rbound=NULL;
	prootU->Qcr=NULL;
	prootU->relK = NULL;
	pTU = prootU;


	for(i=0; i<Ns; i++)		//	Мягкие ограничения
	{
		pTU->pnext = CreateT(Nq, Nj, start_minq, start_maxq, Rabs);
		pTU = pTU->pnext;
		pTU->id = 0;
		pTU->relK = pTK;
		pTU->NCons = i;
		pTU->isSoftCons = 1;
		pTU->slbound = malloc(Nq*sizeof(double));
		pTU->srbound = malloc(Nq*sizeof(double));
		for(indexq=0; indexq<Nq; indexq++)
		{
			pTU->slbound[indexq] = start_minq[indexq] + 0*Rabs[indexq];
			pTU->srbound[indexq] = start_maxq[indexq] - 0*Rabs[indexq];
		}
	}

	for(; i<Nj; i++)		//	Жесткие основные
	{
		pTU->pnext = CreateT(Nq, Nj, start_minq, start_maxq, Rabs);
		pTU = pTU->pnext;
		pTU->id = 0;
		pTU->relK = pTK;
		pTU->NCons = i;
		pTU->isSoftCons = 0;
		pTU->slbound = pTU->lbound;
		pTU->srbound = pTU->rbound;
	}	

	//ppT = pTK->relC;	//	test
	//while(*ppT++)
	//	;

	//for(; i<Nj+Nz*2; i++)		//	Жесткие - границы на z
	//{
	//	pTU->pnext = CreateT(Nq, Nj, start_minq, start_maxq, Rabs);
	//	pTU = pTU->pnext;
	//	pTU->NCons = i;
	//	pTU->isSoftCons = 0;
	//}

///************************************
	if (C_VERSION==0) strcpy(ss1,"6340304.28073996");     /*Borland*/
	if (C_VERSION==1) strcpy(ss1,"6340776.11626546");     /*Visual*/
//****************************
	if(ScreenOutput)
	{
		foutput=fopen("output.txt","wt");
		if(foutput==NULL) 
		{
			printf("\nERROR: can\'t create output.txt");
			exit(1014);
		}
		fprintf(foutput,"\tGetting starting critical points:");
		foutputt=fopen("outputt.txt","wt");
		if(foutput==NULL) 
		{
			printf("\nERROR: can\'t create output.txt");
			exit(1015);
		}
	}
	iStopClause=0;


	XU = U = 0;
	bContinue=TRUE;
	nIteration=0;
	nIterationUbelowZero=0;
	nIterationLaboveZero=0;
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//-------------------------------
	nCrit=1;
	Sl=malloc(nCrit*Nq*sizeof(double));
	if(Sl==NULL) exit(1016);

	for(indexq=0;indexq<Nq;indexq++)
	{
		Sl[indexq]=0.51;
	}

	//for(pT=prootU->pnext; pT; pT=pT->pnext)
	//	AddQcr(pT, Sl);

	free(Sl);
	nCrit = Nt;

	Sl=malloc(nCrit*Nq*sizeof(double));
	if(Sl==NULL) exit(1017);
	for(indexq=0;indexq<Nq*nCrit;indexq++)
	{
		Sl[indexq]=0.51;
	}


	bContinue=1;
	iStopClause=0;
	U=0;			//нужно ли это?

//end of initialization
//start	
	
	t_start=clock();

	while(bContinue)
	{
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//		присвоение управляющим переменным для точек аппроксимации
//		тех значений, которые найдены на этапе создания первоначального
//		разбиения в соответствии с попаданием точек аппроксимации 
//		в ту или иную подобласть

//создание границ для "вероятностных" областей


	if(ScreenOutput)
	{
		printf("\n\n\nIteration #%d",nIteration+1);
		printf("\n<<<CALCULATING UPPER BOUND>>>\n");
	}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//		поиск верхней оценки критерия F (МВА)

	ifail=mvadzo(foutput,
				 prootU,
				 prootK,
				 &Sl,
				 &nCrit,
				 &U,
				 Nd,
				 gamma,
				 Nap,
				 pQap,
				 Na,
				 alpha,
				 ScreenOutput);

	XU=U;

		if((U<0)&&(nIterationUbelowZero==0)) nIterationUbelowZero=nIteration;
//		if(ifail>1) exit(1018);
//		if(ifail!=0) exit(1019);

		//getch();

		XL=L=0;
//====================================================================================

		if(ScreenOutput)
		{
			printf("\n\n\n<<<CALCULATING LOWER BOUND>>>\n");
		}

	XL=L;
//	XL=XU;



			t_end=clock();
			total_time=((double)(t_end - t_start)/CLOCKS_PER_SEC);

			itroutp = fopen("itroutp.dat", "at");
			fprintf(itroutp, "\nItrtn %d - fc = %lf\tV1 = %lf\tV2 = %lf\tTotal time = %lf\tifail = %d\n", nIteration, XU, d[0], d[1], total_time, ifail);
			pTK = prootK;
			while(pTK=pTK->pnext)
			{
				fprintf(itroutp, "\nBounds\tleft\tright\n");
				for(iq=0; iq<Nq; iq++)
					fprintf(itroutp, "Q%d\t%lf\t%lf\n", iq, pTK->lbound[iq], pTK->rbound[iq]);
			}
			pTK = prootK;
			while(pTK=pTK->pnext)
			{
				fprintf(itroutp, "\nB\tZ0\tZ1\n");
				for(iq=0; iq<Nq+1; iq++)
					fprintf(itroutp, "Q%d\t%lf\t%lf\n", iq, pTK->b[0][iq], pTK->b[1][iq]);
			}
			fclose(itroutp);


				//	Cr P
			log = fopen("log.dat", "at");
			fprintf(log, "\tBEFORE SPLIT\n");
			pT = prootU;
			while(pT=pT->pnext)
			{
				fprintf(log, "id%d ", pT->id);
				if(pT->isSoftCons)	fprintf(log, "(soft)");
				else	fprintf(log, "(hard)");
				fprintf(log, " j=%d\n", pT->NCons);
				fprintf(log, "bij:\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", pT->relK->b[0][0], pT->relK->b[0][1], pT->relK->b[0][2], pT->relK->b[0][3], pT->relK->b[0][4], pT->relK->b[0][5]);
				fprintf(log, "bij:\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", pT->relK->b[1][0], pT->relK->b[1][1], pT->relK->b[1][2], pT->relK->b[1][3], pT->relK->b[1][4], pT->relK->b[1][5]);
				fprintf(log, "Bounds:\n");
				for(i=0; i<Nq; i++)
					fprintf(log, "\t%lf\t%lf\n", pT->slbound[i], pT->srbound[i]);
				fprintf(log, "Crit points:\t(%d)\n",pT->NQcr);
				for(i=0; i<pT->NQcr*Nq; i++)
				{
					fprintf(log, "\t%.20lf", pT->slbound[i%Nq] + (pT->srbound[i%Nq]-pT->slbound[i%Nq])*(pT->Qcr[i]));
					if(i%Nq)	fprintf(log, "\n");
				}
			}
			fclose(log);

				ifailC=MaximizeAllCriteria(foutputt,prootK,U,&d,Nd,gamma,ScreenOutput);


			{



//if (nIteration!=1)


	//if(nIteration==2)
	//{
	//	prootK->pnext->isActive = 0;
	//	prootK->pnext->pnext->isActive=1;
	//	prootK->pnext->pnext->pnext->isActive = 0;
	//}


					//	Если не закоментировано, то дробление по всем подобластям критерия
				pTK = prootK;	while(pTK=pTK->pnext)	pTK->isActive = 1;

				//pT = prootU;
				//while(pT=pT->pnext)
				//	pT->isActive = 0;

				pTK = prootK;	pspt = prootU;
				while(pTK=pTK->pnext)
				{
					if(pTK->isActive && getSplitDepth(pTK->r)==splitDepth)	{pTK->isActive = 0;	bContinue = 0;}
					if(pTK->isActive)
					{
						pTK2 = SplitK(pTK, Rabs);	Nap++;
						ptmpK = prootK;
						while((ptmpK=ptmpK->pnext)->pnext);
						ptmpK->pnext = pTK2;

						pT = prootU;
						while(pT=pT->pnext)
							if(isRelated((pTK->id-1)/2, pT->id))
							{
								pT->isActive = 0;
								if((pTK->id-1)/2==pT->id)
								{
									pT2 = Split(pT, Rabs);	Nt++;
									if(pT->isSoftCons)	Ns++;
										//	Устанавливаем соответствия областей критерия и ограничений
									pT->relK = pTK;
									pT2->relK = pTK2;
									ptmp = prootU;	while((ptmp=ptmp->pnext)->pnext);
									ptmp->pnext = pT2;
								}
								else if(isRelated(pTK->id, pT->id))		pT->relK = pTK;
								else if(isRelated(pTK2->id, pT->id))	pT->relK = pTK2;
								else
								{
Nt=Nt;
								}
							}
					}
				}


				pT = prootU;
				while(pT=pT->pnext)
					if(pT->isActive)
					{
						pT->isActive = 0;
						pT2 = Split(pT, Rabs);	Nt++;
						if(pT->isSoftCons)	Ns++;
						pT2->relK = pT->relK;
						ptmp = prootU;	while((ptmp=ptmp->pnext)->pnext);
						ptmp->pnext = pT2;
					};


				free(pQap);
				pQap = malloc(Nap*Nq*sizeof(double));				
				pTK=prootK;	i=0;
				while(pTK=pTK->pnext)
				{
					for(j=0; j<Nq; j++)
						pQap[i*Nq+j] = (pTK->lbound[j] + pTK->rbound[j]) /2;
					i++;
				}
			
			}



		if((L>0)&&(nIterationLaboveZero==0)) nIterationLaboveZero=nIteration;

		ul=fabs((XU-XL)/XU);


			if(ScreenOutput)
			{
				printf("\n-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-");
				printf("\nIteration #%d",nIteration+1);
				printf("\nRegions: %d",Nt);
				printf("\n           XU = %.20f\n           XL = %.20f",XU,XL);
				printf("\n            U = %.20f\n            L = %.20f",U,L);
				printf("\n   |XU-XL/XU| = %.20f\n          eps = %f",ul,__epsUL);
				printf("\n        XU-XL = %.20f\n          eps = %f",XU-XL,__epsULabs);
			
				t_end=clock();
				total_time=((double)(t_end - t_start)/CLOCKS_PER_SEC);
				printf("\n time elapsed = %.4f sec",total_time);
			}




			log = fopen("log.dat", "at");
			fprintf(log, "\tAFTER SPLIT\n");
			pT = prootU;
			while(pT=pT->pnext)
			{
				fprintf(log, "id%d ", pT->id);
				if(pT->isSoftCons)	fprintf(log, "(soft)");
				else	fprintf(log, "(hard)");
				fprintf(log, " j=%d\n", pT->NCons);
				fprintf(log, "bij(0,1,2):\t%lf\t%lf\t%lf\n", pT->relK->b[0][0], pT->relK->b[0][1], pT->relK->b[0][2]);
				fprintf(log, "Bounds:\n");
				for(i=0; i<Nq; i++)
					fprintf(log, "\t%lf\t%lf\n", pT->slbound[i], pT->srbound[i]);
				fprintf(log, "Crit points:\t(%d)\n",pT->NQcr);
				for(i=0; i<pT->NQcr*Nq; i++)
				{
					fprintf(log, "\t%.20lf", pT->slbound[i%Nq] + (pT->srbound[i%Nq]-pT->slbound[i%Nq])*(pT->Qcr[i]));
					if(i%Nq)	fprintf(log, "\n");
				}
			}
			fclose(log);

				//	log
			log = fopen("log.dat", "at");
			fprintf(log, "\t%.20lf\n", XU);
			fclose(log);


//			if(ul<EPS_UL) 

 		//	if((ul<__epsUL)||((U-L)<__epsULabs))
			//{
			//	bContinue=0;
			//	if(iStopClause==0) iStopClause=1;
			//}

			//if(ScreenOutput)
			//{
			//	printf("\n-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-");
			//	if(U<L)
			//	{
			//		printf("\n!!! U<L  !!!");
			//		//getch();
			//	}
			//}


		nIteration++;	Nit++;
		//bContinue=TRUE;
//		getch();
//		if(nIteration==1) 
//		if(nIteration==2) 
//		if(nIteration==3) 
//		if(nIteration==4) 
//		if(nIteration==5) 
//		if(nIteration==6) 
		//if(nIteration==10) 
		//{
		//	bContinue=FALSE;
		//	if(iStopClause==0) iStopClause=3;
		//}
	}//while(bContinue)
	t_end=clock();
	total_time=((double)(t_end - t_start)/CLOCKS_PER_SEC);

//		печать финальных результатов решения задачи
	printf("\n\n\n-=Summary=-");
	printf("\n%d iteration(s) was/were made",nIteration);
	printf("\nTotal Time Elapsed = %.4f sec",total_time);
	printf("\n");
	printf("\ngamma=%.5f",gamma);
	printf("\nalpha=%.5f %.5f",alpha[0],alpha[1]);
	printf("\n");
	printf("\nU");
	printf("\nF=%.5f",U);
	printf("\nNt=%d, Nap=%d",Nt,Nap);
	//for(pT=prootU->pnext; pT; pT=pT->pnext)
	//	printf("\nd=%.5f %.5f",d[0],pT->z[0]);
	printf("\n");
	printf("\nL");
	printf("\nF=%.5f",L);
	printf("\n");
	printf("\nnfunM=%d nfungM=%d",nfunM,nfungM);
	printf("\nnfunU=%d nfungU=%d",nfunU,nfungU);
	printf("\nnfunL=%d nfungL=%d",nfunL,nfungL);


//	if(nIterationUbelowZero>0)
//		printf("\nU<0 at %d iteration",nIterationUbelowZero);
//	if(nIterationLaboveZero>0)
//		printf("\nL>0 at %d iteration",nIterationLaboveZero);
/*	switch(iStopClause)
	{
		case 1:
			printf("\n\nDifference between Upper and Lower bounds is satisfactory.\n");
			break;
		case 2:
			printf("\n\nThere are no active regions or their sizes are too small.\n");
			break;
		case 3:
			printf("\n\nMax allowed iteration number is reached.\n");
			break;
		case 4:
			printf("\n\nStructure is unflexible.\n");
			break;
		default:
			printf("\n\nUnknown stop reason.\n");
	}
*/
//	getch();

/*	if(__Task!=2)
		for(i=0;i<Nap;i++)
			for(indexz=0;indexz<2;indexz++)
				printf("\nZap[%d][%d]= %.10f",i,indexz,pZap[i*2+indexz]);
*/
//		освобождение памяти
	//DeleteList(&proot);
	DeleteList(&prootU);
	DeleteListCr(&prootK);
	//DeleteList(&prootU);	//	Границы областей ограничений верхней оценки ссылаются на критерий
	put_skip(0);

	free(S);
	free(pQap);
	free(Sl);
	//free(alpha);
	//exit(1022);
	_getch();
}
