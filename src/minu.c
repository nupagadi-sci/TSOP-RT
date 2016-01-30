#include "..\headers\my.h"

void SolveMinUdzo(									// возвращает ifail ропуда
					   FILE *F,						// файл отчета
					   struct T *proot,				// области для поиска
					   struct TK *prootK,			// области критерия
					   double *U,					// 
					   double *pQap,
					   double *alpha,
					   int iPrintFlag,			// печать на экран промежуточной информации
					   int ifail[]
					   )
{
	struct T *pT, *pT2, *pTs;
	struct TK *pTK;
	double *pmyx,*pmyxmin,*pmyxmax,*pmykomod,sig[_Nq],kk,FF;
	int iz, iq, it, id, indext, indexd, indexq, indexa, N;
	int i_komod, ij, total_krit, iCrit;
	int Nkomod,	//	Размер массива komod
		Nx,		//	Количество поисковых
		Nf;		//	Количество ограничений в пользователькой функции

		//	Переменные ROPUD
	static double *x,fc,*f,*xmin,*xmax,*multX,*multF,*lam,*ko,*komod,*teta;
	static int *mzamp; static int n,m,me;
	int can,parm;
	int parmget,nAl,mAl;

	FILE *itroutp;
	double spread;
	int i,j;

	FILE *log;
	
		//	----------------
	ROPUD_pT = proot;	//	Передача указателя обрабатываемой области ограничения
	ROPUD_pTK = prootK;	//	и области критерия
		//	----------------

	Nx = Nd +
			2*Ns*Nq +	//	плавающие левая и правая граница областей
				Nap*Nz*(Nq+1);	//	b; Nap соответствует числу областей критерия
	pmyx=(double*)malloc(Nx * sizeof(double));
	if(pmyx==NULL) exit(2);
	pmyxmin=(double*)malloc(Nx * sizeof(double));
	if(pmyxmin==NULL) exit(2);
	pmyxmax=(double*)malloc(Nx * sizeof(double));
	if(pmyxmax==NULL) exit(2);
	

		//	Инициализация total_krit
	total_krit = 0;
	pT = proot;	while(pT=pT->pnext)	total_krit+=pT->NQcr;

	ij = 0;
	for(id=0; id<Nd; id++,ij++)
	{
		pmyx[ij] = d[id];9;
		pmyxmin[ij] = dmin[id];
		pmyxmax[ij] = dmax[id];
	}

			
	//spread = 1;
	//switch(Nit)
	//{
	//case 0:
	//	//pmyx[0] = 5.7999749488811432*1.3;
	//	//pmyx[1] = 7.2757090972328990*1.3;
	//	pmyx[0] = 5.674554*1.0;
	//	pmyx[1] = 7.247084*1.0;
	//	//pmyx[0] = 5.714615*1.05;
	//	//pmyx[1] = 8.343656*1.05;
	//	break;
	//case 1:
	//case 2:
	//	pmyx[0] = 5.674554*1.0;
	//	pmyx[1] = 7.247084*1.0;
	//	//pmyx[0] = 5.714615*1.05;
	//	//pmyx[1] = 8.343656*1.05;
	//	spread = 0.2;
	//	break;
	//case 3:
	//case 4:
	//	pmyx[0] = 5.674554*0.8;
	//	pmyx[1] = 7.247084*0.8;
	//	//pmyx[0] = 5.714615*1.05;
	//	//pmyx[1] = 8.343656*1.05;
	//	spread = 0.02;
	//	break;
	//case 5:
	//case 6:
	//	pmyx[0] = 5.674554*0.8;
	//	pmyx[1] = 7.247084*0.8;
	//	//pmyx[0] = 5.714615*1.05;
	//	//pmyx[1] = 8.343656*1.05;
	//	spread = 0.15;
	//	break;
	//case 7:
	//case 8:
	//	pmyx[0] = 5.674554*0.8;
	//	pmyx[1] = 7.247084*0.8;
	//	//pmyx[0] = 5.714615*1.05;
	//	//pmyx[1] = 8.343656*1.05;
	//	spread = 0.1;
	//	break;
	//}
	//if(Nit>0)	
	//{
	//	spread *= 1;
	//	pmyx[0] = d[0];	pmyx[1] = d[1];
	//	pmyxmin[0] = varSpread(d[0], spread*1.5, dmin[0], dmax[0], 0);
	//	pmyxmax[0] = varSpread(d[0], spread*0.5, dmin[0], dmax[0], 1);
	//	pmyxmin[1] = varSpread(d[1], spread*1.5, dmin[1], dmax[1], 0);
	//	pmyxmax[1] = varSpread(d[1], spread*0.5, dmin[1], dmax[1], 1);
	//}

		//itroutp = fopen("itroutp.dat", "at");
		//fprintf(itroutp, "\nStart d\n");
		//for(id=0; id<Nd; id++)
		//	fprintf(itroutp, "\td[%d]= %lf \tdmin[%d]= %lf \tdmax[%d]= %lf \n", id, pmyx[id], id, pmyxmin[id], id, pmyxmax[id]);
		//fclose(itroutp);


	spread = 0.05;
	pT = proot;
	while(pT=pT->pnext)
		if(pT->isSoftCons)
		{
			for(iq=0; iq<Nq; iq++,ij++)
			{
				pmyx[ij] = pT->lbound[iq]+(pT->rbound[iq]-pT->lbound[iq])*0.25;pT->slbound[iq];
				pmyx[ij+Nq] = pT->rbound[iq]-(pT->rbound[iq]-pT->lbound[iq])*0.25;pT->slbound[iq];
				pmyxmin[ij] = pT->lbound[iq];
				pmyxmax[ij] = pT->rbound[iq];

				pmyxmin[ij+Nq] = pT->lbound[iq];
				pmyxmax[ij+Nq] = pT->rbound[iq];

					//	Добавить другие старты для 5й мягкой, т.к. альфа у нее отличается от других

				switch(Nit)
				{
				case 1:
				case 2:
					spread = 0.03;
					break;
				case 3:
				case 4:
					spread = 0.05;
					break;
				case 5:
				case 6:
					spread = 0.1;
					break;
				case 7:
				case 8:
					spread = 0.01;
					break;
				}
				if(Nit>0)	
				{
					spread *= 2;
					pmyx[ij] = pT->slbound[iq]; pmyx[ij+Nq] = pT->srbound[iq];
					pmyxmin[ij] = varSpread(pT->slbound[iq], spread, pT->lbound[iq], pT->rbound[iq], 0);
					pmyxmax[ij] = varSpread(pT->slbound[iq], spread, pT->lbound[iq], pT->rbound[iq], 1);

					pmyxmin[ij+Nq] = varSpread(pT->srbound[iq], spread, pT->lbound[iq], pT->rbound[iq], 0);
					pmyxmax[ij+Nq] = varSpread(pT->srbound[iq], spread, pT->lbound[iq], pT->rbound[iq], 1);
				}
			}

			if(pT->isLocked)	//	Схлопнулась => нет кр. точек => снова разъедется
			{
				ij-=Nq;
				for(iq=0; iq<Nq; iq++,ij++)
				{
					pmyx[ij] =	start_q[iq];
					pmyxmin[ij] = start_q[iq];
					pmyxmax[ij] = start_q[iq];

					pmyx[ij+Nq] = start_q[iq];
					pmyxmin[ij+Nq] = start_q[iq];
					pmyxmax[ij+Nq] = start_q[iq];
				}
			}

			ij+=Nq;
		};


		//	В силу специфики хранения областей для ограничений, b берем раз в Nj областей
	pTK = prootK;
	while(pTK=pTK->pnext)
	{		
		if(!is_b_min_inheritance && Nit == 0)
			reset_b(&pTK->b);
		switch(Nit)
		{
		case 0:
		case 1:
		case 2:
			spread = 0.07;
			break;
		case 3:
		case 4:
			spread = 0.05;
			break;
		case 5:
		case 6:
			spread = 0.03;
			break;
		case 7:
		case 8:
			spread = 0.1;
			break;
		}
		for(iz=0; iz<Nz; iz++)			//	Берем b
		{
			for(iq=0; iq<Nq+1; iq++, ij++)
			{
				pmyx[ij]	= pTK->b[iz][iq];0;
				pmyxmin[ij]	= bmin;0;
				pmyxmax[ij]	= bmax;0;
				if(Nit>0)	
				{
					spread *= 1;
					pmyx[ij] = pTK->b[iz][iq];
					pmyxmin[ij] = varSpread(pTK->b[iz][iq], spread, bmin, bmax, 0);
					pmyxmax[ij] = varSpread(pTK->b[iz][iq], spread, bmin, bmax, 1);
				}
			}
		}



	}

	Nkomod = 0;7+Nt+total_krit*Nq+Nap*Nq+Na+Nt+Nap+Nap*Nq+Nap*Nq;
	//pmykomod = (double*)malloc(Nkomod*sizeof(double));
	//if(!pmykomod)	exit(2);
//в общем случае
		//	m Число - ограничений
	Nf = 1+Na			//	Ограничения "Эф круглое";
			+Ns*Nq			//	Правая граница больше левой
				+total_krit 	//	Критические точки (основное ограничение)
					+Nap*Nz*2;//		//	Ограничения на z в апроксимационных точках
						//+Nap		//	Положительность критерия (ограничение отсутствует)
	//pmykomod[1]=0;				// me
	//pmykomod[2]=Nd;
	//pmykomod[3]=Nt;
	//pmykomod[4]=Nap;
	//pmykomod[5]=gamma;
	//pmykomod[6]=Na;
	//i_komod=7;

//	pT = proot;	while(pT=pT->pnext)	pmykomod[i_komod++]=pT->NQcr;
//
//	pT = proot;
//	while(pT=pT->pnext)
//	{
//		iCrit = pT->NQcr;
//		for(indexq=0;indexq<Nq*iCrit;indexq++)
//			pmykomod[i_komod++]=pT->Qcr[indexq];
//	}
//
//	for(indexq=0;indexq<Nq*Nap;indexq++)
//		pmykomod[i_komod++]=pQap[indexq];
//	
////вероятности
//	for(indexa=0;indexa<Na;indexa++)
//		pmykomod[i_komod++] = alpha[indexa];
//
////номера ограничений, для которых действуют области
//	pT = proot;	while(pT=pT->pnext)	pmykomod[i_komod++]=pT->NCons;
//
////коэффициенты областей критерия
//	pTK = prootK;	while(pTK=pTK->pnext)	pmykomod[i_komod++]=pTK->a;
//
//
//	pTK = prootK;
//	while(pTK=pTK->pnext)
//		for(indext=0;indext<Nq;indext++)
//			pmykomod[i_komod++]=pTK->ai[indext];


	CreateDatFile("test.dat","XY", Nx, Nf, 0,Nkomod,0,3,pmyx,pmyxmin,pmyxmax,NULL);

	//free(pmykomod);

	can=1;
	parm=2; 
	parmget=1;
	getname("","test","dat","rez",can,parm);
	getmem(&x,&xmin,&xmax,&f,&multX,&multF,&lam,&mzamp,&ko,&komod,&teta,can,&nAl,&mAl);
	callropud1(can,parm,parmget,&x,&fc,&f,&xmin,&xmax,&multX,&multF,&lam,&mzamp,&ko,&komod,&teta,&n,&m,&me);

	nfunU+=(int)ko[31];
	nfungU+=(int)ko[36];

	if(iPrintFlag)
	{
		fprintf(F,"\n\tSolving minU problem:");
		fprintf(F,"\nNumber of function calculations is %.0f, number of direction calculations is %.0f.",ko[31],ko[36]);
	}
	
//	*U=x[0];
	*U=fc;

		// log
	//log = fopen("log.dat", "at");
	//fprintf(log, "%.20lf\n", fc);
	//fclose(log);
		
//	ij=1;
	ij=0;

	for(indexd=0; indexd<Nd; indexd++)
		d[indexd] = x[ij++];

	pT = proot;
	while(pT=pT->pnext)
		if(pT->isSoftCons)
		{
			for(indexq=0;indexq<Nq;indexq++)
				pT->slbound[indexq] = x[ij++];

			for(indexq=0;indexq<Nq;indexq++)
				pT->srbound[indexq] = x[ij++];
		};

		//	Копирование полученных поисковых b
	for(pTK=prootK->pnext; pTK; pTK=pTK->pnext)
		for(iz=0; iz<Nz; iz++)
			for(iq=0; iq<Nq+1; iq++, ij++)
				pTK->b[iz][iq] = x[ij];


	ifail[0]=(int)ko[75];
	ifail[1]=(int)ko[171];

//	printf("\nU=%.20f",x[0]);
//	if(iPrintFlag)
//	{
//
//		printf("\n");
//		printf("\n======================================================");
//		printf("\n");

		itroutp = fopen("itroutp.dat", "at");
		fprintf(itroutp, "ifail=%d ifailF=%d\n", ifail[0], ifail[1]);
		fclose(itroutp);

		printf("\nifail=%d ifailF=%d",ifail[0],ifail[1]);
//		printf("\nnfunM=%d nfungM=%d",nfunM,nfungM);
//		printf("\nnfunU=%d nfungU=%d",nfunU,nfungU);
//		printf("\n");
//		printf("\nF=%f",*U);
//		for(indexd=0;indexd<Nd;indexd++)
//			printf("\ndmin[%1d]=%f d[%1d]=%f dmax[%1d]=%f",indexd,dmin[indexd],indexd,(d)[indexd],indexd,dmax[indexd]);
//		for(iz=0; iz<Nz; iz++)
//			printf("\nzmin[%1d]=%f zmax[%1d]=%f",indexd,zmin[iz],indexd,zmax[iz]);
//
//		for(pT=proot->pnext,indext=0; pT; pT=pT->pnext)
//		{
//			printf("\n[cons %1d]\n\tb[0][0]=%f b[0][1]=%f b[0][2]=%f",indext++,pT->relK->b[0][0],pT->relK->b[0][1],pT->relK->b[0][2]);
//			for(iq=0; iq<pT->NQcr; iq++)
//			{
////				set_z(pT, pT->Qcr+iq*Nq);
//				printf("\n\t");
//				for(iz=0; iz<Nq; iz++)
//					printf("Qcr%d[%d]=%f ",iq,iz,pT->Qcr[iq*Nq+iz]);
////				printf(" -");
////				for(iz=0; iz<Nz; iz++)
////					printf(" z[%d]=%f",0,pT->z[iz]);
//			}
//		}
//
//
//		printf("\n");
//
//		//printf("\nApproximation points");
//		//for(indexa=0;indexa<Nap;indexa++)
//		//{
//		//	printf("\n#%1d weight=%f",indexa+1,a[indexa]);
//		//	printf("\n%f  %f  ",pQap[indexa*Nq+0],pQap[indexa*Nq+1]);
//		//}
//
//		printf("\n    left         sleft         sright       right");
//		printf("\n======================================================");
//
//		pT = proot;
//		indext = 1;
//		while(pT=pT->pnext)
//			if(pT->isSoftCons)
//			{
//				printf("\n  Region #%1d",indext);
//				printf("\n constraint number = %1d",pT->NCons);
//				//printf("\n alpha = %f",pT->alpha);
//				printf("\n%f     %f     %f     %f",pT->lbound[0],pT->slbound[0],pT->srbound[0],pT->rbound[0]);
//				printf("\n%f    %f    %f    %f",pT->lbound[1],pT->slbound[1],pT->srbound[1],pT->rbound[1]);
//				printf("\n");
//				indext++;
//			}
//
//		printf("\n======================================================");
//		
//	}

//	if((ifail==1)&&(ifailF==0))
//	{
//#ifdef DEBUG
//		printf("\n\n!!! ifail=1 & ifailF=0, small change in search variables");
//#endif
//		ifail=1;
//	}
//	if((ifail==1)&&(ifailF==1))
//	{
//#ifdef DEBUG
//		printf("\n\n!!! ifail=1 & ifailF=1, small change in search variables");
//#endif
//		ifail=1;
//	}
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	free(pmyx);
	free(pmyxmin);
	free(pmyxmax);
	put_end(can);
}



