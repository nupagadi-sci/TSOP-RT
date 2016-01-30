#include "..\headers\my.h"
#define NFILE 101  /* the number of files */
static FILE *pFile[NFILE];

int CreateDatFile(char cdf_sz_name[], 
				  char cdf_sz_number[], 
				  int cdf_n, 
				  int cdf_m, 
				  int cdf_me, 
				  int cdf_nmod, 
				  int cdf_nteta, 
				  int cdf_prTask, 
				  double cdf_x[], 
				  double cdf_xmin[], 
				  double cdf_xmax[], 
				  double cdf_komod[])
/*int CreateDatFile(int flag,char cdf_sz_name[], char cdf_sz_number[], int cdf_n, 
				  int cdf_m, int cdf_me, int cdf_nmod, int cdf_nteta, 
				  int cdf_prTask, double cdf_x[], double cdf_xmin[], 
				  double cdf_xmax[], double cdf_komod[])*/
/*
*	Функция создает файл задания.
*	Возвращает 0 в случае успеха.
*/
{
	int cdf_i;
	FILE *F;
	char cdf_sz_fname[40];
	strcpy(cdf_sz_fname,cdf_sz_name);
	F=fopen(cdf_sz_fname,"w+t");
	if (F==NULL)
		return(1);
	else
	{
//		fprintf(F,"nAl=%d mAl=%d nmodAl=%d ntetaAl=%d\n;",cdf_n,cdf_m,cdf_nmod,cdf_nteta);
		fprintf(F,"nAl=%d mAl=%d nmodAl=%d\n;",cdf_n,cdf_m,cdf_nmod);
//		fprintf(F,"\nn=%d m=%d me=%d nmod=%d nteta=%d; \nprTask=%d\n",cdf_n,cdf_m,cdf_me,cdf_nmod,cdf_nteta,cdf_prTask);
		fprintf(F,"\nn=%d m=%d me=%d nmod=%d; \nprTask=%d\n",cdf_n,cdf_m,cdf_me,cdf_nmod,cdf_prTask);
		fprintf(F,"\nmaxIt=1000 zamp=0 prDerM=0 nDirMax=256");
//		fprintf(F,"\nprMet3=1 prMet2=0 prMet1=0");
		//fprintf(F,"\nepsF=1E-6 epsZ=1E-4 epsG=1E-6 epsFcgFg=1E-6");
		//fprintf(F,"\nepsF=1E-5 epsZ=1E-7 epsG=1E-5 epsFcgFg=1E-7");
		//fprintf(F,"\nepsF=1E-5 epsZ=1E-5 epsG=1E-5 epsFcgFg=1E-5");
		//fprintf(F,"\nepsF=1E-7 epsZ=1E-7 epsG=1E-5 epsFcgFg=1E-7");
		//fprintf(F,"\nepsF=1E-4 epsZ=6E-4 epsG=1E-4 epsFcgFg=1E-4");
//		fprintf(F,"\nepsF=1E-7 epsZ=1E-4 epsG=1E-7 epsFcgFg=1E-7");	//	!!!
		fprintf(F,"\nepsF=1E-7 epsZ=1E-7 epsG=1E-7 epsFcgFg=1E-7");
//		fprintf(F,"\niprint=1\niprGetR=1 iprRop=1 iptTras=1");
		fprintf(F,"\niprint=3\niprGetR=2 iprRop=2 iptTras=2");
//		fprintf(F,"\niprint=4\niprGetR=1 iprRop=1 iptTras=2");
//		fprintf(F,"\niprint=4\niprGetR=2 iprRop=2 iptTras=2");
		fprintf(F,"\nnumber=\"%s\"\n;",cdf_sz_number);
		fprintf(F,"\nx0=");
		for(cdf_i=0;cdf_i<cdf_n;cdf_i++)
			fprintf(F," %.16f",cdf_x[cdf_i]);
		fprintf(F,"\nxmin=");
		for(cdf_i=0;cdf_i<cdf_n;cdf_i++)
			fprintf(F," %.16f",cdf_xmin[cdf_i]);
		fprintf(F,"\nxmax=");
		for(cdf_i=0;cdf_i<cdf_n;cdf_i++)
			fprintf(F," %.16f",cdf_xmax[cdf_i]);
//		if(flag==0)
		{
			if(cdf_nmod>0)
			{
				fprintf(F,"\nkomod=");
				for(cdf_i=0;cdf_i<cdf_nmod;cdf_i++)
				{
//!!!!!!!!!!!
//					if (cdf_i>=(4+_Nd+_Nt+cdf_komod[3]))
//						cdf_komod[cdf_i]=(double)(cdf_i+1);
					fprintf(F," %.16f",cdf_komod[cdf_i]);
				}
			}
		}
/*		else
		{
			if(cdf_nmod>0)
			{
				fprintf(F,"\nkomod=");
				for(cdf_i=0;cdf_i<cdf_nmod;cdf_i++)
				{
				}
			}
		}
*/		fprintf(F,"\n;");
		fclose(F);

/*		
		printf("nAl=%d mAl=%d nmodAl=%d ntetaAl=%d\n;",cdf_n,cdf_m,cdf_nmod,cdf_nteta);
		printf("\nn=%d m=%d me=%d nmod=%d nteta=%d prTask=%d\n;",cdf_n,cdf_m,cdf_me,cdf_nmod,cdf_nteta,cdf_prTask);
		printf("\nmaxIt=100 zamp=0");
		printf("\niprint=4\niprGetR=1 iprRop=1 iptTras=2");
		printf("\nnumber=\"%s\"\n;",cdf_sz_number);
		printf("\nx0=");
		for(cdf_i=0;cdf_i<cdf_n;cdf_i++)
			printf(" %f",cdf_x[cdf_i]);
		printf("\nxmin=");
		for(cdf_i=0;cdf_i<cdf_n;cdf_i++)
			printf(" %f",cdf_xmin[cdf_i]);
		printf("\nxmax=");
		for(cdf_i=0;cdf_i<cdf_n;cdf_i++)
			printf(" %f",cdf_xmax[cdf_i]);
		printf("\nnmod=%d\n;",cdf_nmod);
		if(cdf_nmod>0)
		{
			printf("\nkomod=");
			for(cdf_i=0;cdf_i<cdf_nmod;cdf_i++)
				printf(" %f",cdf_komod[cdf_i]);
		}
*/	//	fprintf(F,"\n;");

	}/* fopen */
	return(0);
}