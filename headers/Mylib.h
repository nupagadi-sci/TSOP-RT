#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <conio.h>
//#include <alloc.h>
#include <process.h>
#include <math.h>
#include <time.h>
#define _CHAR     -1
#define _INT      -2
#define _LONG     -4
#define _DOUBLE   -8
#define STOP      exit(0)
#define BEGIN     {
#define END       }

void alloc1(int type,...);
void alloc1L(int type,...);
void alloc2(int type,...);
void alloc2s(int type,...);
void corrAll2(DOUBLE **arr,int n,int m);
void dealloc(void *point1,...);
char *d_to_s(double value,int ndig);
//void exe_fort(char name[],...);
void   get_list(int nBuf,...);
void   get_data(int nBuf,...);
void   get_datat(int nBuf,...);
void   get_start(int nBuf, char *);
void   get_end(int nBuf);
void   put_end(int nFile);
void   put_list(int nFile,...);
void   put_listL(int nFile,...);
void   put_lists(int nFile,...);
void   put_skip(int n);
void   put_start(int nFile, char *);
void   put_start1(int nFile, char *);
