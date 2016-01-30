#include "..\headers\my.h"

extern unsigned _stklen=8096;
int nfunU=0,nfungU=0,nfunM=0,nfungM=0,nfunL=0,nfungL=0;

int main()
{
	if (C_VERSION==0)
		printf("the stack length is %u \n",_stklen);
	putCan();
	start();
	exit(0);
}