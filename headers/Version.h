#if C_VERSION == 0
#define DOUBLE double huge
#define INT int huge
#define HUGE huge
#define CALLOC farcalloc
#define FREE farfree
#define CHAR char huge
#include <alloc.h>
#endif

#if C_VERSION == 1
#define DOUBLE double
#define INT int 
#define HUGE   
#define CALLOC calloc
#define FREE free
#define CHAR char
#include <malloc.h>
#endif
