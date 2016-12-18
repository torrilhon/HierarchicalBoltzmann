/*************************************************************************\
 defs.h  - header with definitions to distinguish PARDISO and MKL linear
           solver routines and index integer types.
\*************************************************************************/
#ifdef PARDISO
#define INT int
#endif

#ifdef MKL_ILP64
#include <mkl.h>
#define INT MKL_INT
#endif
