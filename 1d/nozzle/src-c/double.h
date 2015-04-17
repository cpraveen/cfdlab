//#ifndef MYDOUBLE_H_
//#define MYDOUBLE_H_ 1

#if defined DIF0

#define REALQ double
#define REALX double
#define REALA double

#elif defined DIFQ

#define REALQ adouble
#define REALX double
#define REALA adouble

#include "adolc.h"

#elif defined DIFX

#define REALQ double
#define REALX adouble
#define REALA adouble

#include "adolc.h"

#endif
//#endif
