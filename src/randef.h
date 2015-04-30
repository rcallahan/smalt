#ifdef __cplusplus
extern "C"
{
#endif

#ifndef RANDEF_H
#define RANDEF_H

#include <stdlib.h>
#include <time.h>

#ifdef _WIN32 /* Windows */

#define RANSEED(s) if ((s) <= 0) {srand(time(0));} else {srand((s));}
#define RANDRAW_UNIFORM_1() (((double) rand())/(RAND_MAX + 1))

#else /* POSIX */

#define RANSEED(s) if ((s) <= 0) {srand48(time(0));} else {srand48((s));}
#define RANDRAW_UNIFORM_1() drand48()

#endif

#endif
#ifdef __cplusplus
}
#endif
