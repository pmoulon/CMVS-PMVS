/* rbd 14-Dec-99 for Win32 */
#include "f2c.h"
#if !defined _WIN32
#include <sys/times.h>
#endif
#include <sys/types.h>
#include <time.h>

#ifndef CLK_TCK
#define CLK_TCK 60
#endif

doublereal second_()
{
#if defined _WIN32
  clock_t rusage;

  rusage = clock();
  return (doublereal)(rusage) / CLOCKS_PER_SEC;
#else
  struct tms rusage;

  times(&rusage);
  return (doublereal)(rusage.tms_utime) / CLK_TCK;
#endif
} /* second_ */

