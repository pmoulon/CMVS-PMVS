#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_drotm (const int N, double *X, const int incX, double *Y,
             const int incY, const double *P)
{
#define BASE double
#include "source_rotm.h"
#undef BASE
}
