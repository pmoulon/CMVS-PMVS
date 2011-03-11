#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_dswap (const int N, double *X, const int incX, double *Y,
             const int incY)
{
#define BASE double
#include "source_swap_r.h"
#undef BASE
}
