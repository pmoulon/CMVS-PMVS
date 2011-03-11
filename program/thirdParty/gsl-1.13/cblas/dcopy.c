#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_dcopy (const int N, const double *X, const int incX, double *Y,
             const int incY)
{
#define BASE double
#include "source_copy_r.h"
#undef BASE
}
