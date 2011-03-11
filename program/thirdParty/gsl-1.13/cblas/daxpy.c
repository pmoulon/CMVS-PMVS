#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_daxpy (const int N, const double alpha, const double *X, const int incX,
             double *Y, const int incY)
{
#define BASE double
#include "source_axpy_r.h"
#undef BASE
}
