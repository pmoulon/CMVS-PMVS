#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_dscal (const int N, const double alpha, double *X, const int incX)
{
#define BASE double
#include "source_scal_r.h"
#undef BASE
}
