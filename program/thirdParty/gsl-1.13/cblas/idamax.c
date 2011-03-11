#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

CBLAS_INDEX
cblas_idamax (const int N, const double *X, const int incX)
{
#define BASE double
#include "source_iamax_r.h"
#undef BASE
}
