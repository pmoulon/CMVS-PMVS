#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

double
cblas_dznrm2 (const int N, const void *X, const int incX)
{
#define BASE double
#include "source_nrm2_c.h"
#undef BASE
}
