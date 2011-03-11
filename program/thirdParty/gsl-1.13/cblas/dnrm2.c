#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

double
cblas_dnrm2 (const int N, const double *X, const int incX)
{
#define BASE double
#include "source_nrm2_r.h"
#undef BASE
}
