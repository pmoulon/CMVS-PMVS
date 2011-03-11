#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

float
cblas_snrm2 (const int N, const float *X, const int incX)
{
#define BASE float
#include "source_nrm2_r.h"
#undef BASE
}
