#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

float
cblas_scnrm2 (const int N, const void *X, const int incX)
{
#define BASE float
#include "source_nrm2_c.h"
#undef BASE
}
