#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

double
cblas_dzasum (const int N, const void *X, const int incX)
{
#define BASE double
#include "source_asum_c.h"
#undef BASE
}
