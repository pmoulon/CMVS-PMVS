#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_zdscal (const int N, const double alpha, void *X, const int incX)
{
#define BASE double
#include "source_scal_c_s.h"
#undef BASE
}
