#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_drot (const int N, double *X, const int incX, double *Y, const int incY,
            const double c, const double s)
{
#define BASE double
#include "source_rot.h"
#undef BASE
}
