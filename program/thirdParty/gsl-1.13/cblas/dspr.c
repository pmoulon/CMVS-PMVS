#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_dspr (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
            const int N, const double alpha, const double *X, const int incX,
            double *Ap)
{
#define BASE double
#include "source_spr.h"
#undef BASE
}
