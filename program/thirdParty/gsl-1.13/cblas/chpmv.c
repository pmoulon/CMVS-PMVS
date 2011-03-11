#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_chpmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
             const int N, const void *alpha, const void *Ap, const void *X,
             const int incX, const void *beta, void *Y, const int incY)
{
#define BASE float
#include "source_hpmv.h"
#undef BASE
}
