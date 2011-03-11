#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_dsyr2 (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
             const int N, const double alpha, const double *X, const int incX,
             const double *Y, const int incY, double *A, const int lda)
{
#define BASE double
#include "source_syr2.h"
#undef BASE
}
