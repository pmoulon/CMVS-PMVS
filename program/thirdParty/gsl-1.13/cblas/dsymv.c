#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_dsymv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
             const int N, const double alpha, const double *A, const int lda,
             const double *X, const int incX, const double beta, double *Y,
             const int incY)
{
#define BASE double
#include "source_symv.h"
#undef BASE
}
