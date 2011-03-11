#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_zher2 (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
             const int N, const void *alpha, const void *X, const int incX,
             const void *Y, const int incY, void *A, const int lda)
{
#define BASE double
#include "source_her2.h"
#undef BASE
}
