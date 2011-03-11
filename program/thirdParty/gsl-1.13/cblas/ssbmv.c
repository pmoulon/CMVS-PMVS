#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_ssbmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
             const int N, const int K, const float alpha, const float *A,
             const int lda, const float *X, const int incX, const float beta,
             float *Y, const int incY)
{
#define BASE float
#include "source_sbmv.h"
#undef BASE
}
