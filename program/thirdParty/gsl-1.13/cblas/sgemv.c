#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_sgemv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA,
             const int M, const int N, const float alpha, const float *A,
             const int lda, const float *X, const int incX, const float beta,
             float *Y, const int incY)
{
#define BASE float
#include "source_gemv_r.h"
#undef BASE
}
