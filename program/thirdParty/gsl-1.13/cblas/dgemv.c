#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_dgemv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA,
             const int M, const int N, const double alpha, const double *A,
             const int lda, const double *X, const int incX,
             const double beta, double *Y, const int incY)
{
#define BASE double
#include "source_gemv_r.h"
#undef BASE
}
