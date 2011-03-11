#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_dgemm (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
             const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
             const int K, const double alpha, const double *A, const int lda,
             const double *B, const int ldb, const double beta, double *C,
             const int ldc)
{
#define BASE double
#include "source_gemm_r.h"
#undef BASE
}
