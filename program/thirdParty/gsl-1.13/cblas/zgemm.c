#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_zgemm (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
             const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
             const int K, const void *alpha, const void *A, const int lda,
             const void *B, const int ldb, const void *beta, void *C,
             const int ldc)
{
#define BASE double
#include "source_gemm_c.h"
#undef BASE
}
