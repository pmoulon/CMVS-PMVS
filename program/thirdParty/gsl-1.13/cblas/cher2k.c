#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_cher2k (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
              const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
              const void *alpha, const void *A, const int lda, const void *B,
              const int ldb, const float beta, void *C, const int ldc)
{
#define BASE float
#include "source_her2k.h"
#undef BASE
}
