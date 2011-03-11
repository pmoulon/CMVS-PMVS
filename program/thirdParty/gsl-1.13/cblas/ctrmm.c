#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_ctrmm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
             const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
             const enum CBLAS_DIAG Diag, const int M, const int N,
             const void *alpha, const void *A, const int lda, void *B,
             const int ldb)
{
#define BASE float
#include "source_trmm_c.h"
#undef BASE
}
